#!/usr/bin/env python3
"""Full cylindrical 3D Fourier/MAC evolution diagnostic.

Finite radial no-slip cylinder and periodic z; NOT a whole-space certificate.
The weighted pressure projection is compatible with divergence. Viscosity is
a nonnegative gradient quadratic form, treated by constrained CN solves.
Nonlinearity uses a weighted-adjoint gather/Lamb/scatter construction.
All retained angular modes evolve; angular and axial products are dealiased.
Axis regularity, reconstruction, domain and residual errors remain unbounded.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.fft import fft,ifft,fftfreq,next_fast_len,set_workers
from scipy.sparse import diags,bmat,csr_matrix
from scipy.sparse.linalg import splu
from scipy.interpolate import CubicSpline
from scipy.signal import resample
from initial_full_field import fields


class MAC:
    def __init__(self,nr=192,nz=512,R=8.,L=16.,mapping=4.,J=4,nu=.001,dealias='padding'):
        self.nr,self.nz,self.R,self.L,self.J,self.nu=nr,nz,R,L,J,nu
        self.dealias=dealias;self.nzpad=next_fast_len((3*nz+1)//2)
        self.rf=R*np.sinh(mapping*np.arange(nr+1)/nr)/np.sinh(mapping)
        self.r=(self.rf[:-1]+self.rf[1:])/2
        self.V=(self.rf[1:]**2-self.rf[:-1]**2)/2
        self.dr=np.diff(self.rf);self.dc=np.diff(self.r)
        self.W=self.rf[1:-1]*self.dc
        self.mass=np.concatenate((self.W,self.V,self.V))
        self.nv=3*nr-1;self.z=(np.arange(nz)-nz//2)*L/nz
        self.ki=np.rint(fftfreq(nz)*nz).astype(int);self.kz=2*np.pi*self.ki/L
        self.keep=abs(self.ki)<(nz/3 if dealias=='filter' else nz/2)
        self.modes=4*np.arange(J+1);self.nphi=next_fast_len(3*J+1)
        self.slr=slice(0,nr-1);self.slt=slice(nr-1,2*nr-1);self.slz=slice(2*nr-1,3*nr-1)
        # Cell velocity gather and radial derivative of face velocity.
        rows=[];cols=[];ivals=[];dvals=[]
        for i in range(nr):
            if i>0:rows.append(i);cols.append(i-1);ivals.append(.5);dvals.append(-1/self.dr[i])
            if i<nr-1:rows.append(i);cols.append(i);ivals.append(.5);dvals.append(1/self.dr[i])
        self.I=csr_matrix((ivals,(rows,cols)),shape=(nr,nr-1))
        self.Df=csr_matrix((dvals,(rows,cols)),shape=(nr,nr-1))
        # Interior face gradients plus the outer zero-value wall gradient.
        rr=[];cc=[];vv=[]
        for i in range(nr-1):
            rr.extend((i,i));cc.extend((i,i+1));vv.extend((-1/self.dc[i],1/self.dc[i]))
        rr.append(nr-1);cc.append(nr-1);vv.append(-1/(R-self.r[-1]))
        self.E=csr_matrix((vv,(rr,cc)),shape=(nr,nr))
        self.Wg=np.r_[self.W,R*(R-self.r[-1])]
        self.Dc_even=(self.I@self.E[:-1]).tolil()
        self.Dc_even[-1,-1]-=.5/(R-self.r[-1])
        self.Dc_even=self.Dc_even.tocsr()
        self.Dc_odd=self.Dc_even.tolil()
        self.Dc_odd[0,0]+=.5/self.r[0]
        self.Dc_odd=self.Dc_odd.tocsr()
        # These rows include the radial geometric flux factors exactly.
        self.Divr=self._divr()
        self.Grad=diags(1/self.dc)@csr_matrix((np.tile([-1.,1.],nr-1),
            (np.repeat(np.arange(nr-1),2),np.column_stack((np.arange(nr-1),np.arange(1,nr))).ravel())),shape=(nr-1,nr))
        self.K=[];self.poisson_cache=[];self.stokes_cache={}
        V=diags(self.V);Mrad=diags(self.W)
        radial=self.E.conj().T@diags(self.Wg)@self.E
        for m in self.modes:
            H=diags(self.V/self.r**2)
            Krr=self.Df.T@V@self.Df+(m*m+1)*self.I.T@H@self.I
            Ktt=radial+(m*m+1)*H;Kzz=radial+m*m*H
            Krt=2j*m*self.I.T@H
            K=bmat([[Krr,Krt,None],[Krt.conj().T,Ktt,None],[None,None,Kzz]],format='csc')
            self.K.append(K)
            # Factor -DG; boundary radial face degrees of freedom are removed.
            lower=-self.rf[1:-1]/(self.V[1:]*self.dc)
            upper=-self.rf[1:-1]/(self.V[:-1]*self.dc)
            diagonal=m*m/self.r**2
            diagonal=diagonal.copy();diagonal[1:]-=lower;diagonal[:-1]-=upper
            cp=np.empty((nr,nz));inv=np.empty((nr,nz))
            den=diagonal[0]+self.kz**2;inv[0]=1/den;cp[0]=upper[0]*inv[0]
            for i in range(1,nr):
                den=diagonal[i]+self.kz**2-lower[i-1]*cp[i-1]
                if i==nr-1 and m==0:den[0]=1.
                inv[i]=1/den;cp[i]=upper[i]*inv[i] if i<nr-1 else 0
            self.poisson_cache.append((lower,cp,inv))
        xx=self.r[:4]**2;scale=xx[-1]
        self.axis_weights=np.linalg.solve(np.array([(xx/scale)**j for j in range(4)]),[1.,0,0,0])

    def _divr(self):
        rows=[];cols=[];vals=[]
        for i in range(self.nr):
            if i>0:rows.append(i);cols.append(i-1);vals.append(-self.rf[i]/self.V[i])
            if i<self.nr-1:rows.append(i);cols.append(i);vals.append(self.rf[i+1]/self.V[i])
        return csr_matrix((vals,(rows,cols)),shape=(self.nr,self.nr-1))

    def div_hat(self,a,j):
        return self.Divr@a[self.slr]+1j*self.modes[j]*a[self.slt]/self.r[:,None]+1j*self.kz*a[self.slz]

    def grad_hat(self,p,j):
        return np.concatenate((self.Grad@p,1j*self.modes[j]*p/self.r[:,None],1j*self.kz*p),axis=0)

    def project_hat(self,a,j):
        rhs=self.div_hat(a,j);lower,cp,inv=self.poisson_cache[j]
        q=np.empty_like(rhs);q[0]=rhs[0]*inv[0]
        for i in range(1,self.nr):
            val=rhs[i]-lower[i-1]*q[i-1]
            if i==self.nr-1 and j==0:val[0]=0
            q[i]=val*inv[i]
        for i in range(self.nr-2,-1,-1):q[i]-=cp[i]*q[i+1]
        result=a+self.grad_hat(q,j)
        if j==0:result[self.slr,0]=0
        result[:,~self.keep]=0
        return result

    def project(self,U):
        out=np.empty_like(U)
        for j in range(self.J+1):
            out[j]=ifft(self.project_hat(fft(U[j],axis=-1),j),axis=-1)
        out[0]=out[0].real
        return out

    def gather(self,U):
        return np.stack((self.I@U[self.slr],U[self.slt],U[self.slz]))

    def scatter(self,F):
        return np.concatenate(((self.I.T@(self.V[:,None]*F[0]))/self.W[:,None],F[1],F[2]))

    def curl_cells(self,a,j):
        m=self.modes[j];u=self.gather(a)
        kz=2*np.pi*fftfreq(a.shape[-1],d=self.L/a.shape[-1])
        dz=ifft(1j*kz*fft(u,axis=-1),axis=-1)
        return np.stack((1j*m*u[2]/self.r[:,None]-dz[1],
            dz[0]-self.Dc_even@u[2],
            self.Dc_odd@u[1]+u[1]/self.r[:,None]-1j*m*u[0]/self.r[:,None]))

    def nonlinear(self,U,project=True):
        state=resample(U,self.nzpad,axis=-1) if self.dealias=='padding' else U
        nwork=state.shape[-1]
        uc=np.stack([self.gather(a) for a in state]);wc=np.stack([self.curl_cells(a,j) for j,a in enumerate(state)])
        F=np.empty_like(state)
        cell_force=np.zeros_like(uc)
        for i in range(0,self.nr,24):
            sl=slice(i,min(i+24,self.nr));shape=(3,sl.stop-sl.start,nwork,self.nphi)
            us=np.zeros(shape,complex);ws=np.zeros(shape,complex)
            for j in range(self.J+1):
                us[...,j]=uc[j,:,sl];ws[...,j]=wc[j,:,sl]
                if j:us[...,-j]=uc[j,:,sl].conj();ws[...,-j]=wc[j,:,sl].conj()
            up=ifft(us,axis=-1,norm='forward').real;wp=ifft(ws,axis=-1,norm='forward').real
            cross=np.stack((up[1]*wp[2]-up[2]*wp[1],up[2]*wp[0]-up[0]*wp[2],up[0]*wp[1]-up[1]*wp[0]))
            coeff=fft(cross,axis=-1,norm='forward')
            for j in range(self.J+1):cell_force[j,:,sl]=coeff[...,j]
        for j in range(self.J+1):F[j]=self.scatter(cell_force[j])
        if self.dealias=='padding':F=resample(F,self.nz,axis=-1)
        return self.project(F) if project else F

    def _stokes_factor(self,j,kindex,alpha):
        key=(j,kindex,alpha)
        if key in self.stokes_cache:return self.stokes_cache[key]
        n=self.nr;m=self.modes[j];k=2*np.pi*kindex/self.L
        D=bmat([[self.Divr,diags(1j*m/self.r),1j*k*diags(np.ones(n))]],format='csc')
        BM=-D.conj().T@diags(self.V)
        A=diags(self.mass*(1+alpha*k*k))+alpha*self.K[j]
        full=bmat([[A,BM],[BM.conj().T,None]],format='csc')
        perm=[]
        for i in range(n):
            if i>0:perm.append(i-1)
            perm.extend((n-1+i,2*n-1+i))
            if not(j==0 and kindex==0 and i==n-1):perm.append(self.nv+i)
        perm=np.array(perm)
        lu=splu(full[perm][:,perm],permc_spec='NATURAL')
        value=(lu,perm)
        self.stokes_cache[key]=value
        return value

    def diffuse(self,U,h):
        # Constrained Crank--Nicolson for P(-M^-1 K), alpha=h*nu/2.
        alpha=h*self.nu/2;out=np.zeros_like(U)
        for j in range(self.J+1):
            ah=fft(U[j],axis=-1)
            rhs=self.mass[:,None]*ah-alpha*(self.K[j]@ah+self.mass[:,None]*self.kz**2*ah)
            vh=np.zeros_like(ah)
            for col,kindex in enumerate(self.ki):
                if not self.keep[col]:continue
                lu,perm=self._stokes_factor(j,abs(int(kindex)),alpha)
                b=np.zeros(self.nv+self.nr,complex);b[:self.nv]=rhs[:,col]
                if kindex<0:b[self.slz]*=-1
                solution=np.zeros_like(b);solution[perm]=lu.solve(b[perm])
                if kindex<0:solution[self.slz]*=-1
                vh[:,col]=solution[:self.nv]
            out[j]=ifft(vh,axis=-1)
        out[0]=out[0].real
        return out

    def step(self,U,dt):
        U=self.diffuse(U,dt/2)
        f1=self.nonlinear(U);f2=self.nonlinear(U+dt*f1/2)
        f3=self.nonlinear(U+dt*f2/2);f4=self.nonlinear(U+dt*f3)
        U=U+dt*(f1+2*f2+2*f3+f4)/6
        return self.diffuse(U,dt/2)

    def inner(self,U,V):
        factors=np.ones(self.J+1);factors[1:]=2
        return float(2*np.pi*self.L/self.nz*np.real(np.sum(factors[:,None,None]*self.mass[None,:,None]*U.conj()*V)))

    def dissipation(self,U):
        total=0.
        for j,a in enumerate(U):
            ah=fft(a,axis=-1,norm='forward')
            total+=(1 if j==0 else 2)*np.real(np.sum(ah.conj()*(self.K[j]@ah+self.mass[:,None]*self.kz**2*ah)))
        return float(2*np.pi*self.L*total)

    def divergence_norm(self,U):
        total=0.
        for j,a in enumerate(U):
            d=self.div_hat(fft(a,axis=-1,norm='forward'),j)
            total+=(1 if j==0 else 2)*np.sum(self.V[:,None]*abs(d)**2)
        return float(np.sqrt(2*np.pi*self.L*total))

    def observe(self,U,t):
        a=U[0].real;zc=self.nz//2
        omega=float(self.axis_weights@(a[self.slt][:4,zc]/self.r[:4]))
        # Independent local centered sixth-order derivative for the core.
        h=self.L/self.nz;uz=a[self.slz]
        dz=(uz[:,zc+3]-9*uz[:,zc+2]+45*uz[:,zc+1]-45*uz[:,zc-1]+9*uz[:,zc-2]-uz[:,zc-3])/(60*h)
        b=.5*float(self.axis_weights@dz[:4])
        axial_phase=np.exp(1j*self.kz*(4-self.z[0]))
        receiver=np.real(fft(a[self.slt],axis=-1,norm='forward')@axial_phase)
        spline=CubicSpline(self.r,self.r*receiver)
        candidates=[x for x in spline.derivative().roots(extrapolate=False) if .18<x<.3]
        rstar=.21547557533348247
        rm=max(candidates,key=lambda x:float(spline(x))) if candidates else rstar
        fluctuating=U.copy();fluctuating[0]=0
        return dict(time=t,energy=.5*self.inner(U,U),K=.5*self.inner(fluctuating,fluctuating),
            discrete_dissipation=self.dissipation(U),divergence_L2=self.divergence_norm(U),
            core_b=b,core_Omega=omega,core_beta=b/omega,
            mean_G_fixed=float(spline(rstar)),mean_G_branch=float(spline(rm)),mean_G_branch_radius=float(rm),
            scope='Finite-cylinder discrete observables, not whole-space NS endpoint certificates.')

    def initial(self):
        U=np.zeros((self.J+1,self.nv,self.nz),complex)
        for i in range(0,self.nr,24):
            sl=slice(i,min(i+24,self.nr));local=fields(self.r[sl,None],self.z[None,:])
            for j,m in enumerate(self.modes):
                if m in local['u']:
                    U[j,self.nr-1+sl.start:self.nr-1+sl.stop]=local['u'][m][1]
                    U[j,2*self.nr-1+sl.start:2*self.nr-1+sl.stop]=local['u'][m][2]
        for i in range(0,self.nr-1,24):
            sl=slice(i,min(i+24,self.nr-1));local=fields(self.rf[1+sl.start:1+sl.stop,None],self.z[None,:])
            for j,m in enumerate(self.modes):
                if m in local['u']:U[j,sl]=local['u'][m][0]
        projected=self.project(U)
        return projected,dict(sampled_initial_energy=.5*self.inner(U,U),
            sampled_initial_core_Omega=float(np.real(self.axis_weights@(U[0,self.slt][:4,self.nz//2]/self.r[:4]))),
            removed_axial_mode_count=int(np.sum(~self.keep)),
            sampled_initial_divergence_L2=self.divergence_norm(U),
            projection_and_axial_filter_correction_L2=np.sqrt(self.inner(projected-U,projected-U)))


def run(nr=192,nz=512,R=8.,L=16.,mapping=4.,J=4,nu=.001,T=.001,dt=2e-5,output=None,dealias='padding',backend='fast'):
    started=time.time()
    if backend=='fast':
        from fast_full_mac import FastMAC
        grid=FastMAC(nr,nz,R,L,mapping,J,nu,dealias=dealias)
    else:grid=MAC(nr,nz,R,L,mapping,J,nu,dealias=dealias)
    U,initial=grid.initial();history=[grid.observe(U,0.)]
    raw=grid.nonlinear(U,project=False);projected=grid.project(raw)
    initial['unprojected_nonlinear_energy_work']=grid.inner(U,raw)
    initial['projected_nonlinear_energy_work']=grid.inner(U,projected)
    initial['projection_idempotence_L2']=np.sqrt(grid.inner(grid.project(U)-U,grid.project(U)-U))
    steps=max(1,int(np.ceil(T/dt)));dt=T/steps;viscous_loss=0.
    for i in range(steps):
        old_D=grid.dissipation(U);U=grid.step(U,dt);new_D=grid.dissipation(U)
        viscous_loss+=nu*dt*(old_D+new_D)/2
        if (i+1)%max(1,steps//10)==0 or i+1==steps:
            obs=grid.observe(U,(i+1)*dt);obs['energy_balance_defect']=obs['energy']-history[0]['energy']+viscous_loss
            history.append(obs)
            print(json.dumps(dict(step=i+1,steps=steps,seconds=time.time()-started,**obs)),flush=True)
        if not np.isfinite(U).all():raise RuntimeError('Nonfinite evolution; no endpoint accepted')
    result=dict(status='EXPLORATORY FULL 3D DISCRETE EVOLUTION; NOT VALIDATED WHOLE-SPACE NS',
        parameters=dict(nr=nr,nz=nz,R=R,L=L,mapping=mapping,J=J,nu=nu,T=T,dt=dt,steps=steps,dealias=dealias,backend=backend),
        initial_reconstruction=initial,history=history,seconds=time.time()-started,
        remaining_errors=['radial no-slip wall and periodic axial images','initial sampling/projection/filter errors',
            'axis regularity and continuous divergence-free reconstruction','angular/radial/axial truncation',
            'viscous splitting and time stepping','whole-space tails and space-time residual','endpoint inheritance'])
    if output:Path(output).write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name,typ,default in [('nr',int,192),('nz',int,512),('R',float,8.),('L',float,16.),('mapping',float,4.),('J',int,4),('nu',float,.001),('T',float,.001),('dt',float,2e-5)]:ap.add_argument('--'+name,type=typ,default=default)
    ap.add_argument('--output');ap.add_argument('--dealias',choices=['padding','filter'],default='padding')
    ap.add_argument('--backend',choices=['fast','reference'],default='fast');ap.add_argument('--fft-workers',type=int,default=4)
    args=vars(ap.parse_args());workers=args.pop('fft_workers')
    with set_workers(workers):run(**args)

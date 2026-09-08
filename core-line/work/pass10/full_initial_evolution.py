#!/usr/bin/env python3
"""Full 3D initial temporal coefficients on a mapped cylinder.

All initial components, axial envelopes and angular interactions are included.
Pressure uses Fourier-periodic z and a free-radial exterior condition. This
is a numerical approximation to whole-space initial coefficients, not a
validated compact-3D evolution. The periodic images, radial truncation,
source compatibility repair and differentiation errors are recorded/untested.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.fft import fft,ifft,fftfreq
from scipy.special import kve
from scipy.sparse import csr_matrix
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from initial_full_field import fields,cutoff


def differentiation(x,order=1,width=5):
    rows=[];cols=[];vals=[];n=len(x)
    for i in range(n):
        start=min(max(i-width//2,0),n-width)
        js=np.arange(start,start+width)
        scale=max(abs(x[js]-x[i]))
        yy=(x[js]-x[i])/scale
        rhs=np.zeros(width);rhs[order]=[1,1,2][order]
        weights=np.linalg.solve(np.array([yy**k for k in range(width)]),rhs)/scale**order
        rows.extend([i]*width);cols.extend(js);vals.extend(weights)
    return csr_matrix((vals,(rows,cols)),shape=(n,n))


class Cylinder:
    def __init__(self,nr=384,nz=1024,R=8.,L=16.,mapping=4.):
        self.nr,self.nz,self.R,self.L=nr,nz,R,L
        self.rf=R*np.sinh(mapping*np.arange(nr+1)/nr)/np.sinh(mapping)
        self.r=(self.rf[:-1]+self.rf[1:])/2
        self.vol=(self.rf[1:]**2-self.rf[:-1]**2)/2
        self.z=-L/2+L*np.arange(nz)/nz
        self.kz=2*np.pi*fftfreq(nz,d=L/nz)
        self.Dr=differentiation(self.r,1)
        self.Drr=differentiation(self.r,2)
        xx=self.r[:4]**2;scale=xx[-1]
        self.axis_weights=np.linalg.solve(np.array([(xx/scale)**j for j in range(4)]),np.array([1.,0,0,0]))
        self.repairs=[]

    def dr(self,a):
        return self.Dr@a
    def dz(self,a,n=1):
        return ifft(fft(a,axis=-1)*(1j*self.kz)**n,axis=-1)
    def scalar_lap(self,a,m):
        return self.Drr@a+self.dr(a)/self.r[:,None]+self.dz(a,2)-m*m*a/self.r[:,None]**2
    def gradient(self,u,m):
        rr=self.r[:,None]
        return np.stack((np.stack([self.dr(v) for v in u]),
            np.stack(((1j*m*u[0]-u[1])/rr,(1j*m*u[1]+u[0])/rr,1j*m*u[2]/rr)),
            np.stack([self.dz(v) for v in u])),axis=1)
    def vector_lap(self,u,m):
        out=np.stack([self.scalar_lap(v,m) for v in u])
        out[0]-=(u[0]+2j*m*u[1])/self.r[:,None]**2
        out[1]-=(u[1]-2j*m*u[0])/self.r[:,None]**2
        return out
    def pressure_gradient(self,p,m):
        return np.stack((self.dr(p),1j*m*p/self.r[:,None],self.dz(p)))

    def poisson(self,source,m,tag=''):
        """Approximate -Delta_m p=source; Fourier z, radial exterior DtN."""
        n,nz=self.nr,self.nz;r,rf,V=self.r,self.rf,self.vol
        rhs=fft(source,axis=1,norm='forward')
        repair=0.
        if m==0:
            moment=np.sum(V*rhs[:,0])
            bump=cutoff((r/(self.R/2))**2)
            repair=moment/np.sum(V*bump)
            rhs[:,0]-=repair*bump
            self.repairs.append(dict(tag=tag,raw_mean_source_radial_integral=[float(moment.real),float(moment.imag)],
                uniform_z_compatibility_repair_amplitude=[float(repair.real),float(repair.imag)]))
        lower=-rf[1:-1]/(V[1:]*(r[1:]-r[:-1]))
        upper=-rf[1:-1]/(V[:-1]*(r[1:]-r[:-1]))
        diagonal=np.zeros(n)
        diagonal[1:]-=lower
        diagonal[:-1]-=upper
        kk=abs(self.kz)
        kappa=np.zeros(nz)
        nonzero=kk>0
        if np.any(nonzero):
            x=kk[nonzero]*self.R
            kappa[nonzero]=kk[nonzero]*kve(abs(m-1),x)/kve(m,x)+m/self.R
        kappa[~nonzero]=m/self.R
        outer=self.R*kappa/(V[-1]*(1+kappa*(self.R-r[-1])))
        diagonal+=m*m/r**2
        cp=np.empty((n,nz),float);sol=np.empty_like(rhs)
        den=diagonal[0]+self.kz**2
        cp[0]=upper[0]/den;sol[0]=rhs[0]/den
        for i in range(1,n):
            den=diagonal[i]+self.kz**2
            if i==n-1:den=den+outer
            den=den-lower[i-1]*cp[i-1]
            val=rhs[i]-lower[i-1]*sol[i-1]
            if i==n-1 and m==0:
                # Gauge for the compatible zero Fourier mode, no logarithmic tail.
                den[0]=1.;val[0]=0.
            cp[i]=upper[i]/den if i<n-1 else 0
            sol[i]=val/den
        for i in range(n-2,-1,-1):sol[i]-=cp[i]*sol[i+1]
        return ifft(sol,axis=1,norm='forward')

    def norm2(self,coeff):
        return float(2*np.pi*self.L/self.nz*sum((1 if m==0 else 2)*np.sum(self.vol[:,None]*np.sum(abs(a)**2,axis=0)) for m,a in coeff.items()))
    def inner(self,left,right):
        return float(2*np.pi*self.L/self.nz*sum((1 if m==0 else 2)*np.real(np.sum(self.vol[:,None]*np.sum(np.conj(a)*right.get(m,0),axis=0))) for m,a in left.items()))
    def at(self,a,r,z,radial_order=0):
        # Fourier-interpolate z, then local cubic interpolation in radius.
        spec=fft(a,axis=-1,norm='forward')
        values=np.sum(spec*np.exp(1j*self.kz*(z+self.L/2)),axis=-1)
        lo=max(0,np.searchsorted(self.r,r)-4);hi=min(self.nr,lo+8)
        return complex(CubicSpline(self.r[lo:hi],values[lo:hi])(r,nu=radial_order))
    def axis(self,a,z,zorder=0):
        axial=self.axis_weights@a[:4]
        spec=fft(axial,norm='forward')
        return complex(np.sum(spec*(1j*self.kz)**zorder*np.exp(1j*self.kz*(z+self.L/2))))


def signed(coeff):
    result=dict(coeff)
    result.update({-m:np.conj(a) for m,a in coeff.items() if m})
    return result


def run(nr=384,nz=1024,R=8.,L=16.,mapping=4.,nu=.001,T=.001,second=True,output=None,save_fields=None,pressure='fd'):
    start=time.time()
    if pressure=='fd':
        from high_order_pressure import HighOrderCylinder
        grid=HighOrderCylinder(nr,nz,R,L,mapping)
    else:
        grid=Cylinder(nr,nz,R,L,mapping)
    shape=(nr,nz)
    U={m:np.empty((3,*shape),complex) for m in (0,4)}
    F={m:np.empty((3,*shape),complex) for m in (0,4,8)}
    source={m:np.empty(shape,complex) for m in (0,4,8)}
    gradU={m:np.empty((3,3,*shape),complex) for m in (0,4)}
    maximum_divergence=0.
    block=24
    for i in range(0,nr,block):
        jj=slice(i,min(i+block,nr))
        local=fields(grid.r[jj,None],grid.z[None,:])
        for m in U:
            U[m][:,jj]=local['u'][m]
            gradU[m][:,:,jj]=local['gradient'][m]
            maximum_divergence=max(maximum_divergence,float(np.max(abs(local['divergence'][m]))))
        for m in source:
            source[m][jj]=local['source'][m]
            F[m][:,jj]=nu*local['laplacian'].get(m,0)-local['convection'][m]
    P={m:grid.poisson(s,m,'p0') for m,s in source.items()}
    N0={m:F[m]-grid.pressure_gradient(P[m],m) for m in F}
    gN0={m:grid.gradient(a,m) for m,a in N0.items()}
    divN0={m:np.einsum('ii...->...',g) for m,g in gN0.items()}
    divN0_norm=float(np.sqrt(2*np.pi*L/nz*sum((1 if m==0 else 2)*np.sum(grid.vol[:,None]*abs(a)**2) for m,a in divN0.items())))
    rstar=.21547557533348247
    pzz0=grid.axis(P[0],0,zorder=2).real
    # Outside the supports the exact coefficient equals -grad p; the raw
    # numerical field includes periodic/truncation/differentiation errors.
    exterior=(grid.r[:,None]**2+grid.z[None,:]**2)>6.1**2
    div_exterior=max(float(np.max(abs(x[exterior]))) for x in divN0.values())
    K0=.5*grid.norm2({4:U[4]})
    K1=grid.inner({4:U[4]},N0)
    E0=.5*grid.norm2(U)
    E1=grid.inner(U,N0)
    grad_norm2=2*np.pi*L/nz*sum((1 if m==0 else 2)*np.sum(grid.vol[:,None]*np.sum(abs(g)**2,axis=(0,1))) for m,g in gradU.items())
    result=dict(status='EXPLORATORY full-3D initial coefficients on a mapped finite cylinder; not a validated whole-space evolution',
        parameters=dict(nr=nr,nz=nz,R=R,L=L,mapping=mapping,nu=nu,A=1020,lam=.25,T=T,pressure=pressure),
        analytic_initial_divergence_max=maximum_divergence,
        source_compatibility_repairs=grid.repairs,
        initial=dict(energy=E0,fluctuation_energy=K0,fluctuation_energy_derivative=K1,
            fractional_fluctuation_energy_derivative=K1/K0,
            energy_derivative=E1,viscous_energy_derivative=-nu*float(grad_norm2),
            energy_identity_relative_error=(E1+nu*float(grad_norm2))/max(1,abs(E1)),
            pressure_zz_at_core=pzz0,beta_prime=-4-pzz0/2,
            mean_G=1020*rstar*rstar*float(cutoff((rstar/.315)**2)),
            mean_G_prime=rstar*grid.at(N0[0][1],rstar,4).real,
            radial_mean_acceleration=grid.at(N0[0][0],rstar,4).real,
            vertical_mean_acceleration=grid.at(N0[0][2],rstar,4).real,
            axial_mean_strain_derivative=grid.at(grid.dz(N0[0][2]),rstar,4).real,
            mixed_pressure_radial_unit=[float((grid.at(grid.dr(P[4]),rstar,4)/(.25*1020/2)).real),float((grid.at(grid.dr(P[4]),rstar,4)/(.25*1020/2)).imag)],
            N0_L2=float(np.sqrt(grid.norm2(N0))),numerical_div_N0_L2=divN0_norm,
            numerical_div_N0_exterior_max=div_exterior))
    if second:
        # Exact continuous source is 2 tr(grad u0 grad N0), supported where
        # grad u0 is supported. All signed angular interactions are retained.
        gu,gn=signed(gradU),signed(gN0)
        us,ns=signed(U),signed(N0)
        S1={m:np.zeros(shape,complex) for m in (0,4,8,12)}
        C1={m:np.zeros((3,*shape),complex) for m in S1}
        for ma,ga in gu.items():
            for mb,gb in gn.items():
                m=ma+mb
                if m not in S1:continue
                S1[m]+=2*np.einsum('ij...,ji...->...',ga,gb)
                C1[m]+=np.einsum('ij...,j...->i...',ga,ns[mb])+np.einsum('ij...,j...->i...',gb,us[ma])
        P1={m:grid.poisson(s,m,'p1') for m,s in S1.items()}
        # Delta N0 = Delta F + grad(source). This cancels the harmonic tail
        # exactly at the continuous level, avoiding numerical differentiation
        # of its Laplacian at large radius.
        lapN0={m:grid.vector_lap(F[m],m)+grid.pressure_gradient(source[m],m) for m in F}
        N1={m:nu*lapN0.get(m,0)-C1[m]-grid.pressure_gradient(P1[m],m) for m in S1}
        pzz1=grid.axis(P1[0],0,zorder=2).real
        beta1=-4-pzz0/2
        b1_direct=.5*grid.axis(N0[0][2],0,zorder=1).real
        b2_direct=.5*grid.axis(N1[0][2],0,zorder=1).real
        omega1_direct=float(np.real(grid.axis_weights@(N0[0][1,:4,nz//2]/grid.r[:4])))
        omega2_direct=float(np.real(grid.axis_weights@(N1[0][1,:4,nz//2]/grid.r[:4])))
        approximate={m:U.get(m,0)+T*N0.get(m,0)+.5*T*T*N1[m] for m in N1}
        uz=approximate[0][2]
        bT=.5*grid.axis(uz,0,zorder=1).real
        # The axisymmetric azimuthal component divided by r gives Omega.
        omegaT=float(np.real(grid.axis_weights@(approximate[0][1,:4,nz//2]/grid.r[:4])))
        Gt=grid.r*np.real(approximate[0][1,:,int(round((4+L/2)*nz/L))%nz])
        spl=CubicSpline(grid.r,Gt)
        roots=spl.derivative().roots(extrapolate=False)
        candidates=[r for r in roots if .18<r<.3]
        rmax=max(candidates,key=lambda r:float(spl(r))) if candidates else rstar
        result['second_coefficient']=dict(pressure_zz_time_derivative=pzz1,pressure_T=pzz1+32,
            beta_second=-(pzz1+32)/2-10*beta1,
            mean_G_second=rstar*grid.at(N1[0][1],rstar,4).real,
            N1_L2=float(np.sqrt(grid.norm2(N1))))
        result['direct_core_jet_checks']=dict(b1=b1_direct,b2=b2_direct,
            Omega1=omega1_direct,Omega2=omega2_direct,
            beta1=b1_direct-omega1_direct,
            beta2=b2_direct-omega2_direct-2*b1_direct*omega1_direct+2*omega1_direct**2,
            b1_pressure_discrepancy=b1_direct-(2+beta1),
            b2_pressure_discrepancy=b2_direct-(-4*b1_direct-pzz1/2),
            Omega2_ode_discrepancy=omega2_direct-(2*b1_direct+4))
        result['quadratic_field_at_T']=dict(core_b=bT,core_Omega=omegaT,core_beta=bT/omegaT,
            mean_fixed_receiver=float(spl(rstar)),mean_radial_branch_value=float(spl(rmax)),
            mean_radial_branch_radius=float(rmax),fluctuation_energy=.5*grid.norm2({m:a for m,a in approximate.items() if m}),
            energy=.5*grid.norm2(approximate),
            scope='Values of the full quadratic time polynomial only, not actual NS endpoint values or error bounds.')
        angular_u={m:1j*m*a for m,a in U.items() if m}
        optimum=-grid.inner(angular_u,N0)/grid.norm2(angular_u)
        corot=[]
        for rate in [1020.,optimum]:
            c1={m:N0.get(m,0)+1j*m*rate*U.get(m,0) for m in N1}
            c2={m:N1[m]+2j*m*rate*N0.get(m,0)-(m*rate)**2*U.get(m,0) for m in N1}
            approx={m:U.get(m,0)+T*c1[m]+.5*T*T*c2[m] for m in N1}
            corot.append(dict(rotation_rate=rate,first_coefficient_L2=np.sqrt(grid.norm2(c1)),
                second_coefficient_L2=np.sqrt(grid.norm2(c2)),
                fluctuation_energy_at_T=.5*grid.norm2({m:a for m,a in approx.items() if m}),
                total_energy_at_T=.5*grid.norm2(approx)))
        result['corotating_quadratic_diagnostics']=dict(candidates=corot,
            scope='Exact rotation formulas applied to approximate initial jets. Positive-time polynomial values still have no NS residual/error enclosure; finite-cylinder L2 optimum only.')
        if save_fields:
            saved=dict(r=grid.r,z=grid.z,rf=grid.rf)
            for tag,dd in [('U',U),('N0',N0),('N1',N1),('p0',P),('p1',P1)]:
                saved.update({f'{tag}_m{m}':a for m,a in dd.items()})
            np.savez_compressed(save_fields,**saved)
    result['source_compatibility_repairs']=grid.repairs
    result['seconds']=time.time()-start
    result['unresolved_errors']=['periodic axial images','finite-cylinder velocity norms omit exterior tails; later pressure sources need exterior treatment',
        'legacy source compatibility repair and inconsistent subsequent identities' if pressure=='fv' else 'nonzero numerical zero-mode exterior flux',
        'radial/axial quadrature and differentiation; noncommuting gradient/divergence operators',
        'whole-space smooth divergence-free reconstruction','space-time H4 residual and endpoint inheritance']
    if output:Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name,typ,default in [('nr',int,384),('nz',int,1024),('R',float,8.),('L',float,16.),('mapping',float,4.),('nu',float,.001),('T',float,.001)]:ap.add_argument('--'+name,type=typ,default=default)
    ap.add_argument('--first-only',action='store_true');ap.add_argument('--output');ap.add_argument('--save-fields')
    ap.add_argument('--pressure',choices=['fv','fd'],default='fd')
    args=vars(ap.parse_args());args['second']=not args.pop('first_only')
    run(**args)

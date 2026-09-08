#!/usr/bin/env python3
"""Exploratory whole-space axisymmetric initial-pressure quadrature.

No time stepping, interval certification, or claim of a singular solution.
Profiles are C-infinity flat cutoffs. Pressure uses every retained even
Legendre mode, with analytic radial Green functions on linear source panels.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse, json, time
import numpy as np
import sympy as sp
from scipy.special import expit, eval_legendre
from scipy.integrate import simpson


def cutoff(x, n=0):
    x = np.asarray(x)
    t = (x-.25)/.75
    inside = (t > 0) & (t < 1)
    a = np.clip(t, 1.e-7, 1-1.e-7)
    z = 1/a-1/(1-a)
    f = expit(z)
    z1 = -a**-2-(1-a)**-2
    z2 = 2*a**-3-2*(1-a)**-3
    z3 = -6*a**-4-6*(1-a)**-4
    if n == 0:
        return np.where(t <= 0, 1., np.where(t >= 1, 0., f))
    if n == 1:
        ans = f*(1-f)*z1/.75
    elif n == 2:
        ans = f*(1-f)*((1-2*f)*z1*z1+z2)/.75**2
    elif n == 3:
        ans = f*(1-f)*((1-6*f+6*f*f)*z1**3+3*(1-2*f)*z1*z2+z3)/.75**3
    else:
        raise ValueError(n)
    return np.where(inside, ans, 0.)


def profile_functions(sigma, height, distance, swirl_ratio):
    r,z = sp.symbols('r z', positive=True)
    f = sp.Function('f')
    rr = r*r+z*z
    psi = f(rr)
    ss = (r/sigma)**2
    alpha = f(ss)
    beta = f(((z-distance)/height)**2)+f(((z+distance)/height)**2)
    stream = -z*r*r*alpha*beta
    outer_swirl = r*f((r/(sigma*swirl_ratio))**2)*(f(((z-distance)/(height*swirl_ratio))**2)+f(((z+distance)/(height*swirl_ratio))**2))
    fields = [(-r*(psi+z*sp.diff(psi,z)), 0, 2*z*psi+r*z*sp.diff(psi,r)),
              (0, r*(psi+sp.Rational(1,3)*(r*sp.diff(psi,r)+z*sp.diff(psi,z))), 0),
              (-sp.diff(stream,z)/r, 0, sp.diff(stream,r)/r),
              (0, outer_swirl, 0)]
    modules = {'f': lambda x:cutoff(x,0)}
    for n in range(1,4):
        modules['f'+str(n)] = lambda x,n=n:cutoff(x,n)
    def lower(expr):
        replacements = {}
        for atom in expr.atoms(sp.Subs):
            if isinstance(atom.expr,sp.Derivative) and atom.expr.expr.func == f:
                n = sum(k for _,k in atom.expr.variable_count)
                replacements[atom] = sp.Function('f'+str(n))(atom.point[0])
        return expr.xreplace(replacements)
    result=[]
    for field in fields:
        exprs=[]
        for component in field:
            component=sp.sympify(component)
            exprs.extend([component,sp.diff(component,r),sp.diff(component,z),sp.diff(component,r,2),sp.diff(component,r,z),sp.diff(component,z,2)])
        result.append(sp.lambdify((r,z),[lower(e) for e in exprs],modules=[modules,'numpy'],cse=True))
    return result


class Evaluator:
    def __init__(self,nr=600,nmu=256,lmax=128,sigma=.7,height=1.,distance=4.,swirl_ratio=.45):
        start=time.time()
        if distance-height <= 1 or not 0 < swirl_ratio < .5:
            raise ValueError('Require disjoint core/outer supports and swirl strictly inside pump plateaus.')
        if sigma*swirl_ratio > (distance-height*swirl_ratio)/4:
            raise ValueError('Outer swirl must fit inside r <= |z|/4.')
        if (distance+height+.25)**2 <= sigma*sigma+(distance+height)**2:
            raise ValueError('Radial ceiling does not contain meridional support.')
        self.R=(np.arange(nr)+.5)*(distance+height+.25)/nr
        self.mu,self.weights=np.polynomial.legendre.leggauss(nmu)
        self.RR=self.R[:,None]
        self.MM=self.mu[None,:]
        self.sn=np.sqrt(1-self.MM**2)
        self.r=self.RR*self.sn
        self.z=self.RR*self.MM
        self.ells=np.arange(0,lmax+1,2)
        self.PL=np.array([eval_legendre(int(l),self.mu) for l in self.ells])
        self.PL1=np.array([np.zeros_like(self.mu) if l==0 else l*(self.mu*self.PL[j]-eval_legendre(int(l-1),self.mu))/(self.mu**2-1) for j,l in enumerate(self.ells)])
        self.PL2=(2*self.mu*self.PL1-self.ells[:,None]*(self.ells[:,None]+1)*self.PL)/(1-self.mu**2)
        self.project=((2*self.ells+1)[:,None]*self.weights[None,:]*self.PL/2).T
        self.fields=[]
        for fn in profile_functions(sigma,height,distance,swirl_ratio):
            self.fields.append(np.array([np.broadcast_to(a,self.r.shape) for a in fn(self.r,self.z)]).reshape(3,6,nr,nmu))
        self.fields=np.array(self.fields)
        self.metadata=dict(nr=nr,nmu=nmu,lmax=lmax,sigma=sigma,height=height,distance=distance,swirl_ratio=swirl_ratio)
        self.pressures=[self.pressure_source_hzz(self.source(self.combine(np.eye(4)[j])),j==0,j==1) for j in range(4)]
        self.Cp=self.outer_stress(2)
        self.Cv=self.outer_stress(3)
        self.dP=self.outer_stress(2,gradient=True)
        self.dv=self.outer_stress(3,gradient=True)
        self.setup_seconds=time.time()-start

    def combine(self,amplitudes):
        return np.einsum('i,icdrm->cdrm',amplitudes,self.fields)

    def outer_stress(self,j,gradient=False):
        if gradient:
            G=self.gradient(self.fields[j])
            tensor=np.einsum('ikrm,jkrm->ijrm',G,G)
            return 2*self.integrate_outer_tensor(tensor)
        u=self.fields[j,:,0]
        return self.integrate_outer_tensor(np.einsum('irm,jrm->ijrm',u,u))

    def integrate_outer_tensor(self,tensor):
        r=self.r; z=self.z; R=self.RR
        den=4*np.pi*R**9
        Qrr=-(105*r*r*z*z-15*R*R*(z*z+r*r)+3*R**4)/den
        Qtt=(15*z*z-3*R*R)/(4*np.pi*R**7)
        Qzz=-(105*z**4-90*R*R*z*z+9*R**4)/den
        Qrz=-(105*r*z**3-45*R*R*r*z)/den
        stress=Qrr*tensor[0,0]+Qtt*tensor[1,1]+Qzz*tensor[2,2]+Qrz*(tensor[0,2]+tensor[2,0])
        integrand=2*np.pi*R*R*stress
        return float(simpson(integrand@self.weights,x=self.R))

    def gradient(self,U):
        v,w,q=U[:,0]
        vr,wr,qr=U[:,1]
        vz,wz,qz=U[:,2]
        zero=np.zeros_like(v)
        return np.array([[vr,-w/self.r,vz],[wr,v/self.r,wz],[qr,zero,qz]])

    def source(self,U):
        G=self.gradient(U)
        return np.einsum('ijrm,jirm->rm',G,G)

    def pressure_source_hzz(self,g,isS=False,isW=False,contact=None):
        gl=g@self.project
        if contact is None:
            contact = 6. if isS else (-2. if isW else 0.)
        # The source is exactly constant in the affine core for initial pressure.
        if isS or isW:
            gl[self.R<.5,1:]=0.
        return float(-contact/3+2/5*simpson(gl[:,1]/self.R,x=self.R))

    def pressure_hessian(self,g):
        gl=g@self.project
        gl[self.R<.5,1:]=0.  # All input combinations have the same affine ball.
        C4=float(simpson(gl[:,2]/self.R**3,x=self.R)/9)
        l=self.ells.astype(float)
        a=np.zeros_like(gl); b=np.zeros_like(gl)
        a[0,0]=self.R[0]**2*gl[0,0]/3
        for i in range(1,len(self.R)):
            R=self.R[i]; q=self.R[i-1]/R
            slope=(gl[i]-gl[i-1])/(R-self.R[i-1])
            h3=-np.expm1((l+3)*np.log(q))/(l+3)
            h4=-np.expm1((l+4)*np.log(q))/(l+4)
            a[i]=q**(l+1)*a[i-1]+R**2*(gl[i]*h3+slope*R*(h4-h3))
        def moment(n,logq):
            safe=np.where(n==0,1.,n)
            return np.where(n==0,logq,np.expm1(n*logq)/safe)
        for i in range(len(self.R)-2,-1,-1):
            R=self.R[i]; q=self.R[i+1]/R
            slope=(gl[i+1]-gl[i])/(self.R[i+1]-R)
            h2=moment(2-l,np.log(q)); h3=moment(3-l,np.log(q))
            b[i]=q**(-l)*b[i+1]+R**2*(gl[i]*h2+slope*R*(h3-h2))
        pl=(a+b)/(2*l+1)
        pl1=(-(l+1)*a+l*b)/((2*l+1)*self.RR)
        pl2=-gl-2*pl1/self.RR+l*(l+1)*pl/self.RR**2
        pR=pl1@self.PL; pRR=pl2@self.PL
        pm=pl@self.PL1; pRm=pl1@self.PL1; pmm=pl@self.PL2
        R=self.RR; m=self.MM; s=self.sn
        hrr=s*s*pRR-2*m*s*s*pRm/R+m*m*s*s*pmm/R**2+m*m*pR/R+m*(2-3*m*m)*pm/R**2
        hzz=m*m*pRR+2*m*s*s*pRm/R+s**4*pmm/R**2+s*s*pR/R-3*m*s*s*pm/R**2
        hrz=m*s*pRR+s*(1-2*m*m)*pRm/R-m*s**3*pmm/R**2-m*s*pR/R+s*(3*m*m-1)*pm/R**2
        htt=pR/R-m*pm/R**2
        zero=np.zeros_like(g)
        H=np.array([[hrr,zero,hrz],[zero,htt,zero],[hrz,zero,hzz]])
        residual=float(np.max(np.abs(hrr+htt+hzz+g)))
        pr=s*pR-m*s*pm/R
        pz=m*pR+s*s*pm/R
        return H,residual,residual/max(1.,float(np.max(np.abs(g)))),C4,pr,pz

    def evaluate(self,b=1.,omega=1.,c=1.,A=None):
        if A is None:
            A=np.sqrt(((38/7)*b*b+(2/5)*omega*omega-c*c*self.Cp)/self.Cv)
        U=self.combine([b,omega,c,A]); G=self.gradient(U)
        v,w,q=U[:,0]; vr,wr,qr=U[:,1]; vz,wz,qz=U[:,2]
        vrr,wrr,qrr=U[:,3]; vrz,wrz,qrz=U[:,4]; vzz,wzz,qzz=U[:,5]
        r=self.r
        # Cylindrical components of (u dot grad)u and their r,z derivatives.
        ar=v*vr+q*vz-w*w/r
        at=v*wr+q*wz+v*w/r
        az=v*qr+q*qz
        arr=vr*vr+v*vrr+qr*vz+q*vrz-2*w*wr/r+w*w/r**2
        arz=vz*vr+v*vrz+qz*vz+q*vzz-2*w*wz/r
        atr=vr*wr+v*wrr+qr*wz+q*wrz+(vr*w+v*wr)/r-v*w/r**2
        atz=vz*wr+v*wrz+qz*wz+q*wzz+(vz*w+v*wz)/r
        azr=vr*qr+v*qrr+qr*qz+q*qrz
        azz=vz*qr+v*qrz+qz*qz+q*qzz
        zero=np.zeros_like(v)
        Ga=np.array([[arr,-at/r,arz],[atr,ar/r,atz],[azr,zero,azz]])
        source=np.einsum('ijrm,jirm->rm',G,G)
        H,poisson_error,poisson_relative,C4,pr,pz=self.pressure_hessian(source)
        G1=-Ga-H
        gt=2*np.einsum('ijrm,jirm->rm',G,G1)
        # Exact center jet from known quadratic pressure coefficients, independent
        # of small-R numerical Hessians. This also checks approximate neutrality.
        pzz=-(18/7)*b*b+(2/5)*omega*omega-c*c*self.Cp-A*A*self.Cv
        pp=(-6*b*b+2*omega*omega-pzz)/2
        L=np.array([[-b,-omega,0],[omega,-b,0],[0,0,2*b]])
        L1=-L@L-np.diag([pp,pp,pzz])
        gt0=2*np.trace(L@L1)
        gt_l2=gt@self.project[:,1]
        # In the affine ball, the l=2 part of g_t comes exactly from the
        # quartic harmonic of p: (g_t)_2=-72 b C4 R^2. This replaces origin
        # cancellation and includes the unsampled interval [0,R[0]].
        gt_l2[self.R<.5]=-72*b*C4*self.R[self.R<.5]**2
        pzz1=-gt0/3+2/5*simpson(gt_l2/self.R,x=self.R)-(72/5)*b*C4*self.R[0]**2
        # Independent integration-by-parts representation on separated outer
        # support: 2Pi(u_outer,u1)=-2 integral Q_ij u_outer_i u1_j.
        Gcore=self.gradient(self.combine([b,omega,0.,0.]))
        hcore=2*np.einsum('ijrm,jirm->rm',Gcore,G1)
        hc2=hcore@self.project[:,1]
        hc2[self.R<.5]=-72*b*C4*self.R[self.R<.5]**2
        core_part=-gt0/3+2/5*simpson(hc2/self.R,x=self.R)-(72/5)*b*C4*self.R[0]**2
        uouter=c*self.fields[2,:,0]+A*self.fields[3,:,0]
        accel=np.array([-ar-pr,-at,-az-pz])
        outer_part=-2*self.integrate_outer_tensor(np.einsum('irm,jrm->ijrm',uouter,accel))
        pzz1_split=float(core_part+outer_part)
        viscous=self.dP*c*c+self.dv*A*A
        T=float(pzz1+32*b**3)
        return dict(b=b,omega=omega,c=c,A=float(A),pzz=float(pzz),pzz_prime=float(pzz1),T=T,pzz_prime_split=pzz1_split,T_split=pzz1_split+32*b**3,representation_difference=pzz1_split-float(pzz1),core_part=float(core_part),outer_part=float(outer_part),viscous_coefficient=float(viscous),nu_gate_threshold=float(-T/viscous) if T<0 and viscous>0 else None,contact=float(-gt0/3),poisson_trace_error=poisson_error,poisson_trace_relative=poisson_relative)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--nr',type=int,default=600)
    ap.add_argument('--nmu',type=int,default=256)
    ap.add_argument('--lmax',type=int,default=128)
    ap.add_argument('--sigma',type=float,default=.7)
    ap.add_argument('--height',type=float,default=1.)
    ap.add_argument('--distance',type=float,default=4.)
    ap.add_argument('--swirl-ratio',type=float,default=.45)
    ap.add_argument('--c',type=float,nargs='+',default=[0.,.7,1.,2.,4.])
    args=ap.parse_args()
    params=vars(args).copy(); cs=params.pop('c')
    ev=Evaluator(**params)
    print(json.dumps(dict(metadata=ev.metadata,pressures=ev.pressures,Cp=ev.Cp,Cv=ev.Cv,dP=ev.dP,dv=ev.dv,setup_seconds=ev.setup_seconds)),flush=True)
    for c in cs:
        print(json.dumps(ev.evaluate(c=c)),flush=True)


if __name__=='__main__':
    main()

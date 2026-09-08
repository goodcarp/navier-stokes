#!/usr/bin/env python3
"""Exploratory actual initial Gtt from a finite-box mixed-pressure solve.

The formula retains the complete m=4 pressure, but the numerical elliptic
solve and quadrature are not certified. No NS time integration is performed.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,argparse,time
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.fft import dst,idst
from scipy.interpolate import RectBivariateSpline
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pass7'))
from evaluate_outer_E import fields

mp.mp.dps=50
av=mp.mpf(63)/200
def C(r):
    t=(4*(r/av)**2-1)/3
    if t<=0:return mp.mpf(1)
    if t>=1:return mp.mpf(0)
    return 1/(1+mp.exp(-1/t+1/(1-t)))
def gamma_unit(r):return r*r*C(r)
rstar=mp.findroot(lambda r:mp.diff(gamma_unit,r),(mp.mpf('.2'),mp.mpf('.24')))
Gjets=[float(mp.diff(gamma_unit,rstar,j)) for j in range(5)]

def run(nr=600,nz=3600,R=1.2,Z=6.,A=900.,lam=1.):
    start=time.time();dr=R/(nr+1);dz=2*Z/(nz+1)
    r=(np.arange(nr)+1)*dr;z=-Z+(np.arange(nz)+1)*dz
    _,_,s=fields(r[:,None],z[None,:])
    rhs=dst(s*np.sqrt(r[:,None]),type=1,axis=1,norm='ortho')
    kz2=4/dz**2*np.sin(np.pi*(np.arange(nz)+1)/(2*(nz+1)))**2
    off=-1/dr**2;diag=2/dr**2+(16-.25)/r**2
    cp=np.empty((nr,nz));fp=np.empty_like(rhs)
    denom=diag[0]+kz2;cp[0]=off/denom;fp[0]=rhs[0]/denom
    for i in range(1,nr):
        denom=diag[i]+kz2-off*cp[i-1]
        cp[i]=off/denom;fp[i]=(rhs[i]-off*fp[i-1])/denom
    for i in range(nr-2,-1,-1):fp[i]-=cp[i]*fp[i+1]
    p=idst(fp,type=1,axis=1,norm='ortho')/np.sqrt(r[:,None])
    rr=float(rstar);zz=4.
    # Local interpolation uses only a small smooth patch of the computed pressure.
    ir=np.where(abs(r-rr)<.04)[0];iz=np.where(abs(z-zz)<.06)[0]
    assert len(ir)>4 and len(iz)>4
    pr=0j
    for part,factor in [(p.real,1),(p.imag,1j)]:
        spline=RectBivariateSpline(r[ir],z[iz],part[np.ix_(ir,iz)],kx=3,ky=3,s=0)
        pr+=factor*float(spline(rr,zz,dx=1)[0,0])
    _,_,ss=fields(np.array([[rr]]),np.array([[zz]]));sm=ss[0,0]
    m,k,c,nu=4.,20.,1.4,.001
    G=[A*x for x in Gjets];T0=lam*lam*m*k/(2*rr)
    L2=c*c*rr*rr*G[2]-2*c*nu*rr*G[3]+nu*nu*(G[4]-2*G[3]/rr+3*G[2]/rr**2)
    curvature=lam*lam*m*m*G[2]/(2*rr**2)
    visc=nu*T0*((2*m*m+4)/rr**2-2*k*k)
    pressure=lam*lam*A/2*np.real(np.exp(-1j*k*rr)*((k*k*rr+m*m/rr+1j*k)*pr-1j*k*rr*sm))
    Gt=T0+nu*G[2]
    Gtt=L2+curvature+visc+pressure
    return dict(nr=nr,nz=nz,R=R,Z=Z,A=A,lam=lam,rstar=str(rstar),cutoff_at_max=str(C(rstar)),unit_G_jets=Gjets,
        pressure_radial=[float(pr.real),float(pr.imag)],source=[float(sm.real),float(sm.imag)],
        G0=G[0],Gt=Gt,Gtt=float(Gtt),L2=float(L2),curvature=float(curvature),viscous_torque=float(visc),
        pressure_contribution=float(pressure),linearized_time_Gt_zero=float(-Gt/Gtt),
        seconds=time.time()-start,status='EXPLORATORY finite-box pressure only; A argument need not be exactly neutral. No actual finite duration follows.')

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--nr',type=int,default=600);ap.add_argument('--nz',type=int,default=3600)
    ap.add_argument('--R',type=float,default=1.2);ap.add_argument('--Z',type=float,default=6.)
    ap.add_argument('--A',type=float,default=900.);ap.add_argument('--lam',type=float,default=1.)
    print(json.dumps(run(**vars(ap.parse_args()))),flush=True)

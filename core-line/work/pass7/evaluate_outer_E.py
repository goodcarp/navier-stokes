#!/usr/bin/env python3
"""Exploratory m=4 cylindrical pressure solve; finite box, not certified."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,argparse,time
from pathlib import Path
import numpy as np
from scipy.fft import dst,idst
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'outputs/euler-ns-transfer/pass4/experiments'))
from evaluate_pressure_gate import cutoff

def fields(r,z):
    m,k=4.,20.
    x=(r-.21)**2/.14**2
    e=cutoff(x); er=cutoff(x,1)*2*(r-.21)/.14**2
    err=cutoff(x,2)*(2*(r-.21)/.14**2)**2+cutoff(x,1)*2/.14**2
    q=sum(cutoff((z-sign*4)**2/.6**2) for sign in (-1,1))
    qz=sum(cutoff((z-sign*4)**2/.6**2,1)*2*(z-sign*4)/.6**2 for sign in (-1,1))
    a=e*q; ar=er*q; arr=err*q; az=e*qz
    vx=r*r/.315**2
    vz=sum(cutoff((z-sign*4)**2/.45**2) for sign in (-1,1))
    V=r*cutoff(vx)*vz
    Vr=(cutoff(vx)+2*vx*cutoff(vx,1))*vz
    R2=r*r+z*z
    Qt=3*(4*z*z-r*r)/(4*np.pi*R2**3.5)
    Qtr=15*r*(r*r-6*z*z)/(4*np.pi*R2**4.5)
    Qrz=-15*r*z*(4*z*z-3*r*r)/(4*np.pi*R2**4.5)
    phase=np.exp(1j*k*r)
    d=phase*(r*Qtr*(1j*m*(ar/r-a/r**2)-m*k*a/r)+Qrz*1j*m*az/r)
    s=2*phase*(-m*m*Vr*a/r**2+Vr*(ar+1j*k*a)/r+V*(arr+2j*k*ar-k*k*a)/r)
    local=6*np.pi*m*k*a*a*V*Qtr
    return local,d,s

def evaluate(nr,nz,R,Z):
    start=time.time(); dr=R/(nr+1); dz=2*Z/(nz+1)
    r=(np.arange(nr)+1)[:,None]*dr
    z=-Z+(np.arange(nz)+1)[None,:]*dz
    loc,d,s=fields(r,z)
    # p=f/sqrt(r). Cylindrical Fourier operator becomes positive symmetric.
    rhs=dst(s*np.sqrt(r),type=1,axis=1,norm='ortho')
    kz2=4/dz**2*np.sin(np.pi*(np.arange(nz)+1)/(2*(nz+1)))**2
    off=-1/dr**2
    diag=2/dr**2+(16-.25)/r[:,0]**2
    cp=np.empty((nr,nz),dtype=float)
    fp=np.empty_like(rhs)
    denom=diag[0]+kz2
    cp[0]=off/denom; fp[0]=rhs[0]/denom
    for i in range(1,nr):
        denom=diag[i]+kz2-off*cp[i-1]
        cp[i]=off/denom
        fp[i]=(rhs[i]-off*fp[i-1])/denom
    for i in range(nr-2,-1,-1): fp[i]-=cp[i]*fp[i+1]
    p=idst(fp,type=1,axis=1,norm='ortho')/np.sqrt(r)
    Eloc=loc.sum()*dr*dz
    Epressure=-2*np.pi*np.sum(r*np.real(np.conj(d)*p))*dr*dz
    # Exact whole-space energy inequality, with these numerical integrals only exploratory.
    nd=np.sqrt(np.pi*np.sum(r**3*abs(d)**2)*dr*dz)
    ns=np.sqrt(np.pi*np.sum(r**3*abs(s)**2)*dr*dz)
    return dict(nr=nr,nz=nz,R=R,Z=Z,E_local=float(Eloc),E_pressure=float(Epressure),E=float(Eloc+Epressure),pressure_absolute_energy_bound_numeric=float(2*nd*ns/16),d_r_L2=float(nd),s_r_L2=float(ns),seconds=time.time()-start,status='EXPLORATORY finite-box quadrature; not a rigorous enclosure')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--nr',type=int,default=300);ap.add_argument('--nz',type=int,default=1800);ap.add_argument('--R',type=float,default=1.2);ap.add_argument('--Z',type=float,default=6.)
    print(json.dumps(evaluate(**vars(ap.parse_args()))),flush=True)

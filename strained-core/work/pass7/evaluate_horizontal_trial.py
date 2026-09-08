#!/usr/bin/env python3
"""Exploratory separable whole-space pressure trial, no certified bounds."""
import sys,json
import numpy as np
from scipy.integrate import cumulative_trapezoid,simpson
from evaluate_outer_E import cutoff,fields

def run(nr=4000,nz=1600):
    m=4; R=.35
    r=np.linspace(1e-7,R,nr); z=np.linspace(3.4,4.6,nz)
    _,_,S=fields(r[:,None],np.array([[4.]])); S=S[:,0]
    left=cumulative_trapezoid(r**(m+1)*S,r,initial=0)
    right=-cumulative_trapezoid((r**(1-m)*S)[::-1],r[::-1],initial=0)[::-1]
    P=(r**(-m)*left+r**m*right)/(2*m)
    q=cutoff((z-4)**2/.6**2); vq=cutoff((z-4)**2/.45**2)
    q1=cutoff((z-4)**2/.6**2,1)*2*(z-4)/.6**2
    vq1=cutoff((z-4)**2/.45**2,1)*2*(z-4)/.45**2
    q2=cutoff((z-4)**2/.6**2,2)*(2*(z-4)/.6**2)**2+cutoff((z-4)**2/.6**2,1)*2/.6**2
    vq2=cutoff((z-4)**2/.45**2,2)*(2*(z-4)/.45**2)**2+cutoff((z-4)**2/.45**2,1)*2/.45**2
    F=q*vq; F2=q2*vq+2*q1*vq1+q*vq2
    p=P[:,None]*F[None,:]
    loc,d,s=fields(r[:,None],z[None,:])
    def int2(x): return 2*simpson(simpson(x,x=z,axis=1),x=r)
    nd=np.sqrt(np.pi*int2(r[:,None]**3*abs(d)**2))
    radial_norm=simpson(r**3*abs(P)**2,x=r)+abs(left[-1]/(2*m))**2*R**(4-2*m)/(2*m-4)
    ne=np.sqrt(radial_norm*2*np.pi*simpson(F2*F2,x=z))
    press=-2*np.pi*int2(r[:,None]*np.real(np.conj(d)*p))
    return dict(nr=nr,nz=nz,E_local=float(int2(loc)),Epressure_trial=float(press),residual_weighted_L2=float(ne),d_weighted_L2=float(nd),error_bound_numerical=float(2*nd*ne/m**2),E_upper_numerical=float(int2(loc)+press+2*nd*ne/m**2),radial_P_weighted_square=float(radial_norm),axial_Fpp_square=float(2*simpson(F2*F2,x=z)),status='EXPLORATORY quadrature. The energy error formula is exact but its numerical inputs are not enclosed.')

if __name__=='__main__': print(json.dumps(run(*map(int,sys.argv[1:]))),flush=True)

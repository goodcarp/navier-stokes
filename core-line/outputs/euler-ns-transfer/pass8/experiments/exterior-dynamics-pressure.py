#!/usr/bin/env python3
"""Nonsingular whole-space endcap quadrature for unit pure-swirl pressure.

Exploratory numerical evaluation, not interval certification. It integrates
both physical packets and the full angular Newton kernel, without a box or
prescribed boundary pressure. Exact formulas are documented separately.
"""
import argparse,json
import numpy as np
from scipy.special import expit
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

a=.315; h=.45
def cutoff(x,derivative=False):
    x=np.asarray(x,dtype=float);result=np.zeros_like(x)
    if not derivative:result[x<=.25]=1
    active=(x>.25)&(x<1);t=(4*x[active]-1)/3
    p=expit(1/t-1/(1-t))
    result[active]=(-4/3*p*(1-p)*(t**-2+(1-t)**-2)) if derivative else p
    return result

def quadrature(n,lo,hi):
    x,w=leggauss(n);return (lo+hi)/2+(hi-lo)*x/2,(hi-lo)*w/2

def run(n=64,r=None):
    xi=brentq(lambda x:float(cutoff(x)+x*cutoff(x,True)),.25,.625,xtol=1e-14)
    if r is None:r=a*np.sqrt(xi)
    sr,sw=quadrature(n,0,a);theta,tw=quadrature(n,0,2*np.pi)
    tau,vw=quadrature(n,h/2,h)
    u,uw=quadrature(max(n,96),0,1)
    lower=(sr/a)**2
    # H(sr)=(a²/2) integral_{(sr/a)²}^1 chi(x)² dx.
    nodes=lower[:,None]+(1-lower[:,None])*u
    H=(a*a/2)*(1-lower)*np.sum(cutoff(nodes)**2*uw,axis=1)
    weight=H[:,None]*sr[:,None]*sw[:,None]*tw[None,:]
    dx=r-sr[:,None]*np.cos(theta)[None,:]
    planar=r*r+sr[:,None]**2-2*r*sr[:,None]*np.cos(theta)[None,:]
    labels=['delta','delta_r','delta_z','delta_zz','delta_rr']
    contributions={}
    for center in [4.,-4.]:
        total=np.zeros(len(labels))
        for t,wt in zip(tau,vw):
            fp=float(4*t/h**2*cutoff((t/h)**2)*cutoff((t/h)**2,True))
            for sign in [-1,1]:
                d=4-center-sign*t;R2=planar+d*d;factor=1/(4*np.pi)
                Nz=-factor*d/R2**1.5
                Nzr=factor*3*d*dx/R2**2.5
                Nzz=factor*(2*d*d-planar)/R2**2.5
                Nzzz=factor*3*d*(3*planar-2*d*d)/R2**3.5
                Nzrr=factor*3*d*(R2-5*dx*dx)/R2**3.5
                total-=wt*sign*fp*np.array([np.sum(weight*K) for K in [Nz,Nzr,Nzz,Nzzz,Nzrr]])
        contributions[str(center)]=dict(zip(labels,map(float,total)))
    total={key:sum(row[key] for row in contributions.values()) for key in labels}
    C=float(cutoff((r/a)**2));Cr=float(cutoff((r/a)**2,True)*2*r/a**2)
    centrifugal=r*C*C
    p_r=centrifugal+total['delta_r']
    p_rr=C*C+2*r*C*Cr+total['delta_rr']
    return dict(status='EXPLORATORY whole-space endcap quadrature, not an error enclosure',
        n=n,xi_stationary=xi,r=r,unit_swirl_speed=r*C,
        unit_centrifugal=centrifugal,p_v_r=p_r,p_v_rr=p_rr,p_v_zz=total['delta_zz'],p_v_z=total['delta_z'],
        radial_acceleration=centrifugal-p_r,vertical_acceleration=-total['delta_z'],
        radial_acceleration_r=-total['delta_rr'],vertical_acceleration_z=-total['delta_zz'],
        fractional_centrifugal_remainder=(centrifugal-p_r)/centrifugal,
        harmonic_delta_trace_residual=total['delta_rr']+total['delta_r']/r+total['delta_zz'],
        endcap_contributions=contributions)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=64);p.add_argument('--r',type=float)
    args=p.parse_args();print(json.dumps(run(args.n,args.r),indent=2))

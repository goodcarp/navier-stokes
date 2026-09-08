#!/usr/bin/env python3
"""Independent 1D pressure-derivative quadrature for the flat core cutoff.

Analytic angular integration and pressure Green functions. Floating point
quadrature only: no interval certification and no PDE evolution.
"""
import json
import numpy as np
from scipy.integrate import cumulative_simpson,simpson
from scipy.special import expit

def cutoff(s):
 a=(s-.25)/.75
 inside=(a>0)&(a<1)
 t=np.clip(a,1.e-6,1-1.e-6)
 z=1/t-1/(1-t)
 z1=-t**-2-(1-t)**-2
 z2=2*t**-3-2*(1-t)**-3
 z3=-6*t**-4-6*(1-t)**-4
 p=expit(z)
 p1=p*(1-p)*z1/.75
 p2=p*(1-p)*((1-2*p)*z1*z1+z2)/.75**2
 p3=p*(1-p)*((1-6*p+6*p*p)*z1**3+3*(1-2*p)*z1*z2+z3)/.75**3
 return np.where(a<=0,1.,np.where(a>=1,0.,p)),*(np.where(inside,q,0.) for q in (p1,p2,p3))

def evaluate(n):
 # The only transition is [.25,1]; inner pressure contribution is analytic.
 s=np.linspace(.25,1.,n)
 p,p1,p2,p3=cutoff(s)
 cum=lambda q:cumulative_simpson(q,x=s,initial=0.)
 M=cum(s**3.5*p1*p1)
 N=cum(s**3.5*p*p1)
 V=cum(s**4.5*p1*p1)
 cumT=cum(p1*p1);T=cumT[-1]-cumT
 local=(-824*p*p*p1/35-3664*p*p*p2*s/245-64*p*p*p3*s*s/35
        -2848*p*p1*p1*s/105-1216*p*p1*p2*s*s/105
        +128*p*p1*p3*s**3/105+128*p*p2*p2*s**3/105
        +544*p1**3*s*s/105-64*p1*p1*p2*s**3/105)
 nested=(-128*(3*p1-2*p2*s)*M/(735*s**2.5)
         +3072*p2*N/(245*s**2.5)+1024*p2*V/(105*s**2.5)
         -256*(63*p+96*p1*s+20*p2*s*s)*T/735)
 inner=-(72/5)*(32/21)*cumT[-1]*.25
 t300_source=20/7+inner+simpson(local+nested,x=s)
 Q=.25+cum(p); K=N+7*V/9
 I=simpson(s*p*p1*p1,x=s);J=simpson(s*s*p1**3,x=s)
 t300=(340/49-18368*I/735-448*J/105
       -256*simpson(Q*p1*p1,x=s)/105
       +1536*simpson(s**(-3.5)*p1*K,x=s)/49)
 assert abs(t300-t300_source)<1.e-7
 f=p+2*s*p1/3
 twb=44/35-16*simpson(p1*f*f,x=s)/5
 return dict(n=n,t300=float(t300),tWb=float(twb),dS=0,dW=0,core_b1_Omega1=float(t300+twb),source_vs_reduced=float(t300_source-t300))

if __name__=='__main__':
 print(json.dumps([evaluate(n) for n in (2001,4001,8001,16001)],indent=2))

#!/usr/bin/env python3
"""x6 -- how much of the E_tail damage is repairable?  Three tail models, all with the
theorem's own ball geometry (H3) and the exact affine Gaussian:
  (T-a) the theorem AS WRITTEN:  2N*P + P + sum_k r0^-k sqrt(E|Z|^{2k}) sqrt(P),
        P = 2 P(|Z|>d/2) + P(|Z|>=d)
  (T-b) Cauchy-Schwarz repaired (exact truncated moments instead of sqrt(P)), still d/2
  (T-c) affine-only (no coupling step, so no d/2 and no doubling): 2N P(|Z|>=d) + ...
  (T-d) the seat's sec-5 model: A * (1/2) exp(-d^2/(2 Var_normal)), a ONE-DIMENSIONAL
        half-space crossing at threshold d (not any term of Theorem V.4)
and the smallest L at which each is <= 1/L at d = rho0 sin delta."""
import json, math
import numpy as np
from scipy.stats import chi2
from scipy.special import roots_hermitenorm
R={}
TH=4*(1-math.sqrt(2/3)); I2=4/5-16*math.sqrt(6)/135; I4=-4/7+27*math.sqrt(6)/28
sd=math.sin(math.radians(7.5)); phi0=math.radians(30.0)
_hx,_hw=roots_hermitenorm(600); _hw=_hw/_hw.sum()
def P_gt(a,sy,sz):
    if a<=0: return 1.0
    g=math.sqrt(2*sz)*_hx; rem=a*a-g*g
    return float((_hw*np.where(rem<=0,1.0,chi2.sf(np.maximum(rem,0.0)/(2*sy),df=4))).sum())
def trunc_mom(k,a,sy,sz,nq=4000):
    """E[|Z|^k 1_{|Z|>=a}] by the same chi_4 x normal factorisation."""
    g=math.sqrt(2*sz)*_hx; out=0.0
    for gi,wi in zip(g,_hw):
        rem=a*a-gi*gi
        lo=max(rem,0.0)/(2*sy)
        # integrate over chi2_4 variable t from lo to inf of (2 sy t + gi^2)^{k/2} f_4(t)
        tmax=lo+80.0
        t=np.linspace(lo,tmax,nq)
        f=0.25*t*np.exp(-t/2)                       # chi2_4 density
        out+=wi*float(np.trapz((2*sy*t+gi*gi)**(k/2.0)*f,t))
    return out
def tails(L,f=0.0,d=None):
    nu=sd**2; r0=(1+f)*math.sin(phi0); sy=nu*I2/L; sz=nu*I4/L; N=r0/sd
    if d is None: d=sd
    sC=sy/r0**2
    trC=8*sy+2*sz; trC2=4*(2*sy)**2+(2*sz)**2; m2=trC; m4=trC**2+2*trC2
    g=math.sqrt(2*sz)*_hx; a2=2*sy
    m6=float((_hw*(a2**3*192+3*a2**2*24*g*g+3*a2*4*g**4+g**6)).sum())
    Ph=P_gt(d/2,sy,sz); Pf=P_gt(d,sy,sz); P=2*Ph+Pf
    Ta=2*N*P+P+sum(r0**(-k)*math.sqrt([m2,m4,m6][k-1])*math.sqrt(P) for k in (1,2,3))
    Tb=2*N*P+P+sum(r0**(-k)*trunc_mom(k,d/2,sy,sz) for k in (1,2,3))
    Tc=2*N*Pf+Pf+sum(r0**(-k)*trunc_mom(k,d,sy,sz) for k in (1,2,3))
    Vn=2*(sd**2)*(I2*math.sin(phi0)**2+I4*math.cos(phi0)**2)/L
    A=(1+sd)*math.sin(phi0)/sd+1.0
    Td=A*0.5*math.exp(-d*d/(2*Vn))
    return dict(L=L,d=d,s_C=sC,Ta=Ta,Tb=Tb,Tc=Tc,Td=Td,P_half=Ph,P_full=Pf)
print(" L    T-a (as written)  T-b (CS repaired)  T-c (affine only)  T-d (seat sec5)   1/L      eps_bulk")
rows=[]
for L in [10,20,40,80,160,320,640]:
    t=tails(L); rows.append(t)
    print(f"{L:5d}   {t['Ta']:.4e}       {t['Tb']:.4e}        {t['Tc']:.4e}       {t['Td']:.4e}"
          f"    {1.0/L:.5f}  {t['s_C']:.4e}")
R['tail_models']=rows
def Lmin(key,target='1/L'):
    lo,hi=3.0,3e5
    for _ in range(120):
        mid=0.5*(lo+hi); t=tails(mid)
        tg=(1.0/mid) if target=='1/L' else t['s_C']
        if t[key]<=tg: hi=mid
        else: lo=mid
    return hi
for k,name in [('Ta','theorem as written'),('Tb','Cauchy-Schwarz repaired'),
               ('Tc','affine-only (no coupling)'),('Td','seat sec-5 1-D model')]:
    a=Lmin(k,'1/L'); b=Lmin(k,'eps')
    R[f'Lmin_{k}']=dict(vs_1L=a, vs_eps_bulk=b)
    print(f"  L_min ({name:26s}) vs 1/L: {a:9.3f}   vs eps_bulk: {b:9.3f}")
json.dump(R,open('x6_results.json','w'),indent=1)
print("WROTE x6_results.json")

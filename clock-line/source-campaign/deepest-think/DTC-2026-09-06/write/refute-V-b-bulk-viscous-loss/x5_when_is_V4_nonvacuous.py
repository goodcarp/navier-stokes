#!/usr/bin/env python3
"""x5 -- at what L does Theorem V.4 first say anything?  Optimise the inset f and the
polar angle phi0 to MINIMISE the theorem's own total error E_hess+E_4+E_tail, subject to
(H3) d < r0, and ask when that minimum drops below (i) 1 (non-vacuous), (ii) 1/L (the
campaign's budget), (iii) s_C = eps_bulk (the seat's own stated standard)."""
import json, math
import numpy as np
from scipy.stats import chi2
from scipy.special import roots_hermitenorm
R={}
TH=4*(1-math.sqrt(2/3)); I2=4/5-16*math.sqrt(6)/135; I4=-4/7+27*math.sqrt(6)/28
Cw=math.sqrt(1.5)*TH; n=5; sd=math.sin(math.radians(7.5))
_hx,_hw=roots_hermitenorm(300); _hw=_hw/_hw.sum()
def P_gt(a,sy,sz):
    if a<=0: return 1.0
    g=math.sqrt(2*sz)*_hx; rem=a*a-g*g
    return float((_hw*np.where(rem<=0,1.0,chi2.sf(np.maximum(rem,0.0)/(2*sy),df=4))).sum())
def bound(L,f,phi0,K2hat):
    nu=sd**2; tau=TH/L; r0=(1+f)*math.sin(math.radians(phi0))
    sy=nu*I2/L; sz=nu*I4/L; sC=sy/r0**2; N=r0/sd
    d=f
    if d>=r0: return None
    Rm=r0-d
    V1=n*nu*tau*(math.exp(2*Cw)-1)/Cw; EW=0.5*K2hat*math.exp(Cw)*tau*V1
    trC=8*sy+2*sz; trC2=4*(2*sy)**2+(2*sz)**2; m2=trC; m4=trC**2+2*trC2
    g=math.sqrt(2*sz)*_hx; a2=2*sy
    m6=float((_hw*(a2**3*192+3*a2**2*24*g*g+3*a2*4*g**4+g**6)).sum())
    P=2*P_gt(d/2,sy,sz)+P_gt(d,sy,sz)
    Eh=(r0/Rm**2)*EW+4*N*EW/d
    E4=(r0/Rm**5)*m4
    Et=2*N*P+P+sum(r0**(-k)*math.sqrt([m2,m4,m6][k-1])*math.sqrt(P) for k in (1,2,3))
    return dict(s_C=sC,E_hess=Eh,E_4=E4,E_tail=Et,total=Eh+E4+Et,f=f,phi0=phi0,d=d,r0=r0)
def best(L,K2hat):
    bb=None
    for phi0 in [30.,45.,60.,75.]:
        sp0=math.sin(math.radians(phi0)); fmax=0.99*sp0/(1-sp0)
        for f in np.geomspace(1e-3,min(fmax,80.0),260):
            b=bound(L,f,phi0,K2hat)
            if b and (bb is None or b['total']<bb['total']): bb=b
    return bb
print(" L      K2hat=0: min total   (f,phi0)      |  K2hat=1: min total   (f,phi0)     | eps_bulk(f=0,30d) | 1/L")
rows=[]
for L in [10,20,40,80,160,320,640,1280,2560,5120]:
    b0=best(L,0.0); b1=best(L,1.0)
    eb=sd**2*I2/L/(0.5**2)
    rows.append(dict(L=L,min_total_K0=b0['total'],f0=b0['f'],phi00=b0['phi0'],sC0=b0['s_C'],
                     min_total_K1=b1['total'],f1=b1['f'],phi01=b1['phi0'],sC1=b1['s_C'],
                     eps_bulk_f0=eb,target=1.0/L))
    print(f"{L:5d}   {b0['total']:.4e}  ({b0['f']:.3f},{b0['phi0']:.0f})   |  "
          f"{b1['total']:.4e}  ({b1['f']:.3f},{b1['phi0']:.0f})  |  {eb:.4e}  |  {1.0/L:.5f}")
R['minimised_total']=rows
grid=[10,14,20,28,40,56,80,113,160,226,320,453,640,905,1280,1810,2560,3620,5120,7240,10240,20480,40960,81920]
for tname,tgt in [('<1 (non-vacuous)',lambda L,b:1.0),('<1/L',lambda L,b:1.0/L),
                  ('< s_C (eps_bulk)',lambda L,b:b['s_C'])]:
    for K2hat in [0.0,1.0]:
        hit=None
        for L in grid:
            b=best(L,K2hat)
            if b['total']<=tgt(L,b): hit=L; break
        R[f'L_min_{tname}_K2hat{K2hat}']=hit
        print(f"  smallest grid L with min-total {tname:18s} at K2hat={K2hat}:  L = {hit}")
json.dump(R,open('x5_results.json','w'),indent=1)
print("WROTE x5_results.json")

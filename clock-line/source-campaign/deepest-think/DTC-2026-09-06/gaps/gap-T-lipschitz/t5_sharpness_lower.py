"""
t5_sharpness_lower.py -- how far LEMMA T's constant is from the truth.

Witness family (all hypotheses of LEMMA T satisfied):  Phi = T_{lambda'} , lambda' = lambda(1+s).
By Lemma 1 of lower/prove-lagrangian (re-verified in t1/t3):  a[eta0 o T_{l}^{-1}](0) = (M/2) l L,
so   Delta = (M/2)(lambda'-lambda) L   EXACTLY, while
     mu   = sup |T_{lambda'}x - T_lambda x|/|x| = max(|lambda'-lambda|, |lambda'^-2 - lambda^-2|),
     muJ  = |lambda'^2/lambda^2 - 1|.
Hence the BEST possible constant in |Delta| <= C mu M L is >= 1/2 whenever lambda^3 >= 2.
This file computes the witness ratio numerically (own quadrature, not the identity) and
compares with the proved C(lambda,mu).
"""
import json, numpy as np, mpmath as mp
mp.mp.dps=40
M,RHO0,LL=1.0,1.0,6.0

def I1(lam,mu):
    if mu>=min(lam,lam**-2): return float('inf')
    f=lambda p: mp.sin(p)**2*(mp.sqrt(lam**2*mp.sin(p)**2+lam**-4*mp.cos(p)**2)-mu)**-5
    return float(mp.quad(f,[0,mp.pi/2,mp.pi]))

def a_T(lam,delta=0.0,n=200):
    """a[eta0 o T_lam^{-1}](0) by direct quadrature (independent of Lemma 1)."""
    x,w=np.polynomial.legendre.leggauss(n); tot=0.0
    brk=[0.0]+([delta,np.pi-delta] if delta>0 else [])+[np.pi/2,np.pi]
    brk=sorted(set(brk))
    for a_,b_ in zip(brk[:-1],brk[1:]):
        p=0.5*(b_-a_)*x+0.5*(a_+b_); ww=0.5*(b_-a_)*w
        q=np.sqrt(lam**2*np.sin(p)**2+lam**-4*np.cos(p)**2)
        h=np.ones_like(p) if delta<=0 else np.minimum(1.0,np.minimum(p,np.pi-p)/delta)
        tot+=np.sum(ww*np.abs(np.cos(p))*np.sin(p)**2*q**-5*h)
    return 0.75*M*LL*tot

rows=[]
for lam in (1.0,1.25,1.5,2.0,3.0):
    for s in (0.02,0.01,0.005):
        lp=lam*(1+s)
        mu=max(abs(lp-lam),abs(lp**-2-lam**-2)); muJ=abs(lp**2/lam**2-1)
        D=abs(a_T(lp)-a_T(lam))
        rows.append(dict(lam=lam,s=s,lam_prime=lp,mu=mu,muJ=muJ,
                         Delta_quadrature=D, Delta_lemma1=(M/2)*(lp-lam)*LL,
                         true_ratio_Delta_over_muML=D/(mu*M*LL),
                         proved_C=3*(1+muJ)*lam**2*I1(lam,mu)+lam/2*(muJ/mu),
                         hypotheses_muJ_le_mu=bool(muJ<=mu)))
out=dict(witness=rows,
         note='ratio -> 1/2 exactly once lambda^3 >= 2 (then mu = lambda s); '
              'for lambda^3 < 2 the z-component |lambda\'^-2-lambda^-2| ~ 2s/lambda^2 sets mu '
              'and the ratio -> lambda^3/4.',
         best_possible_C_lower_bound=max(r['true_ratio_Delta_over_muML'] for r in rows),
         proved_C_at_mu0={str(l):3*l**2*I1(l,0)+l/2 for l in (1.0,1.25,1.5,2.0,3.0)})
print(json.dumps(out,indent=1))
json.dump(out,open('t5_results.json','w'),indent=1)

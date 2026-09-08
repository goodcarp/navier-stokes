#!/usr/bin/env python3
"""S6 (unregistered follow-up, reported as an additive find, not as a gate):
the inner-edge plateau offset c_edge(phi) is NOT a universal constant.
For a scale-invariant shell profile w(phi) with log rate kappa = -(3/5) g_1,
   a(rho0+, phi) = kappa log(R/rho0) + c_edge(phi),
   c_edge(phi) = - SUM_{l>=2} ((l+2) g_l /((2l+3)(l-1))) C_{l-1}^{3/2}(cos phi).
(1) bang-bang w = -M sgn(z):   c_edge diverges like A log(1/phi) as phi -> 0.  Fit A.
(2) admissible taper w = -M sgn(z) min(1, phi_ax/delta):  does c_edge saturate?"""
import numpy as np, math, json, hashlib, sys, time
sys.path.insert(0,'.')
from kern import geg_table, Nl
t0=time.time(); LMAX=1001
xq, wq = np.polynomial.legendre.leggauss(2500)
TH = 0.5*math.pi*(xq+1.0); WQ = 0.5*math.pi*wq
CT = geg_table(np.cos(TH), LMAX)
def g_of(wfun):
    return (CT @ (WQ*np.sin(TH)**2*wfun(TH)))/Nl(np.arange(LMAX+1))
def c_edge(g, phi):
    C = geg_table(np.array([math.cos(phi)]), LMAX+1)[:,0]
    ls = np.arange(2, LMAX+1)
    terms = -((ls+2)*g[ls]/((2*ls+3)*(ls-1)))*C[ls-1]
    P = np.cumsum(terms); return float(P[-len(P)//4:].mean())
bang = lambda ph: -np.sign(np.cos(ph))
def taper(delta):
    d = math.radians(delta)
    return lambda ph: -np.sign(np.cos(ph))*np.minimum(1.0, np.minimum(ph, math.pi-ph)/d)
FAM = [("bang-bang (inadmissible: eta unbounded on the axis)", bang),
       ("taper delta=15 deg (admissible)", taper(15)),
       ("taper delta= 7.5 deg (admissible)", taper(7.5))]
res={}
for tag, wf in FAM:
    g = g_of(wf); kap = -0.6*g[1]
    phis = [0.25,0.5,1,2,3,5,7.5,10,15,20,30,45,60,90]
    cs = [c_edge(g, math.radians(p)) for p in phis]
    print(f"\n=== {tag} ===\n   log rate kappa = -(3/5) g_1 = {kap:.6f}")
    print("   phi(deg): " + " ".join(f"{p:>8g}" for p in phis))
    print("   c_edge  : " + " ".join(f"{c:>8.4f}" for c in cs))
    sm = [(p,c) for p,c in zip(phis,cs) if p <= 5]
    A, B = np.polyfit([-math.log(math.radians(p)) for p,_ in sm], [c for _,c in sm], 1)
    pred = [A*(-math.log(math.radians(p)))+B for p,_ in sm]
    r = max(abs(pred[i]-sm[i][1]) for i in range(len(sm)))
    print(f"   fit over phi<=5 deg:  c_edge ~ {A:.5f} * log(1/phi) + {B:.5f}   (max residual {r:.1e})")
    res[tag]=dict(kappa=kap, phis=phis, c_edge=cs, logfit_A=A, logfit_B=B, maxres=r)
print("\nreading: A ~ 1/4 for the bang-bang cap (unbounded), A -> 0 for admissible tapers (saturation).")
json.dump(res, open('s6_results.json','w'), indent=1, default=float)
print(f"elapsed {time.time()-t0:.1f}s  SCRIPT-SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")

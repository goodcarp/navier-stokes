#!/usr/bin/env python3
"""G17 — octave model of the multi-scale stack's self-stretching (frozen positions; prediction for the depletion numerics).
Shells k = 0..N-1 at radii rho_k = rho0 2^k, vorticity w_k(t), w_k(0) = M. Exact facts used: an outer shell at rho_j contributes
(ln 2)/2 * w_j of uniform strain at every point inside it (kappa0 = 1/2 per e-fold, one octave = ln 2 e-folds); a shell's own/near-field
strain on itself is c_self * w_k with c_self = 0.2168 (independent exact constant from the sharp-clock refuter, +0.216773 M on the material).
Inviscid, no swirl, positions frozen (shells move outward slowly on the window): d w_k/dt = a_k w_k, a_k = (ln2/2) sum_{j>k} w_j + c_self w_k.
Outputs: T_{3/2}, T_2, T_4, T_8 (first time max_k w_k reaches x M) vs N, the products with M and with log(R/rho0) = N ln 2,
and the coherence count C(t) = a_0(t)/max_k w_k(t) (axis strain per unit cap; a_0 = (ln2/2) sum_j w_j).
PREDICTION registered before the depletion numerics report: T_{3/2}*M*N ln2 ~ const (log-fast first doubling); T_4*M -> O(1) constant
in N (the log is a stock); C(t) collapses to O(1) within ~1 turnover. Numerics falsify; this model is not the PDE."""
import numpy as np, sys, hashlib
from scipy.integrate import solve_ivp
M=1.0; ln2=np.log(2); cs=0.2168
def run(N, targets=(1.5,2,4,8), tmax=4.0):
    def rhs(t,w):
        S=np.cumsum(w[::-1])[::-1]-w      # sum_{j>k} w_j
        a=(ln2/2)*S+cs*w
        return a*w
    hit={}
    def ev(x):
        f=lambda t,w: np.max(w)-x*M; f.terminal=False; return f
    sol=solve_ivp(rhs,(0,tmax),np.full(N,M),rtol=1e-10,atol=1e-12,dense_output=True,max_step=1e-3)
    ts=np.linspace(0,tmax,40001); W=sol.sol(ts); mx=W.max(axis=0)
    for x in targets:
        i=np.argmax(mx>=x*M); hit[x]=ts[i] if mx[i]>=x*M else np.nan
    a0=(ln2/2)*W.sum(axis=0); C=a0/mx
    Cpts={tt: float(np.interp(tt,ts,C)) for tt in (0,0.25,0.5,1,1.5,2)}
    return hit,Cpts
print(f"{'N':>3} {'log(R/r0)':>9} {'T1.5*M':>8} {'T2*M':>8} {'T4*M':>8} {'T8*M':>8} {'T1.5*M*log':>11} {'T4*M*log':>9}   C(0) C(.5) C(1) C(1.5) C(2)")
rows=[]
for N in (2,3,4,5,6,7,8,10,12):
    h,C=run(N); L=N*ln2
    rows.append((N,h))
    print(f"{N:3d} {L:9.3f} {h[1.5]:8.4f} {h[2]:8.4f} {h[4]:8.4f} {h[8]:8.4f} {h[1.5]*L:11.4f} {h[4]*L:9.4f}   {C[0]:.2f} {C[0.5]:.2f} {C[1]:.2f} {C[1.5]:.2f} {C[2]:.2f}")
t15=[r[1][1.5] for r in rows]; t4=[r[1][4] for r in rows]
print(f"\nspread T1.5*M over N=3..12: x{max(t15[1:])/min(t15[1:]):.2f};  spread T4*M over N=3..12: x{max(t4[1:])/min(t4[1:]):.2f}")
print("PREDICTION: depletion holds in the model iff T4*M spread << T1.5*M spread and C(t) -> O(1) within a turnover")
ok = (max(t4[1:])/min(t4[1:]) < 1.5) and (max(t15[1:])/min(t15[1:]) > 2.0)
print('MODEL VERDICT', 'depletion' if ok else 'no depletion', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())

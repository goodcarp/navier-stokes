"""
r3_accumulation_quadrature.py -- REFUTER
Re-derive the loglog accumulation constant of the attempt's SS4 by NUMERICAL QUADRATURE in
v = log(1/r), r = T-s, avoiding the underflow that killed r2 block (4), and WITHOUT using
the symbolic antiderivative the attempt verified.  Also check the log_+ branch caveat:
the antiderivative F(r) = -5c log(1+log K-(1/5)log r) is undefined where 1+log K-(1/5)log r<=0,
which happens for small K -- does the attempt's restriction to the log_+-active branch cover it?
"""
import numpy as np
from scipy.integrate import quad

def bracket(v, K):                       # 1 + log_+( K r^{-1/5} ),  v = log(1/r)
    return 1.0 + max(np.log(K) + v/5.0, 0.0)

print("Integrand in v:  1/bracket(v,K).   int_0^t ||omega|| ds / c_1 = int_{V0}^{V} dv/bracket.")
print("(dr/r = -dv, so the 1/r cancels exactly; V=log(1/(T-t)), V0=log(1/T).)")
print()
for K in [1.0, 1e2, 1e-3]:
    v_act = max(0.0, -5.0*np.log(K))     # log_+ becomes active at v = -5 log K
    print(f"K = {K:g}   (log_+ active for v > {v_act:.3f};  1+logK = {1+np.log(K):+.3f})")
    for V in [10, 30, 100, 1e3, 1e5, 1e8]:
        I = quad(lambda v: 1.0/bracket(v, K), 0.0, V, limit=500)[0]
        LL = np.log(V)                   # loglog(1/(T-t)) with T=1
        print(f"   V=log(1/(T-t))={V:9.0e}   int/c1={I:12.4f}   loglog={LL:8.4f}   ratio={I/LL:8.4f}")
    print()
print("The ratio must tend to 5 (the claimed coefficient 5c_1).  It does, but only ~ 1/log V,")
print("so the o(1) in `(5c_1+o(1))` is a 1/loglog correction, not a small number at any")
print("physically reachable (T-t):  at (T-t)=1e-100 (V=230) the ratio is still ~")
I=quad(lambda v: 1.0/bracket(v,1.0),0.0,230.0,limit=500)[0]
print(f"   {I/np.log(230.0):.4f}  (vs 5).")
print()
print("BRANCH CHECK: for K<e^{-1} the expression 1+logK-(1/5)log r is NEGATIVE on part of")
print("the range, so F(r)=-5c log(.) is undefined there; the attempt restricts to the")
print("log_+-ACTIVE branch, where the bracket is 1+(logK+v/5) > 1 > 0 by definition of")
print("activation.  Verified numerically:")
for K in [1e-3, 1e-8]:
    v_act = -5.0*np.log(K)
    for dv in [1e-9, 1e-3, 1.0]:
        v = v_act + dv
        print(f"   K={K:g}  v=v_act+{dv:g}:  1+logK+v/5 = {1+np.log(K)+v/5:.9f}  (must be >1)")
print("So SS4's restriction is sound; the additive C(E_0,nu,T) absorbs the inactive branch.")

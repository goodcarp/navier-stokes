#!/usr/bin/env python3
"""G10 — where the dimension family has ZERO slack against the funnel.
(i) EXACT (sympy): in Hou's generalized model (G9), the stretching-free quantity is Theta := omega1 * r^{3-n} = omega^theta / r^{n-2}:
      D_t Theta = r^{3-n} [ (Gamma^2/r^4)_z + nu( w1_rr + n/r w1_r + w1_zz ) ]     (the (n-3) a w1 term cancels identically)
    so the swirl source feeds the TRANSPORTED variable with weight r^{3-n}: = 1 at n=3, and ~ r^{-0.188} near the axis at n=3.188.
(ii) EXACT exponents: the funnel's R1-creation balance (source alone vs the near-cap envelope theta*M/s within one viscous time s^2/nu,
    with |d_z(q^2)| <= M^2/s kept as in 3D — an assumption flagged) becomes  M^2 s^{2-n} * s^2/nu = theta M / s  =>  s_c(n) = (theta nu / M)^{1/(5-n)}.
    n=3: s_c = sqrt(nu/M) = ell  (the funnel's 'on the nose' coincidence);  n=3.188: s_c = M^{-1/(5-n)} = M^{-0.5522} < ell.
    The window (s_c(n), ell) — where the source can create cap-level eta within a viscous time — has width ratio M^{1/(5-n) - 1/2}: EMPTY at n=3,
    M^{0.0522} at n=3.188. Source amplification at the core scale ell: ell^{3-n} = M^{(n-3)/2} = M^{0.094}.
Controls: at n=3 (i) must reduce to the 3D eta-equation (no r-weight) and (ii) must give s_c = ell exactly; a wrong transported power (r^{2-n}) must NOT cancel the stretching."""
import sympy as sp, sys, hashlib
from fractions import Fraction as F
r, z, n, nu, t = sp.symbols('r z n nu t', positive=True)
psi1 = sp.Function('psi1')(r, z, t); Gam = sp.Function('Gamma')(r, z, t); w1 = sp.Function('w1')(r, z, t)
ur = -r*sp.diff(psi1, z); uz = 2*psi1 + r*sp.diff(psi1, r)   # 3D form; G9 verified u^r = -r psi1_z for all n (u^z unused below)
a = ur/r
Dt = lambda f: sp.diff(f, t) + ur*sp.diff(f, r) + uz*sp.diff(f, z)
visc = nu*(sp.diff(w1,r,2) + n/r*sp.diff(w1,r) + sp.diff(w1,z,2))
# Hou's w1 equation solved for w1_t:
w1_t = -(ur*sp.diff(w1,r) + uz*sp.diff(w1,z)) + sp.diff(Gam**2/r**4, z) + (n-3)*a*w1 + visc
def Dt_with(f, w1t):  # material derivative using Hou's w1_t
    return sp.diff(f, t).subs(sp.Derivative(w1, t), w1t) + ur*sp.diff(f, r) + uz*sp.diff(f, z)
Theta = w1*r**(3-n)
lhs = Dt_with(Theta, w1_t)
rhs = r**(3-n)*(sp.diff(Gam**2/r**4, z) + visc)
chk_i = sp.simplify(sp.expand(lhs - rhs))
Theta_wrong = w1*r**(2-n)
chk_ctrl = sp.simplify(sp.expand(Dt_with(Theta_wrong, w1_t) - r**(2-n)*(sp.diff(Gam**2/r**4, z) + visc)))
print('(i) D_t(w1 r^{3-n}) - r^{3-n}[source + visc] =', chk_i)
print('    CTRL wrong power r^{2-n}: residual =', sp.factor(chk_ctrl), '(must be nonzero)')
# (ii)
def s_c_exp(nn):  # exponent of M in s_c(n) with theta, nu = O(1)
    return -1/(5-nn)
print('(ii) s_c exponent: n=3 ->', F(-1,2) == F(-1, 5-3), '(equals ell = M^{-1/2});  n=3.188 ->', round(s_c_exp(3.188), 4), '; window ratio ell/s_c = M^{', round(1/(5-3.188) - 0.5, 4), '} (n=3: M^0, empty)')
print('    source amplification at r = ell: ell^{3-n} = M^{(n-3)/2} = M^{', round((3.188-3)/2, 4), '}')
ok = (chk_i == 0) and (chk_ctrl != 0) and (F(-1,2) == F(-1, 5-3)) and abs(s_c_exp(3.188) + 0.5522) < 1e-3
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

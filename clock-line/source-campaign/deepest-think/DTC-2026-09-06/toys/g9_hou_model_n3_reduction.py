#!/usr/bin/env python3
"""G9 — Hou's generalized axisymmetric system (arXiv 2405.10916, eqs (1.1a-c) as fetched from the html) must reduce at n=3 to the
axisymmetric NS identities verified exactly in G1; and the dimension enters the eta-equation as a STRETCHING term (n-3) a eta.
Hou (transcribed):  Gamma_t + u^r Gamma_r + u^z Gamma_z = nu( Gamma_rr + (n-4)/r Gamma_r + (6-2n)/r^2 Gamma + Gamma_zz )
                    w1_t + u^r w1_r + u^z w1_z = (Gamma^2/r^4)_z - (n-3) psi1_z w1 + nu( w1_rr + n/r w1_r + w1_zz )
                    -(d_rr + n/r d_r + d_zz) psi1 = w1 ;  u^r = -(r^{n-2} psi^theta)_z / r^{n-2},  u^z = (r^{n-2} psi^theta)_r / r^{n-2},  psi^theta = r psi1.
Checks: (a) u^r = -r psi1_z for every n (so -(n-3) psi1_z w1 = +(n-3) a w1 with a = u^r/r);  (b) at n=3 the Gamma operator equals
nu(Lap - (2/r) d_r) [G1 SWIRL-Gam] and the w1 operator equals nu L5 with source d_z(q^2), q = Gamma/r^2 [G1 ETA-EQ]; (c) the (6-2n)/r^2 term is
<= 0 for n>3 (the Gamma maximum principle survives: G = O(1) is n-robust); (d) exponent bookkeeping: over a record interval with transport
charge Xi = ∫a_+ dt = (1/10) log M the extra factor on eta is exp((n-3) Xi) = M^{(n-3)/10} (= M^{0.0188} at 3.188)."""
import sympy as sp, sys, hashlib
r, z, n, nu = sp.symbols('r z n nu', positive=True)
psi1 = sp.Function('psi1')(r, z); Gam = sp.Function('Gamma')(r, z); w1 = sp.Function('w1')(r, z)
psith = r*psi1
ur = -sp.diff(r**(n-2)*psith, z)/r**(n-2); uz = sp.diff(r**(n-2)*psith, r)/r**(n-2)
chk_a = sp.simplify(ur + r*sp.diff(psi1, z))            # (a)
gam_op = lambda nn: sp.diff(Gam,r,2) + (nn-4)/r*sp.diff(Gam,r) + (6-2*nn)/r**2*Gam + sp.diff(Gam,z,2)
Lap = sp.diff(Gam,r,2) + sp.diff(Gam,r)/r + sp.diff(Gam,z,2)
chk_b1 = sp.simplify(gam_op(3) - (Lap - 2*sp.diff(Gam,r)/r))                                  # (b) Gamma at n=3
w_op = lambda nn: sp.diff(w1,r,2) + nn/r*sp.diff(w1,r) + sp.diff(w1,z,2)
L5 = sp.diff(w1,r,2) + 3*sp.diff(w1,r)/r + sp.diff(w1,z,2)
chk_b2 = sp.simplify(w_op(3) - L5)
q = Gam/r**2
chk_b3 = sp.simplify(sp.diff(Gam**2/r**4, z) - sp.diff(q**2, z))                              # source = d_z(q^2)
stretch = sp.simplify(-(n-3)*sp.diff(psi1,z)*w1 - (n-3)*(ur/r)*w1)                            # equals (n-3) a w1
coef = sp.Rational(6,1) - 2*sp.Rational(3188,1000)                                            # (c) at n=3.188
# (d)
extra_exp = (sp.Rational(3188,1000)-3)*sp.Rational(1,10)
print('(a) u^r + r psi1_z =', chk_a, '| (b) Gamma op n=3 residual =', chk_b1, '; w1 op n=3 residual =', chk_b2, '; source residual =', chk_b3)
print('    stretching term identity -(n-3)psi1_z w1 == (n-3) a w1 :', stretch == 0)
print('(c) (6-2n) at n=3.188 =', coef, '(<0: Gamma damped near the axis; max principle survives)')
print('(d) extra eta factor over one transport charge: M^{(n-3)/10} = M^{', extra_exp, '} at n=3.188')
# control: a mis-transcription (n-2)/r in the w1 Laplacian would NOT reduce to L5 at n=3
ctrl = sp.simplify((sp.diff(w1,r,2) + (3-2)/r*sp.diff(w1,r) + sp.diff(w1,z,2)) - L5)
print('CTRL wrong coefficient (n-2)/r at n=3 residual =', ctrl, '(must be nonzero)')
ok = (chk_a == 0) and (chk_b1 == 0) and (chk_b2 == 0) and (chk_b3 == 0) and (stretch == 0) and (coef < 0) and (ctrl != 0)
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

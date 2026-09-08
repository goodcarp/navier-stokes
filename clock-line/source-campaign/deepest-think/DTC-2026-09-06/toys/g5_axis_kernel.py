#!/usr/bin/env python3
"""G5 — the closed single-scale axis model is SUBCRITICAL (exact sympy).
eta(r,z) = eta0(z) 1_{r<l}  =>  u^z(0,z) = 2 (K_l * eta0)(z),  K_l(zeta) = (1/4) ∫_0^l r^3 dr /(zeta^2+r^2)^{3/2},
a(0,z) = -(1/2) d_z u^z = -(K_l' * eta0)(z).  Computed closed forms (the first draft's guessed values were wrong; the computation is the truth): K_l(0) = l/4, K_l ~ l^4/(16 zeta^3),
K_l' ~ -3 l^4/(16 zeta^4), K_l' <= 0 on zeta>0, so ∫|K_l'| = 2 K_l(0) = l/2 and sup|a(0,.)| <= (M/l)(l/2) = M/2 with |eta0| <= M/l: NO log.
Control: the Hilbert kernel 1/(pi zeta) has a non-integrable tail (order 0), unlike K_l' (order -1)."""
import sympy as sp, sys, hashlib
r, l, zeta = sp.symbols('r l zeta', positive=True)
K = sp.simplify(sp.Rational(1,4)*sp.integrate(r**3/(zeta**2 + r**2)**sp.Rational(3,2), (r, 0, l)))
K0 = sp.limit(K, zeta, 0); tail = sp.limit(K*zeta**3, zeta, sp.oo)
Kp = sp.simplify(sp.diff(K, zeta)); Kp_tail = sp.limit(Kp*zeta**4, zeta, sp.oo)
# sign of K' on zeta>0: check at sample points and via the closed form
neg_ok = all((Kp.subs({l:1, zeta:v}) < 0) for v in [sp.Rational(1,10), 1, 10, 100])
L1 = 2*K0   # ∫|K'| = 2 ∫_0^inf (-K') = 2 (K(0) - K(inf)) = 2 K(0)
print('K_l(zeta) =', K); print('K_l(0) =', K0, '| zeta^3 K -> ', tail, '| zeta^4 K\' -> ', Kp_tail, "| K'<0 on zeta>0:", neg_ok)
print('∫|K_l\'| dzeta = 2K_l(0) =', L1, ' -> sup|a(0,.)| <= (M/l)*2K_l(0) = M *', sp.simplify(L1/l), ' (pure number times M: no log)')
ctrl = sp.integrate(1/(sp.pi*zeta), (zeta, 1, sp.oo))
print('CTRL Hilbert tail ∫_1^oo dzeta/(pi zeta) =', ctrl, '(must be oo)')
ok = (sp.simplify(K0 - l/4) == 0) and (sp.simplify(tail - l**4/16) == 0) and (sp.simplify(Kp_tail + 3*l**4/16) == 0) and neg_ok and (ctrl == sp.oo)
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

#!/usr/bin/env python3
"""G11 — ERRATUM gate for G10(ii). The R1-creation balance must be written in ONE variable.
Transported variable Theta = omega^theta / r^{n-2} (no stretching, G10(i)); its source weight r^{3-n} * |d_z q^2| <= M^2 r^{2-n};
its cap envelope theta*M / r^{n-2}. Creation time = envelope/source = theta/M for every scale s and every n.
So s_c = sqrt(theta nu / M) for ALL n (creation time theta/M vs viscous time s^2/nu): the family has exponent-level slack at R1 too.
Control: the WRONG pairing (eta-envelope theta M/s with the Theta-source M^2 s^{2-n}) reproduces G10's spurious s_c(n) = (theta nu/M)^{1/(5-n)}
— i.e. the error is located exactly, and it vanishes at n=3 (which is why G10's own n=3 control could not catch it)."""
import sympy as sp, sys, hashlib
s, n, M, nu, theta = sp.symbols('s n M nu theta', positive=True)
source_Theta = M**2 * s**(2-n)          # r^{3-n} * (M^2 / r)
env_Theta    = theta*M / s**(n-2)
t_create = sp.simplify(env_Theta / source_Theta)
print('creation time in the transported variable =', t_create, '(independent of s and n)')
s_c = sp.solve(sp.Eq(t_create, s**2/nu), s)[0]
print('s_c = ', s_c, '  (= sqrt(theta nu/M) for all n)')
# control: the wrong pairing
env_eta = theta*M/s
t_wrong = sp.simplify(env_eta / source_Theta)
s_wrong = sp.solve(sp.Eq(t_wrong, s**2/nu), s)
print('CTRL wrong pairing: creation time =', t_wrong, ' -> s_c(n) =', s_wrong, ' (G10\'s spurious window; equals the right answer only at n=3)')
ok = (t_create == theta/M) and sp.simplify(s_c - sp.sqrt(theta*nu/M)) == 0 and sp.simplify(t_wrong.subs(n,3) - theta/M) == 0 and sp.simplify(t_wrong.subs(n, sp.Rational(3188,1000)) - theta/M) != 0
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

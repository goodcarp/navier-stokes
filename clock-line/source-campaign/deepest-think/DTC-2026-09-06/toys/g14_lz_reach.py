#!/usr/bin/env python3
"""G14 — the reach of Lei–Zhang's L^2 method in the (G, M) plane (method coordinate, not a theorem about NS).
Facts used: (i) their proof uses delta* only through the absorption K0 delta* ||∇Ω||^2 <= (1/2)||∇Ω||^2 (their (3.6)-(3.7); K0 absolute);
(ii) under the caps, |v^theta| <= min(G/r, M r/2) (Gamma(0)=0, |omega^z| <= M, Gamma <= G), crossover r_G = sqrt(2G/M);
(iii) FBC (1.6) for f supported in r < r0 (f(r0)=0): the effective delta* is the sup of ∫|v^theta|^2 f^2 r dr / ∫ f_r^2 r dr.
On the log region r_G < r < r0 with u = ln(r0/r) in (0, L), L = ln(r0/r_G): weight G^2/r^2 gives ∫ G^2 f^2 du / ∫ f_u^2 du,
whose sup over f(u=0)=0 is the one-sided Wirtinger constant (2L/pi)^2 (extremal f = sin(pi u/(2L))).
So delta*_eff >= (4 G^2/pi^2) L^2 with L = (1/2) ln(M r0^2/(2G)) — a LOG SQUARED (matches Lei–Zhang's |ln r|^{-2} criticality and G13(a)).
Absorption fails when K0 delta*_eff >= 1/2  <=>  ln(M r0^2/(2G)) >= pi/(G sqrt(2 K0))  <=>  M >= M_*(G) := (2G/r0^2) exp(pi/(G sqrt(2K0))).
Checks: (a) Wirtinger constant on (0,L) with one Dirichlet end = (2L/pi)^2 (Rayleigh quotient of the extremal, and a non-extremal test function is smaller);
(b) M_*(G) formula from the threshold identity; (c) control: replacing (2L/pi)^2 by the linear guess L (my first reading) gives a different threshold — the log^2 is real."""
import sympy as sp, sys, hashlib
u, L, G, M, r0, K0 = sp.symbols('u L G M r0 K0', positive=True)
f_ext = sp.sin(sp.pi*u/(2*L))
R_ext = sp.simplify(sp.integrate(f_ext**2, (u, 0, L)) / sp.integrate(sp.diff(f_ext,u)**2, (u, 0, L)))
f_test = u*(2*L - u)     # vanishes at 0, zero slope at L: admissible, non-extremal
R_test = sp.simplify(sp.integrate(f_test**2, (u, 0, L)) / sp.integrate(sp.diff(f_test,u)**2, (u, 0, L)))
print('(a) Rayleigh quotient of the extremal =', R_ext, '| of a test function =', R_test, '| test <= extremal:', sp.simplify(R_ext - R_test) >= 0)
Lval = sp.Rational(1,2)*sp.log(M*r0**2/(2*G))
delta_eff = 4*G**2/sp.pi**2 * Lval**2
Mstar = sp.solve(sp.Eq(K0*delta_eff, sp.Rational(1,2)), M)
Mstar = [m for m in Mstar if sp.simplify(m) != 0]
print('(b) delta*_eff =', sp.simplify(delta_eff), ' | M_*(G) =', [sp.simplify(m) for m in Mstar])
target = (2*G/r0**2)*sp.exp(sp.pi/(G*sp.sqrt(2*K0)))
ok_b = any(sp.simplify(m - target) == 0 for m in Mstar)
# (c) control: linear guess
delta_lin = G**2*Lval
M_lin = sp.solve(sp.Eq(K0*delta_lin, sp.Rational(1,2)), M)
print('(c) CTRL linear-in-L guess would give M =', [sp.simplify(m) for m in M_lin], '(differs: exp(1/(G^2 K0)) vs exp(pi/(G sqrt(2K0))))')
ok = (sp.simplify(R_ext - (2*L/sp.pi)**2) == 0) and (sp.simplify(R_ext - R_test) >= 0) and ok_b and all(sp.simplify(m - target) != 0 for m in M_lin)
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

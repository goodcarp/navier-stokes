"""
s2_leray_nu_and_invariance.py -- DTC-2026-09-06 / lower / corollary-blowup-rate

(1) restores the viscosity power in the Leray/Giga L^q lower bound, by the exact
    ny-normalisation w(x,s) = u(x,s/nu)/nu  which solves the nu=1 equations;
(2) checks that all three forms of the corollary are invariant under the NS
    scaling u -> lam u(lam x, lam^2 t) at fixed nu.
All output is produced by this script.
"""
import sympy as sp

nu, q, s, lam, E0, M, C = sp.symbols('nu q s lambda E0 M C', positive=True)

print("="*72)
print("(1) viscosity power in the Leray/Giga lower bound")
print("="*72)
print("  u solves NS with viscosity nu on [0,T).  Put w(x,sigma) = u(x,sigma/nu)/nu.")
print("  Then w solves the nu=1 equations and blows up at Tw = nu*T.")
print("  nu=1 statement (Leray 1934 / Giga 1986):  ||w(sigma)||_q >= c (Tw-sigma)^(-(1-3/q)/2)")
alpha = (1-3/q)/2
# ||w(sigma)||_q = ||u(t)||_q / nu ,  Tw - sigma = nu (T - t)
lhs = sp.Symbol('U', positive=True)/nu           # ||u(t)||_q / nu
rhs = (nu*s)**(-alpha)                           # (nu (T-t))^(-alpha),  s = T-t
Ubound = sp.simplify(nu*rhs)
print("  => ||u(t)||_q >= c * nu * (nu (T-t))^(-(1-3/q)/2) = c *", sp.powsimp(Ubound, force=True))
for qq in [sp.oo, 4, 6, 9]:
    e_nu = sp.simplify((1 - (1-3/q)/2).subs(q, qq))
    e_s  = sp.simplify((-(1-3/q)/2).subs(q, qq))
    print(f"     q={qq}:  nu exponent = {e_nu},  (T-t) exponent = {e_s}")
print("  q=inf:  ||u(t)||_inf >= c nu^(1/2) (T-t)^(-1/2).")

print()
print("="*72)
print("(2) NS-scaling invariance of the three forms of the corollary")
print("="*72)
# u_lam(x,t) = lam u(lam x, lam^2 t):  E -> E/lam, M -> lam^2 M, (T-t) -> (T-t)/lam^2, nu fixed
sub = {E0: E0/lam, M: lam**2*M, s: s/lam**2}
forms = {
 "(a)  M*(T-t)*(1+log+ Re_E)   [Re_E = E^(2/5) M^(1/5)/nu]":
     M*s*(1+sp.log(E0**sp.Rational(2,5)*M**sp.Rational(1,5)/nu)),
 "(b)  M*(T-t)*(1+log+(E0^(2/5) M^(1/5)/nu))":
     M*s*(1+sp.log(E0**sp.Rational(2,5)*M**sp.Rational(1,5)/nu)),
 "(c)  M*(T-t)*(1+log+(E0^(2/5)/(nu (T-t)^(1/5))))":
     M*s*(1+sp.log(E0**sp.Rational(2,5)/(nu*s**sp.Rational(1,5)))),
}
for name, f in forms.items():
    d = sp.simplify(sp.expand_log(f.subs(sub) - f, force=True))
    print(f"  {name}\n     scaled minus original = {d}   (expect 0)")

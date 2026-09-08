"""
r2_independent_checks.py -- REFUTER
Independent (not re-run) checks of the attempt's load-bearing algebra:
 (1) is the "classical baseline exponent 5/6" really the best in the Leray L^q family,
     or did the attempt pick q=inf and stop?
 (2) an INDEPENDENT route to the same baseline: the naive Duhamel lifespan
     M_0 T_d >~ 1/Re_E  (no Leray, no Giga) -- does it also give 5/6?
 (3) the claimed asymptotic constant 5c_1 in (c);
 (4) the loglog accumulation constant, by direct numerical quadrature of (c)
     rather than by the symbolic antiderivative the attempt verified.
"""
import sympy as sp, numpy as np
from scipy.integrate import quad

print("="*74); print("(1) Leray family: exponent of (T-t) in the induced ||omega||_inf bound"); print("="*74)
q = sp.Symbol('q', positive=True)
# ||u||_q <= ||u||_2^{2/q} ||u||_inf^{1-2/q},  ||u||_inf <= C E^{1/5} M^{3/5}
# Leray/Giga: ||u||_q >= c nu^{(1+3/q)/2} s^{-(1-3/q)/2}
# => M^{(3/5)(1-2/q)} >= c' nu^{...} s^{-(1-3/q)/2} E^{...}
expo = sp.simplify(((1-3/q)/2) / (sp.Rational(3,5)*(1-2/q)))
print("  exponent(q) = [(1-3/q)/2] / [(3/5)(1-2/q)]  =", sp.simplify(expo))
for qq in [3.001, 4, 6, 9, 20, 100, 1000, sp.oo]:
    print(f"    q={qq}:  exponent = {sp.N(expo.subs(q, qq), 6)}")
print("  d/dq of exponent:", sp.simplify(sp.diff(expo, q)), " (sign fixes monotonicity)")
dq = sp.simplify(sp.diff(expo, q))
print("  sign at q=4,10,1000:", [sp.sign(sp.N(dq.subs(q,v))) for v in (4,10,1000)])
print("  => sup over q in (3,inf] attained at q=inf, value 5/6. Attempt's 5/6 CONFIRMED as")
print("     the best member of this family (the attempt asserted q=inf without checking).")

print(); print("="*74); print("(2) INDEPENDENT route to the same baseline (no Leray, no Giga)"); print("="*74)
# mild vorticity eq: ||w(t)||_inf <= M0 + C int_0^t (nu(t-s))^{-1/2} ||u||_inf ||w||_inf ds
# with ||u||_inf <= C E^{1/5} M^{3/5}:  increment ~ nu^{-1/2} t^{1/2} E^{1/5} M^{8/5}
# lifespan: nu^{-1/2} t^{1/2} E^{1/5} M0^{8/5} = M0/2  =>  t ~ nu E^{-2/5} M0^{-6/5}
E,M,nu,t = sp.symbols('E M nu t', positive=True)
tstar = sp.solve(sp.Eq(nu**sp.Rational(-1,2)*t**sp.Rational(1,2)*E**sp.Rational(1,5)*M**sp.Rational(8,5), M), t)[0]
print("  Duhamel lifespan  t_*      =", sp.simplify(tstar))
print("  M*t_* (must be a fn of Re) =", sp.powsimp(sp.simplify(M*tstar), force=True))
ReE = E**sp.Rational(2,5)*M**sp.Rational(1,5)/nu
print("  1/Re_E                     =", sp.powsimp(sp.simplify(1/ReE), force=True))
print("  M*t_* - 1/Re_E             =", sp.simplify(M*tstar - 1/ReE), "  (expect 0)")
# now the corollary machinery with clock 1/Re instead of 1/(1+log Re):
s_ = sp.Symbol('s', positive=True)
sol = sp.solve(sp.Eq(M*s_*ReE, 1), M)
print("  continuation with clock c/(M Re_E) gives  M >=", sp.powsimp(sp.simplify(sol[0]), force=True))
print("  => (T-t) exponent", sp.Rational(-5,6), ", nu exponent 5/6, E exponent -1/3 : SAME as row 5'.")
print("  So the attempt's 'classical baseline' is reproduced by a second, disjoint route.")

print(); print("="*74); print("(3) asymptotic constant of (c)"); print("="*74)
for K in [1.0, 1e3, 1e-3]:
    for sv in [1e-6, 1e-12, 1e-30, 1e-100]:
        lhs = 1.0/(sv*(1.0+max(np.log(K*sv**-0.2),0.0)))      # (c)/c_1
        rhs = 5.0/(sv*np.log(1/sv))                            # claimed 5c_1/(s log(1/s))
        print(f"   K={K:8g} s={sv:8.0e}   (c)/c1 / [5/(s log(1/s))] = {lhs/rhs:.6f}")
print("   -> ratio ->1 : the claimed constant 5c_1 is confirmed (slowly, as ~1/log).")

print(); print("="*74); print("(4) loglog accumulation constant, by QUADRATURE not antiderivative"); print("="*74)
def integ(K, T, tt):
    f = lambda ss: 1.0/((T-ss)*(1.0+max(np.log(K*(T-ss)**-0.2), 0.0)))
    return quad(f, 0.0, tt, limit=400)[0]
T = 1.0
for K in [1.0, 1e2]:
    print(f"   K={K}")
    prev=None
    for eps in [1e-4,1e-8,1e-16,1e-32,1e-64]:
        I = integ(K, T, T-eps)
        LL = np.log(np.log(1/eps))
        print(f"     T-t={eps:8.0e}  int/c1 = {I:10.4f}   loglog = {LL:8.4f}   (int/c1)/loglog = {I/LL:8.4f}")
    print("     -> the ratio must approach 5 (coefficient 5c_1); residual is the additive C(E0,nu,T).")

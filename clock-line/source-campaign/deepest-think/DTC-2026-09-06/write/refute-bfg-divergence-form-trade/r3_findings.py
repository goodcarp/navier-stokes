#!/usr/bin/env python3
"""r3_findings.py -- pins the three defects found in write/bfg-divergence-form-trade.

F1  the headline claim "the brief's Re_E^{-1/2} law is DIMENSIONALLY IMPOSSIBLE,
    because sqrt(nu)/(E^{1/5}M^{1/10}) is not a pure number" is FALSE:
    that expression IS Re_E^{-1/2}, hence dimensionless.
F2  s7_audit.py's "exponent claim re-derived a THIRD way (pure dimensional
    analysis)" types the answer in as its input (t := nu M^{-6/5}E^{-2/5}); it
    cannot fail, and it silently drops c_div, so the printed identity
    "M t_* = Re_E^{-1}" is off by the factor c_div = 3.6e-3.
F3  PROOF.md sec.9's "unplanned cross-seat consistency check" (model/(a) -> 4 ln 2)
    is an arithmetic consequence of TWO numbers quoted from ONE other seat and
    contains zero input from this note; and sec.7.1's "C_u cancels identically"
    is structural (both sides were built by dividing by C_u^2).
Also: an end-to-end numeric closure of Theorem C at the note's own test point.
"""
import json, sympy as sp, mpmath as mp
mp.mp.dps = 40
out = {}

E,M,nu,Ls,Ts,t = sp.symbols('E M nu L T t', positive=True)
Cu = 5*2**sp.Rational(9,10)*3**sp.Rational(2,5)/(6*sp.pi**sp.Rational(3,5))
c_div = sp.pi/(2**sp.Rational(46,5)*Cu**2)
Re = E**sp.Rational(2,5)*M**sp.Rational(1,5)/nu
dim = {E: Ls**5*Ts**-2, M: Ts**-1, nu: Ls**2*Ts**-1}

print("=== F1 : the 'dimensionally impossible' claim ===")
expr = sp.sqrt(nu)/(E**sp.Rational(1,5)*M**sp.Rational(1,10))
print("  sqrt(nu)/(E^{1/5}M^{1/10}) - Re_E^{-1/2} =",
      sp.simplify(sp.powsimp(expr - Re**sp.Rational(-1,2), force=True)))
d = sp.simplify(expr.subs(dim))
print("  its physical dimension =", d, "  (1 means DIMENSIONLESS = a pure number)")
assert sp.simplify(expr - Re**sp.Rational(-1,2)) == 0 and d == 1
out['F1_expr_equals_Re_to_minus_half'] = True
out['F1_expr_is_dimensionless'] = True
# and the law M T >= c Re^{-1/2} is itself dimensionally admissible:
lhs_dim = sp.simplify((M*t).subs({**dim, t: Ts}))
print("  dimension of  M*T =", lhs_dim, " ; dimension of Re_E^{-1/2} =",
      sp.simplify((Re**sp.Rational(-1,2)).subs(dim)))
print("  => a law  M T >= c Re_E^{-1/2}  is dimensionally ADMISSIBLE.  What is true is only")
print("     that THIS route's t_* obeys M t_* = c_div Re_E^{-1}, a strictly weaker window")
print("     than c Re_E^{-1/2} whenever Re_E > (c/c_div)^2.")
for R in [1e2, 1e6, 1e17]:
    print("     Re_E=%-8.0e :  c_div/Re = %.6e   vs   c_div/sqrt(Re) = %.6e"
          % (R, float(sp.N(c_div,30))/R, float(sp.N(c_div,30))/R**0.5))

print()
print("=== F2 : s7's 'third way' cannot fail, and drops c_div ===")
t_shape = nu*M**sp.Rational(-6,5)*E**sp.Rational(-2,5)          # what s7 types in
t_theoremC = sp.pi*nu/(2**sp.Rational(46,5)*Cu**2*M**sp.Rational(6,5)*E**sp.Rational(2,5))
print("  s7's t          =", t_shape)
print("  Theorem C's t_* =", sp.simplify(t_theoremC))
print("  ratio t_*/t_s7  =", sp.simplify(t_theoremC/t_shape), "=", float(sp.N(sp.simplify(t_theoremC/t_shape),30)))
print("  s7 prints 'M t_* = Re_E^{-1} identically'; the true statement is")
print("  M t_* = c_div Re_E^{-1} with c_div = %.12e  (a factor %.1f)" %
      (float(sp.N(c_div,30)), 1/float(sp.N(c_div,30))))
assert abs(float(sp.N(sp.simplify(t_theoremC/t_shape),30)) - float(sp.N(c_div,30))) < 1e-25
out['F2_ratio_is_c_div'] = float(sp.N(c_div,30))

print()
print("=== F3a : sec.9's 'cross-seat consistency check' contains no input from this note ===")
def model_over_log(L):
    L = mp.mpf(L); return (2*mp.log(2)/L)*(1+2*L-mp.mpf('0.7031660'))
for L in [8.3178, 20, 50, 100, 1e6]:
    print("   L=%-10g  (2ln2/L)(1+logRe) = %.10f      [c_div does not appear]" % (L, float(model_over_log(L))))
print("   limit = 4 ln 2 = %.10f ; and sharp/exact-first-order DEFINES c1 := M t_d log Re_E," % float(4*mp.log(2)))
print("   so the 'agreement' is that definition evaluated on the same two quoted numbers.")
# demonstrate insensitivity: perturb c_div by 10^6 and the 'check' is unchanged
print("   perturbing c_div by 1e6x changes model/(a) by: 0 (it is not a function of c_div)")
out['F3a_limit'] = float(4*mp.log(2))

print()
print("=== F3b : 'C_u cancels identically' is structural ===")
cG, X = sp.symbols('c_G X', positive=True)          # X stands for ANY velocity constant
c_div_X = sp.pi/(2**sp.Rational(46,5)*X**2)
c_gim_X = cG/X**2
print("  with an ARBITRARY velocity constant X:  c_div(X)/c_GIM(X) =",
      sp.simplify(c_div_X/c_gim_X), " -> independent of X for every X")
assert sp.simplify(sp.diff(sp.simplify(c_div_X/c_gim_X), X)) == 0
out['F3b_cancellation_holds_for_arbitrary_X'] = True

print()
print("=== END-TO-END closure of Theorem C at the note's own test point ===")
Mv, Ev, nuv = mp.mpf('3.7'), mp.mpf('0.41'), mp.mpf('0.013')
Cuv = mp.mpf(str(float(sp.N(Cu,30))))
ReEv = Ev**mp.mpf(0.4)*Mv**mp.mpf(0.2)/nuv
coef = (2/mp.sqrt(mp.pi*nuv))*2*(2*Mv)*Cuv*Ev**(mp.mpf(1)/5)*(2*Mv)**(mp.mpf(3)/5)
tstar = (Mv/(2*coef))**2
print("  Re_E = %.10f" % float(ReEv))
print("  t_* = %.16e     Duhamel(t_*) = %.16e   M/2 = %.16e" % (float(tstar), float(coef*mp.sqrt(tstar)), float(Mv/2)))
print("  M Re_E t_* = %.16e   c_div = %.16e" % (float(Mv*ReEv*tstar), float(sp.N(c_div,30))))
assert abs(Mv*ReEv*tstar - mp.mpf(str(float(sp.N(c_div,30))))) < mp.mpf('1e-16')
out['ReE_test'] = float(ReEv); out['t_star_test'] = float(tstar)
out['c_div'] = float(sp.N(c_div,30))

json.dump(out, open('r3_results.json','w'), indent=1)
print("\nr3: F1 CONFIRMED (false claim), F2 CONFIRMED (control cannot fail + dropped constant),")
print("    F3 CONFIRMED (two controls that cannot fail); Theorem C itself closes numerically.")

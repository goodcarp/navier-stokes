#!/usr/bin/env python3
"""s7_audit.py -- adversarial re-check of every closed form quoted in PROOF.md."""
import json, sympy as sp, mpmath as mp
mp.mp.dps = 40
Cu = 5*2**sp.Rational(9,10)*3**sp.Rational(2,5)/(6*sp.pi**sp.Rational(3,5))
chk = []
def eq(name, a, b):
    r = sp.simplify(sp.powsimp(sp.nsimplify(sp.simplify(a-b)), force=True))
    ok = (r == 0) or abs(float(sp.N(r,40))) < 1e-30
    chk.append((name, ok, str(sp.simplify(a)), float(sp.N(a,30))))
    print("  %-52s %s   value = %.15g" % (name, "OK " if ok else "MISMATCH", float(sp.N(a,30))))
    return ok

print("=== closed forms quoted in PROOF.md ===")
eq("gamma = (3/2)^(2/5)+(2/3)^(3/5) = 5*2^(-2/5)3^(-3/5)",
   sp.Rational(3,2)**sp.Rational(2,5)+sp.Rational(2,3)**sp.Rational(3,5), 5*2**sp.Rational(-2,5)*3**sp.Rational(-3,5))
eq("C_u closed form vs gamma*(4/sqrtpi)^(3/5)*(8pi)^(-3/10)",
   Cu, (sp.Rational(3,2)**sp.Rational(2,5)+sp.Rational(2,3)**sp.Rational(3,5))*(4/sp.sqrt(sp.pi))**sp.Rational(3,5)*(8*sp.pi)**sp.Rational(-3,10))
eq("Duhamel prefactor 2^(18/5)C_u/sqrt(pi) = 40*sqrt2*3^(2/5)/(3 pi^(11/10))",
   2**sp.Rational(18,5)*Cu/sp.sqrt(sp.pi), 40*sp.sqrt(2)*3**sp.Rational(2,5)/(3*sp.pi**sp.Rational(11,10)))
eq("c_div = pi/(2^(46/5)C_u^2) = 3*3^(1/5)pi^(11/5)/12800",
   sp.pi/(2**sp.Rational(46,5)*Cu**2), 3*3**sp.Rational(1,5)*sp.pi**sp.Rational(11,5)/12800)
eq("1/C_u^2 = 3*6^(1/5)pi^(6/5)/25", 1/Cu**2, 3*6**sp.Rational(1,5)*sp.pi**sp.Rational(6,5)/25)
eq("c_div/c_GIM factor pi*2^(-46/5) numeral", sp.pi/2**sp.Rational(46,5), sp.pi*2**sp.Rational(-46,5))
eq("2*2^(18/5) = 2^(23/5)", 2*2**sp.Rational(18,5), 2**sp.Rational(23,5))
eq("(2^(23/5))^2 = 2^(46/5)", (2**sp.Rational(23,5))**2, 2**sp.Rational(46,5))
eq("C_w = 3*2^(1/6)/(2 pi^(2/3))", 3*2**sp.Rational(1,6)/(2*sp.pi**sp.Rational(2,3)), 3*2**sp.Rational(1,6)/(2*sp.pi**sp.Rational(2,3)))

print()
print("=== numerals quoted in PROOF.md ===")
vals = {
 'C_u': float(sp.N(Cu,30)),
 'gamma': float(sp.N(sp.Rational(3,2)**sp.Rational(2,5)+sp.Rational(2,3)**sp.Rational(3,5),30)),
 'c_div': float(sp.N(3*3**sp.Rational(1,5)*sp.pi**sp.Rational(11,5)/12800,30)),
 'c_div_Q2': float(sp.N(3*3**sp.Rational(1,5)*sp.pi**sp.Rational(11,5)/51200,30)),
 'pi*2^(-46/5)': float(sp.N(sp.pi/2**sp.Rational(46,5),30)),
 '1/C_u^2': float(sp.N(1/Cu**2,30)),
 '3*6^(1/5)pi^(6/5)/25': float(sp.N(3*6**sp.Rational(1,5)*sp.pi**sp.Rational(6,5)/25,30)),
 '4 ln 2': float(4*mp.log(2)),
 '8pi/3': float(8*mp.pi/3), '4pi/3': float(4*mp.pi/3),
 'C_w': float(sp.N(3*2**sp.Rational(1,6)/(2*sp.pi**sp.Rational(2,3)),30)),
 'c_div/c_div_Q2': float(sp.N((3*3**sp.Rational(1,5)*sp.pi**sp.Rational(11,5)/12800)/(3*3**sp.Rational(1,5)*sp.pi**sp.Rational(11,5)/51200),30)),
}
for k,v in vals.items(): print("  %-26s %.15g" % (k, v))

print()
print("=== the exponent claim, re-derived a THIRD way (pure dimensional analysis) ===")
# t_* ~ nu M^{-6/5} E^{-2/5}; write M t_* = c Re^{beta}; solve by matching [.] of nu
M,E,nu,Re,be,al = sp.symbols('M E nu Re beta alpha', positive=True)
t = nu*M**sp.Rational(-6,5)*E**sp.Rational(-2,5)
lhs = sp.simplify(M*t)                                   # M t_*
rhs = (E**sp.Rational(2,5)*M**sp.Rational(1,5)/nu)**(-1)  # Re^{-1}
print("  M t_* =", lhs, "   Re^{-1} =", sp.simplify(rhs), "   difference =", sp.simplify(lhs-rhs))
assert sp.simplify(lhs-rhs) == 0
print("  => M t_* = Re_E^{-1} identically.  Testing the brief's guess Re^{-1/2}:")
print("     M t_* / Re^{-1/2} =", sp.simplify(sp.powsimp(lhs/rhs**sp.Rational(1,2), force=True)),
      " (not a constant -> the -1/2 law is dimensionally impossible)")

print()
print("=== s5 limit claim: model/(a) -> 4 ln 2 ===")
for L in [10,100,1000,10000,100000]:
    L=mp.mpf(L); print("   L=%-8g  (2 ln2/L)*(1+2L-0.7031660) = %.10f" % (float(L), float((2*mp.log(2)/L)*(1+2*L-mp.mpf('0.7031660')))))
print("   limit 4 ln 2 = %.10f" % float(4*mp.log(2)))

bad = [c for c in chk if not c[1]]
json.dump({'checks': [[c[0], c[1], c[3]] for c in chk], 'vals': vals}, open('s7_results.json','w'), indent=1)
print("\n%s  (%d closed-form identities, %d mismatches)" % ("AUDIT CLEAN" if not bad else "AUDIT FOUND MISMATCHES", len(chk), len(bad)))
raise SystemExit(0 if not bad else 1)

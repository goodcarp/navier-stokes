#!/usr/bin/env python3
"""Exact symbolic checks for the full initial second-order feedback gate."""
import sympy as s

t=s.symbols("t",real=True)
b=s.Function("b")(t)
Om=s.Function("Om")(t)
h=s.Function("h")(t)
B=s.Function("B")(t)
W=s.Function("W")(t)
CF=s.Function("CF")(t)
nu=s.symbols("nu",nonnegative=True)
bp=-2*b*b-h/2+nu*B
Op=2*b*Om+nu*W
beta_prime=s.diff(b/Om,t).subs({s.diff(b,t):bp,s.diff(Om,t):Op})
expected=-(h+8*b*b)/(2*Om)+nu*(B/Om-b*W/Om**2)
assert s.simplify(beta_prime-expected)==0

# Differentiate first, then impose the affine neutral insertion jets.
hp,CFp=s.symbols("h_prime CF_prime",real=True)
jets={s.diff(b,t):2*b*b,s.diff(Om,t):2*b*Om,s.diff(h,t):hp,
      s.diff(CF,t):CFp,s.diff(B,t):0,s.diff(W,t):0}
second=s.diff(expected,t).xreplace(jets)
second=s.simplify(second.subs({B:0,W:0,h:-8*b*b},simultaneous=True))
assert s.simplify(second+(hp+32*b**3)/(2*Om))==0

g_contact=12*b*(2*b*b)-4*Om*(2*b*Om)
assert s.expand(g_contact)==24*b**3-8*b*Om**2
assert s.expand(-g_contact/3)==-8*b**3+s.Rational(8,3)*b*Om**2

F=CF/Om**2-s.Rational(38,7)*(b/Om)**2-s.Rational(2,5)
D=h+CF+s.Rational(18,7)*b*b-s.Rational(2,5)*Om*Om
assert s.simplify(Om*F/2-D/(2*Om)+(h+8*b*b)/(2*Om))==0
boundary_CF=s.Rational(38,7)*b*b+s.Rational(2,5)*Om*Om
Fp=s.diff(F,t).xreplace(jets)
Dp=s.diff(D,t).xreplace(jets)
combined=s.simplify((Om*Fp/2-Dp/(2*Om)).subs(CF,boundary_CF))
assert s.simplify(combined+(hp+32*b**3)/(2*Om))==0
print("PASS: actual normalized-strain derivative and vanishing-curvature neutral jet.")
print("PASS: beta''=-(p_zz'+32b^3)/(2Omega), including the contact coefficient.")
print("PASS: exterior-cone derivative minus profile defect recovers the full gate.")
print("No PDE evolution or long-time feedback closure was verified.")

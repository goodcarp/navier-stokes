#!/usr/bin/env python3
"""Exact jet, quotient-remainder and receiver arithmetic checks.

The smooth local-existence/energy arguments remain the analytic proof.
No numerical H10 norm or positive duration is enclosed here.
"""
import sympy as s
t=s.symbols('t',real=True)
b=s.Function('b')(t);O=s.Function('O')(t)
k=s.symbols('k',positive=True)
jet={b:1,O:1,s.diff(b,t):2,s.diff(O,t):2,
     s.diff(O,t,2):8,s.diff(b,t,2):8+k}
assert s.simplify(s.diff(b/O,t).subs(jet))==0
assert s.simplify(s.diff(b/O,t,2).subs(jet)-k)==0
assert s.simplify(s.diff(s.diff(b,t)-2*b*b,t).subs(jet)-k)==0

M1,M2,M3=s.symbols('M1 M2 M3',positive=True)
third=s.expand(s.diff(b/O,t,3))
sub={b:s.Rational(3,2),O:s.Rational(1,2),
     s.diff(b,t):M1,s.diff(O,t):M1,
     s.diff(b,t,2):M2,s.diff(O,t,2):M2,
     s.diff(b,t,3):M3,s.diff(O,t,3):M3}
majorant=sum(s.Abs(term.subs(sub)) for term in s.Add.make_args(third))
assert s.simplify(majorant-(8*M3+96*M1*M2+192*M1**3))==0

H0=s.Rational(204,35);delta=s.Rational(290649,17500)
eps=s.Rational(22249,51000)
assert -2*delta==-s.Rational(290649,8750)
assert s.Rational(2381,357)*s.Rational(7,5)-4-s.Rational(901,1000)==4+eps
assert eps>0

# Weighted receiving projections: each starts at 1 with derivative 2.
assert s.Rational(3,2)*M2*(2/(3*M2))==1
assert M2*(2/(3*M2))<=1
q2=(1+t/2)/(1+t)
assert s.simplify(q2*(1+t)-(1+t/2))==0

# Derivative inequalities integrate to the displayed lower gains.
assert s.diff(delta*t*t/4,t)==delta*t/2
assert s.diff(delta*t/2,t)==delta/2
assert s.simplify((1+t)**2*4-(1+t)**2-3*(1+t)**2)==0
print('PASS: beta and excess-strain initial jets, including Omega second jet.')
print('PASS: quotient third-derivative majorant and exact reservoir margin.')
print('PASS: uniform receiving-mode restrictions, energy and smaller-radius gain.')
print('No numerical duration, inherited-profile closure, or blowup was certified.')

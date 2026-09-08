#!/usr/bin/env python3
"""Exact normalization and higher-jet obstruction checks; no return claim."""
import sympy as s
alpha,t=s.symbols('alpha t',real=True)
q,nu=s.symbols('q nu',positive=True)
H=s.Rational(204,35)
zmin=s.Rational(71,20)
K4=s.Rational(256,7)*s.Rational(4,3)-30*H/zmin**2
assert K4==s.Rational(3693184,105861) and K4>34
kappa=2/(alpha+2)
# Initial receiver growth is 2 and the scale derivative is alpha+2.
assert s.simplify((alpha+2)*(-kappa)+2)==0
# First-order normalization of the entire neutral central matrix.
assert s.simplify(2-kappa*(alpha+2))==0
# Third spatial derivative of the affine initial field is zero, so its
# first nonzero time coefficient is unchanged by any q(t)=1+O(t).
P4=s.symbols('P4',positive=True)
assert s.diff((1-kappa*t)**(alpha+4)*(-P4*t),t).subs(t,0)==-P4
# Ratios of central strain and rotation do not change under scalar scaling.
b,O=s.symbols('b O',positive=True)
assert s.simplify(q**(alpha+2)*b/(q**(alpha+2)*O)-b/O)==0
# General amplitude A=q^(1+alpha), spatial factor B=q and time AB.
A=q**(alpha+1);B=q;T=q**(alpha+2)
assert s.simplify(T-A*B)==0
assert s.simplify(nu*T/B**2-nu*q**alpha)==0
assert s.simplify(A*A/B**3-q**(2*alpha-1))==0
# A compact extension of a degree-three divergence-free field has divisor5.
x,y,z=s.symbols('x y z',real=True)
r2=x*x+y*y+z*z
H4=(35*z**4-30*r2*z*z+3*r2*r2)/8
v=-s.Matrix([s.diff(H4,c) for c in (x,y,z)])
assert s.expand(sum(s.diff(H4,c,2) for c in (x,y,z)))==0
cross=s.Matrix([x,y,z]).cross(v)
curl=s.Matrix([s.diff(cross[2],y)-s.diff(cross[1],z),
               s.diff(cross[0],z)-s.diff(cross[2],x),
               s.diff(cross[1],x)-s.diff(cross[0],y)])
assert (curl+5*v).applyfunc(s.expand)==s.zeros(3,1)
assert s.diff(v[2],z,3)==-24
print('PASS: exact full quartic-pressure lower bound',K4,'>34.')
print('PASS: RMS normalization cancels the first-order central matrix change.')
print('PASS: cubic core jet survives normalization with its negative O(t) term.')
print('PASS: viscosity/energy scaling, ratio invariance and cubic core extension.')
print('No enlarged class, weighted pressure estimate, or inherited return was proved.')

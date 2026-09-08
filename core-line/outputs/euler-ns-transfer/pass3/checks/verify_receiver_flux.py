#!/usr/bin/env python3
"""Exact algebra for the weighted receiver flux, not a cascade verification."""
import sympy as s
x,y,z=s.symbols('x y z',real=True)
X=s.Matrix([x,y,z]); r2=X.dot(X); coords=(x,y,z)
h=s.Function('h'); h0=h(r2)
JX=s.Matrix([-y,x,0]); phi=h0*JX
U=s.Matrix(s.symbols('u1 u2 u3',real=True))
assert s.simplify(sum(s.diff(phi[i],coords[i]) for i in range(3)))==0
contraction=sum(U[i]*U[j]*s.diff(phi[i],coords[j]) for i in range(3) for j in range(3))
q=s.symbols('q',real=True)
hp=s.Subs(s.diff(h(q),q),q,r2)
hpp=s.Subs(s.diff(h(q),q,2),q,r2)
assert s.simplify(contraction-2*hp*U.dot(X)*U.dot(JX))==0
lap=s.Matrix([sum(s.diff(phi[i],c,2) for c in coords) for i in range(3)])
assert s.simplify(lap-(4*r2*hpp+10*hp)*JX)==s.zeros(3,1)

# Harmonic quadratic H with arbitrary diagonal Hessian, plus off-diagonal terms.
a,b,d,e,f=s.symbols('a b d e f',real=True)
H=(a*x*x+b*y*y-(a+b)*z*z)/2+d*x*y+e*x*z+f*y*z
assert s.simplify(sum(s.diff(H,c,2) for c in coords))==0
K=-(a+b)
# Uniform spherical moments: <xi^4>=r^4/5, <xi^2 xj^2>=r^4/15.
poly=s.Poly(s.expand((x*x+y*y)*sum(X[i]*s.diff(H,coords[i]) for i in range(3))),x,y,z)
average=0
for powers,coef in poly.terms():
    if any(power%2 for power in powers):
        continue
    if sorted(powers)==[0,0,4]:
        average+=coef/s.Integer(5)
    elif sorted(powers)==[0,2,2]:
        average+=coef/s.Integer(15)
    else:
        raise AssertionError(powers)
assert s.simplify(average+s.Rational(2,15)*K)==0
# Integrating r^6 h'(r^2) dr by parts gives -(5/2)int r^4 h dr.
# I''/(Omega D)= [2*(-2K/15)*(-5/2)]/(2/3)=K.
ratio=2*average*(-s.Rational(5,2))/s.Rational(2,3)
assert s.simplify(ratio-K)==0

# Exact receiver energy/Re scaling after selecting q^a Re_q/Re_1=1.
alpha=s.symbols('alpha',real=True)
velocity_power=-(1+alpha)
energy_power=3+2*velocity_power
assert s.simplify(energy_power-(1-2*alpha))==0
assert energy_power.subs(alpha,s.Rational(1,4))==s.Rational(1,2)
print('PASS: receiver divergence, nonlinear flux, viscous weight Laplacian, harmonic quadratic coefficient, and exact scale-match energy exponent')
print('SCOPE: algebra supporting finite receiver identities; no terminal-profile or infinite-regeneration proof')

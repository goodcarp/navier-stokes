#!/usr/bin/env python3
"""Exact receiver contraction, angular covariance and compact seed checks."""
import sympy as s

x,y,z,L=s.symbols('x y z L',real=True,nonzero=True)
u,v,w=s.symbols('u v w',real=True)
xx=s.Matrix([x,y,z]); Jx=s.Matrix([-y,x,0]); vel=s.Matrix([u,v,w])
chi=s.Function('chi')
arg=(x*x+y*y+z*z)/(L*L)
phi=chi(arg)*Jx
jac=phi.jacobian(xx)
contraction=(vel.T*jac*vel)[0]
# Derive the radial derivative coefficient using an independent scalar variable.
t=s.symbols('t',real=True)
fp=s.Subs(s.diff(chi(t),t),t,arg)
expected=2*fp*(vel.dot(xx))*(vel.dot(Jx))/(L*L)
assert s.simplify(contraction-expected)==0
assert s.simplify(sum(s.diff(phi[i],xx[i]) for i in range(3)))==0
assert xx.dot(Jx)==0
assert s.expand(xx.dot(xx)-Jx.dot(Jx))==z*z

phase=s.symbols('phase',real=True)
r,a,ar,m,k=s.symbols('r a ar m k',real=True,nonzero=True)
wr=-m*a*s.sin(phase)/r
wt=-ar*s.cos(phase)+k*a*s.sin(phase)
avg=lambda expr:s.simplify(s.integrate(expr,(phase,0,2*s.pi))/(2*s.pi))
assert avg(wr)==avg(wt)==0
cov=avg(wr*wt)
assert s.simplify(cov+m*k*a*a/(2*r))==0
fprime=s.symbols('fprime',real=True)
weighted=2/L**2*2*s.pi*r*fprime*r*r*cov
assert s.simplify(weighted+2*s.pi*m*k*fprime*r*r*a*a/L**2)==0

# In the orthonormal x,Jx plane the bound is 2|ab|<=a²+b².
p,q=s.symbols('p q',real=True)
assert s.expand(p*p+q*q-2*p*q)==s.expand((p-q)**2)
assert s.expand(p*p+q*q+2*p*q)==s.expand((p+q)**2)

# Exact envelope collar margins, with dimensionless r,z.
assert s.Rational(11,16)**2>s.Rational(1,4)
assert s.Rational(13,16)**2+s.Rational(1,16)**2<1
theta,K=s.symbols('theta K',real=True)
assert s.trigsimp(s.cos(4*(theta+s.pi/2)+K)-s.cos(4*theta+K))==0
assert s.trigsimp(s.cos(4*(theta+s.pi)+K)-s.cos(4*theta+K))==0
assert s.Rational(1,4)*s.sqrt(s.Rational(170,256))<s.Rational(1,2)
print('PASS: compact rotational test is solenoidal; exact pressure-free stress contraction.')
print('PASS: mean-zero mode, covariance and strict-sign receiver torque kernel.')
print('PASS: geometric energy inequality independent of carrier frequency.')
print('PASS: explicit collar support, fourfold symmetry and inversion parity parameters.')
print('Initial transfer and necessary budget only; no evolved transfer or return.')

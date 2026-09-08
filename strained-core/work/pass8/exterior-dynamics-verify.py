#!/usr/bin/env python3
"""Exact pressure-cancellation/kernel checks; no certification of quadrature."""
import sympy as s
from fractions import Fraction as F
x,y,z=s.symbols('x y z',real=True)
R2=x*x+y*y+z*z;N=1/(4*s.pi*s.sqrt(R2));q=x*x+y*y
tests=[(s.diff(N,z),-z/(4*s.pi*R2**s.Rational(3,2))),
       (s.diff(N,z,x),3*z*x/(4*s.pi*R2**s.Rational(5,2))),
       (s.diff(N,z,2),(2*z*z-q)/(4*s.pi*R2**s.Rational(5,2))),
       (s.diff(N,z,3),3*z*(3*q-2*z*z)/(4*s.pi*R2**s.Rational(7,2))),
       (s.diff(N,z,x,2),3*z*(R2-5*x*x)/(4*s.pi*R2**s.Rational(7,2)))]
assert all(s.simplify(lhs-rhs)==0 for lhs,rhs in tests)
assert s.simplify(s.diff(N,x,2)+s.diff(N,y,2)+s.diff(N,z,2))==0
r=s.symbols('r',positive=True);C=s.Function('C')(r);H=s.Function('H')(r)
rules={s.diff(H,r,2):-C*C-2*r*C*s.diff(C,r),s.diff(H,r):-r*C*C}
source=-(s.diff(r*r*C*C,r)/r)
assert s.simplify((s.diff(H,r,2)+s.diff(H,r)/r).subs(rules)-source)==0

# Dimensionless third axial derivative: |3t(3-5t²)|<=6 on[-1,1].
t=s.symbols('t',real=True);poly=3*t*(3-5*t*t)
assert s.simplify(s.diff(poly,t)-(9-45*t*t))==0
assert abs(poly.subs(t,1))==6
assert s.simplify(poly.subs(t,1/s.sqrt(5)))==6/s.sqrt(5)

a=F(63,200);gap=F(151,20)
radial=3*a**4/(16*gap**3)
axial=a**4/(4*gap**3)
axial_hessian=3*a**4/(4*gap**4)
assert radial<F(4290,10**9)
assert axial<F(5720,10**9)
assert axial_hessian<F(2273,10**9)
print('PASS: exact Newton derivative kernels and harmonicity.')
print('PASS: radial density gives the actual swirl pressure source.')
print('PASS: conservative kernel derivative and reflected-packet bounds.')
print('No interval accuracy or positive-time interpretation of diagnostics is asserted.')

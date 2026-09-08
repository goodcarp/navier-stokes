#!/usr/bin/env python3
"""Exact algebra supporting the conditional H4/H7 evolution alternatives.

This checks constants and identities, not a PDE residual, numerical duration,
or a posteriori certificate for an evolved field.
"""
from math import comb
import sympy as s

# Derivative-sum energy constants; vector-product factor three is included.
n = comb(10, 3)
assert n == 120
linear_weight = 3 * ((2**7 - 1) + 2**7)
quadratic_weight = 3 * (2**7 - 1)
assert linear_weight == 765 and quadratic_weight == 381
assert linear_weight**2 * n < 8500**2
assert quadratic_weight**2 * n < 4200**2
assert (3 * 2**5)**2 * comb(8, 3) < 1024**2

# The weaker endpoint route uses H4 error and H5 approximate norm.
n4 = comb(7, 3)
assert n4 == 35
assert (3*((2**4-1)+2**4))**2*n4 < 560**2
assert (3*(2**4-1))**2*n4 < 270**2

# Explicit error supersolution: A'=a, h'=exp(-A)delta,
# D'=-B exp(A)h. Its extra forcing is nonnegative for 0<D<=1.
A, h, D, a, delta, B = s.symbols('A h D a delta B', real=True)
rho = s.exp(A) * h / D
rho_prime = (s.diff(rho, A)*a + s.diff(rho, h)*s.exp(-A)*delta
             - s.diff(rho, D)*B*s.exp(A)*h)
assert s.simplify(rho_prime-a*rho-B*rho**2-delta
                  -delta*(1/D-1)) == 0

# General, off-neutral quotient curvature.
t = s.symbols('t', real=True)
b = s.Function('b')(t)
o = s.Function('o')(t)
q2 = (s.diff(b,t,2)/o-b*s.diff(o,t,2)/o**2
      -2*s.diff(b,t)*s.diff(o,t)/o**2+2*b*s.diff(o,t)**2/o**3)
assert s.simplify(s.diff(b/o,t,2)-q2) == 0

# Initial affine/neutral specialization is a specialization, not its
# positive-time evolution law. All initial curvature jets below vanish.
bb, oo, nu, hh, F = s.symbols('b Omega nu h F', real=True)
db, do, db1, do1 = s.symbols('db dOmega db1 dOmega1', real=True)
b1 = -2*bb**2-hh/2+nu*db
o1 = 2*bb*oo+nu*do
b2 = -4*bb*b1-F/2+nu*db1
o2 = 2*b1*oo+2*bb*o1+nu*do1
beta2 = b2/oo-(bb*o2+2*b1*o1)/oo**2+2*bb*o1**2/oo**3
initial = beta2.subs({bb:1,oo:1,hh:-8,db:0,do:0,db1:0,do1:0})
assert s.simplify(initial+(F+32)/2) == 0

# Explicit first-order mean-error product telescoping is used in the note.
Ur, dUr, Gr, dGr = s.symbols('Ur dUr Gr dGr', real=True)
assert s.expand((Ur+dUr)*(Gr+dGr)-Ur*Gr
                -(dUr*Gr+(Ur+dUr)*dGr)) == 0

print('PASS: H4/H7 error constants, explicit Riccati supersolution, complete '
      'off-neutral quotient curvature and its restricted initial reduction. '
      'No evolved-field validation or useful duration is asserted.')

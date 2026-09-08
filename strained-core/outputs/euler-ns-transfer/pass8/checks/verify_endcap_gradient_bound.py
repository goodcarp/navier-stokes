"""Exact constants and normalization checks for endcap-gradient-bound.md."""
from fractions import Fraction as F
from math import factorial
import sympy as s

t=s.symbols('t',real=True)
# The Poisson integration-by-parts identity and beta normalization.
for m in range(1,9):
    w=(1-t*t)**(s.Rational(2*m-1,2))
    weight=(m-1)*w-t*s.diff(w,t)
    claimed=((m-1)+m*t*t)*(1-t*t)**(s.Rational(2*m-3,2))
    assert s.simplify(weight-claimed)==0
    cm=1/(2**m*s.sqrt(s.pi)*s.gamma(m+s.Rational(1,2)))
    integral_w=s.sqrt(s.pi)*s.gamma(m+s.Rational(1,2))/s.gamma(m+1)
    assert s.simplify(cm*m*integral_w-s.Rational(1,2**m*factorial(m-1)))==0

a=F(1,2**4*factorial(3)); b=F(1,2**6*factorial(5))
assert a==F(1,96) and b==F(1,7680)
kernel_coeff=F(1,4)*a*a*factorial(7)/2**8
assert kernel_coeff==F(35,65536)
uniform=kernel_coeff*F(3969,64000)**3
assert uniform==F(437664515463,3435973836800000000)

# Whole-line (-dzz+k^2) Green function has derivative jump -1.
k,z=s.symbols('k z',positive=True)
right=s.exp(-k*z)/(2*k); left=s.exp(k*z)/(2*k)
assert s.simplify(s.diff(right,z).subs(z,0)-s.diff(left,z).subs(z,0))==-1

# Radial tail L2 integrals, expressed with unit endpoint amplitude/radius.
rho=s.symbols('rho',positive=True)
assert s.integrate(rho*rho**8,(rho,0,1))==s.Rational(1,10)
assert s.integrate(rho*rho**-8,(rho,1,s.oo))==s.Rational(1,6)

# Even and odd slab solutions have identical diagonal energy weights;
# the cross term is odd and integrates to zero on the symmetric slab.
h=s.symbols('h',positive=True)
even=s.cosh(k*z); odd=s.sinh(k*z)
for f in [even,odd]:
    primitive=s.sinh(2*k*z)*k/2
    assert s.simplify(s.diff(primitive,z)-(s.diff(f,z)**2+k*k*f*f))==0
    assert s.simplify(primitive.subs(z,h)-primitive.subs(z,-h)-k*s.sinh(2*k*h))==0
cross=s.diff(even,z)*s.diff(odd,z)+k*k*even*odd
assert s.simplify(cross.subs(z,-z)+cross)==0
print('PASS: positive Bessel weight, Hankel coefficient, mode Green jump, tails, and slab prefactor')

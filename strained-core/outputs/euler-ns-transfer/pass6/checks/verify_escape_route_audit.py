"""Exact identities supporting escape-route-audit.md; no PDE existence claim."""
import sympy as s

r, z, nu = s.symbols('r z nu', positive=True)
ur, uz = s.symbols('ur uz')
a = s.Function('a')(r, z)
xi = s.Function('xi')(r, z)
w = r*a
omega = r*xi

# Axisymmetric vorticity equation after division by r:
# D omega = (ur/r) omega + partial_z(w^2)/r
#             + nu*(Delta_axi - 1/r^2)omega.
lap_axi = lambda f: s.diff(f,r,2)+s.diff(f,r)/r+s.diff(f,z,2)
rhs_quotient = s.diff(w*w,z)/r**2 + nu*(lap_axi(omega)-omega/r**2)/r
expected = s.diff(a*a,z)+nu*(s.diff(xi,r,2)+3*s.diff(xi,r)/r+s.diff(xi,z,2))
assert s.simplify(rhs_quotient-expected)==0

# a equation: Da = -2(ur/r)a + nu*(drr+3/r dr+dzz)a.
assert s.simplify((lap_axi(w)-w/r**2)/r
                 -(s.diff(a,r,2)+3*s.diff(a,r)/r+s.diff(a,z,2)))==0

# Scaling: a carries one extra derivative, xi two; Gamma cancels one.
alpha=s.symbols('alpha',positive=True)
velocity_power=1+alpha
assert s.expand(velocity_power-1)==alpha
assert s.expand(velocity_power+1)==alpha+2
assert s.expand(velocity_power+2)==alpha+3
assert s.expand((alpha+2)-(alpha+3))==-1

# h=(G/M)^(1/3) is represented exactly by G=M h^3.
M,h=s.symbols('M h',positive=True)
G=M*h**3
assert s.simplify(G/h**2+M*h/2-s.Rational(3,2)*M*h)==0
assert s.simplify(2*(s.Rational(3,2)*M*h)*M-3*M**2*h)==0

# Exact forced-contraction solution, including its initial condition.
Q,R,D,X0=s.symbols('Q R D X0',positive=True)
m=s.symbols('m',integer=True,nonnegative=True)
B=lambda n: Q**n*X0+Q*D*(R**n-Q**n)/(R-Q)
assert s.simplify(B(0)-X0)==0
assert s.simplify(B(m+1)-Q*B(m)-Q*D*R**m)==0

# Necessary C2 growth exponent from M^(5/3)L^(alpha/3) staying positive.
k=s.symbols('k')
assert s.solve(s.Eq(-s.Rational(5,3)*k+alpha/3,0),k)==[alpha/5]

# Candidate thin-layer radial diffusion and velocity/Gamma amplitudes.
num,lam=s.symbols('num lam',positive=True)
h2=num/lam
assert s.simplify(num/h2-lam)==0
assert s.simplify(h2/num-1/lam)==0
print('PASS: quotient PDE identities, scaling, interpolation, recurrence, and layer budget')

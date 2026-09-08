import sympy as s
from itertools import product
r,z,m,k,theta=s.symbols('r z m k theta', real=True)
e,er=s.symbols('e er',real=True)
phase=m*theta+k*r
wr=-m*e*s.sin(phase)/r
wt=-er*s.cos(phase)+k*e*s.sin(phase)
# Reflection: radial component preserved, azimuthal component negated.
assert s.trigsimp(wr.subs(theta,-theta)+wr.subs(k,-k))==0
assert s.trigsimp(-wt.subs(theta,-theta)+wt.subs(k,-k))==0
assert s.trigsimp(wr.subs({m:-m,k:-k}, simultaneous=True)-wr)==0
assert s.trigsimp(wt.subs({m:-m,k:-k}, simultaneous=True)-wt)==0

N=1/(4*s.pi*s.sqrt(r*r+z*z))
K=s.diff(N,z,2)
Qt=-s.diff(K,r)/r
Qrr=-s.diff(K,r,2)
Qrz=-s.diff(K,r,z)
Qtr=15*r*(r*r-6*z*z)/(4*s.pi*(r*r+z*z)**s.Rational(9,2))
assert s.simplify(s.diff(Qt,r)-Qtr)==0
assert s.simplify(Qrr-Qt-r*Qtr)==0
assert s.simplify(Qrz+15*r*z*(4*z*z-3*r*r)/(4*s.pi*(r*r+z*z)**s.Rational(9,2)))==0

# Cylindrical contraction produces the coefficient 3 Qtheta,r.
assert s.simplify(s.diff(Qt,r)+2*(Qrr-Qt)/r-3*Qtr)==0

# Exact split of (w dot grad)v for v=r B e_theta.
B=s.Function('B')(r,z); f=s.Function('f')(r,theta,z)
Gr=B*s.diff(f,r)
Gt=(B/r+s.diff(B,r))*s.diff(f,theta)
assert s.simplify(Gr-s.diff(B*f,r)+s.diff(B,r)*f)==0
assert s.simplify(Gt-s.diff(B*f,theta)/r-s.diff(B,r)*s.diff(f,theta))==0

Cr,Ck,Ce,rho,mm=s.symbols('Cr Ck Ce rho mm',positive=True)
F=Cr/rho+rho*Ck+Ce/(mm**2*rho)
ro=s.sqrt((Cr+Ce/mm**2)/Ck)
assert s.simplify(s.diff(F,rho).subs(rho,ro))==0
assert s.simplify(F.subs(rho,ro)-2*s.sqrt(Ck*(Cr+Ce/mm**2)))==0

# Frequencies 4,12 have no invariant cubic monomial; quadratic mixed terms vanish.
freqs=[4,-4,12,-12]
assert all(sum(t)!=0 for t in product(freqs,repeat=3))
assert all(a+b!=0 for a in [4,-4] for b in [12,-12])
assert all((m//4)%2==1 for m in [4,12])
print('PASS: reflection signs, Q kernels, local E, projection split, carrier optimum, and mode selection')

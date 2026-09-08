#!/usr/bin/env python3
"""Exact potential/operator checks. No evolution or residual norm is certified."""
import json
import sympy as s

r,z=s.symbols('r z',positive=True)
m=s.symbols('m',integer=True)
P=s.Function('P')(r,z);T=s.Function('T')(r,z)

def velocity(P,T,m):
    return [s.diff(P,r,z)+s.I*m*T/r,
            s.I*m*s.diff(P,z)/r-s.diff(T,r),
            -s.diff(P,r,2)-s.diff(P,r)/r+m*m*P/r**2]

u=velocity(P,T,m)
div=s.diff(r*u[0],r)/r+s.I*m*u[1]/r+s.diff(u[2],z)
assert s.simplify(div)==0
omega_z=s.diff(r*u[1],r)/r-s.I*m*u[0]/r
assert s.simplify(omega_z+s.diff(T,r,2)+s.diff(T,r)/r-m*m*T/r**2)==0

# Check axis powers on generic nonquadratic scalar polynomials for signed modes.
for mode in range(-5,6):
    n=abs(mode)
    pp=r**n*(1+r*r+2*r**4)*(1+z+3*z*z)
    tt=r**n*(2+3*r*r+r**6)*(2+z*z)
    ur,ut,uz=velocity(pp,tt,mode)
    for val,power in [(ur+s.I*ut,abs(mode+1)),
                      (ur-s.I*ut,abs(mode-1)),(uz,n)]:
        regular=s.cancel(val/r**power).expand()
        for (exponent,),coefficient in s.Poly(regular,r).terms():
            assert exponent>=0 and exponent%2==0

# Actual initial-profile identities, with arbitrary cubic profile coefficients.
x=s.symbols('x',real=True)
c=s.symbols('c0:4');f=sum(c[j]*x**j for j in range(4))
F=-s.integrate(f,x)  # Any additive constant is a harmless choice here.
ss=r*r+z*z;phi=f.subs(x,ss);phi1=s.diff(f,x).subs(x,ss)
P0=z*F.subs(x,ss)/2
ur,ut,uz=velocity(P0,0,0)
assert s.expand(ur+r*(phi+2*z*z*phi1))==0
assert s.expand(uz-2*z*(phi+r*r*phi1))==0
T0=F.subs(x,ss)/6-ss*phi/3
assert s.expand(-s.diff(T0,r)-r*(phi+2*ss*phi1/3))==0
a,A=s.symbols('a A',positive=True);q=s.Function('q')(z)
Tv=A*a*a*q*F.subs(x,r*r/(a*a))/2
assert s.simplify(-s.diff(Tv,r)-A*r*q*f.subs(x,r*r/(a*a)))==0

# Common-frequency projection and its potential crosswalk.
kap,k=s.symbols('kap k',real=True,nonzero=True)
D=s.Matrix([[kap/2,-kap/2,s.I*k]])
G=s.Matrix([-kap,kap,s.I*k]);W=s.diag(s.Rational(1,2),s.Rational(1,2),1)
K2=kap*kap+k*k
Pi=s.eye(3)+G*D/K2
assert s.simplify((D*G)[0]+K2)==0
assert s.simplify(D*Pi)==s.zeros(1,3)
assert s.simplify(Pi*Pi-Pi)==s.zeros(3)
assert s.simplify(Pi.conjugate().T*W-W*Pi)==s.zeros(3)
pr,pi,tr,ti=s.symbols('pr pi tr ti',real=True)
pc=pr+s.I*pi;tc=tr+s.I*ti
amps=s.Matrix([s.I*kap*(tc-k*pc),s.I*kap*(tc+k*pc),kap*kap*pc])
assert s.simplify(D*amps)==s.zeros(1,1)
energy=(amps.conjugate().T*W*amps)[0]
assert s.expand(energy-kap*kap*(tr*tr+ti*ti+K2*(pr*pr+pi*pi)))==0

# Exact radial integral in the exterior pressure-residual bound.
h,R=s.symbols('h R',positive=True)
for j in range(5):
    tail=h**(-2*j-5)/s.Integer(2*j+5)+2*R*h**(-2*j-6)/s.Integer(2*j+6)\
        +R*R*h**(-2*j-7)/s.Integer(2*j+7)
    assert s.simplify(s.diff(tail,h)+(h+R)**2*h**(-2*j-8))==0

print(json.dumps(dict(status='PASS',checks=[
    'generic cylindrical divergence and vertical-vorticity identities',
    'axis-regular helical powers for signed modes -5 through 5',
    'exact meridional/core-swirl/outer-swirl initial potential identities',
    'common-frequency kinetic-orthogonal Leray projection and potential energy',
    'exterior pressure-residual radial integrals for derivative orders zero through four'],
    scope='Exact algebra only. No evolved coefficient, Galerkin truncation error, finite-energy Bessel sum, or useful-time NS certificate is asserted.'),indent=2))

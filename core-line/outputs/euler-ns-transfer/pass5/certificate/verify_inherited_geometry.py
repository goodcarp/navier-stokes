#!/usr/bin/env python3
"""Exact checks of actual initial swirl-moment and local-strain identities.

These are symbolic differential identities and rational inequalities, not a
Navier--Stokes time integration or a proof of a return map.
"""
import sympy as s

r,z,c,nu,b,A = s.symbols('r z c nu b A', positive=True)
F=s.Function('F')(r,z)
w=r*F
wt=-c*r*s.diff(w,r)+2*c*z*s.diff(w,z)-c*w+nu*(s.diff(w,r,2)+s.diff(w,r)/r+s.diff(w,z,2)-w/r**2)
Ft=-c*r*s.diff(F,r)+2*c*z*s.diff(F,z)-2*c*F+nu*(s.diff(F,r,2)+3*s.diff(F,r)/r+s.diff(F,z,2))
assert s.simplify(wt/r-Ft)==0

# m=r^3 F is the radial/axial angular-momentum density, up to 2*pi.
m=r**3*F
forward=-s.diff(c*r*m,r)-s.diff(-2*c*z*m,z)+nu*(s.diff(m,r,2)-s.diff(3*m/r,r)+s.diff(m,z,2))
assert s.simplify(r**3*Ft-forward)==0

# The same density law holds for the actual nonlinear meridional velocity;
# its sole transport constraint is the three-dimensional divergence identity.
ur,uz=s.Function('ur')(r,z),s.Function('uz')(r,z)
Ft_general=-ur*s.diff(F,r)-uz*s.diff(F,z)-2*ur*F/r+nu*(s.diff(F,r,2)+3*s.diff(F,r)/r+s.diff(F,z,2))
forward_general=-s.diff(ur*m,r)-s.diff(uz*m,z)+nu*(s.diff(m,r,2)-s.diff(3*m/r,r)+s.diff(m,z,2))
divergence=s.diff(ur,r)+ur/r+s.diff(uz,z)
assert s.simplify(r**3*Ft_general-forward_general-r**3*F*divergence)==0

def generator(phi):
    return s.expand(c*r*s.diff(phi,r)-2*c*z*s.diff(phi,z)+nu*(s.diff(phi,r,2)+3*s.diff(phi,r)/r+s.diff(phi,z,2)))
for phi, expected in [(s.Integer(1),0),(z,-2*c*z),(z**2,-4*c*z**2+2*nu),(r**2,2*c*r**2+8*nu)]:
    assert s.simplify(generator(phi)-expected)==0
d,H2,R2=s.symbols('d H2 R2',positive=True)
variance_rate=(-4*c*(H2+d*d)+2*nu)-2*d*(-2*c*d)
assert s.expand(variance_rate-(-4*c*H2+2*nu))==0
aspect_rate=(2*c*R2+8*nu)/(2*R2)-(-2*c*d)/d
assert s.simplify(aspect_rate-(3*c+4*nu/R2))==0
assert s.simplify((-4*c*H2+2*nu)/(2*H2)+2*c-nu/H2)==0

# The moving-axis local gradient satisfies D_t L=-L^2-Hess(p) initially.
pzz=s.symbols('pzz',real=True)
prr=A*A-3*c*c-pzz/2  # Poisson trace for diag(c,c,-2c)+A*J.
L=s.Matrix([[c,-A,0],[A,c,0],[0,0,-2*c]])
H=s.diag(prr,prr,pzz)
Ld=-L*L-H
assert s.simplify(-Ld[2,2]/2-(2*c*c+pzz/2))==0
assert s.simplify(Ld[0,0]-(2*c*c+pzz/2))==0
assert s.simplify((Ld[1,0]-Ld[0,1])/2+2*c*A)==0
assert s.simplify((2*c*c+pzz/2-2*b*c).subs(pzz,4*b*c-4*c*c))==0

# Finite support margins and local normalized-stress margin for the datum.
aa,hh,dd=s.Rational(63,200),s.Rational(9,20),s.Integer(4)
assert dd-hh>s.Rational(5,2)
assert (dd+hh)**2+aa**2<25
assert aa/(dd-hh)<s.Rational(1,4)
gamma=s.Rational(2381,357)
margin=gamma*s.Rational(7,5)-4-s.Rational(901,1000)
assert margin>4
assert 4*s.Rational(7,5)-4*s.Rational(7,5)**2==-s.Rational(56,25)

# Curl of the initial meridional acceleration on the affine pump plateau.
# Pressure cancels by commutation of its mixed derivatives.
W=s.Function('W')(r,z)
ar=-c*c*r+W**2/r
az=-4*c*c*z
assert s.simplify(s.diff(ar,z)-s.diff(az,r)-s.diff(W**2,z)/r)==0
derivative=-s.Rational(16,3)*(aa/(2*hh))*s.sqrt(s.Rational(5,8))
assert s.simplify(derivative+s.Rational(28,15)*s.sqrt(s.Rational(5,8)))==0
assert s.Rational(5,8)>s.Rational(3,4)**2
lower=s.Rational(7,5)*s.Rational(204,35)/s.Rational(279,40000000)
assert lower==s.Rational(108800000,93)>1169892
print('PASS: azimuthal equation and angular-momentum density adjoint.')
print('PASS: mass, center, variance, transverse moment and aspect-ratio rates.')
print('PASS: moving-axis pump/rotation rates and pressure threshold.')
print('PASS: exact initial packet margins and normalized outer-stress rate >',margin)
print('PASS: initial meridional-vorticity generation; coefficient exceeds',lower)
print('Scope: exact initial identities; no quantitative return time or iteration.')

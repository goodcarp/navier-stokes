#!/usr/bin/env python3
"""Exact pressure-source and initial-jet checks for a compact affine extension.

No imported proof code, PDE time evolution, or cone-invariance claim.
"""
import sympy as sp

x, y, z, s, t, b, Om = sp.symbols("x y z s t b Om", real=True)
p, dp, ddp = sp.symbols("p dp ddp", real=True)
coords = (x, y, z)
X = sp.Matrix(coords)
L = sp.Matrix([[-b, -Om, 0], [Om, -b, 0], [0, 0, 2*b]])
radius = x*x+y*y+z*z
q = (X.T*L*X)[0]
U = (p+sp.Rational(2, 3)*radius*dp)*(L*X)-sp.Rational(2, 3)*dp*q*X

def deriv(expr, coordinate):
    return sp.diff(expr, coordinate) + 2*coordinate*(dp*sp.diff(expr,p)+ddp*sp.diff(expr,dp))

Jac = sp.Matrix(3,3,lambda i,j:deriv(U[i],coords[j]))
assert sp.expand(sp.trace(Jac)) == 0
potential = -sp.Rational(1,3)*p*X.cross(L*X)
curl_potential = sp.Matrix([
    deriv(potential[2],y)-deriv(potential[1],z),
    deriv(potential[0],z)-deriv(potential[2],x),
    deriv(potential[1],x)-deriv(potential[0],y),
])
assert (curl_potential-U).applyfunc(sp.expand) == sp.zeros(3,1)
g = sp.expand(sp.trace(Jac*Jac))
assert sp.expand(g.subs({x:0,y:0,z:0,p:1,dp:0,ddp:0})) == 6*b*b-2*Om*Om
assert sp.expand(g-g.subs(Om,-Om)) == 0
assert sp.expand(-y*sp.diff(g,x)+x*sp.diff(g,y)) == 0

gt=0
for (ix,iz), coeff in sp.Poly(g.subs(y,0),x,z).terms():
    assert ix%2==0 and iz%2==0
    gt += coeff*s**((ix+iz)//2)*(1-t)**(ix//2)*t**(iz//2)
gt=sp.expand(gt)
weighted=sp.Poly(sp.expand((3*t-1)*gt),t)
angular=sp.expand(sum(coeff/sp.Integer(2*n[0]+1) for n,coeff in weighted.terms()))
expected_S = sp.Rational(16,105)*s*(-4*s*s*dp*ddp+6*s*p*ddp+2*s*dp*dp+21*p*dp)
expected_W = sp.Rational(16,135)*s*(2*s*ddp+5*dp)*(2*s*dp+3*p)
assert sp.expand(angular-b*b*expected_S-Om*Om*expected_W) == 0

# These are exact integration-by-parts identities under the profile endpoints:
# int p dp=-1/2; int s p ddp=1/2-I; int s^2 dp ddp=-I.
I=sp.symbols("I",nonnegative=True)
radial_S = (-sp.Rational(32,105)*(-I)
            +sp.Rational(16,35)*(sp.Rational(1,2)-I)
            +sp.Rational(16,105)*I
            +sp.Rational(8,5)*(-sp.Rational(1,2)))
radial_W = (sp.Rational(32,135)*(-I)
            +sp.Rational(16,45)*(sp.Rational(1,2)-I)
            +sp.Rational(16,27)*I
            +sp.Rational(8,9)*(-sp.Rational(1,2)))
assert sp.simplify(radial_S) == -sp.Rational(4,7)
assert sp.simplify(radial_W) == -sp.Rational(4,15)
hzz = sp.expand(-(6*b*b-2*Om*Om)/3+b*b*radial_S+Om*Om*radial_W)
assert hzz == -sp.Rational(18,7)*b*b+sp.Rational(2,5)*Om*Om
hperp = sp.expand((-6*b*b+2*Om*Om-hzz)/2)
assert hperp == -sp.Rational(12,7)*b*b+sp.Rational(4,5)*Om*Om
Hcore=sp.diag(hperp,hperp,hzz)
assert sp.simplify(sp.trace(Hcore)+sp.trace(L*L)) == 0

A,Cv,e,f,h,k = sp.symbols("A Cv e f h k",real=True)
E=sp.Matrix([[e,f,h],[f,-e,k],[h,k,0]])
Hv=sp.diag(Cv/2,Cv/2,-Cv)+E
Lp=sp.expand(-L*L-Hcore-A*A*Hv)
bp=A*A*Cv/2-sp.Rational(5,7)*b*b-Om*Om/5
Omp=2*b*Om
expected_Lp=sp.diag(-bp,-bp,2*bp)+Omp*sp.Matrix([[0,-1,0],[1,0,0],[0,0,0]])-A*A*E
assert (Lp-expected_Lp).applyfunc(sp.expand) == sp.zeros(3,3)
beta,R=sp.symbols("beta R",real=True)
beta_prime=sp.expand((bp*Om-b*Omp)/(Om*Om))
assert sp.simplify(beta_prime.subs(b,beta*Om).subs(A*A*Cv,R*Om*Om)
                   -Om*(R/2-sp.Rational(19,7)*beta*beta-sp.Rational(1,5))) == 0

# The full pressure fourth-derivative kernel has no azimuthal/poloidal coupling.
# Differentiating its explicit Cartesian expression avoids a basis shortcut.
N=1/sp.sqrt(radius)  # The harmless constant 1/(4*pi) is omitted here.
Q=sp.Matrix(3,3,lambda i,j:-sp.diff(N,z,z,coords[i],coords[j]))
azimuthal=sp.Matrix([-y,x,0])
radial_component,axial_component=sp.symbols("radial_component axial_component")
poloidal=sp.Matrix([radial_component*x,radial_component*y,axial_component])
assert sp.simplify((azimuthal.T*Q*poloidal)[0]) == 0

print("PASS: vector-potential extension, divergence, symmetry, no b*Omega pressure source.")
print("PASS: exact spherical pressure-source average and both radial cancellations.")
print("Core p_zz =",hzz)
print("Core p_xx = p_yy =",hperp)
print("Initial b' =",bp,"; Omega' =",Omp)
print("PASS: full matrix jet, anisotropy remainder and normalized-strain threshold.")
print("PASS: pressure-kernel azimuthal/poloidal cross contraction vanishes.")
print("All symbolic checks passed. No time evolution or cone invariance was verified.")

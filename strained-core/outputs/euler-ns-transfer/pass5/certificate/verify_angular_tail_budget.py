#!/usr/bin/env python3
"""Exact finite harmonic support of Q times a radial strain extension.

Checks the spherical kernel against Cartesian differentiation and the
vector-harmonic coefficients; no norm quadrature or numerical certificate.
"""
import sympy as s
x,z,R,mu,a,q=s.symbols('x z R mu a q',real=True)
N=1/s.sqrt(x*x+z*z)  # Omits common 1/(4*pi), evaluated in y=0 plane.
Q=s.Matrix([[-s.diff(N,z,z,x,x),-s.diff(N,z,z,x,z)],
            [-s.diff(N,z,z,x,z),-s.diff(N,z,z,z,z)]])
er=s.Matrix([x,z])/s.sqrt(x*x+z*z)
et=s.Matrix([z,-x])/s.sqrt(x*x+z*z)
B=s.Matrix.hstack(er,et)
actual=s.simplify(B.T*Q*B)
rad=x*x+z*z
expected=s.Matrix([[12-36*z*z/rad,-24*x*z/rad],
                   [-24*x*z/rad,21*z*z/rad-9]])/rad**s.Rational(5,2)
assert (actual-expected).applyfunc(s.simplify)==s.zeros(2,2)

sin=s.sqrt(1-mu*mu)
kernel=s.Matrix([[12-36*mu**2,-24*mu*sin],[-24*mu*sin,21*mu**2-9]])
P=s.Matrix([a*(3*mu*mu-1),-q*mu*sin])
F=(kernel*P).applyfunc(s.expand)
rcoef={0:s.Rational(16,5)*(q-3*a),
       2:s.Rational(16,7)*(q-6*a),
       4:-s.Rational(96,35)*(9*a+2*q)}
vcoef={2:s.Rational(16,7)*a,4:(144*a+42*q)/35}
assert s.expand(F[0]-sum(c*s.legendre(l,mu) for l,c in rcoef.items()))==0
assert s.simplify(F[1]+sin*sum(c*s.diff(s.legendre(l,mu),mu) for l,c in vcoef.items()))==0

# Exact gradient-harmonic orthogonality, including the radial l=0 channel.
for ell in [0,2,4,6,8,10,12]:
 radial=s.integrate(s.expand(F[0]*s.legendre(ell,mu)),(mu,-1,1))
 tangent=s.integrate(s.expand(s.simplify(-F[1]*sin)*s.diff(s.legendre(ell,mu),mu)),(mu,-1,1))
 if ell>4:
  assert s.simplify(radial)==0 and s.simplify(tangent)==0
 else:
  assert s.simplify(radial-2*rcoef.get(ell,0)/(2*ell+1))==0
  assert s.simplify(tangent-2*ell*(ell+1)*vcoef.get(ell,0)/(2*ell+1))==0

# Spherical Sobolev identities for the axisymmetric tangent basis grad_S P_l.
for ell in range(1,9):
 Pl=s.legendre(ell,mu)
 D=s.diff((1-mu*mu)*s.diff(Pl,mu),mu)
 assert s.expand(D+ell*(ell+1)*Pl)==0
 normP=s.integrate(Pl*Pl,(mu,-1,1))
 normV=s.integrate(s.expand((1-mu*mu)*s.diff(Pl,mu)**2),(mu,-1,1))
 assert s.simplify(normV-ell*(ell+1)*normP)==0

ell=s.symbols('ell',nonnegative=True)
assert s.expand(ell*(ell+1)+s.Rational(1,4)-(ell+s.Rational(1,2))**2)==0
print('PASS: Cartesian fourth-derivative kernel -> spherical block.')
print('PASS: Q S_eta has radial degrees 0,2,4 and tangent degrees 2,4 only.')
print('PASS: explicit coefficient contraction and angular Sobolev weights.')
print('No interval norms, full-gate sign, or subsequent dynamics were certified.')

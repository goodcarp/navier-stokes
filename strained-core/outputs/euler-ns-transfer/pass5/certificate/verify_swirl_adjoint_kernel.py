#!/usr/bin/env python3
"""Exact whole-space adjoint ODE and plateau swirl-kernel checks.

No numerical pressure solve, interval quadrature, or gate-sign assertion.
"""
import sympy as s
R=s.symbols('R',positive=True)
mu=s.symbols('mu',real=True)
E,Ep,Epp,Eppp,M,H,J=s.symbols('E Ep Epp Eppp M H J',real=True)
# M=int_0^R t E(t) dt; H=int_R^infty E(t)t^-8 dt;
# J=int_R^infty E'(t)t^-3 dt.
def D(f):
 return (s.diff(f,R)+Ep*s.diff(f,E)+Epp*s.diff(f,Ep)+Eppp*s.diff(f,Epp)
         +R*E*s.diff(f,M)-E/R**8*s.diff(f,H)-Ep/R**3*s.diff(f,J))
F={0:s.Rational(16,5)*(Epp/R**3-Ep/R**4),
   2:s.Rational(16,7)*Epp/R**3-s.Rational(64,7)*Ep/R**4,
   4:-s.Rational(192,35)*Epp/R**3-s.Rational(2088,35)*Ep/R**4-72*E/R**5}
psi={0:s.Rational(16,5)*J,
     2:-s.Rational(16,7)*E/R**3,
     4:s.Rational(192,35)*E/R**3+s.Rational(8,5)*M/R**5-80*R**4*H}
for ell in [0,2,4]:
 assert s.simplify(-D(D(psi[ell]))-2*D(psi[ell])/R+ell*(ell+1)*psi[ell]/R**2-F[ell])==0

# Conversion of q_l/R^5 from eta(s), eta'(s), eta''(s).
p,p1,p2=s.symbols('p p1 p2')
q={0:s.Rational(64,5)*R**4*p2,
   2:-s.Rational(32,7)*R**2*(3*p1-2*R**2*p2),
   4:-s.Rational(24,35)*(105*p+190*R**2*p1+32*R**4*p2)}
sub={p:E,p1:Ep/(2*R),p2:(Epp-Ep/R)/(4*R**2)}
for ell in [0,2,4]: assert s.simplify(q[ell].subs(sub)/R**5-F[ell])==0

# Exact radial Green symmetry, with the radial physical measure retained.
pfun=s.Function('p')(R);vfun=s.Function('v')(R);lam=s.symbols('lam')
L=lambda f:-s.diff(f,R,2)-2*s.diff(f,R)/R+lam*f/R**2
boundary=R**2*(pfun*s.diff(vfun,R)-vfun*s.diff(pfun,R))
assert s.simplify(R**2*(vfun*L(pfun)-pfun*L(vfun))-s.diff(boundary,R))==0

# Radial IBP for swirl source; only the compact-support boundary is removed.
W=s.Function('W')(R);Z=s.Function('Z')(R)
g=-s.diff(W,R)/R-(W+Z)/R**2
assert s.simplify(R**2*vfun*g-(R*s.diff(vfun,R)*W-vfun*Z)+s.diff(R*vfun*W,R))==0

# Angular integration-by-parts endpoint is zero because a smooth swirl
# satisfies w^2=(1-mu^2) times a smooth function at the axis.
for ell in [0,2,4]:
 Pl=s.legendre(ell,mu)
 for power in [0,1,2,3]:
  wsq=(1-mu*mu)*mu**(2*power)
  assert s.integrate(s.expand(mu*s.diff(wsq,mu)*Pl+wsq*(Pl+mu*s.diff(Pl,mu))),(mu,-1,1))==0

# All three plateau simplifications and the final derivative-free kernel.
a,b=s.symbols('a b',real=True)
plateau={0:s.symbols('constant_psi0'),2:s.Rational(16,7)/R**3,
         4:s.Rational(36,7)/R**3+a*R**4+b/R**5}
Mplateau=-R**2/2+s.Rational(5,8)*b
Hplateau=-1/(7*R**7)-a/80
assert s.simplify(psi[4].subs({E:-1,M:Mplateau,H:Hplateau})-plateau[4])==0
kernel=sum(plateau[ell]*mu*s.diff(s.legendre(ell,mu),mu)
           -R*s.diff(plateau[ell],R)*s.legendre(ell,mu) for ell in [0,2,4])
expected=(s.Rational(3,14)*(735*mu**4-370*mu**2+11)/R**3
          +s.Rational(3,2)*a*R**4*(5*mu**2-1)
          +s.Rational(15,8)*b*(21*mu**4-14*mu**2+1)/R**5)
assert s.simplify(kernel-expected)==0
print('PASS: all three exact adjoint radial Poisson equations and source constants.')
print('PASS: radial Green symmetry, swirl-source radial/angular IBP signs.')
print('PASS: plateau adjoints and derivative-free w^2 kernel coefficients.')
print('No interval quadrature, numerical gate sign, or later dynamics was certified.')

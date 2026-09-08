#!/usr/bin/env python3
"""Exact checks for outer-pressure-reservoir.md; not a return-map proof."""
import sympy as s
r,z,c,sigma=s.symbols('r z c sigma',positive=True)
t=s.symbols('t',nonnegative=True)
Q=3*(4*z*z-r*r)/(4*s.pi*(r*r+z*z)**s.Rational(7,2))
transport=s.factor((r*s.diff(Q,r)-2*z*s.diff(Q,z)-2*Q)/Q)
expected=8-6*t/(4-t)-21*t/(1+t)
assert s.simplify(transport-expected.subs(t,r*r/(z*z)))==0
gamma=s.Rational(2381,357)
assert gamma==8-s.Rational(6,63)-s.Rational(21,17)
assert gamma>6
# Each subtracted fraction is increasing on the interval, so endpoint bounds suffice.
assert s.simplify(s.diff(6*t/(4-t),t))==24/(t-4)**2
assert s.simplify(s.diff(21*t/(1+t),t))==21/(t+1)**2
assert s.simplify(expected.subs(t,s.Rational(1,16))-gamma)==0

alpha=s.Function('alpha'); beta=s.Function('beta')
Psi=-c*z*r*r*alpha(r*r/(sigma*sigma))*beta(z)
Pr=-s.diff(Psi,z)/r
Pz=s.diff(Psi,r)/r
assert s.simplify(s.diff(r*Pr,r)/r+s.diff(Pz,z))==0
assert s.simplify(Pr-c*r*alpha(r*r/sigma**2)*(beta(z)+z*s.diff(beta(z),z)))==0
arg=s.symbols('arg')
aprime=s.Subs(s.diff(alpha(arg),arg),arg,r*r/sigma**2)
assert s.simplify(Pz+2*c*z*(alpha(r*r/sigma**2)+(r*r/sigma**2)*aprime)*beta(z))==0
# The axial kernel, vertical velocity square and angular volume give coefficient -48.
assert -s.Rational(24,1)/(4*s.pi)*4*2*s.pi==-48

b,Omega,Cv,CP,D,nu=s.symbols('b Omega Cv CP D nu',positive=True)
# Treat CP as a signed symbol when checking the tuning identity.
CP_signed=s.symbols('CP_signed',real=True)
Asq=(s.Rational(38,7)*b*b+s.Rational(2,5)*Omega*Omega-c*c*CP_signed)/Cv
Ctotal=s.expand(Asq*Cv+c*c*CP_signed)
bprime=Ctotal/2-s.Rational(5,7)*b*b-Omega*Omega/5
assert s.simplify(bprime-2*b*b)==0
etaprime=s.simplify(bprime/Omega-b*(2*b*Omega)/(Omega*Omega))
assert etaprime==0
# For CP<0 all terms of the numerator are positive.
assert s.expand(Asq.subs(CP_signed,-CP)*Cv)==s.Rational(38,7)*b*b+s.Rational(2,5)*Omega*Omega+c*c*CP
# Full-pressure neutral-boundary test. This checks the algebra; p_zz' is unknown.
p,pdot=s.symbols('p pdot',real=True)
beta_prime=-(p+8*b*b)/(2*Omega)
beta_second=(s.diff(beta_prime,b)*2*b*b
             +s.diff(beta_prime,Omega)*2*b*Omega
             +s.diff(beta_prime,p)*pdot).subs(p,-8*b*b)
assert s.simplify(beta_second+(pdot+32*b**3)/(2*Omega))==0
J=s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
L=s.diag(-b,-b,2*b)+Omega*J
g1=2*s.trace(L*(2*b*L))
assert s.expand(g1)==24*b**3-8*b*Omega**2
assert s.expand(-g1/3)==-8*b**3+s.Rational(8,3)*b*Omega**2
print('PASS: exact swirl transport factor, cone lower bound, divergence-free pump, pressure-limit coefficient and tuned initial core ratios')
print('PASS: full-pressure neutral ratio condition and distributional contact coefficient')
print('SCOPE: initial identities only; full meridional-pressure response and cone invariance are unproved')

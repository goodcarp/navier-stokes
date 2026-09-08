#!/usr/bin/env python3
"""Exact algebra checks supporting active-core-pressure.md; not a PDE proof."""
import sympy as s

x,y,z,Omega=s.symbols('x y z Omega', real=True)
r2=x*x+y*y+z*z
psi=s.Function('psi')
P=psi(r2)
coords=(x,y,z)
c=s.Matrix([-Omega*y*P,Omega*x*P,0])
div=lambda v: sum(s.diff(v[i],coords[i]) for i in range(3))
adv=s.Matrix([sum(c[j]*s.diff(c[i],coords[j]) for j in range(3)) for i in range(3)])
assert s.simplify(div(c)) == 0
assert s.simplify(adv-s.Matrix([-Omega**2*x*P**2,-Omega**2*y*P**2,0])) == s.zeros(3,1)
source=-Omega**2*(2*P**2+4*(x*x+y*y)*P*s.Subs(s.Derivative(psi(s.Symbol('q')),s.Symbol('q')),s.Symbol('q'),r2))
assert s.simplify(div(adv)-source)==0
mu=s.symbols('mu',real=True)
ang=s.integrate((3*mu**2-1)*(1-mu**2),(mu,-1,1))/2
assert ang == -s.Rational(4,15)
rad=s.Rational(-1,4) # fundamental theorem and psi(0)=1, psi(infinity)=0
core=s.Rational(2,3)-4*ang*rad
assert core==s.Rational(2,5)

# Fourth-derivative horizontal tensor before the common 1/(4 pi) factor.
N=r2**s.Rational(-1,2)
K=s.Matrix(2,2,lambda i,j:s.diff(N,z,z,coords[i],coords[j]))
yp=s.Matrix([x,y])
expected=((105*z*z-15*r2)*(yp*yp.T)+(3*r2*r2-15*r2*z*z)*s.eye(2))/r2**s.Rational(9,2)
assert s.simplify(K-expected)==s.zeros(2)
t=s.symbols('t',nonnegative=True)
radial=(-12+81*t-12*t*t)/(1+t)**2
tangent=(-12-9*t+3*t*t)/(1+t)**2
assert s.simplify(radial-tangent-15*t*(6-t)/(t+1)**2)==0
# On [0,1/16], the first upper-bound numerator increases and remains negative.
n=s.expand(s.cancel((radial+6)*(1+t)**2))
assert n == -6*t*t+93*t-6
assert s.diff(n,t).subs(t,s.Rational(1,16))>0
assert n.subs(t,s.Rational(1,16))<0
assert s.Rational(1,16)<6 # radial is the larger eigenvalue throughout the interval

# Explicit packet is divergence free for any smooth scalar envelope.
chi=s.Function('chi')(x,y,z)
v=s.Matrix([s.diff(chi,y),-s.diff(chi,x),0])
assert div(v)==0

# Rigid-core vorticity and its initial local RHS.
u=s.Matrix([-Omega*y,Omega*x,0])
curl=lambda a:s.Matrix([s.diff(a[2],y)-s.diff(a[1],z),s.diff(a[0],z)-s.diff(a[2],x),s.diff(a[1],x)-s.diff(a[0],y)])
om=curl(u)
assert om==s.Matrix([0,0,2*Omega])
rhs=s.Matrix([sum(om[j]*s.diff(u[i],coords[j])-u[j]*s.diff(om[i],coords[j]) for j in range(3)) for i in range(3)])
assert rhs==s.zeros(3,1)
print('PASS: rotating-core pressure source, contact term/angular coefficient, fourth-derivative tensor, cone negativity, divergence, and initial vorticity identities')
print('SCOPE: exact algebra plus the stated analytic local-existence argument; no numerical evolution or regeneration theorem')

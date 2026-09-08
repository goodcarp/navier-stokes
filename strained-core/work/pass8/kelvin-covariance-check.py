"""Exact Kelvin polarization and stress-budget checks; no compact-flow claim."""
import sympy as s
from fractions import Fraction as F

t,c,Om,k,h,nu,lam=s.symbols('t c Omega k h nu lambda',positive=True)
S=2*Om*h
K=k+S*t; Q=K*K+h*h; Q0=k*k+h*h
J0=(1-s.exp(-2*c*t))/(2*c)
J1=(1-s.exp(-2*c*t)*(1+2*c*t))/(4*c*c)
J2=(1-s.exp(-2*c*t)*(1+2*c*t+2*c*c*t*t))/(4*c**3)
I=Q0*J0+2*k*S*J1+S*S*J2
assert s.simplify(I.subs(t,0))==0
assert s.simplify(s.diff(I,t)-s.exp(-2*c*t)*Q)==0

# Represent the viscous factor by F(t), substituting its exact derivative.
Ft=s.Function('F')(t)
Fdot=-nu*s.exp(-2*c*t)*Q*Ft
derivative=lambda expr:s.diff(expr,t).subs(s.diff(Ft,t),Fdot)
kap=s.exp(-c*t)*s.Matrix([K,h])
aa=lam*Q0*s.exp(-c*t)*Ft*s.Matrix([-h,K])/Q
A=s.Matrix([[c,-Om],[-Om,c]])
J=s.Matrix([[0,-1],[1,0]])
assert s.simplify(kap.dot(aa))==0
assert s.simplify(kap.diff(t)+(A+Om*J)*kap)==s.zeros(2,1)
rhs=-(A+Om*J)*aa+2*kap*(kap.dot(A*aa))/kap.dot(kap)-nu*kap.dot(kap)*aa
assert s.simplify(aa.applyfunc(derivative)-rhs)==s.zeros(2,1)

sigma=-aa[0]*aa[1]/2
energy=aa.dot(aa)/4
assert s.simplify(derivative(energy)+2*c*energy+2*Om*sigma+2*nu*kap.dot(kap)*energy)==0
logstress=-2*c-2*nu*kap.dot(kap)+S*(h*h-3*K*K)/(K*Q)
assert s.simplify(derivative(sigma)-sigma*logstress)==0
inviscid_sigma=sigma.subs({Ft:1,c:0})
budget=lam*lam*Q0/(8*Om)*(1-Q0/Q)
assert s.simplify(s.diff(budget,t)-inviscid_sigma)==0

hlo=F(1000,63); hhi=F(1600,63); kk=20
assert hhi*hhi<3*kk*kk
ratio=2*((kk*kk+hhi*hhi)/(4*kk*kk+hhi*hhi))**2
assert ratio<F(7,16)
assert (kk*kk+hhi*hhi)/(4*kk*kk+hhi*hhi)<F(1,2)
# h+k^2/h is convex, so its maximum on the closed interval is at an endpoint.
endmax=max(hlo+kk*kk/hlo,hhi+kk*kk/hhi)
assert endmax/(4*390*kk)==F(10369,7862400)<F(33,25000)

# Full-datum initial fluctuation-energy drain from the exact envelope plateaus.
drain=2*80*780*F(63,800)*F(9,10)/(F(6764,75)*F(12,5))
assert drain==F(552825,13528)
assert drain+F(14,5)==F(2953517,67640)>F(4366,100)
print('PASS: exact rotating Kelvin equations, polarization, viscosity, stress depletion, and integrated budget')
print('PASS: actual initial fluctuation-energy constants; no positive-time propagation asserted')

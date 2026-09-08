#!/usr/bin/env python3
"""Exact entropy and clock controls; analytical theorems are not executable.

No other implementation is imported. The translated Gaussian is an
auxiliary scalar-kernel test, not finite-energy Navier--Stokes data.
"""
from pathlib import Path
import hashlib
import sympy as s

print('SOURCE_SHA256',hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
nu,t,L,Y0,C=s.symbols('nu t L Y0 C',positive=True)
mu,beta=s.symbols('mu beta',real=True)
r=s.symbols('r',positive=True)

# A mass-one heat density has variance 2nu*t in each of three axes.
# Its displacement from the reference center is mu in the first axis.
variance=nu*t/L**2
D=mu**2/(4*L**2)+s.Rational(3,2)*(variance-1-s.log(variance))
entropy_upper=s.Rational(3,2)*s.log(L**2/(nu*t))+(mu**2+6*nu*t)/(4*L**2)
assert s.simplify(entropy_upper-D-s.Rational(3,2))==0
tilt_gap=s.expand(D+beta**2*L**2-beta*mu)
positive_gap=(mu/(2*L)-beta*L)**2+s.Rational(3,2)*(variance-1-s.log(variance))
assert s.simplify(tilt_gap-positive_gap)==0
assert s.diff(r-1-s.log(r),r)==1-1/r
assert (r-1-s.log(r)).subs(r,1)==0
# Removing the displacement from entropy makes the tilted bound false.
wrong_tilt_gap=(beta**2*L**2-beta*mu).subs({L:1,mu:2,beta:1})
assert wrong_tilt_gap==-1
print('GAUSSIAN_ENTROPY',D)
print('DROPPED_DISPLACEMENT_CONTROL',wrong_tilt_gap)

# Nash differential inequality has the required t^(-3/2) L2-square power.
Y=(Y0**(-s.Rational(2,3))+s.Rational(4,3)*C*nu*t)**(-s.Rational(3,2))
assert s.simplify(s.diff(Y,t)+2*C*nu*Y**s.Rational(5,3))==0
print('NASH_ODE exact; L1-to-L2 power -3/4 and composed power -3/2')

# Time integral of the logarithm: first branch 0<T<=b^2.
T,b=s.symbols('T b',positive=True)
F=T*(s.log(b/s.sqrt(T))+s.Rational(1,2))
assert s.simplify(s.diff(F,T)-s.log(b/s.sqrt(T)))==0
assert s.limit(F,T,0,dir='+')==0
assert s.simplify(F.subs(T,b*b)-b*b/2)==0
clock_majorant=T*(1+s.log(b/s.sqrt(T)))
assert s.simplify(s.diff(clock_majorant,T)-(s.Rational(1,2)+s.log(b/s.sqrt(T))))==0
assert s.integrate(T*T,(T,0,t))==t**3/3

# Quantitative smallness can be made uniform over the Reynolds number.
# If K=100*C>=100 and c=1/[K(1+log K)], then
# log(1+log K)<=1+log K proves C*c*(1+log(1/c))<=2/100.
K=s.symbols('K',positive=True)
c=1/(K*(1+s.log(K)))
assert s.simplify((K/100)*c*2*(1+s.log(K))-s.Rational(2,100))==0
assert 2*(s.Rational(2,100)+s.Rational(1,100)) < s.Rational(1,2)
print('UNIFORM_BOOTSTRAP_SMALLNESS elementary logarithm slack passed')

# Original NS scaling: E has exponent -1, M exponent +2, nu fixed.
ell_power=-s.Rational(1,5)-2*s.Rational(2,5)
assert ell_power==-1
Re_power=-s.Rational(2,5)+2*s.Rational(1,5)
assert Re_power==0
lam,E,M=s.symbols('lam E M',positive=True)
energy_Re=E**s.Rational(2,5)*M**s.Rational(1,5)/nu
scaled_Re=energy_Re.subs({E:E/lam,M:lam**2*M},simultaneous=True)
assert s.simplify(scaled_Re-energy_Re)==0
log_duration=1/(M*(1+s.log(energy_Re)))
scaled_duration=log_duration.subs({E:E/lam,M:lam**2*M},simultaneous=True)
assert s.simplify(scaled_duration/log_duration-lam**-2)==0
assert (-1)-(-s.Rational(11,10))==s.Rational(1,10)
R=s.symbols('R',positive=True)
ratio=R/(1+s.log(R))
assert s.simplify(s.diff(ratio,R)-s.log(R)/(1+s.log(R))**2)==0
assert s.limit(ratio,R,s.oo)==s.oo
x=s.symbols('x',positive=True)
assert s.limit(s.exp(x/10)/(1+x/5),x,s.oo)==s.oo
j=s.symbols('j',integer=True,nonnegative=True)
assert s.summation(s.Rational(1,2)**j,(j,0,s.oo))==2
print('CLOCK_SCALING -2; energy Reynolds number invariant')
print('LARGE_RECORD_PARABOLIC_RATIO grows; geometric delay lower bounds summable')
print('ALL EXACT CONTROLS PASSED; NO GLOBAL NAVIER_STOKES CONCLUSION')

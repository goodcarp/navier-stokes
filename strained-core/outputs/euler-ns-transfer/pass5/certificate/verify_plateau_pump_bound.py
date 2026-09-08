#!/usr/bin/env python3
"""Exact local-plus-pressure swirl multiplier and rational support bound."""
import sympy as s
rho,z=s.symbols('rho z',positive=True)
R=s.symbols('R',positive=True)
t,a,B=s.symbols('t a B',real=True)
r2=rho*rho+z*z
# Common 1/(4*pi) omitted from every tensor/kernel component.
Qtheta=(15*z*z-3*r2)/r2**s.Rational(7,2)
Qrr=(12*r2*r2-105*rho*rho*z*z)/r2**s.Rational(9,2)
Qrz=rho*z*(45*r2-105*z*z)/r2**s.Rational(9,2)
# P=(rho,0,-2z) on the eta=-1 plateau. This includes both the
# centripetal poloidal acceleration and the swirl-advection coupling.
Klocal=(-2*Qrr+4*z*Qrz/rho-rho*s.diff(Qtheta,rho)
        +2*z*s.diff(Qtheta,z)+2*Qtheta)
H=(21*t*t-14*t+1)/(5*t-1)
assert s.simplify(Klocal/Qtheta-(-15*H).subs(t,z*z/r2))==0

# Convert the dR dmu pressure kernel to physical dx before normalizing.
mu=s.symbols('mu',real=True)
Kmeasure=(s.Rational(3,14)*(735*mu**4-370*mu*mu+11)/R**3
          +s.Rational(3,2)*a*R**4*(5*mu*mu-1)
          +s.Rational(15,8)*B*(21*mu**4-14*mu*mu+1)/R**5)
# (1/(2*pi*R^2)) divided by Qtheta=3(5mu^2-1)/(4pi R^5).
pressure_ratio=s.expand(2*R**3*Kmeasure/(3*(5*mu*mu-1)))
expected_pressure=((735*t*t-370*t+11)/(7*(5*t-1))
                   +a*R**7+5*B*H/(4*R*R))
assert s.simplify(pressure_ratio-expected_pressure.subs(t,mu*mu))==0
combined=s.factor(expected_pressure-15*H)
simple=s.Rational(24,7)+(5*B/(4*R*R)-10)*H+a*R**7
assert s.simplify(combined-simple)==0
assert s.simplify(s.diff(H,t)-3*(35*t*t-14*t+3)/(5*t-1)**2)==0
assert s.expand(35*t*t-14*t+3-(35*(t-s.Rational(1,5))**2+s.Rational(8,5)))==0

# Exact support and join constants.
tmin=s.Rational(710**2,710**2+63**2)
Rmin=s.Rational(71,20)
Rmax2=s.Rational(89,20)**2+s.Rational(63,200)**2
amax=-s.Rational(80,7*6**7)
Bmax=s.Integer(5)
assert Rmin>s.Rational(5,2) and Rmax2<25
assert tmin>s.Rational(16,17)
assert H.subs(t,tmin)>0
assert 5*Bmax/(4*Rmin**2)-10<0
# H is increasing on the support; the envelope decreases with t.
# Its R derivative is 7*amax*R^6 - (5*Bmax/(2*R^3))*H < 0.
assert amax<0 and Bmax>0
upper=s.factor(simple.subs({a:amax,B:Bmax,R:Rmin,t:tmin}))
assert upper == -s.Rational(358640696595957938134855826987,23085491673405775693824000000)
assert upper < -s.Rational(31,2)
print('PASS: local centripetal and swirl-advection kernel, including signs.')
print('PASS: pressure-kernel measure conversion and combined multiplier.')
print('PASS: monotonic upper envelope and exact packet/plateau support.')
print('Exact multiplier upper bound:',upper)
print('Decimal display:',float(upper),'; strictly less than -31/2.')
print('This bounds the combined c*A^2 annular-swirl coefficient, not the full gate alone.')

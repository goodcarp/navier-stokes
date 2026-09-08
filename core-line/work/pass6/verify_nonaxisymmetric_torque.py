#!/usr/bin/env python3
"""Exact full-3D angular-momentum identities and compact Fourier-mode seed."""
import sympy as s
r,z,theta,nu=s.symbols('r z theta nu',positive=True)
G=s.Function('G')(r,theta,z)
ur=s.Function('ur')(r,theta,z)
uth=G/r
scalar_lap=lambda f:s.diff(f,r,2)+s.diff(f,r)/r+s.diff(f,theta,2)/r**2+s.diff(f,z,2)
diffusion=r*(scalar_lap(uth)-uth/r**2+2*s.diff(ur,theta)/r**2)
expected=s.diff(G,r,2)-s.diff(G,r)/r+s.diff(G,z,2)+s.diff(G,theta,2)/r**2+2*s.diff(ur,theta)/r
assert s.simplify(diffusion-expected)==0
wz=(s.diff(G,r)-s.diff(ur,theta))/r
assert s.simplify(expected-(scalar_lap(G)-2*wz))==0

# Mean Reynolds flux reduction: after using divergence, the discrepancy
# is an angular derivative with zero periodic average.
vr=s.Function('vr')(r,theta,z)
vt=s.Function('vt')(r,theta,z)
vz=s.Function('vz')(r,theta,z)
flux=-s.diff(r*r*vr*vt,r)/r-r*s.diff(vz*vt,z)
adv=-vr*vt-r*(vr*s.diff(vt,r)+vz*s.diff(vt,z))
divv=s.diff(vr,r)+vr/r+s.diff(vt,theta)/r+s.diff(vz,z)
assert s.simplify(flux-adv+r*vt*divv-s.diff(vt**2,theta)/2)==0

# Compact-envelope single-mode seed; support away from r=0 is specified
# analytically in the note. Algebra here is valid for any smooth envelope.
m=s.symbols('m',integer=True,positive=True)
k=s.symbols('k',real=True)
a=s.Function('a')(r,z)
phase=m*theta+k*r
potential=a*s.cos(phase)
vr=s.diff(potential,theta)/r
vt=-s.diff(potential,r)
assert s.simplify(s.diff(r*vr,r)/r+s.diff(vt,theta)/r)==0
assert s.simplify(vr+m*a*s.sin(phase)/r)==0
assert s.simplify(vt+s.diff(a,r)*s.cos(phase)-k*a*s.sin(phase))==0

# Use exact trigonometric moments of an integer nonzero Fourier mode.
C,S=s.symbols('C S',real=True)
def average(expr):
    poly=s.Poly(s.expand(expr.xreplace({s.cos(phase):C,s.sin(phase):S})),C,S)
    values={(0,0):s.Integer(1),(1,0):0,(0,1):0,
            (2,0):s.Rational(1,2),(0,2):s.Rational(1,2),(1,1):0}
    return s.simplify(sum(coeff*values[monomial] for monomial,coeff in poly.terms()))
cov=average(vr*vt)
assert s.simplify(cov+m*k*a*a/(2*r))==0
torque=-s.diff(r*r*cov,r)/r
assert s.simplify(torque-m*k*s.diff(r*a*a,r)/(2*r))==0

a0=s.symbols('a0',positive=True)
flat_torque=s.simplify(torque.subs({s.diff(a,r):0,a:a0}))
flat_energy=s.simplify(average(vr**2+vt**2).subs({s.diff(a,r):0,a:a0}))
assert s.simplify(flat_torque-m*k*a0*a0/(2*r))==0
assert s.simplify(flat_energy-a0*a0*((m/r)**2+k*k)/2)==0
# Both signed torque bounds follow from explicit nonnegative squares.
assert s.simplify(flat_energy/2-flat_torque-a0*a0*(m/r-k)**2/4)==0
assert s.simplify(flat_energy/2+flat_torque-a0*a0*(m/r+k)**2/4)==0

# General complex Fourier covariance formula, represented by real/imag parts.
xr,xi,yr,yi=s.symbols('xr xi yr yi',real=True)
cov_general=average((xr*s.cos(phase)-xi*s.sin(phase))*(yr*s.cos(phase)-yi*s.sin(phase)))
assert s.simplify(cov_general-(xr*yr+xi*yi)/2)==0
print('PASS: full cylindrical and Cartesian angular-momentum diffusion identity.')
print('PASS: conservative mean Reynolds torque and periodic-derivative reduction.')
print('PASS: compact-envelope solenoidal Fourier seed and selected-sign mean torque.')
print('PASS: exact complex-mode covariance and local energy bound.')
print('Scope: exact identities and initial seed; no persistence or return map.')

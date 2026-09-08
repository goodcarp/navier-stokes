#!/usr/bin/env python3
"""Interval cutoff jets and exact rational bounds for the outer torque target."""
import sys,json
from pathlib import Path
from fractions import Fraction as F
from math import factorial
here=Path(__file__).resolve().parent
for base in (here,here.parent,here.parent.parent):
    for candidate in (base/'pass5',base/'pass5'/'certificate',base):
        if (candidate/'interval_local_integrals.py').exists():
            sys.path.insert(0,str(candidate))
from interval_local_integrals import iv,interval,cutoff_jets,hull,exact_endpoints

iv.dps=35
n=64
jets=[cutoff_jets(interval(F(1,4)+F(3*j,4*n),
                           F(1,4)+F(3*(j+1),4*n)),2) for j in range(n)]
d1=hull([entry[1] for entry in jets]); d2=hull([entry[2] for entry in jets])
assert bool(d1.a>-4) and bool(d1.b<=0)
assert bool(d2.a>-32) and bool(d2.b<32)

a=F(63,200); r0=F(21,100); width=F(7,50)
rmin,rmax=r0-width,r0+width
zmin,zmax=F(17,5),F(23,5)
assert rmin==F(7,100) and rmax==F(7,20)
assert zmin>F(5,2) and zmax*zmax+rmax*rmax<25
assert rmax/zmin<F(1,4)
assert r0-width/2<a/2
assert F(5,8)<F(4,5)**2
assert F(4,5)*a<r0+width/2
# For xi>=5/8, the analytic derivative factor is <=1-(16/3)xi<0.
assert 1-F(16,3)*F(5,8)<0
curvature=4*F(5,8)*(2*4+F(5,8)*32)
assert curvature==70
torque=40/(F(4,5)*a)
assert torque==F(10000,63)

# exp(5/3)>5 by its positive Taylor series, hence log5<5/3.
assert sum((F(5,3)**j)/factorial(j) for j in range(7))>5
radial_derivative_bound=rmax*(2*4/width)*2
assert radial_derivative_bound==40
phase_energy=400*(rmax*rmax-rmin*rmin)/2
assert phase_energy==F(588,25)
radial_energy_bound=F(80,3)+phase_energy+radial_derivative_bound
cw_upper=3*(F(5,17)**5)*F(12,5)*radial_energy_bound
assert cw_upper<F(3,2)
# Positive horizontal Q on the cone: numerator Qrr =12-81t+12t².
assert 12-F(81,16)>0

H=F(204,35); cv_lower=F(113,20000000); cv_upper=F(279,40000000)
assert H-F(3,2)>0
assert (H-F(3,2))/cv_upper>780**2
assert H/cv_lower<1020**2
assert 780*a*a/4>19
nu=F(1,1000)
gain_lower=F(3,4)**2*torque-70*nu*1020
gain_one=torque-70*nu*1020
assert gain_lower==F(626,35)>0
assert gain_one==F(27509,315)>87
result=dict(status='PASS',scope='initial mean maximum and neutral pressure only; full pressure-derivative gate untested',
    cutoff_jets=dict(panels=n,dps=iv.dps,first_derivative_display=str(d1),
                     second_derivative_display=str(d2),
                     first_derivative_exact=exact_endpoints(d1),
                     second_derivative_exact=exact_endpoints(d2)),
    seed_pressure_upper=str(cw_upper),seed_amplitude_interval=['3/4','1'],
    retuned_swirl_amplitude_bounds=['780','1020'],
    mean_maximum_derivative_lower=str(gain_lower),
    unit_seed_derivative_lower=str(gain_one))
print(json.dumps(result,indent=2))

#!/usr/bin/env python3
"""Independent algebra and exact endpoint replay for the new leading family."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as s

here=Path(__file__).resolve().parent
envelope=json.loads((here/"leading-envelope-certificate.json").read_text())
pressure=json.loads((here/"leading-pressure-certificate.json").read_text())
assert envelope["status"]==pressure["status"]=="PASS"
def dyadic(x):
    return (-1 if x["sign"] else 1)*F(int(x["mantissa"]))*F(2)**x["exponent"]
assert dyadic(envelope["unit_mean_torque"]["exact"]["lower"])>1950
assert dyadic(envelope["R1"]["exact_dyadic_bounds"]["upper"])<200000
assert dyadic(envelope["radial_extraction"]["exact_dyadic_bounds"]["lower"])>F(1,50)
assert dyadic(pressure["E_absolute_upper_expression"]["exact"]["upper"])<F(1,3)

# Covariance with a non-flat envelope and either chirality.
r,phase,m,k,lam=s.symbols("r phase m k lam",real=True,positive=True)
e=s.Function("e")(r)
sn,co=s.sin(phase),s.cos(phase)
wr=-lam*m*e*sn/r
wt=-lam*s.diff(e,r)*co+lam*k*e*sn
cov=s.expand(wr*wt).subs({sn*co:0,sn**2:s.Rational(1,2)})
expected=-lam**2*m*k*e**2/(2*r)
assert s.simplify(cov-expected)==0
torque=-s.diff(r*r*expected,r)/r
assert s.simplify(torque-lam**2*m*k*(e**2+2*r*e*s.diff(e,r))/(2*r))==0
# Substitute the leading phase k=-20 after deriving the general identity.
leading=s.simplify(torque.subs({m:4,k:-20})/lam**2)
assert s.simplify(leading+40*(e**2+2*r*e*s.diff(e,r))/r)==0

# Exact support, kernel cone, and horizontal-gradient bound.
rmin,rmax,zmin=F(1,25),F(6,25),F(17,5)
assert rmax/zmin==F(6,85)<F(1,4)
assert zmin>F(5,2)
assert rmax*rmax+F(23,5)**2<25
ximax=(rmax/F(63,200))**2
assert 8*ximax<5
# C in [0,1], C' in [-4,0] implies -8xi <= C+2xi C' <= 1.

R0=rmax*80*2+400*(rmax*rmax-rmin*rmin)/2+16*F(9,5)
Z0,Z1=F(12,5),F(160,3)
assert R0==F(392,5)
Cw=3*F(5,17)**5*R0*Z0
assert Cw==F(1764000,1419857)<F(3,2)
H=F(204,35)
assert H-F(3,2)/16>900**2*F(279,40000000)
assert H<1020**2*F(113,20000000)

Gt=F(1950,64)-F(19,1000)*1020
energy=80*900*F(1,50)*F(9,10)-F(7,5)*R0*Z0-\
       F(1,1000)*(200000*Z0+R0*Z1)
ratio=energy/(R0*Z0/2)
assert Gt==F(8871,800)>11
assert energy==F(205648,375)>548
assert ratio==F(12853,2205)>5

# Whole-pressure coefficient constants and the common family budget.
# The 4/m projection multiplier is 1, and the horizontal gradient bound is 5.
assert F(4,4)*5==5
local=135*80*(rmax**3-rmin**3)/3*Z0/zmin**7
assert local<F(23,1000)
visc=F(6,1000)*F(5,17)**5*(200000*Z0+R0*Z1)
assert visc==F(9078400,1419857)<7
margin=F(290649,8750)-(F(3999,1000)*F(3,2)+7+F(1020,3))/16
assert margin==F(12493177,1120000)>11
assert F(pressure["common_T_negative_margin"])==margin
assert F(pressure["beta_second_lower"])==margin/2
assert F(envelope["mean_derivative_lower"])==Gt
assert F(envelope["fluctuation_energy_fractional_derivative_lower"])==ratio
print("PASS: non-flat leading covariance/torque, exact retuning and energy "
      "margins, full-pressure E<1/3, and common actual initial T<-11.")

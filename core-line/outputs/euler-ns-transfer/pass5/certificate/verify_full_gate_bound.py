#!/usr/bin/env python3
"""Exact aggregation of saved interval endpoints and analytic coarse bounds."""
import json
from fractions import Fraction as F
from pathlib import Path

here = Path(__file__).resolve().parent
def dyadic(value):
    return (-1 if value["sign"] else 1)*F(int(value["mantissa"]))*F(2)**value["exponent"]

radial = json.loads((here/"radial-strain-interval-1024.json").read_text())
mixed = json.loads((here/"core-twb-interval.json").read_text())
local = json.loads((here/"local-interval-results.json").read_text())
C300_upper = dyadic(radial["exact_dyadic_bounds"]["upper"])
tWb_upper = dyadic(mixed["tWb"]["exact_dyadic_bounds"]["upper"])
Cv_lower = dyadic(local["Cv"]["exact_dyadic_bounds"]["lower"])
dv_upper = dyadic(local["dv"]["exact_dyadic_bounds"]["upper"])
assert C300_upper < 66
assert tWb_upper < 4
assert Cv_lower > 0
assert dv_upper < 901*Cv_lower

kcore_upper = -F(87238,35287)
kpump_upper = -F(358640696595957938134855826987,23085491673405775693824000000)
assert kcore_upper < -F(12,5)
assert kpump_upper < -F(31,2)

b = omega = F(1)
c, nu = F(7,5), F(1,1000)
H = F(38,7)*b*b+F(2,5)*omega*omega
assert H == F(204,35)
initial_pzz = -F(18,7)*b*b+F(2,5)*omega*omega-H
assert initial_pzz == -8*b*b
assert -2*b*b-initial_pzz/2 == 2*b*b

upper = 32+66+4+H*(-F(12,5)-c*F(31,2)+nu*901)
assert upper == -F(290649,8750)
assert upper < 0
print("PASS exact saved dyadic endpoints: C300<66, tWb<4, Cv>0, dv/Cv<901.")
print("PASS exact support bounds: kcore<-12/5, kpump<-31/2.")
print("PASS exact neutral tuning H=204/35, pzz=-8 and b'=2.")
print("Full gate strict upper bound:", upper, "=", float(upper))
print("beta'' strict lower bound:", -upper/2, "=", float(-upper/2))

lam = F(1000)
assert lam*nu == 1
assert (lam*b)**3 == lam**3*b**3
assert (lam*nu)*(lam**2) == lam**3*nu
assert -lam**3*upper/(2*lam*omega) == lam**2*(-upper/(2*omega))
print("PASS rescaling to viscosity one: gate times 10^9, beta'' times 10^6.")
print("Aggregation only: analytic identities and interval range enclosures remain explicit dependencies.")

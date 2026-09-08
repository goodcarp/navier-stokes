#!/usr/bin/env python3
"""
u6 -- the l = 1 strain rate kappa of the mollified datum, isolated feature by feature.

kappa(shell) := -(3/5) g_1 = -(3/4) INT_0^pi omega^theta(phi) cos phi sin^2 phi dphi / M
(far-near-kernel-lemma, Consequence A: this is the l=1 interior mode, exactly constant
inside the shell, and a(0,s) = INT kappa dlog rho').  Bang-bang gives exactly 1/2.

This was produced as a by-product of the u4 control  a(0) = kappa * INT Theta dlog rho
(which reproduces the quadrature to 1.1e-8), and it prices something the assembly's
sec 1.3 does not: the EQUATORIAL mollification, not the axis taper, is what moves kappa.
"""
import json, math
from scipy import integrate

SD = math.sin(math.radians(7.5))
OUT = {}

def kappa(T, Z):
    f = lambda p: 0.75*T(p)*Z(p)*math.cos(p)*math.sin(p)**2
    v, _ = integrate.quad(f, 0, math.pi, limit=400, points=[math.pi/2],
                          epsabs=1e-13, epsrel=1e-13)
    return v

one = lambda p: 1.0
sgn = lambda p: (1.0 if math.cos(p) > 0 else (-1.0 if math.cos(p) < 0 else 0.0))
tanh_tap = lambda p: math.tanh(math.sin(p)/SD)
sharp_tap = lambda p: min(1.0, min(p, math.pi-p)/math.radians(7.5))
def tanh_eq(w):
    return lambda p: math.tanh(math.cos(p)/w)
def lin_eq(dm):
    return lambda p: math.copysign(min(1.0, abs(p - math.pi/2)/math.radians(dm)),
                                   math.cos(p))

OUT['bang_bang'] = kappa(one, sgn)
OUT['sharp_taper_7.5_only'] = kappa(sharp_tap, sgn)
OUT['tanh_taper_7.5_only'] = kappa(tanh_tap, sgn)
OUT['eq_tanh_w0.20_only'] = kappa(one, tanh_eq(0.20))
OUT['eq_linear_dm5_only'] = kappa(one, lin_eq(5.0))
OUT['eq_linear_dm11.5_only'] = kappa(one, lin_eq(11.5))
OUT['campaign_D-C (tanh taper 7.5 + tanh eq w=0.20)'] = kappa(tanh_tap, tanh_eq(0.20))
OUT['ASSEMBLY_1.1 (sharp taper 7.5 + linear dm=5)'] = kappa(sharp_tap, lin_eq(5.0))
for k in list(OUT):
    OUT['deficit_1_minus_2kappa__' + k] = 1 - 2*OUT[k]
# the clock floor eps -> (1-2kappa)/(2kappa) of ASSEMBLY (1.1)
for k in ('sharp_taper_7.5_only', 'campaign_D-C (tanh taper 7.5 + tanh eq w=0.20)',
          'ASSEMBLY_1.1 (sharp taper 7.5 + linear dm=5)'):
    OUT['clock_floor_eps__' + k] = (1 - 2*OUT[k])/(2*OUT[k])

with open('u6_results.json', 'w') as f:
    json.dump(OUT, f, indent=1)
for k, v in OUT.items():
    print(f'{k:56s} {v!r}')

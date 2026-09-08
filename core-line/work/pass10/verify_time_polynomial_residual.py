#!/usr/bin/env python3
"""Exact time-polynomial coefficients and H4 validation obstruction.

Does not calculate the full spatial coefficients or validate an evolved
NS field. The initial mean-jet premise is read from its certified endpoint.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import os
import sympy as s

here = Path(__file__).resolve().parent
t = s.symbols('t', real=True)
L0,L1,L2 = s.symbols('L0 L1 L2')
b = {(i,j):s.symbols(f'B{i}{j}') for i in range(3) for j in range(3)}
a1 = L0-b[0,0]
a2 = L1-b[0,1]-b[1,0]
a3 = L2-b[0,2]-b[2,0]-2*b[1,1]

# B is bilinear, not assumed symmetric. L is the linear viscous operator.
c1 = [1,t]
R1 = a1-(L0+t*L1)+sum(c1[i]*c1[j]*b[i,j]
                               for i in range(2) for j in range(2))
assert s.expand(R1-(-t*a2+t*t*b[1,1])) == 0
c2 = [1,t,t*t/2]
R2 = a1+t*a2-(L0+t*L1+t*t*L2/2)+sum(c2[i]*c2[j]*b[i,j]
                                           for i in range(3) for j in range(3))
expected = -t*t*a3/2+t**3*(b[1,2]+b[2,1])/2+t**4*b[2,2]/4
assert s.expand(R2-expected) == 0

# A rational dual test annihilates every quadratic residual coefficient.
z=F(4,5); phi=F(61,64)
assert 0<phi<1
assert phi*z**3/3-(1-z**3)/3 == 0
assert phi*z**2/2-(1-z**2)/2 == F(1,8)

def endpoint(q):
    return (-1 if q['sign'] else 1)*F(int(q['mantissa']))*F(2)**q['exponent']

certificate_candidates=[here.parent/'pass9'/'general-envelope-initial-jet-certificate.json',
    here.parent.parent/'pass9'/'checks'/'general-envelope-initial-jet-certificate.json']
if os.environ.get('PASS10_INITIAL_JET_CERTIFICATE'):
    certificate_path=Path(os.environ['PASS10_INITIAL_JET_CERTIFICATE'])
else:
    certificate_path=next(p for p in certificate_candidates if p.is_file())
certificate=json.loads(certificate_path.read_text())
assert endpoint(certificate['Gtt_actual_full_pressure']['exact']['upper']) < -13000
assert F(certificate['root_interval'][1]) < F(1,4)

# Same explicit A=1020 and lambda=1/4, with the existing pressure bounds.
d_upper=(1020**2*F(279,40000000)+F(3,2)/16-F(204,35))/2
assert 0<d_upper<F(4,5)
critical_time=F(1,16250)
assert 6500*critical_time == F(2,5)
T=F(1,1000)
assert 6500*T*T == F(13,2000)
assert F(2,5)*T == F(1,2500)
assert F(13,2000)>F(1,2500)

# Sobolev product constants used for the practical residual upper bounds.
assert (3*2**4)**2*comb(7,3) < 288**2
assert (3*2**5)**2*comb(8,3) < 720**2

# Exterior radial integral exponent for a derivative of order j.
j,R=s.symbols('j R', positive=True)
tail=R**(-2*j-5)/(2*j+5)
assert s.simplify(s.diff(tail,R)+R**(-2*j-6)) == 0

print('PASS: exact first/quadratic Taylor residuals; rational cancellation-'
      'resistant residual lower bound; certified initial-jet premise; '
      'linear-candidate H4 cone obstruction at T>=1/16250; product constants '
      'and exterior integral. No actual-flow error lower bound or duration '
      'certificate is claimed.')

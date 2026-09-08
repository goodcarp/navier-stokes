#!/usr/bin/env python3
"""Check the quantitative common-family consequence of the actual certificates."""
import json
from pathlib import Path
from fractions import Fraction as F
import sympy as s
here=Path(__file__).resolve().parent
def endpoint(x):return F((-1 if x['sign'] else 1)*int(x['mantissa']))*F(2)**x['exponent']
D=json.loads((here/'D-viscosity-certificate.json').read_text())
E=json.loads((here/'outer-E-certificate.json').read_text())
assert D['status']==E['status']=='PASS'
assert D['viscous_upper_strictly_below_four']
assert D['radial_upper_strictly_below_120000']
assert E['E_upper_strictly_below_one_sixtieth']
assert endpoint(D['viscous_upper_bound_exact']['upper'])<4
assert endpoint(E['E_upper_expression']['exact']['upper'])<F(1,60)
remaining=F(23199,1000)-F(96,5)
assert remaining==F(3999,1000)>0
bracket=remaining*F(3,2)+4+F(1020,60)
margin=F(290649,8750)-bracket
assert margin==F(435297,70000)>0
beta2=margin/2;g0=F(626,35)
assert beta2==F(435297,140000)
assert beta2/4==F(435297,560000)
assert g0/2==F(313,35)
# Initial discrete-symmetry/flat-core scalar ODE consequence.
t,P=s.symbols('t P',real=True)
# b'=Omega'=2, Omega''=8, b''=-8-P/2 when p_zz'=P.
b=1+2*t+(-8-P/2)*t*t/2
om=1+2*t+8*t*t/2
assert s.simplify(s.diff(b/om,t,2).subs(t,0)+(P+32)/2)==0
assert s.simplify(s.diff(s.diff(b,t)-2*b*b,t).subs(t,0)+(P+32)/2)==0
print('PASS: actual certificate endpoints imply T_lambda < -435297/70000 and beta_second >435297/140000 for every lambda in [3/4,1]; local-time constants agree. Persistence uses smooth local existence, not this arithmetic checker.')

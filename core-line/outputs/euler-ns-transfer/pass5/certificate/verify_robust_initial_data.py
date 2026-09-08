#!/usr/bin/env python3
"""Exact initial-jet algebra and rational-amplitude local feedback bounds."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as s

e,T=s.symbols('e T',real=True)
bp=2+e/2
op=s.Integer(2)
bpp=-4*bp-(T-32)/2
opp=2*bp+2*op
betapp=bpp-opp-2*op*(bp-op)
assert s.simplify(betapp+T/2+5*e)==0
assert s.simplify(bpp-4*bp+T/2+4*e)==0
assert s.simplify(opp-(8+e))==0

def dyadic(x):
    return F((-1 if x['sign'] else 1)*int(x['mantissa']))*F(2)**x['exponent']

record=json.loads((Path(__file__).resolve().parent/'local-interval-results.json').read_text())['Cv']['exact_dyadic_bounds']
cvlo,cvhi=dyadic(record['lower']),dyadic(record['upper'])
assert F(113,20000000)<cvlo and cvhi<F(279,40000000)
H0=F(204,35)
lo=1100**2*F(113,20000000)
hi=1100**2*F(279,40000000)
assert lo==F(13673,2000) and hi==F(33759,4000)
assert lo>H0
delta=F(290649,17500)
assert -(102-F(23199,1000)*H0)/2==delta
assert F(23199,2000)-5==F(13199,2000)
assert F(23199,2000)-4==F(15199,2000)
print('PASS exact nonneutral central second derivatives.')
print('PASS rational A=1100: H lies strictly between',lo,'and',hi)
print('beta prime and S strict lower bound:',(lo-H0)/2)
print('beta second derivative strict lower bound:',delta+F(13199,2000)*(lo-H0))
print('S prime strict lower bound:',delta+F(15199,2000)*(lo-H0))
print('Initial local family only; no quantitative perturbation radius or return.')

#!/usr/bin/env python3
"""Independent rational endpoint/constant replay of the saved E certificate."""
import json
from fractions import Fraction as F
from pathlib import Path
import sympy as s

here = Path(__file__).resolve().parent
data = json.loads((here/"outer-E-certificate.json").read_text())

def dyadic(x):
    return (-1 if x["sign"] else 1)*F(int(x["mantissa"]))*F(2)**x["exponent"]

def upper(key):
    return dyadic(data[key]["exact"]["upper"])

assert data["status"] == "PASS"
assert upper("E_upper_expression") < F(1,60)
assert upper("compatibility_upper_expression") < F(5093339,210000)

# Exact cone support and elementary kernel inequalities used by the code.
tmax = F(7,68)**2
assert 0 < tmax < F(4,3)
assert tmax < 6
t = s.symbols("t", nonnegative=True)
f = (6-t)/(1+t)**5
assert s.factor(s.diff(f,t)) == (4*t-31)/(t+1)**6
assert 4*tmax-31 < 0
# On the cone, 0<=6-t<=6, 0<=4-3t<=4, and (1+t)^(-9/2)<=1.
# Thus Qdiff and Qrz constants are 15*6/(4*pi),15*4/(4*pi).
assert F(15*6,4) == F(45,2)
assert F(15*4,4) == 15
# Angular integration converts the two trial-pair constants.
assert 2*F(45,2) == 45
assert 2*15 == 30
# Exact local coefficient 6*pi*15/(4*pi).
assert F(6*15,4) == F(45,2)

# The all-lambda consequence retains the favorable Cw contribution in D.
remaining_Cw = F(23199,1000)-F(96,5)
assert remaining_Cw == F(3999,1000)
bracket = remaining_Cw*F(3,2)+4+F(1020,60)
assert bracket == F(53997,2000) > 0
margin = F(290649,8750)-bracket
assert margin == F(435297,70000) > 6

print("PASS: dyadic E<1/60, separated kernel constants, and all-lambda "
      "initial gate T_lambda < -435297/70000, conditional on the independently "
      "proved D<4-(96/5)Cw and prior family bounds.")

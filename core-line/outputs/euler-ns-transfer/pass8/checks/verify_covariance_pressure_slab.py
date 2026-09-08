#!/usr/bin/env python3
"""Exact endpoint and algebra replay of the harmonic-slab pressure error."""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as s
from mpmath import iv

here=Path(__file__).resolve().parent
cert=next(p for p in [here.parent/"pass7"/"outer-E-certificate.json",here.parents[1]/"pass7/checks/outer-E-certificate.json"] if p.exists())
data=json.loads(cert.read_text())
ep=data["residual_norm_upper_expression"]["exact"]["upper"]
ne=(-1 if ep["sign"] else 1)*F(int(ep["mantissa"]))*F(2)**ep["exponent"]
assert ne<F(85,4)
assert iv.pi.a > iv.mpf(25)/8

# General positive-weight derivative identity, checked explicitly at m=4.
t=s.symbols("t",real=True)
m=4
w=(1-t*t)**s.Rational(7,2)
weight=(m-1)*w-t*s.diff(w,t)
expected=(m-1)*(1-t*t)**s.Rational(7,2)+(2*m-1)*t*t*(1-t*t)**s.Rational(5,2)
assert s.simplify(weight-expected)==0
# By integration by parts, the integral is m times the original weight.
assert s.simplify(weight-m*w+s.diff(t*w,t))==0

# Exact real even/odd harmonic-slab energy integrand and antiderivative.
k,h=s.symbols("k h",positive=True,real=True)
x=s.symbols("x",real=True)
aa,bb=s.symbols("aa bb",real=True)
f=aa*s.cosh(k*x)+bb*s.sinh(k*x)
energy=s.expand(s.diff(f,x)**2+k*k*f*f)
even=s.simplify((energy+energy.subs(x,-x))/2)
expected_even=k*k*(aa*aa+bb*bb)*s.cosh(2*k*x)
assert s.simplify(even-expected_even)==0
assert s.simplify(s.integrate(expected_even,(x,-h,h))
                 -k*s.sinh(2*k*h)*(aa*aa+bb*bb))==0
# Apply separately to real/imaginary parts for complex Hankel amplitudes.

assert F(2*s.factorial(8),2**8*s.factorial(3)**2)==F(35,4)
odd_tail=F(19,16*3**9)
assert 1+odd_tail<F(10001,10000)
square=(F(85,4)**2/F(16))*F(8,25)*F(35,4)*\
       (F(63,200)**2*F(5,8))**3/F(9,20)**9*F(10001,10000)
assert square==F(2380277273927,95551488000)<25

# Exact Green derivative avoids differentiating interval moment enclosures.
r,mm=s.symbols("r mm",positive=True)
C=s.Function("C")(r);b=s.Function("b")(r)
Ip=s.Function("Ip")(r);Im=s.Function("Im")(r)
P=-2*C*b-(mm-1)*r**(-mm)*Ip-(mm+1)*r**mm*Im
dp=s.diff(P,r).subs({s.diff(Ip,r):r**mm*s.diff(C,r)*b,
                    s.diff(Im,r):-r**(-mm)*s.diff(C,r)*b})
target=-2*C*s.diff(b,r)+mm*(mm-1)*r**(-mm-1)*Ip-mm*(mm+1)*r**(mm-1)*Im
assert s.simplify(dp-target)==0
print("PASS: saved residual ne<85/4, positive Bessel derivative weight, "
      "harmonic-slab energy normalization, whole-space |pressure_r-P'|<5, "
      "and exact derivative-reduced Green gradient.")

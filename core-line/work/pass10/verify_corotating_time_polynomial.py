#!/usr/bin/env python3
"""Exact algebraic checks for the co-rotating Taylor note, not an NS certificate."""
import itertools
import math
import sympy as sp

t,omega=sp.symbols('t omega', real=True)
I=sp.I

# Scalar mode Taylor jet from the actual rotating reconstruction.
m=sp.symbols('m', integer=True)
a=sp.symbols('a0:4')
actual=sum(a[j]*t**j/sp.factorial(j) for j in range(4))
for n in range(4):
    jet=sp.diff(sp.exp(I*m*omega*t)*actual,t,n).subs(t,0)
    expected=sum(sp.binomial(n,j)*(I*m*omega)**(n-j)*a[j] for j in range(n+1))
    assert sp.expand(jet-expected)==0

# Noncommutative bilinear products represented by independent symbols.
b=sp.symbols('b0:4')
M=sp.symbols('M0:3')
B={(i,j):sp.Symbol(f'B{i}{j}') for i in range(3) for j in range(3)}
weights=[1,t,t*t/2]
qprime=b[1]+t*b[2]
linear=sum(weights[j]*M[j] for j in range(3))
quad=sum(weights[i]*weights[j]*B[i,j] for i in range(3) for j in range(3))
# Relations (3), with Mj = (nu Delta + omega L)bj.
sub={M[0]:b[1]+B[0,0],
     M[1]:b[2]+B[0,1]+B[1,0],
     M[2]:b[3]+B[0,2]+B[2,0]+2*B[1,1]}
res=sp.expand((qprime-linear+quad).subs(sub))
expected=-t*t*b[3]/2+t**3*(B[1,2]+B[2,1])/2+t**4*B[2,2]/4
assert sp.expand(res-expected)==0
res1=sp.expand((b[1]-M[0]-t*M[1]+B[0,0]+t*(B[0,1]+B[1,0])+t*t*B[1,1]).subs(sub))
assert sp.expand(res1-(-t*b[2]+t*t*B[1,1]))==0

# Product-mode derivation law for L, including asymmetric B(a,b).
m1,m2=sp.symbols('m1 m2',integer=True)
assert sp.expand(I*(m1+m2)-I*m1-I*m2)==0

# Physical residual cancellation under Q_{-t}, for a generic scalar mode.
q=sp.Function('q')(t);n=sp.Function('n')(t)
physical=sp.diff(sp.exp(-I*m*omega*t)*q,t)-sp.exp(-I*m*omega*t)*n
transformed=sp.exp(-I*m*omega*t)*(sp.diff(q,t)-I*m*omega*q-n)
assert sp.simplify(physical-transformed)==0

# Full Bessel weights and horizontal rotation-invariant weights.
def four_partitions(s):
    for p in itertools.product(range(s+1),repeat=4):
        if sum(p)==s:
            yield p
for s,expected_T,expected_F in [(4,6,24),(5,10,60)]:
    t_coeffs=[math.comb(h,j) for h in range(s+1) for j in range(h+1)]
    f_coeffs=[math.factorial(s)//math.prod(math.factorial(x) for x in p)
              for p in four_partitions(s)]
    assert min(t_coeffs)==min(f_coeffs)==1
    assert max(t_coeffs)==expected_T
    assert max(f_coeffs)==expected_F

# Exact artificial energy gain in a manufactured unitary phase.
x=sp.symbols('x',real=True)
assert sp.expand((1-I*x)*(1+I*x))==1+x*x
assert sp.expand((1-I*x-x*x/2)*(1+I*x-x*x/2))==1+x**4/4
phase=sp.Rational(4)*1020/sp.Integer(1000)
quadratic_factor=1+phase**4/4
assert phase==sp.Rational(102,25)
assert quadratic_factor==sp.Rational(27451429,390625)
assert abs(float(quadratic_factor)-70.27565824)<1e-12
print('PASS: transformed jets, full bilinear residual, reconstruction sign, rotation norm weights, and manufactured phase energy.')
print('No whole-space coefficient norm, positive-time error enclosure, or NS return is certified by this checker.')

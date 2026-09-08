#!/usr/bin/env python3
"""Exact amplitude-polynomial reduction and sufficient pressure threshold."""
from fractions import Fraction as F
import sympy as s

A,lam,nu=s.symbols('A lam nu',real=True)
# Rotating by pi/4 negates the pure m=4 seed. Even projection removes all
# odd seed monomials from the full degree-three/degree-two functional.
terms={(i,j):s.Symbol('c%d%d'%(i,j)) for i in range(4) for j in range(4-i)}
poly=sum(c*A**i*lam**j for (i,j),c in terms.items())
even=s.expand((poly+poly.subs(lam,-lam))/2)
coefficient=s.expand(even.coeff(lam,2))
assert s.Poly(coefficient,A).degree()==1
assert s.simplify(even.coeff(lam,1))==0 and s.simplify(even.coeff(lam,3))==0

H=F(204,35); k=F(23199,1000); m0=F(290649,8750)
assert 102-k*H==-m0
Cw,D,E=s.symbols('Cw D E',real=True)
bound=102-s.Rational(k.numerator,k.denominator)*(s.Rational(H.numerator,H.denominator)-lam**2*Cw)+lam**2*(D+A*E)
assert s.simplify(bound+s.Rational(m0.numerator,m0.denominator)-lam**2*(s.Rational(k.numerator,k.denominator)*Cw+D+A*E))==0
threshold=F(16,9)*m0-k*F(3,2)
assert threshold>24
assert -m0+F(9,16)*(k*F(3,2)+threshold)==0
assert F(3,4)**2*F(10000,63)-F(70,1000)*1020==F(626,35)>0
print('PASS: exact rotation/degree selection leaves lambda²(D_nu+A E).')
print('PASS: existing radial bound and retuning yield the correct sufficient inequality.')
print('Threshold for D_nu+A_(3/4)E:',threshold,'=',float(threshold))
print('Same amplitude already passes the mean-maximum test; interaction bound remains unevaluated.')

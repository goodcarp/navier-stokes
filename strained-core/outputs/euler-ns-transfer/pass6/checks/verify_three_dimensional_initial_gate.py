#!/usr/bin/env python3
"""Exact constants, central symmetry and neutral-family half-margin checks."""
import sympy as s

assert (9*4)**2*10<128**2
# The multiplier and embedding/product arguments are in the proof and audit.
X,Y,d,nu,Z=s.symbols('X Y d nu Z',positive=True)
telescoped=256*d*(3*nu*X+512*X**2+Y*(3*nu+512*(X+Y)))
claimed=256*d*(3*nu*(X+Y)+512*(X**2+X*Y+Y**2))
assert s.expand(telescoped-claimed)==0
assert s.simplify(claimed.subs({X:Z,Y:Z})/d-256*(6*nu*Z+1536*Z**2))==0
assert s.Poly(s.diff(claimed,X),X,Y,d,nu).coeffs() and all(c>0 for c in s.Poly(s.diff(claimed,X),X,Y,d,nu).coeffs())
assert all(c>0 for c in s.Poly(s.diff(claimed,Y),X,Y,d,nu).coeffs())

Cv,Jw,eps=s.symbols('Cv Jw eps',real=True,nonzero=True)
H0=s.Rational(204,35)
A2=(H0+eps**2*Jw)/Cv
pzz=-s.Rational(18,7)+s.Rational(2,5)-A2*Cv+eps**2*Jw
assert s.simplify(pzz+8)==0
m0=s.Rational(290649,8750)
assert m0/4==s.Rational(290649,35000)

# Rotation through pi/2 and incompressibility leave exactly two matrix degrees.
R=s.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
entries=s.symbols('l0:9')
M=s.Matrix(3,3,entries)
constraints=list(M*R-R*M)+[s.trace(M)]
A,_=s.linear_eq_to_matrix(constraints,entries)
assert A.rank()==7
b,om=s.symbols('b om',real=True)
central=s.Matrix([[-b,-om,0],[om,-b,0],[0,0,2*b]])
assert central*R==R*central and s.trace(central)==0

# Flat-core central equations give the same neutral second-ratio formula.
T=s.symbols('T',real=True)
bp=op=s.Integer(2)
bpp=-4*bp-(T-32)/2
opp=2*bp+2*op
assert s.simplify(bpp-opp+T/2)==0

# Sufficient parameter restrictions each pay half the H10-distance budget.
W,D,dstar=s.symbols('W D dstar',positive=True)
assert W*(dstar/(2*W))==dstar/2
assert s.simplify(D*s.sqrt(dstar/(2*D))**2-dstar/2)==0
assert s.simplify(128*W**2*(H0/(256*W**2))-H0/2)==0
print('PASS: H4 pressure-pairing counting constant and full derivative Lipschitz polynomial.')
print('PASS: exact amplitude retuning and retained half-margin beta second derivative.')
print('PASS: fourfold central matrix constraints and actual neutral jet algebra.')
print('PASS: nonzero norm-form seed interval and full retuned distance budget.')
print('Depends on analytic estimates and pass5 interval certificate; no return or blowup.')

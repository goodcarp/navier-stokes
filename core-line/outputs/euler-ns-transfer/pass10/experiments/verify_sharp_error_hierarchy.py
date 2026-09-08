#!/usr/bin/env python3
"""Exact checks for sharp-error-hierarchy.md; no evolution is certified."""
import itertools
import json
import math
from fractions import Fraction as F
import numpy as np
import sympy as s

checks=[]
def record(condition,label):
    assert bool(condition),label
    checks.append(label)

# The selected unchanged outer swirl supplies the lower bounds by axial
# Poincare, using only the conservative rational inequality pi>3.
A=1020;r0=F(63,400);d=F(9,10)
lower4=A*A*r0**4*F(3)**9/(2*d**7)
lower5=A*A*r0**4*F(3)**11/(2*d**9)
record(lower4==F(843075135,64) and lower4>3600**2,
       'exact H4 lower bound')
record(lower5==F(2341875375,16) and lower5>12000**2,
       'exact H5 lower bound')
record(F(560)*12000*F(1,1000)==6720,
       'uniform baseline exponent with exact initialization')
record(F(560)*(3600-1)*F(1,1000)>2015,
       'uniform baseline exponent with H4 initial error at most one')

ell,S=s.symbols('ell S',positive=True)
M={j:s.symbols(f'M{j}',nonnegative=True) for j in range(2,6)}
B=s.zeros(5)
for j in range(5):
    B[j,j]=(j+1)*S
    for k in range(2,j+1):
        B[j,j-k+1]+=s.binomial(j,k)*ell**(k-1)*M[k]
    for k in range(1,j+1):
        B[j,j-k]+=s.binomial(j,k)*ell**k*M[k+1]
for j in range(5):
    for h in range(j):
        record(s.expand(B[j,h]-s.binomial(j+1,h)*ell**(j-h)*M[j-h+1])==0,
               f'combined commutator coefficient row {j} column {h}')
    record(B[j,j]==(j+1)*S and all(B[j,h]==0 for h in range(j+1,5)),
           f'linear hierarchy triangularity row {j}')
q=[2**j-1 for j in range(5)]
record(q==[0,1,3,7,15] and sum(a*a for a in q)==284<17**2,
       'nonlinear component weights and scalar constant')
for j in range(1,5):
    for k in range(1,j+1):
        a,b=sorted((k,j-k+1))
        record(a<=2 and b<=4 and F(j)-a-F(3,2)-b==F(-5,2),
               f'nonlinear Sobolev split and scale j={j} k={k}')

def multi(total):
    for a in range(total+1):
        for b in range(total-a+1):
            yield (a,b,total-a-b)
def multiplicity(alpha):
    return math.factorial(sum(alpha))//math.prod(math.factorial(a) for a in alpha)
alphas=[a for j in range(5) for a in multi(j)]
record(len(alphas)==35 and max(map(multiplicity,alphas))==12
       and min(map(multiplicity,alphas))==1,
       'ordered versus derivative-sum norm multiplicities')

# Summing scalar H2 embeddings of an a-th ordered derivative tensor costs
# no more than Q_a^2+Q_(a+1)^2+Q_(a+2)^2. Check each coefficient independently.
for a in range(3):
    for total in range(a,a+3):
        for gamma in multi(total):
            coefficient=sum(multiplicity(alpha) for alpha in multi(a)
                if all(ai<=gi for ai,gi in zip(alpha,gamma)))
            record(coefficient<=multiplicity(gamma),
                   f'tensor H2 embedding multiplicity a={a} gamma={gamma}')
for length in (F(1,100),F(1,10),F(1),F(2)):
    floor=min(F(1),length**8)
    record(all(length**(2*sum(a))*multiplicity(a)>=floor for a in alphas),
           f'H4 endpoint conversion length={length}')

# Exact integer matrices check that rigid rotation is skew on each input
# derivative index and on the velocity-component index. No numerical
# eigenvalue computation enters this test.
J=np.array([[0,-1,0],[1,0,0],[0,0,0]],dtype=np.int64)
D=np.diag([-1,-1,2]).astype(np.int64)
Id=np.eye(3,dtype=np.int64)
for order in range(5):
    slots=order+1
    n=3**slots
    rotation=np.zeros((n,n),dtype=np.int64)
    strain=np.zeros((n,n),dtype=np.int64)
    for slot in range(slots):
        R=np.ones((1,1),dtype=np.int64)
        H=np.ones((1,1),dtype=np.int64)
        for pos in range(slots):
            R=np.kron(R,J if pos==slot else Id)
            H=np.kron(H,D if pos==slot else Id)
        rotation+=R;strain+=H
    record(np.array_equal(rotation+rotation.T,np.zeros_like(rotation)),
           f'rigid rotation skew on tensor order {order}')
    record(np.max(np.abs(np.diag(strain)))==2*slots
           and np.count_nonzero(strain-np.diag(np.diag(strain)))==0,
           f'tensor strain bound on order {order}')

print(json.dumps({'status':'PASS','checks':len(checks),
    'exact_H4_lower_squared':str(lower4),'exact_H5_lower_squared':str(lower5),
    'uniform_exact_initial_exponent_lower':6720,
    'nonlinear_scalar_constant_squared':284,
    'scope':'Exact rational/combinatorial/tensor checks for the proved '
            'error inequality. No actual evolution norm, residual, '
            'comparison-ODE solution or endpoint gain has been enclosed.'},indent=2))

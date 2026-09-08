#!/usr/bin/env python3
"""Manufactured comparison-system checks; no synthetic input is NS data."""
from fractions import Fraction as F
from math import factorial
import json
from validated_error_hierarchy import (
    Options,Bounds,Q,rational,dyadic_up,sqrt_up,coefficient_matrix,
    exponential_tail_order,linear_flow_enclosure,close_radius,certify_slab,
    propagate,reweight_endpoint,ClosureFailure)

checks=[]
def require(condition,name):
    assert condition,name
    checks.append(name)

for x in (F(0),F(1,3),F(5,7),F(1,10**100),F(10**100),F(17,4)):
    up=dyadic_up(x,48)
    require(up>=x and (x==0 or up-x<=x/F(2**46)),f'outward dyadic rounding {x}')
    rt=sqrt_up(x,48)
    require(rt*rt>=x and (x==0 or rt*rt-x<=x/F(2**44)),f'outward rational square root {x}')
try:
    rational(.001)
except TypeError:
    require(True,'ambiguous binary-float API input is rejected')
else:
    raise AssertionError('binary float accepted')
for label,call in (
        ('negative residual rejected by direct Bounds constructor',
         lambda:Bounds(1,0,(0,0,0,0),(-1,0,0,0,0))),
        ('nonpositive duration rejected before zero-error shortcut',
         lambda:Bounds(0,0,(0,0,0,0),(0,0,0,0,0))),
        ('invalid direct Options constructor rejected',
         lambda:Options(precision_bits=0)),
        ('mismatched radius vectors rejected',lambda:close_radius([1,2],[1]))):
    try:
        call()
    except ValueError:
        require(True,label)
    else:
        raise AssertionError(label)

# A one-dimensional exponential enclosure is checked against tight rational
# reference brackets; its validity follows from the positive Taylor proof.
opt=Options(precision_bits=64,exponential_terms_cap=256)
flow=linear_flow_enclosure([[F(1)]],F(1),opt)
require(F(2718281,10**6)<flow['H_lower'][0][0]
        <=flow['H_upper'][0][0]<F(2718282,10**6),'exponential reference bracket')
require(flow['J_lower'][0][0]<=flow['H_upper'][0][0]-1
        and flow['J_upper'][0][0]>=flow['H_lower'][0][0]-1,
        'integrated exponential exact identity is enclosed')

# Independent matrix-polynomial calculation for zero diagonal. B^5=0, so
# exp(Bh) and its integral are finite exact rational expressions.
def multiply(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),F(0))
             for j in range(len(B[0]))] for i in range(len(A))]
B=coefficient_matrix(F(1),F(0),[1,2,3,4]);h=F(2,7)
nil=linear_flow_enclosure(B,h,opt)
power=[[F(int(i==j)) for j in range(5)] for i in range(5)]
H=[[F(0) for _ in range(5)] for _ in range(5)]
J=[[F(0) for _ in range(5)] for _ in range(5)]
for k in range(5):
    for i in range(5):
        for j in range(5):
            H[i][j]+=h**k*power[i][j]/factorial(k)
            J[i][j]+=h**(k+1)*power[i][j]/factorial(k+1)
    power=multiply(power,B)
require(all(v==0 for row in power for v in row),'manufactured B is nilpotent of order at most five')
require(nil['H_lower']==H==nil['H_upper'],'finite path exponential agrees with exact matrix polynomial')
require(nil['J_lower']==J==nil['J_upper'],'finite path integral agrees with exact matrix polynomial')
chain=linear_flow_enclosure(coefficient_matrix(F(1),F(0),[1,0,0,0]),h,opt)
require(chain['H_upper'][4][0]==F(15,2)*h**4,'manufactured high-order linear feed coefficient')

# Nonzero distinct diagonal rates and a nontrivial off-diagonal path have
# independent elementary exponential formulas. A much longer direct Taylor
# enclosure supplies the reference intervals.
def reference_exp(x,N=160):
    term=F(1);partial=term
    for n in range(1,N+1):
        term=term*x/n;partial+=term
    remainder=(term*x/(N+1))/(1-x/F(N+2))
    return partial,partial+remainder
e1=reference_exp(F(1));e2=reference_exp(F(2))
two=linear_flow_enclosure([[F(1),F(0)],[F(1),F(2)]],F(1),opt)
reference_H=(e2[0]-e1[1],e2[1]-e1[0])
reference_J=((e2[0]-1)/2-(e1[1]-1),(e2[1]-1)/2-(e1[0]-1))
require(two['H_lower'][1][0]<=reference_H[0]<=reference_H[1]<=two['H_upper'][1][0],
        'nonzero-rate path encloses independent exact exponential formula')
require(two['J_lower'][1][0]<=reference_J[0]<=reference_J[1]<=two['J_upper'][1][0],
        'nonzero-rate forced path encloses independent exact formula')

# Exact nonlinear trajectory: z=q*a, a'=284*a^2, ell=1, B=delta=0.
a0=F(1,1000)
initial=tuple(a0*q for q in Q)
problem={'ell':'1','initial_error':list(map(str,initial)),
         'slabs':[{'duration':'1/10','S':'0','M':['0']*4,'residual':['0']*5}],
         'options':{'precision_bits':64}}
result=propagate(problem)
exact_a=a0/(1-284*a0*F(1,10))
require(result['complete'] and F(result['validated_duration'])==F(1,10),
        'manufactured nonlinear interval is certified')
require(all(F(v)>=q*exact_a for v,q in zip(result['endpoint_error'],Q)),
        'certified endpoint encloses exact nonlinear solution')
for slab in result['accepted_slabs']:
    aa=list(map(F,slab['linear_initial_and_residual_upper']))
    bb=list(map(F,slab['nonlinear_response_upper']))
    K=F(slab['closure_K'])
    require(sum((a+b*K)**2 for a,b in zip(aa,bb))<=K,'stored exact radius inequality')

# The exact nonlinear solution remains finite at t=1, but the one-slab box
# fails. Automatic bisection closes it while inheriting every prior bound.
problem['slabs'][0]['duration']='1'
problem['options']['max_slab_bisections']=8
adaptive=propagate(problem)
require(adaptive['complete'] and len(adaptive['accepted_slabs'])>1
        and len(adaptive['failed_attempts'])>0,'adaptive subdivision handles a conservative failure')
exact_a=a0/(1-284*a0)
require(all(F(v)>=q*exact_a for v,q in zip(adaptive['endpoint_error'],Q)),
        'subdivided certificate encloses exact inherited nonlinear solution')
slabs=adaptive['accepted_slabs']
require(all(slabs[i]['initial']==slabs[i-1]['endpoint'] for i in range(1,len(slabs))),
        'no reset of inherited component errors')
require(sum(F(row['duration']) for row in slabs)==F(1),'accepted durations cover the complete request')

jump_problem={'ell':'1','initial_error':['0']*5,
    'slabs':[{'duration':'1/100','S':'0','M':['0']*4,'residual':['0']*5},
             {'duration':'1/100','S':'0','M':['0']*4,'residual':['0']*5,
              'jump_error':['1e-6','0','0','0','0']}]}
jumped=propagate(jump_problem)
require(jumped['complete'] and len(jumped['boundary_updates'])==1,
        'representation jump is included once at the original slab boundary')
require(F(jumped['accepted_slabs'][1]['initial'][0])>=F(1,10**6)
        and F(jumped['endpoint_error'][0])>=F(1,10**6),
        'boundary jump is inherited rather than reset')

problem['options']['max_slab_bisections']=0
failed=propagate(problem)
require(not failed['complete'] and F(failed['validated_duration'])==0
        and failed['endpoint_error']==list(map(str,initial)),
        'failed slab makes no unsupported progress')
require(failed['failed_attempts'][0]['reason']=='nonlinear_radius_does_not_close',
        'nonlinear closure failure is explicit, even for a finite exact solution')

# A giant off-diagonal coefficient is handled polynomially, with no
# exp(1e40) anywhere. These are deliberately synthetic comparison bounds.
giant={'ell':'1','initial_error':['1e-60','0','0','0','0'],
       'slabs':[{'duration':'1','S':'0','M':['0','0','0','1e40'],
                 'residual':['0']*5}],
       'options':{'precision_bits':64,'max_slab_bisections':0}}
large=propagate(giant)
require(large['complete'] and len(large['accepted_slabs'])==1,
        'large off-diagonal feed closes without a large exponential')
require(F(1,10**20)<=F(large['endpoint_error'][4])<F(11,10**21),
        'giant coefficient still yields a small manufactured high-order endpoint')
require(large['accepted_slabs'][0]['largest_diagonal_argument']=='0',
        'exponential control depends on the diagonal only')

# Non-square ell receives a rational upper enclosure of ell^(-5/2).
ell=F(2,3)
certificate=certify_slab(initial,Bounds(F(1,1000),F(0),(F(0),)*4,(F(0),)*5),ell,opt)
c=certificate['nonlinear_constant_upper']
require(c*c*ell**5>=1,'non-square verification length is safely enclosed')

# Exponential budget exhaustion is a distinct, honest failure mechanism.
resource={'ell':'1','initial_error':['1e-6']*5,
          'slabs':[{'duration':'1','S':'100','M':['0']*4,'residual':['0']*5}],
          'options':{'exponential_terms_cap':4,'max_slab_bisections':0}}
exhausted=propagate(resource)
require(not exhausted['complete'] and exhausted['failed_attempts'][0]['reason']
        =='exponential_series_budget_exhausted','exponential resource limit is reported as failure')
resource['initial_error']=['0']*5
zero=propagate(resource)
require(zero['complete'] and zero['endpoint_error']==['0']*5,
        'exact zero error and zero residual remain exactly zero')

new=reweight_endpoint([1,2,3,4,5],F(1,10),F(1,5))
require(new==tuple(F(v)*2**j for j,v in enumerate([1,2,3,4,5])),
        'length changes preserve derivative bounds exactly')
try:
    reweight_endpoint([0]*5,-1,-2)
except ValueError:
    require(True,'negative verification lengths are rejected')
else:
    raise AssertionError('negative lengths accepted')

print(json.dumps({'status':'PASS','checks':len(checks),
    'scope':'Manufactured exact/rational comparison-system tests only. '
            'No supplied test bounds describe or validate the selected NS evolution.',
    'nonlinear_manufactured_accepted_slabs':len(adaptive['accepted_slabs']),
    'giant_off_diagonal_manufactured_endpoint':large['endpoint_error'],
    'checked':checks},indent=2))

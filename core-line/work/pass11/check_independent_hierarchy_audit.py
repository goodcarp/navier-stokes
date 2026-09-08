"""Independent exact-rational tests of the comparison-system propagator.

All input systems are synthetic. No fluid approximation is certified here.
The reference matrix series uses a row-sum norm tail, independently of the
production increasing-path/homogeneous-polynomial kernel implementation.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
import random
from validated_error_hierarchy import (Options, Bounds, Q, ClosureFailure,
    linear_flow_enclosure, close_radius, certify_slab, propagate)


checks=[]
def require(condition,name):
    assert condition,name
    checks.append(name)


def multiply(A,B):
    n=len(A)
    return [[sum((A[i][k]*B[k][j] for k in range(n)),F(0))
             for j in range(n)] for i in range(n)]


def matrix_reference(B,h,N=80):
    """Exact positive Taylor polynomial plus independently proved norm tail."""
    n=len(B)
    power=[[F(i==j) for j in range(n)] for i in range(n)]
    H=[[F(0) for _ in range(n)] for _ in range(n)]
    J=[[F(0) for _ in range(n)] for _ in range(n)]
    for degree in range(N+1):
        for i in range(n):
            for j in range(n):
                H[i][j]+=h**degree*power[i][j]/factorial(degree)
                J[i][j]+=h**(degree+1)*power[i][j]/factorial(degree+1)
        power=multiply(power,B)
    x=h*max(map(sum,B))
    tail=(x**(N+1)/factorial(N+1))/(1-x/F(N+2))
    assert x<F(N+2)
    return H,J,tail,h*tail


opt=Options(precision_bits=48,exponential_terms_cap=160)
rng=random.Random(98137)
for case in range(12):
    n=1+case%5
    B=[[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            # Explicit repeated diagonal rates in half the cases.
            B[i][j]=F(1,3) if (i==j and case%2==0) else F(rng.randrange(1,5),8)
    h=F(1+case%3,4)
    got=linear_flow_enclosure(B,h,opt)
    H,J,th,tj=matrix_reference(B,h)
    for i in range(n):
        for j in range(i+1):
            require(got['H_lower'][i][j]<=H[i][j]
                    and H[i][j]+th<=got['H_upper'][i][j],
                    f'path exponential encloses row-norm reference {case}:{i},{j}')
            require(got['J_lower'][i][j]<=J[i][j]
                    and J[i][j]+tj<=got['J_upper'][i][j],
                    f'integrated path encloses row-norm reference {case}:{i},{j}')

# Exact tangent closure, where the discriminant is zero; a strict bootstrap
# would mishandle this valid bounded-Volterra-box case.
K,data=close_radius((F(1,4),),(F(1),),17)
require(data['discriminant']==0 and K==F(1,4)
        and (F(1,4)+K)**2==K,'tangent nonlinear box closes exactly')
K,data=close_radius((F(0),)*5,(F(1000),)*5)
require(K==0,'zero initial-plus-residual radius closes with arbitrary response')
K,data=close_radius((F(1),F(2)),(F(0),F(0)))
require(K==5,'zero nonlinear response uses the exact linear squared norm')

# A finite exact synthetic Riccati solution nevertheless leaves a partially
# validated interval under a deliberately tight bisection budget.
p={'ell':'1','initial_error':[str(F(q,1000)) for q in Q],
   'slabs':[{'duration':'3/2','S':'0','M':['0']*4,'residual':['0']*5}],
   'options':{'max_slab_bisections':1,'precision_bits':48}}
out=propagate(p)
require(not out['complete'] and out['validated_duration']=='3/4'
        and len(out['accepted_slabs'])==1,'second-half failure retains exactly its certified first half')
require(out['endpoint_error']==out['accepted_slabs'][-1]['endpoint'],
        'failed suffix cannot overwrite the inherited accepted endpoint')
require(sum(F(s['duration']) for s in out['accepted_slabs'])==F(out['validated_duration']),
        'partial prefix durations match the reported endpoint time')
exact_a=F(1,1000)/(1-F(284,1000)*F(out['validated_duration']))
require(all(F(x)>=q*exact_a for x,q in zip(out['endpoint_error'],Q)),
        'partial prefix encloses independently solved Riccati system')

# Add one representation jump before a slab that internally bisects. Its
# size produces the same trajectory as the previous test, without duplicate
# jump application at the second half.
p['initial_error']=['0']*5
p['slabs'][0]['jump_error']=[str(F(q,1000)) for q in Q]
jumped=propagate(p)
require(len(jumped['boundary_updates'])==1 and jumped['endpoint_error']==out['endpoint_error'],
        'jump is charged once even when its source slab bisects and partly fails')

def rejects(call,name):
    try:
        call()
    except (ValueError,TypeError):
        require(True,name)
    else:
        raise AssertionError(name)

# Regression of the direct-constructor gaps found during this audit.
rejects(lambda: certify_slab([1,0,0,0,0],
    Bounds(F(1,100),F(0),(F(0),)*4,(F(-1),)+(F(0),)*4),1),
    'direct Bounds constructor rejects negative residual before certification')
rejects(lambda: certify_slab([0]*5,
    Bounds(F(-1),F(0),(F(0),)*4,(F(0),)*5),1),
    'zero-error shortcut cannot accept a negative duration')
rejects(lambda: Bounds(.01,F(0),(F(0),)*4,(F(0),)*5),
    'direct Bounds constructor rejects binary-float duration')
rejects(lambda: close_radius((F(1),F(2)),(F(0),)),
    'radius closure rejects mismatched vector lengths')
rejects(lambda: close_radius((F(1),),(F(-1),)),
    'radius closure rejects a negative nonlinear response')

summary={'status':'PASS','checks':len(checks),
    'scope':'Synthetic comparison systems; exact rational algebra only. No NS norms or residuals were supplied.',
    'matrix_reference':'Independent ordinary matrix Taylor series with row-sum norm remainder',
    'checked':checks}
Path(__file__).with_name('independent-hierarchy-audit-check.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='checked'},indent=2))

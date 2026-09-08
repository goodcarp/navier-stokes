#!/usr/bin/env python3
"""Rational comparison-system propagation for the Pass10 error hierarchy.

The caller must supply valid whole-slab bounds for an actual approximation.
This program does NOT verify that those inputs describe Navier--Stokes.
Every accepted slab has an exact rational supersolution certificate.
No third-party package is required.
"""
from dataclasses import dataclass
from fractions import Fraction
from math import factorial, isqrt, comb
import argparse
import json

Q = (0, 1, 3, 7, 15)
F = Fraction


def rational(value):
    if isinstance(value, (F, int, str)):
        return F(value)
    raise TypeError('Use an integer, fraction/decimal string, or Fraction; floats are not certified inputs.')


def nonnegative(values, length, label):
    out = tuple(rational(v) for v in values)
    if len(out) != length or any(v < 0 for v in out):
        raise ValueError(f'{label} must contain {length} nonnegative rational bounds.')
    return out


def dyadic_up(x, bits=80):
    """Round a nonnegative rational upward to at most bits significant bits."""
    x = rational(x)
    if x < 0 or bits < 8:
        raise ValueError('Positive precision and a nonnegative value are required.')
    if x == 0:
        return F(0)
    exponent = x.numerator.bit_length() - x.denominator.bit_length()
    threshold = F(2**exponent) if exponent >= 0 else F(1, 2**(-exponent))
    if x < threshold:
        exponent -= 1
    shift = bits - 1 - exponent
    quantum = F(1, 2**shift) if shift >= 0 else F(2**(-shift))
    scaled = x / quantum
    upper_integer = (scaled.numerator + scaled.denominator - 1) // scaled.denominator
    answer = upper_integer * quantum
    assert answer >= x
    return answer


def sqrt_up(x, bits=80):
    """A dyadic relative-precision upper bound for sqrt(x), using integers only."""
    x = rational(x)
    if x < 0 or bits < 8:
        raise ValueError('Square root needs a nonnegative argument and at least eight precision bits.')
    if x == 0:
        return F(0)
    exponent = (x.numerator.bit_length() - x.denominator.bit_length()) // 2
    shift = bits - exponent
    quantum = F(1, 2**shift) if shift >= 0 else F(2**(-shift))
    scaled = x / quantum**2
    n = isqrt(scaled.numerator // scaled.denominator)
    if n*n*scaled.denominator < scaled.numerator:
        n += 1
    upper = n * quantum
    assert upper**2 >= x
    return upper


@dataclass(frozen=True)
class Options:
    precision_bits: int = 80
    exponential_terms_cap: int = 512
    radius_bisections: int = 24
    max_slab_bisections: int = 12
    minimum_duration: F = F(0)

    def __post_init__(self):
        names=('precision_bits','exponential_terms_cap','radius_bisections',
               'max_slab_bisections')
        if any(not isinstance(getattr(self,name),int) for name in names):
            raise TypeError('Iteration and precision limits must be integers.')
        minimum=rational(self.minimum_duration)
        if (self.precision_bits < 8 or self.exponential_terms_cap < 1
                or self.radius_bisections < 0 or self.max_slab_bisections < 0
                or minimum < 0):
            raise ValueError('Invalid validation options.')
        object.__setattr__(self,'minimum_duration',minimum)

    @classmethod
    def read(cls, data=None):
        d = dict(data or {})
        if 'minimum_duration' in d:
            d['minimum_duration'] = rational(d['minimum_duration'])
        ans = cls(**d)
        if (ans.precision_bits < 8 or ans.exponential_terms_cap < 1
                or ans.radius_bisections < 0 or ans.max_slab_bisections < 0
                or ans.minimum_duration < 0):
            raise ValueError('Invalid validation options.')
        return ans


@dataclass(frozen=True)
class Bounds:
    duration: F
    strain: F
    derivatives: tuple
    residual: tuple

    def __post_init__(self):
        h,strain=rational(self.duration),rational(self.strain)
        derivatives=nonnegative(self.derivatives,4,'M2 through M5')
        residual=nonnegative(self.residual,5,'residual')
        if h <= 0 or strain < 0:
            raise ValueError('Duration must be positive and S nonnegative.')
        object.__setattr__(self,'duration',h)
        object.__setattr__(self,'strain',strain)
        object.__setattr__(self,'derivatives',derivatives)
        object.__setattr__(self,'residual',residual)

    @classmethod
    def read(cls, data):
        h, strain = rational(data['duration']), rational(data['S'])
        if h <= 0 or strain < 0:
            raise ValueError('Duration must be positive and S nonnegative.')
        return cls(h, strain, nonnegative(data['M'], 4, 'M2 through M5'),
                   nonnegative(data['residual'], 5, 'residual'))

    def halve(self):
        return Bounds(self.duration/2, self.strain, self.derivatives, self.residual)


class ClosureFailure(Exception):
    def __init__(self, reason, details=None):
        super().__init__(reason)
        self.reason = reason
        self.details = details or {}


def coefficient_matrix(ell, strain, derivatives):
    ell, strain = rational(ell), rational(strain)
    derivatives = nonnegative(derivatives, 4, 'derivatives')
    if ell <= 0 or strain < 0:
        raise ValueError('ell must be positive and strain nonnegative.')
    B = [[F(0) for _ in range(5)] for _ in range(5)]
    for j in range(5):
        B[j][j] = (j+1)*strain
        for i in range(j):
            B[j][i] = comb(j+1, i)*ell**(j-i)*derivatives[j-i-1]
    return B


def exponential_tail_order(x, bits, cap):
    """Return N and an upper bound for sum_(n>N) x^n/n!, x>=0."""
    x = rational(x)
    if x < 0:
        raise ValueError('Nonnegative exponential argument required.')
    term = F(1)
    tolerance = F(1, 2**bits)
    for N in range(cap+1):
        next_term = term*x/(N+1)
        ratio = x/(N+2)
        if ratio < 1:
            tail = next_term/(1-ratio)
            if tail <= tolerance:
                return N, tail
        term = next_term
    raise ClosureFailure('exponential_series_budget_exhausted',
                         {'diagonal_argument':x, 'terms_cap':cap})


def convolution_kernel(rates, duration, N, exp_tail):
    """Positive convolution of exp(rate*t), enclosed by rational series.

    For k+1 rates, K(h)=h^k sum_n h_n(rate_i*h)/(n+k)!.
    Its remainder is at most h^k/k! times the exponential-tail bound at
    max(rate_i*h). The caller may use a larger common diagonal argument.
    """
    h = rational(duration)
    exp_tail=rational(exp_tail)
    rates = tuple(rational(r) for r in rates)
    if (not rates or h < 0 or any(r < 0 for r in rates)
            or not isinstance(N,int) or N < 0 or exp_tail < 0):
        raise ValueError('Kernel inputs must be nonnegative, with at least one rate.')
    k = len(rates)-1
    if all(r == 0 for r in rates):
        exact = h**k/factorial(k)
        return exact, exact
    homogeneous = [F(1)]+[F(0)]*N
    for rate in rates:
        x = rate*h
        for n in range(1, N+1):
            homogeneous[n] += x*homogeneous[n-1]
    lower = h**k*sum((homogeneous[n]/factorial(n+k) for n in range(N+1)), F(0))
    upper = lower+h**k*exp_tail/factorial(k)
    return lower, upper


def increasing_paths(i, j):
    if i == j:
        yield (i,)
        return
    middle = tuple(range(i+1, j))
    for mask in range(1 << len(middle)):
        yield (i,)+(tuple(v for n, v in enumerate(middle) if mask & (1 << n)))+(j,)


def linear_flow_enclosure(B, duration, options=None):
    """Entrywise enclosures of exp(B*h) and its time integral.

    B must be nonnegative and lower triangular. Large off-diagonal entries
    appear only in finite path products, never in the exponential tail.
    """
    opt = options or Options()
    h = rational(duration)
    n = len(B)
    B = [[rational(v) for v in row] for row in B]
    if (h < 0 or any(len(row) != n for row in B)
            or any(B[j][i] < 0 or (i > j and B[j][i] != 0)
                   for j in range(n) for i in range(n))):
        raise ValueError('A nonnegative lower triangular matrix is required.')
    rates = [B[j][j] for j in range(n)]
    N, tail = exponential_tail_order(max(rates, default=F(0))*h,
                                      opt.precision_bits, opt.exponential_terms_cap)
    tail = dyadic_up(tail,opt.precision_bits)
    matrices = [[[F(0) for _ in range(n)] for _ in range(n)] for _ in range(4)]
    Hlo,Hhi,Jlo,Jhi = matrices
    cache = {}
    def kernel(indices, integral):
        key = (indices, integral)
        if key not in cache:
            selected = tuple(rates[i] for i in indices)
            if integral:
                selected = (F(0),)+selected
            cache[key] = convolution_kernel(selected, h, N, tail)
        return cache[key]
    for j in range(n):
        for i in range(j+1):
            for path in increasing_paths(i,j):
                product = F(1)
                for left,right in zip(path[:-1],path[1:]):
                    product *= B[right][left]
                if product == 0:
                    continue
                lo,hi = kernel(path,False)
                Hlo[j][i] += product*lo; Hhi[j][i] += product*hi
                lo,hi = kernel(path,True)
                Jlo[j][i] += product*lo; Jhi[j][i] += product*hi
    return {'H_lower':Hlo,'H_upper':Hhi,'J_lower':Jlo,'J_upper':Jhi,
            'series_order':N,'common_exponential_tail_bound':tail,
            'largest_diagonal_argument':max(rates,default=F(0))*h}


def matvec(A,x):
    return tuple(sum((a*b for a,b in zip(row,x)),F(0)) for row in A)


def close_radius(a,b,bisections=24):
    """Exact rational K satisfying ||a+b*K||^2 <= K, or explicit failure."""
    a=tuple(rational(x) for x in a);b=tuple(rational(x) for x in b)
    if (len(a)!=len(b) or not a or any(x<0 for x in a+b)
            or not isinstance(bisections,int) or bisections < 0):
        raise ValueError('Equal nonempty nonnegative vectors and a nonnegative integer iteration limit are required.')
    A = sum(x*x for x in a)
    P = sum(x*y for x,y in zip(a,b))
    C = sum(y*y for y in b)
    if A == 0:
        return F(0), {'A':A,'P':P,'C':C,'discriminant':F(1)}
    t = 1-2*P
    discriminant = t*t-4*A*C
    if t <= 0 or discriminant < 0:
        raise ClosureFailure('nonlinear_radius_does_not_close',
            {'A':A,'P':P,'C':C,'one_minus_2P':t,'discriminant':discriminant})
    if C == 0:
        K = A/t
    else:
        lo,hi = A,2*A/t
        polynomial = lambda x: C*x*x-t*x+A
        assert polynomial(hi) <= 0
        for _ in range(bisections):
            mid=(lo+hi)/2
            if polynomial(mid) <= 0:
                hi=mid
            else:
                lo=mid
        K=hi
    assert sum((x+y*K)**2 for x,y in zip(a,b)) <= K
    return K, {'A':A,'P':P,'C':C,'discriminant':discriminant}


def certify_slab(initial, bounds, ell, options=None):
    opt = options or Options()
    initial = nonnegative(initial,5,'initial_error')
    ell = rational(ell)
    if ell <= 0:
        raise ValueError('ell must be positive.')
    if not any(initial) and not any(bounds.residual):
        return {'duration':bounds.duration,'initial':initial,'endpoint':initial,
                'whole_slab_component_bound':initial,'whole_slab_radius_squared':F(0),
                'closure_K':F(0),'zero_error_shortcut':True}
    B = coefficient_matrix(ell,bounds.strain,bounds.derivatives)
    c = sqrt_up(ell**(-5),opt.precision_bits)
    assert c*c*ell**5 >= 1
    flow = linear_flow_enclosure(B,bounds.duration,opt)
    H,J = flow['H_upper'],flow['J_upper']
    hz,jd,jq = matvec(H,initial),matvec(J,bounds.residual),matvec(J,Q)
    a = tuple(dyadic_up(x+y,opt.precision_bits) for x,y in zip(hz,jd))
    b = tuple(dyadic_up(c*y,opt.precision_bits) for y in jq)
    K,quadratic = close_radius(a,b,opt.radius_bisections)
    unrounded = tuple(x+y*K for x,y in zip(a,b))
    assert sum(x*x for x in unrounded) <= K
    endpoint = tuple(dyadic_up(x,opt.precision_bits) for x in unrounded)
    assert all(x >= y for x,y in zip(endpoint,unrounded))
    return {'duration':bounds.duration,'initial':initial,'endpoint':endpoint,
            'whole_slab_component_bound':endpoint,
            'whole_slab_radius_squared':sum(x*x for x in endpoint),
            'closure_K':K,'unrounded_endpoint':unrounded,
            'linear_initial_and_residual_upper':a,'nonlinear_response_upper':b,
            'quadratic_coefficients':quadratic,'nonlinear_constant_upper':c,
            'linear_series_order':flow['series_order'],
            'largest_diagonal_argument':flow['largest_diagonal_argument'],
            'common_exponential_tail_bound':flow['common_exponential_tail_bound'],
            'zero_error_shortcut':False}


def serialize(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {str(k):serialize(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [serialize(v) for v in value]
    return value


def reweight_endpoint(endpoint, old_ell, new_ell):
    """Exact inheritance when deliberately changing the verification length."""
    values=nonnegative(endpoint,5,'endpoint')
    old_ell,new_ell=rational(old_ell),rational(new_ell)
    if old_ell <= 0 or new_ell <= 0:
        raise ValueError('Both verification lengths must be positive.')
    ratio=new_ell/old_ell
    return tuple(x*ratio**j for j,x in enumerate(values))


def propagate(problem):
    """Validate a prefix, with adaptive subdivision and inherited errors."""
    ell=rational(problem['ell'])
    if ell <= 0:
        raise ValueError('ell must be positive.')
    initial=nonnegative(problem['initial_error'],5,'initial_error')
    slabs=[Bounds.read(row) for row in problem['slabs']]
    jumps=[nonnegative(row.get('jump_error',[0]*5),5,'jump_error')
           for row in problem['slabs']]
    opt=Options.read(problem.get('options'))
    current=initial;time=F(0);accepted=[];attempts=[];boundary_updates=[]

    def advance(bounds,source_index,depth):
        nonlocal current,time
        start=time
        try:
            cert=certify_slab(current,bounds,ell,opt)
        except ClosureFailure as failure:
            attempts.append({'start':start,'duration':bounds.duration,
                             'depth':depth,'reason':failure.reason,
                             'details':failure.details})
            half=bounds.duration/2
            if depth >= opt.max_slab_bisections or half < opt.minimum_duration:
                return False
            return (advance(bounds.halve(),source_index,depth+1)
                    and advance(bounds.halve(),source_index,depth+1))
        cert.update({'start':start,'source_slab_index':source_index,'bisection_depth':depth,
                     'supplied_bounds':{'S':bounds.strain,'M':bounds.derivatives,
                                        'residual':bounds.residual}})
        accepted.append(cert)
        current=cert['endpoint'];time+=bounds.duration
        return True

    complete=True
    for index,bounds in enumerate(slabs):
        if any(jumps[index]):
            before=current
            current=tuple(dyadic_up(x+y,opt.precision_bits)
                          for x,y in zip(current,jumps[index]))
            boundary_updates.append({'time':time,'source_slab_index':index,
                'before':before,'jump_error':jumps[index],'after':current})
        if not advance(bounds,index,0):
            complete=False
            break
    result={'status':'COMPARISON_INTERVAL_VALIDATED' if complete else 'COMPARISON_PREFIX_ONLY',
        'complete':complete,'ell':ell,'initial_error':initial,
        'requested_duration':sum((s.duration for s in slabs),F(0)),
        'validated_duration':time,'endpoint_error':current,
        'endpoint_radius_squared_upper':sum(x*x for x in current),
        'accepted_slabs':accepted,'failed_attempts':attempts,
        'boundary_updates':boundary_updates,
        'scope':'Conditional comparison-system certificate only. Supplied strain, '
                'higher-derivative, residual and initial-error bounds have not been '
                'verified against any Navier--Stokes approximation by this program. '
                'A failed suffix is not certified and is not a blowup claim.'}
    return serialize(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',help='JSON problem with exact rational/decimal strings')
    parser.add_argument('--output')
    args=parser.parse_args()
    with open(args.input) as handle:
        problem=json.load(handle,parse_float=str)
    result=propagate(problem)
    text=json.dumps(result,indent=2)+'\n'
    if args.output:
        with open(args.output,'w') as handle:handle.write(text)
    print(text,end='')


if __name__=='__main__':
    main()

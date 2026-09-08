#!/usr/bin/env python3
"""Implementation checks of the interval cutoff and range quadrature.

The known integral identities are mathematical independent checks; no sampled
comparison is used as an enclosure proof for the fluid integrals.
"""
import json
from fractions import Fraction as F
from interval_local_integrals import iv,interval,cutoff_jets,range_quadrature,exact_endpoints


def rational_endpoint(data):
    value = F((-1 if data['sign'] else 1)*int(data['mantissa']))
    exponent = data['exponent']
    return value*(2**exponent) if exponent >= 0 else value/F(2**(-exponent))


def contains_result(data,value):
    endpoints = data['exact_dyadic_bounds']
    return rational_endpoint(endpoints['lower']) <= value <= rational_endpoint(endpoints['upper'])


def main():
    iv.dps=35
    for x,expected in [(F(0),(F(1),F(0),F(0))),
                       (F(1,4),(F(1),F(0),F(0))),
                       (F(5,8),(F(1,2),F(-8,3),F(0))),
                       (F(1),(F(0),F(0),F(0)))]:
        for enclosure,answer in zip(cutoff_jets(interval(x)),expected):
            assert contains_result({'exact_dyadic_bounds':exact_endpoints(enclosure)},answer)
    # Cutoff symmetry on its transition gives int_0^1 f = 1/4+(3/4)/2 = 5/8.
    grid=[[F(0),F(1,4),F(1,2),F(3,4),F(1)]]
    mass=range_quadrature(lambda x:cutoff_jets(x,0)[0],grid,256,.005)
    derivative=range_quadrature(lambda x:cutoff_jets(x,1)[1],grid,512,.01)
    assert contains_result(mass,F(5,8))
    assert contains_result(derivative,F(-1))
    print(json.dumps(dict(status='PASS',
        checks=['exact plateau and midpoint cutoff jets',
                'known cutoff mass 5/8 enclosed',
                'fundamental-theorem derivative integral -1 enclosed'],
        cutoff_mass=mass,cutoff_derivative_integral=derivative)),flush=True)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Range-integral certificate for the unchanged C-infinity core cutoff.

The analytic coefficient reduction is established separately in
../pass4/core-pressure-derivative.md. This script certifies the remaining
one-dimensional integral using mpmath.iv outward-rounded ranges; it is not
a formal-kernel proof of Navier--Stokes dynamics.
"""
import argparse
import json
from fractions import Fraction as F
from interval_local_integrals import iv, interval, cutoff_jets, range_quadrature, exact_endpoints


def endpoint_fraction(encoded):
    result = F((-1 if encoded['sign'] else 1)*int(encoded['mantissa']))
    exponent = encoded['exponent']
    return result*F(2)**exponent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--panels', type=int, default=1024)
    parser.add_argument('--dps', type=int, default=35)
    args = parser.parse_args()
    iv.dps = args.dps

    def integrand(s):
        p, dp = cutoff_jets(s, 1)
        return dp*(p+interval(F(2,3))*s*dp)**2

    quadrature = range_quadrature(integrand,
        [[F(1,4),F(1,2),F(5,8),F(3,4),F(1)]], args.panels, .05)
    exact = quadrature['exact_dyadic_bounds']
    integral = interval(endpoint_fraction(exact['lower']),
                        endpoint_fraction(exact['upper']))
    twb = interval(F(44,35))-interval(F(16,5))*integral
    bounds = exact_endpoints(twb)
    upper_below_four = endpoint_fraction(bounds['upper']) < 4
    assert upper_below_four, 'This resolution does not certify tWb < 4.'
    print(json.dumps(dict(
        status='CERTIFIED_WITH_MPMATH_IV',
        cutoff='chi(s)=1 for s<=1/4; expit(1/t-1/(1-t)), t=(4s-1)/3, for 1/4<s<1; 0 for s>=1',
        coefficient='44/35-(16/5)*integral chi_prime(s)*(chi(s)+(2*s/3)*chi_prime(s))^2 ds',
        integration_domain=['1/4','1'], dps=args.dps,
        integral=quadrature,
        tWb=dict(interval=str(twb),exact_dyadic_bounds=bounds,
                 strictly_below_four=upper_below_four),
        limitations='Analytic coefficient identity and trusted interval library are assumed; no full pressure gate or dynamical regeneration is certified.'
    ), indent=2))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""One-dimensional range certificate for an upper bound on nu*d_w.

No pressure solve: use d_w<=2 Qmax ||grad w||_L2^2 and separate the
outer chiral seed's radial/axial factors. Constants are rational.
"""
import sys,json,argparse
from pathlib import Path
from fractions import Fraction as F
here=Path(__file__).resolve().parent
for base in (here,here.parent,here.parent.parent):
    for candidate in (base/'pass5',base/'pass5'/'certificate',base):
        if (candidate/'interval_local_integrals.py').exists():
            sys.path.insert(0,str(candidate))
from interval_local_integrals import iv,interval,cutoff_jets,range_quadrature,exact_endpoints

def rational_endpoint(encoded):
    return F((-1 if encoded['sign'] else 1)*int(encoded['mantissa']))*F(2)**encoded['exponent']

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--panels',type=int,default=1024)
    parser.add_argument('--dps',type=int,default=35);args=parser.parse_args();iv.dps=args.dps
    r0,delta=F(21,100),F(7,50)
    def radial(r):
        y=(r-interval(r0))/interval(delta)
        p,p1,p2=cutoff_jets(y*y,2)
        er=2*y*p1/interval(delta)
        err=(2*p1+4*y*y*p2)/interval(delta*delta)
        cosine=err+er/r-(400+16/r**2)*p
        sine=20*(2*er+p/r)
        return r*(cosine**2+sine**2)
    points=[r0-delta,r0-delta*F(3,4),r0-delta/2,r0,
            r0+delta/2,r0+delta*F(3,4),r0+delta]
    result=range_quadrature(radial,[points],args.panels,.1)
    exact=result['exact_dyadic_bounds']
    radial_interval=interval(rational_endpoint(exact['lower']),rational_endpoint(exact['upper']))
    # Axial mass <=12/5, axial derivative energy <=160/3;
    # radial first derivative/phase energy <=6764/75 (pass6 bound).
    gradient_over_pi=interval(F(12,5))*radial_interval+interval(F(6764,75)*F(160,3))
    viscosity=interval(F(6,1000)*F(5,17)**5)*gradient_over_pi
    threshold=F(5093339,210000)
    upper=rational_endpoint(exact_endpoints(viscosity)['upper'])
    radial_upper=rational_endpoint(exact['upper'])
    assert radial_upper<120000, 'Requested coarse radial certificate not reached.'
    assert upper<4, 'Requested viscous coefficient bound not reached.'
    radial_upper=rational_endpoint(exact['upper'])
    output=dict(status='PASS' if upper<4 and radial_upper<120000 else 'INCONCLUSIVE',
        viscous_upper_strictly_below_four=upper<4,
        radial_upper_strictly_below_120000=radial_upper<120000,
        radial_second_derivative_energy=result,
        viscous_upper_bound_interval=str(viscosity),
        viscous_upper_bound_exact=exact_endpoints(viscosity),
        sufficient_threshold=str(threshold),
        below_threshold=upper<threshold,
        radial_upper_below_120000=radial_upper<120000,
        viscous_upper_below_four=upper<4,
        method='one-dimensional outward range quadrature plus analytic energy/kernel bounds',
        scope='upper bound for nu*d_w, not a direct pressure solve; Euler kernel proof separate')
    print(json.dumps(output,indent=2))

if __name__=='__main__':main()

"""Exact kernel-norm and compact-energy bounds for the initial image tail.

Uses the prior certified cutoff derivative bound |chi'|<4, monotonicity,
and Pass8 R0<=392/5, Z0<=12/5. Does not certify retained image quadrature.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import sympy as s
from initial_periodic_images import analytic_tail_bound


def run():
    t=s.symbols('t',real=True)
    rr=105*t*(1-t)-12;zz=105*t*t-90*t+9
    rz_square=225*t*(1-t)*(7*t-3)**2
    det_minus=s.factor((24-rr)*(24-zz)-rz_square)
    det_plus=s.factor((24+rr)*(24+zz)-rz_square)
    assert s.expand(det_minus-180*(1-t)*(t+3))==0
    assert s.expand(det_plus-36*(16-5*(t-1)**2))==0
    assert s.expand((24-rr)-(105*(t-s.Rational(1,2))**2+s.Rational(39,4)))==0
    assert s.expand((24+rr)-(12+105*t*(1-t)))==0
    # For 0<=t<=1 both 2x2 matrices are PSD by positive first pivot and
    # the displayed nonnegative determinants. The theta eigenvalue 3-15t
    # is also in [-24,24]. Hence ||4pi R^5 D_zz D_ij N||op<=24.
    mu,S,P,P1=s.symbols('mu S P P1',real=True)
    speed_square=S*((1-mu**2)*(P+2*S*mu**2*P1)**2
                     +4*mu**2*(P+S*(1-mu**2)*P1)**2)
    angular_average=s.integrate(speed_square,(mu,-1,1))/2
    assert s.simplify(angular_average-(2*S*P*P+s.Rational(8,5)*S*S*P*P1
                                      +s.Rational(8,15)*S**3*P1*P1))==0
    # Compact-support integration cancels the P^2 and P*P1 terms:
    assert 2+s.Rational(8,5)*(-s.Rational(5,4))==0
    # Hence E(S_Phi)=(8pi/15) integral s^(7/2) Phi'(s)^2 ds.
    pi_upper=F(22,7);c=F(7,5);av=F(63,200)
    weighted_derivative=4+F(5,2)**7*F(16,25)*c*c+6**7*F(3,11)*c*c
    ES=F(8,15)*pi_upper*weighted_derivative
    EW=F(484,135)*pi_upper
    Ev=F(9,20)*pi_upper*av**4
    K=pi_upper/F(32)*F(392,5)*F(12,5)
    Ebound=ES+EW+1020**2*Ev+K
    assert Ebound==F(2015718896182553,7560000000)
    assert Ebound<270000
    tails=[dict(period=L,count=N,bound_exact=str(analytic_tail_bound(L,N)))
           for L in (16,32,64) for N in (8,16,32)]
    return dict(status='PASS',kernel_operator_bound='24/(4*pi*R^5)',
        exact_energy_upper=str(Ebound),coarse_energy_upper=270000,
        absolute_tail_bound='6E/[pi L(LN-6)^4] < 540000/[L(LN-6)^4]',
        tails=tails,dependencies=['cutoff monotonicity and |chi derivative|<4 from pass6 cutoff certificate',
            'Pass8 R0<=392/5 and Z0<=12/5; unchanged supports and lambda=1/4'],
        scope='Exact omitted-image bound; retained image integrals still require quadrature enclosures for a full certificate.')


if __name__=='__main__':
    result=run()
    Path(__file__).with_name('periodic-image-tail-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

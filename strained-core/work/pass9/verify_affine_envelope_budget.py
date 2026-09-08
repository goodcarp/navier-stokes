"""Exact algebra and outward-rounded envelope moment for an affine linear model.

This certifies a reduced-model finite-time energy inequality, not NS evolution.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, json, sys
import sympy as s

here = Path(__file__).resolve().parent
for dependency_dir in (here.parent/'pass5', here.parents[1]/'pass5/certificate'):
    if (dependency_dir/'interval_local_integrals.py').is_file():
        sys.path.insert(0, str(dependency_dir))
        break
else:
    raise FileNotFoundError('Expected sibling pass5/interval_local_integrals.py or pass5/certificate/interval_local_integrals.py')
from interval_local_integrals import cutoff_jets, interval, range_quadrature

DEPENDENCIES = dict(
    cutoff='pass6/checks/outer-mean-maximum-certificate.json: -4 < chi derivative <= 0; chi(5/8)=1/2 follows from the cutoff definition.',
    radius='pass8/checks/fixed-circle-acceleration-certificate.json: 1077377876667/5000000000000 < r_star < 269344469167/1250000000000.',
    amplitude='pass8/checks/leading-envelope-certificate.json: 900 < retuned A < 1020.',
    derived_bounds='h=4/r_star lies between 1000/63 and 1600/63; (r_star/(63/200))^2 < 5/8, hence Omega=A*chi((r_star/(63/200))^2) > 450.',
    layout='Certificate paths are relative to the euler-ns-transfer package root; work/ uses the same filenames directly inside each pass directory.',
    scope='These are declared prior certified input bounds. This checker certifies the new envelope moment and rational energy budget; no PDE quadrature is used.')


def endpoint(d):
    return (-1 if d['sign'] else 1)*F(int(d['mantissa']))*F(2)**int(d['exponent'])


def run(panels=2048):
    x, a, H = s.symbols('x a H', real=True)
    N = (x-a)**2+H
    even = lambda f: s.expand((f+f.subs(x,-x))/2)
    gain = s.cancel(even(N**2)/(x*x+H)-(x*x+a*a+H))
    assert s.simplify(gain-a*a*(5*x*x+a*a+H)/(x*x+H)) == 0
    numerator = even(N**2*(x*x-a*x+H+a*a/3))
    quotient = x**4+(2*H+31*a*a/3)*x*x+H*H+7*H*a*a/3+7*a**4
    remainder = a**4*(a*a-16*H)/3
    assert s.expand(numerator-(x*x+H)*quotient-remainder) == 0

    # x is dimensionless displacement from envelope center: r=.14+.1*x.
    # A2 = integral (e_rr)^2 dr, retaining both symmetric transition halves.
    def a2(x):
        _, f1, f2 = cutoff_jets(x*x,2)
        return 2000*(2*f1+4*x*x*f2)**2
    moment = range_quadrature(a2, [[F(1,2),F(3,4),F(1)]],
                              max_panels=panels, relative_width=F(1,4))
    A2upper = endpoint(moment['exact_dyadic_bounds']['upper'])
    assert A2upper < 800000

    # Prior exact cutoff facts: support width .2, plateau width .1,
    # |chi'|<4. Hence A0 in [.1,.2], A1 <= ||e'||inf TV(e) <= 80*2.
    Hmin, Hmax = F(1000,63)**2, F(1600,63)**2
    assert F(400) < 16*Hmin
    assert F(400) < 4*Hmin  # Gain ratio increases with eta^2.
    ratio_min = 1+400/Hmax
    assert F(400,4)*F(1,10)*ratio_min > 16
    E0upper = (160+(400+Hmax)*F(1,5))/4
    assert E0upper < 93
    tauupper = F(7,5000)
    assert 2*F(7,5)*tauupper < F(1,250)
    viscous_upper = F(1,2000)*tauupper*(800000+(2*Hmax+F(31*400,3))*160
                 +(Hmax**2+F(7*400,3)*Hmax+7*400**2)*F(1,5))
    assert viscous_upper < F(3,2)
    net_lower = F(249,250)*(16-F(3,2))-E0upper/250
    assert net_lower > F(3,20)*E0upper
    return dict(status='PASS', scope='exact rotating-affine linear envelope model only',
        envelope_second_derivative_moment=moment,
        A2_coarse_upper=800000, inviscid_additive_gain_lower=16,
        viscosity_loss_upper=str(viscous_upper), E0_upper=str(E0upper),
        net_additive_gain_lower=str(net_lower), relative_gain_lower='3/20',
        evaluation_time='20/(2*Omega*h)', time_upper=str(tauupper),
        dependencies=DEPENDENCIES)


if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--panels', type=int, default=2048)
    args=ap.parse_args()
    result=run(args.panels)
    (here/'affine-envelope-budget-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

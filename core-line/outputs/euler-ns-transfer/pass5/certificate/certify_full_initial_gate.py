#!/usr/bin/env python3
"""Assemble the finite initial-feedback certificate with exact rational bounds.

Analytical identities are documented/audited in the accompanying notes.
--recompute reruns every interval integral used here. The result certifies
one initial NS derivative inequality, not a singularity or a return map.
"""
import argparse,json,subprocess,sys
from pathlib import Path
from fractions import Fraction as F

HERE=Path(__file__).resolve().parent


def number(endpoint):
    return F((-1 if endpoint['sign'] else 1)*int(endpoint['mantissa']))*F(2)**int(endpoint['exponent'])


def bounds(record):
    e=record['exact_dyadic_bounds']
    return number(e['lower']),number(e['upper'])


def run_json(name,*args):
    result=subprocess.run([sys.executable,str(HERE/name),*args],check=True,capture_output=True,text=True)
    return json.loads(result.stdout)


def assemble(recompute=False):
    if recompute:
        strain=run_json('certify_radial_strain.py','--c','7/5','--panels','1024')
        twb=run_json('certify_core_twb.py','--panels','1024')
        cv=run_json('interval_local_integrals.py','--kind','Cv','--panels','2048','--relative-width','0')
        dv=run_json('interval_local_integrals.py','--kind','dv','--panels','2048','--relative-width','0')
    else:
        strain=json.loads((HERE/'radial-strain-interval-1024.json').read_text())
        twb=json.loads((HERE/'core-twb-interval.json').read_text())
        local=json.loads((HERE/'local-interval-results.json').read_text())
        assert local['profile']=={'radial_scale':'63/200','axial_scale':'9/20','distance':'4'}
        cv,dv=local['Cv'],local['dv']
    assert strain['c']=='7/5' and strain['panels_per_transition']==1024
    Cs=bounds(strain); Tw=bounds(twb['tWb']); Cv=bounds(cv); Dv=bounds(dv)
    assert Cs[1]<66
    assert Tw[1]<4
    assert Cv[0]>0 and Dv[0]>0
    assert Dv[1]/Cv[0]<901
    b=F(1); omega=F(1); c=F(7,5); nu=F(1,1000)
    H=F(38,7)*b*b+F(2,5)*omega*omega
    assert H==F(204,35)
    # These are exact analytical bounds from the kernel notes/checkers.
    kcore=F(-12,5)
    kpump=F(-31,2)
    T_upper=32+66+4+H*(kcore+c*kpump+nu*901)
    assert T_upper==F(-290649,8750)<0
    beta_second_lower=-T_upper/(2*omega)
    reservoir_log_lower=F(2381,357)*c-4*b-nu*901
    assert reservoir_log_lower>4
    # Exact tuning of A uses the actual positive integral, not a rounded A.
    neutral=-F(18,7)*b*b+F(2,5)*omega*omega-H
    assert neutral==-8*b*b
    return dict(
        status='INITIAL_FEEDBACK_CERTIFIED_USING_ANALYTICAL_REDUCTIONS_AND_MPMATH_IV',
        interval_integrals_recomputed=recompute,
        parameters=dict(b=str(b),omega=str(omega),c=str(c),nu=str(nu),A_squared='(204/35)/C_v'),
        coarse_certified_bounds=dict(C300_Phi_upper='66',tWb_upper='4',dv_over_Cv_upper='901',
                                    core_kernel_upper=str(kcore),pump_kernel_upper=str(kpump)),
        full_gate_upper_exact=str(T_upper),full_gate_upper_display=float(T_upper),
        beta_second_derivative_lower_exact=str(beta_second_lower),
        reservoir_log_derivative_lower_exact=str(reservoir_log_lower),
        exact_interval_records=dict(strain=strain,tWb=twb,Cv=cv,dv=dv),
        viscosity_one_scaling='U(t,x)=1000 u(1000 t,x), P(t,x)=1000000 p(1000 t,x)',
        limitations=['One initial pressure derivative and local feedback signs only.',
                     'Analytical identities and outward-rounded interval implementation are part of the certificate trust base.',
                     'No formal kernel replay, quantified return map, inherited profile closure, or NS blowup proof.'])


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--recompute',action='store_true')
    args=ap.parse_args()
    print(json.dumps(assemble(args.recompute),indent=2))

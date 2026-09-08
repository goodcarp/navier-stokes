#!/usr/bin/env python3
"""Whole-pressure upper bound for the new leading-envelope initial family."""
import sys,json,argparse
from pathlib import Path
from fractions import Fraction as F
here=Path(__file__).resolve().parent
for p in [here.parent/'pass5',here.parents[1]/'pass5/certificate',here.parents[2]/'outputs/euler-ns-transfer/pass5/certificate']:
    if (p/'interval_local_integrals.py').exists():sys.path.insert(0,str(p));break
from interval_local_integrals import iv,interval,cutoff_jets,range_quadrature,exact_endpoints
def endpoint(x):return F((-1 if x['sign'] else 1)*int(x['mantissa']))*F(2)**x['exponent']
def recover(res):
    x=res['exact_dyadic_bounds'];return interval(endpoint(x['lower']),endpoint(x['upper']))
def pack(x):return dict(display=str(x),exact=exact_endpoints(x))

def run(panels=1024):
    iv.dps=35
    leading=json.loads((here/'leading-envelope-certificate.json').read_text())
    assert leading['status']=='PASS'
    assert leading['parameters']==dict(radial_center='7/50',radial_width='1/10',m=4,k=-20,seed_amplitude_interval=['1/8','1/4'],nu='1/1000')
    assert F(leading['R0_upper'])==F(392,5)
    assert F(leading['Cw_upper'])<F(3,2)
    assert leading['retuned_A_bounds']==['900','1020']
    assert F(leading['mean_derivative_lower'])>11
    assert F(leading['fluctuation_energy_fractional_derivative_lower'])>5
    assert endpoint(leading['R1']['exact_dyadic_bounds']['upper'])<200000
    def radial(r,which):
        y=(r-interval(F(7,50)))/interval(F(1,10))
        e,e1=cutoff_jets(y*y,1);er=20*y*e1
        if which=='A':return r**5*((er-e/r)**2+400*e*e)
        return r**3*e*e
    def axial(z,which):
        y=(z-4)/interval(F(3,5))
        q,q1=cutoff_jets(y*y,1);qz=2*y*q1/interval(F(3,5))
        return 2*(q*q/z**14 if which=='A' else qz*qz/z**12)
    rp=list(map(F,['.04','.065','.09','.14','.19','.215','.24']))
    zp=list(map(F,['3.4','3.55','3.7','4','4.3','4.45','4.6']))
    res={}
    for which in ['A','B']:
        res['radial_'+which]=range_quadrature(lambda r,w=which:radial(r,w),[rp],panels,.05)
        res['axial_'+which]=range_quadrature(lambda z,w=which:axial(z,w),[zp],panels,.05)
    I={key:recover(val) for key,val in res.items()}
    nd=4*iv.sqrt(iv.pi)*(45/(2*iv.pi)*iv.sqrt(I['radial_A']*I['axial_A'])+
                         15/iv.pi*iv.sqrt(I['radial_B']*I['axial_B']))
    # The exact m4 pressure projection energy inequality uses no finite box:
    # |Epress| <= (4/m) ||r div(Qw)|| ||(w.grad)v||.
    # Horizontal v gradient: |V/r|<=1 and |Vr|<=max(1,8xi)<=5.
    xi_max=(F(6,25)/F(63,200))**2
    assert 8*xi_max<5
    R0,Z0=F(392,5),F(12,5)
    W=iv.sqrt(iv.pi*interval(R0*Z0))
    pressure_abs=5*nd*W
    # The leading chirality reverses the local term's sign, so retain its
    # positive upper bound, using (6-t)/(1+t)^(9/2)<=6.
    local_abs=interval(135*80*(F(6,25)**3-F(1,25)**3)/3*Z0/F(17,5)**7)
    Eabs=local_abs+pressure_abs
    Eok=bool(Eabs.b<interval(F(1,3)).a)
    nu_d_bound=F(6,1000)*F(5,17)**5*(200000*Z0+R0*F(160,3))
    assert nu_d_bound<7
    m0=F(290649,8750)
    margin=m0-(F(3999,1000)*F(3,2)+7+F(1020,3))/16
    assert margin>11
    return dict(status='PASS' if Eok else 'INCONCLUSIVE',panels=panels,dps=iv.dps,
        dependency='leading-envelope-certificate.json',range_integrals=res,
        d_weighted_norm_upper_expression=pack(nd),full_cross_pressure_absolute_upper_expression=pack(pressure_abs),
        local_absolute_upper=str(local_abs),E_absolute_upper_expression=pack(Eabs),
        E_absolute_upper_strictly_below_one_third=Eok,
        nu_d_rational_upper=str(nu_d_bound),nu_d_upper_below_seven=True,
        common_T_negative_margin=str(margin),beta_second_lower=str(margin/2),
        scope='Full actual initial pressure-derivative gate for the separate leading-envelope family, including viscosity and the nonlocal strong-swirl cross pressure; no sustained gain or return.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--panels',type=int,default=1024);ap.add_argument('--output');args=ap.parse_args()
    data=run(args.panels)
    if args.output:Path(args.output).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2),flush=True)

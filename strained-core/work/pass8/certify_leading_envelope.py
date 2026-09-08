"""Exact interval certificate for a leading chiral envelope's initial gain."""
import sys,json,argparse,os
from pathlib import Path
from fractions import Fraction as F
here=Path(__file__).resolve().parent
for p in [here.parent/'pass5',here.parents[1]/'pass5/certificate',here.parents[2]/'outputs/euler-ns-transfer/pass5/certificate']:
    if (p/'interval_local_integrals.py').exists():sys.path.insert(0,str(p));break
from interval_local_integrals import iv,interval,cutoff_jets,range_quadrature,exact_endpoints

def endpoint(e):
    x=F(int(e['mantissa']))*F(2)**e['exponent']
    return -x if e['sign'] else x

def certify(panels=1024):
    iv.dps=35
    r0=F(7,50); dr=F(1,10); av=F(63,200)
    root=json.loads((here/'fixed-circle-acceleration-certificate.json').read_text())
    cutoff_path=next(p for p in [here.parent/'pass6/outer-mean-maximum-certificate.json',here.parents[1]/'pass6/checks/outer-mean-maximum-certificate.json'] if p.exists())
    cutoff_certificate=json.loads(cutoff_path.read_text())
    d1=cutoff_certificate['cutoff_jets']['first_derivative_exact']
    assert endpoint(d1['lower'])>-4 and endpoint(d1['upper'])<=0
    rlo,rhi=map(F,root['unique_global_radial_maximum_interval'])
    r=interval(rlo,rhi)
    def seed_jets(r,order):
        y=(r-interval(r0))/interval(dr)
        jets=cutoff_jets(y*y,order)
        e=jets[0]
        er=2*y*jets[1]/interval(dr) if order else None
        err=(2*jets[1]+4*y*y*jets[2])/interval(dr*dr) if order>=2 else None
        return e,er,err
    e,er,_=seed_jets(r,1)
    torque=-40/r*(e*e+2*r*e*er)
    assert torque.a>1950
    G2=root['G_jets'][2]['exact']
    assert endpoint(G2['lower'])>-19

    def r1(r):
        e,er,err=seed_jets(r,2)
        return r*((err+er/r-(400+16/r**2)*e)**2+400*(2*er+e/r)**2)
    def extraction(r):
        e,_,_=seed_jets(r,0)
        C1=cutoff_jets((r/interval(av))**2,1)[1]
        Cr=2*r*C1/interval(av*av)
        return -r*Cr*e*e
    points=sorted(set([r0-dr,r0-3*dr/4,r0-dr/2,r0,av/2,
                       r0+dr/2,F(1,5),F(21,100),F(11,50),r0+dr]))
    R1=range_quadrature(r1,[points],max_panels=panels,relative_width=.1)
    Ir=range_quadrature(extraction,[points],max_panels=panels,relative_width=.1)
    R1_upper=endpoint(R1['exact_dyadic_bounds']['upper'])
    Ir_lower=endpoint(Ir['exact_dyadic_bounds']['lower'])
    assert R1_upper<200000
    assert Ir_lower>F(1,50)

    # Crude exact energy and pressure bounds, retaining axial viscosity.
    # log6<9/5 follows from exp(9/5)>6 using a positive Taylor truncation.
    from math import factorial
    assert sum(F(9,5)**j/factorial(j) for j in range(9))>6
    rmin=F(1,25);rmax=F(6,25)
    R0_upper=rmax*80*2+400*(rmax*rmax-rmin*rmin)/2+16*F(9,5)
    assert R0_upper==F(392,5)
    Z0=F(12,5);Z1=F(160,3);overlap=F(9,10)
    Cw_upper=3*F(5,17)**5*R0_upper*Z0
    assert Cw_upper<F(3,2)
    H=F(204,35);Cvlo=F(113,20000000);Cvhi=F(279,40000000)
    assert H-F(3,2)/16>900**2*Cvhi
    assert H<1020**2*Cvlo
    Gt_lower=F(1950,64)-F(19,1000)*1020
    assert Gt_lower>11
    energy_margin=80*900*F(1,50)*overlap-F(7,5)*R0_upper*Z0\
                  -F(1,1000)*(200000*Z0+R0_upper*Z1)
    fractional=energy_margin/(R0_upper*Z0/2)
    assert energy_margin>548 and fractional>5
    return dict(status='PASS',scope='actual initial mean maximum and net fluctuation-energy gain; exact central pressure neutrality; this subcertificate does not assess the full pressure derivative; companion leading-pressure certificate does',
        parameters=dict(radial_center=str(r0),radial_width=str(dr),m=4,k=-20,
                        seed_amplitude_interval=['1/8','1/4'],nu='1/1000'),
        root_interval=[str(rlo),str(rhi)],root_certificate_dependency='fixed-circle-acceleration-certificate.json',
        cutoff_derivative_certificate_dependency=os.path.relpath(cutoff_path,here),
        seed_value=dict(display=str(e),exact=exact_endpoints(e)),
        seed_derivative=dict(display=str(er),exact=exact_endpoints(er)),
        unit_mean_torque=dict(display=str(torque),exact=exact_endpoints(torque)),
        R1=R1,radial_extraction=Ir,R1_coarse_upper=200000,extraction_coarse_lower='1/50',
        R0_upper=str(R0_upper),Cw_upper=str(Cw_upper),retuned_A_bounds=['900','1020'],
        mean_derivative_lower=str(Gt_lower),
        fluctuation_energy_derivative_lower_over_pi_lambda_squared=str(energy_margin),
        fluctuation_energy_fractional_derivative_lower=str(fractional),
        trust_base='existing certified root interval and cutoff-jet ranges; new one-dimensional outward-rounded range quadratures')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--panels',type=int,default=1024)
    ap.add_argument('--output',default=str(here/'leading-envelope-certificate.json'))
    args=ap.parse_args(); out=certify(args.panels)
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)

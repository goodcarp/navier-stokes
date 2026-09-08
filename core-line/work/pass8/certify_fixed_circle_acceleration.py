#!/usr/bin/env python3
"""Rigorous initial fixed-circle Gtt budget, using a full-space pressure bound.

The |pressure radial error|<5 premise is proved independently by the m4
harmonic slab estimate from the unchanged pass7 residual certificate.
"""
import sys,json,argparse,time
from pathlib import Path
from fractions import Fraction as F
here=Path(__file__).resolve().parent
for p in [here.parent/'pass5',here.parents[1]/'pass5/certificate',here.parents[2]/'outputs/euler-ns-transfer/pass5/certificate']:
    if (p/'interval_local_integrals.py').exists():sys.path.insert(0,str(p));break
from interval_local_integrals import iv,interval,cutoff_jets,exact_endpoints

AV=F(63,200);R0=F(21,100);DELTA=F(7,50)
RLO=F('0.2154755753334');RHI=F('0.2154755753336')

def pack(x):return dict(display=str(x),exact=exact_endpoints(x))
def mesh(points,n):
    total=points[-1]-points[0]
    for a,b in zip(points[:-1],points[1:]):
        N=max(1,int(n*(b-a)/total)+1)
        for i in range(N):yield a+(b-a)*i/N,a+(b-a)*(i+1)/N

def cjets(r,order=2):
    x=(r/interval(AV))**2
    if order<=2:return cutoff_jets(x,order)
    t=(4*x-1)/3
    assert t.a>0 and t.b<1
    p=cutoff_jets(x,0)[0];h=p*(1-p)
    z1=(-t**-2-(1-t)**-2)*interval(F(4,3))
    z2=(2*t**-3-2*(1-t)**-3)*interval(F(4,3))**2
    z3=(-6*t**-4-6*(1-t)**-4)*interval(F(4,3))**3
    z4=(24*t**-5-24*(1-t)**-5)*interval(F(4,3))**4
    d1=h*z1
    d2=h*((1-2*p)*z1*z1+z2)
    d3=h*((1-6*p+6*p*p)*z1**3+3*(1-2*p)*z1*z2+z3)
    d4=h*((1-14*p+36*p*p-24*p**3)*z1**4+
          6*(1-6*p+6*p*p)*z1*z1*z2+
          3*(1-2*p)*z2*z2+4*(1-2*p)*z1*z3+z4)
    return p,d1,d2,d3,d4

def gjets(r,order=2):
    a=interval(AV);d=cjets(r,order);C=d[0]
    values=[r*r*C,2*r*C+2*r**3/a**2*d[1]]
    if order>=2:values.append(2*C+10*r*r/a**2*d[1]+4*r**4/a**4*d[2])
    if order>=3:values.append(24*r/a**2*d[1]+36*r**3/a**4*d[2]+8*r**5/a**6*d[3])
    if order>=4:values.append(24/a**2*d[1]+156*r*r/a**4*d[2]+112*r**4/a**6*d[3]+16*r**6/a**8*d[4])
    return values

def prove_maximum():
    # Prior pass6 restricts every global radial maximum to (a/2,4a/5).
    # Prove signs outside [.2,.23] and strict concavity inside it.
    for lo,hi in mesh([AV/2,F(1,5)],512):
        assert gjets(interval(lo,hi),1)[1].a>0
    for lo,hi in mesh([F(23,100),AV*F(4,5)],512):
        assert gjets(interval(lo,hi),1)[1].b<0
    for lo,hi in mesh([F(1,5),F(23,100)],512):
        assert gjets(interval(lo,hi),2)[2].b<0
    assert gjets(interval(RLO),1)[1].a>0
    assert gjets(interval(RHI),1)[1].b<0

def green_moments(points,power,n):
    vals=[iv.mpf(0),iv.mpf(0)]
    for lo,hi in mesh(points,n):
        r=interval(lo,hi);dr=interval(hi-lo)
        C,C1=cjets(r,1);Cp=2*r/interval(AV)**2*C1
        e=cutoff_jets(((r-interval(R0))/interval(DELTA))**2,0)[0]
        v=r**power*Cp*e
        vals[0]+=dr*v*iv.cos(20*r);vals[1]+=dr*v*iv.sin(20*r)
    return vals

def run(panels=2048):
    start=time.time();iv.dps=40
    prove_maximum()
    r=interval(RLO,RHI);C=cjets(r,0)[0]
    left=green_moments([AV/2,RLO],4,panels)
    right=green_moments([RHI,F(7,25),AV],-4,panels)
    C1=cjets(r,1)[1];Cp=2*r/interval(AV)**2*C1
    partial=interval(0,RHI-RLO)
    co,si=iv.cos(20*r),iv.sin(20*r)
    left=[left[j]+partial*r**4*Cp*trig for j,trig in enumerate([co,si])]
    right=[right[j]+partial*r**-4*Cp*trig for j,trig in enumerate([co,si])]
    # e=1 and e'=0 at the isolated maximum. This exact derivative of the
    # Green formula cancels every C' boundary term.
    Pr=40*C*si+12*r**-5*left[0]-20*r**3*right[0]
    Pi=-40*C*co+12*r**-5*left[1]-20*r**3*right[1]
    m,k=4,20;c=interval(F(7,5));nu=interval(F(1,1000))
    G=gjets(r,4)
    L2=c*c*r*r*G[2]-2*c*nu*r*G[3]+nu*nu*(G[4]-2*G[3]/r+3*G[2]/r**2)
    Kr=k*k*r+m*m/r
    # K=(k²r+m²/r)+ik. Its phase-rotated product with P' is exact.
    product=(Kr*co+k*si)*Pr-(k*co-Kr*si)*Pi
    error_radius=iv.sqrt(Kr*Kr+k*k)*5
    # At a true maximum V0'=-C, so real[-ikr exp(-ikr)s]=-2k²C.
    seed_coef=m*m*G[2]/(2*r*r)+product/2-k*k*C
    seed_upper=seed_coef+error_radius/2
    visc_factor=(2*m*m+4)/(r*r)-2*k*k
    assert L2.b<0 and seed_upper.b<0 and visc_factor.b<0
    Gtt_upper=780*(L2+interval(F(9,16))*seed_upper)
    amplitude=interval(780,1020);lambda2=interval(F(9,16),1)
    T0=lambda2*m*k/(2*r)
    Gt=T0+nu*amplitude*G[2]
    Gtr=amplitude*(-c*r*G[2]+nu*(G[3]-G[2]/r))-T0/r
    # Exact implicit-function correction for the radial maximum at fixed z=4.
    curvature_denominator=-amplitude*G[2]
    moving_correction=Gtr**2/curvature_denominator
    moving_Gtt_upper=Gtt_upper+moving_correction
    radial_speed=c*r-nu*(G[3]/G[2]-1/r)+T0/(r*amplitude*G[2])
    passed=bool(Gtt_upper.b<-250000 and moving_Gtt_upper.b<-240000
                and Gt.a>85 and Gt.b<186 and radial_speed.a>interval(F(1,5)).b)
    output=dict(status='PASS' if passed else 'INCONCLUSIVE',panels=panels,dps=iv.dps,
        unique_global_radial_maximum_interval=[str(RLO),str(RHI)],
        G_jets=[pack(v) for v in G],C=pack(C),
        Green_pressure_radial_real=pack(Pr),Green_pressure_radial_imag=pack(Pi),
        complete_pressure_radial_error_norm_bound='5',
        mean_transport_diffusion_coefficient=pack(L2),
        seed_coefficient_before_error=pack(seed_coef),
        seed_coefficient_upper_expression=pack(seed_upper),
        viscous_torque_factor=pack(visc_factor),
        Gtt_upper_expression=pack(Gtt_upper),Gtt_upper_strictly_below_minus_250000=passed,
        Gt_interval=pack(Gt),Gtr_interval=pack(Gtr),
        moving_radial_maximum_correction=pack(moving_correction),
        moving_radial_maximum_Gtt_upper_expression=pack(moving_Gtt_upper),
        initial_radial_maximum_speed=pack(radial_speed),
        moving_radial_maximum_Gtt_upper_below_minus_240000=bool(moving_Gtt_upper.b<-240000),
        initial_Gt_between_85_and_186=bool(Gt.a>85 and Gt.b<186),
        radial_speed_above_one_fifth=bool(radial_speed.a>interval(F(1,5)).b),
        seconds=time.time()-start,
        scope='Actual initial fixed-circle and fixed-axial-plane radial critical-branch derivatives; uses independently proved whole-space slab error. Not acceleration of an unrestricted global maximum or a bound on positive-time duration.')
    return output

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--panels',type=int,default=2048);ap.add_argument('--output');args=ap.parse_args()
    data=run(args.panels)
    if args.output:Path(args.output).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2),flush=True)

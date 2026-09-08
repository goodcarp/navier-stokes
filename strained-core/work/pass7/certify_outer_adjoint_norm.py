"""Outward-rounded certificate for ||r div(Qw)||_2 of the actual outer seed."""
import sys,json,time,argparse
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here.parent/'pass5'))
from interval_local_integrals import iv,interval,cutoff_jets,range_quadrature

def decode(e):
    v=F(int(e['mantissa']))*F(2)**e['exponent']
    return -v if e['sign'] else v

def certify(panels=24000,relative_width=.1,dps=30):
    iv.dps=dps
    dr=interval(F(7,50)); hz=interval(F(3,5))
    r0=interval(F(21,100)); z0=interval(4)
    @lru_cache(None)
    def cached_jets(key):
        X=iv.mpf(0); X._mpi_=key
        return cutoff_jets(X,1)
    def integrand(x,y):
        e,e1=cached_jets((x*x)._mpi_)
        q,q1=cached_jets((y*y)._mpi_)
        er=2*x*e1/dr; qz=2*y*q1/hz
        r=r0+dr*x; z=z0+hz*y
        R2=r*r+z*z
        denominator=4*iv.pi*R2**(iv.mpf(9)/2)
        Qdiff=15*r*r*(r*r-6*z*z)/denominator
        Qrz=-15*r*z*(4*z*z-3*r*r)/denominator
        imaginary=Qdiff*(er-e/r)*q+Qrz*e*qz
        real=20*Qdiff*e*q
        # Angular integral pi, doubled for the reflected axial packet.
        return 2*iv.pi*dr*hz*r*16*(real**2+imaginary**2)
    part=[F(-1),F(-1,2),F(0),F(1,2),F(1)]
    t=time.time()
    result=range_quadrature(integrand,[part,part],max_panels=panels,
                            relative_width=relative_width)
    upper=decode(result['exact_dyadic_bounds']['upper'])
    lower=decode(result['exact_dyadic_bounds']['lower'])
    target=F(9,2500) # .0036
    result.update(status='PASS' if upper<target**2 else 'VALID_ENCLOSURE_TARGET_NOT_MET',
                  quantity='nd_squared = ||r div(Qw)||_L2(R3)^2',
                  squared_target=str(target**2),norm_target=str(target),
                  target_proved=upper<target**2,
                  norm_display=[float(max(lower,F(0)))**.5,float(upper)**.5],
                  elapsed_seconds=time.time()-t,dps=dps,
                  trust_base='mpmath.iv outward-rounded range quadrature; no numerical convergence assumption',
                  formula='2*pi*int_upper r*m^2*((k*Qdiff*e*q)^2+(Qdiff*(er-e/r)*q+Qrz*e*qz)^2) dr dz',
                  parameters=dict(m=4,k=20,r0='21/100',dr='7/50',z0=4,hz='3/5'))
    weaker=F(773,200000)
    if upper<weaker**2:
        result['weaker_rational_bound_proved']=str(weaker)
        result['weaker_rational_squared_target']=str(weaker**2)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--panels',type=int,default=24000)
    ap.add_argument('--relative-width',type=float,default=.1)
    ap.add_argument('--dps',type=int,default=30)
    ap.add_argument('--output',default=str(here/'outer-adjoint-norm-certificate.json'))
    args=ap.parse_args(); output=args.output; del args.output
    result=certify(**vars(args))
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

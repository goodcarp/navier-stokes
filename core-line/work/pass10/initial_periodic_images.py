"""Axial periodic-image correction to the full initial core p_zz.

Independent source and stress kernels on the compact analytic datum. No
periodic Poisson solve. Image tails have an analytic bound; the quadrature
itself is floating point and NOT an interval enclosure.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,json,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from initial_full_field import fields,_quadratic_composition,_axial_cutoff
from initial_whole_space_integrals import quadrature


def analytic_tail_bound(period,count):
    """Exact upper bound, using the proved E<270000 and pi>3."""
    L=F(period);n=F(count)
    if n*L<=6:raise ValueError('Require period*count>6.')
    return F(540000)/(L*(L*n-6)**4)


def point_cloud(order=72,A=1020.,lam=.25,c=1.4):
    pieces=[]
    R,wR=quadrature(order,[0.,.5,1.,1.25,2.5,5.,6.])
    mu,wmu=np.polynomial.legendre.leggauss(24)
    r=R[:,None]*np.sqrt(1-mu[None,:]**2);z=R[:,None]*mu[None,:]
    weight=2*np.pi*wR[:,None]*R[:,None]**2*wmu[None,:]
    M=fields(r,z,A=0.,lam=0.,c=c)
    U=M['u'][0].real
    pieces.append(dict(r=r.ravel(),z=z.ravel(),weight=weight.ravel(),
        g=M['source'][0].ravel(),rr=(U[0]**2).ravel(),tt=(U[1]**2).ravel(),
        zz=(U[2]**2).ravel(),rz=(U[0]*U[2]).ravel()))
    r,wr=quadrature(order,[0.,.04,.09,.1575,.19,.24,.315])
    z,wz=quadrature(order,[3.4,3.55,3.7,3.775,4.,4.225,4.3,4.45,4.6])
    for start in range(0,len(z),32):
        Rho,Z=np.broadcast_arrays(r[:,None],z[None,start:start+32])
        weight=4*np.pi*Rho*wr[:,None]*wz[None,start:start+32]
        full=fields(Rho,Z,A=A,lam=lam,c=c)
        mean_M=fields(Rho,Z,A=0.,lam=0.,c=c)
        f=_quadratic_composition(Rho,0.,63/200,0)[0]
        qv=_axial_cutoff(Z,9/20)[0]
        V=Rho*f*qv;U4=full['u'][4]
        pieces.append(dict(r=Rho.ravel(),z=Z.ravel(),weight=weight.ravel(),
            g=(full['source'][0]-mean_M['source'][0]).ravel(),
            rr=(2*abs(U4[0])**2).ravel(),
            tt=(A*A*V*V+2*abs(U4[1])**2).ravel(),
            zz=np.zeros(Rho.size),rz=np.zeros(Rho.size)))
    return {key:np.concatenate([p[key] for p in pieces]) for key in pieces[0]}


def kernel_contribution(cloud,distance):
    r,z,w=cloud['r'],cloud['z'],cloud['weight']
    source_kernel=np.zeros_like(r)
    stress_value=np.zeros_like(r)
    for sign in (-1,1):
        Z=z+sign*distance;Q=r*r+Z*Z
        source_kernel+=(2*Z*Z-r*r)/(4*np.pi*Q**2.5)
        Qrr=(12*Q*Q-105*r*r*Z*Z)/(4*np.pi*Q**4.5)
        Qtt=3*(4*Z*Z-r*r)/(4*np.pi*Q**3.5)
        Qzz=-(105*Z**4-90*Q*Z*Z+9*Q*Q)/(4*np.pi*Q**4.5)
        Qrz=-15*r*Z*(4*Z*Z-3*r*r)/(4*np.pi*Q**4.5)
        stress_value-=Qrr*cloud['rr']+Qtt*cloud['tt']+Qzz*cloud['zz']+2*Qrz*cloud['rz']
    # Exact compact-divergence identity integral g=0. The two-image value at
    # source origin is 1/(pi*d^3); removing it retains the n^-5 cancellation.
    source_kernel-=1/(np.pi*distance**3)
    return float(np.sum(w*source_kernel*cloud['g'])),float(np.sum(w*stress_value))


def run(order=72,periods=(16,32,64),counts=(8,16,32)):
    started=time.monotonic();cloud=point_cloud(order)
    monopole=float(np.sum(cloud['weight']*cloud['g']))
    E=float(.5*np.sum(cloud['weight']*(cloud['rr']+cloud['tt']+cloud['zz'])))
    moment=float(np.sum(cloud['weight']*(2*cloud['zz']-cloud['rr']-cloud['tt'])))
    result=[]
    whole_pzz=-8.709035202638178
    for L in periods:
        by_pair=[kernel_contribution(cloud,n*L) for n in range(1,max(counts)+1)]
        for count in counts:
            source=sum(v[0] for v in by_pair[:count]);stress=sum(v[1] for v in by_pair[:count])
            tail=analytic_tail_bound(L,count)
            result.append(dict(period=L,positive_image_pairs=count,
                correction_source=source,correction_stress=stress,
                source_stress_difference=source-stress,
                periodic_pzz_reference=whole_pzz+stress,
                rigorous_image_tail_bound_exact=str(tail),
                rigorous_image_tail_bound_float=float(tail),
                numerical_energy_tail_bound=6*E/(np.pi*L*(L*count-6)**4)))
    return dict(status='NUMERICAL IMAGE QUADRATURE; EXACT ANALYTIC IMAGE-TAIL BOUND, NO QUADRATURE ENCLOSURE',
        order_per_piece=order,initial_parameters=dict(A=1020,lam=.25,c=1.4),
        whole_space_pzz_reference=whole_pzz,energy_from_stress_cloud=E,
        raw_source_monopole_quadrature=monopole,
        axial_stress_moment=moment,
        leading_two_image_coefficient_over_distance_fifth=6*moment/np.pi,
        point_count=len(cloud['r']),results=result,seconds=time.monotonic()-started,
        note='All source sums subtract the exact zero monopole term; all stress sums decay directly as n^-5. Bounds cover omitted exact image integrals, not quadrature error in retained images.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--order',type=int,default=72)
    ap.add_argument('--output',default=str(Path(__file__).with_name('initial-periodic-images.json')))
    args=ap.parse_args();out=run(args.order)
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)

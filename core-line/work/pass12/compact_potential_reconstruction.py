#!/usr/bin/env python3
"""Axis-regular compact potential fit of an inherited MAC increment.

This is an H7, C6 whole-space APPROXIMATION, not a validated NS solution.
The exact original datum is kept separately; only U(T)-U_discrete(0) is fit.
The velocity induced by the potentials, including axial collar terms, is
divergence free analytically. No sampling/fit/residual bound is certified.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
from scipy.fft import fft, ifft
from scipy.linalg import cholesky_banded, cho_solve_banded

HERE = Path(__file__).resolve().parent
PRIOR = next((p for p in (HERE.parent / 'pass11',
                          HERE.parent.parent / 'pass11' / 'experiments',
                          HERE.parent / 'pass11' / 'experiments')
              if (p / 'initial_full_field.py').is_file()), None)
if PRIOR is None:
    raise ImportError('The accompanying pass11 initial-field and MAC source package is required.')
sys.path.insert(0, str(PRIOR))
from initial_full_field import cutoff, fields
from fast_full_mac import FastMAC


class RadialBasis:
    """r^m B(r²), degree nine, with a C8 compact outer zero extension."""
    def __init__(self, R=8., intervals=96, mapping=4., m=0):
        if m < 0 or int(m) != m:
            raise ValueError('Positive modes only; negative modes are conjugates.')
        if intervals < 10:
            raise ValueError('At least ten radial intervals are required.')
        self.R, self.m, self.degree = R, int(m), 9
        breaks = (R*np.sinh(mapping*np.linspace(0, 1, intervals+1))/np.sinh(mapping))**2
        self.knots = np.r_[np.repeat(0.,10),breaks[1:-1],np.repeat(R*R,10)]
        full_n = len(self.knots)-10
        # Remove the nine coefficients whose endpoint vanishing order is <9.
        self.n = full_n-9
        self.spline = BSpline(self.knots, np.eye(full_n)[:,:self.n], 9,
                             extrapolate=False)
        self.scale = np.maximum(np.sqrt((self.knots[:self.n]+self.knots[10:10+self.n])/2),1.e-12)**m

    def matrices(self,r):
        r=np.asarray(r,float)
        if r.ndim != 1 or np.any(r<0):
            raise ValueError('r must be a one-dimensional nonnegative array.')
        inside=r<=self.R
        s=np.minimum(r,self.R)**2
        b=self.spline(s);bs=self.spline(s,nu=1);bss=self.spline(s,nu=2)
        n=self.m
        power=r[:,None]**n/self.scale
        lower=(n*r[:,None]**(n-1)/self.scale) if n else np.zeros_like(power)
        B=power*b
        Br=lower*b+2*r[:,None]*power*bs
        Lap=4*power*((n+1)*bs+r[:,None]**2*bss)
        # m B/r evaluated without a removable division at the axis.
        mBoverR=lower*b
        return tuple(np.where(inside[:,None],a,0.) for a in (B,Br,Lap,mBoverR))


def axial_window(z, inner=6., outer=8.):
    """C-infinity q=1 on |z|<=inner and q=0 on |z|>=outer."""
    if not 0<inner<outer:
        raise ValueError('Require 0 < inner < outer.')
    z=np.asarray(z,float)
    a=.75/(outer*outer-inner*inner)
    arg=.25+a*(z*z-inner*inner)
    return cutoff(arg),2*a*z*cutoff(arg,1)


def lower_banded(A, bandwidth=19):
    out=np.zeros((bandwidth+1,len(A)))
    for d in range(bandwidth+1):
        out[d,:len(A)-d]=np.diag(A,k=-d)
    return out


def fit_mode(g,delta_hat,j,intervals=96):
    """Weighted coupled P/T least squares, retaining every axial phase."""
    m=int(g.modes[j]);basis=RadialBasis(g.R,intervals,g.mapping,m)
    B,Dr,Lap,mB=basis.matrices(g.r)
    Bf,Df,Lf,mBf=basis.matrices(g.rf[1:-1])
    n=basis.n;nv=2*n
    gram0=np.zeros((nv,nv));gram2=np.zeros_like(gram0);cross=np.zeros_like(gram0)
    gram0[::2,::2]=Lap.T@(g.V[:,None]*Lap)
    gram0[1::2,1::2]=mBf.T@(g.W[:,None]*mBf)+Dr.T@(g.V[:,None]*Dr)
    gram2[::2,::2]=Df.T@(g.W[:,None]*Df)+mB.T@(g.V[:,None]*mB)
    H=Df.T@(g.W[:,None]*mBf)+mB.T@(g.V[:,None]*Dr)
    cross[::2,1::2]=H;cross[1::2,::2]=H.T
    b0,b2,b1=map(lower_banded,(gram0,gram2,cross))
    Ur=delta_hat[g.slr];Ut=delta_hat[g.slt];Uz=delta_hat[g.slz]
    rhs=np.empty((nv,g.nz),complex)
    rhs[::2]=-1j*g.kz*(Df.T@(g.W[:,None]*Ur))-g.kz*(mB.T@(g.V[:,None]*Ut))-Lap.T@(g.V[:,None]*Uz)
    rhs[1::2]=-1j*(mBf.T@(g.W[:,None]*Ur))-Dr.T@(g.V[:,None]*Ut)
    coeff=np.empty_like(rhs);max_backward=0.
    for col,k in enumerate(g.kz):
        band=b0+k*b1+k*k*b2
        scale=1/np.sqrt(band[0])
        if not np.isfinite(scale).all():
            raise RuntimeError('Degenerate potential basis; no ridge is inserted silently.')
        eq=band.copy()
        for d in range(len(eq)):
            eq[d,:nv-d]*=scale[d:]*scale[:nv-d]
        factor=cholesky_banded(eq,lower=True,check_finite=False)
        coeff[:,col]=scale*cho_solve_banded((factor,True),scale*rhs[:,col],check_finite=False)
        # Banded solve residual, diagnostic only; this is NOT a continuum error.
        if col in (0,1,g.nz//2,g.nz-1):
            x=coeff[:,col];res=band[0]*x
            for d in range(1,len(band)):
                res[d:]+=band[d,:nv-d]*x[:nv-d]
                res[:nv-d]+=band[d,:nv-d]*x[d:]
            denom=np.linalg.norm(rhs[:,col])+1.e-300
            max_backward=max(max_backward,float(np.linalg.norm(res-rhs[:,col])/denom))
    P,T=coeff[::2],coeff[1::2]
    fitted=np.concatenate((1j*g.kz*(Df@P)+1j*mBf@T,
                           -g.kz*(mB@P)-Dr@T,-Lap@P))
    return basis,P,T,fitted,max_backward


def evaluate_increment(basis,P,T,k,z_origin,L,r,z,window=True):
    """Evaluate a modal increment on product r,z arrays, including joins."""
    r=np.asarray(r,float);z=np.asarray(z,float);k=np.asarray(k,float)
    B,Dr,Lap,mB=basis.matrices(r)
    E=np.exp(1j*k[:,None]*(z[None,:]-z_origin))
    p=P@E;t=T@E;pz=(P*(1j*k))@E
    ur=Dr@pz+1j*mB@t
    ut=1j*mB@pz-Dr@t
    uz=-Lap@p
    if window:
        q,qp=axial_window(z,outer=L/2)
        ur=q*ur+qp*(Dr@p)
        ut=q*ut+1j*qp*(mB@p)
        uz=q*uz
    return np.stack((ur,ut,uz))


def run(args):
    started=time.perf_counter();path=Path(args.checkpoint)
    with np.load(path) as data:
        U=data['U'];p=json.loads(str(data['parameters']))
        oldr=data['r'];oldrf=data['rf'];oldz=data['z']
    g=FastMAC(**{k:p[k] for k in ('nr','nz','R','L','mapping','J','nu')},dealias='padding')
    g.mapping=p['mapping']
    if not(np.array_equal(oldr,g.r) and np.array_equal(oldrf,g.rf) and np.array_equal(oldz,g.z)):
        raise ValueError('Checkpoint geometry does not match its parameters.')
    U0,initial=g.initial();delta=U-U0
    dh=fft(delta,axis=-1,norm='forward');fit=np.empty_like(dh)
    saved={'k':g.kz,'z_origin':g.z[0],'parameters':json.dumps(p),'intervals':args.intervals}
    results=[];collar_norm_sq=0.;offgrid_data=[]
    for j,m in enumerate(g.modes):
        basis,P,T,fit[j],backward=fit_mode(g,dh[j],j,args.intervals)
        saved[f'P_{m}']=P;saved[f'T_{m}']=T
        # Windowed MAC comparison uses the same staggered component locations.
        pc=ifft(P,axis=-1,norm='forward');tc=ifft(T,axis=-1,norm='forward')
        B,Dr,Lap,mB=basis.matrices(g.r);Bf,Df,_,mBf=basis.matrices(g.rf[1:-1])
        q,qp=axial_window(g.z,outer=g.L/2)
        raw=ifft(fit[j],axis=-1,norm='forward')
        joined=q*raw+np.concatenate((qp*(Df@pc),1j*qp*(mB@pc),np.zeros((g.nr,g.nz))))
        err=joined-delta[j];fac=1 if j==0 else 2
        norm=lambda x:float(2*np.pi*g.L/g.nz*fac*np.sum(g.mass[:,None]*abs(x)**2))
        collar_norm_sq+=norm(joined-raw)
        results.append({'mode':int(m),'delta_L2':np.sqrt(norm(delta[j])),
                        'fit_error_L2_unwindowed':np.sqrt(norm(raw-delta[j])),
                        'fit_error_L2_windowed':np.sqrt(norm(err)),
                        'axial_collar_change_L2':np.sqrt(norm(joined-raw)),
                        'max_sampled_normal_equation_relative_residual':backward})
        print(json.dumps(results[-1]),flush=True)
    total_err=np.sqrt(sum(a['fit_error_L2_windowed']**2 for a in results))
    target=Path(args.output);npz=target.with_suffix('.npz')
    np.savez_compressed(npz,**saved)
    report={'status':'COMPLETED H7 WHOLE-SPACE ENDPOINT RECONSTRUCTION; NOT NS VALIDATION',
       'input_checkpoint':str(path),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'parameters':p,'radial_intervals':args.intervals,'radial_degree':9,
       'regularity':'Compact C8/H9 potentials induce C6/H7 velocity; not C-infinity.',
       'base':'Exact unchanged analytic initial datum u0; no fit of u0 is substituted.',
       'increment':'Fit of U(T)-U_discrete(0); all stored complex axial/angular phases retained.',
       'potential_output':str(npz),'potential_sha256':hashlib.sha256(npz.read_bytes()).hexdigest(),
       'modal_results':results,'total_sampled_windowed_fit_error_L2':total_err,
       'axial_collar_change_L2':np.sqrt(collar_norm_sq),
       'sampled_discrete_increment_L2':np.sqrt(g.inner(delta,delta)),
       'elapsed_seconds':time.perf_counter()-started,
       'uncertified':['sampled norm versus continuous norm','radial/axial derivative errors',
          'time-slab evolution and time derivatives','pressure and full projected residual',
          'difference between whole-space NS and finite-cylinder data','endpoint gain and successor class']}
    target.write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--checkpoint',required=True)
    ap.add_argument('--intervals',type=int,default=96);ap.add_argument('--output',required=True)
    a=ap.parse_args();print(json.dumps(run(a),indent=2))

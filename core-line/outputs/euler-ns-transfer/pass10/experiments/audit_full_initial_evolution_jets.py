#!/usr/bin/env python3
"""Read-only return-hook audit of full_initial_evolution.run.

The original solver is not edited. This records separate initial jet terms
and their incompatibilities; all numbers remain exploratory approximations.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import contextlib
import io
import json
import sys
from pathlib import Path
import numpy as np
import full_initial_evolution as evolution


def audit(loc):
    grid=loc['grid'];nu=loc['nu'];T=loc['T']
    U,N0,N1,F,C1=map(loc.get,('U','N0','N1','F','C1'))
    source=loc['source'];nz=grid.nz
    def zjet(a,order=1):return grid.axis(a,0,zorder=order).real
    def omega(a):return float(np.real(grid.axis_weights@(a[:4,nz//2]/grid.r[:4])))
    def local_zjet(a):
        values=grid.axis_weights@a[:4,nz//2-3:nz//2+4]
        return float(np.real(np.dot([-1,9,-45,0,45,-9,1],values)/(60*grid.L/nz)))
    b=[.5*zjet(dd[0][2]) for dd in (U,N0,N1)]
    o=[omega(dd[0][1]) for dd in (U,N0,N1)]
    beta1=b[1]/o[0]-b[0]*o[1]/o[0]**2
    beta2=b[2]/o[0]-b[0]*o[2]/o[0]**2-2*b[1]*o[1]/o[0]**2+2*b[0]*o[1]**2/o[0]**3
    pzz0=loc['pzz0'];pzz1=loc['pzz1']
    beta1pressure=-4-pzz0/2
    beta2pressure=-(pzz1+32)/2-10*beta1pressure
    lapF=grid.vector_lap(F[0],0)
    components={
        'viscous_DeltaF':.5*nu*zjet(lapF[2]),
        'viscous_grad_source':.5*nu*zjet(source[0],2),
        'differentiated_convection':-.5*zjet(C1[0][2]),
        'pressure':-.5*pzz1,
    }
    expected={
        'b1_from_pressure':2+beta1pressure,
        'b2_from_pressure':-4*(2+beta1pressure)-pzz1/2,
        'Omega1_from_affine_core':2.,
        'Omega2_from_pressure':2*(2+beta1pressure)+4,
        'convection_b2_from_direct_b1':-4*b[1],
        'convection_b2_from_pressure':-4*(2+beta1pressure),
    }
    actualT=(b[0]+T*b[1]+T*T*b[2]/2)/(o[0]+T*o[1]+T*T*o[2]/2)
    reported=loc['result']['quadratic_field_at_T']['core_beta']
    consistency=abs(actualT-reported)
    assert consistency<2.e-9,(consistency,actualT,reported)
    # Linearity of the Fourier/axis extraction, independent of whether its
    # discrete derivative approximates the physical derivative accurately.
    assert abs(sum(components.values())-b[2])<2.e-6*max(1,abs(b[2]))
    return {
        'status':'EXPLORATORY DIRECT-JET AUDIT',
        'parameters':loc['result']['parameters'],
        'b_jets_0_1_2':b,'Omega_jets_0_1_2':o,
        'beta1_direct':beta1,'beta2_direct':beta2,
        'beta1_pressure_formula':beta1pressure,'beta2_pressure_formula':beta2pressure,
        'b2_contributions':components,'continuous_affine_targets':expected,
        'C1_z_derivative_spectral':zjet(C1[0][2]),
        'C1_z_derivative_local_centered_6th':local_zjet(C1[0][2]),
        'N0_z_derivative_spectral':zjet(N0[0][2]),
        'N0_z_derivative_local_centered_6th':local_zjet(N0[0][2]),
        'axis_divergence_N0':float(np.real(zjet(loc['divN0'][0],0))),
        'quadratic_quotient_from_direct_jets':actualT,
        'quadratic_quotient_reported':reported,
        'quotient_reconstruction_absolute_error':consistency,
        'original_summary':loc['result'],
        'scope':'Read-only extraction from the original solver. Direct means '
                'the derivative of its reconstructed grid polynomial, not '
                'a certified derivative of the whole-space NS solution. '
                'Local finite differences are an independent diagnostic, '
                'not replacement values or error bounds.'}


def run(nr=384,nz=1024,pressure='fd'):
    captured=[]
    target=evolution.run.__code__
    def profile(frame,event,arg):
        if event=='return' and frame.f_code is target and arg is not None:
            captured.append(audit(frame.f_locals))
    previous=sys.getprofile()
    try:
        sys.setprofile(profile)
        with contextlib.redirect_stdout(io.StringIO()):
            evolution.run(nr=nr,nz=nz,pressure=pressure,second=True)
    finally:
        sys.setprofile(previous)
    assert len(captured)==1
    return captured[0]


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--nr',type=int,default=384)
    parser.add_argument('--nz',type=int,default=1024)
    parser.add_argument('--pressure',choices=['fv','fd'],default='fd')
    parser.add_argument('--output')
    args=parser.parse_args()
    data=run(args.nr,args.nz,args.pressure)
    text=json.dumps(data,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text)
    print(text,end='')

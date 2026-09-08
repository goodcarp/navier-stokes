#!/usr/bin/env python3
"""Checkpointed finite-cylinder experiment; no continuum error certificate."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.fft import fft,ifft,set_workers
from fast_full_mac import FastMAC


def mode_energies(g,U):
    return [float(np.pi*g.L/g.nz*(1 if j==0 else 2)*
        np.sum(g.mass[:,None]*abs(a)**2)) for j,a in enumerate(U)]


def run(args):
    started=time.perf_counter();target=Path(args.output)
    parameters={k:getattr(args,k) for k in ('nr','nz','R','L','mapping','J','nu')}
    g=FastMAC(**parameters,dealias='padding');U,init=g.initial()
    a0=g.observe(U,0);a0['angular_mode_energies']=mode_energies(g,U)
    ah=fft(U,axis=-1);vh=np.empty_like(ah)
    for j in range(g.J+1):
        vh[j]=-g.nu*((g.K[j]@ah[j])/g.mass[:,None]+g.kz**2*ah[j])
    rhs=g.nonlinear(U)+g.project(ifft(vh,axis=-1))
    rhobs=g.observe(rhs,0);fluct=U.copy();fluct[0]=0
    init['projected_initial_derivatives']={
        'Omega_t':rhobs['core_Omega'],'b_t':rhobs['core_b'],
        'beta_t':(rhobs['core_b']*a0['core_Omega']-a0['core_b']*rhobs['core_Omega'])/a0['core_Omega']**2,
        'mean_G_fixed_t':rhobs['mean_G_fixed'],'K_t':g.inner(fluct,rhs),
        'energy_t':g.inner(U,rhs)}
    del ah,vh,rhs,fluct
    steps=max(1,int(np.ceil(args.T/args.dt)));dt=args.T/steps
    result={'status':'RUNNING FINITE-CYLINDER DIAGNOSTIC; NOT WHOLE-SPACE NS VALIDATION',
        'complete':False,'parameters':dict(parameters,T=args.T,dt=dt,steps=steps,dealias='padding',backend='fast',fft_workers=args.fft_workers),
        'source_sha256':{name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
            for name in ('evolve_full_mac.py','fast_full_mac.py','initial_full_field.py','run_resolved_mac.py')},
        'initial_reconstruction':init,'history':[a0],
        'remaining_errors':['finite radial no-slip wall and periodic axial images',
            'initial sampling and weighted projection','continuous divergence and smooth axis reconstruction',
            'angular/radial/axial truncation','time discretization','whole-space tails and space-time residual','endpoint successor class']}
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'initial':a0,'initial_derivatives':init['projected_initial_derivatives']}),flush=True)
    viscous_loss=0.
    old_D=a0['discrete_dissipation']
    for i in range(steps):
        U=g.step(U,dt)
        if not np.isfinite(U).all():raise RuntimeError('Nonfinite trajectory; no endpoint accepted')
        obs=g.observe(U,(i+1)*dt);new_D=obs['discrete_dissipation']
        viscous_loss+=g.nu*dt*(old_D+new_D)/2;old_D=new_D
        obs['energy_balance_defect']=obs['energy']-a0['energy']+viscous_loss
        obs['angular_mode_energies']=mode_energies(g,U)
        result['history'].append(obs);result['seconds']=time.perf_counter()-started
        if i+1==steps:
            result['complete']=True;result['status']='COMPLETED FINITE-CYLINDER DIAGNOSTIC; NOT WHOLE-SPACE NS VALIDATION'
            if args.checkpoint:np.savez_compressed(args.checkpoint,U=U,r=g.r,rf=g.rf,z=g.z,parameters=json.dumps(result['parameters']))
        target.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'step':i+1,'steps':steps,'seconds':result['seconds'],**obs}),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for name,typ,default in [('nr',int,192),('nz',int,2049),('R',float,8.),('L',float,16.),('mapping',float,4.),('J',int,4),('nu',float,.001),('T',float,.001),('dt',float,.0001),('fft-workers',int,4)]:
        ap.add_argument('--'+name,type=typ,default=default)
    ap.add_argument('--output',required=True);ap.add_argument('--checkpoint')
    args=ap.parse_args()
    with set_workers(args.fft_workers):run(args)

#!/usr/bin/env python3
"""Initial representation and generator calibration only, not time evolution."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,json
from pathlib import Path
import numpy as np
from scipy.fft import fft,ifft,set_workers
from fast_full_mac import FastMAC


def probe(nr,nz,R=8.,L=16.,mapping=4.):
    # J=1 retains the unchanged initial field and the entire mean generator.
    # It is not an acceptable nonlinear evolution truncation by itself.
    g=FastMAC(nr=nr,nz=nz,R=R,L=L,mapping=mapping,J=1)
    U,initial=g.initial();obs=g.observe(U,0)
    ah=fft(U,axis=-1);vh=np.empty_like(ah)
    for j in range(g.J+1):vh[j]=-g.nu*((g.K[j]@ah[j])/g.mass[:,None]+g.kz**2*ah[j])
    viscous=g.project(ifft(vh,axis=-1));rhs=g.nonlinear(U)+viscous
    vobs=g.observe(viscous,0);robs=g.observe(rhs,0)
    fluct=U.copy();fluct[0]=0
    return dict(nr=nr,nz=nz,R=R,L=L,mapping=mapping,**initial,**obs,
        discrete_Omega_prime=robs['core_Omega'],viscous_Omega_prime=vobs['core_Omega'],
        discrete_b_prime=robs['core_b'],
        discrete_beta_prime=(robs['core_b']*obs['core_Omega']-obs['core_b']*robs['core_Omega'])/obs['core_Omega']**2,
        discrete_G_prime=robs['mean_G_fixed'],discrete_K_prime=g.inner(fluct,rhs))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--nr',type=int,default=192)
    ap.add_argument('--nz',type=int,nargs='+',default=[513,1025,2049])
    ap.add_argument('--R',type=float,default=8.);ap.add_argument('--L',type=float,default=16.)
    ap.add_argument('--mapping',type=float,default=4.);ap.add_argument('--output',required=True)
    args=ap.parse_args();rows=[]
    with set_workers(4):
        for nz in args.nz:
            row=probe(args.nr,nz,args.R,args.L,args.mapping);rows.append(row);print(json.dumps(row),flush=True)
    Path(args.output).write_text(json.dumps(dict(status='NUMERICAL INITIAL GRID DIAGNOSTIC ONLY',results=rows),indent=2)+'\n')

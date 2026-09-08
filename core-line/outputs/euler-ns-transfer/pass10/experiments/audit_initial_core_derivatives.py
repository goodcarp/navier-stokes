#!/usr/bin/env python3
"""Read-only diagnosis of initial core derivative aliasing; not a certificate.

Only the explicit initial field is evaluated. No pressure solve or evolution
is repeated. In the exact affine core, all reported viscous terms vanish.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import json
import numpy as np
from high_order_pressure import HighOrderCylinder
from initial_full_field import fields


def run(nr=384,nz=1024,R=8.,L=16.,mapping=4.,nu=.001):
    grid=HighOrderCylinder(nr,nz,R,L,mapping)
    F=np.empty((3,nr,nz),complex)
    U=np.empty_like(F)
    source=np.empty((nr,nz),complex)
    for i in range(0,nr,24):
        sl=slice(i,min(i+24,nr))
        f=fields(grid.r[sl,None],grid.z[None,:])
        U[:,sl]=f['u'][0]
        F[:,sl]=nu*f['laplacian'][0]-f['convection'][0]
        source[sl]=f['source'][0]
    lapF=grid.vector_lap(F,0)
    return {'status':'NUMERICAL DIAGNOSTIC ONLY',
        'parameters':dict(nr=nr,nz=nz,R=R,L=L,mapping=mapping,nu=nu),
        'b0_extracted':.5*grid.axis(U[2],0,zorder=1).real,
        'Fz_z_at_core':grid.axis(F[2],0,zorder=1).real,
        'Fz_z_exact':-4.,
        'source_zz_at_core':grid.axis(source,0,zorder=2).real,
        'viscous_source_contribution_to_b2':.5*nu*grid.axis(source,0,zorder=2).real,
        'lapFz_z_at_core':grid.axis(lapF[2],0,zorder=1).real,
        'viscous_lapF_contribution_to_b2':.5*nu*grid.axis(lapF[2],0,zorder=1).real,
        'exact_viscous_contributions':0.,
        'scope':'These are the independently differentiated explicit initial '
                'field terms used by full_initial_evolution.py. They diagnose '
                'spectral differentiation leakage; no pressure or actual '
                'positive-time NS value is certified.'}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--nr',type=int,default=384)
    p.add_argument('--nz',type=int,default=1024)
    p.add_argument('--output')
    a=p.parse_args()
    result=run(a.nr,a.nz)
    text=json.dumps(result,indent=2)+'\n'
    if a.output:
        from pathlib import Path
        Path(a.output).write_text(text)
    print(text,end='')

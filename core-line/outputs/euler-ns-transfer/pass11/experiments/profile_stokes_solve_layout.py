#!/usr/bin/env python3
"""Small-factor-only timing; no initialization or evolution of a large field."""
import json
import time
from pathlib import Path
import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import splu
from evolve_full_mac import MAC


def main():
    rng = np.random.default_rng(163909)
    records = []
    for n in (96, 384):
        # Axial cache is deliberately tiny. Each tested matrix is assembled at
        # its specified physical k and is independent of this cache length.
        g = MAC(nr=n, nz=16, R=8., L=16., mapping=4., J=1, nu=.001)
        j, alpha = 1, (2e-5/2)*g.nu/2
        for ki in (1, 32):
            k = 2*np.pi*ki/g.L
            D = bmat([[g.Divr,diags(1j*g.modes[j]/g.r),
                       1j*k*diags(np.ones(n))]],format='csc')
            BM=-D.conj().T@diags(g.V)
            A=diags(g.mass*(1+alpha*k*k))+alpha*g.K[j]
            full=bmat([[A,BM],[BM.conj().T,None]],format='csc')
            old,right=[],[]
            for i in range(n):
                if i>0:old.append(i-1)
                old.extend((n-1+i,2*n-1+i,g.nv+i))
                right.extend((n-1+i,2*n-1+i))
                if i<n-1:right.append(i)
                right.append(g.nv+i)
            entries={}
            raw=rng.normal(size=(g.nv+n,2))+1j*rng.normal(size=(g.nv+n,2))
            raw[g.nv:]=0
            for name,p in [('left',old),('right',right)]:
                p=np.array(p)
                lu=splu(full[p][:,p],permc_spec='NATURAL')
                b=raw[p].copy(order='F')
                entries[name]=(lu,b)
                lu.solve(b)
            samples={name:[] for name in entries}
            batch_samples={name:[] for name in entries}
            for rep in range(8):
                # Reverse the measured ordering every repetition.
                order=('left','right') if rep%2==0 else ('right','left')
                for name in order:
                    lu,b=entries[name]
                    t0=time.perf_counter()
                    for _ in range(50):
                        lu.solve(b[:,0]);lu.solve(b[:,1])
                    samples[name].append((time.perf_counter()-t0)/50)
                    t0=time.perf_counter()
                    for _ in range(50):
                        lu.solve(b)
                    batch_samples[name].append((time.perf_counter()-t0)/50)
            for name,(lu,b) in entries.items():
                records.append(dict(nr=n,m=4,k_index=ki,ordering=name,
                    dimension=full.shape[0],matrix_nnz=full.nnz,
                    L_nnz=lu.L.nnz,U_nnz=lu.U.nnz,
                    median_two_separate_rhs_seconds=float(np.median(samples[name])),
                    median_two_batched_rhs_seconds=float(np.median(batch_samples[name])),
                    two_separate_rhs_samples=samples[name],
                    two_batched_rhs_samples=batch_samples[name]))
    out=dict(status='PERFORMANCE DIAGNOSTIC ONLY',
             note='No large field initial() or step() was called. Factor-only tests.',
             records=records)
    Path(__file__).with_name('stokes-solve-layout-timing.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if not k.endswith('samples')} for r in records],indent=2))


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Benchmark unchanged nonlinear operator with one and four FFT workers."""
import hashlib
import inspect
import json
import time
from pathlib import Path
from evolve_full_mac import MAC
import numpy as np
from scipy.fft import set_workers


def main():
    source=Path(__file__).with_name('evolve_full_mac.py')
    before=hashlib.sha256(source.read_bytes()).hexdigest()
    method_source=inspect.getsource(MAC.nonlinear)
    g=MAC(nr=96,nz=256,R=8.,L=16.,mapping=4.,J=4,nu=.001,dealias='padding')
    U,initial=g.initial()
    for workers in (1,4):
        with set_workers(workers):
            g.nonlinear(U)
    results={1:[],4:[]}
    outputs={}
    for workers in (1,4,4,1,1,4):
        with set_workers(workers):
            t0=time.perf_counter()
            out=g.nonlinear(U)
            elapsed=time.perf_counter()-t0
        results[workers].append(elapsed)
        outputs[workers]=out
    rel=np.linalg.norm(outputs[1]-outputs[4])/np.linalg.norm(outputs[1])
    assert rel<1e-12
    result=dict(status='PASS: equivalent nonlinear outputs; timing diagnostic only',
        parameters=dict(nr=96,nz=256,J=4,dealias='padding'),
        source_sha256_before=before,
        source_sha256_after=hashlib.sha256(source.read_bytes()).hexdigest(),
        imported_nonlinear_source=method_source,
        samples_seconds=results,
        medians_seconds={k:float(np.median(v)) for k,v in results.items()},
        four_vs_one_speedup=float(np.median(results[1])/np.median(results[4])),
        relative_output_difference=float(rel),
        note='Same fixed initialized state. Alternating workers, warmed FFT calls, '
             'no cProfile. Concurrent load uncontrolled; no PDE error conclusion.')
    Path(__file__).with_name('fft-worker-timing96.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='imported_nonlinear_source'},indent=2))


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Small-case performance profile only; leaves the solver source unchanged."""
import cProfile
import hashlib
import json
import os
from pathlib import Path
import platform
import pstats
import time

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'evolve_full_mac.py'
before_hash=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
from evolve_full_mac import MAC
import numpy as np
import scipy


def stats(profiler,limit=35):
    parsed=pstats.Stats(profiler)
    rows=[]
    for (filename,line,name),(primitive,total,self_time,cumulative,callers) in parsed.stats.items():
        rows.append({'function':name,'file':Path(filename).name,'line':line,
                     'primitive_calls':primitive,'total_calls':total,
                     'self_seconds':self_time,'cumulative_seconds':cumulative})
    return {'total_primitive_calls':parsed.prim_calls,'total_calls':parsed.total_calls,
            'total_profiled_seconds':parsed.total_tt,
            'by_cumulative':sorted(rows,key=lambda d:d['cumulative_seconds'],reverse=True)[:limit],
            'by_self':sorted(rows,key=lambda d:d['self_seconds'],reverse=True)[:limit]}


def run():
    parameters=dict(nr=96,nz=256,R=8.,L=16.,mapping=4.,J=4,nu=.001,dealias='padding')
    dt=2e-5
    start=time.perf_counter();grid=MAC(**parameters)
    constructor_seconds=time.perf_counter()-start
    start=time.perf_counter();U,initial=grid.initial()
    initialization_seconds=time.perf_counter()-start

    # This exactly matches the floating arithmetic used for each dt/2
    # diffusion call. Only factors are warmed; no extra time step is taken.
    alpha=(dt/2)*grid.nu/2
    wave_numbers=sorted(set(abs(int(k)) for k,keep in zip(grid.ki,grid.keep) if keep))
    warm=cProfile.Profile();warm.enable();start=time.perf_counter()
    for j in range(grid.J+1):
        for k in wave_numbers:grid._stokes_factor(j,k,alpha)
    warm_seconds=time.perf_counter()-start;warm.disable()
    cache_before=len(grid.stokes_cache)

    calls={'diffuse':[],'nonlinear':[]}
    for name in calls:
        original=getattr(grid,name)
        def make_wrapper(function,key):
            def wrapped(*args,**kwargs):
                t=time.perf_counter()
                try:return function(*args,**kwargs)
                finally:calls[key].append(time.perf_counter()-t)
            return wrapped
        setattr(grid,name,make_wrapper(original,name))

    steady=cProfile.Profile();steady.enable();start=time.perf_counter()
    endpoint=grid.step(U,dt)
    step_seconds=time.perf_counter()-start;steady.disable()
    require_cache=len(grid.stokes_cache)
    assert require_cache==cache_before,'The timed step unexpectedly created new factors.'
    assert len(calls['diffuse'])==2 and len(calls['nonlinear'])==4
    assert np.isfinite(endpoint).all()
    pools=[]
    try:
        from threadpoolctl import threadpool_info
        pools=threadpool_info()
    except ImportError:
        pass
    return {'status':'PERFORMANCE PROFILE ONLY',
        'parameters':dict(parameters,dt=dt),
        'source_sha256_before_import':before_hash,
        'source_sha256_after_run':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'environment':{'python':platform.python_version(),'numpy':np.__version__,
            'scipy':scipy.__version__,'platform':platform.platform(),
            'OPENBLAS_NUM_THREADS':os.environ.get('OPENBLAS_NUM_THREADS'),
            'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),
            'threadpools':pools},
        'constructor_seconds':constructor_seconds,
        'initialization_seconds':initialization_seconds,
        'unique_stokes_factors':cache_before,
        'factor_warmup_seconds':warm_seconds,
        'steady_one_step_seconds':step_seconds,
        'steady_stage_timings':{name:{'calls':len(times),'seconds':times,
                                    'sum_seconds':sum(times)} for name,times in calls.items()},
        'factor_count_after_step':require_cache,
        'factor_warmup_profile':stats(warm),
        'steady_step_profile':stats(steady),
        'initialization_diagnostics':initial,
        'endpoint_divergence_L2':grid.divergence_norm(endpoint),
        'scope':'One requested 96x256,J4 step after explicit factor warm-up. '
                'No solver source changes and no large-case rerun. cProfile adds '
                'overhead; background process contention was not controlled. '
                'Timing results are not NS validation results.'}


if __name__=='__main__':
    result=run()
    target=HERE/'mac-step-timing96.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','parameters','constructor_seconds',
        'initialization_seconds','unique_stokes_factors','factor_warmup_seconds',
        'steady_one_step_seconds','steady_stage_timings')},indent=2))
    print('Steady top self-time functions:')
    for row in result['steady_step_profile']['by_self'][:12]:
        print(json.dumps(row))
    print('Factor warm-up top self-time functions:')
    for row in result['factor_warmup_profile']['by_self'][:8]:
        print(json.dumps(row))

#!/usr/bin/env python3
"""Independent radial-reduction controls for the exploratory pressure solver.

These check implementation/convergence on finite-harmonic core fields;
they do not certify the outer-pump gate or provide interval error bounds.
"""
import sys
from pathlib import Path
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here.parent/'experiments'))
sys.path.insert(0,str(here))
import numpy as np
from evaluate_pressure_gate import Evaluator
from evaluate_core_pressure_radial import evaluate

reference=evaluate(8001)
errors=[]
for nr in (800,1600):
    ev=Evaluator(nr=nr,nmu=16,lmax=8)
    assert abs(ev.pressures[0]+18/7)<1.e-6
    assert abs(ev.pressures[1]-.4)<1.e-6
    for field in ev.fields:
        G=ev.gradient(field)
        divergence=np.trace(G,axis1=0,axis2=1)
        assert np.max(np.abs(divergence))<1.e-10*max(1.,np.max(np.abs(G)))
    pure=ev.evaluate(b=1.,omega=0.,c=0.,A=0.)
    mixed=ev.evaluate(b=1.,omega=1.,c=0.,A=0.)
    e=max(abs(pure['pzz_prime']-reference['t300']),
          abs(mixed['pzz_prime']-reference['core_b1_Omega1']))
    errors.append(e)
    assert pure['poisson_trace_relative']<1.e-11
    assert mixed['poisson_trace_relative']<1.e-11
    assert abs(mixed['representation_difference'])<1.e-12
    print(f'Core control nr={nr}: t300={pure["pzz_prime"]:.12f}, '
          f't300+tWb={mixed["pzz_prime"]:.12f}, max error={e:.6g}')
assert errors[1]<1.e-4
assert errors[1]<.3*errors[0]
print('PASS: analytic core pressure coefficients, divergence, full-harmonic trace, '
      'and convergence to an independent one-dimensional cubic reduction.')
print('No interval certification of the outer-pump pressure gate.')

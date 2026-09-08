#!/usr/bin/env python3
"""Recompute scoped finite results; no compact 3D evolution is certified."""
import json,os,subprocess,sys
from pathlib import Path
here=Path(__file__).resolve().parent
experiments=here.parent/'experiments' if (here.parent/'experiments').is_dir() else here
programs=[
    ['verify_general_envelope_jet.py'],
    ['certify_general_envelope_initial_jet.py','--panels','2048','--output',str(here/'general-envelope-initial-jet-certificate.json')],
    ['verify_general_envelope_initial_jet_certificate.py'],
    ['verify_evolution_validation_bridge.py'],
    ['verify_affine_envelope_budget.py','--panels','2048'],
    ['verify_fixed_amplitude_initialization.py'],
    ['check_slice_solver.py'],
    ['compare_slice_runs.py','--folder',str(experiments),'--output',str(here/'slice-comparison.json')],
]
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
for args in programs:
    print('RUN '+args[0],flush=True)
    result=subprocess.run([sys.executable,str(here/args[0]),*args[1:]],env=env)
    if result.returncode:raise SystemExit('FAIL '+args[0])
    print('PASS '+args[0],flush=True)
for name in ['general-envelope-initial-jet-certificate.json','affine-envelope-budget-certificate.json','fixed-amplitude-initialization-certificate.json']:
    assert json.loads((here/name).read_text())['status']=='PASS'
print('PASS: 8 finite programs; two interval certificates recomputed, exact fixed-amplitude algebra replayed, and stored nonlinear slice sensitivity checked. Full compact 3D time evolution, useful actual stage duration, inherited return and blowup remain unproved.',flush=True)

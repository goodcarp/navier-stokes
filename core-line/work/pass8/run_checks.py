#!/usr/bin/env python3
"""Recompute three new interval certificates and replay the finite identities."""
import json,subprocess,sys
from pathlib import Path
here=Path(__file__).resolve().parent
commands=[
 ('verify_covariance_evolution.py',[]),
 ('verify_covariance_pressure_slab.py',[]),
 ('certify_fixed_circle_acceleration.py',['--panels','2048','--output',str(here/'fixed-circle-acceleration-certificate.json')]),
 ('verify_fixed_circle_acceleration_audit.py',[]),
 ('verify_actual_fluctuation_energy.py',[]),
 ('kelvin-covariance-check.py',[]),
 ('verify_endcap_gradient_bound.py',[]),
 ('exterior-dynamics-verify.py',[]),
 ('certify_leading_envelope.py',['--panels','1024','--output',str(here/'leading-envelope-certificate.json')]),
 ('certify_leading_pressure.py',['--panels','1024','--output',str(here/'leading-pressure-certificate.json')]),
 ('verify_leading_family_independent_audit.py',[]),
]
for name,args in commands:
    print('RUN '+name,flush=True)
    p=subprocess.run([sys.executable,str(here/name),*args],cwd=here,capture_output=True,text=True)
    print(p.stdout,flush=True)
    if p.stderr:print(p.stderr,flush=True)
    if p.returncode:raise SystemExit(p.returncode)
    if name.startswith('certify_'):
        data=json.loads(Path(args[-1]).read_text())
        if data['status']!='PASS':raise RuntimeError(name+' failed its interval predicates')
    print('PASS '+name,flush=True)
print('PASS: 11 finite programs; all three new interval certificates recomputed. Model diagnostics remain separate from actual initial identities. No NS time integration, formal replay, quantitative stage duration, inherited return, or blowup proof.',flush=True)

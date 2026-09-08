#!/usr/bin/env python3
"""Recompute the two decisive intervals, then replay the finite identities."""
import json,subprocess,sys
from pathlib import Path
here=Path(__file__).resolve().parent
commands=[
 ('verify_D_pressure_kernels.py',[]),
 ('certify_D_viscosity.py',['--panels','1024']),
 ('verify_E_interaction_reduction.py',[]),
 ('certify_outer_E.py',['--panels','2048','--output',str(here/'outer-E-certificate.json')]),
 ('verify_E_certificate_audit.py',[]),
 ('verify_common_family.py',[]),
 ('verify_next_stage_shear.py',[]),
]
for name,args in commands:
    print('RUN '+name,flush=True)
    result=subprocess.run([sys.executable,str(here/name),*args],cwd=here,text=True,capture_output=True)
    if name=='certify_D_viscosity.py' and result.returncode==0:
        data=json.loads(result.stdout)
        if data['status']!='PASS':raise RuntimeError('D interval certificate did not pass')
        (here/'D-viscosity-certificate.json').write_text(result.stdout)
    if name=='certify_outer_E.py' and result.returncode==0:
        data=json.loads((here/'outer-E-certificate.json').read_text())
        if data['status']!='PASS':raise RuntimeError('E interval certificate did not pass')
    print(result.stdout,flush=True)
    if result.stderr:print(result.stderr,flush=True)
    if result.returncode:raise SystemExit(result.returncode)
    print('PASS '+name,flush=True)
print('PASS: 7 finite check programs (six core checks plus one frozen-model diagnostic); both decisive interval calculations recomputed. No PDE time integration, formal proof replay, inherited return, or blowup proof.',flush=True)

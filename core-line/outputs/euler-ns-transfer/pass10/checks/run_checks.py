#!/usr/bin/env python3
"""Run scoped checks on temporary copies, preserving the published archive."""
import json,os,shutil,subprocess,sys,tempfile
from pathlib import Path

root=Path(__file__).resolve().parent.parent
programs=[
    'verify_time_polynomial_residual.py',
    'verify_cylindrical_NS_representation.py',
    'check_initial_full_field.py',
    'verify_full_initial_coefficient_identities.py',
    'check_manufactured_pressure.py',
    'verify_corotating_time_polynomial.py',
    'verify_periodic_image_tail.py',
    'check_full_reference_comparisons.py',
    'verify_sharp_error_hierarchy.py',
]
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',
    PASS10_INITIAL_JET_CERTIFICATE=str(root.parent/'pass9'/'checks'/'general-envelope-initial-jet-certificate.json'))
with tempfile.TemporaryDirectory(prefix='pass10-checks-') as temp:
    folder=Path(temp)/'experiments'
    shutil.copytree(root/'experiments',folder)
    for name in programs:
        print('RUN '+name,flush=True)
        result=subprocess.run([sys.executable,str(folder/name)],env=env,cwd=folder)
        if result.returncode:raise SystemExit('FAIL '+name)
        print('PASS '+name,flush=True)
    manufactured=json.loads((folder/'manufactured-pressure-checks.json').read_text())
    assert all(all(x.values()) for x in manufactured['checks'].values())
    assert json.loads((folder/'periodic-image-tail-certificate.json').read_text())['status']=='PASS'
print(f'PASS: {len(programs)} scoped programs. No full 3D NS endpoint, commuting projection, complete residual bound, inherited return or blowup is certified.',flush=True)

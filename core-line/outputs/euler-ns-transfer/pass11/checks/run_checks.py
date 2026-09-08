#!/usr/bin/env python3
"""Bounded checks on temporary copies; no production trajectory rerun."""
import os,shutil,subprocess,sys,tempfile
from pathlib import Path
root=Path(__file__).resolve().parent.parent
programs=['verify_poloidal_toroidal_design.py','check_shared_hankel_projection.py',
    'check_full_mac.py','check_fast_full_mac.py','verify_validated_error_hierarchy.py',
    'check_independent_hierarchy_audit.py','check_exterior_harmonic_mass.py',
    'check_exterior_mac_blocks.py','check_compare_mac_endpoints.py']
with tempfile.TemporaryDirectory(prefix='pass11-checks-') as tmp:
    folder=Path(tmp)/'experiments';shutil.copytree(root/'experiments',folder)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    for name in programs:
        print('RUN '+name,flush=True)
        r=subprocess.run([sys.executable,str(folder/name)],cwd=folder,env=env)
        if r.returncode:raise SystemExit('FAIL '+name)
        print('PASS '+name,flush=True)
print('PASS: 9 bounded programs. No validated whole-space NS endpoint, successor class, or blowup is certified.',flush=True)

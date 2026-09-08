#!/usr/bin/env python3
"""Check frozen pass manifests and rerun four bounded programs in a temp copy.

The numerical tests are manufactured/reconstruction checks, not an NS error
certificate. Original source and recorded outputs are never executed in place.
"""
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TOP=ROOT/'core-line/outputs/euler-ns-transfer'
SCRIPTS=['check_c4_adapted_reference.py','verify_comoving_gradient_decay.py',
         'check_compact_potential_reconstruction.py','check_continuous_potential_fit.py']


def main():
    manifest_counts={}
    failures=[]
    # Exclude historical top-before-pass12 snapshots: their relative checks/
    # paths describe their original root location, not the provenance folder.
    manifests=[TOP/'SHA256SUMS.txt']
    for n in range(2,13):
        d=TOP/f'pass{n}'
        manifests += [p for p in d.rglob('SHA256SUMS.txt') if 'provenance' not in p.relative_to(d).parts]
    for p in manifests:
        count=0
        for line in p.read_text().splitlines():
            expected,rel=line.split('  ',1)
            target=p.parent/rel
            if not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest()!=expected:
                failures.append(dict(manifest=str(p.relative_to(ROOT)),file=rel))
            count+=1
        manifest_counts[str(p.relative_to(TOP))]=count
    outputs=[]
    with tempfile.TemporaryDirectory(prefix='forced-route-pass12-') as td:
        temp=Path(td)
        for n in (11,12):
            shutil.copytree(TOP/f'pass{n}'/'experiments',temp/f'pass{n}'/'experiments',
                            ignore=shutil.ignore_patterns('__pycache__','*.pyc','*.npz'))
        cwd=temp/'pass12/experiments'
        for name in SCRIPTS:
            proc=subprocess.run([sys.executable,name],cwd=cwd,text=True,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
            outputs.append(dict(script=name,exit_code=proc.returncode,
                                stdout=proc.stdout,stderr=proc.stderr))
            print(name, 'PASS' if proc.returncode==0 else 'FAILED',flush=True)
    passed=not failures and all(o['exit_code']==0 for o in outputs)
    result=dict(status='PASS' if passed else 'FAILED',
        scope='Package/source integrity plus exact finite algebra and manufactured floating '
              'reconstruction tests; no whole-space NS residual, useful stage, S3 or return certificate.',
        manifest_counts=manifest_counts,manifest_failures=failures,programs=outputs)
    (ROOT/'manifests/pass12-checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='programs'},indent=2))
    return 0 if passed else 1


if __name__=='__main__':raise SystemExit(main())

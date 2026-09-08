#!/usr/bin/env python3
"""Package the completed second-pass research notes and reproduce their checks."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

root=Path(__file__).resolve().parents[2]
src=root/'work/pass2'
dst=root/'outputs/euler-ns-transfer/pass2'
dst.mkdir(parents=True,exist_ok=True)
(dst/'checks').mkdir(exist_ok=True)
notes=['REPORT.md','VERIFICATION.md','active-core-pressure.md',
       'local-pressure-independent-audit.md','affine-stretch-gate.md',
       'correction-window.md','correction-window-independent-audit.md',
       'mean-feedback.md','jeong-yoneda-stage.md','source-check.md']
scripts=['verify_active_core_pressure.py','verify_local_pressure_audit.py',
         'verify_affine_stretch_gate.py','verify_correction_window.py']
for name in notes:
    shutil.copyfile(src/name,dst/name)
for name in scripts:
    shutil.copyfile(src/name,dst/'checks'/name)
runner='''#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys
here=Path(__file__).resolve().parent
for name in '''+repr(scripts)+''':
    print("CHECK:",name,flush=True)
    subprocess.run([sys.executable,str(here/name)],check=True)
print("All packaged checks passed; see VERIFICATION.md for scope.")
'''
(dst/'checks/run_checks.py').write_text(runner)
run=subprocess.run([sys.executable,str(dst/'checks/run_checks.py')],
                   text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
(dst/'CHECK_RESULTS.txt').write_text(run.stdout)
if run.returncode:
    print(run.stdout)
    raise SystemExit(run.returncode)

sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
sources={
  'date':'2026-09-08',
  'primary_pdf':{
    'url':'https://arxiv.org/pdf/2001.02333',
    'inspected_version':'arXiv:2001.02333v2; details in jeong-yoneda-stage.md',
    'local_path':str(src/'jeong-yoneda.pdf'),
    'sha256':sha(src/'jeong-yoneda.pdf')},
  'additional_primary_links':[
    'https://www.irphe.fr/~ledizes/web/JFM2006a.pdf',
    'https://arxiv.org/abs/1002.2489',
    'https://www-fourier.univ-grenoble-alpes.fr/~gallay/Handbook.pdf',
    'https://preprint.press.jhu.edu/ajm/sites/default/files/AJM-choi-jeong-FINAL.pdf'],
  'inherited_provenance':{
    'manifest':'../SOURCES.json',
    'sha256':sha(dst.parent/'SOURCES.json'),
    'fluid_lean_commit':'d0124689230b58b4f86e7b90ac59de06404b3b6b'},
  'verification':{'python':sys.version,'formal_kernel_replay':False,
                  'all_four_packaged_scripts_passed':True,
                  'navier_stokes_blowup_proved':False}}
(dst/'SOURCES.json').write_text(json.dumps(sources,indent=2)+'\n')
paths=sorted(p for p in dst.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
(dst/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(dst)}\n' for p in paths))
print(f'Packaged {len(paths)} files with four passing scripts in {dst}')

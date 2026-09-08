#!/usr/bin/env python3
"""Package and check the completed third-pass finite-transfer research."""
from pathlib import Path
import hashlib, json, shutil, subprocess, sys

root=Path(__file__).resolve().parents[2]
src=root/'work/pass3'
dst=root/'outputs/euler-ns-transfer/pass3'
(dst/'checks').mkdir(parents=True,exist_ok=True)
notes=['REPORT.md','VERIFICATION.md','quantitative-core-stage.md',
       'rotational-receiver.md','material-core-transfer.md',
       'receiver-flux-budget.md','terminal-profile-gap.md',
       'receiver-terminal-independent-audit.md','receiver-source-check.md']
scripts=['verify_quantitative_core_stage.py','verify_rotational_receiver.py',
         'verify_material_core_transfer.py','verify_receiver_flux.py']
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
print("All four packaged checks passed. See VERIFICATION.md for their scope.")
'''
(dst/'checks/run_checks.py').write_text(runner)
run=subprocess.run([sys.executable,str(dst/'checks/run_checks.py')],
                   text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
(dst/'CHECK_RESULTS.txt').write_text(run.stdout)
if run.returncode:
    print(run.stdout)
    raise SystemExit(run.returncode)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources={
 'date':'2026-09-08',
 'local_theory_primary_source':'https://terrytao.wordpress.com/2018/10/09/254a-notes-3-local-well-posedness-for-the-euler-equations/',
 'comparison_primary_sources':[
  'https://arxiv.org/abs/2311.01369v3',
  'https://arxiv.org/abs/1402.0290v3'],
 'inherited_manifests':[
  {'path':'../SOURCES.json','sha256':sha(dst.parent/'SOURCES.json')},
  {'path':'../pass2/SOURCES.json','sha256':sha(dst.parent/'pass2/SOURCES.json')}],
 'active_core_lemma':{
  'path':'../pass2/active-core-pressure.md',
  'sha256':sha(dst.parent/'pass2/active-core-pressure.md')},
 'verification':{
  'python':sys.version,'four_packaged_scripts_passed':True,
  'numerical_NS_evolution':False,'formal_kernel_replay':False,
  'finite_transfer_analytical_argument':True,
  'terminal_profile_return_proved':False,'NS_blowup_proved':False}}
(dst/'SOURCES.json').write_text(json.dumps(sources,indent=2)+'\n')
paths=sorted(p for p in dst.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
(dst/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(dst)}\n' for p in paths))
print(f'Packaged {len(paths)} files and four passing checks in {dst}')

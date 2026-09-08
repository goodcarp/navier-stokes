#!/usr/bin/env python3
"""Prepare this pass's new deliverables; never change older pass archives."""
import hashlib,json,shutil
from pathlib import Path

source=Path(__file__).resolve().parent
root=source.parents[1]/'outputs'/'euler-ns-transfer'
target=root/'pass11';experiments=target/'experiments';data=target/'data'
for p in (target,experiments,data,target/'checks'):p.mkdir(parents=True,exist_ok=True)
for p in source.iterdir():
    if p.name in ('package_pass11.py','update_current_state.py'):continue
    if p.suffix=='.md':
        if p.name=='codex-contribution-map.md':
            text=p.read_text().replace('../../outputs/euler-ns-transfer/','../')
            text=text.replace('](validated_error_hierarchy.py)','](experiments/validated_error_hierarchy.py)')
            text=text.replace('](compare_mac_endpoints.py)','](experiments/compare_mac_endpoints.py)')
            start=text.index('**Link mapping.**');end=text.index('\n\n',start)
            text=text[:start]+'**Link mapping.** This canonical copy links to the neighboring pass archives. The original review and its original source-path links are retained in the source workspace.'+text[end:]
            last=text.index('The top-level archived [STATUS]')
            text=text[:last]+'The current [STATUS](../STATUS.md) and [REPORT](../REPORT.md) incorporate pass11. Historical pass reports remain dated evidence. The continuing absence of a **validated whole-space** endpoint remains accurate.\n'
            (target/p.name).write_text(text)
        else:shutil.copy2(p,target/p.name)
    elif p.suffix in ('.py','.json'):shutil.copy2(p,experiments/p.name)
    elif p.suffix=='.png':shutil.copy2(p,target/p.name)
    elif p.suffix=='.npz':
        q=data/p.name
        if not q.exists() or q.stat().st_size!=p.stat().st_size:shutil.copy2(p,q)
(data/'README.md').write_text('''# Complete numerical endpoint arrays

The NPZ files contain final complex angular coefficient arrays and their
finite-cylinder MAC grid. They preserve phase information for independent
endpoint comparisons. They are numerical checkpoints, not continuous
whole-space solutions or formal certificates.

Large arrays are retained in the local disk archive. The unified GitHub
repository may distribute them as release assets rather than ordinary git
objects; consult its asset inventory for exact download locations. Verify
downloads with the separate `SHA256SUMS.txt` in this directory.

The text/code archive manifest intentionally excludes NPZ payloads. Its
inclusion of this directory's checksum manifest binds the expected payloads
without requiring a large-file checkout for the short verification runner.
''')
sources={'checked_utc':'2026-09-08 approximately 13:20–13:23',
    'primary_release_pdf':{'url':'https://cims.nyu.edu/~tristanb/euler.pdf','sha256':'97ef408bff09b4f6ed9f3867734d1eb2245f3f34e6334b28136c84c02d0ae8d8','relation':'byte-identical to previously inspected source; not duplicated in this pass'},
    'proof_repository':{'url':'https://github.com/tristanbuckmaster/fluid_lean','commit':'d0124689230b58b4f86e7b90ac59de06404b3b6b'},
    'scope':'Bounded public-channel status check, not evidence about unpublished work or a new proof audit.',
    'initial_evaluator':{'source':'../pass10/experiments/initial_full_field.py','sha256':hashlib.sha256((source/'initial_full_field.py').read_bytes()).hexdigest(),'relation':'unchanged selected datum'}}
(target/'SOURCES.json').write_text(json.dumps(sources,indent=2)+'\n')
runner='''#!/usr/bin/env python3
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
'''
(target/'checks'/'run_checks.py').write_text(runner)
def manifest(folder,paths):
    rows=[hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(folder).as_posix() for p in sorted(paths)]
    (folder/'SHA256SUMS.txt').write_text('\n'.join(rows)+'\n')
manifest(data,data.glob('*.npz'))
manifest(target,(p for p in target.rglob('*') if p.is_file() and p!=target/'SHA256SUMS.txt' and p.suffix!='.npz'))
print(target)

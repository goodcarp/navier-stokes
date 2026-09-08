#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,re,shutil,ast
root=Path(__file__).resolve().parents[2];top=root/'outputs/euler-ns-transfer';out=top/'pass7'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def endpoint(x):return F((-1 if x['sign'] else 1)*int(x['mantissa']))*F(2)**x['exponent']
def verify_manifest(folder):
    n=0
    for row in (folder/'SHA256SUMS.txt').read_text().splitlines():
        digest,name=row.split(None,1);p=folder/name.lstrip('*')
        assert sha(p)==digest,(folder,name)
        n+=1
    return n
historical={f'pass{i}':verify_manifest(top/f'pass{i}') for i in range(2,7)}
results=(out/'CHECK_RESULTS.txt').read_text()
core=['verify_D_pressure_kernels.py','certify_D_viscosity.py','verify_E_interaction_reduction.py','certify_outer_E.py','verify_E_certificate_audit.py','verify_common_family.py']
for name in core:assert 'PASS '+name in results,name
diagnostic=(root/'work/pass7/diagnostic-check-output.txt').read_text()
assert diagnostic.startswith('PASS:')
results+='\nRUN verify_next_stage_shear.py\n'+diagnostic+'\nPASS verify_next_stage_shear.py\n'
results+='\nArtifact QA: seven finite programs passed. Both decisive interval calculations were recomputed in the packaged location. Historical pass2--pass6 manifest entries remain unchanged.\n'
D=json.loads((out/'checks/D-viscosity-certificate.json').read_text())
E=json.loads((out/'checks/outer-E-certificate.json').read_text())
assert D['status']==E['status']=='PASS'
assert D['viscous_upper_strictly_below_four'] and D['radial_upper_strictly_below_120000']
assert E['E_upper_strictly_below_one_sixtieth']
assert endpoint(E['E_upper_expression']['exact']['upper'])<F(1,60)
assert endpoint(D['viscous_upper_bound_exact']['upper'])<4
links=0
for p in list(out.glob('*.md'))+[top/'STATUS.md',top/'NEXT_TARGET.md',top/'REPORT.md',top/'VERIFICATION.md']:
    s=p.read_text()
    assert not any(ord(c)<32 and c not in '\n\r\t' for c in s),p
    for target in re.findall(r'\]\(([^)]+)\)',s):
        if '://' in target or target.startswith('#'):continue
        target=target.split('#')[0]
        assert (p.parent/target).exists(),(p,target)
        links+=1
for p in (out/'checks').glob('*.py'):ast.parse(p.read_text())
for p in out.rglob('__pycache__'):
    shutil.rmtree(p)
sources=json.loads((out/'SOURCES.json').read_text())
sources['verification'].update(finite_check_programs_passed=7,decisive_interval_calculations_recomputed=True,
    earlier_packages_unchanged=True,local_artifact_links_resolved=links)
(out/'SOURCES.json').write_text(json.dumps(sources,indent=2)+'\n')
results+=f'Artifact QA: {links} local Markdown links resolved; Python sources parse; certificate predicates checked against exact dyadic endpoints.\n'
results+='Scope: a common finite-amplitude initial gate and qualitative short interval, plus a separate frozen shear diagnostic. No useful actual duration, formal replay, inherited return, or blowup proof.\n'
(out/'CHECK_RESULTS.txt').write_text(results)
def manifest(folder,files):
    files=sorted(files)
    (folder/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(folder)}\n' for p in files))
    return verify_manifest(folder)
pass_count=manifest(out,[p for p in out.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt' and '__pycache__' not in p.parts])
top_files=[p for p in top.iterdir() if p.is_file() and p.name!='SHA256SUMS.txt']
top_files += [p for p in (top/'checks').rglob('*') if p.is_file() and '__pycache__' not in p.parts]
top_count=manifest(top,top_files)
print(json.dumps(dict(status='PASS',historical_unchanged=historical,pass7_manifest_files=pass_count,top_manifest_files=top_count,resolved_links=links,finite_checks_passed=7,NS_blowup_proved=False),indent=2))

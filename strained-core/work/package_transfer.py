"""Package original analysis and finite checks; do not redistribute source papers."""
from pathlib import Path
import datetime
import hashlib
import json
import shutil

OUT=Path('outputs/euler-ns-transfer')
CHECKS=OUT/'checks'
CHECKS.mkdir(parents=True,exist_ok=True)
copies={
    'phase_heat_check.py':'phase_heat_check.py',
    'viscosity-check.py':'viscosity_check.py',
    'verify_direction_scout.py':'direction_check.py',
}
for source,dest in copies.items():
    text=(Path('work')/source).read_text()
    text=text.replace('viscosity-audit.md','../REPORT.md').replace('direction-scout.md','../REPORT.md')
    (CHECKS/dest).write_text(text)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

proof=json.loads(Path('work/sources/proof-scope/manifest.json').read_text())
inputs=[]
for filename,url in [
    ('euler.pdf','https://cims.nyu.edu/~tristanb/euler.pdf'),
    ('statement.pdf','https://cims.nyu.edu/~tristanb/statement.pdf'),
]:
    path=Path('work/sources')/filename
    inputs.append({'name':filename,'url':url,'sha256':sha(path),'bytes':path.stat().st_size,
                   'inspected_copy_mtime_utc':datetime.datetime.fromtimestamp(path.stat().st_mtime,datetime.timezone.utc).isoformat()})
archive=Path('~/<core-line-working-tree>/navier-stokes-research-astra')
local=[]
for rel in ['README.md','FINDINGS.md','docs/outputs/RETURN_Deepest-Think-Astra_PASS14.md',
            'docs/outputs/navier-stokes-pass14/NEXT_TARGET.md','docs/outputs/navier-stokes-pass28/NEXT_TARGET.md']:
    local.append({'repository':'goodcarp/navier-stokes-research-astra','path':rel,'sha256':sha(archive/rel)})
record={
    'prepared_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'papers_downloaded':inputs,
    'additional_primary_sources_read':[
        'https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf',
        'https://www.cscamm.umd.edu/publications/AxisymmetricFlow_CS-08-45.pdf'],
    'proof_scope':{'repository':'tristanbuckmaster/fluid_lean',**proof,
                   'method':'Static selected-file reading and textual proposition comparison only.',
                   'main_theorem_proposition_matches_after_whitespace_normalization':True,
                   'full_build':False,'kernel_replay':False,'dependency_closure_audit':False},
    'local_research_inputs':local,
    'limits':'Paper SHA256 identifies the inspected bytes, not their mathematical correctness. Source papers and proof code are not redistributed in this package.'
}
(OUT/'SOURCES.json').write_text(json.dumps(record,indent=2)+'\n')
print('Packaged three original finite-check programs and source provenance.')

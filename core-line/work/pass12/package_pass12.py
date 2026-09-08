#!/usr/bin/env python3
"""Copy finalized pass12 deliverables and preserve historical passes."""
from pathlib import Path
import hashlib,json,shutil,re

ROOT=Path(__file__).resolve().parents[2]
SRC=Path(__file__).resolve().parent
TOP=ROOT/'outputs/euler-ns-transfer'
DEST=TOP/'pass12'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def manifest(directory,files):
    (directory/'SHA256SUMS.txt').write_text(''.join(f'{digest(p)}  {p.relative_to(directory)}\n' for p in sorted(files)))

def package():
    for d in ('','experiments','data','provenance'):(DEST/d).mkdir(parents=True,exist_ok=True)
    for p in SRC.iterdir():
        if not p.is_file():continue
        if p.name=='package_pass12.py':continue
        if p.suffix=='.md':out=DEST/p.name
        elif p.suffix=='.npz':out=DEST/'data'/p.name
        elif p.suffix in ('.py','.json','.txt','.log'):out=DEST/'experiments'/p.name
        else:continue
        shutil.copy2(p,out)
        if p.suffix=='.md':
            def fix_link(match):
                name=match.group(1)
                if '/' not in name and (SRC/name).is_file() and Path(name).suffix in ('.py','.json','.txt','.log'):
                    return ']('+'experiments/'+name+')'
                return match.group(0)
            out.write_text(re.sub(r'\]\(([^)]+)\)',fix_link,out.read_text()))
    if (SRC/'provenance').exists():
        shutil.copytree(SRC/'provenance',DEST/'provenance',dirs_exist_ok=True)
    # Preserve all earlier pass archives. Validate the records as actually scoped.
    counts={}
    for old in sorted(TOP.glob('pass[0-9]*')):
        if old==DEST:continue
        for sha in old.rglob('SHA256SUMS.txt'):
            count=0
            for line in sha.read_text().splitlines():
                expected,name=line.split('  ',1)
                assert digest(sha.parent/name)==expected,(sha,name)
                count+=1
            counts[str(sha.relative_to(TOP))]=count
    manifest(DEST/'data',list((DEST/'data').glob('*.npz')))
    manifest(DEST,[p for p in DEST.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt' and 'data' not in p.relative_to(DEST).parts])
    links=[]
    for md in DEST.glob('*.md'):
        prose=re.sub(r'```.*?```','',md.read_text(),flags=re.S)
        prose=re.sub(r'`[^`]*`','',prose)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('http:','https:','#')):continue
            if not (md.parent/target.split('#')[0]).exists():links.append([md.name,target])
    assert not links,links
    print(json.dumps({'status':'PASS','historical_manifest_counts':counts,
        'pass12_members':len((DEST/'SHA256SUMS.txt').read_text().splitlines()),
        'data_members':len(list((DEST/'data').glob('*.npz'))),
        'top_level_markdown_link_failures':links},indent=2))

if __name__=='__main__':package()

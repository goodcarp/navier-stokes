#!/usr/bin/env python3
"""Append the owner-frozen pass12 without recollecting the parallel campaign.

Run only after the canonical-ready instruction. Existing source records and
older packages are verified before any intentional current-file replacement.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE='a37d92e6053a5c9de03cc078d65634cf01807f19'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--workspace',type=Path,required=True)
    args=ap.parse_args()
    workspace=args.workspace.resolve()
    snapshot=datetime.now(timezone.utc).isoformat()
    upstream=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    manifest_path=ROOT/'manifests/collected-sources.json'
    manifest=json.loads(manifest_path.read_text())
    records={r['path']:r for r in manifest['records']}
    relocated_arrays=[]
    # Git renames do not move ignored local arrays. Move only these verified
    # known assets into the collaborator's new core-line layout.
    relocated_source_versions=[]
    for row in list(records.values()):
        target=ROOT/row['path']
        if not target.exists() and row['path'].startswith('core-line/'):
            old=ROOT/('codex/'+row['path'].split('/',1)[1])
            if old.is_file() and sha(old)==row['sha256']:
                target.parent.mkdir(parents=True,exist_ok=True)
                old.rename(target)
                relocated_arrays.append(dict(old=str(old.relative_to(ROOT)),new=row['path']))
        if not target.is_file():
            raise RuntimeError('Existing collected source missing: '+row['path'])
        if sha(target)!=row['sha256']:
            # A collaborator changed one working contribution map during the
            # neutral-name migration. Preserve its original source version,
            # and identify the current version as upstream Git content.
            original=row['path'].replace('core-line/','codex/',1).replace('clock-line/','claude/',1)
            data=subprocess.check_output(['git','show',
                '3dc241bc8f22abbd10fe6bb25d81566f5d3f0d5e:'+original],cwd=ROOT)
            if hashlib.sha256(data).hexdigest()!=row['sha256']:
                raise RuntimeError('Cannot preserve original collected source: '+row['path'])
            dest=ROOT/'provenance/first-snapshot-3dc241b'/original
            dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(data)
            preserved=dict(row,path=str(dest.relative_to(ROOT)),previous_archive_path=row['path'])
            records[preserved['path']]=preserved
            records[row['path']]=dict(path=row['path'],source='collaborator-upstream',
                source_path=row['path'],git_commit=upstream,bytes=target.stat().st_size,
                sha256=sha(target),storage='git')
            relocated_source_versions.append(dict(path=row['path'],
                original_preserved_at=preserved['path'],original_sha256=row['sha256'],
                upstream_sha256=sha(target)))
    if relocated_source_versions:
        manifest['source_roots']['collaborator-upstream']='git:https://github.com/goodcarp/forced-route-2026-09@'+upstream
    # Preserve the three original baseline documents edited by the intervening
    # collaborator commits, rather than reverting their current text.
    base_path=ROOT/'manifests/preserved-base.json'
    baseline=json.loads(base_path.read_text())
    relocated_originals=[]
    for row in baseline['files']:
        current=ROOT/row['preserved_path']
        if current.is_file() and sha(current)==row['sha256']:continue
        original=row['original_path']
        data=subprocess.check_output(['git','show',BASE+':'+original],cwd=ROOT)
        if hashlib.sha256(data).hexdigest()!=row['sha256']:
            raise RuntimeError('Original baseline hash mismatch: '+original)
        dest=ROOT/'provenance/original-base-a37d92e'/original
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_bytes(data)
        relocated_originals.append(dict(original_path=original,
            previous_preserved_path=row['preserved_path'],preserved_path=str(dest.relative_to(ROOT)),
            sha256=row['sha256']))
        row['preserved_path']=str(dest.relative_to(ROOT))
    base_path.write_text(json.dumps(baseline,indent=2)+'\n')
    historical={str(p.relative_to(ROOT)):sha(p)
        for n in range(2,12)
        for p in (ROOT/f'core-line/outputs/euler-ns-transfer/pass{n}').rglob('*')
        if p.is_file()}
    sources=[]
    for relative in ['work/pass12','outputs/euler-ns-transfer/pass12']:
        for p in sorted((workspace/relative).rglob('*')):
            if p.is_file() and not any(x in ('__pycache__','.DS_Store') for x in p.parts) and p.suffix!='.pyc':
                sources.append(p)
    sources += sorted(p for p in (workspace/'outputs/euler-ns-transfer').iterdir() if p.is_file())
    changed_current=[]
    copied=[]
    for source in sources:
        relative=source.relative_to(workspace)
        target=ROOT/'core-line'/relative
        digest=sha(source)
        path=str(target.relative_to(ROOT))
        if path in records and records[path]['sha256']!=digest:
            changed_current.append(dict(path=path,previous_sha256=records[path]['sha256'],sha256=digest))
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target)
        if sha(target)!=digest or sha(source)!=digest:
            raise RuntimeError('Source changed during frozen copy: '+str(relative))
        row=dict(path=path,source='codex-workspace',source_path=str(relative),
                 bytes=target.stat().st_size,sha256=digest,
                 storage='release-asset-and-local' if target.suffix=='.npz' or target.stat().st_size>=20000000 else 'git',
                 snapshot='pass12-2026-09-08')
        records[path]=row
        copied.append(path)
    for name,digest in historical.items():
        if sha(ROOT/name)!=digest:raise RuntimeError('Historical pass changed: '+name)
    manifest['records']=sorted(records.values(),key=lambda r:r['path'])
    manifest['final_pass12']=True
    manifest['pass12_snapshot_utc']=snapshot
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    ignores=['# Generated caches and separately downloadable large research arrays.',
             '**/__pycache__/','**/*.pyc','**/.DS_Store','**/.lake/']
    ignores+=['/'+r['path'] for r in manifest['records'] if r['storage']!='git']
    (ROOT/'.gitignore').write_text('\n'.join(ignores)+'\n')
    result=dict(status='PASS',scope='Frozen source copy and historical byte preservation only.',
        snapshot_utc=snapshot,upstream_before_addendum=upstream,
        first_snapshot_commit='3dc241bc8f22abbd10fe6bb25d81566f5d3f0d5e',
        naming={'codex':'core-line','claude':'clock-line'},
        copied_files=len(copied),historical_pass2_through11_files_checked=len(historical),
        existing_source_records_checked=len(manifest['records'])-len([p for p in copied if 'pass12' in Path(p).parts]),
        relocated_local_arrays=relocated_arrays,relocated_originals=relocated_originals,
        relocated_source_versions=relocated_source_versions,
        intentional_current_file_changes=changed_current,
        final_source_records=len(manifest['records']))
    (ROOT/'manifests/pass12-integration.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()

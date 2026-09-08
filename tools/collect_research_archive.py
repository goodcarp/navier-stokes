#!/usr/bin/env python3
"""Collect the authorized research trees without modifying their originals.

This is an archival copier, not a mathematical verifier. Pass 11 is included
only with --final-pass11 after its owner has frozen the package.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CACHE_DIRS = {'.git', '.lake', '__pycache__', '.pytest_cache', '.venv', 'node_modules'}
CACHE_FILES = {'.DS_Store'}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def walk(root):
    for base, dirs, files in os.walk(root):
        dirs[:] = sorted(d for d in dirs if d not in CACHE_DIRS)
        for name in sorted(files):
            p = Path(base) / name
            if name in CACHE_FILES or p.suffix in {'.pyc', '.pyo'}:
                continue
            if p.is_symlink():
                raise ValueError('Source symlink requires explicit review: ' + str(p))
            yield p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workspace', type=Path, required=True)
    ap.add_argument('--campaign', type=Path, required=True)
    ap.add_argument('--final-pass11', action='store_true')
    args = ap.parse_args()
    records, omissions = [], []

    def copy_one(source, destination, source_label, source_relative):
        before = sha(source)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        after = sha(destination)
        if before != after or before != sha(source):
            raise RuntimeError('Source changed during snapshot: ' + str(source))
        records.append(dict(path=str(destination.relative_to(REPO)),
            source=source_label, source_path=str(source_relative),
            bytes=destination.stat().st_size, sha256=after,
            storage='release-asset-and-local' if destination.suffix=='.npz'
                or destination.stat().st_size >= 20_000_000 else 'git'))

    def copy_tree(source, destination, label, source_base=None, allow=None):
        source_base = source if source_base is None else source_base
        for p in walk(source):
            if allow and not allow(p):
                omissions.append(dict(source=label,path=str(p.relative_to(source_base)),
                    reason='External fetch scratch; source URLs/hashes are retained in the curated package.'))
                previous=destination/p.relative_to(source)
                if previous.exists():previous.unlink()
                continue
            copy_one(p, destination/p.relative_to(source), label, p.relative_to(source_base))

    # The original 601-file repository is preserved independently of this copy.
    claude_subtrees=('lower','write','gaps','s3close',
        'rebuild/far-near-kernel-lemma','rebuild/refuter-far-near-kernel',
        'sharp/viscous-numerics','sharp/refuter-correctness',
        'sharp/exact-first-order','sharp/literature-short-time',
        'bfg/corner-numerics','bfg/source-adjudication',
        'verify-astra/refuter-A1-entropy-clock','toys')
    for name in claude_subtrees:
        relative=Path('deepest-think/DTC-2026-09-06')/name
        copy_tree(args.campaign/relative,REPO/'clock-line/source-campaign'/relative,
                  'claude-campaign',args.campaign,
                  lambda p:p.suffix not in {'.pdf','.html'})
    relative=Path('external/alpoge-buckmaster-2026-09-08')
    copy_tree(args.campaign/relative,REPO/'clock-line/source-campaign'/relative,
              'claude-campaign',args.campaign)
    relative=Path('scans/A3_SCAN_2026-09-08_ALPOGE_BUCKMASTER.md')
    copy_one(args.campaign/relative,REPO/'clock-line/source-campaign'/relative,
             'claude-campaign',relative)
    relative=Path('sources/pdfs/alpoge-buckmaster-2026-09-08')
    copy_tree(args.campaign/relative,REPO/'clock-line/source-campaign'/relative,
              'claude-campaign',args.campaign,
              lambda p:p.suffix not in {'.pdf','.html'})

    # Working files are retained as provenance, separately from reader-facing
    # authoritative pass packages. No general Codex/Claude account stores enter.
    for p in sorted((args.workspace/'work').iterdir()):
        if p.is_file() and p.suffix in {'.md','.py','.json','.txt'}:
            copy_one(p,REPO/'core-line/work'/p.name,'codex-workspace',p.relative_to(args.workspace))
    for number in range(2,12 if args.final_pass11 else 11):
        source=args.workspace/'work'/('pass'+str(number))
        if not source.exists():continue
        def allow(p):
            return p.suffix not in {'.html','.pdf'} and (number!=11 or
                p.name not in {'fluid-current-commit.json','fluid-current-tree.json'})
        copy_tree(source,REPO/'core-line/work'/source.name,'codex-workspace',args.workspace,allow)
    copy_tree(args.workspace/'work/sources',REPO/'core-line/work/sources',
              'codex-workspace',args.workspace,
              lambda p:p.suffix not in {'.html','.pdf','.png','.jpeg','.jpg'})

    source=args.workspace/'outputs/euler-ns-transfer'
    for p in walk(source):
        rel=p.relative_to(source)
        if not args.final_pass11 and (len(rel.parts)==1 or rel.parts[0]=='pass11'):
            continue
        copy_one(p,REPO/'core-line/outputs/euler-ns-transfer'/rel,
                 'codex-workspace',p.relative_to(args.workspace))

    manifest=dict(schema_version=1,snapshot_utc=datetime.now(timezone.utc).isoformat(),
        final_pass11=args.final_pass11,
        source_roots={'codex-workspace':str(args.workspace),'claude-campaign':str(args.campaign)},
        scope='Forced-route research and the log-clock sources explicitly cited by its assembly; '
              'no unrelated account/session stores, caches or older campaigns.',
        records=sorted(records,key=lambda x:x['path']),omissions=omissions)
    target=REPO/'manifests/collected-sources.json'
    target.write_text(json.dumps(manifest,indent=2)+'\n')
    large=[r for r in records if r['storage']!='git']
    ignored=['# Generated caches and separately downloadable large research arrays.',
             '**/__pycache__/','**/*.pyc','**/.DS_Store','**/.lake/']
    ignored.extend('/'+r['path'] for r in large)
    (REPO/'.gitignore').write_text('\n'.join(ignored)+'\n')
    print(json.dumps(dict(files=len(records),bytes=sum(r['bytes'] for r in records),
        large_files=len(large),large_bytes=sum(r['bytes'] for r in large),
        final_pass11=args.final_pass11,omissions=len(omissions)),indent=2))


if __name__=='__main__':main()

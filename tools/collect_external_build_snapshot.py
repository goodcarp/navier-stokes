#!/usr/bin/env python3
"""Preserve a bounded, read-only snapshot of the referenced external Lean build.

This never launches Lean or modifies the source checkout. Run after the main
collector, supplying the explicitly referenced euler-blowup directory.
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['build.log', 'run_build.sh', 'lean-toolchain', 'lake-manifest.json',
         'lakefile.toml', 'comparator.json', 'scripts/PrintAxioms.lean']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--source', type=Path, required=True)
    args = ap.parse_args()
    source = args.source.resolve()
    stamp = datetime.now(timezone.utc).isoformat()
    inventory = ROOT / 'manifests/collected-sources.json'
    doc = json.loads(inventory.read_text())
    label = 'external-euler-build'
    dest = ROOT / 'clock-line/replication-snapshot-2026-09-08'
    rows = []
    for name in NAMES:
        p = source / name
        data = p.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        target = dest / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        if hashlib.sha256(p.read_bytes()).hexdigest() != digest:
            raise RuntimeError(f'Source changed during bounded copy: {name}')
        rows.append(dict(path=str(target.relative_to(ROOT)), source=label,
                         source_path=name, bytes=len(data), sha256=digest, storage='git'))
    doc['source_roots'][label] = str(source)
    doc['supplemental_snapshot_utc'] = stamp
    doc['records'] = sorted([r for r in doc['records'] if r['source'] != label] + rows,
                            key=lambda r: r['path'])
    inventory.write_text(json.dumps(doc, indent=2) + '\n')
    print(json.dumps(dict(snapshot_utc=stamp, copied=len(rows), records=rows), indent=2))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Write the exact large-data release inventory; does not upload anything."""
import argparse
import hashlib
import json
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--tag',default='research-archive-2026-09-08')
    args=ap.parse_args()
    manifest=json.loads((ROOT/'manifests/collected-sources.json').read_text())
    existing_path=ROOT/'manifests/data-assets.json'
    existing=json.loads(existing_path.read_text()) if existing_path.exists() else {}
    previous={a['sha256']:a for a in existing.get('assets',[])}
    assets=[]
    names=set()
    by_hash={}
    for row in manifest['records']:
        if row['storage']=='git':continue
        if row['sha256'] in by_hash:
            by_hash[row['sha256']]['local_paths'].append(row['path'])
            continue
        path=Path(row['path'])
        pass_name=next((part for part in path.parts if re.fullmatch(r'pass\d+',part)),None)
        prior=previous.get(row['sha256'])
        name=prior['name'] if prior else ((pass_name+'-' if pass_name else '')+path.name)
        if name in names:name=row['sha256'][:12]+'-'+name
        names.add(name)
        url=prior['url'] if prior else ('https://github.com/goodcarp/forced-route-2026-09/'
            'releases/download/'+args.tag+'/'+name)
        release_tag=url.split('/releases/download/',1)[1].split('/',1)[0]
        kind=('Finite-cylinder Fourier/MAC endpoint' if pass_name=='pass11' else
              'Rejected sampled-fit potential coefficients' if path.name.startswith('compact-fit-') else
              'Continuous-L2-fit potential coefficients' if path.name.startswith('continuous-fit-') else
              'Research data')
        item=dict(name=name,path=row['path'],local_paths=[row['path']],bytes=row['bytes'],
            sha256=row['sha256'],url=url,release_tag=release_tag,kind=kind)
        assets.append(item)
        by_hash[row['sha256']]=item
    tags=sorted(set(a['release_tag'] for a in assets))
    result=dict(schema_version=2,repository='goodcarp/forced-route-2026-09',
        release_tags=tags,release_tag_for_new_assets=args.tag,
        scope='Expected byte-identical release assets and local placements; '
          'this inventory does not by itself establish upload completion.',assets=assets)
    (ROOT/'manifests/data-assets.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Full numerical data','',
       'Large numerical arrays stay in the organized local copy and are distributed',
       'as separate GitHub release assets. The original snapshot and later',
       'addenda keep their own immutable asset URLs:',
       '',
       *['- ['+tag+'](https://github.com/goodcarp/forced-route-2026-09/releases/tag/'+tag+')' for tag in tags],
       '',
       'They are excluded from Git blobs; their exact bytes remain part of the',
       '[source inventory](../manifests/collected-sources.json). The release must',
       'exist and its assets must be downloaded before the large-data part of an',
       'integrity check can pass in a fresh clone.','',
       '| Asset | Role | Local placement relative to repository root | Bytes |',
       '|---|---|---|---:|']
    for a in assets:lines.append('| ['+a['name']+']('+a['url']+') | '+
        a['kind']+' | '+'<br>'.join('`'+p+'`' for p in a['local_paths'])+' | '+str(a['bytes'])+' |')
    lines+=['','Use [the data manifest](../manifests/data-assets.json) for SHA-256 values.',
       'Pass11 files are full Fourier/MAC endpoint fields of the finite-cylinder',
       'experiment. Pass12 files are compact potential reconstructions of an',
       'increment from the same endpoint; both rejected sampled fits and corrected',
       'continuous-L2 fits are preserved. Their induced velocities are C6/H7.',
       'None is a validated whole-space NS state. Keep the accompanying',
       'JSON parameters, solver source hashes and comparison report with them.',
       'No interpolated or independently phase-aligned replacement is the original',
       'checkpoint.','',
       'For a private repository, download through an authenticated GitHub session',
       'or `gh release download`, then place each file at its listed local path.',
       'Run `python3 tools/check_archive.py` afterward. The archive checker validates',
       'the stored hashes; it does not execute arrays or certify the PDE.','']
    (ROOT/'docs/data-assets.md').write_text('\n'.join(lines))
    print(json.dumps(dict(assets=len(assets),bytes=sum(a['bytes'] for a in assets)),indent=2))


if __name__=='__main__':main()

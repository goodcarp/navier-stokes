#!/usr/bin/env python3
"""Verify archive preservation, copied source hashes and curated local links.

This checks the archive's integrity, not the truth of archived mathematics.
Historical link problems are inventoried rather than silently repaired.
"""
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT=Path(__file__).resolve().parents[1]


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def links(path):
    text=path.read_text(errors='replace')
    blank=lambda m:'\n'*m.group(0).count('\n')
    for pattern in [r'```.*?```',r'\\\[.*?\\\]',r'\\\(.*?\\\)',
                    r'(?<!\\)\$\$.*?(?<!\\)\$\$',r'`[^`\n]*`',
                    r'(?<!\\)\$[^$\n]*?(?<!\\)\$']:
        text=re.sub(pattern,blank,text,flags=re.S)
    for m in re.finditer(r'!?\[[^\]]*\]\(([^\n]*?)\)',text):
        target=m.group(1).strip()
        if target.startswith('<') and '>' in target:
            target=target[1:target.index('>')]
        else:target=target.split(' "',1)[0].split(" '",1)[0]
        if not target or target.startswith('#'):continue
        parsed=urlsplit(target)
        if parsed.scheme or parsed.netloc:continue
        rel=unquote(parsed.path)
        if not rel:continue
        dest=(path.parent/rel).resolve()
        if not dest.exists():
            yield dict(file=str(path.relative_to(ROOT)),target=target,
                       line=text[:m.start()].count('\n')+1)


def main():
    base=json.loads((ROOT/'manifests/preserved-base.json').read_text())
    copied=json.loads((ROOT/'manifests/collected-sources.json').read_text())
    hash_failures=[]
    for row in base['files']:
        p=ROOT/row['preserved_path']
        if not p.is_file() or digest(p)!=row['sha256']:
            hash_failures.append(dict(kind='preserved-base',path=row['preserved_path']))
    for row in copied['records']:
        p=ROOT/row['path']
        if not p.is_file() or digest(p)!=row['sha256']:
            hash_failures.append(dict(kind='collected-source',path=row['path']))
    curated=[ROOT/'README.md',ROOT/'EXPLAINER.md',ROOT/'clock-line/README.md',ROOT/'core-line/README.md']
    curated+=list((ROOT/'docs').glob('*.md'))
    missing_curated=[failure for p in curated for failure in links(p)]
    historical=[]
    for name in ['theorems','analysis','audit','lean','replication','strained-core','clock-line/source-campaign','core-line/work','core-line/outputs']:
        for p in (ROOT/name).rglob('*.md'):
            if '.lake' not in p.parts:historical.extend(links(p))
    # Report file names only; never emit a credential-like match's contents.
    secret_patterns={
        'private-key':re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
        'github-token':re.compile(rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,})\b'),
        'openai-token':re.compile(rb'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{40,}\b')}
    secret_hits=[]
    for row in copied['records']:
        p=ROOT/row['path']
        if p.suffix.lower() in {'.pdf','.png','.jpg','.jpeg','.npz','.zip','.gz'}:continue
        data=p.read_bytes()
        for name,pattern in secret_patterns.items():
            if pattern.search(data):secret_hits.append(dict(path=row['path'],pattern=name))
    result=dict(status='PASS' if not(hash_failures or missing_curated or secret_hits)
                    and copied['final_pass11'] else 'INCOMPLETE_OR_FAILED',
        scope='File integrity, curated local path existence and narrow credential-pattern scan only; '
              'no mathematical proof, remote-link health or anchor validation.',
        final_pass11=copied['final_pass11'],final_pass12=copied.get('final_pass12',False),
        preserved_base_files=len(base['files']),
        copied_source_files=len(copied['records']),hash_failures=hash_failures,
        curated_link_failures=missing_curated,credential_pattern_hits=secret_hits,
        historical_link_warnings=historical,
        historical_warning_policy='Recorded without changing verbatim source files; '
          'absolute machine paths, historical copies, and external-fetch omissions may be nonportable.')
    (ROOT/'manifests/archive-check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='historical_link_warnings'},indent=2))
    print('Historical link warnings:',len(historical))
    return 0 if result['status']=='PASS' else 1


if __name__=='__main__':raise SystemExit(main())

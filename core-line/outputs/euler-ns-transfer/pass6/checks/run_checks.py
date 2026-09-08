#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys
here=Path(__file__).resolve().parent
programs=sorted(here.glob('verify_*.py'))
for p in programs:
    print('CHECK '+p.name,flush=True)
    r=subprocess.run([sys.executable,str(p)],cwd=here,capture_output=True,text=True)
    print(r.stdout,end='',flush=True)
    if r.stderr: print(r.stderr,end='',flush=True)
    if r.returncode: raise SystemExit(r.returncode)
print(f'PASS: {len(programs)} finite checks. No return or NS blowup certified.',flush=True)

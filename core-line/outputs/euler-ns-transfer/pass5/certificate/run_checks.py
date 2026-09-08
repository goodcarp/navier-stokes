#!/usr/bin/env python3
"""Run the finite symbolic, interval-consistency and exact arithmetic checks."""
from pathlib import Path
import subprocess,sys
here=Path(__file__).resolve().parent
programs=sorted(here.glob('verify_*.py'))
for program in programs:
    print('CHECK '+program.name,flush=True)
    result=subprocess.run([sys.executable,str(program)],cwd=here,capture_output=True,text=True)
    print(result.stdout,end='',flush=True)
    if result.stderr:
        print(result.stderr,end='',flush=True)
    if result.returncode:
        raise SystemExit(result.returncode)
print(f'PASS: {len(programs)} finite checks. No return map or blowup is certified.',flush=True)

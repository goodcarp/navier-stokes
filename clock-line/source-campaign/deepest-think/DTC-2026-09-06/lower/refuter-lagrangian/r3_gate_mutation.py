#!/usr/bin/env python3
"""r3 (FL-043): does check_constants.py return the same verdict for every outcome?
Mutate one stored number at a time in the *_results.json the gate reads, re-run the gate,
record whether it still exits 0.  A number the gate cannot see is a number the gate does not check.
"""
import json, os, shutil, subprocess, sys, copy, itertools

SRC = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/lower/prove-lagrangian"
WORK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gate_mut")
FILES = ["s%d_results.json" % i for i in (1,2,3,4,5,6,7)]

def leaves(o, path=()):
    if isinstance(o, dict):
        for k,v in o.items(): yield from leaves(v, path+(k,))
    elif isinstance(o, list):
        for i,v in enumerate(o): yield from leaves(v, path+(i,))
    elif isinstance(o, (int,float)) and not isinstance(o, bool):
        yield path, o

def setpath(o, path, val):
    cur = o
    for p in path[:-1]: cur = cur[p]
    cur[path[-1]] = val

def run_gate():
    r = subprocess.run([sys.executable, "check_constants.py"], cwd=WORK,
                       capture_output=True, text=True)
    return r.returncode

os.makedirs(WORK, exist_ok=True)
for f in FILES + ["check_constants.py"]:
    shutil.copy(os.path.join(SRC,f), WORK)
assert run_gate()==0, "baseline gate must pass"

orig = {f: json.load(open(os.path.join(SRC,f))) for f in FILES}
unseen, seen, total = [], 0, 0
for f in FILES:
    for path, val in leaves(orig[f]):
        total += 1
        mut = copy.deepcopy(orig[f])
        newv = val*1.5 + 0.37 if val != 0 else 0.37
        setpath(mut, path, newv)
        json.dump(mut, open(os.path.join(WORK,f),"w"))
        rc = run_gate()
        if rc == 0: unseen.append([f, list(map(str,path)), val, newv])
        else: seen += 1
        json.dump(orig[f], open(os.path.join(WORK,f),"w"))   # restore

res = {"numbers_stored": total, "numbers_the_gate_catches": seen,
       "numbers_the_gate_cannot_see": len(unseen),
       "unseen_sample": unseen}
json.dump(res, open(os.path.join(os.path.dirname(WORK),"r3_results.json"),"w"), indent=1)
print(json.dumps({k:v for k,v in res.items() if k!="unseen_sample"}, indent=1))
print("\nfirst 40 unseen:")
for u in unseen[:40]: print("  ", u[0], "/".join(u[1]), u[2], "->", u[3])

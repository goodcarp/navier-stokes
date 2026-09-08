#!/usr/bin/env python3
"""x4 (FL-043) -- can check_constants.py return the same verdict for every outcome?
Two tests:
 (T1) STATIC: find checks whose 'got' expression cannot depend on any computed result.
 (T2) DYNAMIC: mutate one stored number at a time in the *_results.json the gate reads and
      re-run the gate; a number that never flips the exit code is a number the gate cannot see.
"""
import json, os, shutil, subprocess, sys, copy
SRC="~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/lower/prove-lagrangian"
HERE=os.path.dirname(os.path.abspath(__file__))
WORK=os.path.join(HERE,"gate_mut")
FILES=["s%d_results.json"%i for i in (1,2,3,4,5,6,7)]
os.makedirs(WORK,exist_ok=True)
for f in FILES+["check_constants.py"]: shutil.copy(os.path.join(SRC,f),WORK)
def run(): return subprocess.run([sys.executable,"check_constants.py"],cwd=WORK,capture_output=True,text=True)
base=run(); assert base.returncode==0, base.stdout[-500:]
OUT={"baseline_rc":base.returncode,"baseline_last_line":base.stdout.strip().splitlines()[-1]}
OUT["printed_check_lines"]=sum(1 for L in base.stdout.splitlines() if L.startswith(("PASS","FAIL")))

# ---- T1 static: literal-only checks
src=open(os.path.join(SRC,"check_constants.py")).read()
import re
lits=[]
for m in re.finditer(r"ck\((.*?)\)\n", src, re.S):
    body=m.group(1)
    if "s1[" in body or "s2[" in body or "s3[" in body or "s4[" in body or "s5[" in body or "s6[" in body or "s7[" in body:
        if "*0+" in body: lits.append(("MASKED_TO_ZERO", body.strip()))
    else:
        lits.append(("NO_STORED_INPUT", body.strip()))
OUT["static_unfalsifiable_checks"]=lits

# ---- T2 dynamic
def leaves(o,path=()):
    if isinstance(o,dict):
        for k,v in o.items(): yield from leaves(v,path+(k,))
    elif isinstance(o,list):
        for i,v in enumerate(o): yield from leaves(v,path+(i,))
    elif isinstance(o,(int,float)) and not isinstance(o,bool): yield path,o
def setp(o,path,val):
    cur=o
    for p in path[:-1]: cur=cur[p]
    cur[path[-1]]=val
orig={f:json.load(open(os.path.join(SRC,f))) for f in FILES}
unseen=[]; seen=0; total=0
for f in FILES:
    for path,val in leaves(orig[f]):
        total+=1
        mut=copy.deepcopy(orig[f]); setp(mut,path,val*1.5+0.37 if val!=0 else 0.37)
        json.dump(mut,open(os.path.join(WORK,f),"w"))
        rc=run().returncode
        if rc==0: unseen.append([f,"/".join(map(str,path)),val])
        else: seen+=1
        json.dump(orig[f],open(os.path.join(WORK,f),"w"))
OUT["stored_numbers"]=total; OUT["gate_catches"]=seen; OUT["gate_blind_to"]=len(unseen)
OUT["blind_sample"]=unseen[:60]
print(json.dumps({k:v for k,v in OUT.items() if k!="blind_sample"},indent=1))
print("\nSTATIC unfalsifiable:")
for a,b in lits: print("  [%s] %s"%(a," ".join(b.split())[:110]))
print("\nfirst 40 stored numbers the gate cannot see:")
for u in unseen[:40]: print("   %-18s %-60s %s"%(u[0],u[1],u[2]))
json.dump(OUT,open(os.path.join(HERE,"x4_gate_mutation_results.json"),"w"),indent=1,default=str)

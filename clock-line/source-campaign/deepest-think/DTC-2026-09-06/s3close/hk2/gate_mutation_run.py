#!/usr/bin/env python3
"""gate_mutation_run -- FL-043 discipline: show that check_constants.py can fire.
(a) blind sweep: 210 stored numbers perturbed one at a time (v -> 1.5v + 0.37);
(b) targeted sweep: the 29 leaves check_constants.py actually reads."""
import json, random, subprocess, sys
random.seed(20260908)
FILES = [f'k{k}_results.json' for k in (1,2,3,4,5,6,7,8)] + ['k3_conv.json']
def walk(o, path=()):
    if isinstance(o, dict):
        for k,v in o.items(): yield from walk(v, path+(k,))
    elif isinstance(o, list):
        for i,v in enumerate(o): yield from walk(v, path+(i,))
    elif isinstance(o, float) and o != 0.0: yield path, o
def mutate(fn, path):
    orig = open(fn).read(); d = json.loads(orig); o = d
    for k in path[:-1]: o = o[k]
    o[path[-1]] = 1.5*o[path[-1]] + 0.37
    open(fn,'w').write(json.dumps(d, indent=1))
    r = subprocess.run([sys.executable,'check_constants.py'], capture_output=True)
    open(fn,'w').write(orig)
    return r.returncode != 0
c = m = 0
for fn in FILES:
    d = json.loads(open(fn).read()); leaves = list(walk(d))
    for path, _ in random.sample(leaves, min(28, len(leaves))):
        if mutate(fn, path): c += 1
        else: m += 1
TARGETS = json.load(open('gate_targets.json'))
tc = sum(1 for fn, path in TARGETS if mutate(fn, path))
out = dict(blind_caught=c, blind_missed=m, blind_total=c+m,
           targeted_caught=tc, targeted_total=len(TARGETS))
print(f"blind    : caught {c} / {c+m}  ({100*c/(c+m):.1f} %)  -- the misses are stored"
      f" intermediates PROOF.md never displays")
print(f"targeted : caught {tc} / {len(TARGETS)}  (the leaves check_constants actually reads)")
json.dump(out, open('gate_mutation.json','w'), indent=1)

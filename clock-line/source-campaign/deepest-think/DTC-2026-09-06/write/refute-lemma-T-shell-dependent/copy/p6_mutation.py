"""p6 - MEASURED mutation coverage of check_constants.py.

L-candidate from this campaign's ledger: "a gate must fail on the number it certifies."
So: mutate every stored number, one at a time, by v -> 1.5 v + 0.37, and count how many the
gate catches.  Coverage is reported, not asserted -- a number the gate is blind to is a number
PROOF.md may not quote as certified.
"""
import json, os, shutil, subprocess, sys, tempfile, copy

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = ['p1_results.json', 'p2_results.json', 'p3_results.json', 'p7_results.json',
         'p4_results.json', 'p5_results.json']
CHECKER = os.path.join(HERE, 'check_constants.py')

def walk(obj, path=()):
    if isinstance(obj, bool):
        yield path, obj
    elif isinstance(obj, (int, float)):
        yield path, obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk(v, path + (k,))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from walk(v, path + (i,))

def setin(obj, path, val):
    for p in path[:-1]:
        obj = obj[p]
    obj[path[-1]] = val

base = {f: json.load(open(os.path.join(HERE, f))) for f in FILES}
targets = []
for f in FILES:
    for path, val in walk(base[f]):
        targets.append((f, path, val))

caught, blind = 0, []
tmp = tempfile.mkdtemp(prefix='lemmaTp_mut_')
try:
    for f, path, val in targets:
        for g in FILES:
            shutil.copy(os.path.join(HERE, g), os.path.join(tmp, g))
        d = copy.deepcopy(base[f])
        new = (not val) if isinstance(val, bool) else 1.5*val + 0.37
        setin(d, path, new)
        json.dump(d, open(os.path.join(tmp, f), 'w'))
        r = subprocess.run([sys.executable, CHECKER, tmp], capture_output=True, text=True)
        if r.returncode != 0:
            caught += 1
        else:
            blind.append({'file': f, 'path': [str(p) for p in path], 'value': val})
finally:
    shutil.rmtree(tmp, ignore_errors=True)

out = dict(n_numbers=len(targets), n_caught=caught, n_blind=len(blind),
           coverage=caught/len(targets), blind=blind)
json.dump(out, open(os.path.join(HERE, 'p6_results.json'), 'w'), indent=1)
print("stored numbers mutated: %d   caught: %d   blind: %d   coverage: %.4f"
      % (out['n_numbers'], out['n_caught'], out['n_blind'], out['coverage']))
for b in blind:
    print("   BLIND  %s : %s = %r" % (b['file'], '/'.join(b['path']), b['value']))

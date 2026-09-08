#!/usr/bin/env python3
"""
u7 -- blind mutation test of the gate (FL-043 / L-89 discipline).

Every numeric leaf of every results JSON is perturbed by a relative 1e-3 (or an
absolute 1e-6 when the leaf is 0), one at a time, and `check_constants.py` is run
against the mutated copy.  A leaf that survives is a leaf the gate cannot see; the
log names them all, so that "the gate has teeth" is a measurement rather than a claim.
"""
import json, os, shutil, subprocess, sys, tempfile

FILES = ['u1_results.json', 'u2_results.json', 'u3_results.json',
         'u4_results.json', 'u5_results.json', 'u6_results.json']

def leaves(obj, path=()):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from leaves(v, path + (k,))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from leaves(v, path + (i,))
    elif isinstance(obj, bool):
        yield path, obj
    elif isinstance(obj, (int, float)):
        yield path, obj

def setpath(obj, path, val):
    for k in path[:-1]:
        obj = obj[k]
    obj[path[-1]] = val

def run(cwd):
    r = subprocess.run([sys.executable, 'check_constants.py'], cwd=cwd,
                       capture_output=True, text=True)
    return r.returncode

here = os.path.dirname(os.path.abspath(__file__))
caught = 0; blind = []
tmp = tempfile.mkdtemp()
for f in FILES + ['check_constants.py']:
    shutil.copy(os.path.join(here, f), tmp)
base = {f: json.load(open(os.path.join(here, f))) for f in FILES}
if run(tmp) != 0:
    print('BASELINE FAILS -- aborting'); sys.exit(2)

total = 0
for f in FILES:
    for path, val in list(leaves(base[f])):
        total += 1
        obj = json.loads(json.dumps(base[f]))
        if isinstance(val, bool):
            new = (not val)
        elif val == 0:
            new = 1e-6
        else:
            new = val*(1 + 1e-3)
        setpath(obj, path, new)
        with open(os.path.join(tmp, f), 'w') as fh:
            json.dump(obj, fh)
        rc = run(tmp)
        if rc != 0:
            caught += 1
        else:
            blind.append(f + ':' + '/'.join(map(str, path)))
        with open(os.path.join(tmp, f), 'w') as fh:
            json.dump(base[f], fh)

with open(os.path.join(here, 'u7_mutate_log.txt'), 'w') as fh:
    fh.write(f'leaves mutated: {total}\ncaught: {caught}\nblind: {len(blind)}\n\n')
    fh.write('BLIND LEAVES (mutating these does not fail check_constants.py):\n')
    for b in blind:
        fh.write('  ' + b + '\n')
print(f'{total} leaves, {caught} caught, {len(blind)} blind '
      f'({100.0*caught/max(total,1):.1f} % caught)')
json.dump(dict(total=total, caught=caught, blind=len(blind), blind_leaves=blind),
          open(os.path.join(here, 'u7_results.json'), 'w'), indent=1)
shutil.rmtree(tmp)

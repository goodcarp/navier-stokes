"""r4 - enumerate ALL blind mutation leaves of the target seat's gate.
The seat's gate prints only misses[:40]; PROOF.md section 5 characterises all 200 from that
sample.  This script lists every one and tests the characterisation."""
import sys, json, copy, collections
sys.path.insert(0, '../copy')
import os
os.chdir('../copy')
import importlib.util
spec = importlib.util.spec_from_file_location('cc', 'check_constants.py')
# check_constants runs its gate at import; capture by exec with argv stripped
import io, contextlib
src = open('check_constants.py').read().split("n, fails = run(")[0]
ns = {'__name__': 'ccmod'}
exec(compile(src, 'check_constants.py', 'exec'), ns)
run, leaves, setpath, BASE = ns['run'], ns['leaves'], ns['setpath'], ns['BASE']

tot = 0; misses = []
for f, obj in BASE.items():
    for path, val in list(leaves(obj)):
        J2 = copy.deepcopy(BASE)
        setpath(J2[f], path, 1.5 * val + 0.37)
        tot += 1
        try:
            _, fl = run(J2)
        except Exception:
            fl = ['exception']
        if not fl:
            misses.append((f, list(path), val))
print('total leaves', tot, 'blind', len(misses))
c = collections.Counter((m[0], str(m[1][-1])) for m in misses)
for k, v in sorted(c.items(), key=lambda kv: -kv[1]):
    print('%4d  file g%s  leaf key %r' % (v, k[0], k[1]))
claimed = {'secs', 'S_total_no_resid', 'partial_vs_cesaro_adev', 'resid_p',
           'S_lo', 'jump_tail_rem_bound'}
un = [m for m in misses if str(m[1][-1]) not in claimed]
print('\nblind leaves NOT covered by PROOF.md section 5 characterisation:', len(un))
for m in un[:80]:
    print('   g%s %s = %r' % (m[0], m[1], m[2]))
os.chdir('../scripts')
json.dump({'total_leaves': tot, 'blind': len(misses),
           'blind_by_key': {'g%s:%s' % k: v for k, v in c.items()},
           'uncharacterised': [[m[0], m[1], m[2]] for m in un]},
          open('r4_results.json', 'w'), indent=1)

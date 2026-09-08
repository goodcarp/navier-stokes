"""p8 - audit PROOF.md against the JSONs: every decimal number in the prose must be traceable.

A number in PROOF.md is TRACED if some stored number in p1..p7_results.json, or one of a short
list of exact mathematical constants derived from them, agrees with it to the precision it is
printed at.  Numbers that are settings written to their exact value (grid sizes, delta in degrees,
mu amplitudes, phi_c powers of ten, L values, sample counts) are traced through the JSONs too,
since those JSONs store them.

Reported, not asserted -- an untraced number is one PROOF.md must not claim as computed.
"""
import json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda n: json.load(open(os.path.join(HERE, n)))
pool = []

def harvest(o):
    if isinstance(o, bool):
        return
    if isinstance(o, (int, float)):
        pool.append(float(o)); return
    if isinstance(o, str):
        try:
            pool.append(float(o))
        except ValueError:
            pass
        return
    if isinstance(o, dict):
        for k, v in o.items():
            try:
                pool.append(float(k))
            except ValueError:
                pass
            harvest(v)
        return
    if isinstance(o, list):
        for v in o:
            harvest(v)

for f in ('p1_results.json', 'p2_results.json', 'p3_results.json', 'p4_results.json',
          'p5_results.json', 'p6_results.json', 'p7_results.json'):
    harvest(J(f))

# closed-form constants and simple functions of stored numbers that the prose also quotes
p3 = J('p3_results.json'); p4 = J('p4_results.json')
extra = [math.pi, math.pi**3, 3/(8*math.pi**2), 3/(2*math.pi**2),
         15*math.pi/8, 15*math.pi/4, 3*math.pi/4, 3*math.pi/2, 3*math.pi/8,
         math.sqrt(1.5), 2*(math.sqrt(1.5)-1), 4*(math.sqrt(1.5)-1),
         2*(1-math.sqrt(2/3)), 2*math.log(1.5), 1/15, -1/15, 2/3, 3/2, 1/2, 1/3, 1/4,
         p3['C_rel']/p3['shear_offset_rel_times_L']]
for r in p3['L_table']:
    extra += [r['bound_over_a_true_refuter']*r['L'], r['a_ref_derived'], r['L']]
for r in p4['lambda_rows'] + p4['phi_rows'] + p4['taper_rows']:
    extra += [abs(r['delta']), r['bound'], r['slack'], r['mu'], r['muJ'],
              r.get('muJ_vs_JLambda', 0.0)]
    if 'delta_Lambda' in r:
        extra += [abs(r['delta_Lambda']), r['bound_Lambda'], r['slack_Lambda']]
    extra += [r['muJ']*r['L']]
# three cross-instrument quantities the prose quotes, each a simple function of stored numbers
p5 = J('p5_results.json')
dK1 = [r['delta'] for r in p5['K1_rows']]
extra += [(max(dK1) - min(dK1))/min(dK1),
          abs(dK1[0] - p3['shear_offset_over_M'])/p3['shear_offset_over_M'],
          p5['K2_rows'][0]['mu_single'] - 0.385329]
pool += extra
pool = [p for p in pool if math.isfinite(p)]
pool += [abs(p) for p in pool]   # the prose prints signed z-scores with a unicode minus

txt = open(os.path.join(HERE, 'PROOF.md')).read()
# strip fenced code blocks and the file table (they quote formulas and filenames)
body = re.sub(r'```.*?```', ' ', txt, flags=re.S)
# numbers quoted verbatim FROM OTHER SEATS (at source) -- not computed here, listed explicitly
FOREIGN = {'0.3187', '0.7047', '0.3958', '0.4295', '0.385329', '0.3351', '0.3488', '0.9855',
           '0.4997212305210886', '0.4997212305210992', '1.4980', '1.4914', '1.4736', '1.4444',
           '1.3597', '1.2584', '1.0674'}
body = re.sub(r'\((\d\.\d+)\)', ' EQREF ', body)      # equation labels (2.3), (4.1), ...
nums = re.findall(r'(?<![\w.])(\d+\.\d+(?:e[-+]?\d+)?)', body)
untraced = []
for tok in sorted(set(nums)):
    if tok in FOREIGN:
        continue
    v = float(tok)
    mant = tok.split('e')[0]
    dec = len(mant.split('.')[1])
    scale = 10**int(tok.split('e')[1]) if 'e' in tok else 1
    tol = 0.5*10**(-dec)*scale
    if not any(abs(v - q) <= tol*1.0000001 for q in pool):
        untraced.append(tok)

foreign = sorted(FOREIGN & set(nums))
out = dict(n_distinct_numbers=len(set(nums)), n_untraced=len(untraced), untraced=untraced,
           n_foreign_quoted=len(foreign), foreign_quoted=foreign, pool_size=len(pool))
json.dump(out, open(os.path.join(HERE, 'p8_results.json'), 'w'), indent=1)
print("PROOF.md distinct decimal numbers: %d   untraced: %d   (pool %d)"
      % (out['n_distinct_numbers'], out['n_untraced'], out['pool_size']))
print("quoted verbatim from other seats (at source, not computed here): %d" % len(foreign))
for t in untraced:
    print("   UNTRACED", t)

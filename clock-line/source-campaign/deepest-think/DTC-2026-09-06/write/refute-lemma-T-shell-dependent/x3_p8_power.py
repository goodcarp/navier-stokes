"""x3 - FL-043 test of the target seat's DOCUMENT AUDIT (p8_doc_audit.py).
p8 reports '295 distinct decimal numbers, 0 untraced'.  A gate must fail on the number it
certifies, so: corrupt ONE decimal number in PROOF.md at a time and ask whether p8 notices.
Detection rate = p8's real power.  (p8 is run in a scratch copy; the seat's tree is untouched.)
"""
import json, os, re, shutil, subprocess, sys, tempfile, random
SRC = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/write/lemma-T-shell-dependent"
txt = open(os.path.join(SRC,'PROOF.md')).read()
body = re.sub(r'```.*?```', lambda m: ' '*len(m.group(0)), txt, flags=re.S)
body = re.sub(r'\((\d\.\d+)\)', lambda m: ' '*len(m.group(0)), body)
toks = [(m.start(), m.group(1)) for m in re.finditer(r'(?<![\w.])(\d+\.\d+(?:e[-+]?\d+)?)', body)]
seen, cand = set(), []
for pos, t in toks:
    if t in seen: continue
    seen.add(t); cand.append((pos, t))
rng = random.Random(20260906)
sample = cand if len(cand) <= 300 else rng.sample(cand, 300)
tmp = tempfile.mkdtemp(prefix='p8_power_')
for f in os.listdir(SRC):
    if f.endswith('.json') or f.endswith('.py') or f == 'PROOF.md':
        shutil.copy(os.path.join(SRC,f), os.path.join(tmp,f))
def perturb(tok):
    """bump the last printed digit by 3 (mod 10) -- a change strictly bigger than the
    tolerance p8 itself uses (half a unit in the last place)."""
    mant = tok.split('e')[0]; exp = tok[len(mant):]
    last = mant[-1]
    return mant[:-1] + str((int(last)+3) % 10) + exp
caught, missed = 0, []
for pos, tok in sample:
    new = perturb(tok)
    if new == tok: continue
    mutated = txt[:pos] + new + txt[pos+len(tok):]
    open(os.path.join(tmp,'PROOF.md'),'w').write(mutated)
    r = subprocess.run([sys.executable, os.path.join(tmp,'p8_doc_audit.py')],
                       capture_output=True, text=True, cwd=tmp)
    res = json.load(open(os.path.join(tmp,'p8_results.json')))
    if res['n_untraced'] >= 1 and new in res['untraced']:
        caught += 1
    else:
        missed.append({'token': tok, 'mutated': new})
shutil.rmtree(tmp, ignore_errors=True)
out = dict(n_tested=len(sample), n_caught=caught, n_missed=len(missed),
           detection_rate=caught/len(sample), missed=missed[:60])
json.dump(out, open('x3_results.json','w'), indent=1)
print("p8 document-audit power: %d numbers perturbed in the last printed digit, %d caught, "
      "%d missed -> detection rate %.4f" % (out['n_tested'], caught, len(missed), out['detection_rate']))
print("a sample of numbers p8 does NOT notice being wrong:")
for m in missed[:25]: print("   %-24s -> %s   (still 'traced')" % (m['token'], m['mutated']))

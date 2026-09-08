#!/usr/bin/env python3
"""Descriptive follow-up (NO pre-registered threshold attached to anything here).
1. the H_oct predictions, whose constants were also frozen in PREREG.md
2. a combined dataset: my 9(+) corner runs + the parent seat's 13 runs, read from their JSON
3. which single abscissa explains M T_d
4. how far out of reach the pre-registered floor test actually was
"""
import os, json, math, hashlib, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/sharp/viscous-numerics"
C, D = 1.8659, 5.0730
KBAR, SMIN = 8.68, 1.7432907101226691

def mine():
    rs = []
    for f in sorted(os.listdir(HERE)):
        if f.startswith("results_") and f.endswith(".json"):
            rs += [dict(r, src="corner") for r in json.load(open(os.path.join(HERE, f)))["runs"]]
    return [r for r in rs if abs(r["h"] - 0.125) < 1e-9 and r["T32"]]

def theirs():
    out = []
    for t in ("main", "n7", "resweep"):
        p = os.path.join(SRC, "results_%s.json" % t)
        if os.path.exists(p):
            out += [dict(r, src="parent:" + t) for r in json.load(open(p))["runs"]
                    if r["T32"] and abs(r["h"] - 0.125) < 1e-9 and r["kind"] == "shell"]
    return out

M, P = mine(), theirs()
# dedupe: my re-runs of (N,Re0) that the parent also ran -- keep mine, note agreement
kk = lambda r: (r["N"], round(r["Re0"], 6))
dup = [(kk(a), a["T32"], b["T32"]) for a in M for b in P if kk(a) == kk(b)]
print("REPRODUCTION CHECK -- runs I re-ran that the parent seat also ran:")
for k, a, b in dup:
    print("   N=%d Re0=%-6g  mine %.10f  parent %.10f   rel diff %.2e" % (k[0], k[1], a, b, abs(a - b) / b))
ALL = M + [b for b in P if kk(b) not in {kk(a) for a in M}]
print("   combined dataset: %d runs (%d mine, %d parent-only)" % (len(ALL), len(M), len(ALL) - len(M)))

print("\n1. H_oct  (PRE-REGISTERED constants: s* = max(%.4f, %.2f sqrt(nu/M)), Q' = 0.98 pinned / 1.10 viscous)"
      % (SMIN, KBAR))
def hoct(N, nu):
    s = max(SMIN, KBAR * math.sqrt(nu))
    return (0.98 if s <= SMIN * 1.0000001 else 1.10) / math.log(2.0**N / s), s
rr = []
for r in sorted(M, key=lambda r: (r["N"], r["Re0"])):
    p, sp = hoct(r["N"], r["nu"])
    rr.append(100 * (r["T32"] - p) / p)
    print("   N=%d Re0=%-6g  meas %.4f  H_oct %.4f  resid %+6.2f%%   s* pred %.3f  meas %.3f"
          % (r["N"], r["Re0"], r["T32"], p, rr[-1], sp, r["s_star_at_T32"]))
print("   H_oct: max |resid| = %.2f%%,  rms = %.2f%%" % (max(abs(x) for x in rr), float(np.sqrt(np.mean(np.square(rr))))))
hl = [100 * (r["T32"] - C / (r["logReE"] - D)) / (C / (r["logReE"] - D)) for r in M]
print("   H_log (same runs, same frozen-constant discipline): max |resid| = %.2f%%,  rms = %.2f%%"
      % (max(abs(x) for x in hl), float(np.sqrt(np.mean(np.square(hl))))))

print("\n2. IS M T_d A FUNCTION OF log Re_E ALONE?  (no fit involved -- pure iso-Re_E spread)")
groups = {}
for r in ALL:
    groups.setdefault(round(r["logReE"], 3), []).append(r)
worst = 0.0
for lre, gs in sorted(groups.items()):
    if len(gs) < 2: continue
    t = [g["T32"] for g in gs]
    sp = (max(t) - min(t)) / np.mean(t)
    worst = max(worst, sp)
    print("   log Re_E = %7.3f : %s   spread %.2f%%"
          % (lre, "  ".join("(N=%d,Re0=%g,log(R/rho0)=%.3f) %.4f" % (g["N"], g["Re0"], g["N"] * math.log(2), g["T32"]) for g in gs), 100 * sp))
print("   worst iso-Re_E spread = %.2f%%  vs grid uncertainty 2.25%%  ->  %s"
      % (100 * worst, "M T_d is NOT a function of log Re_E alone" if worst > 3 * 0.0225 else "consistent"))

print("\n3. WHICH SINGLE ABSCISSA?  (free 1- or 2-parameter fits over all %d runs)" % len(ALL))
y = np.array([r["T32"] for r in ALL])
X = dict(logReE=np.array([r["logReE"] for r in ALL]),
         logR=np.array([math.log(r["R"]) for r in ALL]),
         logRs=np.array([math.log(r["R"] / r["s_star_at_T32"]) for r in ALL]),
         halfReE=np.array([0.5 * r["logReE"] for r in ALL]))
def fit_cd(xx):
    best = None
    for d in np.arange(-8.0, float(xx.min()) - 0.02, 0.002):
        c = float(np.sum(y / (xx - d)) / np.sum(1.0 / (xx - d) ** 2))
        rms = float(np.sqrt(np.mean((y - c / (xx - d)) ** 2)))
        if best is None or rms < best[2]: best = (c, d, rms)
    return best
for nm, xx in X.items():
    c, d, rms = fit_cd(xx)
    print("   M T_d = c/(%-8s - d)   c=%7.4f  d=%+7.4f  rms=%.5f  (%.2f%% of mean T_d)"
          % (nm, c, d, rms, 100 * rms / y.mean()))
c1 = float(np.sum(y / X["logRs"]) / np.sum(1 / X["logRs"] ** 2))
print("   M T_d = c/log(R/s*)  (one parameter, d=0)   c=%.4f  rms=%.5f (%.2f%%)"
      % (c1, float(np.sqrt(np.mean((y - c1 / X["logRs"]) ** 2))),
         100 * float(np.sqrt(np.mean((y - c1 / X["logRs"]) ** 2))) / y.mean()))
# mixed abscissa
best = None
for w in np.arange(0.0, 1.001, 0.01):
    xx = w * X["logReE"] + (1 - w) * X["logR"]
    c, d, rms = fit_cd(xx)
    if best is None or rms < best[3]: best = (w, c, d, rms)
print("   M T_d = c/(w logRe_E + (1-w) log(R/rho0) - d):  w=%.2f c=%.4f d=%+.4f rms=%.5f"
      % best)
print("   [w = 0 means the octave count alone; w = 1 means log Re_E alone]")

print("\n4. THE VISCOUS FLOOR the maximum self-selects")
for r in sorted(ALL, key=lambda r: (r["N"], r["Re0"])):
    s, sv = r["s_star_at_T32"], math.sqrt(r["nu"])
    print("   N=%d Re0=%-6g sqrt(nu/M)=%.4f  s*=%6.3f  s*/sqrt(nu/M)=%7.2f  %s"
          % (r["N"], r["Re0"], sv, s, s / sv, "VISCOUS-LIMITED" if s > SMIN * 1.05 else "pinned to the datum cut"))
vb = [r for r in ALL if r["s_star_at_T32"] > SMIN * 1.05]
kk2 = [r["s_star_at_T32"] / math.sqrt(r["nu"]) for r in vb]
print("   on the viscous branch (%d runs): s*/sqrt(nu/M) = %.2f +- %.2f"
      % (len(vb), float(np.mean(kk2)), float(np.std(kk2))))
Qv = [r["T32"] * math.log(r["R"] / r["s_star_at_T32"]) for r in vb]
Qp = [r["T32"] * math.log(r["R"] / r["s_star_at_T32"]) for r in ALL if r not in vb]
print("   Q' = M T_d log(R/s*)  : viscous branch %.4f +- %.4f ; pinned branch %.4f +- %.4f ; all %.4f +- %.4f"
      % (np.mean(Qv), np.std(Qv), np.mean(Qp), np.std(Qp),
         np.mean(Qv + Qp), np.std(Qv + Qp)))

print("\n5. HOW FAR OUT OF REACH WAS THE PRE-REGISTERED FLOOR TEST?")
for target in (0.02, 0.015):
    LR = c1 / target
    print("   to see M T_d = %.3f under M T_d = %.4f/log(R/s*) you need log(R/s*) = %.1f, i.e. R/s* = e^%.0f"
          % (target, c1, LR, LR))
    print("      and since s* >= ~%.1f sqrt(nu/M) on the viscous branch, log Re_E >= 2 log(R/s*) + O(1) = %.0f,"
          " i.e. Re_E >= 10^%.0f" % (np.mean(kk2), 2 * LR, 2 * LR / math.log(10)))
print("   largest Re_E anywhere in this seat or the parent seat: %.3g (log %.2f)"
      % (max(r["ReE"] for r in ALL), max(r["logReE"] for r in ALL)))
print("   => the KILL-FLOOR clause as pre-registered CANNOT fire at any reachable resolution. It has no power.")
print("\nSCRIPT-SHA256", hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest())

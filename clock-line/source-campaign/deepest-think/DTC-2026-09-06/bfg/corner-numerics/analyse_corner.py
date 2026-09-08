#!/usr/bin/env python3
"""Apply the PREREG.md verdicts verbatim to whatever corner runs have landed."""
import sys, os, json, math, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/sharp/viscous-numerics"
C, D = 1.8659, 5.0730                      # H_log, FIXED in PREREG, never refitted
A0_FLOOR = 0.02
KILL_FLOOR_AT = 0.015
SIG_FLOOR = 0.02                           # PREREG floor on sigma_rel

def loadmine():
    rs = []
    for f in sorted(os.listdir(HERE)):
        if f.startswith("results_") and f.endswith(".json"):
            rs += json.load(open(os.path.join(HERE, f)))["runs"]
    return rs

def loadparent(tag):
    p = os.path.join(SRC, "results_%s.json" % tag)
    return json.load(open(p))["runs"] if os.path.exists(p) else []

parent = {r["N"]: r for r in loadparent("main") + loadparent("n7")}
mine = loadmine()
prod = [r for r in mine if abs(r["h"] - 0.125) < 1e-9]         # h = rho0/8 production
conv = [r for r in mine if abs(r["h"] - 0.125) >= 1e-9]        # leg D

def hlog(lre): return C / (lre - D)

def key(r): return (r["N"], round(r["Re0"], 6), round(r["h"], 6))

print("=" * 112)
print("CORNER SWEEP -- datum A, h = rho0/8, box 1.5R, Tmax = 3.0.  M = rho0 = 1.")
print("%4s %2s %7s %8s %8s %9s %9s %8s %8s %8s %7s %6s %8s %8s %7s"
      % ("leg", "N", "Re0", "nu", "sqrtnu", "logRe_E", "M*T_d", "H_log", "resid%", "M*T2",
         "s*", "oct*", "MTd*lgRs", "MTd*lgRe", "sup|e|"))
rows = []
for r in sorted(prod, key=lambda r: (r["N"], r["Re0"])):
    T = r["T32"]
    if T is None:
        print("%4s %2d %7g  --- NO T32 within Tmax ---" % (r.get("label", ""), r["N"], r["Re0"]))
        continue
    lre = r["logReE"]; p = hlog(lre); s = r["s_star_at_T32"]
    row = dict(leg=r.get("label", ""), N=r["N"], Re0=r["Re0"], nu=r["nu"], logReE=lre,
               MTd=T, H_log=p, resid=100 * (T - p) / p, MT2=r["T2"], s_star=s,
               Q=T * math.log(r["R"] / s), P=T * lre, eta=r["etamax_ratio"],
               R=r["R"], a_star=r["a_star_at_T32"], steps=r["nsteps"], wall=r["wall_s"])
    rows.append(row)
    print("%4s %2d %7g %8.5g %8.4f %9.4f %8.4f %8.4f %+8.2f %8s %7.3f %6.2f %8.4f %8.4f %7.4f"
          % (row["leg"], r["N"], r["Re0"], r["nu"], math.sqrt(r["nu"]), lre, T, p, row["resid"],
             ("%.4f" % r["T2"]) if r["T2"] else "----", s, math.log2(s), row["Q"], row["P"],
             r["etamax_ratio"]))

# ---- checks -------------------------------------------------------------------
print("\nCHECK 1  E is nu-independent  =>  logRe_E(N,Re0) = logRe_E(N,100) + ln(Re0/100)")
worst = 0.0
for r in prod + conv:
    if r["N"] in parent and abs(r["h"] - 0.125) < 1e-9:
        pred = parent[r["N"]]["logReE"] + math.log(r["Re0"] / 100.0)
        worst = max(worst, abs(pred - r["logReE"]))
print("   max |measured - predicted| over all h=1/8 runs : %.3e   (PREREG demands <= 1e-9)" % worst)

print("\nCHECK 2  maximum principle sup|eta| non-increasing (PREREG voids a run above 1+4.3e-3)")
bad = [(key(r), r["etamax_ratio"]) for r in prod + conv if r["etamax_ratio"] > 1 + 4.3e-3]
print("   max sup|eta| ratio = %.5f ; runs voided : %s"
      % (max(r["etamax_ratio"] for r in prod + conv), bad if bad else "none"))

# ---- leg D : sigma_rel --------------------------------------------------------
print("\nLEG D  resolution check")
sig = SIG_FLOOR
pairs = []
for r in conv:
    base = [q for q in prod if q["N"] == r["N"] and abs(q["Re0"] - r["Re0"]) < 1e-9]
    if not base or r["T32"] is None or base[0]["T32"] is None: continue
    d = abs(r["T32"] - base[0]["T32"]) / base[0]["T32"]
    pairs.append((r["N"], r["Re0"], r["h"], r["T32"], base[0]["T32"], d))
    print("   N=%d Re0=%-6g  h=1/%-4.1f T32=%.5f   vs h=1/8 T32=%.5f   rel diff %.3f%%"
          % (r["N"], r["Re0"], 1 / r["h"], r["T32"], base[0]["T32"], 100 * d))
    sig = max(sig, d)
print("   sigma_rel = max(0.02 floor, measured) = %.4f   =>  tolerance 3*sigma_rel = %.4f (%.2f%%)"
      % (sig, 3 * sig, 300 * sig))

# ---- KILL-LOG ----------------------------------------------------------------
print("\n" + "=" * 112)
print("PRE-REGISTERED VERDICT 1 : KILL-LOG  (H_log : M T_d = %.4f/(log Re_E - %.4f), constants frozen)" % (C, D))
legA = [r for r in rows if r["leg"] == "A"]
k1 = [(r["Re0"], r["resid"]) for r in legA if abs(r["resid"]) > 300 * sig]
print("  K1  leg-A residuals vs H_log (tolerance +-%.2f%%):" % (300 * sig))
for r in legA:
    print("       Re0=%-6g  meas %.4f  H_log %.4f  resid %+7.2f%%   %s"
          % (r["Re0"], r["MTd"], r["H_log"], r["resid"], "EXCEEDS" if abs(r["resid"]) > 300 * sig else "ok"))
print("      K1 fires: %s" % (bool(k1)))

k2 = None
lo = [r for r in legA if r["Re0"] == 1.0]; hi = [r for r in legA if r["Re0"] == 400.0]
if lo and hi:
    meas = lo[0]["MTd"] / hi[0]["MTd"]
    pred = lo[0]["H_log"] / hi[0]["H_log"]
    rel = abs(meas - pred) / pred
    k2 = rel > 3 * sig
    print("  K2  ratio MTd(Re0=1)/MTd(Re0=400) : measured %.4f  H_log %.4f  rel diff %.2f%%  -> fires: %s"
          % (meas, pred, 100 * rel, k2))

print("  K3  iso-Re_E families (|delta log Re_E| <= 0.01) -- H_log demands equality:")
k3 = False
seen = []
for i, a in enumerate(rows):
    for b in rows[i + 1:]:
        if abs(a["logReE"] - b["logReE"]) <= 0.01:
            rel = abs(a["MTd"] - b["MTd"]) / (0.5 * (a["MTd"] + b["MTd"]))
            f = rel > 3 * sig
            k3 = k3 or f
            print("       logRe_E=%.4f : (N=%d,Re0=%g) MTd=%.4f  vs  (N=%d,Re0=%g) MTd=%.4f"
                  "   spread %.2f%%  -> fires: %s"
                  % (a["logReE"], a["N"], a["Re0"], a["MTd"], b["N"], b["Re0"], b["MTd"],
                     100 * rel, f))
KILLLOG = bool(k1) or bool(k2) or k3
print("  >>> KILL-LOG : %s" % ("FIRES" if KILLLOG else "does NOT fire (H_log survives this corner)"))

# ---- KILL-FLOOR ---------------------------------------------------------------
print("\nPRE-REGISTERED VERDICT 2 : KILL-FLOOR (A0 = %.3f ; kills if any converged run gives M T_d < %.3f)"
      % (A0_FLOOR, KILL_FLOOR_AT))
mn = min(rows, key=lambda r: r["MTd"]) if rows else None
if mn:
    print("  smallest M T_d anywhere in this seat : %.4f  at N=%d Re0=%g" % (mn["MTd"], mn["N"], mn["Re0"]))
    print("  >>> KILL-FLOOR : %s" % ("FIRES" if mn["MTd"] < KILL_FLOOR_AT else "does NOT fire"))

# ---- descriptive fits (NOT pre-registered thresholds) -------------------------
print("\nDESCRIPTIVE (no threshold attached): free two-parameter refits over the corner runs")
if len(rows) >= 4:
    x = np.array([r["logReE"] for r in rows]); y = np.array([r["MTd"] for r in rows])
    lo2 = np.array([math.log(r["R"]) for r in rows])
    def fit_cd(xx):
        best = None
        for d in np.arange(-6.0, min(xx) - 0.05, 0.005):
            c = float(np.sum(y / (xx - d)) / np.sum(1.0 / (xx - d) ** 2))
            rms = float(np.sqrt(np.mean((y - c / (xx - d)) ** 2)))
            if best is None or rms < best[2]: best = (c, d, rms)
        return best
    print("   c/(log Re_E - d)      : c=%.4f d=%.4f rms=%.5f" % fit_cd(x))
    print("   c/(log(R/rho0) - d)   : c=%.4f d=%.4f rms=%.5f" % fit_cd(lo2))
    ss = np.array([math.log(r["R"] / r["s_star"]) for r in rows])
    print("   c/log(R/s*)           : c=%.4f rms=%.5f"
          % (float(np.sum(y / ss) / np.sum(1 / ss ** 2)),
             float(np.sqrt(np.mean((y - (np.sum(y / ss) / np.sum(1 / ss ** 2)) / ss) ** 2)))))
    # two-variable: MTd = c/(A*logReE + B*log(R/rho0) - d)
    from itertools import product
    best = None
    for w in np.arange(0.0, 1.0001, 0.02):
        xx = w * x + (1 - w) * (lo2 / max(lo2) * max(x))  # not used; keep simple below
    # honest simple test: which single abscissa explains the corner better
    print("   (the two abscissae are no longer degenerate here: "
          "corr(log Re_E, log R) = %.4f over these runs)" % float(np.corrcoef(x, lo2)[0, 1]))

json.dump(dict(rows=rows, sigma_rel=sig, KILL_LOG=KILLLOG,
               KILL_FLOOR=bool(mn and mn["MTd"] < KILL_FLOOR_AT),
               K1=k1, K2=k2, K3=k3, conv_pairs=pairs),
          open(os.path.join(HERE, "verdicts.json"), "w"), indent=1, default=float)
print("\nSCRIPT-SHA256", hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest())

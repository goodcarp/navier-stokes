#!/usr/bin/env python3
"""Addendum: (a) clean the viscous/pinned classification, (b) the surviving shape --
log Re_E as a CEILING on the octave count rather than the clock's own abscissa."""
import os, json, math, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/sharp/viscous-numerics"
def mine():
    rs = []
    for f in sorted(os.listdir(HERE)):
        if f.startswith("results_") and f.endswith(".json"):
            rs += json.load(open(os.path.join(HERE, f)))["runs"]
    return [r for r in rs if abs(r["h"] - 0.125) < 1e-9 and r["T32"]]
def theirs():
    o = []
    for t in ("main", "n7", "resweep"):
        p = os.path.join(SRC, "results_%s.json" % t)
        if os.path.exists(p):
            o += [r for r in json.load(open(p))["runs"]
                  if r["T32"] and abs(r["h"] - 0.125) < 1e-9 and r["kind"] == "shell"]
    return o
M = mine(); K = {(r["N"], round(r["Re0"], 6)) for r in M}
ALL = M + [r for r in theirs() if (r["N"], round(r["Re0"], 6)) not in K]
ALL.sort(key=lambda r: (r["N"], r["Re0"]))

print("(a) CLASSIFICATION by whether the viscous floor sits above the datum's inner cut")
print("    (datum cut = the shell's inner edge; the discrete argmax sits at s = 1.743--1.865 rho0)")
vb = [r for r in ALL if r["Re0"] <= 25.0]
pb = [r for r in ALL if r["Re0"] > 25.0]
k = [r["s_star_at_T32"] / math.sqrt(r["nu"]) for r in vb]
print("    viscosity-limited (Re0 <= 25, sqrt(nu/M) >= 0.2 rho0): %d runs, s*/sqrt(nu/M) = %.2f +- %.2f, range [%.2f, %.2f]"
      % (len(vb), np.mean(k), np.std(k), min(k), max(k)))
print("    datum-limited     (Re0 >  25): %d runs, s* = %s (all pinned)"
      % (len(pb), sorted({round(r["s_star_at_T32"], 3) for r in pb})))
Qv = [r["T32"] * math.log(r["R"] / r["s_star_at_T32"]) for r in vb]
Qp = [r["T32"] * math.log(r["R"] / r["s_star_at_T32"]) for r in pb]
print("    Q' = M T_d log(R/s*) : viscous %.4f +- %.4f   datum-limited %.4f +- %.4f   all %.4f +- %.4f"
      % (np.mean(Qv), np.std(Qv), np.mean(Qp), np.std(Qp), np.mean(Qv + Qp), np.std(Qv + Qp)))
KB = float(np.mean(k))

print("\n(b) THE ASSIGNED SWEEP, stated as the pre-registration asked")
A = sorted([r for r in M if r["N"] == 6 and r["Re0"] in (1.0, 4.0, 16.0, 100.0, 400.0)], key=lambda r: r["Re0"])
print("    Re0 1 -> 400 at fixed N = 6 : log Re_E %.4f -> %.4f (a factor %.0f in Re_E)"
      % (A[0]["logReE"], A[-1]["logReE"], A[-1]["ReE"] / A[0]["ReE"]))
print("    M T_d %.4f -> %.4f, a factor %.4f.   H_log demands a factor %.4f."
      % (A[0]["T32"], A[-1]["T32"], A[0]["T32"] / A[-1]["T32"],
         (1.8659 / (A[0]["logReE"] - 5.0730)) / (1.8659 / (A[-1]["logReE"] - 5.0730))))
print("    For comparison, the parent seat's N-sweep covers log Re_E %.3f -> %.3f (a factor %.4f in log-range)"
      % (6.4616, 13.3977, 1.0))
print("    and moves M T_d 1.3446 -> 0.2272, a factor 5.918, over Delta log Re_E = 6.936.")
print("    The corner moves Delta log Re_E = %.3f and gets only a factor %.3f."
      % (A[-1]["logReE"] - A[0]["logReE"], A[0]["T32"] / A[-1]["T32"]))
print("    Same log-Re_E budget, different currency: octaves buy ~%.2fx more speed per unit log Re_E than viscosity."
      % ((math.log(5.918) / 6.936) / (math.log(A[0]["T32"] / A[-1]["T32"]) / (A[-1]["logReE"] - A[0]["logReE"]))))

print("\n(c) THE SURVIVING SHAPE: log Re_E is a CEILING on the octave count, not the clock's abscissa")
print("    s* >= %.2f sqrt(nu/M)  =>  log(R/s*) <= (1/2) log(M R^2/nu) - log %.2f"
      % (KB, KB))
print("    and M R^2/nu = Re_E / (ell/R)^2, with ell/R measured = %.5f"
      % np.mean([r["ell"] / r["R"] for r in ALL]))
off = math.log(KB) + math.log(np.mean([r["ell"] / r["R"] for r in ALL]))
print("    =>  log(R/s*) <= (1/2) log Re_E - %.4f  =: L_max(Re_E)" % off)
c1 = float(np.sum([r["T32"] for r in ALL] / np.array([math.log(r["R"] / r["s_star_at_T32"]) for r in ALL]))
           / np.sum(1.0 / np.array([math.log(r["R"] / r["s_star_at_T32"]) for r in ALL]) ** 2))
print("    with M T_d = %.4f / log(R/s*)  this gives the lower envelope" % c1)
print("        M T_d >= %.4f / ((1/2) log Re_E - %.4f) = %.4f / (log Re_E - %.4f)"
      % (c1, off, 2 * c1, 2 * off))
print("    check: does any run fall BELOW that envelope?")
bad = 0
for r in ALL:
    Lm = 0.5 * r["logReE"] - off
    env = c1 / Lm
    meas_L = math.log(r["R"] / r["s_star_at_T32"])
    flag = "BELOW ENVELOPE" if r["T32"] < env * (1 - 0.0675) else ""
    bad += bool(flag)
    print("      N=%d Re0=%-6g  log(R/s*)=%.4f  L_max=%.4f  M T_d=%.4f  envelope=%.4f  %s"
          % (r["N"], r["Re0"], meas_L, Lm, r["T32"], env, flag))
print("    runs below the envelope by more than 3 sigma_rel: %d / %d" % (bad, len(ALL)))
print("\nSCRIPT-SHA256", hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest())

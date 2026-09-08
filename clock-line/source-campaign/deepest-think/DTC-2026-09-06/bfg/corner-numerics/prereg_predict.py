#!/usr/bin/env python3
"""PRE-RUN: compute the pre-registered predictions of the two competing laws for
every run I am about to launch.  Written and executed BEFORE any production run.

Law H_log  (pure logarithmic clock, the fit the parent seat reported):
      M T_d = c / (log Re_E - d),   (c, d) = (1.8659, 5.0730)   FIXED, not refitted.
Law H_oct  (octave / geometric clock):
      M T_d = Q' / log(R / s*),  s* = max(s_min, kbar sqrt(nu/M)),
      kbar = 8.68, s_min = 1.743 (the datum's own inner cut, measured by the parent seat),
      Q' = 0.98 when s* pins to s_min, 1.10 +- 0.10 on the viscosity-limited branch.

log Re_E is NOT assumed: E is computed from the t=0 velocity field, which does not
depend on nu at all (build/solve never see nu), so log Re_E(N, Re0) =
log Re_E(N, 100) + ln(Re0/100) EXACTLY, with log Re_E(N,100) the parent seat's
measured values.  Those are read from their results_*.json, not typed.
"""
import json, math, os, hashlib
S = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/sharp/viscous-numerics"
L = lambda t: json.load(open(os.path.join(S, "results_%s.json" % t)))["runs"]
main = {r["N"]: r for r in L("main") + L("n7")}
C, D = 1.8659, 5.0730
KBAR, SMIN = 8.68, 1.7432907101226691

def logReE(N, Re0):
    return main[N]["logReE"] + math.log(Re0 / 100.0)

def hlog(N, Re0):
    return C / (logReE(N, Re0) - D)

def hoct(N, Re0):
    R = 1.0 * 2.0**N
    nu = 1.0 / Re0
    s = max(SMIN, KBAR * math.sqrt(nu))
    Qp = 0.98 if s <= SMIN * 1.0000001 else 1.10
    return Qp / math.log(R / s), s

print("parent-seat measured log Re_E at Re0=100 :",
      {N: round(main[N]["logReE"], 4) for N in sorted(main)})
print("parent-seat measured M*T32 at Re0=100    :",
      {N: (round(main[N]["T32"], 4) if main[N]["T32"] else None) for N in sorted(main)})
print()
print("H_log check against the six primary points it was fitted to:")
for N in sorted(main):
    if main[N]["T32"] is None: continue
    p = hlog(N, 100.0)
    print("   N=%d  logReE=%7.4f  meas %.4f  H_log %.4f  resid %+6.2f%%"
          % (N, logReE(N, 100.0), main[N]["T32"], p, 100 * (main[N]["T32"] - p) / p))
print()
print("PLAN AND PREDICTIONS (M = rho0 = 1, so t is M*t):")
hdr = "  leg  N  Re0        nu      log(R/rho0)  logRe_E   H_log MTd   H_oct MTd  (s*)"
print(hdr)
PLAN = [("A", 6, 1.0), ("A", 6, 4.0), ("A", 6, 16.0), ("A", 6, 100.0), ("A", 6, 400.0),
        ("B", 6, 25.0), ("B", 7, 6.25), ("B", 5, 100.0),
        ("C", 7, 100.0), ("E", 7, 400.0)]
rows = []
for leg, N, Re0 in PLAN:
    lo, ho = hlog(N, Re0), hoct(N, Re0)
    print("   %s  %d %7.4g %9.5g %10.4f %9.4f %10.4f %10.4f  (%.3f)"
          % (leg, N, Re0, 1.0 / Re0, N * math.log(2), logReE(N, Re0), lo, ho[0], ho[1]))
    rows.append(dict(leg=leg, N=N, Re0=Re0, nu=1.0/Re0, logRrho0=N*math.log(2),
                     logReE=logReE(N, Re0), H_log=lo, H_oct=ho[0], s_star_pred=ho[1]))
print()
print("ISO-Re_E FAMILIES (the degeneracy breakers) -- H_log says the members are EQUAL:")
for tag, mem in (("logRe_E = %.4f" % logReE(5, 100.0), [(5, 100.0), (6, 25.0), (7, 6.25)]),
                 ("logRe_E = %.4f" % logReE(7, 100.0), [(6, 400.0), (7, 100.0)])):
    print("   %s :" % tag)
    for N, Re0 in mem:
        lo, ho = hlog(N, Re0), hoct(N, Re0)
        print("      N=%d Re0=%-6g  log(R/rho0)=%.4f  H_log %.4f   H_oct %.4f (s*=%.3f)"
              % (N, Re0, N * math.log(2), lo, ho[0], ho[1]))
json.dump(rows, open("prereg_predictions.json", "w"), indent=1)
print()
print("SCRIPT-SHA256", hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest())

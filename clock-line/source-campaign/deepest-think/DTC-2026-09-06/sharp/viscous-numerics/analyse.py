#!/usr/bin/env python3
"""Read results_*.json and apply the PREREG gates verbatim."""
import sys, os, json, math, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))

def load(tag):
    fn = os.path.join(HERE, "results_%s.json" % tag)
    return json.load(open(fn))["runs"] if os.path.exists(fn) else []

def fmt(v, w, p, dash="----"):
    return ("%*.*f" % (w, p, v)) if v is not None else ("%*s" % (w, dash))

def table(runs, title):
    print("=" * 118); print(title)
    hdr = ("kind", "N", "h", "Re0", "grid", "R", "E", "l", "Re_E", "logReE",
           "T32*M", "T2*M", "P", "Q", "oct*", "a*", "sup|eta|", "lbl")
    print("%6s %2s %6s %5s %10s %6s %10s %7s %10s %7s %8s %8s %8s %8s %5s %7s %8s %s" % hdr)
    for r in runs:
        T32, T2, s = r["T32"], r["T2"], r["s_star_at_T32"]
        P = T32 * r["logReE"] if T32 else None
        Q = T32 * math.log(r["R"]) if T32 else None
        oc = math.log2(s) if s else None
        print("%6s %2d %6.4f %5g %4dx%-5d %6g %10.4g %7.3f %10.4g %7.3f %s %s %s %s %s %s %8.5f %s"
              % (r["kind"], r["N"], r["h"], r["Re0"], r["grid"][0], r["grid"][1], r["R"],
                 r["E"], r["ell"], r["ReE"], r["logReE"],
                 fmt(T32, 8, 4), fmt(T2, 8, 4), fmt(P, 8, 4), fmt(Q, 8, 4),
                 fmt(oc, 5, 2), fmt(r["a_star_at_T32"], 7, 4), r["etamax_ratio"], r.get("label", "")))

def gates(runs, tag):
    ok = [r for r in runs if r["T32"] is not None and r["N"] >= 2
          and "visc" not in r.get("label", "")]
    if len(ok) < 3:
        print("  [%s] too few points for the fit" % tag); return None
    x = np.log([r["N"] for r in ok]); y = np.log([r["T32"] for r in ok])
    p = np.polyfit(x, y, 1); beta = -p[0]
    Ps = [r["T32"] * r["logReE"] for r in ok]
    Qs = [r["T32"] * math.log(r["R"]) for r in ok]
    big = max(runs, key=lambda r: r["N"])
    small = [r for r in runs if r["N"] == 2]
    bigfail = big["T32"] is None and bool(small) and small[0]["T32"] is not None
    if beta <= 0.3 or bigfail:
        v = "KILL"
    elif beta >= 0.7 and max(Ps) <= 8.0:
        v = "SURVIVE"
    elif beta >= 0.7:
        v = "SURVIVE-WITH-CAVEAT (beta gate met; P exceeds 8 at the smallest N only)"
    else:
        v = "INCONCLUSIVE"
    print("  [%s] fit log(T32*M) = %.4f - %.4f log N   over N = %s" % (tag, p[1], beta, [r["N"] for r in ok]))
    print("  [%s] P = T32*M*logRe_E     : %s   max %.3f" % (tag, ["%.3f" % q for q in Ps], max(Ps)))
    print("  [%s] Q = T32*M*log(R/rho0) : %s" % (tag, ["%.3f" % q for q in Qs]))
    Q2 = [r["T32"] * math.log(r["R"] / r["s_star_at_T32"]) for r in ok]
    print("  [%s] Q'= T32*M*log(R/s*)    : %s   (s* = radius carrying the max at T32)"
          % (tag, ["%.3f" % q for q in Q2]))
    ok2 = [r for r in ok if r["T2"] is not None]
    if ok2:
        print("  [%s] T2*M*log(R/s*) [s* from T32]  : %s   (frozen-strain value ln2/ln1.5 x Q')"
              % (tag, ["%.3f" % (r["T2"] * math.log(r["R"] / r["s_star_at_T32"])) for r in ok2]))
    A = [r["T32"] * r["a_star_at_T32"] for r in ok]
    print("  [%s] T32*M*a*   (frozen-strain value ln 1.5 = 0.4055) : %s"
          % (tag, ["%.3f" % q for q in A]))
    da = [(ok[i+1]["a_star_at_T32"] - ok[i]["a_star_at_T32"]) / (ok[i+1]["N"] - ok[i]["N"])
          for i in range(len(ok) - 1)]
    print("  [%s] d a* / d octave : %s   -> per e-fold %s"
          % (tag, ["%.4f" % q for q in da], ["%.4f" % (q / math.log(2)) for q in da]))
    print("  [%s] PREREG VERDICT: %s" % (tag, v))
    return dict(Qprime=Q2, T32a=A, da_per_octave=da, beta=float(beta), intercept=float(p[1]), P=Ps, Q=Qs, verdict=v,
                N=[r["N"] for r in ok])

summary = {}
for tag, title in (("main", "PRIMARY: datum A (smoothed bang-bang shell, N octaves), Re0=100, h=1/8, lam=1.5"),
                   ("n7", "PRIMARY extension N=7"),
                   ("rings", "VARIANT B, PREREG horizon T_max = max(0.6, 4/N) -- too short for datum B (DEV-2)"),
                   ("ringsT", "VARIANT: datum B (N discrete dyadic ring-pairs), horizon extended to T=3 (DEV-2)"),
                   ("conv", "GRID CONVERGENCE"),
                   ("resweep", "VISCOSITY SWEEP at fixed N=5 (R=32): which scale survives, and what log is left"),
                   ("controls", "CONTROLS / SENSITIVITY")):
    runs = load(tag)
    if runs:
        table(runs, title)
        if tag == "resweep":
            print("  Re0   nu      sqrt(nu/M)   s*      s*/sqrt(nu/M)   log(R/s*)  Q'=T32*log(R/s*)")
            for r in runs:
                if r["T32"] is None: continue
                sv = math.sqrt(r["nu"]); ss = r["s_star_at_T32"]
                print("  %6g %8.4f %10.4f %8.3f %12.2f %11.4f %10.4f"
                      % (r["Re0"], r["nu"], sv, ss, ss / sv, math.log(r["R"] / ss),
                         r["T32"] * math.log(r["R"] / ss)))
        if tag == "ringsT":
            summary[tag] = gates(runs, tag)
mainr = sorted(load("main") + load("n7"), key=lambda r: r["N"])
if mainr:
    table(mainr, "PRIMARY (datum A) full sweep")
    summary["main"] = gates(mainr, "main")
json.dump(summary, open(os.path.join(HERE, "gates.json"), "w"), indent=1, default=float)
print("=" * 118)
print("SCRIPT-SHA256", hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest())

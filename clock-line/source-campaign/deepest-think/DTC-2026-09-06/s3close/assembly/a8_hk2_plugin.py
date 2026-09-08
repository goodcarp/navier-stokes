"""
a8 -- (H-K2) plugged in.

`s3close/hk2/PROOF.md` (filed 2026-09-08, after a1-a7 were run) closes (H-K2) in a corrected
form.  Its final statement, quoted:

  PROVED    K2 <= 66.6622 * M/rho0   on the whole tube N_tau, lambda in [1,3/2], no log(R/rho0)
                                     -- for the datum (D-B) whose radial cutoff is tanh-mollified
                                        at w0 = 0.25 rho0
  MEASURED  K2  =  2.1903 * M/rho0 (tracked point), 3.0202 (global probe, binding at the taper)
  FALSE     for the sharp-radial-edge datum (D-A) at f = 0: K2 = +infinity;
            with an inset, K2 ~= 8 M/(f rho0).

So C_K is no longer symbolic.  This script re-runs the budget at
  C_K = 66.6622  (their proved bound)
  C_K =  3.0202  (their global measured value)
  C_K =  8/f     (the sharp-edge idealisation at inset f)
and reports the moved L_* and eps.
"""
import json, math, os, sys, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); sys.modules[name] = mod
    spec.loader.exec_module(mod); return mod
a4 = load("a4mod8", "a4_subwindows.py")
eps_N, CSTAR, COLUMNS = a4.eps_N, a4.CSTAR, a4.COLUMNS
assemble = a4.assemble
SHIFT = json.load(open("a2_results.json"))["logReE_minus_2L"]

HK2 = {"proved_bound": 66.6622, "measured_global": 3.0202, "measured_tracked": 2.1903,
       "sharp_edge_times_f": 8.0}
OUT = {"hk2_statement": HK2, "logReE_minus_2L": SHIFT, "c_star": CSTAR}

FS = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]

def eps_best(L, N, col, CK):
    best = float("inf")
    for f in FS:
        e, _ = eps_N(L, N, col, f=f, C_K=CK)
        if e < best:
            best = e
    return best

def eps_best_sharp(L, N, col):
    """C_K = 8/f: the sharp-radial-edge idealisation, where the inset buys down K2"""
    best = float("inf")
    for f in FS:
        e, _ = eps_N(L, N, col, f=f, C_K=8.0/f)
        if e < best:
            best = e
    return best

tab = {}
for name, fn in [("C_K=66.6622", lambda L, N, c: eps_best(L, N, c, HK2["proved_bound"])),
                 ("C_K=3.0202",  lambda L, N, c: eps_best(L, N, c, HK2["measured_global"])),
                 ("C_K=8/f",     eps_best_sharp),
                 ("C_K=1 (a2 baseline)", lambda L, N, c: eps_best(L, N, c, 1.0))]:
    for col in ["proved", "measured"]:
        for N in [1, 16]:
            row = {}
            for L in [160.0, 640.0, 2560.0]:
                e = fn(L, N, col)
                row["eps(L=%g)" % L] = None if not math.isfinite(e) else e
            for target in [0.5, 0.1]:
                lo, hi, best = 1.0, 1e60, None
                if fn(hi, N, col) <= target:
                    for _ in range(60):
                        mid = math.sqrt(lo*hi)
                        if fn(mid, N, col) <= target:
                            hi = mid
                        else:
                            lo = mid
                        if hi/lo < 1.0005:
                            break
                    best = hi
                row["L_star_eps<=%g" % target] = best
                row["log_Lambda_star_eps<=%g" % target] = (2*best + SHIFT) if best else None
            tab["%s | %s | N=%d" % (name, col, N)] = row
OUT["table"] = tab
json.dump(OUT, open("a8_results.json", "w"), indent=1, sort_keys=True, default=str)
print("%-26s %-9s %3s %12s %12s %12s %14s %14s" %
      ("C_K", "column", "N", "eps(160)", "eps(640)", "eps(2560)", "L*(eps<=.5)", "L*(eps<=.1)"))
for k, v in tab.items():
    name, col, N = [x.strip() for x in k.split("|")]
    def g(x):
        return float("nan") if v[x] is None else float(v[x])
    print("%-26s %-9s %3s %12.5g %12.5g %12.5g %14.6g %14.6g" %
          (name, col, N.replace("N=", ""), g("eps(L=160)"), g("eps(L=640)"), g("eps(L=2560)"),
           g("L_star_eps<=0.5"), g("L_star_eps<=0.1")))

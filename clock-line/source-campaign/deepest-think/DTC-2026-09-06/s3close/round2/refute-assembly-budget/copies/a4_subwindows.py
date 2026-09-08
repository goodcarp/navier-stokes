"""
a4 -- what a RESTART lemma would be worth.

The whole reach of the assembly is set by one number: the linearised feedback
exponent  p c, where

    p/(ML) = 1.5 + C'/L + (geom_hi)(2 kappa)(dEps_T'/dmu)(1/geom_lo)

is the rate at which the strain functional amplifies a relative perturbation of the
flow map, and c = M L tau is the window in theta.  Grownwall over [0,tau] costs
e^{pc}.  If the argument could be RESTARTED at N intermediate times -- i.e. if the
reference map Lambda were re-derived from the state at each restart, so the map
error resets to zero -- the cost would be N sub-windows of exponent pc/N each:

    amplification(1 window)  = (e^{pc}-1)/p
    amplification(N windows) = N (e^{pc/N}-1)/p

This script prices that.  It does NOT claim the restart is available: doing so needs
a version of Theorem A / Corollary A1 / Theorem Gamma / Lemma 1 for a field that is
NOT the exact strained plateau (BLOCK 2 of ASSEMBLY.md).
"""
import json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
src = open("a2_budget.py").read().split('if __name__ == "__main__":')[0]
G = {"__name__": "a4"}
exec(src, G)

assemble = G["assemble"]; COLUMNS = G["COLUMNS"]; KAPPA = G["KAPPA"]
bootstrap = G["bootstrap"]; RECORD = G["RECORD"]; C_a_proved = G["C_a_proved"]
PHI0 = G["PHI0"]; eps_Tprime = G["eps_Tprime"]; R_H = G["R_H"]; KAPPA_S = G["KAPPA_S"]
SENS = G["SENS_MEASURED_FACTOR"]

CSTAR = math.log(1.5)/KAPPA

def mu_after(L, c, Cprime, C_a, sens, N):
    """majorant for mu at the end of the window when the map error is reset N-1 times"""
    bs = bootstrap(L, c/N, Cprime, C_a, sens=sens)
    return bs

def eps_N(L, N, column, f=1.0, C_K=1.0):
    """assemble with N sub-windows: the map error resets, the strain deficits add."""
    cp_key, ca_kind, use_V1, sens = COLUMNS[column]
    Cprime = RECORD[cp_key]
    C_a = C_a_proved(math.sin(PHI0)) if ca_kind == "far_near" else RECORD["material_offset_M"]
    bs = bootstrap(L, CSTAR/N, Cprime, C_a, sens=sens)
    if not bs["closed"]:
        return float("inf"), bs
    # the viscous budget is over the WHOLE window (viscosity does not reset)
    r_full = assemble(L, CSTAR, column, C_K=C_K, f=f)
    ell = G["ell_loss"](f)
    eps_a = ell/L + bs["eps_Tprime"] + 1.5*C_a/(KAPPA*L)
    eps_v = r_full["eps_v"]
    eps_delta = 1.0 - 2.0*KAPPA
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0:
        return float("inf"), bs
    eps = (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-eps_delta)) - 1.0
    return eps, bs

if __name__ == "__main__":
    OUT = {"c_star": CSTAR, "note": "N = number of sub-windows a restart lemma would allow"}
    tab = {}
    for column in list(COLUMNS):
        for N in [1, 2, 4, 8, 16, 32, 64]:
            # L_* for eps <= 1/2 and <= 0.1, bisection in log L
            row = {}
            for target in [0.5, 0.1]:
                lo, hi = 1.0, 1e60
                best = None
                if eps_N(hi, N, column)[0] <= target:
                    for _ in range(60):
                        mid = math.sqrt(lo*hi)
                        e = min(eps_N(mid, N, column, f=ff)[0] for ff in [0.25, 1.0, 4.0])
                        if e <= target:
                            hi = mid
                        else:
                            lo = mid
                        if hi/lo < 1.0005:
                            break
                    best = hi
                row["L_star_eps<=%g" % target] = best
            bs = bootstrap(1e18, CSTAR/N, RECORD[COLUMNS[column][0]],
                           C_a_proved(math.sin(PHI0)) if COLUMNS[column][1] == "far_near"
                           else RECORD["material_offset_M"], sens=COLUMNS[column][3])
            row["pc_per_subwindow"] = bs["feedback_exponent_pc"]
            row["amplification_total"] = N*(math.exp(bs["feedback_exponent_pc"])-1.0)
            tab["%s_N=%d" % (column, N)] = row
    OUT["table"] = tab
    with open("a4_results.json", "w") as f:
        json.dump(OUT, f, indent=1, sort_keys=True, default=str)
    print("c* = %.7f" % CSTAR)
    print("%-20s %4s %10s %14s %14s %14s" % ("column", "N", "pc/N", "N(e^{pc/N}-1)",
                                             "L*(eps<=.5)", "L*(eps<=.1)"))
    for k, v in tab.items():
        col, N = k.rsplit("_N=", 1)
        print("%-20s %4s %10.4f %14.4g %14.6g %14.6g"
              % (col, N, v["pc_per_subwindow"], v["amplification_total"],
                 v["L_star_eps<=0.5"] if v["L_star_eps<=0.5"] else float('nan'),
                 v["L_star_eps<=0.1"] if v["L_star_eps<=0.1"] else float('nan')))

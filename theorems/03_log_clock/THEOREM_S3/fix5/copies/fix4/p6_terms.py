"""
p6 -- the term table at the L where the proved column is alive, and the two direct
      sensitivity spot checks (ELL_RAMP and C_K) that p4's L_* bisection cannot resolve.

p4's term table is computed at L = 1e3 .. 1e6, where the proved column has no self-consistent
window at all (its L_* is 1.42e6), so every proved row there reads eps = inf.  This file adds
L = 2e6 and L = 1e7 and prints the decomposition.

It also settles a null result honestly: p4's sensitivity block returns the SAME L_* for
ELL_RAMP = 0.25 and 0.5 and for C_K = 1 and 161.7735.  That is not a caching artefact -- both
perturbations move eps by less than the bisection's own resolution.  The direct evaluation
below shows the size of each move at a fixed (L, c).

Outputs -> p6_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
src = open(os.path.join(HERE, "p4_budget.py")).read().split('if __name__ == "__main__":')[0]
NS = {"__file__": os.path.join(HERE, "p4_budget.py"), "__name__": "p4mod"}
exec(compile(src, "p4_budget.py", "exec"), NS)
os.chdir(HERE)
Datum = NS["Datum"]; bind = NS["bind"]; assemble_c = NS["assemble_c"]
self_consistent = NS["self_consistent"]; FS_F4 = NS["FS_F4"]; T3 = NS["T3"]

if __name__ == "__main__":
    SG = float(sys.argv[1]); SHIFT = float(sys.argv[2])
    P1RES = json.load(open(os.path.join(HERE, "p1_results.json")))
    row = P1RES["rows"]["sigma=%g" % SG]
    D = Datum(SG, SHIFT)
    E0, G0 = row["E0"], row["Gfrak0"]
    CAP = 1.0 + row["eps_cap_certified"]
    bind(D, E0, G0)
    OUT = {"sigma": SG, "cap": CAP}
    tt = []
    for L in (2e6, 1e7):
        for col in ("proved", "proved_Ca_measured", "measured"):
            s = self_consistent(D, L, col, E0, G0, CAP, fs=FS_F4)
            if not s["self_consistent"]:
                tt.append({"L": L, "column": col, "eps": float("inf")})
                continue
            A = assemble_c(D, L, s["c"], col, E0, G0, f=s["best_f"])
            tt.append(A)
            print("L=%.0e %-20s f=%g c/c_*=%.6f mu=%.6g muJ=%.4g epsT=%.6g eps_v=%.3g "
                  "eps_a=%.6g eps=%.6f C_R=%.4f C''=%.4f r_h=%.6f"
                  % (L, col, A["f"], A["c"]/D.c_star, A["mu"], A["muJ"], A["eps_Tprime"],
                     A["eps_v"], A["eps_a"], A["eps"], A["C_R"], A["C_pp"], A["r_h"]),
                  flush=True)
    OUT["term_table"] = tt
    # direct sensitivity at a fixed (L, c) in the proved column
    L0 = 2e6
    c0 = self_consistent(D, L0, "proved", E0, G0, CAP, fs=FS_F4)["c"]
    sens = {}
    for tag, kw in [("ELL_RAMP=0.25", dict(ell_ramp=0.25)),
                    ("ELL_RAMP=0.5", dict(ell_ramp=0.5)),
                    ("C_K=1", dict(C_K=1.0)),
                    ("C_K=161.7735", dict(C_K=161.7735)),
                    ("C_K=40.5556 (p3, mollified, lam_cap)", dict(C_K=40.5556))]:
        A = assemble_c(D, L0, c0, "proved", E0, G0, f=4.0, **kw)
        sens[tag] = {"eps": A["eps"], "eps_a": A["eps_a"], "eps_v": A["eps_v"],
                     "eps_ell": A["eps_ell"]}
    sens["delta_eps_ELL_RAMP"] = sens["ELL_RAMP=0.5"]["eps"] - sens["ELL_RAMP=0.25"]["eps"]
    sens["delta_eps_C_K"] = sens["C_K=161.7735"]["eps"] - sens["C_K=1"]["eps"]
    sens["L_and_c"] = {"L": L0, "c": c0, "c_over_c_star": c0/D.c_star}
    OUT["direct_sensitivity"] = sens
    json.dump(OUT, open(os.path.join(HERE, "p6_results.json"), "w"), indent=1,
              sort_keys=True, default=str)
    print(json.dumps(OUT["direct_sensitivity"], indent=1, default=str))

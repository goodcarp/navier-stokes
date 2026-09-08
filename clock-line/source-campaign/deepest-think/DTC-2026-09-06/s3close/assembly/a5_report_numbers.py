"""a5 -- the numbers quoted in the report line of ASSEMBLY.md."""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
import importlib.util, sys
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); sys.modules[name] = mod
    spec.loader.exec_module(mod); return mod
a4 = load("a4mod", "a4_subwindows.py")
G = {"KAPPA": a4.KAPPA, "R_H": a4.R_H, "COLUMNS": a4.COLUMNS,
     "C_E": a4.G["C_E"], "S_REC": a4.G["S_REC"]}
eps_N = a4.eps_N; CSTAR = a4.CSTAR
A2 = json.load(open("a2_results.json")); A4 = json.load(open("a4_results.json"))
SHIFT = A2["logReE_minus_2L"]

OUT = {"c2": 4.0*math.log(1.5), "c_star_window": CSTAR,
       "logReE_minus_2L": SHIFT, "kappa_delta": G["KAPPA"], "r_h": G["R_H"],
       "C_E": G["C_E"], "s": G["S_REC"]}

# eps at L = 160 and 640, every column, N = 1, 4, 16
tab = {}
for col in list(G["COLUMNS"]):
    for N in [1, 4, 16]:
        for L in [160.0, 640.0]:
            e = min(eps_N(L, N, col, f=ff)[0] for ff in [0.25, 0.5, 1.0, 2.0, 4.0])
            tab["%s_N=%d_L=%g" % (col, N, L)] = (None if not np.isfinite(e) else e)
OUT["eps_at_160_640"] = tab

# Lambda_* from every L_*
lam = {}
for k, v in A4["table"].items():
    for t in ["L_star_eps<=0.5", "L_star_eps<=0.1"]:
        Ls = v[t]
        if Ls in (None, "None"):
            lam["%s_%s" % (k, t)] = None
        else:
            Lsf = float(Ls)
            lam["%s_%s" % (k, t)] = {"L_star": Lsf, "log_Lambda_star": 2*Lsf + SHIFT,
                                     "log10_Lambda_star": (2*Lsf + SHIFT)/math.log(10.0)}
OUT["Lambda_star"] = lam

# the assembled theorem's constant and the numerics ratio
A3 = json.load(open("a3_results.json"))
OUT["numerics"] = {"fit_c": A3["fit_c_over_logRsstar"],
                   "fit_rms_pct": A3["fit_rms_pct_of_mean_Td"],
                   "ratio_vs_bound_L": A3["ratio_vs_bound_L"],
                   "ratio_vs_bound_Leff": A3["ratio_vs_bound_Leff"],
                   "L_range": A3["L_range"]}
with open("a5_results.json", "w") as f:
    json.dump(OUT, f, indent=1, sort_keys=True, default=str)
print(json.dumps(OUT, indent=1, sort_keys=True, default=str))

"""
r4 -- the 19-run confrontation of ASSEMBLY.md section 3.5, from my own transcription.

Sources (read directly, transcribed here, never from a3_numerics.py):
  sharp/viscous-numerics/NOTE.md section 4 table  (N = 2,3,4 at Re0 = 100; T32*M and
        'octave of s*' 0.90 / 0.90 / 0.80, so s* = 2^oct)
  sharp/viscous-numerics/NOTE.md section 5 table  (N = 5, Re0 = 1,2,4,16,400,1600; s* given)
  bfg/corner-numerics/NOTE.md section 2 table     (10 rows; s* given)
Datum of every run: omega^theta = -M sin 2phi, whose strain constant is
kappa = -(3/4) int_0^pi w cos(phi) sin^2(phi) dphi / M = 2/5 (recomputed here in closed form
by sympy, independently of the far/near lemma's quadrature).

What the assembled theorem would say about these runs: nothing.  Its datum is the tapered
plateau (kappa_delta ~ 1/2, taper, equatorial mollification), its L_* is >= 4.5e14 (at c_*) and
its eps(L) is +infinity at every L <= 2560 in the proved column.  The 'eps = 0 bound'
log(3/2)/(kappa L) is a quantity the theorem never asserts: eps >= eps_delta > 0 for all L,
and eps = infinity at these L.  So 'runs above the eps = 0 bound' is not a violation of
anything; it is the statement that the true eps(L) at L ~ 1.4..4.9 is between +0.07 and +2.3.
"""
import json, math, os
import numpy as np
import sympy as sp

phi = sp.symbols('phi')
kappa_sin2 = sp.simplify(-sp.Rational(3, 4)*sp.integrate(-sp.sin(2*phi)*sp.cos(phi)*sp.sin(phi)**2, (phi, 0, sp.pi)))
kappa_plateau = sp.simplify(-sp.Rational(3, 4)*(sp.integrate(-sp.cos(phi)*sp.sin(phi)**2, (phi, 0, sp.pi/2))
                                                + sp.integrate(sp.cos(phi)*sp.sin(phi)**2, (phi, sp.pi/2, sp.pi))))
OUT = {"kappa_sin2phi_closed": str(kappa_sin2), "kappa_plateau_closed": str(kappa_plateau)}
kap = float(kappa_sin2)

RUNS = [
    # (src, N, Re0, M T_d, s*)
    ("vn4", 2, 100, 1.3446, 2**0.90), ("vn4", 3, 100, 0.6667, 2**0.90), ("vn4", 4, 100, 0.4477, 2**0.80),
    ("vn5", 5, 1, 0.9567, 9.119), ("vn5", 5, 2, 0.6859, 5.810), ("vn5", 5, 4, 0.5405, 4.169),
    ("vn5", 5, 16, 0.3911, 2.259), ("vn5", 5, 400, 0.3289, 1.743), ("vn5", 5, 1600, 0.3268, 1.743),
    ("cn2", 5, 100, 0.3379, 1.743), ("cn2", 6, 1, 0.5385, 7.871), ("cn2", 6, 4, 0.3812, 3.654),
    ("cn2", 6, 16, 0.3011, 2.259), ("cn2", 6, 25, 0.2896, 2.138), ("cn2", 6, 100, 0.2716, 1.743),
    ("cn2", 6, 400, 0.2659, 1.743), ("cn2", 7, 6.25, 0.2753, 2.895), ("cn2", 7, 100, 0.2272, 1.743),
    ("cn2", 7, 400, 0.2232, 1.743),
]
OUT["n_runs"] = len(RUNS)
rows = []
for src, N, Re0, MTd, ss in RUNS:
    L = N*math.log(2.0); Le = math.log(2.0**N/ss)
    b = math.log(1.5)/(kap*L); be = math.log(1.5)/(kap*Le)
    rows.append({"src": src, "N": N, "Re0": Re0, "MTd": MTd, "L": L, "Leff": Le,
                 "ratio_L": MTd/b, "ratio_Leff": MTd/be, "implied_eps_L": MTd/b-1.0, "Qprime": MTd*Le})
OUT["rows"] = rows
rL = np.array([r["ratio_L"] for r in rows]); rE = np.array([r["ratio_Leff"] for r in rows])
Q = np.array([r["Qprime"] for r in rows])
OUT["ratio_L"] = {"min": float(rL.min()), "max": float(rL.max()), "mean": float(rL.mean()), "n_above_1": int((rL > 1).sum())}
OUT["ratio_Leff"] = {"min": float(rE.min()), "max": float(rE.max()), "mean": float(rE.mean()), "n_above_1": int((rE > 1).sum())}
OUT["Qprime"] = {"mean": float(Q.mean()), "sd_ddof1": float(Q.std(ddof=1)), "sd_ddof0": float(Q.std(ddof=0))}
OUT["eps0_bound_constant_log1.5_over_kappa"] = math.log(1.5)/kap
x = np.array([1.0/r["Leff"] for r in rows]); y = np.array([r["MTd"] for r in rows])
c_fit = float(np.dot(x, y)/np.dot(x, x)); resid = y - c_fit*x
OUT["fit"] = {"c": c_fit, "rms_pct_of_mean_Td": float(np.sqrt((resid**2).mean())/y.mean()*100)}
OUT["L_range"] = [min(r["L"] for r in rows), max(r["L"] for r in rows)]
OUT["assembly_quoted"] = {"ratio_L": [1.068, 3.271, 1.500], "ratio_Leff": [0.938, 1.185, 1.013],
                          "n_above_L": 19, "n_above_Leff": 8, "fit_c": 1.05338, "rms_pct": 8.776,
                          "Qprime": [1.0272, 0.0762]}
# what the theorem asserts at these L: eps(L) from the assembly's own a2 grid, read-only
a2 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assembly", "a2_results.json")))
OUT["assembly_eps_proved_at_L10_cstar"] = [r["eps"] for r in a2["grid"] if r["column"] == "proved" and r["L"] == 10.0 and abs(r["c"]-0.8113826) < 1e-3]
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r4_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print(json.dumps({k: v for k, v in OUT.items() if k != "rows"}, indent=1, default=str))

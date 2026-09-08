"""
a3 -- confrontation of the assembled upper bound with the campaign's 19 viscous runs.

The 19 rows are TRANSCRIBED (they are measurements made by other seats); their
provenance is recorded per row.  Everything computed from them is computed here.

Sources:
  P1 = sharp/viscous-numerics/NOTE.md primary sweep table (Re0 = 100, h = rho0/8)
  P2 = sharp/viscous-numerics/NOTE.md Re0 sweep at N = 5
  C  = bfg/corner-numerics/NOTE.md sec.2 table (all h = rho0/8)
Datum for every run: omega^theta = -M sin(2 phi) on rho0 < |x| < R = 2^N rho0,
rho0 = M = 1;  T_d = T32 = first time sup|omega^theta| = 1.5 M.
"""
import json, math
import numpy as np
from scipy.integrate import quad

RUNS = [
    # src,  N, Re0,     M T_d,   s*        (s* in units of rho0)
    ("P1",  2, 100,    1.3446,  2**0.90),
    ("P1",  3, 100,    0.6667,  2**0.90),
    ("P1",  4, 100,    0.4477,  2**0.80),
    ("C",   5, 100,    0.3379,  1.743),
    ("C",   6, 1,      0.5385,  7.871),
    ("C",   6, 4,      0.3812,  3.654),
    ("C",   6, 16,     0.3011,  2.259),
    ("C",   6, 25,     0.2896,  2.138),
    ("C",   6, 100,    0.2716,  1.743),
    ("C",   6, 400,    0.2659,  1.743),
    ("C",   7, 6.25,   0.2753,  2.895),
    ("C",   7, 100,    0.2272,  1.743),
    ("C",   7, 400,    0.2232,  1.743),
    ("P2",  5, 1,      0.9567,  9.119),
    ("P2",  5, 2,      0.6859,  5.810),
    ("P2",  5, 4,      0.5405,  4.169),
    ("P2",  5, 16,     0.3911,  2.259),
    ("P2",  5, 400,    0.3289,  1.743),
    ("P2",  5, 1600,   0.3268,  1.743),
]

def kappa_of_profile(w, npts=200001):
    """kappa[w] = -(3/4) int_0^pi w(phi) cos phi sin^2 phi dphi / M  (far-near-kernel Consequence A)"""
    phi = np.linspace(0.0, math.pi, npts)
    return -0.75*np.trapz(w(phi)*np.cos(phi)*np.sin(phi)**2, phi)

if __name__ == "__main__":
    OUT = {}
    # kappa for the runs' datum, computed here
    kap_run = kappa_of_profile(lambda p: -np.sin(2*p))
    OUT["kappa_run_datum_minus_sin2phi"] = kap_run
    OUT["kappa_run_exact_2over5"] = 0.4
    OUT["kappa_run_relerr"] = abs(kap_run-0.4)/0.4
    OUT["sup_abs_sin2phi"] = 1.0

    rows = []
    for src, N, Re0, MTd, sstar in RUNS:
        L = N*math.log(2.0)
        Leff = math.log((2.0**N)/sstar)
        bound_L = math.log(1.5)/(kap_run*L)          # eps = 0 bound on the run's own datum
        bound_Leff = math.log(1.5)/(kap_run*Leff)
        rows.append({"src": src, "N": N, "Re0": Re0, "MTd": MTd, "s_star": sstar,
                     "L": L, "log_R_over_sstar": Leff,
                     "MTd_times_L": MTd*L, "MTd_times_Leff": MTd*Leff,
                     "bound_eps0_L": bound_L, "ratio_meas_over_bound_L": MTd/bound_L,
                     "bound_eps0_Leff": bound_Leff,
                     "ratio_meas_over_bound_Leff": MTd/bound_Leff})
    OUT["rows"] = rows

    q = np.array([r["MTd_times_Leff"] for r in rows])
    OUT["Qprime_mean"] = float(q.mean()); OUT["Qprime_std"] = float(q.std(ddof=1))
    # one-parameter fit  M T_d = c / log(R/s*)  (least squares in T_d)
    x = np.array([1.0/r["log_R_over_sstar"] for r in rows])
    y = np.array([r["MTd"] for r in rows])
    c_fit = float(np.dot(x, y)/np.dot(x, x))
    resid = y - c_fit*x
    OUT["fit_c_over_logRsstar"] = c_fit
    OUT["fit_rms_abs"] = float(np.sqrt((resid**2).mean()))
    OUT["fit_rms_pct_of_mean_Td"] = float(np.sqrt((resid**2).mean())/y.mean()*100)
    OUT["record_fit_c"] = 1.0538
    OUT["record_fit_rms_pct"] = 8.8

    rL = np.array([r["ratio_meas_over_bound_L"] for r in rows])
    rE = np.array([r["ratio_meas_over_bound_Leff"] for r in rows])
    OUT["ratio_vs_bound_L"] = {"min": float(rL.min()), "max": float(rL.max()),
                               "mean": float(rL.mean()),
                               "n_above_1": int((rL > 1).sum()), "n": len(rL)}
    OUT["ratio_vs_bound_Leff"] = {"min": float(rE.min()), "max": float(rE.max()),
                                  "mean": float(rE.mean()),
                                  "n_above_1": int((rE > 1).sum()), "n": len(rE)}
    # implied empirical eps at each L, from  M T_d = (log(3/2)/kappa)(1+eps)/L
    OUT["implied_eps_L"] = {("N%d_Re%g" % (r["N"], r["Re0"])): r["ratio_meas_over_bound_L"]-1.0
                            for r in rows}
    OUT["L_range"] = [min(r["L"] for r in rows), max(r["L"] for r in rows)]
    with open("a3_results.json", "w") as f:
        json.dump(OUT, f, indent=1, sort_keys=True)
    print("kappa(run datum) = %.10f  (exact 2/5, rel %.2e)" % (kap_run, OUT["kappa_run_relerr"]))
    print("Q' = M T_d log(R/s*) : mean %.4f  sd %.4f" % (OUT["Qprime_mean"], OUT["Qprime_std"]))
    print("one-parameter fit c = %.4f  rms %.2f%% of mean T_d  (record: 1.0538, 8.8%%)"
          % (c_fit, OUT["fit_rms_pct_of_mean_Td"]))
    print("L range of the runs: %.3f .. %.3f" % tuple(OUT["L_range"]))
    print("ratio measured/(eps=0 bound at L=log(R/rho0)):  min %.3f max %.3f mean %.3f ; %d/%d above 1"
          % (rL.min(), rL.max(), rL.mean(), (rL > 1).sum(), len(rL)))
    print("ratio measured/(eps=0 bound at L=log(R/s*)):    min %.3f max %.3f mean %.3f ; %d/%d above 1"
          % (rE.min(), rE.max(), rE.mean(), (rE > 1).sum(), len(rE)))

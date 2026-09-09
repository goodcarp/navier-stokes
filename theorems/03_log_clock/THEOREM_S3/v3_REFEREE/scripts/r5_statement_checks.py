"""
r5 -- three arithmetic checks on the displayed statement of THEOREM_S3_v3 sec.1.4.

(a) THE ENERGY-REYNOLDS DICTIONARY.  The theorem proves M T_d <= c_* (1 + eps)/L and then
    displays  T(Lambda) <= c_2 (1 + eps)/log Lambda  with c_2 = 2 c_* and
    log Lambda = 2 L + 2.9135781820.  Since c_*(1+eps)/L = c_2(1+eps)/(2L) and
    2 L = log Lambda - 2.9135781820 < log Lambda, the displayed inequality is STRONGER than
    the proved one by the factor log Lambda/(log Lambda - 2.9135781820).  Quantified here.

(b) eps(L) DECREASING TO THE CLOCK FLOOR.  The statement asserts it; sec.3 item 7 records that
    it is tabulated at five L, not proved.  Tabulated here at twenty-five L with fix5/x2's own
    instrument, and the limit compared with (1 - 2 kappa_delta)/(2 kappa_delta).

(c) THE CLOCK IDENTITIES.  c_* = log(3/2)/kappa_delta, c_2 = 2 c_*, lam_max = e^{3c/4}, and the
    floor, recomputed from kappa_delta.

Outputs -> r5_results.json
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(HERE, "rerun")
sys.path.insert(0, RERUN)
os.chdir(RERUN)

import numpy as np                                                            # noqa: E402
import x2_budget as X2                                                        # noqa: E402

P4 = X2.P4
SHIFT = 2.9135781820
SIGMA = 0.02
OUT = {}

# ---------------------------------------------------------------- (c) the clock identities
P1RES = json.load(open(os.path.join(RERUN, "copies", "fix4", "p1_results.json")))
ROW = P1RES["rows"]["sigma=0.02"]
kap = ROW["kappa_delta"]
c_star = math.log(1.5)/kap
c_2 = 2.0*math.log(1.5)/kap
floor = (1.0 - 2.0*kap)/(2.0*kap)
OUT["c_clock_identities"] = {
    "kappa_delta_from_p1": kap,
    "c_star_here": c_star, "c_star_in_the_theorem": 0.8244732108,
    "c_2_here": c_2, "c_2_in_the_theorem": 1.6489464217,
    "lam_max_at_c_star_here": math.exp(0.75*c_star),
    "lam_max_in_the_theorem": 1.8558724485,
    "lam_max_at_the_certified_cap_here": math.exp(0.75*c_star*1.1647266909),
    "lam_max_at_cap_in_the_theorem": 2.0548738638,
    "clock_floor_here": floor, "clock_floor_in_the_theorem": 1.6700567265e-02,
    "c_2_equals_2_c_star": abs(c_2 - 2*c_star) < 1e-15}

# ---------------------------------------------------------------- (a) the dictionary factor
rows = []
for Lst in [1424207.9721, 1554565.5693, 1761.8423, 2050.2770]:
    logL = 2*Lst + SHIFT
    fac = logL/(2*Lst)
    for eps in [0.1647266909, 0.1]:
        rows.append({"L_star": Lst, "log_Lambda_star": logL,
                     "factor_by_which_the_display_is_stronger": fac,
                     "eps_displayed": eps,
                     "eps_that_the_display_actually_needs": (1+eps)*fac - 1.0,
                     "increase_in_eps": (1+eps)*fac - 1.0 - eps})
OUT["a_energy_Reynolds_dictionary"] = {
    "proved": "M T_d <= c_*(1+eps)/L = c_2 (1+eps)/(2L) = c_2 (1+eps)/(log Lambda - 2.9135781820)",
    "displayed": "T(Lambda) <= c_2 (1+eps)/log Lambda",
    "displayed_is_stronger_since": "log Lambda > 2 L",
    "rows": rows}

# ---------------------------------------------------------------- (b) eps(L)
X6 = json.load(open(os.path.join(RERUN, "x6_results.json")))
E0c, G0c = X6["E_0"]["upper"], X6["Gfrak_0"]["upper"]
CAPc = 1.0 + X6["window_cap_certified_PROVED_pad"]["eps_cap"]
D = P4.Datum(SIGMA, SHIFT)
P4.bind(D, E0c, G0c)
import p1_profile as P1                                                       # noqa: E402
D._rh.clear()
D.lip_pad = X6["h_second_bound"]["proved"]*D.prof.h
D.rh_nsub = 20000
X1 = json.load(open(os.path.join(RERUN, "x1_results.json")))
for f in P4.FS_F4:
    X2.VT[f] = X1["repaired"]["f=%g" % f]["vartheta"]
X2.install(4.0*SIGMA, 2.0)

tab = []
prev = None
mono = True
for L in [2e6, 3e6, 5e6, 1e7, 3e7, 1e8, 1e9, 1e10, 1e12, 1e14, 1e16, 1e18, 1e20, 1e24, 1e30]:
    sc = P4.self_consistent(D, L, "proved", E0c, G0c, CAPc, fs=P4.FS_F4)
    e = (sc["c"]/D.c_star - 1.0) if sc["self_consistent"] else None
    tab.append({"L": L, "self_consistent": bool(sc["self_consistent"]), "eps": e,
                "best_f": sc.get("best_f")})
    if e is not None and prev is not None and e > prev + 1e-15:
        mono = False
    if e is not None:
        prev = e
    print("  L = %-10.3g  self-consistent %-5s  eps = %s" % (L, sc["self_consistent"], e),
          flush=True)
X2.restore()
OUT["b_eps_of_L"] = {"table": tab, "monotone_decreasing_over_the_tabulated_range": mono,
                     "clock_floor": floor,
                     "last_eps": tab[-1]["eps"],
                     "last_eps_over_floor": (tab[-1]["eps"]/floor) if tab[-1]["eps"] else None,
                     "note": "a tabulation, at 15 values of L; still not a proof of monotonicity"}

json.dump(OUT, open(os.path.join(HERE, "r5_results.json"), "w"), indent=1, default=str)
print(json.dumps(OUT["c_clock_identities"], indent=1))
print("dictionary factor at log Lambda_* = %.4f : %.10f, eps must rise from %.10f to %.10f"
      % (rows[0]["log_Lambda_star"], rows[0]["factor_by_which_the_display_is_stronger"],
         rows[0]["eps_displayed"], rows[0]["eps_that_the_display_actually_needs"]))
print("eps(L) monotone decreasing over the tabulated range:", mono,
      " last eps %s vs floor %.10e" % (tab[-1]["eps"], floor))
print("wrote r5_results.json")

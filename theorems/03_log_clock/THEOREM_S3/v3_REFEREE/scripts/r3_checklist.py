"""
r3_checklist -- the P1 checklist run against THEOREM_S3_v3.md and fix5, at the actual
parameters, by script.

P1 (SEAT_HEADER_2026-09-08, proposed 2026-09-08 evening):
  every constant recomputed for THIS datum or proved profile-independent (say which);
  every extremum analytic or an outward-rounded enclosure, never a grid value;
  every imported theorem's hypotheses listed and matched to the statement;
  the window of the estimates versus the interval of the conclusion;
  geometric containments (balls, cones, layers) checked by script at the actual parameters.

This file does the geometry, the quantifiers and the profile-independence claims.  The
extremum item is r2_theta_cert.py (vartheta) and the re-runs of x2/x6.

Outputs -> r3_results.json
"""
import json
import math
import os

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
DEG = math.pi/180.0
DELTA, DM, PHI0 = 7.5*DEG, 5.0*DEG, 30.0*DEG
SIGMA = 0.02
EPS_R = 0.25
FS_BUDGET = [4.0, 8.0, 16.0, 32.0]           # p4_budget.FS_F4, the set the budget minimises over
L_STAR_PROVED = 1424207.9721
VARTHETA_STATED = 0.2623138022


def Theta(u, L, er=EPS_R):
    return 0.5*(mp.tanh((u - er)/er) - mp.tanh((u - L + er)/er))


OUT = {}

# ---------------------------------------------------------------- A. the third branch of d
# THEOREM_S3_v3 sec.1.2 prints
#     d := rho_* min( sin(phi_0 - delta - 4 sigma), sin(pi/2 - delta_m - phi_0 - 4 sigma), f/rho_* )
# The third entry has to be dimensionless like the other two.  ASSEMBLY / pmax-h3v sec.B1 has
#     d = min( f, rho_* sin(phi_0 - delta), rho_* sin(pi/2 - delta_m - phi_0) ) rho_0 ,
# i.e. the inner branch is f rho_0, so the nondimensional entry is f rho_0/rho_*, not f/rho_*.
# As printed, rho_* x (f/rho_*) = f, a pure number, which is a length only in the units rho_0 = 1
# that every script uses.
A = {"printed_third_branch": "f/rho_*", "correct_third_branch": "f rho_0/rho_*",
     "printed_value_of_the_third_branch_times_rho_star": "f (dimensionless)",
     "correct_value": "f rho_0",
     "note": "identical in the units rho_0 = 1 that fix5/x1 and x2 use, so no number moves"}
rows = []
for f in FS_BUDGET + [0.05, 0.5, 1.0, 2.0, 100.0, 1e6]:
    rs = 1.0 + f
    b1 = rs*math.sin(PHI0 - DELTA - 4*SIGMA)
    b2 = rs*math.sin(math.pi/2 - DM - PHI0 - 4*SIGMA)
    b3 = f                                            # = f rho_0 at rho_0 = 1
    rows.append({"f": f, "taper": b1, "equatorial": b2, "inner_f_rho0": b3,
                 "binding": ["taper", "equatorial", "inner"][int(min(
                     range(3), key=lambda i: [b1, b2, b3][i]))],
                 "d": min(b1, b2, b3)})
A["branches"] = rows
A["taper_binds_at_every_f_in_the_budget_set"] = all(
    r["binding"] == "taper" for r in rows if r["f"] in FS_BUDGET)
A["taper_binds_for_all_f_ge"] = "f >= 0.4437: f >= (1+f) sin(phi_0-delta-4 sigma) fails, so the "
A["taper_binds_iff"] = ("f sin(phi_0-delta-4sigma)/(1-sin(phi_0-delta-4sigma)) <= f, i.e. "
                        "f >= sin(.)/(1-sin(.)) = %.10f"
                        % (math.sin(PHI0 - DELTA - 4*SIGMA)/(1 - math.sin(PHI0 - DELTA - 4*SIGMA))))
OUT["A_the_third_branch_of_d"] = A

# ---------------------------------------------------------------- B. geometric containments
B = {}
sin_off = math.sin(PHI0 - DELTA - 4*SIGMA)
grows = []
for f in FS_BUDGET:
    rs = 1.0 + f
    d = rs*sin_off
    r_star = rs*math.sin(PHI0)
    z_star = rs*math.cos(PHI0)
    psi_max = math.asin(d/rs)
    grows.append({
        "f": f, "rho_star": rs, "d": d, "r_star": r_star, "z_star": z_star,
        "R_minus_eq_r_star_minus_d": r_star - d,
        "d_lt_r_star_required_by_H3prime": bool(d < r_star),
        "phi_min_on_ball_deg": (PHI0 - psi_max)/DEG,
        "phi_min_minus_delta_over_sigma": (PHI0 - psi_max - DELTA)/SIGMA,
        "phi_max_on_ball_deg": (PHI0 + psi_max)/DEG,
        "clearance_to_the_equatorial_kink_over_sigma":
            (math.pi/2 - DM - (PHI0 + psi_max))/SIGMA,
        "z_min_on_ball": z_star - d, "z_min_positive": bool(z_star - d > 0),
        "rho_min_on_ball": rs - d, "rho_max_on_ball": rs + d,
        "rho_min_over_rho_0": rs - d})
B["per_f"] = grows
B["scale_free_identities"] = {
    "d/rho_star": sin_off,
    "R_minus/rho_star": math.sin(PHI0) - sin_off,
    "phi_min = delta + 4 sigma exactly, for every f": True,
    "phi_max = 2 phi_0 - delta - 4 sigma exactly, for every f": True,
    "d < r_star  <=>  sin(phi_0-delta-4sigma) < sin(phi_0)": bool(sin_off < math.sin(PHI0)),
    "z_min > 0  <=>  cos(phi_0) > sin(phi_0-delta-4sigma)":
        bool(math.cos(PHI0) > sin_off)}
B["equatorial_branch_never_binds_because"] = (
    "sin(pi/2 - delta_m - phi_0 - 4 sigma) = %.10f > sin(phi_0 - delta - 4 sigma) = %.10f"
    % (math.sin(math.pi/2 - DM - PHI0 - 4*SIGMA), sin_off))
OUT["B_geometric_containments"] = B

# ---------------------------------------------------------------- C. the quantifier "f >= 4"
# (H3'_vartheta) with the stated vartheta is FALSE once the tracked point leaves the shell
# plateau, which happens at f of order e^L.  At x = x_*, with A_sigma(phi_0)/N_sigma = 0.98909112
# (fix4 sec.1.1), the k = 0 clause reads  |1 - Theta(log rho_*) A_sigma(phi_0)/N_sigma| <= vartheta.
mp.mp.dps = 40
AN = mp.mpf("0.9890911210")
C = {"A_sigma(phi_0)/N_sigma": float(AN), "vartheta_stated": VARTHETA_STATED}
rowsC = []
for L in [40.0, 1e5, L_STAR_PROVED]:
    for lab, u in [("u = log rho_* = L - 1", mp.mpf(L) - 1), ("u = L - 0.5", mp.mpf(L) - mp.mpf("0.5")),
                   ("u = L", mp.mpf(L)), ("u = L + 1", mp.mpf(L) + 1)]:
        th = Theta(u, mp.mpf(L))
        t0 = abs(1 - th*AN)
        rowsC.append({"L": L, "where": lab, "f_such_that_log(1+f)=u": float(mp.e**u - 1),
                      "Theta": float(th), "vartheta_0_at_x_star": float(t0),
                      "exceeds_the_stated_vartheta": bool(t0 > VARTHETA_STATED),
                      "exceeds_1": bool(t0 > 1)})
C["rows"] = rowsC
# the exact threshold in u
thr = mp.findroot(lambda x: abs(1 - Theta(x, mp.mpf(40.0))*AN) - VARTHETA_STATED,
                  (mp.mpf(38), mp.mpf(41)), solver="bisect", tol=mp.mpf("1e-30"))
C["threshold_u_at_L_40"] = float(thr)
C["threshold_u_minus_L"] = float(thr - 40)
C["reading"] = ("the k = 0 clause of (H3'_vartheta) fails at the stated vartheta as soon as "
                "log(1+f) > L + %.6f; the statement quantifies over every f >= 4 with no upper "
                "bound, while the budget minimises only over f in %s"
                % (float(thr - 40), FS_BUDGET))
OUT["C_the_f_ge_4_quantifier"] = C

# ---------------------------------------------------------------- D. is vartheta L-independent?
# Theta_L(u) - Theta_inf(u) = -(1/2)(tanh((u-L+eps_r)/eps_r) + 1), and
# 0 <= (1/2)(tanh w + 1) <= e^{2w}; every u-derivative carries the same factor because
# d^k/du^k of (1/2)(tanh w + 1) is (1/2) eps_r^{-k} p_k(T) with p_k divisible by (1 - T^2)
# <= 2(1+T) <= 4 e^{2w}.
D = {}
rowsD = []
for f in FS_BUDGET:
    umax = math.log((1.0 + f)*(1.0 + math.sin(PHI0 - DELTA - 4*SIGMA)))
    for L in [40.0, L_STAR_PROVED]:
        w = (umax - L + EPS_R)/EPS_R
        rowsD.append({"f": f, "u_max_on_ball": umax, "L": L,
                      "bound_on_|Theta_L - Theta_inf| and every derivative/(4/eps_r)^k":
                          float(mp.e**(2*mp.mpf(w))*4)})
D["rows"] = rowsD
D["verdict"] = ("vartheta measured at L = 40 (x1's choice) equals vartheta at every L >= 40 to "
                "the displayed precision and far beyond, for f in the budget's set: the outer "
                "tanh contributes below 1e-100.  PROFILE-INDEPENDENT IN L, proved, for bounded f; "
                "the same estimate is what fails in item C when f grows like e^L.")
OUT["D_L_independence_of_vartheta"] = D

# ---------------------------------------------------------------- E. N in the viscous budget
E = {"N_used_by_t3_budget": "r_*/(M sin delta) i.e. ||eta_0||_inf r_*/M with the KINKED sup",
     "kinked_sup_rho|eta_0|/M": 1.0/math.sin(DELTA),
     "mollified_certified_E_0_upper": 7.5523079773,
     "N_is_an_over-estimate": bool(1.0/math.sin(DELTA) > 7.5523079773),
     "ratio": (1.0/math.sin(DELTA))/7.5523079773,
     "E_hess_and_E_tail_are_increasing_in_N":
         "E_hess = (r_*/R_-^2)(1/2)K_2 e^c tau V_1 + 4 N (1/2) K_2 e^c tau V_1/d and "
         "E_tail = 2 N Pfrak + Pfrak + sum_k r_*^{-k}(E|Z^a|^{2k})^{1/2} Pfrak^{1/2}: both "
         "affine increasing in N, so the substitution is conservative",
     "category": "NOT recomputed for the mollified datum; proved conservative"}
OUT["E_the_N_constant_in_the_budget"] = E

# ---------------------------------------------------------------- F. window vs conclusion
F = {"estimates_hold_on": "[0, tau], tau = c/(M L)",
     "conclusion_asserts": "T_d(u_0) <= t_* = c/(M L)",
     "identical": True,
     "c_range_in_the_statement": [0.8244732108, 1.1647266909*0.8244732108],
     "certified_monotone_cap_eps": 0.1647266909,
     "cap_consistent": True,
     "note": ("the conclusion is exactly the right endpoint of the window, so the theorem "
              "concludes at the last instant its estimates cover; a first-crossing at t = tau "
              "itself is inside the closed window, and the dichotomy of section 4 has to be "
              "stated on the CLOSED interval, which is what the continuation step must deliver")}
OUT["F_window_versus_conclusion"] = F

# ---------------------------------------------------------------- G. the 3.5 sigma cliff claim
G = {"claim_in_THEOREM_S3_v3_sec3_item3":
     "at an offset of 3.5 sigma instead of 4 sigma, 2 vartheta = 1.76645 and the hypothesis fails",
     "what_x1_actually_measured": 0.883225,
     "that_is": "a five-parameter GRID maximum, i.e. a LOWER bound on the supremum",
     "a_lower_bound_of_0.883225_is_below_1": True,
     "so_the_hypothesis_is_not_refuted_at_3.5_sigma_by_this_number": True,
     "reading": ("the cliff is asserted with the same SAFETY = 2 convention that is used to "
                 "assert the hypothesis at 4 sigma; the convention cannot both certify a bound "
                 "and refute one.  r2_theta_cert.py settles the 3.5 sigma case with a certified "
                 "upper bound.")}
OUT["G_the_3.5_sigma_cliff"] = G

json.dump(OUT, open(os.path.join(HERE, "r3_results.json"), "w"), indent=1, default=str)
print("A. third branch of d: printed f/rho_*, correct f rho_0/rho_*; taper binds at every f in",
      FS_BUDGET, "->", OUT["A_the_third_branch_of_d"]["taper_binds_at_every_f_in_the_budget_set"])
print("B. containments: d < r_*, z_min > 0, phi_min = delta + 4 sigma exactly, equatorial "
      "clearance %.1f sigma at f = 4" % grows[0]["clearance_to_the_equatorial_kink_over_sigma"])
print("C. (H3'_vartheta) at the stated vartheta fails once log(1+f) > L %+0.6f"
      % OUT["C_the_f_ge_4_quantifier"]["threshold_u_minus_L"])
for r in rowsC[-4:]:
    print("     L = %.4f  %-22s  Theta = %.6g  vartheta_0 = %.6f  > stated: %s"
          % (r["L"], r["where"], r["Theta"], r["vartheta_0_at_x_star"],
             r["exceeds_the_stated_vartheta"]))
print("D. L-independence bound at f = 32, L = L_*: %.3e" % rowsD[-1][
    "bound_on_|Theta_L - Theta_inf| and every derivative/(4/eps_r)^k"])
print("E. N = 1/sin delta = %.10f against the certified E_0 <= %.10f: conservative"
      % (E["kinked_sup_rho|eta_0|/M"], E["mollified_certified_E_0_upper"]))
print("G. the 3.5 sigma 'failure' is a lower bound 0.883225 < 1, so nothing is refuted there")
print("wrote r3_results.json")

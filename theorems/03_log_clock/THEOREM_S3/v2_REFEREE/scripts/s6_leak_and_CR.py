"""
s6 -- (a) STALE-CONSTANT LEAK TEST.  p4_budget rebinds t3_budget's KAPPA, CSTAR, C2, E0, G0,
      SHIFT but NOT its LAM_MAX and R_H (the KINKED datum's).  Poison them with NaN and check
      that no number on the p4 path changes: that proves lam_max and r_h are taken at the
      self-consistent c everywhere.
      (b) C_R's supremum over the angle: p4 uses C_R_sup_fast (160-point coarse scan + a
      41-point refinement).  Compare with the imported T2.C_R_sup at n = 400.
      (c) the enclosure of the grid maxima E_0, Gfrak_0, N_sigma, Psi_sigma(0): pad the fix4
      grid values with the rigorous second-derivative bounds of s1.
Outputs -> s6_results.json
"""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
FIX4 = os.path.join(os.path.dirname(HERE), "copies", "THEOREM_S3", "fix4")
sys.path.insert(0, FIX4); sys.path.insert(0, os.path.join(FIX4, "imported"))
_cwd = os.getcwd()
import p4_budget as P4
os.chdir(_cwd)
T3, T2 = P4.T3, P4.T2

P1RES = json.load(open(os.path.join(FIX4, "p1_results.json")))
ROW = P1RES["rows"]["sigma=0.02"]
E0, G0 = ROW["E0"], ROW["Gfrak0"]
CAP = 1.0 + ROW["eps_cap_certified"]

if __name__ == "__main__":
    R = {}
    D = P4.Datum(0.02, 2.9135781820)
    P4.bind(D, E0, G0)
    c = 1.0669124225*D.c_star
    base = P4.assemble_c(D, 2e6, c, "proved", E0, G0, f=4.0)
    keep = (T3.LAM_MAX, T3.R_H)
    T3.LAM_MAX, T3.R_H = float("nan"), float("nan")
    P4._CC.clear(); P4._BS.clear()
    pois = P4.assemble_c(D, 2e6, c, "proved", E0, G0, f=4.0)
    T3.LAM_MAX, T3.R_H = keep
    R["stale_constant_leak_test"] = {
        "t3_LAM_MAX_kinked": keep[0], "t3_R_H_kinked": keep[1],
        "identical_under_NaN_poisoning": all(
            (isinstance(base[k], str) or (not np.isfinite(base[k]) and not np.isfinite(pois[k]))
             or base[k] == pois[k]) for k in base),
        "eps_base": base["eps"], "eps_poisoned": pois["eps"],
        "lam_max_used": base["lam_max"], "r_h_used": base["r_h"],
        "theorem_sec_1_2_r_h": 0.9033435309,
        "note": "lam_max and r_h are taken at the self-consistent c; the r_h in the "
                "statement's constant table is the one at c_*, not the one the proof uses"}
    print("leak test:", R["stale_constant_leak_test"]["identical_under_NaN_poisoning"],
          " r_h used at the self-consistent window =", base["r_h"], flush=True)

    # (b) C_R supremum: fast search vs the imported full search
    K = P4.column_constants(2e6, c, "proved", E0, G0)
    cg, gh = K["c_G"], K["Ghat"]
    fastv = {}
    for n in (80, 160, 320):
        v, a = P4.C_R_sup_fast(1.5, cg, E0, gh, ca_flat=None, grad_sinphi=True, n=n)
        fastv["fast n=%d" % n] = {"C_R": v, "argmax_deg": a}
    full = {}
    for n in (200, 400, 800):
        v, a = T2.C_R_sup(1.5, cg, E0, gh, ca_flat=None, grad_sinphi=True, n=n)
        full["T2.C_R_sup n=%d" % n] = {"C_R": v, "argmax_deg": a}
    R["C_R_supremum_search"] = {"c_G": cg, "Ghat": gh, "fast": fastv, "full": full,
                                "fix4_quoted_C_R_at_L2e6": 642.7111}
    print("C_R fast/full:", json.dumps({**fastv, **full}), flush=True)

    # (c) enclosures for the grid maxima, from s1's rigorous derivative bounds
    s1 = json.load(open(os.path.join(HERE, "s1_results.json")))
    B = s1["rigorous_bounds"]
    h = math.pi/600000.0
    R["grid_maximum_enclosures"] = {
        "grid_step_fix4": h,
        "N_sigma": {"fix4": 1.0124508488, "rigorous_pad": 0.125*h*h*B["B_A2"],
                    "my_enclosure": [s1["N_sigma"]["lower"], s1["N_sigma"]["upper"]]},
        "E_0": {"fix4": 7.5523076942, "rigorous_pad_on_supW": 0.125*h*h*B["B_W2"],
                "my_enclosure": [s1["E_0"]["lower"], s1["E_0"]["upper"]]},
        "sup_h_second": {"fix4_called_a_PROVED_Lipschitz_bound": 225.9466903722,
                         "my_scan_of_the_same_quantity": s1["sup_abs_h_second"]["scan_value"],
                         "rigorous_a_priori_bound": s1["sup_abs_h_second"][
                             "rigorous_a_priori_bound_B_A2_over_N"],
                         "fix4_pad_on_sup_h_prime": 225.9466903722*h,
                         "rigorous_pad_on_sup_h_prime": s1["sup_abs_h_second"][
                             "rigorous_a_priori_bound_B_A2_over_N"]*h,
                         "sup_h_prime": s1["sup_abs_h_prime_scan"]}}
    json.dump(R, open(os.path.join(HERE, "s6_results.json"), "w"), indent=1, default=str)
    print(json.dumps(R, indent=1, default=str))

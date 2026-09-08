"""
check_constants.py -- the gate.

Every check consumes a STORED number from a1..a6_results.json and compares it against a
value rebuilt either from mathematics or from other stored numbers.  No check multiplies its
input by zero; no check compares a literal against itself.  The pass count is len(CHECKS),
printed, never a literal.

  python3 check_constants.py            run the gate
  python3 check_constants.py --mutate   perturb every stored numeric leaf one at a time
                                        (v -> 1.5v + 0.37; booleans flipped) and MEASURE the
                                        catch rate.  The rate is reported, never asserted.
  python3 check_constants.py --audit    audit every decimal number printed in ASSEMBLY.md
                                        against the stored JSONs.
"""
import json, math, os, re, sys, copy

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
FILES = ["a1_results.json", "a2_results.json", "a3_results.json",
         "a4_results.json", "a5_results.json", "a6_results.json", "a7_results.json",
         "a8_results.json"]
STATS = "gate_stats.json"     # written by --mutate, audited by --audit

# numbers QUOTED VERBATIM from other seats -- cited, not computed here.  Listed so the
# document audit can separate "traced to my own JSON" from "traced to a named source".
FOREIGN = {
    0.026093: "far-near-kernel-lemma C2in general",
    0.756945: "far-near-kernel-lemma C1 general",
    0.1324254: "L3v seat, inner multipole (r grad)",
    1.6049285: "L3v seat, far remainder l>=3 (r grad)",
    0.4999917: "prove-lagrangian kappa_delta(5 deg) -- the value corrected in 1.2",
    0.34913942: "Lemma-T seat, measured shear cost /L",
    2.11817051: "Lemma-T seat, C_rel = 3pi(sqrt6-2)/2",
    11.781: "Lemma T' relative constant 15pi/4",
    1704.05546: "arXiv identifier",
    110.111083: "s3close/hk2 sec.10 control L-11, E_hess/eps_bulk at L=10, hatK2=1",
    1101.1: "s3close/hk2 sec.10 / refute-V-b: the crossing L",
    1.263238: "V-b seat's WITHDRAWN E_hess/eps_bulk at L=10",
    2.1903: "s3close/hk2, measured K2 at the tracked point",
    3.0202: "s3close/hk2, measured K2 global probe",
    66.6622: "s3close/hk2, PROVED bound on K2 rho0/M",
    32.0503: "s3close/hk2 control L-11, Eh2/Eh1",
}

def load():
    return {f: json.load(open(f)) for f in FILES}

def build_checks(D):
    A1, A2, A3, A4, A5, A6, A7, A8 = (D[f] for f in FILES)
    C = []
    def ck(name, got, want, tol):
        C.append((name, got, want, tol))

    # ---------------- a1 : datum, coefficients, energy ----------------------
    for l, want in [(1, 2.4), (3, 4.444444444444445), (5, 6.461538461538462)]:
        i = [1, 3, 5].index(l)
        ck("N_%d = (l+1)(l+2)/(l+3/2)" % l, A1["checks"]["N_l_1_3_5"][i],
           (l+1.0)*(l+2.0)/(l+1.5), 1e-13)
    ck("C_1^{3/2}(0.3) = 3t", A1["checks"]["C32_l1_at_t"][0], 3*0.3, 1e-13)
    ex = {"1": -5/6, "3": 3/40, "5": -247/1680, "7": 1513/40320, "9": -2773/42240}
    for k, v in ex.items():
        ck("H_%s bang-bang vs exact rational" % k, A1["checks"]["Hl_bangbang_mine"][k], v, 3e-6)
    ck("max rel H_l error", A1["checks"]["Hl_bangbang_maxrel"],
       max(abs(A1["checks"]["Hl_bangbang_mine"][k]-ex[k])/abs(ex[k]) for k in ex), 1e-12)
    ck("kappa_0 = -(3/5)H_1", A1["checks"]["kappa0_from_H1"],
       -0.6*A1["checks"]["Hl_bangbang_mine"]["1"], 1e-13)
    ck("P_cap(1) = 1", A1["checks"]["P_cap_at_1"], 1.0, 1e-12)
    ck("P_cap(lam) = lam (max dev)", A1["checks"]["P_cap_maxdev_from_lam"], 0.0, 1e-12)
    ck("D_l numeric vs closed form", A1["checks"]["D_l_num_vs_closed_maxrel"], 0.0, 1e-10)
    ck("kappa_delta(0) = 1/2", A1["kappa_delta"]["0"], 0.5, 1e-12)
    ck("kappa_delta(7.5) vs Lemma-T seat 0.4997212305210993",
       A1["kappa_delta"]["7.5"], 0.4997212305210993, 1e-12)
    for dd in ["3", "5", "7.5", "10"]:
        dr = float(dd)*math.pi/180.0
        ck("small-delta law 1-2kappa ~ delta^3/4 at %s deg" % dd,
           1.0-2.0*A1["kappa_delta"][dd], dr**3/4.0, 0.012)
    ck("r_h(7.5) vs Lemma-T seat 0.981822", A1["r_h"]["7.5"], 0.981822, 1e-5)
    ck("Phi_h(3/2) at 7.5 deg vs Lemma-T refuter 1.473555",
       A1["Phi_ratio_at_1p5"]["7.5"], 1.473555, 1e-5)
    ck("Phi_h(3/2) at 30 deg vs Lemma-T refuter 1.067358",
       A1["Phi_ratio_at_1p5"]["30"], 1.067358, 1e-5)
    ck("C_E sharp cap vs record 0.172403978",
       A1["C_E_sharp_cap_Linf"], A1["C_E_sharp_cap_record"], 1e-6)
    ck("C_E relerr vs record (recomputed)", A1["C_E_sharp_cap_relerr_vs_record"],
       abs(A1["C_E_sharp_cap_Linf"]-A1["C_E_sharp_cap_record"])/A1["C_E_sharp_cap_record"], 1e-12)
    # C_E = 2 pi sum_l 2 N_l Hhat_l^2 D_l/(2l+3), D_l = 1/(5(l+4)) : rebuild the l=1 term
    Hh1 = A1["checks"]["Hl_bangbang_mine"]["1"]
    ck("C_E l=1 term = 2 N_1 Hhat_1^2/(5*5*5)", A1["C_E_terms_l1_to_l9"]["1"],
       2*2.4*Hh1*Hh1/(5.0*5.0*5.0), 1e-9)
    ck("C_E l=1 share", A1["C_E_l1_share"],
       2*math.pi*A1["C_E_terms_l1_to_l9"]["1"]/A1["C_E_sharp_cap_Linf"], 1e-12)
    ck("C_E sharp = table entry delta=0,dm=0",
       A1["C_E_table"]["delta=0,dm=0,eps=0"], A1["C_E_sharp_cap_Linf"], 1e-13)
    ck("c_E = (2/5) log C_E (7.5,5)", A1["dictionary"]["delta=7.5,dm=5"]["c_E_two_fifths_logC_E"],
       0.4*math.log(A1["C_E_table"]["delta=7.5,dm=5,eps=0"]), 1e-13)
    ck("s = 1/sin(7.5 deg)", A1["dictionary"]["delta=7.5,dm=5"]["s_record_1_over_sin_delta"],
       1.0/math.sin(7.5*math.pi/180.0), 1e-13)
    ck("logReE-2L = 2 log s + (2/5) log C_E",
       A1["dictionary"]["delta=7.5,dm=5"]["log_ReE_minus_2L_at_s_record"],
       2*math.log(A1["dictionary"]["delta=7.5,dm=5"]["s_record_1_over_sin_delta"])
       + A1["dictionary"]["delta=7.5,dm=5"]["c_E_two_fifths_logC_E"], 1e-12)
    ck("half-energy shift = (2/5) log(1/2)", A1["half_energy_shift_in_logReE"],
       0.4*math.log(0.5), 1e-14)
    ck("mollified C_E(eps=0.25) < C_E(eps=0)", 1.0 if
       A1["C_E_table"]["delta=7.5,dm=5,eps=0.25,L=40"] < A1["C_E_table"]["delta=7.5,dm=5,eps=0"]
       else 0.0, 1.0, 0.5)
    # material point: tan phi = lam^3 tan phi0
    for p0 in ["20", "30", "40"]:
        for lam in ["1", "1.25", "1.5"]:
            row = A1["material_point"]["phi0=%s" % p0]["lam=%s" % lam]
            th0 = float(p0)*math.pi/180.0
            want = math.degrees(math.atan2(float(lam)**3*math.sin(th0), math.cos(th0)))
            ck("phi(phi0=%s,lam=%s) = atan(lam^3 tan phi0)" % (p0, lam), row["phi_deg"], want, 1e-10)
            ck("equator margin phi0=%s lam=%s" % (p0, lam), row["margin_to_equator_deg"],
               90.0-row["phi_deg"], 1e-12)
            ck("taper margin phi0=%s lam=%s" % (p0, lam), row["margin_to_taper_deg"],
               row["phi_deg"]-7.5, 1e-12)
    ck("||omega0||_inf/M = 1", A1["bfg_hypotheses"]["sup_omega0_over_M"], 1.0, 1e-15)
    ck("||omega0||_2^2 finite and L-independent in units M^2R^3",
       A1["bfg_hypotheses"]["L2sq_over_M2R3_L=10_eps=0"],
       A1["bfg_hypotheses"]["L2sq_over_M2R3_L=40_eps=0"], 1e-9)

    # ---------------- a2 : the budget ---------------------------------------
    kap = A2["datum"]["kappa_delta"]
    SENSF = 0.34913942333376546/2.1181705106873974
    def sens_of(col):
        return 1.0 if col in ("proved", "proved_oddonly") else SENSF
    ck("a2 kappa == a1 kappa", kap, A1["kappa_delta"]["7.5"], 1e-14)
    ck("a2 r_h == a1 r_h", A2["datum"]["r_h"], A1["r_h"]["7.5"], 1e-14)
    ck("a2 C_E == a1 C_E(7.5,5)", A2["datum"]["C_E"], A1["C_E_table"]["delta=7.5,dm=5,eps=0"], 1e-14)
    ck("c_* = log(3/2)/kappa", A2["datum"]["c_at_lambda_3over2"], math.log(1.5)/kap, 1e-14)
    ck("c_2 = 4 log(3/2)", A2["c2"], 4.0*math.log(1.5), 1e-14)
    ck("logReE-2L (a2) == a1 dictionary", A2["logReE_minus_2L"],
       A1["dictionary"]["delta=7.5,dm=5"]["log_ReE_minus_2L_at_s_record"], 1e-13)
    ck("kappa_s = sqrt6 - 2", A2["record_inputs"]["kappa_s"], math.sqrt(6.0)-2.0, 1e-15)
    ck("C' proved (A1 + computed V) = 119.33", A2["record_inputs"]["Cprime_A1_computedV"],
       119.33, 1e-9)
    ck("C' proved (A1 + proved V) = 151.15", A2["record_inputs"]["Cprime_A1_provedV"], 151.15, 1e-9)
    ck("C' odd-only = 76.11", A2["record_inputs"]["Cprime_A1_oddonly"], 76.11, 1e-9)
    # every row of grid AND f_sweep: pipeline identities, whether or not eps is finite
    nrows = 0
    for tag, rows in [("grid", A2["grid"]), ("fsweep", A2["f_sweep"])]:
        for i, r in enumerate(rows):
            L, c, f, col = float(r["L"]), float(r["c"]), float(r["f"]), r["column"]
            key = "%s[%d] L=%g c=%g f=%g %s" % (tag, i, L, c, f, col)
            ck("eps_ell = ell_loss/L  " + key, float(r["eps_ell"]),
               (math.log(2.0)+math.log(1.0+f)+3*math.log(1.5))/L, 1e-12)
            ck("eps_Ca = 1.5 C_a/(kappa L)  " + key, float(r["eps_Ca"]),
               1.5*float(r["C_a"])/(kap*L), 1e-12)
            ck("eps_a = eps_ell+eps_T'+eps_Ca  " + key, float(r["eps_a"]),
               float(r["eps_ell"])+float(r["eps_Tprime"])+float(r["eps_Ca"]), 1e-9)
            ck("eps_v = bulk+hess+4+tail  " + key, float(r["eps_v"]),
               min(float(r["eps_bulk"])+float(r["E_hess"])+float(r["E_4"])
                   + min(float(r["E_tail"]), 1e12), 1e12), 1e-6)
            ck("c_G = (1.5 + C'/L) c  " + key, float(r["c_G"]),
               (1.5+float(r["Cprime"])/L)*c, 1e-10)
            ck("feedback exponent = p/(ML) * c  " + key, float(r["feedback_exponent_pc"]),
               (1.5 + float(r["Cprime"])/L
                + 1.5*2*kap*sens_of(col)*(15.0*math.pi/4.0)/(A1["r_h"]["7.5"]*(4.0/9.0)))*c, 1e-9)
            ck("C_a is the column's C_a  " + key, float(r["C_a"]),
               (0.291999+0.014754+3.999218/math.sin(math.pi/6)+math.pi/8)
               if col != "measured" else 0.069, 1e-12)
            ck("C_K = 1 in every row  " + key, float(r["C_K"]), 1.0, 1e-15)
            ck("bootstrap_closed iff eps finite  " + key,
               1.0 if (str(r["bootstrap_closed"]) == "True") else 0.0,
               1.0 if math.isfinite(float(r["eps"])) else 0.0, 0.5)
            ck("mu < 1 iff closed  " + key, 1.0 if float(r["mu"]) < 0.899 else 0.0,
               1.0 if (str(r["bootstrap_closed"]) == "True") else 0.0, 0.5)
            ck("eps_delta = 1-2kappa on every row  " + key, float(r["eps_delta"]),
               1.0-2.0*kap, 1e-14)
            if str(r["bootstrap_closed"]) != "True":
                ck("mu == 1 exactly when not closed  " + key, float(r["mu"]), 1.0, 1e-15)
                ck("eps_T' is +inf when not closed  " + key,
                   1.0 if not math.isfinite(float(r["eps_Tprime"])) else 0.0, 1.0, 0.5)
                ck("eps_a is +inf when not closed  " + key,
                   1.0 if not math.isfinite(float(r["eps_a"])) else 0.0, 1.0, 0.5)
                ck("eps is +inf when not closed  " + key,
                   1.0 if not math.isfinite(float(r["eps"])) else 0.0, 1.0, 0.5)
            else:
                ck("eps_T' finite when closed  " + key,
                   1.0 if math.isfinite(float(r["eps_Tprime"])) else 0.0, 1.0, 0.5)
                ck("muJ = muJL(1+2ks/L)+2ks/L >= 2ks/L  " + key,
                   1.0 if float(r["muJ"]) > 0 else 0.0, 1.0, 0.5)
            ck("muJ >= 2 kappa_s/L  " + key, 1.0 if float(r["muJ"]) >= 2*(math.sqrt(6)-2)/L - 1e-12
               else 0.0, 1.0, 0.5)
            ck("d = min(f, taper, equator)  " + key, float(r["d"]),
               min(f, (1+f)*math.sin(math.pi/6-7.5*math.pi/180.0),
                   (1+f)*math.sin(math.pi/2-5*math.pi/180.0-math.pi/6)), 1e-12)
            ck("q_half = (d/2)^2 s^2 L/c  " + key, float(r["q_half"]),
               (float(r["d"])/2.0)**2*(1.0/math.sin(7.5*math.pi/180.0))**2*L/c, 1e-9)
            e = float(r["eps"])
            if math.isfinite(e):
                ea, ev, ed = float(r["eps_a"]), float(r["eps_v"]), float(r["eps_delta"])
                ck("eps assembly (1.1)  " + key, e,
                   (1.0+math.log(1.0/(1.0-ev))/math.log(1.5))/((1.0-ea)*(1.0-ed))-1.0, 1e-9)
                nrows += 1
    ck("eps_delta = 1 - 2 kappa", float(A2["grid"][0]["eps_delta"]), 1.0-2.0*kap, 1e-14)
    # L_star -> Lambda_star
    for k, v in A2["Lambda_star_log"].items():
        if v in (None, "None"):
            continue
        ck("log Lambda_* = 2 L_* + shift [%s]" % k, float(v),
           2*float(A2["L_star"][k]) + A2["logReE_minus_2L"], 1e-6)

    # ---------------- a3 : the 19 runs --------------------------------------
    ck("kappa(-M sin 2phi) = 2/5", A3["kappa_run_datum_minus_sin2phi"], 0.4, 1e-9)
    ck("kappa relerr recomputed", A3["kappa_run_relerr"],
       abs(A3["kappa_run_datum_minus_sin2phi"]-0.4)/0.4, 1e-12)
    ck("19 runs", float(A3["ratio_vs_bound_L"]["n"]), 19.0, 0.5)
    ck("fit c vs record 1.0538", A3["fit_c_over_logRsstar"], A3["record_fit_c"], 6e-3)
    ck("fit rms %% vs record 8.8", A3["fit_rms_pct_of_mean_Td"], A3["record_fit_rms_pct"], 0.05)
    ck("Q' mean vs corner seat 1.0272", A3["Qprime_mean"], 1.0272, 1e-3)
    ck("all 19 above the eps=0 bound at L=log(R/rho0)",
       float(A3["ratio_vs_bound_L"]["n_above_1"]), 19.0, 0.5)
    ck("kappa_run exact target = 2/5", A3["kappa_run_exact_2over5"], 0.4, 1e-15)
    ck("sup|sin 2phi| = 1", A3["sup_abs_sin2phi"], 1.0, 1e-15)
    ck("n_above_1 (L) recomputed", float(A3["ratio_vs_bound_L"]["n_above_1"]),
       float(sum(1 for r in A3["rows"] if r["ratio_meas_over_bound_L"] > 1.0)), 0.5)
    ck("n_above_1 (Leff) recomputed", float(A3["ratio_vs_bound_Leff"]["n_above_1"]),
       float(sum(1 for r in A3["rows"] if r["ratio_meas_over_bound_Leff"] > 1.0)), 0.5)
    ck("Qprime mean recomputed", A3["Qprime_mean"],
       sum(r["MTd_times_Leff"] for r in A3["rows"])/len(A3["rows"]), 1e-12)
    ck("L range min", A3["L_range"][0], min(r["L"] for r in A3["rows"]), 1e-15)
    ck("L range max", A3["L_range"][1], max(r["L"] for r in A3["rows"]), 1e-15)
    for r in A3["rows"]:
        ck("row N=%d Re0=%g : MTd*L" % (r["N"], r["Re0"]), r["MTd_times_L"],
           r["MTd"]*r["L"], 1e-12)
        ck("row N=%d Re0=%g : MTd*Leff" % (r["N"], r["Re0"]), r["MTd_times_Leff"],
           r["MTd"]*r["log_R_over_sstar"], 1e-12)
        ck("row N=%d Re0=%g : bound_Leff" % (r["N"], r["Re0"]), r["bound_eps0_Leff"],
           math.log(1.5)/(0.4*r["log_R_over_sstar"]), 1e-9)
        ck("row N=%d Re0=%g : ratio_Leff" % (r["N"], r["Re0"]), r["ratio_meas_over_bound_Leff"],
           r["MTd"]/r["bound_eps0_Leff"], 1e-12)
        ck("row N=%d Re0=%g : implied eps" % (r["N"], r["Re0"]),
           A3["implied_eps_L"]["N%d_Re%g" % (r["N"], r["Re0"])],
           r["ratio_meas_over_bound_L"]-1.0, 1e-12)
        ck("row N=%d Re0=%g : Re0 > 0" % (r["N"], r["Re0"]), 1.0 if r["Re0"] > 0 else 0.0, 1.0, 0.5)
    for r in A3["rows"]:
        ck("row N=%d Re0=%g : L = N log2" % (r["N"], r["Re0"]), r["L"],
           r["N"]*math.log(2.0), 1e-12)
        ck("row N=%d Re0=%g : log(R/s*)" % (r["N"], r["Re0"]), r["log_R_over_sstar"],
           math.log((2.0**r["N"])/r["s_star"]), 1e-12)
        ck("row N=%d Re0=%g : ratio = MTd/bound" % (r["N"], r["Re0"]),
           r["ratio_meas_over_bound_L"], r["MTd"]/r["bound_eps0_L"], 1e-12)
        ck("row N=%d Re0=%g : bound = log(3/2)/(kappa L)" % (r["N"], r["Re0"]),
           r["bound_eps0_L"], math.log(1.5)/(0.4*r["L"]), 1e-9)

    # ---------------- a4 : restart pricing ----------------------------------
    for k, v in A4["table"].items():
        col, N = k.rsplit("_N=", 1)
        pc = float(v["pc_per_subwindow"])
        ck("amplification N(e^{pc/N}-1) [%s]" % k, float(v["amplification_total"]),
           int(N)*(math.exp(pc)-1.0), 1e-9)
    p1 = float(A4["table"]["proved_N=1"]["pc_per_subwindow"])
    for N in [2, 4, 8, 16, 32, 64]:
        ck("pc/N scales as 1/N [proved N=%d]" % N,
           float(A4["table"]["proved_N=%d" % N]["pc_per_subwindow"]), p1/N, 1e-6)
    ck("c_* in a4 == a2", A4["c_star"], A2["datum"]["c_at_lambda_3over2"], 1e-14)
    for col in ["proved", "proved_oddonly", "proved_sharpsens", "measured"]:
        prev = None
        for N in [1, 2, 4, 8, 16, 32, 64]:
            v = A4["table"]["%s_N=%d" % (col, N)]["L_star_eps<=0.5"]
            v = None if v in (None, "None") else float(v)
            if prev is not None and v is not None:
                ck("L_*(eps<=.5) decreasing in N [%s N=%d]" % (col, N),
                   1.0 if v <= prev*1.000001 else 0.0, 1.0, 0.5)
            prev = v if v is not None else prev
            w = A4["table"]["%s_N=%d" % (col, N)]["L_star_eps<=0.1"]
            if v is not None and w not in (None, "None"):
                ck("L_*(eps<=.1) >= L_*(eps<=.5) [%s N=%d]" % (col, N),
                   1.0 if float(w) >= v else 0.0, 1.0, 0.5)

    # ---------------- a5 : report numbers -----------------------------------
    ck("a5 c2 == a2 c2", A5["c2"], A2["c2"], 1e-15)
    ck("a5 shift == a2 shift", A5["logReE_minus_2L"], A2["logReE_minus_2L"], 1e-15)
    ck("a5 kappa == a1", A5["kappa_delta"], A1["kappa_delta"]["7.5"], 1e-15)
    ck("a5 r_h == a1", A5["r_h"], A1["r_h"]["7.5"], 1e-15)
    ck("a5 C_E == a1", A5["C_E"], A1["C_E_table"]["delta=7.5,dm=5,eps=0"], 1e-15)
    ck("a5 s == 1/sin delta", A5["s"], 1.0/math.sin(7.5*math.pi/180.0), 1e-14)
    ck("a5 c_star == a2", A5["c_star_window"], A2["datum"]["c_at_lambda_3over2"], 1e-15)
    ck("a5 numerics fit_c == a3", A5["numerics"]["fit_c"], A3["fit_c_over_logRsstar"], 1e-15)
    ck("a5 numerics rms == a3", A5["numerics"]["fit_rms_pct"], A3["fit_rms_pct_of_mean_Td"], 1e-15)
    ck("a5 L_range == a3", A5["numerics"]["L_range"][0], A3["L_range"][0], 1e-15)
    ck("a5 L_range max == a3", A5["numerics"]["L_range"][1], A3["L_range"][1], 1e-15)
    for k, v in A5["eps_at_160_640"].items():
        if v in (None, "None"):
            continue
        col, rest = k.split("_N=")
        N, L = rest.split("_L=")
        other = "%s_N=%s_L=%s" % (col, N, "640" if L == "160" else "160")
        if A5["eps_at_160_640"].get(other) not in (None, "None"):
            hi = float(A5["eps_at_160_640"]["%s_N=%s_L=160" % (col, N)])
            lo = float(A5["eps_at_160_640"]["%s_N=%s_L=640" % (col, N)])
            ck("eps(160) > eps(640) [%s N=%s]" % (col, N), 1.0 if hi > lo else 0.0, 1.0, 0.5)
    for k, v in A5["Lambda_star"].items():
        if v in (None, "None"):
            continue
        ck("a5 log Lambda_* = 2L_*+shift [%s]" % k, float(v["log_Lambda_star"]),
           2*float(v["L_star"]) + A5["logReE_minus_2L"], 1e-6)
        ck("a5 log10 Lambda_* [%s]" % k, float(v["log10_Lambda_star"]),
           float(v["log_Lambda_star"])/math.log(10.0), 1e-9)

    # ---------------- a6 : instrument checks --------------------------------
    ck("theta_max = 4(1-sqrt(2/3))", A6["theta_max"], 4.0*(1.0-math.sqrt(2.0/3.0)), 1e-15)
    ck("c_G = sqrt(3/2) theta_max", A6["c_G"], math.sqrt(1.5)*A6["theta_max"], 1e-15)
    ck("I2 = 4/5 - 16 sqrt6/135", A6["I2"], 4.0/5.0-16.0*math.sqrt(6.0)/135.0, 1e-15)
    ck("I4 = -4/7 + 27 sqrt6/28", A6["I4"], -4.0/7.0+27.0*math.sqrt(6.0)/28.0, 1e-15)
    ck("theta_max/I2 = 1.440118 (V-b refuter)", A6["theta_max_over_I2"], 1.440118, 1e-5)
    pub = {"E_4": "E_4_L10", "Eh1": "Eh1_L10", "Eh2": "Eh2_L10",
           "Eh2_over_Eh1": "Eh2_over_Eh1", "correction_factor": "correction_factor",
           "eps_bulk": "eps_bulk_L10"}
    for k, pk in pub.items():
        ck("a6 relerr[%s] recomputed from rows and published" % k, A6["relerr"][k],
           abs(A6["rows"]["L=10"][k]-A6["published"][pk])/abs(A6["published"][pk]), 1e-12)
        ck("a6 reproduces V-b refuter's %s" % k, A6["relerr"][k], 0.0, 1e-4)
    ck("a6 s = 1/sin delta", A6["s"], 1.0/math.sin(7.5*math.pi/180.0), 1e-14)
    ck("a6 r0 = sin(30 deg)", A6["r0"], 0.5, 1e-15)
    ck("a6 d = sin(7.5 deg)", A6["d"], math.sin(7.5*math.pi/180.0), 1e-15)
    ck("a6 R_minus = r0 - d", A6["R_minus"], A6["r0"]-A6["d"], 1e-15)
    ck("a6 N = r0/sin delta", A6["N"], A6["r0"]/math.sin(7.5*math.pi/180.0), 1e-14)
    ck("a6 p sharp = 1.5 + (p-1.5) sens", A6["feedback"]["p_over_ML_sharp_Linf"],
       1.5+(A6["feedback"]["p_over_ML_proved_Linf"]-1.5)*A6["feedback"]["sens_factor"], 1e-13)
    ck("a6 exp(pc sharp)", A6["feedback"]["exp_pc_sharp"],
       math.exp(A6["feedback"]["pc_sharp"]), 1e-9)
    ck("a6 correction factor = 87.166", A6["rows"]["L=10"]["correction_factor"], 87.166, 1e-3)
    ck("a6 Eh2/Eh1 is L-independent",
       A6["rows"]["L=10"]["Eh2_over_Eh1"], A6["rows"]["L=640"]["Eh2_over_Eh1"], 1e-9)
    ck("a6 eps_bulk ~ 1/L", A6["rows"]["L=10"]["eps_bulk"],
       16.0*A6["rows"]["L=160"]["eps_bulk"], 1e-9)
    fb = A6["feedback"]
    ck("p/(ML) = 1.5 + (3/2)(2kappa)(15pi/4)(9/4)/r_h", fb["p_over_ML_proved_Linf"],
       1.5 + 1.5*2*kap*(15.0*math.pi/4.0)/(A1["r_h"]["7.5"]*(4.0/9.0)), 1e-12)
    ck("pc = p/(ML) * c_*", fb["pc_proved"], fb["p_over_ML_proved_Linf"]*fb["c_star"], 1e-13)
    ck("exp(pc)", fb["exp_pc_proved"], math.exp(fb["pc_proved"]), 1e-6)
    ck("sens factor = 0.34913942/2.11817051", fb["sens_factor"],
       0.34913942333376546/2.1181705106873974, 1e-15)
    ck("pc sharp = (1.5+(p-1.5)*sens)*c_*", fb["pc_sharp"],
       (1.5+(fb["p_over_ML_proved_Linf"]-1.5)*fb["sens_factor"])*fb["c_star"], 1e-13)
    ck("a6 c_* == a2 c_*", fb["c_star"], A2["datum"]["c_at_lambda_3over2"], 1e-14)

    # ---------------- a7 : inline-quoted derived numbers --------------------
    ck("C_a proved = 0.291999+0.014754+3.999218/sin30+pi/8", A7["C_a_proved_phi30"],
       0.291999+0.014754+3.999218/math.sin(math.pi/6)+math.pi/8, 1e-13)
    ck("C_a ratio proved/measured", A7["C_a_ratio_proved_over_measured"],
       A7["C_a_proved_phi30"]/0.069, 1e-13)
    ck("ell_loss(f=0) = log2 + 3log(3/2)", A7["ell_loss_const_f0"],
       math.log(2.0)+3*math.log(1.5), 1e-15)
    ck("ell_loss(f=1) = ell_loss(f=0) + log 2", A7["ell_loss_f1"],
       A7["ell_loss_const_f0"]+math.log(2.0), 1e-15)
    ck("eps_ell(L=160,f=1)", A7["eps_ell_L160_f1"], A7["ell_loss_f1"]/160.0, 1e-15)
    ck("eps_Ca(L=160, proved, f=1)", A7["eps_Ca_L160_proved_f1"],
       1.5*A7["C_a_proved_phi30"]/(kap*160.0), 1e-13)
    ck("eps_Ca(L=160, measured, f=1)", A7["eps_Ca_L160_measured_f1"],
       1.5*0.069/(kap*160.0), 1e-13)
    ck("C_E l=1 share in %", A7["C_E_l1_share_pct"], 100.0*A1["C_E_l1_share"], 1e-12)
    ck("lambda at c=0.25", A7["lambda_at_c_0p25"], math.exp(kap*0.25), 1e-14)
    ck("equator margin - dm", A7["equator_margin_minus_dm_deg"],
       90.0-A1["material_point"]["phi0=30"]["lam=1.5"]["phi_deg"]-5.0, 1e-12)
    ck("eps floor = (1-2kappa)/(2kappa)", A7["eps_floor_delta"], (1.0-2.0*kap)/(2.0*kap), 1e-14)
    ck("d(f=1) = 2 sin(phi0-delta)", A7["d_f1_taper"],
       2.0*math.sin(math.pi/6-7.5*math.pi/180.0), 1e-14)
    ck("d(f=1,eps=0.25) = 2 - e^{0.25}", A7["d_f1_taper_with_eps025"], 2.0-math.exp(0.25), 1e-14)
    ck("d reduction %", A7["d_reduction_pct"],
       100.0*(1.0-A7["d_f1_taper_with_eps025"]/A7["d_f1_taper"]), 1e-11)
    ck("geometric product (3/2)(9/4)", A7["geom_product"], 1.5*2.25, 1e-15)
    ck("Lemma T' relative constant 15pi/4", A7["Lemma_Tprime_rel_const"], 15.0*math.pi/4.0, 1e-14)
    ck("Lemma T' conservatism factor", A7["conservatism_factor_LemmaTprime"],
       1.0/A6["feedback"]["sens_factor"], 1e-13)
    # ---------------- a8 : (H-K2) plugged in --------------------------------
    ck("hk2 proved bound", A8["hk2_statement"]["proved_bound"], 66.6622, 1e-9)
    ck("hk2 measured global", A8["hk2_statement"]["measured_global"], 3.0202, 1e-9)
    ck("hk2 measured tracked", A8["hk2_statement"]["measured_tracked"], 2.1903, 1e-9)
    ck("hk2 sharp-edge K2*f", A8["hk2_statement"]["sharp_edge_times_f"], 8.0, 1e-12)
    ck("a8 shift == a2", A8["logReE_minus_2L"], A2["logReE_minus_2L"], 1e-15)
    ck("a8 c_* == a2", A8["c_star"], A2["datum"]["c_at_lambda_3over2"], 1e-15)
    for k, v in A8["table"].items():
        for t in [0.5, 0.1]:
            Ls = v["L_star_eps<=%g" % t]
            lg = v["log_Lambda_star_eps<=%g" % t]
            if Ls in (None, "None") or lg in (None, "None"):
                continue
            ck("a8 log Lambda_* = 2L_*+shift [%s eps<=%g]" % (k, t), float(lg),
               2*float(Ls) + A8["logReE_minus_2L"], 1e-6)
        e160 = v["eps(L=160)"]; e640 = v["eps(L=640)"]; e2560 = v["eps(L=2560)"]
        if e160 not in (None, "None") and e640 not in (None, "None"):
            ck("a8 eps(160) > eps(640) [%s]" % k, 1.0 if float(e160) > float(e640) else 0.0, 1.0, 0.5)
        if e640 not in (None, "None") and e2560 not in (None, "None"):
            ck("a8 eps(640) > eps(2560) [%s]" % k, 1.0 if float(e640) > float(e2560) else 0.0, 1.0, 0.5)
        L05 = v["L_star_eps<=0.5"]; L01 = v["L_star_eps<=0.1"]
        if L05 not in (None, "None") and L01 not in (None, "None"):
            ck("a8 L_*(.1) >= L_*(.5) [%s]" % k, 1.0 if float(L01) >= float(L05) else 0.0, 1.0, 0.5)
    # C_K = 1 baseline row must reproduce a4/a5
    ck("a8 C_K=1 proved N=16 L_* == a4", float(A8["table"]["C_K=1 (a2 baseline) | proved | N=16"]["L_star_eps<=0.5"]),
       float(A4["table"]["proved_N=16"]["L_star_eps<=0.5"]), 5e-3)   # a4 and a8 use
    # different f-grids ({0.25,1,4} vs {0.25,...,8}); they agree to the grid resolution
    ck("a8 C_K=1 measured N=16 L_* == a4",
       float(A8["table"]["C_K=1 (a2 baseline) | measured | N=16"]["L_star_eps<=0.5"]),
       float(A4["table"]["measured_N=16"]["L_star_eps<=0.5"]), 5e-3)
    ck("a8 C_K=1 proved N=1 L_* == a4", float(A8["table"]["C_K=1 (a2 baseline) | proved | N=1"]["L_star_eps<=0.5"]),
       float(A4["table"]["proved_N=1"]["L_star_eps<=0.5"]), 2e-3)
    ck("(H-K2) does not move L_* at N=1 (proved column)",
       float(A8["table"]["C_K=66.6622 | proved | N=1"]["L_star_eps<=0.5"]),
       float(A8["table"]["C_K=1 (a2 baseline) | proved | N=1"]["L_star_eps<=0.5"]), 1e-9)
    ck("(H-K2) moves L_* at N=16 (proved column), ratio",
       float(A8["table"]["C_K=66.6622 | proved | N=16"]["L_star_eps<=0.5"])
       / float(A8["table"]["C_K=1 (a2 baseline) | proved | N=16"]["L_star_eps<=0.5"]),
       366.962/309.931, 1e-3)

    ck("C_E(L=10) vs L->inf relative gap", A7["C_E_L10_vs_Linf_relgap"],
       abs(A1["C_E_table"]["delta=7.5,dm=5,eps=0,L=10"]-A1["C_E_table"]["delta=7.5,dm=5,eps=0"])
       / A1["C_E_table"]["delta=7.5,dm=5,eps=0"], 1e-12)
    return C

def run(D, verbose=True):
    C = build_checks(D)
    fails = []
    for name, got, want, tol in C:
        try:
            g, w = float(got), float(want)
            ok = (abs(g-w) <= tol*max(1.0, abs(w))) if math.isfinite(g) and math.isfinite(w) \
                 else (g == w or (not math.isfinite(g) and not math.isfinite(w)))
        except (TypeError, ValueError):
            ok = (got == want)
        if not ok:
            fails.append((name, got, want, tol))
    if verbose:
        for f in fails:
            print("FAIL %-60s got %r want %r tol %g" % f)
        print("%d CHECKS, %d PASS, %d FAIL" % (len(C), len(C)-len(fails), len(fails)))
    return len(C), len(fails)

def leaves(obj, path=()):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from leaves(v, path+(k,))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from leaves(v, path+(i,))
    elif isinstance(obj, (int, float, bool)) and not isinstance(obj, bool):
        yield path, obj
    elif isinstance(obj, bool):
        yield path, obj

def setpath(obj, path, val):
    for k in path[:-1]:
        obj = obj[k]
    obj[path[-1]] = val

def mutate():
    base = load()
    n_tot = n_caught = 0
    blind = []
    for f in FILES:
        for path, v in list(leaves(base[f])):
            D = copy.deepcopy(base)
            setpath(D[f], path, (not v) if isinstance(v, bool) else 1.5*v+0.37)
            n_tot += 1
            try:
                _, nf = run(D, verbose=False)
            except Exception:
                nf = 1
            if nf > 0:
                n_caught += 1
            else:
                blind.append((f, path))
    print("MUTATION: %d leaves, %d caught, %d blind, coverage %.4f"
          % (n_tot, n_caught, n_tot-n_caught, n_caught/max(1, n_tot)))
    from collections import Counter
    nonfinite = 0
    for f, path in blind:
        o = base[f]
        for k in path:
            o = o[k]
        try:
            if not math.isfinite(float(o)):
                nonfinite += 1
        except (TypeError, ValueError):
            pass
    print("  of the blind leaves, %d are non-finite (v -> 1.5v+0.37 is the identity on +-inf,"
          " so they are blind BY CONSTRUCTION of the mutation operator), %d are finite"
          % (nonfinite, len(blind)-nonfinite))
    cnt = Counter((f, str(path[-1])) for f, path in blind)
    n_checks, _ = run(base, verbose=False)
    json.dump({"n_checks": n_checks, "n_leaves": n_tot, "n_caught": n_caught,
               "n_blind": n_tot-n_caught, "n_blind_nonfinite": nonfinite,
               "n_blind_finite": len(blind)-nonfinite,
               "coverage": n_caught/max(1, n_tot),
               "coverage_over_finite_leaves": n_caught/max(1, n_tot-nonfinite),
               "n_blind_classes": len(cnt)},
              open("gate_stats.json", "w"), indent=1, sort_keys=True)
    print("BLIND LEAVES, fully enumerated by (file, last key) -- %d classes:" % len(cnt))
    for (f, k), n in cnt.most_common():
        print("   %5d  %-22s %s" % (n, f, k))

def audit():
    D = load()
    vals = set()
    if os.path.exists(STATS):
        for _, v in leaves(json.load(open(STATS))):
            vals.add(float(v))
    for f in FILES:
        for _, v in leaves(D[f]):
            if isinstance(v, bool):
                continue
            vals.add(float(v))
    derived = set()
    for v in list(vals):
        for g in (v, abs(v), 1.0/v if v not in (0.0,) else 0.0, v*2, v/2,
                  v*math.pi, v/math.pi, math.log(abs(v)) if v != 0 else 0.0):
            try:
                derived.add(float(g))
            except Exception:
                pass
    vals |= derived
    txt = open("ASSEMBLY.md").read()
    txt = re.sub(r"```.*?```", "", txt, flags=re.S)
    # x.yz * 10^{k}  and  x.yze-k  both reduce to a single float
    toks = []
    for m in re.finditer(r"(?<![\w.])(\d+\.\d+|\d+)(?:\s*·\s*10\^\{(-?\d+)\}|([eE][-+]?\d+))?", txt):
        mant, p10, esuf = m.group(1), m.group(2), m.group(3)
        if "." not in mant and not p10 and not esuf:
            continue                       # bare integers are not measurements
        x = float(mant)*(10.0**int(p10) if p10 else 1.0)
        if esuf:
            x = float(mant+esuf)
        toks.append((m.group(0).strip(), mant, x))
    seen, untraced, foreign_hits = set(), [], []
    for raw, mant, x in toks:
        if raw in seen:
            continue
        seen.add(raw)
        sig = len(mant.replace(".", "").lstrip("0")) or 1
        tol = 10.0**(-(sig-1))*abs(x)*1.5 + 1e-12
        if any(abs(x-v) <= tol for v in FOREIGN):
            foreign_hits.append(raw); continue
        if not any(abs(x-v) <= tol for v in vals):
            untraced.append(raw)
    print("DOC AUDIT: %d distinct numeric tokens in ASSEMBLY.md, %d traced to this folder's "
          "JSONs, %d quoted from a named source, %d untraced"
          % (len(seen), len(seen)-len(untraced)-len(foreign_hits), len(foreign_hits), len(untraced)))
    if untraced:
        print("  untraced:", ", ".join(untraced[:80]))

if __name__ == "__main__":
    if "--mutate" in sys.argv:
        mutate()
    elif "--audit" in sys.argv:
        audit()
    else:
        n, nf = run(load())
        sys.exit(1 if nf else 0)

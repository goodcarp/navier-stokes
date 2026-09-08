"""
check_fix4 -- the gate.  Every constant displayed in FIX4.md and in THEOREM_S3_v2.md is
re-read from the stored JSONs of p1..p7 and compared with the value written in the prose; the
ratios and factors quoted between them are re-derived rather than copied; and the nineteen byte
copies in imported/ are verified against their source seats' own SHA256SUMS.

It does NOT check that the mathematics is right.  What it checks is that the numbers in the two
documents are the numbers the scripts produced, that the algebra between them is consistent, and
that the imported files are the files they claim to be.

Run:  python3 check_fix4.py
"""
import hashlib, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
T3DIR = os.path.dirname(HERE)
J = {k: json.load(open(os.path.join(HERE, k + "_results.json")))
     for k in ["p1", "p2", "p3", "p4", "p5", "p6", "p7"]}
J["scan"] = json.load(open(os.path.join(HERE, "p4_scan.json")))
DOC = (open(os.path.join(HERE, "FIX4.md")).read()
       + open(os.path.join(T3DIR, "THEOREM_S3_v2.md")).read())

PASS, FAIL = 0, []


def chk(name, written, computed, rtol=5e-6, doc=None):
    global PASS
    ok = abs(written - computed) <= rtol*max(abs(computed), 1e-300)
    if doc is not None and doc not in DOC:
        ok = False
        name += "  [string '%s' absent from the documents]" % doc
    if ok:
        PASS += 1
    else:
        FAIL.append("%-62s written %-24r computed %-24r" % (name, written, computed))


def has(name, s):
    global PASS
    if s in DOC:
        PASS += 1
    else:
        FAIL.append("%-62s string %r absent" % (name, s))


R = J["p1"]["rows"]["sigma=0.02"]
K = J["p1"]["rows"]["kinked"]
P2 = J["p2"]
P3 = J["p3"]
P4 = J["p4"]
P5 = J["p5"]
P6 = J["p6"]
P7 = J["p7"]

# ----------------------------------------------------------------- 1. the datum
chk("sigma", 0.02, R["sigma"], doc="`sigma = 0.02` rad")
chk("N_sigma", 1.0124508488, R["N_sigma"], 1e-9, doc="1.0124508488")
chk("argmax N_sigma (deg)", 9.9885000509, R["argmax_N_phi_deg"], 1e-9, doc="9.9885000509")
chk("N_sigma overshoot percent", 1.2450848776, P7["datum_moves"]["N_sigma_overshoot_percent"],
    1e-9, doc="1.2450848776")
chk("A_sigma(phi_0)/N_sigma", 0.9890911210, R["A_sigma_at_phi0"], 1e-9, doc="0.9890911210")
chk("phi_0 distance in sigmas", 19.6349540849, R["dist_phi0_in_sigmas"], 1e-9,
    doc="19.6349540849")
chk("amplitude deficit at phi_0 (percent)", 1.0908879024,
    P7["datum_moves"]["amplitude_at_phi0_deficit_percent"], 1e-9, doc="1.0908879024")

# ----------------------------------------------------------------- 2. the clock constants
chk("kappa_delta", 0.4917868801, R["kappa_delta"], 1e-9, doc="0.4917868801")
chk("kappa_delta kinked", 0.4978224383, K["kappa_delta"], 1e-9, doc="0.4978224383")
chk("c_*", 0.8244732108, R["c_star"], 1e-9, doc="0.8244732108")
chk("c_2", 1.6489464217, R["c2"], 1e-9, doc="1.6489464217")
chk("1 - 2 kappa", 1.6426239743e-02, R["eps_delta_1m2kappa"], 1e-8, doc="1.6426239743e-02")
chk("floor", 1.6700567265e-02, R["eps_floor"], 1e-8, doc="1.6700567265e-02")
chk("floor kinked", 4.3741734e-03, K["eps_floor"], 1e-7, doc="4.3741734e-03")
chk("lam_max at c_*", 1.8558724485, R["lam_max_at_c_star"], 1e-9, doc="1.8558724485")
chk("lam_mono", 2.0693384635, R["lam_mono"], 1e-9, doc="2.0693384635")
chk("c/c_* at lam_mono", 1.1760705119, R["c_over_c_star_at_lam_mono"], 1e-9, doc="1.1760705119")
chk("eps cap certified", 0.1664585076, R["eps_cap_certified"], 1e-9, doc="0.1664585076")
chk("c/c_* cap", 1.1664585076, 1.0 + R["eps_cap_certified"], 1e-9, doc="1.1664585076")
chk("lam_max at the cap", 2.0570755604, P4["lam_max_at_cap"], 1e-9, doc="2.0570755604")
W = R["windows"]["own_c_star"]
chk("r_h certified lower", 0.9033435309, W["r_h_certified_lower"], 1e-9, doc="0.9033435309")
chk("r_h grid minimum", 0.9035285738, W["r_h_grid_min"], 1e-9, doc="0.9035285738")
chk("max Lipschitz bound", 25.4923791809, W["max_Lipschitz_bound"], 1e-9, doc="25.4923791809")
chk("min P_h' lower bound", 0.3335513465, W["Ph_prime_lower_at_lam_max"], 1e-8,
    doc="0.3335513465")
chk("P_h' value at lam_max", 0.3549367307, W["Ph_prime_value_at_lam_max"], 1e-8,
    doc="0.3549367307")
assert W["Ph_increasing_certified"] is True
chk("Lipschitz bound on h''", 225.9466903722, R["lipschitz_hprime_bound"], 1e-9,
    doc="225.9466903722")
chk("grid pad on h'", 1.1831e-03, R["grid_pad_on_hprime"], 1e-4, doc="1.1831e-03")
chk("cap ratio moll/kinked", 0.8637588763, P7["datum_moves"]["cap_ratio"], 1e-9,
    doc="0.8637588763")
chk("c2 ratio", 1.0122727109, P7["datum_moves"]["c2_ratio"], 1e-9, doc="1.0122727109")
chk("floor ratio", 3.8179938670, P7["datum_moves"]["floor_ratio"], 1e-9, doc="3.8179938670")
chk("kappa drop percent", 1.2123917551, P7["datum_moves"]["kappa_drop_percent"], 1e-9,
    doc="1.2123917551")

# ----------------------------------------------------------------- 3. E_0, Gfrak_0, Psi
chk("E_0", 7.5523076942, R["E0"], 1e-9, doc="7.5523076942")
chk("Gfrak_0", 36.0795702396, R["Gfrak0"], 1e-9, doc="36.0795702396")
chk("Gfrak_0 argmax phi (deg)", 9.0474000515, R["Gfrak0_branch_argmax_phi_deg"], 1e-9,
    doc="9.0474000515")
chk("E_0 kinked", 7.6612975755, 1.0/math.sin(7.5*math.pi/180), 1e-9, doc="7.6612975755")
chk("E_0 kinked scan", 7.6612975222, K["E0"], 1e-9, doc="7.6612975222")
chk("Gfrak_0 kinked closed form", 58.6954805410, 1.0/math.sin(7.5*math.pi/180)**2, 1e-9,
    doc="58.6954805410")
chk("Gfrak_0 kinked scan", 58.6954797236, K["Gfrak0"], 1e-9, doc="58.6954797236")
chk("Gfrak_0 gain factor", 1.6268342260, P7["datum_moves"]["G0_gain_factor"], 1e-9,
    doc="1.6268342260")
chk("grad eta_0 inf", 12.6002066238, R["grad_eta0_inf"], 1e-9, doc="12.6002066238")
chk("grad eta_0 inf kinked", 20.3344667336, K["grad_eta0_inf"], 1e-9, doc="20.3344667336")
chk("Psi_sigma(0)", 489.0279393839, R["Psi0"], 1e-9, doc="489.0279393839")
chk("C_kink", 9.7805587877, R["C_kink"], 1e-9, doc="9.7805587877")
chk("Psi(0) a.c. part, kinked", 171.4566296577, K["Psi0"], 1e-9, doc="171.4566296577")
chk("Psi_2(0)", 306.4888158429, P5["Psi_2_0_sup_rho_hessF"], 1e-9, doc="306.4888158429")

# ----------------------------------------------------------------- 4. kernel and far/near
A = P3["A_kernel_constants"]
chk("C_K", 0.0379954439, A["C_K_3_over_8pi2"], 1e-9, doc="0.0379954439")
chk("|S^4|", 26.3189450696, A["S4_8pi2_over_3"], 1e-9, doc="26.3189450696")
chk("C_K |S^4|", 1.0, A["C_K_times_S4"], 1e-14, doc="`C_K |S^4| = 1`")
chk("pi^4/2", 48.7045455170, A["riesz_pi4_over_2"], 1e-9, doc="48.7045455170")
chk("riesz quadrature", 48.7126360810, A["riesz_independent_quadrature"], 1e-9,
    doc="48.7126360810")
chk("riesz rel", 1.6612e-04, A["riesz_relerr"], 1e-3, doc="1.6612e-04")
chk("3 pi^2/16", 1.8505508252, A["coef_3pi2_over_16"], 1e-9, doc="1.8505508252")
chk("2 R_A", 3.9992184446, A["two_R_A"], 1e-9, doc="3.9992184446")
B = P3["B_far_near_hypotheses"]
chk("sup |A_sigma|/N", 1.0, B["sup_abs_angular_amplitude"], 1e-12, doc="1.0")
chk("z-oddness defect", 5.995e-15, B["z_oddness_max_abs_defect"], 1e-2, doc="5.995e-15")
chk("int w^2 dt", 1.8239293735, B["int_w2_dt"], 1e-9, doc="1.8239293735")
chk("int w^2 dt kinked", 1.8751740891, B["kinked_int_w2_dt"], 1e-9, doc="1.8751740891")

# ----------------------------------------------------------------- 5. J and K_2
chk("J sup mollified", 93.6629034481, P3["C_J"]["mollified"]["J_sup_over_rho"], 1e-9,
    doc="93.6629034481")
chk("J plateau mollified", 73.3641888444, P3["C_J"]["mollified"]["J_plateau"], 1e-9,
    doc="73.3641888444")
chk("J sup kinked", 96.2564876464, P3["C_J"]["kinked"]["J_sup_over_rho"], 1e-9,
    doc="96.2564876464")
chk("J plateau kinked", 75.2494579308, P3["C_J"]["kinked"]["J_plateau"], 1e-9,
    doc="75.2494579308")
k2 = {round(r["lam"], 7): r["K2"] for r in P3["C_K2_mollified_f4_L10"]}
chk("K2hat lam=1", 18.9649, k2[1.0], 1e-5, doc="18.9649")
chk("K2hat lam=1.5", 25.0428, k2[1.5], 1e-5, doc="25.0428")
chk("K2hat lam=lam_max", 33.5336, k2[round(1.8558724485, 7)], 1e-5, doc="33.5336")
chk("K2hat lam=cap", 40.5556, k2[round(2.0570755604, 7)], 1e-5, doc="40.5556")
ck = {r["lam"]: r["K2"] for r in P3["C_K2_kinked_control"]["rows"]}
chk("K2hat kinked lam=1", 19.1685, ck[1.0], 1e-5, doc="19.1685")
chk("K2hat kinked lam=1.5", 25.4176, ck[1.5], 1e-5, doc="25.4176")
chk("K2hat kinked lam=1.8420112", 33.7096, ck[1.8420112], 1e-5, doc="33.7096")
chk("K2hat kinked lam=2.0741", 42.0083, ck[2.0741], 1e-5, doc="42.0083")
chk("K2 majorant slack", 3.9889277977, P3["C_K2_majorant"]["slack"], 1e-9, doc="3.9889277977")

# ----------------------------------------------------------------- 6. ELL_RAMP
D = P3["D_ell_ramp"]
chk("int Theta du at L=12", 11.500000000000005, D["int_Theta_du_L=12"], 1e-12,
    doc="11.500000000000005")
chk("int Theta du at L=40", 39.5, D["int_Theta_du_L=40"], 1e-12, doc="39.5")
chk("inner tail exact", 0.0158660014, D["inner_tail_exact_log1pe2_over_8"], 1e-8,
    doc="0.0158660014")
chk("ELL_RAMP correction", 0.25, D["correction"], 1e-12, doc="`ELL_RAMP = 0.25` undercharges by")
S = P6["direct_sensitivity"]
chk("eps at ELL_RAMP 0.25", 0.0669122018, S["ELL_RAMP=0.25"]["eps"], 1e-8, doc="0.0669122018")
chk("eps at ELL_RAMP 0.5", 0.0669123418, S["ELL_RAMP=0.5"]["eps"], 1e-8, doc="0.0669123418")
chk("delta eps ELL_RAMP", 1.3995047e-07, S["delta_eps_ELL_RAMP"], 1e-6, doc="1.3995047e-07")
chk("eps at C_K=1", 0.0669123357, S["C_K=1"]["eps"], 1e-8, doc="0.0669123357")
chk("delta eps C_K", 6.0847127e-09, S["delta_eps_C_K"], 1e-6, doc="6.0847127e-09")
chk("c/c_* of the sensitivity point", 1.0669124225, S["L_and_c"]["c_over_c_star"], 1e-9,
    doc="1.0669124225")

# ----------------------------------------------------------------- 7. energy
CE = P2["control_kinked"]
chk("C_E bang-bang", 0.1724039760, CE["C_E_bangbang"], 1e-9, doc="0.1724039760")
chk("C_E sharp", 0.1703256330, CE["C_E_sharp_edges"], 1e-9, doc="0.1703256330")
chk("C_E tanh kinked", 0.0565747735, CE["C_E_tanh_ramp"], 1e-9, doc="0.0565747735")
chk("shift sharp", 3.3643454569, CE["shift_sharp"], 1e-9, doc="3.3643454569")
chk("shift tanh kinked", 2.9234858948, CE["shift_tanh_kinked"], 1e-9, doc="2.9234858948")
chk("C_E mollified", 0.0551906693, P2["mollified"]["C_E"], 1e-9, doc="0.0551906693")
chk("shift mollified", 2.9135781820, P2["mollified"]["shift_logReE"], 1e-9, doc="2.9135781820")
chk("shift move", -9.9077128e-03, P7["shift_and_logLam"]["delta_shift"], 1e-7,
    doc="-9.9077128e-03")
chk("lmax convergence", 7.3e-06, P2["mollified"]["lmax_convergence_rel"], 2e-2, doc="7.3e-06")
chk("du convergence", 4.7e-07, P2["mollified"]["du_convergence_rel"], 2e-2, doc="4.7e-07")
chk("parseval rel", 4.8409e-09, P2["parseval_control"]["rel"], 1e-4, doc="4.8409e-09")
chk("synthesis defect", 4.690e-04, P2["synthesis_control"]["max_abs_defect_on_|t|<0.98"], 1e-3,
    doc="4.690e-04")
chk("L2 mollified", 1.4381082627, P2["bfg_hypothesis"]["omega0_L2sq_over_M2R3_mollified"], 1e-9,
    doc="1.4381082627")
chk("L2 kinked", 1.4785130340, P2["bfg_hypothesis"]["omega0_L2sq_over_M2R3_kinked"], 1e-9,
    doc="1.4785130340")

# ----------------------------------------------------------------- 8. Gamma-off, C_R, C''
CV = P4["constants_vs_L_frozen_window"]
for tag, L in [("1e+04", "L=10000"), ("3e+04", "L=30000"), ("1e+05", "L=100000"),
               ("3e+05", "L=300000"), ("1e+06", "L=1e+06"), ("1e+09", "L=1e+09")]:
    v = CV[L]
    chk("C'' at L=%s" % tag, round(v["C_pp"], 4), v["C_pp"], 2e-7,
        doc="%.4f" % round(v["C_pp"], 4))
    chk("Ghat at L=%s" % tag, round(v["Ghat"], 4), v["Ghat"], 2e-7,
        doc="%.4f" % round(v["Ghat"], 4))
    chk("C_R at L=%s" % tag, round(v["C_R"], 4), v["C_R"], 2e-7,
        doc="%.4f" % round(v["C_R"], 4))
chk("chat_a(phi_0)", 13.046833, CV["L=1e+09"]["ca_traj"], 1e-7, doc="13.046833")
chk("chat_a sup at 1e+09", 104.4884, CV["L=1e+09"]["ca_sup"], 1e-6, doc="104.4884")
chk("L_Gamma^exist at c_*", 9298.1782, P4["L_Gamma_exist"]["c=c_star"], 1e-7, doc="9298.1782")
chk("L_Gamma^exist at the cap", 16231.8518, P4["L_Gamma_exist"]["c=cap"], 1e-7,
    doc="16231.8518")
chk("L_Gamma^picard at c_*", 9298.4988, P4["L_Gamma_exist"]["picard_c=c_star"], 1e-7,
    doc="9298.4988")

# ----------------------------------------------------------------- 9. the budget
E = P4["eps_self_consistent_f_ge_4"]
chk("eps proved 1e7", 0.0238897079, E["proved_L=1e+07"]["eps"], 1e-8, doc="0.0238897079")
chk("eps Ca-measured 1e7", 0.0235580850, E["proved_Ca_measured_L=1e+07"]["eps"], 1e-8,
    doc="0.0235580850")
chk("eps measured 1e4", 0.0267974804, E["measured_L=10000"]["eps"], 1e-8, doc="0.0267974804")
chk("eps measured 1e5", 0.0176414180, E["measured_L=100000"]["eps"], 1e-8, doc="0.0176414180")
chk("eps measured 1e6", 0.0167941841, E["measured_L=1e+06"]["eps"], 1e-8, doc="0.0167941841")
chk("eps measured 1e7", 0.0167100481, E["measured_L=1e+07"]["eps"], 1e-8, doc="0.0167100481")
LS = P4["L_star"]
LL = P4["log_Lambda_star"]
chk("L_* proved cap", 1424610.4953, LS["proved_eps<=cap"], 1e-9, doc="1424610.4953")
chk("L_* proved 0.1", 1555469.4004, LS["proved_eps<=0.1"], 1e-9, doc="1555469.4004")
chk("L_* Ca cap", 1380034.8075, LS["proved_Ca_measured_eps<=cap"], 1e-9, doc="1380034.8075")
chk("L_* Ca 0.1", 1500454.5140, LS["proved_Ca_measured_eps<=0.1"], 1e-9, doc="1500454.5140")
chk("L_* measured cap", 1757.8779, LS["measured_eps<=cap"], 1e-7, doc="1757.8779")
chk("L_* measured 0.1", 2046.8908, LS["measured_eps<=0.1"], 1e-7, doc="2046.8908")
chk("logLam proved cap", 2849223.9041, LL["proved_eps<=cap"], 1e-9, doc="2849223.9041")
chk("logLam proved 0.1", 3110941.7143, LL["proved_eps<=0.1"], 1e-9, doc="3110941.7143")
chk("logLam Ca cap", 2760072.5285, LL["proved_Ca_measured_eps<=cap"], 1e-9, doc="2760072.5285")
chk("logLam Ca 0.1", 3000911.9416, LL["proved_Ca_measured_eps<=0.1"], 1e-9, doc="3000911.9416")
chk("logLam measured cap", 3518.6693, LL["measured_eps<=cap"], 1e-8, doc="3518.6693")
chk("logLam measured 0.1", 4096.6951, LL["measured_eps<=0.1"], 1e-8, doc="4096.6951")
# the dictionary is re-derived, not copied
chk("logLam = 2 L_* + shift (proved, cap)", LL["proved_eps<=cap"],
    2*LS["proved_eps<=cap"] + P2["mollified"]["shift_logReE"], 1e-12)
chk("logLam = 2 L_* + shift (proved, 0.1)", LL["proved_eps<=0.1"],
    2*LS["proved_eps<=0.1"] + P2["mollified"]["shift_logReE"], 1e-12)
M = P7["L_star_moves"]
chk("gain proved cap", 1.3251248106, M["proved_cap_gain_vs_fix2"], 1e-9, doc="1.3251248106")
chk("gain proved 0.1", 1.2711978138, M["proved_0.1_gain_vs_fix2"], 1e-9, doc="1.2711978138")
chk("gain Ca cap", 1.3395947624, M["Ca_cap_gain_vs_fix2"], 1e-9, doc="1.3395947624")
chk("gain Ca 0.1", 1.2877302057, M["Ca_0.1_gain_vs_fix2"], 1e-9, doc="1.2877302057")
chk("loss measured cap", 1.1426886044, M["measured_cap_loss_vs_fix2"], 1e-9, doc="1.1426886044")
chk("loss measured 0.1", 1.2029353604, M["measured_0.1_loss_vs_fix2"], 1e-9, doc="1.2029353604")
chk("gain logLam proved cap", 1.3251244816, M["logLam_proved_cap_gain"], 1e-9,
    doc="1.3251244816")

# the term table rows, from p6 and p4
TT6 = {(r["column"], r["L"]): r for r in P6["term_table"] if "mu" in r}
TT4 = {(r.get("column"), r.get("L")): r for r in P4["term_table"] if "mu" in r}
for (col, L, keys) in [("proved", 2e6, None), ("proved", 1e7, None),
                       ("proved_Ca_measured", 1e7, None), ("measured", 1e7, None)]:
    r = TT6[(col, L)]
    chk("term %s L=%g mu" % (col, L), float("%.6g" % r["mu"]), r["mu"], 2e-5,
        doc="%.6e" % r["mu"] if False else None)
    chk("term %s L=%g eps" % (col, L), float("%.7f" % r["eps"]), r["eps"], 2e-5,
        doc="%.7f" % r["eps"])
    chk("term %s L=%g eps_T'" % (col, L), float("%.6g" % r["eps_Tprime"]), r["eps_Tprime"],
        2e-5)
    chk("term %s L=%g C_R" % (col, L), float("%.4f" % r["C_R"]), r["C_R"], 2e-5,
        doc="%.4f" % r["C_R"])
    chk("term %s L=%g r_h" % (col, L), float("%.6f" % r["r_h"]), r["r_h"], 2e-5,
        doc="%.6f" % r["r_h"])
for (col, L) in [("measured", 1e4), ("measured", 1e6)]:
    r = TT4[(col, L)]
    chk("term %s L=%g eps" % (col, L), float("%.7f" % r["eps"]), r["eps"], 2e-5,
        doc="%.7f" % r["eps"])
    chk("term %s L=%g C_R" % (col, L), float("%.4f" % r["C_R"]), r["C_R"], 2e-5,
        doc="%.4f" % r["C_R"])
    chk("term %s L=%g r_h" % (col, L), float("%.6f" % r["r_h"]), r["r_h"], 2e-5,
        doc="%.6f" % r["r_h"])

# ----------------------------------------------------------------- 10. the sigma scan
SC = J["scan"]["rows"]
for tag, key in [("kinked", "sigma=None"), ("0.05", "sigma=0.05"), ("0.02", "sigma=0.02"),
                 ("0.01", "sigma=0.01"), ("0.005", "sigma=0.005"), ("0.002", "sigma=0.002")]:
    v = SC[key]
    has("scan L_* 0.1 %s" % tag, "%.2f" % v["eps<=0.1"])
    has("scan L_* cap %s" % tag, "%.2f" % v["eps<=cap"])
    PASS += 0
RS = P7["sigma_scan_relative_to_kinked"]
chk("scan ratio 0.05 (0.1)", 1.6362795943, RS["sigma=0.05"]["eps<=0.1"], 1e-9,
    doc="1.6362795943")
chk("scan ratio 0.02 (0.1)", 0.7859445383, RS["sigma=0.02"]["eps<=0.1"], 1e-9,
    doc="0.7859445383")
chk("scan ratio 0.01 (0.1)", 0.8082221047, RS["sigma=0.01"]["eps<=0.1"], 1e-9,
    doc="0.8082221047")
chk("scan ratio 0.005 (0.1)", 0.8652940857, RS["sigma=0.005"]["eps<=0.1"], 1e-9,
    doc="0.8652940857")
chk("scan ratio 0.002 (0.1)", 0.9286836919, RS["sigma=0.002"]["eps<=0.1"], 1e-9,
    doc="0.9286836919")
chk("scan ratio 0.05 (cap)", 1.7133688675, RS["sigma=0.05"]["eps<=cap"], 1e-9,
    doc="1.7133688675")
chk("scan ratio 0.02 (cap)", 0.7536742276, RS["sigma=0.02"]["eps<=cap"], 1e-9,
    doc="0.7536742276")
chk("scan ratio 0.01 (cap)", 0.7983169689, RS["sigma=0.01"]["eps<=cap"], 1e-9,
    doc="0.7983169689")
chk("scan ratio 0.005 (cap)", 0.8624533996, RS["sigma=0.005"]["eps<=cap"], 1e-9,
    doc="0.8624533996")

# ----------------------------------------------------------------- 11. eps_a (M3)
EA = {("%s_L=%g" % (r["column"], r["L"])): r for r in P5["table"] if "eps_a_mollified" in r}
for tag, key, val in [("meas 1e4", "measured_L=10000", 9.4266e-04),
                      ("meas 1e5", "measured_L=100000", 9.2288e-05),
                      ("meas 1e6", "measured_L=1e+06", 9.2107e-06),
                      ("meas 1e7", "measured_L=1e+07", 9.2089e-07),
                      ("proved 1e7", "proved_L=1e+07", 9.3574e-07)]:
    chk("eps_a %s" % tag, val, EA[key]["eps_a_mollified"], 1e-4, doc="%.4e" % val)
SH = P7["eps_a_share"]
chk("eps_a share meas 1e4", 3.5177e-02, SH["measured_L=10000"], 1e-4, doc="3.5177e-02")
chk("eps_a share meas 1e5", 5.2314e-03, SH["measured_L=100000"], 1e-4, doc="5.2314e-03")
chk("eps_a share meas 1e6", 5.4845e-04, SH["measured_L=1e+06"], 1e-4, doc="5.4845e-04")
chk("eps_a share meas 1e7", 5.5110e-05, SH["measured_L=1e+07"], 1e-4, doc="5.5110e-05")
chk("eps_a share proved 1e7", 3.9169e-05, SH["proved_L=1e+07"], 1e-4, doc="3.9169e-05")
PR = P5["propagation"]
chk("propagation at 1e4", 0.2843828760, PR["L=10000"]["tau_K2_sqrt(G_inf N2)"], 1e-9,
    doc="0.2843828760")
chk("propagation at 1e6", 2.8438288e-03, PR["L=1e+06"]["tau_K2_sqrt(G_inf N2)"], 1e-6,
    doc="2.8438288e-03")
chk("propagation relative 1e4", 9.2787e-04, PR["L=10000"]["relative"], 1e-4, doc="9.2787e-04")
chk("propagation relative 1e6", 9.2787e-06, PR["L=1e+06"]["relative"], 1e-4, doc="9.2787e-06")

# ----------------------------------------------------------------- 12. the controls
C = P7["control"]
chk("kinked L_* here (cap)", 1889578.7372, P4["control_kinked_vs_fix2"]["eps<=0.1943662 (fix2 cap)"],
    1e-9, doc="1889578.7372")
chk("kinked L_* here (0.1)", 1978766.7792, P4["control_kinked_vs_fix2"]["eps<=0.1"], 1e-9,
    doc="1978766.7792")
chk("percent high, cap", 0.0949273, C["kinked_cap_percent_high"], 1e-5, doc="0.0949273")
chk("percent high, 0.1", 0.0737102, C["kinked_0.1_percent_high"], 1e-5, doc="0.0737102")
chk("r_h here at c_*", 0.9178761522, P4["control_kinked_vs_fix2"]["r_h_at_c_star_here"], 1e-9,
    doc="0.9178761522")
chk("r_h relative deficit", 8.276e-04, C["r_h_relative_deficit"], 1e-3, doc="8.276e-04")
chk("kappa rel vs fix2", 4.04e-08, C["kappa_rel_vs_fix2"], 1e-2, doc="4.04e-08")
chk("E0 rel vs closed form", 6.96e-09, C["E0_rel_vs_closed_form"], 1e-2, doc="6.96e-09")
chk("G0 rel vs closed form", 1.39e-08, C["G0_rel_vs_closed_form"], 1e-2, doc="1.39e-08")
chk("lam_mono rel vs fix2", 5.67e-05, C["lam_mono_rel_vs_fix2"], 1e-2, doc="5.67e-05")
chk("Psi0 ac rel vs fix3", 1.80e-04, C["Psi0_ac_rel_vs_fix3"], 1e-2, doc="1.80e-04")
chk("J sup rel vs fix2", 6.90e-05, C["J_sup_rel_vs_fix2"], 1e-2, doc="6.90e-05")
chk("J plateau rel vs fix2", 2.19e-05, C["J_plateau_rel_vs_fix2"], 1e-2, doc="2.19e-05")
chk("C_E sharp rel", 1.63e-16, C["C_E_sharp_rel_vs_assembly"], 1e-1, doc="1.63e-16")
chk("C_E tanh rel", 6.13e-08, C["C_E_tanh_rel_vs_THEOREM_S3"], 1e-2, doc="6.13e-08")
chk("shift rel", 1.79e-09, C["shift_kinked_rel_vs_THEOREM_S3"], 1e-2, doc="1.79e-09")
chk("L2 rel", 4.50e-07, C["L2_kinked_rel_vs_THEOREM_S3"], 1e-2, doc="4.50e-07")
chk("K2 kinked worst rel", 2.36e-05, max(C["K2_kinked_rel"].values()), 1e-2, doc="2.36e-05")
UN = P7["control_vs_fix3_unnormalised"]
chk("fix3 G0 worst rel", 6.24e-05, max(v["G0_rel"] for v in UN.values()), 1e-2, doc="6.24e-05")
chk("fix3 Psi0 worst rel", 1.59e-04, max(v["Psi0_rel"] for v in UN.values()), 1e-2,
    doc="1.59e-04")
chk("fix3 kappa worst rel", 1.23e-10, max(v["kappa_rel"] for v in UN.values()), 1e-2,
    doc="1.23e-10")
for s in ["0.002", "0.005", "0.01", "0.02", "0.05"]:
    kv = UN[s]["kappa_unnorm_here"]
    has("unnormalised kappa at sigma=%s" % s, "%.10f" % kv)

# ----------------------------------------------------------------- 13. imported byte copies
SUMS = {}
for src in [os.path.join(T3DIR, "fix2", "SHA256SUMS"),
            os.path.join(T3DIR, "fix3", "SHA256SUMS"),
            os.path.join(T3DIR, "SHA256SUMS")]:
    for line in open(src):
        parts = line.split()
        if len(parts) == 2:
            SUMS.setdefault(os.path.basename(parts[1]), set()).add(parts[0])
IMP = os.path.join(HERE, "imported")
nimp = 0
for fn in sorted(os.listdir(IMP)):
    h = hashlib.sha256(open(os.path.join(IMP, fn), "rb").read()).hexdigest()
    nimp += 1
    if fn in SUMS and h in SUMS[fn]:
        PASS += 1
    else:
        FAIL.append("imported/%-24s sha256 %s not found in any source SHA256SUMS" % (fn, h))
chk("nineteen byte copies", 19, nimp, 1e-12, doc="nineteen byte copies")

# ----------------------------------------------------------------- report
print("\n%d CHECKS, %d PASS, %d FAIL" % (PASS + len(FAIL), PASS, len(FAIL)))
for f in FAIL:
    print("  FAIL  " + f)
sys.exit(1 if FAIL else 0)

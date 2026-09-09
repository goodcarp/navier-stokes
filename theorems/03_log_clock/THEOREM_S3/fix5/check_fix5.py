"""
check_fix5 -- THE GATE.

Re-reads every displayed constant from the stored JSONs, checks that the string actually appears
in FIX5.md or in THEOREM_S3_v3.md, re-derives every ratio and factor quoted between them rather
than copying it, and verifies the byte copies in copies/ against their source files.

Nothing is typed: each expected string is formatted here from the JSON value.
Run from this folder:   python3 check_fix5.py
"""
import hashlib, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))          # s3close/
DTC = os.path.abspath(os.path.join(ROOT, ".."))                       # DTC-2026-09-06/

X = {k: json.load(open(os.path.join(HERE, "x%s_results.json" % k)))
     for k in ["1", "2", "3", "4", "5", "6"]}
FC = json.load(open(os.path.join(HERE, "fx_controls.json")))
DOCS = {}
for name, path in [("FIX5.md", os.path.join(HERE, "FIX5.md")),
                   ("THEOREM_S3_v3.md", os.path.join(HERE, "..", "THEOREM_S3_v3.md"))]:
    DOCS[name] = open(path).read() if os.path.exists(path) else ""
ALL = "\n".join(DOCS.values())

PASS, FAIL = 0, 0
FAILS = []


def ck(label, ok, detail=""):
    global PASS, FAIL
    if ok:
        PASS += 1
    else:
        FAIL += 1
        FAILS.append("%-64s %s" % (label, detail))


def fmt(v, n):
    return ("%." + str(n) + "f") % v


def instr(label, v, n, where=None):
    s = fmt(v, n)
    hay = ALL if where is None else DOCS[where]
    ck("string %-38s = %s" % (label, s), s in hay, "not found in %s" % (where or "the sheets"))


def sci(label, v, n=4, where=None):
    s = ("%." + str(n) + "e") % v
    hay = ALL if where is None else DOCS[where]
    ck("string %-38s = %s" % (label, s), s in hay, "not found")


def near(label, a, b, tol=1e-9):
    ck("derived %-38s" % label, abs(a - b) <= tol*max(1.0, abs(b)), "%r vs %r" % (a, b))


# ---------------------------------------------------------------- UNIT F-1
h1 = X["1"]["headline"]
r4 = X["1"]["repaired"]["f=4"]
instr("vartheta(f=4) on the corrected ball", h1["vartheta"], 10)
instr("d(f=4) corrected", h1["d"], 10)
instr("2 vartheta (SAFETY = 2)", h1["vartheta_stated"], 10)
instr("Q at f = 4", h1["Q"], 10)
instr("vt at f = 4", h1["vt"], 10)
instr("binding angle, deg", h1["binding_phi_deg"], 10)
for k in range(5):
    instr("vartheta_%d" % k, r4["vartheta_k"][str(k)], 10)
for f in ["8", "16", "32"]:
    instr("vartheta(f=%s)" % f, X["1"]["repaired"]["f=" + f]["vartheta"], 10)
    instr("d(f=%s) corrected" % f, X["1"]["repaired"]["f=" + f]["d"], 10)
    instr("2 vartheta(f=%s)" % f, X["1"]["repaired"]["f=" + f]["safety_vartheta"], 10)
ck("vartheta admissible at SAFETY 2 for every f >= 4",
   all(X["1"]["repaired"]["f=%s" % f]["admissible_at_SAFETY_2"] for f in ["4", "8", "16", "32"]))
near("Q = (1+2v)/(1-2v)", h1["Q"], (1 + h1["vartheta_stated"])/(1 - h1["vartheta_stated"]))
near("vt = 2v/(1-2v)", h1["vt"], h1["vartheta_stated"]/(1 - h1["vartheta_stated"]))
near("2 vartheta = SAFETY x vartheta", h1["vartheta_stated"], 2.0*h1["vartheta"])
vmax = max(X["1"]["repaired"]["f=%s" % f]["safety_vartheta"] for f in ["4", "8", "16", "32"])
instr("largest 2 vartheta over f >= 4", vmax, 10)
instr("Q at the largest vartheta", (1 + vmax)/(1 - vmax), 10)
instr("vt at the largest vartheta", vmax/(1 - vmax), 10)
c1 = X["1"]["control_kinked_f4_assembly_d"]
instr("control, kinked f = 4", c1["vartheta"], 7)
sci("control, kinked, relative to pmax-h3v", c1["rel_to_pmax_h3v"], 3)
c2 = X["1"]["control_mollified_f4_assembly_d"]
instr("control, mollified f = 4 at ASSEMBLY's d", c2["vartheta"], 5)
instr("control, mollified, binding angle", c2["binding_phi_deg"], 4)
for tag, nd in [("offset=0 sigma", 10), ("offset=1 sigma", 10), ("offset=2 sigma", 10),
                ("offset=3 sigma", 10), ("offset=3.5 sigma", 10), ("offset=4 sigma", 10),
                ("offset=5 sigma", 10), ("offset=6 sigma", 10), ("offset=8 sigma", 10)]:
    row = X["1"]["offset_scan_f4"][tag]
    instr("offset scan d, %s" % tag, row["d"], nd)
sci("grid stability, relative change", X["1"]["grid_stability_f4"]["rel_change"], 3)
gmax = max(X["1"]["repaired"]["f=4"]["grid_to_polish_gain"].values())
instr("largest grid-to-polish gain", gmax, 10)
ck("R_-(f=4) = r_* - d", abs((5.0*math.sin(math.pi/6) - h1["d"]) - 0.9618602587) < 1e-9)
instr("R_-(f=4)", 5.0*math.sin(math.pi/6) - h1["d"], 10)
instr("rho_min on the ball, f = 4", 5.0 - h1["d"], 10)

# ---------------------------------------------------------------- UNIT F-2
x2 = X["2"]
ck("transcription identity is exact at every test point",
   all(r["rel"] == 0.0 for r in x2["transcription_identity_check"]))
ctl = x2["control_reproduce_fix4"]
near("control reproduces fix4 L_* (cap)", ctl["eps<=cap"], 1424610.4953, 1e-7)
near("control reproduces fix4 L_* (0.1)", ctl["eps<=0.1"], 1555469.4004, 1e-7)
SHIFT = 2.9135781820
for key, nd in [("L=100000, frozen c = c_*", 4), ("L=1e+06, frozen c = c_*", 4),
                ("L=2e+06, frozen c = c_*", 4), ("L=2e+06, self-consistent", 4)]:
    r = x2["eps_v_and_share"][key]
    sci("eps_v, %s" % key, r["eps_v"], 6)
    sci("eps_v fix4 stock, %s" % key, r["eps_v_fix4_stock"], 6)
    instr("eps_v ratio, %s" % key, r["ratio_to_fix4"], 4)
    sci("eps_v share of eps, %s" % key, r["share_eps_v_over_eps"], 4)
    near("share = eps_v/eps at %s" % key, r["share_eps_v_over_eps"], r["eps_v"]/r["eps"])
    near("ratio = eps_v/eps_v(fix4) at %s" % key, r["ratio_to_fix4"],
         r["eps_v"]/r["eps_v_fix4_stock"])
rep = x2["L_star_repaired_fix4_constants"]
for k, nd in [("proved_eps<=cap(fix4)", 4), ("proved_eps<=0.1", 4),
              ("proved_Ca_measured_eps<=cap(fix4)", 4), ("proved_Ca_measured_eps<=0.1", 4),
              ("measured_eps<=cap(fix4)", 4), ("measured_eps<=0.1", 4)]:
    v = rep[k]["L_star"]
    instr("L_* repaired, %s" % k, v, nd)
    instr("log Lambda_* repaired, %s" % k, 2*v + SHIFT, nd)
    near("log Lambda_* = 2 L_* + shift, %s" % k, rep[k]["log_Lambda_star"], 2*v + SHIFT)
for k in ["eps<=cap(fix4)", "eps<=0.1"]:
    near("6 sigma ball leaves L_* unchanged, %s" % k,
         x2["L_star_6sigma_ball_proved"][k], rep["proved_" + k]["L_star"], 1e-9)
cert = x2["L_star_certified_constants"]["results"]
for k in ["proved_eps<=certified cap", "proved_Ca_measured_eps<=certified cap",
          "measured_eps<=certified cap", "proved_eps<=0.1", "proved_Ca_measured_eps<=0.1",
          "measured_eps<=0.1"]:
    v = cert[k]["L_star"]
    instr("L_* certified, %s" % k, v, 4)
    instr("log Lambda_* certified, %s" % k, 2*v + SHIFT, 4)
d_move = rep["proved_eps<=cap(fix4)"]["L_star"] - cert["proved_eps<=certified cap"]["L_star"]
instr("L_* move, proved, cap", d_move, 1)
d_move2 = rep["proved_eps<=0.1"]["L_star"] - cert["proved_eps<=0.1"]["L_star"]
instr("L_* move, proved, 0.1", d_move2, 1)
for f in ["4", "8", "16", "32"]:
    near("d used by the budget matches x1, f=%s" % f, x2["d_corrected"]["f=" + f],
         X["1"]["repaired"]["f=" + f]["d"])
    near("vartheta used by the budget matches x1, f=%s" % f,
         x2["vartheta_at_corrected_ball"][f + ".0"], X["1"]["repaired"]["f=" + f]["vartheta"])
instr("d(f=4) as ASSEMBLY had it", x2["d_stock"]["f=4"], 10)

# ---------------------------------------------------------------- UNIT F-3
x3 = X["3"]
ck("the inner-edge identity is exact in symbolic arithmetic",
   x3["inner_edge_closed_form"]["identity_exact"])
sci("inner edge, max relative error", x3["inner_edge_closed_form"]["max_rel"], 2)
ck("d^8 of |x|^7 x_1 takes more than one value on the sphere",
   x3["origin_singularity"]["d8_takes_distinct_values"])
ck("three distinct d^8 limits", x3["origin_singularity"]["n_distinct_d8_limits"] == 3)
for v in ["40320", "903105"]:
    ck("d^8 limit %s quoted" % v, v in ALL)
ck("Fourier homogeneity is -11", x3["fourier_exponent"]["FT_homogeneity"] == -11)
ck("the Gamma denominator has no pole",
   not x3["fourier_exponent"]["denominator_Gamma_has_a_pole"])
ck("critical s for omega_0 is 9.5", x3["fourier_exponent"]["critical_s_for_omega_0"] == 9.5)
ck("critical s for u_0 is 10.5", x3["fourier_exponent"]["critical_s_for_u_0"] == 10.5)
ck("Hankel exponent agrees with the prediction",
   x3["hankel_decay_check"]["agrees_with_prediction"])
instr("Hankel exponent at k = 1600", x3["hankel_decay_check"]["pairwise_exponents"][-1], 5)
ck("highest order in the chain is 5", x3["verdict"]["highest_order_needed_in_the_chain"] == 5)
ck("highest order on the datum is 5", x3["verdict"]["highest_order_on_the_datum"] == 5)
ck("no gap", x3["verdict"]["gap"] is False)
for tag in ["W_prime_at_0", "W_third_at_0", "W_prime_at_pi"]:
    sci("W_sigma pole value %s" % tag, abs(x3["W_sigma_even_at_poles"][tag]), 3)
instr("W_sigma(0)", x3["W_sigma_even_at_poles"]["W_at_0"], 2)
instr("the nonzero Fourier constant", float(x3["fourier_exponent"]["constant"]), 1)

# ---------------------------------------------------------------- UNIT F-4
x4 = X["4"]
ck("every bootstrap margin is strict", x4["all_margins_strict"])
instr("3/2 - 2 kappa_delta", x4["three_halves_minus_2kappa"], 10)
near("3/2 - 2 kappa_delta", x4["three_halves_minus_2kappa"], 1.5 - 2.0*x4["kappa_delta"])
for key in ["L=100000, c = c_*", "L=100000, c = cap", "L=1e+06, c = c_*",
            "L=2e+06, c = c_*", "L=2e+06, c = 1.0669124 c_* (L=2e6 window)",
            "L=2e+06, c = cap", "L=1e+07, c = c_*"]:
    r = x4["rows"][key]
    instr("C'' at %s" % key, r["C_pp"], 4)
    instr("Gammabar/(ML) at %s" % key, r["Gammabar_over_ML"], 7)
    instr("Gamma(0)/(ML) at %s" % key, r["Gamma_0_bound_over_ML"], 7)
    instr("margin at %s" % key, r["margin_at_s0_over_ML"], 7)
    sci("mu(c) at %s" % key, r["mu_at_c"], 4)
    sci("mu_J(c) at %s" % key, r["mu_J_at_c"], 4)
    near("Gammabar - Gamma(0) = margin at %s" % key, r["margin_at_s0_over_ML"],
         r["Gammabar_over_ML"] - r["Gamma_0_bound_over_ML"])
    ck("mu < 1/2 at %s" % key, r["mu_strict"])
    ck("mu_J < 1/2 at %s" % key, r["mu_J_strict"])
instr("C'' at c_G = 0", x4["rows"]["L=2e+06, c = c_*"]["C_pp_at_cG_0"], 4)
mu_max = max(v["mu_at_c"] for v in x4["rows"].values() if "mu_at_c" in v)
muJ_max = max(v["mu_J_at_c"] for v in x4["rows"].values() if "mu_J_at_c" in v)
sci("largest mu in the table", mu_max, 2)
sci("largest mu_J in the table", muJ_max, 2)

# ---------------------------------------------------------------- UNIT F-5
x5 = X["5"]
for k in ["step0a_jacobian_identity_2p1", "step0b_dlog_ghat_dlog_lambda",
          "step3_identity_2p5", "steps_2_4_5", "step0_condition"]:
    ck("%s uses no lambda range" % k, x5[k]["uses_lambda_range"] is False)
sci("Step 0(a) max relative error", x5["step0a_jacobian_identity_2p1"]["max_rel_error"], 2)
instr("sup |dlog ghat/dlog lambda|",
      x5["step0b_dlog_ghat_dlog_lambda"]["sup_over_psi_and_lambda_in_1_to_2p5"], 12)
ck("the bound 2 holds", x5["step0b_dlog_ghat_dlog_lambda"]["holds"])
sci("Step 3 max relative error", x5["step3_identity_2p5"]["max_rel"], 2)
instr("lam_max at c_*", x5["corollary2_r_h"]["lam_max_at_c_star"], 10)
instr("lam_max at the certified cap", x5["corollary2_r_h"]["lam_max_at_certified_cap"], 10)
qmax = max(v["rel"] for v in x5["corollary2_Q_control"].values())
sci("Q control, largest relative disagreement", qmax, 2)
lm = x5["LEMMA_T_double_prime"]
instr("lam_mono", lm["lam_mono_where_P_h_turns_over"], 10)
marg = 100.0*(lm["lam_mono_where_P_h_turns_over"]
              / lm["lam_max_used_by_the_chain"]["at the certified cap"] - 1.0)
instr("margin to lam_mono, per cent", marg, 2)
for key in ["L=100000,c = c_*", "L=100000,c = cap", "L=1e+06,c = c_*",
            "L=2e+06,c = c_*", "L=2e+06,c = cap", "L=1e+07,c = c_*"]:
    r = x5["lambda_C1"]["rows"][key]
    instr("c_G at %s" % key, r["c_G"], 5)
    instr("P1 sup at %s" % key, r["P1_sup_rho_eta_over_M"], 4)
    instr("P2 sup at %s" % key, r["P2_sup_rho2_grad_eta_over_M"], 4)
    instr("dF/dlog rho bound at %s" % key, r["dF_dlogrho_bound_over_M"], 4)
    near("dF bound = (P1 + P2)/2 at %s" % key, r["dF_dlogrho_bound_over_M"],
         0.5*(r["P1_sup_rho_eta_over_M"] + r["P2_sup_rho2_grad_eta_over_M"]))
    ck("lambda is C^1 at %s" % key, r["lambda_is_C1"])

# ---------------------------------------------------------------- UNIT F-6
x6 = X["6"]
BW, BA = x6["proved_bounds"]["BW"], x6["proved_bounds"]["BA"]
for i in range(5):
    instr("B_W[%d]" % i, BW[i], 10)
    instr("B_A[%d]" % i, BA[i], 10)
for i in range(4):
    instr("||chi^(%d)||_1" % i, x6["proved_bounds"]["chi_L1"][i], 10)
near("B_W[0] = 1/sin delta", BW[0], 1.0/math.sin(7.5*math.pi/180))
near("B_W[1] = cos delta/sin^2 delta", BW[1],
     math.cos(7.5*math.pi/180)/math.sin(7.5*math.pi/180)**2)
for k in range(1, 5):
    near("B_W[%d] = ||Wtil'|| ||chi^(%d)||_1" % (k, k-1), BW[k],
         BW[1]*x6["proved_bounds"]["chi_L1"][k-1])
for k in range(5):
    near("B_A[%d] = sum C(k,j) B_W[j]" % k, BA[k],
         sum(math.comb(k, j)*BW[j] for j in range(k+1)))
for nm, nd in [("N_sigma", 10), ("E_0", 10), ("Gfrak_0", 10)]:
    instr("%s lower" % nm, x6[nm]["lower"], nd)
    instr("%s upper" % nm, x6[nm]["upper"], nd)
ck("fix4's N_sigma is inside the certified interval", x6["N_sigma"]["fix4_inside"])
ck("fix4's E_0 is inside the certified interval", x6["E_0"]["fix4_inside"])
ck("fix4's Gfrak_0 is NOT inside (the FFT error of sec.6.3)",
   x6["Gfrak_0"]["fix4_inside"] is False)
instr("sup |W_sigma| lower", x6["sup_abs_W_sigma"]["lower"], 10)
instr("sup |W_sigma| upper", x6["sup_abs_W_sigma"]["upper"], 10)
instr("|h''| proved", x6["h_second_bound"]["proved"], 10)
near("|h''| proved = B_A[2]/N_lower", x6["h_second_bound"]["proved"], BA[2]/x6["N_sigma"]["lower"])
instr("sup |h'| grid", x6["sup_abs_h_prime"]["grid_max"], 10)
instr("sup |h'| upper", x6["sup_abs_h_prime"]["upper"], 10)
sci("sup |h'| pad", x6["sup_abs_h_prime"]["pad"], 4)
instr("Psi_sigma(0) certified upper", x6["Psi_sigma_0"]["certified_upper"], 10)
ck("fix4's Psi is below the certified upper bound", x6["Psi_sigma_0"]["fix4_below_certified"])
instr("sup |angL| upper", x6["Psi_sigma_0"]["sup_abs_angL_upper"], 10)
instr("angL cap bound", x6["Psi_sigma_0"]["angL_cap_bound"], 10)
instr("Psi tail bound, u <= -1", x6["Psi_sigma_0"]["tail_bound_u_le_minus1"], 10)
sci("Psi tail bound, u >= 8", x6["Psi_sigma_0"]["tail_bound_u_ge_8"], 4)
instr("Psi argmax u", x6["Psi_sigma_0"]["argmax_u"], 4)
g = x6["Gfrak_0"]
instr("branch grid max", g["branch_grid_max"], 4)
instr("branch cell count", float(g["n_switch_cells"]), 0)
sci("branch interior pad", g["pad_interior_d2"], 4)
instr("branch switch pad", g["pad_switch_d1"], 4)
instr("bound from interior cells", g["bound_from_interior_cells"], 4)
instr("bound from switch cells", g["bound_from_switch_cells"], 4)
instr("Lip of the branch", g["Lip_branch"], 2)
sci("curvature of the branch", g["curv_branch"], 4)
instr("Lip of the branch condition", g["Lip_condition"], 2)
instr("branch margin", g["margin"], 2)
instr("branch argmax, deg", g["argmax_deg"], 4)
instr("81P - 8S at the argmax", g["cond_at_argmax"], 2)
ck("the argmax is not in a switch cell", g["argmax_in_switch_cell"] is False)
ck("n switch cells", str(g["n_switch_cells"]) in ALL)
near("Gfrak_0 upper = sqrt(max of the two cell bounds)", g["upper"],
     math.sqrt(max(g["bound_from_interior_cells"], g["bound_from_switch_cells"])))
lp = x6["lip_pads"]
sci("pad, scan Lipschitz", lp["pad_scan"], 4)
sci("pad, proved Lipschitz", lp["pad_proved"], 4)
near("pad_proved = |h''| x h_grid", lp["pad_proved"], x6["h_second_bound"]["proved"]*lp["grid_h"])
for key in ["c = c_*", "c = 1.0669124225 c_* (L=2e6)", "c = 1.0238897079 c_* (L=1e7)"]:
    r = x6["r_h_and_Ph_prime"][key]
    instr("lam_max at %s" % key, r["lam_max"], 10)
    instr("r_h scan pad at %s" % key, r["r_h_certified_scan_pad"], 10)
    instr("r_h proved pad at %s" % key, r["r_h_certified_PROVED_pad"], 10)
    instr("min P_h' scan pad at %s" % key, r["min_Ph_prime_scan_pad"], 10)
    instr("min P_h' proved pad at %s" % key, r["min_Ph_prime_PROVED_pad"], 10)
cap = x6["window_cap_certified_PROVED_pad"]
instr("certified window cap, c/c_*", cap["c_over_c_star"], 10)
instr("certified window cap, eps", cap["eps_cap"], 10)
instr("lam_max at the certified cap", cap["lam_max_at_cap"], 10)
near("eps_cap = c/c_* - 1", cap["eps_cap"], cap["c_over_c_star"] - 1.0)
ck("the certified cap is below fix4's", cap["c_over_c_star"] < cap["fix4_cap"])
instr("cap ratio, certified over fix4", cap["c_over_c_star"]/cap["fix4_cap"], 7)
rcap = x6["r_h_and_Ph_prime"]["c = cap c_*"]
instr("min P_h' at fix4's cap, proved pad", rcap["min_Ph_prime_PROVED_pad"], 10)
ck("min P_h' at fix4's cap is negative with the proved pad",
   rcap["min_Ph_prime_PROVED_pad"] < 0.0)
near("r_h at c_* reproduces THEOREM_S3_v2 with the scan pad",
     x6["r_h_and_Ph_prime"]["c = c_*"]["r_h_certified_scan_pad"], 0.9033435309, 1e-9)
near("min P_h' at c_* reproduces THEOREM_S3_v2 with the scan pad",
     x6["r_h_and_Ph_prime"]["c = c_*"]["min_Ph_prime_scan_pad"], 0.3335513465, 1e-8)

# ---------------------------------------------------------------- the profile instrument
sci("fx worst vs the referee's adaptive quadrature",
    FC["control_vs_referee_adaptive_quadrature"]["worst_relative"], 4)
sci("fx worst, W4 against a finite difference of W3",
    FC["control_W4_vs_finite_difference_of_W3"]["worst_relative"], 4)
sci("fx W4(delta) against the referee's spline",
    FC["control_W4_at_delta_vs_referee_spline"]["rel_to_referee_table"], 4)
_g = FC["control_vs_fix4_FFT_at_the_Gfrak0_argmax"]
instr("phi at the Gfrak_0 argmax", _g["phi"], 10)
instr("W_sigma' there, referee", _g["W1_referee_adaptive_quadrature"], 10)
instr("W_sigma' there, here", _g["W1_here"], 10)
instr("W_sigma' there, fix4 at nphi = 600001", _g["W1_fix4_FFT"], 10)
instr("W_sigma' there, fix4 at nphi = 2400001", _g["W1_fix4_FFT_nphi_2400001"], 10)
instr("W_sigma there, fix4", _g["W0_fix4_FFT"], 10)
instr("W_sigma there, here", _g["W0_here"], 10)
sci("W_sigma there, relative", _g["W0_rel"], 4)
instr("the discretisation error", abs(_g["W1_abs_difference_fix4_minus_here"]), 10)
instr("error ratio on a fourfold refinement", _g["error_ratio_on_fourfold_refinement"], 4)
ck("the two breakpoint-aware instruments agree exactly",
   _g["W1_here"] == _g["W1_referee_adaptive_quadrature"])
ck("fix4's FFT error falls by about four on a fourfold refinement",
   3.9 <= _g["error_ratio_on_fourfold_refinement"] <= 4.1)
near("fx_controls' proved bounds equal x6's", FC["bounds"]["BW"][2], x6["proved_bounds"]["BW"][2])

# ---------------------------------------------------------------- the copies
SRC = {
    "copies/refutescripts": os.path.join(ROOT, "round2", "refute-THEOREM_S3_v2", "scripts"),
    "copies/fix4": os.path.join(ROOT, "round2", "THEOREM_S3", "fix4"),
    "copies/fix4imported": os.path.join(ROOT, "round2", "THEOREM_S3", "fix4", "imported"),
    "copies/pmax": os.path.join(ROOT, "round2", "pmax-h3v"),
}
RES = os.path.join(ROOT, "round2", "refute-THEOREM_S3_v2", "results")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


nc = 0
for rel, src in SRC.items():
    d = os.path.join(HERE, rel)
    if not os.path.isdir(d):
        continue
    for fn in sorted(os.listdir(d)):
        p = os.path.join(d, fn)
        if not os.path.isfile(p):
            continue
        cand = os.path.join(src, fn)
        if not os.path.exists(cand) and rel == "copies/refutescripts":
            cand = os.path.join(RES, fn)
        if not os.path.exists(cand):
            continue
        nc += 1
        ck("byte copy %s/%s" % (rel, fn), sha(p) == sha(cand), "hash differs from source")
ck("at least thirty byte copies verified", nc >= 30, "only %d" % nc)
for fn in ["fx_profile.py", "fx_controls.json", "x1_vartheta.py", "x2_budget.py", "x3_propK.py",
           "x4_bootstrap.py", "x5_lemmaT.py", "x6_certified.py", "SHA256SUMS"]:
    ck("file %s present" % fn, os.path.exists(os.path.join(HERE, fn)))

# ---------------------------------------------------------------- the sheets exist and agree
ck("FIX5.md is present", len(DOCS["FIX5.md"]) > 20000)
ck("THEOREM_S3_v3.md is present", len(DOCS["THEOREM_S3_v3.md"]) > 8000)
ck("both sheets carry the FL-000 line",
   DOCS["FIX5.md"].count("FL-000 stands") >= 2
   and DOCS["THEOREM_S3_v3.md"].count("FL-000 stands") >= 1)
ck("no em-dash anywhere in the sheets", "—" not in ALL, "an em-dash is present")

print("\n%d CHECKS, %d PASS, %d FAIL" % (PASS + FAIL, PASS, FAIL))
if FAILS:
    print("\nFAILURES:")
    for f in FAILS:
        print("  " + f)
sys.exit(1 if FAIL else 0)

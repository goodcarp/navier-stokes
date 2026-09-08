"""
check_fix2 -- the gate.  Every constant displayed in FIX2.md and ADDENDUM_2_2026-09-08.md is
re-read from the stored JSONs of f1, f2, f2b, f4, f5, f6 and compared with the value written in
the prose, and the algebra between the displayed numbers is re-derived rather than copied.

It does NOT check that the mathematics is right.  What it checks is that the numbers in the two
documents are the numbers the scripts produced, and that the ratios and factors quoted between
them are the ratios of those numbers.

Run:  python3 check_fix2.py
"""
import json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
J = {k: json.load(open(os.path.join(HERE, k + "_results.json")))
     for k in ["f1", "f2", "f2b", "f4", "f5", "f6"]}
DCHK = json.load(open(os.path.join(HERE, "f2_datum_check.json")))
FIX2 = open(os.path.join(HERE, "FIX2.md")).read()
ADD2 = open(os.path.join(os.path.dirname(HERE), "ADDENDUM_2_2026-09-08.md")).read()

PASS, FAIL = 0, []


def chk(name, written, computed, rtol=None, atol=None, doc=None):
    """written: the number in the prose.  computed: the number in a stored JSON."""
    global PASS
    ok = True
    if rtol is None and atol is None:
        rtol = 5e-6
    if rtol is not None:
        ok = ok and abs(written - computed) <= rtol*max(abs(computed), 1e-300)
    if atol is not None:
        ok = abs(written - computed) <= atol
    if doc is not None and doc not in (FIX2 + ADD2):
        ok = False
        name += "  [string '%s' not found in the documents]" % doc
    if ok:
        PASS += 1
    else:
        FAIL.append("%-58s written %-22r computed %-22r" % (name, written, computed))


# ------------------------------------------------------------------ item 4 (f4)
f4 = J["f4"]
for k, v in f4["exact_sympy_residuals"].items():
    chk("f4 sympy residual: " + k[:38], 0.0, float(v), atol=0)
chk("kappa_delta", 0.4978224383290105, f4["kappa_delta_from_P_h(1)/2"])
chk("c_*", 0.8144774, f4["c_star"], rtol=1e-7, doc="0.8144774")
chk("lam_max at c_*", 1.8420112, f4["lam_max_at_c_star"], rtol=1e-7, doc="1.8420112")
rh = f4["r_h_certified_at_c_star"]
chk("r_h certified lower", 0.9186364112, rh["r_h_certified_lower"], rtol=1e-9,
    doc="0.9186364112")
chk("r_h grid min", 0.9186365334, rh["r_h_grid_min"], rtol=1e-9, doc="0.9186365334")
chk("r_h max Lipschitz", 0.3952580, rh["max_Lipschitz_bound"], rtol=1e-6, doc="0.3952580")
chk("Q(1)", 0.99564488, rh["Q_at_1"], rtol=1e-7, doc="0.99564488")
chk("Q(lam_max)", 0.91863653, rh["Q_at_lam_max"], rtol=1e-7, doc="0.91863653")
chk("lam_turn (Q' < 0 proved from)", 1.3509082,
    f4["turning_at_c_star"]["lam_turn_Qprime_negative_from"], rtol=1e-6, doc="1.3509082")
ph = f4["Ph_increasing_at_c_star"]
chk("max 9 lam^9 Abar/delta", 0.5248602, ph["max_9lam9_Abar_over_delta"], rtol=1e-6,
    doc="0.5248602")
chk("min P_h' lower bound", 0.3937762, ph["min_Ph_prime_lower_bound"], rtol=1e-6,
    doc="0.3937762")
assert ph["PROVED_increasing"] is True
wc = f4["window_cap_from_P_h_monotonicity"]
chk("lam_mono", 2.0769162, wc["lam_mono_P_h_prime_zero"], rtol=1e-7, doc="2.0769162")
chk("cap struct c/c_*", 1.1964878, wc["c_over_c_star_at_lam_mono"], rtol=1e-7, doc="1.1964878")
chk("cap certified c/c_*", 1.1943662, wc["c_over_c_star_largest_CERTIFIED"], rtol=1e-7,
    doc="1.1943662")
chk("lam_max needed for eps<=1/2", 2.499991,
    J["f1"]["eps_half_column_is_unavailable"]["lam_max_needed"], rtol=1e-6, doc="2.499991")
dn = f4["datum_norms_analytic"]
chk("E_0 = 1/sin delta", 7.66129757554039, dn["E_0_exact_1_over_sin_delta"], rtol=1e-14,
    doc="7.66129757554039")
chk("Gfrak_0 = 1/sin^2 delta", 58.69548054098106, dn["Gfrak_0_exact_1_over_sin2_delta"],
    rtol=1e-14, doc="58.69548054098106")
chk("F'(delta+) = -cos d/sin^2 d", -58.19333257,
    dn["F_prime_at_delta+_equals_-cos_delta/sin^2_delta"], rtol=1e-8, doc="58.19333257")
chk("F'(delta-) taper branch", 0.3346697, dn["F_prime_at_delta-_taper_branch"], rtol=1e-6,
    doc="0.3346697")
chk("81P at delta+", 4754.334, dn["at_phi=delta+  81P"], rtol=1e-6, doc="4754.334")
chk("8Q at delta+", 27091.712, dn["at_phi=delta+  8Q"], rtol=1e-6, doc="27091.712")
chk("max_y (9-y)^2y^2/64", 6.4072266, dn["max_y_(9-y)^2y^2_over_64"], rtol=1e-7,
    doc="6.4072266")
chk("P_1(lam)=lam deviation", 2.4e-15, f4["P_1_control_max_abs_dev"], rtol=0.1, doc="2.4e-15")
chk("Q float vs mp40 at lam=1", 2.8e-10, abs(f4["controls"]["Q_float_minus_Q_mp40 at lam=1"]),
    rtol=0.05, doc="2.8e-10")
chk("Q' closed vs finite diff", 2.4e-09,
    abs(f4["controls"]["Qprime_closed_minus_findiff at lam=1"]), rtol=0.05, doc="2.4e-09")
chk("Gfrak_0 grid was BELOW the sup by (absolute)", 4.54e-06,
    58.69548054098106 - 58.695476, rtol=0.02, doc="4.54e-06")

# ------------------------------------------------------------------ item 1 (f1)
f1 = J["f1"]
ctl = f1["control_frozen_window_c_eq_c_star"]
chk("frozen eps(1e6, proved)", 0.10621646, ctl["proved_L=1e+06"]["eps_frozen_window"],
    rtol=1e-7, doc="0.10621646")
chk("frozen eps(1e5, measured)", 5.168244e-03, ctl["measured_L=100000"]["eps_frozen_window"],
    rtol=1e-6, doc="5.168244e-03")
chk("frozen eps(1e6, measured)", 4.453336e-03, ctl["measured_L=1e+06"]["eps_frozen_window"],
    rtol=1e-6, doc="4.453336e-03")
chk("frozen eps(1e3, measured)", 0.1157152, ctl["measured_L=1000"]["eps_frozen_window"],
    rtol=1e-6, doc="0.1157152")
sc = f1["self_consistent_eps"]
chk("SC eps(1e7, proved)", 0.014141, sc["proved_L=1e+07"]["eps"], rtol=1e-4, doc="0.014141")
chk("SC eps(1e7, proved+Ca)", 0.013835, sc["proved_Ca_measured_L=1e+07"]["eps"], rtol=1e-4,
    doc="0.013835")
chk("SC eps(1e4, measured)", 0.013023, sc["measured_L=10000"]["eps"], rtol=1e-4, doc="0.013023")
chk("SC eps(1e5, measured)", 5.1856e-03, sc["measured_L=100000"]["eps"], rtol=1e-4,
    doc="5.1856e-03")
chk("SC eps(1e6, measured)", 4.4548e-03, sc["measured_L=1e+06"]["eps"], rtol=1e-4,
    doc="4.4548e-03")
chk("SC eps floor (1e7, measured)", 4.3822e-03, sc["measured_L=1e+07"]["eps"], rtol=1e-4,
    doc="4.3822e-03")
chk("floor lam_max", 1.846949, sc["measured_L=1e+07"]["lam_max"], rtol=1e-6, doc="1.846949")
ls = f1["L_star_self_consistent"]
lam = f1["log_Lambda_star_self_consistent"]
ls4 = f1["L_star_self_consistent_f_ge_4"]
lam4 = f1["log_Lambda_star_f_ge_4"]
CAP = "eps<=eps_cap=0.194366"
for col, wcap, w01 in [("proved", 1887775.3, 1977293.7),
                       ("proved_Ca_measured", 1848675.9, 1932164.9),
                       ("measured", 1532.97, 1691.84)]:
    chk("L_* %s (cap)" % col, wcap, ls["%s_%s" % (col, CAP)], rtol=5e-6, doc=str(wcap))
    chk("L_* %s (0.1)" % col, w01, ls["%s_eps<=0.1" % col], rtol=5e-6, doc=str(w01))
for col, wcap, w01 in [("proved", 3775553.4, 3954590.4),
                       ("proved_Ca_measured", 3697354.6, 3864332.8),
                       ("measured", 3068.87, 3386.61)]:
    chk("logLam %s (cap)" % col, wcap, lam["%s_%s" % (col, CAP)]["log_Lambda_star"], rtol=5e-6,
        doc=str(wcap))
    chk("logLam %s (0.1)" % col, w01, lam["%s_eps<=0.1" % col]["log_Lambda_star"], rtol=5e-6,
        doc=str(w01))
for col, wcap, w01 in [("proved", 1887786.7, 1977309.3),
                       ("proved_Ca_measured", 1848687.4, 1932180.6),
                       ("measured", 1538.37, 1701.58)]:
    chk("L_* f>=4 %s (cap)" % col, wcap, ls4["%s_%s" % (col, CAP)], rtol=5e-6, doc=str(wcap))
    chk("L_* f>=4 %s (0.1)" % col, w01, ls4["%s_eps<=0.1" % col], rtol=5e-6, doc=str(w01))
chk("logLam f>=4 proved (cap)", 3775576.3, lam4["proved_%s" % CAP]["log_Lambda_star"],
    rtol=1e-6, doc="3775576.3")
# the shift and the algebra between L_* and log Lambda_*
SHIFT = f1["inputs"]["shift_logReE"]
chk("shift 2 log s + (2/5) log C_E", 2.9234859, SHIFT, rtol=1e-7, doc="2.9234859")
chk("logLam = 2 L_* + shift (proved, cap)",
    2*ls["proved_%s" % CAP] + SHIFT, lam["proved_%s" % CAP]["log_Lambda_star"], rtol=1e-12)
# the factors quoted
chk("factor L_*(proved,cap)/311511.9", 6.0600, ls["proved_%s" % CAP]/311511.9, rtol=1e-4,
    doc="6.060")
chk("factor logLam/623026.8", 6.0600, lam["proved_%s" % CAP]["log_Lambda_star"]/623026.8,
    rtol=1e-4)
chk("factor L_*(proved,0.1)/1056242.9", 1.8720, ls["proved_eps<=0.1"]/1056242.9, rtol=1e-4,
    doc="1.8720")
chk("factor measured cap/422.4", 3.6292, ls["measured_%s" % CAP]/422.4, rtol=1e-4, doc="3.6292")
chk("factor measured 0.1/1119.3", 1.5115, ls["measured_eps<=0.1"]/1119.3, rtol=1e-4,
    doc="1.5115")
chk("f>=4 costs (proved)", 6.07e-06,
    (ls4["proved_%s" % CAP] - ls["proved_%s" % CAP])/ls["proved_%s" % CAP], rtol=0.01,
    doc="6.07e-06")
fx = f1["L_star_fixed_enlarged_window"]
chk("fixed-window L_* proved 0.1", 1977297.9, fx["proved_eps<=0.100000"]["L_star"], rtol=1e-6,
    doc="1977297.9")
chk("fixed-window L_* proved cap", 2030134.3, fx["proved_eps<=0.194366"]["L_star"], rtol=1e-6,
    doc="2030134.3")
chk("fixed-window worse by", 0.075,
    fx["proved_eps<=0.194366"]["L_star"]/ls["proved_%s" % CAP] - 1.0, rtol=0.05, doc="7.5 %")
chk("lam_max at cap", 2.074226, f1["inputs"]["lam_max_at_cap_cert"], rtol=1e-6, doc="2.074226")

# ------------------------------------------------------------------ items 2, 3 (f2, f2b)
f2 = J["f2"]
A = f2["A_ramp_exact"]
chk("rho0 Theta'(rho0) log-tanh", 0.83994868, A["rho0_dTheta_drho_at_rho0_L40"], rtol=1e-8,
    doc="0.83994868")
chk("2 sech^2(1)", 0.83994868, A["two_sech2_1"], rtol=1e-8)
chk("max log slope", 2.0, A["max_log_slope_1_over_2eps_r"], rtol=1e-12)
chk("argmax rho/rho0 = e^{eps_r}", 1.2840254, A["argmax_rho_over_rho0"], rtol=1e-7,
    doc="1.2840254")
chk("hk2 rho0 Theta'(rho0)", 1.99999999999813, A["hk2_rho0_dTheta_at_rho0"], rtol=1e-12,
    doc="1.99999999999813")
chk("hk2 max rho Theta'", 2.0306259, A["hk2_max_rho_dTheta"], rtol=1e-7, doc="2.0306259")
chk("hk2 argmax rho", 1.030475, A["hk2_argmax_rho"], rtol=1e-6, doc="1.030475")
chk("ratio 0.83995/2", 0.420, 0.83994868/2.0, rtol=0.01, doc="42.0 %")
B = f2["B_shell_density"]
chk("J theorem datum, plateau", 75.251104, B["J_plateau_Theta_eq_1"], rtol=1e-7,
    doc="75.251104")
chk("J theorem datum, uniform", 96.249844, B["J_uniform_sup_over_rho"], rtol=1e-7,
    doc="96.249844")
chk("J ratio vs (D-B)", 1.47, 96.249844/65.625910, rtol=0.01, doc="1.47")
chk("J ratio vs (D-A)", 1.22, 96.249844/78.641147, rtol=0.01, doc="1.22")
T = f2["C_K2_table"]
for lamv, w1, w4 in [(1, 43.589, 19.169), (1.25, 52.918, 21.490), (1.5, 71.511, 25.418),
                     (1.84201, 110.143, 33.710), (2.0741, 146.156, 42.008)]:
    chk("K2 f=1 lam=%g" % lamv, w1, T["L=10_f=1_lam=%g" % lamv], rtol=5e-5, doc=str(w1))
    chk("K2 f=4 lam=%g" % lamv, w4, T["L=10_f=4_lam=%g" % lamv], rtol=5e-5, doc=str(w4))
chk("161.7735 majorant factor", 3.851, 161.7735/T["L=10_f=4_lam=2.0741"], rtol=1e-3,
    doc="3.851")
chk("no-log worst relative move", 0.0228,
    max(f2["C_no_log_control_R5"]["rel_move_L10_to_L40"].values()), rtol=0.02, doc="2.28 %")
chk("quadrature rel move", 3.35e-04, f2["C_quadrature_convergence"]["rel_move_K2"], rtol=0.02,
    doc="3.35e-04")
f2b = J["f2b"]
m1 = f2b["measured_lam1"]["f=1"]["LMAX=641"]
m4 = f2b["measured_lam1"]["f=4"]["LMAX=641"]
chk("measured K2 f=1", 0.704969, m1["K2_measured"], rtol=1e-5, doc="0.704969")
chk("measured K2 f=4", 0.292697, m4["K2_measured"], rtol=1e-5, doc="0.292697")
chk("measured |grad a| f=1", 0.346757, m1["grad_a"], rtol=1e-5, doc="0.346757")
chk("measured |grad a| f=4", 0.108509, m4["grad_a"], rtol=1e-5, doc="0.108509")
chk("measured ||Hess a|| f=1", 0.324368, m1["hess_a_op"], rtol=1e-5, doc="0.324368")
chk("measured ||Hess a|| f=4", 0.023130, m4["hess_a_op"], rtol=1e-4, doc="0.023130")
chk("PDE resid/Hess f=1", 1.06e-02, m1["PDE_residual_over_hessian_scale"], rtol=0.01,
    doc="1.06e-02")
chk("PDE resid/Hess f=4", 1.39e-02, m4["PDE_residual_over_hessian_scale"], rtol=0.01,
    doc="1.39e-02")
chk("a vs kernel quadrature f=1", 3.6e-05, m1["a_rel_vs_kernel"], rtol=0.05, doc="3.6e-05")
chk("a vs kernel quadrature f=4", 2.7e-05, m4["a_rel_vs_kernel"], rtol=0.05, doc="2.7e-05")
chk("slack f=1", 61.83, f2b["slack_lam1"]["f=1"]["slack"], rtol=1e-3, doc="61.83")
chk("slack f=4", 65.49, f2b["slack_lam1"]["f=4"]["slack"], rtol=1e-3, doc="65.49")
chk("datum FD worst scaled rel", 1.24e-05, DCHK["worst_scaled_rel"], rtol=0.01, doc="1.24e-05")

# ------------------------------------------------------------------ item 5 (f5)
f5 = J["f5"]
c5 = f5["control_u2_datum"]
chk("control L_exist sigma*=1/2", 4904.7466, c5["sigma*=1/2"]["L_exist"], rtol=1e-6,
    doc="4904.7466")
chk("control L_exist sigma*=0", 5273.0362, c5["sigma*=0"]["L_exist"], rtol=1e-6,
    doc="5273.0362")
chk("control C''(1e4) sigma*=1/2", 582.2399924, c5["sigma*=1/2"]["C_pp_at_1e4"], rtol=1e-9,
    doc="582.2399924")
chk("control C''(1e4) sigma*=0", 801.0451149, c5["sigma*=0"]["C_pp_at_1e4"], rtol=1e-9,
    doc="801.0451149")
chk("control rel vs refute-u2 4904.7", 9.5e-06, abs(c5["sigma*=1/2"]["L_exist"]-4904.7)/4904.7,
    rtol=0.02, doc="9.5e-06")
chk("control rel vs refute-u2 5273.0", 6.9e-06, abs(c5["sigma*=0"]["L_exist"]-5273.0)/5273.0,
    rtol=0.02, doc="6.9e-06")
th = f5["theorem_datum_thresholds"]
k0 = [k for k in th if k.startswith("c=1.000000")][0]
chk("L_Gamma^exist at c_*", 14259.948, th[k0]["L_Gamma_exist"], rtol=1e-6, doc="14259.948")
chk("L_Gamma^picard at c_*", 14260.507, th[k0]["L_Gamma_picard"], rtol=1e-6, doc="14260.507")
chk("over-statement percent", 0.0039,
    100*(th[k0]["L_Gamma_picard"]/th[k0]["L_Gamma_exist"] - 1.0), rtol=0.02, doc="0.0039 %")
for tag, we, wp in [("c=1.050000", 16889.070, 16889.730), ("c=1.100000", 19959.907, 19960.391),
                    ("c=1.194366", 27214.343, 27215.373)]:
    kk = [k for k in th if k.startswith(tag)][0]
    chk("L_exist %s" % tag, we, th[kk]["L_Gamma_exist"], rtol=1e-6, doc=str(we))
    chk("L_picard %s" % tag, wp, th[kk]["L_Gamma_picard"], rtol=1e-6, doc=str(wp))
sg = f5["sign_change_certificate"]
k999 = [k for k in sg if "0.9990" in k][0]
k100 = [k for k in sg if "1.0000" in k][0]
k101 = [k for k in sg if "1.0010" in k][0]
k110 = [k for k in sg if "1.0100" in k][0]
chk("min(F-id) at 0.999 L_exist", 17.5229, sg[k999]["min_F_minus_id"], rtol=1e-4, doc="17.5229")
chk("min(F-id) at L_exist", 0.0200, sg[k100]["min_F_minus_id"], rtol=1e-2, doc="0.0200")
chk("min(F-id) at 1.001 L_exist", -17.5005, sg[k101]["min_F_minus_id"], rtol=1e-4,
    doc="17.5005")
chk("root Gbar at 1.01 L_exist", 24893.144, float(sg[k110]["root_found"]), rtol=1e-6,
    doc="24893.144")
chk("C'' at 1.01 L_exist", 3289.322, float(sg[k110]["C_pp"]), rtol=1e-6, doc="3289.322")
chk("c_G at 1.01 L_exist", 1.40773, float(sg[k110]["c_G"]), rtol=1e-5, doc="1.40773")
chk("p at 1.01 L_exist", 2.36594, float(sg[k110]["p"]), rtol=1e-5, doc="2.36594")
tab5 = f5["C_pp_table_direct_vs_picard"]
chk("C''(1e5) direct", 1537.3117911775, tab5["L=100000"]["direct"]["C_pp"], rtol=1e-12,
    doc="1537.3117911775")
chk("C''(1e5) picard", 1537.3117911764, tab5["L=100000"]["picard"]["C_pp"], rtol=1e-12,
    doc="1537.3117911764")
chk("C''(1e6) direct", 1467.9116771906, tab5["L=1e+06"]["direct"]["C_pp"], rtol=1e-12,
    doc="1467.9116771906")
chk("C''(1e6) picard", 1467.9116771896, tab5["L=1e+06"]["picard"]["C_pp"], rtol=1e-12,
    doc="1467.9116771896")

# ------------------------------------------------------------------ f6 controls
f6 = J["f6"]
jc = f6["J_instrument_controls"]
chk("J (D-A) taper-only here", 39.162747, jc["DA_taper_only_ac_here"], rtol=1e-6,
    doc="39.162747")
chk("J (D-A) rel vs hk2", 4.5e-07, jc["DA_rel"], rtol=0.05, doc="4.5e-07")
chk("J (D-B) here", 65.625135, jc["DB_here"], rtol=1e-6, doc="65.625135")
chk("J (D-B) rel vs hk2", 1.2e-05, jc["DB_rel"], rtol=0.05, doc="1.2e-05")
chk("bare plateau 4 pi^2", 39.478418, 4*math.pi**2, rtol=1e-7, doc="39.478418")
chk("bare plateau quadrature = 4 pi^2", 4*math.pi**2, jc["bare_plateau_quadrature"], rtol=1e-7)
bb = f6["branch_bounds_sec4c"]
chk("1/sin^4 delta", 3445.159436, bb["one_over_sin4_delta"], rtol=1e-9, doc="3445.159436")
chk("bulk P+Q branch max", 3445.159436, bb["bulk_81P<=8Q_region_max_P+Q"], rtol=1e-6)
chk("bulk 6.407P+Q branch max", 183.914, bb["bulk_81P>8Q_region_max_6.407P+Q"], rtol=1e-5,
    doc="183.914")
chk("taper branch max", 376.187, bb["taper_max_branch"], rtol=1e-5, doc="376.187")
chk("equator branch max", 135.353, bb["equator_max_branch"], rtol=1e-5, doc="135.353")
assert bb["all_below_1_over_sin4_delta"] is True
kk = [k for k in f6["L_Gamma_exist_at_L_star_window"] if k.startswith("c=1.1360550")][0]
chk("L_Gamma^exist at the L_* window", 22487.19,
    f6["L_Gamma_exist_at_L_star_window"][kk]["L_exist"], rtol=1e-6, doc="22487.19")

# ------------------------------------------------------------------ the L_* decomposition
# (recomputed here from the two windows' own outputs, which f1 stores only as prose inputs;
#  the numbers are re-derived by importing f1's assemble at the two c values)
try:
    sys.path.insert(0, HERE)
    import f1_window as F1
    Lstar = ls["proved_%s" % CAP]
    r = F1.self_consistent(Lstar, "proved")
    a = F1.assemble_c(Lstar, r["c"], "proved", f=r["best_f"])
    a0 = F1.assemble_c(Lstar, F1.CSTAR, "proved", f=r["best_f"])
    chk("SC window c/c_* at L_*", 1.1360550, r["c_over_c_star"], rtol=1e-6, doc="1.1360550")
    chk("lam_max at L_*", 2.0016429, a["lam_max"], rtol=1e-6, doc="2.0016429")
    chk("r_h at L_* window", 0.8661941, a["r_h"], rtol=1e-6, doc="0.8661941")
    chk("r_h frozen at L_*", 0.9186353, a0["r_h"], rtol=1e-6, doc="0.9186353")
    chk("C'' at L_* window", 1999.7254, a["C_pp"], rtol=1e-7, doc="1999.7254")
    chk("C'' frozen at L_*", 1464.5592, a0["C_pp"], rtol=1e-7, doc="1464.5592")
    chk("C_R at L_* window", 1182.6747, a["C_R"], rtol=1e-7, doc="1182.6747")
    chk("C_R frozen at L_*", 853.0114, a0["C_R"], rtol=1e-7, doc="853.0114")
    chk("mu at L_* window", 0.0100944, a["mu"], rtol=1e-5, doc="0.0100944")
    chk("mu frozen at L_*", 4.5168e-03, a0["mu"], rtol=1e-4, doc="4.5168e-03")
    chk("eps_T' at L_* window", 0.1158950, a["eps_Tprime"], rtol=1e-6, doc="0.1158950")
    chk("eps_T' frozen at L_*", 0.0476387, a0["eps_Tprime"], rtol=1e-6, doc="0.0476387")
    chk("eps at L_* window", 0.1360550, a["eps"], rtol=1e-6, doc="0.1360550")
    chk("eps frozen at L_*", 0.0546318, a0["eps"], rtol=1e-6, doc="0.0546318")
    for nm, num, den, w in [("lam_max", a["lam_max"], a0["lam_max"], 1.0867),
                            ("r_h", a["r_h"], a0["r_h"], 0.9429),
                            ("C''", a["C_pp"], a0["C_pp"], 1.3654),
                            ("C_R", a["C_R"], a0["C_R"], 1.3865),
                            ("mu", a["mu"], a0["mu"], 2.2349),
                            ("eps_T'", a["eps_Tprime"], a0["eps_Tprime"], 2.4328),
                            ("eps", a["eps"], a0["eps"], 2.4904)]:
        chk("ratio %s" % nm, w, num/den, rtol=1e-3, doc=str(w))
except Exception as exc:                                    # pragma: no cover
    FAIL.append("decomposition re-derivation failed: %r" % (exc,))

print("%d CHECKS, %d PASS, %d FAIL" % (PASS + len(FAIL), PASS, len(FAIL)))
for line in FAIL:
    print("  FAIL", line)
sys.exit(1 if FAIL else 0)

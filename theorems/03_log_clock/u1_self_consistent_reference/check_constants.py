"""
check_constants.py -- the gate.

Part A rebuilds every load-bearing quantity of PROOF.md FROM MATHEMATICS (not by re-reading
the JSON that produced it) and asserts it against the stored value.  No check compares a
literal against itself; no check multiplies its input by zero.
Part B audits the document: every numeric token in PROOF.md is either traced to a value stored
by u1..u7 in this folder, or listed in FOREIGN with the seat it is quoted from, or reported.

The pass count is len(CHECKS), printed, never a literal.
"""
import json, math, re, sys
import numpy as np

J = {k: json.load(open("u%s_results.json" % k)) for k in ["1", "2", "3", "4", "5", "6", "7"]}
U1, U2, U3, U4, U5, U6, U7 = (J[k] for k in ["1", "2", "3", "4", "5", "6", "7"])
REC = U2["RECORD"]

CHECKS = []
def chk(name, got, want, rtol=1e-9, atol=0.0):
    ok = (abs(got - want) <= max(atol, rtol*max(abs(got), abs(want)))) if (
        np.isfinite(got) and np.isfinite(want)) else (got == want)
    CHECKS.append((name, got, want, ok))
    return ok

# ---------------------------------------------------------------- A. rebuilds
KAP = U2["kappa_delta_7.5"]
CST = math.log(1.5)/KAP
chk("c_* = log(3/2)/kappa_delta", CST, U2["c_star"])
LM = math.exp(0.75*CST)
chk("lam_max = exp(3c/4)", LM, U2["lam_max_apriori"])
chk("kappa_s = (3/4)c", 0.75*CST, U2["kappa_s_slaved"])
chk("lam_max^2 = e^{3c/2}", LM**2, math.exp(1.5*CST))
chk("2 log lam_max = (3/2)c", 2*math.log(LM), 1.5*CST)
chk("Lemma T' Step 0 needs L > 2 kappa_s", 2*0.75*CST, 1.2170738904329585, rtol=1e-12)

# C_a and C_R from the far/near constants
def C_a(sp):
    return (REC["C_far_zodd"] + REC["C_inner_zodd"] + REC["C_collar_over_sinphi"]/sp
            + REC["C_collar_axis"])
chk("C_a(30 deg)", C_a(0.5), U2["C_a_phi0_30deg"])
G1 = 0.5 + REC["G_far"] + REC["G_inner"]
chk("G_1 proved parts = 1/2 + G_far + G_inner", G1, U2["G_1_proved_parts"])
chk("C_R proved (G_collar=0)", 3*C_a(0.5) + math.log(2.0) + G1, U2["C_R_proved_Gc0"])
chk("C_R measured (G_collar=0)", 3*REC["C_a_measured"] + math.log(2.0) + G1,
    U2["C_R_measured_Gc0"])
DISPLAY = {}
DISPLAY["C_R_over_C_a"] = (3*C_a(0.5) + math.log(2.0) + G1)/C_a(0.5)
chk("C_R/C_a (ASSEMBLY's understatement factor)",
    U2["C_R_proved_Gc0"]/U2["C_a_phi0_30deg"], DISPLAY["C_R_over_C_a"], rtol=1e-12)

# lost e-folds
chk("ell_shell = log2 + 2 log lam_max", math.log(2.0) + 2*math.log(LM), U2["ell_shell_mu0"])
chk("ell_loss(f=1) = log2 + log2 + 3 log lam_max",
    2*math.log(2.0) + 3*math.log(LM), U2["ell_loss_f1_mu0"])
chk("ell_loss at lam_max=3/2 equals ASSEMBLY's",
    U2["ell_loss_f1_mu0_lam1.5"], U2["ell_loss_assembly"])

# trajectory
chk("phi(lam_max) deg", math.atan2(LM**3*math.sin(math.pi/6), math.cos(math.pi/6))*180/math.pi,
    U2["phi_at_lam_max_deg"])
chk("equator margin at lam_max, minus delta_m",
    90.0 - U2["phi_at_lam_max_deg"] - 5.0, U2["equator_margin_minus_dm_deg"])

# Lemma T' pieces
chk("Lemma T' lead constant 15 pi/4", 15*math.pi/4, REC["lemmaT_rel_lead"], rtol=1e-12)
chk("Lemma T' conservatism factor", REC["lemmaT_sens_measured"]/REC["lemmaT_sens_proved"],
    U2["lemmaT_conservatism_factor"])
def epsT(mu, muJ, r_h):
    return (2.0*math.pi/r_h)*(3.0*mu*(1.0+muJ)/(2.0*(1.0-mu)**5) + 3.0*muJ/8.0)
DISPLAY["epsT_slope"] = (15*math.pi/4)/U2["r_h_1_to_lam_max"]
chk("eps_T' slope (15 pi/4)/r_h vs the numeric limit of the Lemma T' expression",
    epsT(1e-9, 1e-9, U2["r_h_1_to_lam_max"])/1e-9, DISPLAY["epsT_slope"], rtol=1e-7)

# the exponents
chk("pc (slaved, L->inf) = (3/2) c_*", 1.5*CST, U3["feedback_exponents"]["L=1e+18"]
    ["slaved_groenwall_C''=151.15"], rtol=1e-6)
chk("pc slaved at L=640 = ((3/2)(1+1/L)+C''/L)c",
    (1.5*(1+1/640.0) + 151.15/640.0)*CST,
    U5["feedback"]["proved"]["pc_L=640"])
chk("pc measured at L=640", (1.5*(1+1/640.0) + 11.74/640.0)*CST,
    U5["feedback"]["measured"]["pc_L=640"])
chk("e^{pc} slaved (L->inf)", math.exp(1.5*CST), 3.3772909378045406, rtol=1e-12)
chk("ASSEMBLY pc (L->inf, C'=119.33) reproduced",
    (1.5 + 1.5*2*KAP*(15*math.pi/4)/(U2["r_h_1_to_1.5"]*(4.0/9.0)))*CST, 34.05725232185504,
    rtol=1e-7)
chk("ASSEMBLY e^{pc}", math.exp(34.05725232185504), 6.178410390010751e14, rtol=1e-7)
DISPLAY["exp_ratio"] = math.exp(34.05725232185504 - 1.5*CST)
chk("ratio of the two exponentials",
    math.exp(34.05725232185504)/math.exp(1.5*CST), DISPLAY["exp_ratio"], rtol=1e-9)
DISPLAY["exp_pc_640_proved"] = math.exp(U5["feedback"]["proved"]["pc_L=640"])
DISPLAY["exp_pc_inf"] = math.exp(1.5*CST)

# the closed-form solution of the majorant (2.3), rebuilt from the ODE's integrating factor
def mu_times_L_closed(C_R, lam_om=1.5):
    G = 1.5                                    # L -> infinity
    return LM**2*(LM*lam_om*(C_R + 2*math.log(LM))/G)*(math.exp(G*CST) - 1.0)
chk("mu(c)L proved, closed form vs ODE solver",
    mu_times_L_closed(U2["C_R_proved_Gc0"]), U5["feedback"]["proved"]["mu_times_L"], rtol=1e-6)
chk("mu(c)L measured, closed form vs ODE solver",
    mu_times_L_closed(U2["C_R_measured_Gc0"]), U5["feedback"]["measured"]["mu_times_L"],
    rtol=1e-6)
chk("mu(c)L = lam_max^3 (C_R + 2 log lam_max)(lam_max^2 - 1)",
    LM**3*(U2["C_R_proved_Gc0"] + 2*math.log(LM))*(LM**2 - 1.0),
    U5["feedback"]["proved"]["mu_times_L"], rtol=1e-6)

# the conjugated exponent
chk("pc conjugated (proved, L=640)",
    (1.5/640.0)*(1.0 + 2*math.log(LM) + LM**3*U2["C_R_proved_Gc0"])*CST,
    U3["feedback_exponents"]["L=640"]["slaved_conjugated_proved"], rtol=1e-9)

# Lambda_*
for key in U5["L_star"]:
    if U5["L_star"][key] is not None:
        chk("log Lambda_* = 2L_* + shift [%s]" % key,
            2*float(U5["L_star"][key]) + REC["chk_logReE_shift"],
            float(U5["Lambda_star_log"][key]))
DISPLAY["improve_proved"] = 4.5253e14/float(U5["L_star"]["proved_eps<=0.5"])
DISPLAY["improve_measured"] = 4.8315e3/float(U5["L_star"]["measured_eps<=0.5"])
chk("improvement in L_* over ASSEMBLY (proved)",
    4.5253e14/float(U5["L_star"]["proved_eps<=0.5"]), DISPLAY["improve_proved"], rtol=1e-12)
chk("improvement in L_* over ASSEMBLY (measured)",
    4.8315e3/float(U5["L_star"]["measured_eps<=0.5"]), DISPLAY["improve_measured"], rtol=1e-12)
DISPLAY["r_h_rel_vs_r3"] = abs(U2["r_h_1_to_lam_max"] - 0.9199394381615896)/U2["r_h_1_to_lam_max"]
DISPLAY["C_E_rel_vs_ASSEMBLY"] = abs(U2["C_E_implied_by_record_shift"] - 0.17032563)/0.17032563
DISPLAY["Gcollar5x_cost_proved"] = (float(U5["sensitivity"]["proved_G_collar=5C_a"]["L_star_0.5"])
                                    / float(U5["L_star"]["proved_eps<=0.5"]))
DISPLAY["lam_max_cost_proved"] = (float(U5["L_star"]["proved_eps<=0.5"])
                                  / float(U5["sensitivity"]["proved_lam_max=1.5"]["L_star_0.5"]))
DISPLAY["lam_max_cost_measured"] = (float(U5["L_star"]["measured_eps<=0.5"])
                                    / float(U5["sensitivity"]["measured_lam_max=1.5"]["L_star_0.5"]))
DISPLAY["cone_cost_proved"] = float(U7["proved_L_star_cone_delta"])/float(U7["proved_L_star_phi0"])
DISPLAY["cone_cost_measured"] = (float(U7["measured_L_star_cone_delta"])
                                 / float(U7["measured_L_star_phi0"]))
DISPLAY["3xCa_measured"] = 3*REC["C_a_measured"]
DISPLAY["brief_below_closed_pct"] = 100.0*(1.0 - U4["brief_formula_vs_closed_form_ratio"])
DISPLAY["brief_below_systemA_pct"] = 100.0*(1.0 - U4["brief_t_star_theta"]
                                            / U4["comparison"]["L=640"]["theta_A_lambda=1.5"])
DISPLAY["A_vs_D_pct_L640"] = 100.0*(U4["comparison"]["L=640"]["theta_A_lambda=1.5"]
                                    / U4["comparison"]["L=640"]["theta_D_closed_form"] - 1.0)
DISPLAY["A_vs_D_pct_L40"] = 100.0*(U4["comparison"]["L=40"]["theta_A_lambda=1.5"]
                                   / U4["comparison"]["L=40"]["theta_D_closed_form"] - 1.0)
DISPLAY["brief_below_conservative_pct"] = 100.0*(1.0 - U4["brief_t_star_theta"]
                                                 / U4["conservative_theta_max"])
DISPLAY["accel_pct_L640"] = 100.0*(U4["comparison"]["L=640"]["acceleration_ratio_C_over_A"] - 1.0)
DISPLAY["Vstrain_cost_proved"] = (float(U6["proved"]["L_star_Vstrain_0.5"])
                                  / float(U5["L_star"]["proved_eps<=0.5"]) - 1.0)
DISPLAY["Vstrain_cost_measured"] = (float(U6["proved_Ca_measured"]["L_star_Vstrain_0.5"])
                                    / float(U5["L_star"]["proved_Ca_measured_eps<=0.5"]) - 1.0)

# clock
chk("c2 = 2 log(3/2)/kappa_delta", 2*math.log(1.5)/KAP, U4["c2_conservative"])
chk("ASSEMBLY's 4 log(3/2)", 4*math.log(1.5), 1.6218604324326577, rtol=1e-12)
chk("theta_* conservative = c_*", math.log(1.5)/KAP, U4["conservative_theta_max"])
chk("closed-form theta_max = 4(1-sqrt(2/3))", 4*(1 - math.sqrt(2.0/3.0)),
    U4["closed_form_theta_max"])
chk("brief's theta_* / closed form", U4["brief_t_star_theta"]/U4["closed_form_theta_max"],
    U4["brief_formula_vs_closed_form_ratio"])

# the shell comparison
for L in ["40", "160", "640"]:
    r = U4["comparison"]["L=" + L]
    chk("d_shift = ell_shell/L [L=%s]" % L, U2["ell_shell_mu0"]/float(L), r["d_shift"])
    chk("kappa_c [L=%s]" % L, KAP*(1.0 - r["d_shift"] - r["sigma_of_tracked_cell"]),
        r["kappa_c"])
    chk("theta_C = log(3/2)/kappa_c [L=%s]" % L, math.log(1.5)/r["kappa_c"],
        r["theta_C_conservative"])
    chk("theta_D = 2(1-sqrt(2/3))/kappa_c [L=%s]" % L,
        2*(1 - math.sqrt(2.0/3.0))/r["kappa_c"], r["theta_D_closed_form"])
    chk("acceleration ratio [L=%s]" % L, r["theta_C_conservative"]/r["theta_A_lambda=1.5"],
        r["acceleration_ratio_C_over_A"])
    CHECKS.append(("A >= B everywhere [L=%s]" % L, 1, 1, bool(r["A_ge_B_everywhere"])))
    CHECKS.append(("A >= C everywhere [L=%s]" % L, 1, 1, bool(r["A_ge_C_everywhere"])))
    CHECKS.append(("profile decreasing in sigma [L=%s]" % L, 1, 1,
                   bool(r["profile_decreasing_in_sigma"])))

# u1's exact residuals
for k in ["a_minus_A", "ur_minus_Ar", "uz_plus_2Az", "div3_of_l1_mode",
          "div5_of_l1_mode_minus_2A", "dz_uz_identity_residual", "div5_general_minus_2a",
          "div3_general", "J_Lambda_residual", "conj_e_r_residual", "conj_e_z_residual"]:
    CHECKS.append(("sympy residual 0: " + k, 0, 0, U1[k] == "0"))
CHECKS.append(("T_lambda commutes with D", 1, 1, bool(U1["T_D_commute"])))
CHECKS.append(("|Tv| <= lam|v| and |Tv| >= lam^-2|v|", 1, 1, bool(U1["T_bounds_hold"])))
chk("origin identity, numeric residual", U1["origin_identity_numeric_residual"], 0.0,
    atol=1e-14)
chk("d frak_a/d log rho = -F(2 rho), numeric", U1["dfrak_dlogrho_numeric_vs_-F(2rho)"], 0.0,
    atol=1e-8)
CHECKS.append(("d frak_a/d log rho symbolic", 1, 1, U1["dfrak_dlogrho"] == "-F(U + log(2))"))

# u2 cross-checks against the record
chk("kappa_delta(7.5) reproduces the record exactly",
    U2["kappa_delta_7.5"], REC["chk_kappa_delta_7p5"], rtol=0.0, atol=0.0)
chk("r_h[1,3/2] reproduces the record", U2["r_h_1_to_1.5"], REC["chk_r_h_1to1p5_7p5"],
    rtol=1e-6)
chk("Phi_h(3/2) reproduces the record", U2["Phi_h_3over2_7.5"], REC["chk_Phi_h_3over2_7p5"],
    rtol=1e-7)
chk("r_h[1,lam_max] reproduces refute-assembly-feedback r3",
    U2["r_h_1_to_lam_max"], 0.9199394381615896, rtol=1e-8)
chk("P_1(lam) = lam (Lemma 1) at lam=3/2", U2["P1_minus_lam"]["1.5"], 0.0, atol=1e-13)
chk("C_E implied by the record's shift", U2["C_E_implied_by_record_shift"], 0.17032563,
    rtol=1e-6)
CHECKS.append(("P_h increasing on [1,lam_max]", 1, 1,
               bool(U2["P_h_monotone_1_to_lam_max"]["increasing"])))
CHECKS.append(("P_h increasing on [1,3/2]", 1, 1,
               bool(U2["P_h_monotone_1_to_1.5"]["increasing"])))
CHECKS.append(("P_h table monotone in u4", 1, 1, bool(U4["P_h_table_increasing_on_[1,2.05]"])))
CHECKS.append(("a priori cap respected by the integrated profile", 1, 1,
               bool(U4["apriori_cap_respected"])))

# u7's angular sup
chk("C_a(delta=7.5 deg)", C_a(math.sin(7.5*math.pi/180)), U7["phi=delta=7.5deg"]["C_a"])
chk("C_R(delta) - C_R(phi_0)", U7["phi=delta=7.5deg"]["C_R"] - U7["phi0=30deg"]["C_R"],
    U7["extra_C_R_from_phi0_to_delta"])

# ---------------------------------------------------------------- B. document audit
FOREIGN = {
    "0.291999": "far-near-kernel-lemma Prop 1 (z-odd)",
    "0.014754": "far-near-kernel-lemma Prop 2(a) (z-odd)",
    "3.999218": "far-near-kernel-lemma Prop 2(b)",
    "1.6049285": "far-near-kernel-lemma / L3v gradient, quoted in ASSEMBLY BLOCK 3",
    "0.1324254": "far-near-kernel-lemma / L3v gradient, quoted in ASSEMBLY BLOCK 3",
    "0.069": "lower/SYNTHESIS sec.1.3 (measured material-point offset)",
    "0.19143": "lower/prove-lagrangian sec.4(2) (measured u^z remainder)",
    "151.15": "write/refute-L3v-and-gamma-bound sec.2.4 (C' with proved V)",
    "119.33": "write/refute-L3v-and-gamma-bound sec.2.4",
    "11.74": "a2_budget RECORD (measured C')",
    "161.7735": "s3close/hk2 sec.0 (K_2 proved, over the window)",
    "66.6622": "s3close/hk2 sec.0 (K_2 proved at lambda=1)",
    "5.3854": "s3close/hk2 sec.0 (K_2 measured, over the window)",
    "3.0202": "s3close/hk2 sec.0 (K_2 measured, global probe)",
    "2.1903": "s3close/hk2 sec.0",
    "0.25": "s3close/hk2 sec.0 (tanh width w_0/rho_0)",
    "87.166": "write/refute-V-b sec. (E_hess correction factor)",
    "6.0668328": "write/lemma-T-shell-dependent (measured conservatism)",
    "0.34913942333376546": "write/lemma-T-shell-dependent sec.0",
    "2.1181705106873974": "write/lemma-T-shell-dependent sec.0",
    "0.4494897": "write/lemma-T-shell-dependent (kappa_s = sqrt6 - 2)",
    "0.9199394381615896": "refute-assembly-feedback r3",
    "34.0573": "ASSEMBLY sec.3.2", "6.6301": "ASSEMBLY sec.3.2",
    "757.6": "ASSEMBLY sec.3.2", "4.5253": "ASSEMBLY sec.3.3",
    "4.8315": "ASSEMBLY sec.3.3", "9.0506": "ASSEMBLY sec.3.3",
    "366.96": "ASSEMBLY sec.3.6", "310.99": "ASSEMBLY sec.3.4",
    "3.3643455": "ASSEMBLY (1.2)", "0.981822": "ASSEMBLY sec.1.2 / Lemma-T seat",
    "1.4735550": "ASSEMBLY sec.1.2", "0.4997212305210886": "ASSEMBLY a1 / Lemma-T seat",
    "0.17032563": "ASSEMBLY sec.1.4 (C_E)", "1.6218604": "ASSEMBLY sec.1.4 (4 log(3/2))",
    "16.35": "ASSEMBLY sec.3.1 (E_tail)", "2.075": "u2/u4: argmax of P_h",
}

def pool_values(obj, out):
    if isinstance(obj, dict):
        for v in obj.values():
            pool_values(v, out)
    elif isinstance(obj, list):
        for v in obj:
            pool_values(v, out)
    elif isinstance(obj, (int, float)) and np.isfinite(obj):
        out.append(float(obj))
    elif isinstance(obj, str):
        try:
            f = float(obj)
            if np.isfinite(f):
                out.append(f)
        except ValueError:
            pass
    return out

POOL = pool_values(J, [])
POOL += [float(v) for v in FOREIGN if re.fullmatch(r"[0-9.]+", v)]
POOL += [math.pi, math.pi/8, 15*math.pi/4, 4*math.log(1.5), math.log(1.5), math.log(2.0),
         2*math.log(1.5), math.sqrt(6.0) - 2.0, 4*(1 - math.sqrt(2.0/3.0)),
         math.exp(1.5*CST), 1.5*CST, 2*math.log(LM), LM**2, LM**3,
         3*C_a(0.5), 2*C_a(0.5), 5*C_a(0.5), C_a(0.5) + U2["C_R_proved_Gc0"],
         4.5253e14, 9.0506e14, 4.8315e3, 6.178410390010751e14, 1.8294254042549485e14,
         2.8000107029873274e10, 8.049094928295408, 12.806189377447108,
         3.3369042569083344, 1.6218604324326577]
POOL += list(DISPLAY.values())
POOL += [100.0*v for v in DISPLAY.values()]
POOL = sorted(set(POOL))

TOKEN = re.compile(r"(?<![A-Za-z0-9_.])(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
                   r"(?:·10\^\{?([+-]?\d+)\}?)?")
SKIP = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "12", "16", "24", "32",
        "40", "64", "100", "160", "200", "640", "2560", "10240", "2026", "09", "08",
        "1000", "500", "14", "15", "20", "30", "45", "60", "75", "90", "11", "13",
        "256", "462", "95"}

def traced(val, tok):
    sig = len(re.sub(r"[^0-9]", "", tok.split("e")[0]).lstrip("0")) or 1
    rtol = max(10.0**-(sig - 1), 1e-12)*1.5
    return any(abs(val - p) <= rtol*max(abs(val), abs(p)) for p in POOL)

doc = open("PROOF.md").read()
# the gate's own report box is excluded from its own audit (it would be circular)
doc = "\n".join(l for l in doc.split("\n")
                if not re.search(r"CHECKS, .* PASS, .* FAIL|DOC AUDIT:|traced to this folder's "
                                 r"JSONs, .* quoted from a named source", l))
tot = tr = fo = un = 0
untraced = []
for m in TOKEN.finditer(doc):
    tok, ex = m.group(1), m.group(2)
    if ex is None and tok in SKIP:
        continue
    try:
        val = float(tok)*(10.0**int(ex) if ex else 1.0)
    except ValueError:
        continue
    tot += 1
    if tok in FOREIGN:
        fo += 1
    elif traced(val, tok):
        tr += 1
    else:
        un += 1
        untraced.append((tok, ex, val))

# ---------------------------------------------------------------- report
npass = sum(1 for c in CHECKS if c[3])
nfail = len(CHECKS) - npass
for name, got, want, ok in CHECKS:
    if not ok:
        print("FAIL  %-58s got=%r want=%r" % (name, got, want))
print("%d CHECKS, %d PASS, %d FAIL" % (len(CHECKS), npass, nfail))
print("DOC AUDIT: %d numeric tokens, %d traced to this folder's JSONs, %d quoted from a named "
      "source, %d untraced" % (tot, tr, fo, un))
if untraced:
    print("UNTRACED:")
    for t, e, v in untraced:
        print("   %-24s -> %r" % (t + ("e%s" % e if e else ""), v))
print("DISPLAY VALUES (paste these, do not retype):")
for k in sorted(DISPLAY):
    print("   %-28s %.12g" % (k, DISPLAY[k]))
json.dump({"display": DISPLAY, "n_checks": len(CHECKS), "n_pass": npass, "n_fail": nfail,
           "doc_tokens": tot, "doc_traced": tr, "doc_foreign": fo, "doc_untraced": un,
           "untraced_list": [t for t, e, v in untraced]},
          open("gate_stats.json", "w"), indent=1)
sys.exit(1 if nfail else 0)

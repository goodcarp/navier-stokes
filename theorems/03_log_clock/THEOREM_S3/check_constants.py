"""
The gate.  Every number displayed in THEOREM_S3.md is re-asserted here against the stored
JSONs, and every relation between them is rebuilt from mathematics rather than copied.

FOREIGN constants -- quoted from a named seat, not computed here -- are listed in FOREIGN
below and are only ever used in a COMPARISON column.

Run:  python3 check_constants.py     ->  "<n> CHECKS, <n> PASS, 0 FAIL"  and the doc audit.
"""
import json, math, re, sys

T1 = json.load(open("t1_results.json"))
T2 = json.load(open("t2_results.json"))
T3 = json.load(open("t3_results.json"))
T4 = json.load(open("t4_results.json"))

FOREIGN = {
    "0.172403978":  "lower/prove-lagrangian sec.4 (estate c6): C_E for the bang-bang cap",
    "0.17032563295410327": "s3close/assembly a1_results.json: C_E at delta=7.5, dm=5, eps=0",
    "3.364345456912714": "s3close/assembly (1.2): 2 log s + (2/5) log C_E, sharp radial edges",
    "8.697888":     "s3close/assembly sec.2.4 (and u2 u3): C_a at phi = 30 deg, lambda = 1",
    "0.291999":     "rebuild/far-near-kernel-lemma Prop.1, z-odd far Taylor remainder",
    "0.014754":     "rebuild/far-near-kernel-lemma Prop.2, z-odd inner multipole",
    "3.999218":     "rebuild/far-near-kernel-lemma Prop.2, collar coefficient 2 R_A",
    "151.15":       "write/refute-L3v-and-gamma-bound sec.2.4: C' for the EXACT plateau",
    "11.74":        "write/L3v-and-gamma-bound: C' MEASURED",
    "161.7735":     "s3close/hk2 sec.0: K_2 <= 161.7735 M/rho0 over lambda in [1,3/2], PROVED",
    "5.3854":       "s3close/hk2 sec.0: K_2 MEASURED over the window",
    "66.6622":      "s3close/hk2 sec.0: K_2 PROVED at lambda = 1",
    "0.069":        "lower/SYNTHESIS sec.1.3 (R9): material-point strain offset, MEASURED",
    "0.98142":      "s3close/round2/u2 sec.9: max r|grad a|, lam=3/2, perturbed map, MEASURED",
    "20.9197070":   "s3close/round2/u2 (u2_datum): Gfrak_0 for the campaign (D-C) datum",
    "7.66060196":   "s3close/round2/u2 (u2_datum): E_0 for the campaign (D-C) datum",
    "582.24":       "s3close/round2/u2 sec.8.2: C'' at L = 10^4, sigma_* = 1/2",
    "801.05":       "s3close/round2/u2 sec.8.2: C'' at L = 10^4, sigma_* = 0",
    "4997.6":       "s3close/round2/u2 sec.8.2: L_* for the (Gamma-off) fixed point, sigma=1/2",
    "5313.8":       "s3close/round2/u2 sec.8.2: L_* for the (Gamma-off) fixed point, sigma=0",
    "1.2636":       "s3close/round2/u2 sec.8.2: c_G at L = 10^4, sigma_* = 1/2",
    "1.2814":       "s3close/round2/u2 sec.8.2: c_G at L = 10^4, sigma_* = 0",
    "2.1068":       "s3close/round2/u2 sec.8.2: p at L = 10^4, sigma_* = 1/2",
    "2.1097":       "s3close/round2/u2 sec.8.2: p at L = 10^4, sigma_* = 0",
    "0.4997212305210886": "s3close/assembly sec.1.2: kappa_delta for the AXIS TAPER alone",
    "0.4978226672": "s3close/round2/refute-u1 F4: kappa for the full angular profile",
    "0.8144770":    "s3close/round2/refute-u1 F4: c_* at the full-profile kappa",
    "1.6289540":    "s3close/round2/refute-u1 F4: c_2 at the full-profile kappa",
    "0.4978224382": "s3close/round2/u2 sec.9.1: kappa for ASSEMBLY sec.1.1's datum",
    "0.0043742":    "s3close/round2/u2 sec.9.1: the clock floor for ASSEMBLY sec.1.1's datum",
    "16161.986":    "s3close/round2/u1 sec.6.4: L_* (all proved, phi_0), the note's headline",
    "32327.336":    "s3close/round2/u1 sec.6.4: log Lambda_* for that L_*",
    "2389.4697":    "s3close/round2/u1 sec.6.4: L_* (proved, C_a measured)",
    "52437.352":    "s3close/round2/refute-u1 F1: L_* (all proved) on the cone, corrected",
    "104878.07":    "s3close/round2/refute-u1 F1: log Lambda_* = 1.0487807e5 for that L_*",
    "9884.164":     "s3close/round2/refute-u1 F1: L_* (all measured) on the cone, corrected",
    "3.2445":       "s3close/round2/refute-u1 F1: the cone-versus-phi_0 factor",
    "100.6122":     "s3close/round2/refute-u1 F3: sup over all phi of the sharpened C_R",
    "29.024165":    "s3close/round2/u1 sec.3.3: C_R at phi_0 = 30 deg, G_collar = 0",
    "98.289489":    "s3close/round2/u1 sec.6.8: C_R at phi = delta (the cone value)",
    "4.5253e14":    "s3close/assembly sec.3.3: L_* (all proved), a priori reference",
    "9.0506e14":    "s3close/assembly sec.3.3: log Lambda_* for that L_*",
    "34.0573":      "s3close/assembly sec.3.2: the a priori feedback exponent p c_*",
    "1.6218604":    "s3close/assembly sec.1.4: c_2 = 4 log(3/2)",
    "1.6227652":    "s3close/round2/u1 sec.5.3: c_2 = 2 log(3/2)/kappa, taper-only kappa",
    "0.8113826":    "s3close/round2/u1 sec.5.3: c_* at the taper-only kappa",
    "1.8377407":    "s3close/round2/u1 sec.1.1: lam_max at the taper-only c_*",
    "0.9199394":    "s3close/round2/u1 sec.1.1: r_h on [1,lam_max], taper-only profile",
    "5.5754e-04":   "s3close/assembly sec.1.3: 1 - 2 kappa for the axis taper alone",
    "5.5785e-04":   "s3close/assembly sec.1.3: the clock floor for the axis taper alone",
    "5.0531e-02":   "s3close/round2/u2 sec.9.1: 1 - 2 kappa for the campaign (D-C) datum",
    "6.0668328":    "the record's own measured Lemma T' conservatism factor",
    "87.166":       "write/refute-V-b: the corrected E_hess sizing factor",
    "1704.05546":   "the BFG preprint identifier",
}

CHECKS = []


def chk(name, got, want, rtol=1e-9, atol=0.0):
    ok = (abs(got - want) <= max(atol, rtol*abs(want))) if (
        isinstance(got, float) and isinstance(want, float)) else (got == want)
    CHECKS.append((name, got, want, ok))
    return ok


# --------------------------------------------------------------- t1: clock and datum
K = T1["clock"]
kap = T1["kappa"]["THEOREM_datum"]
chk("kappa_delta is P_h(1)/2 of the FULL profile", K["kappa_delta"], kap, 0.0)
chk("c_* = log(3/2)/kappa", K["c_star"], math.log(1.5)/kap, 1e-13)
chk("c_2 = 2 log(3/2)/kappa", K["c2"], 2*math.log(1.5)/kap, 1e-13)
chk("c_2 = 2 c_*", K["c2"], 2*K["c_star"], 1e-13)
chk("lam_max = exp(3 c_*/4)", K["lam_max_apriori"], math.exp(0.75*K["c_star"]), 1e-13)
chk("eps_delta = 1 - 2 kappa", K["eps_delta_1m2kappa"], 1 - 2*kap, 1e-12)
chk("floor = (1-2kappa)/(2kappa)", K["eps_floor_L_to_inf"], (1-2*kap)/(2*kap), 1e-12)
chk("c_2(taper only) = 2log(3/2)/kappa_taper", K["c2_taper_only"],
    2*math.log(1.5)/T1["kappa"]["taper_only"], 1e-13)
chk("c_2(D-C) = 2log(3/2)/kappa_DC", K["c2_campaign_DC"],
    2*math.log(1.5)/T1["kappa"]["campaign_DC"], 1e-13)
chk("floor(D-C) = (1-2k)/(2k)", K["floor_campaign_DC"],
    (1-2*T1["kappa"]["campaign_DC"])/(2*T1["kappa"]["campaign_DC"]), 1e-12)
chk("kappa(taper only) vs assembly 0.4997212305210886",
    T1["kappa"]["taper_only"], 0.4997212305210886, 1e-9)
chk("kappa(full) vs u2 sec.9.1's 0.4978224382", kap, 0.4978224382, 1e-9)
chk("c_2 vs refute-u1 F4's 1.6289540", K["c2"], 1.6289540, 1e-6)
chk("floor vs u2 sec.9.1's 0.0043742", K["eps_floor_L_to_inf"], 0.0043742, 1e-4)
chk("kappa(bang-bang) = 1/2 exactly", T1["kappa"]["bangbang"], 0.5, 1e-13)
chk("P_1(1.5)/1.5 = 1 (control)", T1["control_P1"]["P_1(1.5)/1.5"], 1.0, 1e-12)
chk("P_1(lam_max)/lam_max = 1 (control)", T1["control_P1"]["P_1(lam_max)/lam_max"], 1.0, 1e-12)
chk("P_h increasing on [1, lam_max] (min increment > 0)",
    T1["P_h_monotone"]["min_increment_on_[1,lam_max]"] > 0, True)
chk("P_h turns over above lam_max",
    T1["P_h_monotone"]["first_decrease_at_lam"] > T1["P_h_monotone"]["lam_max"], True)
DN = T1["datum_norms"]["THEOREM_datum"]
chk("E_0 = 1/sin delta for the piecewise-linear taper",
    DN["E0_sup_rho_eta"], 1.0/math.sin(7.5*math.pi/180), 1e-7)
chk("Gfrak_0 is set by the taper kink |d/dphi (1/sin phi)| at phi = delta",
    DN["sup_abs_Fprime"], math.cos(7.5*math.pi/180)/math.sin(7.5*math.pi/180)**2, 1e-6)
chk("E_0(D-C) reproduces u2's 7.66060196",
    T1["datum_norms"]["campaign_DC"]["E0_sup_rho_eta"], 7.66060196, 1e-8)
chk("Gfrak_0(D-C) reproduces u2's 20.9197070",
    T1["datum_norms"]["campaign_DC"]["G0_sup_rho2_grad_eta"], 20.9197070, 1e-6)
chk("Gfrak_0(theorem datum)/Gfrak_0(D-C) = 2.8057",
    DN["G0_sup_rho2_grad_eta"]/T1["datum_norms"]["campaign_DC"]["G0_sup_rho2_grad_eta"],
    2.8057495, 1e-6)
E = T1["energy"]
chk("C_E bang-bang control vs record 0.172403978", E["C_E_bangbang_control"], 0.172403978, 3e-8)
chk("C_E sharp (7.5,5) vs assembly",
    E["C_E_sharp_edges_delta7.5_dm5"], 0.17032563295410327, 1e-12)
chk("shift_sharp = 2 log s + (2/5) log C_E_sharp", E["shift_sharp"],
    2*math.log(E["s_record"]) + 0.4*math.log(E["C_E_sharp_edges_delta7.5_dm5"]), 1e-13)
chk("shift_sharp reproduces assembly (1.2)", E["shift_sharp"], 3.364345456912714, 1e-13)
chk("shift(theorem datum) = 2 log s + (2/5) log C_E_tanh", E["shift_THEOREM_datum"],
    2*math.log(E["s_record"]) + 0.4*math.log(E["C_E_THEOREM_datum_tanh_ramp_eps_r0.25"]), 1e-13)
chk("the radial ramp costs 0.44086 in the shift",
    E["shift_sharp"] - E["shift_THEOREM_datum"], 0.44085956, 1e-6)
chk("s_record = 1/sin delta", E["s_record"], 1.0/math.sin(7.5*math.pi/180), 1e-13)
chk("the D_l recursion agrees with the sharp closed form to 5e-5",
    E["D_l_recursion_vs_closed_form_maxrel_sharp_du2e-5"] < 6e-5, True)
chk("C_E(tanh) is converged in du to 5e-7", E["C_E_tanh_du_convergence_rel"] < 1e-6, True)
chk("C_E(tanh) is converged in l_max to 1e-5", E["C_E_tanh_lmax_convergence_rel"] < 1e-5, True)
chk("omega_0 in L^2 (finite)",
    math.isfinite(T1["bfg_hypothesis"]["omega0_L2sq_over_M2R3_THEOREM_datum"]), True)

# --------------------------------------------------------------- t2: kernel constants
KC = T2["kernel_constants"]
chk("C_K = 3/(8 pi^2)", KC["C_K_3_over_8pi2"], 3.0/(8*math.pi**2), 1e-15)
chk("|S^4| = 8 pi^2/3", KC["S4"], 8*math.pi**2/3, 1e-15)
chk("C_K |S^4| = 1 exactly", KC["C_K_times_S4"], 1.0, 1e-15)
chk("Riesz = pi^4/2", KC["riesz_pi4_over_2"], math.pi**4/2, 1e-15)
chk("Riesz quadrature agrees to 2e-4", KC["riesz_relerr"] < 2e-4, True)
chk("3 pi^2/16 = C_K pi^4/2", KC["coef_3pi2_over_16"],
    KC["C_K_3_over_8pi2"]*KC["riesz_pi4_over_2"], 1e-15)
chk("3 pi^2/16 numeric", KC["coef_3pi2_over_16"], 3*math.pi**2/16, 1e-15)
chk("R_A = (2^5 - 2^-5)^{1/5}", KC["R_A"], (2.0**5 - 2.0**-5)**0.2, 1e-15)
chk("2 R_A reproduces far/near's 3.999218", KC["two_R_A"], 3.999218, 1e-6)

# --------------------------------------------------------------- t2: (Gamma-off) control
CTL = T2["control_reproduce_u2"]
for tag, want in [("sigma*=1/2", 582.24), ("sigma*=0 (R^5)", 801.05)]:
    chk("C''(10^4) reproduces u2's %s" % tag, CTL[tag]["C_pp"], want, 1e-5)
chk("c_G(10^4, sigma=1/2) reproduces u2's 1.2636", CTL["sigma*=1/2"]["c_G"], 1.2636, 1e-4)
chk("c_G(10^4, sigma=0) reproduces u2's 1.2814", CTL["sigma*=0 (R^5)"]["c_G"], 1.2814, 1e-4)
chk("p(10^4, sigma=1/2) reproduces u2's 2.1068", CTL["sigma*=1/2"]["p"], 2.1068, 1e-4)
chk("p(10^4, sigma=0) reproduces u2's 2.1097", CTL["sigma*=0 (R^5)"]["p"], 2.1097, 1e-4)
chk("Chat_a(lam=1, sigma=1/2) reproduces assembly's C_a = 8.697888",
    CTL["Chat_a_lam1_sigma_half_at_c_G_eq_c"], 8.697888, 1e-6)
chk("L_Gamma_*(u2 datum, sigma=0) within 1% of u2's 5313.8",
    abs(CTL["L_gamma_star_u2_datum_sigma0"] - 5313.8)/5313.8 < 0.01, True)
chk("L_Gamma_*(u2 datum, sigma=1/2) within 2% of u2's 4997.6",
    abs(CTL["L_gamma_star_u2_datum_sigma_half"] - 4997.6)/4997.6 < 0.02, True)
chk("C'' = 2 Chat_a + Ghat + lambda (sigma=0)", CTL["sigma*=0 (R^5)"]["C_pp"],
    2*CTL["sigma*=0 (R^5)"]["Chat_a"] + CTL["sigma*=0 (R^5)"]["Ghat"] + 1.5, 1e-12)
chk("c_G = (lambda + C''/L) c (sigma=0)", CTL["sigma*=0 (R^5)"]["c_G"],
    (1.5 + CTL["sigma*=0 (R^5)"]["C_pp"]/1e4)*2*math.log(1.5), 1e-9)
chk("Ghat = (3 pi^2/16) e^{p c_G} Gfrak_0 (sigma=0)", CTL["sigma*=0 (R^5)"]["Ghat"],
    (3*math.pi**2/16)*math.exp(CTL["sigma*=0 (R^5)"]["p"]*CTL["sigma*=0 (R^5)"]["c_G"])
    * 20.9197070, 1e-9)

# --------------------------------------------------------------- t2: C_R
CR = T2["C_R"]
chk("C_R sup over ALL phi equals the sup over the cone (the argmax is interior)",
    CR["proved_supALLphi"]["C_R"], CR["proved_cone_phi>=delta"]["C_R"], 1e-8)
chk("(A-cone) costs nothing: sup over the value at phi_0 is under 1.01",
    CR["proved_supALLphi"]["C_R"]/CR["C_R_at_phi0_30deg_proved"] < 1.01, True)
chk("Chat_a(lam=1, phi_0) = assembly's 8.697888",
    CR["Chat_a_lam1_phi0_control_vs_assembly_8.697888"], 8.697888, 1e-6)
chk("Chat_a(phi_0, lam=3/2) = 13.046833 (u2 sec.7 table)",
    CR["Chat_a_at_phi0_proved"], 13.046833, 1e-6)
chk("ell_ramp = eps_r = 0.25", T2["ell_ramp_outer_edge"], 0.25, 1.5e-4)
chk("L_Gamma_*(theorem datum) exceeds L_Gamma_*(D-C) because Gfrak_0 is 2.8x larger",
    T2["L_gamma_star_theorem_datum"] > T2["L_gamma_star_if_G0_were_DC"], True)

# --------------------------------------------------------------- t3: the budget
D = T3["datum"]
chk("t3 reads t1's kappa", D["kappa_delta"], kap, 0.0)
chk("t3 reads t1's shift", D["logReE_shift_theorem_datum"], E["shift_THEOREM_datum"], 0.0)
LS, LAMS = T3["L_star"], T3["log_Lambda_star"]
for k in LS:
    chk("log Lambda_* = 2 L_* + shift  [%s]" % k,
        LAMS[k]["log_Lambda_star_theorem_datum"],
        2*LS[k] + D["logReE_shift_theorem_datum"], 1e-12)
chk("L_*(proved, eps<=1/2) is in [1e4, 1e6]", 1e4 <= LS["proved_eps<=0.5"] <= 1e6, True)
chk("L_*(proved, eps<=0.1) exceeds L_*(proved, eps<=1/2)",
    LS["proved_eps<=0.1"] > LS["proved_eps<=0.5"], True)
chk("the measured column beats the proved one by 737x",
    LS["proved_eps<=0.5"]/LS["measured_eps<=0.5"], 737.42, 1e-3)
chk("replacing C_a by 0.069 alone moves L_* by 2.8 percent",
    1 - LS["proved_Ca_measured_eps<=0.5"]/LS["proved_eps<=0.5"], 0.0282754, 1e-3)
chk("L_* beats ASSEMBLY's 4.5253e14 by more than 1e9",
    4.5253e14/LS["proved_eps<=0.5"] > 1e9, True)
chk("L_* is 5.94x refute-u1's corrected 52437.352",
    LS["proved_eps<=0.5"]/52437.352, 5.9407, 1e-4)
DOM = T3["dominance_at_Lstar_proved"]
chk("the r|grad a| term is 75.6 percent of C_R",
    DOM["C_R_pieces"]["term_C_share"], 0.75622, 1e-4)
chk("C''(L_*) = 2 Chat_a + Ghat + lambda", DOM["C_pp"],
    2*T3["column_constants"]["proved_L=300000"]["ca_sup"] + DOM["Ghat"] + 1.5, 2e-3)
SEN = T3["sensitivity"]
chk("the drive is LINEAR in C_R: doubling C_R doubles L_*",
    SEN["C_R x 2"]["L_star_0.5_proved"]/SEN["C_R x 1"]["L_star_0.5_proved"], 1.9910, 2e-3)
chk("and times 5 gives 4.96, not an exponential",
    SEN["C_R x 5"]["L_star_0.5_proved"]/SEN["C_R x 1"]["L_star_0.5_proved"], 4.9640, 2e-3)
chk("(H-K2) is now worth 1.4e-6 of L_*",
    SEN["C_K=hk2_proved_window"]["L_star_0.5_proved"]/SEN["C_K=1"]["L_star_0.5_proved"] - 1.0,
    1.355e-6, 5e-2)
chk("the self-consistent cap lam_max = 3/2 is worth 1.8161x",
    LS["proved_eps<=0.5"]/SEN["lam_max=3/2 (self-consistent cap)"]["L_star_0.5"], 1.81607, 1e-4)
chk("(Gamma-off)'s own fixed point does not exist at L = 10^4",
    T3["column_constants"]["proved_L=10000"] is None, True)
chk("and does exist at L = 10^5", T3["column_constants"]["proved_L=100000"] is not None, True)
chk("L_Gamma_* is between 10^4 and 10^5", 1e4 < T3["L_gamma_star"] < 1e5, True)
chk("t3's L_Gamma_* equals t2's", T3["L_gamma_star"], T2["L_gamma_star_theorem_datum"], 1e-9)

lm = D["lam_max"]
chk("lam_max^2 = e^{(3/2) c_*} exactly", lm**2, math.exp(1.5*D["c_star"]), 1e-13)
row = [r for r in T3["term_table_f1"] if r["column"] == "proved" and r["L"] == 1e6][0]
Gfp = 1.5 + row["C_pp"]/1e6
mu_closed = lm**3*(row["C_R"] + 2*1.5*math.log(lm))*(math.exp(Gfp*D["c_star"]) - 1.0)/Gfp
chk("mu(c) L closed form vs the ODE solver at L = 10^6", row["mu"]*1e6/mu_closed, 1.0, 1e-6)
CRp = T3["column_constants"]["proved_L=1e+09"]["C_R"]
chk("mu(c) L at L -> infinity = lam_max^3 (C_R + 2 lam_om log lam_max)(lam_max^2-1)/(3/2)",
    lm**3*(CRp + 2*1.5*math.log(lm))*(lm**2 - 1.0)/1.5, 8499.6195, 1e-6)
chk("eps_T' = (2pi/r_h)[3mu(1+muJ)/(2(1-mu)^5) + 3muJ/8] at L = 10^6", row["eps_Tprime"],
    (2*math.pi/D["r_h_1_to_lam_max"])*(3*row["mu"]*(1+row["muJ"])/(2*(1-row["mu"])**5)
                                       + 3*row["muJ"]/8), 1e-12)
chk("eps = (1 + log(1/(1-eps_v))/log(3/2))/((1-eps_a)(1-eps_delta)) - 1 at L = 10^6",
    row["eps"], (1 + math.log(1/(1-row["eps_v"]))/math.log(1.5))
    / ((1-row["eps_a"])*(1-row["eps_delta"])) - 1, 1e-12)
chk("eps_a = ell/L + eps_T' + ca_traj/(kappa L) at L = 10^6", row["eps_a"],
    row["eps_ell"] + row["eps_Tprime"] + row["eps_Ca"], 1e-12)
chk("eps_v = eps_bulk + E_hess + E_4 + E_tail at L = 10^6", row["eps_v"],
    row["eps_bulk"] + row["E_hess"] + row["E_4"] + row["E_tail"], 1e-12)
chk("eps_Ca = ca_traj/(kappa L) at L = 10^6", row["eps_Ca"], row["ca_traj"]/(kap*1e6), 1e-12)
chk("the viscous column is numerically irrelevant (eps_v < 1e-7 at L = 10^6)",
    row["eps_v"] < 1e-7, True)
chk("E_tail has vanished at L = 10^6", row["E_tail"], 0.0, 0.0, atol=1e-300)
chk("t_* = c_2 (1+eps)/(2 M L): c_2/2 = c_*", D["c2"]/2, D["c_star"], 1e-13)

# --------------------------------------------------------------- t4: derived quantities
chk("t4: C_E bang-bang relerr", T4["C_E_bangbang_relerr_vs_record"],
    abs(E["C_E_bangbang_control"] - 0.172403978)/0.172403978, 1e-12)
chk("t4: shift difference", T4["shift_difference_sharp_minus_tanh"],
    E["shift_sharp"] - E["shift_THEOREM_datum"], 1e-12)
chk("t4: Gfrak_0 ratio", T4["G0_ratio_theorem_over_DC"],
    DN["G0_sup_rho2_grad_eta"]/T1["datum_norms"]["campaign_DC"]["G0_sup_rho2_grad_eta"], 1e-12)
chk("t4: L_* against refute-u1's cone value", T4["Lstar_over_refute_u1_cone"],
    LS["proved_eps<=0.5"]/52437.352, 1e-12)
chk("t4: proved over measured", T4["proved_over_measured"],
    LS["proved_eps<=0.5"]/LS["measured_eps<=0.5"], 1e-12)
chk("t4: percent moved by the measured C_a", T4["percent_Ca_measured"],
    100.0*(1.0 - LS["proved_Ca_measured_eps<=0.5"]/LS["proved_eps<=0.5"]), 1e-12)
chk("t4: Ghat share of C_R, percent", T4["Ghat_share_of_C_R_percent"],
    100.0*DOM["C_R_pieces"]["term_C_share"], 1e-12)
chk("t4: gain from lam_max = 3/2", T4["gain_lam_max_three_halves"],
    LS["proved_eps<=0.5"]/SEN["lam_max=3/2 (self-consistent cap)"]["L_star_0.5"], 1e-12)
chk("t4: (H-K2) is worth", T4["hk2_worth_relative"],
    SEN["C_K=hk2_proved_window"]["L_star_0.5_proved"]
    / SEN["C_K=1"]["L_star_0.5_proved"] - 1.0, 1e-12)
chk("t4: L_star at the C_R-frozen base", T4["L_star_C_R_x1_base"],
    SEN["C_R x 1"]["L_star_0.5_proved"], 1e-12)
chk("t4: C_R scaling is linear (x1/2)", T4["ratio_C_R_half_over_x1"],
    SEN["C_R x 0.5"]["L_star_0.5_proved"]/SEN["C_R x 1"]["L_star_0.5_proved"], 1e-12)
chk("t4: C_R scaling is linear (x2)", T4["ratio_C_R_double_over_x1"],
    SEN["C_R x 2"]["L_star_0.5_proved"]/SEN["C_R x 1"]["L_star_0.5_proved"], 1e-12)
chk("t4: C_R scaling is linear (x5)", T4["ratio_C_R_times5_over_x1"],
    SEN["C_R x 5"]["L_star_0.5_proved"]/SEN["C_R x 1"]["L_star_0.5_proved"], 1e-12)
chk("t4: mu(c) L at L -> infinity", T4["mu_times_L_at_L_infinity"], 8499.6195, 1e-6)

# --------------------------------------------------------------- report
npass = sum(1 for c in CHECKS if c[3])
nfail = len(CHECKS) - npass
for name, got, want, ok in CHECKS:
    if not ok:
        print("FAIL: %-70s got %r want %r" % (name, got, want))
print("%d CHECKS, %d PASS, %d FAIL" % (len(CHECKS), npass, nfail))


# --------------------------------------------------------------- document audit
def flat(o, out):
    if isinstance(o, dict):
        for v in o.values():
            flat(v, out)
    elif isinstance(o, list):
        for v in o:
            flat(v, out)
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        out.append(float(o))
    return out


def ulp_of(tok):
    """half a unit in the last displayed place: the only honest matching tolerance"""
    mant, _, ex = tok.lower().partition("e")
    d = len(mant.split(".")[1]) if "." in mant else 0
    return 0.51*10.0**(-d)*10.0**(int(ex) if ex else 0)


POOL = [p for p in flat(T1, []) + flat(T2, []) + flat(T3, []) + flat(T4, [])
        if math.isfinite(p)]
try:
    doc = open("THEOREM_S3.md").read()
except OSError:
    print("DOC AUDIT: THEOREM_S3.md not present yet")
    sys.exit(0 if nfail == 0 else 1)
doc = doc.split("## GATE")[0]
doc = re.sub(r"§\s*\d+(\.\d+)*", " ", doc)         # section references
doc = re.sub(r"(?m)^#{1,6}\s+\d+(\.\d+)*", " ", doc)   # numbered headings
doc = re.sub(r"arXiv:\S+", " ", doc)                    # the BFG preprint identifier
doc = re.sub(r"sha256 `[0-9a-f.]+`", " ", doc)          # provenance hashes
TOKRE = r"(?<![\w.])\d+\.\d+(?:[eE][-+]?\d+)?|(?<![\w.])\d+[eE][-+]?\d+"
toks = sorted(set(re.findall(TOKRE, doc)))
traced, foreign, untraced = [], [], []
for t in toks:
    if t in FOREIGN:
        foreign.append(t)
        continue
    v, tol = float(t), ulp_of(t)
    (traced if any(abs(p - v) <= tol for p in POOL) else untraced).append(t)
print("DOC AUDIT: %d distinct numeric tokens, %d traced to this folder's JSONs, "
      "%d quoted from a named source, %d untraced"
      % (len(toks), len(traced), len(foreign), len(untraced)))
if untraced:
    print("  untraced:", " ".join(untraced))
sys.exit(0 if (nfail == 0 and not untraced) else 1)

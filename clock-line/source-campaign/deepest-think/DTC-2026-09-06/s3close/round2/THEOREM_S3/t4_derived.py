"""
t4 -- the derived quantities quoted inline in THEOREM_S3.md's prose (ratios, percentages,
differences), computed from t1/t2/t3's stored JSONs so that every one of them is traceable.

Outputs -> t4_results.json
"""
import json, math

T1 = json.load(open("t1_results.json"))
T2 = json.load(open("t2_results.json"))
T3 = json.load(open("t3_results.json"))

E = T1["energy"]
LS = T3["L_star"]
SEN = T3["sensitivity"]
DOM = T3["dominance_at_Lstar_proved"]

OUT = {}

# energy
OUT["C_E_bangbang_relerr_vs_record"] = abs(E["C_E_bangbang_control"] - 0.172403978)/0.172403978
OUT["shift_difference_sharp_minus_tanh"] = E["shift_sharp"] - E["shift_THEOREM_datum"]
OUT["log_Lambda_star_proved_half_with_sharp_shift"] = (
    T3["log_Lambda_star"]["proved_eps<=0.5"]["log_Lambda_star_if_sharp_edges"])

# datum
OUT["G0_ratio_theorem_over_DC"] = (
    T1["datum_norms"]["THEOREM_datum"]["G0_sup_rho2_grad_eta"]
    / T1["datum_norms"]["campaign_DC"]["G0_sup_rho2_grad_eta"])

# the (Gamma-off) control
CTL = T2["control_reproduce_u2"]
OUT["L_gamma_star_u2_relerr_vs_5313.8"] = abs(CTL["L_gamma_star_u2_datum_sigma0"]
                                              - 5313.8)/5313.8

# reach, against the earlier notes
OUT["Lstar_over_refute_u1_cone"] = LS["proved_eps<=0.5"]/52437.352
OUT["assembly_Lstar_over_this_Lstar"] = 4.5253e14/LS["proved_eps<=0.5"]
OUT["proved_over_measured"] = LS["proved_eps<=0.5"]/LS["measured_eps<=0.5"]

# what a change buys, as a factor on L_*(proved, eps <= 1/2)
base = LS["proved_eps<=0.5"]
OUT["factor_Ca_measured"] = LS["proved_Ca_measured_eps<=0.5"]/base
OUT["percent_Ca_measured"] = 100.0*(1.0 - LS["proved_Ca_measured_eps<=0.5"]/base)
OUT["factor_C_R_half"] = SEN["C_R x 0.5"]["L_star_0.5_proved"]/base
OUT["factor_C_R_double"] = SEN["C_R x 2"]["L_star_0.5_proved"]/base
OUT["factor_C_R_times5"] = SEN["C_R x 5"]["L_star_0.5_proved"]/base
OUT["ratio_C_R_half_over_x1"] = (SEN["C_R x 0.5"]["L_star_0.5_proved"]
                                 / SEN["C_R x 1"]["L_star_0.5_proved"])
OUT["L_star_C_R_x1_base"] = SEN["C_R x 1"]["L_star_0.5_proved"]
OUT["ratio_C_R_double_over_x1"] = (SEN["C_R x 2"]["L_star_0.5_proved"]
                                   / SEN["C_R x 1"]["L_star_0.5_proved"])
OUT["ratio_C_R_times5_over_x1"] = (SEN["C_R x 5"]["L_star_0.5_proved"]
                                   / SEN["C_R x 1"]["L_star_0.5_proved"])
OUT["factor_lam_max_three_halves"] = SEN["lam_max=3/2 (self-consistent cap)"]["L_star_0.5"]/base
OUT["gain_lam_max_three_halves"] = base/SEN["lam_max=3/2 (self-consistent cap)"]["L_star_0.5"]
OUT["factor_C_K_one"] = SEN["C_K=1"]["L_star_0.5_proved"]/base
OUT["hk2_worth_relative"] = (SEN["C_K=hk2_proved_window"]["L_star_0.5_proved"]
                             / SEN["C_K=1"]["L_star_0.5_proved"] - 1.0)

# the drive decomposition
OUT["Ghat_share_of_C_R_percent"] = 100.0*DOM["C_R_pieces"]["term_C_share"]

# the closed form of mu(c) L at L -> infinity
D = T3["datum"]
lm = D["lam_max"]
CRinf = T3["column_constants"]["proved_L=1e+09"]["C_R"]
OUT["mu_times_L_at_L_infinity"] = lm**3*(CRinf + 2*1.5*math.log(lm))*(lm**2 - 1.0)/1.5

with open("t4_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True)
print(json.dumps(OUT, indent=1, sort_keys=True))

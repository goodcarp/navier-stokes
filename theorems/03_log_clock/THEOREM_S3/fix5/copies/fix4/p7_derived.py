"""
p7 -- every ratio, percentage and difference quoted in the prose of FIX4.md and
      THEOREM_S3_v2.md, derived here from the stored JSONs so that none is typed.

Outputs -> p7_results.json
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
J = {k: json.load(open(os.path.join(HERE, k + "_results.json")))
     for k in ["p1", "p2", "p3", "p4", "p5", "p6"]}
SCAN = json.load(open(os.path.join(HERE, "p4_scan.json")))

R = J["p1"]["rows"]["sigma=0.02"]
K = J["p1"]["rows"]["kinked"]
P4 = J["p4"]
O = {}

# fix2 / ADDENDUM_3's headline, the number this sheet is measured against (FOREIGN, quoted)
FIX2 = {"L_star_proved_cap": 1887786.7127, "L_star_proved_0.1": 1977309.3012,
        "L_star_Ca_cap": 1848687.4, "L_star_Ca_0.1": 1932180.6,
        "L_star_meas_cap": 1538.37, "L_star_meas_0.1": 1701.58,
        "logLam_proved_cap": 3775576.3488, "logLam_proved_0.1": 3954621.5258,
        "kappa": 0.4978224182, "c_star": 0.8144774, "c2": 1.6289547,
        "floor": 4.3741737e-03, "E0": 7.66129757554039, "G0": 58.69548054098106,
        "shift": 2.9234859, "eps_cap": 0.1943662078602755, "lam_mono": 2.0769162,
        "r_h": 0.9186364112, "J_sup": 96.249844, "J_plateau": 75.251104,
        "K2hat_f4": {"1.0": 19.169, "1.5": 25.418, "1.8420112": 33.710, "2.0741": 42.008}}
O["FOREIGN_fix2_addendum3"] = FIX2

O["datum_moves"] = {
    "kappa_ratio_moll_over_kinked": R["kappa_delta"]/K["kappa_delta"],
    "kappa_drop_percent": 100.0*(1.0 - R["kappa_delta"]/K["kappa_delta"]),
    "floor_ratio": R["eps_floor"]/K["eps_floor"],
    "c2_ratio": R["c2"]/K["c2"],
    "E0_ratio": R["E0"]/K["E0"],
    "G0_ratio": R["Gfrak0"]/K["Gfrak0"],
    "G0_gain_factor": K["Gfrak0"]/R["Gfrak0"],
    "grad_inf_ratio": R["grad_eta0_inf"]/K["grad_eta0_inf"],
    "N_sigma_overshoot_percent": 100.0*(R["N_sigma"] - 1.0),
    "amplitude_at_phi0_deficit_percent": 100.0*(1.0 - R["A_sigma_at_phi0"]),
    "cap_ratio": R["eps_cap_certified"]/K["eps_cap_certified"],
}

O["L_star_moves"] = {
    "proved_cap_gain_vs_fix2": FIX2["L_star_proved_cap"]/P4["L_star"]["proved_eps<=cap"],
    "proved_0.1_gain_vs_fix2": FIX2["L_star_proved_0.1"]/P4["L_star"]["proved_eps<=0.1"],
    "Ca_cap_gain_vs_fix2": FIX2["L_star_Ca_cap"]/P4["L_star"]["proved_Ca_measured_eps<=cap"],
    "Ca_0.1_gain_vs_fix2": FIX2["L_star_Ca_0.1"]/P4["L_star"]["proved_Ca_measured_eps<=0.1"],
    "measured_cap_loss_vs_fix2": P4["L_star"]["measured_eps<=cap"]/FIX2["L_star_meas_cap"],
    "measured_0.1_loss_vs_fix2": P4["L_star"]["measured_eps<=0.1"]/FIX2["L_star_meas_0.1"],
    "logLam_proved_cap_gain": FIX2["logLam_proved_cap"]/P4["log_Lambda_star"]["proved_eps<=cap"],
}

CTL = P4["control_kinked_vs_fix2"]
O["control"] = {
    "kinked_cap_here_over_fix2": CTL["eps<=0.1943662 (fix2 cap)"]/FIX2["L_star_proved_cap"],
    "kinked_cap_percent_high": 100.0*(CTL["eps<=0.1943662 (fix2 cap)"]/FIX2["L_star_proved_cap"] - 1.0),
    "kinked_0.1_percent_high": 100.0*(CTL["eps<=0.1"]/FIX2["L_star_proved_0.1"] - 1.0),
    "r_h_here_over_fix2": CTL["r_h_at_c_star_here"]/FIX2["r_h"],
    "r_h_relative_deficit": 1.0 - CTL["r_h_at_c_star_here"]/FIX2["r_h"],
    "kappa_rel_vs_fix2": abs(K["kappa_delta"] - FIX2["kappa"])/FIX2["kappa"],
    "E0_rel_vs_closed_form": abs(K["E0"] - FIX2["E0"])/FIX2["E0"],
    "G0_rel_vs_closed_form": abs(K["Gfrak0"] - FIX2["G0"])/FIX2["G0"],
    "lam_mono_rel_vs_fix2": abs(K["lam_mono"] - FIX2["lam_mono"])/FIX2["lam_mono"],
    "J_sup_rel_vs_fix2": abs(J["p3"]["C_J"]["kinked"]["J_sup_over_rho"] - FIX2["J_sup"])/FIX2["J_sup"],
    "J_plateau_rel_vs_fix2": abs(J["p3"]["C_J"]["kinked"]["J_plateau"] - FIX2["J_plateau"])/FIX2["J_plateau"],
    "K2_kinked_rel": {str(r["lam"]): abs(r["K2"] - FIX2["K2hat_f4"][str(r["lam"])])
                      / FIX2["K2hat_f4"][str(r["lam"])]
                      for r in J["p3"]["C_K2_kinked_control"]["rows"]},
    "C_E_sharp_rel_vs_assembly": abs(J["p2"]["control_kinked"]["C_E_sharp_edges"]
                                     - 0.17032563295410327)/0.17032563295410327,
    "C_E_tanh_rel_vs_THEOREM_S3": abs(J["p2"]["control_kinked"]["C_E_tanh_ramp"]
                                      - 0.05657477)/0.05657477,
    "shift_kinked_rel_vs_THEOREM_S3": abs(J["p2"]["control_kinked"]["shift_tanh_kinked"]
                                          - FIX2["shift"])/FIX2["shift"],
    "L2_kinked_rel_vs_THEOREM_S3": abs(J["p2"]["bfg_hypothesis"]["omega0_L2sq_over_M2R3_kinked"]
                                       - 1.4785137)/1.4785137,
    "Psi0_ac_rel_vs_fix3": abs(K["Psi0"]*1.0 - 171.4257)/171.4257,
}

# fix3's table is UNNORMALISED: multiply back by N_sigma to compare
FIX3 = {"0.05": {"G0": 25.2888, "Psi0": 244.07, "kappa": 0.4983882748},
        "0.02": {"G0": 36.5308, "Psi0": 495.16, "kappa": 0.4979100241},
        "0.01": {"G0": 43.9502, "Psi0": 901.42, "kappa": 0.4978442142},
        "0.005": {"G0": 49.5547, "Psi0": 1709.06, "kappa": 0.4978278606},
        "0.002": {"G0": 54.2132, "Psi0": 4130.11, "kappa": 0.4978232887}}
O["FOREIGN_fix3_unnormalised"] = FIX3
un = {}
for s, d in FIX3.items():
    r = J["p1"]["rows"]["sigma=%g" % float(s)]
    N = r["N_sigma"]
    un[s] = {"G0_unnorm_here": r["Gfrak0_scan"]*N, "G0_rel": abs(r["Gfrak0_scan"]*N - d["G0"])/d["G0"],
             "Psi0_unnorm_here": r["Psi0"]*N, "Psi0_rel": abs(r["Psi0"]*N - d["Psi0"])/d["Psi0"],
             "kappa_unnorm_here": J["p1"]["unnormalised_control_vs_fix3"]["sigma=%g" % float(s)]["kappa_unnormalised"],
             "kappa_rel": abs(J["p1"]["unnormalised_control_vs_fix3"]["sigma=%g" % float(s)]["kappa_unnormalised"] - d["kappa"])/d["kappa"]}
O["control_vs_fix3_unnormalised"] = un

O["sigma_scan_relative_to_kinked"] = {
    k: {"eps<=0.1": (v["eps<=0.1"]/SCAN["rows"]["sigma=None"]["eps<=0.1"]) if v["eps<=0.1"] not in (None, "None") else None,
        "eps<=cap": (v["eps<=cap"]/SCAN["rows"]["sigma=None"]["eps<=cap"]) if v["eps<=cap"] not in (None, "None") else None}
    for k, v in SCAN["rows"].items()}

O["eps_a_share"] = {
    r["column"] + "_L=%g" % r["L"]: (r["eps_a_mollified"]/r["eps_budget"])
    for r in J["p5"]["table"] if "eps_a_mollified" in r}

O["K2_slack"] = {"majorant_over_max": J["p3"]["C_K2_majorant"]["slack"],
                 "max_lam_range": J["p3"]["C_K2_majorant"]["max_over_lam_mollified"],
                 "mollified_over_kinked_at_lam_1":
                     J["p3"]["C_K2_mollified_f4_L10"][0]["K2"]
                     / J["p3"]["C_K2_kinked_control"]["rows"][0]["K2"]}

O["shift_and_logLam"] = {
    "delta_shift": J["p2"]["mollified"]["shift_logReE"] - J["p2"]["control_kinked"]["shift_tanh_kinked"],
    "logLam_proved_cap_check": 2*P4["L_star"]["proved_eps<=cap"] + J["p2"]["mollified"]["shift_logReE"],
    "logLam_proved_0.1_check": 2*P4["L_star"]["proved_eps<=0.1"] + J["p2"]["mollified"]["shift_logReE"],
}

O["direct_sensitivity"] = J["p6"]["direct_sensitivity"]

json.dump(O, open(os.path.join(HERE, "p7_results.json"), "w"), indent=1, sort_keys=True,
          default=str)
print(json.dumps(O, indent=1, sort_keys=True, default=str))

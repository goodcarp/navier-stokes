"""a7 -- quantities quoted inline in ASSEMBLY.md's prose, computed so they are traceable."""
import json, math
A1 = json.load(open("a1_results.json")); A2 = json.load(open("a2_results.json"))
kap = A1["kappa_delta"]["7.5"]
O = {}
O["ell_loss_const_f0"] = math.log(2.0) + 3.0*math.log(1.5)
O["ell_loss_f1"] = math.log(2.0) + math.log(2.0) + 3.0*math.log(1.5)
O["eps_ell_L160_f1"] = O["ell_loss_f1"]/160.0
O["C_a_proved_phi30"] = 0.291999 + 0.014754 + 3.999218/math.sin(math.pi/6) + math.pi/8
O["eps_Ca_L160_proved_f1"] = 1.5*O["C_a_proved_phi30"]/(kap*160.0)
O["eps_Ca_L160_measured_f1"] = 1.5*0.069/(kap*160.0)
O["C_E_l1_share_pct"] = 100.0*A1["C_E_l1_share"]
O["lambda_at_c_0p25"] = math.exp(kap*0.25)
O["equator_margin_minus_dm_deg"] = (90.0 - A1["material_point"]["phi0=30"]["lam=1.5"]["phi_deg"]) - 5.0
O["eps_floor_delta"] = (1.0-2.0*kap)/(2.0*kap)
O["d_f1_taper"] = 2.0*math.sin(math.pi/6 - 7.5*math.pi/180.0)
O["d_f1_taper_with_eps025"] = 2.0 - math.exp(0.25)
O["d_reduction_pct"] = 100.0*(1.0 - (2.0-math.exp(0.25))/(2.0*math.sin(math.pi/6-7.5*math.pi/180.0)))
O["C_a_ratio_proved_over_measured"] = O["C_a_proved_phi30"]/0.069
O["geom_product"] = 1.5*(9.0/4.0)
O["Lemma_Tprime_rel_const"] = 15.0*math.pi/4.0
O["conservatism_factor_LemmaTprime"] = 2.1181705106873974/0.34913942333376546
O["C_E_L10_vs_Linf_relgap"] = abs(A1["C_E_table"]["delta=7.5,dm=5,eps=0,L=10"]
                                  - A1["C_E_table"]["delta=7.5,dm=5,eps=0"]) \
                              / A1["C_E_table"]["delta=7.5,dm=5,eps=0"]
json.dump(O, open("a7_results.json", "w"), indent=1, sort_keys=True)
print(json.dumps(O, indent=1, sort_keys=True))

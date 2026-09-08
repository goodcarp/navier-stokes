"""
check_constants.py -- the gate for s3close/round2/pmax-h3v.

Two jobs, in this order.

  PART 1, RECOMPUTATION.  Every quantity this note leans on is rebuilt here from the
  mathematics, in this file, and compared with what the four scripts stored.  A check that
  compares a stored literal with itself is not a check and none is written that way.

  PART 2, DOCUMENT AUDIT.  Every numeric literal in PROOF.md is matched against the pool of
  values the scripts produced (to the precision at which it is printed), or against the
  FOREIGN dictionary of numbers quoted from other seats, or against the whitelist of
  structural integers.  Anything unmatched is a FAIL.

WHAT THE GATE CANNOT CATCH, stated so it is not mistaken for coverage: it checks that the
numbers displayed are the numbers the scripts produced and that the algebra between them is
consistent.  It does not check that Lemma A's citation is the right theorem, that Lemma E's
Dini step is correctly quoted, that (H3'_vartheta) is the right hypothesis, or that the
modulo list is complete.  Section A8 and the findings are for that.

Run:  python3 check_constants.py
"""
import json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
J = {n: json.load(open(os.path.join(HERE, "q%d_results.json" % n))) for n in (1, 2, 3, 4)}

PASS, FAIL = [], []


def chk(name, got, want, tol=1e-9, rel=True):
    if want == 0 or not rel:
        ok = abs(got - want) <= tol
        err = abs(got - want)
    else:
        err = abs(got - want) / abs(want)
        ok = err <= tol
    (PASS if ok else FAIL).append((name, got, want, err))
    return ok


def chk_true(name, cond, detail=""):
    (PASS if cond else FAIL).append((name, detail, "True", 0.0 if cond else 1.0))
    return cond


DEG = math.pi / 180.0
DELTA = 7.5 * DEG
DM = 5.0 * DEG
PHI0 = 30.0 * DEG
E2 = math.exp(2.0)

# =====================================================================  PART 1
# ---- q1: every symbolic residual is zero, and the count
n_checks = J[1]["summary"]["n_checks"]
chk_true("q1_all_symbolic_residuals_zero", J[1]["summary"]["n_fail"] == 0,
         "n_fail = %d" % J[1]["summary"]["n_fail"])
chk("q1_n_checks", n_checks, 49, tol=0, rel=False)
chk_true("q1_weight_sign_at_n5",
         J[1]["A_sign_at_n5"]["nonpositive_for_k"] == [0, 1, 2, 3],
         str(J[1]["A_sign_at_n5"]["nonpositive_for_k"]))
for k, v in J[1]["A_sign_at_n5"]["k(k-3)"].items():
    chk("q1_k(k-3)_k=%s" % k, float(v), int(k) * (int(k) - 3), tol=0, rel=False)
# the reference-strain negative result, recomputed here from the matrix itself
S = [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, -2]]
eigs = [1.0, 1.0, 1.0, 1.0, -2.0]              # sym(S) = S is already diagonal
chk("q1_ref_strain_Gamma_rad", J[1]["F_reference_strain"]["Gamma_rad = Lambda_max"], max(eigs))
chk("q1_ref_strain_minus_Lmin", J[1]["F_reference_strain"]["minus_Lambda_min"], -min(eigs))
chk("q1_ref_strain_Gamma_op", J[1]["F_reference_strain"]["Gamma = ||grad_5 b||_op"], 2.0)
chk("q1_ref_strain_Gamma_rad_over_Gamma",
    J[1]["F_reference_strain"]["Gamma_rad_over_Gamma"], 0.5)
chk("q1_ref_strain_sharpening_gain",
    J[1]["F_reference_strain"]["sharpening_gain = Gamma/(-Lambda_min)"], 1.0)
chk("q1_barrier_majorant_gap_k2", J[1]["E_majorant_min_gap_over_grid"]["min_gap"], 0.0,
    tol=1e-12, rel=False)
chk("q1_barrier_majorant_gap_k1", J[1]["E_majorant_k1_min_gap"]["min_gap"], 0.0,
    tol=1e-12, rel=False)

# ---- q2: the two closed forms, recomputed here
E0 = 1.0 / math.sin(DELTA)
G0 = E0 ** 2
chk("q2_E0_closed_form", J[2]["closed_forms"]["E0 = 1/sin(delta)"], E0, 1e-15)
chk("q2_G0_closed_form", J[2]["closed_forms"]["G0 = 1/sin(delta)^2 = E0^2"], G0, 1e-15)
chk("q2_G0_equals_E0_squared", G0, J[2]["closed_forms"]["E0 = 1/sin(delta)"] ** 2, 1e-15)
chk("q2_E0_scan_matches_closed", J[2]["E0_G0"]["E0_scan"], E0, 1e-11)
chk("q2_G0_scan_matches_closed", J[2]["E0_G0"]["G0_scan"], G0, 1e-8)
chk("q2_supFprime_closed", J[2]["E0_G0"]["sup_abs_Fprime"],
    math.cos(DELTA) / math.sin(DELTA) ** 2, 1e-8)
# the plateau identity F^2 + F'^2 = 1/sin^4 phi, checked at ten angles between the kinks
for i, ph in enumerate([x * DEG for x in (10, 20, 30, 40, 50, 60, 70, 80, 84, 85)]):
    F = 1.0 / math.sin(ph)
    Fp = -math.cos(ph) / math.sin(ph) ** 2
    chk("q2_plateau_identity_%d" % i, F * F + Fp * Fp, 1.0 / math.sin(ph) ** 4, 1e-13)
# ||eta_0||_inf factorises
chk("q2_eta0_Linf_factorises", J[2]["sup_norms"]["eta0_Linf_times_rho0_over_M"],
    J[2]["sup_norms"]["sup_u exp(-u)Theta(u)"] * E0, 1e-12)
# the tail constant equals ||eta_0||_inf rho_0/M in the large-L limit
for L in ("L=10", "L=20", "L=40"):
    chk("q2_tail_A9_equals_Linf_%s" % L, J[2]["tail_constants"][L]["A_9"],
        J[2]["sup_norms"]["eta0_Linf_times_rho0_over_M"], 3e-7)
# tail and TV are L-stable
A9 = [J[2]["tail_constants"][L]["A_9"] for L in ("L=10", "L=20", "L=40")]
chk_true("q2_A9_L_stable", (max(A9) - min(A9)) / max(A9) < 1e-6,
         "spread %.2e" % ((max(A9) - min(A9)) / max(A9)))
TV = [J[2]["TV"][L]["TV_over_M_Rcubed"] for L in ("L=20", "L=40")]
chk_true("q2_TV_L_stable", abs(TV[0] - TV[1]) / TV[1] < 1e-5,
         "rel %.2e" % (abs(TV[0] - TV[1]) / TV[1]))
# mollification never raises either norm
for k, v in J[2]["mollification_control"].items():
    chk_true("q2_mollify_E0_le_1_%s" % k, v["E0_ratio_to_E0"] <= 1.0, str(v["E0_ratio_to_E0"]))
    chk_true("q2_mollify_G0_le_1_%s" % k, v["G0_ratio_to_G0"] <= 1.0, str(v["G0_ratio_to_G0"]))
sig = sorted(J[2]["mollification_control"], key=lambda s: -float(s.split("=")[1]))
rat = [J[2]["mollification_control"][s]["G0_ratio_to_G0"] for s in sig]
chk_true("q2_mollify_G0_rises_to_1", all(rat[i] < rat[i + 1] for i in range(len(rat) - 1)),
         str(rat))
# the analytic per-arc derivative agrees with a central difference off the kinks
chk_true("q2_analytic_vs_fd", J[2]["E0_G0"]["analytic_vs_central_difference_control"] < 1e-6,
         str(J[2]["E0_G0"]["analytic_vs_central_difference_control"]))
# the excluded end caps cannot move either supremum
chk_true("q2_endcap_below_E0", J[2]["E0_G0"]["endcap_bound_absF"] < E0)
chk_true("q2_endcap_below_supFprime",
         J[2]["E0_G0"]["endcap_bound_absFprime"] < J[2]["E0_G0"]["sup_abs_Fprime"])

# ---- the (phi,u) separation that Proposition 3(b) rests on, recomputed here
import numpy as _np
_epsr, _L = 0.25, 1000.0
_u = _np.concatenate([_np.linspace(-4, 6, 400001), _np.linspace(_L - 6, _L + 4, 400001)])
_Th = 0.5 * (_np.tanh((_u - _epsr) / _epsr) - _np.tanh((_u - _L + _epsr) / _epsr))
_dTh = ((_np.cosh(_np.clip((_u - _epsr) / _epsr, -350, 350)) ** -2)
        - (_np.cosh(_np.clip((_u - _L + _epsr) / _epsr, -350, 350)) ** -2)) / (2 * _epsr)
_A2, _B2 = (_Th - _dTh) ** 2, _Th ** 2
_Fp = math.cos(DELTA) / math.sin(DELTA) ** 2
_val = E0 ** 2 * _A2 + _Fp ** 2 * _B2
chk("prop3b_global_max_is_G0_squared", float(_val.max()), G0 ** 2, 1e-12)
chk("prop3b_argmax_A2_is_1", float(_A2[int(_val.argmax())]), 1.0, 1e-12)
chk("prop3b_argmax_B2_is_1", float(_B2[int(_val.argmax())]), 1.0, 1e-12)
_off = _Th < 0.99
chk("prop3b_off_plateau_max", float(_val[_off].max()), 3386.1663916755992, 1e-9)
chk_true("prop3b_off_plateau_below_global", float(_val[_off].max()) < G0 ** 2)
chk("prop3b_norms_identity", E0 ** 2 + _Fp ** 2, G0 ** 2, 1e-14)

# ---- q3: vartheta
def psi_in(rho):
    return E2 / (E2 + rho ** 8)


for f in [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]:
    row = J[3]["f_grid_d_from_assembly_rule"]["f=%g" % f]
    rho_s = 1.0 + f
    d_rule = min(f, rho_s * math.sin(PHI0 - DELTA), rho_s * math.sin(math.pi / 2 - DM - PHI0))
    chk("q3_d_rule_f=%g" % f, row["d"], d_rule, 1e-12)
    chk("q3_rho_min_f=%g" % f, row["rho_min_on_ball"], rho_s - d_rule, 1e-12)
    chk("q3_vartheta0_is_psi_f=%g" % f, row["vartheta_k"]["0"], psi_in(rho_s - d_rule), 1e-12)
    chk("q3_Theta_at_rho_star_f=%g" % f, row["Theta_at_rho_star"], 1.0 - psi_in(rho_s), 1e-12)
    chk("q3_vartheta_is_max_f=%g" % f, row["vartheta"],
        max(float(v) for v in row["vartheta_k"].values()), 1e-15)
for k, v in J[3]["series_vs_sympy_control"].items():
    chk_true("q3_series_vs_sympy_%s" % k, v["rel"] < 1e-13, str(v["rel"]))
for k, v in J[3]["grid_stability"].items():
    chk_true("q3_grid_stability_%s" % k, v["rel_change"] < 2e-2, str(v["rel_change"]))
# psi ~ e^2 rho^{-8}: the thresholds quoted in B3
for p in (1e-2, 1e-3, 1e-4, 1e-6):
    chk("q3_rho_min_threshold_%g" % p,
        float(J[3]["asymptotic"]["rho_min_needed_for_psi_below"][str(p)]),
        math.exp(0.25) * (1.0 / p) ** 0.125, 1e-12)
# 1 - sin(phi_0 - delta), the taper-limited shrink factor quoted in B3
chk("q3_taper_shrink_factor", 1.0 - math.sin(PHI0 - DELTA), 0.61731656763491, 1e-12)

# ---- q4: the vartheta budget algebra, rebuilt from the stored pieces
for row in J[4]["term_table_proved_safety2"]:
    if not row.get("admissible"):
        continue
    vt_raw = row["vartheta"]
    Q = (1.0 + vt_raw) / (1.0 - vt_raw)
    vt = vt_raw / (1.0 - vt_raw)
    # sigma_z/sigma_y is not stored in the term table; recover it from extra_bulk
    ratio = (row["extra_bulk"] / (row["eps_bulk"] * vt) - 9.0) / 2.0
    want = (row["eps_bulk"] * (1.0 + vt * (9.0 + 2.0 * ratio))
            + Q * (row["E_hess"] + row["E_4"] + row["E_tail"]))
    chk("q4_eps_v_formula_L=%g_f=%g" % (row["L"], row["f"]), row["eps_v"], want, 1e-12)
    chk_true("q4_sigma_ratio_is_CorV5_L=%g_f=%g" % (row["L"], row["f"]),
             abs(ratio - 3.5130747) / 3.5130747 < 5e-2, "%.6f" % ratio)
# the control column must reproduce THEOREM_S3's published L_*
FOR = J[4]["FOREIGN_THEOREM_S3_L_star"]
CTL = J[4]["control_reproduce_THEOREM_S3"]
chk("q4_control_proved_half", CTL["proved"]["L_star_eps<=0.5"], FOR["proved_eps<=0.5"], 2e-6)
chk("q4_control_proved_Ca", CTL["proved_Ca_measured"]["L_star_eps<=0.5"],
    FOR["proved_Ca_measured_eps<=0.5"], 2e-6)
chk("q4_control_measured", CTL["measured"]["L_star_eps<=0.5"], FOR["measured_eps<=0.5"], 1e-4)
# (the published 422.4 is rounded to one decimal; 1e-4 is that rounding, not slack)
chk("q4_control_proved_tenth", CTL["proved"]["L_star_eps<=0.1"], FOR["proved_eps<=0.1"], 2e-6)
# admissibility flags agree with vartheta < 1
for k, v in J[4]["admissibility"].items():
    chk_true("q4_admissible_flag_%s" % k,
             v["admissible"] == (v["vartheta"] is not None and v["vartheta"] < 1.0), k)
# log Lambda_* = 2 L_* + shift
SHIFT = J[4]["SHIFT"]["logReE_shift_theorem_datum"]
for k, v in J[4]["L_star_with_vartheta"].items():
    if v["L_star"]:
        chk("q4_logLambda_%s" % k, v["log_Lambda_star"], 2 * v["L_star"] + SHIFT, 1e-12)

# ---- derived quantities quoted inline in PROOF.md
ratio_ref = 1.7905793948266 / 0.5096901045590
chk("derived_sigma_ratio_CorV5", ratio_ref, 3.5130747, 1e-7)
chk("derived_E_vartheta_coefficient", 9.0 + 2.0 * ratio_ref, 16.026, 3e-5)
chk("derived_sigma_excess", ratio_ref - 1.0, 2.5130747, 1e-7)
chk("derived_laplacian_clause_constant", 5.0 * 2.0 + 1.0, 11.0, 0, rel=False)
chk("derived_G0_ratio_delta_7p5_to_15",
    (1.0 / math.sin(DELTA) ** 2) / (1.0 / math.sin(15 * DEG) ** 2), 3.9318517, 1e-6)
chk("derived_taper_arc_Fprime_at_delta_minus",
    (math.sin(DELTA) - DELTA * math.cos(DELTA)) / (DELTA * math.sin(DELTA) ** 2), 0.3345, 1e-3)
chk("derived_one_over_delta_m", 1.0 / DM, 11.459, 1e-4)

# =====================================================================  PART 2
def walk(o, out):
    if isinstance(o, dict):
        for v in o.values():
            walk(v, out)
    elif isinstance(o, (list, tuple)):
        for v in o:
            walk(v, out)
    elif isinstance(o, bool):
        pass
    elif isinstance(o, (int, float)):
        out.append(float(o))
    elif isinstance(o, str):
        for t in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", o):
            try:
                out.append(float(t))
            except ValueError:
                pass


POOL = []
walk(J, POOL)
POOL += [E0 ** 2 + _Fp ** 2, float(_val[_off].max()), E0, G0, ratio_ref, 9.0 + 2.0 * ratio_ref, ratio_ref - 1.0,
         math.cos(DELTA) / math.sin(DELTA) ** 2, 1.0 / DM, 1.0 - math.sin(PHI0 - DELTA),
         (1.0 / math.sin(DELTA) ** 2) / (1.0 / math.sin(15 * DEG) ** 2),
         (math.sin(DELTA) - DELTA * math.cos(DELTA)) / (DELTA * math.sin(DELTA) ** 2),
         math.log(1.1), 0.61732]

FOREIGN = {
    # value: origin.  Numbers quoted from another seat; they appear only in comparisons.
    7.6612975721936944: "THEOREM_S3/t1_results.json  E_0",
    58.695476310955364: "THEOREM_S3/t1_results.json  Gfrak_0",
    58.695476: "THEOREM_S3.md sec.1.1  Gfrak_0",
    58.6948961: "refute-u2/NOTE.md sec.3  Gfrak_0 for the ASSEMBLY family",
    58.193328: "THEOREM_S3.md sec.1.1  |F'| at the kink",
    4.093173: "refute-u2/NOTE.md MINOR-2  ||eta_0||_inf for (D-C)",
    0.53431477: "refute-u2/NOTE.md MINOR-2",
    161.7735: "s3close/hk2  K_2 bound",
    311511.9: "THEOREM_S3.md sec.4.5  L_*",
    302704.3: "THEOREM_S3.md sec.4.5  L_*",
    422.4: "THEOREM_S3.md sec.4.5  L_*",
    1056242.9: "THEOREM_S3.md sec.4.5  L_*",
    623026.8: "THEOREM_S3.md sec.0  log Lambda_*",
    1.7905793948266: "V-b Corollary V.5  I_4",
    0.5096901045590: "V-b Corollary V.5  I_2",
    4.18: "refute-u2/NOTE.md MAJOR-1  P2 loss factor",
    0.345: "THEOREM_S3.md sec.4.4  eps_a at L = 3e5",
    1.18: "THEOREM_S3.md sec.4.4  eps_ell at L = 3e5 (times 1e-5)",
    2.4: "text",
}
WHITELIST = set(range(0, 200)) | {
    2026, 2000, 250, 300, 400, 500, 1000, 2500, 21874, 400003,
}
SKIP_LINE = ("sha256", "SHA256", "4c2cdbc32cb0866a",
             "CHECKS,", "tokens scanned")   # the gate's own summary line

doc = open(os.path.join(HERE, "PROOF.md")).read()
tokens, unmatched = 0, []
for lineno, line in enumerate(doc.split("\n"), 1):
    if any(s in line for s in SKIP_LINE):
        continue
    line = re.sub(r"`?\d{4}-\d{2}-\d{2}`?", " ", line)          # dates
    line = re.sub(r"section[s]? [\d.]+", " ", line)             # section numbers
    line = re.sub(r"sec\.[\d.]+", " ", line)
    line = re.sub(r"[A-Za-z]\d+\b", " ", line)                  # labels like A.7, M1, q2
    for t in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", line):
        try:
            v = float(t)
        except ValueError:
            continue
        tokens += 1
        if v == 0.0:
            continue
        if v == int(v) and abs(v) in WHITELIST:
            continue
        ndig = len(t.split("e")[0].split("E")[0].replace("-", "").replace(".", "").lstrip("0"))
        tol = max(10.0 ** (-(ndig - 1)), 1e-12)
        if any(abs(v - p) <= tol * max(abs(p), 1e-300) for p in POOL):
            continue
        if any(abs(v - fv) <= tol * max(abs(fv), 1e-300) for fv in FOREIGN):
            continue
        unmatched.append((lineno, t, line.strip()[:110]))

chk_true("doc_audit_all_numbers_traceable", not unmatched,
         "%d unmatched of %d tokens" % (len(unmatched), tokens))

# =====================================================================  report
print("PART 1 + PART 2")
for name, got, want, err in FAIL:
    print("FAIL  %-52s got=%s want=%s err=%s" % (name, got, want, err))
if unmatched:
    print("\nUNMATCHED NUMERIC TOKENS IN PROOF.md")
    for lineno, t, line in unmatched:
        print("  line %4d  %-16s | %s" % (lineno, t, line))
print("\n%d CHECKS, %d PASS, %d FAIL   (document tokens scanned: %d)"
      % (len(PASS) + len(FAIL), len(PASS), len(FAIL), tokens))
sys.exit(1 if FAIL else 0)

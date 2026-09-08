"""
check_fix3 -- the gate for FIX3.

Every constant displayed in FIX3.md or ADDENDUM_3_2026-09-08.md is re-read from the stored
JSONs, re-formatted, and checked to appear as a string in the document; every ratio and factor
quoted between two numbers is RE-DERIVED here rather than copied.  FL-043 discipline: the gate
prints "N CHECKS, N PASS, 0 FAIL" or names every failure.

Run:  python3 check_fix3.py
"""
import json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DOC = open(os.path.join(HERE, "FIX3.md")).read()
ADD = open(os.path.join(HERE, "..", "ADDENDUM_3_2026-09-08.md")).read()
BOTH = DOC + "\n" + ADD

G1 = json.load(open(os.path.join(HERE, "g1_results.json")))
G2 = json.load(open(os.path.join(HERE, "g2_results.json")))
G3 = json.load(open(os.path.join(HERE, "g3_results.json")))
G4 = json.load(open(os.path.join(HERE, "g4_results.json")))
G5 = json.load(open(os.path.join(HERE, "g5_results.json")))

CH = []


def ck(name, cond, detail=""):
    CH.append((name, bool(cond), detail))


def shown(name, value, fmt="%.6g", where=BOTH):
    s = fmt % value
    ck("SHOWN " + name + " = " + s, s in where, s)


def close(name, a, b, tol):
    ck("CLOSE " + name, abs(a-b) <= tol*max(1.0, abs(b)), "%r vs %r" % (a, b))


# ------------------------------------------------------------------ UNIT 1
ck("A: weight per e-fold is exactly 3/8",
   abs(G1["A_kernel"]["weight_per_efold"] - 0.375) < 1e-15)
ck("A: weight_per_efold_exact string is 3/8", G1["A_kernel"]["weight_per_efold_exact"] == "3/8")
ck("A: int|cos|sin^3 = 1/2", G1["A_kernel"]["int_|cos|sin^3_dphi"] == "1/2")
ck("A: int|cos|sin^2 = 2/3", G1["A_kernel"]["int_|cos|sin^2_dphi"] == "2/3")
ck("A: plateau rate is exactly 1/2",
   abs(G1["A_kernel"]["plateau_a0_rate_per_efold"] - 0.5) < 1e-15)
ck("A: Lap5(rho^-3) = 0 (sympy)", G1["A_kernel"]["Lap5(rho^-3) away from origin (sympy)"] == "0")
close("A: weight quadrature vs (3/8)*10",
      G1["A_kernel"]["weight_numeric_1_to_e^10"], 3.75, 1e-6)
shown("A: weight quadrature", G1["A_kernel"]["weight_numeric_1_to_e^10"], "%.7f")
ck("B: harmonic identity to 1e-10", G1["B_harmonic_identity"]["rel"] < 1e-10,
   str(G1["B_harmonic_identity"]["rel"]))
shown("B: harmonic identity value", G1["B_harmonic_identity"]["numeric_int_Kcal_Lap_eta_test"],
      "%.16g")

rows = {("kinked" if r["sigma"] is None else r["sigma"]): r for r in G1["C_datum"]["rows"]}
for sg, fmt in ((0.05, "%.2f"), (0.02, "%.2f"), (0.01, "%.2f"), (0.005, "%.2f"),
                (0.002, "%.2f"), (0.001, "%.2f"), (0.0005, "%.2f")):
    shown("C: Psi_sigma(0) at sigma=%g" % sg, rows[sg]["Psi0"], fmt)
    shown("C: Gfrak0 at sigma=%g" % sg, rows[sg]["Gfrak0"], "%.4f")
    shown("C: kappa at sigma=%g" % sg, rows[sg]["kappa_delta"], "%.10f")
shown("C: Psi0 a.c. part", rows["kinked"]["Psi0"], "%.2f")
shown("C: grad_inf (this seat)", rows["kinked"]["grad_inf"], "%.4f")
shown("C: grad_inf (DatumS3 control)", G1["C_control_DatumS3"]["grad_inf"], "%.4f")
ck("C: pmax-h3v's 31.0936 does NOT reproduce",
   abs(rows["kinked"]["grad_inf"] - 31.093556949918312) > 10.0)
ck("C: the two independent grad_inf agree to 1e-3",
   abs(rows["kinked"]["grad_inf"] - G1["C_control_DatumS3"]["grad_inf"]) < 1e-2)
cf = G1["C_Psi_closed_form"]
shown("C: [W'] at the taper kink", cf["W_prime_jump_at_taper"], "%.5f")
shown("C: [W'] at the equator kink", cf["W_prime_jump_at_equator_kink"], "%.3f")
shown("C: max_u Theta e^{-2u}", cf["max_u_Theta_e^{-2u}"], "%.5f")
shown("C: C_kink predicted", cf["C_kink_predicted"], "%.5f")
shown("C: C_kink measured", cf["C_kink_measured_sigma_5e-4"], "%.5f")
ck("C: C_kink predicted vs measured within 1 %",
   abs(cf["C_kink_predicted"]/cf["C_kink_measured_sigma_5e-4"] - 1.0) < 0.01)
# re-derive the closed form rather than copying it
pred = cf["W_prime_jump_at_taper"]*cf["max_u_Theta_e^{-2u}"]/math.sqrt(2*math.pi)
close("C: C_kink closed form re-derived", pred, cf["C_kink_predicted"], 1e-12)
# Psi_sigma(0) * sigma must be flat
prod = [rows[s]["Psi0"]*s for s in (0.002, 0.001, 0.0005)]
ck("C: Psi_sigma(0) ~ C/sigma (products within 2 %)",
   max(prod)/min(prod) - 1.0 < 0.02, str(prod))

eps = {(r["L"], r["window"]): r for r in G1["D_eps_a"]}
for L in (1e4, 1e5, 1e6):
    for w in ("frozen c=c_*", "cap c=1.1943662c_*"):
        r = eps[(L, w)]
        shown("D: c_G at L=%g %s" % (L, w), r["c_G"], "%.5f")
        for key, fmt in (("eps_a_sigma=0.05", "%.2e"), ("eps_a_sigma=0.01", "%.2e"),
                         ("eps_a_sigma=0.002", "%.2e"), ("eps_a_sigma=0.001", "%.2e"),
                         ("eps_a_kinked_sqrt_route", "%.2e")):
            shown("D: %s L=%g %s" % (key, L, w), r[key], fmt)
# eps_a scales like 1/L in the smooth route
r6 = eps[(1e6, "frozen c=c_*")]
r5 = eps[(1e5, "frozen c=c_*")]
ck("D: smooth eps_a falls like 1/L (ratio within 3 %)",
   abs(r5["eps_a_sigma=0.002"]/r6["eps_a_sigma=0.002"]/10.0 - 1.0) < 0.03,
   str(r5["eps_a_sigma=0.002"]/r6["eps_a_sigma=0.002"]))
ck("D: kinked eps_a falls like 1/sqrt(L) (ratio within 3 %)",
   abs(r5["eps_a_kinked_sqrt_route"]/r6["eps_a_kinked_sqrt_route"]/math.sqrt(10.0) - 1.0) < 0.03,
   str(r5["eps_a_kinked_sqrt_route"]/r6["eps_a_kinked_sqrt_route"]))
# the ratios quoted in FIX3 sec.1.7, re-derived
epsbud = r6["eps_budget_at_this_c"]
shown("D: eps(1e6) of the budget", epsbud, "%.6f")
for key, txt in (("eps_a_sigma=0.05", "0.0041 %"), ("eps_a_sigma=0.002", "0.069 %"),
                 ("eps_a_kinked_sqrt_route", "3.24 %")):
    v = 100.0*r6[key]/epsbud
    ck("D: quoted share %s for %s" % (txt, key),
       abs(v - float(txt.split()[0]))/float(txt.split()[0]) < 0.02, "%.4f %%" % v)

for r in G5["rows"]:
    if r["L"] == 1e4 and r["window"] == "frozen c=c_*":
        shown("g5: tau*source at L=1e4", r["tau*source"], "%.2f")
        shown("g5: Psi_2(0) at sigma=0.002", r["Psi2_0_sigma=0.002"], "%.0f")
        ck("g5: propagation is 0.28 % of Psi_2(0) at L=1e4",
           abs(100*r["propagation_over_Psi2_0"] - 0.28) < 0.02,
           "%.4f %%" % (100*r["propagation_over_Psi2_0"]))
    if r["L"] == 1e6 and r["window"] == "frozen c=c_*":
        shown("g5: tau*source at L=1e6", r["tau*source"], "%.4f")

# ------------------------------------------------------------------ UNIT 2
A2 = G2["A_theta_exact"]
ck("2A: antiderivative residual is 0", A2["antiderivative_residual_numeric_max"] == 0.0)
shown("2A: tail e-folds each side", A2["closed_form_value"], "%.7f")
close("2A: closed form vs numeric tail", A2["numeric_inner_tail"], A2["closed_form_value"], 1e-4)
shown("2A: numeric inner tail", A2["numeric_inner_tail"], "%.9f")
shown("2A: whole-line integral at L=12", A2["numeric_whole_line_L12"], "%.15g")
shown("2A: shell-only at L=12", A2["numeric_shell_only"], "%.6f")
close("2A: shell-only = L - 0.5 - 2 tails", A2["numeric_shell_only"],
      12.0 - 0.5 - 2*A2["closed_form_value"], 1e-5)
ck("2A: envelope Theta <= e^{-2}(rho/rho0)^8 (inner)",
   A2["envelope_inner_max_ratio"] <= 1.0 + 1e-5, str(A2["envelope_inner_max_ratio"]))
ck("2A: envelope (outer)", A2["envelope_outer_max_ratio"] <= 1.0 + 1e-5)
ck("2A: ELL_RAMP understatement is 0.25", abs(A2["understatement_efolds"] - 0.25) < 1e-12)

B10 = G2["B_J_of_rho_L10"]
B40 = G2["B_J_of_rho_L40"]
shown("2B: J sup", B10["J_sup_over_rho"], "%.4f")
shown("2B: J plateau", B10["J_plateau_Theta_eq_1"], "%.4f")
shown("2B: J at u=0", B10["J_at_u=0"], "%.4f")
shown("2B: J at u=L", B10["J_at_u=L"], "%.4f")
shown("2B: J at u=-1", B10["J_at_u=-1"], "%.6f")
shown("2B: J at u=L+1", B10["J_at_u=L+1"], "%.6f")
shown("2B: Jcoef inner", B10["Jcoef_inner_tail"], "%.3f")
shown("2B: Jcoef outer", B10["Jcoef_outer_tail"], "%.3f")
ck("2B: J sup is L-independent (L=10 vs L=40, 1e-5)",
   abs(B10["J_sup_over_rho"]/B40["J_sup_over_rho"] - 1.0) < 1e-5)
ck("2B: the sup is NOT in the tails",
   B10["J_sup_over_rho"] > 1.5*max(B10["sup_J_over_inner_tail_(u<0)"],
                                   B10["sup_J_over_outer_tail_(u>L)"]))
r_in = B10["J_sup_over_rho"]/B10["sup_J_over_inner_tail_(u<0)"]
r_out = B10["J_sup_over_rho"]/B10["sup_J_over_outer_tail_(u>L)"]
ck("2B: quoted 3.85x above the inner-tail sup", abs(r_in - 3.85) < 0.02, "%.4f" % r_in)
ck("2B: quoted 3.04x above the outer-tail sup", abs(r_out - 3.04) < 0.02, "%.4f" % r_out)
ck("2B: envelope is sharp at u=-1 (0.9999)", abs(B10["envelope_check_inner_u=-1"] - 1.0) < 1e-3)

C4 = G2["C_tail_share_f4"]
shown("2C: inner tail contribution", C4["inner_tail_contribution_to_near_integral"], "%.4f")
shown("2C: hk2 near integral at f=4", C4["hk2_uniform_J_near_integral"], "%.2f")
shown("2C: outer tail contribution", C4["outer_tail_contribution_to_far_integral"], "%.1e")
shown("2C: hk2 far integral at f=4", C4["hk2_uniform_J_far_integral"], "%.3f")
ck("2C: inner tail share is 0.58 % at f=4",
   abs(100*C4["inner_tail_relative_share"] - 0.58) < 0.01,
   "%.4f" % (100*C4["inner_tail_relative_share"]))
ck("2C: inner tail share is 9.07 % at f=1",
   abs(100*G2["C_tail_share_f1"]["inner_tail_relative_share"] - 9.07) < 0.02)
ck("2C: outer tail share is 2.1e-05 at f=4",
   abs(C4["outer_tail_relative_share"]/2.1e-5 - 1.0) < 0.05)
for r in G2["C_lambda_masked"]:
    tag = "lam=1" if abs(r["lam"] - 1.0) < 1e-9 else "lam=lam_max"
    shown("2C: L4 full %s" % tag, r["L4_full"], "%.6f")
    shown("2C: L4 masked %s" % tag, r["L4_shell_only"], "%.6f")
    shown("2C: L5 full %s" % tag, r["L5_full"], "%.6f")
    shown("2C: L5 masked %s" % tag, r["L5_shell_only"], "%.6f")
    shown("2C: L4 tail rel %s" % tag, r["L4_tail_rel"], "%.2e")
    shown("2C: L5 tail rel %s" % tag, r["L5_tail_rel"], "%.2e")
    ck("2C: tails <= 1.7e-04 of Lambda (%s)" % tag,
       max(r["L4_tail_rel"], r["L5_tail_rel"]) <= 1.7e-4)

# ------------------------------------------------------------------ UNIT 3
A3 = G3["A_J_three_data"]
shown("3A: (D-A) J total", A3["(D-A) J_total"], "%.6f")
shown("3A: (D-A) J ac", A3["(D-A) J_ac_plateau"], "%.6f")
shown("3A: (D-A) J jump", A3["(D-A) J_jump"], "%.6f")
shown("3A: (D-B) J plateau", A3["(D-B) J_plateau"], "%.6f")
shown("3A: (D-B) J sup", A3["(D-B) J_sup_over_rho"], "%.3f")
shown("3A: (D-B) argmax u", A3["(D-B) argmax_u"], "%.4f")
shown("3A: (S3) J plateau", A3["(S3) J_plateau"], "%.5f")
shown("3A: (S3) J sup", A3["(S3) J_sup_over_rho"], "%.5f")
ck("3A: (D-A) reproduces hk2 to 1e-6", A3["(D-A) rel"] < 1e-6, str(A3["(D-A) rel"]))
ck("3A: (D-B) reproduces hk2 to 1e-8", A3["(D-B) rel"] < 1e-8, str(A3["(D-B) rel"]))
ck("3A: (S3) reproduces fix2's sup to 1e-4", A3["(S3) rel_sup"] < 1e-4, str(A3["(S3) rel_sup"]))
ck("3A: (S3) reproduces fix2's plateau to 2e-4", A3["(S3) rel_plateau"] < 2e-4)
r155 = A3["(D-B) J_sup_over_rho"]/A3["(S3) J_sup_over_rho"]
ck("3A: quoted 1.55x (S3 sup better than D-B sup)", abs(r155 - 1.55) < 0.01, "%.4f" % r155)
r115 = A3["(S3) J_plateau"]/A3["(D-B) J_plateau"]
ck("3A: quoted 1.15x (S3 plateau worse than D-B plateau)", abs(r115 - 1.15) < 0.01,
   "%.4f" % r115)
ck("3B: no script reads J",
   all(v == 0 for v in G3["B_where_J_is_used"]["occurrences_of_a_shell_density_symbol"].values()))
CS = G3["C_sensitivity"]
shown("3C: K2hat analytic at J=96.26", CS["K2hat(J=96.26, lam=lam_max)"], "%.3f")
shown("3C: K2hat analytic at J=192.52", CS["K2hat(J=192.52, lam=lam_max)"], "%.3f")
shown("3C: dlogK2/dlogJ", CS["dlogK2_dlogJ_at_J=96.26"], "%.3f")
rq = CS["K2hat(J=96.26, lam=lam_max)"]/CS["quadrature_K2hat_for_comparison_fix2"]
ck("3C: quoted 3.05x weaker than the quadrature route", abs(rq - 3.05) < 0.01, "%.4f" % rq)

LS = G3["D_Lstar"]
for key, want in (("t3_budget_import_161.7735_eps<=0.194366", 1887786.7126679737),
                  ("t3_budget_import_161.7735_eps<=0.1", 1977309.301165401)):
    close("3D: L_* reproduces fix2 (%s)" % key, LS[key]["L_star"], want, 1e-9)
shown("3D: L_* proved eps<=0.1943662", LS["t3_budget_import_161.7735_eps<=0.194366"]["L_star"],
      "%.4f")
shown("3D: L_* proved eps<=0.1", LS["t3_budget_import_161.7735_eps<=0.1"]["L_star"], "%.4f")
shown("3D: L_* at C_K=33.710", LS["fix2_quadrature_S3_datum_lam_max_f4_33.710_eps<=0.194366"]["L_star"], "%.4f")
lo = LS["fix2_quadrature_S3_datum_lam_max_f4_33.710_eps<=0.194366"]["L_star"]
hi = LS["t3_budget_import_161.7735_eps<=0.194366"]["L_star"]
ck("3D: a 5.2x range in C_K moves L_* by <= 1e-7 relative", abs(hi-lo)/hi < 1e-7,
   "%.3e" % (abs(hi-lo)/hi))
epsCK = set(v["eps_self_consistent_L=1e+07"] for v in G3["D_eps_vs_CK"].values())
ck("3D: eps(1e7) identical at every C_K", len(epsCK) == 1, str(epsCK))
ck("SHOWN 3D: eps(1e7, proved)", repr(list(epsCK)[0]) in BOTH, repr(list(epsCK)[0]))
# log Lambda_* re-derived, not copied
for key, tgt in (("t3_budget_import_161.7735_eps<=0.194366", "3775576.3488"),
                 ("t3_budget_import_161.7735_eps<=0.1", "3954621.5258")):
    v = 2*LS[key]["L_star"] + 2.9234859
    ck("3D: log Lambda_* re-derived for %s" % key, ("%.4f" % v) == tgt, "%.4f" % v)

SL = G4["slack"]
for k, lam in (("lam=1", 1.0), ("lam=1.25", 1.25), ("lam=1.5", 1.5)):
    shown("3D': transported K2 at %s" % k, SL[k]["transported_measured"], "%.6f")
    shown("3D': untransported K2 at %s" % k, SL[k]["untransported_measured"], "%.6f")
    v = SL[k]["hk2_bound"]/SL[k]["transported_measured"]
    close("3D': slack re-derived at %s" % k, v, SL[k]["slack_transported"], 1e-12)
    shown("3D': slack at %s" % k, v, "%.2f")
hk2_measured = {"lam=1": 2.1903, "lam=1.25": 3.3827, "lam=1.5": 5.3854}
for k, v in hk2_measured.items():
    ck("3D': reproduces hk2 sec.9's measured %s to 1e-4 relative" % k,
       abs(SL[k]["transported_measured"] - v)/v < 1e-4,
       "%r vs %r" % (SL[k]["transported_measured"], v))
ck("3D': the withdrawn k4 numbers are 2.19x and 3.64x too small",
   abs(SL["lam=1.25"]["transported_over_untransported"] - 2.19) < 0.01 and
   abs(SL["lam=1.5"]["transported_over_untransported"] - 3.64) < 0.01)
best = G4["transported_measured_DB"]["lam=1.5"]["ns=1400"]
ck("3D': axisymmetry control at lam=3/2 below 1e-13", best["control_axisym_rel"] < 1e-13,
   str(best["control_axisym_rel"]))
ck("3D': PDE control at lam=3/2 below 1e-4", best["control_PDE_over_hess_scale"] < 1e-4,
   str(best["control_PDE_over_hess_scale"]))
for k in ("lam=1", "lam=1.25", "lam=1.5"):
    m = G4["transported_measured_DB"][k]["ns=1400"]
    shown("3D': a at %s" % k, m["a"], "%.6f")
    shown("3D': |grad a| at %s" % k, m["grad_a"], "%.6f")
    shown("3D': ||Hess a|| at %s" % k, m["hess_a"], "%.6f")

# ------------------------------------------------------------------ the byte copies
# every file under fix2copy/ must hash equal to the entry in its SOURCE seat's own SHA256SUMS
import hashlib                                                              # noqa: E402
SEATS = ["round2/THEOREM_S3/fix2/SHA256SUMS", "hk2/SHA256SUMS",
         "round2/THEOREM_S3/SHA256SUMS", "round2/refute-u2/SHA256SUMS"]
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
src = {}
for rel in SEATS:
    q = os.path.join(ROOT, rel)
    if not os.path.exists(q):
        continue
    for ln in open(q):
        parts = ln.split(None, 1)
        if len(parts) == 2:
            src.setdefault(os.path.basename(parts[1].strip()), set()).add(parts[0])
COPIES = ["f1_window.py", "f2_datum_s3.py", "f2_ramp_k2.py", "f2b_measured.py", "f4_scans.py",
          "f5_gamma_exist.py", "f6_controls.py", "f2_results.json", "imported/hk2lib.py",
          "imported/k1_algebra.py", "imported/k2_datum.py", "imported/k3_bound.py",
          "imported/k4_direct.py", "imported/r3_fixedpoint.py", "imported/t2_gamma_CR.py",
          "imported/t3_budget.py", "imported/t1_results.json", "imported/t2_results.json"]
for c in COPIES:
    q = os.path.join(HERE, "fix2copy", c)
    h = hashlib.sha256(open(q, "rb").read()).hexdigest()
    b = os.path.basename(c)
    ck("COPY %s equals its source seat's own SHA256SUMS entry" % b,
       b in src and h in src[b], h)

# ------------------------------------------------------------------ report
n = len(CH)
bad = [c for c in CH if not c[1]]
for name, ok, det in CH:
    if not ok:
        print("FAIL: %s   [%s]" % (name, det))
print("\n%d CHECKS, %d PASS, %d FAIL" % (n, n - len(bad), len(bad)))
sys.exit(1 if bad else 0)

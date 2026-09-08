"""
p3 -- WHAT DOES AND DOES NOT DEPEND ON THE ANGULAR PROFILE:
      the exact kernel constants, the far/near constants, the shell density J, K_2 for the
      mollified datum, and the ELL_RAMP e-fold count.

A. EXACT KERNEL CONSTANTS.  C_K |S^4| = 1, the Riesz composition pi^4/2 in R^5,
   3 pi^2/16 = C_K * pi^4/2, 2 R_A.  These are properties of K(w) = -d_z G_5(w) and of R^5,
   not of the datum; recomputed here and cross-checked by an independent quadrature (L-98).

B. THE FAR/NEAR CONSTANTS (rebuild/far-near-kernel-lemma sec.2-3).  C1 = 0.291999 (far),
   C2in = 0.014754 (inner), C_collar = 3.999218/sin phi + pi/8.  Their hypotheses are
   (i) |omega^theta| <= M pointwise on the relevant region, (ii) omega^theta odd in z (only
   odd l survive), (iii) the Parseval bound  sum_l g_l^2 N_l = int_{-1}^1 w(t)^2 dt <= 2 M^2.
   All three are checked here for the MOLLIFIED, amplitude-normalised profile.  None of them
   mentions the taper, so the constants are unchanged; this file proves the hypotheses rather
   than asserting the conclusion.

C. THE SHELL DENSITY J and K_2.  fix2's f2_ramp_k2 machinery (hk2's k3_bound.lambda_bulk /
   assemble / sups_on_ball_B, UNCHANGED) is run on the mollified datum, obtained by
   subclassing fix2's DatumS3 and overriding W(phi, k) by the mollified profile of p1.
   The question is only whether the budget's imported majorant K_2 <= 161.7735 M/rho_0 still
   dominates over lambda in [1, lam_max].

D. ELL_RAMP.  int Theta du = L - 2 eps_r = L - 0.5 EXACTLY (fix3 sec.2.2(a)), so the tanh
   ramp costs 0.5 e-folds and t3_budget.py's ELL_RAMP = 0.25 undercharges by 0.25/L.

Outputs -> p3_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "imported"))
import p1_profile as P1                                    # noqa: E402
import t2_gamma_CR as T2                                   # noqa: E402
import k3_bound as K3                                      # noqa: E402  (hk2, unchanged)
from f2_datum_s3 import DatumS3, StrainedS3, geometry_s3   # noqa: E402  (fix2, unchanged)

DEG = math.pi/180.0
DELTA, DM, EPS_R = 7.5*DEG, 5.0*DEG, 0.25
S3SPH = 2*math.pi**2
SIGMA = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
OUT = {"sigma": SIGMA}

# ---------------------------------------------------------------- A. kernel constants
OUT["A_kernel_constants"] = {
    "C_K_3_over_8pi2": T2.C_K, "S4_8pi2_over_3": T2.S4, "C_K_times_S4": T2.CK_S4,
    "riesz_pi4_over_2": T2.RIESZ, "riesz_independent_quadrature": T2.riesz_check(),
    "coef_3pi2_over_16": T2.COEF_G, "R_A": T2.R_A, "two_R_A": T2.TWO_RA,
    "note": "kernel and dimension only; no dependence on the angular profile"}
OUT["A_kernel_constants"]["riesz_relerr"] = abs(
    OUT["A_kernel_constants"]["riesz_independent_quadrature"] - T2.RIESZ)/T2.RIESZ

# ---------------------------------------------------------------- B. far/near hypotheses
prof = P1.Profile(SIGMA)
phi = prof.phi
a0 = prof.a0                                   # A_sigma/N_sigma, the angular amplitude
OUT["B_far_near_hypotheses"] = {}
# (i) sup |omega^theta|/M
OUT["B_far_near_hypotheses"]["sup_abs_angular_amplitude"] = float(np.abs(a0).max())
# (ii) z-oddness: A_sigma(pi - phi) = -A_sigma(phi)
ref = np.interp(math.pi - phi, phi, a0)
OUT["B_far_near_hypotheses"]["z_oddness_max_abs_defect"] = float(np.max(np.abs(a0 + ref)))
# (iii) Parseval:  int_{-1}^{1} w(t)^2 dt  with t = cos phi
t = np.cos(phi)
order = np.argsort(t)
OUT["B_far_near_hypotheses"]["int_w2_dt"] = float(np.trapz((a0**2)[order], t[order]))
OUT["B_far_near_hypotheses"]["parseval_cap_2M2"] = 2.0
OUT["B_far_near_hypotheses"]["kinked_int_w2_dt"] = None
pk = P1.Profile(None)
tk = np.cos(pk.phi); ok = np.argsort(tk)
OUT["B_far_near_hypotheses"]["kinked_int_w2_dt"] = float(np.trapz((pk.a0**2)[ok], tk[ok]))
OUT["B_far_near_hypotheses"]["constants_unchanged"] = {
    "C_far": T2.C_FAR, "C_inner": T2.C_INNER, "C_collar_over_sin_phi": 3.999218,
    "C_collar_axis_pi_over_8": math.pi/8.0}

# ---------------------------------------------------------------- C. J and K_2
class DatumMoll(DatumS3):
    """fix2's DatumS3 with the angular factor replaced by p1's mollified, normalised W."""
    def __init__(self, prof, **kw):
        DatumS3.__init__(self, **kw)
        self.p = prof

    def W(self, phi, k=0):
        ph = np.abs(np.asarray(phi, dtype=float))
        arr = (self.p.w0, self.p.w1, self.p.w2)[k]
        return np.interp(ph, self.p.phi, arr)

def J_of(D, nphi=200001, nu=2001):
    ph = np.linspace(1e-9, math.pi-1e-9, nphi)
    W0 = D.W(ph, 0)
    W1 = D.W(ph, 1)
    s3 = np.sin(ph)**3
    us = np.linspace(-2.0, D.L + 2.0, nu)
    best, arg = 0.0, None
    for uu in us:
        Th = D.Th_u(uu, 0); Tu = D.Th_u(uu, 1)
        val = S3SPH*float(np.trapz(np.sqrt((Tu-Th)**2*W0**2 + Th**2*W1**2)*s3, ph))
        if val > best:
            best, arg = val, float(uu)
    plateau = S3SPH*float(np.trapz(np.sqrt(W0**2 + W1**2)*s3, ph))
    return {"J_sup_over_rho": best, "argmax_u": arg, "J_plateau": plateau}

def k2_bound(D0, lam, f, nd=8, NS=250, NT=100, NC=100, frac=0.95):
    DL = StrainedS3(D0, lam)
    g = geometry_s3(30.0, f, D0.delta_deg, D0.dm_deg, lam)
    dmax = frac*min(g["d_layers"], 0.95*g["r"])
    best = None
    for d in np.linspace(0.05*dmax, dmax, nd):
        Ls = K3.lambda_bulk(DL.grad_norm_fast, (g["r"], g["z"]), d, 3*D0.R, NS, NT, NC)
        sg, sh, sb = K3.sups_on_ball_B(DL, g["r"], g["z"], d)
        ee = [float(v[0]) for v in DL.eta_rz(np.array([g["r"]]), np.array([g["z"]]))]
        A = K3.assemble(sg, sh, sb, d, Ls[4], Ls[5], g["r"], ee[0], ee[1], ee[2])
        A.update(d=d, L4=Ls[4], L5=Ls[5], sup_grad=sg, sup_hess=sh, sup_grad_bdry=sb)
        if best is None or A["K2"] < best["K2"]:
            best = A
    best.update(lam=lam, f=f, r=g["r"], z=g["z"], d_layers=g["d_layers"], L=D0.L)
    return best

if __name__ == "__main__":
    P1RES = json.load(open(os.path.join(HERE, "p1_results.json")))
    row = P1RES["rows"]["sigma=%g" % SIGMA]
    cs = row["c_star"]
    LAMS = [1.0, 1.5, math.exp(0.75*cs), math.exp(0.75*(1.0+row["eps_cap_certified"])*cs)]
    D0 = DatumMoll(prof, delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=10.0)
    Dk = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=10.0)
    OUT["C_J"] = {"mollified": J_of(D0), "kinked": J_of(Dk),
                  "fix2_kinked_sup": 96.249844, "fix2_kinked_plateau": 75.251104,
                  "hk2_DA": 78.641147, "hk2_DB": 65.625910}
    print("J mollified:", OUT["C_J"]["mollified"], flush=True)
    print("J kinked   :", OUT["C_J"]["kinked"], flush=True)
    rows = []
    for lam in LAMS:
        r = k2_bound(D0, lam, 4.0)
        rows.append(r)
        print("K2hat(mollified, f=4, lam=%.6f) = %.4f   (d*=%.4f)" % (lam, r["K2"], r["d"]),
              flush=True)
    OUT["C_K2_mollified_f4_L10"] = rows
    # control: the same call on fix2's kinked datum must reproduce fix2's table
    ctl = []
    for lam in [1.0, 1.5, 1.8420112, 2.0741]:
        r = k2_bound(Dk, lam, 4.0)
        ctl.append({"lam": lam, "K2": r["K2"]})
        print("  control kinked f=4 lam=%.6f: K2hat = %.4f" % (lam, r["K2"]), flush=True)
    OUT["C_K2_kinked_control"] = {"rows": ctl,
                                  "fix2_quoted_f4": {"1.0": 19.169, "1.5": 25.418,
                                                     "1.8420112": 33.710, "2.0741": 42.008}}
    OUT["C_K2_majorant"] = {
        "imported_majorant_161.7735": 161.7735,
        "max_over_lam_mollified": max(r["K2"] for r in rows),
        "slack": 161.7735/max(r["K2"] for r in rows)}

    # ---------------------------------------------------------------- D. ELL_RAMP
    er = EPS_R
    for L in (12.0, 40.0):
        u = np.linspace(-60.0, L + 60.0, 4000001)
        Th = 0.5*(np.tanh((u-er)/er) - np.tanh((u-L+er)/er))
        OUT.setdefault("D_ell_ramp", {})["int_Theta_du_L=%g" % L] = float(np.trapz(Th, u))
        OUT["D_ell_ramp"]["L_minus_2eps_r_L=%g" % L] = L - 2*er
    OUT["D_ell_ramp"]["ramp_cost_e_folds"] = 2*er
    OUT["D_ell_ramp"]["t3_budget_ELL_RAMP"] = 0.25
    OUT["D_ell_ramp"]["correction"] = 2*er - 0.25
    OUT["D_ell_ramp"]["inner_tail_exact_log1pe2_over_8"] = math.log(1+math.exp(-2))/8.0

    json.dump(OUT, open(os.path.join(HERE, "p3_results.json"), "w"), indent=1,
              sort_keys=True, default=str)
    print(json.dumps({k: v for k, v in OUT.items() if k != "C_K2_mollified_f4_L10"},
                     indent=1, default=str))

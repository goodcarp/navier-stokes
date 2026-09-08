"""
u2 -- every constant the slaved-reference argument needs, computed here.

RECORD holds the constants QUOTED from other seats (statements of refereed theorems).
Everything else is computed in this file.

  P_h(lam) = 3 int_0^1 h(arcsin v) v^2 (A v^2 + B)^{-5/2} dv ,  A = lam^2 - lam^-4, B = lam^-4
             (lower/prove-lagrangian sec.2 (2.1); re-derived weight, not imported)
  kappa_h  = P_h(1)/2 ,  r_h(lam_max) = inf_{[1,lam_max]} P_h(lam)/lam

  lam_max  = exp(3c/4)   :  frak_a <= (M_s/2) L <= (3/4) M L  and  tau = c/(ML)
  kappa_s  = (3/4) c     :  |d log lam / d log rho| = |int_0^s d frak_a/d log rho| <= (3/4) c / L
"""
import json, math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import mpmath as mp

DEG = math.pi/180.0

RECORD = {
    # rebuild/far-near-kernel-lemma/NOTE.md sec.2, sec.3 (z-odd constants; PROVED there)
    "C_far_zodd":            0.291999,     # Prop 1, z-odd
    "C_inner_zodd":          0.014754,     # Prop 2(a), z-odd
    "C_collar_over_sinphi":  3.999218,     # Prop 2(b)
    "C_collar_axis":         math.pi/8.0,  # Prop 2(b)
    # the same lemma's GRADIENT constants, quoted in ASSEMBLY sec.2.6 BLOCK 3 (computed by
    # the L3v seat): far and inner parts of r|grad a|.  The COLLAR gradient is NOT proved
    # anywhere -- it is carried here as the free constant G_collar.
    "G_far":                 1.6049285,
    "G_inner":               0.1324254,
    # lower/SYNTHESIS.md sec.1.3 (refuter-lagrangian-b): MEASURED material-point strain offset
    "C_a_measured":          0.069,
    # lower/prove-lagrangian sec.4(2) / s5: MEASURED sup |u^z + 2 A z|/(M rho) for the datum
    "Cz_measured_t0":        0.19143,
    # write/lemma-T-shell-dependent/PROOF.md: Lemma T' relative constant and its measured
    # conservatism factor (the latter is a MEASUREMENT of that seat, not a theorem)
    "lemmaT_rel_lead":       15.0*math.pi/4.0,
    "lemmaT_sens_proved":    0.34913942333376546,
    "lemmaT_sens_measured":  2.1181705106873974,
    # write/L3v-and-gamma-bound + refute-*: Theorem Gamma's O(M) constant
    "Cpp_proved":            151.15,   # Cor A1 + proved V
    "Cpp_computedV":         119.33,
    "Cpp_oddonly":            76.11,
    "Cpp_measured":           11.74,
    # s3close/hk2/PROOF.md sec.0
    "K2_proved_window":      161.7735,
    "K2_proved_lam1":         66.6622,
    "K2_measured_window":      5.3854,
    "K2_measured_global":      3.0202,
    # cross-check targets (record values this file must reproduce with its own instrument)
    "chk_kappa_delta_7p5":   0.4997212305210886,
    "chk_r_h_1to1p5_7p5":    0.981822,
    "chk_Phi_h_3over2_7p5":  1.4735550,
    "chk_logReE_shift":      3.3643455,
}

# ------------------------------------------------------------------ P_h and friends
def h_taper(v, delta):
    """h_delta(phi) = min(1, phi_ax/delta) at phi = arcsin v, v = sin phi in [0,1]."""
    return min(1.0, math.asin(min(1.0, v))/delta)

def h_eq(v, dm):
    """min(1, |phi - pi/2|/delta_m) = min(1, arccos(v)/delta_m)."""
    return min(1.0, math.acos(min(1.0, v))/dm)

def P_h(lam, delta, dm=None, tol=1e-13):
    A = lam**2 - lam**-4
    B = lam**-4
    def f(v):
        w = h_taper(v, delta)
        if dm is not None:
            w *= h_eq(v, dm)
        return w*v*v*(A*v*v + B)**-2.5
    # split at the two kinks: v = sin(delta) and v = cos(dm)
    pts = [0.0, math.sin(delta)] + ([math.cos(dm)] if dm is not None else []) + [1.0]
    pts = sorted(set(p for p in pts if 0.0 <= p <= 1.0))
    tot = 0.0
    for lo, hi in zip(pts[:-1], pts[1:]):
        if hi > lo:
            tot += quad(f, lo, hi, epsabs=tol, epsrel=tol, limit=400)[0]
    return 3.0*tot

def P_1(lam):
    """h == 1: the exact value is lam (prove-lagrangian Lemma 1); computed here to check."""
    A = lam**2 - lam**-4; B = lam**-4
    return 3.0*quad(lambda v: v*v*(A*v*v + B)**-2.5, 0.0, 1.0,
                    epsabs=1e-14, epsrel=1e-14, limit=400)[0]

def r_h(lam_max, delta, dm=None, n=241):
    lams = np.linspace(1.0, lam_max, n)
    return float(min(P_h(l, delta, dm)/l for l in lams))

OUT = {"RECORD": RECORD}

DELTA = 7.5*DEG
DM = 5.0*DEG
PHI0 = 30.0*DEG

# ---- exactness check of Lemma 1's P_1(lam) = lam ------------------------------------
OUT["P1_minus_lam"] = {str(x): P_1(x) - x for x in [1.0, 1.25, 1.5, 1.8377, 2.0, 3.0]}

# ---- kappa and r_h ------------------------------------------------------------------
kappa_delta = 0.5*P_h(1.0, DELTA)
kappa_full = 0.5*P_h(1.0, DELTA, DM)
OUT["kappa_delta_7.5"] = kappa_delta
OUT["kappa_delta_7.5_minus_record"] = kappa_delta - RECORD["chk_kappa_delta_7p5"]
OUT["kappa_full_7.5_5.0"] = kappa_full
OUT["kappa_delta_table"] = {str(d): 0.5*P_h(1.0, d*DEG) for d in [3, 5, 7.5, 10, 15, 20, 30]}

c_star = math.log(1.5)/kappa_delta
OUT["c_star"] = c_star
lam_max = math.exp(0.75*c_star)
OUT["lam_max_apriori"] = lam_max          # exp(3 c/4)
OUT["kappa_s_slaved"] = 0.75*c_star       # sup L |rho lam'/lam|

OUT["r_h_1_to_1.5"] = r_h(1.5, DELTA)
OUT["r_h_1_to_1.5_minus_record"] = r_h(1.5, DELTA) - RECORD["chk_r_h_1to1p5_7p5"]
OUT["r_h_1_to_lam_max"] = r_h(lam_max, DELTA)
OUT["r_h_1_to_lam_max_full"] = r_h(lam_max, DELTA, DM)
OUT["r_h_1_to_1.5_full"] = r_h(1.5, DELTA, DM)
OUT["Phi_h_3over2_7.5"] = P_h(1.5, DELTA)/P_h(1.0, DELTA)
OUT["Phi_h_3over2_minus_record"] = P_h(1.5, DELTA)/P_h(1.0, DELTA) - RECORD["chk_Phi_h_3over2_7p5"]

# ---- monotonicity of P_h (needed for the comparison principle) ----------------------
def mono(lo, hi, delta, dm=None, n=401):
    ls = np.linspace(lo, hi, n)
    vs = np.array([P_h(l, delta, dm) for l in ls])
    d = np.diff(vs)
    return {"increasing": bool((d > 0).all()), "min_increment": float(d.min()),
            "P_lo": float(vs[0]), "P_hi": float(vs[-1])}
OUT["P_h_monotone_1_to_1.5"] = mono(1.0, 1.5, DELTA)
OUT["P_h_monotone_1_to_lam_max"] = mono(1.0, lam_max, DELTA)
OUT["P_h_full_monotone_1_to_lam_max"] = mono(1.0, lam_max, DELTA, DM)

# ---- C_a (far/near kernel lemma, z-odd, at polar angle phi) -------------------------
def C_a(sinphi):
    return (RECORD["C_far_zodd"] + RECORD["C_inner_zodd"]
            + RECORD["C_collar_over_sinphi"]/sinphi + RECORD["C_collar_axis"])
OUT["C_a_phi0_30deg"] = C_a(math.sin(PHI0))
OUT["C_a_table"] = {str(d): C_a(math.sin(d*DEG)) for d in [10, 30, 45, 62.833, 74.4, 90]}

# ---- C_R: the FULL velocity remainder |b(X) - frak_a(|X|) D X| <= C_R M_s |X| -------
# |R_y| <= C_a M r <= C_a M |X|
# |R_z| <= [2 C_a + log(1/sin phi) + G_1] M |z| ,  G_1 = 1/2 + G_far + G_inner + G_collar
def C_R(sinphi, G_collar):
    G1 = 0.5 + RECORD["G_far"] + RECORD["G_inner"] + G_collar
    return 3.0*C_a(sinphi) + math.log(1.0/sinphi) + G1
def C_R_measured(G_collar):
    """the same object with the MEASURED value constants (R9's 0.069 for C_a)."""
    G1 = 0.5 + RECORD["G_far"] + RECORD["G_inner"] + G_collar
    return 3.0*RECORD["C_a_measured"] + math.log(2.0) + G1
OUT["G_1_proved_parts"] = 0.5 + RECORD["G_far"] + RECORD["G_inner"]
OUT["C_R_proved_Gc0"] = C_R(math.sin(PHI0), 0.0)
OUT["C_R_proved_Gc_value"] = C_R(math.sin(PHI0), C_a(math.sin(PHI0)))
OUT["C_R_proved_Gc_5x"] = C_R(math.sin(PHI0), 5.0*C_a(math.sin(PHI0)))
OUT["C_R_measured_Gc0"] = C_R_measured(0.0)
OUT["C_R_assembly_convention"] = C_a(math.sin(PHI0))   # what ASSEMBLY (2.1) actually charges

# ---- lost e-folds -------------------------------------------------------------------
def ell_loss(f, lam_mx, mu):
    """label cut so that the images of the outer sub-shell sit outside 2|X(s)|:
       (1-mu) lam^-2 rho'' >= 2 (1+mu) lam rho_* ,  rho_* = (1+f) rho_0 ."""
    return math.log(2.0) + math.log(1.0+f) + 3.0*math.log(lam_mx) + math.log((1.0+mu)/(1.0-mu))
def ell_shell(lam_mx, mu):
    """the same for a general shell rho: (1-mu) lam^-2 rho'' >= 2 rho."""
    return math.log(2.0) + 2.0*math.log(lam_mx) + math.log(1.0/(1.0-mu))
OUT["ell_loss_f1_mu0"] = ell_loss(1.0, lam_max, 0.0)
OUT["ell_loss_f1_mu0_lam1.5"] = ell_loss(1.0, 1.5, 0.0)
OUT["ell_shell_mu0"] = ell_shell(lam_max, 0.0)
OUT["ell_loss_assembly"] = math.log(2.0) + math.log(2.0) + 3.0*math.log(1.5)   # ASSEMBLY sec.2.4, f=1

# ---- trajectory margins -------------------------------------------------------------
def phi_of_lambda(phi0, lam):
    return math.atan2(lam**3*math.sin(phi0), math.cos(phi0))
OUT["phi_at_lam_1.5_deg"] = phi_of_lambda(PHI0, 1.5)/DEG
OUT["phi_at_lam_max_deg"] = phi_of_lambda(PHI0, lam_max)/DEG
OUT["equator_margin_lam_max_deg"] = 90.0 - phi_of_lambda(PHI0, lam_max)/DEG
OUT["equator_margin_minus_dm_deg"] = 90.0 - phi_of_lambda(PHI0, lam_max)/DEG - 5.0

# ---- Lemma T' sensitivity -----------------------------------------------------------
OUT["lemmaT_conservatism_factor"] = RECORD["lemmaT_sens_measured"]/RECORD["lemmaT_sens_proved"]

# ---- the energy/Reynolds dictionary shift (independent re-derivation of 3.3643455) ---
# log Re_E = 2L + 2 log s + (2/5) log C_E ,  s = 1/sin(delta) .
# C_E for the bang-bang cap is a record number; here only the s-part is recomputed and the
# record shift is checked for internal consistency.
S_REC = 1.0/math.sin(DELTA)
OUT["s_record"] = S_REC
OUT["C_E_implied_by_record_shift"] = math.exp((RECORD["chk_logReE_shift"] - 2*math.log(S_REC))*2.5)

with open("u2_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
for k in sorted(OUT):
    if k != "RECORD":
        print("%-34s %s" % (k, OUT[k]))

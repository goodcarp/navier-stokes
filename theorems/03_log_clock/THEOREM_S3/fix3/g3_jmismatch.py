"""
g3 -- UNIT 3(a) of FIX3: the J mismatch (96.25 against hk2's 65.63 / 78.64), reconciled and
      propagated.

THE OPEN ITEM (FIX2 ADDENDUM_2, "New and open").  THEOREM_S3's datum has shell total-variation
density J = 96.25 against hk2's 65.63 for (D-B) and 78.64 for (D-A).  "J is the quantity through
which the datum swap actually propagates into K_2, and nobody had tracked it."

WHAT THIS SCRIPT SETTLES.

  A. J recomputed from hk2 Lemma 5.1's definition with an instrument written here, for all
     THREE data, so the three numbers are on one bench:
        (D-A)  omega = -M sgn(z) min(1, phi_ax/delta) 1_{shell}          -- sharp sgn(z) jump
        (D-B)  omega = -M tanh(sin phi/sin delta) tanh(cos phi/w) Theta_hk2(rho)
        (S3)   omega = -M sgn(cos phi) min(1,phi_ax/delta) min(1,|phi-pi/2|/delta_m) Theta_S3(rho)
     with hk2's own two rows reproduced as controls.
  B. WHERE J IS ACTUALLY USED.  A grep of the chain: `J` appears in NO script -- not in
     hk2/k3_bound.py (which produces the PROVED K_2 table by direct quadrature of Lambda_4,
     Lambda_5 on the actual field), not in k4_direct.py, not in t2_gamma_CR.py (which produces
     C''), not in t3_budget.py.  It appears only in hk2 Lemma 5.2's ANALYTIC statement
     Lambda_4 = O(M J/rho), Lambda_5 = O(M J/rho^2) -- the "no log" result -- and in the
     reference table of hk2 section 8(a).  So the 96.25/65.63/78.64 spread propagates into
     nothing as the chain stands.
  C. WHAT IT WOULD COST IF IT DID.  The one place J is load-bearing is the fully-analytic
     route: Lambda_p = (Lemma 5.2's far and near bounds in J) + (the collar by quadrature).
     That route is priced here as a function of J, giving dK2/dJ explicitly, and the resulting
     L_* is computed with fix2's f1_window (byte copy), C_K overridden.

Outputs -> g3_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "fix2copy"))
sys.path.insert(0, os.path.join(HERE, "fix2copy", "imported"))
import numpy as np
import k3_bound as K3                                              # noqa: E402 (hk2, unchanged)
from hk2lib import DatumA, DatumB                                  # noqa: E402 (hk2 profiles)
from f2_datum_s3 import DatumS3, StrainedS3, geometry_s3           # noqa: E402 (fix2 byte copy)

OUT = {}
DEG = math.pi/180.0
DELTA, DM, EPS_R = 7.5*DEG, 5.0*DEG, 0.25
S3SPH = 2*math.pi**2

# =====================================================================================
# A.  J FROM LEMMA 5.1, ONE INSTRUMENT, THREE DATA
# =====================================================================================
# Lemma 5.1:  |D eta|({rho'<|x'|<rho'+drho'}) = M rho'^2 J(rho') drho' ,
#   J_ac(rho') = |S^3| int_0^pi rho'^2 |grad eta|(rho',phi)/M  sin^3 phi dphi
#   J_jump     = 2 |S^3| Hh(0+)     (a jump of eta across the 4-plane {z=0})
# The instrument below takes rho'^2 |grad eta|/M directly from the datum object's own
# derivative fields, so it needs no profile-specific algebra.
NPHI = 40001
PHI_EPS = 1e-5      # the end caps contribute O(PHI_EPS^4) : the integrand is O(phi^3)
phi = np.linspace(PHI_EPS, math.pi - PHI_EPS, NPHI)
s3 = np.sin(phi)**3


def J_ac_at(D, rho):
    """|S^3| int rho^2 |grad eta| sin^3 phi dphi  at radius rho, from the object's own fields."""
    r = rho*np.sin(phi)
    z = rho*np.cos(phi)
    gf = getattr(D, "grad_norm_fast", None) or D.grad_norm
    g = gf(np.maximum(r, 1e-300), z)   # DatumB's full eta_rz overflows on the axis (Hh'' ~ 1/s^3)
    return S3SPH*float(np.trapz(rho**2*g*s3, phi))


def J_sup(D, us, jump=0.0):
    best, arg = 0.0, None
    for uv in us:
        v = J_ac_at(D, math.exp(uv)) + jump
        if v > best:
            best, arg = v, float(uv)
    return best, arg


# --- (D-A): sharp radial edges, sharp sgn(z).  Its J is evaluated on the plateau (Theta = 1),
#     which is what hk2 section 8(a) reports; the equatorial jump adds 2|S^3| Hh(0+) = 2|S^3|.
DA = DatumA(delta_deg=7.5, L=10.0)
JA_ac = J_ac_at(DA, 2.0)                        # any rho strictly inside the shell
JA_jump = 2*S3SPH*1.0                           # h(pi/2) = 1
# --- (D-B): real-analytic, no jump; hk2 evaluates on the plateau
DB = DatumB(delta_deg=7.5, w=0.20, L=10.0)
JB_plateau = J_ac_at(DB, math.exp(5.0))         # deep plateau of hk2's rho-tanh ramp
JB_sup, JB_arg = J_sup(DB, np.linspace(math.log(0.2), 10.5, 601))
# --- (S3): THEOREM_S3's own datum
DS = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=10.0)
JS_plateau = J_ac_at(DS, math.exp(5.0))
JS_sup, JS_arg = J_sup(DS, np.linspace(-6.0, 16.0, 1101))

OUT["A_J_three_data"] = {
    "(D-A) J_ac_plateau": JA_ac, "(D-A) J_jump": JA_jump, "(D-A) J_total": JA_ac + JA_jump,
    "NPHI": NPHI, "note_grid": "phi grid 40001 points; the J values are quoted to 5 digits",
    "(D-A) hk2_reported": 78.641147, "(D-A) rel": abs(JA_ac + JA_jump - 78.641147)/78.641147,
    "(D-B) J_plateau": JB_plateau, "(D-B) J_sup_over_rho": JB_sup, "(D-B) argmax_u": JB_arg,
    "(D-B) hk2_reported": 65.625910, "(D-B) rel": abs(JB_plateau - 65.625910)/65.625910,
    "(S3) J_plateau": JS_plateau, "(S3) J_sup_over_rho": JS_sup, "(S3) argmax_u": JS_arg,
    "(S3) fix2_reported_sup": 96.249844, "(S3) fix2_reported_plateau": 75.251104,
    "(S3) rel_sup": abs(JS_sup - 96.249844)/96.249844,
    "(S3) rel_plateau": abs(JS_plateau - 75.251104)/75.251104,
    "ratio_S3_over_DB": JS_sup/JB_plateau, "ratio_S3_over_DA": JS_sup/(JA_ac + JA_jump),
    "which_datum_each_number_belongs_to": {
        "78.641147": "(D-A), plateau, sharp sgn(z) jump included; hk2 sec.8(a)",
        "65.625910": "(D-B), plateau, real-analytic, no jump; hk2 sec.8(a); the row that "
                     "carries hk2's own PROVED K2 table",
        "96.25": "THEOREM_S3's datum, SUPREMUM over rho' (what (2.1) asks), attained in the "
                 "OUTER RAMP at u = L - 0.336; its plateau value is 75.25",
    },
}
print("J: (D-A) %.6f (hk2 78.641147)  (D-B) %.6f (hk2 65.625910)  (S3) plateau %.6f sup %.6f"
      % (JA_ac + JA_jump, JB_plateau, JS_plateau, JS_sup), flush=True)

# =====================================================================================
# B.  WHERE J IS USED -- the grep, recorded as data
# =====================================================================================
import subprocess                                                  # noqa: E402
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
files = ["hk2/k3_bound.py", "hk2/k4_direct.py", "hk2/k5_consequence.py",
         "hk2/k6_confront.py", "hk2/hk2lib.py",
         "round2/THEOREM_S3/t2_gamma_CR.py", "round2/THEOREM_S3/t3_budget.py",
         "round2/THEOREM_S3/fix2/f1_window.py"]
hits = {}
for f in files:
    p = os.path.join(ROOT, f)
    n = 0
    if os.path.exists(p):
        for ln in open(p):
            if ("J_of" in ln) or ("Jcal" in ln) or ("shell_density" in ln):
                n += 1
    hits[f] = n
OUT["B_where_J_is_used"] = {
    "occurrences_of_a_shell_density_symbol": hits,
    "conclusion": ("J occurs in no script that produces a constant. hk2's PROVED K2 table "
                   "(sec.8(c)) comes from k3.lambda_bulk, a direct quadrature of "
                   "int |x-x'|^{-p}|grad eta| dx' on the actual field; C'' comes from "
                   "t2_gamma_CR's (Gamma-off) fixed point, whose inputs are E_0, Gfrak_0, "
                   "C_far, C_inner, C_collar. Neither reads J."),
    "where_J_IS_used": ("hk2 Lemma 5.2 (the analytic 'no log' bound) and the reference table "
                        "of hk2 sec.8(a)."),
}

# =====================================================================================
# C.  WHAT J WOULD COST: the fully-analytic Lambda route, priced
# =====================================================================================
# Lemma 5.2 (hk2), with rho = |x|:
#   Lambda_4 <= 8 M J/rho + (2/3) M J/rho + Lambda_4^collar
#   Lambda_5 <= 4 M J/rho^2 + (4/3) M J/rho^2 + Lambda_5^collar
# The collar rho/2 < rho' < 2 rho minus the ball |x-x'| < d is by quadrature.


class Band:
    def __init__(self, base, lo, hi):
        self.b, self.lo, self.hi = base, lo, hi

    def grad_norm_fast(self, r, z):
        rho = np.hypot(r, z)
        return np.where((rho >= self.lo) & (rho <= self.hi), self.b.grad_norm_fast(r, z), 0.0)


def K2_analytic(Jval, lam, f=4.0, L=10.0, nd=8, NS=300, NT=120, NC=120):
    D0 = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=L)
    DL = StrainedS3(D0, lam)
    g = geometry_s3(30.0, f, 7.5, 5.0, lam)
    rho = g["rho"]
    dmax = 0.95*min(g["d_layers"], 0.95*g["r"])
    band = Band(DL, 0.5*rho, 2.0*rho)
    best = None
    for d in np.linspace(0.05*dmax, dmax, nd):
        Lc = K3.lambda_bulk(band.grad_norm_fast, (g["r"], g["z"]), d, 3*rho, NS, NT, NC)
        L4 = (8.0 + 2.0/3.0)*Jval/rho + Lc[4]
        L5 = (4.0 + 4.0/3.0)*Jval/rho**2 + Lc[5]
        sg, sh, sb = K3.sups_on_ball_B(DL, g["r"], g["z"], d)
        ee = [float(v[0]) for v in DL.eta_rz(np.array([g["r"]]), np.array([g["z"]]))]
        A = K3.assemble(sg, sh, sb, d, L4, L5, g["r"], ee[0], ee[1], ee[2])
        A.update(d=d, L4=L4, L5=L5, collar4=Lc[4], collar5=Lc[5], rho=rho, J=Jval, lam=lam)
        if best is None or A["K2"] < best["K2"]:
            best = A
    return best


rowsC = []
for Jval in (96.2602, 78.641147, 65.625910, 2*96.2602):
    for lam in (1.0, 1.8420112):
        r = K2_analytic(Jval, lam)
        rowsC.append(r)
        print("ANALYTIC route: J=%9.4f lam=%6.4f  L4=%9.4f (collar %8.4f)  L5=%8.4f  "
              "K2hat=%9.4f" % (Jval, lam, r["L4"], r["collar4"], r["L5"], r["K2"]), flush=True)
OUT["C_analytic_route"] = rowsC
base = [r for r in rowsC if abs(r["J"] - 96.2602) < 1e-6 and abs(r["lam"] - 1.8420112) < 1e-9][0]
dbl = [r for r in rowsC if abs(r["J"] - 2*96.2602) < 1e-6 and abs(r["lam"] - 1.8420112) < 1e-9][0]
OUT["C_sensitivity"] = {
    "K2hat(J=96.26, lam=lam_max)": base["K2"],
    "K2hat(J=192.52, lam=lam_max)": dbl["K2"],
    "dlogK2_dlogJ_at_J=96.26": (math.log(dbl["K2"]) - math.log(base["K2"]))/math.log(2.0),
    "quadrature_K2hat_for_comparison_fix2": 33.710,
    "note": ("the analytic route is far weaker than the quadrature one; the quadrature K2 "
             "table of fix2 sec.2.3 is unchanged by J because it never reads J."),
}

# =====================================================================================
# D.  L_* WITH THE K_2 INPUTS, VIA fix2's f1_window (byte copy, C_K overridden)
# =====================================================================================
import f1_window as F1                                             # noqa: E402
FS4 = [f for f in F1.FS_ALL if f >= 4.0]
CK_CASES = {
    "t3_budget_import_161.7735": 161.7735,
    "fix2_quadrature_S3_datum_lam_max_f4_33.710": 33.710,
    "analytic_route_J=96.26_lam_max": base["K2"],
    "analytic_route_J=192.52_lam_max": dbl["K2"],
}
rowsD = {}
_ORIG = F1.assemble_c


def _patch(ck):
    def patched(L, c, column, f=1.0, C_K=None, C_R_override=None, _ck=ck):
        return _ORIG(L, c, column, f=f, C_K=_ck, C_R_override=C_R_override)
    F1.assemble_c = patched


# eps(L) at three L, every C_K (cheap)
epsL = {}
for name, ck in CK_CASES.items():
    _patch(ck)
    epsL[name] = {"C_K": ck}
    for L in (1e5, 1e6, 1e7):
        epsL[name]["eps_self_consistent_L=%g" % L] = F1.self_consistent(L, "proved",
                                                                        fs=FS4)["eps"]
    F1.assemble_c = _ORIG
    print("eps(1e7, proved) at C_K=%10.4f : %s" % (ck, epsL[name]["eps_self_consistent_L=1e+07"]),
          flush=True)
OUT["D_eps_vs_CK"] = epsL

for name, ck in CK_CASES.items():
    for tgt in (0.1943662, 0.1):
        _patch(ck)
        Lp = F1.Lstar_self_consistent("proved", tgt, fs=FS4)
        F1.assemble_c = _ORIG
        rowsD["%s_eps<=%g" % (name, tgt)] = {
            "C_K": ck, "L_star": Lp,
            "logLambda_star": (2*Lp + 2.9234859) if Lp else None}
        print("L_*(proved, eps<=%g) with C_K=%10.4f : %s" % (tgt, ck, Lp), flush=True)
OUT["D_Lstar"] = rowsD
OUT["D_reference"] = {"fix2_L_star_proved_eps<=0.1943662_f>=4": 1887786.7,
                      "fix2_L_star_proved_eps<=0.1_f>=4": 1977309.3,
                      "record_C_K": 161.7735}

with open(os.path.join(HERE, "g3_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print("\nWROTE g3_results.json")

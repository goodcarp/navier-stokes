"""
g2 -- UNIT 2 of FIX3: (Theta-tail), item M7 of THEOREM_S3 section 3.

WHAT M7 SAYS.  hk2's hypothesis (D1) asks for  supp eta(.,t) subset {rho_-(t) <= rho <= rho_+(t)}.
THEOREM_S3's datum has tanh ramps in u = log rho, so

    psi_in(rho)  = e^2/(e^2 + (rho/rho_0)^8)          (Theta = 1 - psi_in - psi_out)
    psi_out(rho) = 1/(1 + e^{-8(u-L)-2})              (its outer mirror)

and supp omega_0 = R^3 with O(rho^8) at the origin and O(rho^{-8}) at infinity.  The tails are
"exponentially small in L and nobody has written the estimate".

WHAT THIS SCRIPT ESTABLISHES.

  A. Theta and its tails, EXACTLY.
       int_{-infinity}^{infinity} Theta(u) du = L - 2 eps_r      (exact, = L - 0.5)
       int_{-infinity}^{0}        Theta(u) du = (eps_r/2) log(1 + e^{-2}) = log(1+e^{-2})/8
     and the same at the outer end.  So the tails are worth 0.0158660 e-folds EACH, and the
     ramps cost 2 eps_r = 0.5 e-folds in total.  t3_budget.py:65 charges ELL_RAMP = 0.25.
  B. the local shell density J(rho') in the tails, with the closed-form envelope
       J(rho') <= Theta(rho') * Jcoef_in   (inner, Theta_u = 8 Theta)
       J(rho') <= Theta(rho') * Jcoef_out  (outer, Theta_u = -8 Theta)
     and Theta <= e^{-2} (rho/rho_0)^8 inside, <= e^{-2}(R/rho)^8 outside -- so the tails LOWER
     the local density and the supremum that (2.1) asks for is attained in the outer RAMP.
  C. the tails' share of Lambda_4 and Lambda_5, analytically (closed form) and by hk2's own
     k3.lambda_bulk with the field masked outside the shell.
  D. the corrected hypothesis (D1') and what changes.

Outputs -> g2_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "fix2copy"))
sys.path.insert(0, os.path.join(HERE, "fix2copy", "imported"))
import numpy as np
import sympy as sp
import k3_bound as K3                                            # noqa: E402 (hk2, unchanged)
from f2_datum_s3 import DatumS3, StrainedS3, geometry_s3         # noqa: E402 (fix2 byte copy)

OUT = {}
DEG = math.pi/180.0
DELTA, DM, EPS_R = 7.5*DEG, 5.0*DEG, 0.25
S3SPH = 2*math.pi**2

# =====================================================================================
# A.  Theta AND ITS TAILS, EXACTLY
# =====================================================================================
u, a, e = sp.symbols('u a e', positive=True)
sig = (1 + sp.tanh((u - a)/e))/2
prim = (e/2)*sp.log(1 + sp.exp(2*(u - a)/e))
res = sp.simplify(sp.diff(prim, u) - sig)
_res_num = max(abs(float(res.subs({u: uv, a: 0.25, e: 0.25}).evalf()))
               for uv in (0.01, 0.1, 0.25, 1.0, 3.0))
tail_exact = sp.simplify(prim.subs(u, 0) - sp.limit(prim, u, -sp.oo))
OUT["A_theta_exact"] = {
    "antiderivative_residual (sympy)": str(res),
    "antiderivative_residual_numeric_max": _res_num,
    "int_{-inf}^{0} (1+tanh((u-a)/e))/2 du": str(sp.simplify(tail_exact)),
    "at a = eps_r = e": str(sp.simplify(tail_exact.subs(a, e))),
    "value_eps_r_0.25": float(sp.simplify(tail_exact.subs({a: sp.Rational(1, 4),
                                                           e: sp.Rational(1, 4)}))),
    "closed_form": "log(1 + e^{-2})/8",
    "closed_form_value": math.log(1.0 + math.exp(-2.0))/8.0,
    "int_Theta_du_whole_line": "L - 2 eps_r  (exact)",
    "ramp_loss_efolds_true": 2*EPS_R,
    "t3_budget_ELL_RAMP": 0.25,
    "understatement_efolds": 2*EPS_R - 0.25,
}
# numeric control at L = 12, and the shell-only integral
Lc = 12.0
ug = np.linspace(-25.0, Lc + 25.0, 8000001)
Th = 0.5*(np.tanh((ug - EPS_R)/EPS_R) - np.tanh((ug - Lc + EPS_R)/EPS_R))
OUT["A_theta_exact"]["numeric_whole_line_L12"] = float(np.trapz(Th, ug))
m_in = ug < 0.0
m_out = ug > Lc
m_sh = (ug >= 0.0) & (ug <= Lc)
OUT["A_theta_exact"]["numeric_inner_tail"] = float(np.trapz(Th[m_in], ug[m_in]))
OUT["A_theta_exact"]["numeric_outer_tail"] = float(np.trapz(Th[m_out], ug[m_out]))
OUT["A_theta_exact"]["numeric_shell_only"] = float(np.trapz(Th[m_sh], ug[m_sh]))
# the envelope bounds Theta <= e^{-2} e^{8u} (u<0) and <= e^{-2} e^{-8(u-L)} (u>L)
uu = np.linspace(-8.0, 0.0, 200001)
Thn = 0.5*(np.tanh((uu - EPS_R)/EPS_R) - np.tanh((uu - Lc + EPS_R)/EPS_R))
_msk = Thn > 1e-10                       # below this Theta is float cancellation noise
OUT["A_theta_exact"]["envelope_inner_max_ratio"] = float(
    np.max(Thn[_msk]/(math.exp(-2.0)*np.exp(8.0*uu[_msk]))))
OUT["A_theta_exact"]["envelope_exact"] = ("Theta = sigma(8u-2) at the inner end and sigma(x) "
                                          "<= e^x, so Theta <= e^{-2} (rho/rho_0)^8 EXACTLY")
uu2 = np.linspace(Lc, Lc + 8.0, 200001)
Thn2 = 0.5*(np.tanh((uu2 - EPS_R)/EPS_R) - np.tanh((uu2 - Lc + EPS_R)/EPS_R))
_msk2 = Thn2 > 1e-10
OUT["A_theta_exact"]["envelope_outer_max_ratio"] = float(
    np.max(Thn2[_msk2]/(math.exp(-2.0)*np.exp(-8.0*(uu2[_msk2] - Lc)))))

# =====================================================================================
# B.  THE LOCAL SHELL DENSITY J(rho') AND ITS TAILS
# =====================================================================================
# hk2 Lemma 5.1 with the ramp kept:
#   |D eta|({rho'<|x'|<rho'+drho'}) = M rho'^2 J(rho') drho' ,
#   J(rho') = |S^3| int_0^pi sqrt( (Theta_u - Theta)^2 W^2 + Theta^2 W_phi^2 ) sin^3 phi dphi .
# The theorem's datum has no jump surface (A is C^1 across the equator), so J = J_ac.


def A_pieces(phi):
    s_sign = np.where(np.cos(phi) >= 0.0, 1.0, -1.0)
    pax = np.minimum(phi, math.pi - phi)
    dpax = np.where(phi < math.pi/2, 1.0, -1.0)
    eqd = np.abs(phi - math.pi/2)
    deqd = np.where(phi > math.pi/2, 1.0, -1.0)
    g1 = np.minimum(1.0, pax/DELTA)
    dg1 = np.where(pax < DELTA, dpax/DELTA, 0.0)
    g2 = np.minimum(1.0, eqd/DM)
    dg2 = np.where(eqd < DM, deqd/DM, 0.0)
    return s_sign*g1*g2, s_sign*(dg1*g2 + g1*dg2)


def Th_u(uv, k=0, L=10.0):
    aa = (uv - EPS_R)/EPS_R
    bb = (uv - L + EPS_R)/EPS_R
    if k == 0:
        return 0.5*(np.tanh(aa) - np.tanh(bb))
    ca = 1.0/np.cosh(np.clip(aa, -350, 350))**2
    cb = 1.0/np.cosh(np.clip(bb, -350, 350))**2
    return 0.5*(ca - cb)/EPS_R


NPHI = 400001
phi = np.linspace(1e-9, math.pi - 1e-9, NPHI)
s_, c_ = np.sin(phi), np.cos(phi)
A0, A1 = A_pieces(phi)
W0 = A0/s_
W1 = (A1*s_ - A0*c_)/s_**2
s3 = s_**3
I_W = float(np.trapz(np.abs(W0)*s3, phi))                     # int |W| sin^3
I_Wp = float(np.trapz(np.abs(W1)*s3, phi))


def J_of_u(uv, L=10.0):
    T0 = Th_u(uv, 0, L)
    T1 = Th_u(uv, 1, L)
    return S3SPH*float(np.trapz(np.sqrt((T1 - T0)**2*W0**2 + T0**2*W1**2)*s3, phi))


# the closed-form tail coefficients:  in the inner tail Theta_u = 8 Theta + O(Theta^2),
# in the outer tail Theta_u = -8 Theta + O(Theta^2), so J <= Theta * Jcoef with
Jcoef_in = S3SPH*float(np.trapz(np.sqrt(49.0*W0**2 + W1**2)*s3, phi))
Jcoef_out = S3SPH*float(np.trapz(np.sqrt(81.0*W0**2 + W1**2)*s3, phi))
J_plateau = S3SPH*float(np.trapz(np.sqrt(W0**2 + W1**2)*s3, phi))

for LL in (10.0, 40.0):
    us = np.linspace(-6.0, LL + 6.0, 6001)
    vals = [J_of_u(x, LL) for x in us]
    i = int(np.argmax(vals))
    key = "B_J_of_rho_L%g" % LL
    OUT[key] = {
        "J_sup_over_rho": float(vals[i]), "argmax_u": float(us[i]),
        "argmax_u_minus_L": float(us[i] - LL),
        "J_plateau_Theta_eq_1": J_plateau,
        "Jcoef_inner_tail": Jcoef_in, "Jcoef_outer_tail": Jcoef_out,
        "J_at_u=0": J_of_u(0.0, LL), "J_at_u=-1": J_of_u(-1.0, LL),
        "J_at_u=-2": J_of_u(-2.0, LL), "J_at_u=L": J_of_u(LL, LL),
        "J_at_u=L+1": J_of_u(LL + 1.0, LL), "J_at_u=L+2": J_of_u(LL + 2.0, LL),
        "sup_J_over_inner_tail_(u<0)": max(J_of_u(x, LL) for x in np.linspace(-6.0, 0.0, 121)),
        "sup_J_over_outer_tail_(u>L)": max(J_of_u(x, LL)
                                           for x in np.linspace(LL, LL + 6.0, 121)),
        "envelope_check_inner_u=-1": J_of_u(-1.0, LL)/(math.exp(-2.0)*math.exp(-8.0)*Jcoef_in),
        "envelope_check_outer_u=L+1": J_of_u(LL + 1.0, LL)/(math.exp(-2.0)*math.exp(-8.0)
                                                            * Jcoef_out),
    }
    print("L=%g: sup J = %.6f at u = L%+.3f ; J(0)=%.4f J(L)=%.4f ; tail sups %.4g / %.4g"
          % (LL, OUT[key]["J_sup_over_rho"], OUT[key]["argmax_u_minus_L"],
             OUT[key]["J_at_u=0"], OUT[key]["J_at_u=L"],
             OUT[key]["sup_J_over_inner_tail_(u<0)"],
             OUT[key]["sup_J_over_outer_tail_(u>L)"]), flush=True)

# =====================================================================================
# C.  THE TAILS' SHARE OF Lambda_4 AND Lambda_5
# =====================================================================================
# hk2 Lemma 5.2 integrates (2.1) over rho' in (0, rho/2) and (2 rho, infinity) ALREADY:
#     int_{rho'>2rho} |x-x'|^{-4} d|D eta| <= 16 M J int_{2rho}^inf rho'^{-2} drho' = 8 M J/rho
#     int_{rho'<rho/2} |x-x'|^{-4} d|D eta| <= 16 rho^{-4} M J (rho/2)^3/3
# with a J that is UNIFORM in rho'.  So the support hypothesis is not used by Lemma 5.2 at all;
# what is used is that (2.1) holds at every rho'.  The closed-form tail shares are:
#   inner tail (rho' < rho_0):  int_0^{rho_0} J(rho') rho'^2 drho'
#                                  <= Jcoef_in e^{-2} rho_0^3/11
#   outer tail (rho' > R):      int_R^inf J(rho') rho'^{-2} drho'  <= Jcoef_out e^{-2}/(9 R)
rho0 = 1.0
for f in (1.0, 4.0):
    rho_star = 1.0 + f
    R = math.exp(10.0)
    near_tail = Jcoef_in*math.exp(-2.0)*rho0**3/11.0
    near_hk2 = OUT["B_J_of_rho_L10"]["J_sup_over_rho"]*(rho_star/2.0)**3/3.0
    far_tail = Jcoef_out*math.exp(-2.0)/(9.0*R)
    far_hk2 = OUT["B_J_of_rho_L10"]["J_sup_over_rho"]/(2.0*rho_star)
    OUT["C_tail_share_f%g" % f] = {
        "rho_star": rho_star,
        "inner_tail_contribution_to_near_integral": near_tail,
        "hk2_uniform_J_near_integral": near_hk2,
        "inner_tail_relative_share": near_tail/near_hk2,
        "outer_tail_contribution_to_far_integral": far_tail,
        "hk2_uniform_J_far_integral": far_hk2,
        "outer_tail_relative_share": far_tail/far_hk2,
        "outer_share_closed_form": "(2/9) e^{-2} (Jcoef_out/J) (rho_*/R) = O(e^{-L})",
    }
    print("f=%g: inner tail share of hk2's near bound = %.4e ; outer tail share of the far "
          "bound = %.4e" % (f, near_tail/near_hk2, far_tail/far_hk2), flush=True)

# --- and the measured Lambda_4, Lambda_5 with the tails masked off -------------------


class Masked:
    """the same field with |grad eta| zeroed outside rho_lo < rho < rho_hi (labels), i.e. the
       compactly-supported idealisation (D1) asks for."""
    def __init__(self, base, lo, hi):
        self.b = base
        self.lo = lo
        self.hi = hi

    def grad_norm_fast(self, r, z):
        rho = np.hypot(r, z)
        g = self.b.grad_norm_fast(r, z)
        return np.where((rho >= self.lo) & (rho <= self.hi), g, 0.0)


D0 = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=10.0)
lam_list = [1.0, 1.8420112]
rowsC = []
for lam in lam_list:
    DL = StrainedS3(D0, lam)
    g = geometry_s3(30.0, 4.0, 7.5, 5.0, lam)
    d = 0.5*min(g["d_layers"], 0.95*g["r"])
    full = K3.lambda_bulk(DL.grad_norm_fast, (g["r"], g["z"]), d, 3*D0.R, 300, 120, 120)
    # the transported shell is T_lam({rho_0 < rho < R}); mask in the LABEL radius
    mk = Masked(DL, lam**-2*1.0, lam*D0.R)
    shell = K3.lambda_bulk(mk.grad_norm_fast, (g["r"], g["z"]), d, 3*D0.R, 300, 120, 120)
    row = {"lam": lam, "d": d, "L4_full": full[4], "L4_shell_only": shell[4],
           "L5_full": full[5], "L5_shell_only": shell[5],
           "L4_tail_rel": (full[4]-shell[4])/full[4],
           "L5_tail_rel": (full[5]-shell[5])/full[5]}
    rowsC.append(row)
    print("lam=%6.4f: L4 full %.6f vs shell-only %.6f (tails %.3e) ; L5 %.6f vs %.6f "
          "(tails %.3e)" % (lam, full[4], shell[4], row["L4_tail_rel"],
                            full[5], shell[5], row["L5_tail_rel"]), flush=True)
OUT["C_lambda_masked"] = rowsC

# =====================================================================================
# D.  THE CORRECTED HYPOTHESIS (D1')
# =====================================================================================
OUT["D_D1prime"] = {
    "old_(D1)": ("|omega^theta(.,t)| <= M_t <= (3/2)M and eta(.,t) supported in "
                 "{rho_-(t) <= rho <= rho_+(t)}"),
    "new_(D1')": ("|omega^theta(.,t)| <= M_t <= (3/2)M on R^5, and the shell total-variation "
                  "density (2.1) holds for EVERY rho' in (0, infinity) with one finite J. "
                  "No support statement is used anywhere in hk2 sections 4-5."),
    "why_it_is_enough": ("Lemma 5.2's four integrals already run over (2 rho, infinity) and "
                         "(0, rho/2); the support hypothesis is never invoked in their proof. "
                         "The collar rho/2 < rho' < 2 rho is evaluated by quadrature on the "
                         "actual field (k3.lambda_bulk), which integrates the tails as they "
                         "are. So (D1') changes no constant."),
    "far_near_lemma": ("far-near-kernel-lemma Propositions 1 and 2 need |omega^theta| <= M on "
                       "{rho' >= 2 rho} and on {rho' <= rho/2}, which holds everywhere; the "
                       "shell support enters only the e-fold COUNT a_far(0) = int Phi[w](0) "
                       "dlog rho', where the tails add exactly log(1+e^{-2})/8 = 0.0158660 "
                       "e-folds at each end, with the SAME sign as the drive."),
    "efold_accounting": {
        "int_Theta_dlogrho_whole_line": "L - 2 eps_r = L - 0.5",
        "each_tail": math.log(1.0 + math.exp(-2.0))/8.0,
        "shell_only": "L - 2 eps_r - 2 log(1+e^{-2})/8 = L - 0.5317320",
        "budget_charges": 0.25,
        "extra_efolds_owed": 2*EPS_R - 0.25,
        "cost_in_eps_a": {"L=1e4": 0.25/1e4, "L=1e5": 0.25/1e5, "L=1e6": 0.25/1e6},
    },
}

with open(os.path.join(HERE, "g2_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print("\nWROTE g2_results.json")
print(json.dumps(OUT["A_theta_exact"], indent=1, sort_keys=True, default=str))

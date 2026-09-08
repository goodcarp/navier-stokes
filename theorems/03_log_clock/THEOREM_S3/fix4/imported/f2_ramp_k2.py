"""
f2 -- ITEMS 2 and 3.

ITEM 2 (the ramp slope).  THEOREM_S3 sec.1.1 says, of
        Theta(rho) = (1/2)[ tanh((u-eps_r)/eps_r) - tanh((u-L+eps_r)/eps_r) ] ,  u = log(rho/rho_0),
    "Written in the log variable, rho Theta'(rho_0) = 1/(2 eps_r) = 2, which is exactly hk2's
     rho_0/(2 w_0) at w_0 = 0.25 rho_0: the inner slope agrees with hk2 (D-B) to the digit."

    Both halves of that sentence are wrong, in different ways, and they are separated here:

      rho Theta'(rho) = dTheta/du  EXACTLY (chain rule).  Its MAXIMUM over rho is
      1/(2 eps_r) = 2, attained at u = eps_r, i.e. at rho = rho_0 e^{eps_r} = 1.2840 rho_0,
      NOT at rho_0.  At rho_0 itself
            rho_0 Theta'(rho_0) = [ sech^2(1) - sech^2((L-eps_r)/eps_r) ] / (2 eps_r)
                                -> 2 sech^2(1) = 0.83994868  as L -> infinity,
      which is 42 % of the claimed 2.

      hk2's (D-B) ramp is tanh in RHO, not in log rho:
            Theta_hk2(rho) = (1/2)[ tanh((rho-rho_0)/w_0) - tanh((rho-R)/w_1) ] ,
            rho Theta_hk2'(rho) = rho[ sech^2((rho-rho_0)/w_0)/w_0 - ... ]/2 ,
            rho_0 Theta_hk2'(rho_0) = rho_0/(2 w_0) - (exp. small) = 2 - eps  at w_0 = rho_0/4.
      So the two profiles agree in the maximum log-slope (2) and disagree at rho_0
      (0.83995 against 2.0).  Neither coincidence is what K_2 depends on.

    WHAT K_2 ACTUALLY DEPENDS ON (hk2 sec.4-5, sec.8): three things, no point slope among them.
      (i)   the SHELL TOTAL-VARIATION DENSITY  J  of hypothesis (2.1),
                |D eta|({rho'<|x'|<rho'+drho'}) <= M rho'^2 J drho' ,
            an integral over the whole angular profile and a supremum over rho' -- the radial
            ramp enters only through the factor sup_rho sqrt((Theta_u-Theta)^2 W^2 + Theta^2 W_phi^2);
      (ii)  the local suprema on the ball B(x,d):  sup|grad eta|, sup||Hess eta||_F, sup_{dB}|grad eta|;
      (iii) the geometry: the distance d from the tracked point to the kink set.
    In (i)-(iii) the radial ramp matters through Theta_u and Theta_uu on the SHELL, and through
    whether the tracked point sits inside the ramp at all.  hk2's own reason for demanding a
    ramp is (2.1)/(D3): a SHARP radial edge makes eta jump, a is then a single-layer potential
    and ||Hess a|| diverges like 1/dist.  Both mollifications kill that; the exponent of the
    tanh is irrelevant to it.  ADDENDUM_1 puts the tracked point at f >= 4, i.e. rho_* = 5 rho_0,
    where BOTH ramps are flat to 10^{-5}, so the near-ball sups are the plateau's and the whole
    radial-ramp question reduces to J.

ITEM 3 (the lambda range).  t3_budget.py:44 imports K_2 <= 161.7735 as a bound over
    lambda in [1, 3/2] while the clock's a priori cap is lambda_max = e^{3c/4} = 1.8420112 (and,
    with item 1's corrected window, up to 2.0741).  Also hk2's constant is for hk2's own datum.
    Both are repaired here by re-running hk2's instrument on THEOREM_S3's datum, over the full
    lambda range, with everything else unchanged.

INSTRUMENT.  imported/hk2lib.py, imported/k3_bound.py, imported/k4_direct.py are BYTE COPIES of
s3close/hk2 (hashes in SHA256SUMS, verified against s3close/hk2/SHA256SUMS).  k3's
lambda_bulk(), assemble() and sups_on_ball_B() are used UNCHANGED; only the datum object is
swapped, for f2_datum_s3.DatumS3 / StrainedS3.  Nothing in s3close/hk2 is written to.

Outputs -> f2_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "imported"))
sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import k3_bound as K3                                   # noqa: E402  (hk2's instrument)
from hk2lib import DatumB, StrainedB, geometry          # noqa: E402  (hk2's own datum)
from f2_datum_s3 import DatumS3, StrainedS3, geometry_s3  # noqa: E402

OUT = {}
DEG = math.pi/180.0
DELTA, DM, EPS_R = 7.5*DEG, 5.0*DEG, 0.25
S3SPH = 2*math.pi**2

# =============================================================== A. the ramp slope, exactly
u, e, Lsym, rho, w0 = sp.symbols('u e L rho w0', positive=True)
Th_log = (sp.tanh((u - e)/e) - sp.tanh((u - Lsym + e)/e))/2
Th_of_rho = Th_log.subs(u, sp.log(rho))
res_chain = sp.simplify(rho*sp.diff(Th_of_rho, rho) - sp.diff(Th_log, u).subs(u, sp.log(rho)))
dTh = sp.diff(Th_log, u)
dTh_inf = sp.limit(dTh, Lsym, sp.oo)                     # single-edge limit
crit = sp.simplify(sp.diff(dTh_inf, u))
OUT["A_ramp_exact"] = {
    "rho dTheta/drho - dTheta/du (sympy)": str(res_chain),
    "dTheta/du, L->infinity": str(sp.simplify(dTh_inf)),
    "d/du of that, zero at u = eps_r": str(sp.simplify(crit.subs(u, e))),
    "max_log_slope_1_over_2eps_r": 1.0/(2*EPS_R),
    "argmax_u_equals_eps_r": EPS_R,
    "argmax_rho_over_rho0": math.exp(EPS_R),
}
def rho_dTheta(uv, L=40.0, er=EPS_R):
    return 0.5*(1.0/math.cosh((uv-er)/er)**2 - 1.0/math.cosh((uv-L+er)/er)**2)/er
OUT["A_ramp_exact"]["rho0_dTheta_drho_at_rho0_L40"] = rho_dTheta(0.0, 40.0)
OUT["A_ramp_exact"]["two_sech2_1"] = 2.0/math.cosh(1.0)**2
OUT["A_ramp_exact"]["value_at_u=eps_r_L40"] = rho_dTheta(EPS_R, 40.0)
OUT["A_ramp_exact"]["THEOREM_S3_claim_at_rho0"] = 2.0
# hk2's ramp, in rho
def rho_dTheta_hk2(r_, w0v=0.25, R=math.exp(10.0), w1f=0.10):
    w1 = w1f*R
    return 0.5*r_*(1.0/math.cosh((r_-1.0)/w0v)**2/w0v - 1.0/math.cosh((r_-R)/w1)**2/w1)
gr = np.linspace(0.5, 3.0, 500001)
vals = np.array([rho_dTheta_hk2(x) for x in gr])
OUT["A_ramp_exact"]["hk2_rho0_dTheta_at_rho0"] = rho_dTheta_hk2(1.0)
OUT["A_ramp_exact"]["hk2_max_rho_dTheta"] = float(vals.max())
OUT["A_ramp_exact"]["hk2_argmax_rho"] = float(gr[int(np.argmax(vals))])
OUT["A_ramp_exact"]["verdict"] = (
    "max log-slope agrees (2 = 1/(2 eps_r) = rho_0/(2 w_0)); the value AT rho_0 does not "
    "(0.83995 for the log-tanh ramp, 2.0 for hk2's rho-tanh ramp).  THEOREM_S3 sec.1.1's "
    "'inner slope agrees with hk2 (D-B) to the digit' compares the max of one with the "
    "point value of the other.")

# =============================================================== B. the shell density J
def J_theorem(D, nphi=400001, nu=4001):
    """J = |S^3| sup_rho int_0^pi sqrt((Theta_u-Theta)^2 W^2 + Theta^2 W_phi^2) sin^3 phi dphi ,
       the uniform-in-rho' constant hypothesis (2.1) asks for (hk2 Lemma 5.1 with the ramp kept)."""
    phi = np.linspace(1e-9, math.pi-1e-9, nphi)
    W0 = D.W(phi, 0)
    W1 = D.W(phi, 1)
    s3 = np.sin(phi)**3
    IW = np.trapz(W0**2*s3, phi)          # placeholder, real work below
    us = np.linspace(-2.0, D.L + 2.0, nu)
    best, arg = 0.0, None
    for uu in us:
        Th = D.Th_u(uu, 0)
        Tu = D.Th_u(uu, 1)
        integ = np.sqrt((Tu - Th)**2*W0**2 + Th**2*W1**2)*s3
        val = S3SPH*np.trapz(integ, phi)
        if val > best:
            best, arg = float(val), float(uu)
    plateau = S3SPH*float(np.trapz(np.sqrt(W0**2 + W1**2)*s3, phi))
    return {"J_uniform_sup_over_rho": best, "argmax_u": arg, "J_plateau_Theta_eq_1": plateau,
            "unused": float(IW)}

DS3 = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=40.0)
OUT["B_shell_density"] = J_theorem(DS3)
OUT["B_shell_density"]["hk2_J_DA_delta7.5_total"] = 78.641147
OUT["B_shell_density"]["hk2_J_DB_delta7.5_w0.20"] = 65.625910
OUT["B_shell_density"]["hk2_J_bare_plateau_ac"] = 2*math.pi**2*(math.sqrt(2) + math.asinh(1))

# =============================================================== C. the PROVED K_2 bound
def k2_bound(D0, lam, f, nd=8, NS=250, NT=100, NC=100, frac=0.95):
    """hk2's k3 assembly, unchanged, on the theorem's datum transported by T_lam."""
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
    best.update(lam=lam, f=f, r=g["r"], z=g["z"], rho=g["rho"], phi_deg=g["phi_deg"],
                d_layers=g["d_layers"], dmax=dmax, L=D0.L)
    return best

LAMS = [1.0, 1.25, 1.5, 1.8420112, 2.0741]     # 1.8420112 = e^{3c_*/4}; 2.0741 = e^{3 c_cap/4}
if __name__ == "__main__":
    rows = []
    for Lv in (10.0, 40.0):
        D0 = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=Lv)
        for f in (1.0, 4.0):
            for lam in LAMS:
                r = k2_bound(D0, lam, f)
                rows.append(r)
                print("L=%4.1f f=%3.1f lam=%6.4f: d*=%.4f |grad a|<=%.4f ||Hess a||<=%.4f "
                      "||Hess u^z||<=%.4f  K2hat <= %.4f"
                      % (Lv, f, lam, r["d"], r["grad_a"], r["hess_a"], r["hess_uz"], r["K2"]))
    OUT["C_K2_bound_theorem_datum"] = rows
    tab = {}
    for r in rows:
        tab["L=%g_f=%g_lam=%g" % (r["L"], r["f"], r["lam"])] = r["K2"]
    OUT["C_K2_table"] = tab
    OUT["C_no_log_control_R5"] = {
        "rel_move_L10_to_L40": {
            "f=%g_lam=%g" % (f, lam):
                abs(tab["L=40_f=%g_lam=%g" % (f, lam)] - tab["L=10_f=%g_lam=%g" % (f, lam)])
                / tab["L=10_f=%g_lam=%g" % (f, lam)]
            for f in (1.0, 4.0) for lam in LAMS}}
    OUT["C_hk2_quoted_for_its_own_datum"] = {"lam=1": 66.6622, "lam=1.25": 107.6738,
                                             "lam=1.5": 161.7735, "f": 0.0,
                                             "used_in_t3_budget_as": 161.7735,
                                             "t3_import_range": "[1, 3/2]"}

    # quadrature convergence control (R2-style)
    D0 = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=10.0)
    conv = {}
    for NS, NT, NC in [(250, 100, 100), (450, 180, 180)]:
        r = k2_bound(D0, 1.5, 4.0, nd=6, NS=NS, NT=NT, NC=NC)
        conv["NS=%d" % NS] = {"K2": r["K2"], "L4": r["L4"], "L5": r["L5"], "d": r["d"]}
    conv["rel_move_K2"] = abs(conv["NS=450"]["K2"] - conv["NS=250"]["K2"])/conv["NS=250"]["K2"]
    OUT["C_quadrature_convergence"] = conv

    # ============================================================ D. the MEASURED K_2 at lam = 1
    # k4's zonal-harmonic instrument needs Hh(t) = W(arccos t); the strained field eta_lam is
    # NOT of the separable form G(rho)H(t), so this instrument measures at lam = 1 only.
    import k4_direct as K4
    class _Adapt(DatumS3):
        def Hh(self, t, k=0):
            return self.W(np.arccos(np.clip(t, -1.0, 1.0)), k)
    meas = {}
    D0 = _Adapt(delta_deg=7.5, dm_deg=5.0, eps_r=EPS_R, L=10.0)
    for LMAX in (81, 161, 321):
        S = K4.Series(D0, LMAX=LMAX)
        g = geometry_s3(30.0, 4.0, 7.5, 5.0, 1.0)
        r_, z_ = g["r"], g["z"]
        A = S.a_field(r_, z_)
        ee = [float(v[0]) for v in D0.eta_rz(np.array([r_]), np.array([z_]))]
        lap = A[3] + 3*A[1]/r_ + A[5]
        k2m, arg = K4.K2_exact(r_, A[1], A[2], A[3], A[4], A[5], ee[0], ee[1], ee[2], ne=121)
        ha = np.linalg.eigvalsh(np.array([[A[3], A[4]], [A[4], A[5]]]))
        meas["LMAX=%d" % LMAX] = {"a": A[0], "grad_a": math.hypot(A[1], A[2]),
                                  "hess_a": max(abs(ha).max(), abs(A[1]/r_)),
                                  "K2_measured": k2m,
                                  "control_Lap5a_minus_dz_eta_rel":
                                      abs(lap - ee[2])/max(abs(ee[2]), 1e-30),
                                  "r": r_, "z": z_}
        print("MEASURED LMAX=%d: K2 = %.6f  (Lap5 a vs d_z eta rel = %.2e)"
              % (LMAX, k2m, meas["LMAX=%d" % LMAX]["control_Lap5a_minus_dz_eta_rel"]))
    OUT["D_measured_lam1_f4"] = meas
    kb = tab["L=10_f=4_lam=1"]
    OUT["D_slack_lam1_f4"] = {"bound": kb, "measured_LMAX321": meas["LMAX=321"]["K2_measured"],
                              "slack": kb/meas["LMAX=321"]["K2_measured"]}

    with open(os.path.join(HERE, "f2_results.json"), "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
    print(json.dumps({k: OUT[k] for k in OUT if not k.startswith("C_K2_bound")},
                     indent=1, sort_keys=True, default=str))

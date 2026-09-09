"""
x2 -- UNIT F-2.  The budget run with (H3'_vartheta) actually carried, and with the corrected
Theorem V.4' ball radius propagated everywhere d enters.

WHAT WAS WRONG.  pmax-h3v sec.B5 states the corrected viscous line

    eps_v(vartheta) = eps_bulk [ 1 + vt (9 + 2 sigma_z/sigma_y) ] + Q ( E_hess + E_4 + E_tail ) ,
    Q := (1+vartheta)/(1-vartheta) ,   vt := vartheta/(1-vartheta) ,

and supplies q4_budget_theta.py to carry it.  fix4 imports `t3_budget.viscous_budget` unchanged,
which returns eps_v = eps_bulk + E_hess + E_4 + E_tail: the vartheta = 0 line, i.e. Theorem V.4
and not V.4'.  Row B4 of THEOREM_S3_v2 claims V.4' while row F4 quotes V.4.

WHAT IS DONE HERE.  `viscous_budget_theta` is `t3_budget.viscous_budget` transcribed with two
changes and nothing else: the ball radius d is exposed (so the F-1 repair
d = rho_* sin(phi_0 - delta - 4 sigma) can be propagated into E_hess's r_*/R_-^2 and 4 N/d terms,
into E_4's r_*/R_-^5, and into E_tail's three Gaussian exit probabilities), and eps_v is
assembled by Corollary V.5' above instead of by the vartheta = 0 sum.  The transcription is
checked against the imported function at d = ASSEMBLY's d and vartheta = 0 before use.  Every
other part of the budget -- fix4's p4 window search, bootstrap, column constants, C_R, the
(Gamma-off) root -- is imported and untouched.

vartheta(f) at the corrected ball comes from x1_results.json, multiplied by SAFETY = 2
(pmax-h3v sec.B5's convention for a grid lower bound on a supremum).

The certified-constants column at the end is UNIT F-6's last step: E_0, Gfrak_0, r_h and the
window cap replaced by x6's certified values, to see whether L_* moves.

Outputs -> x2_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIX4 = os.path.join(HERE, "copies", "fix4")
sys.path.insert(0, HERE)
sys.path.insert(0, FIX4)
sys.path.insert(0, os.path.join(HERE, "copies", "fix4imported"))
_cwd = os.getcwd()
import p4_budget as P4                                            # noqa: E402  (chdirs itself)
os.chdir(_cwd)
T3 = P4.T3
DEG = math.pi/180.0
SIGMA = 0.02
DELTA, DM, PHI0 = 7.5*DEG, 5.0*DEG, 30.0*DEG
SHIFT = 2.9135781820
X1 = json.load(open(os.path.join(HERE, "x1_results.json")))
X6 = json.load(open(os.path.join(HERE, "x6_results.json")))
SAFETY = 2.0

VT = {4.0: X1["repaired"]["f=4"]["vartheta"], 8.0: X1["repaired"]["f=8"]["vartheta"],
      16.0: X1["repaired"]["f=16"]["vartheta"], 32.0: X1["repaired"]["f=32"]["vartheta"]}

# ---------------------------------------------------------------- the corrected ball radius
def d_of(f, offset):
    rho_s = 1.0 + f
    return min(rho_s*math.sin(PHI0 - DELTA - offset),
               rho_s*math.sin(math.pi/2 - DM - PHI0 - offset), f)

# ---------------------------------------------------------------- viscous budget, d and theta
def viscous_budget_theta(L, c, C_K, f, offset=0.0, vartheta=0.0, phi0=PHI0, s=T3.S_REC,
                         delta=DELTA, c_G=None, Cprime=0.0, use_V1_for_Z=True):
    """t3_budget.viscous_budget with d exposed and eps_v assembled by Corollary V.5'."""
    if c_G is None:
        c_G = (1.5 + Cprime/L)*c
    n = T3.RECORD["n_dim"]
    r0 = (1.0 + f)*math.sin(phi0)
    nutau = c/(s*s*L)
    sy = nutau*T3.J_pow(c, -2.0)/c
    sz = nutau*T3.J_pow(c, 4.0)/c
    d = max(d_of(f, offset), 1e-12)
    R_minus = r0 - d
    N = r0/math.sin(delta)              # ||eta_0||_inf r_*/M ; the kinked 1/sin delta is used,
                                        # which is 7.6613 against the datum's true 7.5523:
                                        # conservative, and left as t3_budget has it
    eps_bulk = sy/(r0*r0)
    V1 = n*nutau*(math.exp(2*c_G) - 1.0)/c_G
    tau = c/L
    EW = 0.5*C_K*math.exp(c_G)*tau*V1
    if R_minus <= 0 or d <= 1e-11 or vartheta >= 1.0:
        return {"eps_v": float("inf"), "d": d, "R_minus": R_minus, "admissible": False}
    E_hess = EW*(r0/(R_minus**2)) + 4.0*N*EW/d
    trC = 8*sy + 2*sz
    trC2 = 4*(2*sy)**2 + (2*sz)**2
    trC3 = 4*(2*sy)**3 + (2*sz)**3
    E4 = (r0/(R_minus**5))*(trC**2 + 2*trC2)
    P_a_half = T3.gauss_tail_aniso(d/2.0, sy, sz)
    P_a_full = T3.gauss_tail_aniso(d, sy, sz)
    if use_V1_for_Z:
        q = (d/2.0)**2/nutau
        P_Z_half = min(1.0, 2*n*math.exp(-math.exp(-2*c_G)*q/(4.0*n)))
    else:
        P_Z_half = P_a_half
    Pp = P_Z_half + P_a_half + P_a_full
    m2 = trC
    m4 = trC**2 + 2*trC2
    m6 = trC**3 + 6*trC*trC2 + 8*trC3
    E_tail = 2*N*Pp + Pp + math.sqrt(Pp)*(math.sqrt(m2)/r0 + math.sqrt(m4)/r0**2
                                          + math.sqrt(m6)/r0**3)
    E_tail = min(E_tail, 1e12)
    Q = (1.0 + vartheta)/(1.0 - vartheta)
    vt = vartheta/(1.0 - vartheta)
    ratio = sz/sy
    E_theta = eps_bulk*vt*(9.0 + 2.0*ratio)
    eps_v = eps_bulk + E_theta + Q*(E_hess + E4 + E_tail)
    return {"eps_bulk": eps_bulk, "E_hess": E_hess, "E_4": E4, "E_tail": E_tail,
            "E_vartheta": E_theta, "Q": Q, "vt": vt, "sigma_z_over_sigma_y": ratio,
            "d": d, "R_minus": R_minus, "N": N, "sigma_y": sy, "sigma_z": sz, "c_G": c_G,
            "vartheta": vartheta, "admissible": True,
            "eps_v_theta0_stock_form": eps_bulk + E_hess + E4 + E_tail,
            "eps_v": min(eps_v, 1e12)}

def install(offset, safety):
    """monkey-patch t3_budget.viscous_budget inside fix4's p4, per f."""
    def vb(L, c, C_K, f, **kw):
        for k in ("phi0", "s", "delta"):
            kw.pop(k, None)
        th = safety*VT.get(f, 1.0) if safety > 0 else 0.0
        return viscous_budget_theta(L, c, C_K, f, offset=offset, vartheta=th, **kw)
    T3.viscous_budget = vb
    P4._CC.clear(); P4._BS.clear()

_ORIG_VB = T3.viscous_budget

def restore():
    T3.viscous_budget = _ORIG_VB
    P4._CC.clear(); P4._BS.clear()

if __name__ == "__main__":
    P1RES = json.load(open(os.path.join(FIX4, "p1_results.json")))
    ROW = P1RES["rows"]["sigma=0.02"]
    E0, G0 = ROW["E0"], ROW["Gfrak0"]
    CAP_FIX4 = 1.0 + ROW["eps_cap_certified"]
    D = P4.Datum(SIGMA, SHIFT)
    P4.bind(D, E0, G0)
    OUT = {"vartheta_at_corrected_ball": VT, "SAFETY": SAFETY,
           "offset": 4.0*SIGMA, "offset_over_sigma": 4.0,
           "E0_fix4": E0, "Gfrak0_fix4": G0, "cap_fix4": CAP_FIX4,
           "d_corrected": {"f=%g" % f: d_of(f, 4.0*SIGMA) for f in P4.FS_F4},
           "d_stock": {"f=%g" % f: d_of(f, 0.0) for f in P4.FS_F4}}

    # ---- transcription identity check: at stock d and vartheta = 0 the two must agree exactly
    ident = []
    for (L, c, f) in [(1e6, 0.88, 4.0), (1e4, 0.83, 8.0), (2e6, 0.879, 4.0), (1e7, 0.845, 16.0)]:
        a = _ORIG_VB(L, c, 161.7735, f, c_G=1.24*c/0.83)
        b = viscous_budget_theta(L, c, 161.7735, f, offset=0.0, vartheta=0.0, c_G=1.24*c/0.83)
        ident.append({"L": L, "c": c, "f": f, "eps_v_imported": a["eps_v"],
                      "eps_v_transcribed": b["eps_v"],
                      "rel": abs(a["eps_v"] - b["eps_v"])/max(a["eps_v"], 1e-300)})
    OUT["transcription_identity_check"] = ident
    print("identity (stock d, vartheta = 0):", ["%.2e" % r["rel"] for r in ident], flush=True)

    # ---- CONTROL: reproduce fix4's L_* with the stock instrument
    restore()
    ctl = {}
    for tgt, tag in [(CAP_FIX4 - 1.0, "eps<=cap"), (0.1, "eps<=0.1")]:
        ctl[tag] = P4.Lstar(D, "proved", tgt, E0, G0, CAP_FIX4, fs=P4.FS_F4)
    ctl["fix4_quoted"] = {"eps<=cap": 1424610.4953, "eps<=0.1": 1555469.4004}
    OUT["control_reproduce_fix4"] = ctl
    print("CONTROL reproduce fix4: %s" % ctl, flush=True)

    # ---- eps_v and its share, three L, proved column, with the repair carried
    install(4.0*SIGMA, SAFETY)
    share = {}
    for L in [1e5, 1e6, 2e6]:
        sc = P4.self_consistent(D, L, "proved", E0, G0, CAP_FIX4, fs=P4.FS_F4)
        for tag, c in [("frozen c = c_*", D.c_star)] + ([("self-consistent", sc["c"])]
                                                        if sc["self_consistent"] else []):
            f_best = sc["best_f"] if (sc["self_consistent"] and tag == "self-consistent") else 4.0
            r = P4.assemble_c(D, L, c, "proved", E0, G0, f=f_best)
            vb = viscous_budget_theta(L, c, T3.RECORD["K2_proved_window"], f_best,
                                      offset=4.0*SIGMA, vartheta=SAFETY*VT[f_best],
                                      c_G=P4.column_constants(L, c, "proved", E0, G0)["c_G"])
            v0 = viscous_budget_theta(L, c, T3.RECORD["K2_proved_window"], f_best,
                                      offset=0.0, vartheta=0.0,
                                      c_G=P4.column_constants(L, c, "proved", E0, G0)["c_G"])
            share["L=%g, %s" % (L, tag)] = {
                "L": L, "c_over_c_star": c/D.c_star, "f": f_best,
                "self_consistent": bool(sc["self_consistent"]),
                "vartheta_used": SAFETY*VT[f_best], "Q": vb["Q"], "vt": vb["vt"],
                "sigma_z_over_sigma_y": vb["sigma_z_over_sigma_y"],
                "d_corrected": vb["d"], "d_stock": v0["d"],
                "eps_bulk": vb["eps_bulk"], "E_vartheta": vb["E_vartheta"],
                "E_hess": vb["E_hess"], "E_4": vb["E_4"], "E_tail": vb["E_tail"],
                "eps_v": vb["eps_v"], "eps_v_fix4_stock": v0["eps_v"],
                "ratio_to_fix4": vb["eps_v"]/max(v0["eps_v"], 1e-300),
                "eps": r["eps"], "eps_Tprime": r["eps_Tprime"], "eps_a": r["eps_a"],
                "share_eps_v_over_eps": (vb["eps_v"]/r["eps"]) if np.isfinite(r["eps"])
                                        and r["eps"] > 0 else None}
            print("  L=%.0e %-16s f=%g  eps_v = %.6e (fix4 %.6e, x%.4f)  eps = %s  share = %s"
                  % (L, tag, f_best, vb["eps_v"], v0["eps_v"],
                     vb["eps_v"]/max(v0["eps_v"], 1e-300), r["eps"],
                     share["L=%g, %s" % (L, tag)]["share_eps_v_over_eps"]), flush=True)
    OUT["eps_v_and_share"] = share
    json.dump(OUT, open(os.path.join(HERE, "x2_results.json"), "w"), indent=1, default=str)

    # ---- L_*, three columns, two targets, with the repair carried
    res = {}
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for tgt, tag in [(CAP_FIX4 - 1.0, "eps<=cap(fix4)"), (0.1, "eps<=0.1")]:
            v = P4.Lstar(D, col, tgt, E0, G0, CAP_FIX4, fs=P4.FS_F4)
            res["%s_%s" % (col, tag)] = {"L_star": v,
                                         "log_Lambda_star": (2*v + SHIFT) if v else None}
            print("  L_star REPAIRED %-34s %s" % (col + " " + tag, v), flush=True)
    OUT["L_star_repaired_fix4_constants"] = res
    json.dump(OUT, open(os.path.join(HERE, "x2_results.json"), "w"), indent=1, default=str)

    # ---- robustness: the ball at 6 sigma instead of 4 sigma
    VT6 = {4.0: 0.0884151, 8.0: 0.0884151, 16.0: 0.0884151, 32.0: 0.0884151}
    _sv = dict(VT)
    VT.update(VT6)
    install(6.0*SIGMA, SAFETY)
    r6 = {}
    for tgt, tag in [(CAP_FIX4 - 1.0, "eps<=cap(fix4)"), (0.1, "eps<=0.1")]:
        r6[tag] = P4.Lstar(D, "proved", tgt, E0, G0, CAP_FIX4, fs=P4.FS_F4)
        print("  L_star 6-sigma ball, proved, %-20s %s" % (tag, r6[tag]), flush=True)
    OUT["L_star_6sigma_ball_proved"] = r6
    VT.clear(); VT.update(_sv)

    # ---- UNIT F-6's last step: the same run on x6's CERTIFIED constants
    E0c = X6["E_0"]["upper"]
    G0c = X6["Gfrak_0"]["upper"]
    CAPc = 1.0 + X6["window_cap_certified_PROVED_pad"]["eps_cap"]
    install(4.0*SIGMA, SAFETY)
    P4.bind(D, E0c, G0c)
    # r_h with the PROVED Lipschitz pad
    import p1_profile as P1
    D._rh.clear()
    D.lip_pad = X6["h_second_bound"]["proved"]*D.prof.h
    D.rh_nsub = 20000
    cert = {}
    for tgt, tag in [(CAPc - 1.0, "eps<=certified cap"), (0.1, "eps<=0.1")]:
        for col in ["proved", "proved_Ca_measured", "measured"]:
            v = P4.Lstar(D, col, tgt, E0c, G0c, CAPc, fs=P4.FS_F4)
            cert["%s_%s" % (col, tag)] = {"L_star": v,
                                          "log_Lambda_star": (2*v + SHIFT) if v else None}
            print("  L_star CERTIFIED %-38s %s" % (col + " " + tag, v), flush=True)
    OUT["L_star_certified_constants"] = {
        "E0_used": E0c, "Gfrak0_used": G0c, "cap_used": CAPc,
        "lip_pad_used": D.lip_pad, "rh_nsub": D.rh_nsub, "results": cert}
    json.dump(OUT, open(os.path.join(HERE, "x2_results.json"), "w"), indent=1, default=str)
    restore()
    print("wrote x2_results.json")

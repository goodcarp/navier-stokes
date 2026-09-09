"""
f1 -- ITEM 1: TIME-WINDOW COVERAGE.

THE DEFECT.  THEOREM_S3 sec.1.3 assumes every bootstrap estimate on [0, tau], tau = c_*/(M L)
(theta = M L t in [0, c_*]), and concludes  T_d <= t_* = tau(1+eps) > tau.  The budget is
computed at c = c_* throughout (t3_budget.py:24, and eps_at() passes CSTAR).  So the estimates
are not established on the interval on which the conclusion is drawn.

THE REPAIR.  Run the whole budget on the window actually used, theta in [0, c], and impose
SELF-CONSISTENCY

        c  >=  c_* ( 1 + eps(c, L) )                                                  (SC)

with eps(c,L) the budget's own output AT that c.  Everything that depends on the window is
re-evaluated at c:

    lam_max = e^{3c/4}         (P3: the a priori cap over the window)
    r_h     = inf_{[1,lam_max]} P_h/lam                (f4's certified lower bound at that lam_max)
    ell_loss contains 3 log lam_max
    the majorant ODE integrates theta over [0, c]
    c_G     = Gbar c / L  in (Gamma-off)  (f5's direct root, not the Picard iteration)
    J_pow(c, .), the Gaussian widths sigma_y, sigma_z and the whole viscous budget

eps(c,L) is increasing in c, so the smallest admissible window is the fixed point of (SC), and
that fixed point's eps is the honest eps(L).  Reported alongside: the cruder repair the brief
also allows -- majorise by a slightly larger FIXED c = c_*(1+eps_target) and check that the
budget at that c returns eps <= eps_target (then the doubling is reached inside the window
where the estimates hold).

A STRUCTURAL CAP FALLS OUT (f4 sec.2b).  A2/C6 need P_h increasing on [1, lam_max] for the
quasimonotone comparison.  P_h' vanishes at lam_mono = 2.0769162, so the window may not exceed
c_mono = (4/3) log lam_mono = 1.19649 c_*, i.e. eps <= 0.19649 (certified: 0.19437).
THE eps <= 1/2 COLUMN OF THEOREM_S3 IS THEREFORE NOT AVAILABLE: it needs lam_max = 2.499991.

IMPORTS.  t3_budget.py, t2_gamma_CR.py and t1_results.json are byte copies in ./imported/
(hashes in SHA256SUMS); the viscous half (viscous_budget, gauss_tail_aniso, J_pow, eps_Tprime,
RECORD) is used UNCHANGED, exactly as THEOREM_S3 uses it.  bootstrap() is used unchanged too --
it already takes lam_mx as an argument.  Only the c-dependence that THEOREM_S3 froze is undone.

Outputs -> f1_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "imported"))
os.chdir(os.path.join(HERE, "imported"))          # t3_budget reads t1_results.json from CWD
import numpy as np
import t3_budget as T3                            # noqa: E402
import t2_gamma_CR as T2                          # noqa: E402
os.chdir(HERE)
sys.path.insert(0, HERE)
import f4_scans as F4                             # noqa: E402
import f5_gamma_exist as F5                       # noqa: E402

KAPPA = T3.KAPPA
CSTAR = T3.CSTAR
LAM_OM = T3.LAM_OM
E0 = 1.0/math.sin(7.5*math.pi/180)                # PROVED exactly (f4); t1's grid value 7.6612976
G0 = 1.0/math.sin(7.5*math.pi/180)**2             # PROVED exactly (f4); t1's grid value 58.695476
SHIFT = T3.SHIFT
SHIFT_SHARP = T3.SHIFT_SHARP
ELL_RAMP = T3.ELL_RAMP
FS_ALL = T3.FS
FS_F4 = [f for f in T3.FS if f >= 4.0]            # ADDENDUM_1: the theorem must be stated at f >= 4

_RH = {}
def r_h_of(lam_mx, nsub=5000):
    key = round(lam_mx, 10)
    if key not in _RH:
        _RH[key] = F4.certified_r_h(lam_mx, nsub=nsub)["r_h_certified_lower"]
    return _RH[key]

_CC = {}
def column_constants(L, c, column):
    key = (round(L, 6), round(c, 12), column)
    if key in _CC:
        return _CC[key]
    route, ca_kind, gh_kind, k2key, sens, use_V1 = T3.COLUMNS[column]
    if route == "u2_fixedpoint":
        fp = F5.smallest_root(L, c, LAM_OM, E0, G0, sigma_star=0.0)   # DIRECT root (item 5)
        if fp is None:
            _CC[key] = None
            return None
        Cpp, c_G, Ghat = fp["C_pp"], fp["c_G"], fp["Ghat"]
    else:
        Cpp = T3.RECORD["Cprime_measured"]
        c_G = (1.5 + Cpp/L)*c
        Ghat = T3.RECORD["Ghat_measured"]
    if gh_kind == "measured":
        Ghat = T3.RECORD["Ghat_measured"]
    ca_flat = T3.RECORD["material_offset_M"] if ca_kind == "offset" else None
    C_R, argphi = T2.C_R_sup(LAM_OM, c_G, E0, Ghat, ca_flat=ca_flat,
                             grad_sinphi=(gh_kind == "proved"), phi_lo_deg=0.0, n=400)
    ca_sup = ca_flat if ca_flat is not None else T2.chat_a(0.0, LAM_OM, c_G, E0)
    ca_traj = ca_flat if ca_flat is not None else T2.chat_a(math.sin(30*math.pi/180),
                                                            LAM_OM, c_G, E0)
    out = {"C_pp": Cpp, "c_G": c_G, "Ghat": Ghat, "C_R": C_R, "C_R_argmax_phi_deg": argphi,
           "ca_sup": ca_sup, "ca_traj": ca_traj, "C_K": T3.RECORD[k2key], "sens": sens,
           "use_V1": use_V1}
    _CC[key] = out
    return out

def ell_loss_c(f, mu, lam_mx):
    return (math.log(2.0) + math.log(1.0+f) + 3.0*math.log(lam_mx)
            + math.log((1.0+mu)/max(1.0-mu, 1e-12)) + ELL_RAMP)

_BS = {}
def bootstrap_cached(L, c, column, C_R_override=None):
    key = (round(L, 6), round(c, 12), column, C_R_override)
    if key not in _BS:
        K = column_constants(L, c, column)
        lam_mx = math.exp(0.75*c)
        C_R = K["C_R"] if C_R_override is None else C_R_override
        _BS[key] = (T3.bootstrap(L, c, K["C_pp"], C_R, K["ca_sup"], C1=0.0, lam_mx=lam_mx), C_R)
    return _BS[key]


def assemble_c(L, c, column, f=1.0, C_K=None, C_R_override=None):
    """THEOREM_S3's assemble(), with lam_max, r_h and ell_loss re-read at the window's own c."""
    K = column_constants(L, c, column)
    if K is None:
        return {"L": L, "c": c, "column": column, "f": f, "eps": float("inf"),
                "gamma_off_exists": False}
    lam_mx = math.exp(0.75*c)
    r_h = r_h_of(lam_mx)
    bs, C_R = bootstrap_cached(L, c, column, C_R_override)
    CK = K["C_K"] if C_K is None else C_K
    vb = T3.viscous_budget(L, c, CK, f, c_G=K["c_G"], use_V1_for_Z=K["use_V1"])
    epsT = T3.eps_Tprime(bs["mu"], bs["muJ"], r_h, K["sens"])
    ell = ell_loss_c(f, min(bs["mu"], 0.999) if np.isfinite(bs["mu"]) else 0.999, lam_mx)
    eps_a = ell/L + epsT + K["ca_traj"]/(KAPPA*L)
    eps_v = vb["eps_v"]
    eps_delta = 1.0 - 2.0*KAPPA
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        eps = float("inf")
    else:
        eps = (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-eps_delta)) - 1.0
    return {"L": L, "c": c, "column": column, "f": f, "gamma_off_exists": True,
            "lam_max": lam_mx, "r_h": r_h, "eps_ell": ell/L, "eps_Tprime": epsT,
            "eps_Ca": K["ca_traj"]/(KAPPA*L), "eps_a": eps_a, "eps_v": eps_v,
            "eps_delta": eps_delta, "mu": bs["mu"], "muJ": bs["muJ"], "c_G": K["c_G"],
            "C_pp": K["C_pp"], "Ghat": K["Ghat"], "C_R": C_R, "eps": eps,
            "E_tail": vb["E_tail"], "eps_bulk": vb["eps_bulk"], "E_hess": vb["E_hess"],
            "E_4": vb["E_4"]}

def eps_at(L, c, col, fs=None, **kw):
    fs = FS_ALL if fs is None else fs
    best, arg = float("inf"), None
    for f in fs:
        r = assemble_c(L, c, col, f=f, **kw)
        if r["eps"] < best:
            best, arg = r["eps"], f
    return best, arg

# ------------------------------------------------------------------ the self-consistent window
CAP_STRUCT = 1.1964878432367183      # c/c_* at lam_mono (f4)
CAP_CERT = 1.1943662078602755        # c/c_* at the largest CERTIFIED P_h-monotone window (f4)

def self_consistent(L, col, fs=None, cap=CAP_CERT, n=20, nbis=20):
    """smallest c in [c_*, cap*c_*] with c >= c_*(1+eps(c,L)); returns eps = c/c_*-1 there."""
    los, his = 1.0, cap
    grid = [los + (his-los)*i/n for i in range(n+1)]
    ok = None
    for x in grid:
        c = x*CSTAR
        e, f = eps_at(L, c, col, fs=fs)
        if math.isfinite(e) and CSTAR*(1.0+e) <= c:
            ok = x
            break
    if ok is None:
        return {"eps": float("inf"), "c_over_c_star": None, "best_f": None,
                "self_consistent": False}
    lo, hi = ok - (his-los)/n, ok
    for _ in range(nbis):
        mid = 0.5*(lo+hi)
        c = mid*CSTAR
        e, f = eps_at(L, c, col, fs=fs)
        if math.isfinite(e) and CSTAR*(1.0+e) <= c:
            hi = mid
        else:
            lo = mid
    c = hi*CSTAR
    e, f = eps_at(L, c, col, fs=fs)
    return {"eps": hi - 1.0, "eps_budget_at_that_c": e, "c_over_c_star": hi, "c": c,
            "lam_max": math.exp(0.75*c), "best_f": f, "self_consistent": True}

def Lstar_self_consistent(col, target, fs=None, cap=CAP_CERT):
    lo, hi = 2.0, 1e12
    if self_consistent(hi, col, fs=fs, cap=cap)["eps"] > target:
        return None
    for _ in range(45):
        mid = math.sqrt(lo*hi)
        if self_consistent(mid, col, fs=fs, cap=cap)["eps"] <= target:
            hi = mid
        else:
            lo = mid
        if hi/lo < 1.0000001:
            break
    return hi

def Lstar_fixed_c(col, eps_target, fs=None):
    """the cruder repair: fixed enlarged window c = c_*(1+eps_target); require eps(c) <= target."""
    c = CSTAR*(1.0 + eps_target)
    lo, hi = 2.0, 1e12
    if eps_at(hi, c, col, fs=fs)[0] > eps_target:
        return None, c
    for _ in range(45):
        mid = math.sqrt(lo*hi)
        if eps_at(mid, c, col, fs=fs)[0] <= eps_target:
            hi = mid
        else:
            lo = mid
        if hi/lo < 1.00001:
            break
    return hi, c


def _save(O):
    with open(os.path.join(HERE, "f1_results.json"), "w") as fh:
        json.dump(O, fh, indent=1, sort_keys=True, default=str)


if __name__ == "__main__":
    OUT = {"inputs": {"kappa_delta": KAPPA, "c_star": CSTAR, "E0_exact": E0, "G0_exact": G0,
                      "shift_logReE": SHIFT, "shift_sharp": SHIFT_SHARP,
                      "cap_struct_c_over_cstar": CAP_STRUCT, "cap_cert_c_over_cstar": CAP_CERT,
                      "lam_max_at_c_star": math.exp(0.75*CSTAR),
                      "lam_max_at_cap_cert": math.exp(0.75*CAP_CERT*CSTAR)}}

    # ---- control: at c = c_* with the frozen window, reproduce THEOREM_S3's eps(L) table
    ctl = {}
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for L in [1e3, 1e4, 1e5, 1e6]:
            e, f = eps_at(L, CSTAR, col)
            ctl["%s_L=%g" % (col, L)] = {"eps_frozen_window": e, "best_f": f}
    ctl["THEOREM_S3_quoted"] = {"proved_L=1e6": 0.10622, "proved_Ca_measured_L=1e6": 0.10276,
                                "measured_L=1e3": 0.11572, "measured_L=1e4": 0.012553,
                                "measured_L=1e5": 5.1682e-3, "measured_L=1e6": 4.4533e-3}
    OUT["control_frozen_window_c_eq_c_star"] = ctl
    _save(OUT)

    # ---- the corrected eps(L): self-consistent window
    sc = {}
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for L in [1e3, 1e4, 1e5, 1e6, 1e7]:
            sc["%s_L=%g" % (col, L)] = self_consistent(L, col)
    OUT["self_consistent_eps"] = sc
    _save(OUT)

    # ---- L_* , three columns, corrected
    ls, lam = {}, {}
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for tgt, tag in [(0.1, "eps<=0.1"), (CAP_CERT-1.0, "eps<=eps_cap=0.194366")]:
            k = "%s_%s" % (col, tag)
            v = Lstar_self_consistent(col, tgt)
            ls[k] = v
            lam[k] = {"log_Lambda_star": (2*v + SHIFT) if v else None}
            print("  L_star", k, v, flush=True)
            OUT["L_star_self_consistent"] = ls
            OUT["log_Lambda_star_self_consistent"] = lam
            _save(OUT)
    OUT["L_star_self_consistent"] = ls
    OUT["log_Lambda_star_self_consistent"] = lam

    # ---- the same with ADDENDUM_1's f >= 4 constraint
    ls4, lam4, sc4 = {}, {}, {}
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for L in [1e5, 1e6, 1e7]:
            sc4["%s_L=%g" % (col, L)] = self_consistent(L, col, fs=FS_F4)
        for tgt, tag in [(0.1, "eps<=0.1"), (CAP_CERT-1.0, "eps<=eps_cap=0.194366")]:
            k = "%s_%s" % (col, tag)
            v = Lstar_self_consistent(col, tgt, fs=FS_F4)
            ls4[k] = v
            lam4[k] = {"log_Lambda_star": (2*v + SHIFT) if v else None}
            print("  L_star f>=4", k, v, flush=True)
            OUT["self_consistent_eps_f_ge_4"] = sc4
            OUT["L_star_self_consistent_f_ge_4"] = ls4
            OUT["log_Lambda_star_f_ge_4"] = lam4
            _save(OUT)
    OUT["self_consistent_eps_f_ge_4"] = sc4
    OUT["L_star_self_consistent_f_ge_4"] = ls4
    OUT["log_Lambda_star_f_ge_4"] = lam4

    # ---- the cruder repair: a fixed enlarged window
    fx = {}
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for tgt in [0.1, CAP_CERT-1.0]:
            v, c = Lstar_fixed_c(col, tgt)
            fx["%s_eps<=%.6f" % (col, tgt)] = {"L_star": v, "c": c, "c_over_c_star": c/CSTAR,
                                               "lam_max": math.exp(0.75*c),
                                               "log_Lambda_star": (2*v + SHIFT) if v else None}
    OUT["L_star_fixed_enlarged_window"] = fx

    # ---- what the OLD eps <= 1/2 column would have needed, and why it is unavailable
    OUT["eps_half_column_is_unavailable"] = {
        "c_needed_over_c_star": 1.5, "lam_max_needed": math.exp(0.75*1.5*CSTAR),
        "lam_mono": 2.0769162051407735,
        "reason": "P_h decreases beyond lam_mono, so A2/C6's quasimonotonicity fails on "
                  "[1, lam_max]; the eps<=1/2 window is not admissible at all.",
        "THEOREM_S3_L_star_eps_half_proved": 311511.9}

    # ---- Gamma-off thresholds at the corrected windows (item 5 feeding item 1)
    OUT["L_Gamma_exist_at_windows"] = {}
    for x in [1.0, 1.05, 1.1, CAP_CERT]:
        c = x*CSTAR
        OUT["L_Gamma_exist_at_windows"]["c=%.6f c_*" % x] = {
            "L_exist": F5.L_exist(c, LAM_OM, E0, G0, 0.0),
            "L_picard": F5.L_picard(c, LAM_OM, E0, G0, 0.0)}

    with open(os.path.join(HERE, "f1_results.json"), "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
    print(json.dumps(OUT, indent=1, sort_keys=True, default=str))

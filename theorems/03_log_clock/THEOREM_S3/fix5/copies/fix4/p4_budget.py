"""
p4 -- THE BUDGET for the sigma-mollified datum: eps(L), L_*, log Lambda_*, three columns.

PROVENANCE.  This is fix2's `f1_window.py` (self-consistent window, item 1 of FIX2) with ONE
substitution: every datum constant is the MOLLIFIED datum's, from p1.

    kappa_delta, c_*, lam_max = e^{3c/4}, r_h on [1, lam_max], E_0, Gfrak_0, and the
    log Re_E shift (p2) are rebound; t3_budget's viscous half (viscous_budget,
    gauss_tail_aniso, J_pow, eps_Tprime, RECORD) and its bootstrap() are used UNCHANGED,
    and t2_gamma_CR's chat_a / C_R_sup and f5_gamma_exist's DIRECT (Gamma-off) root are
    used unchanged.  t3_budget's module-level KAPPA, CSTAR, LAM_MAX, R_H, E0, G0 and SHIFT
    are rebound to the mollified values before any call: that rebinding IS the substitution
    under test (L-14's inverse convention, declared).

The self-consistency condition is fix2's, unchanged:

        (SC)     c  >=  c_* ( 1 + eps(c, L) ) ,        lam_max = e^{3c/4} ,

and the admissible window is capped by the largest lam on which P_h is certified increasing
(p1's `eps_cap_certified`, which MOVES with sigma).

Outputs -> p4_results.json  (or p4_scan.json in --scan mode)
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "imported"))
os.chdir(os.path.join(HERE, "imported"))          # t3_budget reads t1_results.json from CWD
import t3_budget as T3                            # noqa: E402
import t2_gamma_CR as T2                          # noqa: E402
import f5_gamma_exist as F5                       # noqa: E402
os.chdir(HERE)
sys.path.insert(0, HERE)
import p1_profile as P1                           # noqa: E402

LAM_OM = T3.LAM_OM
FS_ALL = T3.FS
FS_F4 = [f for f in T3.FS if f >= 4.0]            # ADDENDUM_1: the theorem is stated at f >= 4
SHIFT_KINKED = T3.SHIFT                           # 2.9234859, kinked datum (control only)


class Datum:
    """the mollified datum's constants, and the r_h oracle at any window."""
    def __init__(self, sigma, shift, nphi=600001, rh_nsub=5000):
        self.sigma = sigma
        self.prof, self.inst = P1.build(sigma, nphi=nphi)
        self.kappa = 0.5*self.inst.P_h(1.0)
        self.c_star = math.log(1.5)/self.kappa
        self.E0 = P1.datum_norms(self.prof, nu=401)["E0_scan"] if False else None
        self.rh_nsub = rh_nsub
        self.lip_pad = P1.lipschitz_hpp(self.prof)*(math.pi/(nphi-1))
        self.shift = shift
        self._rh = {}

    def r_h(self, lam_mx):
        k = round(lam_mx, 9)
        if k not in self._rh:
            self._rh[k] = P1.certified_r_h(self.inst, self.prof, lam_mx,
                                           nsub=self.rh_nsub, lip_pad=self.lip_pad)
        return self._rh[k]["r_h_certified_lower"]


def bind(D, E0, G0):
    T3.KAPPA = D.kappa
    T3.CSTAR = D.c_star
    T3.C2 = 2.0*math.log(1.5)/D.kappa
    T3.E0 = E0
    T3.G0 = G0
    T3.SHIFT = D.shift


def C_R_sup_fast(lam_om, c_G, E0, Ghat, ca_flat=None, grad_sinphi=True, n=160,
                 nq_c=2001, nq_f=20001):
    """T2.C_R_of_phi (the instrument, UNCHANGED) with a cheaper search: a coarse scan then a
       41-point refinement of the bracket.  Controlled against T2.C_R_sup in the gate."""
    phis = np.linspace(1e-4, math.pi/2.0, n)
    vals = [T2.C_R_of_phi(p_, lam_om, c_G, E0, Ghat, ca_flat, grad_sinphi, nq=nq_c)
            for p_ in phis]
    k = int(np.argmax(vals))
    a = phis[max(k-1, 0)]
    b = phis[min(k+1, n-1)]
    ph2 = np.linspace(a, b, 41)
    v2 = [T2.C_R_of_phi(p_, lam_om, c_G, E0, Ghat, ca_flat, grad_sinphi, nq=nq_f)
          for p_ in ph2]
    k2 = int(np.argmax(v2))
    return float(v2[k2]), float(ph2[k2]*180/math.pi)


_CC = {}
def column_constants(L, c, column, E0, G0, nphi_CR=160):
    key = (round(L, 6), round(c, 12), column)
    if key in _CC:
        return _CC[key]
    route, ca_kind, gh_kind, k2key, sens, use_V1 = T3.COLUMNS[column]
    if route == "u2_fixedpoint":
        fp = F5.smallest_root(L, c, LAM_OM, E0, G0, sigma_star=0.0)
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
    C_R, argphi = C_R_sup_fast(LAM_OM, c_G, E0, Ghat, ca_flat=ca_flat,
                               grad_sinphi=(gh_kind == "proved"), n=nphi_CR)
    ca_sup = ca_flat if ca_flat is not None else T2.chat_a(0.0, LAM_OM, c_G, E0)
    ca_traj = ca_flat if ca_flat is not None else T2.chat_a(math.sin(30*math.pi/180),
                                                            LAM_OM, c_G, E0)
    out = {"C_pp": Cpp, "c_G": c_G, "Ghat": Ghat, "C_R": C_R, "C_R_argmax_phi_deg": argphi,
           "ca_sup": ca_sup, "ca_traj": ca_traj, "C_K": T3.RECORD[k2key], "sens": sens,
           "use_V1": use_V1}
    _CC[key] = out
    return out


def ell_loss_c(f, mu, lam_mx, ell_ramp):
    return (math.log(2.0) + math.log(1.0+f) + 3.0*math.log(lam_mx)
            + math.log((1.0+mu)/max(1.0-mu, 1e-12)) + ell_ramp)


_BS = {}
def assemble_c(D, L, c, column, E0, G0, f=1.0, C_K=None, C_R_override=None,
               ell_ramp=0.5, eps_a_extra=0.0, nphi_CR=160):
    K = column_constants(L, c, column, E0, G0, nphi_CR=nphi_CR)
    if K is None:
        return {"L": L, "c": c, "column": column, "f": f, "eps": float("inf"),
                "gamma_off_exists": False}
    lam_mx = math.exp(0.75*c)
    r_h = D.r_h(lam_mx)
    C_R = K["C_R"] if C_R_override is None else C_R_override
    bkey = (round(L, 6), round(c, 12), column, C_R_override)
    if bkey not in _BS:
        _BS[bkey] = T3.bootstrap(L, c, K["C_pp"], C_R, K["ca_sup"], C1=0.0, lam_mx=lam_mx)
    bs = _BS[bkey]
    CK = K["C_K"] if C_K is None else C_K
    vb = T3.viscous_budget(L, c, CK, f, c_G=K["c_G"], use_V1_for_Z=K["use_V1"])
    epsT = T3.eps_Tprime(bs["mu"], bs["muJ"], r_h, K["sens"])
    ell = ell_loss_c(f, min(bs["mu"], 0.999) if np.isfinite(bs["mu"]) else 0.999, lam_mx,
                     ell_ramp)
    eps_a = ell/L + epsT + K["ca_traj"]/(D.kappa*L) + eps_a_extra
    eps_v = vb["eps_v"]
    eps_delta = 1.0 - 2.0*D.kappa
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        eps = float("inf")
    else:
        eps = (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-eps_delta)) - 1.0
    return {"L": L, "c": c, "column": column, "f": f, "gamma_off_exists": True,
            "lam_max": lam_mx, "r_h": r_h, "eps_ell": ell/L, "eps_Tprime": epsT,
            "eps_Ca": K["ca_traj"]/(D.kappa*L), "eps_a": eps_a, "eps_v": eps_v,
            "eps_delta": eps_delta, "mu": bs["mu"], "muJ": bs["muJ"], "c_G": K["c_G"],
            "C_pp": K["C_pp"], "Ghat": K["Ghat"], "C_R": C_R, "eps": eps,
            "E_tail": vb["E_tail"], "eps_bulk": vb["eps_bulk"], "E_hess": vb["E_hess"],
            "E_4": vb["E_4"], "ca_traj": K["ca_traj"], "ca_sup": K["ca_sup"], "C_K": CK}


def eps_at(D, L, c, col, E0, G0, fs=None, **kw):
    fs = FS_ALL if fs is None else fs
    best, arg = float("inf"), None
    for f in fs:
        r = assemble_c(D, L, c, col, E0, G0, f=f, **kw)
        if r["eps"] < best:
            best, arg = r["eps"], f
    return best, arg


def self_consistent(D, L, col, E0, G0, cap, fs=None, n=16, nbis=16, **kw):
    los, his = 1.0, cap
    ok = None
    for i in range(n+1):
        x = los + (his-los)*i/n
        e, _ = eps_at(D, L, x*D.c_star, col, E0, G0, fs=fs, **kw)
        if math.isfinite(e) and D.c_star*(1.0+e) <= x*D.c_star:
            ok = x
            break
    if ok is None:
        return {"eps": float("inf"), "c_over_c_star": None, "best_f": None,
                "self_consistent": False}
    lo, hi = ok - (his-los)/n, ok
    for _ in range(nbis):
        mid = 0.5*(lo+hi)
        e, _ = eps_at(D, L, mid*D.c_star, col, E0, G0, fs=fs, **kw)
        if math.isfinite(e) and (1.0+e) <= mid:
            hi = mid
        else:
            lo = mid
    e, f = eps_at(D, L, hi*D.c_star, col, E0, G0, fs=fs, **kw)
    return {"eps": hi - 1.0, "eps_budget_at_that_c": e, "c_over_c_star": hi,
            "c": hi*D.c_star, "lam_max": math.exp(0.75*hi*D.c_star), "best_f": f,
            "self_consistent": True}


def Lstar(D, col, target, E0, G0, cap, fs=None, tol=1e-6, **kw):
    lo, hi = 2.0, 1e12
    if self_consistent(D, hi, col, E0, G0, cap, fs=fs, **kw)["eps"] > target:
        return None
    for _ in range(60):
        mid = math.sqrt(lo*hi)
        if self_consistent(D, mid, col, E0, G0, cap, fs=fs, **kw)["eps"] <= target:
            hi = mid
        else:
            lo = mid
        if hi/lo < 1.0 + tol:
            break
    return hi


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--sigma", type=float, default=None)
    ap.add_argument("--shift", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    P1RES = json.load(open(os.path.join(HERE, "p1_results.json")))

    def datum_for(sg, shift):
        key = "kinked" if sg is None else ("sigma=%g" % sg)
        row = P1RES["rows"][key]
        D = Datum(sg, shift)
        E0 = row["E0"]
        G0 = row["Gfrak0"]
        cap = 1.0 + row["eps_cap_certified"]
        bind(D, E0, G0)
        return D, E0, G0, cap, row

    if a.scan:
        OUT = {"mode": "scan", "note": "coarse settings (n=10,nbis=10,tol=1e-3,nphi_CR=120); "
                                       "used only to choose sigma, never quoted as a constant"}
        rows = {}
        for sg in [None] + P1.SIGMAS:
            D, E0, G0, cap, row = datum_for(sg, SHIFT_KINKED)
            _CC.clear(); _BS.clear()
            r = {"sigma": sg, "kappa": D.kappa, "c_star": D.c_star, "E0": E0, "G0": G0,
                 "eps_floor": row["eps_floor"], "cap": cap}
            for tgt, tag in [(0.1, "eps<=0.1"), (cap-1.0, "eps<=cap")]:
                if tgt <= row["eps_floor"]:
                    r[tag] = None
                    continue
                v = Lstar(D, "proved", tgt, E0, G0, cap, fs=FS_F4, tol=1e-3,
                          nphi_CR=120)
                r[tag] = v
                print("  sigma=%-7s %-12s L_* = %s" % (str(sg), tag, v), flush=True)
            rows["sigma=%s" % sg] = r
            OUT["rows"] = rows
            json.dump(OUT, open(os.path.join(HERE, "p4_scan.json"), "w"), indent=1,
                      sort_keys=True, default=str)
        print(json.dumps(OUT, indent=1, default=str))
        sys.exit(0)

    SG = a.sigma
    SHIFT = a.shift
    D, E0, G0, CAP, row = datum_for(SG, SHIFT)
    OUT = {"sigma": SG, "kappa_delta": D.kappa, "c_star": D.c_star,
           "c2": 2.0*math.log(1.5)/D.kappa, "E0": E0, "Gfrak0": G0,
           "eps_floor": row["eps_floor"], "cap_c_over_c_star": CAP,
           "eps_cap_certified": CAP - 1.0, "shift_logReE": SHIFT,
           "shift_logReE_kinked": SHIFT_KINKED, "ELL_RAMP_used": 0.5,
           "ELL_RAMP_in_t3_budget": T3.ELL_RAMP,
           "lam_max_at_c_star": math.exp(0.75*D.c_star),
           "lam_max_at_cap": math.exp(0.75*CAP*D.c_star)}

    # ---- eps(L), self-consistent window, f >= 4
    sc = {}
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for L in [1e3, 1e4, 1e5, 1e6, 1e7]:
            sc["%s_L=%g" % (col, L)] = self_consistent(D, L, col, E0, G0, CAP, fs=FS_F4)
            print("  eps  %-20s L=%.0e -> %s" % (col, L, sc["%s_L=%g" % (col, L)]["eps"]),
                  flush=True)
    OUT["eps_self_consistent_f_ge_4"] = sc
    json.dump(OUT, open(os.path.join(HERE, a.out or "p4_results.json"), "w"), indent=1,
              sort_keys=True, default=str)

    # ---- L_* and log Lambda_*
    ls, lam = {}, {}
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for tgt, tag in [(CAP-1.0, "eps<=cap"), (0.1, "eps<=0.1")]:
            k = "%s_%s" % (col, tag)
            v = Lstar(D, col, tgt, E0, G0, CAP, fs=FS_F4)
            ls[k] = v
            lam[k] = (2*v + SHIFT) if v else None
            print("  L_star %-34s %s   logLam=%s" % (k, v, lam[k]), flush=True)
            OUT["L_star"] = ls
            OUT["log_Lambda_star"] = lam
            json.dump(OUT, open(os.path.join(HERE, a.out or "p4_results.json"), "w"),
                      indent=1, sort_keys=True, default=str)
    OUT["L_star"] = ls
    OUT["log_Lambda_star"] = lam

    # ---- the term table at the self-consistent window, f = 4
    tt = []
    for col in ["proved", "proved_Ca_measured", "measured"]:
        for L in [1e3, 1e4, 1e5, 1e6]:
            s = self_consistent(D, L, col, E0, G0, CAP, fs=FS_F4)
            if not s["self_consistent"]:
                tt.append({"L": L, "column": col, "eps": float("inf")})
                continue
            tt.append(assemble_c(D, L, s["c"], col, E0, G0, f=s["best_f"]))
    OUT["term_table"] = tt

    # ---- constants against L in the proved column, at the frozen window c = c_*
    ct = {}
    for L in [1e4, 3e4, 1e5, 3e5, 1e6, 1e9]:
        K = column_constants(L, D.c_star, "proved", E0, G0)
        ct["L=%g" % L] = K
    OUT["constants_vs_L_frozen_window"] = ct
    OUT["L_Gamma_exist"] = {
        "c=c_star": F5.L_exist(D.c_star, LAM_OM, E0, G0, 0.0),
        "c=cap": F5.L_exist(CAP*D.c_star, LAM_OM, E0, G0, 0.0),
        "picard_c=c_star": F5.L_picard(D.c_star, LAM_OM, E0, G0, 0.0)
        if hasattr(F5, "L_picard") else None}

    # ---- sensitivity: what the ELL_RAMP correction and eps_a cost
    sens = {}
    for tag, kw in [("ELL_RAMP=0.25 (t3's value)", dict(ell_ramp=0.25)),
                    ("ELL_RAMP=0.5 (corrected)", dict(ell_ramp=0.5))]:
        sens[tag] = Lstar(D, "proved", 0.1, E0, G0, CAP, fs=FS_F4, tol=1e-5, **kw)
    for ck, tag in [(1.0, "C_K=1"), (161.7735, "C_K=hk2 proved window")]:
        sens[tag] = Lstar(D, "proved", 0.1, E0, G0, CAP, fs=FS_F4, tol=1e-5, C_K=ck)
    OUT["sensitivity"] = sens

    # ---- CONTROL (L-11 / L-14): the same instrument on the KINKED datum, at fix2's own
    # window cap and fix2's ELL_RAMP, must reproduce fix2's L_* = 1887786.7 / 1977309.3
    _CC.clear(); _BS.clear()
    Dk, E0k, G0k, CAPk, rowk = datum_for(None, SHIFT_KINKED)
    CAP_FIX2 = 1.1943662078602755
    ctl = {}
    for tgt, tag, cap in [(CAP_FIX2-1.0, "eps<=0.1943662 (fix2 cap)", CAP_FIX2),
                          (0.1, "eps<=0.1", CAP_FIX2)]:
        v = Lstar(Dk, "proved", tgt, E0k, G0k, cap, fs=FS_F4, ell_ramp=0.25)
        ctl[tag] = v
        print("  CONTROL kinked %-28s L_* = %s" % (tag, v), flush=True)
    ctl["fix2_quoted"] = {"eps<=0.1943662": 1887786.7127, "eps<=0.1": 1977309.3012}
    ctl["r_h_at_c_star_here"] = Dk.r_h(math.exp(0.75*Dk.c_star))
    ctl["fix2_r_h_enclosure"] = [0.9186364112, 0.9186365334]
    ctl["kappa_here"] = Dk.kappa
    ctl["fix2_kappa"] = 0.4978224182
    OUT["control_kinked_vs_fix2"] = ctl
    _CC.clear(); _BS.clear()
    bind(D, E0, G0)

    json.dump(OUT, open(os.path.join(HERE, a.out or "p4_results.json"), "w"), indent=1,
              sort_keys=True, default=str)
    print(json.dumps({k: v for k, v in OUT.items() if k not in ("term_table",)},
                     indent=1, default=str))

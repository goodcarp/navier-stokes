"""
t3 -- THE BUDGET of the assembled theorem, with the surviving constants.

PROVENANCE.  This file is  s3close/round2/u1/u5_budget.py  (sha256
4c2cdbc32cb0866a894f6e1ba8853a2e5e81ba47677d3c6bea084cd40854c7b0, itself a copy of
s3close/assembly/a2_budget.py) with the following replaced, and NOTHING else:

  * the datum constants        -- read from t1_results.json (computed in this folder for the
                                  theorem's own datum: ASSEMBLY sec.1.1 angular profile with
                                  hk2's tanh radial ramp), not from u1's u2_results.json.
                                  In particular kappa_delta is the FULL profile's (refute-u1 F4,
                                  u2 sec.9.1), not the axis taper's.
  * bootstrap()                -- C'' from U2's (Gamma-off) fixed point in place of C' = 151.15,
                                  C_1 = 0 in place of C_1 = 1 (U2 sec.7(a)), the drive written
                                  against M rather than M_s, and C_R from t2 (U2's proved
                                  r|grad a| in place of U1's unproved G_collar, the geometry of
                                  refute-u1 F3 kept, the supremum over EVERY angle).
  * COLUMNS, assemble()        -- the three columns the brief asks for.

CARRIED OVER UNCHANGED from u5_budget.py / a2_budget.py, so that the viscous half of the
budget is literally the assembly's own instrument and the comparison is apples to apples:
    gauss_tail_aniso(), viscous_budget(), J_pow(), eps_Tprime(), and the RECORD block.

Window: tau = c/(M L), theta = M L t in [0, c], c = c_* = log(3/2)/kappa_delta.
Outputs -> t3_results.json
"""
import json, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import ncx2

import t2_gamma_CR as T2

# --------------------------------------------------------------- record inputs (unchanged)
RECORD = {
    "lemmaT_abs": "|a[eta0.Phi^-1](0) - a_ref| <= pi M A [3 mu(1+muJ)/(2(1-mu)^5) + 3 muJ/8]",
    "Cprime_A1_provedV":   151.15,      # write/refute-L3v: the EXACT-PLATEAU route (superseded)
    "Cprime_measured":      11.74,
    "C_far_zodd": 0.291999,
    "C_inner_zodd": 0.014754,
    "C_collar_over_sinphi": 3.999218,
    "C_collar_axis": math.pi/8.0,
    "material_offset_M": 0.069,         # R9, MEASURED
    "K2_proved_window": 161.7735,       # s3close/hk2, PROVED over lambda in [1,3/2]
    "K2_measured_window": 5.3854,       # s3close/hk2, MEASURED over the window
    "Ghat_measured": 0.98142,           # u2 sec.9, max r|grad a|, lam=3/2, C^1-perturbed map
    "n_dim": 5,
}

# --------------------------------------------------------------- datum constants (t1)
T1 = json.load(open("t1_results.json"))
KAPPA = T1["clock"]["kappa_delta"]                 # FULL angular profile (taper AND equator)
CSTAR = T1["clock"]["c_star"]
C2 = T1["clock"]["c2"]
LAM_MAX = T1["clock"]["lam_max_apriori"]
R_H = T1["r_h"]["theorem_datum_1_to_lam_max"]
E0 = T1["datum_norms"]["THEOREM_datum"]["E0_sup_rho_eta"]
G0 = T1["datum_norms"]["THEOREM_datum"]["G0_sup_rho2_grad_eta"]
SHIFT = T1["energy"]["shift_THEOREM_datum"]        # 2 log s + (2/5) log C_E, tanh radial ramp
SHIFT_SHARP = T1["energy"]["shift_sharp"]          # ASSEMBLY (1.2), sharp radial edges
DEG = math.pi/180.0
DELTA, DM, PHI0 = 7.5*DEG, 5.0*DEG, 30.0*DEG
S_REC = 1.0/math.sin(DELTA)
LAM_OM = 1.5                                       # M_s/M <= 3/2, the contradiction hypothesis
ELL_RAMP = 0.25                                    # e-folds lost at the outer tanh edge (= eps_r)
SENS_MEASURED_FACTOR = 0.34913942333376546/2.1181705106873974   # record: 1/6.0668

# --------------------------------------------------------------- helpers (unchanged)
def J_pow(c, p):
    k = KAPPA*p
    return c if abs(k) < 1e-14 else (math.exp(k*c)-1.0)/k

def eps_Tprime(mu, muJ, r_h=None, sens=1.0):
    r_h = R_H if r_h is None else r_h
    if mu >= 0.999:
        return float("inf")
    return sens*(2.0*math.pi/r_h)*(3.0*mu*(1.0+muJ)/(2.0*(1.0-mu)**5) + 3.0*muJ/8.0)

def ell_loss(f, mu=0.0):
    return (math.log(2.0) + math.log(1.0+f) + 3.0*math.log(LAM_MAX)
            + math.log((1.0+mu)/max(1.0-mu, 1e-12)) + ELL_RAMP)

# --------------------------------------------------------------- viscous budget (UNCHANGED)
def gauss_tail_aniso(Rad, sy, sz, nz=4000):
    if Rad <= 0:
        return 1.0
    zs = np.linspace(-8.0*math.sqrt(2*sz), 8.0*math.sqrt(2*sz), nz+1)
    w = np.exp(-zs**2/(4*sz))/math.sqrt(4*math.pi*sz)
    rem = Rad**2 - zs**2
    p = np.where(rem <= 0, 1.0, ncx2.sf(np.maximum(rem, 0.0)/(2*sy), 4, 0.0))
    return float(np.trapz(p*w, zs))

def viscous_budget(L, c, C_K, f, phi0=PHI0, s=S_REC, delta=DELTA, c_G=None,
                   Cprime=0.0, use_V1_for_Z=True):
    if c_G is None:
        c_G = (1.5 + Cprime/L)*c
    n = RECORD["n_dim"]
    r0 = (1.0+f)*math.sin(phi0)
    nutau = c/(s*s*L)
    sy = nutau*J_pow(c, -2.0)/c
    sz = nutau*J_pow(c,  4.0)/c
    rho_s = 1.0+f
    d_inner = f
    d_taper = rho_s*math.sin(phi0-delta)
    d_eq = rho_s*math.sin(math.pi/2 - DM - phi0)
    d = max(min(d_inner, d_taper, d_eq), 1e-12)
    R_minus = r0 - d
    N = r0/math.sin(delta)
    eps_bulk = sy/(r0*r0)
    K2 = C_K
    V1 = n*nutau*(math.exp(2*c_G)-1.0)/c_G
    tau = c/L
    EW = 0.5*K2*math.exp(c_G)*tau*V1
    if R_minus <= 0 or d <= 1e-11:
        return {"eps_bulk": float("inf"), "E_hess": float("inf"), "E_4": float("inf"),
                "E_tail": float("inf"), "d": d, "R_minus": R_minus, "N": N,
                "sigma_y": sy, "sigma_z": sz, "q_half": 0.0, "c_G": c_G,
                "eps_v": float("inf")}
    E_hess = EW*(r0/(R_minus**2)) + 4.0*N*EW/d
    trC = 8*sy + 2*sz
    trC2 = 4*(2*sy)**2 + (2*sz)**2
    trC3 = 4*(2*sy)**3 + (2*sz)**3
    E4 = (r0/(R_minus**5))*(trC**2 + 2*trC2)
    P_a_half = gauss_tail_aniso(d/2.0, sy, sz)
    P_a_full = gauss_tail_aniso(d, sy, sz)
    if use_V1_for_Z:
        q = (d/2.0)**2/nutau
        P_Z_half = min(1.0, 2*n*math.exp(-math.exp(-2*c_G)*q/(4.0*n)))
    else:
        P_Z_half = P_a_half
    P = P_Z_half + P_a_half + P_a_full
    m2 = trC
    m4 = trC**2 + 2*trC2
    m6 = trC**3 + 6*trC*trC2 + 8*trC3
    E_tail = 2*N*P + P + math.sqrt(P)*(math.sqrt(m2)/r0 + math.sqrt(m4)/r0**2
                                       + math.sqrt(m6)/r0**3)
    return {"eps_bulk": eps_bulk, "E_hess": E_hess, "E_4": E4, "E_tail": min(E_tail, 1e12),
            "d": d, "R_minus": R_minus, "N": N, "sigma_y": sy, "sigma_z": sz,
            "q_half": (d/2.0)**2/nutau, "c_G": c_G,
            "eps_v": min(eps_bulk + E_hess + E4 + min(E_tail, 1e12), 1e12)}

# --------------------------------------------------------------- bootstrap (REPLACED)
def Lrad_from_m(m, lam_mx):
    lo = lam_mx**-2 - m
    if lo <= 1e-9:
        return float("inf")
    return max(math.log(lam_mx + m), math.log(1.0/lo))

def bootstrap(L, c, Cpp, C_R, ca_sup, C1=0.0, lam_mx=LAM_MAX, lam_om=LAM_OM):
    """
      G          = (3/2)(1 + C_1/L) + C''/L                       (Gamma-off; C_1 = 0, U2 7(a))
      dm/dtheta  = G m + lam_mx (C_R + 2 lam_om log lam_mx)/L ,    m(0) = 0
      dw/dtheta  = (2/L)[ ca_sup + lam_om Lrad(m)/2 ] ,            w(0) = 0
      mu = lam_mx^2 m ,   mu_J = e^w - 1
    Everything is against M, so C_R and ca_sup are absolute constants (t2), not multiples of M_s.
    """
    G = 1.5*(1.0 + C1/L) + Cpp/L
    drive = lam_mx*(C_R + 2.0*lam_om*math.log(lam_mx))/L

    def rhs(th, y):
        m, w = y
        Lr = Lrad_from_m(m, lam_mx)
        if not np.isfinite(Lr):
            return [0.0, 0.0]
        return [G*m + drive, (2.0/L)*(ca_sup + 0.5*lam_om*Lr)]

    def ev(th, y):
        return y[0] - (lam_mx**-2 - 1e-6)
    ev.terminal = True
    ev.direction = 1
    sol = solve_ivp(rhs, (0.0, c), [0.0, 0.0], rtol=1e-10, atol=1e-14, events=ev,
                    max_step=c/200.0)
    m, w = float(sol.y[0, -1]), float(sol.y[1, -1])
    blew = (sol.t[-1] < c - 1e-12) or (not sol.success)
    mu = lam_mx**2*m if not blew else float("inf")
    muJ = (math.exp(w) - 1.0) if (not blew and w < 50) else float("inf")
    return {"m": m, "mu": mu, "muJ": muJ, "c_G": G*c, "Gamma_over_ML": G,
            "closed": (not blew) and mu < 1.0, "feedback_exponent_pc": G*c, "drive": drive}

# --------------------------------------------------------------- the three columns
COLUMNS = {
    #  name                : Gamma-off route, C_a kind, Ghat kind, K2 key, Lemma T' sensitivity
    "proved":              ("u2_fixedpoint", "far_near", "proved",   "K2_proved_window",   1.0,  True),
    "proved_Ca_measured":  ("u2_fixedpoint", "offset",   "proved",   "K2_proved_window",   1.0,  True),
    "measured":            ("L3v_measured",  "offset",   "measured", "K2_measured_window",
                            SENS_MEASURED_FACTOR, False),
}

_CACHE = {}

def column_constants(L, c, column):
    key = (round(L, 9), round(c, 12), column)
    if key in _CACHE:
        return _CACHE[key]
    route, ca_kind, gh_kind, k2key, sens, use_V1 = COLUMNS[column]
    if route == "u2_fixedpoint":
        fp = T2.gamma_off(L, c, LAM_OM, E0, G0, sigma_star=0.0)
        if fp is None:
            _CACHE[key] = None
            return None
        Cpp, c_G, Ghat = fp["C_pp"], fp["c_G"], fp["Ghat"]
    else:
        Cpp = RECORD["Cprime_measured"]
        c_G = (1.5 + Cpp/L)*c
        Ghat = RECORD["Ghat_measured"]
    if gh_kind == "measured":
        Ghat = RECORD["Ghat_measured"]
    ca_flat = RECORD["material_offset_M"] if ca_kind == "offset" else None
    C_R, argphi = T2.C_R_sup(LAM_OM, c_G, E0, Ghat, ca_flat=ca_flat,
                             grad_sinphi=(gh_kind == "proved"), phi_lo_deg=0.0, n=400)
    ca_sup = ca_flat if ca_flat is not None else T2.chat_a(0.0, LAM_OM, c_G, E0)
    ca_traj = ca_flat if ca_flat is not None else T2.chat_a(math.sin(PHI0), LAM_OM, c_G, E0)
    out = {"C_pp": Cpp, "c_G": c_G, "Ghat": Ghat, "C_R": C_R, "C_R_argmax_phi_deg": argphi,
           "ca_sup": ca_sup, "ca_traj": ca_traj, "C_K": RECORD[k2key], "sens": sens,
           "use_V1": use_V1}
    _CACHE[key] = out
    return out

def assemble(L, c, column, C_K=None, f=0.0, lam_mx=LAM_MAX, C_R_override=None):
    K = column_constants(L, c, column)
    if K is None:                                # (Gamma-off)'s own fixed point does not exist
        return {"L": L, "c": c, "column": column, "f": f, "eps": float("inf"),
                "gamma_off_exists": False}
    C_R = K["C_R"] if C_R_override is None else C_R_override
    bs = bootstrap(L, c, K["C_pp"], C_R, K["ca_sup"], C1=0.0, lam_mx=lam_mx)
    CK = K["C_K"] if C_K is None else C_K
    vb = viscous_budget(L, c, CK, f, c_G=K["c_G"], use_V1_for_Z=K["use_V1"])
    epsT = eps_Tprime(bs["mu"], bs["muJ"], R_H, K["sens"])
    ell = ell_loss(f, min(bs["mu"], 0.999) if np.isfinite(bs["mu"]) else 0.999)
    eps_a = ell/L + epsT + K["ca_traj"]/(KAPPA*L)
    eps_v = vb["eps_v"]
    eps_delta = 1.0 - 2.0*KAPPA
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        eps = float("inf")
    else:
        eps = (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-eps_delta)) - 1.0
    return {"L": L, "c": c, "column": column, "C_K": CK, "f": f, "gamma_off_exists": True,
            "eps_ell": ell/L, "eps_Tprime": epsT, "eps_Ca": K["ca_traj"]/(KAPPA*L),
            "eps_a": eps_a, "eps_bulk": vb["eps_bulk"], "E_hess": vb["E_hess"],
            "E_4": vb["E_4"], "E_tail": vb["E_tail"], "eps_v": eps_v, "eps_delta": eps_delta,
            "mu": bs["mu"], "muJ": bs["muJ"], "c_G": K["c_G"], "C_pp": K["C_pp"],
            "Ghat": K["Ghat"], "C_R": C_R, "ca_traj": K["ca_traj"], "ca_sup": K["ca_sup"],
            "feedback_exponent_pc": bs["feedback_exponent_pc"], "d": vb["d"],
            "eps": eps, "bootstrap_closed": bs["closed"]}

FS = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]

def eps_at(L, col, **kw):
    best, arg = float("inf"), None
    for f in FS:
        r = assemble(L, CSTAR, col, f=f, **kw)
        if r["eps"] < best:
            best, arg = r["eps"], f
    return best, arg

def Lstar_for(col, target, **kw):
    lo, hi = 2.0, 1e12
    if eps_at(hi, col, **kw)[0] > target:
        return None
    for _ in range(300):
        mid = math.sqrt(lo*hi)
        if eps_at(mid, col, **kw)[0] <= target:
            hi = mid
        else:
            lo = mid
        if hi/lo < 1.0000001:
            break
    return hi

if __name__ == "__main__":
    OUT = {"datum": {"kappa_delta": KAPPA, "c_star": CSTAR, "c2": C2, "lam_max": LAM_MAX,
                     "r_h_1_to_lam_max": R_H, "E0": E0, "G0": G0, "s_record": S_REC,
                     "logReE_shift_theorem_datum": SHIFT,
                     "logReE_shift_sharp_edges": SHIFT_SHARP,
                     "ell_ramp": ELL_RAMP, "eps_delta": 1.0-2.0*KAPPA,
                     "eps_floor": (1.0-2.0*KAPPA)/(2.0*KAPPA)},
           "record_inputs": RECORD}

    Ls = [1e3, 1e4, 1e5, 1e6]
    # ---- constants per column, per L
    cc = {}
    for col in COLUMNS:
        for L in Ls + [3e4, 3e5, 1e9]:
            cc["%s_L=%g" % (col, L)] = column_constants(L, CSTAR, col)
    OUT["column_constants"] = cc
    OUT["L_gamma_star"] = T2.L_gamma_star(CSTAR, LAM_OM, E0, G0, 0.0)

    # ---- the term table (f = 1)
    tt = []
    for col in COLUMNS:
        for L in Ls + [3e5]:
            tt.append(assemble(L, CSTAR, col, f=1.0))
    OUT["term_table_f1"] = tt

    # ---- eps(L), best f
    el = {}
    for col in COLUMNS:
        for L in Ls:
            e, f = eps_at(L, col)
            el["%s_L=%g" % (col, L)] = {"eps": e, "best_f": f}
    OUT["eps_at_L"] = el

    # ---- L_* and log Lambda_*
    ls, lam = {}, {}
    for col in COLUMNS:
        for target in [0.5, 0.1]:
            k = "%s_eps<=%g" % (col, target)
            ls[k] = Lstar_for(col, target)
            lam[k] = {"log_Lambda_star_theorem_datum": (2*ls[k] + SHIFT) if ls[k] else None,
                      "log_Lambda_star_if_sharp_edges": (2*ls[k] + SHIFT_SHARP) if ls[k] else None}
    OUT["L_star"] = ls
    OUT["log_Lambda_star"] = lam

    # ---- what dominates: the drive decomposition at L = L_*(proved, 1/2)
    Lp = ls["proved_eps<=0.5"]
    if Lp:
        K = column_constants(Lp, CSTAR, "proved")
        cG, Gh = K["c_G"], K["Ghat"]
        pieces = {}
        for tag, kw in [("full", {}),
                        ("Ghat_set_to_zero", {"Ghat": 0.0}),
                        ("ca_far_near_only", {"ca_flat": 0.069})]:
            gh = kw.get("Ghat", Gh)
            caf = kw.get("ca_flat", None)
            v, ar = T2.C_R_sup(LAM_OM, cG, E0, gh, ca_flat=caf, grad_sinphi=True, n=400)
            pieces[tag] = {"C_R": v, "argmax_phi_deg": ar}
        pieces["term_C_share"] = 1.0 - pieces["Ghat_set_to_zero"]["C_R"]/pieces["full"]["C_R"]
        OUT["dominance_at_Lstar_proved"] = {"L_star": Lp, "c_G": cG, "Ghat": Gh,
                                            "C_pp": K["C_pp"], "C_R_pieces": pieces,
                                            "G0": G0, "G0_of_campaign_DC": 20.9197070}

    # ---- sensitivities the theorem's modulo list needs priced
    sens = {}
    sens["lam_max=3/2 (self-consistent cap)"] = {
        "L_star_0.5": Lstar_for("proved", 0.5, lam_mx=1.5),
        "L_star_0.5_measured": Lstar_for("measured", 0.5, lam_mx=1.5)}
    for ck, tag in [(1.0, "C_K=1"), (5.3854, "C_K=hk2_measured"), (66.6622, "C_K=hk2_proved_lam1"),
                    (161.7735, "C_K=hk2_proved_window")]:
        sens[tag] = {"L_star_0.5_proved": Lstar_for("proved", 0.5, C_K=ck)}
    # C_R sensitivity: what a factor in the drive costs (it is linear, not exponential)
    Kref = column_constants(1e6, CSTAR, "proved")
    for mult in [0.5, 1.0, 2.0, 5.0]:
        sens["C_R x %g" % mult] = {
            "L_star_0.5_proved": Lstar_for("proved", 0.5, C_R_override=mult*Kref["C_R"])}
    OUT["sensitivity"] = sens

    with open("t3_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)

    print("kappa=%.13f c*=%.10f c2=%.10f lam_max=%.10f r_h=%.10f"
          % (KAPPA, CSTAR, C2, LAM_MAX, R_H))
    print("E0=%.8f G0=%.8f shift(theorem datum)=%.8f shift(sharp)=%.8f"
          % (E0, G0, SHIFT, SHIFT_SHARP))
    print("L_gamma_star =", OUT["L_gamma_star"])
    print("\neps at c=c_*, best f:")
    for k in sorted(el):
        print("  %-28s eps=%12.6g  f=%s" % (k, el[k]["eps"], el[k]["best_f"]))
    print("\nL_star:", json.dumps(ls, indent=1))
    print("log Lambda_star:", json.dumps(lam, indent=1))
    print("\ndominance:", json.dumps(OUT.get("dominance_at_Lstar_proved", {}), indent=1))
    print("\nsensitivity:", json.dumps(sens, indent=1))

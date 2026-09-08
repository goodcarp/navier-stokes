"""
u5 -- the error budget with the SLAVED reference.

This file is a COPY of  s3close/assembly/a2_budget.py  (sha256 96c1ee0b10be43c3...),
with exactly three things replaced:

  * bootstrap()  -- the majorant system.  The a priori reference is replaced by the true
                    field's own exterior l=1 strain, so the "l=1 rate error" term
                    kappa * eps_T'  DISAPPEARS from the driver of (2.1)-(2.2).  What is left
                    is  d m/dtheta = G m + lam_max lam_om (C_R + 2 log lam_max)/L  and a
                    Jacobian equation that compares J_Phi with lambda^2 directly (so Lemma T'
                    Corollary 3's composition, and its 2 kappa_s/L, are not needed either).
  * COLUMNS      -- the three columns the brief asks for.
  * assemble()   -- eps_a now uses r_h over [1, exp(3c/4)], ell_loss with lam_max = exp(3c/4),
                    and eps_T' evaluated ONCE at the end state.

viscous_budget(), gauss_tail_aniso(), J_pow(), C_a_proved(), eps_Tprime() and the RECORD
block are carried over UNCHANGED from a2_budget.py; the datum constants are re-read from
u2_results.json (computed in this folder) instead of a1_results.json.

Window:  tau = c/(M L),  theta = M L t in [0,c].
"""
import json, math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.stats import ncx2

# --------------------------------------------------------------- record inputs
RECORD = {
    # write/lemma-T-shell-dependent/PROOF.md  (refereed: STANDS)
    "lemmaT_abs":  "|a[eta0.Phi^-1](0) - a_ref| <= pi M A [3 mu(1+muJ)/(2(1-mu)^5) + 3 muJ/8]",
    "kappa_s": math.sqrt(6.0) - 2.0,                 # sup |rho lam'/lam| * L, integro-ODE profile
    # write/L3v-and-gamma-bound + refute-*  (Theorem A / Cor A1 PROVED; C' corrected)
    "Cprime_A1_computedV": 119.33,
    "Cprime_A1_provedV":   151.15,
    "Cprime_A1_oddonly":    76.11,
    "Cprime_computed_L3v":  33.45,
    "Cprime_measured_L3v":  18.5,
    "Cprime_measured":      11.74,
    # rebuild/far-near-kernel-lemma  (z-odd constants)
    "C_far_zodd": 0.291999,
    "C_inner_zodd": 0.014754,
    "C_collar_over_sinphi": 3.999218,
    "C_collar_axis": math.pi/8.0,
    # lower/refuter-lagrangian-b: material-point strain offset, MEASURED
    "material_offset_M": 0.069,
    # write/V-b + refute-V-b : Theorem V.4 and its corrected sizing
    "V4_terms": "E_hess + E_4 + E_tail as in Theorem V.4",
    "n_dim": 5,
}

# --------------------------------------------------------------- datum constants
U2 = json.load(open("u2_results.json"))
DEG = math.pi/180.0
DELTA = 7.5*DEG
DM = 5.0*DEG
KAPPA = U2["kappa_delta_7.5"]                 # computed in u2 (= record to 0.0)
R_H_15 = U2["r_h_1_to_1.5"]                   # r_h on [1,3/2]      (a2's r_h)
R_H = U2["r_h_1_to_lam_max"]                  # r_h on [1,exp(3c/4)] -- the slaved range
S_REC = 1.0/math.sin(DELTA)          # rho0 = s sqrt(nu/M) with the record's normalisation
LOGREE_SHIFT = U2["RECORD"]["chk_logReE_shift"]   # 2 log s + (2/5) log C_E, ASSEMBLY (1.2)
KAPPA_S = 0.75*U2["c_star"]                   # slaved reference: sup L|rho lam'/lam| = 3c/4
PHI0 = 30.0*DEG
CSTAR = U2["c_star"]
LAM_MAX = U2["lam_max_apriori"]               # exp(3 c_*/4) = 1.8377...  (a2 used 3/2)
LAM_OM = 1.5                                  # M_s/M <= 3/2  (contradiction hypothesis)
C_R_PROVED = U2["C_R_proved_Gc0"]             # G_collar = 0
C_R_MEASURED = U2["C_R_measured_Gc0"]
RECORD["Cprime_A1_provedV_brief"] = 151.15
RECORD["K2_proved_window"] = U2["RECORD"]["K2_proved_window"]
RECORD["K2_measured_window"] = U2["RECORD"]["K2_measured_window"]

# --------------------------------------------------------------- helpers
def lam_frozen(theta):
    return math.exp(KAPPA*theta)

def J_pow(c, p):
    """int_0^c lambda(theta)^p dtheta for the frozen model"""
    k = KAPPA*p
    return c if abs(k) < 1e-14 else (math.exp(k*c)-1.0)/k

def ell_loss(f, mu=0.0):
    """e-folds lost from L so that the image of the outer sub-shell is outside 2|X(s)|:
       (1-mu) lam^-2 rho'' >= 2 (1+mu) lam (1+f) rho_0 ."""
    return (math.log(2.0) + math.log(1.0+f) + 3.0*math.log(LAM_MAX)
            + math.log((1.0+mu)/max(1.0-mu, 1e-12)))

def C_a_proved(sinphi):
    return (RECORD["C_far_zodd"] + RECORD["C_inner_zodd"]
            + RECORD["C_collar_over_sinphi"]/sinphi + RECORD["C_collar_axis"])

def phi_of_lambda(phi0, lam):
    return math.atan2(lam**3*math.sin(phi0), math.cos(phi0))

def g_geom(phi, lam):
    return math.sqrt(lam**2*math.sin(phi)**2 + lam**(-4)*math.cos(phi)**2)

# --------------------------------------------------------------- Lemma T' sensitivity
SENS_MEASURED_FACTOR = 0.34913942333376546/2.1181705106873974   # record: proved/measured = 6.0668

def eps_Tprime(mu, muJ, r_h=R_H, sens=1.0):
    """sens = 1 : Lemma T' as proved.  sens = SENS_MEASURED_FACTOR : the same expression
    scaled down by the record's own measured conservatism factor 6.0668 (NOT a theorem)."""
    if mu >= 0.999:
        return float("inf")
    return sens*(2.0*math.pi/r_h)*(3.0*mu*(1.0+muJ)/(2.0*(1.0-mu)**5) + 3.0*muJ/8.0)

# --------------------------------------------------------------- bootstrap majorant ODE
#  *** REPLACED (this is the whole change).  The reference rate is the TRUE field's own
#      exterior l=1 strain at the label radius, so the l=1 rate error vanishes identically
#      and eps_T' does not appear on the right-hand side. ***
def Lrad_from_m(m, lam_mx):
    """|log(|Phi|/|x|)| with |Phi| in [lam^-2|x| - m|x|,  lam|x| + m|x|]."""
    lo = lam_mx**-2 - m
    if lo <= 1e-9:
        return float("inf")
    return max(math.log(lam_mx + m), math.log(1.0/lo))

def bootstrap(L, c, Cpp, C_R, C_a, C1=1.0, lam_mx=LAM_MAX, lam_om=LAM_OM):
    """
      G           = (3/2)(1 + C_1/L) + C''/L                     (hypothesis (Gamma-off))
      d m/dtheta  = G m + lam_mx * lam_om * (C_R + 2 log lam_mx)/L ,      m(0) = 0
      d w/dtheta  = (2 lam_om/L) [ C_a + (1/2) Lrad(m) ] ,                w(0) = 0
      mu          = lam_mx^2 m ,   mu_J = e^w - 1
    """
    G = 1.5*(1.0 + C1/L) + Cpp/L
    c_G = G*c
    drive = lam_mx*lam_om*(C_R + 2.0*math.log(lam_mx))/L

    def rhs(th, y):
        m, w = y
        Lr = Lrad_from_m(m, lam_mx)
        if not np.isfinite(Lr):
            return [0.0, 0.0]
        return [G*m + drive, (2.0*lam_om/L)*(C_a + 0.5*Lr)]

    def ev(th, y):
        return y[0] - (lam_mx**-2 - 1e-6)
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, (0.0, c), [0.0, 0.0], rtol=1e-10, atol=1e-14, events=ev,
                    max_step=c/200.0)
    m, w = float(sol.y[0, -1]), float(sol.y[1, -1])
    blew = (sol.t[-1] < c - 1e-12) or (not sol.success)
    mu = lam_mx**2*m if not blew else float("inf")
    muJ = (math.exp(w) - 1.0) if (not blew and w < 50) else float("inf")
    return {"m": m, "mu": mu, "muJ": muJ, "c_G": c_G, "Gamma_over_ML": G,
            "closed": (not blew) and mu < 1.0,
            "feedback_exponent_pc": G*c, "p_over_ML": G, "drive": drive}

# --------------------------------------------------------------- viscous budget (Thm V.4)
def gauss_tail_aniso(Rad, sy, sz, nz=4000):
    """P(|Z|>Rad) for Z ~ N(0, diag(2 sy I4, 2 sz)) in R^5:
       |Z|^2 = 2 sy X + z^2, X ~ chi^2_4, z ~ N(0,2 sz)."""
    if Rad <= 0:
        return 1.0
    zs = np.linspace(-8.0*math.sqrt(2*sz), 8.0*math.sqrt(2*sz), nz+1)
    w = np.exp(-zs**2/(4*sz))/math.sqrt(4*math.pi*sz)
    rem = Rad**2 - zs**2
    p = np.where(rem <= 0, 1.0, ncx2.sf(np.maximum(rem, 0.0)/(2*sy), 4, 0.0))
    return float(np.trapz(p*w, zs))

def viscous_budget(L, c, C_K, f, phi0=PHI0, s=S_REC, delta=DELTA, c_G=None,
                   Cprime=0.0, use_V1_for_Z=True):
    """
    All lengths in units of rho0 = 1;  nu = M/s^2 (since rho0 = s sqrt(nu/M)),
    tau = c/(M L)  =>  nu tau = c/(s^2 L).
    """
    if c_G is None:
        c_G = (1.5 + Cprime/L)*c
    n = RECORD["n_dim"]
    r0 = (1.0+f)*math.sin(phi0)                 # in units of rho0
    nutau = c/(s*s*L)                            # nu*tau
    sy = nutau*J_pow(c, -2.0)/c                  # nu int lam^-2 dt  (dt = dtheta/(ML))
    sz = nutau*J_pow(c,  4.0)/c
    # ---- distance from the tracked point to the nearest datum discontinuity ----
    rho_s = 1.0+f
    d_inner = f                                   # to the sphere |x| = rho0
    d_taper = rho_s*math.sin(phi0-delta)          # to the taper cone phi = delta
    d_eq    = rho_s*math.sin(math.pi/2 - DM - phi0)   # to the equatorial mollification layer
    d = max(min(d_inner, d_taper, d_eq), 1e-12)
    R_minus = r0 - d
    N = r0/math.sin(delta)                        # ||eta0||_inf r0 / M
    eps_bulk = sy/(r0*r0)
    # ---- E_hess (Theorem V.4, refuter's form) ----
    K2 = C_K                                      # K2 = C_K M/rho0, rho0 = 1
    V1 = n*nutau*(math.exp(2*c_G)-1.0)/c_G
    tau = c/L                                     # M tau, with M = 1
    EW = 0.5*K2*math.exp(c_G)*tau*V1
    if R_minus <= 0 or d <= 1e-11:
        return {"eps_bulk": float("inf"), "E_hess": float("inf"), "E_4": float("inf"),
                "E_tail": float("inf"), "d": d, "R_minus": R_minus, "N": N,
                "sigma_y": sy, "sigma_z": sz, "q_half": 0.0, "c_G": c_G,
                "eps_v": float("inf")}
    E_hess = EW*(r0/(R_minus**2)) + 4.0*N*EW/d
    # ---- E_4 ----
    trC = 8*sy + 2*sz
    trC2 = 4*(2*sy)**2 + (2*sz)**2
    trC3 = 4*(2*sy)**3 + (2*sz)**3
    E4 = (r0/(R_minus**5))*(trC**2 + 2*trC2) if R_minus > 0 else float("inf")
    # ---- E_tail ----
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

# --------------------------------------------------------------- full assembly
#  *** REPLACED: the three columns the brief asks for. ***
COLUMNS = {
    # name : (C'' key, C_a kind, C_R, use Thm V.1 for P(|Z|>d/2), Lemma T' sensitivity, K2 key)
    "proved":             ("Cprime_A1_provedV_brief", "far_near", "proved",   True,  1.0,
                           "K2_proved_window"),
    "proved_Ca_measured": ("Cprime_A1_provedV_brief", "offset",   "measured", True,  1.0,
                           "K2_proved_window"),
    "measured":           ("Cprime_measured",         "offset",   "measured", False,
                           SENS_MEASURED_FACTOR, "K2_measured_window"),
}

def assemble(L, c, column, C_K=None, f=0.0, C1=1.0, G_collar=0.0, lam_mx=LAM_MAX):
    cp_key, ca_kind, cr_kind, use_V1, sens, k2key = COLUMNS[column]
    Cpp = RECORD[cp_key]
    C_a = C_a_proved(math.sin(PHI0)) if ca_kind == "far_near" else RECORD["material_offset_M"]
    C_R = (C_R_PROVED if cr_kind == "proved" else C_R_MEASURED) + G_collar
    if C_K is None:
        C_K = RECORD[k2key]
    bs = bootstrap(L, c, Cpp, C_R, C_a, C1=C1, lam_mx=lam_mx)
    vb = viscous_budget(L, c, C_K, f, c_G=bs["c_G"], Cprime=Cpp, use_V1_for_Z=use_V1)
    r_h_used = r_h_of(lam_mx)
    epsT = eps_Tprime(bs["mu"], bs["muJ"], r_h_used, sens)
    ell = ell_loss(f, min(bs["mu"], 0.999) if np.isfinite(bs["mu"]) else 0.999)
    eps_a = ell/L + epsT + LAM_OM*C_a/(KAPPA*L)
    eps_v = vb["eps_v"]
    eps_delta = 1.0 - 2.0*KAPPA
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        eps = float("inf")
    else:
        eps = (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-eps_delta)) - 1.0
    return {"L": L, "c": c, "column": column, "C_K": C_K, "f": f,
            "eps_ell": ell/L, "eps_Tprime": epsT,
            "eps_Ca": LAM_OM*C_a/(KAPPA*L), "eps_a": eps_a,
            "eps_bulk": vb["eps_bulk"], "E_hess": vb["E_hess"], "E_4": vb["E_4"],
            "E_tail": vb["E_tail"], "eps_v": eps_v, "eps_delta": eps_delta,
            "mu": bs["mu"], "muJ": bs["muJ"], "c_G": bs["c_G"], "r_h": r_h_used,
            "feedback_exponent_pc": bs["feedback_exponent_pc"],
            "Cprime": Cpp, "C_a": C_a, "C_R": C_R, "d": vb["d"], "q_half": vb["q_half"],
            "eps": eps, "bootstrap_closed": bs["closed"]}

def r_h_of(lam_mx):
    """r_h on [1, lam_mx]; u2 stores the two values this file uses."""
    if abs(lam_mx - LAM_MAX) < 1e-12:
        return R_H
    if abs(lam_mx - 1.5) < 1e-12:
        return R_H_15
    raise ValueError("r_h not tabulated for lam_max = %r" % lam_mx)

if __name__ == "__main__":
    cstar = CSTAR
    OUT = {"record_inputs": {k: v for k, v in RECORD.items()},
           "datum": {"kappa_delta": KAPPA, "r_h_slaved": R_H, "r_h_1to1.5": R_H_15,
                     "s_record": S_REC, "logReE_shift": LOGREE_SHIFT,
                     "delta_deg": 7.5, "dm_deg": 5.0, "phi0_deg": 30.0,
                     "lam_max": LAM_MAX, "kappa_s_slaved": KAPPA_S,
                     "C_R_proved": C_R_PROVED, "C_R_measured": C_R_MEASURED,
                     "c_at_lambda_3over2": cstar}}
    Ls = [40.0, 160.0, 640.0, 2560.0]
    FS = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]

    # ---- the grid at c = c_* (and the a2 c-values, for shape) ----
    grid = []
    for col in list(COLUMNS):
        for L in Ls + [10.0]:
            for c in [0.25, 0.5, 1.0, cstar]:
                grid.append(assemble(L, c, col, f=1.0))
    OUT["grid"] = grid

    # ---- eps(L) at c = c_*, best f ----
    def eps_at(L, col, **kw):
        best, arg = float("inf"), None
        for f in FS:
            r = assemble(L, cstar, col, f=f, **kw)
            if r["eps"] < best:
                best, arg = r["eps"], f
        return best, arg
    epsL = {}
    for col in list(COLUMNS):
        for L in Ls:
            e, f = eps_at(L, col)
            epsL["%s_L=%g" % (col, L)] = {"eps": e, "best_f": f}
    OUT["eps_at_L"] = epsL

    # ---- L_* ----
    def Lstar_for(col, target, **kw):
        lo, hi = 2.0, 1e18
        if eps_at(hi, col, **kw)[0] > target:
            return None
        for _ in range(200):
            mid = math.sqrt(lo*hi)
            if eps_at(mid, col, **kw)[0] <= target:
                hi = mid
            else:
                lo = mid
            if hi/lo < 1.0000001:
                break
        return hi
    Lstar, Lam = {}, {}
    for col in list(COLUMNS):
        for target in [0.5, 0.1]:
            k = "%s_eps<=%g" % (col, target)
            Lstar[k] = Lstar_for(col, target)
            Lam[k] = (2*Lstar[k] + LOGREE_SHIFT) if Lstar[k] else None
    OUT["L_star"] = Lstar
    OUT["Lambda_star_log"] = Lam
    OUT["logReE_minus_2L"] = LOGREE_SHIFT
    OUT["c2"] = 2.0*math.log(1.5)/KAPPA

    # ---- feedback exponent, mu(c) L ----
    fb = {}
    for col in list(COLUMNS):
        r = assemble(1e9, cstar, col, f=1.0)
        fb[col] = {"pc": r["feedback_exponent_pc"], "mu_times_L": r["mu"]*1e9,
                   "muJ_times_L": r["muJ"]*1e9, "C_R": r["C_R"], "C_a": r["C_a"],
                   "Cpp": r["Cprime"], "r_h": r["r_h"]}
        for L in Ls:
            fb[col]["pc_L=%g" % L] = assemble(L, cstar, col, f=1.0)["feedback_exponent_pc"]
    OUT["feedback"] = fb

    # ---- sensitivity: the unproved collar-gradient constant, and lam_max ----
    sens = {}
    Ca30 = C_a_proved(math.sin(PHI0))
    for gc, tag in [(0.0, "G_collar=0"), (Ca30, "G_collar=C_a"), (5*Ca30, "G_collar=5C_a")]:
        for col in ["proved", "measured"]:
            sens["%s_%s" % (col, tag)] = {
                "L_star_0.5": Lstar_for(col, 0.5, G_collar=gc),
                "eps_640": eps_at(640.0, col, G_collar=gc)[0]}
    for col in ["proved", "measured"]:
        sens["%s_lam_max=1.5" % col] = {
            "L_star_0.5": Lstar_for(col, 0.5, lam_mx=1.5),
            "eps_640": eps_at(640.0, col, lam_mx=1.5)[0]}
    OUT["sensitivity"] = sens

    # ---- C_K sweep (the hk2 datum) ----
    ck = {}
    for C_K, tag in [(1.0, "1"), (3.0202, "hk2_measured_global"),
                     (5.3854, "hk2_measured_window"), (66.6622, "hk2_proved_lam1"),
                     (161.7735, "hk2_proved_window")]:
        for col in list(COLUMNS):
            ck["%s_C_K=%s" % (col, tag)] = {
                "eps_640": eps_at(640.0, col, C_K=C_K)[0],
                "L_star_0.5": Lstar_for(col, 0.5, C_K=C_K)}
    OUT["C_K_sweep"] = ck

    with open("u5_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)

    print("kappa=%.13f r_h[1,1.5]=%.9f r_h[1,%.6f]=%.9f c*=%.7f"
          % (KAPPA, R_H_15, LAM_MAX, R_H, cstar))
    print("C_R proved=%.6f measured=%.6f ; logReE-2L=%.7f ; c2=%.7f"
          % (C_R_PROVED, C_R_MEASURED, LOGREE_SHIFT, OUT["c2"]))
    print("\nfeedback:")
    print(json.dumps(fb, indent=1))
    print("\neps at c=c_*, best f:")
    for k in sorted(epsL):
        print("  %-28s eps=%12.6g  f=%s" % (k, epsL[k]["eps"], epsL[k]["best_f"]))
    print("\nL_star:", json.dumps(Lstar, indent=1))
    print("log Lambda_star:", json.dumps(Lam, indent=1))
    print("\nsensitivity:", json.dumps(sens, indent=1))

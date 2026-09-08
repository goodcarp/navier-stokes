"""
a2 -- the error budget for the assembled (S3) theorem, and the bootstrap majorant ODE.

Everything is computed here.  Constants imported FROM THE RECORD are listed in
RECORD below with their file of origin; they are the *statements* of the refereed
theorems this assembly consumes, not numbers this seat could re-derive cheaply.

Window:  tau = c/(M L),  theta = M L t in [0,c].
Frozen (conservative) model: a >= kappa M L, lambda(theta) = exp(kappa theta).
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
A1 = json.load(open("a1_results.json"))
DEG = math.pi/180.0
DELTA = 7.5*DEG
DM = 5.0*DEG
KAPPA = A1["kappa_delta"]["7.5"]
R_H = A1["r_h"]["7.5"]
C_E = A1["C_E_table"]["delta=7.5,dm=5,eps=0"]
S_REC = 1.0/math.sin(DELTA)          # rho0 = s sqrt(nu/M) with the record's normalisation
KAPPA_S = RECORD["kappa_s"]
PHI0 = 30.0*DEG
LAM_MAX = 1.5

# --------------------------------------------------------------- helpers
def lam_frozen(theta):
    return math.exp(KAPPA*theta)

def J_pow(c, p):
    """int_0^c lambda(theta)^p dtheta for the frozen model"""
    k = KAPPA*p
    return c if abs(k) < 1e-14 else (math.exp(k*c)-1.0)/k

def ell_loss(f):
    """e-folds lost from L: one octave for the far/near split, the inset, and 3 log lambda"""
    return math.log(2.0) + math.log(1.0+f) + 3.0*math.log(LAM_MAX)

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
def bootstrap(L, c, Cprime, C_a, kappa=KAPPA, r_h=R_H, kappa_s=KAPPA_S,
              lamfac_grad=1.5, lamfac_len=1.5, geom_lo=4.0/9.0, geom_hi=1.5, sens=1.0):
    """
    Majorant system on theta in [0,c] (theta = M L t):

      Gammabar/(M L) = lamfac_grad + Cprime/L        (Theorem Gamma: Gamma <= 2a(0,t)+C'M,
                                                      2a(0,t) <= lamfac_grad * M L on the window)
      dm/dtheta   = (Gammabar/(ML)) m + geom_hi*(2*(kappa*epsT + kappa*c_G/L)
                                                 + C_a*lam/L)
      dmuJL/dth   = 2*[ kappa*epsT + kappa*c_G/L + C_a*lam/L + (Gammabar/(ML))*mu ]
      mu   = m/geom_lo
      muJ  = muJL*(1+2 kappa_s/L) + 2 kappa_s/L
      epsT = eps_Tprime(mu,muJ)

    m := sup_x |Phi(x)-Lambda(x)|/|x| ;  mu := sup |Phi-Lambda|/|Lambda| <= m/geom_lo.
    Returns the state at theta = c, plus a flag if mu left the region mu < 1.
    """
    GoML = lamfac_grad + Cprime/L
    c_G = GoML*c                        # int_0^tau Gamma dt

    def rhs(th, y):
        m, muJL = y
        mu = min(m/geom_lo, 0.9)
        muJ = muJL*(1.0+2.0*kappa_s/L) + 2.0*kappa_s/L
        epsT = eps_Tprime(mu, muJ, r_h, sens)
        lam = min(lam_frozen(th), LAM_MAX)
        drive = 2.0*(kappa*epsT + kappa*c_G/L) + C_a*lam/L
        dm = GoML*m + geom_hi*drive
        dmuJ = 2.0*(kappa*epsT + kappa*c_G/L + C_a*lam/L + GoML*mu)
        return [dm, dmuJ]

    def ev(th, y):
        return y[0]/geom_lo - 0.9
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, (0.0, c), [0.0, 0.0], rtol=1e-7, atol=1e-14, events=ev,
                    max_step=c/50.0)
    m, muJL = sol.y[0, -1], sol.y[1, -1]
    mu = m/geom_lo
    if sol.t[-1] < c - 1e-12:
        mu = 1.0
    muJ = muJL*(1.0+2.0*kappa_s/L) + 2.0*kappa_s/L
    blew = (not sol.success) or (mu >= 0.899) or (not np.isfinite(mu))
    epsT = eps_Tprime(mu, muJ, r_h, sens) if not blew else float("inf")
    # the linearised feedback exponent, the single most informative number here
    p_over_ML = lamfac_grad + Cprime/L + geom_hi*2.0*kappa*sens*(15.0*math.pi/4.0)/(r_h*geom_lo)
    return {"m": m, "mu": mu, "muJ": muJ, "eps_Tprime": epsT, "c_G": c_G,
            "Gamma_over_ML": GoML, "closed": not blew,
            "feedback_exponent_pc": p_over_ML*c, "p_over_ML": p_over_ML}

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
COLUMNS = {
    # name              : (C', C_a, use Thm V.1 for P(|Z|>d/2), Lemma T' sensitivity factor)
    "proved":            ("Cprime_A1_computedV", "far_near", True,  1.0),
    "proved_oddonly":    ("Cprime_A1_oddonly",   "far_near", True,  1.0),
    "proved_sharpsens":  ("Cprime_A1_oddonly",   "far_near", True,  SENS_MEASURED_FACTOR),
    "measured":          ("Cprime_measured",     "offset",   False, SENS_MEASURED_FACTOR),
}

def assemble(L, c, column, C_K=1.0, f=0.0):
    cp_key, ca_kind, use_V1, sens = COLUMNS[column]
    Cprime = RECORD[cp_key]
    C_a = C_a_proved(math.sin(PHI0)) if ca_kind == "far_near" else RECORD["material_offset_M"]
    bs = bootstrap(L, c, Cprime, C_a, sens=sens)
    vb = viscous_budget(L, c, C_K, f, c_G=bs["c_G"], Cprime=Cprime, use_V1_for_Z=use_V1)
    ell = ell_loss(f)
    eps_a = ell/L + bs["eps_Tprime"] + LAM_MAX*C_a/(KAPPA*L)
    eps_v = vb["eps_v"]
    eps_delta = 1.0 - 2.0*KAPPA
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        eps = float("inf")
    else:
        eps = (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-eps_delta)) - 1.0
    return {"L": L, "c": c, "column": column, "C_K": C_K, "f": f,
            "eps_ell": ell/L, "eps_Tprime": bs["eps_Tprime"],
            "eps_Ca": LAM_MAX*C_a/(KAPPA*L), "eps_a": eps_a,
            "eps_bulk": vb["eps_bulk"], "E_hess": vb["E_hess"], "E_4": vb["E_4"],
            "E_tail": vb["E_tail"], "eps_v": eps_v, "eps_delta": eps_delta,
            "mu": bs["mu"], "muJ": bs["muJ"], "c_G": bs["c_G"],
            "feedback_exponent_pc": bs["feedback_exponent_pc"],
            "Cprime": Cprime, "C_a": C_a, "d": vb["d"], "q_half": vb["q_half"],
            "eps": eps, "bootstrap_closed": bs["closed"]}

if __name__ == "__main__":
    OUT = {"record_inputs": {k: v for k, v in RECORD.items()},
           "datum": {"kappa_delta": KAPPA, "r_h": R_H, "C_E": C_E, "s_record": S_REC,
                     "delta_deg": 7.5, "dm_deg": 5.0, "phi0_deg": 30.0,
                     "c_at_lambda_3over2": math.log(1.5)/KAPPA}}
    Ls = [10.0, 40.0, 160.0, 640.0, 2560.0]
    cs = [0.25, 0.5, 1.0, math.log(1.5)/KAPPA]
    grid = []
    for col in list(COLUMNS):
        for L in Ls:
            for c in cs:
                grid.append(assemble(L, c, col, f=1.0))
    OUT["grid"] = grid
    # f-sweep at the doubling window
    cstar = math.log(1.5)/KAPPA
    fs = []
    for col in ["proved", "measured"]:
        for L in Ls:
            for f in [0.25, 0.5, 1.0, 2.0, 4.0]:
                fs.append(assemble(L, cstar, col, f=f))
    OUT["f_sweep"] = fs
    # L_* : smallest L (log-scan) with eps <= 1/2 and <= 0.1, at c = c* and best f
    def eps_at(L, col, C_K=1.0):
        return min(assemble(L, cstar, col, C_K=C_K, f=f)["eps"] for f in
                   [0.05,0.1,0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
    Lstar = {}
    for col in list(COLUMNS):
        for target in [0.5, 0.1]:
            lo, hi = 1.0, 1e60
            if eps_at(hi, col) > target:
                Lstar["%s_eps<=%g" % (col, target)] = None
                continue
            for _ in range(40):
                mid = math.sqrt(lo*hi)
                if eps_at(mid, col) <= target:
                    hi = mid
                else:
                    lo = mid
                if hi/lo < 1.0001:
                    break
            Lstar["%s_eps<=%g" % (col, target)] = hi
    OUT["L_star"] = Lstar
    # Lambda_* = exp(2 L_* + 2 log s + (2/5) log C_E)
    shift = 2*math.log(S_REC) + 0.4*math.log(C_E)
    OUT["logReE_minus_2L"] = shift
    OUT["Lambda_star_log"] = {k: (2*v + shift if v else None) for k, v in Lstar.items()}
    OUT["c2"] = 4.0*math.log(1.5)
    # self-consistent window c = c*(1+eps(c)) by fixed-point iteration
    sc = {}
    for col in list(COLUMNS):
        for L in [1e4, 1e9, 1e15, 1e18, 1e21]:
            c_it, ok = cstar, False
            for _ in range(25):
                e = min(assemble(L, c_it, col, f=f)["eps"] for f in [0.25, 1.0, 4.0])
                if not np.isfinite(e):
                    break
                cn = cstar*(1.0+e)
                if abs(cn-c_it) < 1e-10*max(1.0, c_it):
                    c_it, ok = cn, True
                    break
                c_it = 0.5*c_it + 0.5*cn
            sc["%s_L=%g" % (col, L)] = {"c_fixed": c_it if ok else None,
                                        "converged": ok,
                                        "eps": (c_it/cstar - 1.0) if ok else None}
    OUT["self_consistent_window"] = sc
    with open("a2_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
    print("kappa=%.10f r_h=%.6f C_E=%.9f s=%.6f c*=%.7f" % (KAPPA, R_H, C_E, S_REC, cstar))
    print("logReE - 2L = %.6f ; c2 = %.7f" % (shift, 4*math.log(1.5)))
    print("\n--- eps at c = c* (doubling window), f = 0 ---")
    for col in list(COLUMNS):
        for L in Ls + [1e6, 1e12, 1e18]:
            r = assemble(L, cstar, col, f=1.0)
            print("%-18s L=%9.3g  mu=%9.3g muJ=%9.3g epsT=%9.3g eps_a=%9.3g eps_v=%9.3g  eps=%9.4g  closed=%s"
                  % (col, L, r["mu"], r["muJ"], r["eps_Tprime"], r["eps_a"], r["eps_v"], r["eps"], r["bootstrap_closed"]))
    print("\nL_star:", json.dumps(Lstar, indent=1))
    print("log Lambda_star:", json.dumps(OUT["Lambda_star_log"], indent=1))

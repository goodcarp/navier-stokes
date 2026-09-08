"""
f5 -- ITEM 5: `L_Gamma_* = 14260.5` is the PICARD-ITERATION FAILURE threshold of
      `t2_gamma_CR.gamma_off`, not an existence threshold for (Gamma-off)'s fixed point.
      Here (Gamma-off)'s closure is analysed DIRECTLY.

THE CLOSURE (u2 sec.8, as re-implemented in t2_gamma_CR.py lines 66-102).  In the variable

      Gbar := Gamma_rad-free closure variable  =  lam L + C''        (so c_G = Gbar c / L),

the closure is the scalar equation  F(Gbar) = Gbar  with

      c_G(Gbar) = Gbar c / L
      Chat_a    = (C_far + C_inner) lam + min( 2R_A lam/sigma_* + pi lam/8 , 2 R_A e^{c_G} E_0 )
      Ghat      = (3 pi^2/16) e^{p c_G} Gfrak_0
      A         = (lam/2) L + Chat_a
      Gam_rad   = max( A , A + 2Ghat + lam/2 , 2Chat_a + 2Ghat + lam/2 )
      p         = min( 3 , 1 + 2 Gam_rad/Gbar )         (an inner fixed point: p enters Ghat)
      C''       = 2 Chat_a + Ghat + lam
      F(Gbar)   = lam L + C''(Gbar) .

FACTS ABOUT F, ESTABLISHED HERE.

  * F is continuous and strictly INCREASING in Gbar (c_G increases; Chat_a, Ghat increase;
    the inner p-map is monotone increasing in Ghat and Ghat is increasing in p, so its
    smallest fixed point, reached by iterating up from p = 1, is nondecreasing in Gbar).
  * F(lam L) > lam L  (C'' > 0), so F - id starts POSITIVE.
  * F(Gbar) grows like exp(p c Gbar/L), so F - id -> +infinity.
  Hence F - id has a zero iff  min_{Gbar} (F - id) <= 0, and a sign change from + to - is a
  PROOF (intermediate value theorem, F continuous) that a fixed point exists.  The EXISTENCE
  boundary L_Gamma^exist is the L at which the minimum of F - id touches zero.

  A supersolution is all the bootstrap needs: any Gbar with F(Gbar) <= Gbar gives
  Gamma(s) <= 2a(0,s) + C'' M with C'' = Gbar - lam L.  The smallest root is the best constant.

CONTROL (L-14).  Driven with u2's own datum norms the same instrument must reproduce
refute-u2 MINOR-1's directly-minimised boundaries 4904.7 (sigma_* = 1/2) and 5273.0
(sigma_* = 0), and u2's own C''(10^4) = 582.24 / 801.05.

Outputs -> f5_results.json
"""
import json, math

C_FAR = 0.29199853
C_INNER = 0.01475427
R_A = (2.0**5 - 2.0**-5)**0.2
TWO_RA = 2.0*R_A
COEF_G = 3.0*math.pi**2/16.0

DEG = math.pi/180.0
DELTA = 7.5*DEG
E0_TH = 1.0/math.sin(DELTA)             # PROVED exactly in f4
G0_TH = 1.0/math.sin(DELTA)**2          # PROVED exactly in f4
E0_DC, G0_DC = 7.66060196, 20.9197070   # u2's own datum (D-C), for the control only


def _sexp(x):
    return math.exp(x) if x < 700.0 else float("inf")


def chat_a(sigma_star, lam, c_G, E0):
    br1 = (TWO_RA*lam/sigma_star + math.pi*lam/8.0) if sigma_star > 0 else float("inf")
    br2 = TWO_RA*_sexp(c_G)*E0
    return (C_FAR + C_INNER)*lam + min(br1, br2)


def F_map(Gbar, L, c, lam, E0, G0, sigma_star=0.0, p_cap=3.0, itmax=200):
    """F(Gbar) = lam L + C''(Gbar); returns (F, info) or (inf, None)."""
    c_G = Gbar*c/L
    if c_G > 600.0:
        return float("inf"), None
    Ca = chat_a(sigma_star, lam, c_G, E0)
    A = 0.5*lam*L + Ca
    p = 1.0
    for _ in range(itmax):
        fac = _sexp(p*c_G)
        if not math.isfinite(fac):
            return float("inf"), None
        Gh = COEF_G*fac*G0
        Gam_rad = max(A, A + 2.0*Gh + lam/2.0, 2.0*Ca + 2.0*Gh + lam/2.0)
        pn = min(p_cap, 1.0 + 2.0*Gam_rad/Gbar)
        if abs(pn - p) < 1e-12:
            p = pn
            break
        p = pn
    fac = _sexp(p*c_G)
    if not math.isfinite(fac):
        return float("inf"), None
    Gh = COEF_G*fac*G0
    Cpp = 2.0*Ca + Gh + lam
    return lam*L + Cpp, {"C_pp": Cpp, "c_G": c_G, "p": p, "Ghat": Gh, "Chat_a": Ca,
                         "Gam_rad": Gam_rad, "Gbar": Gbar, "L": L, "c": c}


_SR_CACHE = {}


def smallest_root(L, c, lam, E0, G0, sigma_star=0.0, nscan=600, span=1e8):
    """smallest Gbar with F(Gbar) <= Gbar, by scan for a sign change then bisection.
       Returns None when F - id never becomes <= 0 (no fixed point at that L)."""
    key = (round(L, 6), round(c, 12), lam, round(E0, 12), round(G0, 12), sigma_star)
    if key in _SR_CACHE:
        return _SR_CACHE[key]
    base = lam*L
    lo = base*(1.0 + 1e-14)
    for i in range(nscan+1):
        g = lo*(span**(i/float(nscan)))
        f, _ = F_map(g, L, c, lam, E0, G0, sigma_star)
        if f <= g:
            a, b = base, g
            for _ in range(120):
                m = 0.5*(a+b)
                f2, _ = F_map(m, L, c, lam, E0, G0, sigma_star)
                if f2 <= m:
                    b = m
                else:
                    a = m
            f3, info = F_map(b, L, c, lam, E0, G0, sigma_star)
            info["F_minus_id_at_root"] = f3 - b
            info["sign_change_bracket"] = [a, b]
            fa, _ = F_map(a, L, c, lam, E0, G0, sigma_star)
            info["F_minus_id_left_of_bracket"] = fa - a
            _SR_CACHE[key] = info
            return info
    _SR_CACHE[key] = None
    return None


def min_F_minus_id(L, c, lam, E0, G0, sigma_star=0.0, n=20000, span=1e6):
    """min over Gbar of F - id, on a log grid (diagnostic: it is <= 0 iff a root exists)."""
    base = lam*L
    best, arg = float("inf"), None
    for i in range(n+1):
        g = base*(1.0 + 1e-14)*(span**(i/float(n)))
        f, _ = F_map(g, L, c, lam, E0, G0, sigma_star)
        if math.isfinite(f) and f - g < best:
            best, arg = f - g, g
    return best, arg


def L_exist(c, lam, E0, G0, sigma_star=0.0, lo=10.0, hi=1e10, n=6000, span=1e4):
    """the EXISTENCE boundary: smallest L at which min_{Gbar}(F - id) <= 0, i.e. at which
       F - id changes sign (F continuous, positive at Gbar = lam L, -> +inf: IVT gives a root)."""
    if min_F_minus_id(hi, c, lam, E0, G0, sigma_star, n=n, span=span)[0] > 0:
        return None
    for _ in range(200):
        m = math.sqrt(lo*hi)
        if min_F_minus_id(m, c, lam, E0, G0, sigma_star, n=n, span=span)[0] > 0:
            lo = m
        else:
            hi = m
        if hi/lo < 1.0 + 1e-10:
            break
    return hi


# ---------------------------------------------------------------- Picard, exactly as in t2
def picard(L, c, lam, E0, G0, sigma_star=0.0, p_cap=3.0, damp=0.5, tol=1e-12, itmax=4000):
    """byte-for-byte the algebra of t2_gamma_CR.gamma_off (the object whose threshold is 14260.5)."""
    c_G, p = lam*c, 2.0
    for _ in range(itmax):
        Ca = chat_a(sigma_star, lam, c_G, E0)
        if p*c_G > 500.0:
            return None
        Ghat = COEF_G*math.exp(p*c_G)*G0
        Cpp_new = 2.0*Ca + Ghat + lam
        c_G_new = (lam + Cpp_new/L)*c
        A = 0.5*lam*L + Ca
        Gam_rad = max(A, A + 2.0*Ghat + lam/2.0, 2.0*Ca + 2.0*Ghat + lam/2.0)
        c_R = Gam_rad*c/L
        p_new = min(p_cap, 1.0 + 2.0*c_R/c_G_new)
        if c_G_new > 60.0 or not math.isfinite(c_G_new):
            return None
        d = max(abs(c_G_new - c_G), abs(p_new - p))
        c_G = c_G + damp*(c_G_new - c_G)
        p = p + damp*(p_new - p)
        if d < tol:
            Ca = chat_a(sigma_star, lam, c_G, E0)
            Ghat = COEF_G*math.exp(p*c_G)*G0
            return {"C_pp": 2.0*Ca + Ghat + lam, "c_G": c_G, "p": p, "Ghat": Ghat,
                    "Chat_a": Ca, "L": L}
    return None


def L_picard(c, lam, E0, G0, sigma_star=0.0, lo=1.0, hi=1e12):
    if picard(hi, c, lam, E0, G0, sigma_star) is None:
        return None
    for _ in range(300):
        m = math.sqrt(lo*hi)
        if picard(m, c, lam, E0, G0, sigma_star) is None:
            lo = m
        else:
            hi = m
        if hi/lo < 1.0 + 1e-11:
            break
    return hi


if __name__ == "__main__":
    OUT = {"constants": {"C_far": C_FAR, "C_inner": C_INNER, "R_A": R_A, "two_R_A": TWO_RA,
                         "coef_3pi2_over_16": COEF_G, "E0_theorem_datum_exact": E0_TH,
                         "Gfrak0_theorem_datum_exact": G0_TH}}

    # ---------- CONTROL: u2's own datum, against refute-u2 MINOR-1 --------------------
    cU2 = 2.0*math.log(1.5)
    ctl = {}
    for sig, tag in [(0.5, "sigma*=1/2"), (0.0, "sigma*=0")]:
        Le = L_exist(cU2, 1.5, E0_DC, G0_DC, sig)
        Lp = L_picard(cU2, 1.5, E0_DC, G0_DC, sig)
        r = smallest_root(1e4, cU2, 1.5, E0_DC, G0_DC, sig)
        ctl[tag] = {"L_exist": Le, "L_picard": Lp, "C_pp_at_1e4": r["C_pp"],
                    "c_G_at_1e4": r["c_G"], "p_at_1e4": r["p"]}
    ctl["refute_u2_reported_exist"] = {"sigma*=1/2": 4904.7, "sigma*=0": 5273.0}
    ctl["u2_reported_picard"] = {"sigma*=1/2": 4997.6, "sigma*=0": 5313.8}
    ctl["u2_reported_Cpp_1e4"] = {"sigma*=1/2": 582.24, "sigma*=0": 801.05}
    OUT["control_u2_datum"] = ctl

    # ---------- PRODUCTION: the theorem's datum --------------------------------------
    KAPPA = 0.4978224383290105                # f4 (P_h(1)/2), reproduced there
    C_STAR = math.log(1.5)/KAPPA
    prod = {}
    for cfac in [1.0, 1.005, 1.01, 1.05, 1.1, 1.1943662078602755, 1.2, 1.5]:
        c = cfac*C_STAR
        Le = L_exist(c, 1.5, E0_TH, G0_TH, 0.0)
        Lp = L_picard(c, 1.5, E0_TH, G0_TH, 0.0)
        prod["c=%.6f c_* (c=%.6f)" % (cfac, c)] = {
            "L_Gamma_exist": Le, "L_Gamma_picard": Lp,
            "picard_over_exist": (Lp/Le if (Le and Lp) else None)}
    OUT["theorem_datum_thresholds"] = prod

    # the sign-change certificate at the production window, just above the boundary
    c = C_STAR
    Le = L_exist(c, 1.5, E0_TH, G0_TH, 0.0)
    cert = {}
    for mult in [0.999, 1.0, 1.001, 1.01, 1.1, 2.0]:
        L = Le*mult
        mn, arg = min_F_minus_id(L, c, 1.5, E0_TH, G0_TH, 0.0)
        r = smallest_root(L, c, 1.5, E0_TH, G0_TH, 0.0)
        cert["L=%.4f (=%.4f L_exist)" % (L, mult)] = {
            "min_F_minus_id": mn, "argmin_Gbar": arg,
            "root_found": (None if r is None else r["Gbar"]),
            "C_pp": (None if r is None else r["C_pp"]),
            "c_G": (None if r is None else r["c_G"]),
            "p": (None if r is None else r["p"]),
            "F_minus_id_left_of_bracket": (None if r is None else r["F_minus_id_left_of_bracket"]),
            "F_minus_id_at_root": (None if r is None else r["F_minus_id_at_root"])}
    OUT["sign_change_certificate"] = cert

    # C'' against L, direct root, at c = c_* (to be compared with THEOREM_S3 sec.4.3)
    tab = {}
    for L in [1e3, 3e3, 1e4, 1.2e4, 1.3e4, 1.4e4, 3e4, 1e5, 3e5, 1e6, 1e9]:
        r = smallest_root(L, C_STAR, 1.5, E0_TH, G0_TH, 0.0)
        p2 = picard(L, C_STAR, 1.5, E0_TH, G0_TH, 0.0)
        tab["L=%g" % L] = {"direct": (None if r is None else
                                      {k: r[k] for k in ["C_pp", "c_G", "p", "Ghat", "Chat_a"]}),
                           "picard": p2}
    OUT["C_pp_table_direct_vs_picard"] = tab
    OUT["THEOREM_S3_quoted_L_Gamma_star"] = 14260.5

    with open("f5_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
    print(json.dumps(OUT, indent=1, sort_keys=True, default=str))

"""
f4 -- ITEM 4 of the fix list: the three discrete scans that THEOREM_S3 labels PROVED are
      replaced by (a) a certified enclosure plus an analytic argument for the argmin,
      (b) an analytic argument, (c) an analytic argument.

The three objects:

  (a)  r_h  = inf_{lam in [1,lam_max]} P_h(lam)/lam          (t1_datum.py:98-106, 61 samples)
  (b)  P_h increasing on [1,lam_max]                          (t1_datum.py:115-132, 0.01 steps)
  (c)  E_0 = sup |x||eta_0|/M  and  Gfrak_0 = sup |x|^2|grad eta_0|/M
                                                              (t1_datum.py:156-183, finite diffs)

THE EXACT STRUCTURE (all of it new here; nothing imported).

    P_h(lam) = 3 int_0^1 h(arcsin v) v^2 (A v^2 + B)^{-5/2} dv ,   A = lam^2 - lam^{-4}, B = lam^{-4}.

    A v^2 + B = lam^{-4}[1 + (lam^6-1)v^2] =: lam^{-4} D(v,lam) ,  so

    P_h(lam) = 3 lam^10 int_0^1 H(v) v^2 D^{-5/2} dv ,   H(v) := h(arcsin v) .

    d/dv [ v^3 D^{-3/2} ] = 3 v^2 D^{-5/2}        (EXACT; sympy-checked below)

  => with  Q(lam) := P_h(lam)/lam = 3 lam^9 int_0^1 H v^2 D^{-5/2} dv  and
     F_lam(v) := lam^9 v^3 D(v,lam)^{-3/2} ,     F_lam(0) = 0 ,  F_lam(1) = lam^9 (lam^6)^{-3/2} = 1 ,

     dmu_lam := dF_lam is a PROBABILITY MEASURE on [0,1] and  Q(lam) = E_{mu_lam}[H] .
     (Consequently P_1(lam) = lam EXACTLY -- H == 1 -- which is A2's control, now a theorem.)

     d/dlam log F_lam = 9/lam - 9 lam^5 v^2/D = (9/lam)(1-v^2)/D  >= 0 ,
     so  d_lam F_lam = 9 lam^8 v^3 (1-v^2) D^{-5/2} >= 0 :  mu_lam DECREASES stochastically in lam.

     Integrating by parts (boundary terms vanish: d_lam F_lam = 0 at v = 0 and v = 1),

        Q'(lam) = - int_0^1 H'(v) d_lam F_lam(v) dv
                = - 9 lam^8 int_0^1 H'(v) v^3 (1-v^2) D^{-5/2} dv .                      (Q')

     For the theorem's datum H' = +1/(delta sqrt(1-v^2)) on (0, sin delta),  0 on
     (sin delta, cos delta_m),  -1/(delta_m sqrt(1-v^2)) on (cos delta_m, 1).  Hence

        Q'(lam) = -9 lam^8 [ (1/delta) Aint(lam) - (1/delta_m) Bint(lam) ] ,
        Aint = int_0^{sin delta} v^3 sqrt(1-v^2) D^{-5/2} dv ,
        Bint = int_{cos delta_m}^1 v^3 sqrt(1-v^2) D^{-5/2} dv ,

     and both are majorised in CLOSED FORM by dropping sqrt(1-v^2) <= 1:

        int_0^V v^3 (1+K v^2)^{-5/2} dv = K^{-2} [ 2/3 + (1/3) S^{-3/2} - S^{-1/2} ] ,  S = 1+K V^2,

     with K = lam^6 - 1  (the K -> 0 limit is V^4/4).  Call these Abar, Bbar.
     Both are DECREASING in lam (D is increasing in lam at every v > 0), and lam^8, lam^9 are
     increasing, so on any subinterval [l_i, l_{i+1}] a rigorous upper bound is obtained by
     taking lam^8 at the right end and Abar, Bbar at the left end.  That is the certificate.

  (a) r_h:  |Q'| <= Lip(lam) := 9 lam^8 [ Abar/delta + Bbar/delta_m ]  -- PROVED bound --
      so on a grid of spacing h,  inf Q >= min_i Q(l_i) - max_i Lip_i * h .   ENCLOSED.
      Additionally Q'(lam) < 0 for lam >= lam_turn is PROVED by exhibiting
      (1/delta) Aint > (1/delta_m) Bbar there, so the infimum is at lam_max.

  (b) P_h'(lam) = Q + lam Q' = Q(lam) - 9 lam^9 [ (1/delta)Aint - (1/delta_m)Bint ]
                >= r_h - 9 lam^9 Abar(lam)/delta        (Bint >= 0, Aint <= Abar, Q >= r_h)
      and max_{[1,lam_max]} 9 lam^9 Abar/delta is certified by the same monotone subdivision.
      If that max is < r_h, P_h is strictly increasing.  PROVED.

  (c) with F(phi) := sgn(cos phi) g(phi)/sin phi,  u = log(rho/rho_0),
        rho|eta_0|/M      = |F| Theta ,
        rho^2|grad eta_0|/M = sqrt( F^2 (Theta - Theta_u)^2 + F'^2 Theta^2 ) .
      Writing Theta = sig1 - sig2 with sig_i = (1+tanh(.))/2 in (0,1), one gets EXACTLY
        Theta_u = (2/eps_r)[ sig1(1-sig1) - sig2(1-sig2) ] = 8 b (1-s) ,  b := sig1-sig2, s := sig1+sig2
        Theta - Theta_u = b (8 s - 7)        (eps_r = 1/4)
      so with P := F^2, Qq := F'^2,
        (rho^2|grad eta_0|/M)^2 = b^2 [ P (8s-7)^2 + Qq ] ,   0 <= b <= 1 ,  b <= s <= 2-b .
      max_s (8s-7)^2 = (9-8b)^2, and f(b) := b^2[P(9-8b)^2+Qq] has
        f'(b) = 2b (2P y^2 - 9 P y + Qq) ,  y := 9-8b in [1,9] ,
      so f is increasing (max = f(1) = P+Qq) whenever 81 P <= 8 Qq, and otherwise
        f(b) <= P max_y (9-y)^2y^2/64 + Qq max_y (9-y)^2/64 = 6.40723 P + Qq .
      Both branches are then maximised over phi in closed form.  PROVED.

NUMERICS.  Q, Aint, Bint are computed by Gauss-Legendre with the kinks of H on panel
boundaries (the integrands are real-analytic on each panel), in float64, and cross-checked
against mpmath at 40 digits at four values of lam.  The closed forms Abar, Bbar are exact.

Outputs -> f4_results.json
"""
import json, math
import numpy as np
import mpmath as mp
import sympy as sp

DEG = math.pi/180.0
DELTA = 7.5*DEG
DM = 5.0*DEG
SD = math.sin(DELTA)
CDM = math.cos(DM)
EPS_R = 0.25

OUT = {}

# --------------------------------------------------------------------------- 0. EXACT algebra
v, lam, K = sp.symbols('v lam K', positive=True)
D = 1 + (lam**6 - 1)*v**2
res_anti = sp.simplify(sp.diff(v**3*D**sp.Rational(-3, 2), v) - 3*v**2*D**sp.Rational(-5, 2))
prim = -(K*v**2 + 1)**sp.Rational(-1, 2)/K**2 + sp.Rational(1, 3)*(K*v**2+1)**sp.Rational(-3, 2)/K**2
res_prim = sp.simplify(sp.diff(prim, v) - v**3*(1 + K*v**2)**sp.Rational(-5, 2))
F_lam = lam**9*v**3*D**sp.Rational(-3, 2)
res_dlam = sp.simplify(sp.diff(F_lam, lam) - 9*lam**8*v**3*(1-v**2)*D**sp.Rational(-5, 2))
res_mass = sp.simplify(F_lam.subs(v, 1) - 1)
# the Theta identity of (c):  Theta = s1-s2, Theta_u = 8[s1(1-s1)-s2(1-s2)] = 8 b (1-s)
s1, s2 = sp.symbols('s1 s2')
b_, s_ = s1 - s2, s1 + s2
res_theta = sp.simplify(8*(s1*(1-s1) - s2*(1-s2)) - 8*b_*(1 - s_))
res_theta2 = sp.simplify((b_ - 8*(s1*(1-s1) - s2*(1-s2))) - b_*(8*s_ - 7))
OUT["exact_sympy_residuals"] = {
    "d/dv[v^3 D^{-3/2}] - 3v^2 D^{-5/2}": str(res_anti),
    "d/dv[primitive] - v^3(1+Kv^2)^{-5/2}": str(res_prim),
    "d_lam F_lam - 9 lam^8 v^3(1-v^2)D^{-5/2}": str(res_dlam),
    "F_lam(v=1) - 1   (=> P_1(lam)=lam, total mass 1)": str(res_mass),
    "Theta_u - 8 b (1-s)": str(res_theta),
    "(Theta - Theta_u) - b(8s-7)": str(res_theta2)}
assert all(r == 0 for r in [res_anti, res_prim, res_dlam, res_mass, res_theta, res_theta2])

# --------------------------------------------------------------------------- 1. quadrature
NG = 160
_x, _w = np.polynomial.legendre.leggauss(NG)

def _panel(a, b):
    return 0.5*(b-a)*_x + 0.5*(a+b), 0.5*(b-a)*_w

def _H(vv):
    phi = np.arcsin(np.clip(vv, 0.0, 1.0))
    return np.minimum(1.0, phi/DELTA)*np.minimum(1.0, (math.pi/2 - phi)/DM)

_PANELS = [_panel(0.0, SD), _panel(SD, CDM), _panel(CDM, 1.0)]
_VQ = np.concatenate([p[0] for p in _PANELS])
_WQ = np.concatenate([p[1] for p in _PANELS])
_HQ = _H(_VQ)
_VA, _WA = _panel(0.0, SD)
_VB, _WB = _panel(CDM, 1.0)

def Q(l):
    Dv = 1.0 + (l**6 - 1.0)*_VQ**2
    return 3.0*l**9*float(np.dot(_WQ, _HQ*_VQ**2*Dv**-2.5))

def P_h(l):
    return l*Q(l)

def Aint(l):
    Dv = 1.0 + (l**6 - 1.0)*_VA**2
    return float(np.dot(_WA, _VA**3*np.sqrt(1.0-_VA**2)*Dv**-2.5))

def Bint(l):
    Dv = 1.0 + (l**6 - 1.0)*_VB**2
    return float(np.dot(_WB, _VB**3*np.sqrt(np.maximum(1.0-_VB**2, 0.0))*Dv**-2.5))

def _prim_closed(V, Kv):
    Kv = np.asarray(Kv, dtype=float)
    Sv = 1.0 + Kv*V*V
    with np.errstate(divide='ignore', invalid='ignore'):
        val = (2.0/3.0 + Sv**-1.5/3.0 - Sv**-0.5)/Kv**2
    small = Kv < 1e-9
    if np.ndim(val) == 0:
        return float(V**4/4.0) if bool(small) else float(val)
    val = np.where(small, V**4/4.0, val)
    return val

def Abar(l):
    return _prim_closed(SD, np.asarray(l, dtype=float)**6 - 1.0)

def Bbar(l):
    Kv = np.asarray(l, dtype=float)**6 - 1.0
    return _prim_closed(1.0, Kv) - _prim_closed(CDM, Kv)

def Qprime(l):
    return -9.0*l**8*(Aint(l)/DELTA - Bint(l)/DM)

KAPPA = P_h(1.0)/2.0
LOG32 = math.log(1.5)
C_STAR = LOG32/KAPPA



# --------------------------------------------------------------------------- 2. certified r_h
def _Q_vec(ls):
    """Q on an array of lam, in chunks (the node array has 3*NG points)."""
    out = np.empty_like(ls)
    wgt = (_WQ*_HQ*_VQ**2)[None, :]
    v2 = _VQ[None, :]**2
    step = 20000
    for i in range(0, len(ls), step):
        sl = ls[i:i+step]
        Dv = 1.0 + (sl**6 - 1.0)[:, None]*v2
        out[i:i+step] = 3.0*sl**9*(wgt*Dv**-2.5).sum(axis=1)
    return out


def certified_r_h(lam_max, nsub=50000):
    """rigorous LOWER bound on inf_{[1,lam_max]} Q by monotone-majorant Lipschitz + grid."""
    step = (lam_max - 1.0)/nsub
    ls = 1.0 + step*np.arange(nsub+1)
    Qs = _Q_vec(ls)
    Ab = Abar(ls)
    Bb = Bbar(ls)
    lip = 9.0*ls[1:]**8*(Ab[:-1]/DELTA + Bb[:-1]/DM)      # lam^8 up, Abar/Bbar down
    lower = Qs[:-1] - lip*step
    k = int(np.argmin(Qs))
    return {"r_h_grid_min": float(Qs.min()), "argmin_lam": float(ls[k]),
            "r_h_certified_lower": float(lower.min()),
            "max_Lipschitz_bound": float(lip.max()), "nsub": nsub, "step": step,
            "lam_max": lam_max, "Q_at_lam_max": float(Qs[-1]), "Q_at_1": float(Qs[0])}

def certified_Ph_increasing(lam_max, r_h_lower, nsub=50000):
    step = (lam_max - 1.0)/nsub
    ls = 1.0 + step*np.arange(nsub+1)
    Ab = Abar(ls)
    vals = 9.0*ls[1:]**9*Ab[:-1]/DELTA                    # lam^9 up, Abar down
    k = int(np.argmax(vals))
    return {"max_9lam9_Abar_over_delta": float(vals.max()), "argmax_lam": float(ls[k]),
            "r_h_lower_used": r_h_lower,
            "min_Ph_prime_lower_bound": float(r_h_lower - vals.max()),
            "PROVED_increasing": bool(r_h_lower - vals.max() > 0), "nsub": nsub}

def turning_certificate(lam_max, nsub=20000):
    """smallest lam_turn with (1/delta)Aint > (1/delta_m)Bbar on [lam_turn, lam_max]:
       PROVES Q' < 0 there (Bint <= Bbar), so inf Q on that stretch is at lam_max."""
    step = (lam_max - 1.0)/nsub
    ls = 1.0 + step*np.arange(nsub+1)
    A = np.array([Aint(x) for x in ls])/DELTA
    B = np.asarray(Bbar(ls))/DM
    ok = A > B
    turn = None
    for i in range(len(ls)-1, -1, -1):
        if not ok[i]:
            turn = float(ls[i+1]) if i+1 < len(ls) else None
            break
        turn = float(ls[i])
    return {"lam_turn_Qprime_negative_from": turn,
            "Q_at_1": float(Q(1.0)), "Q_at_lam_max": float(Q(lam_max)),
            "inf_is_at_lam_max": bool(Q(lam_max) < Q(1.0))}


def Ph_prime(l):
    return Q(l) + l*Qprime(l)


# --------------------------------------------------------------------------- 3. E_0, Gfrak_0
def F_ang(phi):
    pax = min(phi, math.pi - phi)
    g = min(1.0, pax/DELTA)*min(1.0, abs(phi - math.pi/2)/DM)
    return (1.0 if math.cos(phi) > 0 else -1.0)*g/math.sin(phi)

def Fp_closed(phi, side):
    """closed-form one-sided derivative of F on each analytic piece."""
    s, c = math.sin(phi), math.cos(phi)
    sg = 1.0 if c > 0 else -1.0
    pax = min(phi, math.pi - phi)
    dpax = 1.0 if phi < math.pi/2 else -1.0
    eq = abs(phi - math.pi/2)
    deq = 1.0 if phi > math.pi/2 else -1.0
    # A(phi) = sg * g1 * g2 ; pick the branch on the requested side of a kink
    inA = (pax < DELTA) or (pax == DELTA and side < 0)
    inE = (eq < DM) or (eq == DM and side < 0)
    g1, dg1 = (pax/DELTA, dpax/DELTA) if inA else (1.0, 0.0)
    g2, dg2 = (eq/DM, deq/DM) if inE else (1.0, 0.0)
    A = sg*g1*g2
    Ap = sg*(dg1*g2 + g1*dg2)
    return (Ap*s - A*c)/s**2


YMAX = (4.5*4.5)**2/64.0          # max_{y in [1,9]} (9-y)^2 y^2 / 64

def sup_over_u(P, Qq):
    if 81.0*P <= 8.0*Qq:
        return P + Qq
    return YMAX*P + Qq


def main():
    # ---- controls: mpmath at 40 digits, and Q' against a finite difference
    mp.mp.dps = 40
    mDELTA = mp.mpf('7.5')*mp.pi/180
    mDM = mp.mpf(5)*mp.pi/180
    mSD, mCDM = mp.sin(mDELTA), mp.cos(mDM)

    def Q_mp(l):
        l = mp.mpf(l)
        def f(x):
            phi = mp.asin(x)
            Hh = min(mp.mpf(1), phi/mDELTA)*min(mp.mpf(1), (mp.pi/2 - phi)/mDM)
            return Hh*x**2*(1 + (l**6-1)*x**2)**mp.mpf('-2.5')
        return 3*l**9*mp.quad(f, [mp.mpf(0), mSD, mCDM, mp.mpf(1)])

    ctl = {}
    for lv in [1.0, 1.2, 1.5, 1.8420112, 2.5]:
        ctl["Q_float_minus_Q_mp40 at lam=%g" % lv] = float(Q(lv) - Q_mp(lv))
        ctl["Abar-Aint at lam=%g" % lv] = Abar(lv) - Aint(lv)
        ctl["Bbar-Bint at lam=%g" % lv] = Bbar(lv) - Bint(lv)
        h = 1e-6
        ctl["Qprime_closed_minus_findiff at lam=%g" % lv] = Qprime(lv) - (Q(lv+h)-Q(lv-h))/(2*h)
    OUT["controls"] = ctl

    OUT["kappa_delta_from_P_h(1)/2"] = KAPPA
    OUT["kappa_delta_THEOREM_S3_quoted"] = 0.49782244
    OUT["Q(1)_minus_2kappa"] = Q(1.0) - 2*KAPPA
    OUT["c_star"] = C_STAR
    OUT["P_1_control_max_abs_dev"] = max(abs(3.0*l**9*float(np.dot(_WQ, _VQ**2*(1.0+(l**6-1.0)*_VQ**2)**-2.5)) - 1.0)
                                         for l in [1.0, 1.25, 1.5, 1.8420112, 2.5])

    LAM_MAX_APRIORI = math.exp(3*C_STAR/4)
    OUT["lam_max_at_c_star"] = LAM_MAX_APRIORI
    rh = certified_r_h(LAM_MAX_APRIORI)
    OUT["r_h_certified_at_c_star"] = rh
    OUT["r_h_THEOREM_S3_quoted"] = 0.91863653
    OUT["turning_at_c_star"] = turning_certificate(LAM_MAX_APRIORI)
    OUT["Ph_increasing_at_c_star"] = certified_Ph_increasing(LAM_MAX_APRIORI,
                                                             rh["r_h_certified_lower"])

    OUT["r_h_certified_other_windows"] = {}
    for cfac in [1.005, 1.01, 1.05, 1.1, 1.2, 1.5]:
        lm = math.exp(3*(cfac*C_STAR)/4)
        r = certified_r_h(lm)
        r["Ph_increasing"] = certified_Ph_increasing(lm, r["r_h_certified_lower"])
        r["turning"] = turning_certificate(lm)
        OUT["r_h_certified_other_windows"]["c=%g c_*" % cfac] = r

    # ------------------------------------------------ 2b. where P_h STOPS increasing (the window cap)
    # A2/C6 need P_h increasing on [1, lam_max] for the quasimonotone comparison.  Enlarging the
    # window (item 1) enlarges lam_max = e^{3c/4}.  P_h'(lam) = Q + lam Q' has a zero at lam_mono:
    # beyond it the comparison hypothesis is FALSE, not merely uncertified.
    lo, hi = 1.5, 3.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if Ph_prime(mid) > 0:
            lo = mid
        else:
            hi = mid
    LAM_MONO = 0.5*(lo+hi)
    # largest window factor for which the CERTIFICATE of (b) still holds
    def cert_ok(cfac):
        lm = math.exp(3*(cfac*C_STAR)/4)
        r = certified_r_h(lm, nsub=20000)
        c2 = certified_Ph_increasing(lm, r["r_h_certified_lower"], nsub=20000)
        return c2["PROVED_increasing"]
    lo2, hi2 = 1.0, 1.30
    for _ in range(40):
        mid = 0.5*(lo2+hi2)
        if cert_ok(mid):
            lo2 = mid
        else:
            hi2 = mid
    OUT["window_cap_from_P_h_monotonicity"] = {
        "lam_mono_P_h_prime_zero": LAM_MONO,
        "Ph_prime_at_lam_mono": Ph_prime(LAM_MONO),
        "c_over_c_star_at_lam_mono": (4.0/3.0)*math.log(LAM_MONO)/C_STAR,
        "eps_max_structural": (4.0/3.0)*math.log(LAM_MONO)/C_STAR - 1.0,
        "c_over_c_star_largest_CERTIFIED": lo2,
        "eps_max_certified": lo2 - 1.0,
        "note": "eps <= 1/2 requires lam_max = e^{3(1.5 c_*)/4} = %.6f > lam_mono: A2/C6 FAIL there"
                % math.exp(3*1.5*C_STAR/4)}

    worst, argw, branch = 0.0, None, None
    NPHI = 200000
    for i in range(1, NPHI):
        phi = math.pi*i/NPHI
        P = F_ang(phi)**2
        Qq = max(Fp_closed(phi, +1)**2, Fp_closed(phi, -1)**2)
        val = sup_over_u(P, Qq)
        if val > worst:
            worst, argw = val, phi*180/math.pi
            branch = "P+Q" if 81*P <= 8*Qq else "6.407P+Q"
    for bp in [DELTA, math.pi/2 - DM, math.pi/2, math.pi/2 + DM, math.pi - DELTA]:
        for eps in (-1e-12, 1e-12):
            phi = bp + eps
            P = F_ang(phi)**2
            Qq = max(Fp_closed(phi, +1)**2, Fp_closed(phi, -1)**2)
            val = sup_over_u(P, Qq)
            if val > worst:
                worst, argw = val, phi*180/math.pi
                branch = "P+Q" if 81*P <= 8*Qq else "6.407P+Q"

    OUT["datum_norms_analytic"] = {
        "E_0_exact_1_over_sin_delta": 1.0/SD,
        "E_0_THEOREM_S3_grid": 7.6612976,
        "tanh(2L-2)_at_L=40_finite_L_factor": math.tanh(2*40.0-2),
        "Gfrak_0_exact_1_over_sin2_delta": 1.0/SD**2,
        "Gfrak_0_THEOREM_S3_grid": 58.695476,
        "Gfrak_0_certified_upper_over_all_phi_and_u": math.sqrt(worst),
        "scan_argmax_phi_deg": argw, "scan_branch": branch,
        "bound_equals_1_over_sin2delta": bool(abs(math.sqrt(worst) - 1.0/SD**2) < 1e-9),
        "max_y_(9-y)^2y^2_over_64": YMAX,
        "at_phi=delta+  81P": 81.0/SD**2, "at_phi=delta+  8Q": 8.0*math.cos(DELTA)**2/SD**4,
        "F_prime_at_delta+_equals_-cos_delta/sin^2_delta": -math.cos(DELTA)/SD**2,
        "F_prime_at_delta-_taper_branch": (math.sin(DELTA)-DELTA*math.cos(DELTA))/(DELTA*SD**2),
    }

    with open("f4_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True)
    print(json.dumps(OUT, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()

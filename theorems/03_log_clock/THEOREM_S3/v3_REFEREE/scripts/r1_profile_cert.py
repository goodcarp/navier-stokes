"""
r1_profile_cert -- CERTIFIED enclosures of W_sigma^{(k)}(phi), k = 0..4, and certified magnitude
bounds on |W_sigma^{(m)}|, m = 0..5, valid on any phi-interval inside the Theorem V.4' ball's
angular range.

WHY.  fix5/x1 evaluates the mollified profile by a Gauss-Legendre convolution and then takes a
GRID maximum of the Taylor-remainder quotients; the theorem states vartheta = 2 x (that grid
maximum).  A grid maximum is a LOWER bound on a supremum and twice a lower bound is not an upper
bound.  To certify an upper bound one needs enclosures, and the only object in the chain that is
not elementary is the mollified profile.  This file supplies them.

THE DECOMPOSITION.  Wtil is continuous, 2 pi periodic, even at 0 and pi, and piecewise
real-analytic with corners only at the kink images.  In the window that matters here, namely
t in [-0.125, 1.25], there is exactly ONE corner, at t = delta:

    t in (-delta, delta) :   Wtil(t) = t/(delta sin t)                 =: Q(t)   (even, analytic)
    t in (delta, pi/2 - delta_m) : Wtil(t) = 1/sin t                   =: P(t)

Distributionally, with J_i := [Wtil^{(i)}](delta) = P^{(i)}(delta) - Q^{(i)}(delta),

    Wtil^{(m)} = {Wtil^{(m)}}  +  sum_{i=1}^{m-1} J_i delta_{t=delta}^{(m-1-i)} ,

so, since W_sigma = Wtil * chi_sigma and (g * chi)^{(m)} = g^{(m)} * chi,

    W_sigma^{(m)}(phi) = sum_{i=1}^{m-1} J_i chi_sigma^{(m-1-i)}(phi - delta)
                         + int {Wtil^{(m)}}(t) chi_sigma(phi - t) dt .

Both pieces are enclosed:

  * the jump sum is closed form;
  * the integral is banded on a fixed t-grid of width sigma/32.  On the P side the derivatives
    of 1/sin alternate in sign with monotone magnitude on (0, pi/2) -- proved below from the
    closed forms, whose numerators are polynomials in cos t with non-negative coefficients -- so
    the band's min and max are its endpoint values.  On the Q side a Cauchy estimate on the disc
    |t| <= 1/2 bounds every derivative, giving a rigorous band pad.  Gaussian band masses come
    from the error function.  The tail outside [-0.125, 1.25] is bounded by
    sup|{Wtil^{(m)}}| x Phi(-16.7), which is below 1e-50 at every order used.

A phi-INTERVAL enclosure is then (value at the midpoint) +/- (halfwidth) x (bound on the next
order), with the bound taken over that same interval.

CONTROL (L-98).  Every enclosure is checked to contain fix5/fx_profile.py's Gauss-Legendre
value, an instrument built on a different principle (fixed-node quadrature in the integration
variable).  A containment failure is reported, not hidden.
"""
import math
import os
import sys

import numpy as np
import sympy as sp
from scipy.special import ndtr

HERE = os.path.dirname(os.path.abspath(__file__))

DEG = math.pi/180.0
DELTA = 7.5*DEG
DM = 5.0*DEG
PHI0 = 30.0*DEG
SIGMA = 0.02
WMAX = 1.0/math.sin(DELTA)

TLO, THI = -0.125, 1.25
WB = SIGMA/32.0                     # band width

# ------------------------------------------------------------------ closed forms, order 0..7
_t = sp.symbols('t', positive=True)
_P = 1/sp.sin(_t)
_Q = _t/(sp.Rational(1, 1)*sp.sin(_t))          # divided by delta below
P_EXPR = [sp.simplify(sp.diff(_P, _t, m)) for m in range(8)]
Q_EXPR = [sp.diff(_Q, _t, m) for m in range(8)]
P_FUN = [sp.lambdify(_t, e, "numpy") for e in P_EXPR]
Q_FUN = [sp.lambdify(_t, e, "numpy") for e in Q_EXPR]


def Pd(m, t):
    return np.asarray(P_FUN[m](np.asarray(t, dtype=float)), dtype=float)


def Qd(m, t):
    """derivatives of Q(t) = t/(delta sin t); the sympy expression carries delta = 1"""
    t = np.asarray(t, dtype=float)
    small = np.abs(t) < 1e-7
    tt = np.where(small, 1e-7, t)
    v = np.asarray(Q_FUN[m](tt), dtype=float)/DELTA
    if np.any(small):
        # the function is even and analytic at 0; use a tiny offset value, the error is O(1e-14)
        v = np.where(small, np.asarray(Q_FUN[m](np.full(tt.shape, 1e-7)), dtype=float)/DELTA, v)
    return v


def sign_alternation_certificate():
    """P^{(m)}(t) = eps_m N_m(cos t)/sin^{m+1} t with N_m having non-negative coefficients and
    eps_m = (-1)^m; hence |P^{(m)}| is strictly decreasing on (0, pi/2).  Checked symbolically."""
    c = sp.symbols('c')
    rows = []
    for m in range(8):
        e = sp.simplify(P_EXPR[m]*sp.sin(_t)**(m+1))
        e = sp.simplify(e.rewrite(sp.cos))
        e = sp.simplify(sp.expand_trig(e))
        # substitute cos t -> c using sin^2 = 1 - c^2 wherever it survives
        e2 = sp.simplify(e.subs(sp.sin(_t)**2, 1 - sp.cos(_t)**2).subs(sp.cos(_t), c))
        e2 = sp.expand(e2)
        pol = sp.Poly(e2*(-1)**m, c) if e2.free_symbols else None
        ok = None
        if pol is not None:
            ok = all(float(co) >= 0 for co in pol.all_coeffs())
        rows.append({"m": m, "numerator_times_(-1)^m": str(sp.simplify(e2*(-1)**m)),
                     "all_coefficients_nonnegative": ok})
    return rows


# ------------------------------------------------------------------ the Gaussian and its tails
def chi(y, k=0, sg=SIGMA):
    y = np.asarray(y, dtype=float)
    ch = np.exp(-0.5*(y/sg)**2)/(sg*math.sqrt(2.0*math.pi))
    s2 = sg*sg
    if k == 0:
        return ch
    if k == 1:
        return -(y/s2)*ch
    if k == 2:
        return (y*y/s2 - 1.0)/s2*ch
    if k == 3:
        return (y/s2**2)*(3.0 - y*y/s2)*ch
    if k == 4:
        return (3.0/s2**2 - 6.0*y*y/s2**3 + y**4/s2**4)*ch
    raise ValueError(k)


_HERM_ROOTS = {0: [], 1: [0.0], 2: [-1.0, 1.0], 3: [-math.sqrt(3.0), 0.0, math.sqrt(3.0)],
               4: [-math.sqrt(3.0 + math.sqrt(6.0)), -math.sqrt(3.0 - math.sqrt(6.0)),
                   math.sqrt(3.0 - math.sqrt(6.0)), math.sqrt(3.0 + math.sqrt(6.0))],
               5: [-2.856970, -1.355626, 0.0, 1.355626, 2.856970]}


def chi_absmax_on(ylo, yhi, k):
    """max |chi^{(k)}(y)| for y in [ylo, yhi]; exact by including the critical points, which are
    the roots of the Hermite polynomial He_{k+1}(y/sigma)."""
    ylo = np.asarray(ylo, dtype=float)
    yhi = np.asarray(yhi, dtype=float)
    m = np.maximum(np.abs(chi(ylo, k)), np.abs(chi(yhi, k)))
    for r in _HERM_ROOTS[k+1]:
        y = r*SIGMA
        inside = (ylo <= y) & (y <= yhi)
        if np.any(inside):
            m = np.where(inside, np.maximum(m, np.abs(chi(np.full(ylo.shape, y), k))), m)
    return m


# ------------------------------------------------------------------ the fixed t-band structure
def _build_bands():
    nl = int(math.ceil((DELTA - TLO)/WB))
    left = DELTA - WB*np.arange(nl, -1, -1)          # ends exactly at DELTA
    nr = int(math.ceil((THI - DELTA)/WB))
    right = DELTA + WB*np.arange(1, nr+1)
    edges = np.concatenate([left, right])
    return edges


EDGES = _build_bands()
TL = EDGES[:-1]
TR = EDGES[1:]
TMID = 0.5*(TL + TR)
NB = TL.size
IS_LEFT = TR <= DELTA + 1e-15

# Cauchy bound for Q on |t| <= DELTA: Q analytic on |t| < pi, take R = 1/2
_R = 0.5
_MQ = _R/(DELTA*abs(math.sin(_R)))            # >= sup_{|t| = R} |t/(delta sin t)| (min |sin| on
_MQ = max(_MQ, _R/(DELTA*math.sinh(_R)))      # the circle is at the real axis for this radius)
QCAUCHY = [math.factorial(m)*_MQ/(_R - DELTA)**(m+1) for m in range(9)]

# band minima and maxima of {Wtil^{(m)}}
GMIN = {}
GMAX = {}
for _m in range(7):
    vlo = np.empty(NB)
    vhi = np.empty(NB)
    # right (P) side: sign constant, |.| decreasing, so endpoints bracket
    pr_l = Pd(_m, np.maximum(TL, DELTA))
    pr_r = Pd(_m, np.maximum(TR, DELTA))
    vlo_p = np.minimum(pr_l, pr_r)
    vhi_p = np.maximum(pr_l, pr_r)
    # left (Q) side: subgrid plus a Cauchy pad
    NSUB = 9
    frac = np.linspace(0.0, 1.0, NSUB)
    ts = TL[:, None] + (TR - TL)[:, None]*frac[None, :]
    ts = np.minimum(ts, DELTA)
    qv = Qd(_m, ts)
    pad = 0.5*(WB/(NSUB - 1))*QCAUCHY[_m+1]
    vlo_q = qv.min(axis=1) - pad
    vhi_q = qv.max(axis=1) + pad
    GMIN[_m] = np.where(IS_LEFT, vlo_q, vlo_p)
    GMAX[_m] = np.where(IS_LEFT, vhi_q, vhi_p)

GABS = {m: np.maximum(np.abs(GMIN[m]), np.abs(GMAX[m])) for m in range(7)}

# the jumps J_i = P^{(i)}(delta) - Q^{(i)}(delta)
JUMP = {i: float(Pd(i, np.array([DELTA]))[0] - Qd(i, np.array([DELTA]))[0]) for i in range(1, 6)}

# tail bound outside [TLO, THI]:  sup |{Wtil^{(m)}}| there times the Gaussian mass beyond 16.7 s
TAIL_SUP = {m: float(np.max(np.abs(Pd(m, np.linspace(DELTA, math.pi/2 - DM, 20001)))))
            for m in range(7)}


def _mass(phi, tl, tr):
    """int_{tl}^{tr} chi_sigma(phi - t) dt, vectorised over phi (column) and bands (row)"""
    return ndtr((phi[:, None] - tl[None, :])/SIGMA) - ndtr((phi[:, None] - tr[None, :])/SIGMA)


def W_value_enclosure(phi, kmax=4):
    """certified enclosure of W_sigma^{(k)}(phi) at the POINTS phi, k = 0..kmax"""
    phi = np.atleast_1d(np.asarray(phi, dtype=float))
    mass = _mass(phi, TL, TR)                          # (nphi, NB), >= 0
    out = []
    for k in range(kmax+1):
        lo = (np.where(GMIN[k][None, :] >= 0, GMIN[k][None, :], GMIN[k][None, :])*mass).sum(axis=1)
        hi = (GMAX[k][None, :]*mass).sum(axis=1)
        lo = (GMIN[k][None, :]*mass).sum(axis=1)
        # the jump terms
        for i in range(1, k):
            v = JUMP[i]*chi(phi - DELTA, k-1-i)
            lo = lo + v
            hi = hi + v
        tail = TAIL_SUP[k]*1e-50
        out.append((lo - tail - 1e-12*np.abs(lo) - 1e-14,
                    hi + tail + 1e-12*np.abs(hi) + 1e-14))
    return out


def W_absbound(phi_lo, phi_hi, mmax=5):
    """certified upper bound on |W_sigma^{(m)}(phi)| for phi in [phi_lo, phi_hi], m = 0..mmax"""
    phi_lo = np.atleast_1d(np.asarray(phi_lo, dtype=float))
    phi_hi = np.atleast_1d(np.asarray(phi_hi, dtype=float))
    # per band, the largest Gaussian mass over phi in the interval: the mass is unimodal in phi
    # with its maximum where the band is centred on phi
    phistar = np.clip(TMID[None, :], phi_lo[:, None], phi_hi[:, None])
    mass = (ndtr((phistar - TL[None, :])/SIGMA) - ndtr((phistar - TR[None, :])/SIGMA))
    out = []
    for m in range(mmax+1):
        b = (GABS[m][None, :]*mass).sum(axis=1)
        for i in range(1, m):
            b = b + abs(JUMP[i])*chi_absmax_on(phi_lo - DELTA, phi_hi - DELTA, m-1-i)
        out.append(b*(1.0 + 1e-12) + 1e-12 + TAIL_SUP[m]*1e-50)
    return out


def W_interval_enclosure(phi_lo, phi_hi, kmax=4):
    """enclosure of W_sigma^{(k)} over the whole interval [phi_lo, phi_hi], k = 0..kmax"""
    phi_lo = np.atleast_1d(np.asarray(phi_lo, dtype=float))
    phi_hi = np.atleast_1d(np.asarray(phi_hi, dtype=float))
    mid = 0.5*(phi_lo + phi_hi)
    om = 0.5*(phi_hi - phi_lo)
    val = W_value_enclosure(mid, kmax)
    bnd = W_absbound(phi_lo, phi_hi, kmax+1)
    out = []
    for k in range(kmax+1):
        lo, hi = val[k]
        out.append((lo - om*bnd[k+1], hi + om*bnd[k+1]))
    return out


if __name__ == "__main__":
    import json
    sys.path.insert(0, os.path.join(HERE, "rerun"))
    import fx_profile as FX                                  # the seat's own instrument (control)
    CONV = FX.Conv(SIGMA)

    OUT = {"what": "certified enclosures of W_sigma^{(k)} on the Theorem V.4' ball's phi range",
           "sigma": SIGMA, "band_width": WB, "n_bands": int(NB),
           "t_window": [TLO, THI], "jumps_J_i": JUMP,
           "Q_Cauchy_bounds": QCAUCHY[:7],
           "sign_alternation_certificate": sign_alternation_certificate()}

    ANG = np.array([DELTA + 4*SIGMA, DELTA + 4.5*SIGMA, DELTA + 5*SIGMA, DELTA + 6*SIGMA,
                    0.25, 0.30, 0.4, PHI0, 0.6, 0.75, 0.83627])
    enc = W_value_enclosure(ANG, 4)
    gl = CONV.derivs(ANG, kmax=4)
    rows = []
    ok = True
    for i, p in enumerate(ANG):
        r = {"phi": float(p), "phi_minus_delta_over_sigma": float((p - DELTA)/SIGMA)}
        for k in range(5):
            lo, hi = enc[k][0][i], enc[k][1][i]
            v = float(gl[k][i])
            c = bool(lo <= v <= hi)
            ok = ok and c
            r["k=%d" % k] = {"lo": float(lo), "hi": float(hi), "GL_value": v, "contains": c,
                             "rel_width": float((hi-lo)/max(abs(v), 1e-300))}
        rows.append(r)
    OUT["control_containment_of_the_GL_instrument"] = {"all_contained": ok, "rows": rows}

    bnd = W_absbound(np.array([DELTA + 4*SIGMA]), np.array([0.83627]), 5)
    OUT["absolute_bounds_over_the_whole_4sigma_range"] = {
        "m=%d" % m: float(bnd[m][0]) for m in range(6)}
    bl = W_absbound(np.array([DELTA + 4*SIGMA]), np.array([DELTA + 4.01*SIGMA]), 5)
    OUT["absolute_bounds_at_the_inner_edge"] = {"m=%d" % m: float(bl[m][0]) for m in range(6)}

    json.dump(OUT, open(os.path.join(HERE, "r1_profile_results.json"), "w"), indent=1,
              default=str)
    print("bands %d, window [%g, %g], band width %.4e" % (NB, TLO, THI, WB))
    print("jumps:", {k: "%.6f" % v for k, v in JUMP.items()})
    print("sign alternation certificate:",
          [(r["m"], r["all_coefficients_nonnegative"]) for r in OUT["sign_alternation_certificate"]])
    print("GL values contained in the certified enclosures at every test angle:", ok)
    for i, p in enumerate(ANG[:4]):
        print("  phi = delta + %.2f sigma :" % ((p - DELTA)/SIGMA),
              " ".join("k%d rel width %.2e" % (k, rows[i]["k=%d" % k]["rel_width"])
                       for k in range(5)))
    print("bounds over the whole 4 sigma range:", OUT["absolute_bounds_over_the_whole_4sigma_range"])
    print("wrote r1_profile_results.json")

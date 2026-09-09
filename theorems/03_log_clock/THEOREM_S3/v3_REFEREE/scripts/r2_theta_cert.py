"""
r2_theta_cert -- a CERTIFIED UPPER BOUND on vartheta, the (H3'_vartheta) constant of
THEOREM_S3_v3 sec.1.4, on the Theorem V.4' ball, at offsets 4 sigma, 5 sigma and 6 sigma.

WHAT IS BEING BOUNDED.  With M = rho_0 = 1, eta_P = -M/r,

    e := eta_0 - eta_P = 1/r - Theta(u) W_sigma(phi)/(rho N_sigma) ,   u = log rho ,

    vartheta_k := sup_{x in B(x_*,d)} sup_{|v| = 1} r(x)^{k+1} |d_v^k e(x)| / k!  ,  k = 0..4 ,
    vartheta   := max_k vartheta_k .

fix5/x1 reports a five-parameter GRID MAXIMUM of the same quantity and THEOREM_S3_v3 states
vartheta = 2 x that number.  A grid maximum is a LOWER bound on a supremum, and twice a lower
bound is not an upper bound.  This file computes an upper bound.

THREE REDUCTIONS, THEN INTERVAL ARITHMETIC.

(1)  e depends on x only through rho = |x| and the polar angle phi; along the line x + s v it
     depends on the direction only through a = <v, Yhat> and c = <v, e_z>, because
     rho(s)^2 = rho^2 + 2 s (r a + z c) + s^2 (this uses |v| = 1) and z(s) = z + s c.

(2)  RESCALE s = rho tau.  Then rho(s) = rho Pp(tau), r(s) = rho Rr(tau), z(s) = rho Z(tau) with

         Pp^2 = 1 + 2 p tau + tau^2 ,  p = a sin phi + c cos phi ,
         Z    = cos phi + c tau ,      Rr^2 = sin^2 phi + 2 a sin phi tau + (1 - c^2) tau^2 ,

     all INDEPENDENT of rho, and

         e(s) = (1/rho) E(tau) ,   E = 1/Rr - Theta(log rho + log Pp) W_sigma(Phi)/(Pp N) ,
         vartheta_k = sup sin^{k+1}(phi) | [tau^k] E | .

     rho survives only inside Theta, whose deviation from 1 on the ball is below 5e-5 at every
     admissible L; so the whole rho extent of the ball is carried in ONE interval and the search
     is one-dimensional, in phi.  (Enlarging the point set is always safe: it is a superset.)

(3)  The tau-coefficient of order k is a polynomial of degree <= k in (a, c), carried EXACTLY as
     a polynomial with interval coefficients, so the direction is never gridded in the main pass.
     Its maximum over the closed unit disc is bounded first by
     sum |coefficient| x max_{a^2+c^2<=1} |a|^i |c|^j , and then sharpened, on the boxes that
     bind, by interval evaluation over a grid of (a, c) cells.

Every elementary function is evaluated in outward-rounded interval arithmetic (r0_iv.py); the
mollified profile's derivatives come from the certified enclosures of r1_profile_cert.py.  The
result is an upper bound valid on the whole ball, so no SAFETY factor is needed or used.

Outputs -> r2_results.json
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r0_iv as R                                                          # noqa: E402
import r1_profile_cert as PC                                              # noqa: E402

DEG = math.pi/180.0
DELTA, DM, PHI0 = 7.5*DEG, 5.0*DEG, 30.0*DEG
SIGMA = 0.02
EPS_R = 0.25
# the certified N_sigma enclosure of fix5/x6, reproduced by re-running x6 here
N_LO, N_HI = 1.0124508487, 1.0124511502


def theta_derivs_iv(U, L):
    """interval enclosures of Theta^{(k)}(u), k = 0..4, for u in the interval U"""
    shape = U.lo.shape
    tot = [R.const(0.0, shape) for _ in range(5)]
    for lo_off, sgn in ((-EPS_R, 1.0), (-(L - EPS_R), -1.0)):
        arg = R.smul(1.0/EPS_R, R.IV(U.lo + lo_off, U.hi + lo_off))
        T = R.tanh(arg)
        one = R.const(1.0, shape)
        T2 = R.mul(T, T)
        om = R.sub(one, T2)
        d = [T, om, R.smul(-2.0, R.mul(T, om)),
             R.smul(-2.0, R.mul(om, R.sub(one, R.smul(3.0, T2)))),
             R.smul(8.0, R.mul(T, R.mul(om, R.sub(R.const(2.0, shape), R.smul(3.0, T2)))))]
        for k in range(5):
            tot[k] = R.add(tot[k], R.smul(0.5*sgn/EPS_R**k, d[k]))
    return tot


def series_E(phi1, phi2, u_lo, u_hi, L=40.0):
    """the tau-series of E on the box phi in [phi1,phi2], log rho in [u_lo,u_hi]"""
    shape = phi1.shape
    PHI = R.IV(phi1, phi2)
    sp_ = R.sin_mono(PHI)
    cp_ = R.cos_mono(PHI)

    P2 = R.szero()
    P2[0] = R.pconst(R.const(1.0, shape))
    P2[1][R.IDX[(1, 0)]] = R.smul(2.0, sp_)
    P2[1][R.IDX[(0, 1)]] = R.smul(2.0, cp_)
    P2[2] = R.pconst(R.const(1.0, shape))

    Z = R.szero()
    Z[0] = R.pconst(cp_)
    Z[1][R.IDX[(0, 1)]] = R.const(1.0, shape)

    R2 = R.szero()
    R2[0] = R.pconst(R.ipow(sp_, 2))
    R2[1][R.IDX[(1, 0)]] = R.smul(2.0, sp_)
    R2[2][R.IDX[(0, 0)]] = R.const(1.0, shape)
    R2[2][R.IDX[(0, 2)]] = R.const(-1.0, shape)

    Rr = R.s_sqrt(R2)
    Pp = R.s_sqrt(P2)
    lP = R.s_log(Pp)
    U = [list(p) for p in lP]
    U[0] = R.pconst(R.add(R.IV(u_lo, u_hi), lP[0][0]))
    Phi = R.s_arctan(R.smul_series(Rr, R.s_inv(Z)))

    Td = theta_derivs_iv(U[0][0], L)
    PHb = Phi[0][0]
    enc = PC.W_interval_enclosure(PHb.lo, PHb.hi, 4)
    Wd = [R.IV(enc[k][0], enc[k][1]) for k in range(5)]

    Th = R.scompose(Td, U)
    Wc = R.scompose(Wd, Phi)
    invN = R.IV(np.full(shape, 1.0/N_HI), np.full(shape, 1.0/N_LO))
    term = R.smul_series(Th, R.smul_series(Wc, R.s_inv(Pp)))
    ri = R.s_inv(Rr)
    E = [R.psub(ri[k], R.pmulIV(invN, term[k])) for k in range(5)]
    return E, sp_


def bound_mono(E, sp_):
    return np.array([np.nextafter(sp_.hi**(k+1)*R.pmag_disc(E[k]), np.inf) for k in range(5)])


_ACCELLS = {}


def ac_cells(n):
    if n in _ACCELLS:
        return _ACCELLS[n]
    e = np.linspace(-1.0, 1.0, n+1)
    A1, C1 = np.meshgrid(e[:-1], e[:-1], indexing="ij")
    A2, C2 = np.meshgrid(e[1:], e[1:], indexing="ij")
    a1, a2, c1, c2 = A1.ravel(), A2.ravel(), C1.ravel(), C2.ravel()
    near = (np.minimum(a1**2, a2**2)*(np.sign(a1) == np.sign(a2))
            + np.minimum(c1**2, c2**2)*(np.sign(c1) == np.sign(c2)))
    keep = near <= 1.0
    _ACCELLS[n] = (a1[keep], a2[keep], c1[keep], c2[keep])
    return _ACCELLS[n]


def bound_grid(E, sp_, n=48):
    """sharpen the direction bound by interval evaluation over (a, c) cells"""
    a1, a2, c1, c2 = ac_cells(n)
    shape = sp_.lo.shape
    out = np.zeros((5,) + shape)
    one = R.IV(np.ones(1), np.ones(1))
    for i in range(a1.size):
        A = R.IV(np.array([a1[i]]), np.array([a2[i]]))
        C = R.IV(np.array([c1[i]]), np.array([c2[i]]))
        Ap = [one]
        Cp = [one]
        for _ in range(4):
            Ap.append(R.mul(Ap[-1], A))
            Cp.append(R.mul(Cp[-1], C))
        for k in range(5):
            tot = None
            for idx in range(R.NMON):
                co = E[k][idx]
                if co is None:
                    continue
                ii, jj = R.MON[idx]
                t = co
                if ii:
                    t = R.mul(t, Ap[ii])
                if jj:
                    t = R.mul(t, Cp[jj])
                tot = t if tot is None else R.add(tot, t)
            if tot is not None:
                out[k] = np.maximum(out[k], tot.mag())
    return np.array([np.nextafter(sp_.hi**(k+1)*out[k], np.inf) for k in range(5)])


def geo_edges(lo, hi, n, beta):
    """a two-sided geometrically graded partition of [lo, hi]: fine at both edges, where the
    smoothing layer (inner) and the largest sin phi (outer) live."""
    t = np.linspace(0.0, 1.0, n + 1)
    g = (np.exp(beta*t) - 1.0)/(math.exp(beta) - 1.0)
    mid = 0.5*(lo + hi)
    left = lo + (mid - lo)*g
    right = hi - (hi - mid)*g[::-1] - 0.0
    right = hi - (hi - mid)*g
    e = np.unique(np.concatenate([left, right[::-1], [lo, hi]]))
    return e


def _eval_boxes(p1, p2, u_lo, u_hi, L, nac, chunk=4000):
    vals = np.empty((5, p1.size))
    for i0 in range(0, p1.size, chunk):
        sl = slice(i0, i0 + chunk)
        E, sp_ = series_E(p1[sl], p2[sl], np.full(p2[sl].shape, u_lo),
                          np.full(p2[sl].shape, u_hi), L)
        m = bound_mono(E, sp_)
        g = bound_grid(E, sp_, nac)
        vals[:, sl] = np.minimum(m, g)
    return vals


def certify(f, offset, n_side=2000, beta=7.0, nac=64, L=40.0, refine_rounds=9, ntop=1200,
            verbose=True, u_override=None):
    rho_s = 1.0 + f
    d = min(rho_s*math.sin(PHI0 - DELTA - offset),
            rho_s*math.sin(math.pi/2 - DM - PHI0 - offset), f)
    psi_max = math.asin(min(d/rho_s, 1.0))
    phi_min, phi_max = PHI0 - psi_max, PHI0 + psi_max
    u_lo = math.log(rho_s - d)
    u_hi = math.log(rho_s + d)
    if u_override is not None:
        u_lo, u_hi = u_override

    edges = geo_edges(phi_min, phi_max, n_side, beta)
    p1, p2 = edges[:-1].copy(), edges[1:].copy()
    vals = _eval_boxes(p1, p2, u_lo, u_hi, L, nac)
    hist = []
    for it in range(refine_rounds + 1):
        per_box = vals.max(axis=0)
        best = float(per_box.max())
        hist.append({"round": it, "n_boxes": int(p1.size), "bound": best,
                     "per_k": [float(v) for v in vals.max(axis=1)],
                     "argmax_phi_deg": float(0.5*(p1[int(np.argmax(per_box))]
                                                  + p2[int(np.argmax(per_box))])/DEG)})
        if verbose:
            print("    round %d  boxes %6d  bound %.6f  per k %s  argmax phi %.4f deg"
                  % (it, p1.size, best, " ".join("%.4f" % v for v in vals.max(axis=1)),
                     hist[-1]["argmax_phi_deg"]), flush=True)
        if it == refine_rounds:
            break
        idx = np.argsort(per_box)[-ntop:]
        pm = 0.5*(p1[idx] + p2[idx])
        np1 = np.concatenate([p1[idx], pm])
        np2 = np.concatenate([pm, p2[idx]])
        nv = _eval_boxes(np1, np2, u_lo, u_hi, L, nac)
        keep = np.setdiff1d(np.arange(p1.size), idx)
        p1 = np.concatenate([p1[keep], np1])
        p2 = np.concatenate([p2[keep], np2])
        vals = np.concatenate([vals[:, keep], nv], axis=1)

    return {"f": f, "offset_over_sigma": offset/SIGMA, "d": d,
            "phi_min_deg": phi_min/DEG, "phi_max_deg": phi_max/DEG,
            "rho_min": rho_s - d, "rho_max": rho_s + d,
            "r_star": rho_s*math.sin(PHI0), "R_minus": rho_s*math.sin(PHI0) - d,
            "certified_vartheta_upper": hist[-1]["bound"],
            "certified_per_k": hist[-1]["per_k"],
            "n_boxes_final": hist[-1]["n_boxes"], "n_direction_cells": int(nac*nac),
            "history": hist}


if __name__ == "__main__":
    t0 = time.time()
    OUT = {"what": "certified UPPER bound on vartheta for (H3'_vartheta) on the Theorem V.4' "
                   "ball; outward-rounded interval arithmetic over phi boxes, rho carried whole, "
                   "direction carried exactly as a polynomial in (a, c)",
           "N_sigma_enclosure": [N_LO, N_HI], "L_used": 40.0,
           "grid_lower_bounds_from_fix5_x1": {}}
    x1 = json.load(open(os.path.join(HERE, "rerun", "x1_results.json")))
    for k, v in x1["repaired"].items():
        OUT["grid_lower_bounds_from_fix5_x1"][k] = {"vartheta_grid": v["vartheta"],
                                                    "stated_2x": v["safety_vartheta"]}
    OUT["offset_scan_grid_lower_bounds"] = {k: v["vartheta"]
                                            for k, v in x1["offset_scan_f4"].items()}

    res = {}
    for off_m in [4.0, 5.0, 6.0]:
        print("  offset %.1f sigma, f = 4" % off_m, flush=True)
        res["f=4,offset=%gsigma" % off_m] = certify(4.0, off_m*SIGMA)
    for f in [8.0, 16.0, 32.0]:
        print("  offset 4.0 sigma, f = %g" % f, flush=True)
        res["f=%g,offset=4sigma" % f] = certify(f, 4.0*SIGMA)
    OUT["certified"] = res
    top4 = max(res["f=%g,offset=4sigma" % f]["certified_vartheta_upper"]
               for f in [8.0, 16.0, 32.0])
    top4 = max(top4, res["f=4,offset=4sigma"]["certified_vartheta_upper"])
    OUT["headline"] = {
        "vartheta_certified_upper_over_f_at_4sigma": top4,
        "vartheta_certified_upper_at_5sigma_f4": res["f=4,offset=5sigma"][
            "certified_vartheta_upper"],
        "vartheta_certified_upper_at_6sigma_f4": res["f=4,offset=6sigma"][
            "certified_vartheta_upper"],
        "vartheta_stated_in_THEOREM_S3_v3": 0.2623138022,
        "vartheta_grid_lower_bound_x1_max_over_f": max(
            v["vartheta_grid"] for v in OUT["grid_lower_bounds_from_fix5_x1"].values()),
        "seconds": time.time() - t0}
    h = OUT["headline"]
    h["certified_below_one"] = bool(top4 < 1.0)
    h["certified_below_the_stated_value"] = bool(top4 <= 0.2623138022)
    h["Q_at_certified"] = (1.0 + top4)/(1.0 - top4) if top4 < 1 else None
    h["vt_at_certified"] = top4/(1.0 - top4) if top4 < 1 else None
    json.dump(OUT, open(os.path.join(HERE, "r2_results.json"), "w"), indent=1, default=str)
    print(json.dumps(OUT["headline"], indent=1))
    print("wrote r2_results.json in %.1f s" % (time.time() - t0))

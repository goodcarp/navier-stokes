"""
g4 -- UNIT 3(b) of FIX3: the honest measured K_2 on the TRANSPORTED field at lambda = 3/2.

THE OPEN ITEM (FIX2 ADDENDUM_2, "Also open").  hk2's quoted measured slack over the window
compares a TRANSPORTED bound with an UNTRANSPORTED measurement: hk2's k4 zonal instrument
evaluates eta_0 at the moved point, which hk2 sec.14 item 7 itself declares wrong for the
bound, while k3's bound uses eta_lam = eta_0 o T_lam^{-1}.  fix2 could not repair it because
eta_lam is not of the separable form G(rho)H(t) the zonal series needs.

THE INSTRUMENT HERE.  A direct 5-D kernel quadrature that never needs separability, because it
puts every derivative on eta and none on the kernel:

    a(x)          = int K(w) eta(x-w) dw ,        K(w) = 3 w_z/(8 pi^2 |w|^5)
    d_i a(x)      = int K(w) (d_i eta)(x-w) dw
    d_i d_j a(x)  = int K(w) (d_i d_j eta)(x-w) dw

all absolutely convergent for (D-B) (real-analytic, |grad^2 eta| = O(rho^{-3}) at infinity and
bounded near x), since |K| ~ |w|^{-4} and dw = |w|^4 d|w| dOmega_4.  In the axisymmetric frame
at the source point x-w, with mu := yhat . yhat' = (r - s sinTheta cos chi)/r' ,

    a_r  = int K  eta_r  mu                       a_z  = int K  eta_z
    a_rr = int K [ eta_rr mu^2 + (eta_r/r')(1-mu^2) ]
    a_rz = int K  eta_rz mu                       a_zz = int K  eta_zz
    a_perp := <yhat_perp . Hess eta . yhat_perp>  = int K [ eta_rr q + (eta_r/r')(1-q) ] ,
        q = s^2 sin^2(Theta) sin^2(chi) / (3 r'^2)      (the S^2 average of (yhat_perp.yhat')^2)

CONTROLS, all independent of the object measured:
  (i)   a_perp must equal a_r/r  (axisymmetry of a);
  (ii)  Lap5 a = a_rr + 3 a_r/r + a_zz must equal d_z eta(x)  (the PDE);
  (iii) at lambda = 1 the whole vector must reproduce hk2's k4 zonal series.

Outputs -> g4_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "fix2copy"))
sys.path.insert(0, os.path.join(HERE, "fix2copy", "imported"))
import numpy as np
import k4_direct as K4                                            # noqa: E402 (hk2, unchanged)
from hk2lib import DatumB, StrainedB, geometry                    # noqa: E402 (hk2, unchanged)

OUT = {}
CK = 3.0/(8*math.pi**2)


def a_and_derivs(D, r, z, ns=900, nT=200, nc=200, smin=1e-6, smax=None):
    """(a, a_r, a_z, a_rr, a_rz, a_zz, a_perp) by direct 5-D kernel quadrature."""
    if smax is None:
        smax = 3*D.R
    xg, wg = np.polynomial.legendre.leggauss(ns)
    lo, hi = math.log(smin), math.log(smax)
    u = 0.5*(hi-lo)*xg + 0.5*(hi+lo)
    wu = 0.5*(hi-lo)*wg
    s = np.exp(u)
    T, wT = np.polynomial.legendre.leggauss(nT)
    T = 0.5*math.pi*(T+1)
    wT = 0.5*math.pi*wT
    X, wX = np.polynomial.legendre.leggauss(nc)
    X = 0.5*math.pi*(X+1)
    wX = 0.5*math.pi*wX
    sT, cT = np.sin(T), np.cos(T)
    sX, cX = np.sin(X), np.cos(X)
    WA = (wT*sT**3)[:, None]*(wX*4*math.pi*sX**2)[None, :]
    KA = WA*cT[:, None]*CK                        # K(w) dOmega_4 factor (the s^4/|w|^5 gone)
    acc = np.zeros(7)
    ST = sT[:, None]
    CT = cT[:, None]
    CXr = cX[None, :]
    SXr = sX[None, :]
    for si, wsi in zip(s, wu*s):
        rp2 = r*r + (si*ST)**2 - 2*r*si*ST*CXr
        rp = np.sqrt(np.maximum(rp2, 1e-300))
        zp = z - si*CT*np.ones_like(CXr)
        e, er, ez, err, erz, ezz = D.eta_rz(rp, zp)
        mu = (r - si*ST*CXr)/rp
        q = (si*ST*SXr)**2/(3.0*rp2)
        er_over = er/rp
        acc[0] += wsi*np.sum(KA*e)
        acc[1] += wsi*np.sum(KA*er*mu)
        acc[2] += wsi*np.sum(KA*ez)
        acc[3] += wsi*np.sum(KA*(err*mu*mu + er_over*(1.0-mu*mu)))
        acc[4] += wsi*np.sum(KA*erz*mu)
        acc[5] += wsi*np.sum(KA*ezz)
        acc[6] += wsi*np.sum(KA*(err*q + er_over*(1.0-q)))
    return acc


def measure(D, Dfield, r, z, ns=900, nT=200, nc=200):
    A = a_and_derivs(Dfield, r, z, ns, nT, nc)
    a, a_r, a_z, a_rr, a_rz, a_zz, a_perp = A
    ee = [float(v[0]) for v in Dfield.eta_rz(np.array([r]), np.array([z]))]
    lap = a_rr + 3.0*a_r/r + a_zz
    k2m, arg = K4.K2_exact(r, a_r, a_z, a_rr, a_rz, a_zz, ee[0], ee[1], ee[2], ne=121)
    ha = np.linalg.eigvalsh(np.array([[a_rr, a_rz], [a_rz, a_zz]]))
    return {"a": a, "a_r": a_r, "a_z": a_z, "a_rr": a_rr, "a_rz": a_rz, "a_zz": a_zz,
            "a_perp": a_perp, "a_r_over_r": a_r/r,
            "control_axisym_rel": abs(a_perp - a_r/r)/max(abs(a_r/r), 1e-300),
            "grad_a": math.hypot(a_r, a_z),
            "hess_a": float(max(abs(ha).max(), abs(a_r/r))),
            "Lap5_a": lap, "dz_eta": ee[2],
            "control_PDE_abs": abs(lap - ee[2]),
            "control_PDE_over_hess_scale": abs(lap - ee[2])/max(abs(a_rr), abs(a_zz),
                                                                abs(a_r/r), 1e-300),
            "K2_measured": k2m, "eta": ee[0], "eta_r": ee[1], "eta_z": ee[2]}


if __name__ == "__main__":
    DB = DatumB(delta_deg=7.5, w=0.20, L=10.0)
    # hk2's own window and tracked point: phi0 = 30 deg, f = 0
    rows = {}
    for lam in (1.0, 1.25, 1.5):
        g = geometry(30.0, 0.0, 7.5, lam)
        DL = StrainedB(DB, lam) if lam != 1.0 else DB
        per = {}
        for (ns, nT, nc) in ((600, 160, 160), (900, 200, 200), (1400, 260, 260)):
            m = measure(DB, DL, g["r"], g["z"], ns, nT, nc)
            per["ns=%d" % ns] = m
            print("lam=%4.2f (%4d,%3d,%3d): a=%12.6f |grad a|=%10.6f ||Hess a||=%10.6f "
                  "K2=%9.5f | axisym %.2e  PDE/scale %.2e"
                  % (lam, ns, nT, nc, m["a"], m["grad_a"], m["hess_a"], m["K2_measured"],
                     m["control_axisym_rel"], m["control_PDE_over_hess_scale"]), flush=True)
        rows["lam=%g" % lam] = per
    OUT["transported_measured_DB"] = rows

    # hk2's own untransported k4 numbers, for the comparison the open item is about
    hk2_k4 = json.load(open(os.path.join(HERE, "..", "..", "..", "hk2", "k4_results.json")))
    unt = {("lam=%g" % r["lam"]): {"K2": r["K2"], "a": r["a"], "grad_a": r["grad_a"],
                                   "hess_a": r["hess_a"]} for r in hk2_k4["measured"]}
    OUT["hk2_k4_untransported"] = unt
    bound = {"lam=1": 66.6622, "lam=1.25": 107.6738, "lam=1.5": 161.7735}
    slack = {}
    for k in ("lam=1", "lam=1.25", "lam=1.5"):
        best = rows[k]["ns=1400"]["K2_measured"]
        slack[k] = {"hk2_bound": bound[k],
                    "transported_measured": best,
                    "slack_transported": bound[k]/best,
                    "untransported_measured": unt[k]["K2"],
                    "slack_as_quoted": bound[k]/unt[k]["K2"],
                    "transported_over_untransported": best/unt[k]["K2"]}
        print("%s: bound %.4f ; transported %.5f (slack %.2f x) ; untransported %.5f "
              "(slack %.2f x)" % (k, bound[k], best, bound[k]/best, unt[k]["K2"],
                                  bound[k]/unt[k]["K2"]), flush=True)
    OUT["slack"] = slack
    OUT["record_K2_measured_window"] = {"t3_budget_RECORD": 5.3854,
                                        "u1_sec6.6_label": "hk2 measured, over the window"}
    with open(os.path.join(HERE, "g4_results.json"), "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
    print("\nWROTE g4_results.json")

"""
f2b -- the MEASURED K_2 for THEOREM_S3's datum at lambda = 1, and the slack against f2's
       proved bound.

Two changes to hk2's k4 instrument, both forced by THEOREM_S3's datum and both declared:

  1. hk2's k4.hhat() integrates the Gegenbauer coefficients with ONE 4001-point Gauss rule on
     [-1,1].  That is right for (D-B), which is real-analytic.  THEOREM_S3's angular profile
     has kinks at t = cos delta, sin delta_m, -sin delta_m, -cos delta, so a single Gauss rule
     mis-resolves the high-l coefficients, which is exactly what the second derivatives of a
     depend on.  Here the coefficients are computed with the kinks on PANEL BOUNDARIES.

  2. k4's PDE control divides by |d_z eta|.  At the tracked point of ADDENDUM_1's f >= 4 the
     datum is in the deep plateau, where d_z eta = O(10^-5) while the Hessian of a is O(10^-2):
     that normalisation reports a 4e-4 absolute residual as a "relative error" of 39.  The
     control is renormalised here by the scale of the object being measured.

  lambda > 1 is NOT measurable with this instrument: eta_lam = eta_0 o T_lam^{-1} is not of the
  separable form G(rho) H(t) that the zonal series needs.  hk2's own measured column has the
  same limitation -- k4 evaluates the UNTRANSPORTED eta_0 at the moved point, which hk2 sec.14
  item 7 declares wrong for the bound.  So the slack below is stated at lambda = 1 only.

Outputs -> f2b_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "imported"))
sys.path.insert(0, HERE)
import numpy as np                                    # noqa: E402
from scipy.special import eval_gegenbauer             # noqa: E402
import k4_direct as K4                                # noqa: E402
from f2_datum_s3 import DatumS3, geometry_s3          # noqa: E402

DEG = math.pi/180.0
DELTA, DM = 7.5*DEG, 5.0*DEG


class DatumS3Zonal(DatumS3):
    def Hh(self, t, k=0):
        return self.W(np.arccos(np.clip(t, -1.0, 1.0)), k)


def hhat_panels(D, LMAX, npan=800):
    """Gegenbauer coefficients of Hh(t) with the kinks of the profile on panel boundaries."""
    kinks = sorted({-1.0, 1.0, math.cos(D.delta), -math.cos(D.delta),
                    math.sin(D.dm), -math.sin(D.dm)})
    x0, w0 = np.polynomial.legendre.leggauss(npan)
    X, W = [], []
    for a, b in zip(kinks[:-1], kinks[1:]):
        X.append(0.5*(b-a)*x0 + 0.5*(a+b))
        W.append(0.5*(b-a)*w0)
    X = np.concatenate(X)
    W = np.concatenate(W)
    H = D.Hh(X, 0)
    ww = W*(1 - X*X)
    out = np.zeros(LMAX+1)
    for l in range(LMAX+1):
        Nl = (l+1)*(l+2)/(l+1.5)
        out[l] = np.sum(ww*H*eval_gegenbauer(l, 1.5, X))/Nl
    return out


def a_kernel(D, x_rz, ns=900, nT=200, nc=200):
    """a(x) = (K * eta)(x) by direct 5-D quadrature (k4 control 2, datum swapped)."""
    r, z = x_rz
    xg, wg = np.polynomial.legendre.leggauss(ns)
    smin, smax = 1e-6, 3*D.R
    u = 0.5*(math.log(smax)-math.log(smin))*xg + 0.5*(math.log(smax)+math.log(smin))
    wu = 0.5*(math.log(smax)-math.log(smin))*wg
    s = np.exp(u)
    T, wT = np.polynomial.legendre.leggauss(nT)
    T = 0.5*math.pi*(T+1); wT = 0.5*math.pi*wT
    X, wX = np.polynomial.legendre.leggauss(nc)
    X = 0.5*math.pi*(X+1); wX = 0.5*math.pi*wX
    sT, cT, cX = np.sin(T), np.cos(T), np.cos(X)
    WA = (wT*sT**3)[:, None]*(wX*4*math.pi*np.sin(X)**2)[None, :]
    tot = 0.0
    for si, wsi in zip(s, wu*s):
        rp = np.sqrt(np.maximum(r*r + (si*sT[:, None])**2 - 2*r*si*sT[:, None]*cX[None, :], 1e-300))
        zp = z - si*cT[:, None]*np.ones_like(cX)[None, :]
        et = D.eta_rz(rp, zp)[0]
        tot += wsi*(3.0/(8*math.pi**2))*np.sum(WA*cT[:, None]*et)
    return float(tot)


if __name__ == "__main__":
    OUT = {}
    D = DatumS3Zonal(delta_deg=7.5, dm_deg=5.0, eps_r=0.25, L=10.0)
    rows = {}
    for f in (1.0, 4.0):
        g = geometry_s3(30.0, f, 7.5, 5.0, 1.0)
        r_, z_ = g["r"], g["z"]
        ee = [float(v[0]) for v in D.eta_rz(np.array([r_]), np.array([z_]))]
        ak = a_kernel(D, (r_, z_))
        per_f = {}
        for LMAX in (81, 161, 321, 641):
            S = K4.Series(D, LMAX=LMAX)
            S.hh = hhat_panels(D, LMAX+2)
            S.cacheP = {}
            A = S.a_field(r_, z_)
            lap = A[3] + 3*A[1]/r_ + A[5]
            scale = max(abs(A[3]), abs(A[5]), abs(A[1]/r_), 1e-300)
            k2m, arg = K4.K2_exact(r_, A[1], A[2], A[3], A[4], A[5], ee[0], ee[1], ee[2], ne=161)
            ha = np.linalg.eigvalsh(np.array([[A[3], A[4]], [A[4], A[5]]]))
            per_f["LMAX=%d" % LMAX] = {
                "a": A[0], "a_kernel_quadrature": ak, "a_rel_vs_kernel": abs(A[0]-ak)/abs(ak),
                "grad_a": math.hypot(A[1], A[2]),
                "hess_a_op": float(max(abs(ha).max(), abs(A[1]/r_))),
                "K2_measured": k2m,
                "PDE_residual_abs": abs(lap - ee[2]),
                "PDE_residual_over_hessian_scale": abs(lap - ee[2])/scale,
                "dz_eta": ee[2]}
            print("f=%g LMAX=%4d: K2=%.6f |grad a|=%.6f ||Hess a||=%.6f  "
                  "PDE resid/Hess = %.2e   a rel vs kernel = %.2e"
                  % (f, LMAX, k2m, math.hypot(A[1], A[2]),
                     per_f["LMAX=%d" % LMAX]["hess_a_op"],
                     per_f["LMAX=%d" % LMAX]["PDE_residual_over_hessian_scale"],
                     per_f["LMAX=%d" % LMAX]["a_rel_vs_kernel"]), flush=True)
        rows["f=%g" % f] = per_f
    OUT["measured_lam1"] = rows
    F2 = json.load(open(os.path.join(HERE, "f2_results.json")))
    tab = F2["C_K2_table"]
    OUT["slack_lam1"] = {}
    for f in (1.0, 4.0):
        b = tab["L=10_f=%g_lam=1" % f]
        m = rows["f=%g" % f]["LMAX=641"]["K2_measured"]
        OUT["slack_lam1"]["f=%g" % f] = {"bound": b, "measured": m, "slack": b/m}
    with open(os.path.join(HERE, "f2b_results.json"), "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
    print(json.dumps(OUT, indent=1, sort_keys=True, default=str))

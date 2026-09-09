"""
fx_profile -- the sigma-mollified angular profile, evaluated by a fixed-node Gauss-Legendre
convolution, with rigorous a priori bounds on every derivative used downstream.

WHY A THIRD INSTRUMENT.  fix4/p1_profile.py evaluates W_sigma by an FFT on a 600001-point
uniform grid and then takes GRID maxima; the referee's scripts/s1_profile_indep.py evaluates it
by scipy adaptive quadrature at one phi at a time.  Neither gives W_sigma''' or W_sigma''''
(needed for the (H3'_vartheta) search at k = 3, 4), and the FFT route gives no certified pad.
This file evaluates

    W_sigma^{(k)}(phi) = int Wtil(t) chi_sigma^{(k)}(phi - t) dt ,     k = 0..4,

with the integration variable t, so the breakpoints of Wtil sit at FIXED abscissae independent
of phi and one panel structure serves every phi.  Panels are the intervals between the kink
images, subdivided to width <= sigma/2, with 24-point Gauss-Legendre on each: the integrand is
real-analytic on every panel, so the quadrature is spectrally accurate, and the same node set is
reused for all five k and all phi (a matrix product).

CONTROLS (L-14).  Reproduce (i) the referee's adaptive-quadrature Ws, Ws1, Ws2, Ws3, and
(ii) fix4's FFT W0, W1, W2 on its own grid, and (iii) the closed-form identity
W_sigma''''(delta) = [W'(delta+) - W'(delta-)] chi_sigma''(0).

RIGOROUS BOUNDS.  ||W_sigma^{(k)}||_inf <= ||Wtil'||_inf ||chi_sigma^{(k-1)}||_1 for k >= 1, and
||W_sigma||_inf <= ||Wtil||_inf, with the L^1 norms of the Gaussian derivatives in closed form
(each is the total variation of the previous derivative, evaluated at its critical points).
These are the numbers the certified pads of x6 use.
"""
import math
import numpy as np

DEG = math.pi/180.0
DELTA = 7.5*DEG
DM = 5.0*DEG
PHI0 = 30.0*DEG
EPS_R = 0.25
SIGMA_DEFAULT = 0.02

# ---------------------------------------------------------------- the kinked profile, exact
def W_exact_scalar(phi):
    """W = A/sin phi on (0,pi), A = sgn(cos)min(1,phi_ax/delta)min(1,|phi-pi/2|/dm)."""
    if phi <= 0.0:
        return 1.0/DELTA
    if phi >= math.pi:
        return -1.0/DELTA
    if phi > math.pi/2:
        return -W_exact_scalar(math.pi - phi)
    s = math.sin(phi)
    if phi < DELTA:
        return (phi/(DELTA*s)) if phi > 1e-8 else (1.0 + phi*phi/6.0)/DELTA
    if phi > math.pi/2 - DM:
        return (math.pi/2 - phi)/(DM*s)
    return 1.0/s

def W_exact(phi):
    phi = np.asarray(phi, dtype=float)
    out = np.empty_like(phi)
    flat = phi.ravel()
    o = out.ravel()
    for i, p in enumerate(flat):
        o[i] = W_exact_scalar(float(p))
    return out

def Wp_exact_scalar(phi):
    """W' on (0,pi), a.e.  W(pi-p) = -W(p) => W'(pi-p) = +W'(p)."""
    if phi > math.pi/2:
        return Wp_exact_scalar(math.pi - phi)
    if phi <= 0.0:
        return 0.0
    s, c = math.sin(phi), math.cos(phi)
    if phi < DELTA:
        num = (phi**3/3.0*(1.0 - phi*phi/10.0)) if phi < 1e-3 else (s - phi*c)
        return num/(DELTA*s*s)
    if phi > math.pi/2 - DM:
        t = math.pi/2 - phi
        return (-s - t*c)/(DM*s*s)
    return -c/(s*s)

def Wtil_scalar(y):
    """even reflection at 0 and at pi, 2 pi periodic."""
    y = math.fmod(y, 2.0*math.pi)
    if y < 0.0:
        y += 2.0*math.pi
    if y > math.pi:
        y = 2.0*math.pi - y
    return W_exact_scalar(y)

WMAX = 1.0/math.sin(DELTA)                      # sup |Wtil|            = 7.66129757554039
WPMAX = math.cos(DELTA)/math.sin(DELTA)**2      # sup |Wtil'|           = 58.19333257...

# ---------------------------------------------------------------- Gaussian and its L^1 norms
def chi_derivs(y, sg, kmax=4):
    """chi_sigma^{(k)}(y), k = 0..kmax, vectorised."""
    ch = np.exp(-0.5*(y/sg)**2)/(sg*math.sqrt(2.0*math.pi))
    s2 = sg*sg
    out = [ch]
    if kmax >= 1:
        out.append(-(y/s2)*ch)
    if kmax >= 2:
        out.append((y*y/s2 - 1.0)/s2*ch)
    if kmax >= 3:
        out.append((y/s2**2)*(3.0 - y*y/s2)*ch)
    if kmax >= 4:
        out.append((3.0/s2**2 - 6.0*y*y/s2**3 + y**4/s2**4)*ch)
    return out

def chi_L1(k, sg):
    """||chi_sigma^{(k)}||_1, closed form (total variation of chi^{(k-1)})."""
    rt = math.sqrt(2.0*math.pi)
    if k == 0:
        return 1.0
    if k == 1:                       # 2 chi(0)
        return 2.0/(sg*rt)
    if k == 2:                       # -4 chi'(sigma)
        return 4.0*math.exp(-0.5)/(sg*sg*rt)
    if k == 3:                       # 2|chi''(0)| + 4 chi''(sqrt3 sigma)
        return (2.0 + 8.0*math.exp(-1.5))/(sg**3*rt)
    if k == 4:                       # TV of chi''' : critical points of chi''' are the roots
        r = np.roots([1.0, 0.0, -6.0, 0.0, 3.0])      # y/sg roots of chi'''' = 0
        r = np.sort(np.real(r[np.abs(np.imag(r)) < 1e-12]))
        pts = [-np.inf] + [float(v)*sg for v in r] + [np.inf]
        tv = 0.0
        prev = 0.0
        for p in pts[1:-1]:
            v = float(chi_derivs(np.array([p]), sg, 3)[3][0])
            tv += abs(v - prev)
            prev = v
        tv += abs(0.0 - prev)
        return tv
    raise ValueError(k)

def W_bounds(sg, kmax=4):
    """rigorous sup bounds on |W_sigma^{(k)}| and on |A_sigma^{(k)}| = |(W_sigma sin)^{(k)}|."""
    BW = [WMAX] + [WPMAX*chi_L1(k-1, sg) for k in range(1, kmax+1)]
    BA = []
    for k in range(kmax+1):
        # (W sin)^{(k)} = sum_j C(k,j) W^{(j)} (sin)^{(k-j)} , |(sin)^{(m)}| <= 1
        BA.append(sum(math.comb(k, j)*BW[j] for j in range(k+1)))
    return {"BW": BW, "BA": BA,
            "chi_L1": [chi_L1(k, sg) for k in range(kmax+1)]}

# ---------------------------------------------------------------- the panel evaluator
class Conv:
    """W_sigma^{(k)}(phi) = sum_j w_j Wtil(t_j) chi^{(k)}(phi - t_j) on a fixed panel set."""
    def __init__(self, sigma=SIGMA_DEFAULT, lo=None, hi=None, gl=24, trunc=16.0):
        self.sigma = sigma
        T = trunc*sigma
        self.lo = (-T) if lo is None else lo
        self.hi = (math.pi + T) if hi is None else hi
        brk = sorted(set([b for b in
                          [-math.pi/2 - DM, -math.pi/2 + DM, -DELTA, DELTA,
                           math.pi/2 - DM, math.pi/2 + DM, math.pi - DELTA, math.pi + DELTA,
                           3*math.pi/2 - DM, 3*math.pi/2 + DM]
                          if self.lo < b < self.hi]))
        edges = [self.lo] + brk + [self.hi]
        wmax = 0.5*sigma
        xs, ws = np.polynomial.legendre.leggauss(gl)
        T_, W_ = [], []
        for a, b in zip(edges[:-1], edges[1:]):
            n = max(1, int(math.ceil((b - a)/wmax)))
            e = np.linspace(a, b, n+1)
            for u, v in zip(e[:-1], e[1:]):
                T_.append(0.5*(v-u)*xs + 0.5*(u+v))
                W_.append(0.5*(v-u)*ws)
        self.t = np.concatenate(T_)
        self.w = np.concatenate(W_)
        self.f = np.array([Wtil_scalar(float(x)) for x in self.t])*self.w
        self.trunc = T

    def derivs(self, phi, kmax=4, chunk=400):
        """returns list of arrays W_sigma^{(k)}(phi), k = 0..kmax."""
        phi = np.atleast_1d(np.asarray(phi, dtype=float))
        out = [np.zeros(phi.shape) for _ in range(kmax+1)]
        for i0 in range(0, phi.size, chunk):
            ph = phi[i0:i0+chunk]
            y = ph[:, None] - self.t[None, :]
            near = np.abs(y) <= self.trunc
            ch = chi_derivs(np.where(near, y, 0.0), self.sigma, kmax)
            for k in range(kmax+1):
                out[k][i0:i0+chunk] = (np.where(near, ch[k], 0.0)*self.f[None, :]).sum(axis=1)
        return out

def A_derivs_from_W(phi, Wd, kmax=4):
    """A_sigma^{(k)} = (W_sigma sin)^{(k)} by Leibniz."""
    s = np.sin(phi)
    c = np.cos(phi)
    sd = [s, c, -s, -c, s]
    return [sum(math.comb(k, j)*Wd[j]*sd[k-j] for j in range(k+1)) for k in range(kmax+1)]


# ---------------------------------------------------------------- controls (L-98)
if __name__ == "__main__":
    import json, os, sys, warnings
    warnings.filterwarnings("ignore")
    HERE = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(HERE, "copies", "refutescripts"))
    sys.path.insert(0, os.path.join(HERE, "copies", "fix4"))
    import s1_profile_indep as S1                       # referee's adaptive quadrature
    import p1_profile as P1                             # fix4's FFT grid
    C = Conv(SIGMA_DEFAULT)
    OUT = {"instrument": "fixed-node Gauss-Legendre convolution in the integration variable t, "
                         "panels split at the kink images, 24 points per panel, panel width "
                         "<= sigma/2",
           "n_nodes": int(C.t.size), "sigma": SIGMA_DEFAULT,
           "bounds": {k: v for k, v in W_bounds(SIGMA_DEFAULT).items()}}

    ANG = [0.05, 0.1309, 0.2109, 0.3, 0.5236, 0.9, 1.2, 1.4835, 2.0, 3.0]
    rows = []
    D = C.derivs(np.array(ANG), kmax=3)
    for i, p in enumerate(ANG):
        ref = [S1.Ws(p), S1.Ws1(p), S1.Ws2(p), S1.Ws3(p)]
        rows.append({"phi": p,
                     "rel": [abs(D[k][i] - ref[k])/max(abs(ref[k]), 1e-300) for k in range(4)]})
    worst = max(max(r["rel"]) for r in rows)
    OUT["control_vs_referee_adaptive_quadrature"] = {
        "angles": ANG, "rows": rows, "worst_relative": worst,
        "orders": "k = 0, 1, 2, 3", "n_angles": len(ANG)}

    fd = []
    for p in [0.2109, 0.3, 0.5236, 0.9, DELTA]:
        h = 1e-4
        x = np.array([p - 2*h, p - h, p + h, p + 2*h])
        Dx = C.derivs(x, kmax=3)
        val = (-Dx[3][3] + 8*Dx[3][2] - 8*Dx[3][1] + Dx[3][0])/(12*h)
        direct = float(C.derivs(np.array([p]), kmax=4)[4][0])
        fd.append({"phi": p, "W4_direct": direct, "W4_from_FD_of_W3": float(val),
                   "rel": abs(direct - val)/max(abs(val), 1e-300)})
    OUT["control_W4_vs_finite_difference_of_W3"] = {
        "rows": fd, "worst_relative": max(r["rel"] for r in fd)}

    s3 = json.load(open(os.path.join(HERE, "copies", "refutescripts", "s3_results.json")))
    mine = float(C.derivs(np.array([DELTA]), kmax=4)[4][0])
    OUT["control_W4_at_delta_vs_referee_spline"] = {
        "here": mine, "referee_spline_table": s3["W4_at_delta_table"],
        "referee_leading_prediction_Delta_times_chi2_0": s3["W4_at_delta_predicted"],
        "rel_to_referee_table": abs(mine - s3["W4_at_delta_table"])/abs(s3["W4_at_delta_table"])}

    prof = P1.Profile(SIGMA_DEFAULT, nphi=600001)
    br = np.where(81.0*prof.w0**2 <= 8.0*prof.w1**2, prof.w0**2 + prof.w1**2,
                  P1.YMAX*prof.w0**2 + prof.w1**2)
    k = int(np.argmax(br))
    ph = float(prof.phi[k])
    Dh = C.derivs(np.array([ph]), kmax=1)
    OUT["control_vs_fix4_FFT_at_the_Gfrak0_argmax"] = {
        "phi": ph, "W0_fix4_FFT": float(prof.W0[k]), "W0_here": float(Dh[0][0]),
        "W0_rel": abs(float(prof.W0[k]) - float(Dh[0][0]))/abs(float(Dh[0][0])),
        "W1_fix4_FFT": float(prof.W1[k]), "W1_here": float(Dh[1][0]),
        "W1_referee_adaptive_quadrature": float(S1.Ws1(ph)),
        "W1_abs_difference_fix4_minus_here": float(prof.W1[k]) - float(Dh[1][0]),
        "W1_rel": abs(float(prof.W1[k]) - float(Dh[1][0]))/abs(float(Dh[1][0])),
        "note": "W_sigma is unaffected (Wtil is continuous); W_sigma' carries a first-order "
                "error because Wtil' jumps at the kink cones and the FFT is a Riemann sum "
                "through the jump"}
    prof2 = P1.Profile(SIGMA_DEFAULT, nphi=2400001)
    OUT["control_vs_fix4_FFT_at_the_Gfrak0_argmax"]["W1_fix4_FFT_nphi_2400001"] = \
        float(np.interp(ph, prof2.phi, prof2.W1))
    OUT["control_vs_fix4_FFT_at_the_Gfrak0_argmax"]["error_ratio_on_fourfold_refinement"] = (
        (float(prof.W1[k]) - float(Dh[1][0]))
        / (float(np.interp(ph, prof2.phi, prof2.W1)) - float(Dh[1][0])))

    json.dump(OUT, open(os.path.join(HERE, "fx_controls.json"), "w"), indent=1, default=str)
    print("worst relative against the referee's adaptive quadrature, k = 0..3 : %.3e" % worst)
    print("worst relative, W4 direct against a finite difference of W3        : %.3e"
          % OUT["control_W4_vs_finite_difference_of_W3"]["worst_relative"])
    print("W4(delta) here %.4f, referee's spline %.4f, relative %.3e"
          % (mine, s3["W4_at_delta_table"],
             OUT["control_W4_at_delta_vs_referee_spline"]["rel_to_referee_table"]))
    g = OUT["control_vs_fix4_FFT_at_the_Gfrak0_argmax"]
    print("at the Gfrak_0 argmax phi = %.8f : W0 rel %.2e, W1 rel %.2e, "
          "fix4 error ratio on a fourfold refinement %.4f"
          % (g["phi"], g["W0_rel"], g["W1_rel"], g["error_ratio_on_fourfold_refinement"]))
    print("wrote fx_controls.json")

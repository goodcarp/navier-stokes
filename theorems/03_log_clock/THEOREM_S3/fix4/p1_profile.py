"""
p1 -- THE SIGMA-MOLLIFIED ANGULAR PROFILE, and every datum constant that depends on it.

THE DEFINITION (fix3 sec.1.6 / pmax-h3v Proposition 3, stated exactly).

  Kinked profile (THEOREM_S3 sec.1.1):
      A(phi) = sgn(cos phi) * min(1, phi_ax/delta) * min(1, |phi - pi/2|/delta_m) ,
      W(phi) = A(phi)/sin phi ,          eta_0 = -(M/rho) Theta(u) W(phi) .

  W is EVEN about phi = 0 and about phi = pi (A is odd there and so is sin phi), so the
  geometrically correct extension of W to the line is the even reflection

      Wtil(phi) = W(|phi|) for |phi| <= pi ,  Wtil(2 pi - phi) = W(phi) ,  2 pi periodic.

  With the Gaussian chi_sigma(y) = exp(-y^2/(2 sigma^2))/(sigma sqrt(2 pi)),

      W_sigma := Wtil * chi_sigma   (restricted to [0, pi]) ,     A_sigma := W_sigma sin phi ,
      N_sigma := sup_{[0,pi]} |A_sigma| ,

  and THE MOLLIFIED DATUM of this sheet is

      omega_0^theta := -M Theta(rho) A_sigma(phi)/N_sigma ,
      eta_0         := omega_0^theta/r = -(M/rho) Theta(u) W_sigma(phi)/N_sigma .

  W_sigma is real-analytic on [0, pi] (Gaussian convolution of a bounded function),
  W_sigma'(0) = W_sigma'(pi) = 0 by evenness, so eta_0 is C^infinity in the angle including
  the axis and the four kink cones are gone.  ||omega_0^theta||_inf = M sup Theta by
  construction.

WHY THE NORMALISATION N_sigma IS NOT OPTIONAL.
  On the bulk delta < phi < pi/2 - delta_m the kinked W is exactly 1/sin phi, which is CONVEX,
  so W_sigma >= W there by Jensen and A_sigma = W_sigma sin phi >= 1: mollifying W OVERSHOOTS
  the amplitude.  The theorem's normalisation is ||omega_0^theta||_inf = M (T_d is the first
  time ||omega||_inf reaches (3/2)M, and lam_om = M_s/M <= 3/2), so the overshoot must be
  divided out.  fix3 sec.1.6's table is the UNNORMALISED profile; its
  kappa_delta(0.05) = 0.4983883 > 0.4978224 is the overshoot, not a gain (sec. "unnormalised"
  in the output reproduces that table and prints N_sigma beside it).

CONTROLS.  With sigma = None the same code returns the kinked profile and must reproduce
fix2 f4's closed forms E_0 = 1/sin delta, Gfrak_0 = 1/sin^2 delta, kappa_delta = 0.4978224182,
and fix3 g1's unnormalised sigma table.

Outputs -> p1_results.json.  Imported as a module by p3, p4, p5.
"""
import json, math, os
import numpy as np
from scipy.signal import fftconvolve
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
DEG = math.pi/180.0
DELTA, DM = 7.5*DEG, 5.0*DEG
SD, CDM = math.sin(DELTA), math.cos(DM)
EPS_R = 0.25
PHI0 = 30.0*DEG
LOG32 = math.log(1.5)
SIGMAS = [0.05, 0.02, 0.01, 0.005, 0.002]
CAP_CERT_FIX2 = 1.1943662078602755      # fix2 sec.1.3, for reference only

# ------------------------------------------------------------------ the kinked profile
# W and W' are evaluated on [0, pi/2] by the STABLE closed form on each analytic piece and
# extended by the exact symmetry  W(pi - phi) = -W(phi),  W'(pi - phi) = +W'(phi)
# (A is odd about pi/2 and sin is even about it).  The naive (A' s - A c)/s^2 loses all its
# digits within 1e-4 of either pole: at phi = pi - 1e-9 it returns 935.56 instead of the
# true 0, and that spurious spike then leaks into every convolution.
def _WWp_half(phi):
    """W, W' on (0, pi/2], phi an array."""
    s, c = np.sin(phi), np.cos(phi)
    W = np.empty_like(phi)
    Wp = np.empty_like(phi)
    tap = phi < DELTA
    eqz = phi > (math.pi/2 - DM)
    bulk = ~tap & ~eqz
    # taper: A = phi/delta
    ph = phi[tap]
    sm = ph < 1e-3
    ss, cc = np.sin(ph), np.cos(ph)
    num = np.where(sm, ph**3/3.0*(1.0 - ph**2/10.0), ss - ph*cc)      # sin - phi cos
    Wt = np.where(sm, (1.0 + ph**2/6.0 + 7.0*ph**4/360.0)/DELTA, ph/(DELTA*ss))
    W[tap] = Wt
    Wp[tap] = num/(DELTA*ss**2)
    # bulk: A = 1
    W[bulk] = 1.0/s[bulk]
    Wp[bulk] = -c[bulk]/s[bulk]**2
    # equatorial: A = (pi/2 - phi)/delta_m
    ph = phi[eqz]
    t = math.pi/2 - ph
    ss, cc = np.sin(ph), np.cos(ph)
    W[eqz] = t/(DM*ss)
    Wp[eqz] = (-ss - t*cc)/(DM*ss**2)
    return W, Wp

def _WWp(phi):
    phi = np.asarray(phi, dtype=float)
    half = np.minimum(phi, math.pi - phi)
    left = phi <= math.pi/2
    W, Wp = _WWp_half(np.where(left, phi, math.pi - phi))
    return np.where(left, W, -W), Wp

def A_kink(phi):
    pax = np.minimum(phi, math.pi - phi)
    eq = np.abs(phi - math.pi/2)
    return np.sign(np.cos(phi))*np.minimum(1.0, pax/DELTA)*np.minimum(1.0, eq/DM)

def A_kink_p(phi):
    pax = np.minimum(phi, math.pi - phi)
    dpax = np.where(phi < math.pi/2, 1.0, -1.0)
    eq = np.abs(phi - math.pi/2)
    deq = np.where(phi > math.pi/2, 1.0, -1.0)
    g1 = np.minimum(1.0, pax/DELTA)
    dg1 = np.where(pax < DELTA, dpax/DELTA, 0.0)
    g2 = np.minimum(1.0, eq/DM)
    dg2 = np.where(eq < DM, deq/DM, 0.0)
    return np.sign(np.cos(phi))*(dg1*g2 + g1*dg2)

def W_kink(phi):
    return _WWp(phi)[0]

def Wp_kink(phi):
    return _WWp(phi)[1]

# ------------------------------------------------------------------ the mollified profile
class Profile:
    """W_sigma, W_sigma', W_sigma'' on a uniform phi grid; the normalised A_sigma and h_sigma."""
    def __init__(self, sigma, nphi=600001, pad=1e-9):
        self.sigma = sigma
        self.nphi = nphi
        self.phi = np.linspace(pad, math.pi - pad, nphi)
        h = self.phi[1] - self.phi[0]
        self.h = h
        s, c = np.sin(self.phi), np.cos(self.phi)
        if sigma is None:
            self.W0 = W_kink(self.phi)
            self.W1 = Wp_kink(self.phi)
            self.W2 = (-2.0*A_kink_p(self.phi)*c/s**2 + A_kink(self.phi)/s
                       + 2.0*A_kink(self.phi)*c*c/s**3)         # a.c. part (A'' = 0)
        else:
            nker = int(math.ceil(8.0*sigma/h))
            lo = self.phi[0] - h*np.arange(nker, 0, -1)
            hi = self.phi[-1] + h*np.arange(1, nker+1)
            full = np.concatenate([lo, self.phi, hi])
            ref = np.where(full < 0.0, -full, np.where(full > math.pi, 2*math.pi - full, full))
            sgn = np.where((full < 0.0) | (full > math.pi), -1.0, 1.0)   # W even => W' odd
            ref = np.clip(ref, 1e-14, math.pi - 1e-14)
            W = W_kink(ref)
            Wp = Wp_kink(ref)*sgn
            y = h*np.arange(-nker, nker+1)
            chi = np.exp(-0.5*(y/sigma)**2)
            chi /= chi.sum()*h
            chip = -(y/sigma**2)*chi
            sl = slice(nker, nker + nphi)
            self.W0 = (fftconvolve(W, chi, mode='same')*h)[sl]
            self.W1 = (fftconvolve(Wp, chi, mode='same')*h)[sl]
            self.W2 = (fftconvolve(Wp, chip, mode='same')*h)[sl]
        self.A0 = self.W0*s
        self.N = float(np.abs(self.A0).max())
        self.argN_deg = float(self.phi[int(np.argmax(np.abs(self.A0)))]/DEG)
        self.w0 = self.W0/self.N
        self.w1 = self.W1/self.N
        self.w2 = self.W2/self.N
        self.a0 = self.A0/self.N
        self.hprof = np.abs(self.a0)
        self.hprof_p = np.sign(self.a0)*(self.w1*s + self.w0*c)
        if sigma is None:
            self.angL = (A_kink_p(self.phi)*c*s + A_kink(self.phi)*(s*s - c*c))/s**3/self.N
        else:
            self.angL = self.w2 + 3.0*(c/s)*self.w1

    def interp(self, arr, x):
        return np.interp(x, self.phi, arr)

# ------------------------------------------------------------------ P_h, kappa, r_h
def _prim_closed(V, Kv):
    """int_0^V v^3 (1+K v^2)^{-5/2} dv = K^{-2}[2/3 + (1/3)S^{-3/2} - S^{-1/2}], S = 1+K V^2."""
    Kv = np.asarray(Kv, dtype=float)
    Sv = 1.0 + Kv*V*V
    with np.errstate(divide='ignore', invalid='ignore'):
        val = (2.0/3.0 + Sv**-1.5/3.0 - Sv**-0.5)/Kv**2
    small = Kv < 1e-9
    if np.ndim(val) == 0:
        return float(V**4/4.0) if bool(small) else float(val)
    return np.where(small, V**4/4.0, val)

NG = 160
_x, _w = np.polynomial.legendre.leggauss(NG)
def _panel(a, b):
    return 0.5*(b-a)*_x + 0.5*(a+b), 0.5*(b-a)*_w
_PAN = [_panel(0.0, SD), _panel(SD, CDM), _panel(CDM, 1.0)]
_VQ = np.concatenate([p[0] for p in _PAN])
_WQ = np.concatenate([p[1] for p in _PAN])

class PhInstrument:
    """Q(lam) = P_h(lam)/lam = E_{mu_lam}[H], H(v) = h_sigma(arcsin v)  (fix2 sec.4)."""
    def __init__(self, prof):
        self.p = prof
        ph = np.arcsin(np.clip(_VQ, 0.0, 1.0))
        self.H = prof.interp(prof.hprof, ph)
        self.wH = _WQ*self.H*_VQ**2
        self.v2 = _VQ**2

    def Q(self, l):
        D = 1.0 + (l**6 - 1.0)*self.v2
        return 3.0*l**9*float(np.dot(self.wH, D**-2.5))

    def Qvec(self, ls):
        D = 1.0 + (ls[:, None]**6 - 1.0)*self.v2[None, :]
        return 3.0*ls**9*(self.wH[None, :]*D**-2.5).sum(axis=1)

    def P_h(self, l):
        return l*self.Q(l)

def lipschitz_hpp(prof):
    """rigorous bound on |h_sigma''| = |(w2 sin + 2 w1 cos - w0 sin)| away from the single
       sign change of A_sigma at phi = pi/2 (a corner of h, where h' jumps sign)."""
    s, c = np.sin(prof.phi), np.cos(prof.phi)
    m = (prof.phi > 1e-3) & (prof.phi < math.pi - 1e-3)
    # on the two end caps the profile is exactly h = phi_ax/delta (kinked) or its Gaussian
    # average, and h'' -> 0 there; the caps are excluded because the closed form for W'' in
    # the (rho,phi) variables loses its digits at the poles.
    return float(np.max(np.abs((prof.w2*s + 2.0*prof.w1*c - prof.w0*s)[m])))

_RHCACHE = {}
def certified_r_h(inst, prof, lam_max, nsub=5000, lip_pad=None):
    """rigorous LOWER bound on inf_{[1,lam_max]} Q by Lipschitz + grid.
       |Q'(lam)| <= 9 lam^8 sup|h'| int_0^1 v^3 D^{-5/2} dv, closed form, decreasing in lam."""
    key = (id(prof), round(lam_max, 10), nsub)
    if key in _RHCACHE:
        return _RHCACHE[key]
    if lip_pad is None:
        lip_pad = lipschitz_hpp(prof)*(math.pi/(prof.nphi - 1))
    ls = np.linspace(1.0, lam_max, nsub+1)
    Qs = inst.Qvec(ls)
    Ibar = _prim_closed(1.0, ls**6 - 1.0)
    sup_hp = float(np.max(np.abs(prof.hprof_p))) + lip_pad
    step = (lam_max - 1.0)/nsub
    lip = 9.0*ls[1:]**8*sup_hp*Ibar[:-1]
    lower = Qs[:-1] - lip*step
    k = int(np.argmin(Qs))
    out = {"r_h_grid_min": float(Qs.min()), "argmin_lam": float(ls[k]),
           "r_h_certified_lower": float(lower.min()), "max_Lipschitz_bound": float(lip.max()),
           "sup_abs_hprime": sup_hp, "nsub": nsub, "step": step, "lam_max": lam_max,
           "Q_at_1": float(Qs[0]), "Q_at_lam_max": float(Qs[-1])}
    _RHCACHE[key] = out
    return out

def _cell_bounds(prof, ncell):
    edges = np.linspace(0.0, 1.0, ncell+1)
    ph_e = np.arcsin(np.clip(edges, 0.0, 1.0))
    idx = np.clip(np.searchsorted(prof.phi, ph_e), 0, prof.nphi-1)
    Mi = np.empty(ncell); mi = np.empty(ncell)
    hp = prof.hprof_p
    for i in range(ncell):
        a, b = idx[i], max(idx[i+1], idx[i]+1)
        seg = hp[a:b+1]
        Mi[i] = seg.max(); mi[i] = seg.min()
    return edges, Mi, mi

def Ph_prime_lower(prof, lam, r_h_lower, lip_pad, ncell=4000, cache={}):
    """P_h' = Q + lam Q' >= r_h_lower - 9 lam^9 sum_i max(M_i + pad, 0) I_i(lam) ,
       I_i = int_cell v^3 D^{-5/2} dv (closed form) >= int_cell v^3 sqrt(1-v^2) D^{-5/2} dv."""
    key = (id(prof), ncell)
    if key not in cache:
        cache[key] = _cell_bounds(prof, ncell)
    edges, Mi, mi = cache[key]
    K = lam**6 - 1.0
    I = _prim_closed(edges[1:], K) - _prim_closed(edges[:-1], K)
    pos = np.maximum(Mi + lip_pad, 0.0)
    return float(r_h_lower - 9.0*lam**9*float(np.dot(pos, I)))

def Ph_prime_value(inst, prof, lam, ncell=20000):
    edges = np.linspace(0.0, 1.0, ncell+1)
    vm = 0.5*(edges[1:] + edges[:-1])
    D = 1.0 + (lam**6 - 1.0)*vm**2
    hp = prof.interp(prof.hprof_p, np.arcsin(np.clip(vm, 0.0, 1.0)))
    w = vm**3*np.sqrt(np.maximum(1.0 - vm**2, 0.0))*D**-2.5
    return inst.Q(lam) - 9.0*lam**9*float(np.dot(hp*w, np.diff(edges)))

# ------------------------------------------------------------------ Theta and the u-scan
def Theta_u(u, k=0, L=40.0, er=EPS_R):
    a = (u - er)/er
    b = (u - L + er)/er
    if k == 0:
        return 0.5*(np.tanh(a) - np.tanh(b))
    ca = 1.0/np.cosh(np.clip(a, -350, 350))**2
    cb = 1.0/np.cosh(np.clip(b, -350, 350))**2
    if k == 1:
        return 0.5*(ca - cb)/er
    ta, tb = np.tanh(a), np.tanh(b)
    if k == 2:
        return (tb*cb - ta*ca)/er**2
    raise ValueError

YMAX = (4.5*4.5)**2/64.0          # max_{y in [1,9]} (9-y)^2 y^2/64   (fix2 sec.4(c))

def _hull_max(px, py, ax, ay):
    """max over the sampled phi of ax*px + ay*py, computed on the convex hull of (px,py)."""
    return float(np.max(ax*px + ay*py))

def _hull(px, py):
    pts = np.column_stack([px, py])
    hull = ConvexHull(pts, qhull_options="Qx")
    return pts[hull.vertices]

def datum_norms(prof, L=40.0, nu=1601):
    """E_0, Gfrak_0 (certified branch bound and direct scan), ||grad eta_0||_inf, Psi_sigma(0).
       The u-scan uses convex hulls in phi so that each u costs O(#hull vertices)."""
    phi, s, c = prof.phi, np.sin(prof.phi), np.cos(prof.phi)
    W0, W1 = prof.w0, prof.w1
    P, S = W0**2, W1**2
    branch = np.where(81.0*P <= 8.0*S, P + S, YMAX*P + S)
    kb = int(np.argmax(branch))
    G0_cert = math.sqrt(float(branch.max()))
    HG = _hull(P, S)                       # for sqrt(a^2 P + b^2 S)
    mask = (phi > 1e-3) & (phi < math.pi - 1e-3)
    HP = _hull(W0[mask], prof.angL[mask])  # for |alpha W0 + beta angL|
    supW = float(np.abs(W0).max())
    uu = np.concatenate([np.linspace(-6.0, 6.0, nu),
                         np.linspace(6.0, L - 6.0, 101),
                         np.linspace(L - 6.0, L + 6.0, nu)])
    supE = 0.0; supG = 0.0; argG = None; supPsi = 0.0; argPsi = None; gradinf = 0.0
    for u in uu:
        T0 = float(Theta_u(u, 0, L)); T1 = float(Theta_u(u, 1, L)); T2 = float(Theta_u(u, 2, L))
        e2 = math.exp(-2.0*u)
        supE = max(supE, supW*T0)
        g2 = float(np.max((T0 - T1)**2*HG[:, 0] + T0**2*HG[:, 1]))
        m = math.sqrt(max(g2, 0.0))
        if m > supG:
            supG, argG = m, float(u)
        gradinf = max(gradinf, m*e2)
        rad = T2 + T1 - 2.0*T0
        v = float(np.max(np.abs(-rad*HP[:, 0] - T0*HP[:, 1])))*e2
        if v > supPsi:
            supPsi, argPsi = v, float(u)
    return {"E0_scan": supE, "Gfrak0_scan": supG, "Gfrak0_argmax_u": argG,
            "Gfrak0_certified_branch_bound": G0_cert,
            "Gfrak0_branch_argmax_phi_deg": float(phi[kb]/DEG),
            "grad_eta0_inf": gradinf, "Psi0": supPsi, "Psi0_argmax_u": argPsi,
            "sup_abs_W": supW, "sup_abs_Wp": float(np.abs(W1).max())}


def build(sigma, nphi=600001):
    prof = Profile(sigma, nphi=nphi)
    inst = PhInstrument(prof)
    return prof, inst

def summarise(sigma, nphi=600001, nsub=20000, verbose=True):
    prof, inst = build(sigma, nphi=nphi)
    kappa = 0.5*inst.P_h(1.0)
    c_star = LOG32/kappa
    lip_h = lipschitz_hpp(prof)
    lip_pad = lip_h*(math.pi/(nphi - 1))
    rec = {"sigma": sigma, "N_sigma": prof.N, "argmax_N_phi_deg": prof.argN_deg,
           "kappa_delta": kappa, "c_star": c_star, "c2": 2.0*LOG32/kappa,
           "eps_delta_1m2kappa": 1.0 - 2.0*kappa,
           "eps_floor": (1.0 - 2.0*kappa)/(2.0*kappa),
           "lam_max_at_c_star": math.exp(0.75*c_star),
           "sup_abs_hprime": float(np.max(np.abs(prof.hprof_p))),
           "lipschitz_hprime_bound": lip_h, "grid_pad_on_hprime": lip_pad}
    rec.update(datum_norms(prof))
    rec["E0"] = rec["E0_scan"]
    rec["Gfrak0"] = rec["Gfrak0_certified_branch_bound"]
    rec["C_kink"] = (sigma*rec["Psi0"]) if sigma is not None else None
    rec["A_sigma_at_phi0"] = float(prof.interp(prof.a0, np.array([PHI0]))[0])
    rec["dist_phi0_to_nearest_kink_rad"] = min(PHI0 - DELTA, math.pi/2 - DM - PHI0)
    rec["dist_phi0_in_sigmas"] = (rec["dist_phi0_to_nearest_kink_rad"]/sigma) if sigma else None
    rec["windows"] = {}
    for tag, x in [("own_c_star", 1.0), ("own_cap_fix2", CAP_CERT_FIX2)]:
        lm = math.exp(0.75*x*c_star)
        rh = certified_r_h(inst, prof, lm, nsub=nsub, lip_pad=lip_pad)
        rh["Ph_prime_lower_at_lam_max"] = Ph_prime_lower(prof, lm, rh["r_h_certified_lower"],
                                                         lip_pad)
        rh["Ph_increasing_certified"] = bool(rh["Ph_prime_lower_at_lam_max"] > 0.0)
        rh["Ph_prime_value_at_lam_max"] = Ph_prime_value(inst, prof, lm)
        rec["windows"][tag] = rh
    lo_, hi_ = 1.5, 3.2
    for _ in range(60):
        mid = 0.5*(lo_ + hi_)
        if Ph_prime_value(inst, prof, mid) > 0:
            lo_ = mid
        else:
            hi_ = mid
    rec["lam_mono"] = 0.5*(lo_ + hi_)
    rec["c_over_c_star_at_lam_mono"] = (4.0/3.0)*math.log(rec["lam_mono"])/c_star
    # the largest CERTIFIED monotone window
    def ok(x):
        lm = math.exp(0.75*x*c_star)
        rh = certified_r_h(inst, prof, lm, nsub=5000, lip_pad=lip_pad)
        return Ph_prime_lower(prof, lm, rh["r_h_certified_lower"], lip_pad) > 0
    a_, b_ = 1.0, 1.35
    for _ in range(30):
        m_ = 0.5*(a_ + b_)
        if ok(m_):
            a_ = m_
        else:
            b_ = m_
    rec["c_over_c_star_largest_CERTIFIED"] = a_
    rec["eps_cap_certified"] = a_ - 1.0
    if verbose:
        print("sigma=%-7s N=%.7f kappa=%.9f floor=%.6e E0=%.5f G0=%.5f Psi0=%.2f "
              "Ckink=%.4f cap=%.7f"
              % (str(sigma), prof.N, kappa, rec["eps_floor"], rec["E0"], rec["Gfrak0"],
                 rec["Psi0"], rec["C_kink"] if rec["C_kink"] else float('nan'),
                 rec["eps_cap_certified"]), flush=True)
    return rec


if __name__ == "__main__":
    OUT = {"definition": {
        "delta_deg": 7.5, "delta_m_deg": 5.0, "eps_r": EPS_R, "phi0_deg": 30.0,
        "mollifier": "Gaussian chi_sigma, std sigma; W even-reflected about phi=0 and phi=pi",
        "normalisation": "omega_0^theta = -M Theta A_sigma/N_sigma, N_sigma = sup|W_sigma sin phi|",
        "sigmas": SIGMAS, "nphi": 600001}}
    rows = {}
    r0 = summarise(None)
    rows["kinked"] = r0
    OUT["control_kinked"] = {
        "kappa_delta": r0["kappa_delta"], "fix2_f4_value": 0.4978224182,
        "E0": r0["E0"], "one_over_sin_delta": 1.0/SD,
        "Gfrak0_branch_bound": r0["Gfrak0"], "one_over_sin2_delta": 1.0/SD**2,
        "N_must_be_1": r0["N_sigma"], "grad_eta0_inf": r0["grad_eta0_inf"],
        "fix3_g1_grad_inf": 20.3339, "Psi0_ac_part": r0["Psi0"], "fix3_g1_Psi0_ac": 171.4257,
        "r_h_at_c_star": r0["windows"]["own_c_star"]["r_h_certified_lower"],
        "fix2_r_h_enclosure": [0.9186364112, 0.9186365334],
        "lam_mono": r0["lam_mono"], "fix2_lam_mono": 2.0769162,
        "eps_cap_certified": r0["eps_cap_certified"], "fix2_eps_cap": 0.1943662078602755}
    json.dump(OUT, open(os.path.join(HERE, "p1_results.json"), "w"), indent=1,
              sort_keys=True, default=str)
    for sg in SIGMAS:
        rows["sigma=%g" % sg] = summarise(sg)
        OUT["rows"] = rows
        json.dump(OUT, open(os.path.join(HERE, "p1_results.json"), "w"), indent=1,
                  sort_keys=True, default=str)
    # the unnormalised comparison, reproducing fix3 g1 sec.1.6's table
    unn = {}
    for sg in SIGMAS:
        pr = Profile(sg)
        v = np.linspace(0.0, 1.0, 400001)
        ph = np.arcsin(np.clip(v, 0.0, 1.0))
        hA = np.abs(pr.interp(pr.A0, ph))
        unn["sigma=%g" % sg] = {"kappa_unnormalised": 0.5*3.0*float(np.trapz(hA*v**2, v)),
                                "N_sigma": pr.N, "Gfrak0_unnormalised_scan": None}
    unn["fix3_g1_quoted_kappa"] = {"0.05": 0.4983882748, "0.02": 0.4979100241,
                                   "0.01": 0.4978442142, "0.005": 0.4978278606,
                                   "0.002": 0.4978232887}
    OUT["unnormalised_control_vs_fix3"] = unn
    OUT["rows"] = rows
    json.dump(OUT, open(os.path.join(HERE, "p1_results.json"), "w"), indent=1,
              sort_keys=True, default=str)
    print(json.dumps({k: v for k, v in OUT.items() if k != "rows"}, indent=1, default=str))

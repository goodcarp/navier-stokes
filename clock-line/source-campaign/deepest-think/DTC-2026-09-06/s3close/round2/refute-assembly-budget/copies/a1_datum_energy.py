"""
a1 -- the datum, its strain constant, its zonal coefficients, its energy, and the
      L <-> log Re_E dictionary.

Seat: s3close/assembly (DTC-2026-09-06 continuation).  Everything here is computed;
nothing is transcribed from another seat except where a row is explicitly labelled
'record' and used only as a comparison target.

Datum (mollified delta-tapered plateau family):

    omega^theta_0(rho,phi) = -M * Theta_{eps}(rho) * g_{delta,delta_m}(phi)
    g(phi) = sgn(cos phi) * min(1, phi_ax/delta) * min(1, |phi - pi/2|/delta_m)
    phi_ax = min(phi, pi-phi)
    Theta_eps(rho) = min(1, log(rho/rho0)/eps, log(R/rho)/eps)      (eps e-folds of ramp)

    support rho0 < |x| < R,  L = log(R/rho0),  rho0 = s sqrt(nu/M).

eta = omega^theta/r, and in the 5D lift  rho*eta = W(t) = -Theta*g(phi)*M/sqrt(1-t^2),
t = cos phi.  Zonal expansion  W(t) = sum_l H_l C_l^{3/2}(t), weight (1-t^2),
N_l = int_{-1}^1 (C_l^{3/2})^2 (1-t^2) dt = (l+1)(l+2)/(l+3/2).

Energy identity derived here (not quoted):  with -Delta_5 psi1 = eta,
    ||u||_2^2 = (1/pi) int psi1 eta dx5 ,   dx5 = 2 pi^2 rho^4 sin^3(phi) drho dphi,
    Green fn of  Psi'' + (4/rho)Psi' - l(l+3)Psi/rho^2 = -f  is
    G_l(rho,rho') = rho_<^l rho_>^{-(l+3)}/(2l+3)  against the weight rho'^3 drho',
hence
    ||u_0||_2^2 = 2 pi R^5 M^2 sum_l N_l Hhat_l^2 * 2 D_l /(2l+3),
    D_l = int int_{x<y} Theta(x)Theta(y) e^{(l+4)x - (l-1)y} dx dy,  x = log(rho/R).
For the sharp shell and L -> infinity, D_l = 1/(5(l+4)) exactly for every l >= 1.
"""
import json, math
import numpy as np
from scipy.special import eval_gegenbauer
from scipy.integrate import quad

OUT = {}

# ---------------------------------------------------------------- profile pieces
def h_taper(phi, delta):
    """axis taper, delta in radians; 1 in the bulk"""
    phi_ax = min(phi, math.pi - phi)
    return min(1.0, phi_ax/delta) if delta > 0 else 1.0

def eq_moll(phi, dm):
    """equatorial mollification factor (|.|), 1 in the bulk"""
    if dm <= 0:
        return 1.0
    return min(1.0, abs(phi - math.pi/2)/dm)

def gprof(phi, delta, dm):
    """|angular profile| on (0,pi/2]; sign handled separately (odd in z)"""
    return h_taper(phi, delta)*eq_moll(phi, dm)

# ---------------------------------------------------------------- N_l, C_l
def N_l(l):
    return (l+1.0)*(l+2.0)/(l+1.5)

def C32(l, t):
    return eval_gegenbauer(l, 1.5, t)

# ------------------------------------------------------------------ H_l
# ---- vectorised composite Gauss-Legendre on [a,b] -------------------------------
_GL = {}
def gl_nodes(n):
    if n not in _GL:
        _GL[n] = np.polynomial.legendre.leggauss(n)
    return _GL[n]

def comp_gauss(f, a, b, npan=200, n=20):
    """f must accept a 1-D numpy array"""
    xg, wg = gl_nodes(n)
    edges = np.linspace(a, b, npan+1)
    lo = edges[:-1][:, None]; hi = edges[1:][:, None]
    mid = 0.5*(lo+hi); half = 0.5*(hi-lo)
    pts = (mid + half*xg[None, :]).ravel()
    wts = (half*wg[None, :]).ravel()
    return float(np.dot(f(pts), wts))

def _gvec(t, delta, dm):
    t = np.clip(np.asarray(t, dtype=float), -1.0, 1.0)
    phi = np.arccos(t)
    phi_ax = np.minimum(phi, np.pi-phi)
    h = np.ones_like(phi) if delta <= 0 else np.minimum(1.0, phi_ax/delta)
    m = np.ones_like(phi) if dm <= 0 else np.minimum(1.0, np.abs(phi-np.pi/2)/dm)
    return h*m

_HCACHE = {}
def H_l(l, delta, dm, M=1.0, npan=200):
    """H_l = -(2M/N_l) int_0^1 g(arccos t) sqrt(1-t^2) C_l(t) dt   (l odd)"""
    if l % 2 == 0:
        return 0.0
    key = (l, round(delta, 12), round(dm, 12))
    if key in _HCACHE:
        return _HCACHE[key]
    brk = {0.0, 1.0}
    if delta > 0: brk.add(math.cos(delta))
    if dm > 0:    brk.add(math.sin(dm))
    brk = sorted(b for b in brk if 0.0 <= b <= 1.0)
    tot = 0.0
    for a, b in zip(brk[:-1], brk[1:]):
        if b <= a: continue
        np_here = max(20, int(npan*(1 + l/40.0)))
        tot += comp_gauss(lambda t: _gvec(t, delta, dm)*np.sqrt(np.maximum(0.0, 1.0-t**2))
                                    * C32(l, t), a, b, npan=np_here, n=20)
    val = -(2.0*M/N_l(l))*tot
    _HCACHE[key] = val
    return val

# ---------------------------------------------------------------- P_h(lambda), kappa, r_h
def P_h(lam, delta, dm):
    """P_h(lam) = 3 int_0^1 h(arcsin v) v^2 (A v^2 + B)^{-5/2} dv, A = lam^2-lam^-4, B = lam^-4"""
    A = lam**2 - lam**(-4)
    B = lam**(-4)
    def f(v):
        phi = math.asin(min(1.0, v))
        return gprof(phi, delta, dm)*v*v*(A*v*v + B)**(-2.5)
    brk = sorted({0.0, 1.0, math.sin(delta) if delta > 0 else 0.0,
                  math.cos(dm) if dm > 0 else 1.0})
    brk = [b for b in brk if 0.0 <= b <= 1.0]
    tot = 0.0
    for a, b in zip(brk[:-1], brk[1:]):
        if b > a:
            val, _ = quad(f, a, b, limit=200, epsabs=1e-13, epsrel=1e-13)
            tot += val
    return 3.0*tot

# ---------------------------------------------------------------- D_l for a radial ramp
def Theta(x, L, eps):
    """x = log(rho/R) in [-L,0]; ramp of eps e-folds at each end"""
    if eps <= 0:
        return 1.0
    return max(0.0, min(1.0, (x+L)/eps, (-x)/eps))

def Theta_a(a, L, eps):
    """same ramp in the variable a = -x = log(R/rho) in [0,L]"""
    if eps <= 0:
        return np.ones_like(np.asarray(a, dtype=float))
    a = np.asarray(a, dtype=float)
    return np.clip(np.minimum((L-a)/eps, a/eps), 0.0, 1.0)

def D_l(l, L, eps, npan=60):
    """
    D_l = int_0^L db Theta(b) e^{-5b} int_0^{L-b} Theta(b+w) e^{-(l+4)w} dw
    (a = -x, b = -y, a > b; exponent (l+4)x-(l-1)y = -(l+4)(a-b)-5b).  All
    exponentials <= 1, so this is overflow-free at every L and l.
    """
    def inner(bvals):
        out = np.empty_like(bvals)
        for i, b in enumerate(bvals):
            wmax = min(L-b, 60.0/(l+4.0))
            if wmax <= 0:
                out[i] = 0.0
                continue
            brk = sorted({0.0, wmax} | ({L-eps-b} if eps > 0 and 0 < L-eps-b < wmax else set())
                                     | ({eps-b} if eps > 0 and 0 < eps-b < wmax else set()))
            tot = 0.0
            for u0, u1 in zip(brk[:-1], brk[1:]):
                tot += comp_gauss(lambda w: Theta_a(b+w, L, eps)*np.exp(-(l+4.0)*w),
                                  u0, u1, npan=npan, n=20)
            out[i] = tot
        return out
    brk = sorted({0.0, L} | ({eps, L-eps} if eps > 0 else set()))
    tot = 0.0
    for u0, u1 in zip(brk[:-1], brk[1:]):
        if u1 <= u0:
            continue
        tot += comp_gauss(lambda b: Theta_a(b, L, eps)*np.exp(-5.0*b)*inner(b),
                          u0, u1, npan=max(8, min(npan, 40)), n=12)
    return tot

def D_l_sharp_exact(l, L):
    """closed form for Theta == 1 on [-L,0]"""
    if abs(l-1.0) < 1e-12:
        return (1.0/5.0)*((1.0-math.exp(-5*L))/5.0 - L*math.exp(-5*L))
    return (1.0/(l+4.0))*((1.0-math.exp(-5*L))/5.0
                          - (math.exp(-5*L)-math.exp(-(l+4.0)*L))/(l-1.0))

# ---------------------------------------------------------------- energy
def C_E(delta, dm, eps, L, lmax=61, sharp_limit=False):
    """C_E = ||u0||_2^2/(M^2 R^5)"""
    tot = 0.0
    terms = {}
    for l in range(1, lmax+1, 2):
        Hh = H_l(l, delta, dm)
        if abs(Hh) < 1e-18:
            continue
        D = 1.0/(5.0*(l+4.0)) if sharp_limit else (
            D_l_sharp_exact(l, L) if eps <= 0 else D_l(l, L, eps))
        term = N_l(l)*Hh*Hh*2.0*D/(2*l+3.0)
        terms[l] = term
        tot += term
    return 2.0*math.pi*tot, terms

# ---------------------------------------------------------------- L2 norm of omega0
def omega_L2sq(delta, dm, eps, L, rho0=1.0):
    """ ||omega_0||_{L^2(R^3)}^2 = 2 pi int |omega|^2 r dr dz = 2 pi R^3 * (radial) * (angular) """
    R = rho0*math.exp(L)
    # int |omega|^2 dx = 2 pi int int omega^2 rho sin(phi) * rho drho dphi
    ang, _ = quad(lambda p: gprof(p, delta, dm)**2*math.sin(p), 0.0, math.pi/2,
                  limit=400, epsabs=1e-13)
    ang *= 2.0                                    # both hemispheres
    rad, _ = quad(lambda x: Theta(x, L, eps)**2*math.exp(3*x), -L, 0.0,
                  limit=400, epsabs=1e-14)        # rho = R e^x, rho^2 drho = R^3 e^{3x} dx
    return 2.0*math.pi*(R**3)*rad*ang

# ================================================================= run
if __name__ == "__main__":
    deg = math.pi/180.0

    # -- 0.  instrument checks -------------------------------------------------
    chk = {}
    chk["N_l_1_3_5"] = [N_l(1), N_l(3), N_l(5)]
    chk["C32_l1_at_t"] = [float(C32(1, 0.3)), 3*0.3]
    # exact bang-bang coefficients (record: H1=-5/6, H3=3/40, H5=-247/1680,
    # H7=1513/40320, H9=-2773/42240) -- recomputed here from scratch
    exact = {1: -5/6, 3: 3/40, 5: -247/1680, 7: 1513/40320, 9: -2773/42240}
    mine = {l: H_l(l, 0.0, 0.0) for l in exact}
    chk["Hl_bangbang_mine"] = mine
    chk["Hl_bangbang_exact"] = exact
    chk["Hl_bangbang_maxrel"] = max(abs(mine[l]-exact[l])/abs(exact[l]) for l in exact)
    chk["kappa0_from_H1"] = -0.6*mine[1]
    # P_h(1) = 1 and P_1(lam) = lam for the bare cap
    chk["P_cap_at_1"] = P_h(1.0, 0.0, 0.0)
    chk["P_cap_maxdev_from_lam"] = max(abs(P_h(lam, 0.0, 0.0)-lam)
                                       for lam in [1.0, 1.1, 1.25, 1.4, 1.5, 2.0, 3.0])
    # D_l numeric vs closed form (sharp)
    chk["D_l_num_vs_closed_maxrel"] = max(
        abs(D_l(l, 8.0, 0.0)-D_l_sharp_exact(l, 8.0))/D_l_sharp_exact(l, 8.0)
        for l in [1, 3, 5, 9, 21])
    OUT["checks"] = chk

    # -- 1.  strain constants kappa_delta, r_h --------------------------------
    kap = {}
    rh = {}
    for d in [0.0, 3.0, 5.0, 7.5, 10.0, 15.0, 30.0]:
        dd = d*deg
        kap["%g" % d] = 0.5*P_h(1.0, dd, 0.0)
        lams = np.linspace(1.0, 1.5, 51)
        rh["%g" % d] = float(min(P_h(l, dd, 0.0)/l for l in lams))
    OUT["kappa_delta"] = kap
    OUT["r_h"] = rh
    OUT["Phi_ratio_at_1p5"] = {"%g" % d: P_h(1.5, d*deg, 0.0)/P_h(1.0, d*deg, 0.0)
                               for d in [3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0]}

    # -- 2.  C_E : sharp cap, L -> infinity  (target: record 0.172403978) ------
    CE_sharp, terms = C_E(0.0, 0.0, 0.0, 1.0, lmax=201, sharp_limit=True)
    OUT["C_E_sharp_cap_Linf"] = CE_sharp
    OUT["C_E_sharp_cap_record"] = 0.172403978
    OUT["C_E_sharp_cap_relerr_vs_record"] = abs(CE_sharp-0.172403978)/0.172403978
    OUT["C_E_l1_share"] = 2*math.pi*terms[1]/CE_sharp
    OUT["C_E_terms_l1_to_l9"] = {str(l): terms[l] for l in [1, 3, 5, 7, 9]}

    # -- 3.  C_E for the tapered / mollified family ---------------------------
    tab = {}
    for d in [0.0, 5.0, 7.5, 15.0]:
        for dm in [0.0, 5.0, 10.0]:
            v, _ = C_E(d*deg, dm*deg, 0.0, 1.0, lmax=201, sharp_limit=True)
            tab["delta=%g,dm=%g,eps=0" % (d, dm)] = v
    for eps in [0.0, 0.1, 0.25, 0.5]:
        for L in [10.0, 40.0, 160.0]:
            v, _ = C_E(7.5*deg, 5.0*deg, eps, L, lmax=121)
            tab["delta=7.5,dm=5,eps=%g,L=%g" % (eps, L)] = v
    OUT["C_E_table"] = tab

    # -- 4.  the L <-> log Re_E dictionary ------------------------------------
    # Re_E = E0^{2/5} M^{1/5}/nu with E0 = ||u||_2^2 ;  nu = M rho0^2/s^2
    # => log Re_E = 2L + 2 log s + (2/5) log C_E     (exact, up to the O(L e^{-5L}) in C_E)
    dic = {}
    for d, dm in [(7.5, 5.0), (7.5, 0.0), (0.0, 0.0)]:
        CEv, _ = C_E(d*deg, dm*deg, 0.0, 1.0, lmax=201, sharp_limit=True)
        s = 1.0/math.sin(d*deg) if d > 0 else None
        dic["delta=%g,dm=%g" % (d, dm)] = {
            "C_E": CEv,
            "c_E_two_fifths_logC_E": 0.4*math.log(CEv),
            "s_record_1_over_sin_delta": s,
            "log_ReE_minus_2L_at_s_record": (2*math.log(s) + 0.4*math.log(CEv)) if s else None,
        }
    OUT["dictionary"] = dic
    # half-energy convention (the brief's E0 = 1/2 int |u|^2): shifts log Re_E by (2/5)log(1/2)
    OUT["half_energy_shift_in_logReE"] = 0.4*math.log(0.5)

    # -- 5.  BFG hypothesis check: omega0 in L^2 cap L^inf --------------------
    bf = {}
    bf["sup_omega0_over_M"] = 1.0     # attained in the bulk by construction
    for L in [10.0, 40.0]:
        bf["L2sq_over_M2R3_L=%g_eps=0" % L] = omega_L2sq(7.5*deg, 5.0*deg, 0.0, L)/(math.exp(L)**3)
        bf["L2sq_over_M2R3_L=%g_eps=0.25" % L] = omega_L2sq(7.5*deg, 5.0*deg, 0.25, L)/(math.exp(L)**3)
    OUT["bfg_hypotheses"] = bf

    # -- 6.  material point trajectory and its margins ------------------------
    mp = {}
    for phi0_deg in [20.0, 30.0, 40.0]:
        phi0 = phi0_deg*deg
        row = {}
        for lam in [1.0, 1.25, 1.5]:
            phit = math.atan2(lam**3*math.sin(phi0), math.cos(phi0))
            row["lam=%g" % lam] = {
                "phi_deg": phit/deg,
                "margin_to_equator_deg": 90.0 - phit/deg,
                "margin_to_taper_deg": phit/deg - 7.5,
            }
        mp["phi0=%g" % phi0_deg] = row
    OUT["material_point"] = mp

    with open("a1_results.json", "w") as f:
        json.dump(OUT, f, indent=1, sort_keys=True)
    print(json.dumps(OUT, indent=1, sort_keys=True)[:6000])

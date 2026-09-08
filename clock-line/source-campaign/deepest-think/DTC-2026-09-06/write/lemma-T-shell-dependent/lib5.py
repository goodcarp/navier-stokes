"""lib5 - own 2-D reduction of the 5-D strain functional at the origin, and the maps.

Every map used is SO(4)-equivariant, so it acts on the meridian half-plane only, and the 5-D
Jacobian of a meridian map (r,z) -> (r',z') is

        J_5 = (r'/r)^3 * det d(r',z')/d(r,z)                                   (*)

(the 3 directions of S^3 orthogonal to y inside R^4 are scaled by r'/r).  Verified against the
exact det D(Lambda) of p1 in p4.

All radii are carried as u = log(rho) so that L = 400 does not overflow.

a[eta_0 o Phi^-1](0) = INT_S Kcal(Phi x) eta_0(x) J_Phi(x) dx_5
     = (3M/4) L INT_0^1 dtau INT_0^pi sgn(cos phi) h(phi) cos(phi_f) e^{-4 (u_f - u)} J_Phi
                                      sin^2(phi) dphi
with Kcal(x) = -(3/(8 pi^2)) x_z/|x|^5, eta_0 = -M sgn(z) h(phi)/r, dx_5 = 2pi^2 rho^4 sin^3 drho dphi.
"""
import numpy as np

M = 1.0

def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b - a)*x + 0.5*(a + b), 0.5*(b - a)*w

def panels(pts, n):
    xs, ws = [], []
    for a, b in zip(pts[:-1], pts[1:]):
        x, w = gl(a, b, n); xs.append(x); ws.append(w)
    return np.concatenate(xs), np.concatenate(ws)

# ---------------------------------------------------------------- the integro-ODE profile
KTH = 2*(1 - np.sqrt(2.0/3.0))
def lam_sigma(sig):   return (1 - (1 - sig)*KTH/2)**-2
def dloglam_dsig(sig): return -KTH/(1 - (1 - sig)*KTH/2)

# ---------------------------------------------------------------- stages
# a stage maps (u, phi) -> (u', phi') and returns its own J_5 factor.
def stage_Lambda(u, phi, L):
    """Lambda(x) = T_{lambda(rho)} x with lambda from the integro-ODE.  u = log(rho/rho0)."""
    sig = u/L
    lam = lam_sigma(sig)
    D = dloglam_dsig(sig)/L                      # dlog lambda / dlog rho
    s, c = np.sin(phi), np.cos(phi)
    rp, zp = lam*s, lam**-2*c                    # (r',z')/rho
    g = np.hypot(rp, zp)
    up = u + np.log(g)
    php = np.arctan2(rp, zp)                     # in (0, pi)
    J = lam**2*(1 + D*(s**2 - 2*c**2))
    return up, php, J, lam, D

def stage_ripple(u, phi, L, mu_a, k=1.0):
    """radial ripple  Psi(w) = s(|w|) w ,  s = 1 + mu_a sin(k pi u / L)."""
    s_ = 1 + mu_a*np.sin(k*np.pi*u/L)
    ds = mu_a*(k*np.pi/L)*np.cos(k*np.pi*u/L)    # ds/dlog|w|
    J = s_**4*(s_ + ds)
    return u + np.log(s_), phi, J

def stage_shear(u, phi, beta_amp):
    """angular shear  phi -> phi + beta_amp sin^2(phi)  (SO(4)-equivariant, C^inf)."""
    php = phi + beta_amp*np.sin(phi)**2
    dbp = beta_amp*np.sin(2*phi)
    J = (np.sin(php)/np.sin(phi))**3*(1 + dbp)
    return u, php, J

def apply_map(u, phi, L, spec):
    """spec: list of ('lambda',) / ('ripple', mu_a, k) / ('shear', beta) applied left to right."""
    Jtot = np.ones_like(u); uu, pp = u, phi
    lam = D = None
    for st in spec:
        if st[0] == 'lambda':
            uu, pp, J, lam, D = stage_Lambda(uu, pp, L)
        elif st[0] == 'ripple':
            uu, pp, J = stage_ripple(uu, pp, L, st[1], st[2])
        elif st[0] == 'shear':
            uu, pp, J = stage_shear(uu, pp, st[1])
        else:
            raise ValueError(st)
        Jtot = Jtot*J
    return uu, pp, Jtot, lam, D

# ---------------------------------------------------------------- profiles
def h_cap(phi):    return np.ones_like(phi)
def h_taper(phi, deg=7.5):
    d = np.deg2rad(deg); pa = np.minimum(phi, np.pi - phi)
    return np.minimum(1.0, pa/d)

# ---------------------------------------------------------------- the functional
def a_of_map(L, spec, h=h_cap, n_phi=200, n_tau=200, taper_deg=None):
    pts = [0.0, np.pi/2, np.pi]
    if taper_deg is not None:
        d = np.deg2rad(taper_deg); pts += [d, np.pi - d]
    ph, wph = panels(np.array(sorted(set(pts))), n_phi)
    ta, wta = gl(0.0, 1.0, n_tau)
    PH, TA = np.meshgrid(ph, ta, indexing='ij')
    WP, WT = np.meshgrid(wph, wta, indexing='ij')
    U = TA*L
    uf, pf, J, lam, D = apply_map(U, PH, L, spec)
    integ = np.sign(np.cos(PH))*h(PH)*np.cos(pf)*np.exp(-4*(uf - U))*J*np.sin(PH)**2
    return 0.75*M*L*float(np.sum(integ*WP*WT))

def hyp_constants(L, spec, n_u=2001, n_phi=1501):
    """mu  = sup |Phi x - Lambda x| / |Lambda x|
       muJ = sup |J_Phi - lambda(rho)^2| / lambda(rho)^2      (the SHEAR-FREE reference weight)"""
    u = np.linspace(0, L, n_u); ph = np.linspace(1e-9, np.pi - 1e-9, n_phi)
    U, PH = np.meshgrid(u, ph, indexing='ij')
    uL, pL, JL, lam, D = apply_map(U, PH, L, [('lambda',)])
    uf, pf, JF, _, _ = apply_map(U, PH, L, spec)
    q = np.exp(uf - uL)
    mu = float(np.max(np.sqrt(np.maximum(1 + q**2 - 2*q*np.cos(pf - pL), 0.0))))
    muJ = float(np.max(np.abs(JF/lam**2 - 1)))
    muJ_vs_JL = float(np.max(np.abs(JF/JL - 1)))
    return mu, muJ, muJ_vs_JL

def lemmaT_prime_bound(L, mu, muJ, lam_bar):
    """pi M A [ 3 mu (1+muJ)/(2(1-mu)^5) + 3 muJ/8 ],  A = L * lam_bar."""
    if mu >= 1:
        return float('inf')
    A = L*lam_bar
    return np.pi*M*A*(3*mu*(1 + muJ)/(2*(1 - mu)**5) + 3*muJ/8.0)

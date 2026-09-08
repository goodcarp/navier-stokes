"""Shared quadrature for the axis-strain functional a[eta0 o Phi^-1](0) in the 5D lift.

a[eta0 o Phi^-1](0) = int_S K(Phi(x)) eta0(x) J_Phi(x) dx_5 ,  K(w) = -(3/8pi^2) w_z/|w|^5.
With dx_5 = 2 pi^2 rho^4 sin^3(phi) drho dphi, eta0 = omega^theta/r, r = rho sin phi:

    a = -(3/4) int_{rho0}^{R} int_0^pi  (Z_Phi/|Phi|^5) omega^theta(phi) J_Phi rho^3 sin^2(phi) dphi drho
      = -(3/4) L int_0^1 int_0^pi (Z_Phi/|Phi|^5) omega^theta J_Phi rho^4 sin^2(phi) dphi dtau,
      rho = rho0 exp(L tau).

Maps are SO(4)-equivariant, given in the (r,z) half-plane, ALWAYS as Psi o T_lambda.
Each map object supplies: image (Rm,Zm), 5-D Jacobian J5, and the image-relative displacement
|Phi - T_lambda x| / |T_lambda x|.
"""
import numpy as np

# ---------- Gauss-Legendre on a list of panels ----------
def gl_nodes(panels, n):
    xs, ws = [], []
    x0, w0 = np.polynomial.legendre.leggauss(n)
    for a, b in zip(panels[:-1], panels[1:]):
        xs.append(0.5*(b-a)*x0 + 0.5*(a+b)); ws.append(0.5*(b-a)*w0)
    return np.concatenate(xs), np.concatenate(ws)

def geometric_panels(lo, hi, n_dec):
    """geometric refinement from hi down to lo (lo>0)"""
    return np.concatenate(([0.0], np.geomspace(lo, hi, n_dec)))

# ---------- data ----------
def omega_bangbang(phi, M=1.0):
    return -M*np.sign(np.cos(phi))

def omega_taper(phi, delta, M=1.0):
    phiax = np.minimum(phi, np.pi-phi)
    return -M*np.sign(np.cos(phi))*np.minimum(1.0, phiax/delta)

# ---------- maps ----------
class MapT:
    """Phi = T_lambda"""
    name = 'T_lambda'
    def __init__(self, lam): self.lam = lam
    def __call__(self, r, z):
        lam = self.lam
        ru, zu = lam*r, z/lam**2
        vr = np.hypot(ru, zu)
        J5 = np.full_like(r, lam**2)
        return ru, zu, J5, np.zeros_like(r)   # Rm, Zm, J5, relative displacement

class MapA:
    """Phi = s(varrho) * T_lambda x ,  s = 1 + mu sin(k pi log(varrho/rho0)/L)  (radial ripple).
    k = 1 (half period, the s3 map A) is NOT a null direction; k = 2 (full period) is -- see s5."""
    def __init__(self, lam, mu, rho0, L, k=1):
        self.lam, self.mu, self.rho0, self.L, self.k = lam, mu, rho0, L, k
        self.name = f'A_radial_ripple_k={k}'
    def __call__(self, r, z):
        lam, mu, L, k = self.lam, self.mu, self.L, self.k
        ru, zu = lam*r, z/lam**2
        vr = np.hypot(ru, zu)
        th = k*np.pi*np.log(vr/self.rho0)/L
        s  = 1 + mu*np.sin(th)
        vs = mu*(k*np.pi/L)*np.cos(th)          # varrho * s'(varrho)
        J5 = lam**2 * s**4 * (s + vs)
        return s*ru, s*zu, J5, np.abs(s-1.0)

class MapAng:
    """Phi = R_beta o T_lambda : rotate in the (r,z) half-plane by beta(phi_u), varrho unchanged."""
    def __init__(self, lam, beta, dbeta, name):
        self.lam, self.beta, self.dbeta, self.name = lam, beta, dbeta, name
    def __call__(self, r, z):
        lam = self.lam
        ru, zu = lam*r, z/lam**2
        vr  = np.hypot(ru, zu)
        phu = np.arctan2(ru, zu)                 # in (0, pi)
        b   = self.beta(phu); db = self.dbeta(phu)
        ph2 = phu + b
        J5  = lam**2 * (np.sin(ph2)/np.sin(phu))**3 * (1.0 + db)
        return vr*np.sin(ph2), vr*np.cos(ph2), J5, 2*np.abs(np.sin(b/2))

# ---------- the functional ----------
def a_of_map(mp, omega, lam, rho0, R, phi_panels, n_phi, n_tau, homogeneous):
    L = np.log(R/rho0)
    ph, wph = gl_nodes(phi_panels, n_phi)
    if homogeneous:
        tau, wtau = np.array([0.5]), np.array([1.0])
    else:
        tau, wtau = gl_nodes([0.0, 1.0], n_tau)
    rho = rho0*np.exp(L*tau)
    PH, RHO = np.meshgrid(ph, rho, indexing='ij')
    WPH, WT  = np.meshgrid(wph, wtau, indexing='ij')
    r = RHO*np.sin(PH); z = RHO*np.cos(PH)
    Rm, Zm, J5, rel = mp(r, z)
    nrm = np.hypot(Rm, Zm)
    integ = -(3.0/4.0)*(Zm/nrm**5)*omega(PH)*J5*RHO**4*np.sin(PH)**2
    a = L*np.sum(integ*WPH*WT)
    muJ = np.max(np.abs(J5/lam**2 - 1.0))
    return a, np.max(rel), muJ

def bound(lam, M, L, mu, muJ):
    return np.pi*lam*M*L*(3*mu*(1+muJ)/(2*(1-mu)**5) + 3*muJ/8.0)

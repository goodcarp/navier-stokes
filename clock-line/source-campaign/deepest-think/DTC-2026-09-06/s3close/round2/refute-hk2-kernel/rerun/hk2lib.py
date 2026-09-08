"""hk2lib -- the datum families, the flow geometry, and the exact eta-derivative fields.

Units throughout: rho0 = 1, M = 1.  x = (y,z), y in R^4, r = |y|, rho = |x|, t = z/rho.

(D-A)  IDEALISED (the analytic datum of V-b sec 1 / prove-lagrangian):
         omega^theta = -M sgn(z) h_delta(phi) 1_{rho0<rho<R},  h_delta = min(1, phi_ax/delta).
       In the plateau bulk (delta<phi<pi-delta, rho0<rho<R, z!=0) this is EXACTLY
         eta = eta_P = -M sgn(z)/r ,  |grad eta_P| = M/r^2,  ||Hess eta_P||_F = sqrt(7) M/r^3.

(D-B)  CAMPAIGN (lower/prove-duhamel/dcommon.py, write/V-b .../b2_gaussian_exact.py):
         omega^theta = -M tanh(sin phi/sin delta) tanh(cos phi/w) Theta(rho),
         Theta(rho) = (1/2)[tanh((rho-rho0)/w0) - tanh((rho-R)/w1)],  w0=.25 rho0, w1=.10 R.
       Record defaults: delta = 15 deg (dcommon) or 7.5 deg (V-b sizing), w = 0.20.
"""
import math
import numpy as np

# ------------------------------------------------------------------ (D-B) fields
class DatumB:
    """the campaign's mollified datum; all derivatives analytic (autodiff-free, hand-coded,
    checked against finite differences by k2)."""
    def __init__(self, delta_deg=7.5, w=0.20, R=None, L=10.0, w0=0.25, w1frac=0.10, M=1.0):
        self.sd = math.sin(math.radians(delta_deg)); self.w = w; self.M = M
        self.R = R if R is not None else math.exp(L)
        self.w0 = w0; self.w1 = w1frac*self.R
        self.delta_deg = delta_deg; self.L = math.log(self.R)

    # radial profile
    def Th(self, rho, k=0):
        a = (rho-1.0)/self.w0; b = (rho-self.R)/self.w1
        if k == 0: return 0.5*(np.tanh(a)-np.tanh(b))
        ca = 1.0/np.cosh(a)**2; cb = 1.0/np.cosh(b)**2
        if k == 1: return 0.5*(ca/self.w0 - cb/self.w1)
        ta = np.tanh(a); tb = np.tanh(b)
        if k == 2: return 0.5*(-2*ta*ca/self.w0**2 + 2*tb*cb/self.w1**2)
        raise ValueError

    # angular profile  A(t) = tanh(sqrt(1-t^2)/sd) * tanh(t/w)   (omega^theta = -M A Theta)
    # eta = -M A(t) Theta(rho) / (rho sqrt(1-t^2)) = -M Hh(t) Theta(rho)/rho
    def Hh(self, t, k=0):
        s = np.sqrt(np.maximum(1.0-t*t, 1e-300))
        u = s/self.sd
        # tanh(u)/s  = (1/sd) * tanh(u)/u  -- regular as s->0
        f = np.tanh(u)/np.maximum(u, 1e-300)/self.sd     # = tanh(s/sd)/s
        g = np.tanh(t/self.w)
        if k == 0: return f*g
        # d/dt : s' = -t/s ;  f = F(s) with F(s)=tanh(s/sd)/s
        Fp = (1.0/(self.sd*np.cosh(u)**2) - np.tanh(u)/s)/s      # dF/ds
        gp = 1.0/(self.w*np.cosh(t/self.w)**2)
        if k == 1: return Fp*(-t/s)*g + f*gp
        # second derivative
        Fpp = ((-2*np.tanh(u)/(self.sd**2*np.cosh(u)**2)) - 2*Fp - (-np.tanh(u)/s + 1/(self.sd*np.cosh(u)**2))/s)/s
        # recompute cleanly: F(s) = T(s)/s with T = tanh(s/sd)
        T   = np.tanh(u); Tp = 1.0/(self.sd*np.cosh(u)**2); Tpp = -2*T*Tp/self.sd
        F   = T/s; Fp = (Tp*s - T)/s**2; Fpp = (Tpp*s**2 - 2*Tp*s + 2*T)/s**3
        gpp = -2*np.tanh(t/self.w)/(self.w**2*np.cosh(t/self.w)**2)
        sp  = -t/s; spp = -(1.0/s + t*t/s**3)
        if k == 2: return (Fpp*sp**2 + Fp*spp)*g + 2*Fp*sp*gp + F*gpp
        raise ValueError

    # eta and its (rho,t) derivatives ; then (r,z) derivatives
    def eta_rz(self, r, z):
        """returns eta, eta_r, eta_z, eta_rr, eta_rz, eta_zz  (Cartesian-consistent, axisym)."""
        rho = np.hypot(r, z); t = z/rho
        s = r/rho
        Th0 = self.Th(rho,0); Th1 = self.Th(rho,1); Th2 = self.Th(rho,2)
        H0 = self.Hh(t,0);   H1 = self.Hh(t,1);   H2 = self.Hh(t,2)
        M = self.M
        # eta = -M H(t) Theta(rho)/rho =: G(rho) H(t) with G = -M Theta/rho
        G0 = -M*Th0/rho
        G1 = -M*(Th1/rho - Th0/rho**2)
        G2 = -M*(Th2/rho - 2*Th1/rho**2 + 2*Th0/rho**3)
        # chain rule rho(r,z), t(r,z)
        rr_ = s; rz_ = t                                  # d rho/dr, d rho/dz
        tr_ = -s*t/rho; tz_ = (1-t*t)/rho
        rrr = t*t/rho; rrz = -s*t/rho; rzz = s*s/rho
        trr = t*(2*s*s - t*t)/rho**2
        trz = s*(2*t*t - s*s)/rho**2
        tzz = -3*s*s*t/rho**2
        eta   = G0*H0
        eta_r = G1*rr_*H0 + G0*H1*tr_
        eta_z = G1*rz_*H0 + G0*H1*tz_
        eta_rr = G2*rr_**2*H0 + G1*rrr*H0 + 2*G1*rr_*H1*tr_ + G0*H2*tr_**2 + G0*H1*trr
        eta_rz = G2*rr_*rz_*H0 + G1*rrz*H0 + G1*(rr_*H1*tz_ + rz_*H1*tr_) + G0*H2*tr_*tz_ + G0*H1*trz
        eta_zz = G2*rz_**2*H0 + G1*rzz*H0 + 2*G1*rz_*H1*tz_ + G0*H2*tz_**2 + G0*H1*tzz
        return eta, eta_r, eta_z, eta_rr, eta_rz, eta_zz

    def grad_norm(self, r, z):
        e, er, ez, *_ = self.eta_rz(r, z)
        return np.sqrt(er*er + ez*ez)

    # ---- lean gradient magnitude (used by the kernel quadratures) ----
    def grad_norm_fast(self, r, z):
        """|grad eta| via |grad f|^2 = (d_rho f)^2 + (1-t^2)(d_t f)^2/rho^2 ; eta = G(rho) H(t)."""
        rho = np.sqrt(r*r + z*z); t = z/rho
        s = np.sqrt(np.maximum(1.0-t*t, 1e-300))
        a = (rho-1.0)/self.w0; b = (rho-self.R)/self.w1
        Th0 = 0.5*(np.tanh(a)-np.tanh(b))
        Th1 = 0.5*((1.0-np.tanh(a)**2)/self.w0 - (1.0-np.tanh(b)**2)/self.w1)
        G0 = -self.M*Th0/rho
        G1 = -self.M*(Th1/rho - Th0/rho**2)
        u = s/self.sd; T = np.tanh(u); Tp = (1.0-T*T)/self.sd     # dT/ds
        F = T/s; Fp = (Tp*s - T)/(s*s)
        g = np.tanh(t/self.w); gp = (1.0-g*g)/self.w
        H0 = F*g
        H1 = Fp*(-t/s)*g + F*gp
        return np.sqrt((G1*H0)**2 + (1.0-t*t)*(G0*H1)**2/(rho*rho))


    def hess_F(self, r, z):
        """Frobenius norm of the 5D Hessian of eta: block [[e_rr,e_rz],[e_rz,e_zz]] + (e_r/r) I3."""
        e, er, ez, err, erz, ezz = self.eta_rz(r, z)
        return np.sqrt(err**2 + 2*erz**2 + ezz**2 + 3*(er/r)**2)

    def hess_op(self, r, z):
        e, er, ez, err, erz, ezz = self.eta_rz(r, z)
        tr = err+ezz; disc = np.sqrt(np.maximum((err-ezz)**2 + 4*erz**2, 0.0))
        lam = np.maximum(np.abs(0.5*(tr+disc)), np.abs(0.5*(tr-disc)))
        return np.maximum(lam, np.abs(er/r))

# ------------------------------------------------------------------ (D-A) fields
class DatumA:
    """idealised: sharp sgn(z), corner taper h=min(1,phi_ax/delta), sharp radial edges."""
    def __init__(self, delta_deg=7.5, R=None, L=10.0, M=1.0):
        self.delta = math.radians(delta_deg); self.R = R if R is not None else math.exp(L)
        self.M = M; self.delta_deg = delta_deg; self.L = math.log(self.R)
    def h(self, phi):
        pax = np.minimum(phi, math.pi-phi)
        return np.minimum(1.0, pax/self.delta)
    def eta(self, r, z):
        rho = np.hypot(r,z); phi = np.arctan2(r, z)
        inside = (rho >= 1.0) & (rho <= self.R)
        return np.where(inside, -self.M*np.sign(z)*self.h(phi)/np.maximum(r,1e-300), 0.0)
    def grad_norm(self, r, z):
        """|grad eta| on the smooth pieces.  eta = -M sgn(z) Hh(phi)/rho, Hh = h(phi)/sin phi,
        so |grad eta| = (M/rho^2) sqrt(Hh^2 + Hh'^2)  (radial and angular parts)."""
        rho = np.hypot(r,z); phi = np.arctan2(r,z); s = np.sin(phi); c = np.cos(phi)
        pax = np.minimum(phi, math.pi-phi)
        h  = np.minimum(1.0, pax/self.delta)
        hp = np.where(pax < self.delta, np.where(phi < math.pi/2, 1.0, -1.0)/self.delta, 0.0)
        Hh  = h/s
        Hhp = (hp*s - h*c)/s**2
        inside = (rho >= 1.0) & (rho <= self.R)
        return np.where(inside, (self.M/rho**2)*np.sqrt(Hh**2 + Hhp**2), 0.0)

    def Hh_of_phi(self, phi):
        s = np.sin(phi); pax = np.minimum(phi, math.pi-phi)
        h = np.minimum(1.0, pax/self.delta)
        hp = np.where(pax < self.delta, np.where(phi < math.pi/2, 1.0, -1.0)/self.delta, 0.0)
        return h/s, (hp*s - h*np.cos(phi))/s**2

# ------------------------------------------------------------------ geometry / window
THETA_MAX = 4.0*(1.0-math.sqrt(2.0/3.0))            # lam(theta_max) = 3/2  (prove-lagrangian 4.1)
def lam_of_theta(th, kap=0.5): return (1.0-kap*th/2.0)**-2

def strain(x_rz, lam):
    """T_lam (r,z) -> (lam r, lam^-2 z)."""
    r,z = x_rz; return (lam*r, z/lam**2)

def dist_to_cone(r, z, ang):
    """distance in the (r,z) half plane from (r,z) to the cone phi = ang about the +z axis
    (equal to the R^5 distance, both sets being axisymmetric)."""
    n = (math.cos(ang), -math.sin(ang))              # unit normal to the ray (sin,cos)
    return abs(r*n[0] + z*n[1])

def dist_to_equator(r, z): return abs(z)

def geometry(phi0_deg=30.0, f=0.0, delta_deg=7.5, lam=1.0, R=None):
    """the tracked material point at label (rho=(1+f)rho0, phi0), pushed by T_lam, and its
    distances to the three layer families.  Taper cone angle transforms as tan -> lam^3 tan."""
    phi0 = math.radians(phi0_deg); rho0 = 1.0+f
    r0, z0 = rho0*math.sin(phi0), rho0*math.cos(phi0)
    r, z = strain((r0,z0), lam)
    rho = math.hypot(r,z)
    dlt = math.radians(delta_deg)
    dlt_lam = math.atan(lam**3*math.tan(dlt))        # the strained taper cone
    d_taper = dist_to_cone(r, z, dlt_lam)
    d_eq    = dist_to_equator(r, z)
    # inner radial edge: the material sphere rho=rho0 becomes the ellipsoid T_lam(S_1);
    # dist(T_lam x, T_lam S) >= sigma_min(T_lam) dist(x,S) = min(lam,lam^-2) * f
    d_in    = min(lam, lam**-2)*f
    return dict(r=r, z=z, rho=rho, phi_deg=math.degrees(math.atan2(r,z)),
                d_taper=d_taper, d_eq=d_eq, d_in_lower=d_in,
                taper_cone_deg=math.degrees(dlt_lam), lam=lam)


class StrainedB:
    """the campaign datum transported by the pure axisymmetric strain:  eta_lam = eta_0 o T_lam^{-1},
    T_lam(y,z) = (lam y, lam^-2 z),  so  eta_lam(r,z) = eta_0(r/lam, lam^2 z).
    (This is the inviscid transport of the scalar eta; omega^theta_lam = lam omega^theta_0 o T^{-1},
    reproducing L3v sec 1's ||omega^theta_lam||_inf = lam M.)  All derivatives by exact chain rule."""
    def __init__(self, base, lam):
        self.b = base; self.lam = float(lam); self.M = base.M; self.R = base.R
    def _pull(self, r, z):
        return r/self.lam, (self.lam**2)*z
    def eta_rz(self, r, z):
        lam = self.lam
        e, er, ez, err, erz, ezz = self.b.eta_rz(*self._pull(r, z))
        return (e, er/lam, ez*lam**2, err/lam**2, erz*lam, ezz*lam**4)
    def grad_norm(self, r, z):
        _, er, ez, *_ = self.eta_rz(r, z)
        return np.sqrt(er*er + ez*ez)
    def grad_norm_fast(self, r, z):
        lam = self.lam
        rp, zp = self._pull(r, z)
        rho = np.sqrt(rp*rp + zp*zp); t = zp/rho
        s = np.sqrt(np.maximum(1.0-t*t, 1e-300))
        B = self.b
        a = (rho-1.0)/B.w0; bb = (rho-B.R)/B.w1
        Th0 = 0.5*(np.tanh(a)-np.tanh(bb))
        Th1 = 0.5*((1.0-np.tanh(a)**2)/B.w0 - (1.0-np.tanh(bb)**2)/B.w1)
        G0 = -B.M*Th0/rho; G1 = -B.M*(Th1/rho - Th0/rho**2)
        u = s/B.sd; T = np.tanh(u); Tp = (1.0-T*T)/B.sd
        F = T/s; Fp = (Tp*s - T)/(s*s)
        g = np.tanh(t/B.w); gp = (1.0-g*g)/B.w
        H0 = F*g; H1 = Fp*(-t/s)*g + F*gp
        # gradient in the PULLED-BACK variables, then push forward
        er0 = G1*(rp/rho)*H0 + G0*H1*(-(rp/rho)*t/rho)
        ez0 = G1*t*H0 + G0*H1*(1.0-t*t)/rho
        return np.sqrt((er0/lam)**2 + (ez0*lam**2)**2)
    def hess_F(self, r, z):
        e, er, ez, err, erz, ezz = self.eta_rz(r, z)
        return np.sqrt(err**2 + 2*erz**2 + ezz**2 + 3*(er/np.maximum(r,1e-300))**2)

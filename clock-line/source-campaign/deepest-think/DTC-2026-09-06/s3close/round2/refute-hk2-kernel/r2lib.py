"""r2lib -- the refuter's own instrument for a, grad a, Hess a and the majorants Lambda_4,
Lambda_5, written from scratch (no code imported from hk2/).

Representation used (hk2 Lemma 4.2, re-derived in NOTE.md sec 2):
  a(x)      = int K(w) eta(x-w) dw ,                        K(w) = C_K w_z/|w|^5, C_K = 3/(8 pi^2)
  grad a(x) = int K(w) g(x-w) dw ,                          g = D eta (bulk part + surface parts)
  Hess a(x) = int_{|w|<d} K(w) Hess eta(x-w) dw            (near, needs eta smooth on B(x,d))
            + int_{|w|>d} gradK(w) (x) g(x-w) dw            (far)
            + oint_{|w|=d} K(w) what (x) g(x-w) dS(w)       (surface of the ball)
with gradK(w) = C_K (e_z - 5 what_z what)/|w|^5.

Coordinates centred at x = (r yhat, z):  w = s what,  what = sin(T)(cos(X) yhat + sin(X) nhat) + cos(T) zhat,
nhat in S^2 (the unit sphere of yhat^perp in R^4), integrated out:  dw = s^4 ds sin^3(T) dT 4 pi sin^2(X) dX.
Source point x' = x - w:  r' = sqrt(r^2 + s^2 sin^2 T - 2 r s sinT cosX),  z' = z - s cosT ,
rhat' = p yhat + q nhat,  p = (r - s sinT cosX)/r',  q = -s sinT sinX/r'.
Tensor components are reduced to the (yhat, zhat) block plus a transverse multiple of P_perp
(the nhat (x) nhat average is P_perp/3).  The transverse part of Hess a must equal a_r/r,
which is used as a consistency control.

Sharp surfaces in D eta (a jump [eta] across a surface S with unit normal n) contribute
[eta] n dH^4|_S to g; they are integrated in surface coordinates.

Bulk discontinuities of g (the sign flip across the equator for sgn(z) data, the sphere
rho' = 1 for a sharp radial edge) are handled by splitting the s-integration of every ray at
the crossing points, so no Gauss panel straddles a jump.
"""
import math
import numpy as np
import sympy as sp

C_K = 3.0/(8*math.pi**2)
C_GK = 4*C_K

# ---------------------------------------------------------------- polar -> (r,z) chain rule (by hand; checked by FD in r2_test)
def polar_to_rz(r, z, f, f_rho, f_phi, f_rhorho, f_rhophi, f_phiphi):
    """f(rho,phi), rho = sqrt(r^2+z^2), phi = atan2(r,z):  return f, f_r, f_z, f_rr, f_rz, f_zz.
    rho_r = s, rho_z = c, phi_r = c/rho, phi_z = -s/rho, rho_rr = c^2/rho, rho_rz = -sc/rho, rho_zz = s^2/rho,
    phi_rr = -2sc/rho^2, phi_rz = (s^2-c^2)/rho^2, phi_zz = 2sc/rho^2   (s = sin phi, c = cos phi)."""
    rho = np.sqrt(r*r + z*z); rho = np.maximum(rho, 1e-300); s = r/rho; c = z/rho
    f_r = f_rho*s + f_phi*c/rho
    f_z = f_rho*c - f_phi*s/rho
    f_rr = f_rhorho*s*s + 2*f_rhophi*s*c/rho + f_phiphi*c*c/rho**2 + f_rho*c*c/rho - 2*f_phi*s*c/rho**2
    f_rz = f_rhorho*s*c + f_rhophi*(c*c - s*s)/rho - f_phiphi*s*c/rho**2 - f_rho*s*c/rho + f_phi*(s*s - c*c)/rho**2
    f_zz = f_rhorho*c*c - 2*f_rhophi*s*c/rho + f_phiphi*s*s/rho**2 + f_rho*s*s/rho + 2*f_phi*s*c/rho**2
    return f, f_r, f_z, f_rr, f_rz, f_zz

# ---------------------------------------------------------------- radial profiles
class RadialTanh:
    """Theta(rho) = (1/2)[tanh((rho-1)/w0) - tanh((rho-R)/w1)]  (dcommon.py)."""
    def __init__(self, R, w0=0.25, w1frac=0.10):
        self.R = R; self.w0 = w0; self.w1 = w1frac*R
    def __call__(self, rho):
        a = np.clip((rho-1.0)/self.w0, -300, 300); b = np.clip((rho-self.R)/self.w1, -300, 300)
        ta = np.tanh(a); tb = np.tanh(b)
        T0 = 0.5*(ta - tb)
        sa = 1.0 - ta*ta; sb = 1.0 - tb*tb
        T1 = 0.5*(sa/self.w0 - sb/self.w1)
        T2 = 0.5*(-2*ta*sa/self.w0**2 + 2*tb*sb/self.w1**2)
        return T0, T1, T2
    support = None

class RadialSharp:
    """Theta = 1 on [1, R], 0 elsewhere (sharp edges -> surface parts handled separately)."""
    def __init__(self, R): self.R = R
    def __call__(self, rho):
        inside = (rho > 1.0) & (rho < self.R)
        z = np.zeros_like(rho)
        return np.where(inside, 1.0, 0.0), z, z

# ---------------------------------------------------------------- angular profiles A(phi), eta = -M A(phi) Theta(rho)/rho
class AngularDB:
    """A(phi) = tanh(sin phi/sin delta) tanh(cos phi/w) / sin phi   (D-B), analytic."""
    def __init__(self, delta_deg, w):
        sd = math.sin(math.radians(delta_deg))
        p = sp.symbols('phi', positive=True)
        A = sp.tanh(sp.sin(p)/sd)*sp.tanh(sp.cos(p)/w)/sp.sin(p)
        self.f = [sp.lambdify(p, sp.diff(A, p, k), 'numpy') for k in range(3)]
        self.w = w; self.sd = sd
    def __call__(self, phi):
        phi = np.clip(phi, 1e-9, math.pi-1e-9)
        return self.f[0](phi), self.f[1](phi), self.f[2](phi)
    jump0 = 0.0        # no equatorial jump

class AngularSharp:
    """A(phi) = sgn(cos phi) h(phi)/sin phi, h = min(1, phi_ax/delta)  (D-A angular profile)."""
    def __init__(self, delta_deg):
        self.delta = math.radians(delta_deg)
    def __call__(self, phi):
        phi = np.clip(phi, 1e-12, math.pi-1e-12)
        pax = np.minimum(phi, math.pi-phi)
        s = np.sin(phi); c = np.cos(phi)
        sg = np.where(phi < math.pi/2, 1.0, -1.0)
        tap = pax < self.delta
        h  = np.where(tap, pax/self.delta, 1.0)
        hp = np.where(tap, sg/self.delta, 0.0)       # dh/dphi (phi_ax = phi for phi<pi/2, pi-phi else)
        A0 = h/s
        A1 = (hp*s - h*c)/s**2
        # do A2 cleanly:  A = h/s ; A' = h'/s - h c/s^2 ; A'' = h''/s - 2 h' c/s^2 - h (1/s) + 2 h c^2/s^3
        #   (using (c/s^2)' = -1/s - 2c^2/s^3 ... let us just differentiate: d/dphi[-h c/s^2] = -h' c/s^2 + h/s + 2 h c^2/s^3)
        A2 = -hp*c/s**2 - hp*c/s**2 + h/s + 2*h*c**2/s**3
        return sg*A0, sg*A1, sg*A2
    @property
    def jump0(self):     # [A] across phi = pi/2 from z<0 to z>0 : A(pi/2^-) - A(pi/2^+) = 1 - (-1) = 2
        return 2.0

class Datum:
    """eta(x) = -M A(phi) Theta(rho)/rho.  Surfaces: equatorial plane if A has a jump; spheres if Theta sharp."""
    def __init__(self, ang, rad, M=1.0, lam=1.0):
        self.ang = ang; self.rad = rad; self.M = M; self.lam = lam
    # --- pull-back for the strained field eta_lam(r,z) = eta_0(r/lam, lam^2 z)
    def _pull(self, r, z): return r/self.lam, self.lam**2*z
    def bulk(self, r, z):
        lam = self.lam
        r0, z0 = self._pull(r, z)
        rho = np.sqrt(r0*r0 + z0*z0); rho = np.maximum(rho, 1e-300); phi = np.arctan2(r0, z0)
        T0, T1, T2 = self.rad(rho); A0, A1, A2 = self.ang(phi)
        G0 = -self.M*T0/rho; G1 = -self.M*(T1/rho - T0/rho**2); G2 = -self.M*(T2/rho - 2*T1/rho**2 + 2*T0/rho**3)
        f = G0*A0; f_rho = G1*A0; f_phi = G0*A1; f_rr_ = G2*A0; f_rp = G1*A1; f_pp = G0*A2
        e, er, ez, err, erz, ezz = polar_to_rz(r0, z0, f, f_rho, f_phi, f_rr_, f_rp, f_pp)
        return e, er/lam, ez*lam**2, err/lam**2, erz*lam, ezz*lam**4
    def gradnorm(self, r, z):
        _, er, ez, *_ = self.bulk(r, z); return np.sqrt(er*er + ez*ez)
    def hessF(self, r, z):
        _, er, ez, err, erz, ezz = self.bulk(r, z)
        return np.sqrt(err**2 + 2*erz**2 + ezz**2 + 3*(er/np.maximum(r, 1e-300))**2)
    # --- surfaces
    def plane_jump(self, rp):
        """[eta] across {z=0} in the +z direction, at radius r' (in the strained frame)."""
        if self.ang.jump0 == 0.0: return None
        r0 = rp/self.lam
        T0, _, _ = self.rad(r0)
        # eta(0+) - eta(0-) = -M Theta/rho [A(pi/2^-) - A(pi/2^+)] = -M Theta jump0 / r'
        return -self.M*T0*self.ang.jump0/np.maximum(r0, 1e-300)
    def sphere_jumps(self):
        """list of (radius, function phi' -> [eta] = eta(rad+) - eta(rad-)) for sharp radial edges (lam = 1 only)."""
        if not isinstance(self.rad, RadialSharp): return []
        assert self.lam == 1.0
        A = self.ang
        inner = (1.0, lambda phi: -self.M*A(phi)[0]/1.0)            # eta(1+) - 0
        outer = (self.rad.R, lambda phi: +self.M*A(phi)[0]/self.rad.R)   # 0 - eta(R-)
        return [inner, outer]
    def has_plane_flip(self):
        return self.ang.jump0 != 0.0

# ---------------------------------------------------------------- quadrature helpers
def gl(n):
    return np.polynomial.legendre.leggauss(n)
def panel_nodes(lo, hi, n):
    """GL nodes on [lo,hi] where lo, hi are arrays (per-direction panels). returns (nodes, weights) shape (n, *lo.shape)."""
    x, w = gl(n)
    x = x.reshape((n,) + (1,)*np.ndim(lo)); w = w.reshape((n,) + (1,)*np.ndim(lo))
    return 0.5*(hi-lo)*x + 0.5*(hi+lo), 0.5*(hi-lo)*w
def graded(a, b, c, nmin_width, ratio=2.0):
    """panel edges on [a,b] refined geometrically toward c in [a,b] down to width nmin_width."""
    edges = [a, b]
    # left side
    lo, hi = a, c
    while hi - lo > nmin_width and c > a:
        mid = c - (c - lo)/ratio if (c - lo)/ratio > nmin_width else lo
        if mid <= lo: break
        edges.append(mid); lo = mid
    lo, hi = c, b
    while hi - lo > nmin_width and c < b:
        mid = c + (hi - c)/ratio if (hi - c)/ratio > nmin_width else hi
        if mid >= hi: break
        edges.append(mid); hi = mid
    if a < c < b: edges.append(c)
    return np.array(sorted(set(edges)))

# ---------------------------------------------------------------- the instrument
class Instrument:
    def __init__(self, D, r, z, nT=160, nX=120, smax=None, extra_bps=(), npan_log=10):
        self.D = D; self.r = r; self.z = z; self.extra_bps = list(extra_bps); self.npan_log = npan_log
        self.nT = nT; self.nX = nX
        self.smax = smax if smax is not None else 4.0*D.rad.R
        # angular nodes; split T at the equator direction only through the s-breakpoints, so plain GL here
        xT, wT = gl(nT); self.T = 0.5*math.pi*(xT+1); self.wT = 0.5*math.pi*wT
        xX, wX = gl(nX); self.X = 0.5*math.pi*(xX+1); self.wX = 0.5*math.pi*wX
        self.WA = (self.wT*np.sin(self.T)**3)[:, None]*(self.wX*4*math.pi*np.sin(self.X)**2)[None, :]
        self.sT = np.sin(self.T)[:, None]*np.ones((1, nX)); self.cT = np.cos(self.T)[:, None]*np.ones((1, nX))
        self.cX = np.cos(self.X)[None, :]*np.ones((nT, 1)); self.sX = np.sin(self.X)[None, :]*np.ones((nT, 1))

    # source point and frame quantities for a block of s (shape (ns, nT, nX))
    def _src(self, s):
        r, z = self.r, self.z
        sT = self.sT[None]; cT = self.cT[None]; cX = self.cX[None]; sX = self.sX[None]
        s = s[..., None, None] if s.ndim == 1 else s
        rp = np.sqrt(np.maximum(r*r + (s*sT)**2 - 2*r*s*sT*cX, 1e-300))
        zp = z - s*cT
        p = (r - s*sT*cX)/rp; q = -s*sT*sX/rp
        return rp, zp, p, q

    def _breakpoints(self):
        """per-direction s-values where the bulk g jumps: equator crossing, sphere crossings."""
        r, z = self.r, self.z
        sT, cT, cX = self.sT, self.cT, self.cX
        bps = []
        if self.D.has_plane_flip():
            with np.errstate(divide='ignore', invalid='ignore'):
                s_eq = z/cT
            bps.append(np.where((s_eq > 0) & np.isfinite(s_eq), s_eq, np.nan))
        for rad, _ in self.D.sphere_jumps():
            xw = r*sT*cX + z*cT; rho2 = r*r + z*z
            disc = xw*xw - rho2 + rad*rad
            with np.errstate(invalid='ignore'):
                sq = np.sqrt(disc)
            s1 = np.where(disc >= 0, xw - sq, np.nan); s2 = np.where(disc >= 0, xw + sq, np.nan)
            bps.append(np.where(s1 > 0, s1, np.nan)); bps.append(np.where(s2 > 0, s2, np.nan))
        for v in self.extra_bps:
            bps.append(np.full((self.nT, self.nX), float(v)))
        return bps

    def _ray_panels(self, d, npan_log=None, n_per=16, n_far=12):
        if npan_log is None: npan_log = self.npan_log
        """build s-nodes (shape (N, nT, nX)) and weights for s in [d, smax] with per-direction breakpoints."""
        bps = self._breakpoints()
        edges = [np.full((self.nT, self.nX), d)]
        for b in bps:
            edges.append(np.where(np.isnan(b) | (b <= d), d, np.minimum(b, self.smax)))
        E = np.sort(np.stack(edges, 0), axis=0)                       # (nb+1, nT, nX)
        nodes = []; wts = []
        # short panels between consecutive breakpoints, geometric sub-panels toward each edge (kernel is peaked at small s)
        for k in range(E.shape[0]-1):
            lo, hi = E[k], E[k+1]
            for j in range(4):    # 4 geometric sub-panels
                a = lo + (hi-lo)*(0.0 if j == 0 else 2.0**(j-4)); b = lo + (hi-lo)*(2.0**(j-3))
                if j == 3: b = hi
                x, w = panel_nodes(a, b, n_per); nodes.append(x); wts.append(w)
        # the far panel [E_last, smax] in log(s)
        lo = E[-1]; hi = np.full_like(lo, self.smax)
        ulo = np.log(lo); uhi = np.log(hi)
        for j in range(npan_log):
            a = ulo + (uhi-ulo)*j/npan_log; b = ulo + (uhi-ulo)*(j+1)/npan_log
            u, wu = panel_nodes(a, b, n_far)
            nodes.append(np.exp(u)); wts.append(wu*np.exp(u))
        return np.concatenate(nodes, 0), np.concatenate(wts, 0)

    def bulk_far_and_grad(self, d, chunk=24):
        """far Hessian term (|w|>d), the s>d part of grad a, of a, and the majorants Lambda_p^{bulk}(d)."""
        S, W = self._ray_panels(d)
        H = np.zeros((2, 2)); mu = 0.0; ga = np.zeros(2); a = 0.0; L4 = L5 = 0.0; L5K = 0.0
        WA = self.WA[None]
        for i0 in range(0, S.shape[0], chunk):
            s = S[i0:i0+chunk]; w = W[i0:i0+chunk]
            rp, zp, p, q = self._src(s)
            e, er, ez, *_ = self.D.bulk(rp, zp)
            sT = self.sT[None]; cT = self.cT[None]; cX = self.cX[None]; sX = self.sX[None]
            gy = er*p; gz = ez; gn = er*q
            Ky = -5*cT*sT*cX; Kz = 1 - 5*cT*cT; Kn = -5*cT*sT*sX      # gradK / (C_K/s^5)
            wt = WA*w/s                                                # s^4 ds / s^5
            H[0, 0] += np.sum(wt*Ky*gy); H[0, 1] += np.sum(wt*Ky*gz)
            H[1, 0] += np.sum(wt*Kz*gy); H[1, 1] += np.sum(wt*Kz*gz)
            mu += np.sum(wt*Kn*gn)/3.0
            wa = WA*w                                                  # s^4 ds / s^4 (for K)
            ga[0] += np.sum(wa*cT*gy); ga[1] += np.sum(wa*cT*gz)
            a += np.sum(wa*cT*e)
            gn_abs = np.sqrt(er*er + ez*ez)
            L4 += np.sum(wa*gn_abs); L5 += np.sum(wt*gn_abs)
            L5K += np.sum(wt*np.sqrt(1 + 15*cT*cT)*gn_abs)
        return dict(H=C_K*H, mu=C_K*mu, ga=C_K*ga, a=C_K*a, L4=L4, L5=L5, L5_gradK=C_K*L5K)

    def bulk_near(self, d, ns=48):
        """near Hessian term (|w|<d), and the s<d parts of grad a and a."""
        x, w = gl(ns); s = 0.5*d*(x+1); ws = 0.5*d*w
        rp, zp, p, q = self._src(s)
        e, er, ez, err, erz, ezz = self.D.bulk(rp, zp)
        WA = self.WA[None]; cT = self.cT[None]
        wt = WA*ws[:, None, None]
        Hyy = err*p*p + (er/rp)*(1-p*p); Hyz = erz*p; Hzz = ezz
        Hnn = (err*q*q + (er/rp)*(3-q*q))/3.0        # transverse: (1/3)[err q^2 + (er/r')(3 - q^2)]
        H = np.array([[np.sum(wt*cT*Hyy), np.sum(wt*cT*Hyz)], [np.sum(wt*cT*Hyz), np.sum(wt*cT*Hzz)]])
        mu = np.sum(wt*cT*Hnn)
        ga = np.array([np.sum(wt*cT*er*p), np.sum(wt*cT*ez)])
        a = np.sum(wt*cT*e)
        return dict(H=C_K*H, mu=C_K*mu, ga=C_K*ga, a=C_K*a)

    def ball_surface(self, d):
        rp, zp, p, q = self._src(np.array([d]))
        e, er, ez, *_ = self.D.bulk(rp, zp)
        WA = self.WA[None]; sT = self.sT[None]; cT = self.cT[None]; cX = self.cX[None]; sX = self.sX[None]
        gy = er*p; gz = ez; gn = er*q
        H = np.array([[np.sum(WA*cT*sT*cX*gy), np.sum(WA*cT*sT*cX*gz)], [np.sum(WA*cT*cT*gy), np.sum(WA*cT*cT*gz)]])
        mu = np.sum(WA*cT*sT*sX*gn)/3.0
        return dict(H=C_K*H, mu=C_K*mu)

    # ---- sharp surfaces ------------------------------------------------------------------
    def plane_terms(self, nr=24, nX=None):
        """equatorial plane {z'=0}: far Hessian part, grad a part, a part (none: eta has no delta), majorants."""
        D = self.D; r, z = self.r, self.z
        if not D.has_plane_flip(): return None
        R = D.rad.R*D.lam
        # r' from 0 to 3R, graded toward r' = r (closest point), log-spaced beyond 3
        edges = np.concatenate([graded(0.0, min(3.0, 3*R), r, max(abs(z), 1e-3)/4), np.exp(np.linspace(math.log(3.0), math.log(3*R), 40))[1:]]) if R > 3 else graded(0.0, 3*R, r, abs(z)/4)
        edges = np.unique(edges)
        xX, wX = gl(self.nX*2); X = 0.5*math.pi*(xX+1); wX = 0.5*math.pi*wX
        Hyz = Hzz = 0.0; az = 0.0; L4 = L5 = 0.0
        for k in range(len(edges)-1):
            rp, wr = panel_nodes(np.array(edges[k]), np.array(edges[k+1]), nr)
            rp = rp[:, None]; wr = wr[:, None]
            J = D.plane_jump(rp)
            w2 = r*r + rp*rp - 2*r*rp*np.cos(X)[None, :] + z*z
            meas = wr*4*math.pi*rp**3*(wX*np.sin(X)**2)[None, :]
            wy = r - rp*np.cos(X)[None, :]
            Ky = -5*z*wy/w2**3.5; Kz = (1 - 5*z*z/w2)/w2**2.5
            Hyz += np.sum(meas*Ky*J); Hzz += np.sum(meas*Kz*J)
            az += np.sum(meas*(z/w2**2.5)*J)
            L4 += np.sum(meas*np.abs(J)/w2**2); L5 += np.sum(meas*np.abs(J)/w2**2.5)
        H = C_K*np.array([[0.0, Hyz], [0.0, Hzz]])
        return dict(H=H, ga=C_K*np.array([0.0, az]), L4=L4, L5=L5)

    def sphere_terms(self, nph=20, nX=None):
        D = self.D; r, z = self.r, self.z
        out = dict(H=np.zeros((2, 2)), mu=0.0, ga=np.zeros(2), L4=0.0, L5=0.0)
        phi_x = math.atan2(r, z); rho_x = math.hypot(r, z)
        for rad, jump in D.sphere_jumps():
            dist = abs(rho_x - rad)
            edges = graded(0.0, math.pi, phi_x, max(dist, 1e-4)/(4*max(rad, 1e-9)))
            # also split at the equator and the taper corners
            edges = np.unique(np.concatenate([edges, [math.pi/2, D.ang.delta if hasattr(D.ang, 'delta') else math.pi/2,
                                                      math.pi - (D.ang.delta if hasattr(D.ang, 'delta') else math.pi/2)]]))
            xe = graded(0.0, math.pi, 0.0, max(dist, 1e-4)/(4*max(rad, 1e-9)))
            for k in range(len(edges)-1):
                ph, wph = panel_nodes(np.array(edges[k]), np.array(edges[k+1]), nph)
                for j in range(len(xe)-1):
                    X, wX = panel_nodes(np.array(xe[j]), np.array(xe[j+1]), 16)
                    ph2 = ph[:, None]; wph2 = wph[:, None]; X2 = X[None, :]; wX2 = wX[None, :]
                    J = jump(ph2)*np.ones_like(X2)
                    sp_, cp_ = np.sin(ph2), np.cos(ph2)
                    wy = r - rad*sp_*np.cos(X2); wn = -rad*sp_*np.sin(X2); wz = z - rad*cp_
                    w2 = wy*wy + wn*wn + wz*wz
                    meas = wph2*wX2*rad**4*sp_**3*4*math.pi*np.sin(X2)**2
                    Ky = -5*wz*wy/w2**3.5; Kz = (1 - 5*wz*wz/w2)/w2**2.5; Kn = -5*wz*wn/w2**3.5
                    ny = sp_*np.cos(X2); nz = cp_; nn = sp_*np.sin(X2)
                    out['H'] += C_K*np.array([[np.sum(meas*Ky*ny*J), np.sum(meas*Ky*nz*J)],
                                              [np.sum(meas*Kz*ny*J), np.sum(meas*Kz*nz*J)]])
                    out['mu'] += C_K*np.sum(meas*Kn*nn*J)/3.0
                    Kw = wz/w2**2.5
                    out['ga'] += C_K*np.array([np.sum(meas*Kw*ny*J), np.sum(meas*Kw*nz*J)])
                    out['L4'] += np.sum(meas*np.abs(J)/w2**2); out['L5'] += np.sum(meas*np.abs(J)/w2**2.5)
        return out

    # ---- assemble ---------------------------------------------------------------------------
    def measure(self, d):
        F = self.bulk_far_and_grad(d); N = self.bulk_near(d); Sb = self.ball_surface(d)
        H = F['H'] + N['H'] + Sb['H']; mu = F['mu'] + N['mu'] + Sb['mu']
        ga = F['ga'] + N['ga']; a = F['a'] + N['a']
        L4 = F['L4']; L5 = F['L5']; L5K = F['L5_gradK']
        P = self.plane_terms()
        if P is not None:
            H = H + P['H']; ga = ga + P['ga']; L4 += P['L4']; L5 += P['L5']; L5K += C_GK*P['L5']
        Sp = self.sphere_terms()
        H = H + Sp['H']; mu += Sp['mu']; ga = ga + Sp['ga']; L4 += Sp['L4']; L5 += Sp['L5']; L5K += C_GK*Sp['L5']
        Hs = 0.5*(H + H.T)
        asym = abs(H[0, 1] - H[1, 0])/max(abs(Hs[0, 1]), 1e-12)
        ev = np.linalg.eigvalsh(Hs)
        mu_check = ga[0]/self.r
        hess_op = max(np.max(np.abs(ev)), abs(mu_check))
        return dict(a=a, a_r=ga[0], a_z=ga[1], grad_a=float(np.hypot(*ga)),
                    a_rr=Hs[0, 0], a_rz=Hs[0, 1], a_zz=Hs[1, 1], hess_asym=asym,
                    mu_integrated=mu, mu_from_a_r=mu_check, hess_op=hess_op,
                    L4=L4, L5=L5, L5_gradK=L5K, d=d)

# ---------------------------------------------------------------- the tensor of Lemma 3.2 and Cor 3.3
def K2_from_fields(r, a_r, a_z, a_rr, a_rz, a_zz, eta, eta_r, eta_z, n=400, seed=0):
    """sup_{|e|=1} ||d_e grad_5 b||_op from the exact second-derivative tensor; e in span(e_r, e_perp, e_z)."""
    Ha = np.zeros((5, 5)); Ha[0, 0] = a_rr; Ha[0, 4] = Ha[4, 0] = a_rz; Ha[4, 4] = a_zz
    for i in (1, 2, 3): Ha[i, i] = a_r/r
    ga = np.zeros(5); ga[0] = a_r; ga[4] = a_z
    q = a_z - eta; qr = a_rz - eta_r; qz = a_zz - eta_z
    Huz = np.zeros((5, 5)); Huz[0, 0] = q + r*qr
    for i in (1, 2, 3): Huz[i, i] = q
    Huz[0, 4] = Huz[4, 0] = r*qz; Huz[4, 4] = -2*a_z - r*a_rz
    y = np.zeros(5); y[0] = r
    def val(e):
        N = np.zeros((5, 5)); He = Ha@e
        for i in range(4):
            N[i, :] = He*y[i] + ga*e[i]; N[i, i] += ga@e
        N[4, :] = Huz@e
        return np.linalg.norm(N, 2)
    best = 0.0
    th = np.linspace(0, math.pi, n); ph = np.linspace(0, math.pi, n//2)
    for t in th:
        for p in ph:
            e = np.array([math.sin(t)*math.cos(p), math.sin(t)*math.sin(p), 0, 0, math.cos(t)])
            v = val(e)
            if v > best: best = v
    return best

def cor33_bound(r, grad_a, hess_a, eta, eta_r, eta_z):
    """hk2 Corollary 3.3, re-implemented."""
    q = grad_a + abs(eta); qr = hess_a + abs(eta_r); qz = hess_a + abs(eta_z)
    S55 = 2*grad_a + r*hess_a
    Huz = max(q, q + r*qr + r*qz, r*qz + S55)
    return r*hess_a + 2*grad_a + Huz, Huz

def ball_sups(D, r0, z0, d, n=301):
    g = np.linspace(-d, d, n)
    RR, ZZ = np.meshgrid(r0+g, z0+g, indexing='ij')
    m = (RR-r0)**2 + (ZZ-z0)**2 <= d*d
    RR = np.where(m, RR, r0); ZZ = np.where(m, ZZ, z0); RR = np.maximum(RR, 1e-6)
    gn = D.gradnorm(RR, ZZ); hf = D.hessF(RR, ZZ)
    th = np.linspace(0, 2*math.pi, 4001)
    gb = D.gradnorm(np.maximum(r0+d*np.cos(th), 1e-6), z0+d*np.sin(th))
    return float(np.max(gn)), float(np.max(hf)), float(np.max(gb))

def hk2_bound(r, d, sup_grad_B, sup_hessF_B, sup_grad_dB, L4, L5, eta, eta_r, eta_z):
    ga = d*sup_grad_B + C_K*L4
    Ha = d*sup_hessF_B + sup_grad_dB + C_GK*L5
    K2, Huz = cor33_bound(r, ga, Ha, eta, eta_r, eta_z)
    return dict(grad_a=ga, hess_a=Ha, hess_uz=Huz, K2=K2)

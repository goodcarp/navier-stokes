"""hb_lib -- shared machinery for the (H-K2) proof, seat s3close/hk2-b.

Contents
  * DatumP : the SHARP model datum  omega = -M sgn(z) h_delta(phi) 1_{rho0<rho<R}
  * DatumB : the campaign's mollified datum (lower/prove-duhamel/dcommon.py docstring)
             omega = -M tanh(sin phi/sin delta) tanh(cos phi/w) Theta(rho)
    both with a strain parameter lam:  eta_lam(r,z) = eta_0(r/lam, lam^2 z).
  * assemble_T : the EXACT 5x5x5 tensor  T[i][j][k] = d_k d_j b_i  from the scalars,
                 verified component-by-component in exact arithmetic by k1_reduction.py.
  * K2_of_T    : max_{|e|=1} || T(.,.,e) ||_op  (the Lipschitz constant of grad_5 b).
  * PolarQuad  : the (s,beta,gamma) quadrature that returns a, grad a, grad^2 a at a point.

Every derivative of the datum is produced by sympy and lambdified: no finite differences.
"""
import math
import numpy as np
import sympy as sp

C_K     = 3.0/(8.0*math.pi**2)        # sup_{|w|=1} |K(w)|,        K = 3 w_z/(8 pi^2 |w|^5)
C_GRADK = 3.0/(2.0*math.pi**2)        # sup_{|w|=1} |grad K(w)|
S4      = 8.0*math.pi**2/3.0          # |S^4|
DELTA_C = 0.2                         # int_{|w|<d} grad K dw = (1/5) e_z


# ----------------------------------------------------------------- data
class _Datum:
    """base: builds lambdified eta and its first and second (r,z) derivatives."""
    def __init__(self, expr, r, z, M=1.0, rho0=1.0, R=None, name=""):
        self.name = name; self.M = M; self.rho0 = rho0; self.R = R
        d = {}
        d['e']  = expr
        d['r']  = sp.diff(expr, r)
        d['z']  = sp.diff(expr, z)
        d['rr'] = sp.diff(expr, r, 2)
        d['rz'] = sp.diff(expr, r, z)
        d['zz'] = sp.diff(expr, z, 2)
        self.f = {k: sp.lambdify((r, z), v, 'numpy') for k, v in d.items()}

    def eta(self, r, z):   return self.f['e'](r, z)
    def der(self, r, z):   return self.f['r'](r, z), self.f['z'](r, z)
    def hess(self, r, z):  return (self.f['rr'](r, z), self.f['rz'](r, z), self.f['zz'](r, z))

    def grad_norm(self, r, z):
        er, ez = self.der(r, z)
        return np.sqrt(er**2 + ez**2)

    def hess_opnorm(self, r, z):
        """|| grad^2_5 eta ||_op at a point, using the SO(4) structure:
           eigenvalues = {d_r eta / r  (x3)} U spec[[d_rr, d_rz],[d_rz, d_zz]]."""
        err, erz, ezz = self.hess(r, z)
        er, _ = self.der(r, z)
        tr = err + ezz; dt = err*ezz - erz*erz
        disc = np.sqrt(np.maximum(tr*tr/4.0 - dt, 0.0))
        l1 = np.abs(tr/2.0 + disc); l2 = np.abs(tr/2.0 - disc)
        return np.maximum(np.maximum(l1, l2), np.abs(er/np.maximum(r, 1e-300)))

    def omega(self, r, z):
        return r*self.eta(r, z)

    def omega_grad(self, r, z):
        """grad_5 (r eta) at (r,0,0,0,z):  (eta + r d_r eta, r d_z eta)."""
        er, ez = self.der(r, z)
        return self.eta(r, z) + r*er, r*ez


def make_datumB(M=1.0, rho0=1.0, R=None, L=10.0, delta_deg=7.5, w=0.12, lam=1.0):
    """campaign datum, strained by lam.  R = rho0 e^L unless given."""
    if R is None: R = rho0*math.exp(L)
    r, z = sp.symbols('r z', positive=True)
    rr = r/lam; zz = lam**2*z                      # T_lam^{-1}
    sd = math.sin(math.radians(delta_deg))
    w0 = 0.25*rho0; w1 = 0.10*R
    s = sp.sqrt(rr**2 + zz**2)
    sph = rr/s; cph = zz/s
    x = sph/sd
    tox = sp.tanh(x)/x
    Th = sp.Rational(1, 2)*(sp.tanh((s - rho0)/w0) - sp.tanh((s - R)/w1))
    eta = -(M/(s*sd))*tox*sp.tanh(cph/w)*Th
    d = _Datum(eta, r, z, M, rho0, R, name=f"B(lam={lam})")
    d.delta = math.radians(delta_deg); d.sd = sd; d.wm = w; d.lam = lam
    return d


def make_datumP(M=1.0, rho0=1.0, R=None, L=10.0, delta_deg=7.5, lam=1.0, smooth=None):
    """SHARP model datum: eta = -M sgn(z) h_delta(phi)/r on rho0<rho<R, strained by lam.
       Used only where the plateau branch is needed in closed form (the near ball)."""
    if R is None: R = rho0*math.exp(L)
    class P:
        pass
    p = P(); p.M = M; p.rho0 = rho0; p.R = R; p.lam = lam
    p.delta = math.radians(delta_deg)
    return p


# ------------------------------------------------- the exact tensor  d_k d_j b_i
def assemble_T(r, z, a_r, a_z, a_rr, a_rz, a_zz, om, om_r, om_z):
    """Exact  T[i,j,k] = d_k d_j b_i  at the point x = (r,0,0,0,z), from the reduction
       verified in k1_reduction.py.  All arguments are the (r,z)-frame scalars."""
    X = np.array([r, 0.0, 0.0, 0.0, z])
    ga = np.array([a_r, 0.0, 0.0, 0.0, a_z])                      # grad a
    H = np.zeros((5, 5))                                          # grad^2 a
    H[0, 0] = a_rr
    H[1, 1] = H[2, 2] = H[3, 3] = a_r/r
    H[0, 4] = H[4, 0] = a_rz
    H[4, 4] = a_zz
    gom = np.array([om_r, 0.0, 0.0, 0.0, om_z])                   # grad omega^theta
    T = np.zeros((5, 5, 5))
    for i in range(5):
        for j in range(5):
            for k in range(5):
                if i < 4:
                    v = X[i]*H[k, j]
                    if i == j: v += ga[k]
                    if i == k: v += ga[j]
                else:
                    if j < 4 and k < 4:
                        dd = ((1.0 if j == k else 0.0)/r - X[j]*X[k]/r**3)
                        v = (ga[4] if j == k else 0.0) + X[j]*H[4, k] - dd*om - (X[j]/r)*gom[k]
                    elif j < 4 and k == 4:
                        v = X[j]*H[4, 4] - (X[j]/r)*gom[4]
                    elif j == 4 and k < 4:
                        v = -3*ga[k] - sum(X[m]*H[m, k] for m in range(4))
                    else:
                        v = -sum(X[m]*H[m, 4] for m in range(4)) - 2*ga[4]
                T[i, j, k] = v
    return T


def K2_of_T(T, ngrid=4000, seed=11):
    """max_{|e|=1} || T(.,.,e) ||_op  by random search + power refinement."""
    rng = np.random.default_rng(seed)
    E = rng.normal(size=(ngrid, 5)); E /= np.linalg.norm(E, axis=1, keepdims=True)
    E = np.vstack([E, np.eye(5), -np.eye(5)])
    M = np.einsum('ijk,nk->nij', T, E)
    s = np.linalg.norm(M, ord=2, axis=(1, 2))
    i = int(np.argmax(s)); best = float(s[i]); e = E[i].copy()
    # local refinement: gradient of sigma_max wrt e is u^T dM v
    for _ in range(400):
        Me = np.einsum('ijk,k->ij', T, e)
        U, S, Vt = np.linalg.svd(Me)
        g = np.einsum('i,ijk,j->k', U[:, 0], T, Vt[0, :])
        g -= (g @ e)*e
        n = np.linalg.norm(g)
        if n < 1e-14: break
        e2 = e + 0.05*g/n; e2 /= np.linalg.norm(e2)
        v2 = np.linalg.norm(np.einsum('ijk,k->ij', T, e2), ord=2)
        if v2 > best: best, e = v2, e2
        else: break
    return best, e


def elementary_norms(r, z, ngrid=3000):
    """nu_m = max_{|e|=1} ||T^{(m)}(.,.,e)||_op for each unit scalar slot m.
       slots: a_r, a_z, a_rr, a_rz, a_zz, a_r/r, omega/r, om_r, om_z."""
    slots = ['a_r', 'a_z', 'a_rr', 'a_rz', 'a_zz', 'a_r_over_r', 'om_over_r', 'om_r', 'om_z']
    out = {}
    for m in slots:
        kw = dict(a_r=0.0, a_z=0.0, a_rr=0.0, a_rz=0.0, a_zz=0.0, om=0.0, om_r=0.0, om_z=0.0)
        # a_r feeds BOTH grad a and the a_r/r Hessian entry; split them into two slots by
        # building the tensor with an explicit override.
        T = _assemble_slot(r, z, m)
        out[m] = K2_of_T(T, ngrid=ngrid)[0]
    return out


def _assemble_slot(r, z, slot):
    """T for a unit value of one scalar slot, all others zero (a_r and a_r/r decoupled)."""
    X = np.array([r, 0.0, 0.0, 0.0, z])
    ga = np.zeros(5); H = np.zeros((5, 5)); gom = np.zeros(5); om = 0.0
    if slot == 'a_r':          ga[0] = 1.0
    elif slot == 'a_z':        ga[4] = 1.0
    elif slot == 'a_rr':       H[0, 0] = 1.0
    elif slot == 'a_rz':       H[0, 4] = H[4, 0] = 1.0
    elif slot == 'a_zz':       H[4, 4] = 1.0
    elif slot == 'a_r_over_r': H[1, 1] = H[2, 2] = H[3, 3] = 1.0
    elif slot == 'om_over_r':  om = r          # so that om/r = 1
    elif slot == 'om_r':       gom[0] = 1.0
    elif slot == 'om_z':       gom[4] = 1.0
    else: raise ValueError(slot)
    T = np.zeros((5, 5, 5))
    for i in range(5):
        for j in range(5):
            for k in range(5):
                if i < 4:
                    v = X[i]*H[k, j]
                    if i == j: v += ga[k]
                    if i == k: v += ga[j]
                else:
                    if j < 4 and k < 4:
                        dd = ((1.0 if j == k else 0.0)/r - X[j]*X[k]/r**3)
                        v = (ga[4] if j == k else 0.0) + X[j]*H[4, k] - dd*om - (X[j]/r)*gom[k]
                    elif j < 4 and k == 4:
                        v = X[j]*H[4, 4] - (X[j]/r)*gom[4]
                    elif j == 4 and k < 4:
                        v = -3*ga[k] - sum(X[m]*H[m, k] for m in range(4))
                    else:
                        v = -sum(X[m]*H[m, 4] for m in range(4)) - 2*ga[4]
                T[i, j, k] = v
    return T


# ------------------------------------------------------------ polar quadrature
class PolarQuad:
    """a, grad a, grad^2 a at x = (r*,0,0,0,z*) for an axisymmetric source eta.

       a(x)   = int_0^inf int_{S^4} K(sig) eta(x - s sig) dOmega_4 ds
       grad a = int_0^inf int_{S^4} K(sig) grad eta(x - s sig) dOmega_4 ds
       d_k d_j a = int_{s<dbar}(1/s) int d_k K(sig)[ (d_j eta)(x-s sig) - (d_j eta)(x) ]
                 + int_{s>dbar}(1/s) int d_k K(sig) (d_j eta)(x - s sig)
                 + (1/5) delta_{k,z} (d_j eta)(x)
       (the sign of the last term is fixed by the control Lap_5 a = d_z eta; the minus
        sign fails that control by exactly 2/5, which is how the sign was settled)
       with sig = (sin b cos g, sin b sin g m, cos b),  m in S^2,
       dOmega_4 = sin^3 b sin^2 g db dg dOmega_2(m).
    """
    def __init__(self, nb=80, ng=80, ns=64, nfar=140):
        self.b, self.wb = np.polynomial.legendre.leggauss(nb)
        self.b = 0.5*math.pi*(self.b + 1.0); self.wb = 0.5*math.pi*self.wb
        self.g, self.wg = np.polynomial.legendre.leggauss(ng)
        self.g = 0.5*math.pi*(self.g + 1.0); self.wg = 0.5*math.pi*self.wg
        self.ns, self.nfar = ns, nfar
        self.vg, self.wv = np.polynomial.legendre.leggauss(ns)
        self.vg = 0.5*(self.vg + 1.0); self.wv = 0.5*self.wv
        self.ug, self.wu = np.polynomial.legendre.leggauss(nfar)

    def _angles(self):
        B, G = np.meshgrid(self.b, self.g, indexing='ij')
        WB, WG = np.meshgrid(self.wb, self.wg, indexing='ij')
        W = WB*WG*np.sin(B)**3*np.sin(G)**2
        return B, G, W

    def eval(self, dat, rs, zs, dbar, smax):
        B, G, W = self._angles()
        sb, cb, sg, cg = np.sin(B), np.cos(B), np.sin(G), np.cos(G)
        er_x, ez_x = dat.der(rs, zs)                       # grad eta at x, (r,z) frame

        FOURPI = 4.0*math.pi
        P2 = FOURPI/3.0                                    # int_{S^2} m_1^2 dOmega_2

        def geom(s):
            y1 = rs - s*sb*cg
            yq = -s*sb*sg
            rp = np.sqrt(y1*y1 + yq*yq)
            zp = zs - s*cb
            return rp, zp, y1/np.maximum(rp, 1e-300), yq/np.maximum(rp, 1e-300)

        # ---- a and grad a : s from 0 to smax, weight ds (integrand regular)
        # substitution s = smax * t^2 concentrates nodes near 0 where eta is largest? use
        # two panels: [0,dbar] and [dbar,smax] with log spacing on the second.
        a_val = 0.0; ga_r = 0.0; ga_z = 0.0
        for lo, hi, logsub in ((0.0, dbar, False), (dbar, smax, True)):
            if logsub:
                u = 0.5*(math.log(hi) - math.log(lo))*self.ug + 0.5*(math.log(hi) + math.log(lo))
                wu = 0.5*(math.log(hi) - math.log(lo))*self.wu
                svals = np.exp(u); sw = wu*svals               # ds = s du
            else:
                svals = lo + (hi - lo)*self.vg; sw = (hi - lo)*self.wv
            for s, ws in zip(svals, sw):
                rp, zp, p, q = geom(s)
                e = dat.eta(rp, zp); erp, ezp = dat.der(rp, zp)
                wgt = W*ws
                a_val += C_K*FOURPI*np.sum(wgt*cb*e)
                ga_r  += C_K*FOURPI*np.sum(wgt*cb*erp*p)
                ga_z  += C_K*FOURPI*np.sum(wgt*cb*ezp)

        # ---- grad^2 a
        d55 = 0.0; d15 = 0.0; d11 = 0.0; d22 = 0.0; d51 = 0.0
        for lo, hi, near in ((0.0, dbar, True), (dbar, smax, False)):
            if near:
                svals = lo + (hi - lo)*self.vg; sw = (hi - lo)*self.wv
            else:
                u = 0.5*(math.log(hi) - math.log(lo))*self.ug + 0.5*(math.log(hi) + math.log(lo))
                wu = 0.5*(math.log(hi) - math.log(lo))*self.wu
                svals = np.exp(u); sw = wu*svals
            for s, ws in zip(svals, sw):
                rp, zp, p, q = geom(s)
                erp, ezp = dat.der(rp, zp)
                gy1 = erp*p                                  # e_1 component of grad eta(y)
                gy5 = ezp                                    # e_5 component
                if near:
                    b1 = gy1 - er_x; b5 = gy5 - ez_x
                else:
                    b1 = gy1; b5 = gy5
                wgt = W*ws/s
                # d_5 K = C_K (1 - 5 cos^2 b) ; d_1 K = -5 C_K cos b sin b cos g
                d55 += C_K*FOURPI*np.sum(wgt*(1.0 - 5.0*cb*cb)*b5)
                d15 += C_K*FOURPI*np.sum(wgt*(-5.0*cb*sb*cg)*b5)
                d51 += C_K*FOURPI*np.sum(wgt*(1.0 - 5.0*cb*cb)*b1)
                d11 += C_K*FOURPI*np.sum(wgt*(-5.0*cb*sb*cg)*b1)
                # (k,j) = (2,2):  d_2 K = -5 C_K cos b sin b sin g m_1 ; (grad eta)_2 = erp*q*m_1
                d22 += C_K*P2*np.sum(wgt*(-5.0*cb*sb*sg)*erp*q)
        d55 += DELTA_C*ez_x
        d51 += DELTA_C*er_x
        return dict(a=a_val, a_r=ga_r, a_z=ga_z, a_zz=d55, a_rz=d15, a_rz_sym=d51,
                    a_rr=d11, a_r_over_r=d22)

    def frakE(self, dat, rs, zs, smax):
        """E(x) = int |x-y|^{-4} |grad eta(y)| dy = int_0^smax int_{S^4} |grad eta(x-s sig)| dOmega ds"""
        B, G, W = self._angles()
        sb, cb, sg, cg = np.sin(B), np.cos(B), np.sin(G), np.cos(G)
        FOURPI = 4.0*math.pi
        tot = 0.0
        for lo, hi, logs in ((0.0, 1e-2*smax, False), (1e-2*smax, smax, True)):
            if logs:
                u = 0.5*(math.log(hi)-math.log(lo))*self.ug + 0.5*(math.log(hi)+math.log(lo))
                wu = 0.5*(math.log(hi)-math.log(lo))*self.wu
                sv = np.exp(u); sw = wu*sv
            else:
                sv = lo + (hi-lo)*self.vg; sw = (hi-lo)*self.wv
            for s, ws in zip(sv, sw):
                y1 = rs - s*sb*cg; yq = -s*sb*sg
                rp = np.sqrt(y1*y1 + yq*yq); zp = zs - s*cb
                gn = dat.grad_norm(rp, zp)
                tot += FOURPI*np.sum(W*ws*gn)
        return float(tot)

    def frakF(self, dat, rs, zs, dbar, smax):
        """F(x,dbar) = int_{|x-y|>dbar} |x-y|^{-5}|grad eta(y)| dy = int_dbar^smax (1/s) int ... ds"""
        B, G, W = self._angles()
        sb, cb, sg, cg = np.sin(B), np.cos(B), np.sin(G), np.cos(G)
        FOURPI = 4.0*math.pi
        u = 0.5*(math.log(smax)-math.log(dbar))*self.ug + 0.5*(math.log(smax)+math.log(dbar))
        wu = 0.5*(math.log(smax)-math.log(dbar))*self.wu
        sv = np.exp(u)
        tot = 0.0
        for s, ws in zip(sv, wu):           # (1/s) ds = du
            y1 = rs - s*sb*cg; yq = -s*sb*sg
            rp = np.sqrt(y1*y1 + yq*yq); zp = zs - s*cb
            gn = dat.grad_norm(rp, zp)
            tot += FOURPI*np.sum(W*ws*gn)
        return float(tot)



    def eval_full(self, dat, rs, zs, dbar, smax):
        """Same decomposition as eval(), but returning the per-shell angular integrals
              J_j(s)    = int_{S^4} K(sig) (d_j eta)(x - s sig) dOmega_4
              I_kj(s)   = int_{S^4} d_k K(sig) (d_j eta)(x - s sig) dOmega_4
        so that BOTH the exact value (signed s-sum) and a bound that keeps the ANGULAR
        cancellation but drops only the s-cancellation (sum of |.|) can be formed.
        Returns a dict with 'exact' and 'absbound' entries for
        a_r, a_z, a_rr, a_rz, a_zz, a_r_over_r."""
        B, G, W = self._angles()
        sb, cb, sg, cg = np.sin(B), np.cos(B), np.sin(G), np.cos(G)
        er_x, ez_x = dat.der(rs, zs)
        er_x = float(er_x); ez_x = float(ez_x)
        FOURPI = 4.0*math.pi
        P2 = FOURPI/3.0
        w5 = 1.0 - 5.0*cb*cb                 # d_5 K / C_K
        w1 = -5.0*cb*sb*cg                   # d_1 K / C_K
        w2 = -5.0*cb*sb*sg                   # d_2 K / (C_K m_1)

        def panels():
            out = []
            for lo, hi, near in ((0.0, dbar, True), (dbar, smax, False)):
                if near:
                    sv = lo + (hi - lo)*self.vg; sw = (hi - lo)*self.wv
                else:
                    u = 0.5*(math.log(hi)-math.log(lo))*self.ug + 0.5*(math.log(hi)+math.log(lo))
                    wu = 0.5*(math.log(hi)-math.log(lo))*self.wu
                    sv = np.exp(u); sw = wu*sv
                out.append((sv, sw, near))
            return out

        ex = dict(a_r=0.0, a_z=0.0, a_rr=0.0, a_rz=0.0, a_zz=0.0, a_r_over_r=0.0)
        ab = dict(a_r=0.0, a_z=0.0, a_rr=0.0, a_rz=0.0, a_zz=0.0, a_r_over_r=0.0)
        for sv, sw, near in panels():
            for s, ws in zip(sv, sw):
                y1 = rs - s*sb*cg; yq = -s*sb*sg
                rp = np.sqrt(y1*y1 + yq*yq); zp = zs - s*cb
                erp, ezp = dat.der(rp, zp)
                p = y1/np.maximum(rp, 1e-300); q = yq/np.maximum(rp, 1e-300)
                g1 = erp*p; g5 = ezp
                # grad a  (no subtraction anywhere; weight ds)
                J1 = C_K*FOURPI*np.sum(W*cb*g1)
                J5 = C_K*FOURPI*np.sum(W*cb*g5)
                ex['a_r'] += ws*J1;  ab['a_r'] += ws*abs(J1)
                ex['a_z'] += ws*J5;  ab['a_z'] += ws*abs(J5)
                # grad^2 a : weight ds/s ; subtract grad eta(x) on the near panel only
                b1 = g1 - er_x if near else g1
                b5 = g5 - ez_x if near else g5
                I11 = C_K*FOURPI*np.sum(W*w1*b1)
                I15 = C_K*FOURPI*np.sum(W*w1*b5)
                I55 = C_K*FOURPI*np.sum(W*w5*b5)
                I22 = C_K*P2*np.sum(W*w2*erp*q)
                f = ws/s
                ex['a_rr'] += f*I11; ex['a_rz'] += f*I15
                ex['a_zz'] += f*I55; ex['a_r_over_r'] += f*I22
                if not near:
                    ab['a_rr'] += f*abs(I11); ab['a_rz'] += f*abs(I15)
                    ab['a_zz'] += f*abs(I55); ab['a_r_over_r'] += f*abs(I22)
        ex['a_zz'] += DELTA_C*ez_x
        ab['a_zz'] = ab['a_zz'] + abs(DELTA_C*ez_x)
        return dict(exact=ex, farabs=ab)


    def shell_integrals(self, dat, rs, zs, svals):
        """J_j(s) and I_kj(s) on a given s-grid (raw, no subtraction).
             J_j(s)  = int_{S^4} K(sig) (d_j eta)(x - s sig) dOmega_4
             I_kj(s) = int_{S^4} d_k K(sig) (d_j eta)(x - s sig) dOmega_4
           and the SUBTRACTED versions Is_kj(s) with (d_j eta)(x) removed (for the near panel)."""
        B, G, W = self._angles()
        sb, cb, sg, cg = np.sin(B), np.cos(B), np.sin(G), np.cos(G)
        er_x, ez_x = float(dat.der(rs, zs)[0]), float(dat.der(rs, zs)[1])
        FOURPI = 4.0*math.pi; P2 = FOURPI/3.0
        w5 = 1.0 - 5.0*cb*cb; w1 = -5.0*cb*sb*cg; w2 = -5.0*cb*sb*sg
        n = len(svals)
        out = {k: np.empty(n) for k in ('J1', 'J5', 'I11', 'I15', 'I55', 'I22',
                                        'Is11', 'Is15', 'Is55')}
        for i, s in enumerate(svals):
            y1 = rs - s*sb*cg; yq = -s*sb*sg
            rp = np.sqrt(y1*y1 + yq*yq); zp = zs - s*cb
            erp, ezp = dat.der(rp, zp)
            p = y1/np.maximum(rp, 1e-300); q = yq/np.maximum(rp, 1e-300)
            g1 = erp*p; g5 = ezp
            out['J1'][i] = C_K*FOURPI*np.sum(W*cb*g1)
            out['J5'][i] = C_K*FOURPI*np.sum(W*cb*g5)
            out['I11'][i] = C_K*FOURPI*np.sum(W*w1*g1)
            out['I15'][i] = C_K*FOURPI*np.sum(W*w1*g5)
            out['I55'][i] = C_K*FOURPI*np.sum(W*w5*g5)
            out['I22'][i] = C_K*P2*np.sum(W*w2*erp*q)
            out['Is11'][i] = C_K*FOURPI*np.sum(W*w1*(g1 - er_x))
            out['Is15'][i] = C_K*FOURPI*np.sum(W*w1*(g5 - ez_x))
            out['Is55'][i] = C_K*FOURPI*np.sum(W*w5*(g5 - ez_x))
        out['er_x'] = er_x; out['ez_x'] = ez_x
        return out

    def gfun(self, dat, rs, zs, svals):
        """g(s) = int_{S^4} |grad_5 eta(x - s sig)| dOmega_4(sig)   on a given s-grid."""
        B, G, W = self._angles()
        sb, cb, sg, cg = np.sin(B), np.cos(B), np.sin(G), np.cos(G)
        FOURPI = 4.0*math.pi
        out = np.empty(len(svals))
        for i, s in enumerate(svals):
            y1 = rs - s*sb*cg; yq = -s*sb*sg
            rp = np.sqrt(y1*y1 + yq*yq); zp = zs - s*cb
            out[i] = FOURPI*np.sum(W*dat.grad_norm(rp, zp))
        return out



def Lambda2(dat, rs, zs, dbar, n=161):
    """sup over the 5D ball B(x,dbar) of ||grad^2_5 eta||_op.
       The SO(4)-orbit of the ball projects onto the 2D disc (r'-r)^2+(z'-z)^2 <= dbar^2."""
    t = np.linspace(-1.0, 1.0, n)
    RR, ZZ = np.meshgrid(rs + dbar*t, zs + dbar*t, indexing='ij')
    msk = ((RR-rs)**2 + (ZZ-zs)**2 <= dbar*dbar) & (RR > 1e-9)
    v = dat.hess_opnorm(RR, ZZ)
    return float(np.max(np.where(msk, v, -np.inf)))

def EF_from_g(svals, g):
    """E = int_0^inf g ds ;  F(dbar) = int_dbar^inf g ds/s.
       svals is a log-uniform grid; trapezoid in log s for F, and in s for E,
       with the s<s0 head handled by g ~ g(0) (E) and by the exact 1/s tail (F)."""
    ls = np.log(svals)
    # E: int g ds = int g s dlog s
    Ecum = np.concatenate([[0.0], np.cumsum(0.5*(g[1:]*svals[1:] + g[:-1]*svals[:-1])*np.diff(ls))])
    Fcum = np.concatenate([[0.0], np.cumsum(0.5*(g[1:] + g[:-1])*np.diff(ls))])
    Etot = Ecum[-1] + g[0]*svals[0]          # head int_0^{s0} g ds ~ g(0) s0
    def F(dbar):
        return float(np.interp(np.log(dbar), ls, Fcum[-1] - Fcum))
    return Etot, F

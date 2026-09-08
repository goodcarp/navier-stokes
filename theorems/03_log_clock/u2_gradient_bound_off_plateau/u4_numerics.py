#!/usr/bin/env python3
"""
u4 -- the numerical confrontation for BLOCK 3.

Instrument (built here, shares no code with hk2/k3, k4, k9 or with L3v's g-series).
Polar coordinates centred on the EVALUATION point kill the kernel singularity
exactly, so nothing is regularised and no principal value appears:

    w = x - x',  what = (n sin al, cos al),  n = cos be e_1 + sin be e_perp in S^3,
    dw = xi^4 dxi sin^3(al) dal * 4 pi sin^2(be) dbe ,
    K(w) = 3 w_z/(8 pi^2 |w|^5) = 3 cos(al)/(8 pi^2 xi^4)
    =>  K(w) dw = (3/(2 pi)) cos(al) sin^3(al) sin^2(be) dxi dal dbe   (no xi !)

    a(x)       = (3/2pi) INT cos al sin^3 al sin^2 be * eta(x') dxi dal dbe
    d_j a(x)   = (3/2pi) INT cos al sin^3 al sin^2 be * (d_j eta)(x') dxi dal dbe

with x' = (r', z'),  r' = sqrt(r^2 - 2 r xi sin al cos be + xi^2 sin^2 al),
z' = z - xi cos al, and (d_{y1} eta)(x') = ((r - xi sin al cos be)/r') d_{r'} eta.

Fields:
  (E) the exact strained plateau            eta_lam = eta_0 o T_lam^{-1}
  (P) the perturbed-map field               eta     = eta_0 o Phi^{-1},
      Phi = Lam o (I + psi),  Lam = T_lam,  psi axisymmetric, z-parity preserving.

Measured: Gamma = ||grad_5 b||_op = max(|a|, sigma_max([[a+Q,P],[P-om,-2a-Q]])),
Q = r d_r a, P = r d_z a, om = r eta.  Compared with 2 a(0) + C'' M.
"""
import json, math, sys
import numpy as np

# ------------------------------------------------------------------ the datum
DELTA = math.radians(7.5); SD = math.sin(DELTA)
W_EQ = 0.20
EPS_R = 0.25
M = 1.0; RHO0 = 1.0

def theta_ramp(rho, L):
    R = math.exp(L)
    a1 = (np.log(rho) - EPS_R)/EPS_R
    a2 = (np.log(rho/R) + EPS_R)/EPS_R
    return 0.5*(np.tanh(a1) - np.tanh(a2))

def theta_ramp_p(rho, L):
    """d Theta/d rho"""
    R = math.exp(L)
    a1 = (np.log(rho) - EPS_R)/EPS_R
    a2 = (np.log(rho/R) + EPS_R)/EPS_R
    return 0.5*(1.0/np.cosh(a1)**2 - 1.0/np.cosh(a2)**2)/(EPS_R*rho)

def eta0(rp, zp):
    """eta_0(r,z) = -(M/rho) h(phi) Theta(rho),  h = tanh(s/sd) tanh(c/w)/s."""
    rho = np.hypot(rp, zp)
    rho = np.where(rho < 1e-300, 1e-300, rho)
    s = rp/rho; c = zp/rho
    s = np.where(s < 1e-300, 1e-300, s)
    h = np.tanh(s/SD)*np.tanh(c/W_EQ)/s
    return -(M/rho)*h*theta_ramp(rho, eta0.L)

def eta0_grad(rp, zp):
    """(d_r eta_0, d_z eta_0) by exact differentiation."""
    rho = np.hypot(rp, zp); rho = np.where(rho < 1e-300, 1e-300, rho)
    s = rp/rho; c = zp/rho
    s = np.where(s < 1e-300, 1e-300, s)
    T = np.tanh(s/SD); Z = np.tanh(c/W_EQ)
    h = T*Z/s
    # dh/ds and dh/dc at fixed (s,c) on the unit circle -> use dh/dphi
    Tp = (1.0/np.cosh(s/SD)**2)/SD          # dT/ds
    Zp = (1.0/np.cosh(c/W_EQ)**2)/W_EQ      # dZ/dc
    # phi derivative:  ds/dphi = c, dc/dphi = -s
    dh_dphi = (Tp*c*Z + T*Zp*(-s))/s - T*Z*c/s**2
    Th = theta_ramp(rho, eta0.L); Thp = theta_ramp_p(rho, eta0.L)
    # eta = -(1/rho) h(phi) Theta(rho)
    d_rho = -(Thp/rho - Th/rho**2)*h
    d_phi = -(Th/rho)*dh_dphi
    # (r,z) <- (rho,phi):  d_r = s d_rho + (c/rho) d_phi ; d_z = c d_rho - (s/rho) d_phi
    return s*d_rho + (c/rho)*d_phi, c*d_rho - (s/rho)*d_phi
eta0.L = 10.0

# ------------------------------------------------------- the map perturbation
class Pert:
    """psi(alpha) = (psi_r(r,z) yhat, psi_z(r,z)); psi_r even in z, psi_z odd in z."""
    def __init__(self, eps, L):
        self.eps = eps; self.L = L
        self.k = 2*math.pi
    def _uc(self, r, z):
        rho = np.hypot(r, z); rho = np.where(rho < 1e-300, 1e-300, rho)
        return rho, r/rho, z/rho, np.log(rho)/self.L
    def val(self, r, z):
        rho, s, c, u = self._uc(r, z)
        pr = self.eps*rho*(1 - c*c)*np.sin(self.k*u)
        pz = self.eps*rho*c*(1 - c*c)*np.cos(self.k*u)
        return pr, pz
    def jac(self, r, z):
        """analytic d(psi_r,psi_z)/d(r,z)."""
        rho, sn, c, u = self._uc(r, z)
        k = self.k; e = self.eps
        A = 1 - c*c; Ap = -2*c
        B = c*(1 - c*c); Bp = 1 - 3*c*c
        S = np.sin(k*u); C = np.cos(k*u)
        j11 = e*(sn*A*S - sn*c*Ap*S + sn*A*k*C/self.L)
        j12 = e*(c*A*S + sn*sn*Ap*S + c*A*k*C/self.L)
        j21 = e*(sn*B*C - sn*c*Bp*C - sn*B*k*S/self.L)
        j22 = e*(c*B*C + sn*sn*Bp*C - c*B*k*S/self.L)
        return j11, j12, j21, j22

def pert_norms(pert, L, n=600):
    """sup |psi|/|alpha| and sup ||Dpsi||_op over the shell (rho in [1, e^L])."""
    rho = np.exp(np.linspace(0, L, n)); ph = np.linspace(1e-6, math.pi-1e-6, n)
    RR, PP = np.meshgrid(rho, ph, indexing='ij')
    r = RR*np.sin(PP); z = RR*np.cos(PP)
    pr, pz = pert.val(r, z)
    rel = np.max(np.hypot(pr, pz)/RR)
    j11, j12, j21, j22 = pert.jac(r, z)
    # in R^5 the y-block of Dpsi has eigenvalue psi_r/r on the 3 transverse directions
    tr3 = np.max(np.abs(pr/np.maximum(r, 1e-12)))
    Jm = np.stack([np.stack([j11, j12], -1), np.stack([j21, j22], -1)], -2)
    sv = np.linalg.svd(Jm.reshape(-1, 2, 2), compute_uv=False)[:, 0]
    return float(rel), float(max(sv.max(), tr3))

# ------------------------------------------------------------ the field maps
def field_exact(lam):
    def f(rp, zp, want_grad):
        a_r = rp/lam; a_z = zp*lam*lam
        if not want_grad:
            return eta0(a_r, a_z), None, None
        gr, gz = eta0_grad(a_r, a_z)
        return eta0(a_r, a_z), gr/lam, gz*lam*lam
    return f

def field_perturbed(lam, pert, nit=14):
    def f(rp, zp, want_grad):
        b_r = rp/lam; b_z = zp*lam*lam                 # beta = T_lam^{-1} x
        ar = b_r.copy(); az = b_z.copy()
        for _ in range(nit):                            # alpha = beta - psi(alpha)
            pr, pz = pert.val(ar, az)
            nr = b_r - pr; nz = b_z - pz
            ar, az = nr, nz
        val = eta0(ar, az)
        if not want_grad:
            return val, None, None
        gr, gz = eta0_grad(ar, az)
        j11, j12, j21, j22 = pert.jac(ar, az)
        # DPhi = T_lam (I + Dpsi); grad eta = (DPhi)^{-T} grad eta_0
        # (r,z)-block: (I+Dpsi)^{-T} then T_lam^{-1} = diag(1/lam, lam^2)
        A11 = 1+j11; A12 = j12; A21 = j21; A22 = 1+j22
        det = A11*A22 - A12*A21
        # (I+Dpsi)^{-T} = [[A22, -A21],[-A12, A11]]/det
        vr = (A22*gr - A21*gz)/det
        vz = (-A12*gr + A11*gz)/det
        return val, vr/lam, vz*lam*lam
    return f

# ------------------------------------------------------------ the quadrature
def gl_panels(edges, n):
    xs, ws = np.polynomial.legendre.leggauss(n)
    X = []; W = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        X.append(0.5*(hi-lo)*xs + 0.5*(hi+lo)); W.append(0.5*(hi-lo)*ws)
    return np.concatenate(X), np.concatenate(W)

def make_nodes(L, nal=200, nbe=32, nxi=12, dlog=0.75):
    al, wal = gl_panels(np.linspace(0, math.pi, 9), nal//8)
    be, wbe = gl_panels(np.linspace(0, math.pi, 5), nbe//4)
    R = math.exp(L)
    lo = -8.0; hi = math.log(4*R)
    edges = [0.0] + list(np.exp(np.arange(lo, hi + dlog, dlog)))
    xi, wxi = gl_panels(np.array(edges), nxi)
    return (al, wal, be, wbe, xi, wxi)

def evaluate(field, x, nodes, chunk=400000):
    """return a, d_r a, d_z a at x = (r, z)."""
    r, z = x
    al, wal, be, wbe, xi, wxi = nodes
    ca = np.cos(al); sa = np.sin(al)
    cb = np.cos(be); sb = np.sin(be)
    Wang = (ca*sa**3*wal)[:, None]*(sb**2*wbe)[None, :]      # (nal, nbe)
    Wang = Wang.ravel()
    SA = np.repeat(sa, len(be)); CA = np.repeat(ca, len(be))
    CB = np.tile(cb, len(al))
    acc = np.zeros(3)
    nA = len(Wang)
    step = max(1, chunk//max(nA, 1))
    for i0 in range(0, len(xi), step):
        X = xi[i0:i0+step]; WX = wxi[i0:i0+step]
        dy = X[:, None]*SA[None, :]
        y1 = r - dy*CB[None, :]
        rp = np.sqrt(np.maximum(r*r - 2*r*dy*CB[None, :] + dy*dy, 0.0))
        zp = z - X[:, None]*CA[None, :]
        val, gr, gz = field(rp, zp, True)
        cosdir = np.where(rp > 1e-300, y1/np.maximum(rp, 1e-300), 0.0)
        w = WX[:, None]*Wang[None, :]
        acc[0] += float(np.sum(w*val))
        acc[1] += float(np.sum(w*(cosdir*gr)))
        acc[2] += float(np.sum(w*gz))
    return acc*(3.0/(2.0*math.pi))

def gamma_op(a, dra, dza, om, r):
    Q = r*dra; P = r*dza
    B = np.array([[a+Q, P], [P-om, -2*a-Q]])
    sv = np.linalg.svd(B, compute_uv=False)
    return float(max(abs(a), sv[0]))

# ------------------------------------------------------------------- driver
def run():
    OUT = {}
    from scipy import integrate as _it
    # control: a(0) = kappa_prof * INT Theta dlog rho  exactly (l=1 interior mode)
    kap, _ = _it.quad(lambda p: 0.75*np.tanh(np.sin(p)/SD)*np.tanh(np.cos(p)/W_EQ)
                      * np.cos(p)*np.sin(p)**2, 0, math.pi, limit=300)
    OUT['kappa_profile'] = kap
    for L in (10.0, 40.0):
        eta0.L = L
        nodes = make_nodes(L)
        nodes_fine = make_nodes(L, nal=400, nbe=64, nxi=16, dlog=0.5)
        R = math.exp(L)
        # slab sample: 2 rho0 <= rho <= R/2, sin phi >= 1/2
        rhos = np.exp(np.linspace(math.log(2.0), math.log(R/2), 4))
        phis = [math.radians(d) for d in (30, 60, 90)]
        pts = [(rho*math.sin(p), rho*math.cos(p)) for rho in rhos for p in phis]
        for lam in (1.0, 1.5):
            cases = {'exact': field_exact(lam)}
            for eps_t in (0.02, 0.05):
                pert = Pert(1.0, L)
                rel, dn = pert_norms(pert, L)
                pert.eps = eps_t/max(rel, dn)
                rel2, dn2 = pert_norms(pert, L)
                cases[f'pert{eps_t}'] = field_perturbed(lam, pert)
                OUT[f'pertnorm_L{int(L)}_eps{eps_t}'] = dict(sup_rel=rel2, sup_Dpsi=dn2,
                                                             mu_est=rel2*lam**3)
            for nm, fld in cases.items():
                a0 = evaluate(fld, (0.0, 0.0), nodes)[0]
                if nm == 'exact' and lam == 1.0:
                    thint, _ = _it.quad(lambda lg: theta_ramp(math.exp(lg), L),
                                        -5.0, L + 5.0, limit=400)
                    OUT[f'a0_control_L{int(L)}'] = dict(quadrature=a0,
                                                        exact_l1=kap*thint,
                                                        rel=abs(a0-kap*thint)/abs(kap*thint))
                rows = []
                for (r, z) in pts:
                    a, dra, dza = evaluate(fld, (r, z), nodes)
                    om = r*fld(np.array([r]), np.array([z]), False)[0][0]
                    g = gamma_op(a, dra, dza, om, r)
                    rows.append(dict(r=r, z=z, rho=math.hypot(r, z),
                                     sinphi=r/math.hypot(r, z),
                                     a=a, r_grad_a=r*math.hypot(dra, dza),
                                     omega=om, Gamma=g))
                gmax = max(row['Gamma'] for row in rows)
                OUT[f'L{int(L)}_lam{lam}_{nm}'] = dict(
                    a0=a0, a0_over_half_ML=a0/(0.5*L),
                    Gamma_max_measured=gmax,
                    max_r_grad_a=max(row['r_grad_a'] for row in rows),
                    max_abs_a=max(abs(row['a']) for row in rows),
                    max_abs_omega=max(abs(row['omega']) for row in rows),
                    rows=rows)
            # convergence control on one point, exact field
            pc = pts[len(pts)//2]
            c1 = evaluate(cases['exact'], pc, nodes)
            c2 = evaluate(cases['exact'], pc, nodes_fine)
            OUT[f'conv_L{int(L)}_lam{lam}'] = dict(
                coarse=list(c1), fine=list(c2),
                rel=[abs(c1[i]-c2[i])/max(abs(c2[i]), 1e-30) for i in range(3)])
    with open('u4_results.json', 'w') as f:
        json.dump(OUT, f, indent=1)
    for k in sorted(OUT):
        v = OUT[k]
        if isinstance(v, dict) and 'rows' in v:
            print(f"{k:28s} a0={v['a0']:.5f} ({v['a0_over_half_ML']:.5f} x ML/2)  "
                  f"Gamma_max={v['Gamma_max_measured']:.5f}  max r|grad a|={v['max_r_grad_a']:.5f} "
                  f" max|om|={v['max_abs_omega']:.5f}")
        else:
            print(f"{k:28s} {v}")

if __name__ == '__main__':
    run()

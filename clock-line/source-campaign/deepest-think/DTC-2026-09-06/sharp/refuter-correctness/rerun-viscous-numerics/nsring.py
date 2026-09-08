#!/usr/bin/env python3
"""
sharp/viscous-numerics : does the ACTUAL viscous axisymmetric no-swirl NS evolution of a
dyadic N-ring datum double sup|omega^theta| in time ~ 1/(M log Re_E)?

Exact system solved (no model, no closure):
    eta = omega^theta / r,   L5 = d_rr + (3/r) d_r + d_zz,   L5 psi1 = -eta,
    Psi = r^2 psi1,   u^r = -(1/r) d_z Psi = r a,  a = -d_z psi1,
    u^z = (1/r) d_r Psi = 2 psi1 + r d_r psi1,     D_t eta = nu L5 eta.

Discretisation follows the estate's validated persistence-toy (psi1 on NODES, eta on CELL
CENTRES, face volume fluxes = differences of Psi so the transport velocity is EXACTLY
divergence free), with the sparse-LU Poisson replaced by a separable DST-I(z) +
tridiagonal(r) direct solver so that multi-octave grids are affordable.

Numerics falsify; they never prove.
"""
import sys, os, json, time, math, hashlib
import numpy as np
from scipy.fft import dst, idst

# ---------------------------------------------------------------- grid ------------

class Grid:
    def __init__(self, rmax, zmin, zmax, h):
        self.h = float(h)
        self.Nr = int(round(rmax / h))
        self.Nz = int(round((zmax - zmin) / h))
        self.z0 = float(zmin)
        self.Rn = np.arange(self.Nr + 1) * self.h
        self.Zn = self.z0 + np.arange(self.Nz + 1) * self.h
        self.Rc = (np.arange(self.Nr) + 0.5) * self.h
        self.Zc = self.z0 + (np.arange(self.Nz) + 0.5) * self.h
    def cell_mesh(self):
        return np.meshgrid(self.Rc, self.Zc, indexing="ij")
    def node_mesh(self):
        return np.meshgrid(self.Rn, self.Zn, indexing="ij")

# ------------------------------------------------- separable Poisson (L5 on nodes) --

class Poisson:
    """L5 psi1 = -eta on nodes i=0..Nr, j=0..Nz with Dirichlet on i=Nr, j=0, j=Nz.
    Axis row uses the regular limit L5 f = 4 f_rr + f_zz, f even in r  ->  8(f1-f0)/h^2.
    bad=True replaces the r^3 weight by r^1 (plain 3D Laplacian) -- the operator control."""
    def __init__(self, g, bad=False):
        self.g = g
        Nr, Nz, h = g.Nr, g.Nz, g.h
        p = 1.0 if bad else 3.0
        lo = np.zeros(Nr); di = np.zeros(Nr); up = np.zeros(Nr)
        c0 = (4.0 if bad else 8.0) / h**2
        di[0] = -c0; up[0] = c0
        i = np.arange(1, Nr)
        r = i * h; rp = (i + 0.5) * h; rm = (i - 0.5) * h
        cp = rp**p / (r**p * h**2); cm = rm**p / (r**p * h**2)
        lo[1:] = cm; di[1:] = -(cp + cm); up[1:] = cp
        self.lo, self.di, self.up = lo, di, up
        K = Nz - 1                                   # z-modes
        k = np.arange(1, Nz)
        self.mu = (2.0 - 2.0 * np.cos(np.pi * k / Nz)) / h**2      # (K,)
        # precompute Thomas factors for  (T_r - mu_k I) x = b
        d = di[:, None] - self.mu[None, :]                          # (Nr,K)
        cpr = np.zeros((Nr, K)); den = np.zeros((Nr, K))
        den[0] = 1.0 / d[0]; cpr[0] = up[0] * den[0]
        for ii in range(1, Nr):
            den[ii] = 1.0 / (d[ii] - lo[ii] * cpr[ii - 1])
            cpr[ii] = up[ii] * den[ii]
        self.cpr, self.den = cpr, den

    def solve(self, eta_node, g_bot=None, g_top=None, g_right=None):
        """eta_node: (Nr+1,Nz+1) values of eta at nodes. Returns psi1 on nodes."""
        g = self.g; Nr, Nz, h = g.Nr, g.Nz, g.h
        b = -eta_node[:Nr, 1:Nz].copy()                             # (Nr, Nz-1)
        if g_bot is not None: b[:, 0] -= g_bot[:Nr] / h**2
        if g_top is not None: b[:, -1] -= g_top[:Nr] / h**2
        if g_right is not None: b[Nr - 1, :] -= self.up[Nr - 1] * g_right[1:Nz]
        bh = dst(b, type=1, axis=1)
        y = np.empty_like(bh); x = np.empty_like(bh)
        y[0] = bh[0] * self.den[0]
        for ii in range(1, Nr):
            y[ii] = (bh[ii] - self.lo[ii] * y[ii - 1]) * self.den[ii]
        x[Nr - 1] = y[Nr - 1]
        for ii in range(Nr - 2, -1, -1):
            x[ii] = y[ii] - self.cpr[ii] * x[ii + 1]
        sol = idst(x, type=1, axis=1)
        psi = np.zeros((Nr + 1, Nz + 1))
        psi[:Nr, 1:Nz] = sol
        if g_bot is not None: psi[:Nr, 0] = g_bot[:Nr]
        if g_top is not None: psi[:Nr, Nz] = g_top[:Nr]
        if g_right is not None: psi[Nr, :] = g_right
        return psi

# --------------------------------------------------------------- operators ---------

def cell_to_node(g, eta, zsym):
    """average cell field (Nr,Nz) to nodes (Nr+1,Nz+1); even in r; zsym='odd' -> eta odd in z."""
    Nr, Nz = g.Nr, g.Nz
    pad = np.zeros((Nr + 2, Nz + 2))
    pad[1:-1, 1:-1] = eta
    pad[0, 1:-1] = eta[0]                       # even across the axis
    pad[-1, 1:-1] = 0.0
    if zsym == "odd":
        pad[1:-1, 0] = -eta[:, 0]               # odd across z = 0
    return 0.25 * (pad[:-1, :-1] + pad[1:, :-1] + pad[:-1, 1:] + pad[1:, 1:])

def face_fluxes(g, psi1):
    h = g.h
    RR, _ = g.node_mesh()
    Psi = RR**2 * psi1
    Fr = -(Psi[:, 1:] - Psi[:, :-1]) / h        # (Nr+1, Nz)  = r u^r at r-faces
    Fz = (Psi[1:, :] - Psi[:-1, :]) / h         # (Nr, Nz+1)  = r u^z at z-faces
    return Fr, Fz

def _upw(qm2, qm1, qp1, qp2, F):
    left = qm1 + (1.0 / 6.0) * (qm1 - qm2) + (1.0 / 3.0) * (qp1 - qm1)
    right = qp1 - (1.0 / 6.0) * (qp2 - qp1) - (1.0 / 3.0) * (qp1 - qm1)
    return np.where(F >= 0.0, left, right)

def advect_rhs(g, eta, Fr, Fz, zsym):
    Nr, Nz, h = g.Nr, g.Nz, g.h
    E = np.zeros((Nr + 4, Nz + 4))
    E[2:-2, 2:-2] = eta
    E[1, 2:-2] = eta[0]; E[0, 2:-2] = eta[1]                  # even across r=0
    if zsym == "odd":
        E[2:-2, 1] = -eta[:, 0]; E[2:-2, 0] = -eta[:, 1]      # odd across z=0
    qm2 = E[0:Nr + 1, 2:-2]; qm1 = E[1:Nr + 2, 2:-2]
    qp1 = E[2:Nr + 3, 2:-2]; qp2 = E[3:Nr + 4, 2:-2]
    fr = Fr * _upw(qm2, qm1, qp1, qp2, Fr)
    qm2 = E[2:-2, 0:Nz + 1]; qm1 = E[2:-2, 1:Nz + 2]
    qp1 = E[2:-2, 2:Nz + 3]; qp2 = E[2:-2, 3:Nz + 4]
    fz = Fz * _upw(qm2, qm1, qp1, qp2, Fz)
    div = (fr[1:, :] - fr[:-1, :]) / h + (fz[:, 1:] - fz[:, :-1]) / h
    return -div / g.Rc[:, None]

def lap5(g, eta, zsym):
    """conservative L5 on cell centres: (1/r^3) d_r(r^3 d_r) + d_zz, no flux at r=0,
    eta = 0 outside the box; zsym='odd' -> eta odd across z=0."""
    Nr, Nz, h = g.Nr, g.Nz, g.h
    r = g.Rc[:, None]
    rp = (np.arange(Nr) + 1.0)[:, None] * h
    rm = (np.arange(Nr) + 0.0)[:, None] * h
    ep = np.zeros_like(eta); em = np.zeros_like(eta)
    ep[:-1] = eta[1:] - eta[:-1]; ep[-1] = -eta[-1]           # eta=0 outside r=rmax
    em[1:] = eta[1:] - eta[:-1]; em[0] = 0.0                  # rm=0 kills it anyway
    out = (rp**3 * ep - rm**3 * em) / (r**3 * h**2)
    zz = np.zeros_like(eta)
    zz[:, 1:-1] = (eta[:, 2:] - 2 * eta[:, 1:-1] + eta[:, :-2]) / h**2
    gb = -eta[:, 0] if zsym == "odd" else 0.0
    zz[:, 0] = (eta[:, 1] - 2 * eta[:, 0] + gb) / h**2
    zz[:, -1] = (0.0 - 2 * eta[:, -1] + eta[:, -2]) / h**2
    return out + zz

def node_velocity(g, psi1):
    h = g.h
    RR, _ = g.node_mesh()
    a = np.zeros_like(psi1)
    a[:, 1:-1] = -(psi1[:, 2:] - psi1[:, :-2]) / (2 * h)
    a[:, 0] = -(psi1[:, 1] - psi1[:, 0]) / h
    a[:, -1] = -(psi1[:, -1] - psi1[:, -2]) / h
    ur = RR * a
    dpr = np.zeros_like(psi1)
    dpr[1:-1, :] = (psi1[2:, :] - psi1[:-2, :]) / (2 * h)
    dpr[0, :] = 0.0
    dpr[-1, :] = (psi1[-1, :] - psi1[-2, :]) / h
    uz = 2.0 * psi1 + RR * dpr
    return ur, uz, a

def energy(g, psi1, zsym):
    """E = ||u||_2^2 (no 1/2) over R^3 = int (ur^2+uz^2) 2 pi r dr dz, trapezoid on nodes."""
    ur, uz, _ = node_velocity(g, psi1)
    RR, _ = g.node_mesh()
    w = np.ones_like(RR); w[0, :] *= 0.5; w[-1, :] *= 0.5; w[:, 0] *= 0.5; w[:, -1] *= 0.5
    E = float(np.sum(w * (ur**2 + uz**2) * 2 * np.pi * RR) * g.h**2)
    return 2.0 * E if zsym == "odd" else E

def interp2(g, F, r, z):
    h = g.h
    rr = abs(r)
    i = min(max(int(rr / h), 0), g.Nr - 1)
    j = min(max(int((z - g.z0) / h), 0), g.Nz - 1)
    x = (rr - g.Rn[i]) / h; y = (z - g.Zn[j]) / h
    return float((1 - x) * (1 - y) * F[i, j] + x * (1 - y) * F[i + 1, j]
                 + (1 - x) * y * F[i, j + 1] + x * y * F[i + 1, j + 1])

# ------------------------------------------------------------------- data ----------

def datum_shell(g, N, rho0=1.0, sub=3):
    """eta_0 = -(2 z / s^2) * Theta(s), Theta = (1/2)[tanh((s-rho0)/w0) - tanh((s-R)/w1)]."""
    R = rho0 * 2.0**N
    w0 = 0.25 * rho0; w1 = 0.10 * R
    RR, ZZ = g.cell_mesh(); h = g.h
    off = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros_like(RR)
    for dr in off:
        for dz in off:
            r = RR + dr * h; z = ZZ + dz * h
            s2 = r * r + z * z; s = np.sqrt(s2)
            Th = 0.5 * (np.tanh((s - rho0) / w0) - np.tanh((s - R) / w1))
            acc += -(2.0 * z / np.maximum(s2, 1e-30)) * Th
    return acc / sub**2, R

def datum_rings(g, N, rho0=1.0, sub=3, gam=0.20):
    """N Gaussian ring-pairs at s_k = rho0 2^k on the 45-degree cone, self-similar cores."""
    R = rho0 * 2.0**(N - 1)
    RR, ZZ = g.cell_mesh(); h = g.h
    off = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros_like(RR)
    c = 1.0 / math.sqrt(2.0)
    for k in range(N):
        sk = rho0 * 2.0**k; sig = gam * sk
        rk = sk * c; zk = sk * c
        Ck = math.sqrt(2.0) / sk
        for dr in off:
            for dz in off:
                r = RR + dr * h; z = ZZ + dz * h
                acc += -Ck * np.exp(-((r - rk)**2 + (z - zk)**2) / (2 * sig**2))
    return acc / sub**2, R

def hill_eta(g, rho=1.0, A=1.0, sub=4, smooth=0.0):
    RR, ZZ = g.cell_mesh(); h = g.h
    off = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros_like(RR)
    for dr in off:
        for dz in off:
            s = np.sqrt((RR + dr * h)**2 + (ZZ + dz * h)**2)
            acc += 0.5 * (1.0 - np.tanh((s - rho) / smooth)) if smooth > 0 else (s < rho).astype(float)
    return A * acc / sub**2

def hill_psi1_exact(g, rho=1.0, A=1.0):
    U = 2.0 * A * rho**2 / 15.0
    RR, ZZ = g.node_mesh()
    S = np.sqrt(RR**2 + ZZ**2)
    inside = 0.75 * U * (1.0 - (S / rho)**2) + 0.5 * U
    outside = 0.5 * U * (rho / np.maximum(S, 1e-12))**3
    return np.where(S < rho, inside, outside), U

def sha_self():
    return hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()

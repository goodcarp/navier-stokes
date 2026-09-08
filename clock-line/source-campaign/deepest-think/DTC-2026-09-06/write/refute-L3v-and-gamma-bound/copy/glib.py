"""Shared library: the strained delta-tapered plateau, its zonal coefficients H_l(lambda),
and composite Gauss-Legendre panels.  No constant is typed from memory."""
import numpy as np
from scipy.special import roots_legendre

_GL = {}
def gl(p):
    if p not in _GL:
        _GL[p] = roots_legendre(p)
    return _GL[p]

def comp_gauss(edges, npanel, p=16):
    """composite p-point Gauss-Legendre: npanel panels on each [edges[i],edges[i+1]]"""
    xs, ws = gl(p)
    X, W = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        e = np.linspace(a, b, npanel + 1)
        for u, v in zip(e[:-1], e[1:]):
            X.append(0.5 * (v - u) * xs + 0.5 * (u + v))
            W.append(0.5 * (v - u) * ws)
    return np.concatenate(X), np.concatenate(W)

def phitilde(theta, lam):
    """polar angle of T_lambda^{-1}x when x has polar angle theta (0<theta<pi)"""
    return np.arctan2(np.sin(theta) / lam, lam**2 * np.cos(theta))

def dphitilde_dtheta(theta, lam):
    return lam**-3 / (np.cos(theta)**2 + lam**-6 * np.sin(theta)**2)

def hdelta(phi, delta):
    if delta <= 0:
        return np.ones_like(phi)
    return np.minimum(1.0, np.minimum(phi, np.pi - phi) / delta)

def theta_delta(lam, delta):
    """the polar angle theta at which phitilde(theta)=delta (the strained taper corner)"""
    return np.arctan(lam**3 * np.tan(delta)) if delta > 0 else 0.0

def W_lambda(theta, lam, delta, M=1.0):
    """rho * eta_lambda at polar angle theta  =  -lambda M sgn(cos th) h(phitilde)/sin th"""
    return -lam * M * np.sign(np.cos(theta)) * hdelta(phitilde(theta, lam), delta) / np.sin(theta)

def omega_lambda(theta, lam, delta, M=1.0):
    """omega^theta of the strained field at polar angle theta"""
    return -lam * M * np.sign(np.cos(theta)) * hdelta(phitilde(theta, lam), delta)

def H_direct(lmax, lam, delta, M=1.0, npanel=None, p=16):
    """instrument I1.  H_l = -(2 lam M/N_l) int_0^{pi/2} h(phitilde) C_l^{3/2}(cos th) sin^2 th dth
       for odd l; 0 for even l."""
    if npanel is None:
        npanel = max(60, int(1.5 * lmax / p) + 20)
    thd = theta_delta(lam, delta)
    edges = [0.0, thd, np.pi / 2] if thd > 0 else [0.0, np.pi / 2]
    th, w = comp_gauss(edges, npanel, p)
    ct = np.cos(th)
    amp = w * hdelta(phitilde(th, lam), delta) * np.sin(th)**2
    H = np.zeros(lmax + 1)
    Cm1 = np.ones_like(ct); C = 3.0 * ct
    for l in range(1, lmax + 1):
        if l % 2 == 1:
            Nl = (l + 1) * (l + 2) / (l + 1.5)
            H[l] = -(2 * lam * M / Nl) * float(amp @ C)
        Cm1, C = C, ((2 * l + 3) * ct * C - (l + 2) * Cm1) / (l + 1)
    return H

_p0 = {0: 1.0}
def Pn0(n):
    """P_n(0) exactly by the two-term recurrence P_n(0) = -(n-1)/n P_{n-2}(0)"""
    if n % 2: return 0.0
    if n not in _p0:
        _p0[n] = Pn0(n - 2) * (-(n - 1.0) / n)
    return _p0[n]

def H_byparts(lmax, lam, delta, M=1.0, npanel=None, p=16):
    """instrument I2 (the representation the proof uses):
       H_l = -(1/N_l)[ [phi](0) P_{l+1}(0) + int phi'_cl P_{l+1} dt ],  phi = W (1-t^2),
       [phi](0) = -2 lam M,  int phi'_cl P_{l+1} dt = -2 int_0^{pi/2} P_{l+1}(cos th) D(th) dth,
       D(th) = d/dth[ -lam M h(phitilde(th)) sin th ]."""
    if npanel is None:
        npanel = max(60, int(1.5 * lmax / p) + 20)
    thd = theta_delta(lam, delta)
    edges = [0.0, thd, np.pi / 2] if thd > 0 else [0.0, np.pi / 2]
    th, w = comp_gauss(edges, npanel, p)
    ct = np.cos(th)
    pt = phitilde(th, lam)
    hp = np.where(pt < delta, 1.0 / delta, 0.0) if delta > 0 else np.zeros_like(pt)
    D = -lam * M * (hp * dphitilde_dtheta(th, lam) * np.sin(th) + hdelta(pt, delta) * np.cos(th))
    H = np.zeros(lmax + 1)
    Pm1 = np.ones_like(ct); P = ct.copy()
    for n in range(1, lmax + 1):
        Pm1, P = P, ((2 * n + 1) * ct * P - n * Pm1) / (n + 1)   # now P = P_{n+1}
        if n % 2 == 1:
            Nl = (n + 1) * (n + 2) / (n + 1.5)
            smooth = -2.0 * float((w * D) @ P)
            H[n] = -(1.0 / Nl) * ((-2 * lam * M) * Pn0(n + 1) + smooth)
    return H

def term(l, Hl):
    """the L3v summand |H_l| ||C_l^{3/2}||_inf / (l(l+3)-4)"""
    return abs(Hl) * ((l + 1) * (l + 2) / 2.0) / (l * (l + 3) - 4.0)

#!/usr/bin/env python3
"""ns5d.py -- 5D-lift Biot-Savart for axisymmetric NO-SWIRL 3D incompressible flow.

DERIVATION (all of it re-derived here, nothing quoted):
  Stokes stream function Psi:  u^r = -(1/r) d_z Psi ,  u^z = (1/r) d_r Psi ,
  omega^theta = d_z u^r - d_r u^z  =>  d_rr Psi - (1/r) d_r Psi + d_zz Psi = -r omega^theta.
  Put eta := omega^theta / r  and  Psi = r^2 psi.  Then the LHS becomes r^2 (d_rr + (3/r) d_r + d_zz) psi,
  so                      Delta_5 psi = -eta,     Delta_5 = d_rr + (3/r) d_r + d_zz,
  i.e. psi is the Newtonian potential in R^5 = R^4_r x R_z of the source eta (radial in the R^4 factor).
  G_5(x) = 1/(8 pi^2 |x|^3)   (since G_d = 1/((d-2)|S^{d-1}| |x|^{d-2}), |S^4| = 8 pi^2/3).
  Velocity:  u^r = -(1/r) d_z (r^2 psi) = -r d_z psi   =>   a := u^r/r = -d_z psi.
             u^z =  (1/r) d_r (r^2 psi) = 2 psi + r d_r psi.
  With d^4 x' = r'^3 dr' sin^2(th) dth dOmega_2  (|S^2| = 4 pi):
      psi(r,z) = (1/(2 pi)) int int eta(r',z') r'^3 J3 dr' dz',   J3 = int_0^pi sin^2 th (A - B cos th)^{-3/2} dth
      a(r,z)   = (3/(2 pi)) int int eta(r',z') r'^3 (z - z') J5 dr' dz',  J5 = int_0^pi sin^2 th (A-B cos th)^{-5/2} dth
  with A = r^2 + r'^2 + (z-z')^2,  B = 2 r r'.

CLOSED FORM for the theta integrals (derived, then checked numerically):
  int_{-1}^{1} sqrt(1-u^2) (A - B u)^{-p} du, u = 1-2w:
      = 4 (A-B)^{-p} int_0^1 sqrt(w(1-w)) (1 + m w)^{-p} dw,  m = 2B/(A-B)
      = 4 (A-B)^{-p} B(3/2,3/2) 2F1(p, 3/2; 3; -m) = (pi/2) (A-B)^{-p} 2F1(p,3/2;3;-m).
  Pfaff  2F1(p,3/2;3;-m) = (1+m)^{-3/2} 2F1(3-p, 3/2; 3; m/(1+m))   and  1+m = (A+B)/(A-B), so
      J5 = (pi/2) 2F1(1/2, 3/2; 3; k2) / [ (A-B) (A+B)^{3/2} ]
      J3 = (pi/2) 2F1(3/2, 3/2; 3; k2) / (A+B)^{3/2}
  with the classical elliptic modulus   k2 = 2B/(A+B) = 4 r r' / ((r+r')^2 + (z-z')^2)  in [0,1].
  A - B = (r-r')^2 + (z-z')^2  is the only singular factor, and 2F1(1/2,3/2;3;.) is FINITE at k2=1
  (c-a-b = 1 > 0; value 16/(3 pi)), so J5 ~ 1/(3 r^3 t^2) at separation t -> 0 -- exactly the flat-space
  codimension-2 behaviour.  2F1(3/2,3/2;3;.) has the expected log divergence at k2 = 1 (ring self-energy).

ENERGY:  E := ||u||_2^2 (no 1/2) = int u.(curl A), A = (Psi/r) e_theta  =>  E = 2 pi int int omega^theta Psi dr dz
         = 2 pi int int eta psi r^3 dr dz = (1/pi) int_{R^5} eta psi = (1/pi) ||grad psi||^2_{L^2(R^5)}.
"""
import numpy as np
from scipy.special import hyp2f1

# ---------------------------------------------------------------- kernels
def J5(r, z, rp, zp):
    dz = z - zp
    sm = (r + rp)**2 + dz*dz          # A + B
    df = (r - rp)**2 + dz*dz          # A - B
    k2 = 4.0*r*rp/sm
    k2 = np.clip(k2, 0.0, 1.0)
    return 0.5*np.pi*hyp2f1(0.5, 1.5, 3.0, k2)/(df*sm**1.5)

# 2F1(3/2,3/2;3;x) has c-a-b = 0, hence a LOG singularity at x = 1 (the ring self-energy log).
# scipy loses it below 1-x ~ 1e-10 and returns inf at x = 1, so splice in the exact leading behaviour
#   2F1(a,b;a+b;x) = (G(c)/(G(a)G(b))) [ -log(1-x) + 2psi(1)-psi(a)-psi(b) ] + O((1-x)log(1-x)),
# with a=b=3/2, c=3:  G(3)/G(3/2)^2 = 8/pi,  2psi(1)-2psi(3/2) = 4 log 2 - 4.
_C0 = 8.0/np.pi
_K0 = 4.0*np.log(2.0) - 4.0
def _F3(one_minus_x):
    om = np.asarray(one_minus_x, float)
    om = np.clip(om, 1e-300, 1.0)
    out = np.empty_like(om)
    big = om > 1e-9
    out[big] = hyp2f1(1.5, 1.5, 3.0, 1.0 - om[big])
    sml = ~big
    out[sml] = _C0*(-np.log(om[sml]) + _K0)
    return out

def J3(r, z, rp, zp):
    dz = z - zp
    sm = (r + rp)**2 + dz*dz          # A + B
    df = (r - rp)**2 + dz*dz          # A - B  ( = 1 - k2, times sm )
    return 0.5*np.pi*_F3(df/sm)/sm**1.5

# ---------------------------------------------------------------- quadrature builders
def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b - a)*x + 0.5*(a + b), 0.5*(b - a)*w

def t_panels(d, sig, nt, ext=7.0):
    """panel edges in t (distance from the evaluation point) resolving a Gaussian of width sig at distance d."""
    hi = d + ext*sig
    lo = d - ext*sig
    mids = [d + f*sig for f in (-5, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 5)]
    edges = [e for e in mids if 0.0 < e < hi]
    if lo > 0.0:
        edges = [0.0, lo] + edges
    else:
        # evaluation point sits inside the core: geometric ladder into t = 0 (integrand ~ t^2 log t)
        lad = [sig*4.0**(-j) for j in range(1, 7)]
        edges = [0.0] + lad + edges
    edges = sorted(set([0.0] + edges + [hi]))
    ts, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        x, w = gl(nt, a, b)
        ts.append(x); ws.append(w)
    return np.concatenate(ts), np.concatenate(ws)

# ---------------------------------------------------------------- ring datum
def make_ring(rk, zk, eta0, C):
    C = np.asarray(C, float)
    return dict(r=float(rk), z=float(zk), eta0=float(eta0), C=C, Cinv=np.linalg.inv(C),
                sig=float(np.sqrt(np.max(np.linalg.eigvalsh(C)))))

def eta_of(rr, zz, ring):
    dr = rr - ring['r']; dz = zz - ring['z']
    Ci = ring['Cinv']
    q = Ci[0, 0]*dr*dr + 2.0*Ci[0, 1]*dr*dz + Ci[1, 1]*dz*dz
    return ring['eta0']*np.exp(-0.5*q)

def _polar_grid(r0, z0, ring, nt, nal):
    d = np.hypot(ring['r'] - r0, ring['z'] - z0)
    t, wt = t_panels(d, ring['sig'], nt)
    al = 2.0*np.pi*np.arange(nal)/nal
    wal = np.full(nal, 2.0*np.pi/nal)
    T, AL = np.meshgrid(t, al, indexing='ij')
    W = np.outer(wt, wal)
    rp = r0 + T*np.cos(AL)
    zp = z0 + T*np.sin(AL)
    return T, AL, W, rp, zp

def a_of_ring(r0, z0, ring, nt=28, nal=384):
    """contribution of one Gaussian eta-ring to a = u^r/r at (r0,z0), polar quadrature about (r0,z0)."""
    T, AL, W, rp, zp = _polar_grid(r0, z0, ring, nt, nal)
    good = rp > 0.0
    e = np.where(good, eta_of(rp, zp, ring), 0.0)
    dz = -T*np.sin(AL)                                   # z0 - z'
    K = np.where(good, J5(r0, z0, np.where(good, rp, 1.0), zp), 0.0)
    integ = (3.0/(2.0*np.pi))*e*np.where(good, rp, 0.0)**3*dz*K*T
    return float(np.sum(integ*W))

def psi_of_ring(r0, z0, ring, nt=28, nal=384):
    T, AL, W, rp, zp = _polar_grid(r0, z0, ring, nt, nal)
    good = rp > 0.0
    e = np.where(good, eta_of(rp, zp, ring), 0.0)
    K = np.where(good, J3(r0, z0, np.where(good, rp, 1.0), zp), 0.0)
    integ = (1.0/(2.0*np.pi))*e*np.where(good, rp, 0.0)**3*K*T
    return float(np.sum(integ*W))

def a_field(r0, z0, rings, **kw):
    return sum(a_of_ring(r0, z0, R, **kw) for R in rings)

def psi_field(r0, z0, rings, **kw):
    return sum(psi_of_ring(r0, z0, R, **kw) for R in rings)

def uz_field(r0, z0, rings, h=None, **kw):
    """u^z = 2 psi + r d_r psi, 4th-order central differences in r (psi is smooth)."""
    if h is None:
        h = 1e-3*max(r0, 1e-12)
    p = [psi_field(r0 + k*h, z0, rings, **kw) for k in (-2, -1, 1, 2)]
    dpsi = (p[0] - 8*p[1] + 8*p[2] - p[3])/(12*h)
    return 2.0*psi_field(r0, z0, rings, **kw) + r0*dpsi

def ur_field(r0, z0, rings, **kw):
    return r0*a_field(r0, z0, rings, **kw)

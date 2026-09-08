#!/usr/bin/env python3
"""b2 -- the EXACT instrument for the affine (pure axisymmetric strain) case.

For b(x,t) = a(t) diag(1,1,1,1,-2) x the backward reverse-time diffusion of Theorem V.1
Step 0 has Z_tau = Y_tau - x_0 EXACTLY Gaussian, mean 0, covariance
        C_tau = diag( 2 sigma_y I_4 , 2 sigma_z ),
        sigma_y = nu \int_0^tau lam(u)^{-2} du ,  sigma_z = nu \int_0^tau lam(u)^{4} du .
Hence  eta(X(tau),tau) = E[ eta_0(x_0 + Z) ]  can be evaluated to quadrature precision:
the 4 y-directions collapse to a noncentral-chi density with 4 dof, the z-direction to a
1D normal.  No PDE solve, no Monte-Carlo error, no numerical diffusion.

Tests
  T0  Monte-Carlo of the backward SDE vs the closed-form covariance  (checks the algebra)
  T1  pure plateau eta_P = -M/r      : relative loss vs s := sigma_y/r0^2 = nu int dt/r(t)^2
  T2  the campaign's smooth taper datum at s0 = 1.5, 2, 3 rho0
  T3  sigma_z-independence of the bulk loss (a sharp falsifiable prediction)
  T4  the third/fourth-order remainder scaling
"""
import json, math
import numpy as np
from scipy.special import ive
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite_e import hermegauss

RES = {}
rng = np.random.default_rng(20260906)

# ---------------------------------------------------------------- the datum (self-contained)
def datum_eta(s, cphi, sphi, M=1.0, rho0=1.0, R=16.0, delta_deg=7.5, w=0.12):
    """the campaign's admissible taper datum, eta = omega^theta / r, as an analytic formula:
       omega^theta = -M tanh(sin phi / sin delta) tanh(cos phi / w) Theta(s)
       eta         = -(M/(s sin delta)) (tanh x / x)|_{x=sin phi/sin delta} tanh(cos phi/w) Theta(s)
       Theta(s) = 1/2 [ tanh((s-rho0)/w0) - tanh((s-R)/w1) ],  w0 = .25 rho0, w1 = .10 R """
    sd = math.sin(math.radians(delta_deg))
    w0 = 0.25*rho0; w1 = 0.10*R
    x = np.maximum(sphi/sd, 1e-12)
    tox = np.tanh(x)/x
    Th = 0.5*(np.tanh((s-rho0)/w0) - np.tanh((s-R)/w1))
    return -(M/(s*sd))*tox*np.tanh(cphi/w)*Th

def plateau_eta(s, cphi, sphi, M=1.0, **kw):
    """the bare plateau eta_P = -M sgn(z)/r  (r = s sin phi)."""
    return -M*np.sign(cphi)/(s*sphi)

# ---------------------------------------------------------------- exact Gaussian average
def gauss_average(fun, r0, z0, sig_y2, sig_z2, nR=900, nZ=241, nrad=12.0):
    """E[ f(|y0+G_y|, z0+G_z) ] with G_y ~ N(0, sig_y2 I_4), G_z ~ N(0, sig_z2).
       f is called as f(s, cos phi, sin phi) with s = sqrt(R^2+z^2)."""
    sy = math.sqrt(sig_y2); sz = math.sqrt(sig_z2)
    Rlo = max(1e-12, r0 - nrad*sy); Rhi = r0 + nrad*sy
    xg, wg = leggauss(nR)
    Rv = 0.5*(Rhi-Rlo)*xg + 0.5*(Rhi+Rlo); wR = 0.5*(Rhi-Rlo)*wg
    # noncentral chi (k=4): p(R) = (R^2/(sig^2 r0)) exp(-(R-r0)^2/(2 sig^2)) * ive(1, R r0/sig^2)
    arg = Rv*r0/sig_y2
    pR = (Rv**2/(sig_y2*r0))*np.exp(-(Rv-r0)**2/(2*sig_y2))*ive(1, arg)
    hx, hw = hermegauss(nZ)                       # weight exp(-x^2/2), sum hw = sqrt(2 pi)
    Zv = z0 + sz*hx; wZ = hw/math.sqrt(2*math.pi)
    RR, ZZ = np.meshgrid(Rv, Zv, indexing='ij')
    SS = np.sqrt(RR**2 + ZZ**2)
    val = fun(SS, ZZ/SS, RR/SS)
    tot = (wR[:, None]*wZ[None, :]*val).sum()
    mass = (wR*pR).sum()*1.0
    tot = (wR[:, None]*pR[:, None]*wZ[None, :]*val).sum()
    return float(tot), float(mass)

# ---------------------------------------------------------------- window / strain model
def window(kappa=0.5, theta=None, ntab=4001):
    """prove-lagrangian (4.1) at the tracked shell sigma=0: lam(theta) = (1-kappa theta/2)^{-2},
       theta = M L t.  Returns tabulated theta grid and lam."""
    if theta is None:
        theta = 2*(1-math.sqrt(2.0/3.0))/kappa      # accelerated doubling window, lam = 3/2
    th = np.linspace(0.0, theta, ntab)
    lam = (1 - kappa*th/2)**-2
    return th, lam, theta

# ---------------------------------------------------------------- T0: covariance check
def T0(nu=2e-3, tau=0.4, npath=200000, nstep=800, a0=1.0):
    """Monte-Carlo the backward SDE with the affine strain; compare Cov(Z_tau) with the
       closed form C_tau = diag(2 nu int lam^{-2}, ..., 2 nu int lam^4)."""
    dt = tau/nstep
    # a(t) constant = a0  ->  lam(u) = exp(a0 u)
    lamf = lambda u: np.exp(a0*u)
    u = np.linspace(0, tau, 20001)
    sy = nu*np.trapz(lamf(u)**-2, u); sz = nu*np.trapz(lamf(u)**4, u)
    # reverse characteristic: x_s = X(tau - s);  Z_s = Y_s - x_s
    Z = np.zeros((npath, 5))
    S = np.diag([a0, a0, a0, a0, -2*a0])
    for k in range(nstep):
        # dZ = -S Z ds + sqrt(2 nu) dW   (affine: A_s = S exactly, no dependence on Z)
        Z = Z - (Z @ S.T)*dt + math.sqrt(2*nu*dt)*rng.standard_normal((npath, 5))
    cov = (Z.T @ Z)/npath
    mean = Z.mean(axis=0)
    pred = np.diag([2*sy]*4 + [2*sz])
    return dict(nu=nu, tau=tau, a0=a0, npath=npath, nstep=nstep,
                sigma_y=float(sy), sigma_z=float(sz),
                cov_diag_mc=[float(cov[i, i]) for i in range(5)],
                cov_diag_pred=[float(pred[i, i]) for i in range(5)],
                rel_err_diag=[float(abs(cov[i, i]-pred[i, i])/pred[i, i]) for i in range(5)],
                max_offdiag_over_trace=float(np.abs(cov-np.diag(np.diag(cov))).max()/np.trace(cov)),
                mean_over_sd=[float(mean[i]/math.sqrt(pred[i, i]/npath)) for i in range(5)])

RES['T0_covariance'] = T0()
print("T0", json.dumps(RES['T0_covariance'], indent=1))

# ---------------------------------------------------------------- T1: pure plateau
def T1():
    rows = []
    r0 = 1.0; z0 = 3.0                     # far from the equator so the sgn jump is invisible
    for s in [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]:
        sig_y2 = 2*s*r0**2                 # sigma_y = s r0^2  ->  variance 2 sigma_y
        sig_z2 = sig_y2                    # (T3 varies this)
        val, mass = gauss_average(plateau_eta, r0, z0, sig_y2, sig_z2, nR=1200, nZ=201, nrad=14.0)
        ref = plateau_eta(np.array([math.hypot(r0, z0)]), np.array([z0/math.hypot(r0, z0)]),
                          np.array([r0/math.hypot(r0, z0)]))[0]
        rel_loss = 1.0 - val/ref
        rows.append(dict(s=s, mass=mass, rel_loss=rel_loss, rel_loss_over_s=rel_loss/s,
                         resid_over_s2=(rel_loss - s)/s**2))
    return rows
RES['T1_pure_plateau'] = T1()
print("\nT1  s, rel_loss/s, (rel_loss-s)/s^2")
for x in RES['T1_pure_plateau']:
    print(f"   {x['s']:9.2e}  mass={x['mass']:.12f}  {x['rel_loss_over_s']:.9f}  {x['resid_over_s2']:.6f}")

# ---------------------------------------------------------------- T3: sigma_z independence
def T3():
    rows = []
    r0 = 1.0; z0 = 3.0; s = 1e-3
    sig_y2 = 2*s*r0**2
    for fac in [0.1, 1.0, 5.0, 25.0, 100.0]:
        val, mass = gauss_average(plateau_eta, r0, z0, sig_y2, sig_y2*fac, nR=1200, nZ=241, nrad=14.0)
        ref = -1.0/r0
        rows.append(dict(sigz_over_sigy=fac, rel_loss=1.0-val/ref))
    return rows
RES['T3_sigma_z_independence'] = T3()
print("\nT3 sigma_z/sigma_y, rel_loss  (theorem: independent)")
for x in RES['T3_sigma_z_independence']:
    print(f"   {x['sigz_over_sigy']:7.2f}  {x['rel_loss']:.12e}")

json.dump(RES, open('b2_results.json', 'w'), indent=1)
print("\nWROTE b2_results.json")

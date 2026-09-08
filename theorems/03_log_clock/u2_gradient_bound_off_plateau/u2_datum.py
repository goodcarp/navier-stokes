#!/usr/bin/env python3
"""
u2 -- the datum constants that feed the two weighted maximum principles.

Datum (D-C): the campaign datum of hk2 (D-B) with the radial ramp written in the
LOG variable (so both radial edges are mollified on the same relative width), which
is the scale-invariant version of hk2's tanh ramp and agrees with it at the inner
edge (rho*Theta'(rho_0) = 1/(2 eps_r) = 2 = rho_0/(2 w_0) at w_0 = 0.25 rho_0):

    omega^theta_0 = -M * T(phi) * Z(phi) * Theta(rho) ,
    T = tanh(sin phi / sin delta)          (axis taper, delta = 7.5 deg)
    Z = tanh(cos phi / w)                  (equatorial mollifier, w = 0.20)
    Theta(rho) = 1/2 [ tanh((log(rho/rho_0) - eps_r)/eps_r)
                       - tanh((log(rho/R) + eps_r)/eps_r) ] ,   eps_r = 0.25

    eta_0 = omega^theta_0/(rho sin phi) = -(M/rho) h(phi) Theta(rho),
    h(phi) := T(phi) Z(phi)/sin phi .

Constants produced (all dimensionless, M = rho_0 = 1):
    E0 = sup rho |eta_0|                 -> feeds P1 (the collar bound)
    G0 = sup rho^2 |grad eta_0|          -> feeds P2 (the gradient bound)
    F0 = rho_0^2 sup |grad eta_0|        -> the flat companion of P2
    Jshell = the shell total-variation density of hk2 Lemma 5.1 (for comparison only)
Also: the same four for the strained field eta_lambda = eta_0 o T_lambda^{-1}, as a
falsification control on the growth factors the maximum principles allow.

Run: python3 u2_datum.py
"""
import json, math
import numpy as np
from scipy import integrate, optimize

M = 1.0; RHO0 = 1.0
DELTA = math.radians(7.5)
W_EQ = 0.20
EPS_R = 0.25

OUT = {'delta_deg': 7.5, 'w_eq': W_EQ, 'eps_r': EPS_R}

def Theta(rho, L):
    R = math.exp(L)
    return 0.5*(math.tanh((math.log(rho) - EPS_R)/EPS_R)
                - math.tanh((math.log(rho/R) + EPS_R)/EPS_R))

def rhoThetap(rho, L):
    """rho * dTheta/drho"""
    R = math.exp(L)
    a1 = (math.log(rho) - EPS_R)/EPS_R
    a2 = (math.log(rho/R) + EPS_R)/EPS_R
    return 0.5*(1.0/math.cosh(a1)**2 - 1.0/math.cosh(a2)**2)/EPS_R

import mpmath as mp
mp.mp.dps = 50
_SD = mp.tanh  # placeholder to keep linters quiet
SD = mp.sin(mp.mpf(DELTA))
WQ = mp.mpf(W_EQ)

def hfun_mp(phi):
    """h(phi) = tanh(sin phi/sin delta) tanh(cos phi/w)/sin phi, 50-digit."""
    p = mp.mpf(phi); s = mp.sin(p); c = mp.cos(p)
    if s == 0:
        return mp.tanh(c/WQ)/SD
    return mp.tanh(s/SD)*mp.tanh(c/WQ)/s

def hprime_mp(phi):
    """analytic d h/d phi, 50-digit (no cancellation at 50 dps)."""
    p = mp.mpf(phi); s = mp.sin(p); c = mp.cos(p)
    T = mp.tanh(s/SD); Z = mp.tanh(c/WQ)
    Tp = mp.sech(s/SD)**2*(c/SD)
    Zp = -mp.sech(c/WQ)**2*(s/WQ)
    if s == 0:
        return mp.mpf(0)
    return (Tp*Z + T*Zp)/s - T*Z*c/s**2

def hfun(phi):
    return float(hfun_mp(phi))

def hprime(phi):
    return float(hprime_mp(phi))

def scan(L, n_rho=4001, n_phi=4001):
    """sup of rho|eta_0|, rho^2|grad eta_0|, |grad eta_0| over the (rho,phi) half-plane."""
    R = math.exp(L)
    lr = np.linspace(math.log(0.2), math.log(R) + 1.0, n_rho)
    rhos = np.exp(lr)
    phis = np.linspace(1e-9, math.pi - 1e-9, n_phi)
    hs = np.array([hfun(p) for p in phis])
    hps = np.array([hprime(p) for p in phis])
    Th = np.array([Theta(r, L) for r in rhos])
    rTh = np.array([rhoThetap(r, L) for r in rhos])
    E = np.max(np.abs(np.outer(Th, hs)))
    A = np.abs(np.outer(rTh - Th, hs))          # rho^2 |d_rho eta| part
    B = np.abs(np.outer(Th, hps))               # rho^2 |ang| part
    Gmat = np.sqrt(A**2 + B**2)
    G = Gmat.max()
    iG = np.unravel_index(Gmat.argmax(), Gmat.shape)
    F = (Gmat/np.outer(rhos**2, np.ones_like(phis))).max()
    iF = np.unravel_index((Gmat/np.outer(rhos**2, np.ones_like(phis))).argmax(), Gmat.shape)
    return dict(E0=float(E), G0=float(G), F0=float(F),
                G0_at=(float(rhos[iG[0]]), math.degrees(phis[iG[1]])),
                F0_at=(float(rhos[iF[0]]), math.degrees(phis[iF[1]])))

for L in (10.0, 40.0):
    OUT[f'datum_L{int(L)}'] = scan(L)

# sup |h| in closed-ish form: attained on the axis, = tanh(1/w)/sin(delta)
OUT['E0_axis_closed'] = math.tanh(1.0/W_EQ)/math.sin(DELTA)
OUT['one_over_sin_delta'] = 1.0/math.sin(DELTA)

# shell total-variation density of hk2 Lemma 5.1 (for comparison with the pointwise route)
S3 = 2.0*math.pi**2
def Jshell():
    f = lambda t: math.sqrt(hfun(math.acos(t))**2 + (1-t*t)*dhdt(t)**2)*(1-t*t)
    def dhdt(t, hh=1e-6):
        return (hfun(math.acos(min(1-1e-12, t+hh))) - hfun(math.acos(max(-1+1e-12, t-hh))))/(2*hh)
    v, _ = integrate.quad(f, -1+1e-9, 1-1e-9, limit=400)
    return S3*v
OUT['J_shell_ac_campaign'] = Jshell()
OUT['note_J'] = 'hk2 (D-B) delta=7.5 w=0.20 reports 65.625910 for the a.c. shell density'

# ---- control: the strained field, exact.  eta_lambda = eta_0 o T_lambda^{-1}
# grad eta_lambda = diag(lam^{-1} I4, lam^2) grad eta_0 (at the label point), and
# |T_lambda alpha| in [lam^{-2}|alpha|, lam|alpha|].  So the exact growth factors are
#   sup rho|eta_lam| <= lam * E0 ,   sup rho^2|grad eta_lam| <= lam^2 * lam^2 * G0
# The maximum principles allow e^{c_G} and e^{3 c_G} with c_G >= 2 log lam.
for lam in (1.0, 1.25, 1.5):
    cG = 2*math.log(lam) if lam > 1 else 0.0
    OUT[f'strain_lam{lam}'] = dict(
        exact_rho_eta_factor=lam,
        maxprin_rho_eta_factor=math.exp(cG),
        exact_rho2_grad_factor=lam**4,
        maxprin_rho2_grad_factor=math.exp(3*cG))

with open('u2_results.json', 'w') as f:
    json.dump(OUT, f, indent=1)
for k, v in OUT.items():
    print(f'{k:26s} {v}')

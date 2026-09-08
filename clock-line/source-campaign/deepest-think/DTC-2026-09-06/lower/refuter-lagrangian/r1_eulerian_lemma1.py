#!/usr/bin/env python3
"""r1: independent (Eulerian, no Jacobian, no material change of variables) check of LEMMA 1.

Derivation used here (independent of the attempt's Lagrangian route):
  3D Biot-Savart, x -> 0 limit:  a(0) = (3/4) \int\int w(r,z) (-z) r^2 (r^2+z^2)^{-5/2} dr dz
  (derived from u = (1/4pi)\int w x (x-y)/|x-y|^3 dV with the azimuthal integral expanded to O(r_x))
  Deformed field, EULERIAN:  w_lam(r,z) = lam * w_0(r/lam, lam^2 z)   [eta materially conserved]
Everything below is a quadrature over EULERIAN (rho,phi); the map T_lam appears only inside
the argument of w_0, never as a Jacobian.
"""
import numpy as np, json
from scipy import integrate

M = 1.0

def h_taper(phi, delta):
    if delta <= 0: return 1.0
    pax = min(phi, np.pi - phi)
    return min(1.0, pax/delta)

def h_sin2(phi, _):
    # profile w = -M sin(2 phi):  encode as w = -M * g(phi), g = sin2phi (signed, no sgn z factor)
    return np.sin(2*phi)

def a_eulerian(lam, L, prof='flat', delta=0.0, nphi=200001, rho_in=1.0):
    """a_lam(0)/(M L) by Eulerian quadrature.  ZERO use of the 5D Jacobian."""
    R = rho_in*np.exp(L)
    phi = np.linspace(1e-12, np.pi-1e-12, nphi)
    # Eulerian point at angle phi, radius rho -> material label (r/lam, lam^2 z)
    # label angle phi0 = atan2(rho sin phi / lam, lam^2 rho cos phi)
    phi0 = np.arctan2(np.sin(phi)/lam, lam**2*np.cos(phi))
    # g(phi) = |label|/rho ; the label radius interval (rho_in,R) maps to a rho-interval of
    # log-length exactly L for EVERY phi, so \int drho/rho = L identically.  We do not use that:
    g = np.sqrt(np.sin(phi)**2/lam**2 + lam**4*np.cos(phi)**2)
    loglen = np.log(R/g) - np.log(rho_in/g)      # computed, not assumed
    if prof == 'flat':
        gprof = np.sign(np.cos(phi0))*np.array([h_taper(p, delta) for p in phi0])
    elif prof == 'sin2':
        gprof = np.sin(2*phi0)
    w_lam = -M*lam*gprof                          # w_lam = lam * w_0(label)
    integrand = 0.75*w_lam*(-np.cos(phi))*np.sin(phi)**2*loglen
    return np.trapz(integrand, phi)/(M*L)

def P_lagrangian(lam, prof='flat', delta=0.0, n=400001):
    """the attempt's own formula, for comparison: P_h(lam)=3\int_0^1 h(v) v^2 (Av^2+B)^{-5/2}dv"""
    A = lam**2 - lam**-4; B = lam**-4
    v = np.linspace(0.0, 1.0, n)
    if prof == 'flat':
        phi0 = np.arcsin(np.clip(v,0,1))
        hv = np.array([h_taper(p, delta) for p in phi0])
        return 3*np.trapz(hv*v**2*(A*v**2+B)**-2.5, v)
    else:
        # w = -M sin2phi  =>  eta = -2M cos phi / rho ; general reduction
        phi = np.linspace(1e-12, np.pi/2-1e-12, n)
        Q = lam**2*np.sin(phi)**2 + lam**-4*np.cos(phi)**2
        return (3.0/2)*np.trapz(2*2*np.cos(phi)**2*np.sin(phi)**3*Q**-2.5, phi)

out = {}
lams = [1.0,1.1,1.25,1.5,1.75,2.0,3.0,5.0]
out['flat_eulerian_vs_lam'] = []
for lam in lams:
    for L in (5.0, 20.0):
        v = a_eulerian(lam, L)*2      # a/(ML) = P/2  ->  P
        out['flat_eulerian_vs_lam'].append([lam, L, v, abs(v-lam)/lam])
out['flat_lagrangian_vs_lam'] = [[lam, P_lagrangian(lam), abs(P_lagrangian(lam)-lam)/lam] for lam in lams]

# tapered: Eulerian vs the attempt's Lagrangian P_h
out['taper_eulerian_vs_lagrangian'] = []
for ddeg in (5.0,7.5,10.0,15.0,30.0):
    d = np.deg2rad(ddeg)
    for lam in (1.0,1.25,1.5):
        e = a_eulerian(lam, 20.0, delta=d)*2
        l = P_lagrangian(lam, delta=d)
        out['taper_eulerian_vs_lagrangian'].append([ddeg, lam, e, l, abs(e-l)/max(l,1e-30)])

# kappa reproduction (lam=1)
out['kappa'] = {'flat': a_eulerian(1.0,20.0),
                'taper_30d': a_eulerian(1.0,20.0,delta=np.deg2rad(30)),
                'taper_15d': a_eulerian(1.0,20.0,delta=np.deg2rad(15)),
                'taper_7.5d': a_eulerian(1.0,20.0,delta=np.deg2rad(7.5)),
                'sin2phi': a_eulerian(1.0,20.0,prof='sin2')}
out['P_sin2_over_P_flat_at_lam1'] = P_lagrangian(1.0,prof='sin2')

# the identity behind Lemma 1, stated Eulerianly: log-length of the deformed shell is L at every phi
phi = np.linspace(1e-9, np.pi-1e-9, 9)
for lam in (1.5,3.0):
    g = np.sqrt(np.sin(phi)**2/lam**2 + lam**4*np.cos(phi)**2)
    out.setdefault('loglength_invariance',[]).append([lam, float(np.max(np.abs(np.log(np.exp(7.0)/g)-np.log(1.0/g)-7.0)))])
print(json.dumps(out, indent=1))
json.dump(out, open('r1_results.json','w'), indent=1)

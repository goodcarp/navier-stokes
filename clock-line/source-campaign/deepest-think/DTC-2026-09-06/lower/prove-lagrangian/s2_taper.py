#!/usr/bin/env python3
"""s2_taper.py -- the axis taper: admissibility, kappa loss, and the strain-monotonicity margin.

Datum:  omega^theta_0 = -M sgn(z) h_delta(phi),  h_delta(phi) = min(1, phi_ax/delta),
        phi_ax = min(phi, pi-phi).  Then |omega_0|_inf = M and eta = omega^theta/r obeys
        |eta| <= M/(rho_0 sin delta) -- ADMISSIBLE (the bang-bang itself is not).

Two things are checked here:
 (T1) EXACT deficit bound.  P_{h_delta}(lam) = lam - D(lam,delta) with
      0 <= D(lam,delta) = 3 int_0^{w}(1-h) v^2 (A v^2+B)^{-5/2} dv
                       <= 3 int_0^{w} v^2 (A v^2+B)^{-5/2} dv = w^3 / ( B (A w^2+B)^{3/2} ),
      w = sin delta, A = lam^2-lam^-4, B = lam^-4   (same antiderivative as Lemma 1).
 (T2) MONOTONICITY MARGIN.  The theorem only needs Phi_h(lam) := P_h(lam)/P_h(1) >= 1 for
      lam in [1, 3/2].  Find the largest delta for which that holds with margin.
"""
import numpy as np, json
OUT = {}

def h_taper(delta):
    def h(phi):
        pax = np.minimum(phi, np.pi-phi)
        return np.minimum(1.0, pax/delta)
    return h

def P_h(h, lam, n=400001):
    A = lam**2 - lam**-4; B = lam**-4
    phi = np.linspace(0.0, np.pi/2, n); v = np.sin(phi)
    return 3.0*np.trapz(h(phi)*v**2*(A*v**2+B)**-2.5*np.cos(phi), phi)

def D_bound(lam, delta):
    w = np.sin(delta); A = lam**2 - lam**-4; B = lam**-4
    return w**3/(B*(A*w**2+B)**1.5)

lams = np.linspace(1.0, 1.5, 26)
tab = {}
for dd in [30.0, 20.0, 15.0, 10.0, 7.5, 5.0, 3.0]:
    d = np.deg2rad(dd); h = h_taper(d)
    P1 = P_h(h, 1.0)
    P = np.array([P_h(h, L) for L in lams])
    Phi = P/P1
    Dtrue = np.array([L - P_h(h, L) for L in lams])
    Dbnd  = np.array([D_bound(L, d) for L in lams])
    tab[f'{dd}deg'] = dict(kappa=0.5*P1,
                           Phi_min=float(Phi.min()), Phi_at_1_5=float(Phi[-1]),
                           D_true_at_1_5=float(Dtrue[-1]), D_bound_at_1_5=float(Dbnd[-1]),
                           bound_is_upper=bool(np.all(Dbnd >= Dtrue - 1e-12)),
                           Phi_monotone_ge_1=bool(np.all(Phi >= 1.0 - 1e-12)))
OUT['taper_table'] = tab
OUT['lams'] = lams.tolist()

# the conservative fully-explicit margin at delta = 7.5 deg
d = np.deg2rad(7.5)
OUT['explicit_margin_7.5deg'] = {
    'kappa_delta': 0.5*P_h(h_taper(d), 1.0),
    'P_h(1)_lower_bound_1_minus_w3': 1 - np.sin(d)**3,
    'lam_minus_Dbound_at_1.5': 1.5 - D_bound(1.5, d),
    'guaranteed_Phi_lower_at_1.5': (1.5 - D_bound(1.5, d))/1.0,
}
print(json.dumps(OUT, indent=1, default=str))
json.dump(OUT, open('s2_results.json','w'), indent=1, default=str)

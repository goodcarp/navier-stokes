"""g9 — the proved jump-discontinuity bound  |H_l(lambda)| <= sqrt(2/pi)(2 lam M + V_lam) l^{-3/2},
uniform on lambda in [1,3/2], and its verification against the computed coefficients.

Route (all steps elementary):
  C_l^{3/2} = P_{l+1}'  (g1),  so with phi(t) := W(t)(1-t^2) = -lam M sgn(t) h(phitilde) sin(phi),
      int_{-1}^1 W C_l^{3/2}(1-t^2) dt = int phi P_{l+1}' dt = -[phi](0) P_{l+1}(0) - int phi'_cl P_{l+1} dt
  (phi vanishes at t = +-1, and its only jump is at t = 0, of size [phi](0) = -2 lam M).
  Szego's bound  |P_n(cos th)| < sqrt(2/(pi n sin th))  then gives
      |H_l| <= (1/N_l) sqrt(2/(pi(l+1))) ( 2 lam M + V_lam ),
      V_lam := int_{-1}^1 |phi'_cl(t)| (1-t^2)^{-1/4} dt = 2 int_0^{pi/2} |phi'(cos th)| sin^{1/2}th dth,
  and  l^{3/2}(l+3/2)/((l+1)(l+2)sqrt(l+1)) <= 1 for l >= 3  (checked), whence the display.
"""
import json
import numpy as np
import glib

M = 1.0
out = {}

# ---- Szego bound, numerically falsified over n <= 600 --------------------------------
x = np.cos(np.linspace(1e-9, np.pi - 1e-9, 400001))
w = (1 - x**2)**0.25
P0, P1 = np.ones_like(x), x.copy()
worst_n, worst_nh = 0.0, 0.0
for n in range(1, 601):
    P0, P1 = P1, ((2 * n + 1) * x * P1 - n * P0) / (n + 1)
    v = np.max(np.abs(P1) * w)           # P1 is now P_{n+1}
    worst_n = max(worst_n, v / np.sqrt(2 / (np.pi * (n + 1))))
    worst_nh = max(worst_nh, v / np.sqrt(2 / (np.pi * (n + 1.5))))
out['szego_worst_ratio_n'] = float(worst_n)
out['szego_worst_ratio_n_plus_half'] = float(worst_nh)

# ---- the prefactor sup --------------------------------------------------------------
ll = np.arange(3, 10**7, 2, dtype=float)
pref = ll**1.5 * (ll + 1.5) / ((ll + 1) * (ll + 2) * np.sqrt(ll + 1))
out['prefactor_sup_l_ge_3'] = float(np.max(pref))
out['prefactor_at_l3'] = float(pref[0])

# ---- V_lambda -----------------------------------------------------------------------
def phiprime(th, lam, delta):
    """phi'(t) for t = cos th in (0,1):  phi(t) = -lam M h(phitilde) sin th,
       phi'(t) = lam M [ h'(phitilde) dphitilde/dphi + h(phitilde) t/s ]"""
    pt = glib.phitilde(th, lam)
    hp = np.where(pt < delta, 1.0 / delta, 0.0) if delta > 0 else np.zeros_like(pt)
    return lam * M * (hp * glib.dphitilde_dtheta(th, lam) + glib.hdelta(pt, delta) * np.cos(th) / np.sin(th))

def V(lam, delta, n=200001):
    """V = 2 int_0^{pi/2} |phi'(cos th)| sqrt(sin th) dth.  For delta = 0 the integrand has an
       integrable th^{-1/2} singularity at the axis, so a graded mesh (th ~ j^4) is used there."""
    thd = glib.theta_delta(lam, delta)
    edges = [0.0, thd, np.pi / 4, np.pi / 2] if 0 < thd < np.pi / 4 else \
            ([0.0, np.pi / 4, np.pi / 2] if thd == 0 else [0.0, np.pi / 4, thd, np.pi / 2])
    edges = sorted(set([e for e in edges if e <= np.pi / 2]))
    tot = 0.0
    for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        if a == 0.0:
            x = np.linspace(0.0, 1.0, n)
            th = a + (b - a) * x**4
            th[0] = 1e-16
        else:
            th = np.linspace(a, b, n)
        y = np.abs(phiprime(th, lam, delta)) * np.sqrt(np.sin(th))
        tot += np.trapz(y, th)
    return 2.0 * tot

# finite-difference check of phiprime
for lam, dd in [(1.0, 7.5), (1.5, 7.5), (1.25, 30.0)]:
    delta = np.deg2rad(dd)
    th = np.array([0.4, 0.8, 1.2, 1.5])
    t = np.cos(th); h = 1e-6
    def PH(tv):
        thv = np.arccos(tv)
        return -lam * M * glib.hdelta(glib.phitilde(thv, lam), delta) * np.sin(thv)
    fd = (PH(t + h) - PH(t - h)) / (2 * h)
    an = phiprime(th, lam, delta)
    out['phiprime_fd_check_lam%.2f_d%.1f' % (lam, dd)] = float(np.max(np.abs(fd - an) / np.abs(an)))

rows = {}
for dd in [0.0, 7.5, 15.0, 30.0]:
    delta = np.deg2rad(dd)
    for lam in [1.0, 1.1, 1.25, 1.4, 1.5]:
        Vv = V(lam, delta)
        C = np.sqrt(2 / np.pi) * (2 * lam * M + Vv)
        H = glib.H_direct(4001, lam, delta, npanel=400)
        ls = np.arange(3, 4002, 2)
        ratio = np.max(np.abs(H[ls]) * ls**1.5) / C
        rows['lam%.2f_del%.1f' % (lam, dd)] = {
            'V': float(Vv), 'C_bound': float(C),
            'max_l32_absH': float(np.max(np.abs(H[ls]) * ls**1.5)),
            'bound_over_truth': float(1.0 / ratio),
            'l32_absH_at_4001': float(abs(H[4001]) * 4001**1.5),
            'analytic_V_upper': float(4 * lam * M + (4.0 / 3.0) * lam**2.5 * M * np.sqrt(max(delta, 1e-16)) * 2)
            if dd > 0 else float(4 * lam * M),
        }
        print('%4.1f %.2f  V=%.4f  C=%.4f  max l^1.5|H_l|=%.4f  bound/truth=%.2f'
              % (dd, lam, Vv, C, rows['lam%.2f_del%.1f' % (lam, dd)]['max_l32_absH'],
                 1.0 / ratio), flush=True)
out['rows'] = rows

# ---- the resulting uniform L3v bound ------------------------------------------------
from scipy.special import zeta
tailsum = float(zeta(1.5) - 1 - 2**-1.5)                      # sum_{l>=3} l^{-3/2}
ratio_sup = float(np.max((ll + 1) * (ll + 2) / (2 * (ll + 4) * (ll - 1))))
out['sum_l_ge_3_l^-1.5'] = tailsum
out['C_ratio_sup'] = ratio_sup
out['uniform_L3v_bound'] = {k: float(v['C_bound'] * tailsum * ratio_sup) for k, v in rows.items()}
print(json.dumps({'sum': tailsum, 'ratio_sup': ratio_sup,
                  'szego_n': out['szego_worst_ratio_n'],
                  'prefactor_sup': out['prefactor_sup_l_ge_3']}, indent=1))
with open('g9_results.json', 'w') as f:
    json.dump(out, f, indent=1)

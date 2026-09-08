"""r5 - independent (light) re-computation of the L3v sum
        S(lam,delta) = sum_{l>=3 odd} |H_l| ||C_l^{3/2}||_inf / ((l+4)(l-1)).
Own instrument: composite 32-point Gauss-Legendre (the seat uses 16-point), 800 panels per
sub-interval, panels broken at the strained taper corner; partial sums to l = 6001; analytic
two-sided tail from the exact jump coefficients (Wallis-bounded remainder) plus a MEASURED
(not fitted) residual envelope.  QUADPACK cross-check of H_l at l = 3, 51.
"""
import json
import numpy as np
from scipy.integrate import quad

LMAX = 6001
MM = 600000

def phit(th, lam): return np.arctan2(np.sin(th) / lam, lam**2 * np.cos(th))
def hdel(p, dl):  return np.ones_like(p) if dl <= 0 else np.minimum(1.0, np.minimum(p, np.pi - p) / dl)
def theta_delta(lam, dl): return float(np.arctan(lam**3 * np.tan(dl))) if dl > 0 else 0.0

def gauss_panels(edges, npanel, p=32):
    xs, ws = np.polynomial.legendre.leggauss(p)
    X, W = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        e = np.linspace(a, b, npanel + 1)
        for uu, vv in zip(e[:-1], e[1:]):
            X.append(0.5 * (vv - uu) * xs + 0.5 * (uu + vv)); W.append(0.5 * (vv - uu) * ws)
    return np.concatenate(X), np.concatenate(W)

def Hs(lam, dl, lmax=LMAX, npanel=800, p=32):
    thd = theta_delta(lam, dl)
    edges = [0.0, thd, np.pi / 2] if thd > 0 else [0.0, np.pi / 2]
    th, w = gauss_panels(edges, npanel, p)
    ct = np.cos(th)
    amp = w * hdel(phit(th, lam), dl) * np.sin(th)**2
    H = np.zeros(lmax + 2)
    a, b = np.ones_like(ct), 3.0 * ct
    for l in range(1, lmax + 1):
        if l % 2 == 1:
            H[l] = -(2 * lam * (l + 1.5) / ((l + 1) * (l + 2))) * float(amp @ b)
        a, b = b, ((2 * l + 3) * ct * b - (l + 2) * a) / (l + 1)
    return H

P0 = np.zeros(LMAX + 2); P0[0] = 1.0
for n in range(2, LMAX + 2, 2): P0[n] = P0[n - 2] * (-(n - 1.0) / n)
def Hjump(l, lam): return (2 * lam * (l + 1.5) / ((l + 1) * (l + 2))) * P0[l + 1]
def term(l, H):    return abs(H) * ((l + 1) * (l + 2) / 2.0) / ((l + 4.0) * (l - 1.0))

mm = np.arange(1, MM + 1, dtype=np.float64)
cbin = np.exp(np.cumsum(np.log((2 * mm - 1) / (2 * mm))))
def jump_tail(lam, L0):
    m0 = (L0 + 1) // 2
    l = 2.0 * mm[m0:] - 1.0
    s = float(lam * np.sum((l + 1.5) * cbin[m0:] / ((l + 4.0) * (l - 1.0))))
    rem = float(lam / np.sqrt(np.pi) / np.sqrt(float(MM) - 1.0))
    return s, rem

out, rows, qp = {}, {}, {}
for dd in [0.0, 7.5, 15.0, 30.0]:
    dl = np.deg2rad(dd)
    for lam in [1.0, 1.25, 1.5]:
        H = Hs(lam, dl)
        ls = np.arange(3, LMAX + 1, 2)
        partial = float(sum(term(l, H[l]) for l in ls))
        Rl = np.array([abs(H[l] - Hjump(l, lam)) for l in ls])
        sel = ls >= 1001
        Cenv = float(np.max(Rl[sel] * ls[sel]**2.0))
        resid_tail = float(Cenv * (5.0 / 14.0) / (2.0 * LMAX))
        jt, jrem = jump_tail(lam, LMAX)
        S_lo, S_hi = partial + jt - jrem - resid_tail, partial + jt + jrem + resid_tail
        rows['lam%.2f_del%.1f' % (lam, dd)] = {'partial': partial, 'jump_tail': jt,
            'jump_tail_rem': jrem, 'resid_env_C': Cenv, 'resid_tail_bound': resid_tail,
            'S_lo': S_lo, 'S_hi': S_hi, 'S_mid': 0.5 * (S_lo + S_hi), 'halfwidth': 0.5 * (S_hi - S_lo)}
        print('%4.1f %.2f  S in [%.5f, %.5f]  mid %.5f' % (dd, lam, S_lo, S_hi, 0.5 * (S_lo + S_hi)), flush=True)
        if lam in (1.0, 1.5) and dd in (0.0, 7.5):
            def Cl(l, x):
                a, b = 1.0, 3.0 * x
                for n in range(1, l): a, b = b, ((2 * n + 3) * x * b - (n + 2) * a) / (n + 1)
                return b
            for l in [3, 51]:
                thd = theta_delta(lam, dl)
                pts = [0.0, thd, np.pi / 2] if thd > 0 else [0.0, np.pi / 2]
                tot = 0.0
                for aa, bb in zip(pts[:-1], pts[1:]):
                    v, _ = quad(lambda th: (1.0 if dl <= 0 else min(1.0, min(phit(th, lam), np.pi - phit(th, lam)) / dl))
                                * Cl(l, np.cos(th)) * np.sin(th)**2, aa, bb, limit=300, epsabs=1e-12, epsrel=1e-12)
                    tot += v
                hq = -(2 * lam * (l + 1.5) / ((l + 1) * (l + 2))) * tot
                qp['lam%.2f_del%.1f_l%d' % (lam, dd, l)] = float(abs(H[l] - hq) / abs(hq))
out['rows'] = rows; out['quadpack_reldiff'] = qp
out['max_quadpack_reldiff'] = max(qp.values())
seat = json.load(open('../copy/g5_results.json'))['L3v_sums']
cmp = {}
for k, v in rows.items():
    if k in seat:
        cmp[k] = {'seat_S': seat[k]['S'], 'mine_lo': v['S_lo'], 'mine_hi': v['S_hi'],
                  'inside': bool(v['S_lo'] <= seat[k]['S'] <= v['S_hi']),
                  'seat_minus_my_mid': seat[k]['S'] - v['S_mid']}
out['vs_seat'] = cmp
out['all_inside'] = bool(all(v['inside'] for v in cmp.values()))
print(); print('max QUADPACK reldiff', out['max_quadpack_reldiff'])
for k, v in cmp.items():
    print('%-18s seat %.4f  mine [%.4f, %.4f]  inside=%s  seat-mid %+.5f'
          % (k, v['seat_S'], v['mine_lo'], v['mine_hi'], v['inside'], v['seat_minus_my_mid']))
print('all seat values inside my two-sided interval:', out['all_inside'])
json.dump(out, open('r5_results.json', 'w'), indent=1)

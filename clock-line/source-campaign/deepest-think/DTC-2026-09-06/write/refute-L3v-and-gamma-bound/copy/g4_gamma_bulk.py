"""g4 — the bulk angular quantities that enter the Gamma bound.

The bulk (homogeneous-degree-1) deviatoric stream function is  psi_dev = rho f(t),
    f(t) = sum_{l>=3, odd} c_l C_l^{3/2}(t),   c_l = H_l/((l+4)(l-1)),
and, since  d_z[rho^gamma C_l(t)] = rho^{gamma-1}[gamma t C_l + (1-t^2) C_l'],  with gamma = 1,
    a_dev(t) = -(t f + (1-t^2) f') = -sum_l c_l [ (l+4) t C_l(t) - (l+1) C_{l+1}(t) ]
using  (1-t^2) C_n' = (n+3) t C_n - (n+1) C_{n+1}   (verified symbolically here).
The third quantity is obtained from the exact identity proved in g3,
    (1-t^2) a'_dev = 3 f + 3 t a_dev + (1-t^2) W_dev ,
and cross-checked against a finite-difference of the a_dev series away from the jump/corner.
"""
import json, time
import numpy as np
import sympy as sp
import glib

M = 1.0
LMAX = 4001
out = {'LMAX': LMAX}
t0 = time.time()

# ---- the Gegenbauer derivative identity used above ---------------------------------
tsym = sp.symbols('t')
resid = []
for n in range(0, 12):
    lhs = (1 - tsym**2) * sp.diff(sp.gegenbauer(n, sp.Rational(3, 2), tsym), tsym)
    rhs = (n + 3) * tsym * sp.gegenbauer(n, sp.Rational(3, 2), tsym) - (n + 1) * sp.gegenbauer(n + 1, sp.Rational(3, 2), tsym)
    resid.append(sp.simplify(sp.expand(lhs - rhs)))
out['gegenbauer_deriv_identity_residuals'] = [str(x) for x in resid]
assert all(x == 0 for x in resid)

def bulk(lam, delta, tgrid, lmax=LMAX):
    H = glib.H_direct(lmax, lam, delta, npanel=max(60, int(1.5 * lmax / 16) + 20))
    ls = np.arange(3, lmax + 1, 2)
    c = np.zeros(lmax + 3)
    c[ls] = H[ls] / ((ls + 4.0) * (ls - 1.0))
    f = np.zeros_like(tgrid); ad = np.zeros_like(tgrid)
    fs, ads = [], []
    Cm1, C0 = np.ones_like(tgrid), 3.0 * tgrid       # C_0, C_1
    Cs = [Cm1, C0]
    prev, cur = Cm1, C0
    for l in range(1, lmax + 1):
        nxt = ((2 * l + 3) * tgrid * cur - (l + 2) * prev) / (l + 1)   # C_{l+1}
        L = l
        if L >= 3 and L % 2 == 1:
            f += c[L] * cur                                    # cur = C_L
            ad += -c[L] * ((L + 4) * tgrid * cur - (L + 1) * nxt)
        prev, cur = cur, nxt
        fs.append(f.copy()); ads.append(ad.copy())
    fs = np.array(fs); ads = np.array(ads)
    k = fs.shape[0]
    return fs[-1], ads[-1], fs[k // 2:].mean(0), ads[k // 2:].mean(0), H

res = {}
for dd in [0.0, 7.5, 15.0, 30.0]:
    delta = np.deg2rad(dd)
    for lam in [1.0, 1.25, 1.5]:
        th = np.linspace(1e-4, np.pi / 2, 4001)     # northern half; equator is an endpoint
        tg = np.cos(th)
        fp, ap, fc, ac, H = bulk(lam, delta, tg)
        om = glib.omega_lambda(th, lam, delta)
        W = glib.W_lambda(th, lam, delta)
        Wdev = W - 3 * H[1] * tg
        ident = 3 * fc + 3 * tg * ac + (1 - tg**2) * Wdev
        # finite-difference cross-check of the identity, away from the corner and the equator
        d = np.gradient(ac, tg)
        fd = (1 - tg**2) * d
        thd = glib.theta_delta(lam, delta)
        ok = (np.abs(th - np.pi / 2) > 0.05) & (np.abs(th - thd) > 0.05) & (th > 0.05)
        res['lam%.2f_del%.1f' % (lam, dd)] = {
            'sup_f': float(np.max(np.abs(fc))),
            'sup_adev': float(np.max(np.abs(ac))),
            'sup_s2_adevprime': float(np.max(np.abs(ident))),
            'identity_vs_finite_difference_maxdiff': float(np.max(np.abs(ident - fd)[ok])),
            'identity_vs_fd_scale': float(np.max(np.abs(ident)[ok])),
            'sup_s2_Wdev': float(np.max(np.abs((1 - tg**2) * Wdev))),
            'sup_omega': float(np.max(np.abs(om))),
            'partial_vs_cesaro_f': float(np.max(np.abs(fp - fc))),
            'partial_vs_cesaro_adev': float(np.max(np.abs(ap - ac))),
            'H1': float(H[1]), 'kappa': float(-0.6 * H[1]),
        }
        r = res['lam%.2f_del%.1f' % (lam, dd)]
        print('%4.1f %.2f  sup|f| %.4f  sup|a_dev| %.4f  sup s^2|a_dev\'| %.4f   id-vs-fd %.2e (scale %.3f)'
              % (dd, lam, r['sup_f'], r['sup_adev'], r['sup_s2_adevprime'],
                 r['identity_vs_finite_difference_maxdiff'], r['identity_vs_fd_scale']), flush=True)

out['cases'] = res
out['secs'] = time.time() - t0
with open('g4_results.json', 'w') as fjs:
    json.dump(out, fjs, indent=1)
print('secs', out['secs'])

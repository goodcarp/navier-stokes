"""g8 — assembly of the Gamma bound, plus the two external verifications.

(1) closed-form representation of a_dev and the two kernel constants 3pi/4 and 2/3,
(2) shell-boundary (edge) correction sums at one octave,
(3) reproduction of far-near-kernel-lemma's inner-edge profile c_edge(phi) at lambda = 1
    (which is also the refuter's material-point offset at lambda = 1),
(4) the assembled constant C' in   Gamma <= 2 a(0) + C' M,   proved and measured columns.
"""
import json
import numpy as np
import sympy as sp
import glib

M = 1.0
LMAX = 4001
out = {}

# ---------------------------------------------------------------- 0. kernel monotonicity, proved
u = sp.symbols('u', positive=True)
for p, name in [(sp.Rational(1, 2), 'p_half'), (sp.Integer(1), 'p_one')]:
    tau = sp.symbols('tau')
    J = sp.integrate((1 - tau**2)**p, (tau, u, 1))
    D = (1 - u**2)**(p + 1) - 3 * u * J
    out['kernel_monotone_D_prime_' + name] = str(sp.simplify(sp.diff(D, u)))
    out['kernel_monotone_D_at_1_' + name] = str(sp.simplify(sp.limit(D, u, 1)))
out['K_f_sup_exact'] = float(3 * sp.pi / 4)
out['K_omega_sup_exact'] = 2.0 / 3.0

# ---------------------------------------------------------------- 1. series + closed form
def bulk(lam, delta, tgrid, lmax=LMAX):
    H = glib.H_direct(lmax, lam, delta, npanel=max(60, int(1.5 * lmax / 16) + 20))
    c = np.zeros(lmax + 3)
    ls = np.arange(3, lmax + 1, 2)
    c[ls] = H[ls] / ((ls + 4.0) * (ls - 1.0))
    f = np.zeros_like(tgrid); ad = np.zeros_like(tgrid)
    fs, ads = [], []
    prev, cur = np.ones_like(tgrid), 3.0 * tgrid
    for l in range(1, lmax + 1):
        nxt = ((2 * l + 3) * tgrid * cur - (l + 2) * prev) / (l + 1)
        if l >= 3 and l % 2 == 1:
            f += c[l] * cur
            ad += -c[l] * ((l + 4) * tgrid * cur - (l + 1) * nxt)
        prev, cur = cur, nxt
        fs.append(f.copy()); ads.append(ad.copy())
    fs, ads = np.array(fs), np.array(ads)
    k = fs.shape[0]
    return fs[k // 2:].mean(0), ads[k // 2:].mean(0), H, c

def edge_sums(H, c, lmax=LMAX):
    """the two shell-boundary corrections, evaluated at one octave in (factors 2^-(l+4), 2^-(l-1)),
       as sup-norm bounds on their contribution to a and to r|grad a|."""
    a_in = a_out = g_in = g_out = 0.0
    for l in range(3, min(lmax, 200) + 1, 2):
        h = abs(H[l])
        supC = lambda n: (n + 1) * (n + 2) / 2.0
        from math import comb
        supdC = lambda n: 0.0 if n == 0 else 3.0 * comb(n + 3, n - 1)
        # inner boundary term:  -H_l rho_in^{l+4} rho^{-(l+3)} C_l /((l+4)(2l+3))
        #   a-contribution  = -(that)_z = -H_l (l+1) (rho_in/rho)^{l+4} C_{l+1}/((l+4)(2l+3))
        w = 2.0 ** (-(l + 4))
        a_in += h * (l + 1) * supC(l + 1) * w / ((l + 4) * (2 * l + 3))
        n = l + 1                                     # exterior harmonic rho^{-(n+3)}C_n
        g_in += h * (l + 1) * ((n + 3) * supC(n) + supdC(n)) * w / ((l + 4) * (2 * l + 3))
        # outer boundary term:  -H_l rho^l R^{1-l} C_l /((l-1)(2l+3))
        w2 = 2.0 ** (-(l - 1))
        a_out += h * (l + 2) * supC(l - 1) * w2 / ((l - 1) * (2 * l + 3))
        n = l - 1                                     # interior harmonic rho^n C_n
        g_out += h * (l + 2) * (n * supC(n) + supdC(n)) * w2 / ((l - 1) * (2 * l + 3))
    return a_in, a_out, g_in, g_out

# ---------------------------------------------------------------- 2. c_edge at lambda = 1
def c_edge(delta, phis, lmax=LMAX):
    """a(rho_0+, phi) - kappa M log(R/rho_0) in the limit R/rho_0 -> infinity"""
    tg = np.cos(phis)
    fc, ac, H, c = bulk(1.0, delta, tg, lmax)
    s = np.zeros_like(tg)
    prev, cur = np.ones_like(tg), 3.0 * tg
    for l in range(1, lmax + 1):
        nxt = ((2 * l + 3) * tg * cur - (l + 2) * prev) / (l + 1)
        if l >= 3 and l % 2 == 1:
            s += -H[l] * (l + 1) * nxt / ((l + 4) * (2 * l + 3))    # inner-boundary a-contribution
        prev, cur = cur, nxt
    return ac + s, H

phis = np.deg2rad(np.array([1, 2, 3, 5, 7.5, 10, 15, 20, 30, 45, 60, 75, 90.0]))
ce_bb, Hbb = c_edge(0.0, phis)
ce_15, _ = c_edge(np.deg2rad(15.0), phis)
ce_75, _ = c_edge(np.deg2rad(7.5), phis)
out['c_edge_phis_deg'] = [float(np.rad2deg(p)) for p in phis]
out['c_edge_bangbang'] = [float(v) for v in ce_bb]
out['c_edge_taper_15deg'] = [float(v) for v in ce_15]
out['c_edge_taper_7.5deg'] = [float(v) for v in ce_75]
out['c_edge_bangbang_at_10deg'] = float(ce_bb[5])
out['c_edge_bangbang_at_30deg'] = float(ce_bb[8])

# ---------------------------------------------------------------- 3. the assembly
g6 = json.load(open('g6_results.json'))
g5 = json.load(open('g5_results.json'))
rows = {}
for dd in [0.0, 7.5, 15.0, 30.0]:
    delta = np.deg2rad(dd)
    for lam in [1.0, 1.25, 1.5]:
        th = np.linspace(1e-4, np.pi / 2, 4001)
        tg = np.cos(th)
        fc, ac, H, c = bulk(lam, delta, tg)
        om = glib.omega_lambda(th, lam, delta)
        W = glib.W_lambda(th, lam, delta)
        Wdev = W - 3 * H[1] * tg
        s2ad = 3 * fc + 3 * tg * ac + (1 - tg**2) * Wdev
        # closed-form a_dev
        thd = glib.theta_delta(lam, delta)
        integ = 3 * fc * np.sqrt(1 - tg**2) + (1 - tg**2) * om - 3 * H[1] * tg * (1 - tg**2)**1.5
        y = integ * np.sin(th)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(th))])
        a_closed = np.where(np.abs(tg) < 0.999, -(1 - tg**2)**-1.5 * cum, np.nan)
        mask = np.abs(tg) < 0.99
        kap = -0.6 * H[1]
        om_eff = om - 3 * H[1] * tg * np.sqrt(1 - tg**2)
        a_in, a_out, g_in, g_out = edge_sums(H, c)
        S_L3v = g5['L3v_sums']['lam%.2f_del%.1f' % (lam, dd)]['S'] if 'lam%.2f_del%.1f' % (lam, dd) in g5['L3v_sums'] else float('nan')
        meas = {
            'kappa': float(kap),
            'sup_f': float(np.max(np.abs(fc))),
            'sup_adev': float(np.max(np.abs(ac))),
            'sup_r_grad_adev': float(np.max(np.abs(s2ad))),
            'sup_omega_eff': float(np.max(np.abs(om_eff))),
            'sup_s2_Wdev': float(np.max(np.abs((1 - tg**2) * Wdev))),
            'closed_vs_series_maxdiff': float(np.nanmax(np.abs(a_closed - ac)[mask])),
            'edge_a': float(a_in + a_out), 'edge_grad': float(g_in + g_out),
            'L3v_sum': S_L3v,
        }
        pr = {}
        # proved envelopes for the two data norms (no measurement used):
        #   |omega_eff| <= lam M + (3/2)|H_1| ,  |(1-t^2)W_dev| <= lam M + 3|H_1| max_t t(1-t^2)
        pr['sup_omega_eff_env'] = lam * M + 1.5 * abs(H[1])
        pr['sup_s2_Wdev_env'] = lam * M + 3 * abs(H[1]) * (2.0 / (3 * np.sqrt(3)))
        pr['adev'] = out['K_f_sup_exact'] * meas['sup_f'] + out['K_omega_sup_exact'] * meas['sup_omega_eff']
        pr['adev_full'] = out['K_f_sup_exact'] * S_L3v + out['K_omega_sup_exact'] * pr['sup_omega_eff_env']
        pr['r_grad_adev'] = 3 * meas['sup_f'] + 3 * pr['adev'] + meas['sup_s2_Wdev']
        pr['r_grad_adev_full'] = 3 * S_L3v + 3 * pr['adev_full'] + pr['sup_s2_Wdev_env']
        pr['a1_offset'] = 0.8 * kap
        pr['r_grad_a1'] = float(g6['sup_r_grad_a1_over_kappaM_slab']) * kap
        pr['collar_ellip_a'] = g6['C_far_a_l3'] * lam + g6['C_in_a_l1'] * lam + kap * 3 * np.log(lam) if lam > 1 else 0.0
        pr['collar_ellip_grad'] = (g6['C_far_grad_l3'] + g6['C_in_grad_l1']) * lam if lam > 1 else 0.0
        Cp_proved = (2 * (pr['a1_offset'] + pr['adev'] + meas['edge_a'] + pr['collar_ellip_a'])
                     + (pr['r_grad_a1'] + pr['r_grad_adev'] + meas['edge_grad'] + pr['collar_ellip_grad'])
                     + lam)
        Cp_meas = (2 * (pr['a1_offset'] + meas['sup_adev'] + meas['edge_a'] + pr['collar_ellip_a'])
                   + (pr['r_grad_a1'] + meas['sup_r_grad_adev'] + meas['edge_grad'] + pr['collar_ellip_grad'])
                   + lam)
        Cp_full = (2 * (pr['a1_offset'] + pr['adev_full'] + meas['edge_a'] + pr['collar_ellip_a'])
                   + (pr['r_grad_a1'] + pr['r_grad_adev_full'] + meas['edge_grad'] + pr['collar_ellip_grad'])
                   + lam)
        rows['lam%.2f_del%.1f' % (lam, dd)] = {'measured': meas, 'proved_pieces': pr,
                                               'C_prime_fully_proved': float(Cp_full),
                                               'C_prime_proved': float(Cp_proved),
                                               'C_prime_measured_bulk': float(Cp_meas)}
        print('%4.1f %.2f  kap %.4f  sup|f| %.4f  |a_dev| m %.4f p %.4f  r|grad a_dev| m %.4f p %.4f'
              '  edge a %.4f g %.4f   C_proved %.2f  C_meas %.2f  clo-ser %.1e'
              % (dd, lam, kap, meas['sup_f'], meas['sup_adev'], pr['adev'],
                 meas['sup_r_grad_adev'], pr['r_grad_adev'], meas['edge_a'], meas['edge_grad'],
                 Cp_proved, Cp_meas, meas['closed_vs_series_maxdiff']), '  C_full %.2f' % Cp_full, flush=True)

out['rows'] = rows
# crossover: C'/(2 a(0)) = C'/(2 kappa M L) < 0.1
out['crossover_L_for_10pct_full'] = {k: float(v['C_prime_fully_proved'] / (2 * v['measured']['kappa'] * 0.1))
                                     for k, v in rows.items()}
out['crossover_L_for_10pct'] = {k: float(v['C_prime_proved'] / (2 * v['measured']['kappa'] * 0.1))
                                for k, v in rows.items()}
out['crossover_L_for_10pct_measured'] = {k: float(v['C_prime_measured_bulk'] / (2 * v['measured']['kappa'] * 0.1))
                                         for k, v in rows.items()}
print(json.dumps({k: out[k] for k in ['c_edge_bangbang_at_10deg', 'c_edge_bangbang_at_30deg',
                                      'K_f_sup_exact', 'K_omega_sup_exact']}, indent=1))
print('c_edge bangbang:', [round(v, 4) for v in out['c_edge_bangbang']])
print('c_edge 15deg   :', [round(v, 4) for v in out['c_edge_taper_15deg']])
with open('g8_results.json', 'w') as f:
    json.dump(out, f, indent=1)

"""r1 - independent re-derivation of the Part A foundations and Theorem A.
Nothing is imported from the target seat; every object is rebuilt here.
Instruments: sympy (exact), scipy QUADPACK adaptive quadrature (H_l, V_lambda),
mpmath 30-digit adaptive quadrature (cross-check at small l).
"""
import json
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.integrate import quad

mp.mp.dps = 30
out = {}

# ---------- 1. C_l^{3/2} = P_{l+1}', N_l, sup norms ----------------------------------
t = sp.symbols('t')
res_C, res_N, res_sup = [], [], []
for l in range(1, 11):
    C = sp.gegenbauer(l, sp.Rational(3, 2), t)
    res_C.append(sp.simplify(sp.expand(C - sp.diff(sp.legendre(l + 1, t), t))))
    N = sp.integrate(sp.expand(C**2 * (1 - t**2)), (t, -1, 1))
    res_N.append(sp.simplify(N - sp.Rational((l + 1) * (l + 2), 1) / (l + sp.Rational(3, 2))))
    res_sup.append(sp.simplify(C.subs(t, 1) - sp.Rational((l + 1) * (l + 2), 2)))
out['gegenbauer_is_Pprime_residuals'] = [str(r) for r in res_C]
out['N_l_residuals'] = [str(r) for r in res_N]
out['supC_at_1_residuals'] = [str(r) for r in res_sup]
assert all(r == 0 for r in res_C + res_N + res_sup), 'A1 identities FAILED'
print('A1 identities: residual 0 for l<=10', flush=True)

tt = np.linspace(-1, 1, 200001)
worst = 0.0
Cm1, Cc = np.ones_like(tt), 3.0 * tt
for l in range(1, 60):
    worst = max(worst, np.max(np.abs(Cc)) / ((l + 1) * (l + 2) / 2.0))
    Cm1, Cc = Cc, ((2 * l + 3) * tt * Cc - (l + 2) * Cm1) / (l + 1)
out['supC_ratio_max_over_l_lt_60'] = float(worst)

# ---------- 2. prefactor inequality, PROVED symbolically ------------------------------
x = sp.symbols('x', positive=True)
diff = sp.expand((x + 1)**3 * (x + 2)**2 - x**3 * (x + sp.Rational(3, 2))**2)
poly = sp.Poly(diff, x)
out['prefactor_poly'] = str(diff)
out['prefactor_all_coeffs_nonneg'] = bool(all(c >= 0 for c in poly.all_coeffs()))
print('prefactor poly', diff, 'coeffs all >= 0:', out['prefactor_all_coeffs_nonneg'], flush=True)

# ---------- 3. the strained tapered profile ------------------------------------------
def phit(th, lam):
    return np.arctan2(np.sin(th) / lam, lam**2 * np.cos(th))

def hdel(p, dl):
    return 1.0 if dl <= 0 else min(1.0, min(p, np.pi - p) / dl)

def theta_delta(lam, dl):
    return float(np.arctan(lam**3 * np.tan(dl))) if dl > 0 else 0.0

def Cl32_scalar(l, ct):
    a, b = 1.0, 3.0 * ct
    for n in range(1, l):
        a, b = b, ((2 * n + 3) * ct * b - (n + 2) * a) / (n + 1)
    return b

def H_quad(l, lam, dl):
    """H_l = -(2 lam/N_l) int_0^{pi/2} h(phit) C_l(cos th) sin^2 th dth  (odd l), QUADPACK."""
    N = (l + 1) * (l + 2) / (l + 1.5)
    thd = theta_delta(lam, dl)
    pts = [0.0, thd, np.pi / 2] if thd > 0 else [0.0, np.pi / 2]
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        v, _ = quad(lambda th: hdel(phit(th, lam), dl) * Cl32_scalar(l, np.cos(th)) * np.sin(th)**2,
                    a, b, limit=400, epsabs=1e-13, epsrel=1e-13)
        tot += v
    return -(2 * lam / N) * tot

# exact rationals, derived symbolically here
sym = {}
for l in [1, 3, 5, 7, 9]:
    C = sp.gegenbauer(l, sp.Rational(3, 2), t)
    N = sp.Rational((l + 1) * (l + 2), 1) / (l + sp.Rational(3, 2))
    I = sp.integrate(C * sp.sqrt(1 - t**2), (t, 0, 1))
    sym[l] = sp.nsimplify(sp.simplify(-2 * I / N))
target = {1: sp.Rational(-5, 6), 3: sp.Rational(3, 40), 5: sp.Rational(-247, 1680),
          7: sp.Rational(1513, 40320), 9: sp.Rational(-2773, 42240)}
out['exact_Hl_symbolic'] = {str(l): str(v) for l, v in sym.items()}
out['exact_Hl_matches_target'] = {str(l): bool(sp.simplify(sym[l] - target[l]) == 0) for l in target}
out['quad_vs_exact_Hl'] = {str(l): float(abs(H_quad(l, 1.0, 0.0) - float(sym[l]))) for l in target}
print('exact H_l match:', out['exact_Hl_matches_target'], flush=True)
print('quad vs exact:', out['quad_vs_exact_Hl'], flush=True)

# mpmath cross-check at a strained tapered case
def H_mp(l, lam, dl):
    N = mp.mpf((l + 1) * (l + 2)) / (l + mp.mpf('1.5'))
    thd = mp.atan(lam**3 * mp.tan(dl)) if dl > 0 else mp.mpf(0)
    pts = [mp.mpf(0), thd, mp.pi / 2] if thd > 0 else [mp.mpf(0), mp.pi / 2]
    def f(th):
        p = mp.atan2(mp.sin(th) / lam, lam**2 * mp.cos(th))
        h = mp.mpf(1) if dl <= 0 else min(mp.mpf(1), min(p, mp.pi - p) / dl)
        a, b = mp.mpf(1), 3 * mp.cos(th)
        for n in range(1, l):
            a, b = b, ((2 * n + 3) * mp.cos(th) * b - (n + 2) * a) / (n + 1)
        return h * b * mp.sin(th)**2
    return -(2 * lam / N) * mp.quad(f, pts)
mpchk = {}
for l in [3, 9, 21, 51]:
    a = H_quad(l, 1.5, np.deg2rad(7.5)); b = float(H_mp(l, mp.mpf('1.5'), mp.mpf('7.5') * mp.pi / 180))
    mpchk[str(l)] = {'quad': a, 'mpmath': b, 'reldiff': abs(a - b) / abs(b)}
out['mpmath_vs_quad_lam1.5_d7.5'] = mpchk
print('mpmath vs quad:', {k: '%.2e' % v['reldiff'] for k, v in mpchk.items()}, flush=True)

# ---------- 4. V_lambda ---------------------------------------------------------------
def Vlam(lam, dl):
    """V = 2 int_0^{pi/2} |Phi'(cos th)| sqrt(sin th) dth ; substitute th = u^2 near 0."""
    thd = theta_delta(lam, dl)
    def dphit_dth(th):
        return lam**-3 / (np.cos(th)**2 + lam**-6 * np.sin(th)**2)
    def g(th):
        p = phit(th, lam)
        hp = (1.0 / dl) if (dl > 0 and p < dl) else 0.0
        return abs(lam * (hp * dphit_dth(th) + hdel(p, dl) * np.cos(th) / np.sin(th))) * np.sqrt(np.sin(th))
    brk = [0.0] + ([np.sqrt(thd)] if thd > 0 else []) + [np.sqrt(np.pi / 2)]
    tot = 0.0
    for a, b in zip(brk[:-1], brk[1:]):
        v, _ = quad(lambda u: g(max(u, 1e-14)**2) * 2 * u, a, b, limit=400, epsabs=1e-12, epsrel=1e-12)
        tot += v
    return 2 * tot

rows = {}
LS = [3, 5, 7, 9, 11, 15, 21, 31, 51, 81, 121, 201, 301]
for dd in [0.0, 7.5, 15.0, 30.0]:
    dl = np.deg2rad(dd)
    for lam in [1.0, 1.1, 1.25, 1.4, 1.5]:
        V = Vlam(lam, dl)
        thd = theta_delta(lam, dl)
        Vb = 2 * lam * (2 + np.sqrt(np.sin(thd)))
        C = np.sqrt(2 / np.pi) * (2 * lam + V)
        worst, wl = 0.0, 0
        for l in LS:
            r = abs(H_quad(l, lam, dl)) * l**1.5 / C
            if r > worst: worst, wl = r, l
        rows['lam%.2f_del%.1f' % (lam, dd)] = {
            'V': float(V), 'V_analytic_bound_2lam(2+sqrt(sin th_d))': float(Vb),
            'V_bound_holds': bool(V <= Vb), 'V_le_6lam': bool(V <= 6 * lam),
            'C_bound': float(C), 'thmA_worst_ratio': float(worst), 'thmA_worst_l': wl,
            'thmA_holds': bool(worst < 1.0), 'theta_delta_deg': float(np.rad2deg(thd))}
        print('%4.1f %.2f V=%.5f (<= %.5f  %s) C=%.5f  thmA worst %.4f at l=%d'
              % (dd, lam, V, Vb, V <= Vb, C, worst, wl), flush=True)
out['rows'] = rows
out['thmA_holds_everywhere'] = bool(all(v['thmA_holds'] for v in rows.values()))
out['V_bound_holds_everywhere'] = bool(all(v['V_bound_holds'] for v in rows.values()))
out['thmA_worst_ratio_overall'] = float(max(v['thmA_worst_ratio'] for v in rows.values()))

# ---------- 5. Corollary A1 pieces ----------------------------------------------------
ls = np.arange(3, 10**7, 1, dtype=float)
out['supC_over_denom'] = float(np.max((ls + 1) * (ls + 2) / (2 * (ls + 4) * (ls - 1))))
out['zeta32_tail_all_l'] = float(mp.zeta(1.5) - 1 - mp.mpf(2)**mp.mpf('-1.5'))
odd = np.arange(3, 20000001, 2, dtype=float)
out['sum_odd_l_ge3'] = float(np.sum(odd**-1.5) + 1.0 / np.sqrt(20000001.0))
out['A1_coefficient_all_l'] = out['supC_over_denom'] * out['zeta32_tail_all_l']
out['A1_coefficient_odd_only'] = out['supC_over_denom'] * out['sum_odd_l_ge3']
out['A1_per_lamM'] = out['A1_coefficient_all_l'] * float(np.sqrt(2 / np.pi)) * 8
out['A1_per_lamM_odd_only'] = out['A1_coefficient_odd_only'] * float(np.sqrt(2 / np.pi)) * 8
print(json.dumps({k: out[k] for k in ['prefactor_all_coeffs_nonneg', 'supC_over_denom',
                                      'zeta32_tail_all_l', 'sum_odd_l_ge3', 'A1_coefficient_all_l',
                                      'A1_coefficient_odd_only', 'A1_per_lamM', 'A1_per_lamM_odd_only',
                                      'thmA_holds_everywhere', 'V_bound_holds_everywhere',
                                      'thmA_worst_ratio_overall', 'supC_ratio_max_over_l_lt_60']}, indent=1))
json.dump(out, open('r1_results.json', 'w'), indent=1)

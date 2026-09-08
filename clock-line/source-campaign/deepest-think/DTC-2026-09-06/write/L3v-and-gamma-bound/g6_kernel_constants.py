"""g6 — the explicit constants of the Gamma bound.

Independent re-derivation (no constant imported from a neighbouring seat):
  * the shell functional  Phi_int, Phi_ext  from the radial Green function of Delta_5,
  * the far / inner Cauchy-Schwarz+Parseval constants for a_dev and for r|grad a_dev|,
  * the collar bound, with the delta-taper repair that removes the 1/sin(phi) blow-up,
  * the exact l = 1 field of a spherical shell and its r|grad a|,
  * the edge-collar (geometric) constants at one octave.
"""
import json
import numpy as np
import sympy as sp

out = {}
t = sp.symbols('t', real=True)

# ---------------------------------------------------------------- 0. the two z-derivative identities
r, zz = sp.symbols('r z_', positive=True)
rho = sp.sqrt(r**2 + zz**2)
res1, res2 = [], []
for l in range(1, 11):
    res1.append(sp.simplify(sp.expand(sp.simplify(
        sp.diff(rho**l * sp.gegenbauer(l, sp.Rational(3, 2), zz / rho), zz)
        - (l + 2) * rho**(l - 1) * sp.gegenbauer(l - 1, sp.Rational(3, 2), zz / rho)))))
for l in range(0, 10):
    res2.append(sp.simplify(sp.expand(sp.simplify(
        sp.diff(rho**(-l - 3) * sp.gegenbauer(l, sp.Rational(3, 2), zz / rho), zz)
        + (l + 1) * rho**(-l - 4) * sp.gegenbauer(l + 1, sp.Rational(3, 2), zz / rho)))))
out['I1_residuals'] = [str(x) for x in res1]
out['I2_residuals'] = [str(x) for x in res2]

# ---------------------------------------------------------------- 1. norms used in the sums
def supC(n):            # ||C_n^{3/2}||_inf = C_n^{3/2}(1)
    return (n + 1) * (n + 2) / 2.0
def supdC(n):           # ||dC_n^{3/2}/dt||_inf = 3 C_{n-1}^{5/2}(1) = 3*binom(n+3,4)... computed
    if n == 0:
        return 0.0
    from math import comb
    return 3.0 * comb(n + 3, 4) / comb(4, 4) if False else 3.0 * comb(n - 1 + 4, n - 1)
# check supdC against sympy
chk = []
tt = np.linspace(-1, 1, 40001)
for n in range(1, 12):
    f = sp.lambdify(t, sp.diff(sp.gegenbauer(n, sp.Rational(3, 2), t), t), 'numpy')
    chk.append(float(np.max(np.abs(f(tt))) / supdC(n)))
out['supdC_check_ratio'] = chk
def N(n):
    return (n + 1) * (n + 2) / (n + 1.5)

LMX = 400
def C_far(lmin, kind):
    """sqrt(2) * int_0^{1/2} Q(u) du/u  for the interior (far-shell) expansion.
       kind='a'   : Q^2 = sum ((l+2)/(2l+3))^2 supC(l-1)^2 u^{2(l-1)} / N_l
       kind='grad': same with supC(l-1) replaced by the r|grad| envelope of the solid harmonic
                    of degree n=l-1:  n*supC(n) + supdC(n)   (r|grad F| <= rho^{n-1}(n||C||+||sC'||))"""
    u = np.linspace(1e-9, 0.5, 200001)
    Q2 = np.zeros_like(u)
    for l in range(lmin, LMX + 1, 2):
        n = l - 1
        amp = supC(n) if kind == 'a' else (n * supC(n) + supdC(n))
        Q2 += ((l + 2) / (2 * l + 3.0))**2 * amp**2 * u**(2 * n) / N(l)
    Q = np.sqrt(Q2)
    return float(np.sqrt(2) * np.trapz(Q / u, u))

def C_in(lmin, kind):
    """sqrt(2) * int_0^{1/2} S(v) dv/v for the exterior (inner-shell) expansion,
       S^2 = sum ((l+1)/(2l+3))^2 amp^2 v^{2(l+4)}/N_l, amp = supC(l+1) (or the grad envelope
       of the exterior harmonic rho^{-(n+3)}C_n with n=l+1: (n+3)supC(n)+supdC(n))."""
    v = np.linspace(1e-9, 0.5, 200001)
    S2 = np.zeros_like(v)
    for l in range(lmin, LMX + 1, 2):
        n = l + 1
        amp = supC(n) if kind == 'a' else ((n + 3) * supC(n) + supdC(n))
        S2 += ((l + 1) / (2 * l + 3.0))**2 * amp**2 * v**(2 * (l + 4)) / N(l)
    S = np.sqrt(S2)
    return float(np.sqrt(2) * np.trapz(S / v, v))

out['C_far_a_l3'] = C_far(3, 'a')      # deviatoric far remainder for a
out['C_far_a_l1'] = C_far(1, 'a')      # cross-check vs far-near-kernel-lemma C1 (z-odd) = 0.291999
out['C_far_grad_l3'] = C_far(3, 'grad')
out['C_in_a_l3'] = C_in(3, 'a')
out['C_in_a_l1'] = C_in(1, 'a')        # cross-check vs C2in (z-odd) = 0.014754
out['C_in_grad_l1'] = C_in(1, 'grad')
out['C_in_grad_l3'] = C_in(3, 'grad')

# ---------------------------------------------------------------- 2. the collar bound
S4 = 8 * np.pi**2 / 3.0                         # |S^4|
RA = (2.0**5 - 2.0**-5)**0.2                    # bathtub radius / rho
out['S4'] = float(S4); out['R_A_over_rho'] = float(RA)
# generic (no taper):   |eta| <= M/r' ,  1/r' <= 2/r on {r'>=r/2};  plus the r'<r/2 piece
out['collar_generic_1_over_s'] = float((3 / (8 * np.pi**2)) * 2 * S4 * RA)      # coefficient of 1/sin(phi)
out['collar_generic_axis_piece'] = float(np.pi / 8)
# tapered:  |eta| <= lam M (pi/2)/(delta rho') <= lam M pi/(delta rho)  on the collar
def collar_taper(delta):
    return float((3 / (8 * np.pi**2)) * (np.pi / delta) * S4 * RA)
out['collar_taper_coeff_over_delta'] = float(collar_taper(1.0))
out['collar_taper_7.5deg'] = collar_taper(np.deg2rad(7.5))
out['collar_taper_15deg'] = collar_taper(np.deg2rad(15.0))
out['collar_taper_30deg'] = collar_taper(np.deg2rad(30.0))
# crossover angle where 4/s = c/delta
out['collar_crossover_sin_phi_7.5deg'] = float(out['collar_generic_1_over_s'] / collar_taper(np.deg2rad(7.5)))

# ---------------------------------------------------------------- 3. exact l=1 field of a spherical shell
R, r0, kap, Mm = sp.symbols('R rho_0 kappa M', positive=True)
H1v = -sp.Rational(5, 3) * kap * Mm                      # kappa = -(3/5) H_1
psi1r = (H1v / 5) * rho * (sp.log(R) - sp.log(rho)) + (H1v / 25) * (rho - r0**5 / rho**4)
psi1 = 3 * (zz / rho) * psi1r
a1k = sp.simplify(-sp.diff(psi1, zz))
out['a1_in_kappa'] = str(a1k)
g = r * sp.sqrt(sp.diff(a1k, r)**2 + sp.diff(a1k, zz)**2)
gg = sp.lambdify((r, zz, R, r0, kap, Mm), g, 'numpy')
aa = sp.lambdify((r, zz, R, r0, kap, Mm), a1k, 'numpy')
Rv, r0v, kapv = 4096.0, 1.0, 0.5
lg = np.linspace(np.log(2 * r0v), np.log(Rv / 2), 400)
ph = np.linspace(1e-4, np.pi - 1e-4, 400)
LG, PH = np.meshgrid(lg, ph)
RH = np.exp(LG)
rv, zv = RH * np.sin(PH), RH * np.cos(PH)
G1 = gg(rv, zv, Rv, r0v, kapv, 1.0)
A1 = aa(rv, zv, Rv, r0v, kapv, 1.0)
out['sup_r_grad_a1_over_kappaM_slab'] = float(np.max(np.abs(G1)) / kapv)
out['sup_a1_minus_kappaMlog_over_kappaM'] = float(np.max(np.abs(A1 - kapv * np.log(Rv / RH))) / kapv)
out['a1_at_origin_check'] = float(kapv * np.log(Rv / r0v))
# a1 on the whole shell (no slab restriction), for comparison
lg2 = np.linspace(np.log(r0v), np.log(Rv), 400)
LG2, PH2 = np.meshgrid(lg2, ph); RH2 = np.exp(LG2)
G2 = gg(RH2 * np.sin(PH2), RH2 * np.cos(PH2), Rv, r0v, kapv, 1.0)
out['sup_r_grad_a1_over_kappaM_fullshell'] = float(np.max(np.abs(G2)) / kapv)

# ---------------------------------------------------------------- 4. edge-collar geometric constants
# the two edge collars are at log-distance >= log 2 from any slab point; their contribution to
# a and to r|grad a| is bounded by the same far/inner sums restricted to u,v <= 1/2, i.e. by
# C_far_*_l1 and C_in_*_l1 above, times the source amplitude.
out['edge_bound_a_l3_plus_inner_l1'] = float(out['C_far_a_l3'] + out['C_in_a_l1'])
out['edge_bound_grad'] = float(C_far(1, 'grad') + out['C_in_grad_l1'])
out['outer_edge_l1_logwidth_3logl_at_1.5'] = float(3*np.log(1.5))
out['C_far_grad_l1'] = float(C_far(1, 'grad'))

print(json.dumps(out, indent=1))
with open('g6_results.json', 'w') as f:
    json.dump(out, f, indent=1)

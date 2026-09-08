"""p4 - numerical verification of LEMMA T' on maps Phi that are NOT Lambda.

Reference value (the quantity step P of prove-lagrangian consumes):
        a_ref  =  (M/2) INT_{rho0}^R P_h(lambda(rho)) dlog rho      ( = (M/2) L lam_bar for h == 1 )
Hypotheses measured, never assumed:
        mu   = sup |Phi x - Lambda x| / |Lambda x|
        muJ  = sup |J_Phi - lambda(rho)^2| / lambda(rho)^2
Bound:  |a[eta_0 o Phi^-1](0) - a_ref| <= pi M A [ 3 mu (1+muJ)/(2(1-mu)^5) + 3 muJ/8 ],  A = L lam_bar.
"""
import json, numpy as np, lib5 as G

out = {}
lam_bar = float(np.sqrt(1.5))
KAPPA_S = 2*(np.sqrt(1.5) - 1)            # = sqrt(6)-2, derived in p3 (sympy exact)

def bound_vs_Lambda(L, mu, muJ_L):
    """(T'_Lambda): the literal form of the task brief -- difference against
    a[eta_0 o Lambda^-1](0), hypotheses stated against J_Lambda."""
    if mu >= 1:
        return float('inf')
    return (np.pi*L*lam_bar*(1 + 2*KAPPA_S/L)
            * (3*mu*(1 + muJ_L)/(2*(1 - mu)**5) + 3*muJ_L/8.0))
Ls = [8.317766166719343, 20.0, 50.0, 100.0, 400.0]

# ---------------------------------------------------------- 0. instrument cross-checks
xc = {}
# (i) a single T_lambda must give (M/2) lambda L  ->  use a constant-lambda stage
def a_single(L, lam, n_phi=400, n_tau=4):
    ph, wph = G.panels(np.array([0.0, np.pi/2, np.pi]), n_phi)
    s, c = np.sin(ph), np.cos(ph)
    g = np.hypot(lam*s, lam**-2*c)
    return 0.75*L*float(np.sum(np.sign(c)*(lam**-2*c/g)/g**4*lam**2*s**2*wph))
xc['a_single_over_ML'] = {str(l): a_single(10.0, l)/10.0 for l in (0.7, 1.0, 1.25, 1.5, 2.0)}
xc['a_single_max_err'] = max(abs(a_single(10.0, l)/10.0 - l/2) for l in (0.7, 1.0, 1.25, 1.5, 2.0))
assert xc['a_single_max_err'] < 1e-12, xc['a_single_max_err']
# (ii) det D Lambda of lib5 against p1's closed form, at random points
rng = np.random.default_rng(20260906)
uu = rng.uniform(0, 30, 500); pp = rng.uniform(1e-3, np.pi - 1e-3, 500)
_, _, JL, lam, D = G.apply_map(uu, pp, 30.0, [('lambda',)])
cl = lam**2*(1 + D*(np.sin(pp)**2 - 2*np.cos(pp)**2))
xc['detDLambda_max_rel'] = float(np.max(np.abs(JL/cl - 1)))
assert xc['detDLambda_max_rel'] < 1e-15
# (iii) Step-3 shell integral:  INT |eta_0| lambda^2 |Lambda x|^-4 dx_5 = pi^3 M INT lambda dlog rho
def step3(L, n_phi=300, n_tau=300):
    ph, wph = G.panels(np.array([0.0, np.pi/2, np.pi]), n_phi)
    ta, wta = G.gl(0.0, 1.0, n_tau)
    PH, TA = np.meshgrid(ph, ta, indexing='ij'); WP, WT = np.meshgrid(wph, wta, indexing='ij')
    lam = G.lam_sigma(TA)
    g = np.hypot(lam*np.sin(PH), lam**-2*np.cos(PH))
    return 2*np.pi**2*L*float(np.sum(lam**2*np.sin(PH)**2/g**4*WP*WT))
xc['step3'] = []
for L in Ls:
    lhs = step3(L); rhs = np.pi**3*L*lam_bar
    xc['step3'].append(dict(L=L, lhs=lhs, rhs=rhs, rel=abs(lhs - rhs)/rhs))
xc['step3_max_rel'] = max(d['rel'] for d in xc['step3'])
assert xc['step3_max_rel'] < 1e-13, xc['step3_max_rel']
# (iv) taper cross-check of the instrument at lambda = 1
def a_taper1(L, deg=7.5, n_phi=400):
    d = np.deg2rad(deg)
    ph, wph = G.panels(np.array(sorted({0.0, d, np.pi/2, np.pi - d, np.pi})), n_phi)
    c = np.cos(ph)
    return 0.75*L*float(np.sum(np.sign(c)*G.h_taper(ph, deg)*c*np.sin(ph)**2*wph))
xc['kappa_delta_7p5'] = a_taper1(10.0)/10.0
out['cross_checks'] = xc

# ---------------------------------------------------------- 1. Phi = Lambda  (mu = 0)
rows = []
for L in Ls:
    a_true = G.a_of_map(L, [('lambda',)])
    a_ref = 0.5*L*lam_bar
    mu, muJ, muJvsJL = G.hyp_constants(L, [('lambda',)])
    b = G.lemmaT_prime_bound(L, mu, muJ, lam_bar)
    rows.append(dict(case='Phi = Lambda', L=L, mu=mu, muJ=muJ, muJ_vs_JLambda=muJvsJL,
                     a_ref=a_ref, a_phi=a_true, delta=a_true - a_ref, bound=b,
                     holds=bool(b >= abs(a_true - a_ref)),
                     slack=b/abs(a_true - a_ref), rel_err=abs(a_true - a_ref)/a_ref,
                     rel_err_times_L=abs(a_true - a_ref)/a_ref*L,
                     bound_rel=b/a_ref, bound_rel_times_L=b/a_ref*L))
out['lambda_rows'] = rows

# ---------------------------------------------------------- 2. Phi != Lambda
cases = []
for L in (8.317766166719343, 50.0):
    for name, spec in [
        ('Lambda + ripple(mu_a=0.02,k=1)',  [('lambda',), ('ripple', 0.02, 1.0)]),
        ('Lambda + ripple(mu_a=0.05,k=3)',  [('lambda',), ('ripple', 0.05, 3.0)]),
        ('Lambda + shear(beta=0.02)',       [('lambda',), ('shear', 0.02)]),
        ('Lambda + shear(beta=0.08)',       [('lambda',), ('shear', 0.08)]),
        ('Lambda + ripple + shear',         [('lambda',), ('ripple', 0.03, 2.0), ('shear', 0.04)]),
        ('shear + Lambda (pre-composed)',   [('shear', 0.05), ('lambda',)]),
    ]:
        a_phi = G.a_of_map(L, spec)
        a_ref = 0.5*L*lam_bar
        a_Lam = G.a_of_map(L, [('lambda',)])
        mu, muJ, muJvsJL = G.hyp_constants(L, spec)
        b = G.lemmaT_prime_bound(L, mu, muJ, lam_bar)
        d = a_phi - a_ref
        bL = bound_vs_Lambda(L, mu, muJvsJL)
        dL = a_phi - a_Lam
        cases.append(dict(case=name, L=L, mu=mu, muJ=muJ, muJ_vs_JLambda=muJvsJL,
                          a_ref=a_ref, a_phi=a_phi, a_Lambda=a_Lam,
                          delta=d, bound=b, holds=bool(b >= abs(d)),
                          slack=(b/abs(d) if d != 0 else float('inf')),
                          delta_Lambda=dL, bound_Lambda=bL,
                          holds_Lambda=bool(bL >= abs(dL)),
                          slack_Lambda=(bL/abs(dL) if dL != 0 else float('inf'))))
out['phi_rows'] = cases
out['all_hold_Lambda'] = all(r['holds_Lambda'] for r in cases)
out['worst_slack_Lambda'] = min(r['slack_Lambda'] for r in cases)
assert out['all_hold_Lambda']
out['all_hold'] = all(r['holds'] for r in rows + cases)
out['worst_slack'] = min(r['slack'] for r in rows + cases)
assert out['all_hold'], [r for r in rows + cases if not r['holds']]

# ---------------------------------------------------------- 3. tapered datum
tap = []
for L in (8.317766166719343, 50.0):
    for name, spec in [('Phi = Lambda', [('lambda',)]),
                       ('Lambda + shear(beta=0.05)', [('lambda',), ('shear', 0.05)])]:
        a_phi = G.a_of_map(L, spec, h=lambda p: G.h_taper(p, 7.5), taper_deg=7.5)
        # a_ref for the taper: (M/2) INT P_h(lambda) dlog rho, computed with the SAME instrument
        ph, wph = G.panels(np.array(sorted({0.0, np.deg2rad(7.5), np.pi/2,
                                            np.pi - np.deg2rad(7.5), np.pi})), 400)
        ta, wta = G.gl(0.0, 1.0, 300)
        PH, TA = np.meshgrid(ph, ta, indexing='ij'); WP, WT = np.meshgrid(wph, wta, indexing='ij')
        lam = G.lam_sigma(TA)
        g = np.hypot(lam*np.sin(PH), lam**-2*np.cos(PH))
        a_ref = 0.75*L*float(np.sum(np.sign(np.cos(PH))*G.h_taper(PH, 7.5)
                                    * (lam**-2*np.cos(PH)/g)/g**4*lam**2*np.sin(PH)**2*WP*WT))
        mu, muJ, _ = G.hyp_constants(L, spec)
        b = G.lemmaT_prime_bound(L, mu, muJ, lam_bar)   # A uses lam_bar (h <= 1 -> same majorant)
        d = a_phi - a_ref
        tap.append(dict(case=name, L=L, mu=mu, muJ=muJ, a_ref=a_ref, a_phi=a_phi, delta=d,
                        bound=b, holds=bool(b >= abs(d)),
                        slack=(b/abs(d) if d != 0 else float('inf'))))
out['taper_rows'] = tap
assert all(r['holds'] for r in tap)

# ---------------------------------------------------------- 4. quadrature convergence
conv = []
for L in (8.317766166719343, 400.0):
    v = [G.a_of_map(L, [('lambda',), ('shear', 0.04)], n_phi=n, n_tau=n) for n in (100, 200, 400)]
    conv.append(dict(L=L, vals=v, rel_200_400=abs(v[2] - v[1])/abs(v[2])))
out['convergence'] = conv
out['conv_max'] = max(c['rel_200_400'] for c in conv)
assert out['conv_max'] < 1e-10, out['conv_max']

json.dump(out, open('p4_results.json', 'w'), indent=1)
print("cross-checks: a[T_lam](0)/(ML) max err vs lam/2 =", xc['a_single_max_err'])
print("              det D Lambda vs p1 closed form, max rel =", xc['detDLambda_max_rel'])
print("              STEP 3  INT|eta_0|lam^2|Lam x|^-4 = pi^3 M INT lam dlog rho, max rel =",
      xc['step3_max_rel'])
print("              kappa_delta(7.5deg) =", repr(xc['kappa_delta_7p5']))
print("\nPhi = Lambda   (mu = 0 exactly; only the Jacobian shear enters)")
print(f"{'L':>9s} {'mu':>10s} {'muJ':>10s} {'muJ*L':>10s} {'|Delta|/a_ref*L':>16s} {'bound/a_ref*L':>14s} {'slack':>8s}")
for r in rows:
    print(f"{r['L']:9.3f} {r['mu']:10.2e} {r['muJ']:10.6f} {r['muJ']*r['L']:10.7f} "
          f"{r['rel_err_times_L']:16.6f} {r['bound_rel_times_L']:14.6f} {r['slack']:8.3f}")
print("\nPhi != Lambda")
print(f"{'L':>7s}  {'case':34s} {'mu':>9s} {'muJ':>9s} {'|Delta|':>10s} {'bound':>11s} {'slack':>8s} holds")
for r in cases:
    print(f"{r['L']:7.2f}  {r['case']:34s} {r['mu']:9.5f} {r['muJ']:9.5f} {abs(r['delta']):10.5f} "
          f"{r['bound']:11.4f} {r['slack']:8.2f} {r['holds']}")
print("\n(T'_Lambda): the literal brief form, difference against a[eta_0 o Lambda^-1](0)")
print(f"{'L':>7s}  {'case':34s} {'mu':>9s} {'muJ vs JLam':>12s} {'|Delta_Lam|':>12s} "
      f"{'bound':>11s} {'slack':>9s} holds")
for r in cases:
    print(f"{r['L']:7.2f}  {r['case']:34s} {r['mu']:9.5f} {r['muJ_vs_JLambda']:12.5f} "
          f"{abs(r['delta_Lambda']):12.5f} {r['bound_Lambda']:11.4f} {r['slack_Lambda']:9.2f} "
          f"{r['holds_Lambda']}")
print("   worst slack:", out['worst_slack_Lambda'])

print("\ntapered datum (7.5 deg)")
for r in tap:
    print(f"{r['L']:7.2f}  {r['case']:28s} mu={r['mu']:.5f} muJ={r['muJ']:.5f} "
          f"|Delta|={abs(r['delta']):.5f} bound={r['bound']:.4f} slack={r['slack']:.2f} {r['holds']}")
print("\nquadrature convergence (200 vs 400 nodes):", [(c['L'], f"{c['rel_200_400']:.2e}") for c in conv])
print("ALL ROWS HOLD:", out['all_hold'], "  worst slack:", out['worst_slack'])

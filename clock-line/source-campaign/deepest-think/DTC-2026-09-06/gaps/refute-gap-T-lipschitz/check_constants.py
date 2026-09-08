"""check_constants.py - re-asserts every number displayed in my NOTE.md against my result JSONs."""
import json, numpy as np
r1=json.load(open('r1_results.json')); r2=json.load(open('r2_results.json'))
r3=json.load(open('r3_results.json')); r4=json.load(open('r4_results.json'))
def ok(n,g,w,t):
    assert abs(g-w)<=t, f"{n}: {g} vs {w} (tol {t})"; print(f"  OK  {n:52s} {g!r}")

print("r1 - Lemma T's own algebra, re-derived independently")
for k in ('gradK_norm_residual','step_a_residual','step_b_residual','step_c_residual',
          'step_d_residual','antiderivative_residual','per_shell_core_residual','a_integrand_residual'):
    assert r1[k]=='0', (k, r1[k]); print(f"  OK  {k:52s} 0")
assert r1['J_of_lambda']=='pi/(2*lambda)' and r1['per_shell_core']=='lambda/3'
assert r1['gradK_sup_symbolic']=='3/(2*pi**2)'
ok('grad K sup constant', r1['gradK_sup_num'], 3/(2*np.pi**2), 1e-15)
ok('K sup constant',      r1['K_sup_num'],     3/(8*np.pi**2), 1e-15)
ok('angular identity, max mpmath error', r1['angular_identity_max_err'], 0.0, 1e-15)
ok('det D Lambda, max FD relative error', r1['detDLambda_max_rel_err_fd'], 0.0, 1e-7)

print("r2 - the bridge")
ok('kappa*theta terminal', r2['kappa_theta_terminal'], 0.36700683814454793, 1e-12)
ok('log-mean of lambda over the shell', r2['lam_logmean'], np.sqrt(1.5), 1e-10)
mus  = [r['mu_single']  for r in r2['rows']]
ok('mu(H2) vs best single T_lambda, min over L', min(mus), 0.3853290, 1e-6)
ok('mu(H2) vs best single T_lambda, max over L', max(mus), 0.3853290, 1e-6)
ok('spread of mu(H2) across L = 8.3 .. 400',     max(mus)-min(mus), 0.0, 1e-12)
ok('lambdabar*', r2['rows'][0]['lambar_star'], 1.177, 1e-9)
ok('muJ(H3) vs single, min over L', min(r['muJ_single'] for r in r2['rows']), 0.6278128, 1e-6)
ok('Lemma-T relative bound, single T_lambda, min over L',
   min(r['bound_single_rel'] for r in r2['rows']), 66.1116, 1e-3)
for r in r2['rows']:
    assert r['bound_single_rel'] > 60, r
    ok(f"muJ_local * L  (L={r['L']:g})", r['muJ_local']*r['L'], 0.8990, 2e-3)
    ok(f"Lemma-T rel bound, shell-local ref, * L  (L={r['L']:g})", r['bound_local_rel']*r['L'], 2.08, 0.06)
    ok(f"true rel error of the shell-resolved ref * L  (L={r['L']:g})",
       r['rel_err_shell_ref']*r['L'], 0.343, 0.01)
for r in r2['rows']:
    assert 0.17 < r['rel_err_lam1'] < 0.23 and 0.17 < r['rel_err_lam15'] < 0.23, r
ok('worst single-lambda-at-endpoint relative error',
   max(max(r['rel_err_lam1'], r['rel_err_lam15']) for r in r2['rows']), 0.2237, 1e-3)
for r in r2['repair']:
    ok(f"per-shell repair, sum(bound)/a * L (L={r['L']:g})", r['rel']*r['L'], 7.4, 1.6)
ok('quadrature convergence (200 vs 400 nodes)',
   max(c['rel_200_400'] for c in r2['convergence']), 0.0, 1e-12)
x = r2['cross_check']
ok('a[T_1](0)/(M L)      (estate: 1/2)',   x['a_single_lam1_over_ML'][0], 0.5, 1e-12)
ok('a[T_1.5](0)/(M L)    (estate: 3/4)',   x['a_single_lam_15_over_ML'], 0.75, 1e-12)
ok('kappa_delta(7.5 deg) (estate: 0.4997212)', x['kappa_delta'], 0.4997212305210886, 1e-12)

print("r3 - the shear between shells is angular, not radial")
ok('max polar-angle shift, deg', r3['max_polar_angle_shift_deg'], 17.16108, 1e-4)
ok('lambda^3 at sigma=0 (tan-ratio)', r3['lambda_cubed_ratio_extremes'][0], 3.375, 1e-9)
ok('sup angular part / |T_lambdabar x|', r3['sup_rel_angular_part'], 0.3175481, 1e-6)
ok('max angular fraction of the displacement', r3['max_angular_fraction'], 1.0, 1e-9)
ok('log-periodic null direction, diffeo rows', r3['max_rel_diffeo'], 0.0, 1e-13)

print("r4 - the seat's control IS a diffeomorphism")
assert r4['all_monotone'] and r4['all_onto'] and not r4['any_claim_violation']
for r in r4['rows']:
    ok(f"min(1+beta') at mu_b={r['mu_beta']}, phi_c={r['phi_c']:g}",
       r['min_1_plus_dbeta'], 1-r['mu_beta'], 1e-9)
print("\nALL ASSERTIONS PASSED")

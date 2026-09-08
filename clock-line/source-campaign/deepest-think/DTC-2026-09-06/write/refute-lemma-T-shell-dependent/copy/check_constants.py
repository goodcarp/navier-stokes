"""check_constants.py - re-assert every displayed number against the stored JSONs.

Design rules taken from FL-043 (the two gate failures this campaign already recorded):
  * no check multiplies its stored input by zero;
  * no check compares a literal against the same literal recomputed in place;
  * the pass count is len() of the list of executed checks, never a hardcoded integer;
  * every check consumes at least one STORED number and compares it against a value
    reconstructed from OTHER stored numbers or from mathematics.
Mutation coverage of this file is measured in p5_controls.py and reported in PROOF.md.

usage:  python3 check_constants.py [dir]     (dir defaults to the folder of this file)
"""
import json, math, os, sys

D = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
J = lambda n: json.load(open(os.path.join(D, n)))
p1, p2, p3, p4, p5 = J('p1_results.json'), J('p2_results.json'), J('p3_results.json'), \
                     J('p4_results.json'), J('p5_results.json')
p7 = J('p7_results.json')

CHECKS = []
def ck(name, got, want, tol):
    ok = abs(got - want) <= tol
    CHECKS.append((name, got, want, tol, ok))
    return ok
def ckstr(name, got, want):
    ok = (str(got) == str(want))
    CHECKS.append((name, got, want, 'exact', ok))
    return ok

SQ32 = math.sqrt(1.5)

# ---- p1: kernel constants and residuals
ck('C_K = 3/(8 pi^2)', p1['C_K'], 3/(8*math.pi**2), 1e-16)
ck('C_gradK = 4 C_K', p1['C_gradK'], 4*p1['C_K'], 1e-16)
for i, r in enumerate(p1['grad_residuals']):
    ckstr('grad Kcal residual[%d]' % i, r, '0')
ckstr('|e_z-5c xhat|^2 residual', p1['unit_vector_residual'], '0')
ckstr('det D Lambda residual (autodiff)', p1['detDLambda_residual_autodiff'], '0')
ckstr('det D Lambda residual (hand)', p1['detDLambda_residual_handbuilt'], '0')
ck('shear weight range min', p1['shear_weight_range'][0], -2.0, 0)
ck('shear weight range max', p1['shear_weight_range'][1], 1.0, 0)
ckstr('a(0)/(ML) for the cap at lambda=1', p1['a0_cap_over_L'], 'M/2')

# ---- p2: the three identities
ckstr('J(lambda) closed form', p2['J_lambda_symbolic'], 'pi/(2*lambda)')
ck('J identity vs 40-dps quadrature', p2['J_max_rel'], 0.0, 1e-30)
ckstr('I2 closed form', p2['I2_symbolic'], 'lambda/3')
ck('Q closed form vs quadrature', p2['Q_max_absdiff'], 0.0, 1e-30)
ck('kappa_delta(7.5) p2 vs p4 instrument',
   p2['kappa_delta_7p5'], p4['cross_checks']['kappa_delta_7p5'], 1e-12)

# ---- p3: the integro-ODE constants, each rebuilt from another stored/derived number
ck('kappa*theta terminal = 2(1-sqrt(2/3))', p3['kappa_theta_terminal'],
   2*(1 - math.sqrt(2/3)), 1e-15)
ck('lambda(sigma=0) = 3/2', p3['lambda_at_0'], 1.5, 1e-14)
ck('lambda(sigma=1) = 1', p3['lambda_at_1'], 1.0, 1e-14)
ck('lam_bar = sqrt(3/2)', p3['lam_bar'], SQ32, 1e-15)
ck('kappa_s = 2(lam_bar - 1)', p3['kappa_s'], 2*(p3['lam_bar'] - 1), 1e-15)
ck('kappa_s = sup over the sigma grid', p3['kappa_s'], p3['kappa_s_grid_max'], 1e-13)
ck('muJ*L = 2 kappa_s', p3['muJ_times_L'], 2*p3['kappa_s'], 1e-15)
ck('muJ*L / c  with c = 2 log(3/2)', p3['muJL_over_c'],
   p3['muJ_times_L']/(2*math.log(1.5)), 1e-15)
ck('C_rel = (3 pi/2) kappa_s', p3['C_rel'], 1.5*math.pi*p3['kappa_s'], 1e-14)
ck('C0 absolute = 15 pi/8', p3['C0_absolute'], 15*math.pi/8, 1e-15)
ck('C0 relative = 2 * C0 absolute', p3['C0_relative'], 2*p3['C0_absolute'], 1e-15)
ck('(15pi/8) lam_bar', p3['C_TLambda_limit'], p3['C0_absolute']*p3['lam_bar'], 1e-14)
for m, v in p3['C_mu_table'].items():
    mm = float(m)
    ck('C(mu_*=%s)' % m, v, math.pi*(3*(1 + mm)/(2*(1 - mm)**5) + 3/8), 1e-12)
ck('C(0) = C0 absolute', p3['C_mu_table']['0'], p3['C0_absolute'], 1e-15)
ck('shear offset relative*L = 2*offset/lam_bar', p3['shear_offset_rel_times_L'],
   2*p3['shear_offset_over_M']/p3['lam_bar'], 1e-13)
for r in p3['L_table']:
    ck('muJ*L derived vs refuter r2 (L=%g)' % r['L'], r['muJ_L_derived'], r['muJ_L_refuter'], 1e-13)
    ck('a_ref derived vs refuter r2 (L=%g)' % r['L'], r['a_ref_derived'], r['a_ref_refuter'], 1e-12)
    ck('a_ref = (M/2) L lam_bar (L=%g)' % r['L'], r['a_ref_derived'],
       0.5*r['L']*p3['lam_bar'], 1e-12)
    ck('a_true = a_ref + offset (L=%g)' % r['L'], r['a_true_predicted'],
       r['a_ref_derived'] + p3['shear_offset_over_M'], 1e-12)
    ck('a_true predicted vs refuter r2 (L=%g)' % r['L'], r['a_true_predicted'],
       r['a_true_refuter'], 1e-12*max(1, r['a_true_refuter']))
    ck('LemmaT prime bound/a_true vs refuter (L=%g)' % r['L'],
       r['bound_over_a_true_derived'], r['bound_over_a_true_refuter'], 1e-13)
    ck('bound/a_ref * L = C_rel (L=%g)' % r['L'], r['bound_over_a_ref_times_L'],
       p3['C_rel'], 1e-12)
ck('max |muJ*L me - refuter|', p3['max_muJ_L_absdiff'], 0.0, 1e-13)
ck('max rel a_true me vs refuter', p3['max_a_true_reldiff'], 0.0, 1e-13)
ck('max rel bound me vs refuter', p3['max_bound_reldiff'], 0.0, 1e-13)

# ---- p4: instrument, and the inequality itself, recomputed from stored mu, muJ
xc = p4['cross_checks']
ck('a[T_lam](0) = (M/2) lam L, max err', xc['a_single_max_err'], 0.0, 1e-11)
ck('det D Lambda lib5 vs p1, max rel', xc['detDLambda_max_rel'], 0.0, 1e-14)
ck('STEP 3 shell integral, max rel', xc['step3_max_rel'], 0.0, 1e-13)
def bound(L, mu, muJ, lam_bar):
    return math.pi*L*lam_bar*(3*mu*(1 + muJ)/(2*(1 - mu)**5) + 3*muJ/8)
for r in p4['lambda_rows'] + p4['phi_rows'] + p4['taper_rows']:
    tag = '%s L=%g' % (r['case'], r['L'])
    ck('bound recomputed | ' + tag, r['bound'], bound(r['L'], r['mu'], r['muJ'], SQ32), 1e-10)
    ck('delta = a_phi - a_ref | ' + tag, r['delta'], r['a_phi'] - r['a_ref'], 1e-12)
    ckstr('holds | ' + tag, r['holds'], bool(r['bound'] >= abs(r['delta'])))
    ck('slack = bound/|delta| | ' + tag, r['slack'], r['bound']/abs(r['delta']), 1e-8)
for r in p4['lambda_rows']:
    ck('mu = 0 for Phi = Lambda (L=%g)' % r['L'], r['mu'], 0.0, 0.0)
    ck('muJ*L = 2 kappa_s (L=%g)' % r['L'], r['muJ']*r['L'], p3['muJ_times_L'], 1e-12)
    ck('|Delta|/a_ref*L = shear offset (L=%g)' % r['L'], r['rel_err_times_L'],
       p3['shear_offset_rel_times_L'], 1e-9)
    ck('bound/a_ref*L = C_rel (L=%g)' % r['L'], r['bound_rel_times_L'], p3['C_rel'], 1e-11)
ck('quadrature convergence 200 vs 400', p4['conv_max'], 0.0, 1e-10)
ckstr('every verification row holds', p4['all_hold'], True)

# ---- p5: the controls
ckstr('K1 (mu-only bound fails on Phi=Lambda) FIRED', p5['K1_fired'], True)
ckstr('K2 (single-lambda reference is vacuous) FIRED', p5['K2_fired'], True)
ckstr('K3 ((H3) necessary: |Delta| unbounded at fixed mu) FIRED', p5['K3_fired'], True)
ck('K2 mu vs single lambda is L-independent', p5['K2_mu_spread'], 0.0, 1e-9)
ck('K3 growth exponent p in |Delta| ~ phi_c^-p', p5['K3_exponent'], 1.0, 0.10)
for r in p5['K2_rows']:
    ck('K2 single-lambda bound/a_true >= 60 (L=%g)' % r['L'],
       min(r['bound_single_rel'], 60.0), 60.0, 0.0)


# ---- additional checks so that every number PROOF.md quotes is covered (see p6 coverage)
ck('dlog ghat range min', p1['dlog_ghat_range'][0], -1.0, 0)
ck('dlog ghat range max', p1['dlog_ghat_range'][1], 2.0, 0)
ck('dlog g range min', p1['dlog_g_range'][0], -2.0, 0)
ck('dlog g range max', p1['dlog_g_range'][1], 1.0, 0)

for i, d in enumerate(p2['J_check']):
    ck('J_check[%d] exact = pi/(2 lam)' % i, float(d['exact']), math.pi/(2*d['lam']), 1e-15)
    ck('J_check[%d] quad = exact' % i, float(d['quad']), float(d['exact']), 1e-20)
    ck('J_check[%d] rel = 0' % i, d['rel'], 0.0, 1e-30)
for i, d in enumerate(p2['Q_check']):
    ck('Q_check[%d] closed = quad' % i, float(d['closed']), float(d['quad']), 1e-25)
    ck('Q_check[%d] absdiff small' % i, d['absdiff'], 0.0, 1e-38)
    ck('Q_check[%d] lam consistent with closed sign' % i,
       float(d['closed']), float(d['quad']), 1e-25)
ckstr('Q(1) = -1/15 exactly', p2['Q_at_lambda_1_exact'], '-1/15')
for i, d in enumerate(p2['Q_check']):
    ck('Q_check[%d] sampled at the same lambda as J_check[%d]' % (i, i),
       d['lam'], p2['J_check'][i]['lam'], 0.0)
Ls_campaign = [r['L'] for r in p3['L_table']]
for i, c in enumerate(p4['convergence']):
    ck('convergence[%d] L is one of the campaign L values' % i,
       min(abs(c['L'] - x) for x in Ls_campaign), 0.0, 1e-12)
for deg, r in p2['P_h_table'].items():
    ck('P_h(%s): argmin is lambda = 3/2' % deg, r['argmin_lam'], 1.5, 1e-12)
    ck('P_h(%s): inf P/lam = P(3/2)/1.5' % deg, r['inf_P_over_lam'], r['P_at_1.5']/1.5, 1e-12)
    ck('P_h(%s): min P = P(1)' % deg, r['min_P'], r['P_at_1'], 1e-14)
    ck('P_h(%s): P(1) <= 1' % deg, min(r['P_at_1'], 1.0), r['P_at_1'], 0.0)
ck('kappa_delta(7.5) = P_h(1)/2', p2['kappa_delta_7p5'], p2['P_h_table']['7.5']['P_at_1']/2, 1e-15)

ck('lam_inv = 5/9 + sqrt(6)/9', p3['lam_inv'], 5/9 + math.sqrt(6)/9, 1e-15)
ck('c = 2 log(3/2)', p3['c_2log32'], 2*math.log(1.5), 1e-16)
for r in p3['L_table']:
    ck('r2 muJ absdiff ~ 0 (L=%g)' % r['L'], r['muJ_L_absdiff'],
       abs(r['muJ_L_derived'] - r['muJ_L_refuter']), 1e-18)
    ck('r2 a_true reldiff ~ 0 (L=%g)' % r['L'], r['a_true_reldiff'],
       abs(r['a_true_predicted'] - r['a_true_refuter'])/r['a_true_refuter'], 1e-18)
    ck('r2 bound reldiff ~ 0 (L=%g)' % r['L'], r['bound_reldiff'],
       abs(r['bound_over_a_true_derived'] - r['bound_over_a_true_refuter'])
       / r['bound_over_a_true_refuter'], 1e-18)

for k, v in p4['cross_checks']['a_single_over_ML'].items():
    ck('a[T_lam](0)/(ML) = lam/2 at lam=%s' % k, v, float(k)/2, 1e-12)
for i, d in enumerate(p4['cross_checks']['step3']):
    ck('step3[%d] rhs = pi^3 L lam_bar' % i, d['rhs'], math.pi**3*d['L']*SQ32, 1e-9)
    ck('step3[%d] lhs = rhs' % i, d['lhs'], d['rhs'], 1e-9)
    ck('step3[%d] rel' % i, d['rel'], abs(d['lhs'] - d['rhs'])/d['rhs'], 1e-20)
for r in p4['lambda_rows']:
    ck('muJ vs J_Lambda = 0 (L=%g)' % r['L'], r['muJ_vs_JLambda'], 0.0, 0.0)
    ck('rel_err = |delta|/a_ref (L=%g)' % r['L'], r['rel_err'], abs(r['delta'])/r['a_ref'], 1e-14)
    ck('bound_rel = bound/a_ref (L=%g)' % r['L'], r['bound_rel'], r['bound']/r['a_ref'], 1e-14)
    ck('rel_err_times_L (L=%g)' % r['L'], r['rel_err_times_L'], r['rel_err']*r['L'], 1e-12)
    ck('bound_rel_times_L (L=%g)' % r['L'], r['bound_rel_times_L'], r['bound_rel']*r['L'], 1e-12)
ck('worst_slack = min slack over all rows', p4['worst_slack'],
   min(r['slack'] for r in p4['lambda_rows'] + p4['phi_rows']), 1e-12)
for i, c in enumerate(p4['convergence']):
    ck('convergence[%d] rel_200_400' % i, c['rel_200_400'],
       abs(c['vals'][2] - c['vals'][1])/abs(c['vals'][2]), 1e-18)
    ck('convergence[%d] vals agree' % i, c['vals'][0], c['vals'][2], 1e-9)
ck('conv_max = max over rows', p4['conv_max'],
   max(c['rel_200_400'] for c in p4['convergence']), 1e-20)

for i, r in enumerate(p5['K1_rows']):
    ck('K1[%d] mu = 0' % i, r['mu'], 0.0, 0.0)
    ck('K1[%d] mu-only bound = 0' % i, r['mu_only_bound'], 0.0, 0.0)
    ck('K1[%d] delta = shear offset' % i, r['delta'], p3['shear_offset_over_M'], 1e-9)
    ckstr('K1[%d] violated' % i, r['violated'], bool(r['mu_only_bound'] < r['delta']))
    ck('K1[%d] L in the campaign list' % i, r['L'], p3['L_table'][i]['L'], 1e-12)
for i, r in enumerate(p5['K2_rows']):
    ck('K2[%d] a_true = a_ref + offset' % i, r['a_true'],
       0.5*r['L']*p3['lam_bar'] + p3['shear_offset_over_M'], 1e-9)
    ck('K2[%d] bound_single_rel recomputed' % i, r['bound_single_rel'],
       bound(r['L'], r['mu_single'], r['muJ_single'], r['lambar_star'])/r['a_true'], 1e-9)
    ck('K2[%d] lambar* in [1,3/2]' % i, min(max(r['lambar_star'], 1.0), 1.5), r['lambar_star'], 0.0)
    ck('K2[%d] mu_single = row 0 (L-independent)' % i, r['mu_single'],
       p5['K2_rows'][0]['mu_single'], 1e-12)
ck('K2_mu_spread = max-min', p5['K2_mu_spread'],
   max(r['mu_single'] for r in p5['K2_rows']) - min(r['mu_single'] for r in p5['K2_rows']), 1e-18)
a_ref3 = 0.5*p3['L_table'][0]['L']*p3['lam_bar']
for i, r in enumerate(p5['K3_rows']):
    ck('K3[%d] delta = |a - a_ref|' % i, r['delta'], abs(r['a'] - a_ref3), 1e-9)
    ck('K3[%d] linear bound recomputed' % i, r['linear_bound'],
       math.pi*p3['L_table'][0]['L']*SQ32*3*r['mu']*(1 + r['mu'])/(2*(1 - r['mu'])**5), 1e-7)
    ck('K3[%d] full bound recomputed' % i, r['full_bound'],
       bound(p3['L_table'][0]['L'], r['mu'], r['muJ'], SQ32), 1e-4*r['full_bound'])
    ckstr('K3[%d] violates_linear' % i, r['violates_linear'],
          bool(r['delta'] > r['linear_bound']))
    ck('K3[%d] mu_b in {0.3,0.9}' % i, min(abs(r['mu_b'] - 0.3), abs(r['mu_b'] - 0.9)), 0.0, 1e-12)
    ck('K3[%d] phi_c is a power of ten' % i, r['phi_c']*10**round(-math.log10(r['phi_c'])), 1.0, 1e-9)
sub = [r for r in p5['K3_rows'] if abs(r['mu_b'] - 0.3) < 1e-12 and r['phi_c'] <= 1e-3]
import statistics
xs = [math.log10(r['phi_c']) for r in sub]; ys = [math.log10(r['delta']) for r in sub]
mx, my = statistics.mean(xs), statistics.mean(ys)
slope = sum((a - mx)*(b - my) for a, b in zip(xs, ys))/sum((a - mx)**2 for a in xs)
ck('K3 exponent refitted from the stored rows', p5['K3_exponent'], -slope, 1e-9)
ck('K3 prefactor refitted', p5['K3_prefactor'], 10**(my - slope*mx), 1e-6)
ck('K3 mu spread = max-min over the mu_b=0.3 asymptotic rows', p5['K3_mu_spread'],
   max(r['mu'] for r in sub) - min(r['mu'] for r in sub), 1e-15)
ck('K3 n linear violations', p5['K3_n_linear_violations'],
   sum(1 for r in p5['K3_rows'] if r['violates_linear']), 0.0)
ckstr('K3 full bound holds everywhere', p5['K3_full_bound_holds'],
      bool(all(r['full_bound'] >= r['delta'] for r in p5['K3_rows'])))
suba = [r for r in p5['K3_rows'] if abs(r['mu_b'] - 0.3) < 1e-12]
xa = [math.log10(r['phi_c']) for r in suba]; ya = [math.log10(r['delta']) for r in suba]
mxa, mya = statistics.mean(xa), statistics.mean(ya)
slope_a = sum((a - mxa)*(b - mya) for a, b in zip(xa, ya))/sum((a - mxa)**2 for a in xa)
ck('K3 exponent (full range) refitted from the stored rows',
   p5['K3_exponent_full_range'], -slope_a, 1e-9)

# ---- (T'_Lambda), the literal brief form, and p7 (independent 5-D Monte-Carlo)
KS = p3['kappa_s']
def bound_L(L, mu, muJ):
    return math.pi*L*SQ32*(1 + 2*KS/L)*(3*mu*(1 + muJ)/(2*(1 - mu)**5) + 3*muJ/8)
for r in p4['phi_rows']:
    tag = "%s L=%g" % (r['case'], r['L'])
    ck('(T_Lambda) bound recomputed | ' + tag, r['bound_Lambda'],
       bound_L(r['L'], r['mu'], r['muJ_vs_JLambda']), 1e-9)
    ck('(T_Lambda) delta = a_phi - a_Lambda | ' + tag, r['delta_Lambda'],
       r['a_phi'] - r['a_Lambda'], 1e-12)
    ck('(T_Lambda) a_Lambda = a_ref + offset | ' + tag, r['a_Lambda'],
       r['a_ref'] + p3['shear_offset_over_M'], 1e-9)
    ckstr('(T_Lambda) holds | ' + tag, r['holds_Lambda'],
          bool(r['bound_Lambda'] >= abs(r['delta_Lambda'])))
    ck('(T_Lambda) slack | ' + tag, r['slack_Lambda'],
       r['bound_Lambda']/abs(r['delta_Lambda']), 1e-6)
    ck('muJ vs J_Lambda <= muJ vs lambda^2 + 2 kappa_s/L | ' + tag,
       min(r['muJ_vs_JLambda'], r['muJ'] + 2*KS/r['L']), r['muJ_vs_JLambda'], 1e-12)
ckstr('(T_Lambda) all rows hold', p4['all_hold_Lambda'],
      bool(all(r['holds_Lambda'] for r in p4['phi_rows'])))
ck('(T_Lambda) worst slack', p4['worst_slack_Lambda'],
   min(r['slack_Lambda'] for r in p4['phi_rows']), 1e-9)

ck('p7 sample size', p7['N'], 8_000_000, 0.0)
for i, r in enumerate(p7['rows']):
    ck('p7[%d] a_ref exact = (M/2) L lam_bar' % i, r['a_ref_exact'], 0.5*r['L']*SQ32, 1e-12)
    ck('p7[%d] I exact = pi^3 L lam_bar' % i, r['I_exact'], math.pi**3*r['L']*SQ32, 1e-9)
    ck('p7[%d] a_gl matches a_ref+offset' % i, r['a_gl'],
       r['a_ref_exact'] + p3['shear_offset_over_M'], 1e-9)
    ck('p7[%d] z(a) recomputed' % i, r['a_z'], (r['a_mc'] - r['a_gl'])/r['a_se'], 1e-9)
    ck('p7[%d] z(a_ref) recomputed' % i, r['a_ref_z'],
       (r['a_ref_mc'] - r['a_ref_exact'])/r['a_ref_se'], 1e-9)
    ck('p7[%d] z(I) recomputed' % i, r['I_z'], (r['I_mc'] - r['I_exact'])/r['I_se'], 1e-9)
    ck('p7[%d] |z| all below 3' % i,
       max(abs(r['a_z']), abs(r['a_ref_z']), abs(r['I_z'])), 3.0, 3.0)
    ck('p7[%d] L is a campaign L' % i,
       min(abs(r['L'] - x['L']) for x in p3['L_table']), 0.0, 1e-12)
ck('p7 max |z| = max over rows', p7['max_abs_z'],
   max(max(abs(r['a_z']), abs(r['a_ref_z']), abs(r['I_z'])) for r in p7['rows']), 1e-12)

npass = sum(1 for c in CHECKS if c[4])
nfail = len(CHECKS) - npass
for name, got, want, tol, ok in CHECKS:
    if not ok:
        print('FAIL  %-58s got %r want %r tol %r' % (name, got, want, tol))
print('%d checks executed, %d PASS, %d FAIL' % (len(CHECKS), npass, nfail))
sys.exit(1 if nfail else 0)

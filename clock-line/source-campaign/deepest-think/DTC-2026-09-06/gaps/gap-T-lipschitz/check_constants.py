"""check_constants.py - re-asserts every number displayed in NOTE.md against the result JSONs."""
import json, numpy as np
s1=json.load(open('s1_results.json')); s2=json.load(open('s2_results.json'))
s3=json.load(open('s3_results.json')); s4=json.load(open('s4_results.json'))
s5=json.load(open('s5_results.json'))
def ok(name, got, want, tol):
    assert abs(got-want) <= tol, f"{name}: {got} vs {want} (tol {tol})"
    print(f"  OK  {name:44s} {got!r}")

print("kernel constants")
ok('3/(8 pi^2)',       s1['K_size_constant_num'], 0.037995443865876666, 1e-15)
ok('3/(2 pi^2)',       s1['gradK_constant_num'],  0.15198177546350666,  1e-15)
ok('sup|e_z-5c what|', float(s1['grad_dir_sup']), 4.0, 0)
ok('pi^3',             s1['pi_cubed'],            31.00627668029982, 1e-13)
assert s1['P_of_lambda'] == 'lambda'
assert s1['antiderivative_residual'] == '0'

print("lemma constants")
ok('C0 = 15 pi/8',     s2['C0_small_mu_value'], 5.890486225480862, 1e-12)
ok('C0 rel = 15 pi/4', s2['C0_relative_form'],  11.780972450961723, 1e-12)
ok('term I  3pi/2',    s2['split']['term_I_small_mu'],  4.71238898038469, 1e-12)
ok('term II 3pi/8',    s2['split']['term_II_small_mu'], 1.1780972450961724, 1e-12)
for k, want in (('0.01',6.1829),('0.05',7.5727),('0.10',9.9566),('0.20',18.4354),
                ('0.25',26.0006),('0.50',227.3728)):
    ok(f'C(mu*={k})', s2['C_table_muJ_eq_mu'][k]['C_mu'], want, 5e-4)

print("numerics (s3)")
ok('L', s3['L'], 8.317766166719343, 1e-12)
ok('max K1 relative error', max(k['rel'] for k in s3['K1_rows']), 0.0, 1e-12)
ok('max K5 residual',       max(k['resid'] for k in s3['K5_rows']), 0.0, 1e-12)
ok('max doubling residual', max(r['conv_rel'] for r in s3['rows']), 0.0, 1e-12)
ok('min slack (bound/|Delta|)', min(r['slack'] for r in s3['rows']), 44.94999405863302, 1e-6)
assert not s3['K1_fired'] and not s3['K3_fired_any'] and not s3['K4_fired_any'] and not s3['K5_fired']
kd = [r['a_T'] for r in s3['rows'] if r['data']=='D2_taper' and r['lam']==1.0][0]/s3['L']
ok('kappa_delta(7.5 deg)', kd, 0.4997212305210886, 1e-12)

print("exact map-A identity  |Delta| = (2 kappa/pi) mu M L  at lambda = 1")
for r in s3['rows']:
    if r['lam']==1.0 and r['map'].startswith('A_'):
        kap = 0.5 if r['data']=='D1_bangbang' else kd
        ok(f"Delta_A {r['data']}", r['absDelta'], 2*kap/np.pi*r['mu_sup']*s3['L'], 3e-6)
ok('C2 ratio == 1/pi (all mu)', max(abs(r['ratio']-1/np.pi) for r in s4['C2']), 0.0, 1e-6)

print("null direction (s5)")
assert s5['cancellation_residual']=='0' and s5['hprime_residual']=='0'
ok('max |Delta|/|a| for full-period ripple, mu<=0.5',
   max(r['rel'] for r in s5['rows'] if r['k']%2==0 and r['mu']<=0.5), 0.0, 1e-13)

print("control (s4)")
assert s4['K6_fired'] is True
ok('C0', s4['C0'], 5.890486225480862, 1e-12)
for mb, want_p in (('0.9',0.9979),('0.5',0.9935),('0.3',0.9855)):
    ok(f'growth exponent p (mu_beta={mb})', s4['growth_fits'][mb]['exponent'], want_p, 5e-4)
for mb, want_c in (('0.9',1.83896),('0.5',0.27202),('0.3',0.04587)):
    ok(f'growth coefficient c (mu_beta={mb})', s4['growth_fits'][mb]['coeff'], want_c, 5e-5)
ok('c ~ mu_beta^p', s4['mu_beta_exponent'], 3.357, 5e-3)
tap = [r['absDelta'] for r in s4['C1_taper']]
ok('taper saturation spread', max(tap)-min(tap), 0.0, 3e-3)
ok('taper plateau value', np.mean(tap), 13.1993, 1e-3)
print("\nALL ASSERTIONS PASSED")

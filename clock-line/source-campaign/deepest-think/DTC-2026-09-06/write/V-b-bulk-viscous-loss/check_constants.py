#!/usr/bin/env python3
"""check_constants.py -- re-assert every number displayed in PROOF.md against the results
files.  Prints FAILURES: none, or lists them.  Nothing is typed from memory: each assertion
compares a literal that appears in PROOF.md with the value produced by b1..b5."""
import json, math, os
H = os.path.dirname(os.path.abspath(__file__))
def J(n): return json.load(open(os.path.join(H, n)))
B1, B2, B3, B4, B5 = J('b1_results.json'), J('b2_results.json'), J('b3_results.json'), J('b4_results.json'), J('b5_results.json')
F = []
def eq(name, got, want, tol=5e-6, rel=True):
    d = abs(got-want)/max(abs(want), 1e-300) if rel else abs(got-want)
    if not (d <= tol): F.append(f"{name}: got {got!r} want {want!r} (dev {d:.3e})")
def ident(name, got, want):
    if str(got).strip() != str(want).strip(): F.append(f"{name}: got {got!r} want {want!r}")

# ---- b1 exact identities
ident('A1', B1['A1_r3_Lap5_inv_r'], '-1')
ident('A2', B1['A2_cart_R3_Lap5_inv_r'], '-1')
ident('A3', B1['A3_bulk_operator_residual'], '0')
ident('B_trace', B1['B_hessian_trace_times_r3'], '-1')
ident('B_op', B1['B_hessian_opnorm_times_r3'], '2')
ident('B_eig', sorted(B1['B_hessian_eigenvalues_at_r'].items()), sorted({'0':1,'-1/r**3':3,'2/r**3':1}.items()))
for k, v in B1['C_legendre_residuals'].items(): ident(f'C_leg_{k}', v, '0')
ident('D_coeffs', B1['D_Delta_y_k_of_inv_r_coefficients'], ['-1','-3','-45','-1575'])
ident('E_ident', B1['E_identification_residual'], '0')
ident('E_sigz', B1['E_sigma_z_absent'], '0')
ident('F_aff2', B1['F_V_affine'], ['0','0'])
ident('F_aff5', B1['F_V_affine_n5'], ['0','0','0','0','0'])
ident('F_nl_lin', B1['F_V_nonlinear_linearised_in_eps'], ['-2*epsilon','0'])
ident('F_n1', B1['F_V_n1_residual_vs_minus_Xpp_over_Xp3'], '0')

# ---- b2 the exact instrument
T0 = B2['T0_covariance']
eq('T0_worst_rel_err', max(T0['rel_err_diag']), 3.634299786245147e-3, 1e-9)
eq('T0_mc_se', math.sqrt(2.0/T0['npath']), 3.1622776601683794e-3, 1e-9)
eq('T0_sigma_ratio', T0['sigma_z']/T0['sigma_y'], 3.5894, 1e-4)
eq('T0_worst_over_se', max(T0['rel_err_diag'])/math.sqrt(2.0/T0['npath']), 1.1493, 1e-3)
eq('T0_maxmean_sd', max(abs(v) for v in T0['mean_over_sd']), 2.0128497619252963, 1e-9)
eq('T0_offdiag', T0['max_offdiag_over_trace'], 7.426100900084578e-4, 1e-9)
want_over_s = [1.015825634, 1.004569342, 1.001507566, 1.000450677, 1.000150075, 1.000045005, 1.000014996]
want_res    = [1.582563, 1.523114, 1.507566, 1.502255, 1.500746, 1.500176, 1.499638]
for i, row in enumerate(B2['T1_pure_plateau']):
    eq(f'T1_over_s_{i}', row['rel_loss_over_s'], want_over_s[i], 1e-8)
    eq(f'T1_res_{i}', row['resid_over_s2'], want_res[i], 1e-5)
z = [r['rel_loss'] for r in B2['T3_sigma_z_independence']]
eq('T3_spread', max(z)/min(z)-1.0, 3.3e-8, 0.2)
for i, v in enumerate(z[:4]): eq(f'T3_val_{i}', v, 1.001507566422e-03, 1e-11)

# ---- b3 the campaign confrontation
want_enh = {1.25:4.3507, 1.5:2.0361, 2.0:1.0411, 2.5:1.0047, 3.0:1.0036, 4.0:1.0036, 6.0:1.0037}
for row in B3['T1_exact_bulk_rate']: eq(f"enh_{row['s0']}", row['enhancement'], want_enh[row['s0']], 1e-4)
def g(run, s0, key): return [x for x in B3[f'T3_exact_gauss_{run}'] if x['s0'] == s0][0][key]
for run, s0, w in [('N4',1.5,3.015), ('N4',2.0,2.106), ('N4',3.0,1.078),
                   ('N4-fine',1.5,1.518), ('N4-fine',2.0,1.072), ('N4-fine',3.0,0.842),
                   ('N4-Re400',3.0,1.149), ('N4-Re25',3.0,0.750)]:
    eq(f'meas_over_exact_{run}_{s0}', g(run, s0, 'meas_over_exact'), w, 6e-4)
eq('bulk_ratio_1p5', g('N4',1.5,'exact_gauss_rel_loss')/g('N4',1.5,'s'), 4.2478, 1e-4)
eq('h_ref_1p5', B3['T4_ratios']['h_refinement_1p5'], 3.8975, 1e-4)
eq('h_ref_2p0', B3['T4_ratios']['h_refinement_2p0'], 15.3862, 1e-4)
ex = B3['T4_ratios']['nu_16x_1p5_excess']; me = B3['T4_ratios']['nu_16x_1p5_meas']
eq('excess_spread', ex[0]/ex[2], 1.8170, 1e-3)
eq('phys_spread', g('N4-Re25',1.5,'exact_gauss_rel_loss')/g('N4-Re400',1.5,'exact_gauss_rel_loss'), 21.6053, 1e-4)

# ---- b4 the necessity of grad^2 u
lk = [r['loss_times_kappa'] for r in B4['S1_kappa_sweep']]
eq('S1_spread', max(lk)/min(lk)-1.0, 1.705442549245e-4, 1e-6)
eq('S4_spread', B4['S4_loss_kappa_product_spread']['max_over_min'], 1.0001705442549245, 1e-12)
for r in B4['S1_kappa_sweep']: eq(f"S1_ratio_{r['kappa']}", r['meas_over_pred'], 1.0, 3e-4)
for r in B4['S2_nu_sweep']:    eq(f"S2_ratio_{r['nu']}", r['meas_over_pred'], 1.0, 1e-3)
res = [r['loss'] for r in B4['S3_resolution']]
eq('S3_spread', (max(res)-min(res))/max(res), 1.1272431201105165e-09, 1e-9)

# ---- b5 sizing
eq('theta_max', B5['theta_max'], 0.7340136762890959, 1e-12)
eq('I2', B5['I2'], 0.5096901045590307, 1e-12); ident('I2_exact', B5['I2_exact'], '4/5 - 16*sqrt(6)/135')
eq('I4', B5['I4'], 1.790579394826636, 1e-12);  ident('I4_exact', B5['I4_exact'], '-4/7 + 27*sqrt(6)/28')
eq('c_window', B5['c_window'], 0.8989794855663562, 1e-12)
eq('C_tau', B5['C_tau_int'], 0.8109302162163288, 1e-12)
eq('Gamma0tau', B5['c_frozen_Gamma0'], 0.7340136762890959, 1e-12)
eq('sigz_over_sigy', B5['sigma_z_over_sigma_y'], 3.513074667941207, 1e-12)
eq('Lmin_1L', B5['L_min_one_dissipation_length_vs_1L'], 24.60320563224208, 1e-9)
eq('Lmin_bulk', B5['L_min_one_dissipation_length_vs_eps_bulk'], 48.33674519170008, 1e-9)
eq('hess_prefactor', B5['hess_prefactor'], 34.420007109996334, 1e-12)
for L, w in [('10',0.7916165916770723), ('40',3.1664663667082893),
             ('160',12.665865466833157), ('640',50.66346186733263)]:
    eq(f'K2hat_{L}', B5['K2hat_threshold'][L], w, 1e-12)
for row in B5['bulk_table']:
    eq(f"epsL_{row['L']}", row['eps_over_target'], 0.03473454, 1e-6)
    eq(f"ratio_{row['L']}", row['ratio'], 1.4401, 1e-4)
want_f = {10:0.25972, 20:0.19228, 40:0.14189, 80:0.10441, 160:0.07663}
for row in B5['inset_table_affine']:
    eq(f"f_{row['L']}", row['f_required'], want_f[row['L']], 1e-4)
eq('c2_L10', [r for r in B5['inset_table_affine'] if r['L']==10][0]['c2'], 1.50272, 1e-5)
want_AP = {5:1.1390, 10:0.48678, 20:0.088903, 40:2.9654e-3, 80:3.2993e-6, 160:4.0840e-12}
for row in B5['brief_one_dissipation_length']:
    eq(f"AP_{row['L']}", row['tail'], want_AP[row['L']], 1e-4)

# ---- remaining displayed rows
for L, w in [(20,1.48105),(40,1.47291),(80,1.46985),(160,1.46871)]:
    eq(f'c2_L{L}', [r for r in B5['inset_table_affine'] if r['L']==L][0]['c2'], w, 1e-5)
want_hess = {1.0:[1.2632,0.3158,0.0790,0.0197], 3.0:[3.7897,0.9474,0.2369,0.0592],
             10.0:[12.6324,3.1581,0.7895,0.1974]}
for row in B5['hess_table']:
    for i, L in enumerate([10,40,160,640]):
        eq(f"hess_{row['K2_hat']}_{L}", row[f'L{L}'], want_hess[row['K2_hat']][i], 2.5e-3)
want_t3 = {('N4',1.5):(2.348e-3,9.975e-3,3.008e-2), ('N4',2.0):(1.407e-3,1.564e-3,3.294e-3),
           ('N4',3.0):(6.805e-4,6.855e-4,7.389e-4),
           ('N4-fine',1.5):(2.348e-3,9.949e-3,1.511e-2), ('N4-fine',2.0):(1.407e-3,1.563e-3,1.675e-3),
           ('N4-fine',3.0):(6.800e-4,6.850e-4,5.765e-4),
           ('N4-Re400',3.0):(1.682e-4,1.693e-4,1.946e-4), ('N4-Re25',3.0):(2.897e-3,2.934e-3,2.201e-3)}
for (run,s0),(a,b,c) in want_t3.items():
    eq(f'sC_{run}_{s0}', g(run,s0,'s'), a, 5e-4)
    eq(f'ex_{run}_{s0}', g(run,s0,'exact_gauss_rel_loss'), b, 5e-4)
    eq(f'me_{run}_{s0}', g(run,s0,'meas_loss'), c, 5e-4)
want_eps = {10:3.4735e-3,20:1.7367e-3,40:8.6836e-4,80:4.3418e-4,160:2.1709e-4,320:1.0855e-4,640:5.4273e-5}
want_lag = {10:5.0022e-3,20:2.5011e-3,40:1.2505e-3,80:6.2527e-4,160:3.1264e-4,320:1.5632e-4,640:7.8159e-5}
for row in B5['bulk_table']:
    eq(f"epsb_{row['L']}", row['eps_bulk_exact'], want_eps[row['L']], 6e-5)
    eq(f"lag_{row['L']}", row['eps_lagrangian'], want_lag[row['L']], 6e-5)

eq('T5_datum', B3['T5_datum_agreement']['max_rel_diff_bulk'], 8.807e-4, 1e-4)
eq('T5_npts', B3['T5_datum_agreement']['npts'], 11376, 0, rel=False)

print("FAILURES:", "none" if not F else "")
for f in F: print("  ", f)
raise SystemExit(1 if F else 0)

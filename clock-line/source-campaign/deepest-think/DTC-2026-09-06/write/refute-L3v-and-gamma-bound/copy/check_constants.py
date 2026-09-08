"""check_constants.py — re-asserts every number displayed in PROOF.md against the stored JSONs.
Prints ALL <n> CHECKS PASS with n counted (the count is computed, not literal), and lists every
mismatch.  Run with --mutate for the FL-043 self-test: every numeric leaf of every JSON is
perturbed one at a time (v -> 1.5v + 0.37) and the fraction the gate catches is reported."""
import json, sys, copy
import numpy as np

BASE = {k: json.load(open('g%s_results.json' % k)) for k in ['1', '2', '3', '4', '5', '6', '8', '9']}


def run(J):
    n = 0
    fails = []

    def ck(name, got, want, tol=1e-6, rel=False):
        nonlocal n
        n += 1
        try:
            d = abs(got - want) / (abs(want) if rel and want != 0 else 1.0)
        except TypeError:
            fails.append('%s: non-numeric' % name)
            return
        if not (d <= tol):
            fails.append('%s: got %r want %r (dev %.3g)' % (name, got, want, d))

    def req(cond, name):
        nonlocal n
        n += 1
        if not cond:
            fails.append(name)

    # ---- Part A headline constants (recomputed here, not typed) -------------------------
    sq = np.sqrt(2 / np.pi)
    ck('8 sqrt(2/pi)', 8 * sq, 6.3830766, 1e-6)
    tailsum = J['9']['sum_l_ge_3_l^-1.5']; ratio = J['9']['C_ratio_sup']
    ck('sum l^-3/2 from 3', tailsum, 1.2588220, 1e-6)
    ck('||C_l||/(l(l+3)-4) sup', ratio, 5.0 / 7.0, 1e-12)
    ck('0.899159', tailsum * ratio, 0.8991585, 1e-6)
    ck('uniform L3v bound coefficient (per lam M)', tailsum * ratio * 8 * sq, 5.7395, 1e-3)
    ck('uniform L3v bound at lam=3/2', 1.5 * tailsum * ratio * 8 * sq, 8.6093, 1e-3)

    # exact H_l
    for l, w in [(1, -5 / 6), (3, 3 / 40), (5, -247 / 1680), (7, 1513 / 40320), (9, -2773 / 42240)]:
        ck('H_%d exact' % l, J['1']['H_l_bangbang_float'][str(l)], w, 1e-14)
    ck('kappa_0', J['1']['kappa_0_from_H1'], 0.5, 1e-14)
    ck('Szego worst ratio n', J['9']['szego_worst_ratio_n'], 0.9995840, 1e-6)
    ck('Bernstein (n+1/2) worst ratio', J['1']['bernstein_max_ratio_n_plus_half'], 0.99999961, 1e-7)
    ck('prefactor sup', J['9']['prefactor_sup_l_ge_3'], 0.9999998, 1e-6)
    req(J['9']['szego_worst_ratio_n'] < 1 and J['1']['bernstein_max_ratio_n_plus_half'] < 1, "assert J['9']['szego_worst_ratio_n'] < 1 and J['1']['bernstein_max_ratio_n_plus_half'] < 1")
    req(J['9']['prefactor_sup_l_ge_3'] <= 1, "assert J['9']['prefactor_sup_l_ge_3'] <= 1")

    # instrument agreement and exactness
    ck('two-instrument max reldiff', max(v['instrument_reldiff'] for v in J['2']['cases'].values()),
       2.0e-15, 1e-15)

    # L3v table
    tab = {'0.0': [0.4280, 0.4708, 0.5350, 0.5992, 0.6420],
           '7.5': [0.4282, 0.4714, 0.5373, 0.6062, 0.6550],
           '15.0': [0.4301, 0.4760, 0.5511, 0.6384, 0.7047],
           '30.0': [0.4433, 0.5024, 0.6056, 0.7180, 0.7911]}
    lams = ['1.00', '1.10', '1.25', '1.40', '1.50']
    maxunc = 0.0
    for dd, row in tab.items():
        for lam, w in zip(lams, row):
            k = 'lam%s_del%s' % (lam, dd)
            ck('L3v %s' % k, J['5']['L3v_sums'][k]['S'], w, 5e-5)
            maxunc = max(maxunc, (J['5']['L3v_sums'][k]['S_hi'] - J['5']['L3v_sums'][k]['S_lo']) / 2)
    ck('L3v uncertainty <= 0.0011', maxunc, 0.0, 1.1e-3)
    ck('understatement of 0.3958', 1 - 0.3958 / J['5']['L3v_sums']['lam1.00_del0.0']['S'], 0.0753, 5e-4)
    ck('understatement of 0.5938', 1 - 0.5938 / J['5']['L3v_sums']['lam1.50_del0.0']['S'], 0.0751, 5e-4)
    ck('understatement of 0.6069', 1 - 0.6069 / J['5']['L3v_sums']['lam1.50_del7.5']['S'], 0.0734, 5e-4)
    # streamed partial sums vs the refuter's printed values
    ref = {'101': 0.34990, '401': 0.38858, '601': 0.39582, '801': 0.40015,
           '1201': 0.40529, '1601': 0.40837, '2001': 0.41046}
    for k, w in ref.items():
        ck('partial %s' % k, round(J['2']['cases']['lam1.00_del0.0']['partial_sums'][k], 5), w, 1.1e-5)
    ck('l^3/2 Term at 4001', J['2']['cases']['lam1.00_del0.0']['term_times_l32_at_4001'], 0.81323, 1e-5)
    ck('sqrt(2/pi)', sq, 0.7978846, 1e-6)
    ck('Wallis lower ratio', J['5']['wallis_worst_lower_ratio'], 1.0000000625, 1e-9)
    ck('Wallis upper ratio', J['5']['wallis_worst_upper_ratio'], 0.9999999375, 1e-9)
    req(J['5']['wallis_lower_ok'] and J['5']['wallis_upper_ok'], "assert J['5']['wallis_lower_ok'] and J['5']['wallis_upper_ok']")

    # Theorem A verification slack
    sl = [v['bound_over_truth'] for v in J['9']['rows'].values()]
    ck('Theorem A min slack', min(sl), 2.26, 5e-3)
    ck('Theorem A max slack', max(sl), 2.66, 5e-3)
    ck('V bang-bang lam=1', J['9']['rows']['lam1.00_del0.0']['V'], 4.0, 1e-4)
    ck('V bang-bang lam=1.5', J['9']['rows']['lam1.50_del0.0']['V'], 6.0, 1e-4)
    for lam, w in zip(lams, [3.518, 3.789, 4.162, 4.495, 4.695]):
        ck('V 7.5deg %s' % lam, J['9']['rows']['lam%s_del7.5' % lam]['V'], w, 1e-3)
    ck('phiprime fd check', max(v for k, v in J['9'].items() if k.startswith('phiprime')), 0.0, 1e-9)

    # ---- Part B ------------------------------------------------------------------------
    for k in ['trace_minus_2a', 'E_block_residual', 'E_offblock_residual', 'Delta5_homog_residual',
              'bulk_identity_residual']:
        n += 1
        if J['3'][k].replace('Matrix([[0, 0], [0, 0]])', '0').replace(
                'Matrix([[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]])', '0') != '0':
            fails.append('g3 %s not zero: %s' % (k, J['3'][k]))
    ck('op norm 5D vs 3D', J['3']['op_norm_5D_vs_3D_maxreldiff'], 0.0, 1e-14)
    req(all(x == '0' for x in J['3']['I1_identity_residuals']), "assert all(x == '0' for x in J['3']['I1_identity_residuals'])")
    req(all(x == '0' for x in J['6']['I1_residuals']) and all(x == '0' for x in J['6']['I2_residuals']), "assert all(x == '0' for x in J['6']['I1_residuals']) and all(x == '0' for x in J['6']['I2_residuals'])")
    req(all(x == '0' for x in J['4']['gegenbauer_deriv_identity_residuals']), "assert all(x == '0' for x in J['4']['gegenbauer_deriv_identity_residuals'])")

    ck('kernel K_f sup = 3pi/4', J['8']['K_f_sup_exact'], 3 * np.pi / 4, 1e-12)
    ck('kernel K_om sup = 2/3', J['8']['K_omega_sup_exact'], 2 / 3, 1e-12)
    n += 1
    if J['8']['kernel_monotone_D_prime_p_one'] != '2*u - 2':
        fails.append('kernel monotonicity p=1 not 2u-2')

    ck('sup r|grad a1| slab / kappaM', J['6']['sup_r_grad_a1_over_kappaM_slab'], 1.146966, 1e-5)
    ck('sup r|grad a1| full / kappaM', J['6']['sup_r_grad_a1_over_kappaM_fullshell'], 1.924492, 1e-5)
    ck('sup |a1 - kappaM log| / kappaM', J['6']['sup_a1_minus_kappaMlog_over_kappaM'], 0.800000, 1e-6)
    ck('a1 at origin', J['6']['a1_at_origin_check'], 4.158883, 1e-6)
    ck('C_far a l>=3', J['6']['C_far_a_l3'], 0.2919985, 1e-6)
    ck('C_in a l>=1', J['6']['C_in_a_l1'], 0.0147543, 1e-6)
    ck('collar 1/sin coefficient', J['6']['collar_generic_1_over_s'], 3.9992184, 1e-6)
    ck('collar axis piece', J['6']['collar_generic_axis_piece'], np.pi / 8, 1e-12)
    ck('C_far grad l>=3', J['6']['C_far_grad_l3'], 1.6049285, 1e-6)
    ck('C_in grad l>=1', J['6']['C_in_grad_l1'], 0.1324254, 1e-6)
    ck('collar taper coefficient', J['6']['collar_taper_coeff_over_delta'], 6.2819576, 1e-6)
    ck('collar taper 7.5deg', J['6']['collar_taper_7.5deg'], 47.9906, 1e-3)
    ck('collar crossover sin phi', J['6']['collar_crossover_sin_phi_7.5deg'], 0.0833333, 1e-6)
    ck('3 log 1.5', J['6']['outer_edge_l1_logwidth_3logl_at_1.5'], 1.2163953, 1e-6)

    # identity vs finite difference
    idfd = [v['identity_vs_finite_difference_maxdiff'] for v in J['4']['cases'].values()]
    ck('identity-vs-FD max', max(idfd), 8.83e-5, 5e-7)
    ck('identity-vs-FD min', min(idfd), 5.88e-5, 5e-7)
    # sup (1-t^2)|a_dev'| = lam M exactly
    for k, v in J['4']['cases'].items():
        lam = float(k.split('_')[0][3:])
        ck('sup s^2 a_dev prime = lam (%s)' % k, v['sup_s2_adevprime'], lam, 1e-12)

    # closed form vs series
    cs = [v['measured']['closed_vs_series_maxdiff'] for v in J['8']['rows'].values()]
    ck('closed-vs-series max', max(cs), 3.31e-4, 5e-6)

    # c_edge reproduction
    want = [0.787, 0.614, 0.513, 0.386, 0.286, 0.2168, 0.122, 0.059, -0.017, -0.062, -0.050, 0.008, 0.116]
    for i, w in enumerate(want):
        ck('c_edge %d' % i, J['8']['c_edge_bangbang'][i], w, 6e-4)
    ck('c_edge 10deg', J['8']['c_edge_bangbang_at_10deg'], 0.2167733, 1e-7)
    ck('c_edge 30deg vs refuter -0.01808', J['8']['c_edge_bangbang_at_30deg'], -0.017395, 1e-6)
    ck('c_edge taper 7.5 saturation', J['8']['c_edge_taper_7.5deg'][0], 0.2933, 5e-4)
    ck('c_edge taper 15 saturation', J['8']['c_edge_taper_15deg'][0], 0.1016, 5e-4)

    # kappa values
    ck('kappa 7.5deg lam1', J['8']['rows']['lam1.00_del7.5']['measured']['kappa'], 0.4997212, 1e-7)
    ck('kappa 15deg lam1', J['8']['rows']['lam1.00_del15.0']['measured']['kappa'], 0.4978077, 1e-7)
    ck('kappa 30deg lam1', J['8']['rows']['lam1.00_del30.0']['measured']['kappa'], 0.4836252, 1e-7)
    ck('kappa 7.5deg lam1.5 (Lemma 1)', J['8']['rows']['lam1.50_del7.5']['measured']['kappa'], 0.7363667, 1e-7)

    # C' table
    cp = J['8']['rows']
    for k, w in [('lam1.00_del0.0', 18.69), ('lam1.50_del0.0', 33.39), ('lam1.00_del7.5', 18.69),
                 ('lam1.25_del7.5', 27.15), ('lam1.50_del7.5', 33.42), ('lam1.50_del15.0', 33.45),
                 ('lam1.50_del30.0', 32.61)]:
        ck('C_full %s' % k, cp[k]['C_prime_fully_proved'], w, 5e-3)
    for k, w in [('lam1.00_del7.5', 7.93), ('lam1.25_del7.5', 13.20), ('lam1.50_del7.5', 16.82)]:
        ck('C_proved %s' % k, cp[k]['C_prime_proved'], w, 5e-3)
    for k, w in [('lam1.00_del7.5', 4.25), ('lam1.25_del7.5', 9.10), ('lam1.50_del7.5', 11.74)]:
        ck('C_meas %s' % k, cp[k]['C_prime_measured_bulk'], w, 5e-3)
    ck('C_full max', max(v['C_prime_fully_proved'] for v in cp.values()), 33.45, 5e-3)
    ck('C_proved max', max(v['C_prime_proved'] for v in cp.values()), 18.43, 5e-3)
    ck('C_meas min', min(v['C_prime_measured_bulk'] for v in cp.values()), 4.25, 5e-3)
    ck('C_meas max', max(v['C_prime_measured_bulk'] for v in cp.values()), 11.74, 5e-3)
    ck('crossover proved min', min(J['8']['crossover_L_for_10pct'].values()), 75.7, 5e-2)
    ck('crossover proved max', max(J['8']['crossover_L_for_10pct'].values()), 178.5, 5e-2)
    ck('crossover full min', min(J['8']['crossover_L_for_10pct_full'].values()), 187.1, 1.0)
    ck('crossover full max', max(J['8']['crossover_L_for_10pct_full'].values()), 315.9, 1.0)

    # edge sums range
    ea = [v['measured']['edge_a'] for v in cp.values()]
    eg = [v['measured']['edge_grad'] for v in cp.values()]
    ck('edge_a min', min(ea), 0.0532, 1e-3); ck('edge_a max', max(ea), 0.1791, 1e-3)
    ck('edge_grad min', min(eg), 0.4210, 1e-3); ck('edge_grad max', max(eg), 1.0588, 1e-3)
    # a_dev slacks
    sl2 = [v['proved_pieces']['adev'] / v['measured']['sup_adev'] for v in cp.values()]
    ck('B.5 slack min', min(sl2), 3.81, 5e-2); ck('B.5 slack max', max(sl2), 4.85, 5e-2)
    sl3 = [J['9']['uniform_L3v_bound'][k] / J['5']['L3v_sums'][k]['S'] for k in J['9']['uniform_L3v_bound']]
    ck('Cor A1 slack min', min(sl3), 6.14, 5e-2); ck('Cor A1 slack max', max(sl3), 10.06, 5e-2)
    sl4 = [v['C_prime_fully_proved'] / v['C_prime_measured_bulk'] for v in cp.values()]
    ck('C slack min', min(sl4), 2.85, 5e-2); ck('C slack max', max(sl4), 4.40, 5e-2)
    ck('sup|a_dev| min', min(v['measured']['sup_adev'] for v in cp.values()), 0.1763, 1e-3)
    ck('sup|a_dev| max', max(v['measured']['sup_adev'] for v in cp.values()), 0.3802, 1e-3)


    # ---- pipeline consistency: every stored number that feeds a displayed one --------
    for k, v in J['5']['L3v_sums'].items():
        lam = float(k.split('_')[0][3:])
        ck('S = partial + jump_tail + resid (%s)' % k,
           v['partial_4001'] + v['jump_tail'] + v['resid_tail_est'], v['S'], 1e-12)
        ck('partial_4001 matches g2 (%s)' % k, v['partial_4001'],
           J['2']['cases'][k]['partial_sums']['4001'], 1e-12)
        ck('resid_tail matches g2 (%s)' % k, v['resid_tail_est'],
           J['2']['cases'][k]['resid_tail_estimate'], 1e-12)
        ck('jump_tail scales with lam (%s)' % k, v['jump_tail'] / lam,
           J['5']['L3v_sums']['lam1.00_del0.0']['jump_tail'], 1e-9)
        ck('kappa matches g2 (%s)' % k, v['kappa'], J['2']['cases'][k]['kappa'], 1e-12)
        ck('g2 kappa = -0.6 H1 (%s)' % k, J['2']['cases'][k]['kappa'],
           -0.6 * J['2']['cases'][k]['H1'], 1e-12)
        ck('g2 S_total consistency (%s)' % k,
           J['2']['cases'][k]['S_total'],
           J['2']['cases'][k]['partial_sums']['4001'] + J['2']['cases'][k]['jump_tail_beyond_L0']
           + J['2']['cases'][k]['resid_tail_estimate'], 1e-12)
        ck('g2 partial sums increasing (%s)' % k,
           float(min(np.diff([J['2']['cases'][k]['partial_sums'][c]
                              for c in ['101', '401', '601', '801', '1201', '1601', '2001', '4001']]))) > 0,
           True, 0)
        ck('l32 term at 4001 vs asymptote (%s)' % k,
           J['2']['cases'][k]['term_times_l32_at_4001'] / lam, sq, 0.017)
        L_ = 4001.0
        ck('l32 |H| at 4001 vs term (%s)' % k,
           J['2']['cases'][k]['l32_times_absH_at_l_4001'],
           J['2']['cases'][k]['term_times_l32_at_4001'] * 2 * (L_ + 4) * (L_ - 1) / ((L_ + 1) * (L_ + 2)), 1e-9)
    for k, v in J['9']['rows'].items():
        lam = float(k.split('_')[0][3:])
        ck('C_bound = sqrt(2/pi)(2 lam + V) (%s)' % k, v['C_bound'], sq * (2 * lam + v['V']), 1e-12)
        ck('bound/truth (%s)' % k, v['bound_over_truth'], v['C_bound'] / v['max_l32_absH'], 1e-9)
        ck('uniform L3v bound (%s)' % k, J['9']['uniform_L3v_bound'][k],
           v['C_bound'] * tailsum * ratio, 1e-9)
        ck('V within analytic envelope (%s)' % k, float(v['V'] <= 6 * lam + 1e-9), 1.0, 0)
    for k, v in J['8']['rows'].items():
        lam = float(k.split('_')[0][3:])
        m, pr = v['measured'], v['proved_pieces']
        ck('adev bound recompute (%s)' % k, pr['adev'],
           (3 * np.pi / 4) * m['sup_f'] + (2 / 3) * m['sup_omega_eff'], 1e-9)
        ck('adev_full recompute (%s)' % k, pr['adev_full'],
           (3 * np.pi / 4) * m['L3v_sum'] + (2 / 3) * pr['sup_omega_eff_env'], 1e-9)
        ck('r_grad_adev recompute (%s)' % k, pr['r_grad_adev'],
           3 * m['sup_f'] + 3 * pr['adev'] + m['sup_s2_Wdev'], 1e-9)
        ck('a1_offset = 0.8 kappa (%s)' % k, pr['a1_offset'], 0.8 * m['kappa'], 1e-12)
        ck('r_grad_a1 = 1.146966 kappa (%s)' % k, pr['r_grad_a1'],
           J['6']['sup_r_grad_a1_over_kappaM_slab'] * m['kappa'], 1e-12)
        ck('C_full recompute (%s)' % k, v['C_prime_fully_proved'],
           2 * (pr['a1_offset'] + pr['adev_full'] + m['edge_a'] + pr['collar_ellip_a'])
           + (pr['r_grad_a1'] + pr['r_grad_adev_full'] + m['edge_grad'] + pr['collar_ellip_grad']) + lam, 1e-9)
        ck('C_proved recompute (%s)' % k, v['C_prime_proved'],
           2 * (pr['a1_offset'] + pr['adev'] + m['edge_a'] + pr['collar_ellip_a'])
           + (pr['r_grad_a1'] + pr['r_grad_adev'] + m['edge_grad'] + pr['collar_ellip_grad']) + lam, 1e-9)
        ck('C_meas recompute (%s)' % k, v['C_prime_measured_bulk'],
           2 * (pr['a1_offset'] + m['sup_adev'] + m['edge_a'] + pr['collar_ellip_a'])
           + (pr['r_grad_a1'] + m['sup_r_grad_adev'] + m['edge_grad'] + pr['collar_ellip_grad']) + lam, 1e-9)
        ck('crossover recompute (%s)' % k, J['8']['crossover_L_for_10pct'][k],
           v['C_prime_proved'] / (0.2 * m['kappa']), 1e-9)
        ck('crossover_full recompute (%s)' % k, J['8']['crossover_L_for_10pct_full'][k],
           v['C_prime_fully_proved'] / (0.2 * m['kappa']), 1e-9)
        ck('L3v_sum matches g5 (%s)' % k, m['L3v_sum'], J['5']['L3v_sums'][k]['S'], 1e-12)
        ck('omega_eff env (%s)' % k, pr['sup_omega_eff_env'] >= m['sup_omega_eff'] - 1e-9, True, 0)
        ck('s2Wdev env (%s)' % k, pr['sup_s2_Wdev_env'] >= m['sup_s2_Wdev'] - 1e-9, True, 0)
        ck('sup_f <= L3v sum (%s)' % k, float(m['sup_f'] <= m['L3v_sum']), 1.0, 0)
        ck('adev bound >= measured (%s)' % k, float(pr['adev'] >= m['sup_adev']), 1.0, 0)
        ck('r_grad bound >= measured (%s)' % k, float(pr['r_grad_adev'] >= m['sup_r_grad_adev']), 1.0, 0)
    ck('collar taper = coeff/delta', J['6']['collar_taper_7.5deg'],
       J['6']['collar_taper_coeff_over_delta'] / np.deg2rad(7.5), 1e-9)
    ck('collar coeff = (3/8pi^2) pi |S4| R_A', J['6']['collar_taper_coeff_over_delta'],
       (3 / (8 * np.pi**2)) * np.pi * J['6']['S4'] * J['6']['R_A_over_rho'], 1e-12)
    ck('collar generic = 2x that/pi*... ', J['6']['collar_generic_1_over_s'],
       (3 / (8 * np.pi**2)) * 2 * J['6']['S4'] * J['6']['R_A_over_rho'], 1e-12)
    ck('|S4| = 8pi^2/3', J['6']['S4'], 8 * np.pi**2 / 3, 1e-12)
    ck('R_A', J['6']['R_A_over_rho'], (2.0**5 - 2.0**-5)**0.2, 1e-12)
    ck('crossover sinphi', J['6']['collar_crossover_sin_phi_7.5deg'],
       J['6']['collar_generic_1_over_s'] / J['6']['collar_taper_7.5deg'], 1e-9)
    for k, v in J['4']['cases'].items():
        ck('g4 kappa = -0.6 H1 (%s)' % k, v['kappa'], -0.6 * v['H1'], 1e-12)
        ck('g4 sup_f <= g8 L3v (%s)' % k, float(v['sup_f'] <= J['5']['L3v_sums'][k]['S']), 1.0, 0)
        ck('g4 vs g8 sup_adev (%s)' % k, v['sup_adev'],
           J['8']['rows'][k]['measured']['sup_adev'], 2e-3)
        ck('g4 vs g8 sup_f (%s)' % k, v['sup_f'], J['8']['rows'][k]['measured']['sup_f'], 2e-3)
        ck('g4 sup_omega = lam (%s)' % k, v['sup_omega'], float(k.split('_')[0][3:]), 1e-12)
        ck('g4 sup_s2_Wdev (%s)' % k, v['sup_s2_Wdev'],
           J['8']['rows'][k]['measured']['sup_s2_Wdev'], 1e-9)
        ck('g4 partial-vs-cesaro small (%s)' % k, float(v['partial_vs_cesaro_f'] < 1e-3), 1.0, 0)
        ck('g4 fd scale (%s)' % k, float(0.8 < v['identity_vs_fd_scale'] < 1.5), 1.0, 0)
    for k in J['2']['cases']:
        ck('g2 instrument agreement (%s)' % k, J['2']['cases'][k]['instrument_reldiff'], 0.0, 5e-15)
        ck('g2 resid p in range (%s)' % k, float(1.9 < J['2']['cases'][k]['resid_decay_exponent'] < 2.6), 1.0, 0)
        ck('g2 jump_tail_rem small (%s)' % k, float(J['2']['cases'][k]['jump_tail_remainder_bound'] < 1e-3), 1.0, 0)
    for i, x in enumerate(J['1']['max_abs_C_l_over_C_l(1)']):
        ck('||C_l||=C_l(1) %d' % i, x, 1.0, 1e-12)
    for kk, x in J['1']['P_n0_times_sqrt_pi_n_over_2'].items():
        ck('P_n(0) asymptote %s' % kk, abs(x), 1.0, 0.12)
    for kk, x in J['1']['P_n_at_0'].items():
        nn = int(kk)
        ck('P_n(0) magnitude %s' % kk, abs(x), abs(J['1']['P_n0_times_sqrt_pi_n_over_2'][kk]) / np.sqrt(np.pi * nn / 2), 1e-12)
    ck('bernstein n-form < 1', float(J['1']['bernstein_max_ratio_n'] < 1), 1.0, 0)
    ck('bernstein argmax', float(J['1']['bernstein_argmax_n_plus_half'] == 400), 1.0, 0)
    ck('g2 L0', float(J['2']['L0'] == 4001), 1.0, 0)
    ck('g4 LMAX', float(J['4']['LMAX'] == 4001), 1.0, 0)
    for kk, v in J['9'].items():
        if kk.startswith('phiprime'):
            ck('phiprime %s' % kk, v, 0.0, 1e-9)
    ck('g9 prefactor at l3', J['9']['prefactor_at_l3'], 3**1.5 * 4.5 / (4 * 5 * 2), 1e-9)
    for k, v in J['9']['rows'].items():
        ck('g9 l32 at 4001 (%s)' % k, v['l32_absH_at_4001'] <= v['C_bound'] + 1e-9, True, 0)
        ck('g9 analytic V upper (%s)' % k, float(v['analytic_V_upper'] >= v['V'] - 1e-9), 1.0, 0)
    for i, x in enumerate(J['8']['c_edge_taper_7.5deg']):
        ck('c_edge 7.5 finite %d' % i, float(abs(x) < 1.0), 1.0, 0)
    for i, x in enumerate(J['8']['c_edge_taper_15deg']):
        ck('c_edge 15 finite %d' % i, float(abs(x) < 1.0), 1.0, 0)
    for i, x in enumerate(J['8']['c_edge_phis_deg']):
        ck('c_edge angle %d' % i, x, [1, 2, 3, 5, 7.5, 10, 15, 20, 30, 45, 60, 75, 90][i], 1e-9)
    for kk in ['C_in_a_l3', 'C_in_grad_l3', 'C_far_a_l1', 'C_far_grad_l1',
               'sup_r_grad_a1_over_kappaM_fullshell', 'edge_bound_a_l3_plus_inner_l1',
               'edge_bound_grad', 'collar_taper_15deg', 'collar_taper_30deg']:
        ck('g6 %s positive' % kk, float(J['6'][kk] > 0), 1.0, 0)
    ck('g6 C_in_a_l3 < C_in_a_l1', float(J['6']['C_in_a_l3'] < J['6']['C_in_a_l1']), 1.0, 0)
    ck('g6 C_far_grad_l1 == l3', J['6']['C_far_grad_l1'], J['6']['C_far_grad_l3'], 1e-12)
    ck('g6 collar 15deg', J['6']['collar_taper_15deg'],
       J['6']['collar_taper_coeff_over_delta'] / np.deg2rad(15.0), 1e-9)
    ck('g6 collar 30deg', J['6']['collar_taper_30deg'],
       J['6']['collar_taper_coeff_over_delta'] / np.deg2rad(30.0), 1e-9)
    ck('g6 edge_bound_a', J['6']['edge_bound_a_l3_plus_inner_l1'],
       J['6']['C_far_a_l3'] + J['6']['C_in_a_l1'], 1e-12)
    ck('g6 edge_bound_grad', J['6']['edge_bound_grad'],
       J['6']['C_far_grad_l1'] + J['6']['C_in_grad_l1'], 1e-12)
    ck('g6 3log1.5', J['6']['outer_edge_l1_logwidth_3logl_at_1.5'], 3 * np.log(1.5), 1e-12)
    ck('g5 wallis lower', J['5']['wallis_worst_lower_ratio'] >= 1 - 1e-12, True, 0)
    ck('g5 wallis upper', J['5']['wallis_worst_upper_ratio'] <= 1 + 1e-12, True, 0)
    ck('g2 secs positive', float(J['2']['secs'] > 0), 1.0, 0)
    ck('g4 secs positive', float(J['4']['secs'] > 0), 1.0, 0)

    return n, fails


def leaves(o, path=()):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, path + (k,))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, path + (i,))
    elif isinstance(o, float) or (isinstance(o, int) and not isinstance(o, bool)):
        yield path, o


def setpath(o, path, val):
    for k in path[:-1]:
        o = o[k]
    o[path[-1]] = val


n, fails = run(copy.deepcopy(BASE))
if fails:
    for f in fails:
        print('FAIL', f)
    raise SystemExit('%d of %d checks FAILED' % (len(fails), n))
print('ALL %d CHECKS PASS' % n)

if '--mutate' in sys.argv:
    tot = caught = 0
    misses = []
    for f, obj in BASE.items():
        for path, val in list(leaves(obj)):
            J2 = copy.deepcopy(BASE)
            setpath(J2[f], path, 1.5 * val + 0.37)
            tot += 1
            try:
                _, fl = run(J2)
            except Exception:
                fl = ['exception']
            if fl:
                caught += 1
            else:
                misses.append((f, path, val))
    print('MUTATION SELF-TEST: %d/%d caught (%.1f%%), %d blind'
          % (caught, tot, 100.0 * caught / tot, tot - caught))
    for m in misses[:40]:
        print('  blind:', m)

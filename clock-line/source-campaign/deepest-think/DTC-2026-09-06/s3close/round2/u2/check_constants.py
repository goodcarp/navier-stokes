#!/usr/bin/env python3
"""
check_constants.py -- re-assert every number displayed in PROOF.md against the
JSON results produced by u1..u5.  FL-043 discipline: a displayed constant that is
not checked here is a defect.  Exit code 1 on any failure.

Usage:  python3 check_constants.py
"""
import json, math, sys, os

FAIL = []
N = [0]

def load(p):
    with open(p) as f:
        return json.load(f)

def chk(name, got, want, rtol=1e-9, atol=0.0):
    N[0] += 1
    try:
        g = float(got); w = float(want)
    except (TypeError, ValueError):
        if got != want:
            FAIL.append(f'{name}: {got!r} != {want!r}')
        return
    if not (abs(g - w) <= max(atol, rtol*abs(w))):
        FAIL.append(f'{name}: {g!r} != {w!r} (rel {abs(g-w)/max(abs(w),1e-300):.3e})')

def zero(name, s):
    N[0] += 1
    t = str(s).replace(' ', '')
    if t not in ('0', 'Matrix([[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]])'):
        FAIL.append(f'{name}: residual not zero: {s!r}')

u1 = load('u1_results.json'); u2 = load('u2_results.json')
u3 = load('u3_results.json'); u5 = load('u5_results.json')
u4 = load('u4_results.json') if os.path.exists('u4_results.json') else None

# ---- sec 2 and sec 5: the exact algebra -------------------------------------
zero('A1 div5b-2a', u1['A1_div5b_minus_2a_residual'])
zero('A2 grad5b structure', u1['A2_grad5b_structure_residual'])
zero('C advective split', u1['C_advective_split_residual'])
zero('C grad commutes Lap5', u1['C_grad_commutes_with_Lap5_residual'])
for k in ('D_lap5_rho1_g_residual', 'D_lap5_rho2_g_residual',
          'D_k1_rearranged_residual', 'D_k2_rearranged_residual'):
    zero(k, u1[k])
chk('A3 worst ratio', u1['A3_worst_ratio_op_over_bound'], 0.9989318482494669, 1e-12)
N[0] += 1
if u1['A3_violated']:
    FAIL.append('A3: envelope 2|a|+r|grad a|+|omega| violated')
chk('b(0) r', u1['B_b_at_origin'][0], 0.0, atol=0)
chk('b(0) z', u1['B_b_at_origin'][1], 0.0, atol=0)

# ---- sec 6: kernel constants -------------------------------------------------
chk('C_K', u1['F_C_K'], 3.0/(8*math.pi**2), 1e-15)
chk('|S^4|', u1['F_S4'], 8*math.pi**2/3, 1e-15)
chk('C_K |S^4| = 1', u1['F_C_K_times_S4'], 1.0, 1e-14)
chk('3 pi^2/16', u1['F_3pi2_over_16'], 1.8505508252042546, 1e-14)
chk('C_K pi^4/2 = 3pi^2/16', u1['F_C_K_times_pi4_over_2'], u1['F_3pi2_over_16'], 1e-13)
chk('Riesz closed', u1['E_riesz_closed_pi4_over_2'], 48.70454551700121, 1e-13)
chk('Riesz classical formula', u1['E_riesz_classical_formula'], 48.70454551700121, 1e-13)
chk('Riesz numeric', u1['E_riesz_numeric'], 48.70454551715526, 1e-12)
chk('Riesz rel err', u1['E_riesz_rel_err'], 3.1630060173002208e-12, 1e-6)
chk('R_A/rho', u1['F_R_A_over_rho'], 1.9996092223226414, 1e-14)
chk('collar coeff', u1['F_collar_coeff_4RA_over_2'], 3.9992184446452828, 1e-14)

# ---- sec 1/3/5: datum constants ---------------------------------------------
chk('E0', u2['E0_axis_closed'], 7.66060196237754, 1e-12)
chk('E0 (scan)', u2['datum_L10']['E0'], 7.66060196237754, 1e-9)
chk('1/sin delta', u2['one_over_sin_delta'], 7.66129757554039, 1e-12)
chk('G0', u2['datum_L10']['G0'], 20.919707021704003, 1e-9)
chk('G0 at L=40', u2['datum_L40']['G0'], 20.919707021704003, 1e-9)
chk('F0', u2['datum_L10']['F0'], 8.494612736878054, 1e-8)
chk('J_ac reproduces hk2 65.625910', u2['J_shell_ac_campaign'], 65.625910, 1e-6)
chk('datum delta_deg', u2['delta_deg'], 7.5, 1e-12)
chk('datum w_eq', u2['w_eq'], 0.20, 1e-12)
chk('datum eps_r', u2['eps_r'], 0.25, 1e-12)
chk('G0 attained at phi (deg)', u2['datum_L10']['G0_at'][1], 6.660, 1e-4)

# ---- sec 7: far/near constants ----------------------------------------------
chk('C1 far z-odd', u3['C1_far_zodd'], 0.29199853254882013, 1e-10)
chk('C1 vs record 0.291999', u3['C1_far_zodd'], 0.291999, 5e-6)
chk('C2 inner z-odd', u3['C2_inner_zodd'], 0.014754271582665845, 1e-10)
chk('C2 vs record 0.014754', u3['C2_inner_zodd'], 0.014754, 5e-5)
chk('collar coeff (u3)', u3['collar_coef'], 3.9992184446452828, 1e-13)
chk('collar vs record 3.999218', u3['collar_coef'], 3.999218, 5e-7)
chk('pi/8', u3['collar_axis_piece_pi_over_8'], math.pi/8, 1e-15)
chk('sigma crossover', u3['sigma_crossover_lam1.5_cG1.2164'], 0.05834931948483575, 1e-10)
chk('phi crossover deg', math.degrees(math.asin(u3['sigma_crossover_lam1.5_cG1.2164'])),
    3.345, 1e-3)
chk('c window 2log(3/2)', u3['c_window_frozen'], 0.8109302162163288, 1e-14)
chk('u3 copy of E0', u3['E0'], 7.66060196237754, 1e-12)
chk('u3 copy of G0', u3['G0'], 20.919707021704003, 1e-9)
chk('u3 copy of 3pi^2/16', u3['riesz_C_K_pi4_2'], 1.8505508252042546, 1e-14)
chk('record C1 quoted', u3['record_far_near']['C1'], 0.291999, 1e-12)
chk('record C2in quoted', u3['record_far_near']['C2in'], 0.014754, 1e-12)
chk('record collar quoted', u3['record_far_near']['collar'], 3.999218, 1e-12)
chk('record axis quoted', u3['record_far_near']['axis'], math.pi/8, 1e-14)
chk('sigma crossover lam1', u3['sigma_crossover_lam1.0_cG0.8109'],
    0.05834931948483575, 1e-10)

# Chat_a table
chk('Chat_a lam1 sig1/2 = ASSEMBLY 8.697888',
    u3['Linf_lam1.0_p2']['Chat_a'], 8.697888775120775, 1e-11)
chk('Chat_a lam1 vs ASSEMBLY 8.697888', u3['Linf_lam1.0_p2']['Chat_a'], 8.697888, 1e-6)
chk('Chat_a lam1.25 sig1/2', u3['Linf_lam1.25_p2']['Chat_a'], 10.8724, 1e-4)
chk('Chat_a lam1.5 sig1/2', u3['Linf_lam1.5_p2']['Chat_a'], 13.046833162681162, 1e-11)
chk('Chat_a lam1 axis', u3['Linf_lam1.0_axis_p2']['Chat_a'], 69.238699, 1e-8)
chk('Chat_a lam1 sin7', u3['Linf_lam1.0_sin7_p2']['Chat_a'], 31.338654, 1e-7)
chk('Chat_a lam1.25 sin7', u3['Linf_lam1.25_sin7_p2']['Chat_a'], 39.173318, 1e-7)
chk('Chat_a lam1.5 sin7', u3['Linf_lam1.5_sin7_p2']['Chat_a'], 47.007982, 1e-7)
chk('Chat_a lam1.25 axis', u3['Linf_lam1.25_axis_p2']['Chat_a'], 84.8075, 1e-4)
chk('Chat_a lam1.5 axis', u3['Linf_lam1.5_axis_p2']['Chat_a'], 103.858, 1e-4)

# sec 8.1 table
for lam, p3, p2, mp, cpp in [(1.0, 440.9651, 195.9845, 96.0324, 214.3802),
                             (1.25, 810.1045, 293.9767, 183.7737, 316.9714),
                             (1.5, 1488.257, 440.9651, 324.1093, 468.5587)]:
    chk(f'Ghat p3 lam{lam}', u3[f'Linf_lam{lam}_p3']['Ghat'], p3, 1e-5)
    chk(f'Ghat p2 lam{lam}', u3[f'Linf_lam{lam}_p2']['Ghat'], p2, 1e-5)
    chk(f'Ghat map lam{lam}', u3[f'Linf_lam{lam}_map_mu0.05']['Ghat'], mp, 1e-5)
    chk(f'C2prime p2 lam{lam}', u3[f'Linf_lam{lam}_p2']['C2prime'], cpp, 1e-5)
    for tg in ('p3', 'p2', 'map_mu0.05'):
        chk(f'c_G Linf lam{lam} {tg}', u3[f'Linf_lam{lam}_{tg}']['c_G'],
            lam*0.8109302162163288, 1e-12)
for lam, cpp in [(1.0, 335.4619), (1.25, 464.8417), (1.5, 650.1811)]:
    chk(f'C2prime p2 axis lam{lam}', u3[f'Linf_lam{lam}_axis_p2']['C2prime'], cpp, 1e-4)

# sec 8.2 fixed-point table
def row(mode, mu, lam, sig):
    for r in u3['table']:
        if (r['mode'] == mode and abs(r['mu']-mu) < 1e-12 and abs(r['lam']-lam) < 1e-12
                and abs(r['sigma_star']-sig) < 1e-12):
            return r
    raise KeyError((mode, mu, lam, sig))

H = 0.5
for mode, mu, lam, Ls, C4, cg, pp in [
        ('proved', 0.0, 1.0, 2189.7, 232.821, 0.8298, 2.0629),
        ('proved', 0.0, 1.25, 3277.2, 361.4774, 1.0430, 2.0790),
        ('proved', 0.0, 1.5, 4997.6, 582.2400, 1.2636, 2.1068),
        ('p3', 0.0, 1.5, 9908.9, 3604.6869, 1.5087, None),
        ('map', 0.0, 1.5, 670.0, 329.5321, 1.2431, 2.0591),
        ('map', 0.02, 1.5, 696.2, 342.0495, 1.2441, 2.0615),
        ('map', 0.05, 1.5, 736.5, 361.3406, 1.2457, 2.0652)]:
    r = row(mode, mu, lam, H)
    chk(f'Lstar {mode} mu{mu} lam{lam} sig1/2', r['Lstar'], Ls, 1e-4)
    chk(f'C2prime(1e4) {mode} mu{mu} lam{lam}', r['L10000']['C2prime'], C4, 1e-6)
    chk(f'c_G(1e4) {mode} mu{mu} lam{lam}', r['L10000']['c_G'], cg, 1e-4)
    if pp is not None:
        chk(f'p(1e4) {mode} mu{mu} lam{lam}', r['L10000']['p'], pp, 1e-4)

for mode, mu, lam, Ls, C4 in [('proved', 0.0, 1.0, 2435.2, 362.8466),
                              ('proved', 0.0, 1.25, 3577.2, 526.7639),
                              ('proved', 0.0, 1.5, 5313.8, 801.0451),
                              ('map', 0.05, 1.5, 1172.3, 557.8960)]:
    r = row(mode, mu, lam, 0.0)
    chk(f'Lstar axis {mode} mu{mu} lam{lam}', r['Lstar'], Ls, 1e-4)
    chk(f'C2prime(1e4) axis {mode} mu{mu} lam{lam}', r['L10000']['C2prime'], C4, 1e-6)

sd = math.sin(math.radians(7.5))
r = row('proved', 0.0, 1.5, sd)
chk('Lstar sin7.5 proved lam1.5', r['Lstar'], 5088.2, 1e-4)
chk('C2prime(1e4) sin7.5 proved', r['L10000']['C2prime'], 658.118, 1e-5)
r = row('map', 0.05, 1.5, sd)
chk('Lstar sin7.5 map0.05', r['Lstar'], 788.2, 1e-4)
chk('C2prime(1e4) sin7.5 map0.05', r['L10000']['C2prime'], 431.1578, 1e-6)

# 8.3 mu pricing at L -> infinity
for m, w in (('0.00', 321.5704), ('0.02', 333.4470), ('0.05', 351.7030)):
    chk(f'C2prime Linf map mu{m} lam1.5', u3[f'Linf_lam1.5_map_mu{m}']['C2prime'], w, 1e-6)
chk('mu pricing pct on C2prime',
    100*(u3['Linf_lam1.5_map_mu0.05']['C2prime']/u3['Linf_lam1.5_map_mu0.00']['C2prime']-1),
    9.37, 1e-3)
chk('mu pricing pct on Lstar',
    100*(row('map', 0.05, 1.5, 0.5)['Lstar']/row('map', 0.0, 1.5, 0.5)['Lstar'] - 1),
    9.9183, 1e-4)

# ---- sec 4/5: u5 ------------------------------------------------------------
chk('Lambda_max closed form err', u5['a_closed_form_max_abs_rel_err'],
    1.023743022700939e-15, 1e-6)
chk('envelope violations', u5['b_envelope_violations'], 0, atol=0)
chk('envelope worst ratio', u5['b_envelope_worst_ratio'], 0.9984150986286499, 1e-10)
chk('Gamma_rad/Gamma strain', u5['c_strain_Gamma_rad_over_Gamma'], 0.5, 1e-14)
for lam in (1.25, 1.5):
    d = u5[f'c_strain_lam{lam}']
    chk(f'P1 sharp lam{lam}', d['P1_rel'], 0.0, atol=1e-15)
    chk(f'P2 sharp lam{lam}', d['P2_rel'], 0.0, atol=1e-15)
    chk(f'P1 exact lam{lam}', d['P1_exact'], lam, 1e-14)
    chk(f'P2 exact lam{lam}', d['P2_exact'], lam**4, 1e-14)
    chk(f'P1 maxprinciple lam{lam}', d['P1_maxprinciple'], lam, 1e-14)
    chk(f'P2 maxprinciple lam{lam}', d['P2_maxprinciple'], lam**4, 1e-14)
    chk(f'c_G strain lam{lam}', d['c_G'], 2*math.log(lam), 1e-14)
    chk(f'crude P1 lam{lam}', d['crude_P1_eGamma'], lam**2, 1e-13)
    chk(f'||T_lam|| lam{lam}', u5[f'd_Tlam_op_lam{lam}'], lam, 1e-13)
    chk(f'||T_lam^-1|| lam{lam}', u5[f'd_Tlam_inv_op_lam{lam}'], lam**2, 1e-13)
chk('crude P2 at lam1.5', u5['c_strain_lam1.5']['crude_P2_e3Gamma'], 11.390625, 1e-12)

# ---- derived numbers quoted in the prose ------------------------------------
chk('ratio 16.55 (bound vs hk2 TV bound)', 1.8505508252042546*20.919707021704003
    / (4.6795*0.5), 16.55, 1e-3)
chk('predicted 16.78', 2*26.318945069571622*20.919707021704003/65.62590995517674,
    16.78, 1e-3)
chk('|S4| G0 / J  = 8.39', 26.318945069571622*20.919707021704003/65.62590995517674,
    8.39, 1e-3)
_cg = 1.5*2*math.log(1.5)
_pen = 26.318945069571622*u2['datum_L10']['G0']/u2['J_shell_ac_campaign']
chk('e^{4c_G} = 129.8', math.exp(4*_cg), 129.8, 1e-3)
chk('e^{5c_G} = 437.89', math.exp(5*_cg), 437.89, 1e-4)
chk('TV-route penalty p=3', math.exp(4*_cg)/_pen, 15.47, 1e-3)
chk('TV-route penalty p=2', math.exp(5*_cg)/_pen, 52.19, 1e-3)

# ---- sec 9: the confrontation (only if u4 has finished) ----------------------
if u4 is not None:
    chk('a0 control rel (L=10)', u4['a0_control_L10']['rel'], 0.0, atol=2e-3)
    chk('a0 control rel (L=40)', u4['a0_control_L40']['rel'], 0.0, atol=2e-3)
    for k, v in u4.items():
        if k.startswith('conv_'):
            N[0] += 1
            if max(v['rel'][:1]) > 1e-4 or max(v['rel']) > 5e-2:
                FAIL.append(f'{k}: convergence {v["rel"]}')
    # R1 / R3 / R4 controls
    chk('R1 a0 control L10 rel', u4['a0_control_L10']['rel'], 1.0936332280292203e-08, 1e-6)
    chk('R1 a0 control L40 rel', u4['a0_control_L40']['rel'], 2.6302605883077273e-09, 1e-6)
    chk('R3 kappa_profile (D-C)', u4['kappa_profile'], 0.4747343553883626, 1e-9)
    # the measurement table
    tabl = [(10, 1.0, 'exact', 4.50998, 8.49649, -0.52347, 0.56703),
            (10, 1.0, 'pert0.02', 4.50886, 8.50075, -0.51697, 0.57478),
            (10, 1.0, 'pert0.05', 4.50634, 8.50532, -0.50736, 0.58633),
            (10, 1.5, 'exact', 6.69585, 12.57542, -0.81629, 0.96387),
            (10, 1.5, 'pert0.02', 6.69520, 12.57948, -0.81092, 0.97093),
            (10, 1.5, 'pert0.05', 6.69356, 12.58392, -0.80321, 0.98142),
            (40, 1.0, 'exact', 18.75201, 36.98055, -0.52347, 0.56705),
            (40, 1.0, 'pert0.05', 18.74756, 36.97815, -0.51696, 0.59040),
            (40, 1.5, 'exact', 27.84065, 54.86501, -0.81629, 0.96387),
            (40, 1.5, 'pert0.05', 27.83646, 54.86332, -0.80961, 0.96977)]
    for L, lam, nm, a0, gm, gm2, rga in tabl:
        d = u4[f'L{L}_lam{lam}_{nm}']
        chk(f'a0 L{L} lam{lam} {nm}', d['a0'], a0, 1e-5)
        chk(f'Gmax L{L} lam{lam} {nm}', d['Gamma_max_measured'], gm, 1e-5)
        chk(f'Gmax-2a0 L{L} lam{lam} {nm}', d['Gamma_max_measured']-2*d['a0'], gm2, 1e-4)
        chk(f'r|grad a| L{L} lam{lam} {nm}', d['max_r_grad_a'], rga, 1e-4)
        N[0] += 1
        if d['max_r_grad_a'] >= lam:
            FAIL.append(f'r|grad a| >= lam M at L{L} lam{lam} {nm}')
    # R4: L-independence of r|grad a|
    for lam, want in ((1.0, 3.6e-5), (1.5, 1.3e-8)):
        r10 = u4[f'L10_lam{lam}_exact']['max_r_grad_a']
        r40 = u4[f'L40_lam{lam}_exact']['max_r_grad_a']
        chk(f'R4 L-independence lam{lam}', abs(r10-r40)/r10, 0.0, atol=want)
    # perturbation percentages
    for L, lam, e2, e5 in ((10, 1.0, 1.37, 3.40), (10, 1.5, 0.73, 1.82),
                           (40, 1.0, None, 4.12), (40, 1.5, None, 0.61)):
        b = u4[f'L{L}_lam{lam}_exact']['max_r_grad_a']
        if e2 is not None:
            chk(f'pct 0.02 L{L} lam{lam}',
                100*(u4[f'L{L}_lam{lam}_pert0.02']['max_r_grad_a']/b - 1), e2, 5e-3)
        chk(f'pct 0.05 L{L} lam{lam}',
            100*(u4[f'L{L}_lam{lam}_pert0.05']['max_r_grad_a']/b - 1), e5, 5e-3)
    # perturbation norms and mu
    for L, eps, rel, mu1, mu15 in ((10, 0.02, 0.0166542125868058, 0.01665, 0.05621),
                                   (10, 0.05, 0.04163553146701449, 0.04164, 0.14052),
                                   (40, 0.02, 0.017292463430167545, None, None),
                                   (40, 0.05, 0.043231158575418864, None, None)):
        d = u4[f'pertnorm_L{L}_eps{eps}']
        chk(f'sup|psi|/|a| L{L} eps{eps}', d['sup_rel'], rel, 1e-9)
        chk(f'||Dpsi|| L{L} eps{eps}', d['sup_Dpsi'], eps, 1e-9)
        if mu1 is not None:
            chk(f'mu lam1 L{L} eps{eps}', d['sup_rel'], mu1, 1e-3)
            chk(f'mu lam1.5 L{L} eps{eps}', d['sup_rel']*1.5**3, mu15, 1e-3)
    # slack table
    G1 = u3['Linf_lam1.0_p2']['Ghat']; G15 = u3['Linf_lam1.5_p2']['Ghat']
    chk('slack r|grad a| lam1', G1/u4['L10_lam1.0_exact']['max_r_grad_a'], 345.6, 1e-3)
    chk('slack r|grad a| lam1.5', G15/u4['L10_lam1.5_exact']['max_r_grad_a'], 457.5, 1e-3)
    C1p = u3['Linf_lam1.0_p2']['C2prime']; C15p = u3['Linf_lam1.5_p2']['C2prime']
    b1 = 2*u4['L40_lam1.0_exact']['a0'] + C1p
    b15 = 2*u4['L40_lam1.5_exact']['a0'] + C15p
    chk('bound Gamma L40 lam1', b1, 251.88, 1e-4)
    chk('bound Gamma L40 lam1.5', b15, 524.24, 1e-4)
    chk('slack Gamma L40 lam1', b1/u4['L40_lam1.0_exact']['Gamma_max_measured'], 6.81, 1e-3)
    chk('slack Gamma L40 lam1.5', b15/u4['L40_lam1.5_exact']['Gamma_max_measured'], 9.55, 1e-3)
    chk('decomposition 11.39 x 8.39 x 4.79', 11.390625*8.39*4.79, 457.5, 5e-3)
    # convergence (R2): the displayed envelope, exactly
    convs = {k: v for k, v in u4.items() if k.startswith('conv_')}
    N[0] += 1
    if len(convs) != 4:
        FAIL.append(f'expected 4 convergence rows, got {len(convs)}')
    ra = [v['rel'][0] for v in convs.values()]
    rr = [v['rel'][1] for v in convs.values()]
    rz = [v['rel'][2] for v in convs.values()]
    chk('R2 max rel on a', max(ra), 1.9200822805008024e-06, 1e-9)
    chk('R2 min rel on a', min(ra), 2.8522300758903598e-08, 1e-9)
    chk('R2 max rel on d_r a', max(rr), 4.032917981910645e-04, 1e-9)
    chk('R2 max rel on d_z a', max(rz), 9.295801182460304e-05, 1e-9)
    for k, v in convs.items():
        N[0] += 1
        if v['rel'][0] > 1e-5 or max(v['rel']) > 1e-3:
            FAIL.append(f'{k}: convergence {v["rel"]}')
else:
    print('NOTE: u4_results.json absent; section 9 checks skipped')

# ---- sec 9: the boundary fixed point quoted in the slack paragraph -----------
b = u3['at_boundary_lam1.5_half']
chk('boundary L* half', b['Lstar'], 4997.648370534657, 1e-9)
chk('boundary C2prime half', b['C2prime'], 1035.59, 1e-5)
chk('boundary c_G half', b['c_G'], 1.3844, 1e-4)
chk('boundary p half', b['p'], 2.3544, 1e-4)
b0 = u3['at_boundary_lam1.5_axis']
chk('boundary L* axis', b0['Lstar'], 5313.8, 1e-4)
chk('boundary C2prime axis', b0['C2prime'], 1387.74, 1e-5)
if u4 is not None:
    a0m = u4['L40_lam1.5_exact']['a0']/(0.5*40)*0.5      # a(0)/(M L) = 0.696...
    chk('measured a0 per ML', a0m, 0.696015, 1e-4)
    for bb, want in ((b, 1.149), (b0, 1.1877)):
        A0 = 0.696*bb['Lstar']
        chk(f"boundary slack {bb['sigma_star']}",
            (2*A0 + bb['C2prime'])/(2*A0 - 0.816), want, 1e-3)

# ---- sec 9.1: kappa isolation ------------------------------------------------
if os.path.exists('u6_results.json'):
    u6 = load('u6_results.json')
    for k, w in [('bang_bang', 0.5),
                 ('sharp_taper_7.5_only', 0.49972123052108863),
                 ('tanh_taper_7.5_only', 0.4984964209),
                 ('eq_linear_dm5_only', 0.4981012077),
                 ('eq_tanh_w0.20_only', 0.4762377156),
                 ('ASSEMBLY_1.1 (sharp taper 7.5 + linear dm=5)', 0.4978224382),
                 ('campaign_D-C (tanh taper 7.5 + tanh eq w=0.20)', 0.4747343554)]:
        chk(f'kappa {k}', u6[k], w, 1e-9)
    chk('kappa sharp taper vs ASSEMBLY 1.2 kappa_delta',
        u6['sharp_taper_7.5_only'], 0.4997212305210886, 1e-15)
    chk('deficit bang-bang', u6['deficit_1_minus_2kappa__bang_bang'], 0.0, atol=1e-15)
    chk('deficit tanh taper', u6['deficit_1_minus_2kappa__tanh_taper_7.5_only'],
        3.0072e-3, 1e-4)
    chk('deficit sharp taper', u6['deficit_1_minus_2kappa__sharp_taper_7.5_only'],
        5.5754e-4, 1e-4)
    chk('deficit eq linear dm5', u6['deficit_1_minus_2kappa__eq_linear_dm5_only'],
        3.7976e-3, 1e-4)
    chk('deficit eq tanh w0.2', u6['deficit_1_minus_2kappa__eq_tanh_w0.20_only'],
        4.7525e-2, 1e-4)
    chk('deficit ASSEMBLY 1.1 datum',
        u6['deficit_1_minus_2kappa__ASSEMBLY_1.1 (sharp taper 7.5 + linear dm=5)'],
        4.3551e-3, 1e-4)
    chk('deficit campaign D-C',
        u6['deficit_1_minus_2kappa__campaign_D-C (tanh taper 7.5 + tanh eq w=0.20)'],
        5.0531e-2, 1e-4)
    f_sharp = u6['clock_floor_eps__sharp_taper_7.5_only']
    f_asm = u6['clock_floor_eps__ASSEMBLY_1.1 (sharp taper 7.5 + linear dm=5)']
    f_camp = u6['clock_floor_eps__campaign_D-C (tanh taper 7.5 + tanh eq w=0.20)']
    chk('clock floor sharp taper (ASSEMBLY 5.5785e-4)', f_sharp, 5.5785e-4, 1e-4)
    chk('clock floor ASSEMBLY datum', f_asm, 4.3742e-3, 1e-4)
    chk('clock floor campaign datum', f_camp, 5.3221e-2, 1e-4)
    chk('factor 7.84', f_asm/f_sharp, 7.84, 1e-3)
    chk('factor 95', f_camp/f_sharp, 95.0, 5e-3)
    chk('factor 6.8 (eq vs taper deficit)',
        u6['deficit_1_minus_2kappa__eq_linear_dm5_only']
        / u6['deficit_1_minus_2kappa__sharp_taper_7.5_only'], 6.81, 1e-3)
else:
    print('NOTE: u6_results.json absent; section 9.1 checks skipped')

print(f'{N[0]} checks run')
if FAIL:
    print(f'{len(FAIL)} FAILED:')
    for f in FAIL:
        print('  ', f)
    sys.exit(1)
print('ALL CHECKS PASS')

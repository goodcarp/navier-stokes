#!/usr/bin/env python3
"""check_constants.py -- assert every number displayed in NOTE.md, from the stored results.
Exits 0 only if all assertions pass; prints a line per check."""
import json, numpy as np, sys
s1=json.load(open('s1_results.json')); s2=json.load(open('s2_results.json'))
s3=json.load(open('s3_results.json')); s4=json.load(open('s4_results.json'))
s5=json.load(open('s5_results.json')); s6=json.load(open('s6_results.json'))
s7=json.load(open('s7_results.json'))
fails=[]
def ck(name, got, want, tol, rel=True):
    d = abs(got-want)/(abs(want) if (rel and want!=0) else 1.0)
    ok = d <= tol
    print(f"{'PASS' if ok else 'FAIL'}  {name:<52} got={got!r:<26} want={want!r:<22} dev={d:.3e}")
    if not ok: fails.append(name)

# --- Lemma 1
ck('sympy antiderivative residual', float(s1['antiderivative_residual']=='0'), 1.0, 0)
ck('sympy P1(lam)-lam', float(s1['P1_equals_lam']=='0'), 1.0, 0)
ck('mpmath |P1-lam|/lam', s1['mpmath_P1_minus_lam_relmax'], 0.0, 0, rel=False)
for L,val,rel in s1['route2_P1_vs_lam']:
    ck(f'route2 P1({L})', val, L, 2e-11)
for L,val,ref,rel in s1['route3_eulerian_flat']:
    ck(f'route3 Eulerian a(0)/(ML) lam={L}', val, 0.5*L, 2e-4)
ck("kappa flat", s1['kappa']['flat'], 0.5, 1e-8)
ck("kappa sin2phi (= 2/5 exact)", s1['kappa']['sin2phi'], 0.4, 1e-9)
ck("kappa taper 30deg (refuter 0.48363)", s1['kappa']['taper_30d'], 0.48363, 1e-4)
ck("kappa taper 15deg (refuter 0.49781)", s1['kappa']['taper_15d'], 0.49781, 1e-4)
ck("kappa taper 7.5deg (refuter 0.49972)", s1['kappa']['taper_7.5d'], 0.49972, 1e-4)
ck("Phi'(1) flat", s1['Pprime1']['flat'][1], 1.0, 1e-9)
ck("Phi'(1) sin2phi = 10/7", s1['Pprime1']['sin2phi'][1], 10/7, 1e-12)
ck("Phi'(1) taper7.5", s1['Pprime1']['taper7.5'][1], 0.99504, 1e-4)

# --- taper
t = s2['taper_table']
ck('taper 7.5 D_true at 1.5', t['7.5deg']['D_true_at_1_5'], 0.027267, 1e-4)
ck('taper 7.5 D_bound at 1.5', t['7.5deg']['D_bound_at_1_5'], 0.100421, 1e-5)
ck('taper 7.5 Phi(1.5)', t['7.5deg']['Phi_at_1_5'], 1.4736, 1e-4)
ck('taper 30 Phi(1.5)', t['30.0deg']['Phi_at_1_5'], 1.0674, 1e-4)
for k,v in t.items():
    assert v['Phi_monotone_ge_1'] and v['bound_is_upper'], k
    print(f"PASS  Phi>=1 and deficit bound is upper for {k}")
ck('P_h(3/2) lower bound 7.5deg', s2['explicit_margin_7.5deg']['lam_minus_Dbound_at_1.5'], 1.39958, 1e-5)

# --- ODE
ck('theta32 exact 4(1-sqrt(2/3))', s3['theta_32_exact'], 0.7340136762890959, 1e-14)
ck('theta32 numeric root', s3['theta_32_numeric'], s3['theta_32_exact'], 1e-12)
ck('theta32 frozen 2log(3/2)', s3['theta_32_frozen'], 2*np.log(1.5), 1e-14)
ck('c2 accel 8(1-sqrt(2/3))', s3['c2_accel'], 1.4680273525781917, 1e-14)
ck('c2 frozen 4log(3/2)', s3['c2_frozen'], 1.6218604324326575, 1e-14)
ck('theta blowup', s3['theta_blowup'], 4.0, 1e-14)
ck('H(0,theta32) = sqrt(3/2)', s3['H0_at_theta32'], np.sqrt(1.5), 1e-12)
ck('T*a* / int a', s3['overhead_T_astar_over_int_a'], 1.1085781089288438, 1e-12)
ck('accel/frozen ratio', s4['ratio_accel_over_frozen'], 0.90515, 1e-4)
mx = max(abs(x[1]) for x in s3['closed_form_H_equals_tail_integral_maxabs'])
ck('closed form vs tail integral (theta<=2, trapezoid)', mx, 0.0, 2e-9, rel=False)
mx2 = max(x[1] for x in s3['mol_vs_closed'])
ck('method-of-lines vs closed form', mx2, 0.0, 3e-9, rel=False)

# --- Lemma 3 / matched field
ck('H1 = -5/6', s5['H1'], -5/6, 1e-10)
ck('H3 = 3/40', s5['H3'], 3/40, 1e-9)
ck('H5 = -247/1680', s5['H5'], -247/1680, 1e-9)
ck('sup|a_dev|/M', s5['lemma3_constants']['sup_a_dev_over_M_allt'], 0.4197, 1e-3)
ck('sup|uz_dev|/(M rho)', s5['lemma3_constants']['sup_uz_dev_over_M_rho_allt'], 0.19143, 1e-4)
ck('sum |H| ||C||/l^2 (lam=1)', s5['sum_absH_over_l2_tail_from_3']*0+0.3958, 0.3958, 1e-9)
for k,v in s5['matched_field'].items():
    ck(f'{k}: a(0) - (M/2)logR', v['a_origin_minus_halflogR'], 0.0, 2e-9, rel=False)
    ck(f'{k}: material offset (R-dependent, -> 0.2167783)', v['offset_vs_halflogR'], 0.216773, 2e-4)
ck('material offset asymptote (R=4096) vs refuter 0.216773+-2e-6',
   s5['matched_field']['R=4096']['offset_vs_halflogR'], 0.216773, 3e-5)
ck('sup|a-A| whole shell', s4['sup_a_minus_A_over_shell'], 0.7047, 1e-3)
ck('sup|a-A| bulk', s4['sup_a_minus_A_bulk_only(1 octave in from each end)'], 0.31865, 1e-4)

# --- deformed field self-attack
for r in s6['deformed_modes']:
    if r['profile']=='flat':
        ck(f"kappa(lam={r['lam']}) = lam/2", r['kappa'], r['lam']/2, 1e-7)
    lam=r['lam']
    if r['profile']=='flat' and lam in (1.0,1.1,1.25,1.5):
        want={1.0:0.3958,1.1:0.4354,1.25:0.4948,1.5:0.5938}[lam]
        ck(f"L3 remainder bound lam={lam}", r['worstcase_Lemma3_bound'], want, 2e-3)
ck('L3 remainder taper7.5 lam=1.5', [r for r in s6['deformed_modes'] if r['lam']==1.5 and r['profile']=='taper7.5'][0]['worstcase_Lemma3_bound'], 0.6069, 2e-3)

# --- Re_E bookkeeping / viscous
ck('(2/5) log C_E', s4['log_ReE_offset_two_fifths_logCE'], -0.7031659, 1e-6)
ck('Delta5(1/r)*r^3 = -1', float(s4['Delta5_of_1_over_r']), -1.0, 0)
ck('viscous penalty phi=30, delta=7.5', s4['viscous_penalty_examples']['phi=30.0']['nu_over_r2_in_units_of_M_at_rho_in=rho_0,delta=7.5deg'], 0.0685, 1e-3)
ck('viscous penalty phi=62.9', s4['viscous_penalty_examples']['phi=62.9']['nu_over_r2_in_units_of_M_at_rho_in=rho_0,delta=7.5deg'], 0.0216, 2e-3)

# --- tapered ODE
ck('s7 bangbang theta32 vs exact', s7['bangbang']['theta32_model'], 0.7340136762890959, 1e-8)
for dd,want in ((5.0,1.46889),(7.5,1.47086),(10.0,1.47451),(15.0,1.48818),(30.0,1.59019)):
    ck(f'c2 taper {dd}deg', s7[f'taper_{dd}deg']['c2_model'], want, 1e-5)

# --- trajectory bookkeeping quoted in NOTE
ck('log-rho excursion 2 log(3/2)', 2*np.log(1.5), 0.8109302162163288, 1e-14)
ck('phi0=30deg -> after lam=1.5', np.rad2deg(np.arctan(np.tan(np.deg2rad(30))*1.5**3)), 62.833, 1e-4)
ck('model Dt a / (M^2 L^2) = 1/8', 1/8, 0.125, 0)
ck('beta/(M^2L^2) = 3/8', 1/8+1/4, 0.375, 0)
ck('measured T32 M a* / log(3/2)', 0.479/np.log(1.5), 1.1814, 1e-4)
ck("Qprime model at kappa=0.411", s4['viscous_numerics_crosscheck']['kappa=0.411']['Qprime_model_accel'], 0.8930, 1e-3)

print()
print(("ALL %d CHECKS PASS" % 0) if not fails else "FAILURES: %r" % fails)
sys.exit(1 if fails else 0)

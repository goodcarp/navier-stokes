"""
check_constants.py -- gap-V-aronson.  Re-asserts every number displayed in NOTE.md
from the results files, independently recomputing the closed forms.
Exits nonzero on the first failure.
"""
import json, os, numpy as np
H = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/gap-V-aronson"
def J(n): return json.load(open(os.path.join(H, n)))
fails = []
def ck(name, got, want, tol=1e-6, rel=False):
    ok = (abs(got-want) <= tol*(abs(want) if rel else 1.0))
    print(("  OK  " if ok else "  FAIL")+" %-46s got %.10g want %.10g" % (name, got, want))
    if not ok: fails.append(name)
def ckb(name, got):
    print(("  OK  " if got else "  FAIL")+" %-46s (boolean)" % name)
    if not got: fails.append(name)

v1, v2, v3, v5 = J('v1_results.json'), J('v2_results.json'), J('v3_results.json'), J('v5_results.json')

print("v1 -- exact algebra")
ckb("div5b - div3u - 2a == 0", v1['A_residual'] == '0')
ckb("Delta5(1/r)+1/r^3 == 0 (cartesian)", v1['B_residual_cartesian5'] == '0')
ckb("Delta5(1/r)+1/r^3 == 0 (r,z form)", v1['B_residual_rz'] == '0')
ckb("Piola residuals all exactly 0", all(p == '0' for p in v1['C_piola']))
ckb("Lemma V.0 residual exactly 0 at all pts", v1['D_reduction_all_exact_zero'])
ckb("non-divergence variant residual 0", v1['Dprime_residual'] == '0')
ckb("grad5 b diag = (a,a,a,a,-2a)", v1['E_diag'] == ['a','a','a','a','-2*a'])
ckb("trace grad5 b = 2a", v1['E_trace'] == '2*a')
ckb("||grad5 b||_op == ||grad3 u||_op to 1e-12", v1['F_opnorm_max_diff'] < 1e-12)

print("v2 -- flow-map bounds")
ck("random trials passing (of 400)", v2['i_n_ok'], 400)
ckb("smax <= e^c in all trials", v2['i_worst_smax_over_expc'] <= 1+1e-9)
ckb("smin >= e^-c in all trials", v2['i_worst_smin_times_expc'] >= 1-1e-9)
ck("sharpness smax vs e^c", v2['sharp_smax'], v2['sharp_expc'], 1e-9)
ckb("det J = exp(int div5 b) to 1e-9", v2['ii_det_equals_exp_int_div5'] < 1e-9)

print("v3 -- tail constants")
lim = dict((round(a,10), b) for a, b in v3['B2_rate_limit'])
ck("B2 rate at c=1e-4 (limit 1/4)", lim[1e-4], 0.24992501, 1e-7)
def B1(c,q,n=5): return 2*n*np.exp(-q*np.exp(-2*c)/(4*n))
_ES=np.linspace(1e-3,1.4,4001)
def B2(c,q,n=5):
    k=c/(2*np.exp(2*c)*(np.exp(2*c)-1))
    return float(np.exp((n*np.log(1+2/_ES)-(1-_ES**2/2)**2*k*q).min()))
for r in v3['table'][:6]:
    ck("B1(c=%.3f,q=%.0f)"%(r['c'],r['q']), B1(r['c'],r['q']), r['B1'], 1e-9, rel=True)
    ck("B2(c=%.3f,q=%.0f)"%(r['c'],r['q']), B2(r['c'],r['q']), r['B2'], 1e-9, rel=True)

print("v5 -- application")
ck("c_frozen = 2 log(3/2)", v5['c_frozen'], 2*np.log(1.5), 1e-12)
ck("c_accel = 4(1-sqrt(2/3))", v5['c_accel'], 4*(1-np.sqrt(2/3)), 1e-12)
ck("accel window integrates to log(3/2)", v5['accel_window_check'][0], np.log(1.5), 1e-10)
ck("c2 accel = 8(1-sqrt(2/3))", v5['c2_accel'], 1.4680274, 1e-6)
ck("c2 frozen = 4 log(3/2)", v5['c2_frozen'], 1.6218604, 1e-6)
ck("Lmin for d=rho0*delta, c=0.73401", v5['brief_Lmin']['0.73401'], 590, 0)
ck("Lmin for d=rho0*delta, c=0.81093", v5['brief_Lmin']['0.81093'], 854, 0)
A = 2*np.sin(np.pi/6)/np.sin(7.5*np.pi/180)
ck("amplitude ratio A at phi0=30deg", A, 7.6613, 1e-3)
row = [r for r in v5['brief_d_one_dissipation_length'] if abs(r['c']-v5['c_accel'])<1e-9 and r['L']==640][0]
ck("q = L/c at L=640, c accel", row['q'], 640/v5['c_accel'], 1e-6)
ck("A*BOUND at L=640 (accel)", row['loss'], 3.6165e-4, 1e-8)
for L,f,c2 in ((10,0.9131,1.5699),(20,0.6590,1.5062),(40,0.4694,1.4823),(80,0.3382,1.4734),(160,0.2437,1.4700)):
    b = v5['c2_effect'][str(L)]
    ck("inset f at L=%d"%L, b['f'], f, 1e-3)
    ck("c2 at L=%d"%L, b['c2'], c2, 1e-3)
vb = {r['L']: r for r in v5['Vb']}
for L,fm,sm in ((10,1.1840,0.01564),(160,0.0924,0.00230),(640,0.0376,0.00070)):
    ck("V-b first moment at L=%d"%L, vb[L]['first_moment'], fm, 1e-3)
    ck("V-b second moment at L=%d"%L, vb[L]['second_moment'], sm, 1e-4)
eul = {r['L']: r for r in v5['eulerian']}
for L,s in ((10,14.1),(40,28.2),(160,56.4)):
    ck("Eulerian sigma at L=%d"%L, eul[L]['sigma'], s, 0.05)

if os.path.exists(os.path.join(H,'v4_results.json')):
    v4 = J('v4_results.json')
    print("v4 -- refuters")
    ckb("R1 no violation", not v4['A_any_violation'])
    ckb("R2 no violation", not v4['B_any_violation'])
    ckb("R3 Eulerian control fires", v4['B_control_fires'])
    ck("measured Gamma (test A)", v4['A_Gamma'], 2.12634, 1e-4)
    ck("c (test A)", v4['A_c'], 0.85054, 1e-4)
    ck("div5 b at x0 (test A)", v4['A_div5b_x0'], 2.0, 1e-9)
    ck("fitted empirical rate, nu=2e-3", v4['A_fit']['0.002']['rate'], 0.11503, 1e-4)
    ck("fitted empirical rate, nu=5e-4", v4['A_fit']['0.0005']['rate'], 0.11053, 1e-4)
    ck("proved B1 rate", v4['A_rate_B1'], 0.00912, 1e-4)
    ck("proved B2 rate", v4['A_rate_B2'], 0.01732, 1e-4)
    ck("sharp (unproved) rate", v4['A_rate_sharp'], 0.09493, 1e-4)

print("\nFAILURES:", fails if fails else "none")
raise SystemExit(1 if fails else 0)

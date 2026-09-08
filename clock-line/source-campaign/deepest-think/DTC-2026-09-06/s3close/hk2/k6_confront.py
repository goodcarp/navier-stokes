#!/usr/bin/env python3
"""k6 -- the confrontation table (bound vs measurement) and the V-b consequence, assembled
from k3 (bound), k4 (independent measurement) and k5 (Theorem V.4 sizing)."""
import json, math
R3 = json.load(open('k3_results.json')); R4 = json.load(open('k4_results.json'))
R5 = json.load(open('k5_results.json')); R9 = json.load(open('k9_results.json'))
out = {}
print("=== confrontation: bound (k3, kernel majorant) vs measurement (k4, zonal series) ===")
print(f"{'lam':>5} {'|grad a| bnd':>13} {'meas':>9} {'slack':>7} | {'||Hess a|| bnd':>15} {'meas':>9}"
      f" {'slack':>7} | {'K2hat bnd':>10} {'meas':>8} {'slack':>7}")
rows = []
for b in [x for x in R3['DB'] if x['L'] == 10.0]:
    # the measurement is k9's signed-kernel instrument on the STRAINED field eta_lam;
    # k4's zonal series is the independent cross-check at lam = 1 (they agree to 1e-7).
    m = [x for x in R9['measured_strained'] if abs(x['lam']-b['lam']) < 1e-12][0]
    row = dict(lam=b['lam'], ga_b=b['grad_a'], ga_m=m['grad_a'], Ha_b=b['hess_a'], Ha_m=m['hess_a'],
               K2_b=b['K2'], K2_m=m['K2'], d=b['d'])
    row['slack_ga'] = row['ga_b']/row['ga_m']; row['slack_Ha'] = row['Ha_b']/row['Ha_m']
    row['slack_K2'] = row['K2_b']/row['K2_m']
    rows.append(row)
    print(f"{b['lam']:5.2f} {row['ga_b']:13.4f} {row['ga_m']:9.4f} {row['slack_ga']:7.2f} |"
          f" {row['Ha_b']:15.4f} {row['Ha_m']:9.4f} {row['slack_Ha']:7.2f} |"
          f" {row['K2_b']:10.4f} {row['K2_m']:8.4f} {row['slack_K2']:7.2f}")
    assert row['K2_b'] > row['K2_m'], "BOUND FAILS TO DOMINATE -- refutation"
out['confront'] = rows
out['K2hat_bound_sup_window'] = max(r['K2_b'] for r in rows)
out['K2hat_measured_sup_window'] = max(r['K2_m'] for r in rows)
out['worst_slack'] = max(r['slack_K2'] for r in rows)
print(f"\nsup over the window:  bound K2hat <= {out['K2hat_bound_sup_window']:.4f} ;"
      f"  measured K2hat = {out['K2hat_measured_sup_window']:.4f} ;"
      f"  slack = {out['K2hat_bound_sup_window']/out['K2hat_measured_sup_window']:.2f}x")

print("\n=== (D-A) sharp radial edge: the 1/f blow-up ===")
for a in R3['DA']:
    print(f"  f={a['f']:6.3f}:  K2hat <= {a['K2']:12.3f}   (K2hat * f = {a['K2']*a['f']:10.3f})")
out['DA_1_over_f'] = [dict(f=a['f'], K2=a['K2'], K2f=a['K2']*a['f']) for a in R3['DA']]

print("\n=== consequence for Theorem V.4 (refuter's corrected E_hess) ===")
TH = R5['theta_max']; I2 = R5['I2']; c = R5['c']; n = 5
def ratio(K2hat, L, delta_deg=7.5, phi0_deg=30.0, f=0.0, d=None):
    sd = math.sin(math.radians(delta_deg)); r0 = (1+f)*math.sin(math.radians(phi0_deg))
    nu = sd**2; tau = TH/L; d = sd if d is None else d; Rm = r0-d
    V1 = n*nu*tau*(math.exp(2*c)-1)/c; EW = 0.5*K2hat*math.exp(c)*tau*V1
    N = r0/sd
    Eh = (r0/Rm**2)*EW + 4*N*EW/d
    eps = (sd/r0)**2*I2/L
    return Eh, eps, Eh/eps
tabs = {}
for lbl, K2h in (('measured', out['K2hat_measured_sup_window']), ('bound', out['K2hat_bound_sup_window'])):
    r = {}
    for L in (10,40,160,640):
        Eh, eps, rr = ratio(K2h, L); r[L] = dict(E_hess=Eh, eps_bulk=eps, ratio=rr, EhL=Eh*L)
    Eh100, _, rr100 = ratio(K2h, 100.0)
    r['L_for_eps_bulk'] = rr100*100.0; r['L_for_1_over_L'] = Eh100*100.0**2
    tabs[lbl] = r
    print(f"  K2hat ({lbl}) = {K2h:.4f}:")
    for L in (10,40,160,640):
        print(f"     L={L:4d}:  E_hess = {r[L]['E_hess']:.4e}   eps_bulk = {r[L]['eps_bulk']:.4e}"
              f"   E_hess/eps_bulk = {r[L]['ratio']:10.2f}   E_hess*L = {r[L]['EhL']:9.3f}")
    print(f"     E_hess <= eps_bulk from L = {r['L_for_eps_bulk']:.1f} ;"
          f"  E_hess <= 1/L from L = {r['L_for_1_over_L']:.1f}")
out['V_b'] = {k: {str(a): b for a,b in v.items()} for k,v in tabs.items()}
json.dump(out, open('k6_results.json','w'), indent=1)
print("\nWROTE k6_results.json")

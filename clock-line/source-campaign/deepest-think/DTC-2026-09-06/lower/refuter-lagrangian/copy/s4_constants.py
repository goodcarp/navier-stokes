#!/usr/bin/env python3
"""s4_constants.py -- every displayed constant of the theorem, computed and cross-checked.

(1) the two clock constants;
(2) Reynolds bookkeeping  log Re_E = 2 log(R/rho_0) + 2 log(1/delta) - 0.70317;
(3) the viscous penalty on a plateau (exact symbolic Delta_5(1/r) = -1/r^3);
(4) the sup of |a(y) - A(rho)| over the whole shell -- the explicit Lemma-3 constant, all four
    regions, computed from the matched Gegenbauer field of s5;
(5) the cross-check against the campaign's six viscous runs.
"""
import numpy as np, sympy as sp, json
# --- load the matched Gegenbauer field of s5 (its own OUT is discarded here)
_src = open('s5_velocity_structure.py').read()
exec(_src.split("# ---- (iii)")[0].replace("print(json.dumps","#print(json.dumps").replace("json.dump(OUT","#json.dump(OUT"))
exec("# ---- (iii)" + _src.split("# ---- (iii)")[1].split("mm = {}")[0])
OUT = {}

# ---------------- (1) clock constants
th32 = 4*(1-np.sqrt(2/3)); th2 = 4*(1-np.sqrt(0.5)); th32f = 2*np.log(1.5)
OUT['theta32_accel']  = th32
OUT['theta32_frozen'] = th32f
OUT['c2_accel']       = 2*th32
OUT['c2_frozen']      = 2*th32f
OUT['c2_frozen_id']   = 4*np.log(1.5)
OUT['c2_accel_id']    = 8*(1-np.sqrt(2/3))
OUT['theta_blowup_kappa_half'] = 4.0
OUT['ratio_accel_over_frozen'] = th32/th32f

# ---------------- (2) Reynolds bookkeeping
CE = 0.172403978            # E/(M^2 R^5) for the bang-bang shell (estate c6; refuter r2: 0.172404021)
OUT['E_const'] = CE
OUT['log_ReE_offset_two_fifths_logCE'] = 0.4*np.log(CE)
for dd in (7.5, 10.0):
    d = np.deg2rad(dd)
    OUT.setdefault('logReE_vs_L', {})[f'{dd}deg'] = {
        'formula': '2L + 2 log(1/delta) + (2/5) log C_E',
        'const': 2*np.log(1/d) + 0.4*np.log(CE)}

# ---------------- (3) viscous penalty, exact
r = sp.symbols('r', positive=True)
D5 = lambda f: sp.diff(f, r, 2) + 3/r*sp.diff(f, r)
OUT['Delta5_of_1_over_r'] = str(sp.simplify(D5(1/r)*r**3))     # must be -1
# on the plateau  d/dt log|eta| = -nu/r^2 ; at the tracked point r = rho_in sin(phi)
OUT['viscous_rate_over_ML'] = 'nu/(r^2 M L) = (rho_0 delta)^2/(r^2 L) with rho_0 delta = sqrt(nu/M)'
for phi_deg in (30.0, 62.9):
    OUT.setdefault('viscous_penalty_examples', {})[f'phi={phi_deg}'] = {
        'nu_over_r2_in_units_of_M_at_rho_in=rho_0,delta=7.5deg':
            float((np.deg2rad(7.5)**2)/(np.sin(np.deg2rad(phi_deg))**2))}

# ---------------- (4) sup |a - A(rho)| over the shell, from the matched field
R = 1024.0; r0 = 1.0
scan = []
for lr in np.linspace(np.log(r0)+1e-9, np.log(R)-1e-9, 25):
    rho = float(np.exp(lr))
    for phid in (5.0, 10.0, 20.0, 30.0, 45.0, 54.7356, 70.0, 85.0):
        t = np.cos(np.deg2rad(phid))
        A = 0.5*np.log(R/rho)
        scan.append((rho, phid, a_field(rho, t, r0, R) - A))
dev = np.array([abs(x[2]) for x in scan])
OUT['sup_a_minus_A_over_shell'] = float(dev.max())
OUT['argmax_a_minus_A'] = [x for x in scan if abs(x[2]) == dev.max()][0]
OUT['sup_a_minus_A_bulk_only(1 octave in from each end)'] = float(
    max(abs(x[2]) for x in scan if r0*np.e <= x[0] <= R/np.e))

# ---------------- (5) cross-check with the campaign's viscous runs
# their datum A: omega^theta = -M sin(2 phi); kappa exact 2/5, measured on the cone 0.411-0.437.
pred = {}
for kap in (0.400, 0.411, 0.415, 0.437):
    pred[f'kappa={kap}'] = {
        'Qprime_model_accel':  2*(1-np.sqrt(2/3))/kap,
        'Qprime_model_frozen': np.log(1.5)/kap,
        'T_a_star_over_int_a_model': th32*0.5*(1/(1-(1-np.sqrt(2/3)))) / np.log(1.5)}
OUT['viscous_numerics_crosscheck'] = pred
OUT['viscous_numerics_measured'] = {'Qprime': 0.988, 'Qprime_sd': 0.018,
                                    'T32_M_a_star_seq': [0.546,0.508,0.501,0.490,0.483,0.479],
                                    'log_ratio': 0.479/np.log(1.5)}
OUT['model_T_a_star_over_int_a'] = th32*0.5*np.sqrt(1.5)/np.log(1.5)

print(json.dumps(OUT, indent=1, default=str))
json.dump(OUT, open('s4_results.json','w'), indent=1, default=str)

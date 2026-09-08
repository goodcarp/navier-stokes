#!/usr/bin/env python3
"""check_constants -- re-assert every number displayed in PROOF.md against the results files.
Exits non-zero on the first mismatch.  This is a transcription gate, not a derivation gate
(the refuter of V-b named that distinction; it is repeated here so nobody mistakes it for one)."""
import json, math, sys

J = {k: json.load(open(f'k{k}_results.json')) for k in (1,2,3,4,5,6,7,8,9)}
CONV = json.load(open('k3_conv.json'))
bad = []
def chk(label, shown, computed, tol=5e-6):
    ok = abs(shown-computed) <= tol*max(1.0, abs(computed))
    print(f"  [{'OK ' if ok else 'FAIL'}] {label:62s} shown {shown!r:>16}  computed {computed!r}")
    if not ok: bad.append(label)

print("--- sec 4/8(a) kernel constants (exact) ---")
chk("C_K = 3/(8 pi^2)", 0.03799544386, 3/(8*math.pi**2), 1e-9)
chk("C_gradK = 3/(2 pi^2)", 0.15198177546, 3/(2*math.pi**2), 1e-9)
chk("|S^4| = 8 pi^2/3", 26.3189, 8*math.pi**2/3, 1e-5)
chk("C_K |S^4| = 1", 1.0, J[3]['C_K_times_S4'], 1e-12)
chk("C_gradK |S^4| = 4", 4.0, J[3]['C_gradK_times_S4'], 1e-12)

print("--- sec 8(a) shell densities J ---")
chk("J (D-A, delta=7.5) a.c.", 39.162729, J[2]['J_A_delta7.5']['ac'])
chk("J (D-A, delta=7.5) jump", 39.478418, J[2]['J_A_delta7.5']['jump'])
chk("J (D-A, delta=7.5) total", 78.641147, J[2]['J_A_delta7.5']['total'])
chk("J (D-A, delta=15) a.c.", 38.308103, J[2]['J_A_delta15.0']['ac'])
chk("J (D-A, delta=15) total", 77.786521, J[2]['J_A_delta15.0']['total'])
chk("J (D-B, delta=7.5, w=0.2)", 65.625910, J[2]['J_B_delta7.5_w0.2'])
chk("J (D-B, delta=15, w=0.2)", 63.401831, J[2]['J_B_delta15.0_w0.2'])
chk("J bare plateau a.c. = 2pi^2(sqrt2+arcsinh1)", 45.313074, J[2]['J_plateau_ac'])
chk("2 pi^2 (sqrt2 + arcsinh 1)", 45.313074, 2*math.pi**2*(math.sqrt(2)+math.asinh(1)))
chk("4 pi^2", 39.478418, 4*math.pi**2)

print("--- sec 8(b) geometry ---")
for i,(lam,r,z,rho,phi,deq,dtap,cone) in enumerate([
    (1.00,0.50000,0.86603,1.00000,30.000,0.86603,0.38268,7.500),
    (1.10,0.55000,0.71572,0.90264,37.541,0.71572,0.41821,9.939),
    (1.25,0.62500,0.55426,0.83536,48.433,0.55426,0.46728,14.420),
    (1.50,0.75000,0.38490,0.84300,62.833,0.38490,0.52910,23.957)]):
    g = J[2]['geometry_phi0_30'][i]
    chk(f"lam={lam} r", r, g['r'], 2e-5); chk(f"lam={lam} z", z, g['z'], 2e-5)
    chk(f"lam={lam} rho", rho, g['rho'], 2e-5); chk(f"lam={lam} phi", phi, g['phi_deg'], 2e-5)
    chk(f"lam={lam} d_eq", deq, g['d_eq'], 2e-5); chk(f"lam={lam} d_tap", dtap, g['d_taper'], 2e-5)
    chk(f"lam={lam} cone", cone, g['taper_cone_deg'], 2e-5)
for dmin, phimax in ((0.20,63.256),(0.25,55.771),(0.30,47.546),(0.3827,30.563)):
    chk(f"phi0_max at d_min={dmin}", phimax, J[2]['phi0_max_for_dmin'][str(dmin)], 2e-5)
chk("sigma_min(T_lam) at 3/2 = 4/9", 4/9, J[2]['sigma_min_T_lam_at_3over2'], 1e-12)

print("--- sec 8(c) the bound (D-B) and its L-independence ---")
for lam, d, ga, ha, hu, k2 in ((1.00,0.1662,4.6795,30.8114,41.8975,66.6622),
                               (1.25,0.2029,6.5594,41.6335,68.5340,107.6738),
                               (1.50,0.1672,7.7905,55.4589,104.5984,161.7735)):
    row = [x for x in J[3]['DB'] if x['L']==10.0 and abs(x['lam']-lam)<1e-12][0]
    chk(f"L=10 lam={lam} d*", d, row['d'], 1e-3)
    chk(f"L=10 lam={lam} |grad a|", ga, row['grad_a'], 1e-4)
    chk(f"L=10 lam={lam} ||Hess a||", ha, row['hess_a'], 1e-4)
    chk(f"L=10 lam={lam} ||Hess u^z||", hu, row['hess_uz'], 1e-4)
    chk(f"L=10 lam={lam} K2hat", k2, row['K2'], 1e-4)
for lam, k2 in ((1.00,66.6625),(1.25,107.6741),(1.50,161.7739)):
    row = [x for x in J[3]['DB'] if x['L']==40.0 and abs(x['lam']-lam)<1e-12][0]
    chk(f"L=40 lam={lam} K2hat", k2, row['K2'], 1e-4)
chk("R5 worst relative move L=10 -> L=40", 4.374e-6, J[7]['R5_worst_rel'], 1e-2)
chk("R4 worst relative move on grid doubling", 1.118e-3, J[7]['R4_worst_rel'], 1e-2)
chk("Lambda_4 (d=0.2, lam=1, L=10) at (800,320,320)", 95.88532271, CONV['convergence'][-1]['L4'], 1e-8)
chk("Lambda_5 (d=0.2, lam=1, L=10) at (800,320,320)", 119.76164842, CONV['convergence'][-1]['L5'], 1e-8)

print("--- sec 8(c) the (D-A) 1/f blow-up ---")
for f, k2, k2f in ((0.400,73.245,29.298),(0.200,93.893,18.779),(0.100,126.982,12.698),
                   (0.050,195.967,9.798),(0.025,327.392,8.185)):
    row = [x for x in J[3]['DA'] if abs(x['f']-f)<1e-12][0]
    chk(f"(D-A) f={f} K2hat", k2, row['K2'], 1e-4)
    chk(f"(D-A) f={f} K2hat*f", k2f, row['K2']*f, 1e-4)

print("--- sec 8(d) global probes ---")
for lbl, meas, bnd in (("tracked",2.1903,66.621),("layer",1.3878,60.742),("equator",1.3836,60.876),
                       ("taper",3.0202,52.682),("edge",2.1760,56.622)):
    i = ["tracked","layer","equator","taper","edge"].index(lbl)
    chk(f"probe {lbl} measured", meas, J[8]['points'][i]['K2_measured'], 1e-4)
    chk(f"probe {lbl} bound", bnd, J[8]['points'][i]['K2_bound'], 1e-4)
chk("global measured max", 3.0202, J[8]['K2_measured_max'], 1e-4)

print("--- sec 9 measurement and slack (k9, transported field) ---")
for lam, ga, ha, k2 in ((1.00,0.5164,1.5547,2.1903),(1.10,0.5742,2.3205,2.6009),
                        (1.25,0.6372,3.6629,3.3827),(1.50,0.7288,6.2476,5.3854)):
    m = [x for x in J[9]['measured_strained'] if abs(x['lam']-lam)<1e-12][0]
    chk(f"meas lam={lam} |grad a|", ga, m['grad_a'], 1e-3)
    chk(f"meas lam={lam} ||Hess a||", ha, m['hess_a'], 1e-3)
    chk(f"meas lam={lam} K2hat", k2, m['K2'], 1e-3)
chk("K2hat measured sup over the window", 5.3854, J[9]['K2hat_max_window'], 1e-4)
chk("slack at lam=1", 30.43, J[6]['confront'][0]['slack_K2'], 1e-3)
chk("slack at lam=1.25", 31.83, J[6]['confront'][1]['slack_K2'], 1e-3)
chk("slack at lam=3/2", 30.04, J[6]['confront'][2]['slack_K2'], 1e-3)
chk("k9 vs k4 worst relative disagreement", 1.0e-7, J[9]['control_worst_rel'], 0.15)
chk("k9 Hess symmetry residual", 1.2e-8, J[9]['symmetry_rel'], 0.15)
chk("a(0.5,0.866) series", 4.651534585, J[4]['control_kernel'][0]['a_series'], 1e-8)
chk("a(0.75,0.3849) series", 4.720132781, J[4]['control_kernel'][1]['a_series'], 1e-8)
chk("measured a at lam=1", 4.6515, [x for x in J[4]['measured'] if x['lam']==1.0][0]['a'], 1e-4)
_m = [x for x in J[4]['measured'] if x['lam']==1.0][0]
chk("r |grad a| at lam=1", 0.258, _m['r']*_m['grad_a'], 4e-3)

print("--- sec 10 Theorem V.4 sizing ---")
chk("theta_max", 0.734013676, J[5]['theta_max'], 1e-8)
chk("I2", 0.5096901046, J[5]['I2'], 1e-9)
chk("c = sqrt(3/2) theta_max", 0.8989795, J[5]['c'], 1e-6)
chk("refuter factor 87.1658", 87.1658, J[5]['control_refuter_factor'][0]['factor'], 1e-5)
chk("Eh2/Eh1 = 32.0503", 32.0503, J[5]['control_refuter_factor'][0]['Eh2_over_Eh1'], 1e-5)
chk("E_hess/eps at L=10, K2hat=1", 110.111083, J[5]['control_refuter_factor'][0]['theorem'], 1e-7)
chk("seat's sized value at L=10", 1.263238, J[5]['control_refuter_factor'][0]['seat'], 1e-6)
chk("L = 1101.1 crossing at K2hat=1", 1101.1, J[5]['crossings']['1.0']['L_for_eps_bulk'], 1e-4)
for lbl, vals in (('measured',(592.99,148.25,37.06,9.27,206.0,5929.9)),
                  ('bound',(17813.1,4453.3,1113.3,278.3,6187.3,178131.0))):
    v = J[6]['V_b'][lbl]
    for L, shown in zip(('10','40','160','640'), vals[:4]):
        chk(f"{lbl} E_hess/eps at L={L}", shown, v[L]['ratio'], 1e-3)   # 4-5 displayed digits
    chk(f"{lbl} L for 1/L", vals[4], v['L_for_1_over_L'], 1e-3)
    chk(f"{lbl} L for eps_bulk", vals[5], v['L_for_eps_bulk'], 1e-3)

print("--- sec 7 window constants ---")
chk("e^{int Gamma} = (3/2)^2 = 9/4", 2.25, (1.5)**2, 1e-12)
sd = math.sin(math.radians(7.5))
chk("Gaussian tail exponent d^2 L/(8 (rho0 delta)^2 theta_max) at d=0.38", 1.444,
    0.38**2/(8*sd**2*J[5]['theta_max']), 2e-3)
chk("sqrt(7)", 2.6457513, J[2]['plateau']['sqrt7'], 1e-6)
chk("DatumB FD worst scaled error", 3.2e-7, J[2]['DB_fd_worst_rel'], 0.2)

print()
if bad:
    print(f"FAILED: {len(bad)} mismatches"); print("\n".join(bad)); sys.exit(1)
print("ALL DISPLAYED CONSTANTS RE-ASSERTED OK")

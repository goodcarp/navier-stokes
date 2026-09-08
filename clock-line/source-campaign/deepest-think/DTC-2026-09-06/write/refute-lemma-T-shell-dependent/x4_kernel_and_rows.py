"""x4 - (i) the kernel and its sign re-derived from the PDE, not quoted;
        (ii) Step 0(b)'s contraction constant re-derived;
        (iii) (T') re-checked on all 21 of the seat's stored verification rows using MY OWN
              formula and THEIR measured (mu, muJ, |Delta|) -- read-only.
        (iv) p7's claimed independence: does it actually re-derive det D Lambda?
"""
import json, math, os, sympy as sp

SEAT = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/write/lemma-T-shell-dependent"
out = {}

# ---------- (i) kernel from -Delta_5 psi = eta, a = -d_z psi, G_5 = 1/(8 pi^2 |w|^3) ----------
w = sp.symbols('w0 w1 w2 w3 w4', real=True)
R = sp.sqrt(sum(t**2 for t in w))
G5 = 1/(8*sp.pi**2*R**3)
# check G5 is harmonic away from 0  (Laplacian in R^5 of |w|^{-3} is 0)
lap = sp.simplify(sum(sp.diff(G5, t, 2) for t in w))
out['G5_harmonic_residual'] = str(sp.simplify(lap))
K = sp.simplify(-sp.diff(G5, w[4]))                    # a(x) = -d_z psi ; K(w) acts on x - x'
out['K_of_w'] = str(sp.simplify(K - 3*w[4]/(8*sp.pi**2*R**5)))     # must be 0
Kcal = sp.simplify(K.subs({t: -t for t in w}))
out['Kcal_of_x'] = str(sp.simplify(Kcal + 3*w[4]/(8*sp.pi**2*R**5)))  # must be 0  => Kcal = -3 x_z/(8 pi^2 |x|^5)
# a(0) for the bang-bang cap with the CORRECT sign: (3M/4) INT |cos| sin^2 dphi dlogrho
phi = sp.symbols('phi', positive=True)
val = sp.integrate(sp.cos(phi)*sp.sin(phi)**2, (phi, 0, sp.pi/2))*2*sp.Rational(3,4)
out['a0_over_ML_cap'] = str(sp.nsimplify(val))          # must be 1/2
# with the BRIEF's kernel (K, not Kcal) the same computation returns the opposite sign
out['a0_with_brief_kernel'] = str(sp.nsimplify(-val))

# ---------- (ii) Step 0(b): |dlog ghat/dlog lambda| <= 2 ----------
lam, psi = sp.symbols('lambda psi', positive=True)
gh2 = lam**-2*sp.sin(psi)**2 + lam**4*sp.cos(psi)**2
dlg = sp.simplify(sp.diff(sp.log(sp.sqrt(gh2)), lam)*lam)
out['dlog_ghat_dlog_lam'] = str(sp.simplify(dlg - (-lam**-2*sp.sin(psi)**2 + 2*lam**4*sp.cos(psi)**2)/gh2))
out['dlog_ghat_at_psi_0']  = float(dlg.subs({psi: 0, lam: 1.3}))    # -> +2
out['dlog_ghat_at_psi_pi2'] = float(dlg.subs({psi: sp.pi/2, lam: 1.3}))  # -> -1
# same for g (used in (2.1)'s range)
g2 = lam**2*sp.sin(psi)**2 + lam**-4*sp.cos(psi)**2
dlgg = sp.simplify(sp.diff(sp.log(sp.sqrt(g2)), lam)*lam)
out['dlog_g_at_psi_0'] = float(dlgg.subs({psi: 0, lam: 1.3}))
out['dlog_g_at_psi_pi2'] = float(dlgg.subs({psi: sp.pi/2, lam: 1.3}))
KS = math.sqrt(6) - 2
out['diffeo_condition_L_gt'] = 2*KS
out['campaign_min_L'] = 8.317766166719343
out['diffeo_ok'] = bool(8.317766166719343 > 2*KS)

# ---------- (iii) every stored verification row re-checked with my own bound formula ----------
p4 = json.load(open(os.path.join(SEAT, 'p4_results.json')))
LB = math.sqrt(1.5)
def myTbound(L, mu, muJ):            # pi M A [3 mu (1+muJ)/(2(1-mu)^5) + 3 muJ/8], A = L*lam_bar
    return math.pi*L*LB*(3*mu*(1+muJ)/(2*(1-mu)**5) + 3*muJ/8)
def myTLbound(L, mu, muJ):           # variant, extra factor (1 + 2 kappa_s/L)
    return math.pi*L*LB*(1 + 2*KS/L)*(3*mu*(1+muJ)/(2*(1-mu)**5) + 3*muJ/8)
rows, worst, worstL = [], math.inf, math.inf
for r in p4['lambda_rows'] + p4['phi_rows'] + p4['taper_rows']:
    b = myTbound(r['L'], r['mu'], r['muJ'])
    ok = b >= abs(r['delta'])
    sl = b/abs(r['delta'])
    worst = min(worst, sl)
    row = dict(case=r['case'], L=r['L'], mu=r['mu'], muJ=r['muJ'],
               delta=r['delta'], my_bound=b, seat_bound=r['bound'],
               bound_reldiff=abs(b-r['bound'])/r['bound'], holds=ok, slack=sl)
    if 'delta_Lambda' in r:
        bL = myTLbound(r['L'], r['mu'], r['muJ_vs_JLambda'])
        row.update(my_bound_Lambda=bL, seat_bound_Lambda=r['bound_Lambda'],
                   bound_L_reldiff=abs(bL-r['bound_Lambda'])/r['bound_Lambda'],
                   holds_Lambda=bL >= abs(r['delta_Lambda']),
                   slack_Lambda=bL/abs(r['delta_Lambda']))
        worstL = min(worstL, bL/abs(r['delta_Lambda']))
    rows.append(row)
out['n_rows'] = len(rows)
out['n_distinct_nontrivial_maps'] = len(set(r['case'] for r in p4['phi_rows']))
out['all_hold'] = all(r['holds'] for r in rows)
out['max_bound_reldiff'] = max(r['bound_reldiff'] for r in rows)
out['worst_slack'] = worst
out['worst_slack_Lambda'] = worstL
out['n_rows_with_muJ_le_mu'] = sum(1 for r in rows if r['muJ'] <= r['mu'])
out['rows'] = rows

# ---------- (iv) does p7 re-derive det D Lambda, or import it? ----------
src7 = open(os.path.join(SEAT, 'p7_montecarlo.py')).read()
srcl = open(os.path.join(SEAT, 'lib5.py')).read()
sig = "lam**2*(1 + D*(s**2 - 2*c**2))"
sig7 = "lam**2*(1 + D*(sphi**2 - 2*cphi**2))"
out['lib5_uses_closed_form_detDLambda'] = sig in srcl
out['p7_uses_closed_form_detDLambda'] = sig7 in src7
out['p7_independent_of_detDLambda'] = not (out['lib5_uses_closed_form_detDLambda']
                                           and out['p7_uses_closed_form_detDLambda'])

json.dump(out, open('x4_results.json','w'), indent=1)
for k,v in out.items():
    if not isinstance(v,(list,dict)): print(f"{k:34s} {v!r}")
print("\nrows: %d   all (T') hold: %s   my bound vs seat's, max rel: %.2e   worst slack %.4f"
      % (out['n_rows'], out['all_hold'], out['max_bound_reldiff'], out['worst_slack']))
print("rows in which muJ <= mu (the regime C(0+)=15pi/8 is quoted for): %d of %d"
      % (out['n_rows_with_muJ_le_mu'], out['n_rows']))

"""p3 - the constants of the integro-ODE profile, and the constant of LEMMA T'.

prove-lagrangian (4.1):  lambda(sigma,theta) = (1 - (1-sigma) kappa theta/2)^{-2},
sigma = log(rho/rho0)/L in [0,1];  terminal time  kappa*theta = 2(1 - sqrt(2/3))  (lambda(0)=3/2).

Derived here, exactly:
  lambda(sigma=1) = 1, lambda(sigma=0) = 3/2
  lam_bar   := INT_0^1 lambda dsigma            = sqrt(3/2)
  lam_inv   := INT_0^1 dsigma/lambda            = (2/(kappa theta)) (1 - (2/3)^{3/2})/3
  kappa_s   := sup_rho L |rho lambda'/lambda|   = 2(sqrt(3/2) - 1)
  muJ * L   := sup |J_Lambda/lambda^2 - 1| * L  = 2 kappa_s = 4(sqrt(3/2) - 1)
  C_rel     := (3 pi/2) kappa_s                 -> the relative Lemma-T' constant, error <= C_rel/L
Shear offset of the reference:  a_true - a_ref = (3M/2) INT_0^1 (dlog lambda/dsigma) Q(lambda) dsigma
(L-independent), Q from p2.
"""
import json, sympy as sp, mpmath as mp

mp.mp.dps = 40
out = {}
sig = sp.Symbol('sigma')
kth = 2*(1 - sp.sqrt(sp.Rational(2, 3)))
out['kappa_theta_terminal_exact'] = str(sp.simplify(kth))
out['kappa_theta_terminal'] = float(kth)

lam_s = (1 - (1 - sig)*kth/2)**-2
out['lambda_at_0'] = float(lam_s.subs(sig, 0)); out['lambda_at_1'] = float(lam_s.subs(sig, 1))
assert sp.simplify(lam_s.subs(sig, 0) - sp.Rational(3, 2)) == 0
assert sp.simplify(lam_s.subs(sig, 1) - 1) == 0

lam_bar = sp.simplify(sp.integrate(lam_s, (sig, 0, 1)))
out['lam_bar_exact'] = str(lam_bar); out['lam_bar'] = float(lam_bar)
assert sp.simplify(lam_bar - sp.sqrt(sp.Rational(3, 2))) == 0

lam_inv = sp.simplify(sp.integrate(1/lam_s, (sig, 0, 1)))
out['lam_inv_exact'] = str(sp.nsimplify(sp.simplify(lam_inv)))
out['lam_inv'] = float(lam_inv)

dlog = sp.simplify(sp.diff(sp.log(lam_s), sig))          # = -kth/(1-(1-sigma)kth/2)
out['dlog_lambda_dsigma'] = str(dlog)
sup_dlog = sp.simplify(sp.Abs(dlog.subs(sig, 0)))        # monotone in sigma; max at sigma=0
kappa_s = sp.simplify(sup_dlog)
out['kappa_s_exact'] = str(sp.simplify(sp.nsimplify(kappa_s)))
out['kappa_s'] = float(kappa_s)
assert sp.simplify(kappa_s - 2*(sp.sqrt(sp.Rational(3, 2)) - 1)) == 0
# verify it is the sup over [0,1]
grid = [float(sp.Abs(dlog).subs(sig, sp.Rational(k, 2000))) for k in range(2001)]
out['kappa_s_grid_max'] = max(grid)
assert abs(max(grid) - float(kappa_s)) < 1e-13

muJL = sp.simplify(2*kappa_s)
out['muJ_times_L_exact'] = str(sp.nsimplify(muJL)); out['muJ_times_L'] = float(muJL)
c_window = 2*sp.log(sp.Rational(3, 2))
out['c_2log32'] = float(c_window)
out['muJL_over_c'] = float(muJL/c_window)

# the two Lemma-T' shape constants and the C(mu_*) table (recomputed here, not copied)
out['C0_absolute'] = float(sp.Rational(15, 8)*sp.pi)          # 15 pi/8
out['C0_relative'] = float(sp.Rational(15, 4)*sp.pi)          # 15 pi/4
mus = sp.Symbol('m')
Cfun = sp.pi*(3*(1 + mus)/(2*(1 - mus)**5) + sp.Rational(3, 8))
out['C_mu_table'] = {str(m): float(Cfun.subs(mus, sp.Rational(str(m))))
                     for m in ('0', '0.01', '0.05', '0.10', '0.20', '0.25', '0.50')}
out['C_TLambda_limit_exact'] = str(sp.nsimplify(sp.Rational(15, 8)*sp.pi*lam_bar))
out['C_TLambda_limit'] = float(sp.Rational(15, 8)*sp.pi*lam_bar)

C_rel = sp.simplify(sp.Rational(3, 2)*sp.pi*kappa_s)
out['C_rel_exact'] = str(sp.nsimplify(C_rel)); out['C_rel'] = float(C_rel)

# --------------------------------------------------- shear offset of the reference value
def Q(l):
    l = mp.mpf(l)
    f = lambda p: mp.cos(p)*mp.sin(p)**2*(mp.sin(p)**2 - 2*mp.cos(p)**2) \
                  / (l**2*mp.sin(p)**2 + l**-4*mp.cos(p)**2)**mp.mpf(2.5)
    return mp.quad(f, [0, mp.pi/2])
lam_f = sp.lambdify(sig, lam_s, 'mpmath')
dlog_f = sp.lambdify(sig, dlog, 'mpmath')
shear_int = mp.quad(lambda s: dlog_f(s)*Q(lam_f(s)), [0, 1])
offset = mp.mpf(3)/2*shear_int                       # = (a_true - a_ref)/M , L-independent
out['shear_offset_over_M'] = float(offset)
out['shear_offset_rel_times_L'] = float(offset/(mp.mpf(1)/2*mp.mpf(str(float(lam_bar)))))
# relative offset = offset/( (1/2) L lam_bar )  ->  (offset*2/lam_bar)/L

# --------------------------------------------------- the L table, and the refuter's r2 numbers
ref = json.load(open('~/Desktop/Solve Navier Stokes/campaign/deepest-think/'
                     'DTC-2026-09-06/gaps/refute-gap-T-lipschitz/r2_results.json'))
rows = []
for r in ref['rows']:
    L = r['L']
    a_ref = 0.5*L*float(lam_bar)
    a_true_pred = a_ref + float(offset)
    bound_abs = mp.pi*float(lam_bar)*L*(3*float(muJL)/L/8.0)      # mu = 0, muJ = muJL/L
    rows.append(dict(
        L=L,
        muJ_L_derived=float(muJL), muJ_L_refuter=r['muJ_local']*L,
        muJ_L_absdiff=abs(float(muJL) - r['muJ_local']*L),
        a_ref_derived=a_ref, a_ref_refuter=r['a_ref_shell'],
        a_true_predicted=a_true_pred, a_true_refuter=r['a_true'],
        a_true_reldiff=abs(a_true_pred - r['a_true'])/r['a_true'],
        bound_over_a_ref_times_L=float(bound_abs)/a_ref*L,
        bound_over_a_true_derived=float(bound_abs)/a_true_pred,
        bound_over_a_true_refuter=r['bound_local_rel'],
        bound_reldiff=abs(float(bound_abs)/a_true_pred - r['bound_local_rel'])/r['bound_local_rel'],
    ))
out['L_table'] = rows
out['max_muJ_L_absdiff'] = max(r['muJ_L_absdiff'] for r in rows)
out['max_a_true_reldiff'] = max(r['a_true_reldiff'] for r in rows)
out['max_bound_reldiff'] = max(r['bound_reldiff'] for r in rows)

json.dump(out, open('p3_results.json', 'w'), indent=1)
print("kappa*theta (terminal) =", out['kappa_theta_terminal_exact'], "=", out['kappa_theta_terminal'])
print("lam_bar = INT lambda dsigma =", out['lam_bar_exact'], "=", repr(out['lam_bar']))
print("lam_inv = INT dsigma/lambda =", out['lam_inv_exact'], "=", repr(out['lam_inv']))
print("kappa_s = sup L|rho lam'/lam| =", out['kappa_s_exact'], "=", repr(out['kappa_s']))
print("muJ*L   = 2 kappa_s          =", out['muJ_times_L_exact'], "=", repr(out['muJ_times_L']))
print("        = ", repr(out['muJL_over_c']), "* c   with c = 2 log(3/2) =", repr(out['c_2log32']))
print("C_rel   = (3pi/2) kappa_s    =", out['C_rel_exact'], "=", repr(out['C_rel']))
print("C0 = 15pi/8 =", repr(out['C0_absolute']), "  15pi/4 =", repr(out['C0_relative']),
      "  (15pi/8) lam_bar =", out['C_TLambda_limit_exact'], "=", repr(out['C_TLambda_limit']))
print("C(mu_*) table:", {k: round(v, 6) for k, v in out['C_mu_table'].items()})
print("shear offset (a_true - a_ref)/M =", repr(out['shear_offset_over_M']),
      "  -> relative =", repr(out['shear_offset_rel_times_L']), "/L")
print()
print(f"{'L':>9s} {'muJ*L me':>10s} {'muJ*L ref':>10s} {'a_true me':>12s} {'a_true ref':>12s}"
      f" {'bnd/a_true me':>14s} {'bnd/a_true ref':>15s} {'C_rel meas':>11s}")
for r in rows:
    print(f"{r['L']:9.3f} {r['muJ_L_derived']:10.7f} {r['muJ_L_refuter']:10.7f} "
          f"{r['a_true_predicted']:12.7f} {r['a_true_refuter']:12.7f} "
          f"{r['bound_over_a_true_derived']:14.7f} {r['bound_over_a_true_refuter']:15.7f} "
          f"{r['bound_over_a_true_refuter']*r['L']:11.5f}")
print("\nmax |muJ*L me - ref| =", out['max_muJ_L_absdiff'],
      "  max rel a_true =", out['max_a_true_reldiff'],
      "  max rel bound =", out['max_bound_reldiff'])

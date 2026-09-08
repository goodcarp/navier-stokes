"""p5 - the registered controls.  A gate that cannot fail certifies nothing (FL-043).

K1  The mu-only bound  |Delta| <= C mu M A  must FAIL on Phi = Lambda, where mu = 0 exactly and
    Delta != 0.  If it does not fire, the muJ term of LEMMA T' is decoration.
K2  Against the best SINGLE T_lambdabar the hypothesis constant mu must be Theta(1), L-independent,
    and the resulting single-lambda bound must be vacuous (>= 60x the value).  Own instrument;
    this re-runs the refuter's central claim without importing its numbers.
K3  At FIXED mu, |Delta| must be unbounded as the Jacobian is allowed to blow up in a shrinking
    axis layer -- so (H3) cannot be dropped from LEMMA T' either.  Registered growth: phi_c^-1.
"""
import json, numpy as np, lib5 as G

out = {}
lam_bar = float(np.sqrt(1.5))
Ls = [8.317766166719343, 20.0, 50.0, 100.0, 400.0]

# ------------------------------------------------------------------ K1
k1 = []
for L in Ls:
    a_phi = G.a_of_map(L, [('lambda',)])
    a_ref = 0.5*L*lam_bar
    mu, muJ, _ = G.hyp_constants(L, [('lambda',)])
    d = abs(a_phi - a_ref)
    mu_only = np.pi*L*lam_bar*(3*mu*(1 + muJ)/(2*(1 - mu)**5))     # muJ term DELETED
    k1.append(dict(L=L, mu=mu, delta=d, mu_only_bound=mu_only, violated=bool(mu_only < d)))
out['K1_rows'] = k1
out['K1_fired'] = all(r['violated'] for r in k1)

# ------------------------------------------------------------------ K2
def mu_vs_single(L, lb, n_u=801, n_phi=601):
    u = np.linspace(0, L, n_u); ph = np.linspace(1e-9, np.pi - 1e-9, n_phi)
    U, PH = np.meshgrid(u, ph, indexing='ij')
    uL, pL, JL, lam, D = G.apply_map(U, PH, L, [('lambda',)])
    s, c = np.sin(PH), np.cos(PH)
    rb, zb = lb*s, lb**-2*c
    ub = U + np.log(np.hypot(rb, zb)); pb = np.arctan2(rb, zb)
    q = np.exp(uL - ub)
    mu = float(np.max(np.sqrt(np.maximum(1 + q**2 - 2*q*np.cos(pL - pb), 0.0))))
    muJ = float(np.max(np.abs(JL/lb**2 - 1)))
    return mu, muJ

k2 = []
for L in Ls:
    grid = np.linspace(1.0, 1.5, 251)
    mus = [mu_vs_single(L, lb, 201, 201)[0] for lb in grid]
    lb_star = float(grid[int(np.argmin(mus))])
    mu_s, muJ_s = mu_vs_single(L, lb_star)
    a_true = G.a_of_map(L, [('lambda',)])
    b = G.lemmaT_prime_bound(L, mu_s, muJ_s, lb_star)     # single-lambda Lemma T, A = lb_star * L
    k2.append(dict(L=L, lambar_star=lb_star, mu_single=mu_s, muJ_single=muJ_s,
                   a_true=a_true, bound_single_rel=b/a_true))
out['K2_rows'] = k2
out['K2_mu_spread'] = max(r['mu_single'] for r in k2) - min(r['mu_single'] for r in k2)
out['K2_fired'] = bool(all(r['bound_single_rel'] > 60 for r in k2)
                       and out['K2_mu_spread'] < 1e-9)

# ------------------------------------------------------------------ K3
def m_ramp(t):
    t = np.clip(t, 0.0, 1.0); return t**2*(3 - 2*t)
def dm_ramp(t):
    return np.where((t > 0) & (t < 1), 6*t*(1 - t), 0.0)

def a_with_axis_ramp(L, mu_b, phi_c, n_phi=200, n_tau=120):
    pts = sorted({0.0, phi_c, 2*phi_c, 10*phi_c, np.pi/2,
                  np.pi - 10*phi_c, np.pi - 2*phi_c, np.pi - phi_c, np.pi})
    pts = [p for p in pts if 0.0 <= p <= np.pi]
    ph, wph = G.panels(np.array(pts), n_phi)
    ta, wta = G.gl(0.0, 1.0, n_tau)
    PH, TA = np.meshgrid(ph, ta, indexing='ij'); WP, WT = np.meshgrid(wph, wta, indexing='ij')
    U = TA*L
    uL, pL, JL, lam, D = G.apply_map(U, PH, L, [('lambda',)])
    t1, t2 = pL/phi_c, (np.pi - pL)/phi_c
    beta = mu_b*m_ramp(t1)*m_ramp(t2)*np.cos(pL)
    dbeta = mu_b*(dm_ramp(t1)/phi_c*m_ramp(t2)*np.cos(pL)
                  - m_ramp(t1)*dm_ramp(t2)/phi_c*np.cos(pL)
                  - m_ramp(t1)*m_ramp(t2)*np.sin(pL))
    pf = pL + beta
    Jf = JL*(np.sin(pf)/np.sin(pL))**3*(1 + dbeta)
    integ = np.sign(np.cos(PH))*np.cos(pf)*np.exp(-4*(uL - U))*Jf*np.sin(PH)**2
    a = 0.75*L*float(np.sum(integ*WP*WT))
    mu = float(np.max(2*np.abs(np.sin(beta/2))))
    muJ = float(np.max(np.abs(Jf/lam**2 - 1)))
    return a, mu, muJ

L3 = 8.317766166719343
a_ref3 = 0.5*L3*lam_bar
k3 = []
for mu_b in (0.3, 0.9):
    for phi_c in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
        a, mu, muJ = a_with_axis_ramp(L3, mu_b, phi_c)
        lin = np.pi*L3*lam_bar*3*mu*(1 + mu)/(2*(1 - mu)**5)
        k3.append(dict(mu_b=mu_b, phi_c=phi_c, mu=mu, muJ=muJ, a=a, delta=abs(a - a_ref3),
                       linear_bound=lin, violates_linear=bool(abs(a - a_ref3) > lin),
                       full_bound=G.lemmaT_prime_bound(L3, mu, muJ, lam_bar)))
out['K3_rows'] = k3
sub_all = [r for r in k3 if r['mu_b'] == 0.3]
sub = [r for r in sub_all if r['phi_c'] <= 1e-3]           # asymptotic subrange
x = np.log10([r['phi_c'] for r in sub]); y = np.log10([r['delta'] for r in sub])
p, c0 = np.polyfit(x, y, 1)
out['K3_exponent'] = float(-p)
out['K3_prefactor'] = float(10**c0)
xa = np.log10([r['phi_c'] for r in sub_all]); ya = np.log10([r['delta'] for r in sub_all])
out['K3_exponent_full_range'] = float(-np.polyfit(xa, ya, 1)[0])
out['K3_mu_spread'] = float(max(r['mu'] for r in sub) - min(r['mu'] for r in sub))
out['K3_full_bound_holds'] = bool(all(r['full_bound'] >= r['delta'] for r in k3))
out['K3_n_linear_violations'] = sum(1 for r in k3 if r['violates_linear'])
out['K3_fired'] = bool(out['K3_n_linear_violations'] > 0 and out['K3_mu_spread'] < 0.05
                       and out['K3_full_bound_holds'])

json.dump(out, open('p5_results.json', 'w'), indent=1)
print("K1  the mu-only bound on Phi = Lambda   (mu = 0, Delta != 0)")
for r in k1:
    print(f"   L={r['L']:8.3f}  mu={r['mu']:.1e}  |Delta|={r['delta']:.6f}  "
          f"mu-only bound={r['mu_only_bound']:.6f}  VIOLATED={r['violated']}")
print("   K1 FIRED:", out['K1_fired'])
print("\nK2  best single T_lambdabar (own instrument)")
print(f"   {'L':>9s} {'lambar*':>8s} {'mu':>10s} {'muJ':>9s} {'LemT bound/a_true':>18s}")
for r in k2:
    print(f"   {r['L']:9.3f} {r['lambar_star']:8.4f} {r['mu_single']:10.6f} "
          f"{r['muJ_single']:9.4f} {r['bound_single_rel']:18.3f}")
print("   mu spread over L =", out['K2_mu_spread'], "  K2 FIRED:", out['K2_fired'])
print("\nK3  axis ramp at fixed mu, L =", L3)
for r in k3:
    print(f"   mu_b={r['mu_b']:.1f} phi_c={r['phi_c']:.0e}  mu={r['mu']:.5f}  muJ={r['muJ']:.3e}  "
          f"|Delta|={r['delta']:14.4f}  mu-only bound={r['linear_bound']:12.4f}  "
          f"violates={r['violates_linear']}  full bound={r['full_bound']:.4e}")
print(f"   fitted |Delta| ~ {out['K3_prefactor']:.4g} phi_c^-{out['K3_exponent']:.4f}"
      f"   mu spread={out['K3_mu_spread']:.2e}   K3 FIRED: {out['K3_fired']}"
      f"   full bound holds: {out['K3_full_bound_holds']}")

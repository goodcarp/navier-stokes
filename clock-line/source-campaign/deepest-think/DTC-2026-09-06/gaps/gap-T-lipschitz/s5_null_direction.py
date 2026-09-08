"""s5 - ADDITIVE FIND: the exact null direction of a[.](0) among radial rescalings.

For Psi(u) = s(|u|) u composed after T_lambda, the integrand factor is EXACT:

    (Z_Phi/|Phi|^5) J_Phi = z_u lambda^2 (s + varrho s'(varrho)) / varrho^5 ,   varrho=|u|,

the s^4 cancelling identically.  Writing h(varrho) := varrho (s(varrho)-1), s + varrho s' - 1 = h'(varrho),

    a_Phi(0) - a_{T_lambda}(0)
      = -(3/4) int_0^pi (cos phi sin^2 phi / g(phi)^5) omega(phi)
                 * [ s(R g) - s(rho0 g) + int_{rho0 g}^{R g} (s-1) d varrho / varrho ] dphi ,
      g(phi) = sqrt(lambda^2 sin^2 phi + lambda^-4 cos^2 phi).

=> if s is log-periodic with period L/k (k a positive integer) and has zero log-mean over a period,
   then the bracket VANISHES IDENTICALLY and a is EXACTLY unchanged, for every lambda, every datum,
   and every amplitude mu (however large).  Radial rearrangement is a null direction of the axis strain.
"""
import json, numpy as np, sympy as sp
from lib_quad import *

M, rho0, R = 1.0, 1.0, 4096.0
L = np.log(R/rho0)

# symbolic: the s^4 cancellation and h' identity
vr, s_, lam_ = sp.symbols('varrho s lambda', positive=True)
sfun = sp.Function('s')
expr = (sfun(vr)*sp.Symbol('z_u'))/(sfun(vr)*vr)**5 * lam_**2*sfun(vr)**4*(sfun(vr) + vr*sp.diff(sfun(vr), vr))
expr = sp.simplify(expr)
target = sp.Symbol('z_u')*lam_**2*(sfun(vr) + vr*sp.diff(sfun(vr), vr))/vr**5
res1 = sp.simplify(expr - target)
h = vr*(sfun(vr)-1)
res2 = sp.simplify(sp.diff(h, vr) - (sfun(vr) + vr*sp.diff(sfun(vr), vr) - 1))

out = dict(cancellation_residual=str(res1), hprime_residual=str(res2))
assert res1 == 0 and res2 == 0

# numeric: k = 2 (one full log-period over the shell) must give Delta == 0 for ANY mu
P = np.array([0.0, np.pi/2, np.pi])
om = lambda p: omega_bangbang(p, M)
rows = []
for lam in (1.0, 1.5):
    aT,_,_ = a_of_map(MapT(lam), om, lam, rho0, R, P, 120, 120, True)
    for mu in (0.05, 0.2, 0.5, 0.9, 3.0):
        for k in (1, 2, 4):
            aP,_,muJ = a_of_map(MapA(lam, mu, rho0, L, k=k), om, lam, rho0, R, P, 160, 160, False)
            rows.append(dict(lam=lam, mu=mu, k=k, a_T=aT, a_Phi=aP, Delta=aP-aT,
                             rel=abs(aP-aT)/abs(aT)))
out['rows'] = rows
out['null_max_rel_k_even'] = max(r['rel'] for r in rows if r['k'] % 2 == 0)
json.dump(out, open('s5_results.json','w'), indent=1)
print("symbolic residuals:", out['cancellation_residual'], out['hprime_residual'])
print(f"{'lam':>5s} {'mu':>5s} {'k':>3s} {'a_T':>11s} {'a_Phi':>11s} {'Delta':>13s} {'rel':>10s}")
for r in rows:
    print(f"{r['lam']:5.2f} {r['mu']:5.2f} {r['k']:3d} {r['a_T']:11.7f} {r['a_Phi']:11.7f} {r['Delta']:13.3e} {r['rel']:10.2e}")
print(f"\nmax relative Delta over all even-k (full-period) rows = {out['null_max_rel_k_even']:.3e}")

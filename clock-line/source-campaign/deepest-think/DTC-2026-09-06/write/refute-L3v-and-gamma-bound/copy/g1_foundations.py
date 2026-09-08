"""g1 — foundations for L3v and the Gamma bound.  Exact algebra (sympy) + a numerical
falsification pass on the two classical inequalities the proof uses.

Everything printed here is computed; nothing is typed from memory.
"""
import json
import numpy as np
import sympy as sp

out = {}
t = sp.symbols('t')
rho, z, r, lam, M = sp.symbols('rho z r lambda M', positive=True)

# ---------------------------------------------------------------- 1. C_l^{3/2} = P_{l+1}'
res = []
for l in range(0, 15):
    C = sp.gegenbauer(l, sp.Rational(3, 2), t)
    P = sp.legendre(l + 1, t)
    res.append(sp.simplify(sp.expand(C - sp.diff(P, t))))
out['C_l_32_equals_dP_lp1'] = [str(x) for x in res]
assert all(x == 0 for x in res)

# ---------------------------------------------------------------- 2. norms
Nl = []
sup = []
for l in range(0, 15):
    C = sp.gegenbauer(l, sp.Rational(3, 2), t)
    n = sp.integrate(C**2 * (1 - t**2), (t, -1, 1))
    Nl.append(sp.nsimplify(sp.simplify(n - sp.Rational((l+1)*(l+2), 1)/(l + sp.Rational(3, 2)))))
    sup.append(sp.simplify(C.subs(t, 1) - sp.Rational((l+1)*(l+2), 2)))
out['N_l_residual'] = [str(x) for x in Nl]
out['sup_C_l_residual'] = [str(x) for x in sup]
assert all(x == 0 for x in Nl) and all(x == 0 for x in sup)

# ||C_l||_inf = C_l(1): verify the max of |C_l| on [-1,1] is at t=1 (l<=14)
maxcheck = []
tt = np.linspace(-1, 1, 20001)
for l in range(1, 15):
    f = sp.lambdify(t, sp.gegenbauer(l, sp.Rational(3, 2), t), 'numpy')
    v = np.max(np.abs(f(tt)))
    maxcheck.append(float(v / ((l + 1) * (l + 2) / 2.0)))
out['max_abs_C_l_over_C_l(1)'] = maxcheck

# ---------------------------------------------------------------- 3. exact H_l, bang-bang, lambda=1
# eta_0 = rho^{-1} W(t), W = -M sgn(t)/sqrt(1-t^2);  H_l = (1/N_l) int W C_l (1-t^2) dt
Hl_exact = {}
for l in [1, 3, 5, 7, 9]:
    C = sp.gegenbauer(l, sp.Rational(3, 2), t)
    I = -2 * sp.integrate(C * sp.sqrt(1 - t**2), (t, 0, 1))   # = int_{-1}^1 W C (1-t^2) dt , M=1
    Nl_ = sp.Rational((l+1)*(l+2), 1) / (l + sp.Rational(3, 2))
    Hl_exact[l] = sp.nsimplify(sp.simplify(I / Nl_))
out['H_l_bangbang_exact'] = {str(k): str(v) for k, v in Hl_exact.items()}
out['H_l_bangbang_float'] = {str(k): float(v) for k, v in Hl_exact.items()}
out['kappa_0_from_H1'] = float(-sp.Rational(3, 5) * Hl_exact[1])

# ---------------------------------------------------------------- 4. radial Green solution
# psi_l'' + (4/rho) psi_l' - l(l+3) psi_l/rho^2 = -H_l/rho  has particular solution
# H_l rho /(l(l+3)-4) for l != 1;  and for l = 1, psi = (H_1/5) rho log(R/rho) + ...
ll = sp.symbols('l', positive=True, integer=True)
psip = M * rho / (ll * (ll + 3) - 4)
lhs = sp.simplify(sp.diff(psip, rho, 2) + 4 / rho * sp.diff(psip, rho) - ll * (ll + 3) * psip / rho**2)
out['radial_particular_residual'] = str(sp.simplify(lhs + M / rho))
assert sp.simplify(lhs + M / rho) == 0

R = sp.symbols('R', positive=True)
psi1 = M * rho * sp.log(R / rho) / 5
lhs1 = sp.simplify(sp.diff(psi1, rho, 2) + 4 / rho * sp.diff(psi1, rho) - 1 * (1 + 3) * psi1 / rho**2)
out['radial_l1_residual'] = str(sp.simplify(lhs1 + M / rho))
assert sp.simplify(lhs1 + M / rho) == 0

# full shell solution coefficients (exact), eta_l(rho) = H_l/rho on (rho0,R)
rho0 = sp.symbols('rho_0', positive=True)
Hl = sp.symbols('H_l')
for l in [3, 5]:
    inner = rho**(-(l + 3)) * sp.integrate(rr**(l) * rr**4 * (Hl / rr), (rr, rho0, rho)) if False else None
rr = sp.symbols('rr', positive=True)
def psi_shell(l):
    inner = rho**(-(l + 3)) * sp.integrate(rr**(l + 3) * Hl, (rr, rho0, rho))
    outer = rho**l * sp.integrate(rr**(-l) * Hl, (rr, rho, R))
    return sp.simplify((inner + outer) / (2 * l + 3))
p3 = psi_shell(3)
out['psi_3_shell'] = str(sp.simplify(sp.expand(p3)))
bulk3 = Hl * rho / (3 * 6 - 4)
out['psi_3_bulk_residual_at_infinite_shell'] = str(sp.simplify(sp.limit(sp.limit(p3, rho0, 0), R, sp.oo) - bulk3))

p1 = sp.simplify((rho**(-4) * sp.integrate(rr**4 * Hl, (rr, rho0, rho)) + rho * sp.integrate(Hl / rr, (rr, rho, R))) / 5)
out['psi_1_shell'] = str(sp.simplify(sp.expand(p1)))

# ---------------------------------------------------------------- 5. Legendre sup bound (Szego 7.33.3 form)
# test  |P_n(t)| <= sqrt(2/(pi*(n+1/2))) * (1-t^2)^{-1/4}   and the weaker n-version
from numpy.polynomial import legendre as L
def Pn_vals(nmax, x):
    """all P_n(x), n=0..nmax, via recurrence; x is a 1-D array"""
    out_ = np.empty((nmax + 1, x.size))
    out_[0] = 1.0
    if nmax >= 1:
        out_[1] = x
    for n in range(1, nmax):
        out_[n + 1] = ((2 * n + 1) * x * out_[n] - n * out_[n - 1]) / (n + 1)
    return out_

x = np.cos(np.linspace(1e-9, np.pi - 1e-9, 200001))
nmax = 400
P = Pn_vals(nmax, x)
w = (1 - x**2)**0.25
ratio_np = []
ratio_n = []
for n in range(1, nmax + 1):
    ratio_np.append(float(np.max(np.abs(P[n]) * w) / np.sqrt(2.0 / (np.pi * (n + 0.5)))))
    ratio_n.append(float(np.max(np.abs(P[n]) * w) / np.sqrt(2.0 / (np.pi * n))))
out['bernstein_max_ratio_n_plus_half'] = max(ratio_np)
out['bernstein_max_ratio_n'] = max(ratio_n)
out['bernstein_argmax_n_plus_half'] = int(np.argmax(ratio_np) + 1)

# P_n(0)
Pn0 = {}
for n in [2, 4, 10, 100, 1000]:
    Pn0[str(n)] = float(Pn_vals(n, np.array([0.0]))[n, 0])
out['P_n_at_0'] = Pn0
out['P_n0_times_sqrt_pi_n_over_2'] = {k: float(v * np.sqrt(np.pi * int(k) / 2)) for k, v in Pn0.items()}

print(json.dumps(out, indent=1)[:4000])
with open('g1_results.json', 'w') as f:
    json.dump(out, f, indent=1)

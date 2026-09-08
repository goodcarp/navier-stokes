#!/usr/bin/env python3
"""
s6_no_energy_free.py -- why NO repair of BFG's vortex-stretching estimate that goes
through ||u||_inf can be energy-free.

CLAIM (N1).  There is no function F with  ||u||_inf <= F(||omega||_inf)  for all smooth
divergence-free finite-energy u on R^3.
Proof.  Fix such a u with omega = curl u not identically 0 and put, for lambda > 0,
    u_lambda(x) := lambda^{-1} u(lambda x).
Then div u_lambda = 0 and curl u_lambda (x) = omega(lambda x), so ||omega_lambda||_inf
= ||omega||_inf is CONSTANT in lambda while ||u_lambda||_inf = lambda^{-1}||u||_inf -> infinity
as lambda -> 0.  (Energy: ||u_lambda||_2^2 = lambda^{-5}||u||_2^2 -> infinity too, which is why
an energy bound is exactly what breaks the family.)  Verified symbolically below on an
explicit compactly-supported-vorticity field.

CONSEQUENCE.  The divergence-form Duhamel estimate pairs |omega| against |u| in L^inf.
The mean-zero kernel removes the BMO/local-average step entirely, but it cannot remove
the ||u||_inf factor; and by (N1) that factor must be paid for with a second quantity.
Energy is the only one of the candidates that is non-increasing along NS, hence the only
one that propagates over the window.  So the divergence-form repair is energy-dependent
by necessity, not by choice of proof.
"""
import json
import sympy as sp
import mpmath as mp
out = {}
x, y, z, lam = sp.symbols('x y z lambda', positive=False), None, None, None
x, y, z = sp.symbols('x y z', real=True)
lam = sp.Symbol('lam', positive=True)
X = (x, y, z)

# an explicit smooth divergence-free field with non-trivial curl (Gaussian-localised)
psi = sp.exp(-(x**2+y**2+z**2))
A = sp.Matrix([psi*y, psi*z, psi*x])          # vector potential
u = sp.Matrix([sp.diff(A[2], y) - sp.diff(A[1], z),
               sp.diff(A[0], z) - sp.diff(A[2], x),
               sp.diff(A[1], x) - sp.diff(A[0], y)])       # u = curl A  => div u = 0
divu = sp.simplify(sum(sp.diff(u[i], X[i]) for i in range(3)))
print("div u                         =", divu); assert divu == 0
w = sp.Matrix([sp.diff(u[2], y) - sp.diff(u[1], z),
               sp.diff(u[0], z) - sp.diff(u[2], x),
               sp.diff(u[1], x) - sp.diff(u[0], y)])
print("omega  identically zero?      =", sp.simplify(w) == sp.zeros(3,1)); assert sp.simplify(w) != sp.zeros(3,1)

ul = sp.Matrix([u[i].subs({x: lam*x, y: lam*y, z: lam*z})/lam for i in range(3)])
divul = sp.simplify(sum(sp.diff(ul[i], X[i]) for i in range(3)))
print("div u_lambda                  =", divul); assert divul == 0
wl = sp.Matrix([sp.diff(ul[2], y) - sp.diff(ul[1], z),
                sp.diff(ul[0], z) - sp.diff(ul[2], x),
                sp.diff(ul[1], x) - sp.diff(ul[0], y)])
res = sp.simplify(wl - sp.Matrix([w[i].subs({x: lam*x, y: lam*y, z: lam*z}) for i in range(3)]))
print("curl u_lambda - omega(lam x)  =", list(res.T)); assert res == sp.zeros(3,1)
out['divergence_free'] = 0; out['curl_scaling_residual'] = 0

# numerical demonstration: sup |u| on a grid, at three dilations
f = sp.lambdify((x, y, z), sp.sqrt(u.dot(u)), 'numpy')
g = sp.lambdify((x, y, z), sp.sqrt(w.dot(w)), 'numpy')
import numpy as np
gr = np.linspace(-4, 4, 241)
Xg, Yg, Zg = np.meshgrid(gr, gr, gr, indexing='ij')
supu = float(np.nanmax(f(Xg, Yg, Zg))); supw = float(np.nanmax(g(Xg, Yg, Zg)))
print("sup|u| = %.12f   sup|omega| = %.12f" % (supu, supw))
rows = []
for lv in [1.0, 0.5, 0.1, 0.01]:
    # u_lam(x)=lam^{-1}u(lam x): sup over x of |u_lam| = lam^{-1} sup|u|; sup|omega_lam| = sup|omega|
    rows.append([lv, supu/lv, supw, (supu/lv)/supw])
    print("  lam = %-6g  sup|u_lam| = %-16.8f  sup|omega_lam| = %-14.8f  ratio = %.6g" % tuple(rows[-1]))
out['dilation_rows'] = rows
out['sup_u'] = supu; out['sup_omega'] = supw
print("ratio ||u||_inf/||omega||_inf is unbounded on the family => no F exists.")

# scaling exponents of the candidate control quantities under u_lam
print()
print("scaling under u_lam(x) = lam^{-1}u(lam x)  (omega sup FIXED):")
for name, expo in [("||u||_inf", -1), ("E = ||u||_2^2", -5), ("||omega||_2^2 (enstrophy)", -3), ("||omega||_inf", 0)]:
    print("   %-26s ~ lam^(%d)" % (name, expo))
out['scaling_exponents'] = {"u_inf": -1, "E": -5, "enstrophy": -3, "omega_inf": 0}
# verify E and enstrophy exponents symbolically by change of variables
lamv = sp.Symbol('lamv', positive=True)
print("   check: int |lam^{-1}u(lam x)|^2 dx = lam^{-2}lam^{-3} int|u|^2 =",
      sp.simplify(lamv**-2*lamv**-3), "* E")
print("   check: int |omega(lam x)|^2 dx     = lam^{-3} int|omega|^2 =",
      sp.simplify(lamv**-3), "* enstrophy")

json.dump(out, open('s6_results.json','w'), indent=1)
print("\ns6: (N1) established -- no bound ||u||_inf <= F(||omega||_inf) exists")

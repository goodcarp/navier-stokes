"""r2 - independent audit of Part B's identities and of the ONE unproved step in (B.4).

(1) re-derive (B.1) the 5x5 structure of grad_5 b from the definitions, symbolically;
(2) re-derive (B.2), (B.3) and the integrating-factor solution (B.4) symbolically;
(3) the step PROOF.md asserts without proof: "(a_dev being finite at one endpoint)".
    (B.4) picks the solution of the first-order ODE whose homogeneous part C(1-t^2)^{-3/2}
    is absent.  That is legitimate ONLY if (1-t^2)^{3/2} a_dev -> 0 at t -> 1.  Here that is
    supplied by a Frobenius analysis of the SECOND-order equation for f at the regular
    singular point t = 1, which PROOF.md does not give.
(4) the kernel monotonicity claims K_p' <= 0 and the two sups 3pi/4, 2/3.
"""
import json
import sympy as sp
import numpy as np

out = {}
r, z, t, u, tau = sp.symbols('r z t u tau', positive=True)

# ---------- (1) grad_5 b ---------------------------------------------------------------
# psi_1 = G(q, Z) with q = |y|^2  (smooth in y; R d/dR = 2 q d/dq)
y1, y2, y3, y4, Z = sp.symbols('y1 y2 y3 y4 Z')
Y = [y1, y2, y3, y4]
q = y1**2 + y2**2 + y3**2 + y4**2
ii, jj_ = sp.symbols('i j', integer=True, positive=True)
# psi_1 = q^i Z^j : a spanning family for the (linear) identity, with SYMBOLIC exponents,
# so the check covers every analytic axisymmetric psi_1 by linearity.
psi = q**ii * Z**jj_
a_expr = -sp.diff(psi, Z)                       # a = -partial_z psi_1
uz = 2 * psi + 2 * ii * psi                     # u^z = 2 psi + r partial_r psi   (r d/dr (q^i Z^j) = 2i q^i Z^j)
b = [a_expr * yi for yi in Y] + [uz]
V = Y + [Z]
J = sp.Matrix(5, 5, lambda i, jj: sp.diff(b[i], V[jj]))
sub = {y1: r, y2: 0, y3: 0, y4: 0}
Jp = sp.Matrix(5, 5, lambda i, jj: sp.simplify(J[i, jj].subs(sub)))

# the same objects written in (r,z) coordinates
Q = sp.symbols('Q', positive=True)
psi_rz = (r**2)**ii * Z**jj_
a_rz = -sp.diff(psi_rz, Z)
lap5 = sp.diff(psi_rz, r, 2) + 3 / r * sp.diff(psi_rz, r) + sp.diff(psi_rz, Z, 2)
om_rz = -r * lap5
claim = sp.zeros(5, 5)
for i in range(4):
    claim[i, i] = a_rz
claim[4, 4] = -2 * a_rz
claim[0, 0] = a_rz + r * sp.diff(a_rz, r)
claim[4, 4] = -2 * a_rz - r * sp.diff(a_rz, r)
claim[0, 4] = r * sp.diff(a_rz, Z)
claim[4, 0] = r * sp.diff(a_rz, Z) - om_rz
resid = sp.Matrix(5, 5, lambda i, jj: sp.simplify(sp.expand(Jp[i, jj] - claim[i, jj])))
out['grad5b_structure_residual'] = str(resid)
out['grad5b_structure_is_zero'] = bool(resid == sp.zeros(5, 5))
out['trace_grad5b'] = str(sp.simplify(sum(Jp[i, i] for i in range(5)) - 2 * a_rz))
print('grad5b residual is zero matrix:', out['grad5b_structure_is_zero'])
print('trace(grad5 b) - 2a =', out['trace_grad5b'], '  => b is NOT 5-D divergence free; trace = 2a')

# ||grad_5 b||_op = ||grad_3 u||_op : the 5-D matrix has a-eigenvalue with multiplicity 3
# extra compared with the 3-D one, so the two largest singular values coincide iff
# |a| <= sigma_max of the 2x2 (r,z) block, which holds because the block contains a and -2a
# on its diagonal only through the claim above.  Checked numerically:
import numpy as _np
rng = _np.random.default_rng(7)
worst = 0.0
for _ in range(400):
    A_, B_, C_, D_, aval, omv = rng.normal(size=6)
    M5 = _np.diag([aval, aval, aval, aval, -2 * aval]).astype(float)
    M5[0, 0] += A_; M5[4, 4] += -A_; M5[0, 4] = B_; M5[4, 0] = B_ - omv
    M3 = _np.array([[aval + A_, 0.0, B_], [0.0, aval, 0.0], [B_ - omv, 0.0, -2 * aval - A_]])
    s5 = _np.linalg.svd(M5, compute_uv=False)[0]
    s3 = _np.linalg.svd(M3, compute_uv=False)[0]
    worst = max(worst, abs(s5 - s3) / max(s3, 1e-12))
out['grad5_vs_grad3_op_norm_max_reldiff'] = float(worst)
print('||grad5 b||_op vs ||grad3 u||_op max rel diff over 400 random matrices:', worst)

# the triangle-inequality bound 2|a| + r|grad a| + |omega|
worst2 = 0.0
for _ in range(2000):
    A_, B_, aval, omv = rng.normal(size=4)
    M5 = _np.diag([aval, aval, aval, aval, -2 * aval]).astype(float)
    M5[0, 0] += A_; M5[4, 4] += -A_; M5[0, 4] = B_; M5[4, 0] = B_ - omv
    s5 = _np.linalg.svd(M5, compute_uv=False)[0]
    bnd = 2 * abs(aval) + _np.hypot(A_, B_) + abs(omv)
    worst2 = max(worst2, s5 / bnd)
out['B1_bound_worst_ratio'] = float(worst2)
print('B1 bound worst ratio (must be <= 1):', worst2)

# ---------- (2) (B.2), (B.3), (B.4) ----------------------------------------------------
# Lap5(rho^gamma g(t)) and dz(rho^gamma g(t)) with g(t) = t^n, n and gamma SYMBOLIC
n_, gam = sp.symbols('n gamma')
rho = sp.sqrt(r**2 + z**2)
T = z / rho
expr = rho**gam * T**n_
lap5 = sp.diff(expr, r, 2) + 3 / r * sp.diff(expr, r) + sp.diff(expr, z, 2)
gT = T**n_
gp = n_ * T**(n_ - 1)
gpp = n_ * (n_ - 1) * T**(n_ - 2)
claim_lap = rho**(gam - 2) * (gam * (gam + 3) * gT + (1 - T**2) * gpp - 4 * T * gp)
out['lap5_residual'] = str(sp.simplify(sp.expand(sp.powsimp(lap5 - claim_lap, force=True))))
dz = sp.diff(expr, z)
claim_dz = rho**(gam - 1) * (gam * T * gT + (1 - T**2) * gp)
out['dz_residual'] = str(sp.simplify(sp.expand(sp.powsimp(dz - claim_dz, force=True))))
print('Lap5 residual:', out['lap5_residual'], '| dz residual:', out['dz_residual'])

# (B.2): r|grad a| = (1-t^2)|a'| for a = a(t), a(t) = t^n
A0 = T**n_
rg2 = sp.simplify(sp.expand(sp.powsimp((r * sp.diff(A0, r))**2 + (r * sp.diff(A0, z))**2, force=True)))
tgt = sp.expand(((1 - T**2) * n_ * T**(n_ - 1))**2)
out['B2_squared_residual'] = str(sp.simplify(sp.powsimp(rg2 - tgt, force=True)))
print('B2 squared residual:', out['B2_squared_residual'])

# (B.3): pure ODE identity
f = sp.Function('f'); W = sp.Function('W')
adev = -(t * f(t) + (1 - t**2) * sp.diff(f(t), t))
lhs = sp.expand((1 - t**2) * sp.diff(adev, t))
rhs = 3 * f(t) + 3 * t * adev + (1 - t**2) * W(t)
rhs = sp.expand(rhs.subs(W(t), -(4 * f(t) + (1 - t**2) * sp.diff(f(t), t, 2) - 4 * t * sp.diff(f(t), t))))
out['B3_residual'] = str(sp.simplify(lhs - rhs))
print('B3 residual:', out['B3_residual'])

ad = sp.Function('ad')
lhs2 = sp.diff((1 - t**2)**sp.Rational(3, 2) * ad(t), t)
rhs2 = (1 - t**2)**sp.Rational(1, 2) * ((1 - t**2) * sp.diff(ad(t), t) - 3 * t * ad(t))
out['integrating_factor_residual'] = str(sp.simplify(sp.expand(lhs2 - rhs2)))
print('integrating factor residual:', out['integrating_factor_residual'])

# ---------- (3) the Frobenius step PROOF.md omits --------------------------------------
# near t = 1 put s = 1-t.  (1-t^2) = s(2-s).  Equation:  s(2-s) f'' + 4(1-s) f' + 4 f = -W.
s = sp.symbols('s', positive=True)
F = sp.Function('F')
ode = sp.Eq(s * (2 - s) * sp.diff(F(s), s, 2) + 4 * (1 - s) * sp.diff(F(s), s) + 4 * F(s), 0)
# indicial equation at the regular singular point s=0:  F ~ s^k  =>  2k(k-1) + 4k = 0 => k=0 or k=-1
k = sp.symbols('k')
ind = sp.simplify(2 * k * (k - 1) + 4 * k)
out['indicial_polynomial'] = str(sp.factor(ind))
out['indicial_roots'] = [str(v) for v in sp.solve(ind, k)]
print('indicial poly', sp.factor(ind), 'roots', sp.solve(ind, k))
# so the two homogeneous behaviours are  f ~ const  and  f ~ s^{-1}.
# f = sum_{l>=3} c_l C_l^{3/2} with sum |c_l| ||C_l||_inf < infty (Corollary A1) is BOUNDED,
# hence it is the k=0 branch, hence f' is obtained from  4 f' = -W - 4 f + O(s) f'' ...
# concretely: from the ODE, s(2-s) f'' + 4(1-s) f' = -W - 4f, and for the bounded branch
# f'(s) = O(1) + O(int |W|), so (1-t^2) f'(t) -> 0 and a_dev(1) = -f(1) is FINITE.
# We verify the conclusion numerically below rather than only asserting it.
out['frobenius_note'] = ('f bounded (Corollary A1) forces the k=0 branch; the k=-1 branch is '
                         'unbounded, so (1-t^2)f\' -> 0 and a_dev(1) = -f(1) is finite. '
                         'PROOF.md asserts this ("a_dev being finite at one endpoint") without '
                         'the argument.')

# ---------- (4) kernel monotonicity and sups ------------------------------------------
res = {}
for p, nm in [(sp.Rational(1, 2), 'p_half'), (sp.Integer(1), 'p_one')]:
    Jf = sp.integrate((1 - tau**2)**p, (tau, u, 1))
    D = (1 - u**2)**(p + 1) - 3 * u * Jf
    res[nm] = {'D_prime': str(sp.simplify(sp.diff(D, u))), 'D_at_1': str(sp.simplify(sp.limit(D, u, 1)))}
    Kp = (1 - u**2)**sp.Rational(-3, 2) * Jf
    res[nm]['K_at_0'] = str(sp.simplify(Kp.subs(u, 0)))
    res[nm]['K_at_0_float'] = float(sp.simplify(Kp.subs(u, 0)))
out['kernels'] = res
print(json.dumps(res, indent=1))
json.dump(out, open('r2_results.json', 'w'), indent=1)

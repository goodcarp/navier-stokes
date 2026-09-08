"""g3 — the exact structure of grad_5 b for an axisymmetric no-swirl field, and the
bulk identity that converts the second angular derivative of the strain into zeroth-order data.
All residuals are sympy-exact."""
import json
import sympy as sp

out = {}
y1, y2, y3, y4, z = sp.symbols('y1 y2 y3 y4 z', real=True)
r_ = sp.sqrt(y1**2 + y2**2 + y3**2 + y4**2)
r, zz = sp.symbols('r z_', positive=True)
psi = sp.Function('psi')

# ---- 1. grad_5 b -------------------------------------------------------------------
Psi = psi(r_, z)
dz = sp.diff(Psi, z)
b = [-dz * y1, -dz * y2, -dz * y3, -dz * y4,
     2 * Psi + r_ * sp.diff(Psi, r_.args[0]) if False else 2 * Psi + (y1 * sp.diff(Psi, y1) + y2 * sp.diff(Psi, y2) + y3 * sp.diff(Psi, y3) + y4 * sp.diff(Psi, y4))]
# note: r d_r psi = sum y_i d_{y_i} psi
V = [y1, y2, y3, y4, z]
Jac = sp.Matrix(5, 5, lambda i, j: sp.diff(b[i], V[j]))
sub = {y1: r, y2: 0, y3: 0, y4: 0, z: zz}
J = sp.simplify(Jac.subs(sub))
a_expr = sp.simplify(-sp.diff(psi(r, zz), zz))
out['trace_minus_2a'] = str(sp.simplify(J.trace() - 2 * a_expr))
D = sp.diag(1, 1, 1, 1, -2) * a_expr
E = sp.simplify(J - D)
out['E_matrix'] = str(E.tolist())
# claimed entries
Q = sp.simplify(r * sp.diff(a_expr, r))          # r d_r a
P = sp.simplify(r * sp.diff(a_expr, zz))         # r d_z a
lap5 = sp.diff(psi(r, zz), r, 2) + 3 / r * sp.diff(psi(r, zz), r) + sp.diff(psi(r, zz), zz, 2)
omth = sp.simplify(-r * lap5)                    # omega^theta = r*eta = -r Delta_5 psi
claim = sp.Matrix([[Q, P], [P - omth, -Q]])
Esub = sp.Matrix([[E[0, 0], E[0, 4]], [E[4, 0], E[4, 4]]])
out['E_block_residual'] = str(sp.simplify(Esub - claim))
M0 = sp.Matrix(5, 5, lambda i, j: Esub[0,0] if (i==0 and j==0) else (Esub[0,1] if (i==0 and j==4) else (Esub[1,0] if (i==4 and j==0) else (Esub[1,1] if (i==4 and j==4) else 0))))
out['E_offblock_residual'] = str(sp.simplify(E - M0))

# ---- 2. ||grad_5 b||_op = ||grad_3 u||_op ------------------------------------------
ur = r * a_expr
uz = 2 * psi(r, zz) + r * sp.diff(psi(r, zz), r)
G3 = sp.Matrix([[sp.diff(ur, r), 0, sp.diff(ur, zz)],
                [0, a_expr, 0],
                [sp.diff(uz, r), 0, sp.diff(uz, zz)]])
import numpy as np
test_psi = sp.sin(2 * r) * sp.exp(-zz) + r**3 * zz**2 / 7 + sp.log(1 + r) * sp.cos(zz)
Jn = sp.Matrix(5, 5, lambda i, j: sp.simplify(J[i, j].subs({psi(r, zz): test_psi}).doit()))
G3n = sp.Matrix(3, 3, lambda i, j: sp.simplify(G3[i, j].subs({psi(r, zz): test_psi}).doit()))
Jf = sp.lambdify((r, zz), Jn, 'numpy'); Gf = sp.lambdify((r, zz), G3n, 'numpy')
rng = np.random.default_rng(11)
dd = []
for _ in range(200):
    rv, zv = rng.uniform(0.3, 3), rng.uniform(-3, 3)
    s5 = np.linalg.svd(np.array(Jf(rv, zv), float), compute_uv=False)
    s3 = np.linalg.svd(np.array(Gf(rv, zv), float), compute_uv=False)
    dd.append(abs(s5.max() - s3.max()) / s3.max())
out['op_norm_5D_vs_3D_maxreldiff'] = float(max(dd))

# ---- 3. the derivative identity  d_z[rho^l C_l^{3/2}(z/rho)] = (l+2) rho^{l-1} C_{l-1}
rho = sp.sqrt(r**2 + zz**2)
res = []
for l in range(1, 13):
    lhs = sp.diff(rho**l * sp.gegenbauer(l, sp.Rational(3, 2), zz / rho), zz)
    rhs = (l + 2) * rho**(l - 1) * sp.gegenbauer(l - 1, sp.Rational(3, 2), zz / rho)
    res.append(sp.simplify(sp.expand(sp.simplify(lhs - rhs))))
out['I1_identity_residuals'] = [str(x) for x in res]

# ---- 4. Delta_5 (rho^gamma f(t)) ---------------------------------------------------
t = sp.symbols('t', real=True)
g = sp.Function('g')
gam = sp.symbols('gamma')
expr = rho**gam * g(zz / rho)
lap = sp.simplify(sp.diff(expr, r, 2) + 3 / r * sp.diff(expr, r) + sp.diff(expr, zz, 2))
tt = zz / rho
target = rho**(gam - 2) * (gam * (gam + 3) * g(tt) + (1 - tt**2) * sp.Subs(sp.Derivative(g(t), t, 2), t, tt) - 4 * tt * sp.Subs(sp.Derivative(g(t), t), t, tt))
out['Delta5_homog_residual'] = str(sp.simplify(sp.expand(sp.simplify(lap - target))))

# ---- 5. the bulk identity (1-t^2) a'_dev = 3 f + 3 t a_dev + (1-t^2) W_dev ---------
# with psi_bulk = rho f(t),  W_dev = -[4f + (1-t^2)f'' - 4 t f'],  a_dev = -(t f + (1-t^2) f')
F = sp.Function('f')
fd = F(t)
Wdev = -(4 * fd + (1 - t**2) * sp.diff(fd, t, 2) - 4 * t * sp.diff(fd, t))
adev = -(t * fd + (1 - t**2) * sp.diff(fd, t))
lhs = (1 - t**2) * sp.diff(adev, t)
rhs = 3 * fd + 3 * t * adev + (1 - t**2) * Wdev
out['bulk_identity_residual'] = str(sp.simplify(lhs - rhs))

# ---- 6. r|grad a| = (1-t^2)|a'(t)| for a homogeneous-degree-0 strain ---------------
A = sp.Function('A')
aa = A(zz / rho)
lhs2 = r * sp.sqrt(sp.diff(aa, r)**2 + sp.diff(aa, zz)**2)
tt = zz / rho
rhs2 = (1 - tt**2) * sp.Abs(sp.Subs(sp.Derivative(A(t), t), t, tt))
out['r_grad_a_residual'] = str(sp.simplify(sp.simplify(lhs2**2 - rhs2**2)))

print(json.dumps(out, indent=1))
with open('g3_results.json', 'w') as f:
    json.dump(out, f, indent=1)

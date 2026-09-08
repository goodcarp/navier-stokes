"""
v1_reduction.py -- gap-V-aronson, DTC-2026-09-06.

Exact verification of the algebraic facts the Gaussian-localization lemma rests on.

A. div_5 b = div_3 u + 2a  for the 5D lift of an axisymmetric no-swirl field;
   hence div_5 b = 2a when u is 3D-incompressible.  (sympy, exact)
B. Delta_5 (1/r) = -1/r^3.  (sympy, exact)
C. Piola identity  sum_l d_{alpha_l} ( det(J) (J^{-1})_{l i} ) = 0.  (sympy, exact rationals)
D. THE REDUCTION.  With x = X(alpha), zeta(alpha) := eta(X(alpha)),
       det(J) * (Delta_x eta)(X(alpha))  ==  sum_l d_{alpha_l} ( det(J) G^{lk} d_{alpha_k} zeta ),
       G := (J^T J)^{-1},  J := dX/dalpha.
   => in Lagrangian labels the diffusion operator is EXACTLY divergence form with
   respect to the measure det(J) dalpha: no drift, no zeroth-order term.  The
   compressibility of the lift is entirely the weight det(J).
E. the 5D lift of the uniform axisymmetric strain u = (a r, -2 a z):
   grad_5 b = diag(a,a,a,a,-2a), operator norm 2a, div_5 b = 2a.

Everything is evaluated in EXACT rational arithmetic at random rational points
(no floating point), so a residual printed as 0 is an exact 0.
"""
import json, random, sympy as sp

OUT = {}
HERE = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/gap-V-aronson"

# ---------------- A ----------------
y = sp.symbols('y1:5', real=True); z = sp.Symbol('z', real=True)
r = sp.sqrt(sum(yi**2 for yi in y))
ur = sp.Function('u_r'); uz = sp.Function('u_z')
b = [ur(r, z)*yi/r for yi in y] + [uz(r, z)]
div5 = sum(sp.diff(b[i], v) for i, v in enumerate(list(y)+[z]))
Rr = sp.Symbol('R', positive=True)
# evaluate on the slice y = (R,0,0,0) (legitimate: after differentiation the
# expression depends on y only through r, and every field is a function of (r,z))
div5 = div5.doit().subs({y[0]: Rr, y[1]: 0, y[2]: 0, y[3]: 0})
div5 = sp.simplify(div5.subs(sp.sqrt(Rr**2), Rr))
div3 = sp.Derivative(ur(Rr, z), Rr).doit() + ur(Rr, z)/Rr + sp.Derivative(uz(Rr, z), z).doit()
residA = sp.simplify(div5 - div3 - 2*ur(Rr, z)/Rr)
print("A. div_5 b =", div5)
print("A. div_5 b - div_3 u - 2a  =", residA)
OUT['A_div5'] = str(div5); OUT['A_residual'] = str(residA)

# ---------------- B ----------------
lap5 = sum(sp.diff(1/r, v, 2) for v in list(y)+[z])
residB = sp.simplify(lap5 + 1/r**3)
rr = sp.Symbol('rr', positive=True)
residB2 = sp.simplify(sp.diff(1/rr, rr, 2) + 3/rr*sp.diff(1/rr, rr) + 1/rr**3)
print("B. Delta_5(1/r) + 1/r^3 =", residB, " ; (d_rr+3/r d_r)(1/r)+1/r^3 =", residB2)
OUT['B_residual_cartesian5'] = str(residB); OUT['B_residual_rz'] = str(residB2)

# ---------------- C, D ----------------
a1, a2, a3, t = sp.symbols('a1 a2 a3 t', real=True)
A = [a1, a2, a3]
X = [a1 + t*a2**2 + t*a3**3/5,
     a2*(1 + t*a1) + t*a3/3,
     a3*(1 + t*a1*a2/4) + t*a1**2]           # nonlinear, det J != 1
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
eta = x1**3*x2 + x2**2*x3**2 + x1*x3**3 + 2*x1**2*x3   # polynomial => exact rationals
lap_eta = sum(sp.diff(eta, v, 2) for v in (x1, x2, x3))

J = sp.Matrix(3, 3, lambda i, j: sp.diff(X[i], A[j]))
D = sp.expand(J.det())
adj = J.adjugate()                     # J^{-1} = adj/D  (polynomial entries)
# G = J^{-1} J^{-T} = adj adj^T / D^2 ;  D*G = adj adj^T / D  -> keep as (adj adj^T)/D
S = sp.expand(adj*adj.T)               # so that det(J) G = S / D

sub = {x1: X[0], x2: X[1], x3: X[2]}
zeta = eta.subs(sub, simultaneous=True)
gz = sp.Matrix([sp.diff(zeta, A[k]) for k in range(3)])
flux = (S*gz)                          # = D^2 * G * grad zeta ; we need D*G*grad zeta = flux/D
lhs = D*lap_eta.subs(sub, simultaneous=True)          # det(J) * (Delta eta)(X)
rhs = sum(sp.diff(flux[l]/D, A[l]) for l in range(3))

random.seed(20260906)
def rnd():
    return {a1: sp.Rational(random.randint(-9, 9), 7),
            a2: sp.Rational(random.randint(-9, 9), 5),
            a3: sp.Rational(random.randint(-9, 9), 3),
            t : sp.Rational(random.randint(1, 9), 11)}
piola_max = 0
red = []
for _ in range(8):
    pt = rnd()
    Jn = J.subs(pt); Dn = sp.Rational(sp.nsimplify(D.subs(pt)))
    if Dn == 0: continue
    # Piola at this point needs symbolic differentiation first
    resid = sp.nsimplify(sp.together(sp.expand((lhs - rhs).subs(pt, simultaneous=True))))
    scale = sp.Abs(sp.nsimplify(lhs.subs(pt, simultaneous=True)))
    red.append((sp.simplify(resid), scale))
piola = [sp.expand(sum(sp.diff(D*adj[l, i]/D, A[l]) for l in range(3))) for i in range(3)]
piola = [sp.simplify(p) for p in piola]
print("C. Piola residuals (exact):", piola)
OUT['C_piola'] = [str(p) for p in piola]
print("D. reduction residual (exact) at 8 random rational points:")
allzero = True
for resid, scale in red:
    print("     residual =", resid, "   |LHS| =", scale)
    if sp.simplify(resid) != 0: allzero = False
OUT['D_reduction_all_exact_zero'] = bool(allzero)
OUT['D_reduction_residuals'] = [str(v[0]) for v in red]
OUT['D_reduction_scales'] = [str(v[1]) for v in red]

# D'. the non-divergence ("drift") form, for the record:
#     (Delta eta)(X) = d_l(G^{lk} d_k zeta) + G^{lk}(d_l log D)(d_k zeta)
GG = sp.expand(adj*adj.T)              # = D^2 G
alt = sum(sp.diff(GG[l, :].dot(gz)/D**2, A[l]) for l in range(3)) \
    + sum((GG[l, :].dot(gz)/D**2)*sp.diff(sp.log(D), A[l]) for l in range(3))
pt = rnd()
residDp = sp.simplify((lap_eta.subs(sub, simultaneous=True) - alt).subs(pt, simultaneous=True))
print("D'. non-divergence form residual =", residDp)
OUT['Dprime_residual'] = str(residDp)

# ---------------- E ----------------
aa = sp.Symbol('a', real=True)
bs = [aa*yi for yi in y] + [-2*aa*z]
Gr = sp.Matrix(5, 5, lambda i, j: sp.diff(bs[i], (list(y)+[z])[j]))
print("E. grad_5 b diag =", [Gr[i, i] for i in range(5)], " trace =", sp.simplify(Gr.trace()),
      " (= 2a)  op-norm = 2|a|")
OUT['E_diag'] = [str(Gr[i, i]) for i in range(5)]
OUT['E_trace'] = str(sp.simplify(Gr.trace()))

with open(HERE + "/v1_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, default=str)
print("\nwrote v1_results.json")

# ---------------- F. ||grad_5 b||_op == ||grad_3 u||_op, exactly (independent numeric check) ----------------
import numpy as np
rng = np.random.default_rng(20260906)
def lift_grad(ur, uz, dur_r, dur_z, duz_r, duz_z, yv, zv):
    r = np.linalg.norm(yv); e = yv/r
    G5 = np.zeros((5,5))
    G5[:4,:4] = dur_r*np.outer(e,e) + (ur/r)*(np.eye(4)-np.outer(e,e))
    G5[:4, 4] = dur_z*e
    G5[4, :4] = duz_r*e
    G5[4, 4]  = duz_z
    G3 = np.array([[dur_r, 0, dur_z],[0, ur/r, 0],[duz_r, 0, duz_z]])   # (r, theta, z) frame
    return G5, G3
worst = 0.0
for _ in range(400):
    r = float(rng.uniform(0.2,3)); zv = float(rng.normal())
    yv = rng.normal(size=4); yv = yv/np.linalg.norm(yv)*r
    ur = float(rng.normal()); uz = float(rng.normal())
    dur_r, dur_z, duz_r = (float(rng.normal()) for _ in range(3))
    duz_z = -dur_r - ur/r                      # 3D incompressibility
    G5, G3 = lift_grad(ur,uz,dur_r,dur_z,duz_r,duz_z,yv,zv)
    worst = max(worst, abs(np.linalg.norm(G5,2)-np.linalg.norm(G3,2)))
    # and the divergence identity
    assert abs(np.trace(G5) - 2*ur/r) < 1e-10
print("F. max | ||grad_5 b||_op - ||grad_3 u||_op | over 400 random incompressible jets = %.3e" % worst)
OUT['F_opnorm_max_diff'] = worst
with open(HERE + "/v1_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, default=str)

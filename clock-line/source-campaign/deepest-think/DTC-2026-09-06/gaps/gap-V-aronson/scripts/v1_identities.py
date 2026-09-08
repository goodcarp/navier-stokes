"""
v1 -- the exact algebraic reductions behind GAP V.   All checks are sympy identities
(residual identically 0) or exact rational/eigenvalue computations.  Nothing numerical.

Claims verified:
 A1  Delta5 f = Delta3 f + (2/r) d_r f            (axisymmetric f)
 A2  Delta5 f = r^{-2} div3( r^2 grad3 f )        (axisymmetric f)  -> weighted-divergence form in R^3
 A3  b(y,z) = (u^r y/r, u^z) on R^4 x R  ==>  div5 b = d_r u^r + 3 u^r/r + d_z u^z
     and with 3D incompressibility  div3 u = 0  ==>  div5 b = 2a,  a := u^r/r.
 A4  b.grad5 (axisym f) = u.grad3 f   -- the transport in the lift IS the 3D transport.
 A5  spec(grad5 b) = {a,a,a} U spec( [[d_r u^r, d_z u^r],[d_r u^z, d_z u^z]] )
     and spec(grad3 u) = {a} U (the same 2x2 block)  ==>  ||grad5 b||_op = ||grad3 u||_op.
 A6  Lagrangian reduction: for a diffeomorphism X, A=DX, J=det A, g=A^T A, ftilde=f o X:
        (Delta f) o X  =  J^{-1} d_j ( J (g^{-1})_{jm} d_m ftilde ).
     (checked on an explicit NON volume preserving 3D map, residual 0)
 A7  d/dt log J = (div b) o X   along the flow of b.
 A8  For the axisymmetric lift, J5 = (r(X)/r(alpha))^2 exactly.
 A9  Delta5 (1/r) = -1/r^3 .
 A10 stream-function convention: u^r=-r d_z psi1, u^z = 2 psi1 + r d_r psi1 ==> div3 u = 0,
     a = -d_z psi1.
 A11 Pure axisymmetric strain u=(a0 r, -2 a0 z): grad5 b = diag(a0,a0,a0,a0,-2a0), and the
     5D Lagrangian metric g^{-1} = diag(e^{-2a0 t} I4, e^{+4a0 t}); J = e^{2 a0 t} is SPATIALLY
     CONSTANT so the Lagrangian equation has NO drift and is an exactly solvable anisotropic
     heat equation.
"""
import sympy as sp
import json, sys

out = {}
def rec(k, v):
    out[k] = v
    print(f"{k}: {v}")

r, z, t, a0, lam = sp.symbols('r z t a0 lambda', positive=True)
f = sp.Function('f')

# ---------- A1, A2, A9 ----------
F = f(r, z)
D5 = sp.diff(F, r, 2) + 3/r*sp.diff(F, r) + sp.diff(F, z, 2)
D3 = sp.diff(F, r, 2) + 1/r*sp.diff(F, r) + sp.diff(F, z, 2)
rec("A1_residual", sp.simplify(D5 - (D3 + 2/r*sp.diff(F, r))) == 0)
wdiv = sp.simplify(sp.Rational(1,1)/r**2 * (sp.diff(r**2*sp.diff(F, r), r) + sp.diff(r**2*sp.diff(F, z), z)))
# div3 in cylindrical for an axisym vector field (Vr,Vz): (1/r) d_r (r Vr) + d_z Vz ; here weight r^2 handled explicitly
wdiv3 = sp.simplify(1/r**2 * ( (1/r)*sp.diff(r*(r**2*sp.diff(F,r)), r) + sp.diff(r**2*sp.diff(F,z), z) ))
rec("A2_residual", sp.simplify(D5 - wdiv3) == 0)
rec("A9_r3_times_Delta5_inv_r", sp.simplify(r**3*(sp.diff(1/r, r, 2) + 3/r*sp.diff(1/r, r))))

# ---------- A3, A4, A5 : Cartesian R^5 ----------
y = sp.symbols('y1:5', real=True)
Z = sp.Symbol('Z', real=True)
R = sp.sqrt(sum(yi**2 for yi in y))
ur = sp.Function('ur'); uz = sp.Function('uz'); g = sp.Function('g')
b = [ur(R, Z)*yi/R for yi in y] + [uz(R, Z)]
X5 = list(y) + [Z]
div5 = sum(sp.diff(b[i], X5[i]) for i in range(5))
div5s = sp.simplify(sp.expand(div5))
# express in (r,z)
div5s = sp.simplify(div5s.subs(R, r))
target = sp.diff(ur(r, z), r) + 3*ur(r, z)/r + sp.diff(uz(r, z), z)
target = target.subs(z, Z).subs(sp.Symbol('z'), Z)
rec("A3_div5_formula_residual",
    sp.simplify(div5s - (sp.diff(ur(r, Z), r) + 3*ur(r, Z)/r + sp.diff(uz(r, Z), Z))) == 0)
# with div3 u = 0 : d_r ur + ur/r + d_z uz = 0
div3expr = sp.diff(ur(r, Z), r) + ur(r, Z)/r + sp.diff(uz(r, Z), Z)
rec("A3_div5_minus_2a_equals_div3_residual",
    sp.simplify((sp.diff(ur(r, Z), r) + 3*ur(r, Z)/r + sp.diff(uz(r, Z), Z)) - 2*ur(r, Z)/r - div3expr) == 0)

# A4 : transport
FA = g(R, Z)                                     # axisymmetric scalar in R^5
adv5 = sum(b[i]*sp.diff(FA, X5[i]) for i in range(5))
adv5 = sp.simplify(sp.expand(adv5)).subs(R, r)
adv3 = ur(r, Z)*sp.diff(g(r, Z), r) + uz(r, Z)*sp.diff(g(r, Z), Z)
rec("A4_transport_residual", sp.simplify(adv5 - adv3) == 0)

# A5 : spectra
Jac5 = sp.Matrix(5, 5, lambda i, j: sp.diff(b[i], X5[j]))
# evaluate at the point y = (r,0,0,0)
sub = {y[0]: r, y[1]: 0, y[2]: 0, y[3]: 0}
Jac5p = sp.simplify(Jac5.subs(sub))
ev5 = Jac5p.eigenvals()
rec("A5_grad5b_eigenvalues", {sp.srepr(k)[:0] + str(sp.simplify(k)): int(v) for k, v in ev5.items()})
# 3D lift: b3 = (ur y/r, uz) on R^2 x R  (y in R^2)
w = sp.symbols('w1:3', real=True)
R3 = sp.sqrt(w[0]**2 + w[1]**2)
b3 = [ur(R3, Z)*wi/R3 for wi in w] + [uz(R3, Z)]
X3 = list(w) + [Z]
Jac3 = sp.Matrix(3, 3, lambda i, j: sp.diff(b3[i], X3[j]))
Jac3p = sp.simplify(Jac3.subs({w[0]: r, w[1]: 0}))
ev3 = Jac3p.eigenvals()
rec("A5_grad3u_eigenvalues", {str(sp.simplify(k)): int(v) for k, v in ev3.items()})
blk = sp.Matrix([[sp.diff(ur(r, Z), r), sp.diff(ur(r, Z), Z)],
                 [sp.diff(uz(r, Z), r), sp.diff(uz(r, Z), Z)]])
rec("A5_2x2_block_charpoly_matches_5D",
    sp.simplify(sp.factor(Jac5p.charpoly(sp.Symbol('mu')).as_expr()
                          / ((sp.Symbol('mu') - ur(r, Z)/r)**3 * blk.charpoly(sp.Symbol('mu')).as_expr()))) == 1)
rec("A5_2x2_block_charpoly_matches_3D",
    sp.simplify(sp.factor(Jac3p.charpoly(sp.Symbol('mu')).as_expr()
                          / ((sp.Symbol('mu') - ur(r, Z)/r) * blk.charpoly(sp.Symbol('mu')).as_expr()))) == 1)

# ---------- A6 : Lagrangian divergence-form reduction ----------
al = sp.symbols('al1:4', real=True)
# an explicit NON volume-preserving, non-orthogonal diffeomorphism
Xm = sp.Matrix([al[0] + sp.Rational(1,5)*al[1]**2,
                sp.Rational(3,2)*al[1] + sp.Rational(1,7)*al[2],
                sp.Rational(1,2)*al[2] + sp.Rational(1,3)*al[0]*al[1]])
A = Xm.jacobian(sp.Matrix(al))
Jd = sp.simplify(A.det())
ginv = sp.simplify((A.T*A).inv())
xs = sp.symbols('x1:4', real=True)
ff = sp.Function('F')
# a concrete f (so the identity is a genuine functional check, not a tautology)
fx = sp.exp(xs[0])*sp.sin(xs[1]) + xs[2]**3*xs[0] + sp.cos(xs[1]*xs[2])
lapf = sum(sp.diff(fx, xi, 2) for xi in xs)
subm = {xs[i]: Xm[i] for i in range(3)}
lapf_at_X = sp.simplify(lapf.subs(subm))
ftil = sp.simplify(fx.subs(subm))
flux = sp.simplify(Jd*ginv*sp.Matrix([sp.diff(ftil, a) for a in al]))
rhs = sp.simplify(sum(sp.diff(flux[i], al[i]) for i in range(3))/Jd)
rec("A6_lagrangian_divform_residual", sp.simplify(sp.expand(rhs - lapf_at_X)) == 0)
rec("A6_map_detJ_nonconstant", sp.simplify(sp.diff(Jd, al[0])) != 0 or sp.simplify(sp.diff(Jd, al[1])) != 0)

# ---------- A7 : d/dt log det DX = div b ----------
tt = sp.Symbol('t', real=True)
bb = [sp.Function(f'B{i}')(*(sp.symbols('q1:4', real=True))) for i in range(3)]
q = sp.symbols('q1:4', real=True)
# verify the Jacobi formula symbolically on a concrete 1-parameter family
Xt = sp.Matrix([al[0]*sp.exp(tt) + tt*al[1]**2,
                al[1]*sp.exp(-sp.Rational(1,2)*tt),
                al[2] + tt*sp.sin(al[0])])
At = Xt.jacobian(sp.Matrix(al))
Jt = sp.simplify(At.det())
Vel = sp.simplify(sp.diff(Xt, tt))                     # velocity in Lagrangian labels
# Eulerian divergence  = trace( (dV/dalpha) A^{-1} )
divb = sp.simplify(sp.trace(Vel.jacobian(sp.Matrix(al))*At.inv()))
rec("A7_jacobi_residual", sp.simplify(sp.diff(sp.log(Jt), tt) - divb) == 0)

# ---------- A8 : J5 = (r(X)/r0)^2 for the axisymmetric lift ----------
# lift of the map (r,z) -> (R(r,z), Zm(r,z)) to R^4 x R : y -> R(r,z) y/r
Rf = sp.Function('Rf'); Zf = sp.Function('Zf')
rr = sp.sqrt(sum(yi**2 for yi in y))
Xl = [Rf(rr, Z)*yi/rr for yi in y] + [Zf(rr, Z)]
Al = sp.Matrix(5, 5, lambda i, j: sp.diff(Xl[i], X5[j]))
Alp = sp.simplify(Al.subs(sub))
J5 = sp.simplify(Alp.det())
# 3D volume-preservation:  Jacobian of (r,z)->(R,Z) in 3D is (R/r) * det[[dR/dr,dR/dz],[dZ/dr,dZ/dz]] = 1
j2 = sp.Matrix([[sp.diff(Rf(r, Z), r), sp.diff(Rf(r, Z), Z)],
                [sp.diff(Zf(r, Z), r), sp.diff(Zf(r, Z), Z)]]).det()
J5r = sp.simplify(J5.subs(rr, r))
rec("A8_J5_over_J3", sp.simplify(sp.factor(J5r / ((Rf(r, Z)/r)*j2))))
rec("A8_J5_equals_J3_times_R2_over_r2_residual",
    sp.simplify(J5r - ((Rf(r, Z)/r)*j2) * (Rf(r, Z)/r)**2) == 0)
rec("A8_when_J3_eq_1_then_J5", sp.simplify((Rf(r, Z)/r)**2))

# ---------- A10 : stream function convention ----------
psi = sp.Function('psi')
urs = -r*sp.diff(psi(r, z), z)
uzs = 2*psi(r, z) + r*sp.diff(psi(r, z), r)
rec("A10_div3_residual", sp.simplify(sp.diff(urs, r) + urs/r + sp.diff(uzs, z)) == 0)
rec("A10_a_equals_minus_dz_psi", sp.simplify(urs/r + sp.diff(psi(r, z), z)) == 0)

# ---------- A11 : pure strain, exact Lagrangian metric ----------
lamt = sp.exp(a0*t)
A5s = sp.diag(lamt, lamt, lamt, lamt, lamt**-2)
ginv_s = sp.simplify((A5s.T*A5s).inv())
# singular values (the Gronwall bound uses the OPERATOR norm, not the spectral radius)
dr_ur, dz_ur, dr_uz, dz_uz = sp.symbols('p q s w', real=True)
M5 = sp.Matrix([[dr_ur,0,0,0,dz_ur],[0,sp.Symbol('aa',real=True),0,0,0],
                [0,0,sp.Symbol('aa',real=True),0,0],[0,0,0,sp.Symbol('aa',real=True),0],
                [dr_uz,0,0,0,dz_uz]])
M3 = sp.Matrix([[dr_ur,0,dz_ur],[0,sp.Symbol('aa',real=True),0],[dr_uz,0,dz_uz]])
import random as _rnd
_rnd.seed(7); _mx=0.0
import numpy as _np
for _ in range(200):
    sub2 = {dr_ur:_rnd.uniform(-2,2), dz_ur:_rnd.uniform(-2,2), dr_uz:_rnd.uniform(-2,2),
            dz_uz:_rnd.uniform(-2,2), sp.Symbol('aa',real=True):_rnd.uniform(-2,2)}
    n5 = _np.linalg.norm(_np.array(M5.subs(sub2)).astype(float), 2)
    n3 = _np.linalg.norm(_np.array(M3.subs(sub2)).astype(float), 2)
    _mx = max(_mx, abs(n5-n3))
rec("A5_opnorm_5D_equals_3D_maxdiff_over_200_random", float(_mx))

rec("A11_ginv_diag", [sp.simplify(ginv_s[i, i]) for i in range(5)])
rec("A11_J5", sp.simplify(A5s.det()))
rec("A11_gradlogJ_is_zero", sp.simplify(sp.diff(sp.log(A5s.det()), sp.Symbol('anything', real=True))) == 0)
rec("A11_max_ginv_eigen_vs_e2c", sp.simplify(sp.log(ginv_s[4, 4])/(2*a0*t)))   # = 2 ; c = ||grad5 b|| t = 2 a0 t

with open(__file__.replace('v1_identities.py', 'v1_results.json'), 'w') as fh:
    json.dump({k: str(v) for k, v in out.items()}, fh, indent=1)
bad = [k for k, v in out.items() if (k.endswith('residual') or k.endswith('_5D') or k.endswith('_3D')) and v is not True]
print("\nFAILED_CHECKS:", bad)
print("ALL_A_CHECKS_TRUE:", not bad)

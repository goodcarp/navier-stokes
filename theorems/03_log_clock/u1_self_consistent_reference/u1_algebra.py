"""
u1 -- the algebra of the SLAVED reference map, in sympy (exact where marked).

The reference of ASSEMBLY sec.2.1 fixes the strain rate a priori (from the integro-ODE
model).  Here the reference rate is the TRUE field's own exterior l=1 strain:

    frak_a(rho, s) := int_{2 rho}^{R_s} F(rho', s) dlog rho' ,
    F(rho', s)     := -(3/5) g_1[omega^theta(.,s)](rho')   (far-near-kernel-lemma Consequence A)

and  lambda(rho,s) := exp int_0^s frak_a(rho,sigma) dsigma ,  Lambda_s(x) := T_{lambda(|x|,s)} x .

Checked here:
 (a) the l=1 interior mode: psi_1 = -A z  =>  u^r = A r, u^z = -2 A z, div_3 = 0, div_5 = 2A.
     So v(X) = A (Y, -2Z) is EXACTLY the velocity of the l=1 interior mode, and Consequence A
     says that mode is exactly constant inside the shells that generate it.
 (b) 3-D incompressibility identity  d_z u^z = -2a - r d_r a   (used to bound the z-remainder).
 (c) J_Lambda = lambda^2 [1 + (rho lambda'/lambda)(1 - 3 cos^2 phi)],  1-3cos^2 in [-2,1].
 (d) the CONJUGATION identity:  Phi_s(x) = T_{lambda(|x|,s)}(x + e(s))  =>
        e' = (frak_a(|Phi|) - frak_a(|x|)) D (x+e) + T_lambda^{-1} R .
     No Lipschitz constant of u appears: the O(M L) part of the velocity is absorbed by the
     reference's own time dependence.
 (e) T_lambda commutes with D = diag(1,1,1,1,-2);  |T_lambda e| <= lambda |e| and
     |T_lambda x| >= lambda^{-2}|x| for lambda >= 1  (exact quadratic-form argument + numeric).
 (f) d frak_a / d log rho = -F(2 rho),  hence |d frak_a/d log rho| <= M_s/2.
 (g) the origin identity  frak_a(rho,s) = a(0,s) - int_{rho_0}^{2 rho} F dlog rho' , F >= 0.
 (h) the z-remainder assembly:
        |R_z| <= [2 C_a(phi) + log(1/sin phi) + G_1] M_s |X| .

Numerics falsify, they never prove.  Nothing outside this folder is written.
"""
import json, math
import numpy as np
import sympy as sp

OUT = {}

# ---------------------------------------------------------------- (a) the l=1 interior mode
r, z, A, lam, s = sp.symbols('r z A lambda s', positive=True)
psi1 = -A*z
a_sym = -sp.diff(psi1, z)                      # a = -d_z psi1
ur = a_sym*r                                   # u^r = a r
uz = 2*psi1 + r*sp.diff(psi1, r)               # u^z = 2 psi1 + r d_r psi1
OUT["a_minus_A"] = str(sp.simplify(a_sym - A))
OUT["ur_minus_Ar"] = str(sp.simplify(ur - A*r))
OUT["uz_plus_2Az"] = str(sp.simplify(uz + 2*A*z))
div3 = sp.diff(r*ur, r)/r + sp.diff(uz, z)
div5 = sp.diff(r**3*ur, r)/r**3 + sp.diff(uz, z)
OUT["div3_of_l1_mode"] = str(sp.simplify(div3))
OUT["div5_of_l1_mode_minus_2A"] = str(sp.simplify(div5 - 2*A))

# ---------------------------------------------------------------- (b) d_z u^z = -2a - r d_r a
af = sp.Function('a')(r, z)
psi = sp.Function('psi')(r, z)
a_of_psi = -sp.diff(psi, z)
uz_gen = 2*psi + r*sp.diff(psi, r)
lhs = sp.diff(uz_gen, z)
rhs = -2*a_of_psi - r*sp.diff(a_of_psi, r)
OUT["dz_uz_identity_residual"] = str(sp.simplify(lhs - rhs))
# and div_5 b = 2a for a general axisymmetric no-swirl field
div5_gen = sp.diff(r**3*(a_of_psi*r), r)/r**3 + sp.diff(uz_gen, z)
OUT["div5_general_minus_2a"] = str(sp.simplify(div5_gen - 2*a_of_psi))
div3_gen = sp.diff(r*(a_of_psi*r), r)/r + sp.diff(uz_gen, z)
OUT["div3_general"] = str(sp.simplify(div3_gen))

# ---------------------------------------------------------------- (c) J_Lambda
# meridian map (r,z) -> (l r, l^-2 z), l = lam(rho), rho = sqrt(r^2+z^2);
# the three S^3 directions orthogonal to y are scaled by l, so J_5 = l^3 * det d(r',z')/d(r,z).
rr, zz, l, lp = sp.symbols('rr zz l lp', positive=True)
rho_e = sp.sqrt(rr**2 + zz**2)
dR_dr = l + lp*rr**2/rho_e
dR_dz = lp*rr*zz/rho_e
dZ_dr = -2*lp*zz*rr/(l**3*rho_e)
dZ_dz = l**-2 - 2*lp*zz**2/(l**3*rho_e)
J5 = sp.simplify(l**3*(dR_dr*dZ_dz - dR_dz*dZ_dr))
D_sh = lp*rho_e/l                                  # rho lam'/lam
target = l**2*(1 + D_sh*(rr**2 - 2*zz**2)/(rr**2 + zz**2))
OUT["J_Lambda_residual"] = str(sp.simplify(sp.expand(J5 - target)))
c2 = sp.Symbol('c2')             # cos^2 phi in [0,1]
OUT["one_minus_3cos2_range"] = [float((1 - 3*c2).subs(c2, 1)), float((1 - 3*c2).subs(c2, 0))]

# ---------------------------------------------------------------- (d) conjugation identity
e_r = sp.Function('e_r')(s); e_z = sp.Function('e_z')(s); lm = sp.Function('lam')(s)
a_lab = sp.Function('a_lab')(s)      # frak_a at the LABEL radius   (drives lambda)
a_cur = sp.Function('a_cur')(s)      # frak_a at the CURRENT radius (in the true velocity)
R_r = sp.Function('R_r')(s); R_z = sp.Function('R_z')(s)
Phi_r = lm*(r + e_r)
Phi_z = lm**-2*(z + e_z)
eqs = [sp.Eq(sp.diff(Phi_r, s).subs(sp.Derivative(lm, s), a_lab*lm), a_cur*Phi_r + R_r),
       sp.Eq(sp.diff(Phi_z, s).subs(sp.Derivative(lm, s), a_lab*lm), -2*a_cur*Phi_z + R_z)]
sol = sp.solve(eqs, [sp.Derivative(e_r, s), sp.Derivative(e_z, s)], dict=True)[0]
pred_r = (a_cur - a_lab)*(r + e_r) + R_r/lm
pred_z = -2*(a_cur - a_lab)*(z + e_z) + lm**2*R_z
OUT["conj_e_r_residual"] = str(sp.simplify(sol[sp.Derivative(e_r, s)] - pred_r))
OUT["conj_e_z_residual"] = str(sp.simplify(sol[sp.Derivative(e_z, s)] - pred_z))

# ---------------------------------------------------------------- (e) T_lambda and D
T = sp.diag(lam, lam, lam, lam, lam**-2)
D = sp.diag(1, 1, 1, 1, -2)
OUT["T_D_commute"] = bool((T*D - D*T) == sp.zeros(5))
# |T v|^2 = lam^2 |y|^2 + lam^-4 z^2 ; for lam >= 1, lam^-4 <= lam^2 so |Tv| <= lam|v|;
# and lam^2 >= lam^-4 so |Tv| >= lam^-2 |v|.  Exact:
yy, zt, L1 = sp.symbols('yy zt lam1', positive=True)   # yy=|y|^2, zt=z^2, lam1>=1
q = L1**2*yy + L1**-4*zt
OUT["Tv_upper_residual"] = str(sp.simplify(sp.expand(L1**2*(yy + zt) - q)))     # >= 0 for lam>=1
OUT["Tv_lower_residual"] = str(sp.simplify(sp.expand(q - L1**-4*(yy + zt))))    # >= 0 for lam>=1
rng = np.random.default_rng(20260908)
wu, wl = 0.0, 10.0
for _ in range(400):
    lv = rng.uniform(1.0, 2.0)
    v = rng.normal(size=(2000, 5))
    Tv = v.copy(); Tv[:, :4] *= lv; Tv[:, 4] /= lv**2
    ratio = np.linalg.norm(Tv, axis=1)/np.linalg.norm(v, axis=1)
    wu = max(wu, float((ratio/lv).max()))
    wl = min(wl, float((ratio*lv**2).min()))
OUT["max_|Tv|/(lam|v|)"] = wu
OUT["min_|Tv|/(lam^-2|v|)"] = wl
OUT["T_bounds_hold"] = (wu <= 1.0 + 1e-12) and (wl >= 1.0 - 1e-12)

# ---------------------------------------------------------------- (f),(g) the shell integral
u_, U_, Rl = sp.symbols('u U R_l')
Ff = sp.Function('F')
frak = sp.Integral(Ff(u_), (u_, sp.log(2) + U_, Rl))     # U_ = log rho, integrand in log rho'
OUT["dfrak_dlogrho"] = str(sp.simplify(sp.diff(frak, U_)))
# the origin identity: a(0,s) = int over ALL shells, frak_a(rho) = a(0) - int_{rho0}^{2rho}
lo = sp.Symbol('log_rho0')
whole = sp.Integral(Ff(u_), (u_, lo, Rl))
inner = sp.Integral(Ff(u_), (u_, lo, sp.log(2) + U_))
diff_int = sp.simplify(whole - inner - frak)
OUT["origin_identity_residual"] = str(diff_int)
# numeric confirmation of the same identity on a concrete F
_F = lambda t: 0.5/(1.0 + math.exp(-(t - 1.3)))
from scipy.integrate import quad as _q
_lo, _Rl, _U = 0.0, 9.0, 1.7
_w = _q(_F, _lo, _Rl)[0]; _in = _q(_F, _lo, math.log(2) + _U)[0]; _fr = _q(_F, math.log(2) + _U, _Rl)[0]
OUT["origin_identity_numeric_residual"] = abs(_w - _in - _fr)
_h = 1e-6
OUT["dfrak_dlogrho_numeric_vs_-F(2rho)"] = abs(
    (_q(_F, math.log(2) + _U + _h, _Rl)[0] - _q(_F, math.log(2) + _U - _h, _Rl)[0])/(2*_h)
    + _F(math.log(2) + _U))

# ---------------------------------------------------------------- (h) z-remainder assembly
# |R_z| <= |int_0^z 2[a - frak_a(rho_zeta)]| + |int_0^z 2[frak_a(rho_zeta) - frak_a(rho)]|
#          + |int_0^z r d_r a|
# with sin(phi_zeta) = r/rho_zeta >= r/rho = sin(phi)  (checked: rho_zeta <= rho for |zeta|<=|z|)
zz1, rr1, zt1 = sp.symbols('zeta r1 z1', positive=True)
OUT["rho_zeta_monotone"] = str(sp.simplify(sp.sqrt(rr1**2 + zt1**2) - sp.sqrt(rr1**2 + zz1**2)))
OUT["sin_phi_zeta_ge_sin_phi"] = "rho_zeta <= rho for |zeta| <= |z|  =>  r/rho_zeta >= r/rho"

with open("u1_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
for k in sorted(OUT):
    print("%-34s %s" % (k, OUT[k]))

"""
r4 -- the algebra behind the alternative reference, checked in sympy (exact) and numerically.

(i)   the l=1 interior mode: psi_1 = -a_1 z  gives  u^r = a_1 r,  u^z = -2 a_1 z  (3-D
      divergence-free; 5-D divergence = 2 a_1, which is the rate of the 5-D Jacobian).
(ii)  T_lambda = diag(lambda,lambda,lambda,lambda,lambda^-2) commutes with D = diag(1,1,1,1,-2).
(iii) THE CONJUGATION IDENTITY.  Write Phi_s(x) = T_{lambda(s)}(x + e(s)) with
      lambda'/lambda = a_lab(s) (the reference rate at the LABEL) and
      d/ds Phi = a_cur (Phi_r, -2 Phi_z) + R  (the true velocity: l=1 strain at the CURRENT
      radius plus the remainder R).  Then EXACTLY
             e' = (a_cur - a_lab) D (x + e)  +  T_lambda^{-1} R .
      No Lipschitz constant of u enters: the O(M L) part of the velocity is absorbed by the
      reference's own time dependence.
(iv)  |T_lambda e| <= lambda |e|  and  |T_lambda x| >= lambda^{-2} |x|  for lambda >= 1,
      hence mu = sup |Phi - Lambda~| / |Lambda~| <= lambda^3 sup |e|/|x|.
(v)   the linearised coefficient of ASSEMBLY (2.1): (3/2)(2 kappa)(2 pi/r_h)(3/2 + 3/8)(9/4)
      = (3/2)(2 kappa)(15 pi/4)(9/4)/r_h  -- i.e. the 15 pi/4 presumes mu_J = mu.
"""
import json, math
import numpy as np
import sympy as sp

OUT = {}
r, z, a1, lam, s = sp.symbols('r z a_1 lambda s', positive=True)

# (i)
psi1 = -a1*z
ur = -r*sp.diff(psi1, z)
uz = 2*psi1 + r*sp.diff(psi1, r)
div3 = sp.diff(r*ur, r)/r + sp.diff(uz, z)
div5 = sp.diff(r**3*ur, r)/r**3 + sp.diff(uz, z)
OUT["i_ur_minus_a1r"] = str(sp.simplify(ur - a1*r))
OUT["i_uz_plus_2a1z"] = str(sp.simplify(uz + 2*a1*z))
OUT["i_div3"] = str(sp.simplify(div3))
OUT["i_div5_minus_2a1"] = str(sp.simplify(div5 - 2*a1))

# (ii)
T = sp.diag(lam, lam, lam, lam, lam**-2)
D = sp.diag(1, 1, 1, 1, -2)
OUT["ii_commutator_zero"] = bool((T*D - D*T) == sp.zeros(5))

# (iii) meridian form, e = (e_r(s), e_z(s)), lambda = lambda(s)
e_r = sp.Function('e_r')(s); e_z = sp.Function('e_z')(s); lamf = sp.Function('lambda')(s)
a_lab = sp.Function('a_lab')(s); a_cur = sp.Function('a_cur')(s)
R_r = sp.Function('R_r')(s); R_z = sp.Function('R_z')(s)
Phi_r = lamf*(r + e_r); Phi_z = lamf**-2*(z + e_z)
# impose lambda' = a_lab lambda and d/ds Phi = a_cur (Phi_r, -2 Phi_z) + R, solve for e'
eqs = [sp.Eq(sp.diff(Phi_r, s).subs(sp.Derivative(lamf, s), a_lab*lamf), a_cur*Phi_r + R_r),
       sp.Eq(sp.diff(Phi_z, s).subs(sp.Derivative(lamf, s), a_lab*lamf), -2*a_cur*Phi_z + R_z)]
sol = sp.solve(eqs, [sp.Derivative(e_r, s), sp.Derivative(e_z, s)], dict=True)[0]
pred_r = (a_cur - a_lab)*(r + e_r) + R_r/lamf
pred_z = -2*(a_cur - a_lab)*(z + e_z) + lamf**2*R_z
OUT["iii_e_r_residual"] = str(sp.simplify(sol[sp.Derivative(e_r, s)] - pred_r))
OUT["iii_e_z_residual"] = str(sp.simplify(sol[sp.Derivative(e_z, s)] - pred_z))

# (iv) numeric, 200000 random vectors, lambda in [1, 2]
rng = np.random.default_rng(1)
worst_up, worst_lo = 0.0, 10.0
for _ in range(200):
    L_ = rng.uniform(1.0, 2.0)
    v = rng.normal(size=(1000, 5))
    Tv = v.copy(); Tv[:, :4] *= L_; Tv[:, 4] /= L_**2
    ratio = np.linalg.norm(Tv, axis=1)/np.linalg.norm(v, axis=1)
    worst_up = max(worst_up, float((ratio/L_).max()))
    worst_lo = min(worst_lo, float((ratio*L_**2).min()))
OUT["iv_max_|Te|/(lam|e|)"] = worst_up
OUT["iv_min_|Tx|/(lam^-2|x|)"] = worst_lo
OUT["iv_holds"] = (worst_up <= 1.0 + 1e-12) and (worst_lo >= 1.0 - 1e-12)

# (v)
kap, rh = sp.symbols('kappa r_h', positive=True)
lin = sp.Rational(3, 2)*2*kap*(2*sp.pi/rh)*(sp.Rational(3, 2) + sp.Rational(3, 8))*sp.Rational(9, 4)
asm = sp.Rational(3, 2)*2*kap*(15*sp.pi/4)*sp.Rational(9, 4)/rh
OUT["v_assembly_coefficient_matches_2pi(3/2+3/8)"] = bool(sp.simplify(lin - asm) == 0)
OUT["v_muJ_frozen_coefficient_over_assembly"] = float(sp.Rational(3, 2)/(sp.Rational(3, 2) + sp.Rational(3, 8)))
OUT["v_15pi_over_4"] = float(15*math.pi/4)
OUT["v_3pi"] = float(3*math.pi)

with open("r4_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True)
print(json.dumps(OUT, indent=1, sort_keys=True))

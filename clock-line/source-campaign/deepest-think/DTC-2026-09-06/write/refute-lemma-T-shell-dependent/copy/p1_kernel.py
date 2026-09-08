"""p1 - the 5D strain kernel at the origin, its gradient, the two sharp constants,
and det D(Lambda) for the SHELL-DEPENDENT strain Lambda(x) = T_{lambda(|x|)} x.

Exact (sympy).  Nothing is imported from another seat.

Convention.  far-near-kernel-lemma sec 0:  a(x) = INT K(x-x') eta(x') dx',
K(w) = 3 w_z/(8 pi^2 |w|^5).  At x = 0, writing the integration variable as x,
    a(0) = INT Kcal(x) eta(x) dx ,     Kcal(x) := K(-x) = -(3/(8 pi^2)) x_z/|x|^5 .
Kcal is the kernel used by gap-T-lipschitz.  The task brief quotes '+3 w_z/(8 pi^2 |w|^5)' as
the origin kernel; that is K, not Kcal, and with it a(0) of the campaign datum comes out
NEGATIVE.  Sign fixed and checked in section 4 below.
"""
import json, sympy as sp

out = {}
ph = sp.Symbol('phi', real=True)

# ------------------------------------------------- 1. Kcal and grad Kcal in R^5
X = sp.symbols('x1:6', real=True)                     # x5 == z
rho = sp.sqrt(sum(x**2 for x in X))
Kcal = -sp.Rational(3, 8)/sp.pi**2 * X[4]/rho**5
c = X[4]/rho

grad = [sp.simplify(sp.diff(Kcal, x)) for x in X]
claim = [sp.simplify(-sp.Rational(3, 8)/sp.pi**2 * rho**-5 * ((1 if i == 4 else 0) - 5*c*X[i]/rho))
         for i in range(5)]
res_grad = [sp.simplify(grad[i] - claim[i]) for i in range(5)]
out['grad_residuals'] = [str(r) for r in res_grad]
assert all(r == 0 for r in res_grad), res_grad

vec2 = sum(((1 if i == 4 else 0) - 5*c*X[i]/rho)**2 for i in range(5))
res_vec = sp.simplify(sp.expand(sp.simplify(vec2)) - (1 + 15*c**2))
out['unit_vector_residual'] = str(res_vec)
assert res_vec == 0

C_K, C_gradK = sp.Rational(3, 8)/sp.pi**2, sp.Rational(3, 2)/sp.pi**2
out['C_K_exact'], out['C_K'] = str(C_K), float(C_K)
out['C_gradK_exact'], out['C_gradK'] = str(C_gradK), float(C_gradK)

# ------------------------------------------------- 2. det D(Lambda), exact, two routes
# Only the pointwise pair (lambda(rho), rho lambda'(rho)/lambda(rho)) can enter the determinant,
# and the family lambda(rho) = a rho^m realises every such pair; so verifying on it is general.
a_, m_ = sp.symbols('a m', positive=True)
r_, z_ = sp.symbols('r z', positive=True)
Y = sp.symbols('y1:5', real=True)
ZZ = sp.Symbol('zz', real=True)
rh = sp.sqrt(sum(y**2 for y in Y) + ZZ**2)
lamf = a_*rh**m_
LamMap = [lamf*Y[i] for i in range(4)] + [lamf**-2*ZZ]
vars5 = list(Y) + [ZZ]
Jmat = sp.Matrix(5, 5, lambda i, j: sp.diff(LamMap[i], vars5[j]))
sub = {Y[0]: r_, Y[1]: 0, Y[2]: 0, Y[3]: 0, ZZ: z_}
Jmat_s = sp.simplify(Jmat.subs(sub))                 # block diagonal: lam*I_3 (+) 2x2 in (y1,z)
detJ_s = sp.simplify(Jmat_s.det())
rho_s = sp.sqrt(r_**2 + z_**2)
claim_det = (a_*rho_s**m_)**2 * (1 + m_*(r_**2 - 2*z_**2)/rho_s**2)   # rho lam'/lam = m
res_det = sp.simplify(sp.expand(sp.powsimp(detJ_s - claim_det, force=True)))
out['detDLambda_residual_autodiff'] = str(res_det)
assert res_det == 0, res_det

# route 2: hand-built entries with free values Lv = lambda, Lp = lambda'
Lv, Lp = sp.symbols('Lv Lp', positive=True)
P = sp.Symbol('P', positive=True)                    # rho
Mh = sp.zeros(5, 5)
yv = [r_, 0, 0, 0]
for i in range(4):
    for j in range(4):
        Mh[i, j] = Lv*(1 if i == j else 0) + Lp*yv[i]*yv[j]/P
    Mh[i, 4] = Lp*yv[i]*z_/P
for j in range(4):
    Mh[4, j] = -2*Lv**-3*Lp*z_*yv[j]/P
Mh[4, 4] = Lv**-2 - 2*Lv**-3*Lp*z_**2/P
det_h = sp.simplify(Mh.det().subs(P, sp.sqrt(r_**2+z_**2)))
claim_h = Lv**2*(1 + (sp.sqrt(r_**2+z_**2)*Lp/Lv)*(r_**2-2*z_**2)/(r_**2+z_**2))
out['detDLambda_residual_handbuilt'] = str(sp.simplify(det_h - claim_h))
assert sp.simplify(det_h - claim_h) == 0
out['detDLambda_formula'] = "lambda^2 * (1 + (rho lambda'/lambda) * (sin^2 phi - 2 cos^2 phi))"

# weight range:  sin^2 - 2 cos^2 = 1 - 3 cos^2 in [-2, 1]  =>  |.| <= 2
res_w = sp.simplify(sp.sin(ph)**2 - 2*sp.cos(ph)**2 - (1 - 3*sp.cos(ph)**2))
out['shear_weight_identity_residual'] = str(res_w)
assert res_w == 0
out['shear_weight_range'] = [-2.0, 1.0]

# ------------------------------------------------- 3. two logarithmic derivative ranges
lm, ps = sp.symbols('lam psi', positive=True)
gh2 = lm**-2*sp.sin(ps)**2 + lm**4*sp.cos(ps)**2          # |T_lam^{-1} w|^2/|w|^2
d1 = sp.simplify(lm*sp.diff(sp.log(sp.sqrt(gh2)), lm)
                 - (-lm**-2*sp.sin(ps)**2 + 2*lm**4*sp.cos(ps)**2)/gh2)
out['dlog_ghat_residual'] = str(d1); assert d1 == 0
out['dlog_ghat_range'] = [-1.0, 2.0]     # convex combination of -1 and +2

g2 = lm**2*sp.sin(ph)**2 + lm**-4*sp.cos(ph)**2           # |T_lam x|^2/|x|^2
d2 = sp.simplify(lm*sp.diff(sp.log(sp.sqrt(g2)), lm)
                 - (lm**2*sp.sin(ph)**2 - 2*lm**-4*sp.cos(ph)**2)/g2)
out['dlog_g_residual'] = str(d2); assert d2 == 0
out['dlog_g_range'] = [-2.0, 1.0]        # convex combination of +1 and -2

# ------------------------------------------------- 4. sign sanity: a(0) for the cap at lambda=1
# a(0) = INT Kcal eta_0 dx5,  eta_0 = -M sgn(z)/r,  dx5 = 2 pi^2 rho^4 sin^3(phi) drho dphi.
# radial powers: rho (from x_z) * rho^-5 * rho^-1 (from 1/r) * rho^4 = rho^-1  ->  dlog rho, factor L
M_ = sp.Symbol('M', positive=True)
ang = sp.integrate(2*sp.pi**2*sp.Rational(3, 8)/sp.pi**2*M_*sp.Abs(sp.cos(ph))*sp.sin(ph)**2,
                   (ph, 0, sp.pi))
out['a0_cap_over_L'] = str(sp.simplify(ang))
assert sp.simplify(ang - M_/2) == 0

json.dump(out, open('p1_results.json', 'w'), indent=1)
for k in ['grad_residuals', 'unit_vector_residual', 'C_K', 'C_gradK',
          'detDLambda_residual_autodiff', 'detDLambda_residual_handbuilt', 'detDLambda_formula',
          'shear_weight_identity_residual', 'shear_weight_range',
          'dlog_ghat_residual', 'dlog_ghat_range', 'dlog_g_residual', 'dlog_g_range',
          'a0_cap_over_L']:
    print(f"{k:36s} {out[k]}")

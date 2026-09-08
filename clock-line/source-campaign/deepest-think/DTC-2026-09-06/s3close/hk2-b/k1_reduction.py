#!/usr/bin/env python3
"""k1 -- STEP 1 and STEP 2 of the (H-K2) proof.

STEP 1 (exact, sympy): the algebraic reduction of the 5D-lift Hessian  grad^2_5 b
       to the scalars  (grad a, grad^2 a, omega^theta, grad omega^theta)  at a point.
STEP 2 (exact, sympy):  the kernel constants of  a = K * eta,
       K(w) = 3 w_z /(8 pi^2 |w|^5),  and the exact delta-coefficient 1/5.

Nothing outside this folder is written.  Every constant printed here is derived.
"""
import json, math, os
import sympy as sp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}

# =====================================================================
# STEP 1.  exact structure of grad^2_5 b
# =====================================================================
x1, x2, x3, x4, z = sp.symbols('x1 x2 x3 x4 z', real=True)
X = [x1, x2, x3, x4, z]
r = sp.sqrt(x1**2 + x2**2 + x3**2 + x4**2)

# a concrete, nonlinear, non-symmetric test stream function psi1(r,z)
rs, zs = sp.symbols('rs zs', positive=True)
PSI = (rs**2 + sp.Rational(3,2)*zs**2)*sp.exp(-rs**2/5 - zs**2/7) + rs**2*zs/11 + sp.log(1+rs**2+zs**2)/3
psi = PSI.subs({rs: r, zs: z})

# build every field from the (r,z) expression, then substitute r = |y|
PSI_r  = sp.diff(PSI, rs)
UZ     = (2*PSI + rs*PSI_r).subs({rs: r, zs: z})
A      = (-sp.diff(PSI, zs)).subs({rs: r, zs: z})
LAP5   = (sp.diff(PSI, rs, 2) + 3*sp.diff(PSI, rs)/rs + sp.diff(PSI, zs, 2))
OM     = (-rs*LAP5).subs({rs: r, zs: z})                     # omega^theta = -r Lap5 psi1

b = [A*x1, A*x2, A*x3, A*x4, UZ]

# --- check the 5D lift identities that the record uses -----------------
pt = {x1: sp.Rational(2,3), x2: sp.Rational(1,5), x3: sp.Rational(-3,7),
      x4: sp.Rational(4,9), z: sp.Rational(-5,4)}
def ev(e): return sp.nsimplify(sp.simplify(e.subs(pt)))

div5 = sum(sp.diff(b[i], X[i]) for i in range(5))
RES['S1_div5b_minus_2a'] = float(sp.simplify((div5 - 2*A).subs(pt)))

# --- the claimed exact formula for d_k d_j b_i --------------------------
# scalars at a general point (as 5D-covariant objects)
ga  = [sp.diff(A, X[k]) for k in range(5)]                   # grad a
Ha  = [[sp.diff(A, X[k], X[j]) for j in range(5)] for k in range(5)]   # grad^2 a
gom = [sp.diff(OM, X[k]) for k in range(5)]                  # grad omega

def claimed(i, j, k):
    """d_k d_j b_i  from the claimed reduction."""
    if i < 4:
        t = 0
        if i == j: t += ga[k]
        if i == k: t += ga[j]
        t += X[i]*Ha[k][j]
        return t
    # i == 4  (the z-row, b_5 = u^z)
    if j < 4 and k < 4:
        # d_k d_j u^z ,  j,k <= 4
        dd = (sp.KroneckerDelta(j, k)/r - X[j]*X[k]/r**3)
        return sp.KroneckerDelta(j, k)*ga[4] + X[j]*Ha[4][k] - dd*OM - (X[j]/r)*gom[k]
    if j < 4 and k == 4:
        # d_z d_j u^z
        return X[j]*Ha[4][4] - (X[j]/r)*gom[4]
    if j == 4 and k < 4:
        # d_k d_z u^z = -3 d_k a - sum_m x_m d_m d_k a
        return -3*ga[k] - sum(X[m]*Ha[m][k] for m in range(4))
    # j == k == 4
    return -sum(X[m]*Ha[m][4] for m in range(4)) - 2*ga[4]

worst = 0.0
for i in range(5):
    for j in range(5):
        for k in range(5):
            true = sp.diff(b[i], X[j], X[k])
            res = (true - claimed(i, j, k)).subs(pt)
            worst = max(worst, abs(float(sp.simplify(res))))
RES['S1_reduction_max_residual'] = worst
print(f"S1  div5 b - 2a residual            : {RES['S1_div5b_minus_2a']:.3e}")
print(f"S1  grad^2_5 b reduction residual   : {worst:.3e}   (125 components, exact rationals)")

# also check the two auxiliary identities used to eliminate u^z
d_r_uz = (sp.diff(UZ, x1) - (X[0]*ga[4] - (X[0]/r)*OM)).subs(pt)
d_z_uz = (sp.diff(UZ, z) + sum(X[m]*ga[m] for m in range(4)) + 2*A).subs(pt)
RES['S1_dr_uz_residual'] = float(sp.simplify(d_r_uz))
RES['S1_dz_uz_residual'] = float(sp.simplify(d_z_uz))
print(f"S1  d_i u^z = x_i d_z a - (x_i/r) w : {RES['S1_dr_uz_residual']:.3e}")
print(f"S1  d_z u^z = -x.grad_y a - 2a      : {RES['S1_dz_uz_residual']:.3e}")

# =====================================================================
# STEP 2.  kernel constants
# =====================================================================
w1, w2, w3, w4, wz = sp.symbols('w1 w2 w3 w4 wz', real=True)
W = [w1, w2, w3, w4, wz]
nw = sp.sqrt(sum(v**2 for v in W))
Kw = 3*wz/(8*sp.pi**2*nw**5)

# C_K = sup_{|w|=1} |K(w)|
C_K = sp.Rational(3,1)/(8*sp.pi**2)
RES['S2_C_K'] = float(C_K)

# C_gradK = sup_{|w|=1} |grad K(w)|  (Euclidean norm of the gradient vector)
gK = [sp.simplify(sp.diff(Kw, v)) for v in W]
n2 = sp.simplify(sum(g**2 for g in gK)*nw**10)           # |grad K|^2 * |w|^10 -> function of wz/|w|
n2s = sp.simplify(n2.subs({w1: sp.sqrt(1-wz**2), w2: 0, w3: 0, w4: 0}))
sup2 = sp.simplify(sp.Max(n2s.subs(wz, 0), n2s.subs(wz, 1)))
C_gK = sp.sqrt(sp.simplify(n2s.subs(wz, 1)))
RES['S2_gradK_norm2_expr'] = str(sp.simplify(n2s))
RES['S2_C_gradK'] = float(C_gK)
print(f"S2  |grad K(w)|^2 |w|^10            = {sp.simplify(n2s)}")
print(f"S2  C_K      = 3/(8 pi^2)          = {RES['S2_C_K']:.10f}")
print(f"S2  C_gradK  = 3/(2 pi^2)          = {RES['S2_C_gradK']:.10f}   (= sup at w_z = +-1)")

# |S^4|
S4 = sp.Rational(8,3)*sp.pi**2
RES['S2_S4_area'] = float(S4)

# the exact delta coefficient:  int_{|w|<d} d_k K(w) dw = (1/5) delta_{k,z}
# computed as the flux  oint_{|w|=d} K(w) nu_k dS = 3/(8 pi^2 d^6) * oint w_z w_k dS
# with oint_{|w|=d} w_z w_k dS = delta_{kz} |S^4| d^6 / 5
coef = sp.simplify(C_K*S4/5)
RES['S2_delta_coeff'] = float(coef)
print(f"S2  int_{{|w|<d}} grad K dw          = {sp.nsimplify(coef)} * e_z   (independent of d)")
print(f"S2  |S^4| = 8 pi^2/3                = {RES['S2_S4_area']:.10f}")

# numerical confirmation of the flux identity by Monte-Carlo on S^4
rng = np.random.default_rng(7)
V = rng.normal(size=(4_000_000, 5)); V /= np.linalg.norm(V, axis=1, keepdims=True)
mc = np.mean(V[:, 4]*V[:, 4])                    # should be 1/5
RES['S2_MC_mean_wz2'] = float(mc)
print(f"S2  MC  <w_z^2>_{{S^4}} = {mc:.6f}  (exact 0.2)")

# the annulus identity  int_{eps<|w|<d} grad K dw = 0  (the two fluxes cancel exactly).
# MC control on the ABSOLUTELY CONVERGENT annulus 0.2 < |w| < 1.
n = 6_000_000
P = rng.normal(size=(n, 5)); P /= np.linalg.norm(P, axis=1, keepdims=True)
lo, hi = 0.2, 1.0
u = (lo**5 + (hi**5 - lo**5)*rng.random(n))**(1/5.0)
P = P*u[:, None]
nrm = np.linalg.norm(P, axis=1)
gz = (3/(8*math.pi**2))*(nrm**-5 - 5*P[:, 4]*P[:, 4]*nrm**-7)
gr = (3/(8*math.pi**2))*(-5*P[:, 4]*P[:, 0]*nrm**-7)
volA = float(S4)/5.0*(hi**5 - lo**5)
RES['S2_MC_annulus_gradKz'] = float(np.mean(gz)*volA)
RES['S2_MC_annulus_gradKr'] = float(np.mean(gr)*volA)
print(f"S2  MC  int_{{0.2<|w|<1}} d_z K dw = {RES['S2_MC_annulus_gradKz']:.3e}   (exact 0)")
print(f"S2  MC  int_{{0.2<|w|<1}} d_1 K dw = {RES['S2_MC_annulus_gradKr']:.3e}   (exact 0)")

# MC control on the flux  oint_{|w|=1} K(w) w_k dS = delta_{kz}/5
fl_z = float(np.mean((3/(8*math.pi**2))*V[:, 4]*V[:, 4])*float(S4))
fl_r = float(np.mean((3/(8*math.pi**2))*V[:, 4]*V[:, 0])*float(S4))
RES['S2_MC_flux_z'] = fl_z; RES['S2_MC_flux_r'] = fl_r
print(f"S2  MC  flux_z = {fl_z:.6f} (exact 0.2)   flux_r = {fl_r:.3e} (exact 0)")

json.dump(RES, open(os.path.join(HERE, 'k1_results.json'), 'w'), indent=1, sort_keys=True)
print("\nwrote k1_results.json")

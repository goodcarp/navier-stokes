#!/usr/bin/env python3
"""
u1 -- the exact algebra BLOCK 3 rests on.  Nothing here is quoted; everything is
re-derived in this file with sympy (exact) or checked numerically at high precision.

Items
  A. the 5-D lift:  div_5 b = 2a ;  grad_5 b = a*diag(1,1,1,1,-2) + E ;
     ||grad_5 b||_op <= 2|a| + r|grad a| + |omega^theta|      (this is L3v's (B.1))
  B. b(0) = 0  (so |b(x)| <= Gamma |x| by the mean value theorem)
  C. the transported-gradient equation has NO zeroth-order term:
        D_t (d_j eta) = -(d_j b).grad eta + nu Lap_5 (d_j eta)
     -- div_5 b = 2a does not appear (advective form).
  D. the weighted-Laplacian identities used by the two maximum principles
        Lap_5 (rho^k g) = rho^k Lap g + 2k rho^{k-2} x.grad g + k(k+3) rho^{k-2} g
     and their rearrangements for k = 1, 2.
  E. the Riesz composition constant  INT_{R^5} |x-x'|^{-4}|x'|^{-2} dx' = (pi^4/2)/|x| .
  F. C_K = 3/(8 pi^2), C_K |S^4| = 1, and 3 pi^2 / 16 = C_K * pi^4/2 .

Run:  python3 u1_algebra.py
"""
import json, math
import sympy as sp
import numpy as np
from scipy import integrate

OUT = {}

# ---------------------------------------------------------------- A, B, C, D
y1, y2, y3, y4, z = sp.symbols('y1 y2 y3 y4 z', real=True)
Y = sp.Matrix([y1, y2, y3, y4])
X5 = [y1, y2, y3, y4, z]
q = y1**2 + y2**2 + y3**2 + y4**2          # = r^2
r = sp.sqrt(q)

# an arbitrary axisymmetric no-swirl stream function psi1 = F(q, z);
# use a spanning family F = q**i * z**j with i, j symbolic-free but generic rationals
def lift(F):
    """return (a, b_vector, omega_theta, eta) for psi1 = F(q,z)."""
    a = -sp.diff(F, z)                      # a = u^r / r = -d_z psi1
    uz = 2*F + 2*q*sp.diff(F, q_sym)        # u^z = 2 psi1 + r d_r psi1 = 2F + 2q F_q
    b = [sp.simplify(a*Y[i]) for i in range(4)] + [sp.simplify(uz)]
    lap5 = sum(sp.diff(F_sub, v, 2) for v in X5)
    return a, b, lap5

q_sym = sp.symbols('q', positive=True)
zz = sp.symbols('zz', real=True)

# Build F as a function of (q, z) then substitute q -> |y|^2 to work in R^5.
i_, j_ = sp.Rational(3, 2), 3               # a generic non-polynomial member
Fqz = q_sym**i_ * zz**j_ + q_sym**2*zz + q_sym*zz**5 + zz**3   # generic enough
F_sub = Fqz.subs({q_sym: q, zz: z})

a_expr = -sp.diff(F_sub, z)
uz_expr = 2*F_sub + 2*q*sp.diff(Fqz, q_sym).subs({q_sym: q, zz: z})
b_vec = [sp.simplify(a_expr*Y[i]) for i in range(4)] + [sp.simplify(uz_expr)]

# A1: div_5 b = 2a
div5 = sum(sp.diff(b_vec[k], X5[k]) for k in range(5))
res_div = sp.simplify(div5 - 2*a_expr)
OUT['A1_div5b_minus_2a_residual'] = str(sp.simplify(res_div))

# A2: grad_5 b structure at the frame point y = (r,0,0,0).
frame = {y2: 0, y3: 0, y4: 0}
Gb = sp.Matrix(5, 5, lambda I, J: sp.diff(b_vec[I], X5[J]))
Gb_f = sp.simplify(Gb.subs(frame))
rr = sp.symbols('rr', positive=True)
Gb_f = sp.simplify(Gb_f.subs({y1: rr}))

a_f = sp.simplify(a_expr.subs(frame).subs({y1: rr}))
# d_r a and d_z a at the frame point:  d_r = d_{y1} there
ar_f = sp.simplify(sp.diff(a_expr, y1).subs(frame).subs({y1: rr}))
az_f = sp.simplify(sp.diff(a_expr, z).subs(frame).subs({y1: rr}))
# omega^theta = r * eta, eta = -Lap_5 psi1
eta_expr = -sum(sp.diff(F_sub, v, 2) for v in X5)
eta_f = sp.simplify(eta_expr.subs(frame).subs({y1: rr}))
om_f = sp.simplify(rr*eta_f)

E_pred = sp.zeros(5, 5)
E_pred[0, 0] = rr*ar_f
E_pred[0, 4] = rr*az_f
E_pred[4, 0] = rr*az_f - om_f
E_pred[4, 4] = -rr*ar_f
D_pred = a_f*sp.diag(1, 1, 1, 1, -2)
res_struct = sp.simplify(Gb_f - (D_pred + E_pred))
OUT['A2_grad5b_structure_residual'] = str(res_struct)

# A3: ||grad_5 b||_op <= 2|a| + r|grad a| + |omega^theta| -- numeric, 4000 random matrices
rng = np.random.default_rng(20260908)
worst = 0.0
for _ in range(4000):
    A, P, Q, W = rng.normal(size=4)          # a, r d_z a, r d_r a, omega
    Mx = A*np.diag([1., 1., 1., 1., -2.])
    Mx[0, 0] += Q
    Mx[0, 4] += P
    Mx[4, 0] += P - W
    Mx[4, 4] += -Q
    lhs = np.linalg.norm(Mx, 2)
    rhs = 2*abs(A) + math.hypot(P, Q) + abs(W)
    worst = max(worst, lhs/rhs)
OUT['A3_worst_ratio_op_over_bound'] = worst
OUT['A3_violated'] = bool(worst > 1 + 1e-12)

# B: b(0) = 0 for a z-odd psi1 (psi1 odd in z <=> eta odd in z).
Fodd = (q_sym**sp.Rational(3, 2))*zz + q_sym*zz**3 + zz**5
Fodd_sub = Fodd.subs({q_sym: q, zz: z})
uz_odd = 2*Fodd_sub + 2*q*sp.diff(Fodd, q_sym).subs({q_sym: q, zz: z})
a_odd = -sp.diff(Fodd_sub, z)
b0 = [sp.limit(sp.simplify(a_odd*Y[0]).subs({y2: 0, y3: 0, y4: 0, z: 0}), y1, 0),
      sp.simplify(uz_odd).subs({y1: 0, y2: 0, y3: 0, y4: 0, z: 0})]
OUT['B_b_at_origin'] = [str(sp.simplify(t)) for t in b0]

# C: no zeroth-order term in the gradient equation.
#    d_j( b . grad eta ) = (d_j b).grad eta + b . grad(d_j eta)  -- identity, and
#    d_j Lap_5 = Lap_5 d_j  in Cartesian R^5.  Verified on a generic pair.
etaT = sp.Function('e')(y1, y2, y3, y4, z)
lhs = sp.diff(sum(b_vec[k]*sp.diff(etaT, X5[k]) for k in range(5)), y1)
rhs = (sum(sp.diff(b_vec[k], y1)*sp.diff(etaT, X5[k]) for k in range(5))
       + sum(b_vec[k]*sp.diff(sp.diff(etaT, y1), X5[k]) for k in range(5)))
OUT['C_advective_split_residual'] = str(sp.simplify(lhs - rhs))
com = sp.simplify(sp.diff(sum(sp.diff(etaT, v, 2) for v in X5), y1)
                  - sum(sp.diff(sp.diff(etaT, y1), v, 2) for v in X5))
OUT['C_grad_commutes_with_Lap5_residual'] = str(com)

# D: the weighted Laplacian identities.
g = sp.Function('g')(y1, y2, y3, y4, z)
rho2 = y1**2 + y2**2 + y3**2 + y4**2 + z**2
rho = sp.sqrt(rho2)
def lap5(f):
    return sum(sp.diff(f, v, 2) for v in X5)
for k in (1, 2):
    lhs = lap5(rho**k*g)
    rhs = (rho**k*lap5(g) + 2*k*rho**(k-2)*sum(X5[m]*sp.diff(g, X5[m]) for m in range(5))
           + k*(k+3)*rho**(k-2)*g)
    OUT[f'D_lap5_rho{k}_g_residual'] = str(sp.simplify(lhs - rhs))

# rearranged forms actually used:
#   k=1 :  rho*Lap|eta| = Lap V - (2/rho^2) x.grad V - 2 V/rho^2 ,  V = rho|eta|
#   k=2 :  rho^2*Lap w  = Lap W - (4/rho^2) x.grad W - 2 W/rho^2 ,  W = rho^2 w
V = rho*g
lhsV = sp.simplify(rho*lap5(g))
rhsV = sp.simplify(lap5(V) - 2/rho2*sum(X5[m]*sp.diff(V, X5[m]) for m in range(5)) - 2*V/rho2)
OUT['D_k1_rearranged_residual'] = str(sp.simplify(lhsV - rhsV))
W = rho2*g
lhsW = sp.simplify(rho2*lap5(g))
rhsW = sp.simplify(lap5(W) - 4/rho2*sum(X5[m]*sp.diff(W, X5[m]) for m in range(5)) - 2*W/rho2)
OUT['D_k2_rearranged_residual'] = str(sp.simplify(lhsW - rhsW))

# ---------------------------------------------------------------- E, F
C_K = 3.0/(8.0*math.pi**2)
S4 = 8.0*math.pi**2/3.0
OUT['F_C_K'] = C_K
OUT['F_S4'] = S4
OUT['F_C_K_times_S4'] = C_K*S4
OUT['F_3pi2_over_16'] = 3*math.pi**2/16.0
OUT['F_C_K_times_pi4_over_2'] = C_K*math.pi**4/2.0

# E: Riesz composition in R^5:  INT |e - x'|^{-4} |x'|^{-2} dx' = pi^4/2  at |e| = 1.
# Reduce with the S^4 spherical average: dx' = t^4 dt dOmega_4,
# INT_{S^4} |e - t sigma|^{-4} dOmega_4 = |S^3| INT_0^pi sin^3(th) (1 - 2t cos th + t^2)^{-2} dth
S3 = 2.0*math.pi**2
def sph_avg_m4(t):
    f = lambda th: math.sin(th)**3/(1.0 - 2.0*t*math.cos(th) + t*t)**2
    v, _ = integrate.quad(f, 0.0, math.pi, limit=400, epsabs=1e-13, epsrel=1e-13,
                          points=[0.0, math.pi])
    return S3*v
def riesz_integrand(t):
    return sph_avg_m4(t)*t**4*t**(-2)
I1, _ = integrate.quad(riesz_integrand, 0.0, 1.0, limit=400, epsabs=1e-12, epsrel=1e-12,
                       points=[1.0])
I2, _ = integrate.quad(riesz_integrand, 1.0, np.inf, limit=400, epsabs=1e-12, epsrel=1e-12)
OUT['E_riesz_numeric'] = I1 + I2
OUT['E_riesz_closed_pi4_over_2'] = math.pi**4/2.0
OUT['E_riesz_rel_err'] = abs((I1 + I2) - math.pi**4/2.0)/(math.pi**4/2.0)

# a second, independent check of the same constant by the classical formula
# INT |x-y|^{-al}|y|^{-be} dy = C |x|^{n-al-be},
# C = pi^{n/2} G((n-al)/2)G((n-be)/2)G((al+be-n)/2) / (G(al/2)G(be/2)G(n-(al+be)/2))
from math import gamma as G
n, al, be = 5, 4, 2
Cclass = (math.pi**(n/2)*G((n-al)/2)*G((n-be)/2)*G((al+be-n)/2)
          / (G(al/2)*G(be/2)*G(n-(al+be)/2)))
OUT['E_riesz_classical_formula'] = Cclass

# bathtub constant for the two-octave collar:  R_A = (2^5 - 2^-5)^{1/5}
OUT['F_R_A_over_rho'] = (2.0**5 - 2.0**-5)**0.2
OUT['F_collar_coeff_4RA_over_2'] = 2.0*(2.0**5 - 2.0**-5)**0.2   # = 3.999218...

with open('u1_results.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
for k, v in OUT.items():
    print(f'{k:44s} {v}')

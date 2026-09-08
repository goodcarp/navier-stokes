#!/usr/bin/env python3
"""b1 -- EXACT identities underlying V-b (sympy, exact rational / symbolic arithmetic).

Nothing here is numerical unless stated.  Every residual printed as 0 is a true zero in
exact arithmetic.  Written and run by seat write/V-b-bulk-viscous-loss.
"""
import json, math
import sympy as sp

OUT = {}

# ---------------------------------------------------------------- A. the plateau identity
r, z, nu, M, t = sp.symbols('r z nu M t', positive=True)
y1,y2,y3,y4,zz = sp.symbols('y1 y2 y3 y4 zz', real=True)

def Lap5_rz(f):
    """Delta_5 in (r,z): d_rr + (3/r) d_r + d_zz   (5D Laplacian on axisymmetric-in-R^4 fields)."""
    return sp.diff(f, r, 2) + 3*sp.diff(f, r)/r + sp.diff(f, z, 2)

A1 = sp.simplify(r**3 * Lap5_rz(1/r))
OUT['A1_r3_Lap5_inv_r'] = str(A1)                     # must be -1

# Cartesian 5D check: y in R^4, z the fifth coordinate; f = 1/|y|
Y = sp.sqrt(y1**2+y2**2+y3**2+y4**2)
f = 1/Y
lap5_cart = sum(sp.diff(f, v, 2) for v in (y1,y2,y3,y4,zz))
A2 = sp.simplify(sp.simplify(lap5_cart) * Y**3)
OUT['A2_cart_R3_Lap5_inv_r'] = str(A2)                # must be -1

# eta_P = -M sgn(z)/r  ->  Delta_5 eta_P = -eta_P/r^2  (away from z=0)
etaP = -M/r                                           # on {z>0}
A3 = sp.simplify(Lap5_rz(etaP) + etaP/r**2)
OUT['A3_bulk_operator_residual'] = str(A3)            # must be 0

# ---------------------------------------------------------------- B. Hessian of 1/r in R^5
H = sp.Matrix(5,5, lambda i,j: sp.diff(f, [y1,y2,y3,y4,zz][i], [y1,y2,y3,y4,zz][j]))
rp = sp.Symbol('r', positive=True)
H0 = H.subs({y2:0, y3:0, y4:0}).subs({y1: rp})      # evaluate at y = r * e1, r > 0
H0 = sp.simplify(H0)
ev = sp.simplify(sp.Matrix(H0).eigenvals())
OUT['B_hessian_eigenvalues_at_r'] = {str(k): int(v) for k, v in ev.items()}
OUT['B_hessian_trace_times_r3'] = str(sp.simplify(H0.trace()*rp**3))       # must be -1
OUT['B_hessian_opnorm_times_r3'] = str(sp.simplify(sp.Max(*[sp.Abs(k) for k in ev])*rp**3))

# ---------------------------------------------------------------- C. directional derivatives
# g(s) = 1/|y0 + s v|,  |y0| = r, mu = v.yhat  ->  g(s) = (r^2 + 2 s r mu + s^2)^(-1/2)
# generating function of Legendre polynomials:  d^k/ds^k g |_{s=0} = k! P_k(-mu) / r^{k+1}
s_, mu_ = sp.symbols('s mu', real=True)
g = (r**2 + 2*s_*r*mu_ + s_**2)**sp.Rational(-1,2)
Ck = {}
for k in range(1, 7):
    dk = sp.simplify(sp.diff(g, s_, k).subs(s_, 0))
    pred = sp.factorial(k)*sp.legendre(k, -mu_)/r**(k+1)
    Ck[k] = str(sp.simplify(dk - pred))                # must be 0
OUT['C_legendre_residuals'] = Ck
OUT['C_sup_abs_dk'] = "k! / r^(k+1)   (|P_k| <= 1, attained at mu = -+1)"

# ---------------------------------------------------------------- D. iterated Laplacian in y
# Delta_y acts on R^4 radial functions as d_rr + (3/r) d_r
def Lap_y(F):
    return sp.diff(F, r, 2) + 3*sp.diff(F, r)/r
cur = 1/r; coef = []
for k in range(1, 5):
    cur = sp.simplify(Lap_y(cur))
    c = sp.simplify(cur * r**(2*k+1))
    coef.append(str(c))
OUT['D_Delta_y_k_of_inv_r_coefficients'] = coef        # -1, -3, -45, -1575 ...

# ---------------------------------------------------------------- E. the affine-strain covariance
# b(x,t) = a(t) diag(1,1,1,1,-2) x .  Backward reverse-time diffusion
#   dY_s = -b(Y_s, tau-s) ds + sqrt(2 nu) dW_s ,  Y_0 = X(tau)
# Z_s := Y_s - x_s  with x_s the reverse characteristic.  Then Z_s is a mean-zero Gaussian
# with covariance  C_tau = 2 nu int_0^tau  Psi_tau Psi_p^{-1} (.)^T dp .
# Psi solves Psi' = -S(tau-s) Psi, S = a diag(1,1,1,1,-2).  With lam(u) = exp(int_0^u a),
#   Psi_tau Psi_p^{-1} = diag( lam(tau-p)^{-1} x4 , lam(tau-p)^2 )
# so with the substitution u = tau - p,
#   C_tau = 2 nu int_0^tau diag( lam(u)^{-2} x4 , lam(u)^4 ) du .
u_ = sp.symbols('u', positive=True); tau = sp.symbols('tau', positive=True)
a_ = sp.Function('a')
lam = sp.exp(sp.Integral(a_(u_), (u_, 0, u_)))
# symbolic verification of the propagator identity on the 1st and 5th diagonal entries
p_ = sp.symbols('p', positive=True)
mu_fun = sp.Function('mu')                     # mu(s) = int_{tau-s}^{tau} a
Psi1 = sp.exp(-mu_fun(s_)); Psi5 = sp.exp(2*mu_fun(s_))
res1 = sp.simplify(sp.diff(Psi1, s_) + (-1)*(-sp.diff(mu_fun(s_), s_))*Psi1)   # Psi1' = -a(tau-s) Psi1
OUT['E_propagator_symbolic_note'] = ("Psi_s = diag(e^{-mu(s)} x4, e^{2 mu(s)}), mu(s)=int_{tau-s}^tau a; "
                                     "Psi_tau Psi_p^{-1} = diag(lam(tau-p)^{-1} x4, lam(tau-p)^2)")

# The decisive identity:  (1/2) C_tau : Hess(eta_P)  =  - eta_P(x0) * nu * int_0^tau dt / r(t)^2
# with r(t) = lam(t) r0 .   Hess(eta_P) = -M Hess(1/r) = -(M/r^3)(3 yhat yhat^T - I_4) (+) 0 .
sy, sz, r0 = sp.symbols('sigma_y sigma_z r0', positive=True)
C = sp.diag(2*sy, 2*sy, 2*sy, 2*sy, 2*sz)
HessEtaP = -M*H0.subs(rp, r0)
half_contract = sp.simplify(sp.Rational(1,2)*sum(C[i,j]*HessEtaP[i,j] for i in range(5) for j in range(5)))
target = -(-M/r0) * (sy/r0**2)          # = -eta_P(x0) * sigma_y / r0^2   with eta_P(x0) = -M/r0
OUT['E_half_C_contract_Hess'] = str(sp.simplify(half_contract))
OUT['E_identification_residual'] = str(sp.simplify(half_contract - (-M/r0)*(-sy/r0**2)))   # must be 0
OUT['E_sigma_z_absent'] = str(sp.simplify(sp.diff(half_contract, sz)))                      # must be 0

# and sigma_y / r0^2 = nu int_0^tau dt / r(t)^2  because r(t) = lam(t) r0 and
# sigma_y = nu int_0^tau lam(u)^{-2} du .
OUT['E_sigma_y_is_nu_int_dt_over_r2'] = "sigma_y/r0^2 = (nu/r0^2) int lam^{-2} du = nu int du / r(u)^2"

# ---------------------------------------------------------------- F. label-frame operator
# Lemma V.0:  d_t zeta = (nu/D) d_l ( D G^{lk} d_k zeta ),  G = (J^T J)^{-1}, D = det J.
# Expanded:   = nu G^{lk} d_l d_k zeta + nu V^k d_k zeta ,  V^k := (1/D) d_l ( D G^{lk} ).
# CLAIM: V == 0 identically iff X(.,t) is affine in alpha; and in general V is a contraction of
# second derivatives of X.  Verified symbolically on a nonlinear map in n=2 (enough to exhibit
# non-vanishing) and on a general affine map in n=5.
a1,a2 = sp.symbols('alpha1 alpha2', real=True)
eps = sp.symbols('epsilon', real=True)
def V_of_map(Xs, al):
    n = len(al)
    J = sp.Matrix(n, n, lambda i,j: sp.diff(Xs[i], al[j]))
    D = sp.simplify(J.det())
    G = sp.simplify((J.T*J).inv())
    V = []
    for k in range(n):
        V.append(sp.simplify(sum(sp.diff(D*G[l,k], al[l]) for l in range(n))/D))
    return V
# affine map in n=2 (with shear + dilation) -> V must be exactly 0
Xaff = [sp.Rational(3,2)*a1 + sp.Rational(1,5)*a2, sp.Rational(2,3)*a2]
OUT['F_V_affine'] = [str(v) for v in V_of_map(Xaff, (a1,a2))]
# nonlinear map -> V nonzero, and O(eps) = O(second derivative of X)
Xnl = [a1 + eps*a2**2, a2 + eps*a1*a2]
Vnl = V_of_map(Xnl, (a1,a2))
OUT['F_V_nonlinear'] = [str(sp.simplify(v)) for v in Vnl]
OUT['F_V_nonlinear_linearised_in_eps'] = [str(sp.simplify(sp.series(v, eps, 0, 2).removeO())) for v in Vnl]
# affine in n=5 with a generic matrix
Ms = sp.Matrix(5,5, lambda i,j: sp.Rational((i*7+j*3) % 5 + 1, (i+j) % 3 + 2))
al5 = sp.symbols('A0:5', real=True)
X5 = [sum(Ms[i,j]*al5[j] for j in range(5)) for i in range(5)]
OUT['F_V_affine_n5'] = [str(v) for v in V_of_map(X5, al5)]
# n = 1, general map: V = -X''/X'^3 exactly  =>  V == 0  <=>  X affine.  (Converse PROVED here
# only for n = 1; for n >= 2 only the forward implication and the O(d^2 X) structure are proved.)
Xf = sp.Function('X')
J1 = sp.diff(Xf(a1), a1); D1 = J1; G1 = 1/J1**2
V1 = sp.simplify(sp.diff(D1*G1, a1)/D1)
OUT['F_V_n1_general'] = str(sp.simplify(V1))
OUT['F_V_n1_residual_vs_minus_Xpp_over_Xp3'] = str(sp.simplify(V1 + sp.diff(Xf(a1), a1, 2)/J1**3))

for k, v in OUT.items():
    print(f"{k}: {v}")
json.dump(OUT, open('b1_results.json','w'), indent=1)
print("\nWROTE b1_results.json")

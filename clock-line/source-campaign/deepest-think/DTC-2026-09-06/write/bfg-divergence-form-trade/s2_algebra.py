#!/usr/bin/env python3
"""
s2_algebra.py -- the exact pointwise constant in the divergence-form Duhamel term.

With div u = div omega = 0 the two nonlinear terms of the vorticity equation are

    u_i d_i omega_j - omega_i d_i u_j = d_i A_{ij},   A_{ij} = u_i omega_j - omega_i u_j,

and A = u (x) omega - omega (x) u is ANTISYMMETRIC.  Moving the derivative onto the
heat kernel, the j-th component of the Duhamel integrand is sum_i (d_i G) A_{ij}
= [A^T grad G]_j = -[A grad G]_j, with  A v = (v.omega) u - (v.u) omega.  Hence the
whole vector is controlled by

    Q := sup { |(n.u) omega - (n.omega) u| / (|u||omega|) : |n| = 1, u,omega != 0 }.

CLAIM (A1):  Q = 1 exactly.  (The triangle inequality only gives 2.)
Proof.  a = u/|u|, b = omega/|omega|, c = a.b, alpha = n.a, beta = n.b.  Then
   |(n.u)omega - (n.omega)u|^2/(|u|^2|omega|^2) = alpha^2 + beta^2 - 2 c alpha beta.
For a UNIT n the pair (alpha,beta) ranges over exactly the ellipse
   {(alpha,beta) : alpha^2 - 2 c alpha beta + beta^2 <= 1 - c^2}
(Gram feasibility: the 3x3 Gram matrix of (a,b,n) must be PSD, and in dim >= 3 every
such (alpha,beta) is attainable).  The objective IS that quadratic form, so
   sup = 1 - c^2 <= 1,  attained at c = 0 (u perpendicular omega) on the boundary.
Checked below symbolically (Gram determinant) and by numerical maximisation.
"""
import json
import numpy as np
import sympy as sp

out = {}
rng = np.random.default_rng(20260906)

# ---- 1. antisymmetry and the action of A ---------------------------------
u1,u2,u3,w1,w2,w3,v1,v2,v3 = sp.symbols('u1 u2 u3 w1 w2 w3 v1 v2 v3', real=True)
u = sp.Matrix([u1,u2,u3]); w = sp.Matrix([w1,w2,w3]); v = sp.Matrix([v1,v2,v3])
A = u*w.T - w*u.T
r1 = sp.expand(A + A.T)
print("A + A^T                          :", list(r1)); assert r1 == sp.zeros(3,3)
r2 = sp.expand(A*v - ((v.dot(w))*u - (v.dot(u))*w))
print("A v - [(v.w)u - (v.u)w]          :", list(r2.T)); assert r2 == sp.zeros(3,1)
out['A_antisym_residual'] = 0; out['Av_identity_residual'] = 0

# ---- 2. the divergence form, on explicit polynomial divergence-free fields
x,y,z = sp.symbols('x y z', real=True); X = (x,y,z)
U = sp.Matrix([x**2*y - z**3, -2*x*y**2 + y*z, x*z**2 + 3*y**2*z - x**3])
W = sp.Matrix([y*z**2 + x**3*y, -x*z**2 - x**2*y**2*sp.Rational(3,2), 5*x*y**3 - z*y*z])
# force divergence-free by removing the divergence via a corrector in the 3rd slot
U = sp.Matrix([U[0], U[1], sp.integrate(-sp.diff(U[0],x)-sp.diff(U[1],y), z)])
W = sp.Matrix([W[0], W[1], sp.integrate(-sp.diff(W[0],x)-sp.diff(W[1],y), z)])
dU = sp.expand(sum(sp.diff(U[i],X[i]) for i in range(3)))
dW = sp.expand(sum(sp.diff(W[i],X[i]) for i in range(3)))
print("div U, div W                     :", dU, dW); assert dU == 0 and dW == 0
lhs = [sp.expand(sum(U[i]*sp.diff(W[j],X[i]) for i in range(3)) - sum(W[i]*sp.diff(U[j],X[i]) for i in range(3))) for j in range(3)]
rhs = [sp.expand(sum(sp.diff(U[i]*W[j] - W[i]*U[j], X[i]) for i in range(3))) for j in range(3)]
dres = [sp.expand(lhs[j]-rhs[j]) for j in range(3)]
print("u.grad w - w.grad u - div(A)     :", dres); assert all(d == 0 for d in dres)
out['divergence_form_residual'] = [str(d) for d in dres]

# ---- 3. Q = 1: the Gram/feasibility statement, symbolically ---------------
al, be, c = sp.symbols('alpha beta c', real=True)
Gram = sp.Matrix([[1, c, al],[c, 1, be],[al, be, 1]])       # Gram of (a,b,n), all unit
det = sp.expand(Gram.det())
print("det Gram(a,b,n)                  =", det)
# det >= 0 (PSD) <=> alpha^2 - 2c alpha beta + beta^2 <= 1 - c^2
constraint = sp.expand((1-c**2) - (al**2 - 2*c*al*be + be**2))
print("(1-c^2) - [a^2-2c ab+b^2]        =", constraint, "   det - that =", sp.expand(det - constraint))
assert sp.expand(det - constraint) == 0
objective = sp.expand(al**2 + be**2 - 2*c*al*be)
print("objective - constraint form      =", sp.expand(objective - (al**2 - 2*c*al*be + be**2)))
assert sp.expand(objective - (al**2 - 2*c*al*be + be**2)) == 0
print("=> sup objective = 1 - c^2 <= 1, attained at c = 0.  Q = 1 EXACTLY.")
out['Q_symbolic'] = "sup = 1 - c^2 <= 1"

# ---- 4. Q = 1: vectorised numerical maximisation --------------------------
N = 4000000
uu = rng.normal(size=(N,3)); ww = rng.normal(size=(N,3)); nn = rng.normal(size=(N,3))
nn /= np.linalg.norm(nn, axis=1)[:,None]
num = np.linalg.norm((nn*uu).sum(1)[:,None]*ww - (nn*ww).sum(1)[:,None]*uu, axis=1)
den = np.linalg.norm(uu,axis=1)*np.linalg.norm(ww,axis=1)
vals = num/den
print("Q  random max over %d samples   = %.15f   (claimed exactly 1)" % (N, vals.max()))
assert vals.max() <= 1 + 1e-12
out['Q_random_max'] = float(vals.max())
# exact attainment: u = e1, omega = e2, n = (e1+e2)/sqrt(2)
e1 = np.array([1.,0,0]); e2 = np.array([0,1.,0]); nq = (e1+e2)/np.sqrt(2)
att = np.linalg.norm(np.dot(nq,e1)*e2 - np.dot(nq,e2)*e1)/1.0
print("Q  at u=e1, omega=e2, n=(e1+e2)/sqrt2 = %.15f  -> the sup is attained" % att)
assert abs(att - 1.0) < 1e-15
out['Q_attained'] = float(att); out['Q_exact'] = 1.0; out['Q_triangle'] = 2.0

json.dump(out, open('s2_results.json','w'), indent=1)
print("\ns2: divergence form verified; Q = 1 exactly (Gram argument), attained")

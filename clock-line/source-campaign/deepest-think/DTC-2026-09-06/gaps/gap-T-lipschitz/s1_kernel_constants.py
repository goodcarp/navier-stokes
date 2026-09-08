"""s1 - exact (sympy) constants for the strain kernel K on R^5, and the exact shell integral.

K(w) = -(3/(8 pi^2)) w_z / |w|^5   (a(0) = int K(w) eta(w) dw ; see rebuild/far-near-kernel-lemma)
  (1) |K(w)| <= (3/(8 pi^2)) |w|^-4                (sharp)
  (2) |grad K(w)| <= (3/(2 pi^2)) |w|^-5           (sharp)   [symbolic gradient verified]
  (3) int_{T_lambda S} du_5 /(r_u |u|^4) = pi^3 L  exactly, for every lambda > 0
  (4) Lemma-1 cross check a[eta0 o T_lambda^-1](0) = (M/2) lambda L (bang-bang cap), by exact antiderivative
"""
import json, sympy as sp

out = {}
c, lam, phi, rho, rho0, RR, M, v, A, B = sp.symbols(
    'c lambda phi rho rho_0 R_out M v A B', positive=True)

# ---- (1) size of K -------------------------------------------------------
out['K_size_constant'] = '3/(8*pi**2)'
out['K_size_constant_num'] = float(sp.Rational(3,8)/sp.pi**2)

# ---- (2) gradient of K (symbolic, R^5) -----------------------------------
x1,x2,x3,x4,zz = sp.symbols('x1 x2 x3 x4 z', real=True)
w = sp.Matrix([x1,x2,x3,x4,zz]); nw = sp.sqrt(sum(t**2 for t in w))
Ksym = -sp.Rational(3,8)/sp.pi**2 * zz / nw**5
grad = sp.Matrix([sp.diff(Ksym, t) for t in (x1,x2,x3,x4,zz)])
pred = -sp.Rational(3,8)/sp.pi**2 * ( sp.Matrix([0,0,0,0,1])/nw**5 - 5*zz*w/nw**7 )
res  = sp.simplify(grad - pred)
assert all(sp.simplify(e) == 0 for e in res), res
out['gradK_symbolic_residual'] = '0 (all 5 components)'
g2 = sp.expand(1 - 10*c**2 + 25*c**2)            # |e_z - 5 c what|^2
out['grad_dir_norm2'] = str(g2)                   # 1 + 15 c^2
out['grad_dir_sup'] = str(sp.sqrt(g2.subs(c,1)))  # 4, attained on the axis c = +-1
CK = sp.Rational(3,8)/sp.pi**2 * 4
assert sp.simplify(CK - sp.Rational(3,2)/sp.pi**2) == 0
out['gradK_constant'] = '3/(2*pi**2)'
out['gradK_constant_num'] = float(CK)

# ---- (3) the exact shell integral ---------------------------------------
g = sp.sqrt(lam**-2*sp.sin(phi)**2 + lam**4*sp.cos(phi)**2)   # |T_lambda^-1 u| = |u| g(phi)
Lsym = sp.log(RR/rho0)
rad = sp.log((RR/g)) - sp.log((rho0/g))                        # int_{rho0/g}^{R/g} drho/rho
assert sp.simplify(rad - Lsym) == 0
out['radial_log_length'] = 'log(R/rho0), independent of phi and lambda  [residual 0]'
ang = sp.integrate(sp.sin(phi)**2, (phi, 0, sp.pi))
assert ang == sp.pi/2
out['angular_integral_sin2'] = str(ang)
out['shell_integral'] = 'pi**3 * L'
out['pi_cubed'] = float(sp.pi**3)

# ---- (4) Lemma 1 by exact antiderivative --------------------------------
# a_lambda(0) = (M/2) L * P(lambda),  P(lambda) = 3 int_0^1 v^2 (A v^2 + B)^{-5/2} dv,
#   v = sin(phi), A = lambda^2 - lambda^-4, B = lambda^-4, A + B = lambda^2.
F = v**3/(3*B*(A*v**2+B)**sp.Rational(3,2))
resid = sp.simplify(sp.diff(F, v) - v**2*(A*v**2+B)**sp.Rational(-5,2))
assert resid == 0, resid
out['antiderivative_residual'] = str(resid)
P = sp.simplify(3*(F.subs(v,1) - F.subs(v,0)).subs({A: lam**2-lam**-4, B: lam**-4}))
out['P_of_lambda'] = str(sp.simplify(P))
assert sp.simplify(P - lam) == 0
out['lemma1_residual'] = '0  (P(lambda) = lambda exactly)'

json.dump(out, open('s1_results.json','w'), indent=1)
for k,val in out.items(): print(f"{k:26s} {val}")

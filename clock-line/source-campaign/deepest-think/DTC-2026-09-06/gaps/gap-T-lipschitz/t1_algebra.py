"""
t1_algebra.py -- gap-T-lipschitz.  Exact/symbolic ingredients of LEMMA T.

Establishes, with sympy (exact) plus numerical cross-checks:
  (A) K(w) = (3/(8 pi^2)) w_z |w|^{-5} is homogeneous of degree -4 on R^5,
      and grad K obeys  |grad K(w)| = (3/(8 pi^2)) |w|^{-5} sqrt(1+15 t^2), t = w_z/|w|,
      so  sup_w |w|^5 |grad K(w)| = 4*(3/(8 pi^2)) = 3/(2 pi^2)   [C_K, EXACT]
  (B) the 5D lift Jacobian of an axisymmetric map (r,z) -> (F,G) is
      J5 = (F/r)^3 (F_r G_z - F_z G_r)        [checked against a numerical 5x5 det]
  (C) the shell weight identity  INT_S (M/r) |x|^{-4} dx5 = pi^3 M log(R/rho0)
  (D) I2(lambda) := INT_0^pi |cos p| sin^2 p  q_l(p)^{-5} dp = (2/3) lambda, EXACT,
      where q_l(p)^2 = lambda^2 sin^2 p + lambda^{-4} cos^2 p.
      (D) is equivalent to Lemma 1 of lower/prove-lagrangian: a_lambda(0) = (M/2) lambda L.
Every printed number is computed here.
"""
import json
import numpy as np
import sympy as sp

out = {}

# ---------------------------------------------------------------- (A) kernel
w1,w2,w3,w4,wz = sp.symbols('w1 w2 w3 w4 wz', real=True)
W = sp.Matrix([w1,w2,w3,w4,wz])
nw = sp.sqrt(sum(c**2 for c in W))
K  = sp.Rational(3,1)/(8*sp.pi**2) * wz / nw**5

s = sp.symbols('s', positive=True)
homog = sp.simplify(K.subs({w1:s*w1,w2:s*w2,w3:s*w3,w4:s*w4,wz:s*wz}) - s**(-4)*K)
out['K_homogeneity_residual_deg_-4'] = str(sp.simplify(homog))

gradK = sp.Matrix([sp.diff(K,v) for v in (w1,w2,w3,w4,wz)])
# claimed:  gradK = (3/8pi^2) |w|^{-5} ( e_z - 5 t u ),  u = w/|w|, t = wz/|w|
t_  = wz/nw
U   = W/nw
claim = sp.Rational(3,1)/(8*sp.pi**2)*nw**(-5)*(sp.Matrix([0,0,0,0,1]) - 5*t_*U)
out['gradK_form_residual'] = str(sp.simplify(sp.expand(gradK-claim)))
# |e_z - 5 t u|^2 = 1 + 15 t^2
tt = sp.symbols('t', real=True)
nrm2 = sp.simplify((sp.Matrix([0,0,0,0,1]) - 5*t_*U).dot(sp.Matrix([0,0,0,0,1]) - 5*t_*U))
out['grad_dir_norm2_minus_1_plus_15t2'] = str(sp.simplify(nrm2 - (1+15*t_**2)))
CK = sp.Rational(3,1)/(8*sp.pi**2)*sp.sqrt(1+15*1)          # t^2 = 1 is the max
out['C_K_symbolic'] = str(sp.nsimplify(sp.simplify(CK)))
out['C_K_value'] = float(CK)
assert sp.simplify(CK - 3/(2*sp.pi**2)) == 0

# ---------------------------------------------------------------- (B) Jacobian formula
def J5_numeric(F,G,r,z,h=1e-6):
    """det of the true 5x5 Jacobian of  y -> (F(|y|,z)/|y|) y ,  z -> G(|y|,z),
       evaluated by central differences at a generic 5D point with |y|=r."""
    # place y generically on the 3-sphere of radius r
    d = np.array([0.31,-0.52,0.77,0.19]); d = d/np.linalg.norm(d)
    y0 = r*d; x0 = np.concatenate([y0,[z]])
    def Phi5(x):
        y = x[:4]; zz = x[4]; rr = np.linalg.norm(y)
        return np.concatenate([F(rr,zz)/rr*y, [G(rr,zz)]])
    Jm = np.zeros((5,5))
    for j in range(5):
        e = np.zeros(5); e[j] = h
        Jm[:,j] = (Phi5(x0+e)-Phi5(x0-e))/(2*h)
    return np.linalg.det(Jm)

lam = 1.37
F  = lambda r,z: lam*r + 0.21*r*z/np.hypot(r,z)
G  = lambda r,z: z/lam**2 + 0.13*(z*z-r*r)/np.hypot(r,z)
rs, zs = 0.83, 1.41
# analytic (F/r)^3 (F_r G_z - F_z G_r) by finite differences of F,G in (r,z)
h=1e-6
Fr=(F(rs+h,zs)-F(rs-h,zs))/(2*h); Fz=(F(rs,zs+h)-F(rs,zs-h))/(2*h)
Gr=(G(rs+h,zs)-G(rs-h,zs))/(2*h); Gz=(G(rs,zs+h)-G(rs,zs-h))/(2*h)
J_formula = (F(rs,zs)/rs)**3*(Fr*Gz-Fz*Gr)
J_true    = J5_numeric(F,G,rs,zs)
out['J5_formula'] = J_formula; out['J5_numeric_5x5_det'] = J_true
out['J5_rel_diff'] = abs(J_formula-J_true)/abs(J_true)
# and for T_lambda: J5 = lambda^2
FT = lambda r,z: lam*r; GT = lambda r,z: z/lam**2
out['J5_Tlambda_numeric'] = J5_numeric(FT,GT,rs,zs)
out['J5_Tlambda_exact']   = lam**2

# ---------------------------------------------------------------- (C) shell weight
rho,phi,M,L,rho0,R = sp.symbols('rho phi M L rho0 R', positive=True)
# dx5 = 2 pi^2 rho^4 sin^3(phi) drho dphi ;  r = rho sin phi
integrand = (M/(rho*sp.sin(phi)))*rho**(-4)*2*sp.pi**2*rho**4*sp.sin(phi)**3
val = sp.integrate(sp.integrate(integrand,(phi,0,sp.pi)),(rho,rho0,R))
out['shell_weight_integral'] = str(sp.simplify(val))
out['shell_weight_minus_pi3_M_L'] = str(sp.simplify(val - sp.pi**3*M*sp.log(R/rho0)))

# ---------------------------------------------------------------- (D) I2 = (2/3) lambda^3
lm = sp.symbols('lambda', positive=True)
v  = sp.symbols('v', positive=True)           # v = sin phi
A  = lm**2 - lm**(-4); B = lm**(-4)
# on [0,pi/2]: |cos| sin^2 q^{-5} dphi ; put v = sin phi, dphi = dv/cos phi
#   -> v^2 (A v^2 + B)^{-5/2} dv        (the |cos| cancels the Jacobian)
anti = v**3/(3*B*(A*v**2+B)**sp.Rational(3,2))
res  = sp.simplify(sp.diff(anti,v) - v**2*(A*v**2+B)**sp.Rational(-5,2))
out['I2_antiderivative_residual'] = str(sp.simplify(res))
I2_sym = sp.simplify(2*anti.subs(v,1))        # two hemispheres
out['I2_closed_form'] = str(sp.simplify(I2_sym))
out['I2_minus_two_thirds_lambda'] = str(sp.simplify(I2_sym - sp.Rational(2,3)*lm))
# Lemma-1 consistency:  a_lambda(0) = (3/4) M L I2(lambda) = (M/2) lambda L
out['lemma1_consistency  a_lam/(M L) - lambda/2'] = str(sp.simplify(sp.Rational(3,4)*I2_sym - lm/2))
# and the Term-II mass:  INT |eta0| |K(T_lam x)| dx5 = (3/4) M lambda^{-2} L I2 = (1/2) lambda^{-1} M L
out['termII_mass_over_ML  minus 1/(2 lambda)'] = str(sp.simplify(sp.Rational(3,4)*lm**-2*I2_sym - 1/(2*lm)))

# numerical confirmation of I2 by my own Gauss-Legendre quadrature
def gl(n):
    x,w = np.polynomial.legendre.leggauss(n); return x,w
def I2_num(lamv,n=400):
    x,w = gl(n); tot=0.0
    for a,b in ((0.0,np.pi/2),(np.pi/2,np.pi)):
        p = 0.5*(b-a)*x+0.5*(a+b); ww = 0.5*(b-a)*w
        q = np.sqrt(lamv**2*np.sin(p)**2 + lamv**-4*np.cos(p)**2)
        tot += np.sum(ww*np.abs(np.cos(p))*np.sin(p)**2*q**-5)
    return tot
out['I2_numeric_check'] = {f'{lv}':[I2_num(lv), 2/3*lv, abs(I2_num(lv)-2/3*lv)]
                           for lv in (1.0,1.1,1.25,1.5,2.0,0.8)}

print(json.dumps(out, indent=1, default=str))
json.dump(out, open('t1_results.json','w'), indent=1, default=str)

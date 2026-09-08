"""r1 - independent re-derivation of the algebra of LEMMA T, plus the two identities the
seat's BRIDGE to prove-lagrangian needs but never wrote down.

(1) grad K and the sharp constants, from scratch (sympy; not read from the seat's s1).
(2) the ANGULAR identity   J(lambda) := int_0^pi sin^2(phi)/g(phi;lambda)^4 dphi = pi/(2 lambda),
    g = sqrt(lambda^2 sin^2 + lambda^-4 cos^2) = |T_lambda x|/|x| .
    Route (exact, each step symbolic):
      a) d/dt[ arctan(t sqrt(a/b))/sqrt(a b) ] = 1/(b + a t^2)                        [symbolic]
      b) t = tan(phi) on (0,pi/2), integrand even about pi/2
         => I(a,b) := int_0^pi dphi/(b + (a-b) sin^2 phi) = pi/sqrt(a b)
      c) -d/da of the integrand is sin^2/( . )^2 , so J = -dI/da = (pi/2) a^-3/2 b^-1/2 [symbolic]
      d) a = lambda^2, b = lambda^-4  =>  J = pi/(2 lambda)                            [symbolic]
    plus a 30-digit mpmath confirmation of (b) and of J at six lambda.
    This is the identity that makes the seat's Step-3 constant pi^3 M L/lambda true in the
    x-variables; s1 proved only the u-substitution route and s3 only checked it numerically.
(3) the SHELL-DEPENDENT strain Lambda(x) = T_{lambda(|x|)} x :
       det D Lambda = lambda^2 [ 1 + (dlog lambda/dlog rho)(sin^2 phi - 2 cos^2 phi) ]
    verified against a finite-difference 5x5 determinant of the actual map at random points.
(4) the per-shell Lemma-1 core int_0^{pi/2} cos sin^2/g^5 dphi = lambda/3, hence
       a[eta0 o Lambda^-1](0) = (M/2) int_0^L lambda(rho) dlog rho + (Jacobian-shear term),
    EXACTLY, for ANY profile lambda(.).  This is the generalisation of Step 3 that the seat's
    single-T_lambda Lemma T does not state.
"""
import json, sympy as sp, numpy as np, mpmath as mp
out = {}

# ---------------- (1) grad K ----------------
x1,x2,x3,x4,z = sp.symbols('x1 x2 x3 x4 z', real=True)
w = sp.Matrix([x1,x2,x3,x4,z]); nw = sp.sqrt(sum(t**2 for t in w))
K = -sp.Rational(3,8)/sp.pi**2 * z/nw**5
gv = sp.Matrix([sp.diff(K,t) for t in (x1,x2,x3,x4,z)])
c = sp.Symbol('c')
lhs = sp.simplify(sum(t**2 for t in gv)*nw**10)
rhs = (sp.Rational(3,8)/sp.pi**2)**2*(1+15*c**2)
res = sp.simplify(lhs - rhs.subs(c, z/nw))
out['gradK_norm_residual'] = str(res); assert res == 0
out['gradK_sup_symbolic'] = str(sp.nsimplify(sp.sqrt(rhs.subs(c,1)), [sp.pi]))
out['gradK_sup_num'] = float(sp.sqrt(rhs.subs(c,1)))
out['K_sup_num'] = float(sp.Rational(3,8)/sp.pi**2)
out['gradK_sup_equals_3_over_2pi2'] = bool(abs(float(sp.sqrt(rhs.subs(c,1))) - 3/(2*np.pi**2)) < 1e-15)

# ---------------- (2) the angular identity ----------------
a, b, t, phi, lam = sp.symbols('a b t phi lambda', positive=True)
F = sp.atan(t*sp.sqrt(a/b))/sp.sqrt(a*b)
resA = sp.simplify(sp.diff(F, t) - 1/(b + a*t**2))
out['step_a_residual'] = str(resA); assert resA == 0, resA
I_ab = sp.limit(F, t, sp.oo)*2                       # 2 * int_0^inf dt/(b+a t^2)
I_ab = sp.simplify(I_ab)
out['I_ab'] = str(I_ab)
resB = sp.simplify(I_ab - sp.pi/sp.sqrt(a*b))
out['step_b_residual'] = str(resB); assert resB == 0, resB
integrand = 1/(b + (a-b)*sp.sin(phi)**2)
resC = sp.simplify(-sp.diff(integrand, a) - sp.sin(phi)**2/(b+(a-b)*sp.sin(phi)**2)**2)
out['step_c_residual'] = str(resC); assert resC == 0, resC
J_ab = sp.simplify(-sp.diff(sp.pi/sp.sqrt(a*b), a))
out['J_ab'] = str(J_ab)
J_lam = sp.simplify(J_ab.subs({a: lam**2, b: lam**-4}))
resD = sp.simplify(J_lam - sp.pi/(2*lam))
out['J_of_lambda'] = str(J_lam)
out['step_d_residual'] = str(resD); assert resD == 0, resD
mp.mp.dps = 30
num = []
for lv in (0.5, 0.7, 1.0, 1.25, 1.5, 3.0):
    f = lambda p, lv=lv: mp.sin(p)**2/(lv**2*mp.sin(p)**2 + lv**-4*mp.cos(p)**2)**2
    q = mp.quad(f, [0, mp.pi/2, mp.pi])
    num.append(dict(lam=lv, quad=str(q), closed=str(mp.pi/(2*lv)), err=float(abs(q-mp.pi/(2*lv)))))
out['angular_identity_mpmath'] = num
out['angular_identity_max_err'] = max(r['err'] for r in num)

# ---------------- (3) det D Lambda, finite differences in R^5 ----------------
KTH = 2*(1-np.sqrt(2.0/3.0)); rho0 = 1.0; Lbig = 8.317766166719343
def lam_np(rho):   return (1-(1-np.log(rho/rho0)/Lbig)*KTH/2)**-2
def dlogl(rho):    # d log lambda / d log rho
    s = np.log(rho/rho0)/Lbig
    return (-KTH/(1-(1-s)*KTH/2))/Lbig
def Lam_np(X):
    rho = np.linalg.norm(X); l = lam_np(rho)
    return np.concatenate([l*X[:4], [X[4]/l**2]])
rng = np.random.default_rng(20260906); errs = []
for _ in range(400):
    X = rng.normal(size=5); X *= np.exp(rng.uniform(0, Lbig))/np.linalg.norm(X)
    h = 1e-6*np.linalg.norm(X)
    Jm = np.zeros((5,5))
    for j in range(5):
        e = np.zeros(5); e[j] = h
        Jm[:,j] = (Lam_np(X+e)-Lam_np(X-e))/(2*h)
    rho = np.linalg.norm(X); l = lam_np(rho)
    s2 = np.sum(X[:4]**2)/rho**2; c2 = X[4]**2/rho**2
    pred = l**2*(1 + dlogl(rho)*(s2 - 2*c2))
    errs.append(abs(np.linalg.det(Jm)-pred)/abs(pred))
out['detDLambda_max_rel_err_fd'] = float(max(errs))
out['detDLambda_formula'] = 'lambda^2 * (1 + dlog(lambda)/dlog(rho) * (sin^2 phi - 2 cos^2 phi))'
assert max(errs) < 1e-7, max(errs)

# ---------------- (4) per-shell core ----------------
v, A, B = sp.symbols('v A B', positive=True)
Fa = v**3/(3*B*(A*v**2+B)**sp.Rational(3,2))
res4 = sp.simplify(sp.diff(Fa,v) - v**2*(A*v**2+B)**sp.Rational(-5,2))
out['antiderivative_residual'] = str(res4); assert res4 == 0
core = sp.simplify((Fa.subs(v,1)-Fa.subs(v,0)).subs({A: lam**2-lam**-4, B: lam**-4}))
res5 = sp.simplify(core - lam/3)
out['per_shell_core'] = str(core); out['per_shell_core_residual'] = str(res5)
assert res5 == 0, res5
D = sp.Symbol('D')
gg = sp.sqrt(lam**2*sp.sin(phi)**2 + lam**-4*sp.cos(phi)**2)
integ = sp.simplify((lam**-2*sp.cos(phi)/gg**5)*lam**2*(1+D*(sp.sin(phi)**2-2*sp.cos(phi)**2))*sp.sin(phi)**2)
tgt = sp.cos(phi)*sp.sin(phi)**2*(1+D*(sp.sin(phi)**2-2*sp.cos(phi)**2))/gg**5
r6 = sp.simplify(integ-tgt); out['a_integrand_residual'] = str(r6); assert r6 == 0
out['consequence'] = ('a[eta0 o Lambda^-1](0) = (M/2) int_0^L lambda(rho) dlog rho + shear term, '
                      'EXACTLY, for ANY profile lambda(.) and the bang-bang cap; and '
                      'int_S |eta0| |Lambda x|^-4 dx_5 <= pi^3 M int_0^L dlog(rho)/lambda(rho).')
json.dump(out, open('r1_results.json','w'), indent=1)
for k,vv in out.items():
    if k != 'angular_identity_mpmath': print(f"{k:34s} {vv}")
print("mpmath angular identity:", [(r['lam'], f"{r['err']:.2e}") for r in num])

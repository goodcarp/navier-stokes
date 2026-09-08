"""p2 - the three shell identities Lemma T' needs.

  (I1)  J(lambda) := INT_0^pi sin^2(phi) / g(phi;lambda)^4 dphi = pi/(2 lambda)      EXACT
  (I2)  INT_0^{pi/2} cos(phi) sin^2(phi) / g(phi;lambda)^5 dphi = lambda/3            EXACT
  (I3)  Q(lambda) := INT_0^{pi/2} cos sin^2 (sin^2 - 2 cos^2)/g^5 dphi                closed form

with g(phi;lambda)^2 = lambda^2 sin^2(phi) + lambda^-4 cos^2(phi) = |T_lambda x|^2/|x|^2.
(I1) is what makes the shell integral of Lemma T' work for a rho-DEPENDENT lambda: at each
fixed rho the phi-integral is done with lambda = lambda(rho) and returns pi/(2 lambda(rho)),
so the radial integral is just  INT dlog(rho)/lambda(rho).

Also: P_h(lambda) for the delta-tapered profile, and r_h := inf_{[1,3/2]} P_h/lambda.
"""
import json, sympy as sp, mpmath as mp

out = {}
phi, lam, A, B, v, aa, bb = sp.symbols('phi lambda A B v a b', positive=True)
g2 = lam**2*sp.sin(phi)**2 + lam**-4*sp.cos(phi)**2

# ---------------------------------------------------- (I1) via the exact antiderivative
# d/dphi [ atan(sqrt(a/b) tan phi)/sqrt(a b) ] = 1/(b cos^2 + a sin^2)
F = sp.atan(sp.sqrt(aa/bb)*sp.tan(phi))/sp.sqrt(aa*bb)
res_anti = sp.simplify(sp.diff(F, phi) - 1/(bb*sp.cos(phi)**2 + aa*sp.sin(phi)**2))
out['antiderivative_residual'] = str(res_anti); assert res_anti == 0
# hence  INT_0^pi dphi/(b cos^2 + a sin^2) = pi/sqrt(ab)   (two half-periods)
# differentiate in a:  INT_0^pi sin^2/(b cos^2 + a sin^2)^2 dphi = (pi/2) a^{-3/2} b^{-1/2}
dF = sp.simplify(-sp.diff(sp.pi/sp.sqrt(aa*bb), aa) - sp.pi/2*aa**sp.Rational(-3,2)*bb**sp.Rational(-1,2))
out['da_residual'] = str(sp.simplify(dF)); assert sp.simplify(dF) == 0
# a = lambda^2, b = lambda^-4  =>  (pi/2) lambda^-3 lambda^2 = pi/(2 lambda)
J_claim = sp.simplify((sp.pi/2*aa**sp.Rational(-3,2)*bb**sp.Rational(-1,2)).subs({aa: lam**2, bb: lam**-4}))
out['J_lambda_symbolic'] = str(J_claim)
assert sp.simplify(J_claim - sp.pi/(2*lam)) == 0

# independent 40-digit quadrature of (I1)
mp.mp.dps = 40
def Jnum(l):
    l = mp.mpf(l)
    f = lambda p: mp.sin(p)**2/(l**2*mp.sin(p)**2 + l**-4*mp.cos(p)**2)**2
    return mp.quad(f, [0, mp.pi/2, mp.pi])
lams = [0.5, 1.0, 1.1, 1.2247448713915890, 1.25, 1.5, 3.0]
out['J_check'] = [{'lam': l, 'quad': mp.nstr(Jnum(l), 25),
                   'exact': mp.nstr(mp.pi/(2*mp.mpf(l)), 25),
                   'rel': float(abs(Jnum(l) - mp.pi/(2*mp.mpf(l)))/(mp.pi/(2*mp.mpf(l))))}
                  for l in lams]
out['J_max_rel'] = max(d['rel'] for d in out['J_check'])
assert out['J_max_rel'] < 1e-30, out['J_max_rel']

# ---------------------------------------------------- (I2),(I3) via v = sin(phi)
# g^2 = A v^2 + B  with A = lambda^2 - lambda^-4, B = lambda^-4;  cos phi dphi = dv
Asub, Bsub = lam**2 - lam**-4, lam**-4
I2 = sp.integrate(v**2*(A*v**2 + B)**sp.Rational(-5,2), (v, 0, 1))
I2 = sp.simplify(I2.subs({A: Asub, B: Bsub}))
out['I2_symbolic'] = str(sp.simplify(I2))
assert sp.simplify(I2 - lam/3) == 0

I4 = sp.simplify(sp.integrate(v**4*(A*v**2 + B)**sp.Rational(-5,2), (v, 0, 1)))
Qsym = sp.simplify((3*I4 - 2*sp.integrate(v**2*(A*v**2+B)**sp.Rational(-5,2), (v, 0, 1)))
                   .subs({A: Asub, B: Bsub}))
Qsym = sp.simplify(sp.powsimp(sp.radsimp(Qsym), force=True))
out['Q_symbolic'] = str(Qsym)
Q_lim1 = sp.simplify(sp.limit(Qsym, lam, 1))
out['Q_at_lambda_1_exact'] = str(Q_lim1)
def Qf(l):
    if abs(float(l) - 1.0) < 1e-12:
        return mp.mpf(str(sp.N(Q_lim1, 40)))
    return mp.mpf(str(sp.N(Qsym.subs(lam, sp.Rational(l)), 40)))
def Qnum(l):
    l = mp.mpf(l)
    f = lambda p: mp.cos(p)*mp.sin(p)**2*(mp.sin(p)**2 - 2*mp.cos(p)**2) \
                  /(l**2*mp.sin(p)**2 + l**-4*mp.cos(p)**2)**mp.mpf(2.5)
    return mp.quad(f, [0, mp.pi/2])
out['Q_check'] = [{'lam': l, 'closed': mp.nstr(Qf(l), 20), 'quad': mp.nstr(Qnum(l), 20),
                   'absdiff': float(abs(Qf(l) - Qnum(l)))} for l in lams]
out['Q_max_absdiff'] = max(d['absdiff'] for d in out['Q_check'])
assert out['Q_max_absdiff'] < 1e-30, out['Q_max_absdiff']

# ---------------------------------------------------- P_h for the delta taper (numeric)
def P_h(l, deg):
    mp.mp.dps = 25
    l = mp.mpf(l); d = mp.mpf(deg)*mp.pi/180
    def f(p):
        h = min(mp.mpf(1), p/d)
        return 3*h*mp.cos(p)*mp.sin(p)**2/(l**2*mp.sin(p)**2 + l**-4*mp.cos(p)**2)**mp.mpf(2.5)
    return mp.quad(f, [0, min(d, mp.pi/2), mp.pi/2])
tab = {}
for deg in (3, 5, 7.5, 10, 15, 20, 30):
    row = {}
    ls = [1.0 + k*0.02 for k in range(26)]
    vals = [(l, float(P_h(l, deg))) for l in ls]
    row['P_at_1'] = vals[0][1]; row['P_at_1.5'] = vals[-1][1]
    row['inf_P_over_lam'] = min(pv/l for l, pv in vals)
    row['argmin_lam'] = min(vals, key=lambda t: t[1]/t[0])[0]
    row['min_P'] = min(pv for _, pv in vals)
    tab[str(deg)] = row
out['P_h_table'] = tab
out['kappa_delta_7p5'] = tab['7.5']['P_at_1']/2.0     # a(0)/(ML) at lambda=1 is P_h/2

json.dump(out, open('p2_results.json', 'w'), indent=1)
print("(I1)  INT_0^pi sin^2/g^4 dphi =", out['J_lambda_symbolic'],
      "   max rel err vs 40-dps quadrature =", out['J_max_rel'])
print("(I2)  INT_0^{pi/2} cos sin^2/g^5 dphi =", out['I2_symbolic'])
print("(I3)  Q(lambda) =", out['Q_symbolic'])
print("      max |closed - quad| =", out['Q_max_absdiff'])
print("kappa_delta(7.5 deg) = P_h(1)/2 =", repr(out['kappa_delta_7p5']),
      "   (prove-lagrangian: 0.4997212)")
print("delta  P_h(1)     P_h(3/2)   inf P_h/lambda on [1,3/2]")
for k, r in tab.items():
    print(f"{k:>5s}  {r['P_at_1']:.6f}  {r['P_at_1.5']:.6f}   {r['inf_P_over_lam']:.6f}")

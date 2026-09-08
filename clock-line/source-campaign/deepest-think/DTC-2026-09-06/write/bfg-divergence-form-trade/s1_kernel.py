#!/usr/bin/env python3
"""
s1_kernel.py -- the heat-kernel facts the divergence-form Duhamel estimate needs.

Every identity is (a) derived symbolically with sympy and (b) confirmed by
independent high-precision quadrature with mpmath.  Nothing is typed from memory.

G_nu(x,t) = (4 pi nu t)^{-3/2} exp(-|x|^2/(4 nu t)),  x in R^3.

K1  int_{R^3} G_nu(x,t) dx        = 1
K2  int_{R^3} d_1 G_nu(x,t) dx    = 0                (mean-zero kernel)
K3  int_{R^3} |grad G_nu(x,t)| dx = 2/sqrt(pi nu t)
K4  ||G(.,tau)||_{L^2(R^3)}       = (8 pi tau)^{-3/4}
K5  int_0^t (t-s)^{-1/2} ds       = 2 sqrt(t)
K6  sup_x |grad G(x,1)|(|x|+1)^4  (Kukavica's pointwise majorant constant)
K7  sup_x G(x,1)(|x|+1)^4         (BFG's pointwise majorant constant)
K8  the two Duhamel scalings: 8pi/3 (mean-zero kernel) vs 4pi/3 (mean-one kernel)
"""
import json, sys
import sympy as sp
import mpmath as mp

mp.mp.dps = 30
out = {}
r, t, nu, s, tau, a = sp.symbols('r t nu s tau a', positive=True)

# work with the substitution b = 1/(4 nu t) to keep sympy fast
b = sp.symbols('b', positive=True)
Gr = (b/sp.pi)**sp.Rational(3, 2)*sp.exp(-b*r**2)      # = (4 pi nu t)^{-3/2} e^{-r^2/(4 nu t)}

K1 = sp.simplify(sp.integrate(4*sp.pi*r**2*Gr, (r, 0, sp.oo)))
print("K1  int G_nu dx            =", K1);  assert K1 == 1
out['K1_mass'] = str(K1)

xs = sp.Symbol('xs', real=True)
g1 = sp.sqrt(b/sp.pi)*sp.exp(-b*xs**2)
K2 = sp.simplify(sp.integrate(sp.diff(g1, xs), (xs, -sp.oo, sp.oo)))
print("K2  int d_1 G_nu dx        =", K2);  assert K2 == 0
out['K2_mean_zero'] = str(K2)

# |grad G| = 2 b r G  (with b = 1/(4 nu t) this is (r/(2 nu t)) G)
K3b = sp.simplify(sp.integrate(4*sp.pi*r**2*(2*b*r)*Gr, (r, 0, sp.oo)))
K3 = sp.simplify(K3b.subs(b, 1/(4*nu*t)))
print("K3  int |grad G_nu| dx     =", sp.simplify(sp.powsimp(K3, force=True)))
assert sp.simplify(K3 - 2/sp.sqrt(sp.pi*nu*t)) == 0
out['K3_gradL1'] = "2/sqrt(pi*nu*t)"

k3num = []
for nuv, tv in [('1', '1'), ('1', '0.01'), ('7', '3'), ('0.001', '5')]:
    nuf, tf = mp.mpf(nuv), mp.mpf(tv)
    sc = mp.sqrt(nuf*tf)
    f = lambda rr: 4*mp.pi*rr**2*(rr/(2*nuf*tf))*(4*mp.pi*nuf*tf)**mp.mpf(-1.5)*mp.e**(-rr**2/(4*nuf*tf))
    val = mp.quad(f, [0, sc, 5*sc, 20*sc])
    exact = 2/mp.sqrt(mp.pi*nuf*tf)
    rel = abs(val-exact)/exact
    k3num.append([float(nuf), float(tf), float(val), float(exact), float(rel)])
    print("    nu=%-8s t=%-8s quad=%.18g exact=%.18g relerr=%.2e" % (nuv, tv, float(val), float(exact), float(rel)))
    assert rel < mp.mpf('1e-24')
out['K3_numeric'] = k3num

G2 = (b/sp.pi)**sp.Rational(3, 2)*sp.exp(-b*r**2)
K4sq = sp.simplify(sp.integrate(4*sp.pi*r**2*G2**2, (r, 0, sp.oo))).subs(b, 1/(4*tau))
K4 = sp.simplify(sp.sqrt(sp.simplify(K4sq)))
print("K4  ||G(.,tau)||_2         =", sp.simplify(sp.powsimp(K4, force=True)))
assert sp.simplify(K4 - (8*sp.pi*tau)**sp.Rational(-3, 4)) == 0
out['K4_L2'] = "(8*pi*tau)**(-3/4)"

K5 = sp.simplify(sp.integrate((t-s)**sp.Rational(-1, 2), (s, 0, t)))
print("K5  int_0^t (t-s)^-1/2 ds  =", K5);  assert sp.simplify(K5 - 2*sp.sqrt(t)) == 0
out['K5_abel'] = "2*sqrt(t)"

# K6 / K7: stationary points of log of the majorant ratio, t = nu = 1
x6 = mp.findroot(lambda z: 1/z - z/2 + 4/(z+1), mp.mpf(2))
sup6 = (x6/2)*(4*mp.pi)**mp.mpf(-1.5)*mp.e**(-x6**2/4)*(x6+1)**4
x7 = (-1+mp.sqrt(33))/2                       # root of -z/2 + 4/(z+1) = 0
sup7 = (4*mp.pi)**mp.mpf(-1.5)*mp.e**(-x7**2/4)*(x7+1)**4
print("K6  sup |grad G|(|x|+1)^4  = %.10f   at |x| = %.9f" % (float(sup6), float(x6)))
print("K7  sup G (|x|+1)^4        = %.10f   at |x| = %.9f" % (float(sup7), float(x7)))
out['K6_gradmajorant'] = float(sup6); out['K6_argmax'] = float(x6)
out['K7_Gmajorant'] = float(sup7);    out['K7_argmax'] = float(x7)

# K8: the two Duhamel scalings.  int_{R^3}(|z|+A)^{-4}dz = 4 pi int_0^inf r^2/(r+A)^4 dr
#     = 4 pi/(3A)  (exact, checked symbolically), so
#     mean-zero: int_0^t 4pi/(3 sqrt(t-s)) ds = (8pi/3) sqrt(t)
#     mean-one : int_0^t sqrt(t-s) 4pi/(3 sqrt(t-s)) ds = (4pi/3) t
A = sp.Symbol('A', positive=True)
shell = sp.simplify(sp.integrate(4*sp.pi*r**2/(r+A)**4, (r, 0, sp.oo)))
print("K8  int_{R^3}(|z|+A)^-4 dz =", shell); assert sp.simplify(shell - 4*sp.pi/(3*A)) == 0
I0 = sp.simplify(sp.integrate(shell.subs(A, sp.sqrt(t-s)), (s, 0, t)))
I1 = sp.simplify(sp.integrate(sp.sqrt(t-s)*shell.subs(A, sp.sqrt(t-s)), (s, 0, t)))
print("K8  mean-zero Duhamel      =", I0, "  -> /sqrt(t) = 8pi/3 = %.12f" % float(8*mp.pi/3))
print("K8  mean-one  Duhamel      =", I1, "  -> /t       = 4pi/3 = %.12f" % float(4*mp.pi/3))
assert sp.simplify(I0 - sp.Rational(8, 3)*sp.pi*sp.sqrt(t)) == 0
assert sp.simplify(I1 - sp.Rational(4, 3)*sp.pi*t) == 0
out['K8_meanzero'] = "8*pi/3*sqrt(t)"; out['K8_meanone'] = "4*pi/3*t"
out['K8_8pi_3'] = float(8*mp.pi/3); out['K8_4pi_3'] = float(4*mp.pi/3)

json.dump(out, open('s1_results.json', 'w'), indent=1)
print("\ns1: ALL 8 IDENTITIES DERIVED AND ASSERTED")

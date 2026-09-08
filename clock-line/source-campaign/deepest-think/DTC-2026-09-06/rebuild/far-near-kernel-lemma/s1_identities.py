#!/usr/bin/env python3
"""S1 - the two exact zonal-harmonic identities the far/near split rests on (sympy, exact).

Setting: R^5 = R^4_y x R_z, rho = |x|, t = z/rho.  Zonal harmonics of R^5 are rho^l C_l^{3/2}(t)
(interior) and rho^{-l-3} C_l^{3/2}(t) (exterior).  a = -d_z psi1, so we need d_z of both families.

I1  d_z [ rho^l      C_l^{3/2}(z/rho) ] =  (l+2)  rho^{l-1}   C_{l-1}^{3/2}(z/rho)     l >= 1
    [PREREG registered (2l+1); K2 FIRED on that form at l=2 (12z vs 15z).  The exact constant is
     (l+2) = l*C_l(1)/C_{l-1}(1); the pre-registered (2l+1) came from matching the coefficient of
     t^l in C_l instead of the coefficient of z^l in rho^l C_l(z/rho).  Both forms are tested below.]
I2  d_z [ rho^{-l-3} C_l^{3/2}(z/rho) ] = -(l+1)  rho^{-l-4}  C_{l+1}^{3/2}(z/rho)     l >= 0
I3  N_l := INT_{-1}^{1} C_l^{3/2}(t)^2 (1-t^2) dt = (l+1)(l+2)/(l+3/2)
I4  matching for a unit shell source: alpha_l = g_l/(2l+3)  (interior coefficient)
I5  h_l (Gegenbauer coefficients of the bang-bang source -M sgn(t)/sqrt(1-t^2)), h_1 = -5M/6,
    hence Phi(0) = -(3/5) h_1 = M/2 = kappa_0.
KILLS: K2 (I1), K3 (I2).
"""
import sympy as sp, hashlib, sys
s, z = sp.symbols('s z', positive=True)          # s = |y| (transverse radius), z axial
rho = sp.sqrt(s**2 + z**2); t = z/rho
LMAX = 14
fails = []

print("=== I1: d_z [ rho^l C_l^{3/2}(t) ] = c_l rho^{l-1} C_{l-1}^{3/2}(t):  PREREG c_l=(2l+1) vs REPAIRED c_l=(l+2) ===")
prereg_fired = False
for l in range(1, LMAX+1):
    lhs = sp.diff(rho**l * sp.gegenbauer(l, sp.Rational(3,2), t), z)
    ok_pre = sp.simplify(sp.expand(sp.simplify(lhs - (2*l+1)*rho**(l-1)*sp.gegenbauer(l-1, sp.Rational(3,2), t)))) == 0
    ok_rep = sp.simplify(sp.expand(sp.simplify(lhs - (l+2)  *rho**(l-1)*sp.gegenbauer(l-1, sp.Rational(3,2), t)))) == 0
    print(f"  l={l:2d}   (2l+1): {'OK  ' if ok_pre else 'FAIL'}    (l+2): {'OK' if ok_rep else 'FAIL'}")
    if not ok_pre: prereg_fired = True
    if not ok_rep: fails.append(('I1-repaired', l))
print("  KILL K2 fired on the pre-registered constant (2l+1):", prereg_fired)
print("  repaired constant (l+2) = l*C_l(1)/C_{l-1}(1) holds exactly for every l <= %d" % LMAX)

print("=== I2: d_z [ rho^{-l-3} C_l^{3/2}(t) ] = -(l+1) rho^{-l-4} C_{l+1}^{3/2}(t) ===")
for l in range(0, LMAX+1):
    lhs = sp.diff(rho**(-l-3) * sp.gegenbauer(l, sp.Rational(3,2), t), z)
    rhs = -(l+1) * rho**(-l-4) * sp.gegenbauer(l+1, sp.Rational(3,2), t)
    ok = sp.simplify(sp.expand(sp.simplify(sp.together(lhs - rhs)))) == 0
    print(f"  l={l:2d}  {'OK' if ok else 'FAIL'}")
    if not ok: fails.append(('I2', l))

print("=== I3: N_l closed form ===")
u = sp.symbols('u', real=True)
def N_cf(l): return sp.Rational((l+1)*(l+2),1)/sp.Rational(2*l+3,2)
for l in range(0, 9):
    ex = sp.integrate(sp.gegenbauer(l, sp.Rational(3,2), u)**2*(1-u**2), (u,-1,1))
    ok = sp.simplify(ex - N_cf(l)) == 0
    print(f"  l={l}  exact {ex}  closed form {N_cf(l)}  {'OK' if ok else 'FAIL'}")
    if not ok: fails.append(('I3', l))

print("=== I4: shell matching alpha_l = g_l/(2l+3) ===")
# A_l = alpha rho^l (rho<1), beta rho^{-l-3} (rho>1); continuity + jump [A'] = -g_l
al, be, g = sp.symbols('al be g')
sol = sp.solve([sp.Eq(al, be), sp.Eq(-(l_:=0)*0 + (-(sp.Symbol('L')+3)*be - sp.Symbol('L')*al), -g)], [al, be], dict=True)
L = sp.Symbol('L')
sol = sp.solve([sp.Eq(al, be), sp.Eq(-(L+3)*be - L*al, -g)], [al, be], dict=True)[0]
print("  alpha_l =", sp.simplify(sol[al]), "  (want g/(2L+3))")
ok = sp.simplify(sol[al] - g/(2*L+3)) == 0
print("  ", "OK" if ok else "FAIL")
if not ok: fails.append(('I4', 0))

print("=== I5: bang-bang source coefficients h_l, and Phi(0) = M/2 ===")
def mom(k): return sp.Rational(1,2)*sp.gamma(sp.Rational(k+1,2))*sp.gamma(sp.Rational(3,2))/sp.gamma(sp.Rational(k+1,2)+sp.Rational(3,2))
def h_coef(l):
    C = sp.Poly(sp.expand(sp.gegenbauer(l, sp.Rational(3,2), u)), u)
    return sp.nsimplify(sp.simplify(-2*sum(co*mom(k) for (k,), co in C.terms())/N_cf(l)), rational=True)
H = {l: h_coef(l) for l in [1,3,5,7,9,11,13]}
print("  h_l =", H)
ok = (H[1] == sp.Rational(-5,6))
print(f"  h_1 = {H[1]} (want -5/6) -> Phi(0) = -3/5*h_1 = {sp.Rational(-3,5)*H[1]} (want 1/2)  {'OK' if ok else 'FAIL'}")
if not ok: fails.append(('I5', 1))
# sharpness of kappa_0 = 1/2: Phi(0) = -(3/4) INT_0^pi w cos sin^2 dphi, |.| <= (3/4) M INT |cos| sin^2 = M/2
ph = sp.symbols('ph')
print("  (3/4) INT_0^pi |cos phi| sin^2 phi dphi =", sp.Rational(3,4)*sp.integrate(sp.Abs(sp.cos(ph))*sp.sin(ph)**2,(ph,0,sp.pi)), " (sharp bound on |Phi(0)|/M)")

print("\nFAILURES:", fails if fails else "none")
print('SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(1 if fails else 0)

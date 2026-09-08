#!/usr/bin/env python3
"""
r4_recheck.py -- second refuter pass over write/bfg-divergence-form-trade.

Reads NOTHING from the target folder.  Everything below is derived here from the
definitions.  Purpose: (i) reconfirm the load-bearing constants by my own route,
(ii) pin the two claims I believe are defective (the "dimensionally impossible"
headline, and the vacuous entries inside the folder's own audit), (iii) count how
many of s7_audit.py's "9 closed forms" are capable of failing.
"""
import json
import sympy as sp
import mpmath as mp
mp.mp.dps = 50
OUT = {}

def line(t): print("\n" + t + "\n" + "-"*len(t))

# ============================================================ A. kernel facts
line("A. kernel facts, my own route")
r, nu, t, tau, s = sp.symbols('r nu t tau s', positive=True)
G = (4*sp.pi*nu*t)**sp.Rational(-3,2)*sp.exp(-r**2/(4*nu*t))
mass = sp.simplify(sp.integrate(4*sp.pi*r**2*G, (r, 0, sp.oo)))
gradL1 = sp.simplify(sp.integrate(4*sp.pi*r**2*(r/(2*nu*t))*G, (r, 0, sp.oo)))
print("  int G            =", mass)
print("  int |grad G|     =", sp.simplify(sp.powsimp(gradL1, force=True)))
assert mass == 1
assert sp.simplify(gradL1 - 2/sp.sqrt(sp.pi*nu*t)) == 0
# probabilistic cross-route: E|X|/(2 nu t), X ~ N(0, 2 nu t I_3), E|X| = 2 sigma sqrt(2/pi)
sig = sp.sqrt(2*nu*t)
prob = sp.simplify(2*sig*sp.sqrt(2/sp.pi)/(2*nu*t))
print("  E|X|/(2 nu t)    =", sp.simplify(sp.powsimp(prob, force=True)), " residual",
      sp.simplify(prob - gradL1))
assert sp.simplify(prob - gradL1) == 0
G2 = (4*sp.pi*tau)**sp.Rational(-3,2)*sp.exp(-r**2/(4*tau))
L2 = sp.sqrt(sp.simplify(sp.integrate(4*sp.pi*r**2*G2**2, (r, 0, sp.oo))))
print("  ||G(.,tau)||_2   =", sp.simplify(sp.powsimp(L2, force=True)))
assert sp.simplify(L2 - (8*sp.pi*tau)**sp.Rational(-3,4)) == 0
OUT['gradL1'] = "2/sqrt(pi*nu*t)"; OUT['L2'] = "(8*pi*tau)**(-3/4)"

# ============================================================ B. Q = 1
line("B. Q = sup |(n.u)w - (n.w)u|/(|u||w|), my own route")
c, al, be = sp.symbols('c alpha beta', real=True)
Gram = sp.Matrix([[1,c,al],[c,1,be],[al,be,1]])
det = sp.expand(Gram.det())
obj = al**2 + be**2 - 2*c*al*be
print("  det Gram         =", det)
print("  det + obj - (1-c^2) =", sp.expand(det + obj - (1-c**2)))
assert sp.expand(det + obj - (1-c**2)) == 0
# so obj <= 1-c^2 <= 1 on the feasible (PSD) set.  Independent numeric maximisation:
import numpy as np
rng = np.random.default_rng(4242)
N = 3000000
u = rng.normal(size=(N,3)); w = rng.normal(size=(N,3)); n = rng.normal(size=(N,3))
n /= np.linalg.norm(n,axis=1)[:,None]
val = np.linalg.norm((n*u).sum(1)[:,None]*w - (n*w).sum(1)[:,None]*u, axis=1)/(
      np.linalg.norm(u,axis=1)*np.linalg.norm(w,axis=1))
print("  numeric max over %d samples = %.15f" % (N, val.max()))
assert val.max() <= 1 + 1e-12
OUT['Q'] = 1.0; OUT['Q_numeric_max'] = float(val.max())

# ============================================================ C. Lemma U
line("C. Lemma U constant, my own optimisation")
A_, B_, M, E = sp.symbols('A B M E', positive=True)
f = A_*sp.sqrt(tau) + B_*tau**sp.Rational(-3,4)
ts = sp.solve(sp.diff(f, tau), tau)[0]
fmin = sp.simplify(sp.powsimp(f.subs(tau, ts), force=True))
gam = sp.simplify(sp.powsimp(fmin/(A_**sp.Rational(3,5)*B_**sp.Rational(2,5)), force=True))
print("  tau_*  =", ts)
print("  gamma  =", gam, "=", sp.nsimplify(gam), "=", float(sp.N(gam,30)))
Cu = sp.simplify(sp.powsimp(fmin.subs({A_: 4*M/sp.sqrt(sp.pi),
        B_: (8*sp.pi)**sp.Rational(-3,4)*sp.sqrt(E)})/(E**sp.Rational(1,5)*M**sp.Rational(3,5)),
        force=True))
print("  C_u    =", sp.simplify(Cu), "=", float(sp.N(Cu,30)))
# my own independent bisection on f'(tau)=0 at M=E=1, no closed form used
g  = lambda v: (4/mp.sqrt(mp.pi))*mp.sqrt(v) + (8*mp.pi)**mp.mpf('-0.75')*v**mp.mpf('-0.75')
dg = lambda v: (2/mp.sqrt(mp.pi))/mp.sqrt(v) - mp.mpf('0.75')*(8*mp.pi)**mp.mpf('-0.75')*v**mp.mpf('-1.75')
lo, hi = mp.mpf('1e-6'), mp.mpf('10')
for _ in range(300):
    mid = (lo+hi)/2
    if dg(mid) > 0: hi = mid
    else: lo = mid
tstar_num = (lo+hi)/2
print("  bisection: tau_* = %.18f   f(tau_*) = %.18f" % (float(tstar_num), float(g(tstar_num))))
assert abs(mp.mpf(str(sp.N(Cu,40))) - g(tstar_num)) < mp.mpf('1e-25')
OUT['C_u'] = float(sp.N(Cu,30)); OUT['tau_star'] = float(tstar_num)
OUT['gamma'] = float(sp.N(gam,30))

# ============================================================ D. Theorem C
line("D. Theorem C / c_div, my own bootstrap algebra")
tt, Re = sp.symbols('t Re', positive=True)
Cnu = 2*1/sp.sqrt(sp.pi*nu)              # Q = 1
Duh = Cnu*(2*sp.sqrt(tt))*(2*M)*(Cu*E**sp.Rational(1,5)*(2*M)**sp.Rational(3,5))
tstar = sp.simplify(sp.powsimp(sp.solve(sp.Eq(Duh, M/2), tt)[0], force=True))
print("  Duhamel(t) <=", sp.simplify(sp.powsimp(Duh, force=True)))
print("  t_*        =", tstar)
nu_of_Re = E**sp.Rational(2,5)*M**sp.Rational(1,5)/Re
cdiv = sp.simplify(sp.powsimp(tstar.subs(nu, nu_of_Re)*M*Re, force=True))
print("  c_div = M Re t_* =", sp.nsimplify(cdiv), " = %.18e" % float(sp.N(cdiv,40)))
assert sp.simplify(cdiv - 3*3**sp.Rational(1,5)*sp.pi**sp.Rational(11,5)/12800) == 0
assert sp.simplify(cdiv - sp.pi/(2**sp.Rational(46,5)*Cu**2)) == 0
assert sp.simplify(sp.diff(cdiv, M)) == 0 and sp.simplify(sp.diff(cdiv, E)) == 0
OUT['c_div'] = float(sp.N(cdiv,40))
# fully numeric end-to-end at an independent test point of my own choosing
Mv, Ev, nuv = mp.mpf('0.83'), mp.mpf('11.9'), mp.mpf('0.207')
Cun = mp.mpf(str(sp.N(Cu,40)))
coef = (2/mp.sqrt(mp.pi*nuv))*2*(2*Mv)*Cun*Ev**(mp.mpf(1)/5)*(2*Mv)**(mp.mpf(3)/5)
tsn  = (Mv/(2*coef))**2
ReEn = Ev**(mp.mpf(2)/5)*Mv**(mp.mpf(1)/5)/nuv
print("  independent test point M=%s E=%s nu=%s -> Re_E=%.12f" % (Mv,Ev,nuv,float(ReEn)))
print("  t_*=%.12e  Duhamel(t_*)=%.15f  M/2=%.15f  M Re t_*=%.18e"
      % (float(tsn), float(coef*mp.sqrt(tsn)), float(Mv/2), float(Mv*ReEn*tsn)))
assert abs(coef*mp.sqrt(tsn) - Mv/2) < mp.mpf('1e-40')
assert abs(Mv*ReEn*tsn - mp.mpf(str(sp.N(cdiv,40)))) < mp.mpf('1e-30')
OUT['c_div_endtoend'] = float(Mv*ReEn*tsn)

# ============================================================ E. F1
line("E. F1 -- is Re_E^{-1/2} 'dimensionally impossible'?")
expr = sp.sqrt(nu)/(E**sp.Rational(1,5)*M**sp.Rational(1,10))
ReE  = E**sp.Rational(2,5)*M**sp.Rational(1,5)/nu
print("  sqrt(nu)/(E^{1/5}M^{1/10}) - Re_E^{-1/2} =",
      sp.simplify(sp.powsimp(expr - ReE**sp.Rational(-1,2), force=True)))
assert sp.simplify(sp.powsimp(expr - ReE**sp.Rational(-1,2), force=True)) == 0
L, T = sp.symbols('L T', positive=True)
dim = {E: L**5/T**2, M: 1/T, nu: L**2/T}
for nm, ex in [("Re_E", ReE), ("Re_E^{-1/2}", ReE**sp.Rational(-1,2)),
               ("sqrt(nu)/(E^{1/5}M^{1/10})", expr), ("M*t_*", M*tstar)]:
    d = sp.simplify(sp.powsimp(ex.subs(dim), force=True))
    pure = not ({L, T} & d.free_symbols)
    print("  [%-28s] = %s   -> dimensionless: %s" % (nm, d, pure))
    assert pure
print("  => every one of them is DIMENSIONLESS.  'M t_*/Re_E^{-1/2} is not a pure number'")
print("     is FALSE: that ratio IS Re_E^{-1/2}, a dimensionless group.")
OUT['F1_expr_equals_ReE_minus_half'] = True
OUT['F1_all_dimensionless'] = True
# the numeric mismatch inside the same subsection
print("  M t_* (Theorem C) = c_div/Re_E,   but sec.6.2 displays M t_* = Re_E^{-1} 'identically'.")
print("  ratio of the two displays = 1/c_div = %.10f" % float(1/sp.N(cdiv,30)))
OUT['F1_display_gap'] = float(1/sp.N(cdiv,30))
# and Re^{-1/2} is a STRONGER lower bound than Re^{-1} for Re>1
for R in [1, 1e2, 1e6, 1e17]:
    Rm = mp.mpf(R); cd = mp.mpf(str(sp.N(cdiv,30)))
    print("  Re_E=%-8.0e  c_div/Re = %-14.6e  c_div/sqrt(Re) = %-14.6e  stronger by %.4e"
          % (float(Rm), float(cd/Rm), float(cd/mp.sqrt(Rm)), float(mp.sqrt(Rm))))

# ============================================================ F. F2 / F6
line("F. controls in the target's own audit that cannot fail")
# (i) s7's 'third way' types the answer in, and its input differs from Theorem C's t_*
t_s7 = nu*M**sp.Rational(-6,5)*E**sp.Rational(-2,5)
print("  s7 'third way' input t = nu M^{-6/5}E^{-2/5};  t_*(Theorem C)/t = %s = %.12e"
      % (sp.nsimplify(sp.simplify(sp.powsimp(tstar/t_s7, force=True))),
         float(sp.N(sp.simplify(tstar/t_s7), 30))))
assert sp.simplify(sp.powsimp(tstar/t_s7 - cdiv, force=True)) == 0
print("  -> the two differ by exactly c_div; the s7 section can never see Theorem C's constant.")
OUT['F2_ratio_is_c_div'] = True
# (ii) s7 check #9 compares an expression to ITSELF
lhs9 = 3*2**sp.Rational(1,6)/(2*sp.pi**sp.Rational(2,3))
rhs9 = 3*2**sp.Rational(1,6)/(2*sp.pi**sp.Rational(2,3))
print("  s7 check 9 (C_w): lhs is rhs ->", sp.srepr(lhs9) == sp.srepr(rhs9),
      " (identical sympy trees; the check cannot fail for ANY value)")
OUT['F6_s7_check9_tautology'] = bool(sp.srepr(lhs9) == sp.srepr(rhs9))
# (iii) how many of the 9 can fail?
print("  s7's 9 checks, classified:")
cls = [("gamma reduced form","substantive"),("C_u vs gamma-route","substantive"),
       ("Duhamel prefactor","substantive"),("c_div two forms","substantive"),
       ("1/C_u^2 reduced","substantive"),("pi/2^(46/5) vs pi*2^(-46/5)","notational identity"),
       ("2*2^(18/5)=2^(23/5)","trivial exponent arithmetic"),
       ("(2^(23/5))^2=2^(46/5)","trivial exponent arithmetic"),
       ("C_w = C_w","TAUTOLOGY, cannot fail")]
for a,b in cls: print("    %-32s %s" % (a,b))
OUT['F6_substantive_of_9'] = 5
# (iv) the C_u cancellation
X, cG = sp.symbols('X c_G', positive=True)
ratio = (sp.pi/(2**sp.Rational(46,5)*X**2))/(cG/X**2)
print("  with ANY symbol X in place of C_u:  c_div(X)/c_GIM(X) =", sp.simplify(ratio),
      " d/dX =", sp.simplify(sp.diff(ratio, X)))
assert sp.simplify(sp.diff(ratio, X)) == 0
OUT['F3b_cancellation_structural'] = True
# (v) the cross-seat '4 ln 2'
print("  cross-seat 'check': (2 ln2/L)(1+2L-0.7031660) for L =")
for Lv in [8.3178, 20, 50, 100, 1e6]:
    Lm = mp.mpf(Lv)
    print("    L=%-10g  %.10f" % (Lv, float((2*mp.log(2)/Lm)*(1+2*Lm-mp.mpf('0.7031660')))))
print("    limit = 4 ln 2 = %.10f  -- forced by the two QUOTED formulas alone; c_div absent."
      % float(4*mp.log(2)))
print("    offset check: (2/5)*ln(0.172403978) = %.10f  (the quoted -0.7031660)"
      % float(mp.mpf(2)/5*mp.log(mp.mpf('0.172403978'))))
OUT['F3a_offset_from_Eshell'] = float(mp.mpf(2)/5*mp.log(mp.mpf('0.172403978')))

# ============================================================ G. the log-clock comparison
line("G. clock comparison reproduced by my own arithmetic")
cd = mp.mpf(str(sp.N(cdiv,30)))
for R in [mp.mpf(10)**k for k in [0,2,5,10,20,40]]:
    print("  Re=%-9.0e  (a)/c_1=%-14.6e  (b)=%-14.6e  (a)/(b)=%.6e * c_1"
          % (float(R), float(1/(1+mp.log(R))), float(cd/R), float((1/(1+mp.log(R)))/(cd/R))))
print("  sup_{Re>=1} (1+log Re)/Re = 1 at Re=1  ->  (a)>=(b) for all Re>=1 iff c_1 >= c_div")
for c1v in ['1e-3','1e-6']:
    tgt = cd/mp.mpf(c1v)
    R = mp.findroot(lambda RR: RR/(1+mp.log(RR)) - tgt, mp.mpf(2))
    print("  crossover at c_1=%s : Re_E = %.6g" % (c1v, float(R)))
print("  datum table:")
for Lv in [8.3178, 20, 50, 100]:
    Lm = mp.mpf(str(Lv)); lR = 2*Lm - mp.mpf('0.7031660'); ReEv = mp.e**lR
    print("    L=%-9.4f logRe=%-12.6f Re=%-14.6e (b)=%-14.6e (a)/c_1=%-14.6e model=%-10.6f model/(b)=%.6e"
          % (Lv, float(lR), float(ReEv), float(cd/ReEv), float(1/(1+lR)),
             float(2*mp.log(2)/Lm), float((2*mp.log(2)/Lm)/(cd/ReEv))))

json.dump(OUT, open('r4_results.json','w'), indent=1)
print("\nr4: recheck complete")

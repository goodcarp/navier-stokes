#!/usr/bin/env python3
"""
s4_clock.py -- assemble the divergence-form clock, extract its Reynolds exponent
exactly, and compare it with the log clock, the GIM/Kukavica velocity clock and
BFG's claimed Theorem 8/10.

Inputs (all produced by s1..s3 in this folder, none typed from memory):
   int |grad G_nu(.,t)| dx = 2/sqrt(pi nu t)                       [s1 K3]
   int_0^t (t-s)^{-1/2} ds = 2 sqrt(t)                             [s1 K5]
   Q = 1                                                            [s2]
   ||u||_inf <= C_u E^{1/5} M^{3/5},  C_u = 5*2^{9/10}3^{2/5}/(6 pi^{3/5})  [s3]
   Re_E = E^{2/5} M^{1/5} / nu     (estate convention, E = ||u||_2^2, no 1/2)
"""
import json
import sympy as sp
import mpmath as mp
mp.mp.dps = 40

s1 = json.load(open('s1_results.json'))
s2 = json.load(open('s2_results.json'))
s3 = json.load(open('s3_results.json'))
out = {}

M, E, nu, t, Re, Cu, Q, kap = sp.symbols('M E nu t Re C_u Q kappa', positive=True)
Cu_exact = 5*2**sp.Rational(9,10)*3**sp.Rational(2,5)/(6*sp.pi**sp.Rational(3,5))
assert abs(float(sp.N(Cu_exact,30)) - s3['C_u']) < 1e-14
Qv = sp.Integer(1)
assert s2['Q_exact'] == 1.0

# ---------- the a-priori bound ----------------------------------------------
# ||omega(t)||_inf <= M + C(nu) int_0^t (t-s)^{-1/2}||omega||_inf ||u||_inf ds,
# C(nu) = 2 Q / sqrt(pi nu)
Cnu = 2*Qv/sp.sqrt(sp.pi*nu)
print("C(nu) = 2Q/sqrt(pi nu) =", Cnu)
out['C_nu'] = "2/sqrt(pi*nu)"

# ---------- the bootstrap ----------------------------------------------------
# on [0,T*] assume ||omega||_inf <= K*M with K = 2; then ||u||_inf <= C_u E^{1/5}(K M)^{3/5}
K = sp.Integer(2)
Duh = Cnu*(2*sp.sqrt(t))*(K*M)*(Cu_exact*E**sp.Rational(1,5)*(K*M)**sp.Rational(3,5))
Duh = sp.simplify(sp.powsimp(Duh, force=True))
print("Duhamel(t) <=", Duh)
out['Duhamel_expr'] = str(Duh)

tsol = sp.solve(sp.Eq(Duh, M/2), t)
tstar = sp.simplify(sp.powsimp(tsol[0], force=True))
print("t_* (Duhamel = M/2) =", tstar)
out['t_star'] = str(tstar)

# ---------- exponent extraction (the honest way: solve for the exponents) ----
# claim t_* = c * M^alpha * Re^beta with Re = E^{2/5}M^{1/5}/nu.
al, be = sp.symbols('alpha beta', real=True)
# exponents of t_* in (M, E, nu):
eM = sp.degree(sp.Poly(sp.simplify(sp.log(tstar).rewrite(sp.log).expand(force=True)), ), gen=None) if False else None
# do it robustly: read the exponents off by differentiating log
LT = sp.expand_log(sp.log(tstar), force=True)
expM = sp.simplify(M*sp.diff(LT, M))
expE = sp.simplify(E*sp.diff(LT, E))
expN = sp.simplify(nu*sp.diff(LT, nu))
print("exponents of t_*:  M^(%s)  E^(%s)  nu^(%s)" % (expM, expE, expN))
out['exponents_t_star'] = {'M': str(expM), 'E': str(expE), 'nu': str(expN)}
sols = sp.solve([sp.Eq(expN, -be),
                 sp.Eq(expE, sp.Rational(2,5)*be),
                 sp.Eq(expM, al + sp.Rational(1,5)*be)], [al, be], dict=True)
print("solve overdetermined system (3 eqs, 2 unknowns) ->", sols)
assert len(sols) == 1 and sols[0][al] == -1 and sols[0][be] == -1
out['alpha'] = -1; out['beta'] = -1

# ---------- rewrite in Re_E --------------------------------------------------
nu_of_Re = E**sp.Rational(2,5)*M**sp.Rational(1,5)/Re
tRe = sp.simplify(sp.powsimp(tstar.subs(nu, nu_of_Re), force=True))
print("t_* in terms of (M, Re_E) =", tRe)
c_div = sp.simplify(sp.powsimp(tRe*M*Re, force=True))
print("c_div = M Re_E t_*        =", c_div, " = %.12e" % float(sp.N(c_div, 40)))
assert sp.simplify(tRe - c_div/(M*Re)) == 0
assert sp.simplify(sp.diff(c_div, E)) == 0 and sp.simplify(sp.diff(c_div, M)) == 0
c_div_closed = sp.pi/(2**sp.Rational(46,5)*Cu_exact**2)
print("closed form pi/(2^(46/5) C_u^2) =", sp.simplify(c_div_closed), " = %.12e" % float(sp.N(c_div_closed, 40)))
assert sp.simplify(c_div - c_div_closed) == 0
out['c_div'] = float(sp.N(c_div, 40))
out['c_div_closed_form'] = "pi/(2^(46/5) C_u^2) = 3 pi^(11/5) 2^(-73/10) 3^(-4/5) / 25"

# a fully reduced closed form
c_div_red = sp.nsimplify(sp.simplify(sp.powsimp(sp.radsimp(c_div_closed), force=True)))
print("c_div reduced             =", c_div_red)
out['c_div_reduced'] = str(c_div_red)

# ---------- the same computation without the antisymmetry (Q = 2) -----------
Duh2 = Cnu*2*(2*sp.sqrt(t))*(K*M)*(Cu_exact*E**sp.Rational(1,5)*(K*M)**sp.Rational(3,5))
t2 = sp.simplify(sp.solve(sp.Eq(Duh2, M/2), t)[0].subs(nu, nu_of_Re)*M*Re)
print("c_div if Q=2 (triangle ineq) =", sp.simplify(t2), " = %.12e  (ratio %s)" % (float(sp.N(t2,30)), sp.simplify(c_div/t2)))
out['c_div_Q2'] = float(sp.N(t2, 30)); out['Q_gain_factor'] = float(sp.N(sp.simplify(c_div/t2), 30))

# ---------- GIM / Kukavica velocity clock -----------------------------------
# T >= c_G nu / ||u_0||_inf^2 ; feed the same Lemma U bound
cG = sp.Symbol('c_G', positive=True)
T_gim = cG*nu/(Cu_exact*E**sp.Rational(1,5)*M**sp.Rational(3,5))**2
T_gim_Re = sp.simplify(sp.powsimp(T_gim.subs(nu, nu_of_Re), force=True))
print("T_GIM in (M,Re_E)         =", T_gim_Re)
cgim = sp.simplify(T_gim_Re*M*Re)
print("c_GIM = M Re_E T_GIM      =", cgim)
ratio = sp.simplify(c_div/cgim)
print("c_div / c_GIM             =", ratio, "  ( = pi 2^(-46/5)/c_G ; C_u CANCELS ) = %.12e / c_G" % float(sp.N(sp.pi/2**sp.Rational(46,5),30)))
assert sp.simplify(ratio - sp.pi/(2**sp.Rational(46,5)*cG)) == 0
out['c_div_over_c_GIM'] = "pi*2^(-46/5)/c_G = %.12e / c_G" % float(sp.N(sp.pi/2**sp.Rational(46,5),30))
out['pi_2_46_5'] = float(sp.N(sp.pi/2**sp.Rational(46,5), 30))

# ---------- the three clocks side by side ------------------------------------
c1 = sp.Symbol('c_1', positive=True)
print()
print("CLOCK COMPARISON (dimensionless window  M * T):")
print("  (a) log clock   (S1, Astra)        M T >= c_1/(1 + log_+ Re_E)")
print("  (b) divergence  (this note)        M T >= %.12e / Re_E" % float(sp.N(c_div,30)))
print("  (c) GIM/Kukavica velocity clock    M T >= (c_G/1) / Re_E        [same exponent]")
print("  (d) BFG Thm 8/10 as claimed        M T >= 1/c                   [no Re at all]")
rows = []
for R in [mp.mpf(10)**k for k in [0,2,5,10,20,40]]:
    a = 1/(1+max(mp.log(R), mp.mpf(0)))
    b = mp.mpf(str(sp.N(c_div,30)))/R
    rows.append([float(R), float(a), float(b), float(a/b)])
    print("   Re_E=%-9.0e   (a)/c_1 = %-14.6e   (b) = %-14.6e   ratio (a)/(b) = %.6e /c_1" % tuple(rows[-1]))
out['comparison_rows'] = rows

# crossover: c_1/(1+log Re) = c_div/Re  ->  Re/(1+log Re) = c_div/c_1
print()
print("crossover  Re/(1+log Re) = c_div/c_1 :")
cross = []
for c1v in ['1', '0.1', '0.01', '1e-3', '1e-6']:
    tgt = mp.mpf(str(sp.N(c_div,30)))/mp.mpf(c1v)
    if tgt < mp.mpf('1'):
        # Re/(1+log Re) >= 1 for Re>=1 with equality at Re=1 -> log clock wins for all Re>=1
        cross.append([float(mp.mpf(c1v)), None]); print("   c_1=%-8s -> target %.4e < 1: the log clock is larger for EVERY Re_E >= 1" % (c1v, float(tgt)))
    else:
        R = mp.findroot(lambda RR: RR/(1+mp.log(RR)) - tgt, mp.mpf(2))
        cross.append([float(mp.mpf(c1v)), float(R)]); print("   c_1=%-8s -> crossover at Re_E = %.6g" % (c1v, float(R)))
out['crossover'] = cross

json.dump(out, open('s4_results.json','w'), indent=1)
print("\ns4: clock assembled; Reynolds exponent = -1 EXACTLY (overdetermined solve consistent)")

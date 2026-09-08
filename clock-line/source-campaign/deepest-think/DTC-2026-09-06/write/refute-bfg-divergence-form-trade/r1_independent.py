#!/usr/bin/env python3
"""r1_independent.py -- REFUTER's own route to every load-bearing constant of
write/bfg-divergence-form-trade.  Nothing is read from that folder; every number
is recomputed here from a definition, by a route chosen to differ from theirs.
"""
import json
import numpy as np
import sympy as sp
import mpmath as mp
mp.mp.dps = 50
R = {}
def rec(k, v, note=""):
    R[k] = float(v); print("  %-46s %.18g   %s" % (k, float(v), note)); return v

print("=== A. kernel L^1 of |grad G_nu| : three routes ===")
# route 1 (mine): direct 3-D CARTESIAN quadrature, no radial substitution
def gradL1_cart(nu, t, Ncut=9.0, n=241):
    nu, t = float(nu), float(t)
    a = 4*nu*t
    L = Ncut*np.sqrt(a)
    x = np.linspace(-L, L, n); h = x[1]-x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r2 = X*X+Y*Y+Z*Z
    G = (np.pi*a)**-1.5*np.exp(-r2/a)
    gmag = np.sqrt(r2)/(2*nu*t)*G
    return gmag.sum()*h**3
# route 2 (mine): 1-D radial mpmath quad on the exact radial density
def gradL1_rad(nu, t):
    nu, t = mp.mpf(nu), mp.mpf(t)
    a = 4*nu*t
    f = lambda r: 4*mp.pi*r**2 * (r/(2*nu*t)) * (mp.pi*a)**mp.mpf(-1.5)*mp.e**(-r**2/a)
    return mp.quad(f, [0, mp.sqrt(a), 6*mp.sqrt(a), 40*mp.sqrt(a)])
closed = lambda nu,t: 2/mp.sqrt(mp.pi*mp.mpf(nu)*mp.mpf(t))
for (nu,t) in [(1,1),(1,0.01),(7,3),(0.001,5)]:
    c = closed(nu,t); q = gradL1_rad(nu,t); cart = gradL1_cart(nu,t)
    print("   nu=%-7g t=%-6g closed=%.18g  radial=%.18g (rel %.1e)  cartesian=%.6g (rel %.1e)"
          % (nu,t,float(c),float(q),float(abs(q-c)/c), cart, abs(cart-float(c))/float(c)))
    assert abs(q-c)/c < mp.mpf('1e-30')
    assert abs(cart-float(c))/float(c) < 2e-5
rec('gradL1_nu1_t1', closed(1,1))
# int G = 1 and int d1 G = 0
print("   int G_nu(nu=1,t=1) =", float(mp.quad(lambda r: 4*mp.pi*r**2*(4*mp.pi)**mp.mpf(-1.5)*mp.e**(-r**2/4), [0,2,12,80])))
print("   int d_1 G = 0 by oddness in x_1 (exact)")
rec('G_L2_tau1', (8*mp.pi)**mp.mpf(-0.75), "||G(.,1)||_2 = (8 pi tau)^{-3/4}")
# independent numeric ||G(.,1)||_2
n2 = mp.sqrt(mp.quad(lambda r: 4*mp.pi*r**2*((4*mp.pi)**mp.mpf(-1.5)*mp.e**(-r**2/4))**2, [0,2,12,80]))
print("   ||G(.,1)||_2 numeric =", float(n2), " closed =", float((8*mp.pi)**mp.mpf(-0.75)))
assert abs(n2-(8*mp.pi)**mp.mpf(-0.75)) < mp.mpf('1e-30')

print()
print("=== B. Q = sup |(n.u)w - (n.w)u|/(|u||w|) : brute force on the sphere, my parameterisation ===")
rng = np.random.default_rng(20260906)
best = 0.0
for _ in range(60):
    U = rng.normal(size=(200000,3)); W = rng.normal(size=(200000,3)); N = rng.normal(size=(200000,3))
    U /= np.linalg.norm(U,axis=1,keepdims=True); W /= np.linalg.norm(W,axis=1,keepdims=True)
    N /= np.linalg.norm(N,axis=1,keepdims=True)
    V = (np.einsum('ij,ij->i',N,U))[:,None]*W - (np.einsum('ij,ij->i',N,W))[:,None]*U
    best = max(best, float(np.linalg.norm(V,axis=1).max()))
rec('Q_random_max_12e6', best, "must be <= 1")
assert best <= 1.0 + 1e-12
# exact structured search: objective alpha^2+beta^2-2 c alpha beta on the admissible ellipse
c_s, al, be, lam = sp.symbols('c alpha beta lam', real=True)
obj = al**2+be**2-2*c_s*al*be
con = 1-c_s**2-(al**2-2*c_s*al*be+be**2)
print("   objective - (1-c^2) on the constraint boundary:", sp.simplify(obj-(1-c_s**2)-(-con)))
rec('Q_exact', 1, "sup_n = 1-c^2 <= 1, max over c at c=0")
# attained
u=np.array([1,0,0.]); w=np.array([0,1,0.]); nn=np.array([1,1,0.])/np.sqrt(2)
rec('Q_attained', np.linalg.norm(nn.dot(u)*w-nn.dot(w)*u))

print()
print("=== C. C_u : my own symbolic optimisation, then an independent 1-D numeric min ===")
tau, A, B = sp.symbols('tau A B', positive=True)
f = A*sp.sqrt(tau) + B*tau**sp.Rational(-3,4)
tstar = sp.solve(sp.diff(f,tau), tau)[0]
fmin = sp.simplify(sp.powsimp(f.subs(tau,tstar), force=True))
print("   tau_* =", sp.simplify(tstar))
print("   f(tau_*) =", fmin)
gamma_mine = sp.simplify(fmin/(A**sp.Rational(3,5)*B**sp.Rational(2,5)))
print("   gamma =", sp.nsimplify(sp.simplify(gamma_mine)), "=", float(sp.N(gamma_mine,30)))
rec('gamma', sp.N(gamma_mine,30))
Amine = 4/sp.sqrt(sp.pi); Bmine = (8*sp.pi)**sp.Rational(-3,4)
Cu_mine = sp.simplify(gamma_mine*Amine**sp.Rational(3,5)*Bmine**sp.Rational(2,5))
Cu_mine = sp.simplify(sp.powsimp(Cu_mine, force=True))
print("   C_u (symbolic) =", Cu_mine, "=", float(sp.N(Cu_mine,30)))
rec('C_u', sp.N(Cu_mine,30))
# numeric min, different optimiser (mpmath findroot on f')
Af = 4/mp.sqrt(mp.pi); Bf = (8*mp.pi)**mp.mpf(-0.75)
g = lambda tv: Af*mp.sqrt(tv)+Bf*tv**mp.mpf(-0.75)
fp = lambda tv: Af/(2*mp.sqrt(tv)) - mp.mpf('0.75')*Bf*tv**mp.mpf('-1.75')
lo, hi = mp.mpf('1e-8'), mp.mpf('1e6')
for _ in range(400):
    mid=(lo+hi)/2
    if fp(mid) < 0: lo=mid
    else: hi=mid
ts = (lo+hi)/2
rec('tau_star', ts); rec('C_u_numeric', g(ts))
assert abs(g(ts)-mp.mpf(str(float(sp.N(Cu_mine,30))))) < mp.mpf('1e-14')
# closed form claimed by the note
Cu_claim = 5*2**sp.Rational(9,10)*3**sp.Rational(2,5)/(6*sp.pi**sp.Rational(3,5))
print("   note's closed form 5*2^(9/10)3^(2/5)/(6 pi^(3/5)) - mine =",
      sp.simplify(sp.powsimp(sp.nsimplify(sp.simplify(Cu_claim-Cu_mine)),force=True)))
assert abs(float(sp.N(Cu_claim-Cu_mine,40))) < 1e-35

print()
print("=== D. c_div : my own bootstrap, symbolically, from scratch ===")
M0,E0,nu,t,Cu = sp.symbols('M0 E0 nu t C_u', positive=True)
Q = 1
Cnu = 2*Q/sp.sqrt(sp.pi*nu)
uinf = Cu*E0**sp.Rational(1,5)*(2*M0)**sp.Rational(3,5)
duh = Cnu*(2*sp.sqrt(t))*(2*M0)*uinf              # K5: int_0^t (t-s)^{-1/2} ds = 2 sqrt t
print("   Duhamel(t) =", sp.simplify(duh))
sol = sp.solve(sp.Eq(duh, M0/2), t)
tstar_sym = sp.simplify(sol[0]).subs(Cu, Cu_mine)
tstar_sym = sp.simplify(sp.powsimp(tstar_sym, force=True))
print("   t_* =", tstar_sym)
ReE = E0**sp.Rational(2,5)*M0**sp.Rational(1,5)/nu
cdiv_mine = sp.simplify(sp.powsimp(sp.simplify(M0*ReE*tstar_sym), force=True))
print("   c_div = M0 * Re_E * t_* =", cdiv_mine, "=", float(sp.N(cdiv_mine,30)))
rec('c_div', sp.N(cdiv_mine,30))
cdiv_claim = 3*3**sp.Rational(1,5)*sp.pi**sp.Rational(11,5)/12800
print("   note's 3*3^(1/5)pi^(11/5)/12800 - mine =", sp.simplify(cdiv_claim-cdiv_mine))
assert abs(float(sp.N(cdiv_claim-cdiv_mine,40))) < 1e-35
rec('c_div_over_piCu2', sp.N(sp.pi/(2**sp.Rational(46,5)*Cu_mine**2),30), "= pi/(2^{46/5}C_u^2)")
rec('inv_Cu2', sp.N(1/Cu_mine**2,30))
rec('pi_2_m46_5', sp.N(sp.pi*2**sp.Rational(-46,5),30))
# Q=2 variant
duh2 = (2*2/sp.sqrt(sp.pi*nu))*(2*sp.sqrt(t))*(2*M0)*uinf
t2 = sp.simplify(sp.solve(sp.Eq(duh2,M0/2), t)[0]).subs(Cu,Cu_mine)
rec('c_div_Q2', sp.N(sp.simplify(M0*ReE*t2),30)); rec('ratio_Q1_Q2', sp.N(cdiv_mine/sp.simplify(M0*ReE*t2),30))

print()
print("=== E. THE DIMENSIONAL CLAIM the note makes in its headline ===")
Es,Ms,nus = sp.symbols('E M nu', positive=True)
Re = Es**sp.Rational(2,5)*Ms**sp.Rational(1,5)/nus
expr = sp.sqrt(nus)/(Es**sp.Rational(1,5)*Ms**sp.Rational(1,10))
print("   sqrt(nu)/(E^{1/5}M^{1/10})  minus  Re_E^{-1/2}  =",
      sp.simplify(sp.powsimp(expr - Re**sp.Rational(-1,2), force=True)))
# dimensions: [E]=L^5 T^-2, [M]=T^-1, [nu]=L^2 T^-1
Ls,Ts = sp.symbols('L T', positive=True)
dim = {Es: Ls**5*Ts**-2, Ms: Ts**-1, nus: Ls**2*Ts**-1}
print("   dimension of sqrt(nu)/(E^{1/5}M^{1/10}) =", sp.simplify(expr.subs(dim)))
print("   dimension of Re_E                        =", sp.simplify(Re.subs(dim)))
rec('dim_of_sqrtnu_over_E15M110_is_one', 1 if sp.simplify(expr.subs(dim))==1 else 0,
    "1 == the expression IS dimensionless, i.e. IS a pure number")

print()
print("=== F. the campaign datum, my own arithmetic ===")
rows=[]
for L in [mp.mpf('8.3178'), mp.mpf(20), mp.mpf(50), mp.mpf(100)]:
    logRe = 2*L - mp.mpf('0.7031660'); Re_ = mp.e**logRe
    powwin = mp.mpf(str(float(sp.N(cdiv_mine,30))))/Re_
    logwin = 1/(1+logRe); model = 2*mp.log(2)/L
    rows.append([float(L),float(logRe),float(Re_),float(powwin),float(logwin),float(model),
                 float(model/powwin), float(model/logwin)])
    print("   L=%-9g logRe=%-12.6f Re=%-14.6e pow=%-14.6e log/c1=%-12.6e model=%-10.6f model/pow=%-12.6e model/log=%.6f"
          % tuple(rows[-1]))
R['datum_rows']=rows
print("   limit of model/log as L->inf = 2ln2*2 = 4ln2 =", float(4*mp.log(2)))
print("   NOTE: model/log = (2ln2/L)(1+2L-0.7031660) uses ONLY two numbers quoted from")
print("         sharp/exact-first-order; c_div does not appear.  Zero input from this note.")

print()
print("=== G. crossover and the 'log beats power for Re>=1' claim ===")
h = lambda r: r/(1+mp.log(r))
print("   min of Re/(1+log Re) on [1,inf) :", [ (float(r), float(h(mp.mpf(r)))) for r in [1,1.5,2,10,1e6] ])
for c1 in ['1e-3','1e-6']:
    f = lambda r: h(r) - mp.mpf(str(float(sp.N(cdiv_mine,30))))/mp.mpf(c1)
    r0 = mp.findroot(f, mp.mpf(10))
    print("   c_1=%s -> crossover Re_E = %.6g" % (c1, float(r0)))
    R['crossover_'+c1]=float(r0)

json.dump(R, open('r1_results.json','w'), indent=1, default=float)
print("\nr1: independent recomputation complete")

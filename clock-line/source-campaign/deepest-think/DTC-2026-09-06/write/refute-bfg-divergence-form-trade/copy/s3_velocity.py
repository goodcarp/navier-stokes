#!/usr/bin/env python3
"""
s3_velocity.py -- the energy-dependent velocity bound, with its exact constant.

LEMMA U.  Let u be smooth, divergence free, with u in L^2(R^3) and omega = curl u
in L^infty, and set E := ||u||_2^2 (the estate convention: no 1/2), M := ||omega||_inf.
Split with the (auxiliary, non-dynamical) heat semigroup at scale tau:

    u = e^{tau Delta} u + (u - e^{tau Delta} u),
    (i)  ||e^{tau Delta}u||_inf <= ||G(.,tau)||_2 ||u||_2 = (8 pi tau)^{-3/4} E^{1/2}
    (ii) u - e^{tau Delta}u = int_0^tau curl(G_sigma * omega) dsigma      [Delta u = -curl omega]
         ==> ||u - e^{tau Delta}u||_inf <= M int_0^tau 2/sqrt(pi sigma) dsigma = (4/sqrt(pi)) M sqrt(tau)

    ||u||_inf <= A tau^{1/2} + B tau^{-3/4},   A = (4/sqrt(pi)) M,  B = (8 pi)^{-3/4} E^{1/2}.

Minimising over tau > 0 gives  ||u||_inf <= C_u E^{1/5} M^{3/5}  with C_u derived here.

Also computed, for contrast: the same split against ||omega||_2 instead of ||u||_2,
which gives ||u||_inf <= C_w M^{1/3} ||omega||_2^{2/3} -- a bound that CANNOT be used in a
clock because enstrophy is not non-increasing under NS, while energy is.
"""
import json
import sympy as sp
import mpmath as mp
mp.mp.dps = 40
out = {}

A, B, tau = sp.symbols('A B tau', positive=True)
f = A*sp.sqrt(tau) + B*tau**sp.Rational(-3,4)
crit = sp.solve(sp.diff(f, tau), tau)
tau_star = sp.simplify(crit[0])
print("tau_*            =", tau_star)
fmin = sp.simplify(sp.powsimp(sp.simplify(f.subs(tau, tau_star)), force=True))
print("min_tau f(tau)   =", fmin)

gamma = sp.simplify(sp.powsimp(fmin/(A**sp.Rational(3,5)*B**sp.Rational(2,5)), force=True))
gamma = sp.nsimplify(sp.simplify(gamma))
print("f_min / (A^(3/5) B^(2/5)) =", gamma, " = ", sp.simplify(gamma - ((sp.Rational(3,2))**sp.Rational(2,5) + (sp.Rational(2,3))**sp.Rational(3,5))))
assert sp.simplify(gamma - ((sp.Rational(3,2))**sp.Rational(2,5) + (sp.Rational(2,3))**sp.Rational(3,5))) == 0
gamma_val = mp.mpf(3)/2
gamma_num = (mp.mpf(3)/2)**(mp.mpf(2)/5) + (mp.mpf(2)/3)**(mp.mpf(3)/5)
print("gamma            = (3/2)^(2/5) + (2/3)^(3/5) = %.15f" % float(gamma_num))
out['gamma'] = float(gamma_num)
out['tau_star_formula'] = "(3B/(2A))^(4/5)"

# substitute A = (4/sqrt(pi)) M,  B = (8 pi)^{-3/4} E^{1/2}
M, E = sp.symbols('M E', positive=True)
Aex = 4/sp.sqrt(sp.pi)*M
Bex = (8*sp.pi)**sp.Rational(-3,4)*sp.sqrt(E)
bound = sp.simplify(sp.powsimp(fmin.subs({A: Aex, B: Bex}), force=True))
print("||u||_inf bound  =", bound)
Cu_sym = sp.simplify(sp.powsimp(bound/(E**sp.Rational(1,5)*M**sp.Rational(3,5)), force=True))
print("C_u (symbolic)   =", Cu_sym)
Cu = mp.mpf(str(sp.N(Cu_sym, 40)))
# independent numeric route: minimise numerically at M=E=1
g = lambda tv: (4/mp.sqrt(mp.pi))*mp.sqrt(tv) + (8*mp.pi)**mp.mpf(-0.75)*tv**mp.mpf(-0.75)
# independent 1-D minimisation by golden section (no closed form used)
lo, hi = mp.mpf('1e-6'), mp.mpf('1e3')
phi = (mp.sqrt(5)-1)/2
for _ in range(400):
    x1 = hi - phi*(hi-lo); x2 = lo + phi*(hi-lo)
    if g(x1) < g(x2): hi = x2
    else: lo = x1
tstar = (lo+hi)/2
Cu_num = g(tstar)
print("C_u              = %.15f   (symbolic)" % float(Cu))
print("C_u              = %.15f   (numeric 1-D minimisation, M=E=1, tau_*=%.10f)" % (float(Cu_num), float(tstar)))
assert abs(Cu - Cu_num) < mp.mpf('1e-25')
# closed form check: C_u = gamma * (4/sqrt(pi))^{3/5} * (8 pi)^{-3/10}
Cu_closed = gamma_num*(4/mp.sqrt(mp.pi))**(mp.mpf(3)/5)*(8*mp.pi)**(-mp.mpf(3)/10)
print("C_u closed form  = gamma*(4/sqrt(pi))^(3/5)*(8pi)^(-3/10) = %.15f" % float(Cu_closed))
assert abs(Cu - Cu_closed) < mp.mpf('1e-30')
out['C_u'] = float(Cu)
out['C_u_closed_form'] = "((3/2)^(2/5)+(2/3)^(3/5)) * (4/sqrt(pi))^(3/5) * (8*pi)^(-3/10)"

# dimension / scaling check: u_lam(x,t)=lam u(lam x, lam^2 t)
lam = sp.Symbol('lam', positive=True)
# M -> lam^2 M, E = ||u||_2^2 -> lam^2 * lam^{-3} E = lam^{-1} E, ||u||_inf -> lam ||u||_inf
lhs = lam*(E**sp.Rational(1,5)*M**sp.Rational(3,5))
rhs = (lam**-1*E)**sp.Rational(1,5)*(lam**2*M)**sp.Rational(3,5)
print("scaling residual (E^{1/5}M^{3/5} is exactly of velocity dimension):", sp.simplify(lhs-rhs))
assert sp.simplify(lhs-rhs) == 0
out['scaling_residual'] = 0

# ---- the enstrophy variant (for contrast, NOT used in the clock) ---------
# ||e^{tau D}u||_inf <= ||K_tau||_2 ||omega||_2 with K_tau = e^{tau D} * BiotSavart kernel.
# ||K_tau||_2^2 = (2 pi)^{-3} int |xi|^{-2} e^{-2 tau |xi|^2} dxi = (2 pi)^{-3} 4 pi int_0^inf e^{-2 tau q^2} dq
xi, q = sp.symbols('xi q', positive=True)
K2sq = sp.simplify((2*sp.pi)**-3*4*sp.pi*sp.integrate(sp.exp(-2*tau*q**2), (q, 0, sp.oo)))
K2n = sp.simplify(sp.sqrt(K2sq))
print("||K_tau||_2      =", sp.simplify(sp.powsimp(K2n, force=True)))
Aw = 4/sp.sqrt(sp.pi)*M; Bw = K2n
fw = Aw*sp.sqrt(tau) + Bw*sp.Symbol('W', positive=True)      # W = ||omega||_2
tw = sp.symbols('tw', positive=True)
W = sp.Symbol('W', positive=True)
fw = Aw*sp.sqrt(tw) + K2n.subs(tau, tw)*W
tw_star = sp.solve(sp.diff(fw, tw), tw)[0]
fw_min = sp.simplify(sp.powsimp(sp.simplify(fw.subs(tw, tw_star)), force=True))
Cw = sp.simplify(sp.powsimp(fw_min/(M**sp.Rational(1,3)*W**sp.Rational(2,3)), force=True))
print("||u||_inf <= C_w M^{1/3}||omega||_2^{2/3},  C_w =", Cw, "= %.12f" % float(sp.N(Cw, 30)))
out['C_w_enstrophy_variant'] = float(sp.N(Cw, 30))
out['note_enstrophy'] = "enstrophy is not non-increasing under NS; energy is. This variant cannot propagate."

json.dump(out, open('s3_results.json','w'), indent=1)
print("\ns3: LEMMA U constant derived two independent ways and asserted equal")

#!/usr/bin/env python3
"""G6 — a DYNAMICAL falsifier control for the funnel: Hou's generalized-dimension axisymmetric family (arXiv 2405.10916:
nearly self-similar blowup at fractional dimension n = 3.188, profile solving constant-viscosity equations in that dimension).
Any step of the estate's (K)/(K-core) chain that is insensitive to n cannot be the decisive one (FL-043 vacuity check, one level up).
Exact facts computed here for the n-dimensional axisymmetric no-swirl kernel (lift to R^{n+2}, r radial in R^{n+1}):
  a(0,z) = (|S^n|/|S^{n+1}|) ∫∫ eta r'^n (z-z')/(r'^2+(z-z')^2)^{(n+2)/2} dr' dz'
  bang-bang extremal: a(0) = kappa0(n) M log(R/rho0),  kappa0(n) = (2/n) |S^n|/|S^{n+1}| = (2/n) Gamma((n+2)/2)/(sqrt(pi) Gamma((n+1)/2))
  energy capacity radius: delta(n) = M^{-2/(n+2)};   ell = M^{-1/2} (n-free);   L(n) = log(delta/ell) = ((n-2)/(2(n+2))) log M
  literal-reading deficit (estate: demand (1/10) log M over ONE epoch 1/M vs supply kappa0 M L log M / M) = 1/kappa0(n)  (= 2 at n=3)
Controls: n=3 must reproduce kappa0 = 1/2, delta = M^{-2/5}, L = (1/10) log M, deficit = 2 exactly."""
import sympy as sp, sys, hashlib
n = sp.symbols('n', positive=True)
kappa0 = (2/n)*sp.gamma((n+2)/2)/(sp.sqrt(sp.pi)*sp.gamma((n+1)/2))
delta_exp = -2/(n+2); L_coeff = (n-2)/(2*(n+2)); deficit = 1/kappa0
k3 = sp.simplify(kappa0.subs(n,3)); d3 = sp.simplify(delta_exp.subs(n,3)); L3 = sp.simplify(L_coeff.subs(n,3)); df3 = sp.simplify(deficit.subs(n,3))
print('n=3 :  kappa0 =', k3, ' delta exponent =', d3, ' L coefficient =', L3, ' deficit =', df3)
nH = sp.Rational(3188,1000)
kH = sp.N(kappa0.subs(n,nH), 6); dH = sp.N(delta_exp.subs(n,nH), 6); LH = sp.N(L_coeff.subs(n,nH), 6); dfH = sp.N(deficit.subs(n,nH), 6)
print('n=3.188 (Hou blowup dimension):  kappa0 =', kH, ' delta exponent =', dH, ' L coefficient =', LH, ' deficit =', dfH)
# angular integral identity used: ∫_0^pi sin^{n-1} phi |cos phi| dphi = 2/n
phi = sp.symbols('phi', positive=True)
ang3 = sp.integrate(sp.sin(phi)**2*sp.Abs(sp.cos(phi)), (phi, 0, sp.pi))
print('angular integral n=3: ∫ sin^2 |cos| =', ang3, '(must be 2/3)')
# Hardy margin of the estate's W = r^9 h conjugation: effective dimension d = 2k-2 with k=9 -> 16, coefficient (k-1)(k-3)=48 vs Hardy((d-2)/2)^2=49
d_eff = lambda k: 2*k - 2
hardy = lambda d: sp.Rational((d-2)**2, 4)
print('estate Hardy margin at n=3: coefficient (k-1)(k-3) =', (9-1)*(9-3), ' vs Hardy(d=16) =', hardy(16), ' margin', hardy(16)-(9-1)*(9-3))
# in dimension n the lift is R^{n+2} (radial in R^{n+1}); the same conjugation exponent k gives effective dimension n+1+2(k-1)... report the n-derivative of the margin qualitatively:
print('n-sensitivity: kappa0, delta, L, deficit vary smoothly with n; the one-unit Hardy margins are the n-sensitive steps (a shift of the lift dimension by 0.188 moves a Hardy constant by ~ (d-2)/2 * 0.188 ~ 1.3 units at d=16: the margin 49-48=1 does NOT survive)')
ok = (k3 == sp.Rational(1,2)) and (d3 == sp.Rational(-2,5)) and (L3 == sp.Rational(1,10)) and (df3 == 2) and (ang3 == sp.Rational(2,3)) and (abs(float(kH) - 0.479) < 0.01)
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

#!/usr/bin/env python3
"""G8 — the direct-cap 1/20 clock as an exponent identity, and its dimension-dependence.
Source structure (Sol's frozen memo 085cfd46…, §2–§3, as audited in AUDIT_DIRECTCAP_CORE.md): exponents of M —
  P = M^p (detector mass, p = -7/5), H = M^h (block length, h = -19/20), delta = M^{d} (terminal detector radius, d = -2/5),
  core split ell = M^{-beta}; exterior bound |W_{>=R}| <~ G M P H R^{-3}  (power 3 from |q| <= G/r^2, |r q_z| <= M);
  annular clock  = log(delta/R_eps) = (d - (1+p+h)/3) log M  - eps log M;
  core clock     = (1/2) log( m2(tau*)/m2(0) ),  m2 = ∫ R r^{-2} dY (r^{-2} = the R^4-harmonic weight; drift coefficient 2),
                   m2(tau*) >~ ell^{-1} M^{-2} H^{-1} (bare dr dz measure after the exact r^3 r^{-2} r cancellation),
                   m2(0) <~ P delta^{-2}   =>  core clock = (1/2)[beta - 2 - h - p + 2 d].
Dimension n generalization of ONLY the structurally forced pieces (flagged): the lift is R^{n+2} (radial part R^{n+1}),
harmonic weight r^{-(n-1)}, drift coefficient (n-1); the core integrand weight becomes r^{n-3}; delta(n) = M^{-2/(n+2)} (energy capacity);
the exterior power 3 and the cap power 2 are algebraic (kept); p and h are DETECTOR-DEFINED and NOT re-derived here (held at n=3 values,
and alternatively p(n) = n*d(n) - 1/5 if P = phi*delta^n with phi = M^{-1/5}).
Checks that can fail: (i) n=3 reproduces R_eps = -9/20, both clocks = 1/20, equalizer beta* = 11/20 (the audit's independent 'bonus fact');
(ii) the equalizer formula beta*(n) = [2 + h + p - (n-1)(1+p+h)/3]/(n-2) gives 11/20 at n=3;
(iii) the n=3.188 values are then compared with kappa0(n)L(n) = 0.0553, L(n)/2 = 0.0572, and 1/20."""
from fractions import Fraction as F
import math, sys, hashlib
def kappa0(n): return (2/n)*math.gamma((n+2)/2)/(math.sqrt(math.pi)*math.gamma((n+1)/2))
def clocks(n, p, h, d, beta):
    R_exp = (1 + p + h)/3
    ann = d - R_exp                       # minus eps
    core = ((n-2)*beta - 2 - h - p + (n-1)*d)/(n-1)
    return R_exp, ann, core
def beta_star(n, p, h): return (2 + h + p - (n-1)*(1+p+h)/3)/(n-2)
p3, h3, d3 = F(-7,5), F(-19,20), F(-2,5)
R3, ann3, core3 = clocks(F(3), p3, h3, d3, F(11,20))
b3 = beta_star(F(3), p3, h3)
print('n=3 ledger: R_eps exponent =', R3, '| annular clock =', ann3, '| core clock at beta=11/20 =', core3, '| equalizer beta* =', b3)
ok = (R3 == F(-9,20)) and (ann3 == F(1,20)) and (core3 == F(1,20)) and (b3 == F(11,20))
# exact core formula check against the audit's stated "1/2(beta + 7/20 - 2 gamma)" with gamma = -d = 2/5
beta = F(11,20); gamma = F(2,5)
audit_core = F(1,2)*(beta + F(7,20) - 2*gamma)
print('audit core formula (1/2)(beta+7/20-2gamma) =', audit_core, '| equals structural core:', audit_core == core3)
ok &= (audit_core == core3)
# dimension n = 3.188, two pinning choices for p
n = 3.188; d = -2/(n+2)
for label, p in [('p held at -7/5', -7/5), ('p(n) = n*d(n) - 1/5 (P = phi delta^n, phi = M^{-1/5})', n*d - 0.2)]:
    h = -19/20
    bs = beta_star(n, p, h)
    R, ann, core = clocks(n, p, h, d, bs)
    print(f'n=3.188 [{label}]: d = {d:.4f}, p = {p:.4f}, R_eps = {R:.4f}, beta* = {bs:.4f}, annular = core = {ann:.4f} (core {core:.4f})')
L = (n-2)/(2*(n+2)); print(f'compare: kappa0(n)L(n) = {kappa0(n)*L:.4f}, L(n)/2 = {L/2:.4f}, fixed 1/20 = 0.0500')
# sensitivity: which single exponent moves the clock most (partial derivatives of the annular clock)
print('annular clock partials: d(ann)/d(d) = +1, d(ann)/d(p) = -1/3, d(ann)/d(h) = -1/3  -> the clock is set by delta and the detector exponents, not by the axis kernel')
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

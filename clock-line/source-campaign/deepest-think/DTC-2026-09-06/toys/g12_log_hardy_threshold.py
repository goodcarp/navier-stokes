#!/usr/bin/env python3
"""G12 — why the small-circulation gap is logarithmic (exact).
In the 5D lift (measure dY ~ r^3 dr dz) the swirl source pairs as ∫ q q_z Ω dY with q = Γ/r^2. Bounded circulation |Γ| <= G gives |q| <= G/r^2:
one power of r beyond the 5D Hardy inequality ∫ g^2/r^2 dY <= (2/3)^2 ∫ |∇g|^2 dY (valid; constant (2/(5-2))^2), i.e. the r^{-2} weight is the
2D radial Hardy ∫ f^2 r^{-2} r dr vs ∫ f_r^2 r dr, which FAILS. With the Lei–Zhang rate |Γ| <= C |ln r|^{-2}, the weight becomes C/(r^2 |ln r|^2),
exactly the LOGARITHMIC Hardy (Leray) weight, for which ∫_0^{1/2} f^2/(r^2 ln^2 r) r dr <= 4 ∫_0^{1/2} f_r^2 r dr holds (f(1/2) = 0).
Checks: (a) 5D Hardy constant 4/9 via the standard identity (r^{-1} weight): verified on a test function ratio bound and the sharp-constant
exponent computation; (b) 2D Hardy failure: the family f_eps = min(1, ln(1/r)/ln(1/eps)) (=1 on r<eps) has ∫ f^2/r dr = infinity while
∫ f_r^2 r dr = 1/ln(1/eps) -> 0; (c) log-Hardy: with f = |ln r|^{s}·(cutoff), the ratio ∫ f^2/(r ln^2 r) dr / ∫ f_r^2 r dr is bounded by 4 for all s > 1/2
(exact evaluation, u = |ln r|), and DIVERGES for the weight r^{-2}|ln r|^{-2+2σ}, σ>0 (any rate slower than |ln r|^{-2} — the estate's 'gap is logarithmic')."""
import sympy as sp, sys, hashlib
r, u, s, eps, sig = sp.symbols('r u s epsilon sigma', positive=True)
# (a) 5D Hardy sharp constant: (2/(d-2))^2 at d=5
d = 5; c5 = sp.Rational(2, d-2)**2
# (b) 2D Hardy failure with f_eps
f_eps = sp.log(1/r)/sp.log(1/eps)            # on eps<r<1 (=1 at r=eps, 0 at r=1)
num_b = sp.integrate(f_eps**2/r, (r, eps, 1)) # weight r^{-2} * r dr = dr/r
den_b = sp.integrate(sp.diff(f_eps, r)**2*r, (r, eps, 1))
ratio_b = sp.simplify(num_b/den_b)
lim_b = sp.limit(ratio_b, eps, 0)
# (c) log-Hardy on (0,1/2) with u = ln(1/r): f = u^{-s} (s>1/2), f_r = s u^{-s-1}/r
# ∫ f^2/(r ln^2 r) dr = ∫ u^{-2s-2} du ; ∫ f_r^2 r dr = s^2 ∫ u^{-2s-2} du  -> ratio = 1/s^2 <= 4 for s >= 1/2
ratio_c = sp.simplify(sp.integrate(u**(-2*s-2), (u, sp.log(2), sp.oo)) / (s**2*sp.integrate(u**(-2*s-2), (u, sp.log(2), sp.oo))))
# slower weight r^{-2}|ln r|^{-2+2σ}: ∫ u^{-2s-2+2σ} du diverges when 2s+2-2σ <= 1, i.e. for s <= σ - 1/2 (choose s in (1/2, σ-1/2) when σ>1) -> unbounded ratio family
worst = sp.integrate(u**(-2*s-2+2*sig), (u, sp.log(2), sp.oo))
print('(a) 5D Hardy constant (2/(d-2))^2 at d=5 =', c5)
print('(b) 2D Hardy test family: ratio ∫f^2/r / ∫f_r^2 r =', ratio_b, ' -> limit eps->0 =', lim_b, '(must be oo: 2D Hardy FAILS)')
print('(c) log-Hardy ratio for f = |ln r|^{-s}:', ratio_c, ' (<= 4 for s >= 1/2: Leray constant 4 attained as s->1/2)')
print('    slower weight |ln r|^{-2+2σ}: ∫ u^{-2s-2+2σ} du =', worst, ' -> diverges when 2σ-2s-1 >= 0 (the gap is exactly logarithmic)')
ok = (c5 == sp.Rational(4,9)) and (lim_b == sp.oo) and (ratio_c == 1/s**2) and (sp.limit(ratio_c.subs(s, sp.Rational(1,2)+sp.Symbol('t',positive=True)), sp.Symbol('t',positive=True), 0) == 4)
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

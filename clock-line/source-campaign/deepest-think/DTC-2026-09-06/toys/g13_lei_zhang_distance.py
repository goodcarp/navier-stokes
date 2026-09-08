#!/usr/bin/env python3
"""G13 — exact distance of the record scenario from the Lei–Zhang regularity rate, and the FBC (1.6) form-ratio for the cap-saturating swirl.
(a) Lei–Zhang Cor. 1.3 rate at the core scale ell = M^{-1/2}: |Gamma| <= C |ln ell|^{-2} = 4C/(ln M)^2.  Scenario: Gamma ~ G = O(1) at r ~ ell.
    Distance = G (ln M)^2 / (4C): the scenario exceeds the known-regular rate by (ln M)^2 at its own core scale.
(b) FBC (1.6) form ratio for the cap-saturating swirl |v^theta| = min(G/r, M r/2) (cross-over at r_G = sqrt(2G/M)), test function f = 1 on r < r0
    (dx ~ r dr in the meridional plane; the r-derivative part vanishes inside, so the ratio is the bare weighted mass): 
    ∫_{r<r0} |v^theta|^2 r dr = ∫_0^{r_G} (M^2 r^2/4) r dr + ∫_{r_G}^{r0} (G^2/r^2) r dr = G^2/4 * ... + G^2 ln(r0/r_G)  -> the log(M) term is G^2 * (1/2) ln(M r0^2/(2G)).
    So at fixed circulation G the (1.6) mass grows like (G^2/2) ln M: the FBC's smallness delta* is beaten by exactly one logarithm with coefficient G^2/2.
Controls: at G -> 0 the ratio -> 0 (small swirl passes); the r_G-inner piece is M-independent (G^2/4 exactly)."""
import sympy as sp, sys, hashlib
r, M, G, C, r0 = sp.symbols('r M G C r0', positive=True)
ell = 1/sp.sqrt(M)
rate_at_ell = C*sp.log(1/ell)**(-2)
dist = sp.simplify(G/rate_at_ell)
print('(a) Lei–Zhang rate at ell:', sp.simplify(rate_at_ell), ' | scenario/rate =', dist)
rG = sp.sqrt(2*G/M)
inner = sp.integrate((M**2*r**2/4)*r, (r, 0, rG))
outer = sp.integrate((G**2/r**2)*r, (r, rG, r0))
mass = sp.simplify(inner + outer)
print('(b) inner piece =', sp.simplify(inner), ' | outer piece =', sp.simplify(outer), ' | total =', mass)
lead = sp.simplify(sp.expand_log(mass, force=True))
print('    expanded:', lead)
ctrl0 = sp.limit(mass, G, 0)
ok = (sp.simplify(dist - G*sp.log(M)**2/(4*C)) == 0) and (sp.simplify(inner - G**2/4) == 0) and (ctrl0 == 0) and (sp.simplify(sp.diff(outer, M) - G**2/(2*M)) == 0)
print('CTRL G->0 mass =', ctrl0, '(must be 0); d(outer)/dM = G^2/(2M) (the (G^2/2) ln M growth):', sp.simplify(sp.diff(outer, M)))
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

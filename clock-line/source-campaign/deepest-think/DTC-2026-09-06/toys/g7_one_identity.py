#!/usr/bin/env python3
"""G7 — are the three '1/20' coincidences one identity?  Exact Fractions.
Objects (exponents of M): ell = -1/2 (viscous core), delta = -2/(n+2) (energy capacity), L = delta - ell = log-count coefficient,
kappa0 = axis strain per e-fold per unit M.  Estate numbers at n=3: direct-cap toll per block (1/20) log M; telescope window gamma < 1/20;
localisation ladder stepping by M^{-1/20} (delta=M^{-8/20}, M^{-9/20} penetration, ell=M^{-10/20}, M^{-11/20} movable core, s_min=M^{-12/20}).
Claims checked exactly:
  (a) one epoch (1/M) of full-log extremal strain kappa0*M*L*log M integrates to kappa0*L*log M = (1/20) log M at n=3  -> equals the direct-cap toll VALUE
  (b) sqrt(delta*ell) = M^{-9/20} = delta * M^{-1/20} = ell * M^{+1/20}: the unprotectable-shell radius IS the ladder rung between delta and ell
  (c) the ladder step 1/20 = (delta - ell)/2 = L/2 = kappa0*L  (only because kappa0 = 1/2 exactly at n=3)
Prediction for the registrar (testable against Sol's direct-cap derivation): if the toll is kappa0(n)*L(n) in dimension n, then at n=3.188 it is 0.0553 log M, not 0.05.
Control: at n=4 (kappa0 = (2/4)*Gamma(3)/(sqrt(pi)*Gamma(5/2)) = 1/2 * 2/(sqrt(pi)*3 sqrt(pi)/4) = 4/(3 pi)) the identity (c) must FAIL (kappa0 != 1/2)."""
from fractions import Fraction as F
import math, sys, hashlib
def kappa0(n): return (2/n)*math.gamma((n+2)/2)/(math.sqrt(math.pi)*math.gamma((n+1)/2))
ell = F(-1,2); delta = F(-2,5); L = delta - ell
k3 = F(1,2); assert abs(kappa0(3.0) - 0.5) < 1e-12
a = k3*L; b = (delta+ell)/2; c = L/2
print('(a) one-epoch full-log strain integral coefficient kappa0*L =', a, '| direct-cap toll 1/20:', a == F(1,20))
print('(b) sqrt(delta*ell) exponent =', b, '| = delta - 1/20:', b == delta - F(1,20), '| = ell + 1/20:', b == ell + F(1,20))
print('(c) ladder step L/2 =', c, '| equals kappa0*L:', c == a, '(holds iff kappa0 = 1/2)')
nH = 3.188; LH = (nH-2)/(2*(nH+2)); tollH = kappa0(nH)*LH
print(f'prediction at n=3.188: kappa0 = {kappa0(nH):.4f}, L = {LH:.4f}, kappa0*L = {tollH:.4f} log M  (vs 1/20 = 0.05); ladder half-step L/2 = {LH/2:.4f}')
k4 = kappa0(4.0); print(f'CTRL n=4: kappa0 = {k4:.4f} (exact 4/(3pi) = {4/(3*math.pi):.4f}); kappa0 == 1/2? {abs(k4-0.5)<1e-9} (must be False)')
ok = (a == F(1,20)) and (b == F(-9,20)) and (c == a) and abs(k4 - 4/(3*math.pi)) < 1e-12 and abs(k4-0.5) > 1e-3 and abs(tollH-0.0553) < 5e-4
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

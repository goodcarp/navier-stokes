"""
REFUTER, second numeric pass.  My own first pass (rf1, mpmath with breakpoints [0,1,10,R-1,inf])
LOST accuracy at R=1e32 and broke down at R=1e64 -- J0 fell to 12.61, breaking monotonicity.
This script re-does the same integrals on log-spaced panels and, more importantly, computes the
R -> infinity LIMITS in closed/quadrature form so the whole J0/J1 table can be checked against
something that does not depend on any panel choice.

Limits (h_R(r) = (1 - log(1+r)/log R)_+, w(r) = 1/(r^4+1), c_R = avg_{B_1} h_R):
  J0(R)         -> 4 pi int_0^inf r^2 w dr = 4 pi * (pi/(2 sqrt 2)) = sqrt(2) pi^2   (monotone up)
  c_R           = 1 - A/log R,  A = 3 int_0^1 log(1+r) r^2 dr = 2 log 2 - 5/6
  J1(R)*log R   -> 4 pi int_0^inf |A - log(1+r)| r^2 w dr
"""
import json, math
from mpmath import mp, mpf, quad, pi, log, sqrt, inf
mp.dps = 30
out = {}

A = 2*log(2) - mpf(5)/6
out['A_avg_constant'] = dict(closed='2 log 2 - 5/6', value=str(A),
                             quad=str(3*quad(lambda r: log(1+r)*r*r, [0,1])))
lim_J0 = 4*pi*quad(lambda r: r*r/(r**4+1), [0,1,10,inf])
out['limit_J0'] = dict(quad=str(lim_J0), closed_sqrt2_pi2=str(sqrt(2)*pi**2))
lim_J1logR = 4*pi*quad(lambda r: abs(A - log(1+r))*r*r/(r**4+1),
                       [0, mp.e**A - 1, 1, 10, inf])
out['limit_J1_times_logR'] = str(lim_J1logR)

rows = []
for k in (2, 4, 8, 16, 32, 64):
    R = mpf(10)**k; lg = log(R)
    h = lambda r: max(mpf(0), 1 - log(1+r)/lg)
    c = 1 - A/lg
    # log-spaced panels out to R-1, then the vanishing tail
    br = [mpf(0)] + [mpf(10)**mpf(j) for j in range(0, k+1)] + [R-1, inf]
    br = sorted(set(br))
    J0 = 4*pi*quad(lambda r: h(r)/(r**4+1)*r*r, br)
    J1 = 4*pi*quad(lambda r: abs(h(r)-c)/(r**4+1)*r*r, br)
    rows.append(dict(R=f'1e{k}', avg_B1=str(c), J0=str(J0), J0_logR=str(J0*lg),
                     J1=str(J1), J1_logR=str(J1*lg)))
out['table_logspaced_panels'] = rows
print(json.dumps(out, indent=1))

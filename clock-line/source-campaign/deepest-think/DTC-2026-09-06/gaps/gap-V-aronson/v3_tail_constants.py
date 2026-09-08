"""
v3_tail_constants.py -- gap-V-aronson.

The two explicit tail bounds for  P(|Z_tau| >= d),  Z = (backward stochastic
characteristic) - (backward deterministic characteristic), in dimension n=5.

Both start from  dZ_s = -A_s Z_s ds + sqrt(2 nu) dW_s,  Z_0 = 0,  ||A_s|| <= Gamma,
c := Gamma tau.

(B1) pathwise Gronwall + reflection principle
        |Z_tau| <= e^{c} sqrt(2 nu) sup_{s<=tau} |W_s|
        P(sup_{s<=tau}|W_s| >= R) <= 2n exp(-R^2/(2 n tau))
     =>  P(|Z_tau| >= d) <= 2n exp( - d^2 e^{-2c} / (4 n nu tau) ).

(B2) adapted propagator + exponential supermartingale + eps-net on S^{n-1}
        Z_tau = Phi_tau N_tau,  N_s = sqrt(2 nu) int_0^s Phi_r^{-1} dW_r (adapted),
        ||Phi_tau|| <= e^c, ||Phi_r^{-1}|| <= e^{Gamma r},
        <theta.N>_tau <= V0 := nu tau (e^{2c}-1)/c   for every unit theta,
        N(eps) <= (1+2/eps)^n
     =>  P(|Z_tau| >= d) <= (1+2/eps)^n exp( - (1-eps^2/2)^2 c d^2 / (2 e^{2c}(e^{2c}-1) nu tau) ).

Both are reported; the lemma uses the minimum.  q := d^2/(nu tau).
As c -> 0, (B2) exponent -> q/4, the exact heat-kernel exponent (checked below).
"""
import json, numpy as np

HERE = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/gap-V-aronson"
n = 5
OUT = {}

def B1(c, q, n=5):
    return 2*n*np.exp(-q*np.exp(-2*c)/(4*n))

_ES = np.linspace(1e-3, 1.4, 4001)
def B2(c, q, n=5, eps=None):
    k = c/(2*np.exp(2*c)*(np.exp(2*c)-1))
    if eps is not None:
        return (1+2/eps)**n*np.exp(-(1-eps*eps/2)**2*k*q)
    logv = n*np.log(1+2/_ES) - (1-_ES*_ES/2)**2*k*q
    i = int(np.argmin(logv))
    return float(np.exp(logv[i])), float(_ES[i])

def BOUND(c, q, n=5):
    b1 = B1(c, q, n); b2, e2 = B2(c, q, n)
    return min(b1, b2), b1, b2, e2

# --- sanity: c->0 limit of the B2 rate is 1/4 ---
lim = [(c, c/(2*np.exp(2*c)*(np.exp(2*c)-1))) for c in (1e-4, 1e-3, 1e-2, 0.1)]
print("B2 rate coefficient c/(2 e^{2c}(e^{2c}-1)) as c->0 (exact limit 1/4 = 0.25):")
for c, v in lim: print("   c=%.0e -> %.8f" % (c, v))
OUT['B2_rate_limit'] = lim
# exact 1D heat-kernel comparison at c=0: P(|N(0,2 nu tau)| >= d) = erfc(d/sqrt(4 nu tau))
from math import erfc, sqrt
print("\nsanity, c=0, n=1: true tail erfc(sqrt(q)/2) vs B2 exponent exp(-q/4):")
for q in (4, 9, 16, 25, 36):
    print("   q=%2d   erfc=%.4e   exp(-q/4)=%.4e   ratio=%.4f" % (q, erfc(sqrt(q)/2), np.exp(-q/4), erfc(sqrt(q)/2)/np.exp(-q/4)))
OUT['c0_sanity'] = [dict(q=q, erfc=erfc(sqrt(q)/2), exp=float(np.exp(-q/4))) for q in (4,9,16,25,36)]

# --- the table the NOTE quotes ---
print("\nn=5.  P(|Z_tau| >= d) bounds.   q = d^2/(nu tau)")
print(" %6s %8s | %12s | %12s %6s | %12s" % ("c", "q", "B1", "B2", "eps*", "min"))
rows = []
for c in (0.05, 0.2, 0.4, 0.811, 1.2, 2.0):
    for q in (10, 25, 50, 100, 200, 400):
        b, b1, b2, e2 = BOUND(c, q)
        rows.append(dict(c=c, q=q, B1=b1, B2=b2, eps=e2, best=b))
        print(" %6.3f %8.0f | %12.4e | %12.4e %6.3f | %12.4e" % (c, q, b1, b2, e2, b))
OUT['table'] = rows

# crossover q where B2 beats B1
print("\ncrossover q (B2 < B1):")
cross = {}
for c in (0.05, 0.2, 0.4, 0.811, 1.2, 2.0):
    qq = None
    for q in np.arange(1, 4000, 2.0):
        b1 = B1(c, q); b2, _ = B2(c, q)
        if b2 < b1: qq = float(q); break
    cross[str(c)] = qq
    print("   c=%.3f  q* = %s" % (c, qq))
OUT['crossover_q'] = cross

# --- the inversion the application needs: smallest q with BOUND <= target ---
def q_needed(c, target, n=5):
    lo, hi = 1.0, 1e7
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if BOUND(c, mid, n)[0] <= target: hi = mid
        else: lo = mid
    return hi
print("\nq needed for the tail bound to fall below a target:")
tg = [1e-2, 1e-3, 1e-4]
qn = {}
for c in (0.2, 0.4, 0.811, 1.2):
    qn[str(c)] = {str(t): q_needed(c, t) for t in tg}
    print("   c=%.3f : " % c + "  ".join("target %.0e -> q>=%.1f" % (t, qn[str(c)][str(t)]) for t in tg))
OUT['q_needed'] = qn

with open(HERE+"/v3_results.json","w") as fh: json.dump(OUT, fh, indent=1, default=float)
print("\nwrote v3_results.json")

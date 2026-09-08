"""
r2_sharp_explicit_form.py -- REFUTER2 of corollary-blowup-rate.
Claim under attack: the attempt's form (c) is "THE explicit form" of (b).
It is not.  (b) is  M >= f(M)  with f(M) = c1/(s(1+log_+(K M^(1/5)))) STRICTLY DECREASING in M,
so (b) is EQUIVALENT to  M >= M_fix,  the unique root of M = f(M).  M_fix is the sharp
explicit bound; the attempt's m = c1/(s(1+log_+(K s^(-1/5)))) is a weakening of it.
Blocks:
 (1) monotonicity of f, hence the equivalence, symbolically.
 (2) M_fix by bisection at high precision; ratio M_fix/m over a grid: how much (c) throws away.
 (3) the one-step sharpening (c') that recovers most of it, and its validity check.
 (4) asymptotics: both (c) and (c') have the same 5c1 coefficient, so the headline is safe.
"""
import numpy as np, sympy as sp
from mpmath import mp, mpf, log as mlog, findroot
mp.dps = 40

print("="*78); print("(1) f is strictly decreasing => (b) <=> M >= M_fix"); print("="*78)
M, s, c1, K = sp.symbols('M s c1 K', positive=True)
f = c1/(s*(1+sp.log(K*M**sp.Rational(1,5))))      # log_+ active branch
df = sp.simplify(sp.diff(f, M))
print("  d/dM [ c1/(s(1+log(K M^(1/5)))) ] =", df)
print("  sign: numerator -c1 < 0, denominator 5 s M (1+log(K M^(1/5)))^2 > 0 on the active branch")
print("  => f strictly decreasing; M - f(M) strictly increasing; the solution set of M>=f(M)")
print("     is exactly [M_fix, inf).  So M_fix is the SHARP explicit form of (b).")

def lp(x):
    x = mpf(x)
    return mlog(x) if x > 1 else mpf(0)

def Mfix(Kv, sv, cv):
    g = lambda m: m - cv/(sv*(1+lp(Kv*m**mpf(0.2))))
    m0 = cv/(sv*(1+lp(Kv/sv**mpf(0.2))))     # bracket around the attempt's m (ratio is O(1))
    lo, hi = m0/mpf('1e6'), m0*mpf('1e6')
    for _ in range(2000):
        mid = (lo*hi)**mpf('0.5')
        if g(mid) < 0: lo = mid
        else: hi = mid
        if hi/lo < 1+mpf('1e-30'): break
    return (lo*hi)**mpf('0.5')

print(); print("="*78)
print("(2) how lossy is the attempt's (c)?   ratio M_fix / m,  m = c1/(s(1+log_+(K s^(-1/5))))")
print("="*78)
print("      K        s          c1        m            M_fix        M_fix/m")
worst = 0.0
for Kv in [mpf('1e-4'), mpf(1), mpf('1e4'), mpf('1e12')]:
    for sv in [mpf('1e-2'), mpf('1e-6'), mpf('1e-12'), mpf('1e-30')]:
        for cv in [mpf('0.9'), mpf('0.1'), mpf('1e-4')]:
            mm = cv/(sv*(1+lp(Kv/sv**mpf(0.2))))
            mf = Mfix(Kv, sv, cv)
            r = float(mf/mm)
            worst = max(worst, r)
            if cv == mpf('0.1'):
                print(f"   {float(Kv):8.0e} {float(sv):9.0e} {float(cv):9.0e}  "
                      f"{float(mm):11.4e}  {float(mf):11.4e}   {r:8.5f}")
print(f"   max M_fix/m over the whole grid = {worst:.6f}")
print("   (>1 everywhere: (c) is strictly weaker than the sharp form, but the loss is a")
print("    bounded factor, never an order.  So (c) is CORRECT and LOSSY, not wrong.)")

print(); print("="*78)
print("(3) the one-step sharpening (c'): put M >= c1/s into the log instead of 1/s")
print("     (c')  M >= c1/( s (1 + [log K + (1/5)log(1/s) + (1/5)log c1]_+) )")
print("    valid by the same bootstrap because c1<1 makes the bracket SMALLER, so we must")
print("    check the bootstrap still closes.  Direct numerical check against M_fix:")
print("="*78)
print("      K        s          c1        (c)          (c')         M_fix     (c')<=M_fix?")
ok = True
for Kv in [mpf('1e-4'), mpf(1), mpf('1e4'), mpf('1e12')]:
    for sv in [mpf('1e-2'), mpf('1e-6'), mpf('1e-12'), mpf('1e-30')]:
        for cv in [mpf('0.9'), mpf('0.1'), mpf('1e-4')]:
            mm  = cv/(sv*(1+lp(Kv/sv**mpf(0.2))))
            arg = Kv*(cv/sv)**mpf(0.2)
            mmp = cv/(sv*(1+lp(arg)))
            mf  = Mfix(Kv, sv, cv)
            good = mmp <= mf*(1+mpf('1e-25'))
            ok &= bool(good)
            if cv == mpf('0.1'):
                print(f"   {float(Kv):8.0e} {float(sv):9.0e} {float(cv):9.0e}  "
                      f"{float(mm):11.4e}  {float(mmp):11.4e}  {float(mf):11.4e}   {bool(good)}")
print("   (c') <= M_fix on every grid point:", ok, "  <-- (c') is a valid lower bound")

print(); print("="*78)
print("(4) asymptotic coefficient: is it 5c1 for (c), (c') and M_fix alike?")
print("     effective coefficient  kappa(s) := M * s * log(1/s) / c1   ->  5 ?")
print("="*78)
Kv, cv = mpf(1), mpf('0.1')
print("        s            kappa[(c)]   kappa[(c')]  kappa[M_fix]")
for e in [6, 12, 30, 100, 1000, 10000, 100000]:
    sv = mpf(10)**(-e)
    L  = mlog(1/sv)
    mm  = cv/(sv*(1+lp(Kv/sv**mpf(0.2))))
    mmp = cv/(sv*(1+lp(Kv*(cv/sv)**mpf(0.2))))
    mf  = Mfix(Kv, sv, cv)
    print(f"   1e-{e:<6d}   {float(mm*sv*L/cv):10.5f}   {float(mmp*sv*L/cv):10.5f}   "
          f"{float(mf*sv*L/cv):10.5f}")
print("   all three -> 5, so the headline '(5c1+o(1))/((T-t)log(1/(T-t)))' is right for each;")
print("   but note how slowly: the approach is O(1/log(1/s)).")

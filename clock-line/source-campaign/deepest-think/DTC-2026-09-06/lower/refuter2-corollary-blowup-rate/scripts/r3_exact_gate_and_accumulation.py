"""
r3_exact_gate_and_accumulation.py -- REFUTER2 of corollary-blowup-rate.

(1) REPLACE the attempt's blind 5-parameter sampler by an EXACT decision procedure.
    (b) is  M >= f(M) with f decreasing, so (b) <=> M >= M_fix.  Hence a candidate
    explicit form  M >= m_p  is TRUE iff  m_p <= M_fix everywhere, and FALSE iff
    m_p > M_fix somewhere.  After the reduction  M = mu/s,  kappa = K s^(-1/5),
    everything depends on (kappa, c1) alone for M_fix, and on (kappa, c1, s) for m_p.
(2) Show the exponent-mutant family p < 1/5 is FALSE for EVERY p, exhibit witnesses,
    and compute the |log10 s| at which the witness first appears -- i.e. exactly how far
    outside the attempt's sampling box (s in [1e-12,1e4]) the counterexamples live.
(3) Seed distribution for the uniform sampler on the two mutants nearest the truth.
(4) The Section-4 accumulation constant by high-precision quadrature, the effective
    loglog coefficient, and the (T-t) at which it first reaches 4.5 out of 5.
"""
import numpy as np
from mpmath import mp, mpf, log as mlog, quad, exp as mexp
mp.dps = 60

def lp(x):
    x = mpf(x)
    return mlog(x) if x > 1 else mpf(0)

def mu_fix(kappa, c1):
    """root of mu*(1+log_+(kappa*mu^(1/5))) = c1 ; strictly increasing LHS in mu."""
    g = lambda mu: mu*(1+lp(kappa*mu**mpf(0.2))) - c1
    lo, hi = c1/mpf('1e60'), c1
    for _ in range(400):
        mid = (lo*hi)**mpf('0.5')
        if g(mid) < 0: lo = mid
        else: hi = mid
        if hi/lo < 1+mpf('1e-40'): break
    return (lo*hi)**mpf('0.5')

print("="*78)
print("(1) sanity: the exact procedure reproduces the truth of (c) (p = 1/5)")
print("="*78)
print("     kappa       c1        mu_fix      mu_p(p=1/5)   mu_p<=mu_fix (c) true?")
allok = True
for kexp in range(-6, 25, 3):
    kappa = mpf(10)**kexp
    for c1 in [mpf('0.999'), mpf('0.5'), mpf('1e-3'), mpf('1e-6')]:
        mf = mu_fix(kappa, c1)
        mp_ = c1/(1+lp(kappa))              # p = 1/5 : s^(1/5-p) = 1
        ok = mp_ <= mf*(1+mpf('1e-30'))
        allok &= bool(ok)
        if c1 == mpf('0.5'):
            print(f"   {float(kappa):9.0e} {float(c1):9.0e}  {float(mf):11.4e}  "
                  f"{float(mp_):11.4e}   {bool(ok)}")
print("   (c) verified TRUE on every grid point of the exact procedure:", allok)

print()
print("="*78)
print("(2) the exponent mutants p < 1/5 are FALSE FOR EVERY p -- explicit witnesses")
print("    Take K = s^p exactly.  Then K s^(-p) = 1, so log_+ = 0 and m_p = cc/s (its")
print("    largest possible value), while kappa = K s^(-1/5) = s^(p-1/5) -> infinity.")
print("    A violation exists as soon as mu_fix(kappa,c1) < c1, i.e. as soon as kappa > 1.")
print("="*78)
print("     p        1/5-p     log10(1/s) needed   kappa      mu_fix       m_p*s=c1   violation")
for p in [mpf('0.1'), mpf('0.19'), mpf('0.199'), mpf('0.1999')]:
    d = mpf('0.2') - p
    # choose kappa = 1e4  =>  log10(1/s) = 4/d
    L10 = mpf(4)/d
    s   = mpf(10)**(-L10)
    kappa = s**(p - mpf('0.2'))
    c1 = mpf('0.5')
    mf = mu_fix(kappa, c1)
    viol = mf < c1
    print(f"   {float(p):7.4f}  {float(d):8.1e}   {float(L10):14.4g}   {float(kappa):8.1e}  "
          f"{float(mf):10.4e}   {float(c1):8.4f}   {bool(viol)}")
print()
print("    => the counterexample to the mutant p exists only for  log10(1/s) >~ 4/(1/5-p).")
print("       The attempt's sampler draws s in [1e-12, 1e4], i.e. log10(1/s) <= 12.")
print("       Smallest detectable defect:  1/5-p >= 4/12 = 0.333  -- i.e. NO mutant with")
print("       p in (-0.133, 0.2) can be seen at K=s^p.  The sampler also draws K from")
print("       E0^(2/5)/nu with E0 in [1e-8,1e8], nu in [1e-8,1e4]: log10 K in [-7.2, 11.2].")
print("       Requiring K = s^p with log10(1/s)=4/d gives log10 K = -4p/d, e.g. -796 at")
print("       p=0.199 -- roughly 800 decades outside the box.")
print()
print("    Concrete witness for p = 0.199, written back in the original variables:")
p = mpf('0.199'); d = mpf('0.2')-p; L10 = mpf(4)/d
s = mpf(10)**(-L10); K = s**p; kappa = K*s**(-mpf('0.2')); c1 = mpf('0.5')
mf = mu_fix(kappa, c1)
Mw = mf/s*mpf('1.000000001')          # any M in [M_fix, c1/s) violates the mutant
print(f"      s = T-t = 1e-{float(L10):.0f} ,  K = E0^(2/5)/nu = 1e-{float(-mlog(K)/mlog(10)):.0f} ,"
      f"  c1 = 0.5")
print(f"      M = {float(mf):.6f}/s * (1+1e-9)")
print(f"      (b) holds:   M*s*(1+log_+(K M^(1/5))) = "
      f"{float(Mw*s*(1+lp(K*Mw**mpf(0.2)))):.9f}  >= c1 = 0.5   ->",
      bool(Mw*s*(1+lp(K*Mw**mpf(0.2))) >= c1))
print(f"      mutant m_p = c1/(s(1+log_+(K s^-p))) = c1/s = {float(c1):.4f}/s ;  "
      f"M*s = {float(Mw*s):.6f} < 0.5  -> mutant VIOLATED:", bool(Mw*s < c1))

print()
print("="*78)
print("(3) uniform sampler, seed distribution on the two mutants nearest the truth")
print("="*78)
def logp(x): return np.maximum(np.log(x), 0.0)
N = 2_000_000
for pm in [0.19, 0.199]:
    cnt = []
    for seed in range(101, 121):
        rng = np.random.default_rng(seed)
        E0 = 10.0**rng.uniform(-8, 8, N); nu = 10.0**rng.uniform(-8, 4, N)
        sN = 10.0**rng.uniform(-12, 4, N); c1N = 10.0**rng.uniform(-6, 0.6, N)
        MN = 10.0**rng.uniform(-12, 14, N)
        b = MN*sN*(1.0+logp(E0**0.4*MN**0.2/nu)) >= c1N
        cc = np.minimum(c1N, 1.0)
        ch = MN >= cc/(sN*(1.0+logp(E0**0.4/(nu*sN**pm))))
        cnt.append(int((b & (~ch)).sum()))
    cnt = np.array(cnt)
    print(f"   mutant p={pm}:  N={N}, 20 seeds -> min={cnt.min()} median={np.median(cnt):.1f} "
          f"max={cnt.max()}  zeros={int((cnt==0).sum())}/20")
print("   Both mutants are FALSE (block 2 gives explicit witnesses); the uniform gate")
print("   returns PASS for p=0.199 on every seed.  FL-043: same verdict, different truth.")

print()
print("="*78)
print("(4) Section 4 accumulation: high-precision quadrature vs the claimed closed form")
print("    I(V) = int_{s=e^-V}^{1} c/(r(1+log K+(1/5)log(1/r))) dr,  substitute r=e^-v:")
print("    I(V) = int_0^V c dv/(1+log K+v/5) = 5c[log(1+log K+V/5) - log(1+log K)]")
print("="*78)
c = mpf('0.1'); K = mpf(1)
print("        V=log(1/(T-t))     quadrature        closed form       claimed 5c*loglog")
print("                                                                 = 5c*log V")
for V in [mpf(10), mpf(100), mpf('1e3'), mpf('1e6'), mpf('1e8')]:
    f = lambda v: c/(1+mlog(K)+v/5)
    Q = quad(f, [0, V])
    C = 5*c*(mlog(1+mlog(K)+V/5) - mlog(1+mlog(K)))
    L = 5*c*mlog(V)
    print(f"   {float(V):16.4e}   {float(Q):15.9f}   {float(C):15.9f}   {float(L):15.9f}")
print()
print("    effective loglog coefficient  kappa_acc(V) := I(V)/(c*loglog(1/(T-t)))")
print("                                                = I(V)/(c*log V)   -> 5")
print("        T-t            V           kappa_acc     (headline says 5)")
for e in [6, 12, 100, 1000, 1e4, 1e5, 1e6, 1e7, 1e9]:
    V = mpf(e)*mlog(10)
    I = 5*c*(mlog(1+mlog(K)+V/5) - mlog(1+mlog(K)))
    print(f"   1e-{float(e):<10.0f}  {float(V):10.3e}   {float(I/(c*mlog(V))):9.5f}")
print()
print("    Solving kappa_acc(V) = 4.5 :")
lo, hi = mpf(10), mpf('1e60')
for _ in range(400):
    mid = (lo*hi)**mpf('0.5')
    I = 5*c*(mlog(1+mlog(K)+mid/5) - mlog(1+mlog(K)))
    if I/(c*mlog(mid)) < mpf('4.5'): lo = mid
    else: hi = mid
V45 = (lo*hi)**mpf('0.5')
print(f"      V = {float(V45):.4e}   i.e.  T-t = e^-V = 1e-{float(V45/mlog(10)):.4e}")
print("    So the '5c_1 loglog' coefficient is never within 10% of 5 at any (T-t) that")
print("    could be written down.  The statement is asymptotically correct; the number 5")
print("    is not an approximation to anything reachable.")

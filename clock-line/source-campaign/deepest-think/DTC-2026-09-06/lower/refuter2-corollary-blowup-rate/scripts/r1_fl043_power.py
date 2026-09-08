"""
r1_fl043_power.py -- REFUTER2 of corollary-blowup-rate.
FL-043: could the attempt's block (B) gate return the same verdict for every outcome?
Three probes the attempt (and the first refuter) did not run:
  P1  seed fragility of the attempt's OWN control (the "min(c1,1) guard is load-bearing"
      demonstration, reported as "11 violations").
  P2  detection power as a continuous function of perturbation size (exponent 1/5 -> 1/5+eps),
      i.e. the smallest defect the gate can see at N = 4e6.
  P3  does the gate distinguish TRUE from TRUE-BUT-WEAKER?  If not, it certifies the
      direction of (c) but says nothing about its constant.
  P4  tightness of the sampled ensemble: how many draws are anywhere near the boundary.
Everything printed is produced by this script.
"""
import numpy as np

def logp(x): return np.maximum(np.log(x), 0.0)

def draw(seed, N):
    rng = np.random.default_rng(seed)
    E0 = 10.0**rng.uniform(-8, 8, N)
    nu = 10.0**rng.uniform(-8, 4, N)
    s  = 10.0**rng.uniform(-12, 4, N)
    c1 = 10.0**rng.uniform(-6, 0.6, N)
    M  = 10.0**rng.uniform(-12, 14, N)
    return E0, nu, s, c1, M

def bhold(E0, nu, s, c1, M):
    return M*s*(1.0 + logp(E0**0.4*M**0.2/nu)) >= c1

N = 4_000_000
print("="*78)
print("P1  seed fragility of the attempt's own control")
print("    control = form (c) WITHOUT the min(c1,1) guard; attempt reports 11 violations")
print("    at seed 20260906.  Re-run over 40 seeds, same N =", N)
print("="*78)
counts = []
for k, seed in enumerate([20260906] + list(range(1, 40))):
    E0, nu, s, c1, M = draw(seed, N)
    b = bhold(E0, nu, s, c1, M)
    Re_star = E0**0.4/(nu*s**0.2)
    c_nomin = M >= c1/(s*(1.0 + logp(Re_star)))
    v = int((b & (~c_nomin)).sum())
    counts.append(v)
    if k < 6 or v == 0:
        print(f"   seed {seed:>9}   violations = {v}")
counts = np.array(counts)
print(f"   over {len(counts)} seeds: min={counts.min()}  median={np.median(counts):.1f} "
      f"max={counts.max()}  mean={counts.mean():.2f}")
print(f"   seeds returning ZERO violations (i.e. the control silently PASSES) = "
      f"{int((counts==0).sum())}/{len(counts)}")
print(f"   detection rate per draw = {counts.mean()/N:.3e}")

print()
print("="*78)
print("P2  power curve: exponent 1/5 -> 1/5+eps inside Re_* = E0^(2/5)/(nu s^(1/5))")
print("    (the exponent that actually carries the physics).  seed 20260906, N =", N)
print("="*78)
E0, nu, s, c1, M = draw(20260906, N)
b = bhold(E0, nu, s, c1, M)
cc = np.minimum(c1, 1.0)
print("    eps        exponent      violations   fires?")
for eps in [0.0, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.05, 1/6-0.2, 0.05, -0.05,
            -1e-2, -1e-3, -1e-4]:
    p = 0.2 + eps
    Re_m = E0**0.4/(nu*s**p)
    ch = M >= cc/(s*(1.0 + logp(Re_m)))
    v = int((b & (~ch)).sum())
    print(f"   {eps:+9.4g}   {p:9.6f}   {v:10d}   {'YES' if v else 'no'}")
print("    NOTE: eps<0 shrinks Re_*, shrinks the denominator, STRENGTHENS the claim -> must fire.")
print("          eps>0 grows Re_*, WEAKENS the claim -> a true-but-weaker statement, never fires.")

print()
print("="*78)
print("P3  TRUE vs TRUE-BUT-WEAKER: does the gate see the constant at all?")
print("="*78)
Re_star = E0**0.4/(nu*s**0.2)
for tag, num in [("(c) as stated, cc=min(c1,1)", cc),
                 ("cc/2   (weaker, still true)", cc/2),
                 ("cc/1e6 (absurdly weak, true)", cc/1e6),
                 ("2*cc   (stronger, false?)",   2*cc),
                 ("10*cc  (stronger, false?)",   10*cc)]:
    ch = M >= num/(s*(1.0 + logp(Re_star)))
    v = int((b & (~ch)).sum())
    print(f"   {tag:32s} violations = {v:9d}   {'FIRES' if v else 'PASSES'}")
print("   => every weakening passes. The gate certifies the DIRECTION of (c), not its constant.")

print()
print("="*78)
print("P4  how tight is the ensemble?  m = cc/(s(1+log+ Re_*)) is the (c) threshold.")
print("="*78)
m = cc/(s*(1.0 + logp(Re_star)))
ratio = M[b]/m[b]
print(f"   draws with (b) true                       = {b.sum()}")
for f in [1.0, 1.01, 1.1, 2.0, 10.0, 1e2, 1e4]:
    print(f"   of those, M/m <= {f:8g}                 = {int((ratio<=f).sum()):9d}"
          f"   ({100.0*(ratio<=f).mean():.4f}% )")
print(f"   median M/m = {np.median(ratio):.3e}   geometric mean = "
      f"{np.exp(np.mean(np.log(ratio))):.3e}")
print("   => the ensemble is overwhelmingly far from the boundary; the '0 violations in")
print("      4,000,000' headline is carried by a vanishing fraction of informative draws.")

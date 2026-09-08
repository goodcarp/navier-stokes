"""
r1_fl043_mutants.py -- REFUTER, DTC-2026-09-06/lower/refuter-corollary-blowup-rate
FL-043 on the attempt's block (B): could that gate return "0 violations" for EVERY outcome?
Method: re-run the identical sampler, but replace the CORRECT explicit form (c) by a family
of MUTANTS.  A gate with discriminating power must return >0 violations for wrong mutants
and 0 for the true one.  Numerics falsify, never prove.
"""
import numpy as np
rng = np.random.default_rng(20260906)          # same seed as the attempt's s1
def logp(x): return np.maximum(np.log(x), 0.0)
N = 4_000_000
E0 = 10.0**rng.uniform(-8, 8, N)
nu = 10.0**rng.uniform(-8, 4, N)
s  = 10.0**rng.uniform(-12, 4, N)
c1 = 10.0**rng.uniform(-6, 0.6, N)
M  = 10.0**rng.uniform(-12, 14, N)
cc = np.minimum(c1, 1.0)

# hypothesis (b), exactly as the attempt writes it
b = M*s*(1.0 + logp(E0**0.4*M**0.2/nu)) >= c1
print(f"samples {N};  (b) true = {int(b.sum())}  ({100*b.mean():.2f}%)")
print()
print("mutant                                        (b)&~(mutant)   verdict")
print("-"*74)

def gate(name, c_ok, expect):
    bad = b & (~c_ok)
    v = int(bad.sum())
    ok = ("PASS" if v == 0 else "FIRES")
    flag = "  <-- GATE BLIND" if (expect == "fires" and v == 0) else ""
    print(f"{name:<45} {v:>13}   {ok}{flag}")
    return v

# TRUE statement (c): must be 0
gate("(c) TRUE  exponent 1/5, +1, min(c1,1)",
     M >= cc/(s*(1.0 + logp(E0**0.4/(nu*s**0.2)))), "pass")

# M1: wrong exponent on (T-t) inside the log  (1/5 -> 1/4)
gate("M1  exponent 1/4 in Re_*",
     M >= cc/(s*(1.0 + logp(E0**0.4/(nu*s**0.25)))), "fires")
# M2: wrong exponent 1/6
gate("M2  exponent 1/6 in Re_*",
     M >= cc/(s*(1.0 + logp(E0**0.4/(nu*s**(1/6.))))), "fires")
# M3: wrong energy exponent 2/5 -> 1/5
gate("M3  E0 exponent 1/5 in Re_*",
     M >= cc/(s*(1.0 + logp(E0**0.2/(nu*s**0.2)))), "fires")
# M4: drop the +1 inside the bracket
gate("M4  bracket = log_+ Re_* (no +1)",
     M >= cc/(s*logp(E0**0.4/(nu*s**0.2))), "fires")
# M5: no min guard (the attempt's own positive control)
gate("M5  no min(c1,1) guard",
     M >= c1/(s*(1.0 + logp(E0**0.4/(nu*s**0.2)))), "fires")
# M6: constant inflated by 5 (the asymptotic constant, applied at every t)
gate("M6  numerator 5*cc instead of cc",
     M >= 5*cc/(s*(1.0 + logp(E0**0.4/(nu*s**0.2)))), "fires")
# M7: log_+ replaced by log (can go negative -> denominator <1)
gate("M7  log instead of log_+",
     M >= cc/(s*(1.0 + np.log(E0**0.4/(nu*s**0.2)))), "fires")
# M8: nu in the numerator of Re_*
gate("M8  nu multiplied not divided",
     M >= cc/(s*(1.0 + logp(E0**0.4*nu/s**0.2))), "fires")

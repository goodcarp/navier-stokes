#!/usr/bin/env python3
"""S2 - the explicit constants C1 (far Taylor remainder), C2in (inner multipole), C2collar.
All from the REPAIRED interior identity d_z[rho^l C_l] = (l+2) rho^{l-1} C_{l-1} (S1) and Parseval.

Notation: w = omega^theta on a shell, |w| <= M; G = w/sin(phi) = SUM_l g_l C_l^{3/2}(t);
Parseval: SUM_l g_l^2 N_l = INT_{-1}^{1} w(t)^2 dt <= 2 M^2,  N_l = (l+1)(l+2)/(l+3/2).

Interior (|xi|<1):  Phi(xi) = -SUM_{l>=1} ((l+2) g_l/(2l+3)) |xi|^{l-1} C_{l-1}^{3/2}(t)
Exterior (|xi|>1):  Phi(xi) = +SUM_{l>=0} ((l+1) g_l/(2l+3)) |xi|^{-(l+4)} C_{l+1}^{3/2}(t)
|C_m^{3/2}(t)| <= C_m^{3/2}(1) = (m+1)(m+2)/2 on [-1,1].
"""
import numpy as np, math, hashlib, json
from scipy.integrate import quad
Nl   = lambda l: (l+1)*(l+2)/(l+1.5)
Cm1  = lambda m: (m+1)*(m+2)/2.0          # C_m^{3/2}(1)

def Q(u, odd_only, LMAX=4000):
    ls = np.arange(3, LMAX, 2) if odd_only else np.arange(2, LMAX)
    c  = ((ls+2)/(2*ls+3))**2 * Cm1(ls-1)**2 / Nl(ls)
    with np.errstate(under='ignore'):
        return math.sqrt(float(np.sum(c*u**(2.0*(ls-1)))))
def S(v, odd_only, LMAX=4000):
    ls = np.arange(1, LMAX, 2) if odd_only else np.arange(0, LMAX)
    c  = ((ls+1)/(2*ls+3))**2 * Cm1(ls+1)**2 / Nl(ls)
    with np.errstate(under='ignore'):
        return math.sqrt(float(np.sum(c*v**(2.0*(ls+4)))))

out = {}
for tag, odd in (("general (any bounded omega^theta)", False), ("z-odd (near-cap stack)", True)):
    C1,  e1 = quad(lambda u: Q(u, odd)/u, 0, 0.5, limit=200)
    C2i, e2 = quad(lambda v: S(v, odd)/v, 0, 0.5, limit=200)
    C1 *= math.sqrt(2); C2i *= math.sqrt(2); e1 *= math.sqrt(2); e2 *= math.sqrt(2)
    out[tag] = dict(C1=C1, C1_quaderr=e1, C2in=C2i, C2in_quaderr=e2)
    print(f"{tag}:")
    print(f"   C1   (far Taylor remainder, |a_far(x)-a_far(0)| <= C1 M)      = {C1:.6f}   (quad err {e1:.1e})")
    print(f"   C2in (inner multipole,       |a_inner(x)|        <= C2in M)   = {C2i:.6f}   (quad err {e2:.1e})")
    print(f"   leading small-u behaviour: Q(u)/u -> {Q(1e-4,odd)/1e-4:.6f},  S(v)/v^{4 if not odd else 5} -> {S(1e-2,odd)/1e-2**(4 if not odd else 5):.6f}")

# --- collar constant (rho/2 < rho' < 2 rho), proved by two crude but explicit pieces -------------
# |a_collar| <= (3/(8 pi^2)) INT_A |omega^theta| /(r' |x-x'|^4) dx',  A = {rho/2<|x'|<2rho}, rho=1.
# A1 = A cap {r' >= r/2}: 1/r' <= 2/r; bathtub/rearrangement: INT_{A1}|x-x'|^{-4}dx' <= |S^4| R_A,
#      |A| = (|S^4|/5)(2^5 - 2^-5) rho^5  =>  R_A = (2^5-2^-5)^{1/5} rho.
# A2 = A cap {r' <  r/2}: |x-x'|^2 >= (z-z')^2 + r^2/4 and (1/r')dx' = 2 pi^2 r'^2 dr' dz';
#      INT_0^{r/2} r'^2 dr' INT_R dz'/((z-z')^2+r^2/4)^2 = (r^3/24)(4 pi/r^3) = pi/6.
S4 = 8*math.pi**2/3
R_A = (2**5 - 2**-5)**0.2
A1 = (3/(8*math.pi**2)) * 2 * S4 * R_A          # coefficient of 1/s
A2 = (3/(8*math.pi**2)) * 2*math.pi**2 * math.pi/6
print(f"\ncollar bound  |a_collar| <= ({A1:.6f}/s + {A2:.6f}) M,   s = sin(phi) at the evaluation point")
print(f"   R_A = (2^5-2^-5)^(1/5) = {R_A:.6f} rho ;  A2 = pi/8 = {math.pi/8:.6f}")
for phid in (10, 30, 45, 90):
    s = math.sin(math.radians(phid)); print(f"   phi={phid:3d} deg  s={s:.4f}  ->  bound {A1/s + A2:8.3f} M")
out['collar'] = dict(coef_over_s=A1, const=A2, R_A=R_A)

# --- assembled total ------------------------------------------------------------------------------
print("\nTOTAL:  | a(x) - (M/2) * (e-folds of near-cap weight outside 2|x|) |  <=  C_tot(phi) M")
for tag, odd in (("general", False), ("z-odd", True)):
    C1 = out["general (any bounded omega^theta)" if not odd else "z-odd (near-cap stack)"]['C1']
    C2 = out["general (any bounded omega^theta)" if not odd else "z-odd (near-cap stack)"]['C2in']
    for phid in (10, 45, 90):
        s = math.sin(math.radians(phid))
        print(f"   {tag:8s} phi={phid:3d} deg :  C1={C1:.4f} + collar={A1/s+A2:.3f} + C2in={C2:.5f}  =  {C1+A1/s+A2+C2:.3f}")
json.dump(out, open('s2_constants.json','w'), indent=1)
print('\nSCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())

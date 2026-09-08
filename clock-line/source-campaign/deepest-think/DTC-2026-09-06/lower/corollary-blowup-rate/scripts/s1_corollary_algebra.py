"""
s1_corollary_algebra.py  --  DTC-2026-09-06 / lower / corollary-blowup-rate
Everything printed here is produced by this script. Numerics falsify, never prove;
the symbolic blocks (A),(C),(D) are identities, the numeric blocks (B),(E) are searches
for counterexamples to steps I claim are theorems.
"""
import sympy as sp
import numpy as np

print("="*78); print("(A) NS scaling: invariance of M*T and Re_E"); print("="*78)
lam, E, M, nu, T = sp.symbols('lambda E M nu T', positive=True)
# u_lam(x,t) = lam u(lam x, lam^2 t), nu fixed.
# ||u_lam||_2^2 = lam^2 * lam^-3 * E = E/lam ;  ||curl u_lam||_inf = lam^2 M ; times scale as lam^-2
E_l = E/lam; M_l = lam**2*M; T_l = T/lam**2
ell   = E**sp.Rational(1,5)*M**sp.Rational(-2,5)
ell_l = E_l**sp.Rational(1,5)*M_l**sp.Rational(-2,5)
Re    = M*ell**2/nu
Re_l  = M_l*ell_l**2/nu
print("  ell_lam/ell           =", sp.simplify(ell_l/ell), "   (expect 1/lambda)")
print("  Re_E(lam)/Re_E        =", sp.simplify(Re_l/Re),   "   (expect 1)")
print("  Re_E as E,M,nu        =", sp.simplify(Re), "  (expect E**(2/5)*M**(1/5)/nu)")
print("  (M*T)_lam/(M*T)       =", sp.simplify((M_l*T_l)/(M*T)), "   (expect 1)")

print(); print("="*78)
print("(B) the implicit -> explicit bootstrap:  (b) implies (c).")
print("    (b) M*s*(1+log+(E0^(2/5)M^(1/5)/nu)) >= c1")
print("    (c) M >= cc/(s*(1+log+(E0^(2/5)/(nu*s^(1/5)))))   with cc=min(c1,1)")
print("    Search for a counterexample: (E0,nu,s,c1,M) with (b) true and (c) false.")
print("="*78)
rng = np.random.default_rng(20260906)
def logp(x): return np.maximum(np.log(x), 0.0)
N = 4_000_000
# NB: distinct names from the sympy symbols above (an earlier version of this
# script shadowed nu and M with numpy arrays, which silently turned block (C)
# into a four-million-element sympification).
E0n = 10.0**rng.uniform(-8, 8, N)
nun = 10.0**rng.uniform(-8, 4, N)
sn  = 10.0**rng.uniform(-12, 4, N)
c1n = 10.0**rng.uniform(-6, 0.6, N)     # includes c1 > 1 to exercise the min(c1,1)
Mn  = 10.0**rng.uniform(-12, 14, N)
cc = np.minimum(c1n, 1.0)
b_holds = Mn*sn*(1.0 + logp(E0n**0.4*Mn**0.2/nun)) >= c1n
Re_star = E0n**0.4/(nun*sn**0.2)
c_holds = Mn >= cc/(sn*(1.0 + logp(Re_star)))
bad = b_holds & (~c_holds)
print(f"  samples                 = {N}")
print(f"  (b) true                = {b_holds.sum()}")
print(f"  (b) true and (c) false  = {int(bad.sum())}    <-- must be 0")
if bad.sum():
    i = np.argmax(bad); print("   witness:", E0n[i], nun[i], sn[i], c1n[i], Mn[i])

print(); print("  Same search with cc replaced by c1 (i.e. WITHOUT the min(c1,1) guard):")
c_holds_nomin = Mn >= c1n/(sn*(1.0 + logp(Re_star)))
bad2 = b_holds & (~c_holds_nomin)
print(f"  (b) true and (c-no-min) false = {int(bad2.sum())}   <-- expected > 0: the guard is needed")

print(); print("="*78)
print("(C) the classical comparator: Leray/Giga  ||u(t)||_inf >= c nu^(1/2)(T-t)^(-1/2)")
print("    fed through  ||u||_inf <= C E^(1/5) M^(3/5)   (Astra pass 8, eq. (4))")
print("="*78)
a, b, g = sp.symbols('a b g', real=True)   # M >= const * nu^a * E^b * (T-t)^g
sdt = sp.Symbol('s', positive=True)
# E^(1/5) M^(3/5) >= c nu^(1/2) s^(-1/2)  =>  M >= c' nu^(5/6) E^(-1/3) s^(-5/6)
sol = sp.solve(sp.Eq(sp.Rational(3,5)*a, sp.Rational(1,2)), a)[0]
print("  exponent of nu   =", sol, "  (from (3/5)a = 1/2)")
print("  exponent of E    =", sp.Rational(-1,3))
print("  exponent of (T-t)=", sp.Rational(-5,6))
# check by substitution
Cc = sp.Symbol('C', positive=True)
lhs = (E**sp.Rational(-1,3)*nu**sp.Rational(5,6)*sdt**sp.Rational(-5,6))
print("  check E^(1/5)*(that)^(3/5) =",
      sp.powsimp(sp.expand_power_base(E**sp.Rational(1,5)*lhs**sp.Rational(3,5), force=True), force=True),
      "  (expect nu^(1/2)/s^(1/2))")
print("  => classical route gives exponent 5/6 in (T-t); the corollary gives 1 minus a log.")

print(); print("="*78)
print("(D) accumulation of vorticity implied by the corollary (compare BKM / Ingimarson-Kukavica)")
print("="*78)
r, K, cst = sp.symbols('r K c', positive=True)
# integrand c/(r*(1+log+(K r^(-1/5)))) in r=T-s, on the branch where the log+ is active.
# Verify the claimed antiderivative by differentiation (no sp.integrate call).
integrand = cst/(r*(1+sp.log(K)-sp.Rational(1,5)*sp.log(r)))
F = -5*cst*sp.log(1+sp.log(K)-sp.Rational(1,5)*sp.log(r))
print("  claimed antiderivative F(r) = -5c*log(1 + log K - (1/5)log r)")
print("  F'(r) - integrand           =", sp.simplify(sp.diff(F, r) - integrand), "  (expect 0)")
print("  => int_?^t ||omega|| ds >= 5*c*log( 1 + log K + (1/5)log(1/(T-t)) ) + O(1)")
print("  i.e.  >= 5*c*loglog(1/(T-t)) + O(1).")
print("  Ingimarson-Kukavica (Euler, Thm 2.1): int_0^t ||omega|| >= (1/C) loglog(1/(T*-t)).")

print(); print("="*78)
print("(E) where the corollary overtakes the Leray-derived comparator (both constants set to 1)")
print("="*78)
print("   ours:      M >= c1/(s*(1+log+(E0^(2/5)/(nu s^(1/5)))))")
print("   classical: M >= c2*nu^(5/6) E0^(-1/3) s^(-5/6)")
print("   The crossover LOCATION depends on the two unknown constants c1,c2; what is")
print("   constant-free is the RATIO's growth exponent.  Set c1=c2=1 and tabulate:")
print()
print("     E0     nu        (T-t)      corollary      classical      ratio")
for (E0v, nuv) in [(1.0,1.0),(1.0,1e-2),(1e3,1e-3),(1e-3,1e-1)]:
    for sv in [1e0,1e-3,1e-6,1e-9,1e-12,1e-15]:
        ours = 1.0/(sv*(1.0+logp(E0v**0.4/(nuv*sv**0.2))))
        clas = nuv**(5/6)*E0v**(-1/3)*sv**(-5/6)
        print(f"   {E0v:7g} {nuv:7g}  {sv:9.0e}  {ours:13.4e}  {clas:13.4e}  {ours/clas:9.3e}")
    print()
print("   Fitted slope of log(ratio) vs log(1/s) over s in [1e-15,1e-6]:")
for (E0v, nuv) in [(1.0,1.0),(1e3,1e-3)]:
    ss = 10.0**np.arange(-6,-15.01,-0.25)
    ours = 1.0/(ss*(1.0+logp(E0v**0.4/(nuv*ss**0.2))))
    clas = nuv**(5/6)*E0v**(-1/3)*ss**(-5/6)
    sl = np.polyfit(np.log(1/ss), np.log(ours/clas), 1)[0]
    print(f"     E0={E0v:g} nu={nuv:g}:  slope = {sl:.4f}   (-> 1/6 = {1/6:.4f} as s->0, minus the log)")

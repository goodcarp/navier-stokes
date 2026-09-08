#!/usr/bin/env python3
"""R3: is the bootstrap closure CIRCULAR or does a universal c exist?
Exact where possible (Fraction/sympy), high-precision mpmath elsewhere."""
import sympy as sp
from sympy import Rational as Q
from fractions import Fraction
import mpmath as mp
mp.mp.dps=40

# increment/M0 at T=H=c/(M0*A), A=1+log_+ Re,  a=M0*H=c/A
# = C*a*[1 + (1/2)*log_+(Re*A/c) + a^2]
def Finc(c,C,Re):
    A=1+max(mp.log(Re),mp.mpf(0))
    a=mp.mpf(c)/A
    lg=mp.log(Re*A/c); lg=lg if lg>0 else mp.mpf(0)
    return C*a*(1+lg/2+a**2)

print("=== does sup over Re of the increment stay <= 1/2 for small c? ===")
Res=[mp.mpf(10)**k for k in range(-6,61)]
for C in [1,10,100,1000,10**6]:
    found=None
    k=0
    while k<80:
        c=mp.mpf(2)**(-k)
        m=max(Finc(c,C,R) for R in Res)
        if m<=mp.mpf(1)/2:
            found=(k,m); break
        k+=1
    print(f"  C={C:<8} smallest c=2^-k that works: k={found[0]:<3} c={float(2.0**-found[0]):.3e}  sup_inc={float(found[1]):.4f}")

print()
print("=== is the closure circular? symbolic check that the log REGENERATES ===")
c,A,Re=sp.symbols('c A Re',positive=True)
lhs=sp.log(sp.sqrt(Re*A/c))          # log_+(ell0/sqrt(nu H)) with H=c/(M A)
rhs=sp.Rational(1,2)*(sp.log(Re)+sp.log(A)+sp.log(1/c))
print("  log(ell0/sqrt(nu H)) - (1/2)[log Re + log A + log(1/c)] =",sp.simplify(lhs-rhs))
print("  -> the substitution T=H reproduces (1/2)log Re, which is EXACTLY the")
print("     quantity A=1+log_+Re in the denominator absorbs. Bounded, not circular.")
print("  a*(1/2)*log A <= (c/A)*(1/2)*A = c/2 using log A<=A:",
      "log A<=A for A>=1 ->", all(mp.log(x)<=x for x in [1,1.0001,2,10,1e6,1e30]))
print("  a*(1/2)*log(1/c) = (c/(2A))log(1/c) <= (c/2)log(1/c) -> 0 as c->0:",
      [ (float(cc), float(cc/2*mp.log(1/cc))) for cc in [mp.mpf('0.1'),mp.mpf('0.01'),mp.mpf('1e-6')]])

print()
print("=== GAIN over the standard bounded-velocity baseline ===")
# baseline: mild theory restart time ~ c1*nu/||u||_inf^2, ||u||_inf<=C E^{1/5}M^{3/5}=C M ell
# => T_base ~ nu/(M^2 ell^2) = 1/(M*Re).  ratio H0/T_base = c*Re/(1+log Re)
for R in [mp.mpf(1),mp.mpf(10)**3,mp.mpf(10)**6,mp.mpf(10)**12]:
    print(f"  Re={float(R):.0e}  H0/T_base = c*{float(R/(1+mp.log(R))):.4e}")
print("  (blind re-deriver's independent L^p scan gave the same ceiling: sup_{p>3} of")
print("   the mild-vorticity exponent is Re^-1, attained only as p->inf.)")

print()
print("=== the geometric-record series (corollary) ===")
for R0 in [mp.mpf(1),mp.mpf(10)**3,mp.mpf(10)**6,mp.mpf(10)**12]:
    ssum=mp.nsum(lambda j: mp.mpf(2)**(-j)/(1+max(mp.log(R0)+j*mp.log(2)/5,mp.mpf(0))),[0,mp.inf])
    print(f"  Re0={float(R0):.0e}  sum_j H_j*(M0/c) = {float(ssum):.10f}  (FINITE -> accumulation NOT excluded)")
# and the ODE majorant blows up in finite time
I=mp.quad(lambda Mv: 1/(Mv**2*(1+mp.log(Mv)/5)),[mp.e,mp.inf])
print(f"  int_e^inf dM/(M^2(1+log M/5)) = {float(I):.6f} < inf  -> the clock's own ODE majorant")
print("  blows up in finite time; the theorem CANNOT imply regularity. Self-consistent.")

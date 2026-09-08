#!/usr/bin/env python3
"""R5: last hostile checks -- (a) the package's cap-substitution makes the blind
re-deriver's 'missing monotonicity sub-lemma' unnecessary; (b) enlarging ell
upward is legal; (c) the Duhamel source constant bookkeeping; (d) the pure-HEAT
route cannot reach the log clock (so the advective kernel is load-bearing)."""
import sympy as sp
from sympy import Rational as Q
import mpmath as mp
E,M,M0,nu,tau,T,c=sp.symbols('E M M0 nu tau T c',positive=True)

# (a)/(b) cap substitution
ell=lambda EE,MM: EE**Q(1,5)*MM**Q(-2,5)
print("ell(E0,2*M0)/ell(E0,M0) =",sp.simplify(ell(E,2*M0)/ell(E,M0)),
      "=",float(sp.simplify(ell(E,2*M0)/ell(E,M0))))
print("  -> cap length is SMALLER than ell0, so log_+(ell_cap/sqrt(nu tau)) <= log_+(ell0/sqrt(nu tau)):",
      float(2**sp.Rational(-2,5))<1)
print("  and (7) needs only L>=ell_cap, which L=max(ell_cap,sqrt(nu tau)) satisfies.")
print("  => the package never evaluates the estimate at the INSTANTANEOUS M(s);")
print("     it uses the fixed cap 2*M0 throughout, so the blind re-deriver's item [11]")
print("     (monotonicity of M^2[1+log_+(ell(M)/sqrt(nu tau))+(M tau)^2] in M) is NOT")
print("     a gap in THIS package. Verify the monotonicity anyway:")
f=M**2*(1+sp.log(ell(E,M)/sp.sqrt(nu*tau))+(M*tau)**2)
print("   d/dM f / M =",sp.simplify(sp.diff(f,M)/M),
      "   -> = 8/5 + 2log_+ + 4 M^2 tau^2 > 0:",
      sp.simplify(sp.diff(f,M)/M - (Q(8,5)+2*sp.log(ell(E,M)/sp.sqrt(nu*tau))+4*M**2*tau**2))==0)

# (c) constants bookkeeping in (16)
print()
print("(16) bookkeeping: source <= int |omega| * (kernel avg of |grad u|)")
print("   <= 2M0 * C*(2M0)*[1+log_+(ell0/sqrt(nu(t-s)))+(2M0(t-s))^2]")
print("   = 4C M0^2 [ ... + 4(M0(t-s))^2 ] -> absorbed into a single C. consistent.")

# (d) pure-heat route ceiling: advection Duhamel with heat kernel costs ||u||_inf
print()
print("(d) HEAT-kernel route (BFG-style mild form): the ADVECTION Duhamel term is")
print("    int_0^T ||grad G_{nu(T-s)}||_1 ||u ox omega||_inf ds = C U M sqrt(T/nu).")
print("    With U = C M ell, closing at <= M/2 forces  M ell sqrt(T/nu) <~ 1, i.e.")
print("    T <~ nu/(M^2 ell^2) = 1/(M*Re).  So the pure-heat route CANNOT reach the")
print("    logarithmic clock -- it caps at the standard baseline. The advective kernel")
print("    (which swallows advection and leaves only the stretching pairing) is the")
print("    load-bearing new ingredient of pass 8.")
for R in [mp.mpf(10)**k for k in [0,3,6,12]]:
    print(f"      Re={float(R):.0e}:  (1/(M Re)) vs c/(M(1+log Re)) -> gain {float(R/(1+mp.log(R))):.3e} x c")

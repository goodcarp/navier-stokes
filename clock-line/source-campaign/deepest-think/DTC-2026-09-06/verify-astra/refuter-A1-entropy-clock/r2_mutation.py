#!/usr/bin/env python3
"""R2 / FL-043: can entropy_clock_gate.py's controls FAIL on a false clock?
Re-implement each assert-block of the gate, then substitute clocks the
theorem does NOT prove, and record pass/fail."""
import sympy as sp
from sympy import Rational as Q
nu,E,M,lam=sp.symbols('nu E M lam',positive=True)
Re=E**Q(2,5)*M**Q(1,5)/nu
sub={E:E/lam,M:lam**2*M}
def gate_scaling(clock):
    """The two scaling asserts the gate actually performs on the clock."""
    try:
        ok_inv = sp.simplify(Re.subs(sub,simultaneous=True)-Re)==0
        ok_dur = sp.simplify(clock.subs(sub,simultaneous=True)/clock-lam**-2)==0
        return bool(ok_inv and ok_dur)
    except Exception as e:
        return f"ERR {e}"
cands={
 "CLAIMED  c/(M(1+log Re))      ": 1/(M*(1+sp.log(Re))),
 "FALSE-A  c/M   (no log; BFG-like, STRONGER)": 1/M,
 "FALSE-B  c/(M(1+log Re)^5)    ": 1/(M*(1+sp.log(Re))**5),
 "FALSE-C  c/(M(1+log Re)^(1/2))": 1/(M*(1+sp.log(Re))**Q(1,2)),
 "FALSE-D  c/(M*Re)  (the OLD conservative clock)": 1/(M*Re),
 "FALSE-E  c*Re^7/M (absurdly long)": Re**7/M,
 "FALSE-F  c/(M*exp(Re))        ": sp.exp(-Re)/M,
 "WRONG-G  c/(M^2 ell)  (not lam^-2)": 1/(M**2*E**Q(1,5)*M**Q(-2,5)),
}
print("=== gate scaling controls: does the gate reject a false clock? ===")
for k,v in cands.items():
    print(f"{k:48s} gate-passes={gate_scaling(v)}")

print()
print("=== the gate's UNIFORM_BOOTSTRAP_SMALLNESS block, verbatim algebra ===")
K=sp.symbols('K',positive=True)
c=1/(K*(1+sp.log(K)))
lhs=(K/100)*c*2*(1+sp.log(K))
print("  (K/100)*c*2*(1+log K) =",sp.simplify(lhs),"  gate asserts ==2/100 ->",sp.simplify(lhs-Q(2,100))==0)
print("  NOTE: c cancels identically; this holds for EVERY K>0 and is independent")
print("  of the entropy bound, the kernel bound, and the Duhamel closure.")
# what the closure ACTUALLY needs:  C*a*[1 + (1/2)log_+(Re*A/c) + a^2] <= 1/2, a=c/A, A=1+log_+Re
print("  the quantity the proof must bound, F(c,C,A) = C*(c/A)*[1+(1/2)log(Re*A/c)+(c/A)^2],")
print("  is never formed or evaluated anywhere in entropy_clock_gate.py.")
# does the gate's identity survive replacing c by ANYTHING?
for repl,name in [(1/K,"c=1/K"),(sp.Integer(1),"c=1"),(K,"c=K (absurd, c must be <1)")]:
    l2=(K/100)*repl*2*(1+sp.log(K))
    print(f"    with {name:22s}: expression = {sp.simplify(l2)}  (gate would need ==2/100: {sp.simplify(l2-Q(2,100))==0})")

print()
print("=== the gate's DROPPED_DISPLACEMENT control ===")
L,mu,beta=sp.symbols('L mu beta',positive=True),sp.symbols('mu'),sp.symbols('beta')
L=sp.symbols('L',positive=True)
val=(beta**2*L**2-beta*mu).subs({L:1,mu:2,beta:1})
print("  hard-coded witness beta^2 L^2 - beta mu at (L,mu,beta)=(1,2,1) =",val,"  gate asserts ==-1")
print("  this is a single hand-picked numeric substitution; it is a genuine sign")
print("  witness for 'dropping the displacement breaks the tilt', but it tests an")
print("  algebraic identity, not any step of the NS argument.")

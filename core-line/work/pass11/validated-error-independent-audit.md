# Independent audit of the rational comparison propagator

The current `validated_error_hierarchy.py` and `validated-error-propagation.md` are accepted for their stated conditional comparison-system purpose. The path kernels, positive remainder, nonlinear closure and inherited subdivision are mathematically consistent. Two direct-constructor validation holes found during this audit were fixed by the owning agent and independently retested. No NS approximation norms or residual bounds were available or inferred.

The independent checker `check_independent_hierarchy_audit.py` now passes 161 exact-rational assertions; its record is `independent-hierarchy-audit-check.json`. The original checker also passes its updated 44 checks. The independent script imports the production functions but obtains its linear reference from ordinary matrix powers with a row-sum-norm exponential remainder, rather than the production path or complete-homogeneous recurrence.

## Linear propagation

For a nonnegative lower-triangular matrix, variation of constants assigns one term to every strictly increasing sequence of off-diagonal transitions. Intermediate diagonal evolution is the convolution of the scalar exponentials along that sequence. The implementation's path enumeration has the correct multiplicity, including direct edges and the no-edge diagonal path. Adding a zero rate integrates that convolution once and therefore produces \(J(h)=\int_0^h e^{Bt}dt\).

The coefficient recurrence updates each positive-rate geometric generating factor in ascending polynomial degree. This correctly computes the complete homogeneous polynomials, including repeated or zero rates. For a path with \(k\) edges,

\[
h_n(x_0,\ldots,x_k)\le {n+k\choose k}x^n,
\qquad \frac{{n+k\choose k}}{(n+k)!}=\frac1{k!n!},
\]

where \(x\ge\max x_i\). Thus the implementation's kernel tail \(h^k\,\mathrm{tail}_{\exp}/k!\) is correct. For the first omitted exponential term at degree \(N+1\), all subsequent term ratios are at most \(x/(N+2)\); the geometric majorant has the correct denominator. Only the largest diagonal argument appears in this tail. Large off-diagonal coefficients legitimately enter finite positive products outside it.

The independent checks include 12 nontrivial matrices of sizes one through five, repeated diagonal rates in half the cases, arbitrary positive lower entries, and the integrated kernel. Every production interval contains the tighter independently generated rational Taylor reference interval. These supplement the existing exact nilpotent and distinct-rate tests.

## Nonlinear closure and rounding

Let \(a,b\ge0\) be the outward-rounded linear and nonlinear response bounds. Writing \(A=|a|^2\), \(P=a\cdot b\), \(C=|b|^2\), \(t=1-2P\), the closure polynomial is exactly
\(f(K)=CK^2-tK+A\). For \(A>0,t>0,D=t^2-4AC\ge0\), the proposed \(K=2A/t\) gives \(f(K)=-AD/t^2\le0\). The initial lower bisection endpoint \(A\) lies below or at the smaller root: since \(t\le1\), \(AC\le t^2/4\) places \(A\) left of the polynomial minimum, and \(f(A)=CA^2+2PA\ge0\). Retaining the feasible upper endpoint throughout bisection is therefore sound. The zero-response and zero-radius branches are correct.

A non-strict box inequality is enough. Positive monotone Volterra iteration starts from zero (or the positive linear response), stays below \(a+bK\), and converges by the bounded local Lipschitz/factorial estimate on the finite interval. A strict first-exit bootstrap is unnecessary. In particular the independent tangent example \(a=1/4,b=1,K=1/4,D=0\) must be accepted and is accepted.

The integer square-root algorithm encloses \(\ell^{-5/2}\) in the required direction. Dyadic upward rounding occurs before the quadratic decision, then again for the inherited endpoint. The endpoint radius after the second rounding need not equal or lie below `closure_K`. The implementation and note correctly keep those distinct: `closure_K` controls the unrounded box, while the stored rounded component bounds dominate that box and are safe future initial data.

## Inheritance, jumps and partial failure

Recursive bisection mutates the inherited endpoint and elapsed time only after an accepted child. Python's short-circuit `and` prevents the second child from running if the first failed. If the first child succeeds but the second fails, the accepted first child's endpoint and time remain the certified prefix. The independent synthetic Riccati example forces exactly this situation: a requested duration \(3/2\), at most one bisection, and initial \(q/1000\) returns precisely a validated duration \(3/4\). Its endpoint encloses the exact \(q\alpha(t)\), \(\alpha'=284\alpha^2\), solution.

Representation jumps are applied once before an original slab and are inherited through every internal split. Replacing the preceding example's initial vector with an equal time-zero jump produces exactly the same accepted endpoint and only one boundary-update record. A failure immediately after a jump leaves an endpoint bound at that same time in the **new** representation; consumers should retain the boundary-update metadata when resuming. This is valid zero-time inheritance, not progress through the failed slab.

The fixed verification length and `reweight_endpoint` scaling are correct. Neither subdivision nor reweighting erases prior error. The failure labels distinguish resource exhaustion and conservative box failure; neither is a claim about singularity of the comparison system or NS.

## Concrete bugs found and fixed

Before the fix, `Bounds.read` validated JSON-derived data but the public dataclass constructor bypassed those checks. Two exact reproductions were:

1. `Bounds(Fraction(1,100),0,(0,)*4,(-1,0,0,0,0))` (with `Fraction` imported), passed directly to `certify_slab` with initial `(1,0,0,0,0)`, returned a purported order-zero endpoint below one. The supplied negative norm bound violates the documented premise, but the public API should reject it rather than emit a certificate.
2. A direct `Bounds(-1,0,(0,)*4,(0,)*5)` with zero initial error reached the zero-error shortcut and accepted a negative duration.

The owning agent added validation/canonicalization to `Bounds.__post_init__`, so both cases and a binary-float duration now fail before propagation. `Options.__post_init__` also validates direct construction. `close_radius` now rejects unequal-length or negative vectors instead of relying on `zip` or its internal caller. The convolution helper rejects negative orders/tails. I read these changes and the independent regression tests pass. Valid previously supplied JSON calculations retain the same mathematical result.

The low-level kernel still necessarily relies on its documented premise that its supplied nonnegative `exp_tail` truly bounds the relevant exponential tail. Production `linear_flow_enclosure` supplies that quantity through the proved exact routine. Passing an invented small tail directly to a low-level helper is not an independently certified input.

No remaining blocking correctness issue was found in the supported path. The next missing work is supplying certified whole-space, whole-slab bounds and an approximation reconstruction, not further synthetic propagation examples.

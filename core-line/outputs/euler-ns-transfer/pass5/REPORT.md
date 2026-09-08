# Fifth pass: a certified initial feedback inequality

8 September 2026. **The full initial pressure inequality is now certified for an explicit compact smooth, unforced ordinary Navier–Stokes datum. Same-solution return and blowup remain open.**

Replacing the thin cylindrical pump with a radial annular strain makes the entire initial pressure-feedback test reduce exactly to angular degrees zero, two and four and bounded radial integrals. This reduction includes the nonlocal pressure response and all viscosity terms. It applies to this initial scalar functional; it is not a finite-dimensional evolution model.

For the exact datum with central strain and rotation both one, annular strength 7/5, and viscosity 1/1000, outward-rounded interval integration and exact analytical bounds give

\[
p_{zz}'(0)+32\le-\frac{290649}{8750}<-33,
\qquad (b/\Omega)''(0)\ge\frac{290649}{17500}>16.
\]

The [main proof](CERTIFIED_INITIAL_FEEDBACK.md) specifies the datum, exact amplitude tuning and complete component identity. The [independent reconstruction audit](full-gate-certificate-audit.md) checks that no interaction or viscous contribution was omitted. Multiplying velocity by 1000 and accelerating time by 1000 converts the same construction to viscosity one.

The [finite-time corollary](finite-feedback-corollary.md) proves that one actual local solution has a positive interval with increasing strain/rotation ratio, strain growth faster than 2b², increasing normalized exterior swirl stress, and actual receiving-velocity/Reynolds gain. It gives a sufficient interval in terms of the new datum's H10 norm. That norm has not been numerically enclosed, so this package claims no decimal duration or useful-size stage gain.

A [robustness corollary](robust-initial-data.md) also gives strict favorable initial signs with the direct rational swirl amplitude A=1100, without exact neutral tuning. This is a separate local datum; it does not provide an invariant class.

This supersedes the immediate initial-sign task for a **new radial datum**. It does not certify the older cylindrical candidate's approximate value of −241.35. The old data and reports remain unchanged.

The next mathematical obligation is to carry the actual terminal velocity and inherited exterior geometry into a reusable class with retained gain. The [return calculation](actual-return-next-lemma.md) proves an unavoidable first-order cubic core deformation: p_zzzz(0) >= 3693184/105861 > 34. Receiver rescaling cancels the central matrix's first-order change but cannot erase this higher derivative. The [inherited-geometry calculation](inherited-geometry-first-variation.md) supplies exact packet-moment rates and a necessary aspect-ratio return budget, and retains the packet-induced meridional deformation. The [independent finite-time audit](finite-feedback-independent-audit.md) reviews the local corollary and higher-jet calculation.

Favorable initial derivatives do not establish persistence to a singular time, and choosing smaller receivers does not supply an infinite sequence of compatible stages. A proposed larger return class is stated explicitly; its invariance is unproved.

To reproduce, run `python3 certificate/run_checks.py`, then `python3 certificate/certify_full_initial_gate.py --recompute` from this directory. The latter recomputes all four interval integrals. The saved [certificate](certificate/FULL_INITIAL_GATE_CERTIFICATE.json) contains the exact dyadic endpoints. [Verification and limitations](VERIFICATION.md) describe the trust base: analytical arguments and mpmath interval arithmetic, without a formal proof-kernel replay or an NS blowup proof.

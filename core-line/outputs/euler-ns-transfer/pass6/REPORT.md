# Sixth pass: reject the axisymmetric return and quantify three-dimensional transfer

8 September 2026. **The proposed indefinitely repeating axisymmetric rotating-core class is ruled out. A concrete three-dimensional initial-data family retains favorable local feedback and supplies positive mean torque, but a usable repeated stage remains unproved.**

## The obstruction changes the next task

The actual axisymmetric circulation maximum principle gives

\[
 G(V_{n+1})\le q_n^\alpha G(V_n),\qquad
 G(V)=\|rV_\theta\|_\infty .
\]

Under indefinitely shrinking scales and fixed \(\alpha>0\), this tends to zero. The pass5 return class's positive normalized rotation and uniform local profile bounds instead force \(G\ge g_*>0\). The [circulation proof](circulation-return-obstruction.md) derives the contradiction, a finite return-count bound, and a global ceiling for the rotational receiving Reynolds observable. It uses [Lei–Zhang's established scalar maximum principle](https://msp.org/pjm/2017/289-1/pjm-v289-n1-p06-s.pdf). Finite smooth axisymmetric forcing with a finite weighted torque budget cannot remove this obstruction.

Adding a bounded cubic core coordinate to the old class is therefore insufficient. This is a stronger restriction than the first-order shape defects found in pass5. The earlier finite local certificate and receiving gain remain valid.

The [escape audit](escape-route-audit.md) shows that allowing rotation to vanish does not rescue a class with uniform global C² bounds over complete normalized stages: its meridional profiles converge locally to zero. This conclusion needs those stronger global assumptions. At the critical exponent \(\alpha=0\), a uniform full-stage velocity bound implies the axisymmetric Type I bound excluded by [Chen–Strain–Tsai–Yau, Theorem 1.1](https://arxiv.org/pdf/0709.4230). Neither statement proves regularity for arbitrary axisymmetric flows. A genuine second scale or loss of those uniform bounds remains a separate research possibility.

## A concrete three-dimensional initial transfer

The [full angular-momentum equation](nonaxisymmetric-torque-budget.md) identifies the terms absent under axisymmetry. After azimuthal averaging, direct pressure torque cancels; transfer into mean rotation comes from actual quadratic velocity correlations.

An explicit compact solenoidal fourfold perturbation supplies strictly positive initial torque into a fixed receiver. The [receiver calculation](rotational-receiver-torque.md) proves the exact sign and frequency-independent bound

\[
 |Q_L(v)|\le C_\chi\|v\|_2^2.
\]

Thus a required receiver increment \(D\) from this stress term must pay
\(\int_0^\tau\|v\|_2^2dt\ge2I_LD/C_\chi\).
Increasing carrier frequency cannot provide unlimited mean transfer at fixed fluctuation energy.

The [three-dimensional embedding](three-dimensional-initial-gate.md) retunes the remote swirl exactly and gives an explicit norm-form interval of nonzero perturbation amplitudes for which

\[
 p_{zz}'(0)+32\le-\frac{290649}{17500},\qquad
 (b/\Omega)''(0)\ge\frac{290649}{35000}>0.
\]

These estimates use the complete three-dimensional pressure. The [functional audit](three-dimensional-gate-audit.md) verifies the Sobolev constants and retuning requirement, and the [full-family audit](nonaxisymmetric-family-audit.md) checks the actual seed, central symmetry, pressure tuning and receiver transfer together. No numerical useful-size perturbation amplitude or duration has been certified.

## The remaining decisive test

Initial compatibility and positive local torque do not prove enough sustained torque for a return. The compact receiver seed can initially draw on the existing outer circulation budget; it does not immediately increase that global mean maximum.

A separate [outer-supported family](outer-mean-maximum-target.md) addresses that distinction. At seed amplitudes \(3/4\le\lambda\le1\), exact central-pressure retuning remains possible, and a cutoff interval certificate proves an actual initial mean-maximum growth rate exceeding \(626/35\), including viscosity. **The full central pressure-derivative gate has not been checked at those finite amplitudes.** This is distinct from the small-amplitude family that retains the gate.

The [remaining pressure test](outer-pressure-compatibility-test.md) reduces the full calculation to two new interaction coefficients. At the already torque-admissible amplitude \(\lambda=3/4\), a sufficient condition is

\[
 D_\nu+A_{3/4}E < \frac{5093339}{210000}\approx24.254 .
\]

Neither interaction coefficient has been evaluated or enclosed. The [independent compatibility audit](outer-compatibility-independent-audit.md) checks the expansion and inequality directions. This is the immediate finite-amplitude compatibility test, not an asserted inequality.

The next target is one actual three-dimensional stage with a common amplitude/time window that simultaneously preserves the full feedback margin, supplies the required evolved stress, pays viscosity, and enters a controlled inherited successor state. The endpoint includes the old field, new angular modes, pressure tails, core deformation and outer geometry. No fresh packet or affine-core reset is allowed.

The previous pass's twelve checks and interval certificate remain unchanged. Run python3 checks/run_checks.py for this pass's finite identities. [Verification](VERIFICATION.md) and [sources](SOURCES.json) state the scope. No same-solution return, infinite compatible sequence, singular solution, or formal proof-kernel verification has been obtained.

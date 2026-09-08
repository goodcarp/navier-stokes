# Seventh pass: the same three-dimensional family passes both initial tests

8 September 2026. **The finite-amplitude compatibility gap from pass6 is closed.**
The same compact, smooth three-dimensional Navier–Stokes datum now has certified
favorable central pressure feedback and positive initial growth of the global
azimuthal-mean angular-momentum maximum, including viscosity.
**A useful quantitative stage, inherited return and blowup remain unproved.**

The [common-family theorem](common-family-feedback.md) applies to every seed
amplitude 3/4 <= lambda <= 1, with the earlier exact pressure retuning:

| Quantity for the actual datum, at viscosity 1/1000 | Certified result |
|---|---|
| Central pressure | p_zz(0) = -8 exactly |
| Complete pressure-feedback test | p_zz'(0)+32 < -435297/70000 (below -6.2185) |
| Initial ratio acceleration | (b/Omega)''(0) > 435297/140000 (above 3.109) |
| Mean angular momentum at an initial maximizing circle | time derivative >626/35 (above 17.885) |

These requirements now hold on **one and the same family**. The earlier tiny
inner perturbation is not being combined with a different outer datum.
The fixed-point mean estimate and smooth local existence give simultaneous
strict gains for a common short interval, but no useful numerical interval
length or prescribed finite gain has been enclosed.

## What resolved the full pressure calculation

The [D kernel](D-pressure-kernel.md) includes all core and meridional pressure
responses and proves D_nu < 4-(96/5)Cw. Its viscosity estimate requires one
new one-dimensional outward range integral.

The [E reduction](E-interaction-reduction.md) retains the nonlocal cross
pressure of the strong outer swirl and the seed. Its complete bound is
E<1/60. The proof uses an exact radial Green trial and an energy estimate
for its whole-space residual, with both radial tails included.
The 2048-panel range calculation bounds E above by less than 0.016616141.
That is an upper bound, not an approximation to the coefficient's value.

Keeping the favorable Cw term while substituting these two estimates closes
the inequality for the full amplitude interval. Independent audits check
[D](D-pressure-independent-audit.md), [E and the common interval](E-certificate-independent-audit.md),
and [the actual interval implementation](outer-E-interval-independent-audit.md).

## What remains

The axisymmetric return obstruction from pass6 remains valid.
The new candidate permits initial three-dimensional transfer, but its evolving
correlation, viscosity, pressure, core deformation and outer geometry must
still be controlled together. Positive initial derivatives do not justify
integrating a favorable Riccati inequality to blowup.

The next task is a useful quantitative interval that retains gain and leaves
a controlled inherited endpoint, followed by a genuine return or another
compatible continuation mechanism. There is no packet reset, artificial
prescribed stress, singular solution, or formal proof-kernel replay here.

Run python3 checks/run_checks.py to recompute both decisive interval
calculations and the seven finite checks. [Verification](VERIFICATION.md)
distinguishes the analytic proof, interval evidence and open obligations.

The [next-stage shear diagnostic](next-stage-shear-diagnostic.md) finds favorable initial winding but fast pressure and polarization changes. Its roughly 10^-3 target interval is a model diagnostic, not a certified duration for the actual solution.

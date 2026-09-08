# Current research state

Updated 8 September 2026, after the ninth transfer pass.
**The full Navier–Stokes goal remains open. No singular solution, complete
blowup proof or formal proof-kernel verification has been obtained.**

The preceding eighth pass made concrete progress by constructing and auditing
three compatible initial gains. This ninth pass also makes concrete progress:
it encloses the new seed's actual initial receiver deceleration, evaluates
a resolved nonlinear planar gain pulse, proves a finite gain in an exact
affine linear envelope model, and derives a full-field validation route.
An explicit rational amplitude now removes implicit pressure retuning from
the selected next initial datum.

## Selected next full-field test

Use `u0=M+1020v+(1/4)w_L`, with the unchanged compact leading-envelope
profiles and viscosity `1/1000`. This is an explicit initial datum, not a
regenerated endpoint. The [fixed-amplitude calculation](pass9/fixed-amplitude-initialization.md)
and its independent audit give actual initial bounds:

- Core ratio derivative `beta'(0)>0.02484` and acceleration `beta''(0)>6.09`.
- Receiving mean derivative `G_t(0,r_*,4)>102.49`.
- Fractional nonaxisymmetric energy derivative `K'(0)/K(0)>5.829`.

The amplitude intentionally gives a strict initial core cone. Its correct
off-neutral identity is used; the old neutral shortcut is not propagated.
The pressure-neutral pass8 family remains a valid historical result.

## What the new evolution evidence establishes

The [full-pressure new-envelope jet](pass9/general-envelope-jet-evaluation.md)
retains the steep envelope, viscosity and axial covariance flux. Its fresh
whole-space pressure-gradient error is below 1.3. Uniformly over the stated
independent amplitude ranges, the initial fixed-circle acceleration is
below -13000; tracking its radial maximum at fixed packet-center height
still gives acceleration below -12000. These are initial statements, not
turning-time bounds or acceleration bounds for the global spatial maximum.

In the nonlinear planar model, the selected amplitude shows about 0.175%
gain in the tracked mean branch and 14.6% fluctuation-energy gain at time
0.001. Radial/time refinement supports that observation. A longer A=950
comparison run shows a later torque reversal and loss of fixed-receiver
gain. These are finite-radius slice calculations; they omit the compact
axial pressure, evolving exterior strain, central core and older field.
They are not time evolution of the full compact 3D datum.

The separate affine linear model has a certified finite-time energy gain
above 15% for the entire envelope and a finite positive shear-work budget.
Its energy normalization, local torque and missing nonlinear terms are
kept separate from both the slice and the actual NS solution.

## Current next action

Build and validate one complete compact three-dimensional approximation
from the explicit datum, carrying all mean, axial, core and older-field
effects. The [H4 endpoint bridge](pass9/evolution-validation-bridge.md)
can certify retained mean, energy, core-cone and full-velocity gains without
requiring every intermediate derivative to stay positive. It must also
enclose the actual inherited endpoint in a successor class stable in its
named norm. H7 curvature validation remains an optional stronger route.

No full-field residual, useful actual duration or terminal-class margin
has yet been enclosed. Numerical reconstruction, exterior tails, domain
errors and every time-slab error still need to be carried. A fresh packet
or a reduced-model pulse cannot replace this step.
[NEXT_TARGET.md](NEXT_TARGET.md) retains the full-stage and infinite-assembly
obligations.

Eight packaged finite programs passed, including two recomputed interval
certificates and the fixed-amplitude exact algebra. Pass2–pass8 manifests
remain unchanged. See the [ninth-pass report](pass9/REPORT.md) and
[verification](pass9/VERIFICATION.md) for scope.

The pass6 bounded axisymmetric rotating-core return obstruction remains
valid under its stated assumptions; no general axisymmetric regularity
theorem is asserted. There is no missing authorization or external blocker.
The remaining gap is mathematical, and the full goal remains active.

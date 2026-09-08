# Tenth pass: full-field calibration and usable error estimates

**The full Navier–Stokes goal remains open.** This pass builds the complete
initial 3D field and a pressure/coefficient prototype, proves an obstruction
to the linear Taylor validator, and separates several numerical artifacts
from the actual initial dynamics. It has not produced a validated full 3D
endpoint, inherited return, singular solution, or formal proof-kernel replay.

## What is now concrete

The [full initial evaluator](initial-full-field.md) represents the selected
`u0=M+1020v+(1/4)w_L` without changing its profiles or retuning its amplitude.
It includes the whole radial pump, both axial endcaps and all initial
angular interactions. The [cylindrical representation](cylindrical-NS-representation.md)
derives the complete velocity, pressure, helical diffusion, axis regularity
and projection identities. Independent algebra and Cartesian-difference
checks support the implementation.

The [whole-space integrals](initial-whole-space-integrals.md) provide
independent initial pressure and energy references. They give approximately
`beta'(0)=0.35452`, `K=6.44377`, and `K'/K=69.43544`. These decimal values are
quadrature references, distinct from the already certified pass9 initial
inequalities. The full pump carries most of the initial energy; a planar
slice cannot represent this pressure or its evolution.

## The numerical failures are identified, not counted as gains

The [coefficient findings](full-coefficient-findings.md) document the
production grids and their limits. Independent manufactured pressure tests
show convergence, but a small scalar solve residual does not imply a
commuting projection or correct axis orders. The refined production field
still has a substantial divergence error. A coarse apparent failure of the
core condition is traced mainly to axial differentiation of large outer
convection terms, not to a reliable NS trajectory.

The [image calculation](initial-periodic-images.md) independently accounts
for the axial-domain pressure error: the period-16 correction is about
`+0.48013`; period 32 gives `+0.008383`. After those corrections, the refined
solver differs from the reference by about `5e-5`. An exact kernel/energy
argument bounds the **omitted image tail**; the retained image integrals
and the solver discretization are still numerical, not interval-enclosed.

The large quadratic-polynomial energy increase at time 0.001 is not accepted
as growth. The [co-rotating construction](corotating-time-polynomial.md)
resums rigid angular phase while retaining the full pressure, generated
modes and harmonic tails. It reduces that artifact, but the remaining
spread is large. No global frame removes all differential rotation or
the evolving meridional field.

## What the validation argument must change

The [time-polynomial obstruction](time-polynomial-residual.md) proves that
the linear Taylor field cannot pass the specified H4 core-cone validator
at times `T>=1/16250`. At 0.001 its residual requires a radius above 0.0065,
while its core margin permits less than 0.0004. This is an obstruction to
that approximation and certificate, not a lower bound on actual NS error
or an exclusion of a useful physical stage.

The [sharper error hierarchy](sharp-error-hierarchy.md) addresses a second
problem. An exact lower bound `||u0||H5>12000` makes the original uniform
growth exponent exceed 6720 at time 0.001. A derived five-component error
system instead puts `(j+1)||sym grad v||infinity` on its diagonal and keeps
derivatives two through five in explicit lower-triangular couplings.
Rigid rotation cancels. This is a conditional estimate with checked
constants, not a numerical enclosure of its coefficients or a proof of
physical stability. Smaller verification lengths have explicit endpoint
and nonlinear costs; no favorable choice has yet been certified.

The next approximation must resolve the full viscous evolution and use
an error estimate with usable growth constants. Rotation, mean feedback,
axial motion, pressure tails and inherited errors must all remain in the
same field. A positive initial jet or a reduced-model pulse cannot replace
the required endpoint and successor-class membership.

See [VERIFICATION.md](VERIFICATION.md) for the exact scope of the checkers
and [the current target](../NEXT_TARGET.md) for the remaining finite-stage
and infinite-assembly obligations.

# Eleventh pass: evolving the full field and preserving the evidence

**The Navier–Stokes blowup goal remains open.** This pass replaces the
initial time-polynomial diagnostic with an uninterrupted, fully coupled
three-dimensional numerical evolution on a finite cylinder. It also
implements an exact-rational propagator for the conditional error hierarchy.
Neither component supplies a validated whole-space stage or a usable return.

The selected datum remains `u0=M+1020v+(1/4)w_L`, with viscosity `.001`.
There is no amplitude retuning, replacement packet, phase reset, or affine
core reset. The mean flow, radial and axial velocities, both endcaps, viscous
diffusion and generated angular modes all evolve together.

![Recorded finite-cylinder pulse](trajectory-summary.png)

## The numerical experiment now has a coherent discrete evolution

The [operator audit](full-mac-independent-audit.md) establishes the intended
weighted divergence/gradient adjoint, orthogonal pressure projection,
positive viscous gradient form, constrained Crank–Nicolson step and
energy-neutral Lamb nonlinearity. These are properties of the finite
discretization. The radial boundary is a no-slip wall, and the axial
direction is periodic.

A [separate audit of the faster implementation](fast-full-mac-independent-audit.md)
compares fused Fourier products, cached sparse matrices and signed-frequency
paired solves with the reference. The largest complete-step difference in
its six small-grid comparisons is below `6e-16`. This speedup changes neither
the model nor its boundary conditions. Its raw nonlinear output is already
restricted to retained axial modes; its norm must not be used as a bound on
the complete unresolved nonlinear product.

Conservation alone is insufficient. At radial resolution 192, odd axial
grids 1025 and 2049 give virtually the same initial fluctuation energy but
initial core rotation derivatives approximately `-4.7542` and `1.99684`.
The exact initial identity is `Omega'(0)=2`. Earlier coarse/filter histories
are retained as failed calibration diagnostics. The completed time run
uses 2049 axial points and modes `0,4,8,12,16`.

## What the first completed run actually shows

The baseline uses 192 radial cells, radius 8, axial period 16 and time step
`1e-4`, ending at `T=.001`. These are **finite-cylinder numerical values**:

| Observable | Represented initial field | At T=.001 |
|---|---:|---:|
| Nonaxisymmetric energy K | 6.456398 | 6.887351 |
| Mean circulation at the original receiver | 41.870181 | 41.715110 |
| Mean circulation at the selected shifted receiver | 41.870182 | 41.879439 |
| Core rotation Omega | 1.000000 | 1.001995 |
| Core strain minus rotation, b−Omega | −4.55e-6 | 6.60e-5 |

The energy increase is about 6.67% relative to the same represented initial
field. This is much smaller than the old quadratic-polynomial spike. The
fixed receiver loses mean circulation; the selected receiver moves from
radius `.215471` to `.219724` at axial height 4 and retains a small increase.
The selection is a radial stationary-point search in `.18<r<.3`, not a
proved global maximum or a guaranteed continuously tracked branch.

At the endpoint, energies in modes 8,12,16 are approximately `.036612`,
`.00078567`, and `.00002427`. Generated modes are retained, but the small
last retained mode does not bound omitted modes. Discrete divergence stays
near `4.4e-11`; the integrated energy-balance defect is about `-.001473`.
That defect is small compared with total energy but still has to be compared
with the much smaller gain margins.

The run with half the time step also completes. It gives `K(T)=6.888747`,
a 6.70% increase over its same represented initial value; the selected
mean value is `41.879450`, and `b−Omega=6.5191e-5`. The original receiver
still loses mean circulation. Its energy-balance defect falls to about
`-4.29e-5`.

The [phase-preserving endpoint comparison](experiments/resolved-time-comparison.json)
finds an absolute full-field L2 difference `.0026355` and nonaxisymmetric
L2 difference `.0026214`, or `.0706%` of the finer run's nonaxisymmetric
amplitude. Endpoint K differs by `.0013963`, the selected mean by
`1.0752e-5`, and the core ratio by `8.49e-7`. This supports time-step
consistency of these finite-cylinder diagnostics; it is not a certified
time-error remainder or a spatial/domain convergence result. The
[endpoint criteria](trajectory-diagnostic-criteria.md) retain those separate
requirements and the eventual certified margins.

## The principal unresolved error is visible already at initialization

The discrete projection changes the sampled initial field by about `1.80`
in its L2 norm on the 192-cell grid. A separate 384-by-4097 initial probe
reduces this to `.4263`, gives `K0≈6.444111`, and recovers
`Omega'(0)≈1.99918`. Its initial receiving derivative is about `100.469`
and energy derivative `K'(0)≈441.195`, compared with the independent
whole-space references from pass10. These are improvements, not enclosures.
In particular, the receiving derivative is still below the pass9 exact
initial lower bound `102.495`, so the finer probe has demonstrable remaining
representation/discretization error despite its better core calibration.

The endpoint core margin is only of order `1e-5` to `1e-4`. Spatial errors,
the initial representation change, periodic images, the radial wall, and
smooth-axis reconstruction have not been bounded within it. The earlier
initial image correction at period 16 is substantial; it cannot simply be
multiplied by T to obtain a later domain-error estimate.

## The error propagation is implemented; its fluid inputs are still missing

The [rational slab propagator](validated-error-propagation.md) advances the
five ordered-derivative error components from pass10, inheriting every
previous endpoint and any representation jump. Its linear flow is enclosed
by finite increasing-index paths and positive rational series. Large
off-diagonal derivative coefficients do not enter its exponential argument.
An exact quadratic feasibility test closes each nonlinear slab, with
adaptive subdivision and explicit partial-failure reporting.

Its own checker passes 44 tests; an [independent implementation audit](validated-error-independent-audit.md)
passes 161 assertions. All example coefficients are manufactured. No
actual fluid trajectory has certified whole-slab strain, higher derivative,
residual and initial-error inputs for this propagator. A tested comparison
algorithm is not a certified NS error radius.

## A route to the next approximation is specified

The [poloidal–toroidal design](poloidal-toroidal-evolution-design.md) gives
compact smooth potentials for the unchanged datum and a possible
divergence-free Galerkin representation. A [shared-Hankel comparison](solver-strategy-review.md)
identifies why mixing order-dependent radial frequency grids breaks the
simple helical projector. Finite sums of point-frequency Bessel modes also
fail the required whole-space finite-energy condition.

The [exterior-energy calculation](exterior-harmonic-mass-audit.md) and
[implementable boundary blocks](exterior-mac-block-specification.md) replace
the wall pressure condition by a compatible harmonic-exterior metric and
trace-matched positive viscosity form. Those blocks are independently
checked but are not used in the reported production run. A harmonic
exterior approximation still needs a viscous-tail residual and sufficiently
smooth joins; the true viscous exterior cannot simply be declared irrotational.

The next decisive result remains a smooth whole-space approximation with
enclosed residual and inherited error, retained endpoint gains, and
membership in a named stable successor class. No numerical gain here
establishes that class or an infinite singular chain.

See [VERIFICATION.md](VERIFICATION.md) for reproducibility and exact scope.
The [bounded public-release check](new-release-status.md) found no new
ordinary-NS manuscript in the checked primary channels; it makes no claim
about unpublished work or who will finish first.

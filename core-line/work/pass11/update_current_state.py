#!/usr/bin/env python3
"""Update only current top-level navigation; preserve pass2–pass10 archives."""
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parents[2]/'outputs'/'euler-ns-transfer'
status='''# Current research state

Updated 8 September 2026, after the eleventh transfer pass.
**The full Navier–Stokes goal remains open. No singular solution, complete
blowup proof, validated useful compact 3D stage, or inherited return has been
obtained.**

The [eleventh pass](pass11/REPORT.md) supplies an uninterrupted numerical
evolution of the complete selected field on a finite cylinder, a compatible
discrete pressure/viscosity formulation, independent operator audits, and
an exact-rational comparison-system propagator. These close implementation
gaps, not the whole-space endpoint or infinite-assembly obligations.

## Selected datum and prior established results

Keep `u0=M+1020v+(1/4)w_L`, viscosity `1/1000`, and the unchanged compact
profiles. No new packet, amplitude retuning, phase reset or affine-core reset
is used. The [pass9 initial-sign certificates](pass9/fixed-amplitude-initialization.md)
and correct off-neutral pressure identity remain in force. Whole-space
initial quadratures from pass10 give numerical references
`beta'(0)≈.35452`, `K≈6.44377`, and `K'/K≈69.43544`; these decimals are
distinct from the rational initial-sign certificates.

Earlier actual local finite estimates, reduced models and obstructions
retain their own scopes. The bounded axisymmetric return obstruction does
not assert general axisymmetric regularity. No initial jet, changing-data
family, isolated affine pulse or reset quotient supplies a reusable stage.

## Full numerical evolution and its limits

The new MAC representation has a weighted-compatible divergence/gradient
pair, orthogonal pressure projection, positive viscous quadratic form and
dealiased energy-neutral nonlinear generator. Independent checks compare
its faster implementation with the reference to floating-point roundoff.
The reported runs retain modes `0,4,8,12,16`, mean feedback, meridional and
axial motion and both packet endcaps.

The completed baseline reaches `T=.001` on a 192-by-2049 grid, radial radius 8
and axial period 16. Its fluctuation energy rises from about 6.4564 to 6.8874.
The original receiving circle loses mean circulation; a shifted receiver
retains a small gain. The core cone margin is only about 6.6e-5. These are
finite no-slip/periodic cylinder diagnostics, not actual whole-space gains.
The half-time-step run also completes, with K about 6.88875, core cone
margin 6.519e-5, and the same receiver behavior. The phase-preserving
full-field L2 difference is .0026355; core-ratio values differ by 8.49e-7.
This is a consistency comparison, not an enclosed time-error remainder.

Coarse axial grids can conserve energy while producing the wrong core
derivative. The initial discrete projection itself has L2 change about 1.80
at 192 radial cells, falling to .426 on the 384-cell probe. Smooth Cartesian
axis orders, full spatial and angular truncation, radial-wall and periodic
image errors, and whole-space tails remain unbounded. Positive endpoint
decimals cannot bypass those errors.

The [pass10 image bounds](pass10/initial-periodic-images.md) explain a large
initial periodic-pressure correction and certify only the omitted image
tail. They do not convert a later periodic solution into a whole-space one.
The linear Taylor validator obstruction and quadratic phase artifacts
remain valid under their stated assumptions.

## Validation and the next approximation

The [sharp ordered-derivative hierarchy](pass10/sharp-error-hierarchy.md)
now has a [rational slab propagator](pass11/validated-error-propagation.md).
It encloses the positive triangular linear flow, closes an explicit
quadratic nonlinear inequality, and inherits prior errors across every
subdivision and representation jump. The accompanying examples are
manufactured; no actual whole-slab fluid coefficients, residuals, initial
errors or useful comparison radius have been certified.

[Potential-based initial data](pass11/poloidal-toroidal-evolution-design.md)
and [compatible exterior blocks](pass11/exterior-mac-block-specification.md)
specify improvements to the approximation. They are independently checked
designs, not hidden features of the no-slip production solver. Exterior
vorticity/heat tails and sufficiently smooth joins still need bounds.

The required next result remains a controlled smooth divergence-free,
axis-regular whole-space evolution with retained gains and membership in
a named stable successor class. [NEXT_TARGET.md](NEXT_TARGET.md) retains
the finite-stage and infinite-assembly obligations. A useful finite
endpoint would not itself prove a singular chain.

See [verification scope](pass11/VERIFICATION.md). No new proof-kernel
replay is claimed by this Codex pass. Historical pass2–pass10 archives remain
unchanged. There is no missing authorization or external blocker; the
remaining gaps are mathematical and computational, and the full goal stays
active.
'''
(root/'STATUS.md').write_text(status)
path=root/'REPORT.md';text=path.read_text();start=text.index('**Latest:**');end=text.index('\n\n',start)
text=text[:start]+'''**Latest:** [The eleventh-pass report](pass11/REPORT.md) gives a complete finite-cylinder time evolution, compatible pressure and viscosity operators, independently checked acceleration of the implementation, and an exact-rational propagator for the conditional error hierarchy. The endpoint diagnostics still lack whole-space reconstruction, residual and inherited-return bounds. No validated compact 3D stage or blowup has been obtained; [STATUS.md](STATUS.md) records the current state.'''+text[end:]
path.write_text(text)
path=root/'NEXT_TARGET.md';text=path.read_text();start=text.index('## Immediate target after');end=text.index('## Progress and remaining gap after the third pass')
text=text[:start]+'''## Immediate target after the eleventh pass

Keep the fixed datum `u0=M+1020v+(1/4)w_L`, viscosity `1/1000`, and unchanged
profiles. The initial-sign work is closed. Do not retune the pressure
amplitude or replace the evolved field by a fresh packet or affine core.

**The required result remains a validated full 3D endpoint with retained
gain and inherited geometry.** The [finite-cylinder time evolution](pass11/REPORT.md)
is useful evidence for a modest pulse, with a moving receiver and a small
core margin. It is not a whole-space stage. Its compatible MAC identities
remove the earlier discrete pressure/divergence mismatch but do not enforce
every smooth-axis order or enclose continuous errors.

Use the [diagnostic criteria](pass11/trajectory-diagnostic-criteria.md) to
prioritize independent time, radial, axial and angular comparisons.
Preserve complex phase and compare the complete field as well as the seed;
total energy alone conceals core and torque errors. The initial projection
change and failed coarse axial calibration must remain in the error budget.
Keep radial wall and axial image effects separate. A time-zero image
correction times T is not an evolved domain bound.

Construct a smooth divergence-free, axis-regular whole-space approximation
using the [potential representation](pass11/poloidal-toroidal-evolution-design.md)
or another justified scheme. The [exterior block specification](pass11/exterior-mac-block-specification.md)
supplies a compatible harmonic-exterior discrete model, but smooth joins
and omitted viscous vorticity/heat tails must still be charged. Merely
changing a pressure boundary row or filtering axis values after projection
would not preserve the required coupled identities.

Enclose whole-slab symmetric strain, derivative orders two through five,
all projected-residual components, initial representation and subsequent
jump errors. Supply those actual bounds to the [rational error propagator](pass11/validated-error-propagation.md)
or another justified estimate. Its manufactured tests are not fluid inputs.
Choose a verification length with its explicit nonlinear, pointwise and
successor-H4 costs. Every later slab inherits the preceding error vector.

Apply the [endpoint bridge](pass9/evolution-validation-bridge.md): retained
mean/seed/core and full receiving velocity/Reynolds gains must survive the
radius. Prove that the actual normalized inherited endpoint lies in a named
stable successor class, including old fields, deformations, phases and
tails. A positive finite pulse and a numerical receiver shift do not supply
that class. Positive derivatives throughout the stage are optional stronger
tests, not replacements for the endpoint condition.

Uniform repeated compatibility, admissible all-order physical-force
assembly and a genuinely divergent continuation quantity remain separate
obligations. Smooth forcing is allowed. The [bounded axisymmetric return
obstruction](pass6/circulation-return-obstruction.md) and the integrated
fixed-receiver stress budget retain their stated scopes.

'''+text[end:];path.write_text(text)
# Refresh exactly the prior top-level manifest members; nested passes have
# their own manifests and are deliberately not folded into this one.
manifest=root/'SHA256SUMS.txt';members=[line.split('  ',1)[1] for line in manifest.read_text().splitlines() if line]
manifest.write_text(''.join(hashlib.sha256((root/p).read_bytes()).hexdigest()+'  '+p+'\n' for p in members))
print(root/'STATUS.md')

# Next mathematical target: one viscous return with retained gain

**Open.** This is a precisely posed research target, not an existence theorem or a claim that we are close to a Clay solution.

## Objective

Construct one actual smooth three-dimensional Navier–Stokes stage with fixed positive viscosity and admissible smooth forcing, if needed, that both amplifies an active field and leaves a quantitatively usable next-stage configuration. Build viscosity into the leading evolution from the outset.

The source's continuous shooting and causal selection can guide the argument. Its exactly linear phase cell, off-axis ring geometry and super-separated Euler frequencies cannot be assumed to survive.

## Immediate target after the eleventh pass

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

## Progress and remaining gap after the third pass

The [quantitative core theorem](pass3/quantitative-core-stage.md) establishes one actual finite interval with positive central vorticity and radius-uniform rotational receiving-mode gain. The [receiver lemma](pass3/rotational-receiver.md) converts this to actual weighted RMS velocity and a smaller receiver with increasing Reynolds number. Its selected velocity/Re/energy ratios match the proposed powers at one receiving time. These are single-solution finite estimates, with uniform constants over the stated fixed profile/parameter family.

The missing terminal requirement is now explicit. The output has a first-order symmetric strain, whereas its rotational gain and selected radius change start at second order. [Rescaling cannot remove that profile mismatch](pass3/terminal-profile-gap.md). The next proof must evolve an invariant class containing that strain or supply a genuine return phase, preserving gain and restoring the outer environment. No such class or return has been proved.

The [material circulation](pass3/material-core-transfer.md) remains unchanged through order four. The [fixed receiver gains through flux](pass3/receiver-flux-budget.md) from a changing footprint. A restart must track that supplying reservoir; it cannot infer an unlimited circulation supply from the finite velocity-gain theorem.

## Progress and remaining gap after the second pass

The [active-core pressure lemma](pass2/active-core-pressure.md) now gives explicit compact smooth data for which remote vortex packets cause an initially rotating core to gain axial strain and vorticity at the first nonzero time orders. It is an actual local unforced NS statement at fixed positive viscosity. It does not yet provide a quantified finite gain or usable return.

The [affine model audit](pass2/affine-stretch-gate.md) shows why a restart must restore effective circulation/velocity scale and localization margin, rather than merely stretch the same packet. The [Jeong–Yoneda comparison](pass2/jeong-yoneda-stage.md) supplies an actual viscous gradient-amplification benchmark near the candidate scale exponents, with exact decoupling and domain limitations.

The next coupled-stage argument must evolve the intended mean stress in the velocity. The [mean-feedback calculation](pass2/mean-feedback.md) shows that oscillation alone does not make this stress small. The [correction-window lemma](pass2/correction-window.md) leaves a conditional all-order window, but its uniform estimates for the full viscous residual remain unproved.

## Candidate scales

For a nonaxisymmetric core of size \(L=L_0\rho\), use the provisional budget
\[
\mathrm{Re}=\mathrm{Re}_0\rho^{-1/4},\qquad
U=\frac{\nu\mathrm{Re}_0}{L_0}\rho^{-5/4},\qquad
M=M_0\rho^{-1/16},\qquad h=L/M.
\]
These scales permit viscosity to be smaller than core strain in the principal calculation while keeping geometric-stage energy and dissipation sums finite. They supply no sign, stability, profile or regeneration information. They are an input to a feasibility test, not a prescribed solution.

## The finite-stage lemma must supply all of the following

1. **An actual coupled field.** Choose a smooth divergence-free base plus an active perturbation. Include the perturbation's effects on pressure and the older flow. If a linearized-wave approximation is used, justify a small carrier amplitude or a cancellation of its potentially faster self-interaction, and retain the leading slow mean stress in the coupled evolution.

2. **A quantitative frequency window.** Derive the lower threshold needed to make localization and correction errors small, and the upper threshold allowed by diffusion. Show that the interval is nonempty for the required derivative orders. A phrase such as “take the frequency sufficiently large” is not a completed estimate.

3. **Retained full-field gain.** Prove a lower bound for an actual velocity/vorticity/strain quantity, with its role in the following stage specified. Pay diffusion across growth, transition, return and holding. A principal quotient returning to zero is insufficient, as the exact 0.294-amplitude control demonstrates. An upper propagator bound is not a lower growth estimate.

4. **A usable terminal configuration.** State the target geometric class in a named norm and show that the evolved, fully coupled field enters it. Carry inherited centers, deformations and older fields. Quantify how the receiving core obtains the velocity/circulation scale and restored localization margin needed for its Reynolds number to increase. Do not replace the evolved field by a fresh input profile, silently reset old phases, or rely on a scalar pressure closure.

5. **An admissible physical force.** Bound the vector residual after pressure selection, including every older/newer interaction and temporal join. State mixed derivative bounds with their scale dependence. For the proposed absolute-summability assembly, enforce the all-order condition in the report, or prove another all-order convergence estimate. Smooth forcing is allowed; nonsmooth control disguised as a residual is not admissible.

## Decision table

| Observation | Mathematical consequence |
|---|---|
| Correction threshold exceeds diffusion ceiling | Reject this parameter/profile choice before building an infinite hierarchy |
| Return occurs but retained gain is at most one | Reject it as an amplifying stage; do not count the return as progress toward a cascade |
| Large gradient comes only from changing the initial datum | Another finite family, not same-solution regeneration |
| Wave heating destroys a required phase cell without a replacement estimate | Rebuild the center/gradient argument before using that stage |
| Residual is small only at fixed accuracy or only in a weak norm | No all-order physical-force conclusion follows |
| Positive finite gain but no usable terminal geometry | A finite-stage result only; next-stage regeneration remains open |
| All five obligations hold | A candidate finite-stage lemma worth independent proof audit; infinite compatibility is still separate |

## What comes after one successful stage

An infinite construction would additionally need one fixed smooth initial datum, a causal chain of compatible stages, finite total time, a genuine divergent continuation quantity, and a globally admissible force. Stage energies being summable does not construct that chain. Establish uniform quantitative transfer margins before attempting to formalize it.

Our PASS14 and PASS28 controls remain useful negative tests. They should be used to detect another changing-data family or short isolated episode, not presented as already fulfilling this target.

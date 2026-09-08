# Next mathematical target: one viscous return with retained gain

**Open.** This is a precisely posed research target, not an existence theorem or a claim that we are close to a Clay solution.

## Objective

Construct one actual smooth three-dimensional Navier–Stokes stage with fixed positive viscosity and admissible smooth forcing, if needed, that both amplifies an active field and leaves a quantitatively usable next-stage configuration. Build viscosity into the leading evolution from the outset.

The source's continuous shooting and causal selection can guide the argument. Its exactly linear phase cell, off-axis ring geometry and super-separated Euler frequencies cannot be assumed to survive.

## Immediate target after the ninth pass

**The next test has an explicit initial field.** Choose
`u0=M+1020v+(1/4)w_L`, viscosity `1/1000`, with the unchanged compact
leading envelope centered at 7/50, width 1/10, angular mode four and radial
carrier -20. The [fixed-amplitude theorem](pass9/fixed-amplitude-initialization.md)
preserves initial mean and fluctuation-energy gains and gives a strict
initial core cone. No implicitly retuned amplitude has to be represented.
This field is not a regenerated endpoint of an earlier solution.

**The next required result is a validated full 3D endpoint with retained
gain and inherited geometry.** The [H4 bridge](pass9/evolution-validation-bridge.md)
supplies sufficient residual, error-majorant and endpoint tests. Construct
a smooth, divergence-free whole-space approximation including the actual
axial envelopes, generated vertical motion, evolved mean, central core and
older field. Bound its complete vector residual and initial reconstruction
error. Every subsequent time slab inherits the previous error radius.

Validate prescribed endpoint mean, fluctuation-energy, core-cone and
full-field velocity/Reynolds gains as needed by the successor mechanism.
Prove actual endpoint membership in a named stable successor class after
charging rescaling and coordinate errors. H4 is sufficient only if that
successor lemma is stable in H4. Positive derivative or curvature signs at
every intermediate time are optional stronger conditions, not automatic
requirements of retained gain. Use the H7 alternative if those conditions
or a stronger successor class are needed.

The finite planar simulations select an early time around 0.001 for testing,
not a certified compact-3D duration. At the selected A=1020 they show about
0.175% tracked mean-branch gain and 14.6% finite-radius fluctuation-energy
gain at that time, with refinement checks. The longer A=950 comparison
shows a later torque reversal. Its timing cannot be assigned to the actual
selected solution. The exact affine-envelope finite gain is also a separate
linear model, with a finite shear-work budget.

The actual full-pressure initial-jet certificate now shows deceleration
for the leading seed too, including radial recentering at fixed axial height.
It does not supply an actual turning time or rule out a useful endpoint gain.
Do not substitute those initial jets for a time-dependent error enclosure.

Raw slice data are not already a full-field H4 approximation. They omit
three-dimensional effects and exterior energy tails; residual numerical
mean circulation gives a nonintegrable whole-plane mean-energy tail unless
it is removed in a controlled reconstruction. A small finite-grid energy
balance diagnostic does not bound the complete 3D residual.

The pass8 three-sign result, pass9 explicit-amplitude inequalities and
reduced-model finite pulse are closed diagnostic work. Do not repeat them
as the main target. No affine-core reset, fresh outer seed, prescribed mean
stress or pressure closure is a substitute for the evolved endpoint.
Uniform repeated compatibility and admissible all-order force assembly
remain separate after one actual stage.

The [pass6 circulation theorem](pass6/circulation-return-obstruction.md)
continues to rule out the former bounded axisymmetric rotating-core return.
The conditional escape audit also remains valid under its stated stronger
global assumptions. These do not exclude all second-scale mechanisms or
prove general axisymmetric regularity.

The fixed-receiver stress budget remains |Q_L(v)|<=C_chi ||v||2². A required
increment D supplied by this stress costs integrated perturbation energy
at least 2 I_L D/C_chi. Frequency alone cannot evade that cost.

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

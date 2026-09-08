# Euler stage-return mechanism and the local Navier–Stokes gaps

Date: 2026-09-08. Scope: a focused reading of the supplied paper *Blowup for the Euler Equations with Smooth Forcing*, chiefly §§3, 7, 11 and 12, with coordinate definitions from §2. Source copy: `work/sources/euler.pdf`; extracted text: `work/sources/euler.txt`. This note does not audit the whole proof or its formalization. Paper theorems below are reported as the paper's results; proposed transfer statements are distinguished explicitly.

## Finding

The paper supplies the kind of **one-history iteration** absent from PASS14 and PASS28. Its return mechanism is, however, a **controlled forced construction**. A prescribed base rotation, tiny activated seeds, and higher-order residual corrections jointly prepare each new stage. Neither an autonomous return mechanism for our collars nor a positive-viscosity version follows from the Euler theorem.

The most promising immediately reusable pieces are the weighted amplitude lemma, the continuous shooting argument for a return, and the separation of low-order continuation from derivative estimates. The printed scale schedule cannot simply be reused after adding fixed positive viscosity.

## What is amplified and what is returned

**Reported paper identities — §3, (3.1)–(3.4).** The new layer is transported by the entire older corrected meridional velocity. Its material map begins at the current full-flow point: `X_m(a,tau_m)=Z(tau_m)+a`. The phase is `N_m p_m·A_m`; its covector is `zeta_m=(D_a X_m)^(-T)p_m`. The circulation amplitude `T_m` and normalized principal reduced-vorticity amplitude `OmegaHat_m` solve a material two-component system

\[
\partial_t\binom{T_m}{\widehat\Omega_m}
=\begin{pmatrix}0&\widehat a_m\\\widehat b_m&0\end{pmatrix}
\binom{T_m}{\widehat\Omega_m},
\]

where

\[
\widehat a_m=-\frac{(J\zeta_m)\cdot\nabla\Gamma_{<m}}{\sigma_{m-1}\kappa_m},
\qquad
\widehat b_m=\sigma_{m-1}(2W\Gamma_{<m})\zeta_{m,1}.
\]

These retain the spatially varying Euler metric and coupling. Positive entries and a growing eigenline produce amplification. Positivity alone does not supply the required derivative bounds or the weighted propagator estimate.

**Reported paper return — Lemma 7.4, Lemma 7.5 and Proposition 7.6.** Let `alpha` denote the base rotation rate, and `E_m` the complete older perturbation gradient at this layer's older-flow center. Define

\[
F_m=-\frac{(E_m^T\zeta_m)_1}{\zeta_{m,2}},
\qquad \dot\zeta_{m,1}=(F_m-\alpha)\zeta_{m,2}.
\]

Choosing `alpha=F_m-dot(z_mu)/zeta_m,2` makes the first covector component follow a prescribed profile `z_mu`. It first decreases to zero through a nonnegative transition and then receives a short negative pulse with trial amplitude `mu`. In the rescaled time, the quotient `v=c0 OmegaHat/T` obeys

\[
v'=c_r(1-\eta)-kv^2
\quad\text{then}\quad
v'=-\mu c\Lambda_{\rm ret}\eta_2-kv^2.
\]

Under uniform positive coefficient, support, denominator and entry-ratio bounds, all trials exist through the common endpoint. The endpoint quotient is positive at `mu=0` and negative at `mu=c_p`. Continuity gives a nonempty compact zero set; its least element is selected. No differentiability of this selector is needed or claimed for the infinite construction.

On the selected member, both `zeta_m,1` and the **central principal** `OmegaHat_m` vanish. Switching to `alpha=F_m` keeps them zero and keeps central principal `T_m` constant during holding. This does **not** freeze the complete field, its gradient, off-center amplitudes or correction levels.

## Why this is one history, rather than a succession of unrelated examples

**Reported paper structure — Propositions 7.2–7.3 and 11.1.** Every current trial shares the same already selected past. Feedback is constructed by a Volterra-type contraction: the difference between states produced by trial rates is bounded by the time integral of their difference. The instantaneous rigid rotation cancels from the feedback expression. The triangular layer/level construction has a finite derivative allocation for every requested finite output order.

Continuation uses bounded rotation, positive denominator and symbol margins, and support staying within the fixed chart. It does not start by assuming the high-order estimates it aims to prove. Older-field estimates are uniform over all current trials; the youngest layer's completed holding and subsequent full-flow path are estimated on the selected trial.

At insertion, old material labels, phases and displacements are **inherited**. They are not reset to zero. Complete corrected fields determine the next center and coupling. Two separate estimates are essential:

1. A cone lower bound preserves the **complete** circulation gradient at the full-flow point through holding and the following growth.
2. A contracted transverse-gradient estimate, with older-center displacement terms retained, makes the next layer's growing coefficients positive on its entire insertion collar.

The principal clock `sigma_m^2` is defined from the norm of the principal gradient sum at the selected endpoint. It is not silently identified with the complete gradient. §11 proves their required relation separately.

**Reported paper sequence — §12 and Proposition 12.3.** Fix all profiles, geometry and the derivative-order cost function first. A single sufficiently large `y=log N_1` determines

\[
N_i=N_{i-1}^{Q_i},\quad Q_i=q_0+i,\quad
\ell_i=N_{i-1}^{-4}\;(i\ge2),\quad
J_i=2k_i+8,
\]

with admitted orders `k_i` defined by both `k_i<=Q_i/c_Q` and `K_ho(k_i)<=log N_i`. These orders increase to infinity. The seed is `T_i^seed=-N_i^(-k_i-8)` and the target amplitude is `-N_i^(-7/8)`. The desired logarithmic gain is `(k_i+8-7/8)log N_i`. Angles and growth clocks are chosen causally from previously selected clocks; no future realized field is assumed.

The complete stage window has duration at most `C/sigma_(i-2)` for `i>=2`, so the window lengths sum to a finite terminal time. Flat activations and endpoint profiles match time derivatives when selected histories are joined. The first stage has separate fixed constants and is not justified by fictitious earlier indices.

## Why the force is smooth despite infinitely many stages

**Reported paper estimate — (11.24), (12.14)–(12.17).** The force estimate treats three sources separately:

- Activation of the tiny seed, with its actual short-ramp amplitude.
- The terminal remainder of the finite correction hierarchy; correction depth grows with admitted derivative order.
- Both phase means, including circulation-square and variable-symbol terms. Means have slow physical differentiation rates and cannot be dropped just because an oscillatory inverse exists.

All physical components are recovered using one fixed compact inverse-curl operator, including activation plus terminal remainder together and each compatible mean separately. The result is a bound of the form

\[
\|f_i\|_{C^k_{x,t}([0,T_*))}\le C_kN_i^{-1/2}
\qquad (k\le k_i).
\]

Its constants include physical-coordinate and recovery costs. For a fixed order, sufficiently late terms are summable; the finitely many early terms need a separate all-order lifetime estimate. §12 proves those fixed-layer and base bounds using finite clock moments of every fixed order. This distinction is necessary: eventual admission alone does not control early layers near the terminal time. The common compact off-axis support is also retained. Only the force, not the singular state, is extended through the endpoint.

This is qualitatively different from declaring an arbitrary residual to be a force. Our PASS1 example failed because its residual curl could not extend smoothly through concentration.

## Mapping to the local archive

| Local asset/gap | What the paper supplies | Remaining transfer obligation |
|---|---|---|
| PASS14: large finite amplification across different parents, with shrinking physical intervals | A causal selection of infinitely many compatible stages and a single finite terminal time | Construct such stages for positive-viscosity dynamics or prove a new residual correction scheme; Euler compatibility is insufficient |
| PASS28: one short strain-reinforcing collar episode; initial core irrotational | An active coupled circulation/vorticity amplitude and a precise return of one principal component | Derive a comparable active pair for the collar and identify an admissible mechanism that controls its return; no autonomous collar return is demonstrated |
| PASS28: later pressure and geometric feedback unclosed | Complete older fields, displaced centers and nonlinear material maps remain in every stage | Pay those terms with actual inherited bounds; prescribing fresh collars or a scalar pressure law would recreate the old gap |
| PASS1: nonzero, nonsmooth concentration residual | A hierarchy canceling leading errors and controlling every physical force derivative | Prove force smoothness for any adapted construction, including viscous residuals and all means |
| PASS30–31: an unproved bound near the axis in standard NS | A fixed off-axis Euler cascade | No direct payment of the axis concentration budget; a moving-ring version loses the paper's fixed positive radial and metric margins |

## A reusable lemma and a sharp first gate

**Established ODE calculation, independently checked here.** For a prescribed smooth scalar `d(t)`, consider the common-damping pair

\[
T'=a\Omega-dT,\qquad \Omega'=bT-d\Omega.
\]

Where `T` is nonzero, `v=c0 Omega/T` satisfies exactly the same Riccati equation as without `d`:

\[
v'=c_0b-(a/c_0)v^2.
\]

Multiplying both amplitudes by `exp(integral d)` removes the common damping. Consequently the return quotient and its endpoint signs survive for prescribed coefficients, but

\[
(\log|T|)'=(a/c_0)v-d.
\]

During a nominal Euler hold with `b=Omega=0`, one has `T'=-dT`, not `T'=0`. A successful damped return is therefore not a successful amplifying iteration. If `d`, `a`, or `b` depend on the changed full field, this algebra alone does not preserve the paper's coefficient hypotheses.

**Actionable proposed lemma, not proved for a PDE.** Extend the return comparison to a two-component amplitude system with common damping plus explicit error terms. Require: a uniform nonvanishing denominator/amplitude margin, a quantitative perturbation bound for the quotient sufficient to retain the two endpoint signs, positive net gain after integrating damping through the entire growth/return/hold history, and an inherited next-stage transverse margin. First establish or refute these inputs for one actual viscosity-bearing stage. This identifies an exact gate before building any infinite hierarchy.

**Leading-frequency diagnostic, not an NS construction.** In physical coordinates, the classical axisymmetric diffusion operators for circulation and reduced azimuthal vorticity are respectively

\[
\partial_{zz}+\partial_{rr}-r^{-1}\partial_r,
\qquad
\partial_{zz}+\partial_{rr}+3r^{-1}\partial_r.
\]

In `y=(z,r^2/2)`, these are `partial_11+r^2 partial_22` and that same operator plus `4 partial_2`. Thus both amplitudes have the same leading high-frequency damping

\[
d_m=\nu N_m^2(\zeta_{m,1}^2+r^2\zeta_{m,2}^2),
\]

but lower-order diffusion terms differ and require a separate derivation and error bound. On the paper's fixed off-axis region with nondegenerate covectors, this leading rate is comparable to `nu N_m^2`. Its printed growth clock is comparable to `Lambda_m sigma_(m-2)`, much smaller along its sequence. Any fixed positive viscosity therefore defeats this unchanged leading-amplitude schedule eventually. With the same increasingly separated frequencies, even a fixed positive fractional dissipative order gives the analogous leading obstruction. This does not exclude an entirely different hypodissipative schedule or construction.

For a force-based attempt to preserve exactly the Euler velocity, the exact momentum identity instead requires `f_NS=f_E-nu Delta u`. The paper's smooth bound for `f_E` does not control that added term. Its curl is pressure-independent. The next calculation must address its smoothness through the terminal time, not only smallness on finite regular intervals.

## Transfer verdict

The new paper offers a concrete target for an original next lemma: **a gain-preserving return under unavoidable damping and fully inherited geometry**. The shooting and gluing logic is useful now; the ability to prescribe rotation, activate seeds and correct a residual is where the forced Euler construction has capabilities our unforced NS orbit does not yet possess. Direct transplantation is blocked by that control issue, by leading viscous damping, and by the loss of uniform geometry in any move toward the symmetry axis. None of these observations claims an audit failure in the Euler paper or settles a modified construction.

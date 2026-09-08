# Audit of the auxiliary phase-heat lemma

Date: 2026-09-08. Reviewed `work/phase_heat_lemma.md` and `work/phase_heat_check.py` read-only, and compared the principal variables with equations (2.1)–(2.6) and (3.1)–(3.4) of the supplied Euler paper. Ran `python3 work/phase_heat_check.py` with the system interpreter: exit code 0, all nine recorded checks PASS. No full PDE construction or formalization was audited.

## Verdict

The principal-model algebra, noncommuting-matrix heat factorization, derivative estimate, and reset-with-amplitude-loss example are correct within their stated auxiliary scope. The important limitation is preserved: coefficients and geometry are prescribed in this model, and none of the omitted physical residuals is controlled. A few small formulation changes would make the statement fully precise; no substantive sign or normalization error was found.

## Principal variables and general periodic profiles

At highest phase order, write the perturbation streamfunction as `psi(s,t)`. Then

\[
\xi=N^2\kappa\psi_{ss},\qquad
v_{\rm new}=N J\zeta\,\psi_s
=\frac{J\zeta}{N\kappa}\partial_s^{-1}\xi.
\]

All these equalities refer to the principal phase terms, with slow-label and envelope derivatives omitted. The zero-mean convention for the inverse derivative is essential. Defining `eta=(sigma/N)partial_s^-1 xi` therefore gives the circulation equation

\[
\gamma_t=-\frac{d_{\rm old}}{\sigma\kappa}\eta
\]

at principal order. The linear principal part of the reduced-vorticity source is `N c zeta_1 partial_s gamma`. Applying `(sigma/N)partial_s^-1` gives `eta_t=sigma c zeta_1 gamma` for mean-zero circulation. This verifies both entries of the matrix in the note.

The Euler paper uses a smooth odd periodic `F` that agrees with `s` near zero and a primitive `P` with `P'=F`. Its principal fields are `gamma=T F` and `xi=Omega F'`. Since `F` has zero mean, the chosen primitive of `F'` is exactly `F`, giving `eta=OmegaHat F` with `OmegaHat=sigma Omega/N`. No sinusoidal identity is needed.

The diffusion metric is also correct: the physical principal Laplacian in volume coordinates has symbol

\[
q=\zeta_1^2+r^2\zeta_2^2=r^2\kappa.
\]

For circulation and reduced azimuthal vorticity the two physical diffusion operators differ at lower order but share this symbol. Consequently the common coefficient `d_nu=nu N^2 q` is appropriate for the principal model; using `nu N^2 kappa` without the factor `r^2` would be incorrect before a fixed unit-ring normalization.

An additional useful qualification: heat evolution generally changes the full phase profile. A packet starting with `R_in F` becomes `Phi R_in` times `exp(D partial_s^2)F`. It does not preserve the identity `F(s)=s` near zero. Thus the exact phase-model result cannot retain the paper's linear-phase-core arguments or identify a physical gradient with an unchanged `F'(0)=1`. The current note does not claim that transfer, but this should be explicit in any later stage theorem.

## Factorization and derivative bound

The formula

\[
W(t)=\Phi(t,s)e^{D(t,s)\partial_\theta^2}W(s)
\]

is valid because every component matrix commutes with the scalar phase heat operator. It does not replace `Phi` by an exponential of the integral of `B`, so noncommutativity of `B(t_1)` and `B(t_2)` causes no problem. The script tests an actual noncommuting example with three phase harmonics and verifies the complete PDE residual and initial data symbolically.

The displayed Parseval bound is correct for vector-valued phase-mean-zero data with the Euclidean matrix operator norm. Every nonzero integer frequency obeys `exp(-n^2 D)<=exp(-D)`. Differentiation removes the zero mode when its order is positive. The Duhamel formula with a source is also exact for this phase-independent coefficient problem.

Two precision edits are recommended:

1. Specify `d_nu` continuous and nonnegative for the classical formulation. Alternatively, require local integrability and state a mild solution with the PDE holding almost everywhere. Nonnegativity alone does not ensure the defining integral exists.
2. State explicitly that decay in a full inhomogeneous periodic `H^k` norm requires mean-zero data. Without mean zero, the same decay holds for the positive-order derivative norms or homogeneous seminorms, not for an undifferentiated constant component. This is a clarification of the note's wording, not a failure of its Parseval argument.

For later use, the forcing Duhamel term requires the same mean convention if one wants to apply a spectral-gap estimate to it. A phase mean cannot be assigned the `e^{-D}` factor.

## Return example and verification scope

The unequal-damping quotient is

\[
v'=b+(d_1-d_2)v-av^2,
\]

so equal damping cancels exactly. The logarithmic amplitude identity is also correct wherever `T` does not vanish.

On the displayed interval, `cos t+sin t` is positive. The explicit pair solves the stated damped rotation system and reaches `Omega=0` at `t=pi/4`, retaining

\[
T=\sqrt2e^{-\pi/2}=0.2939861162\ldots<1.
\]

The corresponding undamped retained amplitude is `sqrt(2)>1`. The script checks the formulas and ODE residual symbolically; its final strict-amplitude inequality is checked numerically. For an entirely exact account of that inequality, square it and use `e^pi>1+pi>2`. Then `2e^{-pi}<1`.

The phrase “contracting cycle” should be narrowed to “contracting reset pulse” or “contracting return step.” The example is one finite return pulse, not a closed recurrent orbit or a proven periodic cycle. It establishes exactly the intended falsifier: a success test that records only `Omega=0` can miss lost amplitude.

The script's phase-dependent coefficient counterexample correctly rejects the commutation assumption when `B` depends on phase. It is a commutator example, not specifically a zero-mode example despite the internal variable name. That naming has no effect on the result.

## Scope of the proposed PDE gate

The proposed net logarithmic inequality is a reasonable **design requirement for a chosen principal component or harmonic**. It must not be promoted from an upper propagator bound to a lower growth theorem. A future formulation should define `g_m` as the actual inviscid logarithmic gain of the retained component, or provide a separate lower propagator/cone estimate. It also needs the effect of profile smoothing at the measuring point, residual forcing, lower-order diffusion, inherited geometry, and later holding losses.

For a general periodic profile, a first-harmonic spectral upper bound is enough to demonstrate damping or refute a claimed unconditional growth inference. It is not enough to certify that a particular physical circulation gradient grows. The note presently labels its criterion as a design requirement rather than a sufficient blowup theorem, which is appropriate.

## Smooth forcing is an available route

The archival comparison must not create an artificial requirement that the new construction be unforced. Fefferman's official statements C and D explicitly permit a smooth prescribed force for ordinary three-dimensional Navier–Stokes with positive viscosity. The whole-space force must satisfy the stated rapid decay of every mixed derivative in space and time; the periodic version requires spatial periodicity and the corresponding time decay. [Official problem statement, pp. 1–2](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf), checked 2026-09-08.

Therefore a controlled rotation, activated seed, or correction force is not excluded merely because it is forced. The relevant obligation is that the **total physical force** satisfy the correct admissibility conditions through the proposed singularity and on the required time domain, while the original positive-viscosity equation holds. Forcing can compensate damping in the auxiliary equation; whether its complete physical derivatives remain admissible is a separate hard estimate. The present lemma proves no impossibility for forced C/D constructions.

Our earlier mechanism-transfer note described the missing controls of the archive's unforced orbit. That comparison should be read as a description of the existing route, not a restriction on the newly authorized research. The constructive target may properly be a smooth-forced NS stage with gain, a valid return, summable complete residuals, and ultimately one admissible singular solution.

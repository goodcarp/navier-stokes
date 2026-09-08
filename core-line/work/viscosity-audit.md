# Viscosity audit of the published forced Euler construction

Date: 2026-09-08. Scope: a read-only mathematical comparison with the downloaded 112-page *Blowup for the Euler Equations with Smooth Forcing*. Page numbers below are the manuscript's printed page numbers. This note does not independently establish its Euler theorem, build the Lean project, or make claims about the authors' separate unreleased hypodissipative construction.

**Result:** the exact published velocity cannot be transferred to ordinary Navier–Stokes with admissible smooth forcing. Moreover, in the explicit auxiliary principal-profile model, retaining the published infinite frequency schedule makes every fixed positive order of dissipation eventually dominate the published growth rate. This second conclusion is a statement about a precisely specified model and schedule, not a no-go theorem for other constructions.

## 1. An obstruction for the actual velocity, independent of asymptotics

Theorem 1.1, p. 2, places the circulation `Gamma = r u_phi` in a fixed solid torus with `0 < R < r0/2` and bounds `Gamma`, `u_r`, and `u_z` uniformly on `[0,T*)`. Consequently

\[
 \|u_\phi(t)\|_\infty
 \le (r_0-R)^{-1}\|\Gamma(t)\|_\infty.
\]

Thus the full velocity is uniformly bounded. Proposition 13.3, pp. 107–108, explicitly gives the compact support of the full velocity, this swirl inequality, and the resulting energy bound. The vorticity and circulation gradient blow up while velocity amplitude remains bounded.

Suppose this same velocity were an ordinary Navier–Stokes solution with fixed `nu > 0`, smooth initial data, and a smooth compactly supported force `f_NS` through `T*`. For a smooth preterminal solution, testing

\[
 u_t+(u\cdot\nabla)u+\nabla p=\nu\Delta u+f_{NS}
\]

against `-Delta u`, and using incompressibility to remove pressure, gives

\[
 \frac{d}{dt}\|\nabla u\|_2^2+\nu\|\Delta u\|_2^2
 \le \frac{C}{\nu}\|u\|_\infty^2\|\nabla u\|_2^2
     +\frac{C}{\nu}\|f_{NS}\|_2^2.
\]

The nonlinear estimate is `|<u·grad u, Delta u>| <= ||u||_infty ||grad u||_2 ||Delta u||_2`; Young's inequality absorbs the last factor. On the finite interval all coefficients on the right are integrable. Grönwall bounds `||grad u||_2` uniformly, while the energy estimate bounds `||u||_2`. Standard local strong well-posedness in `H^1` then continues the solution through `T*`, contradicting the claimed vorticity breakdown. This is also a direct instance of the velocity continuation criterion summarized in Fefferman's Clay statement, p. 3.

Allowing a different pressure does not circumvent this argument. The naive residual prescription `f_NS = f_Euler - nu Delta u` cannot have the required smoothness through blowup. Arbitrary gradient changes of pressure/force do not remove the continuation obstruction.

## 2. A separate geometry obstruction

The manuscript's concentration occurs on the axisymmetric lift of a moving meridional point, hence on a circle a fixed distance from the axis: equation (4.2), p. 11, and Theorem 1.1. For ordinary Navier–Stokes in the suitable finite-energy class with smooth forcing, Caffarelli–Kohn–Nirenberg partial regularity excludes an off-axis axisymmetric singularity. A singular point off the axis rotates into a whole singular circle, whose one-dimensional parabolic Hausdorff measure is positive; CKN gives zero measure to the singular set. See the primary author exposition linked below, section 3.2.

This applies to a candidate classical solution from smooth data in the standard energy setting: energy estimates permit continuation as a suitable weak solution, and weak–strong uniqueness identifies it with the preterminal classical solution. It is not a claim that every distributional NS solution is regular away from the axis without suitability or energy hypotheses.

Therefore an ordinary-NS extension cannot preserve the fixed off-axis axisymmetric location of breakdown. Concentration would need to approach the axis, or exact axisymmetry would need to be abandoned. Merely increasing amplitude on the same fixed ring does not evade this separate obstruction. This does not exclude off-axis forcing whose resulting solution later concentrates on the axis.

## 3. Exact ordinary-diffusion operators and the principal symbol

Use the paper's physical volume coordinates, equation (2.1), p. 3:

\[
 y_1=z,\qquad y_2=r^2/2,\qquad
 \xi=\omega_\phi/r,\qquad \Gamma=ru_\phi.
\]

The usual axisymmetric NS circulation diffusion is

\[
 \nu(\partial_r^2-r^{-1}\partial_r+\partial_z^2)\Gamma,
\]

and reduced-vorticity diffusion is

\[
 \nu(\partial_r^2+3r^{-1}\partial_r+\partial_z^2)\xi.
\]

Since `partial_r = r partial_y2`, these become respectively

\[
 \nu D\Gamma,\qquad \nu(D+4\partial_{y_2})\xi,
 \qquad D=\partial_{y_1}^2+2y_2\partial_{y_2}^2.
\]

For the material phase `s_m = N_m p_m·A_m(y,t)` from (3.2), p. 7, with `zeta_m = H_m^T p_m`, the common second-order principal phase symbol is

\[
 q_m=\zeta_{m,1}^2+r^2\zeta_{m,2}^2=r^2\kappa_m,
 \qquad \kappa_m=r^{-2}\zeta_{m,1}^2+\zeta_{m,2}^2.
\]

Both principal equations thus acquire `nu N_m^2 q_m partial_s^2`. The first-order `4 partial_y2` term, derivatives of material amplitudes, localization, and varying coefficients remain additional terms. This is why adding a constant damping term is a principal approximation rather than the full NS equation.

Section 12 uses the magnified coordinates of (2.10), p. 5. There the metric factor is `c_g = 1 + 2 epsilon_geom x_2` and the corresponding symbol is `q_m' = zeta_1^2 + c_g zeta_2^2`; the effective viscosity is `nu_eff = nu/epsilon_geom^2`. The final fixed physical dilation gives only another fixed factor. None changes the frequency exponents below.

On the growth intervals this symbol has a uniform positive lower bound. Section 9, pp. 48–49, fixes a containing region with `0<c_-<=c_g<=c_+` and fixed covector intervals `0<z_-<|zeta_m|<z_+`. Lemma 9.2, pp. 60–61, recovers those intervals on growth/transition/return. Thus at the material center

\[
 q_m'\ge q_*:=\min(1,c_-)z_-^2>0.
\]

The paper's young-label estimates extend such bounds onto the controlled material tubes, but the central lower bound already suffices for the auxiliary central-profile calculation here. An off-axis radius tending to zero or uncontrolled phase compression would invalidate the fixed `q_*` hypothesis; those are changes of scheme.

## 4. Exact statement of the auxiliary model

The Euler principal amplitude matrix in (3.4), p. 7, is

\[
 B_m=\begin{pmatrix}0&\widehat a_m\\\widehat b_m&0\end{pmatrix}.
\]

Freeze the material label, retain the supplied Euler coefficient history `B_m(t)` and positive phase symbol `q_m(t)`, and consider a periodic profile equation with common scalar diffusion. An ordinary-diffusion Fourier mode `k != 0` then obeys

\[
 R_k'=[B_m(t)-\nu_{eff}N_m^2q_m(t)k^2 I]R_k.
\]

At constant coefficients its eigenvalues are `+-sqrt(a_hat b_hat) - nu_eff N_m^2 q_m k^2`. For time-varying coefficients the scalar damping commutes exactly with the matrix evolution, so multiplying the Euler solution by `exp(-integral damping)` is exact for this auxiliary model.

The actual periodic profile `F` equals `s` near zero and is not a sine; its second derivative cannot be replaced globally by `-F`. The Fourier-mode statement concerns each nonzero mode. The common phase heat semigroup treats their full superposition. Local vanishing of `F''` on its linear cell does not imply that global diffusion vanishes or remains absent there during a growth interval.

For any common scalar damping `d(t)`, the projective variable `v=c0 Omega_hat/T` obeys

\[
 v'=c_0\widehat b-(\widehat a/c_0)v^2,
 \qquad (\log|T|)'=(\widehat a/c_0)v-d.
\]

The first equation is unchanged, while the second loses the full damping. Thus the return test (7.14) can look successful while amplitude amplification has failed. The Euler holding conclusion in Lemma 7.4, p. 32, changes from constant `T` to decaying `T` in this model.

For the quantitative comparison below, generalized dissipation means `nu (-Delta)^alpha`, with ordinary viscosity at `alpha=1`. The frozen high-frequency model has damping `nu_eff N_m^(2alpha) q_m^alpha |k|^(2alpha)`, where the magnification factor is now `nu_eff = nu/epsilon_geom^(2alpha)`, up to the other fixed physical dilation. For noninteger alpha this is a principal-symbol model; no claim is made here that the full nonlocal operator preserves the localization or existing correction scheme.

## 5. Quantitative falsifier for the published infinite schedule

Equation (12.1), p. 95, fixes `beta=1/8` and `bar_beta=7/8`. Equations (12.3)–(12.4), pp. 95–96, prescribe

\[
 N_m=N_{m-1}^{Q_m},\quad Q_m=q_0+m,\quad
 \Lambda_m=(k_m+e-\bar\beta)\log N_m,\quad k_m\le Q_m/c_Q,
\]

\[
 s_m=\Lambda_m\sigma_{m-2}/\sigma_{m-1},\qquad
 \Gamma_m^{gr}=\sigma_{m-1}\sin s_m,\qquad
 \sigma_j\asymp N_j^{\beta/2}.
\]

The last clock bounds are the prescribed and realized scalar windows, discussed in section 12.2, p. 97. Lemma 9.2, p. 61, supplies the logarithmic amplitude growth `g_m=(1+O(epsilon_s)) Gamma_m^gr` on growth. The same upper bound is used explicitly in (11.17), p. 89. The scalar windows are uniform in the layer.

Using only `sin s <= s`, for `m>=3`,

\[
 g_m\le C\Lambda_m\sigma_{m-2}
 \le C\Lambda_m N_m^{\beta/(2Q_mQ_{m-1})}.
\]

Every nonzero principal Fourier mode has damping at least `d_m >= nu_eff q_*^alpha N_m^(2alpha)`. Hence

\[
 \boxed{\displaystyle
 \frac{d_m}{g_m}\ge
 \frac{c_\alpha\nu_{eff}}{\Lambda_m}
 N_m^{\,2\alpha-\beta/(2Q_mQ_{m-1})}.}
\]

For ordinary viscosity the exponent is exactly

\[
 2-\frac{1}{16Q_mQ_{m-1}}.
\]

For **every fixed `alpha>0`**, `Q_mQ_{m-1}->infinity`, so eventually the exponent is at least `alpha`. Also `x_m=log N_m=y product_(j=2)^m Q_j` by (12.5), p. 96, while

\[
 \log\Lambda_m\le \log(Q_m/c_Q+8)+\log x_m=o(x_m).
\]

It follows that `d_m/g_m -> infinity` for every fixed `nu_eff>0` and `alpha>0`, within this same published schedule and auxiliary model. The logarithmic amplitude derivative becomes negative on the relevant growing line for all sufficiently large layers. The earliest destroyed step is the positive logarithmic growth and target hit in Lemma 9.2, p. 61, before the return-root argument or all-order force summability.

The Euler target needs logarithmic gain `Lambda_m`. Its growth duration is comparable to `Lambda_m/Gamma_m^gr`, equation (12.18), p. 104. At this duration the dissipative exponent divided by the proposed Euler gain is again comparable to `d_m/Gamma_m^gr`, which diverges. A small but fixed viscosity can admit a finite number of favorable stages; it cannot preserve this infinite schedule.

This is **not** a contradiction of an unreleased hypodissipative claim. That construction can use different scales, symbols, amplitudes, or a different leading balance. The statement proved here is deliberately about the existing Euler schedule under common principal damping.

## 6. What conditional window remains?

For a fixed finite layer with supplied positive growth `g_m` and phase symbol `q_m`, a necessary growing-mode gate is

\[
 \nu_{eff}N_m^{2\alpha}q_m^\alpha |k|^{2\alpha}<g_m.
\]

It is sufficient for exponential growth only in the constant-coefficient two-by-two model along its growing eigenvector. For slowly varying coefficients the exact propagator gain must exceed the integrated damping. Neither version establishes a localized nonlinear PDE construction.

Any serious endpoint redesign must obtain a nonempty compatible window between the lower frequency required by the correction estimates and the upper frequency imposed by this viscous gate. In the published schedule that window eventually closes. Ordinary NS additionally requires escaping both the bounded-velocity and off-axis-axisymmetry obstructions in sections 1–2; changing only the frequency ladder is insufficient.

## Verification and primary sources

`viscosity-check.py` uses SymPy to verify the coordinate-transformed differential operators, phase-symbol identity, damped characteristic polynomial, Riccati cancellation, and exact two-generation frequency exponent. It does not verify the asymptotic limit by numerical sampling; the limit is established by the displayed inequalities. It imports no Lean or downloaded project code.

- [Published Euler manuscript](https://cims.nyu.edu/~tristanb/euler.pdf), especially pp. 2, 3, 5, 7, 32–36, 48–49, 60–63, 89, 95–97, 104, 107–108.
- [Euler formal theorem statement](https://github.com/tristanbuckmaster/fluid_lean/blob/main/euler-blowup/Challenge.lean). Its exported theorem does not state the bounded-velocity feature used here; that feature is in the paper's stronger construction.
- [Fefferman's official Clay statement](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf), pp. 1–4, for admissible forcing, NS velocity breakdown, and CKN partial regularity.
- [Characterization and Regularity for Axisymmetric Solenoidal Vector Fields with Application to Navier–Stokes Equation](https://www.cscamm.umd.edu/publications/AxisymmetricFlow_CS-08-45.pdf), section 3.2, for the off-axis singularity exclusion in the standard axisymmetric NS setting.

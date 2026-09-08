# Forced Euler → Navier–Stokes: the first transfer gate

September 8, 2026. **Status: the direct transfer fails; a redesigned viscous stage remains open.**

**Latest:** [The ninth-pass report](pass9/REPORT.md) adds the leading seed's full-pressure initial deceleration, a resolved finite pulse in a nonlinear planar model, an exact affine-envelope finite-gain theorem and an H4 full-field endpoint-validation route. The selected next initial field now has explicit rational amplitude and an audited strict core cone. No compact 3D persistence interval, inherited return or blowup has been proved; [STATUS.md](STATUS.md) records the current state.

**Continuation:** [The second-pass report](pass2/REPORT.md) adds an explicit compact active-core local lemma for actual NS, a published viscous comparison family, and sharper correction/mean-feedback conditions. The first-pass analysis below remains valid. No NS blowup proof has been obtained.

This pass extracted the relevant construction mechanism, compared it with our PASS14/PASS28 gaps, derived a viscosity test, and checked the finite algebra independently. It did not construct a singular Navier–Stokes solution or independently replay the released Lean proof. None of the elementary lemmas below is claimed to be new to the literature.

## What the release gives us

The Euler manuscript uses a controlled return: rotation steers a wave covector, a pulse resets one principal vorticity amplitude, and inherited geometry plus separate gradient estimates prepare the next stage. A causal selection and one scale sequence assemble infinitely many stages; corrections control all mixed derivatives of the physical force. These are the useful ingredients to compare with our missing same-solution regeneration. See §§7, 11–12 and Proposition 12.3 of the [Euler manuscript](https://cims.nyu.edu/~tristanb/euler.pdf).

Our PASS14 gives large finite amplification across different initial data. PASS28 gives one reinforcing episode with an initially irrotational core. Neither gives repeated amplification on one fixed solution. This release supplies a relevant mechanism to study, but no existing gap in our Navier–Stokes arguments is closed by citation alone.

The main strategic change is to require a **frequency window compatible with viscosity before developing an infinite correction hierarchy**.

## 1. The actual Euler velocity cannot be reused unchanged

The source construction bounds circulation and meridional velocity and keeps the circulation away from the axis. Thus \(u^\phi=\Gamma/r\) is bounded as well, even though vorticity diverges.

For ordinary Navier–Stokes with fixed \(\nu>0\), test the momentum equation against \(-\Delta u\). Young's inequality gives

\[
\frac{d}{dt}\|\nabla u\|_2^2+\nu\|\Delta u\|_2^2
\le \frac C\nu\|u\|_\infty^2\|\nabla u\|_2^2+
\frac C\nu\|f\|_2^2.
\]

For bounded velocity and an admissible smooth force on a finite interval, Grönwall and the energy estimate give a uniform \(H^1\) bound. Local strong theory then continues the solution. Consequently the same velocity cannot solve ordinary NS with such a force through the proposed singular time. Changing pressure cannot evade this argument.

There is a separate geometric restriction: a suitable axisymmetric NS solution cannot have an off-axis singular point. Rotation would produce a singular circle, contradicting Caffarelli–Kohn–Nirenberg's zero one-dimensional parabolic Hausdorff measure conclusion. A redesign must permit unbounded velocity and, if it retains axisymmetry, move concentration toward the axis. See [Fefferman's statement](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf), pp. 1–4, and [Liu–Wang](https://www.cscamm.umd.edu/publications/AxisymmetricFlow_CS-08-45.pdf), §3.2.

Smooth forcing remains permitted in Clay alternatives C/D. Our investigation imposes no additional unforced requirement.

## 2. An exact auxiliary lemma: reset success can conceal amplitude loss

In volume coordinates \(y=(z,r^2/2)\), direct differentiation gives the viscous operators

\[
\nu\mathcal D\Gamma,\qquad \nu(\mathcal D+4\partial_{y_2})\xi,
\quad \mathcal D=\partial_{y_1}^2+r^2\partial_{y_2}^2,
\quad \xi=\omega^\phi/r.
\]

Their principal phase symbol is the same:
\[
q=\zeta_1^2+r^2\zeta_2^2>0.
\]

Freeze the material-label dependence into prescribed time coefficients. For mean-zero periodic phase fields, normalize the pair as
\[
W=\left(\gamma,\frac{\sigma}{N}\partial_\theta^{-1}\xi\right)^T.
\]
The leading system has the form
\[
W_t=B(t)W+d(t)W_{\theta\theta},\qquad d(t)=\nu N^2q(t).
\]

**Lemma.** If \(B,d\) are continuous, independent of phase, and \(d\ge0\), let \(\Phi(t,s)\) solve the component equation and set \(D(t,s)=\int_s^t d\). Then
\[
W(t)=\Phi(t,s)e^{D(t,s)\partial_\theta^2}W(s).
\]

Matrix multiplication and the phase heat operator commute. Differentiating this formula proves the equation and initial condition. Matrices at different times need not commute: their time ordering stays inside \(\Phi\). Fourier mode \(n\) is multiplied by \(e^{-n^2D}\), so Parseval gives, for mean-zero data,
\[
\|\partial_\theta^jW(t)\|_2
\le \|\Phi(t,s)\|e^{-D(t,s)}\|\partial_\theta^jW(s)\|_2.
\]
For \(j\ge1\), differentiation removes the need for a mean-zero assumption. With an added residual \(h\), the solution also contains
\[
\int_s^t\Phi(t,\tau)e^{D(t,\tau)\partial_\theta^2}h(\tau)\,d\tau.
\]

This handles a general periodic carrier, not just a sine. It is an exact auxiliary-model result, **not** a full localized NS reduction. Phase-dependent coefficients, material derivatives, different lower-order viscous terms and nonlinear corrections require separate estimates. Heat also changes an initially linear phase cell; a pointwise gradient lower bound must be rebuilt.

For a single harmonic the pair reads
\[
T'=a\Omega-dT,\qquad \Omega'=bT-d\Omega.
\]
Where \(T\ne0\),
\[
(\Omega/T)'=b-a(\Omega/T)^2,
\qquad (\log|T|)'=a\Omega/T-d.
\]
The return quotient does not see common damping. The retained amplitude does.

**Exact negative control.** Take \(a=1,b=-1,d=2\), starting at \((T,\Omega)=(1,1)\). Then
\[
(T,\Omega)=e^{-2t}(\cos t+\sin t,\cos t-\sin t).
\]
At \(t=\pi/4\), the reset is perfect: \(\Omega=0\). But \(T=\sqrt2e^{-\pi/2}\approx0.294\), compared with \(\sqrt2\approx1.414\) without damping. The strict loss follows from \(e^\pi>1+\pi>2\). A return-only test would accept a contracting pulse.

## 3. The published scale schedule fails the principal damping test

Using the source's equations (12.1)–(12.4), its fixed \(\beta=1/8\), and its growth-clock bounds, the two-generation conversion for \(m\ge3\) is
\[
N_m=N_{m-1}^{Q_m},\quad Q_m=q_0+m,
\qquad g_m\le C\Lambda_mN_m^{1/(16Q_mQ_{m-1})}.
\]
The controlled phase symbol has a fixed positive lower bound on these growth windows.

For dissipation \(\nu(-\Delta)^\alpha\), ordinary viscosity being \(\alpha=1\), the principal Fourier damping therefore satisfies
\[
\boxed{\quad
\frac{d_m}{g_m}\ge
\frac{c_\alpha\nu_{\rm eff}}{\Lambda_m}
N_m^{2\alpha-1/(16Q_mQ_{m-1})}\longrightarrow\infty.
\quad}
\]

Indeed \(Q_mQ_{m-1}\to\infty\), while \(\log\Lambda_m=o(\log N_m)\). For every fixed \(\alpha>0\), the positive power eventually dominates the logarithmic factor. The source's positive-growth/target-hit step is lost in this model before the return or force summability steps can help.

**Scope:** this retains the Euler coefficient history and its scale schedule while adding principal damping. It rules out that unchanged transfer mechanism. It neither analyzes a self-consistent redesigned viscous flow nor contradicts the authors' separate hypodissipative claim. For noninteger \(\alpha\), localization of the full nonlocal operator is an additional issue.

The practical gate is
\[
N_{\rm correction}(m,k)\le N_m\le N_{\rm diffusion}(m).
\]
The left side buys accuracy; the right side prevents diffusion from overwhelming amplification. Both must hold for the orders needed at each stage.

## 4. A candidate redesign survives elementary scaling only

Merely shrinking a controlled thin ring is insufficient. If its radius is \(R\), older-gradient scale \(L\ll R\), carrier wavelength \(h=L/M\), and \(|\nabla\Gamma|\lesssim G_0/L\), the frozen growth calculation gives
\[
\lambda\lesssim G_0/(R^{3/2}L^{1/2}),\qquad
d_\nu\asymp\nu M^2/L^2.
\]
Growth therefore requires
\[
G_0/\nu\gtrsim M^2(R/L)^{3/2}.
\]
At fixed aspect ratio, ring radius cancels. This is a conditional controlled-profile obstruction, not a theorem against all near-axis concentration.

A broader candidate is a nonaxisymmetric nested core. Set \(\rho=L/L_0\) and choose
\[
\mathrm{Re}=\mathrm{Re}_0\rho^{-a},\quad
M=M_0\rho^{-b},\qquad 0<2b<a<1/2.
\]
With \(U=\nu\mathrm{Re}/L\), the optimistic strain/diffusion comparison, energy and stage time scale as
\[
\frac{\nu/h^2}{U/L}\asymp\rho^{a-2b},\quad
E\asymp\rho^{1-2a},\quad
\tau\asymp\rho^{2+a}.
\]
For \(a=1/4,b=1/16\), these powers are respectively \(1/8,1/2,9/4\). Core dissipation over a stage scales as \(\rho^{3/4}\); even the estimate with a comparable-amplitude fine carrier gives the summable power \(5/8\). These budgets do not immediately forbid a cascade on geometric scales.

They do **not** provide an amplifying solution. A comparable-amplitude carrier also introduces potentially larger nonlinear terms; small amplitude or an actual cancellation must be proved. The core must generate the right strain, survive its own evolution, and pass a usable configuration to the next stage.

## 5. The physical-force condition remains severe

Let \(\mathcal R_m\) be the normalized physical vector residual increment of the **assembled** field after choosing pressure, including localization, joining and interactions with older stages. On a fixed stage,
\[
f_m=(U_m^2/L_m)\mathcal R_m.
\]
If \(e_{m,k,j}=\|\partial_y^k\partial_s^j\mathcal R_m\|_\infty\), a direct sufficient summability requirement is
\[
\sum_m\rho_m^{-[3+2a+k+(2+a)j]}e_{m,k,j}<\infty
\qquad\text{for every fixed }k,j.
\]
Carrier derivative costs must be included in \(e\). Smooth temporal joins and a common compact support are also needed. Fixed numerical accuracy on isolated stages does not establish this condition. Solving a principal equation exactly removes only one residual contribution.

Nor is time-dependent shrinking a symmetry. For \(u=a_0(t)V(x/L(t),s(t))\), \(s'=a_0/L\), the normalized residual includes
\[
V_s+(V\cdot\nabla)V+\nabla P-\frac{\nu}{a_0L}\Delta V
+\frac{La_0'}{a_0^2}V-\frac{L'}{a_0}(y\cdot\nabla)V.
\]
The last two terms must be solved or controlled, not discarded.

## Outcome and next dependency

This pass rejects direct reuse of the velocity, fixed ring geometry and published Euler scale schedule for ordinary NS. It establishes useful auxiliary-model calculations and identifies a scaling candidate that is not immediately excluded by elementary budgets.

The next target is **one viscous stage with positive retained gain, a nonempty correction/diffusion frequency window, and an inherited terminal geometry**. Its precise obligations and rejection conditions are in [NEXT_TARGET.md](NEXT_TARGET.md). It is not yet proved.

## Verification and provenance

Three independently prepared symbolic programs were run successfully. Checks cover coordinate operators, phase normalization, damping, a noncommuting time-dependent matrix example, the exact failed-amplification pulse, and scale exponents. Independent analytical reviews found no blocking error; their qualifications about phase-profile distortion, force cross-interactions and nonlinear carrier size are incorporated above.

These are finite algebra checks and mathematical reviews, not a full PDE or Lean verification. The current closing estimate in our original NS archive remains unproved.

The public proof source was pinned to commit `d0124689230b58b4f86e7b90ac59de06404b3b6b`. The exported Challenge and Solution propositions match textually after whitespace normalization. That is not type checking or a dependency audit; the exported proposition also does not certify every geometric or quantitative detail used in this comparison. See [VERIFICATION.md](VERIFICATION.md) and [SOURCES.json](SOURCES.json).

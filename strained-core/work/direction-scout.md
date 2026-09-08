# Direction scout: ordinary viscosity is a principal-scale constraint

Date: 2026-09-08. Scope: one necessary-condition calculation and one precisely posed replacement target. This is not a Navier–Stokes blowup construction.

## Finding

The released Euler construction's ability to choose the new carrier frequency arbitrarily large after fixing the older flow is incompatible with transferring its growing-wave argument unchanged to ordinary Navier–Stokes. Ordinary viscosity acts at order frequency squared, while the frozen Euler growth exponent of §3 is independent of the new carrier frequency. Moving the same thin-ring design toward the axis does not automatically fix this: under controlled dimensionless profiles, it scales the growth and the damping together.

A more plausible *scaling target*, though no existence lemma is proved here, is a nonaxisymmetric nested core with increasing local Reynolds number and unbounded velocity. Its first gate is an amplification estimate that includes viscosity at principal order and retains a nonempty interval of allowed carrier frequencies. That gate must precede importing the infinite-stage bookkeeping.

## 1. Exact reduced viscous operators and the frozen fast-wave test

Use the paper's variables

\[
y=(z,r^2/2),\qquad \Gamma=r u^\phi,\quad \xi=\omega^\phi/r,
\quad A=\operatorname{diag}(r^{-2},1),\quad W=r^{-4}.
\]

For ordinary Navier–Stokes, direct cylindrical differentiation gives

\[
D_v\Gamma=\nu\mathcal D\Gamma+r f^\phi,
\qquad
D_v\xi=\nu(\mathcal D+4\partial_{y_2})\xi
 +W\partial_{y_1}(\Gamma^2)+(\nabla\times f)^\phi/r,
\]

where

\[
\mathcal D=\partial_{y_1}^2+r^2\partial_{y_2}^2.
\]

Thus both scalar equations have the same second-order principal diffusion symbol. For the material phase \(s=Np\cdot A_m(y,t)\), write \(\zeta=H_m^Tp\). The paper uses

\[
\kappa=\zeta^TA\zeta,
\quad d=J\zeta\cdot\nabla_y\Gamma_{<m},
\quad c=2W\Gamma_{<m}.
\]

The viscous principal symbol is

\[
q=\zeta_1^2+r^2\zeta_2^2=r^2\kappa>0.
\]

**Carrier qualification.** The actual paper takes a smooth periodic carrier \(F\) with \(F(s)=s\) near zero, not a single sine wave. Viscosity introduces leading \(F''\) and \(F'''\) terms, so merely putting a scalar damping term into its original fixed-carrier ansatz does not close that ansatz. Instead Fourier-resolve the frozen fast-phase linearization. On its nonzero harmonic \(n\), the two-amplitude system, after normalizing the vorticity amplitude by \(nN\), has inviscid eigenvalues

\[
\pm\lambda,\qquad \lambda^2=-\frac{dc\zeta_1}{\kappa},
\]

and viscous eigenvalues

\[
-\nu n^2N^2q\pm\lambda.
\]

This is an exact calculation for the frozen-coefficient linear fast-wave system. It is a necessary test for retaining that growing-wave mechanism, not a theorem excluding other nonlinear PDE mechanisms or all transient growth.

If \(\lambda>0\), a growing harmonic requires

\[
N^2<\frac{\lambda}{\nu n^2q}.
\]

For the leading material system with time-dependent coefficients, the scalar damping commutes with the matrix propagator:

\[
P_{\nu,n}(t,s)=
\exp\left[-\nu n^2N^2\int_s^tq(\tau)\,d\tau\right]P_{E,n}(t,s).
\]

Once the older fields and the interval are fixed, \(q\) is positive and the right side's Euler factor is independent of \(N\) in the paper's normalization. Sending \(N\to\infty\) therefore destroys this mechanism. A viscous replacement must prove a quantitative overlap

\[
N_{\rm correction}(m,k)\le N_m\le N_{\rm diffusion}(m)
\]

at each stage and for every correction order eventually needed. The original freedom to take the carrier “sufficiently large” cannot simply be inherited.

## 2. Why moving the thin ring toward the axis is insufficient by itself

This section is conditional on a controlled-profile model, not a no-go theorem for all axisymmetric blowup.

Let the ring radius be \(R\), the physical scale of the older circulation gradient be \(L\ll R\), and the new physical carrier wavelength be \(h=L/M\), where \(M\gg1\) is scale separation. Suppose

\[
|\Gamma|\le G_0,\qquad |\nabla_{z,r}\Gamma|\le C G_0/L.
\]

Writing \(k=(\zeta_1,r\zeta_2)\), the Euler growth exponent obeys

\[
\lambda^2
=-\frac{2\Gamma}{r^3}
\frac{k_z^2\Gamma_r-k_zk_r\Gamma_z}{|k|^2},
\quad
\lambda^2\le \frac{2|\Gamma|\,|\nabla_{z,r}\Gamma|}{r^3}.
\]

Consequently, in a thin annulus with \(r\asymp R\),

\[
\lambda\lesssim\frac{G_0}{R^{3/2}L^{1/2}},
\qquad D_{\nu}\asymp\frac{\nu M^2}{L^2}.
\]

Positive frozen growth therefore requires, up to controlled constants,

\[
\boxed{\quad\frac{G_0}{\nu}\gtrsim M^2(R/L)^{3/2}.\quad}
\]

If the aspect ratio \(L/R=\delta\) is held fixed, taking \(R\to0\) cancels out of this inequality. With bounded circulation and fixed \(\delta\), \(M\to\infty\) remains impossible in this controlled-profile mechanism. If the aspect ratio gets thinner, the requirement worsens. A much larger circulation gradient than \(G_0/L\) introduces an additional scale and is a substantive change requiring a new calculation.

Near-axis smoothness also imposes parity/pole conditions: preterminal smooth axisymmetric swirl has \(u^\phi=O(r)\) at the axis and hence \(\Gamma=O(r^2)\) there. Annuli shrinking toward the axis can have constants that deteriorate in time, but replacing the paper's fixed positive-radius coordinate estimates requires explicitly controlling that deterioration. One may not reuse their geometry constants uniformly.

## 3. A scaling-compatible replacement target

Consider nonaxisymmetric cores of size \(L_m\to0\), velocity \(U_m\), and local Reynolds number

\[
\mathrm{Re}_m=U_m L_m/\nu.
\]

For a core strain of order \(U_m/L_m\) feeding a carrier with separation \(M_m=L_m/h_m\), its optimistic growth-versus-diffusion gate is

\[
\mathrm{Re}_m\gtrsim M_m^2.
\]

This is a dimensional necessary condition for such a strain-driven wave, not proof that the strain has the required orientation, duration or nonlinear feedback. If the fine carrier itself has amplitude comparable to \(U_m\), its nonlinear turnover rate may be \(M_m\) times the assumed core strain. A principal linear-growth use therefore additionally requires a small carrier amplitude or a proved cancellation; the energy bookkeeping below does not supply that cancellation.

Fix a reference length \(L_0\), set \(\rho_m=L_m/L_0\), and choose dimensionally consistent candidate scales

\[
\mathrm{Re}_m=\mathrm{Re}_0\rho_m^{-\alpha},
\quad U_m=\frac{\nu\mathrm{Re}_0}{L_0}\rho_m^{-1-\alpha},
\quad M_m=M_0\rho_m^{-\beta},
\qquad 0<2\beta<\alpha<\tfrac12.
\]

Then diffusion relative to core strain scales as

\[
\frac{\nu/h_m^2}{U_m/L_m}
=\frac{M_0^2}{\mathrm{Re}_0}\rho_m^{\alpha-2\beta}\to0.
\]

The characteristic stage lifetime and core energy are

\[
\tau_m\asymp L_m/U_m
=\frac{L_0^2}{\nu\mathrm{Re}_0}\rho_m^{2+\alpha},
\qquad E_m\asymp U_m^2L_m^3
=\nu^2\mathrm{Re}_0^2L_0\rho_m^{1-2\alpha}.
\]

For a geometric sequence of \(\rho_m\), these times and energies are summable. Under a core-gradient estimate \(|\nabla u|\asymp U_m/L_m\), dissipation over a stage is

\[
\nu U_m^2L_m\tau_m\asymp
\nu^2\mathrm{Re}_0L_0\rho_m^{1-\alpha},
\]

also summable. If the fine carrier itself has amplitude comparable to \(U_m\), this estimate acquires \(M_m^2\); its exponent becomes \(1-\alpha-2\beta>0\) under the chosen range. Thus the *scale budget* is not ruled out by these basic energy estimates. No PDE orbit satisfying it is supplied.

For a concrete rational exponent example, \(\alpha=1/4\), \(\beta=1/16\) gives

\[
U_m\sim\rho_m^{-5/4},\quad
\tau_m\sim\rho_m^{9/4},\quad
E_m\sim\rho_m^{1/2},\quad
D_\nu/\mathrm{strain}\sim\rho_m^{1/8}.
\]

This is a Type-II dimensional target, not a time-dependent rescaling of a known solution. Breaking axisymmetry removes the automatic ring of simultaneous singular points; it does not provide amplification or regeneration.

## 4. Smooth forcing is an exceptionally strong second gate

On an individual stage, nondimensionalize using fixed \(L_m,U_m,\tau_m=L_m/U_m\). Let \(\mathcal R_m\) be the normalized **physical vector equation residual increment of the assembled solution after a pressure is chosen**, including interaction with every retained older field, so that

\[
f_m=(U_m^2/L_m)\mathcal R_m.
\]

Any pressure-gradient part must first be removed or accounted for; a scalar/vorticity residual is not interchangeable with this physical force norm. If

\[
e_{m,k,j}=\|\partial_y^k\partial_s^j\mathcal R_m\|_\infty,
\]

then, up to fixed dimensional prefactors,

\[
\|\partial_x^k\partial_t^j f_m\|_\infty
\lesssim
\rho_m^{-[3+2\alpha+k+(2+\alpha)j]}e_{m,k,j}.
\]

The safe direct-series sufficient condition for a force smooth through the accumulation time is

\[
\boxed{\quad
\sum_m\rho_m^{-[3+2\alpha+k+(2+\alpha)j]}e_{m,k,j}<\infty
\quad\text{for every fixed }k,j.\quad}
\]

Smooth temporal activation and endpoint compatibility are needed as well. Growing dimensionless carrier frequencies are already included in \(e_{m,k,j}\), so this condition must not hide their derivative costs. This termwise condition is sufficient, not necessary for every possible cancellation scheme.

A fixed nonzero normalized profile residual scales like \(U_m^2/L_m\), which diverges. A fixed accuracy computation at each stage therefore cannot establish a smooth physical forcing construction. Corrections must improve faster than every required power of the shrinking scale. Solving the principal viscous equation exactly can remove one contribution, but remaining localization, gluing and cross-stage nonlinear residuals still require the same all-order physical bounds.

For comparison, a genuinely time-dependent change \(u(x,t)=a(t)V(x/L(t),s(t))\), \(s'=a/L\), gives the exact normalized residual

\[
V_s+(V\cdot\nabla)V+\nabla P-\frac{\nu}{aL}\Delta V
+\frac{La'}{a^2}V-\frac{L'}{a}(y\cdot\nabla)V.
\]

The last two terms are real PDE terms, not harmless relabeling. Unless canceled by a proved profile evolution (or a legitimate pressure gradient), they enter the physical force multiplied by \(a^2/L\). This excludes claiming a shrinking copy of the Euler construction solves ordinary Navier–Stokes.

## 5. Minimal next lemma, before an infinite cascade

**Viscous amplification with a frequency window.** On one exact smooth older Navier–Stokes stage, construct a localized divergence-free perturbation with its full viscous leading evolution, and prove simultaneously:

1. A prescribed circulation/strain or vorticity gain on an explicitly chosen carrier interval \(N_{\rm correction}\le N\le N_{\rm diffusion}\).
2. Controlled nonlinear feedback on the older flow and a terminal velocity/strain configuration usable by the next stage.
3. A bound on the physical force correction, including mixed derivatives, in a norm that will satisfy the scaled summability condition above.

The first numerical or analytic test should be whether the lower correction threshold fits below the upper diffusion threshold as local Reynolds number increases. Failure is informative and should stop that adaptation. Success on one stage would still leave the same-solution regeneration lemma and infinite-stage compatibility unproved.

## Sources and verification boundaries

- Alpöge/Buckmaster, *Blowup for the Euler Equations with Smooth Forcing*, local extraction `work/sources/euler.txt`, especially Theorem 1.1, equations (2.1)–(2.4), (2.10)–(2.15), and (3.1)–(3.4). Source URL supplied in the research context: https://cims.nyu.edu/~tristanb/euler.pdf . A fresh web-open attempt failed; the calculation used the locally available extraction and should be matched against the PDF before publication. The claims in this memo about viscosity, the thin-ring inequality and the alternative scales are our deductions, not claims made by the Euler paper.
- Liu–Wang, *Characterization and Regularity for Axisymmetric Solenoidal Vector Fields with Application to Navier–Stokes Equation*, equations (1.2)–(1.3) and §2 pole conditions: https://www.cscamm.umd.edu/publications/AxisymmetricFlow_CS-08-45.pdf . Successfully accessed 2026-09-08. Only the cylindrical viscous equations and smoothness/parity constraints are used from this source.

`python3 work/verify_direction_scout.py` passed symbolic checks of both reduced diffusion operators, the physical-coordinate growth identity, the frozen viscous characteristic polynomial, and the rational exponent example. These are algebra checks only. No Lean verification, Navier–Stokes numerical simulation, or claimed singular solution was performed for this scout.

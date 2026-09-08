# Jeong–Yoneda: an actual viscous amplification family and its transfer limits

Status: focused source extraction plus explicit scaling deductions, 8 September 2026. This is not an audit of the complete paper or a new blowup theorem. Primary source: Jeong–Yoneda, *Vortex stretching and enhanced dissipation for the incompressible 3D Navier–Stokes equations*, [arXiv:2001.02333v2](https://arxiv.org/pdf/2001.02333v2), posted 5 June 2020 (PDF title page dated 8 June). Local source: `work/pass2/jeong-yoneda.pdf`. Page references below are printed page numbers.

## 1. What the theorem actually supplies

**Published theorem.** Theorem 1.1, pp. 3–4, supplies smooth, unforced, globally smooth 2.5-dimensional Navier–Stokes solutions at viscosities \(\nu_n\downarrow0\), on
\[
\mathbb T_n^3=(\mathbb R/(2L_n\mathbb Z))^2\times(\mathbb R/(2\mathbb Z)),
\]
with controlled initial kinetic energy and large time-averaged enstrophy relative to that energy. For every fixed \(0<\bar a_0<1\), its displayed conclusion is an asymptotic lower bound of the form
\[
\nu_n^{\bar a_0}\,\frac1\delta\int_0^\delta\|\nabla u_n^{\nu_n}(t)\|_2^2\,dt
\gtrsim\|u_{n,0}\|_2^2.
\]
The source's statement uses a liminf; this notation here suppresses the asymptotic qualifier. When \(\bar a_0\le1/2\), Remark 1.2 permits \(L_n=1\), which is the branch used throughout this note. In that branch \(\|u_{n,0}\|_2^2\asymp1\). Section 1.3.4 also describes a small-exponent branch with uniformly bounded initial \(H^1\) velocity, so the effect is not merely large enstrophy already inserted at time zero.

This is an actual viscous PDE result. Its conclusion is an integrated norm bound for a family of data and viscosities. It does not supply a prescribed terminal profile, a lower bound at every point or time, or consecutive stages of one solution.

## 2. Geometric scales, viscosity and gain

**Source definitions.** Sections 1.5.2, 2.2 and 2.3.5 use
\[
\bar\ell=\ell/L=2^{-n},\qquad
\widetilde\ell=\ell^{1+c\delta},\qquad
\|\omega^L_{n,0}\|_\infty=1.
\]
Here \(L\) is the horizontal domain/base-vortex scale, \(\ell\) smooths the antiparallel base-vorticity configuration, and \(\widetilde\ell\ll\ell\) is the small-vortex scale. The paper allows the positive absolute constants denoted by \(c,C\) to change between estimates. The reference interval is \([0,\delta]\), with \(\delta>0\) fixed independently of \(n\).

Equation (2.37), p. 23, reads
\[
\boxed{\nu_n=\frac c{\delta^4}\,
\bar\ell^{\,4c_0\delta(1-C\delta)}
\ell^{\,4(1+C\delta)}L^{-2}.}
\]
This is the explicit small-viscosity choice used by the proof, independent of the small-vorticity amplitude parameter \(q\). It is a sufficient perturbative choice, not a proved physical saturation threshold or necessary upper limit for amplification. With \(L=1\) and \(\epsilon=\ell=\bar\ell\), write
\[
\nu_n=C_\delta\epsilon^K,\qquad
K=4(1+C\delta)+4c_0\delta(1-C\delta)=4+O(\delta).
\]
The constants are fixed after selecting the paper's sufficiently small \(\delta\); they are not numerical optimization data.

**Source mechanism.** Lemma 2.4, p. 14, proves for the *Euler base flow*, in its specified transported central region,
\[
\partial_1\eta_1(t,x)\ge
\exp\!\left[\frac2\pi t(1-\varepsilon_n-C\delta)\log\frac1\epsilon\right],
\qquad \varepsilon_n\to0.
\]
Small horizontal vorticity is stretched by this deformation. Section 2.3 then bounds the viscous–inviscid differences separately for the base velocity, passive velocity and both vorticities. Equation (2.33) demands that the time-averaged squared vorticity error be much smaller than the Euler small-vorticity squared norm. The accompanying displayed lower bound uses
\[
\mathcal E=\exp\!\left((1+C\delta)\|\partial_1u^L_{n,1}\|_{L^\infty_{t,x}}\right)
\asymp\epsilon^{-(2/\pi)(1+O(\delta))}
\]
and a squared-norm gain of order \(\mathcal E^{2\delta-c\delta^2}/\delta\), hence a positive power \(\epsilon^{-g_\delta}\), \(g_\delta=4\delta/\pi+O(\delta^2)>0\), after fixing sufficiently small \(\delta\). The precise constants should be retained from the original estimates in any future proof.

There is a notation slippage in that displayed line on p. 22: its left-hand side carries a viscous superscript although the immediately preceding target and following Cauchy formula concern the Euler solution. The safe use is the proof's stated Euler lower-bound/viscous-error comparison and Theorem 1.1's integrated viscous conclusion. We do not import a viscous Cauchy formula or a pointwise viscous lower bound.

**Consequence.** The established stage gain is polynomial in \(1/\epsilon\), equivalently a positive but small power of \(1/\nu_n\). This is useful evidence for one actual viscous stretching episode. An integral lower bound implies a large norm at some time, but does not identify a usable output state for a subsequent stage.

## 3. Exact conversion to fixed viscosity

**Scaling deduction, not a new existence theorem.** Take the \(L=1\) solutions, select target viscosity \(\nu>0\), and choose an isotropic domain scale \(\lambda_n>0\). In fixed reference units set
\[
A_n=\frac{\nu}{\nu_n\lambda_n},\qquad
U_n(x,t)=A_nu_n\left(\frac{x}{\lambda_n},\frac{A_n}{\lambda_n}t\right),\qquad
P_n(x,t)=A_n^2p_n\left(\frac{x}{\lambda_n},\frac{A_n}{\lambda_n}t\right).
\]
Then \((U_n,P_n)\) solves unforced NS with viscosity \(\nu\) on \((\mathbb R/(2\lambda_n\mathbb Z))^3\). Direct substitution gives the required relation \(A_n\lambda_n\nu_n=\nu\). Ordinary NS scaling at fixed viscosity alone would not change \(\nu_n\); the amplitude/time conversion in this formula is essential.

The reference interval becomes
\[
\tau_n=\frac{\delta\nu_n\lambda_n^2}{\nu},
\]
and the kinetic energy in **one shrinking periodic domain** becomes
\[
\|U_n(0)\|_2^2=A_n^2\lambda_n^3\|u_n(0)\|_2^2
\asymp\nu^2\nu_n^{-2}\lambda_n.
\]
Relative gradient-norm amplification is preserved because every squared gradient norm acquires the common factor \(A_n^2\lambda_n\), and a time average cancels the time conversion. Reynolds number based on the domain scale and order-one source base velocity is \(\mathrm{Re}_n\asymp A_n\lambda_n/\nu\asymp\nu_n^{-1}\).

Choose \(\lambda_n=\nu_n^p\), \(p>2\), in these reference units. Then
\[
E_n\asymp\lambda_n^{1-2/p},\quad
\tau_n\asymp\lambda_n^{2+1/p},\quad
\mathrm{Re}_n\asymp\lambda_n^{-1/p},\quad
A_n\asymp\lambda_n^{-(1+1/p)}.
\]
The initial ratio of domain scale to small-vortex scale is
\[
M_n=\widetilde\ell_n^{-1}
\asymp\lambda_n^{-(1+c\delta)/(pK)}.
\]
Thus the candidate parameters \(a=1/p\), \(b=(1+c\delta)/(pK)\) obey \(0<2b<a<1/2\) for sufficiently small fixed \(\delta\). In particular \(p=4\) gives
\[
a=1/4,\qquad b=1/16+O(\delta),\qquad
E_n\asymp\lambda_n^{1/2},\qquad
\tau_n\asymp\lambda_n^{9/4}.
\]
This explains a real connection to the earlier proposed exponent window; \(b=1/16\) is a small-\(\delta\) leading value, not an exact extracted exponent. The elementary initial scale ratio \(M_n^2/\mathrm{Re}_n\) tends to zero. It does not replace the paper's stronger viscous-error estimates, and \(M_n\) here is the **initial** scale separation: deformation can make finer scales during a stage.

## 4. Why this does not regenerate velocity, Reynolds number or the next stage

**Exact invariant structure.** Write
\[
u=(v_1,v_2,w)(x_1,x_2,t),\quad
\omega=(\partial_2w,-\partial_1w,\partial_1v_2-\partial_2v_1).
\]
Equation (1.5) is exactly
\[
\partial_tv+v\cdot\nabla_hv+\nabla_hp=\nu\Delta_hv,
\quad\nabla_h\cdot v=0,
\qquad
\partial_tw+v\cdot\nabla_hw=\nu\Delta_hw.
\]
The base \(v\) is autonomous 2D NS. The vertical velocity \(w\) is a passive scalar. There is no term through which the amplified \(\nabla_hw\) changes \(v\) or its pressure. Indeed \((w e_3)\cdot\nabla=0\), and
\[
-\Delta_hp=\sum_{i,j=1}^2\partial_iv_j\,\partial_jv_i
\]
contains no \(w\). The horizontal vorticity stretches, but cannot become the next stronger base-strain generator while this invariant structure is preserved. The global smoothness of these 2.5D solutions is explicitly noted in §1.1.

**Velocity bounds, deduced from the exact equations.** The maximum principle gives \(\|w(t)\|_\infty\le\|w(0)\|_\infty\). The planar scalar vorticity also obeys the maximum principle. For mean-zero planar velocity on a torus of scale \(L\), the integrable periodic Biot–Savart kernel gives
\[
\|v(t)\|_\infty\le C L\|\operatorname{curl}_h v(0)\|_\infty.
\]
An initial constant planar mean, if present, is separately conserved. For this paper's \(L=1\), normalized base, the bound is uniform in \(n\) and time. This is a bound, not a maximum principle for \(v\), but it prevents unbounded base-velocity amplification across the normalized family. Whatever amplitude is chosen for \(w(0)\), the passive scalar cannot increase its own supremum.

The fixed-viscosity rescaling therefore begins with the large velocity amplitude \(A_n\); at \(p=4\), \(A_n\asymp\lambda_n^{-5/4}\). It does not dynamically produce that amplitude from a previous bounded-amplitude state. Large gradient/enstrophy gain is not a demonstrated gain in velocity amplitude or in the Reynolds number needed by a smaller successor. Simply putting all these pre-existing high-amplitude stages into a single initial datum would lose smooth bounded initial data, even if their separate energies were summable.

**Domain limitation.** The energy factor \(\lambda_n\) above is for one shrinking periodic domain. Periodically extending that domain to \(\mathbb R^3\) has infinite total energy. Tiling a fixed torus with the shrinking cells removes the helpful volume factor: there are \(\asymp\lambda_n^{-3}\) cells and energy grows like \(A_n^2\). Cutting out a single cell or localizing in the third direction changes the PDE and destroys the exact 2.5D decoupling. No compact, mutually coupled stage on one fixed \(\mathbb R^3\) or torus follows from this calculation.

## 5. Bounded next target and the role of force

The useful asset is a quantitative actual-NS baseline: a finite viscous gradient-amplification episode survives at \(\nu_n\asymp\epsilon^{4+O(\delta)}\), and its scale arithmetic fits a shrinking finite-energy family after explicit fixed-viscosity conversion. The missing target is a **localized transfer lemma for coupled stages**: show that an incoming smooth state generates a stronger velocity/strain at the next smaller scale, produces a specified usable terminal state, and controls pressure, localization, dissipation and inherited history on one domain. Rotating isolated copies does not establish this lemma; overlapping differently oriented copies introduces interactions absent from the source theorem.

Smooth forcing remains allowed for Clay's forced alternatives C/D. The unforced nature of this source is an asset, not a restriction on our proposed research. Force may assist localization or transfer, but the *total* force must satisfy the required all-order smoothness through the proposed accumulation time (and the appropriate spatial conditions). The scale conversion does not establish that forcing condition or license inserting ever larger fine-scale velocities by an uncontrolled force.

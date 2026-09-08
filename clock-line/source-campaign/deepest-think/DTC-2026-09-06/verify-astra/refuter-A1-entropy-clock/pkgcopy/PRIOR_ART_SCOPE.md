# Pass 8: bounded primary-source scope audit

Date: 2026-09-06. Role: `/root/log_clock_sources`. This is a literature scope check, not independent validation of the proposed proof. No current twin, pulse, wire, or private mailbox mathematics was read.

## Target and search outcome

The proposed statement concerns smooth unforced three-dimensional Navier–Stokes flow on \(\mathbb R^3\), with

\[
E=\|u(t_0)\|_2^2,\qquad M=\|\operatorname{curl}u(t_0)\|_\infty,
\qquad \operatorname{Re}_E=E^{2/5}M^{1/5}/\nu,
\]

and a local cap \(\|\omega\|_\infty\le 3M/2\) on an interval of length
\(c/[M(1+\log_+\operatorname{Re}_E)]\). The proposed mechanism uses the actual scalar advection–diffusion kernel, a density bound, a second moment bound, BMO exponential integrability, and relative entropy.

No exact match for that combination of hypotheses and quantitative conclusion was located in this bounded search. That is a search outcome, not a novelty finding. The five works below are the closest retained primary sources. Searches covered combinations of Navier–Stokes, local existence/lifespan, bounded vorticity, logarithmic energy/Reynolds dependence, heat kernel, John–Nirenberg, and relative entropy. Unrelated compressible, two-dimensional, rotating, and higher-Sobolev-norm lifespan results were not treated as matches.

## 1. Actual divergence-free scalar kernel

Elias Hess-Childs, Renaud Raquépas, Keefer Rowan, *Divergence-free drifts decrease concentration*: [arXiv v1 PDF](https://arxiv.org/pdf/2503.16723v1), [journal record](https://doi.org/10.1016/j.jfa.2025.111314), JFA 290(7), 111314 (2026).

ArXiv Definition 1.6, printed p. 2, explicitly allows time-dependent bounded divergence-free drift in \(\partial_t f=\Delta f-u\cdot\nabla f\). Theorem 1.7, p. 3, gives concentration comparison for positive finite measures \(\mu\preceq\eta\), with \(\eta\) symmetric decreasing. The text expressly includes \(\mu=\delta_y,\eta=\delta_0\).

Corollary 1.10, p. 3, gives all \(L^p\) norm upper bounds and a variance **lower** bound against pure heat flow. Differential entropy is also bounded **below** by heat entropy when the initial measure has some positive moment. A delta satisfies that condition. After diffusivity scaling, the delta-data density bound is \((4\pi\nu t)^{-d/2}\).

This directly supports the density input. It does not state the desired second moment upper bound, BMO pairing, or three-dimensional vorticity lifespan. Applying it to a backward row kernel requires explicitly identifying the reversed drift and its sign. The journal preview renumbers the main results as Theorem 1.5 and Corollary 1.8; the inspected arXiv history lists v1 only.

## 2. A stronger stated bounded-vorticity clock

Zachary Bradshaw, Aseel Farhat, Zoran Grujić, *An algebraic reduction of the “scaling gap” in the Navier–Stokes regularity problem*: [arXiv v4](https://arxiv.org/pdf/1704.05546v4).

Theorem 8, printed p. 8, states a mild vorticity solution for initial \(\omega_0\in L^2\cap L^\infty\) on \(T\ge c^{-1}\|\omega_0\|_\infty^{-1}\), in the unit-viscosity convention. The preceding paragraph explicitly says the \(L^2\) assumption has no quantitative effect. Theorem 10, printed p. 11, states the corresponding analyticity and amplification estimate for any prescribed factor greater than one. Remark 11 mentions a logarithmic correction for localization, without an explicit energy–vorticity formula.

This is a stronger *stated* clock than the present target, not an exact match to its energy dependence or method. The earlier local project audit raised a specific concern about the displayed Picard/BMO estimate; the present source task does not validate that proof or reject the theorem. Therefore this source prevents an unqualified “no earlier stronger claim” statement, but is not being used to certify the new target.

## 3. BMO/Besov local existence and continuation

Hideo Kozono, Takayoshi Ogawa, Yasushi Taniuchi, *Navier–Stokes equations in the Besov space near \(L^\infty\) and BMO*, Kyushu Journal of Mathematics 57 (2003), 303–324: [publisher PDF](https://www.jstage.jst.go.jp/article/kyushujm/57/2/57_2_303/_pdf/-char/en), [DOI](https://doi.org/10.2206/kyushujm.57.303).

Theorem 1, p. 306, constructs a unique mild solution from divergence-free velocity in inhomogeneous \(B^0_{\infty,\infty}\), with logarithmically weighted \(L^\infty\) smoothing. Its following remark gives \(T_*=C_\varepsilon\|a\|_{B^0_{\infty,\infty}}^{-2/(1-\varepsilon)}\) for sufficiently small \(\varepsilon>0\). Theorem 2, p. 307, gives continuation when the time integral of the homogeneous \(\dot B^0_{\infty,\infty}\) vorticity norm is finite.

These are pertinent endpoint existence and continuation results. Their displayed hypotheses and clocks are different from a cap depending only on \(E,M,\nu\). They do not themselves provide the proposed entropy/kernel estimate. Lemmas 2.1–2.2 concern frequency splitting and logarithmic heat smoothing, so a further comparison with their proof could be useful before any novelty assessment.

## 4. Stochastic representation of actual three-dimensional NS

Peter Constantin, Gautam Iyer, *A stochastic Lagrangian representation of the three-dimensional incompressible Navier–Stokes equations*, CPAM 61 (2008), 330–345: [author preprint](https://arxiv.org/pdf/math/0511067), [DOI](https://doi.org/10.1002/cpa.20192).

Theorem 2.2 gives the stochastic Weber representation using the flow \(dX=u\,dt+\sqrt{2\nu}\,dW\), its inverse, expectation, and Leray projection, for regular divergence-free data. Proposition 2.7 states the vorticity representation \(\omega=\mathbb E[((\nabla X)\omega_0)\circ X^{-1}]\). Theorem 2.6 records local existence with time controlled by a higher Hölder norm, independently of viscosity.

This is direct precedent for averaging along noisy trajectories of the actual velocity. Its retained deformation gradient is precisely where three-dimensional stretching enters. The inspected statements do not replace that gradient by an \(E,M\) logarithmic clock, nor do they use the proposed BMO/entropy scalar-kernel pairing.

## 5. Relative entropy with advection and diffusion

Peter Constantin, Gautam Iyer, *Stochastic Lagrangian transport and generalized relative entropies*, Communications in Mathematical Sciences 4 (2006), 767–777: [author preprint](https://arxiv.org/pdf/math/0608797), [DOI](https://doi.org/10.4310/CMS.2006.v4.n4.a5).

The paper considers a linear diffusion–transport operator with variable diffusivity and potential. Theorem 1 (credited there to Michel–Mischler–Perthame) states that \(\int H(f/\rho)\phi\rho\) decreases when \(f,\rho>0\) evolve under the same forward operator, \(\phi\ge0\) under its backward adjoint, and \(H\) is convex. Sections 3–5 supply stochastic representations and a Jensen-based proof.

This is a close conceptual antecedent for relative entropy and an actual advecting flow. Its comparison density solves the same forward operator. A freely selected Gaussian reference in a one-time entropy bound is a different use, so the theorem cannot be substituted without checking that distinction. The inspected result does not assert the target Navier–Stokes lifespan.

## Attribution limits for the proposed proof

The complete kernel-to-vorticity-cap argument still requires a proof. In particular, the backward kernel identity, its moment upper bound, the Gaussian exponential BMO bound with an explicitly controlled mean, the resulting entropy pairing, and the continuation/first-exit closure should each be demonstrated in the main note. Source 1 can supply the sharp scalar density bound, but no source here supplies the entire chain.

The elementary inequality \(\int Fp\le H(p\mid g)+\log\int e^F g\) is an entropy variational inequality. Its use in the proposed argument should not be described as a new principle. Likewise, the bounded search gives no basis for claiming new BMO theory or a new stochastic representation. No conclusion about Clay regularity follows from this literature audit.

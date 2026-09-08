# Angular error budget for the full initial-pressure gate

There are two rigorous reductions. For the original cylindrical pump, the angular tail can be bounded through a small number of nonnegative integral norms, using only the poloidal adjoint weight. For the proposed radial annular strain pump, a stronger selection rule makes the entire pressure tail above degree four vanish exactly. Neither reduction certifies the current numerical sign until the remaining integral and datum uncertainties are enclosed.

## 1. Original pump: remove the large azimuthal weight exactly

Keep the exact datum fixed and compare the exact mixed pressure evaluation with the evaluation using the exact angularly truncated whole-space pressure. Write

\[
 u_o=cP+Av,\qquad g_o=c^2g_P+A^2g_v,\qquad
 e=(I-\mathcal P_{\le L})g_o,\qquad \phi=N*e.
\]

Here `P` is poloidal and `v` azimuthal, both axisymmetric. Their bilinear pressure source vanishes. The core source has degrees at most four, so it can be omitted from `e` for `L>=4`.

The tensor `Q=-D²∂zzN` has no azimuthal–poloidal coupling. Since `∇φ` is poloidal,

\[
 (QAv)\cdot\nabla\phi=0
\]

pointwise. The previously proved degree-four selection rule removes the core part of the pressure error. Therefore the entire angular error is exactly

\[
 \delta T_{\rm ang}=2\int_{\mathbb R^3}F\cdot\nabla\phi\,dx,
 \qquad \boxed{F=cQP.}
 \tag{1}
\]

In particular, using `||Q(cP+Av)||` would unnecessarily include the large swirl amplitude in the adjoint weight. The amplitude `A²` still enters the pressure source `g_o`; it has not disappeared from the error.

For even scalar source parity, let `l*` be the first even integer strictly above `L`, and set

\[
 \lambda_* =\ell_*(\ell_*+1),\qquad \kappa_*=\ell_*+\frac12.
\]

## 2. A sharper weighted pressure energy bound

Write `R=|x|`. Radial Hardy and angular Poincare give, for a potential containing only degrees at least `l*`,

\[
 \int|\nabla\phi|^2\ge
 \left(\lambda_*+\frac14\right)\int\frac{|\phi|^2}{R^2}
 =\kappa_*^2\|\phi/R\|_2^2.
\]

The radial Hardy inequality here follows directly from

\[
 \int_0^\infty (R\partial_R\phi+\phi/2)^2dR
 =\int_0^\infty R^2|\partial_R\phi|^2dR
     -\frac14\int_0^\infty|\phi|^2dR,
\]

at each angle, first for suitable smooth functions and then by energy-space approximation. Compact smooth source gives a decaying Newton solution; the radial boundary terms vanish. The whole-space energy identity yields

\[
 \|\nabla\phi\|_2^2=\int e\phi
       \le\|Re\|_2\|\phi/R\|_2.
\]

Consequently

\[
 \boxed{\|\nabla N*e\|_2\le\kappa_*^{-1}\|Re\|_2,
 \qquad\|(N*e)/R\|_2\le\kappa_*^{-2}\|Re\|_2.}
 \tag{2}
\]

This improves the previous `Rmax/sqrt(lambda*)` factor: the source is measured with its actual radius weight, and the denominator includes radial Hardy. All potential tails are retained; no artificial boundary condition is imposed at the source radius.

The simplest rigorous angular budget, requiring **no computed tail coefficients**, is thus

\[
 \boxed{|\delta T_{\rm ang}|\le
       \frac{2}{\kappa_*}\|cQP\|_2\,\|Rg_o\|_2.}
 \tag{3}
\]

At `L=1280`, `l*=1282` and `kappa*=1282.5`. For example, allocating a rigorous error allowance of `60` to this term requires the rigorously enclosed product in (3) to be below `38475`. The number `60` is an optional budget allocation, not a claim that the present norms pass it.

## 3. An additional tail factor from the adjoint weight

The gradient of a degree-`l` scalar pressure is a combination of the degree-`l` radial and poloidal vector spherical harmonics. Therefore only the matching high-degree **vector-harmonic projection** of `F` pairs with `∇φ` in (1). This is not the same as projecting each Cartesian component at the same degree cutoff.

In spherical coordinates, write

\[
 F=a(R,\mu)e_R+b(R,\mu)e_\theta,
 \qquad \mu=z/R,
\]

and define

\[
 d=\operatorname{div}_{S^2}(b e_\theta)
   =-\partial_\mu\!\left(\sqrt{1-\mu^2}\,b\right),
 \qquad \mathcal A=-\partial_\mu[(1-\mu^2)\partial_\mu].
\]

Set the following squared norms, all with physical volume measure:

\[
\begin{aligned}
 F_0^2&=\int(a^2+b^2),\\
 F_1^2&=\int[(1-\mu^2)a_\mu^2+d^2],\\
 F_2^2&=\int[(\mathcal A a)^2+(1-\mu^2)d_\mu^2],
\end{aligned}
 \tag{4}
\]

and source norms

\[
\begin{aligned}
 G_0^2&=\int R^2g_o^2,\\
 G_1^2&=\int R^2(1-\mu^2)(\partial_\mu g_o)^2,\\
 G_2^2&=\int R^2(\mathcal A g_o)^2.
\end{aligned}
 \tag{5}
\]

For an axisymmetric tangent field the surface curl is zero. Expanding its poloidal component in `∇S Y_l/sqrt(l(l+1))` proves that the three quantities in (4) weight its vector-harmonic squared coefficients respectively by `1,lambda_l,lambda_l²`; the radial component has the same weights. The analogous scalar spectral identities give (5).

Combining those identities with (2) gives nine valid bounds:

\[
 \boxed{|\delta T_{\rm ang}|\le
   \min_{0\le j,k\le2}
   \frac{2F_jG_k}{\kappa_*\lambda_*^{(j+k)/2}}.}
 \tag{6}
\]

One can use whichever rigorously enclosed product is smallest. Enclosing these few norms replaces enclosing every high Legendre coefficient. Optional subtraction of rigorously enclosed energies of a few known low modes sharpens the norms further; such subtraction is not required for validity.

For the explicit profiles, `G0 F2` and `G1 F1` both give order `l*^(−3)` decay using cutoff derivatives only through order three. `G2 F0` has the same spectral decay but the pump source requires the fourth cutoff derivative. Thus the adjoint-weight route can be cheaper to certify. The evaluator's existing three cutoff derivatives are sufficient for the first two choices, after generating the needed angular derivatives analytically.

### Scalar adjoint alternative

Let `h=div F`. Compact support away from the origin justifies integration by parts in (1). Applying the weighted pressure estimate to the high-angular part of `h` gives the alternative exact bound

\[
 \boxed{|\delta T_{\rm ang}|\le
       2\kappa_*^{-2}\|Rh\|_2\,\|Rg_o\|_2.}
 \tag{7}
\]

More source or adjoint angular derivatives give the corresponding factors `lambda*^(−k/2)`. This is useful if the full divergence has cancellation; it is not automatically smaller than (6). In spherical components,
`h=R^(−2)∂R(R²a)+d/R`. The gradient adjoint and this scalar adjoint retain the same exact pressure functional.

## 4. What to enclose rigorously for the original pump

The norms in (4)–(5) are nonnegative two-dimensional integrals. In `(R,mu)` their measure is `2pi R² dR dmu`; the source norms therefore use `2pi R⁴ dR dmu`. They can alternatively be integrated over the natural cylindrical support rectangles using `2pi rho d rho dz`.

For cylindrical derivative implementation, the angular derivative at fixed `R` is

\[
 D_\theta=z\partial_\rho-\rho\partial_z,
\]

and

\[
 \Delta_{S^2}g=D_\theta^2g+(z/\rho)D_\theta g.
\]

Apparent axis denominators have removable limits for the smooth axisymmetric profiles. In the spherical formulas, `b` vanishes like `sqrt(1−mu²)`, so `sqrt(1−mu²)b` is the smooth quantity to differentiate. An interval implementation should preserve that cancellation rather than divide two uncertain vanishing quantities.

A rigorous numerical certificate must enclose:

1. `C_P,C_v` and the exact neutral parameter `A²=(204/35−4 C_P)/C_v`, including positivity of `C_v`.
2. One or more pairs `Fj,Gk` in (6), or the two norms in (3), with the exact datum parameters enclosed.
3. The retained pressure-mode projection/radial errors, including the low core channels `0,2,4` and their contact contribution.
4. The explicit advection, viscosity, outer integral and final arithmetic errors.

For source norms, `g_o=4g_P+A²g_v`. Their squared norms are quadratic polynomials in `A²` whose three coefficients are unit-profile norm/cross integrals. Those few integrals can be enclosed once and reused for the entire tuned amplitude interval. The adjoint norms depend on `c=2` but not on `A`.

At fixed `b=Omega=1,c=2,nu=1/100`, the exact gate itself is affine in `A²`, by the parity-reduced cubic formula from pass4. Tuning uncertainty can therefore be handled directly by an interval in `A²`; an enclosure of its square root is not intrinsically required.

For an arbitrary compact low-mode Poisson source residual `r_low`, (2) with minimum degree zero gives

\[
 |\delta\mathcal M_{\rm outer,low}|
       \le4F_0\|Rr_{\rm low}\|_2.
 \tag{8}
\]

This is only the outer part. Core low-mode errors require their separate exact three-mode functional. After removing modes `0,2,4`, the remaining even residual has degree at least six, and the outer bound improves to `2F0 ||Rr||/6.5`. An angular quadrature error in a retained coefficient is a retained-mode error; the high cutoff `l*=1282` must not be assigned to it.

An explicit total budget is

\[
 T_{.01}\le\widehat T_{.01}+E_{\rm ang}
 +E_{\rm retained,outer}+E_{\rm retained,core}
 +E_{\rm explicit}+E_{\rm tuning}+E_{\rm arithmetic}.
 \tag{9}
\]

A certificate requires the right side to be strictly negative. The exploratory split value `−241.348818` is a useful target margin, but observed grid convergence or source/split agreement is not one of the rigorous error terms in (9). Derivative cutoffs, interval quadrature, radial reconstruction and amplitude tuning remain to be enclosed.

## 5. Stronger branch: a radial annular pump has no high-pressure tail

Now replace the cylindrical pump by

\[
 P=S_\eta,
\]

where `eta(R²)` vanishes near zero, equals `−1` on an annulus containing the swirl support, and vanishes outside a larger radius. This is the same radial divergence-free strain-extension operator used for the core. Put

\[
 a_0=\eta(R^2),\qquad q_0=3a_0+2R^2\eta'(R^2).
\]

Its spherical components are exactly

\[
 P_R=Ra_0(3\mu^2-1),\qquad
 P_\theta=-Rq_0\mu\sqrt{1-\mu^2}.
\]

The radial/polar block of the full pressure tensor is

\[
 Q\big|_{R,\theta}=\frac1{4\pi R^5}
 \begin{pmatrix}
 12-36\mu^2&-24\mu\sqrt{1-\mu^2}\\
 -24\mu\sqrt{1-\mu^2}&21\mu^2-9
 \end{pmatrix}.
 \tag{10}
\]

Thus the exact adjoint weight has components

\[
\begin{aligned}
 F_R={c\over4\pi R^4}
    [-12a_0(3\mu^2-1)^2+24q_0\mu^2(1-\mu^2)],\\
 F_\theta={c\over4\pi R^4}\mu\sqrt{1-\mu^2}
    [-24a_0(3\mu^2-1)-q_0(21\mu^2-9)].
\end{aligned}
 \tag{11}
\]

The radial component is an even polynomial of degree four. The tangent component is `sqrt(1−mu²)` times an odd polynomial of degree three, hence a combination of `∇S P2` and `∇S P4`. Therefore

\[
 \boxed{\int F\cdot\nabla\phi_\ell=0\quad\text{for every }\ell>4.}
 \tag{12}
\]

This is stronger than a norm bound: for this branch, the full pressure functional in the mixed evaluation needs **only degrees `0,2,4`**, regardless of the high angular content of the narrow swirl source. The core already has exactly the same selection rule, and the azimuthal part of the outer weight still pairs to zero.

For direct implementation, write

\[
 F={c\over4\pi R^4}
   \sum_{\ell=0,2,4}
     [r_\ell(R)P_\ell e_R+v_\ell(R)\nabla_{S^2}P_\ell],
\]

with

\[
\begin{array}{c|cc}
 \ell&r_\ell&v_\ell\\ \hline
 0&\frac{16}{5}(q_0-3a_0)&0\\
 2&\frac{16}{7}(q_0-6a_0)&\frac{16}{7}a_0\\
 4&-\frac{96}{35}(9a_0+2q_0)&\frac{144a_0+42q_0}{35}.
\end{array}
\]

If the pressure is `sum p_l(R)P_l(mu)`, its exact contribution to the outer mixed derivative is the one-dimensional contraction

\[
 \boxed{2\int F\cdot\nabla p
 =2c\sum_{\ell=0,2,4}{1\over2\ell+1}
   \int R^{-2}\left[r_\ell p_\ell'
        +{\ell(\ell+1)\over R}v_\ell p_\ell\right]dR.}
 \tag{13}
\]

The integrals are nonsingular because `eta` vanishes near zero. This identity retains the complete pressure rather than declaring the high source modes absent. It permits their exact elimination from this particular scalar observable.

The same radial quadratic identity gives `pzz[P](0)=0` and `dP=2Pi(P,Delta P)=0`, since both the annular profile and its Laplacian-direction profile vanish near zero. Those are initial scalar cancellations, not statements that the pump velocity or later viscous effects vanish.

This structural branch removes `E_ang` from (9) altogether. It does **not** remove angular quadrature error in the three retained source coefficients, their whole-space radial pressure reconstruction, the local nonlinear integrals, tuning, or subsequent dynamical obligations. The old favorable numerical sign must also be recomputed for this changed pump; it cannot be inherited from the cylindrical datum.

`verify_angular_tail_budget.py` checks (10) by Cartesian differentiation of the Newton kernel, then checks the exact harmonic coefficients and their orthogonality. The analytic norm estimates above are proved in this note. The checker does not certify any numerical norms or the gate sign.

# Audit of a compact meridional flow increasing the swirl pressure contribution

Status: independent derivation of exact initial-time identities and a thin-support asymptotic. The calculations below support simultaneous rotational growth, neutral initial strain ratio, and increasing **swirl-only** pressure-supply ratio for a suitably chosen smooth compact datum. They do not determine the derivative of the total outer pressure contribution or prove persistence, profile closure or regeneration.

## Definitions and sign conventions

Use cylindrical coordinates \((r,\theta,z)\), with \(W=w(r,z)e_\theta\) smooth, axisymmetric, compactly supported away from the origin, and contained in
\[
\mathcal C=\{r^2/z^2\le1/16\}.
\]
Let \(N(x)=1/(4\pi|x|)\), and define the smooth pressure tensor on this support by
\[
Q_{ij}(x)=-\partial_{zzij}N(x).
\]
For purely azimuthal \(W\), its pressure coefficient is
\[
C_\theta=-\partial_{zz}p[W](0)
=\int Q_{ij}W_iW_j\,dx
=\int Q_\theta(r,z)w^2\,dx,
\]
where
\[
Q_\theta=\frac{3(4z^2-r^2)}{4\pi(r^2+z^2)^{7/2}}>0\quad\text{on }\mathcal C.
\]
All integrals use the three-dimensional volume element, \(dx=r\,dr\,d\theta\,dz\). The horizontal restriction of \(Q_{ij}\) is positive definite on the narrow cone by the exact kernel estimate in `work/pass2/active-core-pressure.md`.

Here \(W\) denotes the **outer** swirl, not the full azimuthal velocity, which also includes the rotating core. To define its time derivative as an observable of the actual solution, fix a smooth axisymmetric spatial cutoff equal to one on a neighborhood of the outer initial swirl support and zero on a neighborhood of the entire compact core support. Choose its transition region in the initial zero-swirl gap and its support away from the origin, and apply it to the actual evolving azimuthal velocity. At the initial instant the cutoff commutators vanish because the initial swirl and its derivatives are zero in that transition region. Thus its derivative is exactly the outer swirl equation used below. At positive time this localization adds commutator terms; no later-time identity without those terms is asserted. Equivalently one may track the outer initial swirl through the linear swirl equation driven by the actual meridional velocity, but that is a labeled component, not the complete swirl field.

Suppose the initial meridional velocity agrees with
\[
P_c=cr e_r-2cz e_z,\qquad c>0,
\]
on an open neighborhood of the entire support of \(W\). A disjoint compact core contributes no initial velocity there. This is an assumption about the initial datum, not the later meridional flow.

## Exact derivative of the swirl pressure functional

The axisymmetric swirl equation for ordinary NS is
\[
w_t+P_r\partial_rw+P_z\partial_zw+\frac{P_r}{r}w
=\nu\left(\partial_{rr}+r^{-1}\partial_r+\partial_{zz}-r^{-2}\right)w.
\]
There is no azimuthal pressure gradient. On the initial swirl support its nonviscous part is
\(w_t=-crw_r+2czw_z-cw\). Since \(\nabla\cdot P_c=0\), integration by parts gives
\[
(C_\theta')_{\rm nonviscous}
=c\int\big(r\partial_rQ_\theta-2z\partial_zQ_\theta-2Q_\theta\big)w^2\,dx.
\]
Set \(t=r^2/z^2\), using \(t\) only as a dimensionless spatial ratio in this formula. Direct differentiation gives
\[
\boxed{\frac{r\partial_rQ_\theta-2z\partial_zQ_\theta-2Q_\theta}{Q_\theta}
=\mathcal A(t):=8-\frac{6t}{4-t}-\frac{21t}{1+t}.}
\]
Moreover
\[
\mathcal A'(t)=-\frac{24}{(4-t)^2}-\frac{21}{(1+t)^2}<0,
\qquad
\mathcal A(1/16)=\frac{2381}{357}>6.
\]
Thus \(\gamma:=2381/357\) is a valid uniform initial lower multiplier on the stated support cone.

For the viscous term it is essential to keep the **Cartesian tensor** and vector field, rather than treating \(Q_\theta\) as a harmonic scalar. Since \(\Delta Q_{ij}=0\) away from the origin, compact integration by parts gives
\[
2\nu\int Q_{ij}W_i\Delta W_j\,dx
=-2\nu\int Q_{ij}\partial_kW_i\partial_kW_j\,dx.
\]
Define
\[
D_\theta=\int Q_{ij}\partial_kW_i\partial_kW_j\,dx\ge0.
\]
The sign follows because \(W\) and each of its Cartesian derivatives are horizontal and supported in \(\mathcal C\), where the horizontal tensor is positive. This identity includes the cylindrical \(-w/r^2\) term and the derivatives of the azimuthal unit vector. No contact term occurs because the support is separated from the evaluation point.

Combining the two contributions proves exactly, at insertion,
\[
\boxed{C_\theta'
=c\int Q_\theta w^2\mathcal A(r^2/z^2)\,dx-2\nu D_\theta
\ge c\gamma C_\theta-2\nu D_\theta.}
\]
For the strained core, \(\Omega'=2b\Omega\). Therefore
\[
\left(\frac{C_\theta}{\Omega^2}\right)'
\ge\frac{C_\theta}{\Omega^2}
\left(c\gamma-4b-\frac{2\nu D_\theta}{C_\theta}\right)>0
\]
whenever
\[
\boxed{c>\frac{4b+2\nu D_\theta/C_\theta}{\gamma}.}
\]
If \(W=Av\) with fixed swirl profile \(v\), both \(C_\theta\) and \(D_\theta\) scale as \(A^2\); their quotient is independent of the amplitude subsequently chosen for tuning.

## Compact meridional realization and its pressure sign

Let \(\alpha\) be smooth, compactly supported in its nonnegative argument and equal to one near zero. Let \(\beta\) be smooth, even, compactly supported near \(z=\pm d\), vanish on a neighborhood of zero, and equal one on the axial swirl supports. For \(\sigma>0\), take the axisymmetric streamfunction
\[
\Psi=-czr^2\alpha(r^2/\sigma^2)\beta(z).
\]
Using \(P_r=-r^{-1}\partial_z\Psi\), \(P_z=r^{-1}\partial_r\Psi\), one obtains
\[
\boxed{P_r=cr\alpha(r^2/\sigma^2)(\beta+z\beta'),\qquad
P_z=-2cz\big(\alpha(r^2/\sigma^2)+(r^2/\sigma^2)\alpha'(r^2/\sigma^2)\big)\beta.}
\]
This is smooth across the axis, compact, divergence free and odd under full inversion. Choose the swirl support inside open plateaus where \(\alpha=1\), \(\alpha'=0\), \(\beta=1\), \(\beta'=0\), and choose the axial support disjoint from the core. It then realizes the required affine meridional velocity exactly on the swirl support.

Let \(P_1=P_c/c\) and \(C_P=-\partial_{zz}p[P_1](0)\). Keep \(\alpha,\beta\) fixed and let \(\sigma\to0\). With \(s=r/\sigma\) and
\(F(s)=\alpha(s^2)+s^2\alpha'(s^2)\),
\[
(P_1)_r=\sigma s\alpha(s^2)(\beta+z\beta'),\qquad
(P_1)_z=-2zF(s)\beta.
\]
On the axis away from zero,
\[
Q_{zz}(0,z)=-\partial_{zzzz}N(0,z)=-\frac6{\pi|z|^5}.
\]
Hence the vertical-velocity square contributes, using \(2\pi r\,dr\,dz=2\pi\sigma^2s\,ds\,dz\),
\[
-48\sigma^2
\left[\int_0^\infty sF(s)^2\,ds\right]
\left[\int_{\mathbb R}|z|^{-3}\beta(z)^2\,dz\right].
\]
The remaining terms are \(O(\sigma^4)\): the horizontal velocity is \(O(\sigma)\), the radial-vertical kernel component is \(O(\sigma)\), and the vertical kernel differs from its axial value by \(O(\sigma^2)\). The axial support is bounded away from zero, so these bounds are uniform on the fixed rescaled support. Consequently
\[
\boxed{\frac{C_P}{\sigma^2}\longrightarrow
-48\left[\int_0^\infty s(\alpha(s^2)+s^2\alpha'(s^2))^2\,ds\right]
\left[\int_{\mathbb R}|z|^{-3}\beta(z)^2\,dz\right]<0.}
\]
Both integrals are strictly positive for the specified nonzero cutoffs. Thus \(C_P<0\) for sufficiently small \(\sigma\).

## Pressure additivity and simultaneous tuning

Although \(P_c\) and \(W\) overlap, their pressure cross source vanishes for axisymmetric fields. Indeed
\[
(P_c\cdot\nabla)W+(W\cdot\nabla)P_c
=\left(P_r w_r+P_z w_z+\frac{P_r}{r}w\right)e_\theta,
\]
which has zero divergence. The pressure Poisson equation and the decaying normalization therefore give \(p[P_c+W]=p[P_c]+p[W]\). The compact core is disjoint, so its quadratic pressure cross sources also vanish initially. It follows that the **total initial outer** coefficient is exactly
\[
C_{\rm out}=A^2C_v+c^2C_P.
\]
The neutral initial strain-ratio condition from `strained-receiver.md` is
\[
C_{\rm out}=C_*:=\frac{38}{7}b^2+\frac25\Omega^2.
\]
Once geometry is fixed with \(C_P<0\), first choose \(c\) above the swirl-ratio threshold, then choose
\[
\boxed{A^2=\frac{C_*-c^2C_P}{C_v}>0.}
\]
This choice does not disturb that threshold because \(D_\theta/C_\theta\) is amplitude independent. It establishes simultaneously at the initial instant: \(\Omega'>0\), \((b/\Omega)'=0\), and \((C_\theta/\Omega^2)'>0\).

## Limits requiring explicit retention

- The growing ratio is the swirl contribution alone. The meridional pressure contribution also evolves, so its derivative can offset that gain. No sign for \(C_{\rm out}'\) or the total normalized pressure supply follows from these identities.
- Equality of the initial strain-ratio derivative to zero does not give forward invariance or a positive-time neutral ratio. Its next derivative is not controlled here.
- The \(\sigma\to0\) sign argument does not keep an arbitrary fixed swirl support inside the meridional plateau. A compatible swirl family must fit inside the shrinking radial plateau. Its diffusion quotient can grow as \(\sigma\) shrinks, so the required \(c\), compensating swirl amplitude, energy and Sobolev norms need not remain bounded.
- The formula for the favorable multiplier applies where the initial meridional field equals \((cr,-2cz)\). Later support deformation, cutoff effects and evolution of the meridional field require new estimates. Even the affine reference flow increases \(r^2/z^2\) at rate \(6c\), so persistence in the narrow cone cannot be assumed.
- The tensor signs and asymptotic coefficient were independently checked by differentiation and exact algebra. They do not replace the analysis needed for a complete ordinary-NS stage. This construction uses no force, and it imposes no prohibition on smooth forcing in the wider C/D program.

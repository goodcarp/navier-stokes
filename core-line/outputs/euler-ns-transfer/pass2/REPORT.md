# Second pass: an active core and the missing regeneration mechanism

8 September 2026. **Navier–Stokes blowup remains unproved.** This pass adds an explicit local lemma for an actual finite-energy NS solution, connects our scale proposal to an established viscous family, and makes the conditions for a successful next stage more precise.

## 1. A concrete local lemma for the actual equation

Take a compact rotating core, with velocity equal to `Omega(-y,x,0)` near the origin, and disjoint compact horizontal vortex packets above and below it. An explicit nonaxisymmetric choice is given in [active-core-pressure.md](active-core-pressure.md). The entire datum is smooth, compactly supported and divergence free. The force is zero and viscosity is any fixed positive number.

The exact pressure calculation gives

\[
K=A^2 C_v-\frac25\Omega^2,
\qquad
C_v\ge\frac6{4\pi}\int\frac{|v(x)|^2}{|x|^5}\,dx>0.
\]

For an outer-packet amplitude A with K>0, the local smooth NS solution satisfies

\[
\partial_z u_z(t,0)=Kt+o(t),\qquad
\omega_z(t,0)=2\Omega+\Omega Kt^2+o(t^2).
\]

Thus the initially active central vortex begins to acquire stronger axial strain and vorticity. This is a calculation of the actual PDE's initial derivatives, supported by local smooth existence; it is not merely a prescribed strain ODE. The core's own pressure contributes axial compression, and the positive threshold explicitly pays that contribution.

Two independent symbolic derivations agree, including the distributional contact term in the core pressure and the sign of the outer-packet kernel. Their scope is recorded in [local-pressure-independent-audit.md](local-pressure-independent-audit.md). No originality claim is made for this elementary local mechanism.

**What is still missing:** a quantitative large or retained finite gain, a favorable terminal geometry, and a restart that increases the relevant Reynolds number. A large initial acceleration alone does not establish any of these; its useful interval can shrink with the data parameters.

## 2. The proposed scale window has a published viscous counterpart

Jeong–Yoneda's [actual NS amplification theorem](https://arxiv.org/abs/2001.02333) supplies gradient amplification in a family of globally smooth 2.5-dimensional flows. A fixed-viscosity conversion of its quantitative parameters gives, in one shrinking periodic domain,

\[
\mathrm{Re}\asymp\lambda^{-1/4},\quad
M\asymp\lambda^{-1/16+O(\delta)},\quad
E\asymp\lambda^{1/2},\quad
\tau\asymp\lambda^{9/4}.
\]

This matches the leading exponents of our earlier feasibility calculation. The derivation and source equation are in [jeong-yoneda-stage.md](jeong-yoneda-stage.md).

The limitation is exact: their vertical velocity is a passive scalar, and its amplified gradients cannot feed back into the planar base flow. Velocity remains bounded in the normalized source setting. The large velocity in the fixed-viscosity conversion is already present in its initial data, while its favorable energy accounting uses a shrinking domain. This is a useful quantitative benchmark, not a localization or regeneration theorem on one fixed domain.

## 3. Stretching alone spends the resources needed for restart

The exact affine examples in [affine-stretch-gate.md](affine-stretch-gate.md) separate several effects:

- Isotropic transverse stretching can greatly amplify vorticity while the carrier Reynolds number decreases through viscosity.
- A stretched circular vortex retains the total cross-sectional circulation that fixes its transverse Reynolds number.
- Anisotropic waves can increase wavelength Reynolds number, but the exact product of that number and the relevant localization margin decreases. Their gain consumes the initial margin needed to keep a localized wave approximation valid.

These are statements about specified exact models and a stated localization ansatz. They do not rule out a genuinely three-dimensional mechanism that gathers vorticity, changes geometry or restores localization. They identify what such a mechanism has to accomplish.

## 4. An all-order correction budget can fit, conditionally

For the provisional scales `Re~rho^(-1/4)`, `M~rho^(-1/16)`, suppose a new viscous stage actually supplies uniform estimates for its full normalized residual. Under the precise hypothesis in [correction-window.md](correction-window.md), depth

\[
J=38K+87
\]

gives physical mixed derivatives through order K bounded by `A_K rho^2`. A slow order schedule can pay even constants of size `exp[O(K^2 log K)]` while preserving the principal diffusion margin. The scalar correction majorant and endpoint quantifiers have been independently checked.

This is conditional compatibility. The crucial uniform viscous estimates have not been proved. A correction loss `Re^p/M^s` needs `p/s<1/4` for our proposed carrier scale; a polynomial gap with any carrier choice needs `p/s<1/2`. Boundary cases depend on constants.

The leading packet mean is a separate issue. [mean-feedback.md](mean-feedback.md) derives an exact localized-wave stress that stays of leading size as carrier frequency increases. Its pressure-projected divergence can drive useful strain. If an uncancelled mean is instead placed in the force with only fixed algebraic smallness, the stated shrinking-core model fails sufficiently high derivative bounds. Increasing oscillatory correction depth alone does not remove that mean.

## Next target

Retain the intended pressure/Reynolds-stress feedback in the evolving velocity. Prove a finite transfer on localized divergence-free data that increases the next core's usable velocity/circulation scale, restores its localization margin, and leaves the required terminal geometry. The total residual must include the mean, activation, joins and inherited fields, with all-order bounds. The active-core local lemma supplies a favorable initial sign; every sustained-transfer obligation remains open.

The full objective remains a verified NS result. These local and conditional lemmas do not complete it. See [VERIFICATION.md](VERIFICATION.md) for the exact scope of the checks and [../NEXT_TARGET.md](../NEXT_TARGET.md) for the continuing finite-stage requirements.

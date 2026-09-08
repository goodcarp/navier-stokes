# Independent audit: compact affine-core pressure and the next feedback jet

2026-09-08. Reviewed `affine-core-pressure.md`, ran its original checker, and independently recomputed the strain contribution using cylindrical components. Also checked the root's proposed full second derivative of `beta=b/Omega`. No original file was changed.

**Verdict:** the cutoff-independent coefficient `-18/7` and the proposed full-feedback identities are correct. No blocking correction is needed. The exterior functional's derivative cannot determine the next strain-ratio derivative without the evolving core-profile defect.

## 1. Independent strain coefficient and contact term

The stated vector-potential extension simplifies, for its unit poloidal part, to

`u_r=-r[psi+2z^2 psi']`,

`u_z=2z[psi+r^2 psi']`.

Here `psi` is evaluated at `s=r^2+z^2`, and primes are derivatives in `s`. For axisymmetric poloidal velocity the pressure source can be computed independently as

`g_S=(partial_r u_r)^2+(u_r/r)^2+(partial_z u_z)^2+2(partial_z u_r)(partial_r u_z)`.

Inserting these components and averaging with `3mu^2-1`, `mu=z/sqrt(s)`, reproduces exactly

`< (3mu^2-1) g_S >`

`= (16s/105)[-4s^2 psi'psi''+6s psi psi''+2s(psi')^2+21psi psi']`.

The radial measure is `dr/r=ds/(2s)`, which explains the factor one-half in the note. The endpoint identities

`integral psi psi'=-1/2`,

`integral s psi psi''=1/2-I`,

`integral s^2 psi'psi''=-I`, `I=integral s(psi')^2`,

are correct: boundary terms vanish because the cutoff is constant near zero and compactly supported. Their substitution gives `-4/7`, with the coefficient of `I` exactly zero. The local source is `g_S(0)=6`, so the distributional Hessian contact term is `-6/3=-2`. Therefore

`p_zz[strain](0)=-2-4/7=-18/7`.

For swirl, the extension uses the radial profile `f=psi+2s psi'/3`; it has the same endpoint values `f(0)=1`, `f(infinity)=0`. The earlier universal radial-swirl calculation already implies `p_zz=2Omega^2/5`, agreeing with the displayed independent cancellation. The source has no strain–swirl cross term. Taking the pressure trace then gives the transverse coefficient `-12b^2/7+4Omega^2/5`.

The original Cartesian checker and the independent cylindrical checker both pass. The conclusion is specific to this vector-potential extension; it does not imply universal cutoff independence for arbitrary compact affine extensions.

## 2. Full initial dynamics and tuning

At the axisymmetric affine insertion,

`L=diag(-b,-b,2b)+Omega J`,

and the exact gradient equation is `L_t=-L^2-D^2p_0`, because `u_0(0)=0` and `Delta u_0=0` locally. Consequently

`b'=-2b^2-p_zz/2`, `Omega'=2bOmega`.

The computed compact-core pressure and disjoint outer contribution reproduce

`b'=A^2 C_v/2-5b^2/7-Omega^2/5`.

At the tuned insertion `b'=2b^2`, the full pressure satisfies `p_zz=-8b^2` and `beta'=0` for `beta=b/Omega`.

The first two viscosity-curvature observations are justified even though the affine matrix now has a nonzero symmetric part. Locally,

`u_t(0)=-L^2 x-grad p_0`,

`Delta p_0=-tr(L^2)` is constant,

and hence `Delta u_0=Delta u_t(0)=0` throughout the initial affine neighborhood. It follows that both curvature terms in the later `b,Omega` equations and their **first** time derivatives vanish at insertion. Higher curvature derivatives are not covered by this reasoning.

Differentiating the actual ratio equation at the tuned point therefore gives exactly

`beta''(0)=-(p_zz'(0)+32b^3)/(2Omega)`.

No substituted time-dependent compact-core pressure formula is used in this identity.

## 3. The full pressure derivative and its contact term

For the actual full datum, including any outer meridional pump,

`u_t(0)=nu Delta u_0-P div(u_0 tensor u_0)`,

`p_t(0)=N*[2 tr(grad u_0 grad u_t(0))]`, `N=1/(4pi|x|)`.

The Leray projection here retains the full pressure response. Although `u_t` need not be compactly supported, the pressure-source derivative is compactly supported because `grad u_0` is compactly supported. Its Hessian convolution is therefore well-defined with the usual contact/PV interpretation.

At the tuned central matrix,

`2tr(L L_t)=12b b'-4Omega Omega'=24b^3-8bOmega^2`.

Accordingly the contact contribution to `p_zz'(0)` is exactly

`-(24b^3-8bOmega^2)/3`.

This contact contribution is not the whole pressure derivative. The principal-value integral of the full source derivative remains to be evaluated or bounded.

## 4. Why the exterior derivative is insufficient

For an evolving exterior stress functional `C_F(t)`, define a pressure-profile defect by

`D(t)=p_zz(t,0)-[-18b(t)^2/7+2Omega(t)^2/5-C_F(t)]`.

At the specified insertion `D(0)=0`. There is no reason from initial support separation or the central matrix form alone that `D'(0)=0`.

Let

`F(t)=C_F(t)/Omega(t)^2-(38/7)[b(t)/Omega(t)]^2-2/5`.

At the tuned boundary `F(0)=beta'(0)=0`, direct differentiation gives

`beta''(0)=(Omega/2)F'(0)-D'(0)/(2Omega)`.

Thus `F'(0)<0` proves loss of that exterior cone observable, but does not by itself prove `beta''(0)<0`. Likewise a favorable swirl-only `C_theta'` does not fix `C_F'`, `D'`, or `beta''`. The pure-azimuthal reservoir corollary in the source note is correctly restricted; a meridional pump requires its additional evolving contribution to be retained.

## Verification

`python3 work/pass4/derive_affine_core_pressure.py` passed. `python3 work/pass4/verify_affine_pressure_independent.py` independently passed the cylindrical angular coefficient, cutoff cancellation, tuned pressure-source/contact derivative, exact `beta''` formula and defect identity. These verify initial algebraic identities, not the sign of the full cubic pressure integral or forward invariance of a profile class.

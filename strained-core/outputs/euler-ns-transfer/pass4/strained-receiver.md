# First-order receiving gain with an inherited affine strain

Status: exact initial-time identities for one ordinary, unforced NS solution, at every fixed viscosity \(\nu>0\). This expands the local receiving calculation to an initially strained core. It does not assert affine evolution at positive time, an invariant profile class, or regeneration.

## Initial field and exact pressure split

Let
\[
Jx=(-x_2,x_1,0),\qquad D=\operatorname{diag}(-1,-1,2),\qquad
M=bD+\Omega J,\qquad b,\Omega>0.
\]
Assume the smooth divergence-free compact datum equals \(Mx\) on a ball about the origin. It may include the disjoint paired outer packets from the earlier pressure construction. Let \(p_0\) denote its full decaying initial pressure, and write \(H=\nabla^2p_0(0)\). No pressure component is omitted.

Within the affine core,
\[
\Delta u_0=0,\qquad
u_t(0,x)=-M^2x-\nabla p_0(x),\qquad
\operatorname{tr}H=-\operatorname{tr}M^2=2\Omega^2-6b^2.
\]
The symmetric part of \(M^2\) is
\[
C=\operatorname{diag}(b^2-\Omega^2,b^2-\Omega^2,4b^2),
\qquad M^2=C-2b\Omega J.
\]
Consequently the **exact local initial acceleration** splits as
\[
\boxed{u_t(0,x)=2b\Omega Jx+\nabla\psi(x),\qquad
\psi=-\tfrac12 x^TCx-p_0,\qquad\Delta\psi=0.}
\]
The rotational term is pressure independent. The harmonic gradient contains the full pressure response. Its Hessian at the center is \(-C-H\), generally not axisymmetric. Viscosity contributes zero to this first derivative because the initial field is affine on an open neighborhood; viscosity affects subsequent derivatives and the time interval.

## Exact finite-radius receiving coefficients

For a nonnegative, nonzero smooth radial weight \(\chi_r(x)=\chi(|x|^2/r^2)\) supported strictly inside the affine core, put
\[
I_r=\int\chi_r x_1^2\,dx=\int\chi_r x_2^2\,dx=\int\chi_r x_3^2\,dx.
\]
Define two actual weighted velocity coefficients,
\[
\Omega_r(t)=\frac{\int\chi_r u(t)\cdot Jx\,dx}{2I_r},\qquad
b_r(t)=\frac{\int\chi_r u(t)\cdot Dx\,dx}{6I_r}.
\]
Initially \(\Omega_r(0)=\Omega\) and \(b_r(0)=b\), at every receiver radius. The test \(\chi_rJx\) is divergence free; \(\chi_rDx\) generally is not. Hence only the rotational projection automatically annihilates every pressure gradient.

The radial harmonic first-moment identity proved in `work/pass3/rotational-receiver.md` is
\[
\int\chi_r x_i h(x)\,dx=I_r\partial_i h(0)
\]
for every harmonic \(h\) on the receiver ball. Apply it to \(h=\partial_j\psi\). Orthogonality of symmetric and antisymmetric matrices then gives
\[
\boxed{\Omega_r'(0)=2b\Omega,\qquad
b_r'(0)=\tfrac12\psi_{33}(0)=-2b^2-\tfrac12H_{33}=:b'_0.}
\]
In particular the normalized rotational mode has
\[
\boxed{a_r(t):=\Omega_r(t)/\Omega,\qquad a_r(0)=1,\qquad a_r'(0)=2b>0.}
\]
These are identities at finite receiver radii, not just pointwise center formulas. Nonquadratic harmonic pressure terms integrate out of their initial derivatives. All pressure-dependent anisotropy and higher spatial terms remain present in the actual velocity and may affect later times.

Equivalently, the central matrix derivative is \(\nabla u_t(0,0)=-M^2-H\). Its rotational coefficient is \(2b\Omega\), and its axisymmetric strain coefficient, defined by projection onto \(D\), is \(b'_0\). For nonaxisymmetric outer packets the remaining symmetric traceless components need not vanish. Circular paired azimuthal packets preserve the relevant central axisymmetry initially; this does not remove nonaffine higher spatial terms.

## RMS velocity: the exact pressure dependence

Set
\[
E_r(t)=\frac12\int\chi_r|u(t)|^2\,dx,\qquad
V_r(t)=\left(\frac{2E_r(t)}{\int\chi_r}\right)^{1/2}.
\]
The two weighted modes are orthogonal, with squared norms \(2I_r\) and \(6I_r\). Initially
\[
E_r(0)=I_r(\Omega^2+3b^2).
\]
Since the initial field is exactly the sum of these modes, its orthogonal residual is initially zero and contributes no first energy derivative. Therefore
\[
\boxed{E_r'(0)=I_r(4b\Omega^2+6b b'_0),\qquad
\frac{V_r'(0)}{V_r(0)}=
\frac{2b\Omega^2+3b b'_0}{\Omega^2+3b^2}.}
\]
Equivalently,
\[
E_r'(0)=I_r(4b\Omega^2-12b^3-3bH_{33}).
\]
This separates the pressure-free rotational contribution from the pressure-dependent strain contribution. It is exact at every admissible finite radius.

Unlike the purely rotating datum, an increasing rotational projection alone does **not** establish total RMS growth: the initial datum also carries strain energy. For \(b>0\), the necessary and sufficient condition for a positive initial RMS derivative is
\[
b'_0>-\frac23\Omega^2.
\]
The stronger shape-ratio condition below implies \(V_r'(0)/V_r(0)\ge2b\). The general weighted decomposition still supplies
\[
E_r(t)\ge I_r\big(\Omega_r(t)^2+3b_r(t)^2\big),
\]
but it does not supply lower bounds for the two coefficients at later times without a remainder estimate.

## Compact extension and the explicit pressure margin

The independent computation in `work/pass4/affine-core-pressure.md` uses the divergence-free compact extension
\[
u_M=\operatorname{curl}\!\left[-\frac{\vartheta(|x|^2)}3\,x\times(Mx)\right]
=\left(\vartheta+\frac{2|x|^2}{3}\vartheta'\right)Mx
-\frac{2\vartheta'}3(x\cdot Mx)x,
\]
where \(\vartheta=1\) near zero and vanishes outside the core cutoff. With disjoint outer velocity \(Av\), write \(C_v=-\partial_{33}p[v](0)>0\). The pressure calculation gives
\[
H_{33}=-\frac{18}{7}b^2+\frac25\Omega^2-A^2C_v.
\]
Using that coefficient in the receiving formulas yields
\[
\boxed{b'_0=\frac12A^2C_v-\frac57b^2-\frac15\Omega^2,
\qquad \Omega'_0=2b\Omega.}
\]
The cutoff pressure coefficient is imported from the separate exact pressure calculation, rather than inferred from a homogeneous affine-flow model. The formula applies to the chosen compact extension; an arbitrary different compact extension need not have this pressure coefficient.

## Isotropic scale changes cannot repair a degrading strain ratio

For an initially affine core of radius \(L\), the rotational and strain Reynolds observables are
\[
\mathcal R_\Omega=\frac{\Omega L^2}\nu,\qquad
\mathcal R_b=\frac{bL^2}\nu,
\qquad\eta=\frac b\Omega.
\]
An isotropic change from \(L\) to \(qL\), comparing the same receiving coefficients, multiplies both Reynolds observables by \(q^2\) and leaves \(\eta\) unchanged. Any common NS amplitude renormalization also cancels out of this ratio. At the initial time,
\[
\boxed{\eta'_0=\frac{b'_0-2b^2}\Omega.}
\]
Thus prevention of ratio degradation **at first order** requires
\[
\boxed{b'_0\ge2b^2,\quad\text{equivalently }H_{33}\le-8b^2.}
\]
For the displayed compact extension this becomes
\[
\boxed{\frac12A^2C_v\ge\frac{19}{7}b^2+\frac15\Omega^2.}
\]
Writing \(\mathcal P=A^2C_v/\Omega^2\) for the dimensionless outer pressure supply, one obtains
\[
\eta'_0=\Omega\left(\frac{\mathcal P}{2}-\frac{19}{7}\eta^2-\frac15\right).
\]
The notation \(\mathcal P\) here is the scalar pressure-supply ratio called \(R\) in the parent's stage budget, not the Leray projection. Equality is neutral only at first order. A strict positive margin gives a short initial interval of ratio increase for the fixed smooth datum; persistence across a full stage requires estimates for the evolving pressure supply and all other state variables. No invariant cone is proved here.

For clarity, if the receiving radius is selected as \(q(t)=1-\kappa t+O(t^2)\), then
\[
\frac{\mathcal R_{\Omega,\mathrm{new}}}{\mathcal R_{\Omega,\mathrm{old}}}
=1+(2b-2\kappa)t+O(t^2),\qquad
\frac{\mathcal R_{b,\mathrm{new}}}{\mathcal R_{b,\mathrm{old}}}
=1+\left(\frac{b'_0}{b}-2\kappa\right)t+O(t^2).
\]
Rotational Reynolds gain needs \(\kappa<b\) at this order. The difference between the two logarithmic gain rates is \(b'_0/b-2b\), independent of the radius selection. Shrinking the receiver therefore cannot cure a negative initial shape-ratio derivative.

The matrix contractions, generic trace-compatible pressure Hessian formulas, RMS derivative and compact-core substitution were checked by exact SymPy algebra. This check certifies finite identities only. All observations concern fixed Eulerian receivers, which may exchange energy and material with their surroundings. They provide no transported-loop circulation gain, prescribed terminal profile or outer-packet regeneration. Smooth forcing remains available in the wider C/D program; these initial identities assume the unforced equation.

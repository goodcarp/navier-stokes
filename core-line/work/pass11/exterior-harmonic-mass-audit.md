# Exterior harmonic energy as an open radial projection metric

The proposed boundary mass and Robin gradient are correct for an exterior velocity constrained to be a harmonic gradient. They give a compatible way to remove the impermeable-wall pressure condition from the MAC projection. This is an exact continuum exterior-energy identity combined with a specified lumped interior boundary mass; it is not a proof that the evolved NS exterior remains irrotational. Viscosity, trace matching, axial periodicity and whole-space reconstruction still require separate treatment.

The independent checker is `check_exterior_harmonic_mass.py`, with record `exterior-harmonic-mass-check.json`. It checks the boundary algebra symbolically, integrates the exterior mass and full gradient energy numerically, and constructs a small weighted-adjoint MAC projector independently of the evolving solver.

## Exact exterior mode

For a mode \(e^{im\theta+ikz}\), let \(n=|m|\), and write
\(u_{\rm ext}=\nabla[f(r)e^{im\theta+ikz}]\). The scalar equation is

\[
f''+r^{-1}f'-(m^2/r^2+k^2)f=0.
\]

For \(k\ne0\), the decaying solution is proportional to \(K_n(|k|r)\); the defining equation and decay are recorded in [DLMF 10.25](https://dlmf.nist.gov/10.25). Define

\[
\kappa=-\frac{f'(R)}{f(R)}
 =-\frac{|k|K_n'(|k|R)}{K_n(|k|R)}>0,
\qquad q=u_r(R)=f'(R).
\]

For \(k=0,n>0\), use \(f\propto r^{-n}\), so \(\kappa=n/R\). The exterior traces are exactly

\[
f(R)=-q/\kappa,\quad
u_r(R)=q,\quad u_\theta(R)=-\frac{im}{R\kappa}q,
\quad u_z(R)=-\frac{ik}{\kappa}q.                       \tag{1}
\]

For the \((m,k)=(0,0)\) mode, the nonconstant scalar harmonic solution is logarithmic, with \(u_r\propto1/r\) and infinite radial kinetic energy. Its finite-energy radial degree of freedom must therefore be removed: \(q=0\). A constant scalar gauge gives zero velocity. An azimuthal circulation \(u_\theta\propto1/r\), which would require a multivalued potential, also has infinite radial energy per axial period. Do not manufacture a finite zero-mode mass by dividing by \(\kappa=0\).

## Exterior kinetic mass

Multiply the scalar ODE by \(r\overline f\), integrate, and use decay at infinity:

\[
\int_R^\infty r\left(|f'|^2+(m^2/r^2+k^2)|f|^2\right)dr
=-R\overline{f(R)}f'(R)=\frac R\kappa|q|^2.             \tag{2}
\]

This is the **entire three-component** exterior radial L2 norm. It is not just the radial velocity energy, and no separate tangent/axial exterior mass should be added to it. The common angular/axial Parseval factors and the conventional energy factor \(1/2\) are applied afterward. For positive angular modes stored once, their conjugate-mode factor remains two.

Let the last pressure center be \(r_c=R-\delta\). If the MAC radial boundary degree uses the lumped interior half-cell mass \(R\delta\), its full mass becomes

\[
W_R=R\left(\delta+\kappa^{-1}\right).                  \tag{3}
\]

The exterior part of (3) is exact. The interior part is the specified face-centered lumping, consistent with ordinary \(r_{\rm face}\,\Delta r_{\rm dual}\) weights. It is not the exact volume integral \(\int_{R-\delta}^Rr\,dr=R\delta-\delta^2/2\). Replacing it by that exact dual-volume mass would change the discrete Robin denominator; these are two different discretization choices.

## Weighted-adjoint divergence and Robin gradient

The last pressure cell's radial divergence receives \(Rq/V_{\rm last}\). With pressure weights \(V_i\), require the full kinetic adjoint
\(G=-M^{-1}D^*V\). Its boundary row is then

\[
(Gp)_R=-\frac R{W_R}p_{\rm last}
=-\frac{\kappa}{1+\kappa\delta}p_{\rm last}.             \tag{4}
\]

This is precisely the boundary value obtained from the exterior DtN relation \(p_r(R)=-\kappa p(R)\) and linear center-to-boundary reconstruction
\(p(R)=p_{\rm last}+\delta p_r(R)\). Equation (4) is exact for the stated discrete pair; the linear reconstruction has its normal spatial truncation error. The resulting

\[
\Pi=I-G(DG)^{-1}D
\]

is kinetic-orthogonal and idempotent, and removes its own gradients. The zero mode uses the compatible flux constraint and pressure gauge instead of an invertible open-boundary scalar operator.

The outer mass depends on both angular and axial frequency through \(\kappa\). Thus the full kinetic norm is diagonal **after axial Fourier transformation**, not a single unchanged physical-space radial diagonal mass. The adjoint nonlinear scatter must use that same metric. Retaining the previous physical-space mass division for the outer degree would destroy the claimed kinetic adjointness. The dependence is even in \(m,k\), preserving the expected reality symmetries.

## Full exterior viscous stiffness

Harmonic velocity satisfies \(\Delta u_{\rm ext}=0\) pointwise, but its gradient energy is nonzero. Since its Cartesian components are harmonic, Green's identity gives

\[
\int_R^\infty r|\nabla u_{\rm ext}|^2dr
=-R\operatorname{Re}\sum_{a=r,\theta,z}
\overline{u_a(R)}\,\partial_ru_a(R).                    \tag{5}
\]

At \(R\), direct differentiation of (1) and the scalar ODE give

\[
\partial_ru_r=-\left[\frac1R+\frac{m^2/R^2+k^2}{\kappa}\right]q,
\]
\[
\partial_ru_\theta=im\left(\frac1R+\frac1{R^2\kappa}\right)q,
\qquad \partial_ru_z=ikq.
\]

Substitution in (5) yields the positive scalar stiffness

\[
D_{\rm ext}=1+\frac{2R(m^2/R^2+k^2)}\kappa
               +\frac{m^2}{R^2\kappa^2},\qquad
\int_R^\infty r|\nabla u_{\rm ext}|^2dr=D_{\rm ext}|q|^2. \tag{6}
\]

For \(k=0,n>0\), this simplifies to \(D_{\rm ext}=2(n+1)\). If the implementation separately adds \(k^2\) times the **full** kinetic mass, its additional radial/horizontal exterior stiffness must be

\[
D_{\rm ext}-k^2R/\kappa
=1+\frac{2m^2/R+Rk^2}\kappa+\frac{m^2}{R^2\kappa^2};     \tag{7}
\]

adding (6) in that convention would count the axial contribution twice.

Equations (6)–(7) do not justify retaining the old no-slip interior boundary penalties. The interior vector trace must match (1); otherwise a tangential jump creates a distributional vortex sheet and the composite velocity is not H1. A compatible weak stiffness must impose that trace through its interpolation/test space, include the interior derivatives with the new boundary values, and add the exterior block once. Stronger Sobolev residual validation needs correspondingly smooth joins or a separate reconstruction with a quantified defect. One radial mass edit is sufficient for the projection identity but insufficient for a full viscous solver.

## Pressure and dynamics scope

At the initial compact datum, and at radii outside the support of the initial nonlinear forcing, this exterior condition reproduces the correct scalar open-exterior pressure relation. It removes the artificial immediate impermeable-wall condition, while axial periodic images remain a separate error.

For an irrotational exterior, \(u\times\omega=0\). Therefore the scalar removed by a Lamb-form projection is harmonic there, even when the exterior velocity is nonzero. This scalar is the Bernoulli pressure \(p+|u|^2/2\), up to the projection's sign convention. The physical pressure \(p\) itself generally is not harmonic there: its \(-|u|^2/2\) part supplies the nonzero exterior quadratic pressure source. An evolved pressure-curvature observable must include this conversion.

The actual viscous flow normally develops exterior vorticity and does not stay in this harmonic-gradient class. Constraining it to that class defines an approximation space whose omitted diffusion/vorticity tail must be charged. Nor does radial finite energy per axial period provide finite energy after repeating the field along the whole axial line. These are explicit reconstruction and residual obligations, not reasons to discard the compatible boundary improvement.

## Independent checks

The checker tests \((m,k)=(0,0.7),(0,2.3),(1,0),(4,0),(4,0.7),(8,0),(8,2.3)\). Direct integration uses all Cartesian Hessian terms written in cylindrical components, independently of the boundary formula. A small nonuniform seven-cell MAC matrix is assembled from fluxes and (3); its gradient is computed from the weighted adjoint rather than inserting (4). It then checks (4), divergence removal, idempotence, kinetic self-adjointness and complete removal of a random pressure gradient. These tests check the stated construction and its constants; they are not interval certification of special functions or of an evolved field.

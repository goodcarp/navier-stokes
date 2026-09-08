# PASS3: the active core's first vorticity gain is canceled by material area contraction

2026-09-08. This analyzes the actual compact, unforced NS datum in `work/pass2/active-core-pressure.md`, with its fixed finite parameters. It does not change the initial data, assert exact viscous Kelvin conservation, or give a general obstruction to three-dimensional regeneration.

## Finding

The positive central vorticity term `Omega K t^2` does **not** produce circulation gain in a material transverse core. Its entire second-order gain is canceled by the core's area contraction. In fact, the viscous material circulation has vanishing derivatives through order four at the initial time:

`Gamma(t)=Gamma(0)+O(t^5)`.

This stronger statement follows from the initially rigid rotation and harmonic local acceleration; it is a finite Taylor-jet statement, not exact conservation for positive time. A fifth-order viscous term is the first not excluded by the calculation below, and its sign is not established.

The pressure-generated **strain** does have a positive initial gain. Thus the datum supplies an onset of strain that can act on a different receiving perturbation, while it does not yet supply the increasing circulation Reynolds number of the same material vortex packet.

## 1. The material transverse plane is invariant

Write the initial datum as in PASS2:

`u_0=Omega psi(|x|^2)(-y,x,0)+A v_remote`,

with the paired packets centered at `z=±d`. The initial horizontal velocity is even under `z -> -z`, and the vertical velocity is zero, hence odd under that reflection. NS and local uniqueness preserve the vector reflection symmetry. Therefore

`u_z(t,x,y,0)=0`, `partial_z u_perp(t,x,y,0)=0`.

The plane `z=0` is a material plane, and its vorticity is normal to that plane. Full inversion symmetry also keeps the origin fixed as a material particle.

Choose a disk `D_r={z=0,x^2+y^2<r^2}` whose closure lies strictly inside the initial rigid-rotation region. For sufficiently small positive time its material image `D_r(t)` is a smooth planar disk, generally deformed rather than circular. Let `X(t,a)` be its horizontal material map and

`j(t,a)=det D_a X(t,a)>0`.

The full flow is volume preserving. The horizontal area factor instead obeys

`d_t j(t,a)=-lambda(t,X(t,a),0) j(t,a)`,

where `lambda=partial_z u_z`, since `div_perp u_perp=-lambda` on the invariant plane.

## 2. Central deformation and exact second-order cancellation

At the origin PASS2 gives

`lambda(t,0)=K t+O(t^2)`,

`omega_z(t,0)=2 Omega[1+(K/2)t^2]+O(t^3)`,

where `K=A^2 C_v-(2/5)Omega^2>0`. Consequently

`j(t,0)=1-(K/2)t^2+O(t^3)`

and

`omega_z(t,0) j(t,0)=2 Omega+O(t^3)`.

For a direct deformation check, let `L=Omega J` be the horizontal initial rotation matrix and `H=partial_t grad_perp u_perp(0,0)`. The local acceleration is a gradient, so `H` is symmetric, and incompressibility gives `tr H=-K`. The horizontal deformation matrix has expansion

`F_perp(t)=I+tL+(t^2/2)(H+L^2)+O(t^3)`.

Its determinant is `1+(tr H/2)t^2+O(t^3)`, exactly the area factor above. Also

`F_perp(t)^T F_perp(t)=I+t^2 H+O(t^3)`.

Thus positive `K` guarantees contraction of **area**, not necessarily contraction of both principal transverse lengths. An anisotropic packet can have one expanding principal direction while its area decreases.

Define the central area radius `ell_A(t)=r sqrt(j(t,0))`. If the rotational velocity scale is defined by `U_rot=omega_z ell_A/2`, then

`Re_rot(t)=U_rot ell_A/nu=omega_z(t,0) r^2 j(t,0)/(2nu)`.

The second-order Reynolds gain is zero. This is a precisely defined rotational/core-flux scale, not a claim about the full velocity norm on the deformed boundary. Pressure creates an irrotational strain component already at order `t`, and the full velocity norm can reflect that additional component.

For a finite disk, put `K(a)=-partial_zz p_0(a,0)`, which is continuous and positive on sufficiently small disks around the origin. Pointwise in initial labels,

`j(t,a)=1-(K(a)/2)t^2+O(t^3)`,

`omega_z(t,X(t,a),0)=2Omega[1+(K(a)/2)t^2]+O(t^3)`.

The disk area satisfies

`|D_r(t)|=pi r^2-(t^2/2) integral_Dr K(a) da+O(t^3)`.

The same gain/contraction cancellation holds under the integral. Replacing `K(a)` by its central value for a finite disk would introduce an additional spatial approximation, which is unnecessary here.

## 3. Viscous Kelvin identity and the justified vanishing jets

Let `C_r(t)=boundary D_r(t)` and define the actual material circulation

`Gamma_r(t)=oint_Cr(t) u(t)·dx = integral_Dr(t) omega_z(t) dA`.

For smooth unforced NS its exact evolution is

`Gamma_r'(t)=nu oint_Cr(t) Delta u(t)·dx`

`             =nu integral_Dr(t) Delta omega_z(t) dA`.

Equivalently, the pulled-back local flux density obeys

`d_t[omega_z(t,X(t,a),0) j(t,a)]`

`    =nu j(t,a) Delta omega_z(t,X(t,a),0)`.

The Laplacian here is the full three-dimensional Laplacian. Replacing it by a two-dimensional cross-sectional Laplacian would lose the axial diffusion term.

Here is a local jet proof that this right side is `O(t^4)`. Write `u_n=partial_t^n u|_(t=0)`, `omega_n=curl u_n`, and `p_n=partial_t^n p|_(t=0)`. On an open initial core neighborhood,

`u_0=L_3 x`, `L_3=Omega diag_block(J,0)`, `omega_0=b=2Omega e_z`,

with `L_3` constant and skew. We use `R=(L_3 x)·grad`, the infinitesimal rigid-rotation operator. It commutes with the Laplacian on scalar and vector components.

1. `omega_0` is constant and `omega_1=0` identically in this neighborhood. Thus `Delta omega_0=Delta omega_1=0` there.

2. The initial acceleration is

   `u_1=-L_3^2 x-grad p_0`.

   Since `Delta p_0=-tr(L_3^2)=2Omega^2` is constant, `Delta u_1=0`. Moreover `grad u_1` is symmetric. Therefore

   `omega_2=(b·grad)u_1`, and `Delta omega_2=0`.

3. Differentiating the pressure Poisson equation once gives

   `Delta p_1=-2 tr(L_3 grad u_1)=0`,

   because a skew matrix has zero trace pairing with a symmetric matrix. The next velocity jet is

   `u_2=-L_3 u_1-Ru_1-grad p_1+nu Delta u_1`.

   Hence `Delta u_2=0`, using the commutation of `Delta` and `R`. Differentiating the vorticity equation twice yields

   `omega_3=-R omega_2+L_3 omega_2+(b·grad)u_2+nu Delta omega_2`.

   Every term is harmonic locally, so `Delta omega_3=0`.

Thus `partial_t^j Delta omega|_(t=0)=0` locally for `j=0,1,2,3`. Since the solution is smooth up to its initial time, Taylor expansion on a smaller compact core gives

`Delta omega(t,x)=O(t^4)`.

Material transport and the area factor do not change this order. The exact viscous identity consequently implies

`Gamma_r^(j)(0)=0` for `j=1,2,3,4`,

`Gamma_r(t)=2Omega pi r^2+O(t^5)`.

The same statement holds pointwise for pulled-back flux density:

`omega_z(t,X(t,a),0)j(t,a)=2Omega+O(t^5)`.

One can identify the first not-yet-excluded coefficient without claiming its sign:

`Gamma_r^(5)(0)=nu integral_Dr Delta omega_4,z(a,0) da`.

Terms from differentiating the moving surface vanish because all preceding Laplacian jets vanish on a neighborhood. No assertion is made that this fifth derivative is nonzero, positive, or independent of viscosity. These calculations do not justify vanishing of higher circulation jets or exact circulation conservation for any positive time interval.

## 4. What must change in a receiving core

The circulation Reynolds number is

`Re_Gamma(t)=Gamma_r(t)/(2pi nu)=Omega r^2/nu+O(t^5)`.

If instead one defines the mean normal vorticity on the material disk by `omega_bar=Gamma_r/|D_r(t)|`, its area radius by `ell_r=sqrt(|D_r(t)|/pi)`, and its rotational velocity scale by `omega_bar ell_r/2`, the resulting Reynolds number is exactly `Gamma_r/(2pi nu)`. Area contraction cannot increase it independently of circulation.

There is a different variable with a positive onset: a strain-based receiver parameter,

`Re_strain(t)=lambda(t,0) ell_A(t)^2/nu`

`            =(K r^2/nu)t+O(t^2)`.

This can improve the initial growth-versus-diffusion balance for a separately specified wave or seed. It measures the strain furnished by the coupled outer packets, not a gain in the material core's own circulation. Turning it into a regenerative stage requires quantitative persistence and feedback estimates that are not supplied by the second jet.

There is also no contradiction with a gain on a **fixed Eulerian** rotational receiver. In the central infinitesimal-core approximation, write the vorticity factor as `G=1+Kt^2/2+O(t^3)`. A smaller Eulerian radius `q r` can have `q^2 G>1` when `q` is sufficiently close to one. To second order this requires `q^2>1-Kt^2/2`: that receiver shrinks more slowly than the original material area radius. Its enlarged initial material footprint collects neighboring circulation. This observation neither proves a finite prescribed gain nor supplies an infinite reservoir; it explains why a material-loop obstruction and a fixed-Eulerian modal gain can both be correct.

A proposed transfer must therefore explicitly change at least one of the following:

- the receiving loop/surface, so it is not merely the same material circulation packet shrinking under stretching, with the additional intercepted flux and return strands accounted for;
- the role of the receiver, using pressure-generated strain to amplify a new perturbation and proving how the strain/localization configuration is replenished;
- the circulation itself through the actual viscous line/flux integral, with its sign and size proved rather than inferred from positive vorticity.

This identifies the gap in this datum's attempted use as a restart. It does not rule out other three-dimensional mechanisms or establish a sign restriction on general viscous Kelvin circulation.

## Verification

`python3 work/pass3/verify_material_core_transfer.py` verifies the deformation determinant and metric expansions, second-order flux cancellation, the rigid-rotation/Laplacian commutator, skew–symmetric pressure pairing, and harmonic jet identities. The argument above then uses exact NS and material-transport identities. It does not numerically integrate the flow, prove uniform-time amplification, or build a cascade.

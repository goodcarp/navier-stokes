# Inherited packet geometry and pump strength at the certified datum

The certified initial pressure feedback can coexist with positive geometric margins for a short stage. It does not preserve the inserted packet shape: the packet moves toward the origin, expands transversely, and loses local angular velocity relative to the core. Ordinary viscosity further broadens its moment-based aspect ratios. These are actual initial Navier--Stokes identities, rather than derivatives of a manually reset family.

Throughout, use cylindrical radius `r`, axial coordinate `z`, spherical radius `R`, and the actual axisymmetric solution. The initial parameters are

\[
b=\Omega=1,\qquad c=\tfrac75,\qquad \nu=\tfrac1{1000},
\]

with a radial annular strain pump `c S_eta`, `eta=-1` on `5/2<=R<=5`. The exterior swirl packet has

\[
w_0^+(r,z)=A r\,\chi(r^2/a^2)\chi((z-d)^2/h^2),
\quad a=\tfrac{63}{200},\quad h=\tfrac9{20},\quad d=4,
\]

plus its reflected partner. `A` is the positive neutral-tuning amplitude. Both packets lie strictly inside the pump plateau, where the actual initial meridional velocity is exactly

\[
u_{r,0}=cr,\qquad u_{z,0}=-2cz.
\tag{1}
\]

The core velocity vanishes there. The swirl changes azimuthal motion but does not change these initial meridional trajectories.

## 1. Material geometry and finite margins

For material labels initially in the packet, the initial rates are

\[
r'=cr,\quad z'=-2cz,\quad d'=-2cd,
\quad a'=ca,\quad h'=-2ch.
\tag{2}
\]

Here `a,h` in (2) refer to the boundary of the **transported set of initial labels**. They are not support radii of the viscous solution. At `nu>0`, an azimuthal packet generally has immediate diffusive tails, so exact compact packet/core separation cannot be reused at positive time.

Initially the packet satisfies the strict bounds

\[
R\ge d-h=3.55>2.5,\qquad
R^2\le(d+h)^2+a^2<25,\qquad
\frac r{|z|}\le\frac{63}{710}<\frac14.
\tag{3}
\]

Thus there is a positive interval on which the transported packet-label set stays inside the original annulus and cone, by continuity of the smooth flow. This is compatible with the initially positive normalized swirl-stress derivative. It does not supply a numerical stage duration: that requires bounds on velocity/acceleration and on the chosen concentration error. The initial margins alone cannot be extrapolated with a frozen linear flow.

The transported plateau boundaries also deform anisotropically. For a point on an initial sphere,

\[
\frac{R'}R=c(1-3\mu^2),\qquad \mu=z/R.
\tag{4}
\]

The angular variation in (4) cannot be removed by a single isotropic dilation. Following both packet and plateau labels preserves their inclusion for short time, but does not preserve an exactly spherical, affine pump profile.

## 2. Actual viscous packet moments

There is a useful actual-state observable that avoids pretending the support remains compact. Given the actual meridional velocity, the azimuthal equation is linear in the swirl. Evolve the upper initial packet as a tagged solution `w^+` of this equation. The sum of tagged solutions recovers the actual swirl, including the separately tagged core and lower packet. The tags introduce no external force and do not alter the velocity used in transport.

Writing `F=w^+/r`, the equation is

\[
\partial_tF+u_r\partial_rF+u_z\partial_zF
=-2\frac{u_r}{r}F
+\nu\left(\partial_{rr}F+\frac3r\partial_rF+\partial_{zz}F\right).
\tag{5}
\]

For smooth axisymmetric solutions and vanishing moment boundary terms, the positive angular-momentum measure

\[
dM=2\pi r^2w^+\,dr\,dz=2\pi r^3F\,dr\,dz
\]

has constant total mass. The following identities follow by integration by parts; in particular they hold at the compact initial datum. Normalize `dM` to a probability measure and write expectation `E`. Define

\[
D=E[z],\quad \mathsf H^2=E[(z-D)^2],\quad
\mathsf R^2=E[r^2].
\]

For general later meridional velocity, the exact moment equations are

\[
D'=E[u_z],\qquad
(\mathsf R^2)'=2E[ru_r]+8\nu,\qquad
(\mathsf H^2)'=2E[(z-D)u_z]+2\nu.
\tag{6}
\]

The coefficient `8 nu` comes from the four-dimensional radial diffusion operator for `F`, and includes the cylindrical basis contribution. At the inserted datum, (1) applies on all initial tagged mass, so

\[
\boxed{D'=-2cD,\quad
(\mathsf R^2)'=2c\mathsf R^2+8\nu,\quad
(\mathsf H^2)'=-4c\mathsf H^2+2\nu.}
\tag{7}
\]

Here `D(0)=d`. The initial widths are explicitly determined by the unchanged cutoff:

\[
\mathsf R^2(0)=a^2\frac{\int_0^1\xi^2\chi(\xi)d\xi}
{\int_0^1\xi\chi(\xi)d\xi},\qquad
\mathsf H^2(0)=h^2\frac{\int_{-1}^1 y^2\chi(y^2)dy}
{\int_{-1}^1\chi(y^2)dy}.
\]

Consequently two dimensionless shape ratios satisfy

\[
\boxed{\left.\frac d{dt}\log\frac{\mathsf R}{D}\right|_0
=3c+\frac{4\nu}{\mathsf R^2}>0,\qquad
\left.\frac d{dt}\log\frac{\mathsf H}{D}\right|_0
=\frac{\nu}{\mathsf H^2}>0.}
\tag{8}
\]

These ratios are unaffected by any isotropic spatial or velocity normalization. They therefore measure a real inherited change in the packet geometry.

## 3. Relation to the shrinking receiving core

For the receiving-RMS scale match with exponent `alpha in (0,1/2)`, the neutral initial core gives

\[
L(t)=L(0)(1-\kappa t+O(t^2)),\qquad
\kappa=\frac2{\alpha+2}\in(\tfrac45,1).
\]

Combining this with (7), the initial normalized rates are

\[
\begin{aligned}
(\log(D/L))'&=\kappa-2c\in(-2,-\tfrac95),\\
(\log(\mathsf R/L))'&=\kappa+c+\frac{4\nu}{\mathsf R^2}
>\tfrac{11}{5},\\
(\log(\mathsf H/L))'&=\kappa-2c+\frac{\nu}{\mathsf H^2}.
\end{aligned}
\tag{9}
\]

The packet approaches the core in normalized axial distance and broadens in normalized transverse radius. The axial-width rate includes diffusion and is stated without discarding that term. For material label widths, remove the diffusion terms in (9).

## 4. Stress growth does not imply packet rotation or pump renewal

At the axis point initially `(0,d)`, follow its actual material trajectory and define local packet rotation `A_loc` and meridional strain `c_loc=-partial_z u_z/2`. Initially the velocity gradient is `diag(c,c,-2c)+A J`, and all local viscous derivatives vanish because the core of this packet and the pump plateau are affine. The exact gradient equation gives

\[
\boxed{A_{\rm loc}'=-2cA,\qquad
c_{\rm loc}'=2c^2+\tfrac12p_{zz}(0,d).}
\tag{10}
\]

The pressure in (10) is the full actual initial pressure at the exterior point. Its value is not the certified pressure derivative at the origin. The central rotation satisfies `Omega'=2b Omega`, so

\[
\left.\frac d{dt}\log\frac{A_{\rm loc}}\Omega\right|_0=-2(c+b)=-\tfrac{24}{5}.
\tag{11}
\]

To keep the **local pump/core strain ratio** from decreasing at insertion requires

\[
\boxed{p_{zz}(0,d)\ge4bc-4c^2=-\tfrac{56}{25}.}
\tag{12}
\]

No such exterior pressure bound is supplied by the origin feedback certificate. Moreover, a single axis-point condition does not control the pump over the full packet neighborhood.

The exact affine meridional plateau also immediately acquires a nonaffine component on the axial swirl transition. With the convention `omega_theta=partial_z u_r-partial_r u_z`, the actual equation gives there

\[
\left.\partial_t\omega_\theta\right|_0
=\frac1r\partial_z(w_0^2).
\tag{13}
\]

The mixed pressure derivatives cancel in (13), and the initially affine meridional viscous term vanishes. The right side is nonzero in the packet's axial cutoff region. Thus reinserting the original exact affine-plateau formulas at the next time would omit a real inherited meridional vorticity component.

This initial deformation has a large explicit coefficient. At

    r=a/2,  z=d+h sqrt(5/8),

the radial cutoff is one, while the axial cutoff and its derivative are exactly 1/2 and -8/3. Equation (13) consequently gives

    partial_t omega_theta(0,r,z)=-(28/15) A² sqrt(5/8).

Since sqrt(5/8)>3/4 and the certified upper bound is C_v<279/40000000, the neutral amplitude yields

    |partial_t omega_theta(0,r,z)| > (7/5)A²
       > 108800000/93 > 1,169,892.

This is an exact initial derivative at the stated point in the affine pump plateau. It is not a lower bound on later vorticity, an upper bound on the stage lifespan, or a proof that useful finite transfer fails. It quantifies a deformation that a continuation estimate must retain; the large exterior amplitude is not a passive reservoir.

In contrast, the previously audited initial outer-swirl stress calculation gives, for a smooth spatial window equal to one near the packets and zero near the core,

\[
\left.\frac d{dt}\log\frac{C_\theta}{\Omega^2}\right|_0
\ge\gamma c-4b-\nu\frac{d_v}{C_v},\qquad
\gamma=\frac{2381}{357}.
\]

The existing interval certificate `d_v/C_v<901` makes this lower bound strictly larger than `4` at the present parameters. This windowed observable is well defined at positive time; the simple initial formula does not discard the later window-flux terms. Its initial increase is compatible with (11): transport toward the origin increases the pressure weight while the packet's local rotation decays. It is not evidence that the exterior packet's own velocity, circulation density, or strain profile has been regenerated.

## 5. A measurable necessary return condition

For the same inherited upper packet, (6) implies the exact aspect-ratio budget

\[
\log\frac{(\mathsf R/D)(T)}{(\mathsf R/D)(0)}
=\int_0^T\left[
\frac{E[ru_r]}{\mathsf R^2}
-\frac{E[u_z]}D
+\frac{4\nu}{\mathsf R^2}\right]dt,
\tag{14}
\]

as long as `D>0` and these moments exist. To return to a class with `mathsf R/D<=beta_max`, the right side must not exceed `log(beta_max/beta_initial)`. Returning to exactly the same aspect ratio requires zero net budget. At insertion the integrand is strictly positive by (8). Sustained axial compression and transverse expansion cannot repeatedly restore the same shape without an offsetting later transport regime, a proved transfer into a different receiving packet, or a class allowing controlled cumulative shape change. Merely shrinking coordinates leaves (14) unchanged.

The moment condition is necessary rather than sufficient. A useful full return criterion must also control concentration outside the pump/cone window, the inherited meridional vorticity in (13), local pump strength, the complete velocity profile in an appropriate norm, and the pressure feedback margin. A strictly positive initial shape derivative alone does not exclude all later returns: the actual meridional field can change, and no lower bound over every subsequent stage has been proved here.

`verify_inherited_geometry.py` checks the swirl-density adjoint, the exact moment rates, the local matrix identities, the pressure threshold, and the initial rational margins. It does not integrate the PDE or provide a stage duration.

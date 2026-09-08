# Fully three-dimensional angular-momentum and mean-torque budget

Nonaxisymmetry removes the scalar maximum principle for `Gamma=r u_theta` through a pressure torque and a vector-viscosity coupling. Its azimuthal mean has a more restrictive budget: both direct terms average to zero, and transfer into the mean comes from quadratic velocity correlations. A compact, explicitly divergence-free single Fourier mode can give a selected-sign mean torque at insertion. This is an actual initial-data construction, not a later-time return mechanism.

All equations below concern smooth **unforced** incompressible Navier--Stokes with ordinary positive viscosity. This choice analyzes the present unforced datum; it does not forbid a separate smooth-forcing research route. Write

\[
D_t=\partial_t+u_r\partial_r+\frac{u_\theta}{r}\partial_\theta+u_z\partial_z,
\qquad \Gamma=r u_\theta.
\]

In a fully three-dimensional flow, `Gamma` is a pointwise axial angular-momentum component. The circulation around a fixed coordinate circle is `2 pi` times its azimuthal mean, rather than `2 pi Gamma` at an individual angle.

## 1. The full angular equation

The theta-component momentum equation is

\[
D_tu_\theta+\frac{u_ru_\theta}{r}
=-\frac1r\partial_\theta p
+\nu\left(\Delta u_\theta-\frac{u_\theta}{r^2}
+\frac2{r^2}\partial_\theta u_r\right).
\]

Multiplication by `r`, including `D_t r=u_r`, gives the exact identity

\[
\boxed{D_t\Gamma=-\partial_\theta p+
\nu\left(\Gamma_{rr}-\frac1r\Gamma_r+\Gamma_{zz}
+\frac1{r^2}\Gamma_{\theta\theta}
+\frac2r\partial_\theta u_r\right).}
\tag{1}
\]

Equivalently, in Cartesian coordinates `Gamma=x u_y-y u_x`,

\[
D_t\Gamma=-\partial_\theta p+\nu(\Delta\Gamma-2\omega_z),
\quad \omega_z=\frac1r(\Gamma_r-\partial_\theta u_r).
\]

At an interior positive maximum of `Gamma`, the scalar diffusion terms are nonpositive, but `-p_theta+(2 nu/r) partial_theta u_r` has no fixed sign. Under axisymmetry both terms vanish, recovering the scalar maximum-principle argument. The angular scalar diffusion `nu Gamma_thetatheta/r^2` itself remains dissipative; it is the additional coupling to `u_r`, not that scalar term, which can act as a signed source.

Pressure still has zero circulation around every closed loop. For a material loop, the unforced Kelvin identity is `d/dt integral u dot dl = nu integral Delta u dot dl`. Equation (1) does not create a new pressure term in that material-loop identity. The coordinate-circle and pointwise angular-momentum observables have different transport budgets.

## 2. Azimuthal mean and a single Fourier mode

Let an overbar mean `(2 pi)^(-1) integral_0^(2 pi) dtheta`, let `U=bar u`, and write `v=u-U`. The mean meridional field obeys `(r U_r)_r/r+(U_z)_z=0`. Set `G=bar Gamma=r U_theta`. Averaging the conservative transport form of (1) gives

\[
\boxed{\partial_tG+U_rG_r+U_zG_z
=\nu\left(G_{rr}-\frac1rG_r+G_{zz}\right)+\mathcal T,}
\tag{2}
\]

where

\[
\boxed{\mathcal T=-\frac1r\partial_r\left(r^2\overline{v_rv_\theta}\right)
-r\partial_z\overline{v_zv_\theta}.}
\tag{3}
\]

Both `bar p_theta` and `overline(partial_theta u_r)` vanish exactly. Pressure can change the correlations through the subsequent three-dimensional evolution; there is no direct pressure-torque term in (2). The divergence-free constraint and periodicity also give the useful equivalent formula

\[
\mathcal T=-\overline{v_rv_\theta}
-r\overline{v_r\partial_rv_\theta+v_z\partial_zv_\theta}.
\tag{4}
\]

For a single real perturbation

\[
v=\varepsilon\operatorname{Re}\big(V(r,z)e^{im\theta}\big),\qquad m\ne0,
\]

its solenoidal constraint is `(r V_r)_r/r+(im/r)V_theta+(V_z)_z=0`, and

\[
\mathcal T=-\frac{\varepsilon^2}{2r}\partial_r
\left(r^2\operatorname{Re}(V_r\overline{V_\theta})\right)
-\frac{\varepsilon^2r}{2}\partial_z
\operatorname{Re}(V_z\overline{V_\theta}).
\tag{5}
\]

In particular, adding a mean-zero mode to a fixed axisymmetric datum changes the initial mean torque at order `epsilon^2`, with no order-`epsilon` base/perturbation contribution. This does not imply that the pointwise nonaxisymmetric pressure response is quadratic: that response can already have a linear Fourier component. After positive time, the nonlinear flow generally creates other Fourier modes, so (5) is not a closed one-mode evolution model.

For a compact perturbation, (3) also implies `integral T r dr dz=0`. Its positive local torque is balanced by negative torque elsewhere. This conclusion is an exact initial compact-data identity; later global angular-momentum conservation requires the corresponding moment boundary terms to vanish.

## 3. Compact seed with a selected-sign mean torque

Choose integer `m>=1`, nonzero real `k`, and a real smooth compact envelope `a(r,z)` supported in `r0-delta<r<r0+delta`, with `0<delta<r0`. Define the globally smooth vector potential and velocity

\[
\mathcal A=a(r,z)\cos(m\theta+kr)e_z,\qquad v=\nabla\times\mathcal A.
\]

The support avoids the axis, so extension by zero near the axis causes no coordinate singularity. Explicitly, with `phi=m theta+kr`,

\[
v_r=-\frac{ma}{r}\sin\phi,\qquad
v_\theta=-a_r\cos\phi+ka\sin\phi,\qquad v_z=0.
\tag{6}
\]

This is exactly divergence free, has zero azimuthal mean, and is a single `m` mode. A concrete unchanged-cutoff envelope is

\[
a(r,z)=a_0\chi((r-r_0)^2/\delta^2)\chi((z-z_0)^2/h^2).
\]

Its Fourier correlation and mean torque are exactly

\[
\overline{v_rv_\theta}=-\frac{mk}{2r}a^2,\qquad
\boxed{\mathcal T=\frac{mk}{2r}\partial_r(ra^2).}
\tag{7}
\]

On the nonzero flat part of the envelope, `T=mk a0^2/(2r)`. Its sign is selected by `mk`; for the actual perturbation `epsilon v`, multiply (7) by `epsilon^2`. The compensating torque is in the envelope transition. Neither pressure nor viscosity supplies an additional initial *mean* torque for this mode; its mean viscous velocity is zero.

High Fourier or radial-phase frequencies do not give unlimited torque at fixed local velocity energy. On the flat part,

\[
\overline{|v|^2}=\frac{a_0^2}{2}\left[(m/r)^2+k^2\right],\qquad
|\mathcal T|\le\tfrac12\overline{|v|^2}.
\tag{8}
\]

Equality holds in (8) when `|k|=|m|/r`. Holding the vector-potential amplitude fixed while increasing frequency also increases velocity and is not a small perturbation at fixed velocity norm.

## 4. Necessary torque budget for a Gamma-carrying scale return

Suppose the intended receiving component has `U~L^(-1-alpha)`, so `Gamma~L U~L^(-alpha)` with `alpha>0`. This section applies to a scheme that actually requires a growing mean rotational/angular-momentum component. It is not a universal necessary condition on every conceivable three-dimensional blowup mechanism, since a flow can have large mean-zero or meridional velocity.

Let `kappa=-L'/L>0` and write `G(t,L rho,L zeta)=L^(-alpha) H(t,rho,zeta)`. Equation (2) becomes

\[
\partial_tH+\left(\frac{U_{r,z}}L+\kappa(\rho,\zeta)\right)\cdot\nabla H
+\alpha\kappa H
=\frac\nu{L^2}\left(H_{\rho\rho}-\frac1\rho H_\rho+H_{\zeta\zeta}\right)
+L^\alpha\mathcal T(t,L\rho,L\zeta).
\tag{9}
\]

For a fixed positive rescaled peak, (9) requires at least `T>=alpha kappa G` at that peak; a negative diffusion contribution increases the necessary source. More generally, for a smooth bounded decaying mean with positive maximum `M(t)=sup G`, the parabolic comparison bound is

\[
M(T)\le M(0)+\int_0^T\|\mathcal T_+(t)\|_\infty\,dt.
\tag{10}
\]

Thus a claimed increase of the **mean maximum** by `q^(-alpha)` needs integrated positive torque at least `(q^(-alpha)-1) M(0)`. A smaller local receiving core can initially borrow a larger pre-existing exterior mean maximum; applying this multiplier bound to that local value instead of the global maximum would be invalid. An infinite Gamma-growing sequence eventually needs an unbounded actual torque budget or another proved change of mechanism.

There is also an exact smallness estimate. In a region `r<=C L`, if `|v|<=epsilon U` and `|(partial_r,partial_z)v_theta|<=K epsilon U/L`, (4) gives

\[
|\mathcal T|\le(1+CK)\varepsilon^2U^2.
\tag{11}
\]

Consequently, on the turnover time `tau~L/U`, where `kappa~U/L`, sustaining an order-one prescribed mean-growth rate requires a corresponding order-`epsilon^2` budget. For fixed profile constants, arbitrarily small `epsilon` cannot be counted as a leading-order mean torque. A longer stage or larger derivative norm changes this comparison and must be controlled in the actual dynamics.

More precisely, combining the fixed-peak requirement with (11) gives the necessary inequality

\[
\varepsilon^2\ge\frac{\alpha\kappa G}{(1+CK)U^2}
=\frac{\alpha}{1+CK}\left(\frac{\kappa L}{U}\right)
\left(\frac{G}{LU}\right),
\]

even before paying a strictly negative diffusion term. This implication uses a fixed positive rescaled peak as in (9), not an arbitrary transient local value borrowed from an exterior reservoir.

For profiles with ordinary curvature at scale `L`, the viscous cost has scale `nu Gamma/L^2=U^2/Re`, `Re=UL/nu`; the mean torque has scale `epsilon^2 U^2`. This is a scaling comparison, not a proved lower bound on diffusion for every profile. Raising a wavenumber also changes velocity gradients and viscous decay, so an amplitude-small, increasingly oscillatory perturbation is not automatically small in the norm needed for a continuation or return argument.

## 5. Scope of the construction

The seed supplies a concrete nonaxisymmetric direction with positive initial mean torque and an explicit compensating transfer. It does not prove favorable pressure at the core, persistence of the phase correlation, preservation of the receiving profile, or a repeated stage. Equation (2), the quadratic size (5), and the inherited-state evolution must all remain in a prospective return argument. The separate rotational-receiver calculation can test whether this torque has a useful projection into the chosen receiving region.

`verify_nonaxisymmetric_torque.py` checks the coordinate identity, mean-flux reduction, compact-seed solenoidality and covariance, selected-sign torque, and the exact local energy bound. These are symbolic identities, not a numerical evolution or formal-kernel proof.

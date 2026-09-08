# Whole-space exterior pressure and the initial meridional deformation scale

The large outer swirl does not contribute its bare centrifugal force to the meridional evolution. Its actual pressure cancels most of that force in the diagnostic below. The remaining response is still large compared with the initially prescribed pump strain. This note gives exact whole-space pressure formulas, conservative analytic bounds, and explicitly exploratory values. It does not certify a later-time deformation or a stage duration.

## 1. Exact centrifugal cancellation for the unit exterior swirl

Use the unchanged pass 7 datum and cylindrical radius `r`. Write the unit exterior swirl as

\[
v_\theta=r C(r)F(z),\quad C(r)=\chi(r^2/a^2),\quad a=63/200,
\]
\[
F(z)=\chi((z-4)^2/h^2)+\chi((z+4)^2/h^2),\quad h=9/20.
\]

The two axial bumps are disjoint. Put `f=F^2`, and define

\[
H(r)=\int_r^a sC(s)^2ds,\quad \rho(r,z)=H(r)f(z),
\quad \Phi=N*\rho,\quad N=(4\pi|x|)^{-1}.
\]

The exact source is

\[
g_v=-\frac1r\partial_r(v_\theta^2)
=\Delta_h\rho,
\]

where `Delta_h=partial_rr+r^(-1)partial_r` on radial functions. Since `-Delta Phi=rho`, the decaying whole-space pressure is exactly

\[
\boxed{p_v=-\rho-\partial_{zz}\Phi.}
\tag{1}
\]

No boundary pressure or infinite-cylinder approximation is imposed. Set `delta=p_v+rho=-partial_zz Phi`. At a point on the axial plateau `z=4`, `f=1` and all its derivatives vanish locally, so

\[
\begin{aligned}
p_{v,r}&=rC^2+\delta_r,\\
p_{v,rr}&=C^2+2rCC_r+\delta_{rr},\\
p_{v,z}&=\delta_z,\qquad p_{v,zz}=\delta_{zz}.
\end{aligned}
\tag{2}
\]

Therefore the unit swirl's meridional acceleration contribution is

\[
\boxed{a_{v,r}=v_\theta^2/r-p_{v,r}=-\delta_r,
\qquad a_{v,z}=-\delta_z.}
\tag{3}
\]

The bare term `r C^2` has cancelled exactly. It is only the remaining nonlocal axial-localization response that drives (3). The function `delta` is harmonic on the entire source-free axial slab `|z-4|<h/2`, hence

\[
\delta_{rr}+\delta_r/r+\delta_{zz}=0
\tag{4}
\]

there. This gives an independent numerical check of the computed pressure derivatives.

## 2. Nonsingular compact endcap integrals

An integration by parts in the source axial coordinate gives

\[
\boxed{\delta(x)=-\int_{\mathbb R^3}H(s)f'(\zeta)
\partial_zN(x-y)\,dy.}
\tag{5}
\]

The derivative `f'` is supported only in the four endcaps `h/2<=|zeta-4|<=h` and `h/2<=|zeta+4|<=h`. At the evaluation circle `z=4`, every source in (5) has vertical separation at least `h/2=9/40`. All further derivatives needed in (2) can therefore be taken under a nonsingular compact integral.

For target `(r,0,z)`, source `(s cos theta,s sin theta,zeta)`, set

\[
d=z-\zeta,\quad X=r-s\cos\theta,\quad
q=r^2+s^2-2rs\cos\theta,\quad R_y^2=q+d^2.
\]

The kernels used by the evaluator are exactly

\[
\begin{aligned}
N_z&=-d/(4\pi R_y^3),\\
N_{zr}&=3dX/(4\pi R_y^5),\\
N_{zz}&=(2d^2-q)/(4\pi R_y^5),\\
N_{zzz}&=3d(3q-2d^2)/(4\pi R_y^7),\\
N_{zrr}&=3d(R_y^2-5X^2)/(4\pi R_y^7).
\end{aligned}
\tag{6}
\]

Thus `(delta,delta_r,delta_z,delta_zz,delta_rr)` equals minus the integral of `H f'` times the five kernels in (6), with measure `s ds dtheta dzeta`.

One can instead integrate the ring angle analytically. Its exact Newton kernel is

\[
\int_0^{2\pi}N\,d\theta
=\frac{K(k)}{\pi\sqrt{(r+s)^2+d^2}},\qquad
k=\frac{2\sqrt{rs}}{\sqrt{(r+s)^2+d^2}},
\]

where `K(k)=integral_0^(pi/2)(1-k^2 sin^2 t)^(-1/2)dt`. This leaves a two-dimensional nonsingular endcap integral. The present exploratory evaluator directly performs the smooth angular integral; either representation supports interval quadrature without a large pressure solver.

The upper packet is exactly symmetric about `z=4`. It therefore contributes `p_{+,z}(r,4)=0`. Its axial second derivative need not vanish. The small nonzero first axial derivative at this circle comes from the reflected lower packet.

## 3. Conservative bounds on the reflected packet

The radial density mass obeys the exact Fubini identity and bound

\[
\int_{\mathbb R^2}H(|y_h|)dy_h
=\pi\int_0^a s^3C(s)^2ds\le\frac{\pi a^4}{4}.
\]

The total variation of one axial squared bump is two. For the lower packet observed at `z=4`, the gap is `g=8-h=151/20`. From

\[
|N_{zr}|\le\frac3{8\pi R_y^3},\quad
|N_{zz}|\le\frac1{2\pi R_y^3},\quad
|N_{zzz}|\le\frac3{2\pi R_y^4},
\]

one obtains

\[
\begin{aligned}
|\delta_{-,r}|&\le\frac{3a^4}{16g^3}<4.290\cdot10^{-6},\\
|\delta_{-,z}|&\le\frac{a^4}{4g^3}<5.720\cdot10^{-6},\\
|\delta_{-,zz}|&\le\frac{3a^4}{4g^4}<2.273\cdot10^{-6}.
\end{aligned}
\tag{7}
\]

These are rigorous but conservative bounds, independent of the exploratory quadrature. The same formulas with `g=h/2` bound the nearby upper packet; those much broader bounds do not by themselves determine a useful deformation time or the sign of its Hessian.

## 4. Whole-space diagnostic at the maximizing circle

The evaluator solves `chi(xi)+xi chi'(xi)=0` numerically and finds

\[
\xi_*=0.46792364389312446\ldots,\qquad
r_*=0.21547557533348247\ldots.
\]

No uniqueness proof or interval location for that root is asserted here. The earlier rigorous localization of every maximizing radius remains `(a/2,a sqrt(5/8))`.

The nonsingular ring/endcap quadratures at 48 and 96 nodes give the following unit-swirl values. These are **exploratory numerical values, not certified intervals**:

| Quantity at `(r_*,4)` | 96-node value |
|---|---:|
| Bare centrifugal term `r C^2` | `0.168428941855066` |
| Full pressure derivative `p_v,r` | `0.162065114234723` |
| Residual outward acceleration `a_v,r` | `0.006363827620343` |
| `p_v,zz` | `0.030158603663512` |
| `p_v,z` | `0.000000191039833465` |
| `partial_r a_v,r` | `0.000624733708675` |

The 48/96 differences are below `2e-11` for these entries, and the harmonic trace check (4) is below `2e-16`. Agreement is a diagnostic, not an error estimate. The computation includes the whole lower packet and has no exterior pressure boundary or missing radial tail: (5) is an exact compact-source representation of the whole-space pressure.

In this diagnostic, about `96.22%` of the centrifugal force is cancelled. The upper packet contributes `p_+,zz approximately0.0301586992`; the lower correction is approximately`-9.56e-8`. The upper axial first derivative is exactly zero by symmetry.

## 5. What this says about a 10^(-3) stage

For the actual family, multiply this pressure/acceleration contribution by the exact `A_lambda^2`, where `780<A_lambda<1020`. Using the diagnostic unit values, that component would contribute

\[
A_\lambda^2 a_{v,r}\ \text{approximately }3872\text{ to }6621,
\]

and

\[
\partial_z a_{v,z}=-A_\lambda^2p_{v,zz}
\ \text{approximately }-18348\text{ to }-31377.
\]

The initial mean axial strain is `partial_z U_z=-2c=-2.8`. A frozen first-derivative comparison changes that strain by its own initial size on roughly `0.9e-4` to `1.5e-4`. This indicates that the initial affine-pump velocity profile is not a slowly changing input on a `10^(-3)` interval.

Geometric displacement is a different estimate. A frozen second-order trajectory contribution over `10^(-3)` gives a fractional radial change `A^2 a_v,r t^2/(2r_*)` of about `0.9%` to `1.5%`, and an axial-width contribution of about `0.9%` to `1.6%`. These numbers are component-only Taylor diagnostics. They do not prove the actual deformation at that time, because subsequent acceleration, seed interactions, and pressure evolution are not frozen. Rapid change of an initial strain does not by itself prove catastrophic displacement or failure of a useful stage.

## 6. Retaining the rest of the actual datum

The analysis above isolates one exact component, not the total meridional acceleration. For

\[
u_0=M+A_\lambda v+\lambda w,
\]

azimuthal averaging at the exterior plateau gives the exact decomposition

\[
\begin{aligned}
\partial_t\overline{u_r}
&=-c^2r-p_{M,r}
+A_\lambda^2a_{v,r}
+\lambda^2[-\overline{(w\cdot\nabla)w}_r-\overline{p_w}_r],\\
\partial_t\overline{u_z}
&=-4c^2z-p_{M,z}
-A_\lambda^2p_{v,z}-\lambda^2\overline{p_w}_z.
\end{aligned}
\tag{8}
\]

Here the cylindrical radial convective component includes the centripetal basis term. Mean-zero cross modes vanish in (8), the mixed `M/w` pressure source is zero as proved in pass 7, and the initial mean meridional viscous term vanishes on the affine plateau. At the seed's flat envelope, its unprojected radial acceleration contribution is
`m^2/(2r^3)+k^2/(2r)`, but its own whole-space pressure must again be retained before calling that a surviving force.

For a material path in the full nonaxisymmetric flow, further advective and angular correlations also enter; (8) is the initial azimuthal-mean Eulerian equation. The unit-swirl computation alone does not certify the sign or size of every term in (8), nor a positive-time bound. A controlled `10^(-3)` stage needs the full inherited-state evolution and a remainder estimate, not an extrapolation of the large `A^2` term.

Files: `exterior-dynamics-pressure.py`, `exterior-dynamics-pressure-48.json`, `exterior-dynamics-pressure-96.json`, and `exterior-dynamics-verify.py`. The first three provide exploratory evaluation; the checker verifies the exact differential kernels, cancellation identity, and rational remote bounds.

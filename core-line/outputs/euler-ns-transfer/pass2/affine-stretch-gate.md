# PASS2: affine vortex stretching does not supply the increasing-Reynolds restart

2026-09-08. This is a stress test of the proposed nonaxisymmetric shrinking-core scale budget. It establishes exact statements in explicitly defined affine-flow classes and a conditional localization ceiling. It supplies neither a finite-energy singular solution nor a general NS regularity theorem.

## Outcome

An exact volume-preserving strain can increase vorticity and velocity while shrinking transverse length scales. For the circular stretched vortex and for isotropically compressed Kelvin modes, however, it does **not** increase the transverse Reynolds number needed by the earlier frequency-window proposal. Circulation fixes that Reynolds number; viscosity decreases the corresponding Fourier-mode quantity.

Anisotropic Kelvin waves can show increasing Reynolds number if it is measured using their weakly compressed wavelength. The same calculation shows the loss of the transverse localization margin: the gain consumes the packet's initial separation between wavelength and the strongly compressed envelope. It is a finite-stage resource, not an indefinitely repeatable restart.

The useful next target is therefore **regeneration of circulation/three-dimensional geometry together with a restored localization margin**, rather than another demonstration of affine vorticity gain.

## 1. Exact Kelvin-wave evolution on a volume-preserving strain

Let `S(t)` be a smooth symmetric trace-free `3 x 3` matrix. The affine velocity

\[
u_b(x,t)=S(t)x,
\qquad p_b(x,t)=-\tfrac12 x\cdot(S'+S^2)x
\]

solves the unforced NS equations at every viscosity. Let `F'=SF`, `F(0)=I`. Then `det F=1`.

A single real Kelvin wave `v=Re[a(t) exp(i k(t)·x)]` with `k·a=0` has identically zero self-advection. With the pressure mode chosen to enforce incompressibility, the full velocity `u_b+v` is an exact, generally infinite-energy NS solution when

\[
k'=-S^Tk,
\quad
a'=-Sa+2k\frac{k\cdot Sa}{|k|^2}-\nu|k|^2a.
\]

Because the base is irrotational, its wave-vorticity amplitude `b=i k x a` obeys the particularly simple equation

\[
b'=Sb-\nu|k|^2b.
\]

Thus

\[
\boxed{\quad k(t)=F(t)^{-T}k_0,
\quad b(t)=F(t)b_0\exp\!\left[-\nu\int_0^t|F(s)^{-T}k_0|^2ds\right].\quad}
\]

The velocity is recovered by `a=i k x b/|k|^2`. In particular `|a|=|b|/|k|` for real polarization up to the harmless complex phase, since `b·k=0`.

These identities include wavevector deformation, vorticity stretching and enhanced viscous damping. They are not a frozen-wavevector approximation. They rely on the stated symmetric base; a vortical affine base brings additional terms.

## 2. Exact finite-stage gain ceiling under isotropic transverse compression

Take `S=diag(-a/2,-a/2,a)` with constant `a>0`, `b_0` axial, and `k_0` transverse. Write

\[
X=e^{at}\ge1,
\qquad \eta=\nu|k_0|^2/a.
\]

Then

\[
|k|^2=|k_0|^2X,
\quad \frac{|b(t)|}{|b_0|}=X e^{-\eta(X-1)},
\quad \frac{|a(t)|}{|a_0|}=X^{1/2}e^{-\eta(X-1)}.
\]

Consequently the exact maximum vorticity gain over all `t>=0` is

\[
G_\omega(\eta)=
\begin{cases}
\eta^{-1}e^{-1+\eta},&0<\eta<1,\\
1,&\eta\ge1,
\end{cases}
\]

and the maximum velocity gain is

\[
G_u(\eta)=
\begin{cases}
(2\eta)^{-1/2}e^{-1/2+\eta},&0<\eta<1/2,\\
1,&\eta\ge1/2.
\end{cases}
\]

But the Reynolds number measured at inverse carrier frequency `h=1/|k|` satisfies

\[
\boxed{\quad\mathrm{Re}_h(t)
=\frac{|a(t)|}{\nu|k(t)|}
=\mathrm{Re}_h(0)e^{-\eta(X-1)}.\quad}
\]

It never increases. The exact result distinguishes substantial finite vorticity gain from the increasing-Reynolds restart required by `Re_m ~ rho_m^(-alpha)`.

The same ceiling for vorticity follows with `a` replaced by an upper bound `a_max` for an arbitrary nonnegative time-dependent axial strain: write `A=integral a`, `X=e^A`, and use

`integral_0^t e^{A(s)} ds >= (e^{A(t)}-1)/a_max`.

More basically, `Re_h(t)=Re_h(0) exp[-nu |k_0|^2 integral e^A]` remains exact. Unboundedly prescribing the affine strain has not generated that strain from a finite-energy fluid.

## 3. Exact stretched Gaussian vortex and a wider nonaxisymmetric bound

For arbitrary smooth axial strain `a(t)`, consider

\[
u_r=-\tfrac12a(t)r,\qquad u_z=a(t)z,\qquad
u_\theta=\frac{\mathcal C}{2\pi r}(1-e^{-r^2/s(t)}),
\quad
\omega_z=\frac{\mathcal C}{\pi s(t)}e^{-r^2/s(t)}.
\]

This is an exact unforced NS velocity with a suitable pressure if

\[
\boxed{\quad s'=-a(t)s+4\nu.\quad}
\]

The total cross-sectional circulation is exactly `C`, and

\[
s(t)=e^{-A(t)}\left[s(0)+4\nu\int_0^te^{A(\tau)}d\tau\right],
\qquad A(t)=\int_0^ta(\tau)d\tau.
\]

For constant `a`, the width tends to `4 nu/a` and vorticity saturates. For all such strain histories, `||u_theta||_infty sqrt(s)` is the same constant multiple of `|C|`; it follows directly by setting `r=sqrt(s) rho`. Thus its transverse Reynolds number is fixed by `|C|/nu` even during compression. This is not the same assertion as constancy of circulation around every finite viscous material loop: viscous loops can lose or gain circulation in general. Here the exact conserved quantity is the whole cross-sectional integral.

The obstruction is wider than axisymmetry. Consider the exact z-independent axial-vorticity class

\[
u=(B(t)x_\perp+v(x_\perp,t),\,a(t)z),
\quad B=B^T,\quad \operatorname{tr}B=-a,
\quad \nabla_\perp\cdot v=0,
\quad v=K_{2D}*\omega.
\]

The vorticity is `(0,0,omega)`; its cross-section need not be circular or symmetric. The scalar equation is exactly

\[
\partial_t\omega+(Bx_\perp+v)\cdot\nabla\omega
=a\omega+\nu\Delta_\perp\omega,
\]

or, equivalently,

\[
\partial_t\omega+\nabla_\perp\cdot[(Bx_\perp+v)\omega]
=\nu\Delta_\perp\omega.
\]

For smooth sufficiently decaying solutions, `M(t)=||omega(t)||_1 <= M(0)` by the scalar Kato inequality; signed circulation is conserved, and for single-sign vorticity its magnitude equals `M`. A two-dimensional Biot–Savart estimate gives an explicit uniform ceiling. With `Omega=||omega||_infty`, split the kernel at radius `R`:

\[
\|v\|_\infty\le \Omega R+\frac{M}{2\pi R}
\le\sqrt{\frac{2M\Omega}{\pi}}
\quad\text{at }R=\sqrt{M/(2\pi\Omega)}.
\]

Define the vorticity-area length `ell_omega=sqrt(M/Omega)` when `M,Omega>0`. Then

\[
\boxed{\quad
\frac{\|v\|_\infty\ell_\omega}{\nu}
\le\sqrt{2/\pi}\frac{M(0)}{\nu}.
\quad}
\]

This is a proved bound on the induced vortex velocity at the explicitly defined vorticity-area scale. It is not a bound on the externally supplied affine velocity, on every possible definition of a highly anisotropic core size, or on general three-dimensional flows. Merely breaking circular cross-sectional symmetry in this stretched-vortex class therefore does not yield unbounded growth of this core Reynolds number.

## 4. Anisotropic Kelvin-wave escape consumes localization

Let

`S=diag(-a1,-a2,a1+a2)`, with constants `a1>a2>0`,

and choose `k_0=k0 e_y`, `b_0=b0 e_z`, so wave velocity is in the x direction. Write

\[
A_i=a_it,
\quad D=\frac{\nu k_0^2}{2a_2}(e^{2a_2t}-1).
\]

The exact global mode has

\[
|k|=k_0e^{A_2},\qquad
U=U_0e^{A_1-D},\qquad
\mathrm{Re}_h=\mathrm{Re}_{h,0}e^{A_1-A_2-D}.
\]

Thus a claim that **every** anisotropic Kelvin mode has nonincreasing wavelength Reynolds number would be false. There can be finite gain here; in the inviscid equation the gain is arbitrarily large with time.

However, a material envelope of initial x-width `Lx0` compresses to `Lx=Lx0 e^{-A1}`. Its scale separation from the carrier is

\[
\mathcal M(t)=|k(t)|L_x(t)
=k_0L_{x0}e^{-(A_1-A_2)}.
\]

Consequently

\[
\boxed{\quad
\mathrm{Re}_h(t)\mathcal M(t)
=\frac{U_0L_{x0}}{\nu}e^{-D}.
\quad}
\]

If the localized-wave construction needs `M(t)>=M_min>>1`, then its wave Reynolds number is bounded by the initial transverse envelope Reynolds number divided by `M_min`. Its amplification factor cannot exceed `M(0)/M_min`, even before viscous loss.

The localization condition has a direct divergence-free meaning. Multiplying `v_x=U F(k y)` by an x-envelope creates divergence of size `U/Lx`; a usual phase-antiderivative correction in the y-component has size `U/(k Lx)`. Such a correction ceases to be small precisely as `M=kLx` approaches one. This is a conditional ceiling for retaining the small-correction localized Kelvin-wave ansatz, not a theorem that all other localization methods fail. Once this condition fails, the previously neglected velocity component/geometry has become principal and must be analyzed anew. The envelope's own diffusion scale `nu/Lx^2`, absent from the exact infinite plane wave, is a further cost.

This identifies what a restart would have to regenerate: scale separation, transverse envelope Reynolds number, and the appropriate strain orientation. Choosing a new larger carrier by fiat restores separation but lowers its initial carrier Reynolds number; this trade cannot be omitted from a gain calculation.

## 5. Finite energy, circulation and the unresolved restart

All exact global examples above have a spatially unbounded affine background; Kelvin waves additionally have infinite spatial extent, and straight vortices are z-independent. None meets Clay's finite-energy whole-space requirements. Localizing the affine strain requires a divergence-free extension and changes the pressure and evolution. Vortex lines must also close or continue; an axially stretched material segment expands in length. `det F=1` prevents the same positive-volume material core from shrinking in all three directions.

The earlier isotropic nested-core budget `U_m L_m/nu -> infinity` therefore has no realization here by repeatedly stretching the same circulation-carrying packet. A viable escape must provide genuinely three-dimensional regeneration, for example:

1. A newly selected, nonmaterial receiving core intercepts/collects additional aligned vorticity flux from several strands, increasing its effective circulation while maintaining the pressure/strain geometry.
2. Folding and reorientation restore an effective short axial size and sufficient transverse localization margin, with the return strands and their cancellation explicitly included.
3. Viscous transport and any smooth force are included in the circulation budget, rather than treating Kelvin circulation as exactly conserved for arbitrary NS material loops.

These are proposed obligations, not demonstrated mechanisms. The precise finite-stage target is a return map on **localized divergence-free data** that increases the core circulation/induced-velocity scale and restores `M>=M_min`, after accounting for the finite reservoir and every incoming/returning strand. The exact affine amplification formulas above provide a baseline and a falsifier: if a proposed restart is only a change of scale on the same packet, it cannot meet that target.

## Verification and sources

`python3 work/pass2/verify_affine_stretch_gate.py` checks the full Gaussian stretched-vorticity residual, its circulation normalization, Kelvin-mode amplitude/vorticity/viscosity identities, isotropic gain critical points, and anisotropic Reynolds–localization product. This is exact algebra for the models stated here, not verification of a finite-energy NS construction.

Background primary sources checked 2026-09-08:

- Gallay–Maekawa, *Three-dimensional stability of Burgers vortices*: https://arxiv.org/abs/1002.2489 . Establishes the classical Burgers-vortex setting and circulation Reynolds-number convention; no stability result is needed for the calculations above.
- Fabijonas–Holm, *Multi-frequency Craik–Criminale solutions of the Navier–Stokes equations*: https://arxiv.org/abs/nlin/0304049 . Context for exact Kelvin-wave constructions; the restricted symmetric-affine formulas here were derived directly.
- Gallay–Maekawa, *Existence and stability of viscous vortices*: https://www-fourier.univ-grenoble-alpes.fr/~gallay/Handbook.pdf , §1 and §4 for straight-vortex circulation and the Biot–Savart representation. The explicit split-kernel constant and localization ceiling above are our calculations.

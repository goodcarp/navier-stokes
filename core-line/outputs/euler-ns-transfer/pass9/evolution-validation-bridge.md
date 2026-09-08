# A posteriori bridge for one actual leading-envelope NS evolution

**Conditional validation theorem, not a validated evolution.** This note
states two sufficient, computable alternatives that would transfer a
full-field approximation on `[0,T]` to the actual pass8 solution. The
preferred minimal route validates retained endpoint gains in H4; the
optional H7 route also reaches all-time derivative and curvature tests.
No such
approximation, useful duration, or error radius has yet been certified.
The size of the constants below does not establish useful persistence.

Fix one leading-envelope datum: either a pass8 pressure-neutral member
with its exact retuned `A_*`, or the selected explicit datum
`A_*=1020, lambda_*=1/4`. The target is the ordinary unforced NS solution
from that one datum, at `nu=1/1000`. Its amplitude stays fixed during
evolution. The H4 endpoint route and the error estimates below apply to
both choices. The explicit datum has `beta'(0)=d>0`; the optional
curvature integration below retains this additional linear term.

## 0. Preferred minimal route: H4 endpoint validation

The target is a retained finite gain and usable inherited endpoint. It
does **not** require the mean derivative, energy derivative, or core
ratio curvature to stay positive at every intermediate instant.
The following lower-order route avoids imposing those stronger gates.

Use the derivative-sum norm defined below, the full residual (1), and
`E4=||u−v||4`. If `||R||4<=delta4` and `||v||5<=V5`, then

\[
 \boxed{D^+E_4\le560V_5 E_4+270E_4^2+\delta_4.}
 \tag{0a}
\]

There are 35 multiindices of order at most four. Split the Leibniz
terms at `|beta|<=2` versus `|beta|>=3`, putting the lower derivative
factor in L-infinity using H2. This gives the sufficient linear
constant `93 sqrt(35)<560` and quadratic constant `45 sqrt(35)<270`.
Top-order transport cancellation and the viscous sign are the same as
in (2). Formula (3), with `a=560V5`, `4200` replaced by `270`, and
`e0,delta` replaced by their H4 counterparts, gives a certified radius
`rho4` as soon as its denominator stays positive. This also gives actual
H4 continuation through `T`. Initial-amplitude, reconstruction, and
slab-boundary errors must still be carried.

At a chosen endpoint receiver `(r_T,z_T)`, the actual mean satisfies

\[
 G(u(T),r_T,z_T)\ge G(v(T),r_T,z_T)-r_T\rho_4(T).
 \tag{0b}
\]

Let `M0_up` be a certified upper bound for the actual initial mean
maximum, and let `eta_G>0` be the desired gain. The sufficient condition

\[
 G(v(T),r_T,z_T)-r_T\rho_4(T)
       \ge(1+\eta_G)M_{0,\rm up}
 \tag{0c}
\]

proves at least that gain in the actual mean maximum. No intermediate
torque or acceleration sign is required. The approximation's receiver
coordinate must be the one at which its value is enclosed; uncertainty
in that coordinate also incurs the corresponding spatial error.

For actual nonaxisymmetric energy, the contraction of `I−P0` in L2
gives the particularly simple endpoint lower bound

\[
 K(u(T))\ge\tfrac12
    \left(\sqrt{2K(v(T))}-\rho_4(T)\right)_+^2.
 \tag{0d}
\]

Requiring the right side to exceed `(1+eta_K)K0_up` proves a prescribed
energy gain. `K0_up` must enclose the same exact initial datum's energy;
a crude bound may be too wide for a small useful gain. This test is
valid whether the fluctuation energy rises monotonically or not.

H4 also controls the core gradient. At the actual symmetric core,
`|b−b_v|,|Omega−Omega_v|<=rho4`. Therefore

\[
 \Omega_v(T)-\rho_4(T)>0,
\]
\[
 b_v(T)-(1+\eta_\beta)\Omega_v(T)
              -(2+\eta_\beta)\rho_4(T)\ge0
 \tag{0e}
\]

imply `beta(T)>=1+eta_beta`, from the actual initial `beta(0)=1`.
Neither neutrality at time `T` nor positive `beta″` throughout the
interval is assumed. A separate lower bound on `Omega_v−rho4` can
require retained rotational amplitude, rather than only a ratio gain.

If a full velocity/Reynolds gain is required, use a specified
nonnegative receiving weight `chi_L` of mass `M_L>0` and define
`V_L(a)=||sqrt(chi_L)a||2/sqrt(M_L)`. Then

\[
 V_L(u(T))\ge V_L(v(T))-
           \sqrt{\|\chi_L\|_\infty/M_L}\,\rho_4(T).
 \tag{0f}
\]

Multiplying this lower bound by `L/nu` gives a retained Reynolds-number
test for that receiver. This is a full-velocity observable; mean angular
momentum or fluctuation-energy gain need not be substituted for it.

Finally, H4 controls C2 errors and a successor class formulated in the
full H4 field. The analogue of (15) is

\[
 \epsilon_{{\rm end},4}\le
 U^{-1}L^{-3/2}\max(1,L^4)\rho_4(T),
 \tag{0g}
\]

plus actual coordinate/scale parameter errors. A distance to the boundary
of the claimed H4 successor class larger than (0g) proves membership for
the **same inherited actual endpoint**. This is enough only if the
successor-stage lemma is in fact stable in that named class; an H7
requirement cannot be silently replaced by H4. The trajectory/deformation
bounds (14) also work with `rho4`. No old field is removed and no
envelope is rebuilt at the endpoint.

Conditions (0a)–(0g) are conditional tests, not evaluated claims about
the present candidate. They can certify useful retained gains even if
some instantaneous derivatives change sign. A positive denominator for
a tiny time still does not establish the prescribed margins.

## 1. Optional stronger route: a norm reaching the derivative observables

Use the archive's derivative-sum norm

\[
 \|a\|_s^2=\sum_{|\alpha|\le s}\sum_{j=1}^3
                  \|\partial^\alpha a_j\|_{L^2(\mathbb R^3)}^2.
\]

Let `v(t,x)` be a sufficiently smooth, exactly divergence-free approximate
**whole-space velocity**, with its actual mean, fluctuation, and older
field all included. Define the Leray-projected vector residual

\[
 \mathcal R=v_t-\nu\Delta v+\mathbb P(v\cdot\nabla v),
 \quad N(a)=\nu\Delta a-\mathbb P(a\cdot\nabla a),
 \quad E=\|u-v\|_7.
 \tag{1}
\]

Certify an initial error `E(0)<=e0`, a residual bound
`||R(t)||7<=delta(t)`, and `||v(t)||8<=V8(t)` on the **whole interval**.
An exact reduced solution is eligible only after it is embedded in a full
velocity and its full residual is bounded. Discarded vertical motion,
mean stress, localization, pressure projection, and temporal joins belong
in (1). A small scalar ODE residual is not a substitute.

The actual initial error includes uncertainty in the exact pressure
retuning. For example, replacing `A_*` by `A_num` contributes
`|A_num−A_*| ||v_outer||7` to the triangle-inequality upper bound for
initial error, together with all representation errors. The H4 route
uses the analogous term with `||v_outer||4`. The coarse interval
`900<A_*<1020` does not establish a small `e0`.

`H6` suffices for the scalar instantaneous pressure derivative used at
insertion. We use **H7** because the positive-time core-ratio curvature
also involves viscous fifth spatial derivatives evaluated at the core.
An energy-norm error alone controls neither that quantity nor pointwise
mean torque.

## 2. Explicit error majorant

On the actual solution's existence interval, the following elementary
estimate holds in the upper-Dini-derivative sense:

\[
 \boxed{D^+E\le a(t)E+4200E^2+\delta(t),
       \qquad a(t)=8500V_8(t).}
 \tag{2}
\]

To verify the constants, apply all derivatives of order at most seven
to the error equation, retain the favorable viscous energy term, and
cancel each divergence-free top-order transport term. There are
`binomial(10,3)=120` multiindices. Each vector product costs at most a
factor three. In a transport commutator, split at `|beta|<=3` versus
`|beta|>=4`; the lower derivative factor has enough remaining derivatives
for the established scalar `H2 -> L-infinity` bound with constant one.
The commutator weight is at most `2^7−1=127`. The nontransport term
`e dot grad v` has weight at most `128` and uses `||v||8`.
Thus the linear coefficient is at most `765 sqrt(120)<8500`, and the
quadratic coefficient at most `381 sqrt(120)<4200`. The Leray projection
commutes with derivatives and can be removed when paired with a
divergence-free differentiated error.

The large linear coefficient is a safe baseline. It may be replaced by
any **verified** upper bound `a(t)` for the logarithmic growth of the full
linearized viscous NS operator in this norm. Such a replacement must
include unresolved spatial/frequency tails; an eigenvalue calculation
on a finite matrix alone does not prove it.

One entirely explicit supersolution is obtained by setting

\[
 A(t)=\int_0^t a(s)\,ds,\quad
 h(t)=e_0+\int_0^t e^{-A(s)}\delta(s)\,ds,
\]
\[
 D(t)=1-4200\int_0^t e^{A(s)}h(s)\,ds,\qquad
 \boxed{\rho(t)=e^{A(t)}h(t)/D(t).}
 \tag{3}
\]

If the integrals are rigorously enclosed and `D(t)>0` through `T`, then
`E(t)<=rho(t)`. Differentiation verifies that (3) is a supersolution;
the extra term has the favorable sign because `h′>=0` and `0<D<=1`.
Equivalently, a validated scalar Riccati integration of (2) can give
a tighter bound. A finite error bound and the standard Sobolev continuation
criterion extend the actual solution through `T`.

Time subdivision is allowed: the next slab inherits its previous
certified error radius. It is not restarted at zero. A changed numerical
representation at a slab boundary also incurs its measured jump in H7.

These baseline estimates may give an extremely small interval for the
steep, large-amplitude pass8 field. A finite Sobolev norm or a positive
formal denominator at a microscopic time is not a useful-stage result.

## 3. Pointwise mean torque and actual receiving gain

For a field `a`, put `U=P0 a`, `w=a−U`, and, in cylindrical components,

\[
 G=rU_\theta,\quad R=\overline{w_rw_\theta},\quad
 Z=\overline{w_zw_\theta},\quad
 \mathcal T=-2R-r\partial_rR-r\partial_zZ.
 \tag{4}
\]

The approximation must evolve these quantities from its own full field;
they are not prescribed stresses. On a fixed exterior region
`r_-<=r<=r_+`, where `r_->0`, set `Y=||v||7` and `X=Y+rho`.
Conservative explicit consequences of H7 embedding are

\[
 |\mathcal T(u)-\mathcal T(v)|
 \le (8+32r_+)(X+Y)\rho=:\epsilon_T.
 \tag{5}
\]

For example, meridional component derivatives of a field error are
bounded by `2rho`; those of its fluctuation by `4rho`. Product
telescoping gives `|delta R|<=4(X+Y)rho` and
`|delta R_r|,|delta Z_z|<=16(X+Y)rho`, proving (5). Computing narrower
local component intervals from the approximation can improve this bound.

Positive torque alone is not positive `G_t` after insertion. The exact
instantaneous NS functional is

\[
 H_G(a)=-U_rG_r-U_zG_z+
       \nu(G_{rr}-G_r/r+G_{zz})+\mathcal T(a).
 \tag{6}
\]

Let tildes denote quantities of `v`. A directly computable local error
bound for (6) is

\[
\begin{aligned}
 \epsilon_G={}&\rho\big[|\widetilde G_r|+|\widetilde G_z|
 +( |\widetilde U_r|+\rho)(1+2r_+)
 +( |\widetilde U_z|+\rho)2r_+\big]\\
 &+\nu(8r_++2+r_-^{-1})\rho+\epsilon_T.
\end{aligned} \tag{7}
\]

Here `Delta_*G=r(U_theta,rr+U_theta,zz)+U_theta,r−U_theta/r`;
the diffusion bound in (7) retains every term. The first line explicitly
pays the mean radial and axial evolution errors. If, on the chosen
fixed receiving circle and all `0<=t<=T`,

\[
 H_G(v)-\epsilon_G\ge g_*>0,
 \tag{8}
\]

then `G(u,T)>=G(u,0)+g_*T`. Starting at the initial global mean maximum
gives at least this increment in the supremum. For a specified moving
receiver, add its actual prescribed velocity dotted with `grad G` to
(6), and charge the corresponding gradient error; a moving global
maximum is not silently substituted for this receiver.

## 4. Integrated fluctuation energy is easier, but remains coupled

Define

\[
 K(a)=\tfrac12\|(I-P_0)a\|_2^2,\quad
 J(a)=-\int w_iw_j\partial_jU_i\,dx-\nu\|\nabla w\|_2^2.
\]

For the actual solution, `K′=J(u)`. Rotational averaging is an orthogonal
projection in L2 and in the full first-gradient norm. Product
telescoping and the conservative embedding `||grad a||infinity<=3||a||7`
give

\[
 |K(u)-K(v)|\le\tfrac12(X+Y)\rho=: \epsilon_K,
\]
\[
 |J(u)-J(v)|\le
 [3(X^2+XY+Y^2)+\nu(X+Y)]\rho=: \epsilon_J.
 \tag{9}
\]

Thus the verified condition

\[
 J(v)-\epsilon_J\ge\kappa[K(v)+\epsilon_K],\qquad \kappa>0,
 \tag{10}
\]

on the whole interval implies `K(u,T)>=exp(kappa T)K(u,0)`.
Unlike an initial fractional derivative, this is a finite-time gain.
A desired factor `1+eta` requires the actual certified interval to satisfy
`kappa T>=log(1+eta)`. It cannot be inferred from finiteness of constants.

## 5. Full off-neutral core ratio, including viscous curvature

The actual odd/C4 symmetries persist by uniqueness, keeping the origin
stationary and the central gradient in the form
`diag(−b,−b,2b)+Omega J`. It is convenient to enforce these symmetries
exactly on `v` as well and then recompute its full residual. Put

\[
 h=p_{zz}(0),\quad F=h_t=2\Pi(u,N(u)),\quad
 d_b=\tfrac12\Delta\partial_z u_z(0),\quad
 d_\Omega=\tfrac12\Delta(\partial_xu_y-\partial_yu_x)(0),
\]

and let `d_b1,d_Omega1` be the same two linear functionals applied to
`N(u)`. For the actual unforced solution, **at arbitrary positive time**,

\[
 b_1=-2b^2-h/2+\nu d_b,\qquad
 \Omega_1=2b\Omega+\nu d_\Omega,
\]
\[
 b_2=-4bb_1-F/2+\nu d_{b1},\qquad
 \Omega_2=2b_1\Omega+2b\Omega_1+\nu d_{\Omega1},
\]
\[
 \boxed{\beta''=\frac{b_2}{\Omega}
 -\frac{b\Omega_2+2b_1\Omega_1}{\Omega^2}
 +\frac{2b\Omega_1^2}{\Omega^3}.}
 \tag{11}
\]

The initial shortcut `beta″=−(F+32)/2` uses initial affine-core
cancellations and neutrality. It is not the equation to propagate.

Here is an explicit interval recipe for (11). Enclose its ingredients
for `v`, then enlarge those intervals by the following error radii:

| Quantity | Sufficient error radius |
|---|---:|
| `b, Omega` | `rho` |
| `h=Pi(u,u)` | `128(X+Y)rho` |
| `F=2Pi(u,N(u))` | `256rho[3nu(X+Y)+512(X²+XY+Y²)]` |
| `d_b,d_Omega` | `3rho` |
| `d_b1,d_Omega1` | `3 L_N rho`, with `L_N=3nu+1024(X+Y)` |

Indeed `||N(u)−N(v)||5<=L_N rho`, using
`||B(a,b)||5<=1024||a||7||b||7` and `||Delta a||5<=3||a||7`.
The pressure bounds are the previously audited whole-space bilinear
estimates; their multiplier definition retains the pressure contact term.
Interval evaluation of (11), after certifying `Omega>=omega_min>0`,
therefore encloses the actual ratio curvature. The high derivative core
terms need not be computed through a differentiated numerical pressure:

\[
 \Delta N(u)=\nu\Delta^2u-\Delta(u\cdot\nabla u)
                     +\nabla\operatorname{tr}(\nabla u\nabla u).
 \tag{12}
\]

This exact local cancellation includes their fifth-order viscous terms.
The full nonlocal `h,F` still have to be enclosed separately.
If the lower interval endpoint of (11) exceeds `b_*>0` throughout
`[0,T]`, a neutral pass8 datum with `beta(0)=1,beta′(0)=0` satisfies
`beta(T)>=1+b_*T²/2`. For the selected explicit amplitude, `beta(0)=1`
and `beta′(0)=d>0` instead give `beta(T)>=1+dT+b_*T²/2`.
This condition includes changes in the core's spatial curvature and
pressure, rather than assuming it remains affine.

## 6. The approximate time derivative is not the NS functional

When evaluating (6), (9), or (11), apply `N` to the approximate field.
Do not interchange this with differentiating a fitted numerical time
curve without paying its residual. In particular, if
`||R||7<=delta`, then

\[
 |\partial_tG(v)-H_G(v)|\le r_+\delta,
\]
\[
 |\partial_tK(v)-J(v)|\le Y\delta,
\quad
 |\partial_t h(v)-F(v)|\le256Y\delta.
 \tag{13}
\]

Naively differentiating a numerical `beta(t)` twice can additionally
introduce time derivatives of its residual. Formula (11) avoids that
unstated requirement by evaluating the actual NS acceleration
functionals at each approximate field.

## 7. Exterior trajectories, deformations, and the inherited endpoint

To compare actual and approximate material trajectories, let
`X′=u(t,X)`, `Xtilde′=v(t,Xtilde)`, and suppose
`||grad v||infinity<=L`, `||grad²v||infinity<=M2`, either globally or
on a certified tube containing both trajectories and, at each time,
the entire straight segment joining their positions. This segment
condition justifies the mean-value estimates. Their distance and
flow-Jacobian error then satisfy the usable differential bounds

\[
 d_X'\le Ld_X+\rho,
\]
\[
 d_D'\le(L+3\rho)d_D+
       (3\rho+M_2d_X)\|D\widetilde X\|.
 \tag{14}
\]

Enclose the tube itself using these bounds, including the connecting
segments; do not assume the required containment.
Centers and widths can be checked through (14), but a tracked material
image of the initial packet is not the support of a viscous evolving
component. Generated tails and the older field remain in (1). Mean
maximizing circles also need not follow material trajectories.

At time `T`, the same actual solution satisfies
`||u(T)−v(T)||7<=rho(T)`. Suppose a successor class is defined in the
full normalized field `a(x_c+Ly)/U`, with fixed positive `L,U` and the
current coordinate axes. Its inherited norm error is bounded by

\[
 \epsilon_{\rm end}\le
 U^{-1}L^{-3/2}\max(1,L^7)\rho(T).
 \tag{15}
\]

Translation uncertainty adds a separately bounded term from the
approximation's next derivative; uncertain scaling or rotation must
likewise be charged. If the normalized approximation lies strictly
inside the successor class by more than (15) and these parameter errors,
the *actual inherited endpoint* lies in that class. Every old or exterior
field is retained in this whole-space comparison. Selecting new fitting
parameters is allowed; replacing the endpoint by the fitted profile is
not. None of (8), (10), or (11) alone proves this terminal membership.

## What would constitute a useful certificate

A minimal successful finite validation would provide the full
reconstructed field and space-time residual enclosures, a nonblowing-up
H4 error majorant, the required prescribed endpoint gains, and an actual
terminal-class margin surviving (0g). The optional H7 route can additionally
certify positive margins in (8), (10), and the full (11) throughout the
interval and a terminal margin in (15). Those stronger monotonicity and
curvature conditions are not prerequisites for the minimal endpoint route.
Finite-box or spectral computations need their domain, projection,
interpolation, roundoff, time, and tail errors included. Computing
another initial jet or observing floating-point convergence supplies
none of those missing interval obligations.

Smooth forcing remains allowed by the broader target. For a chosen
physical force, subtract its Leray projection in (1); the error estimate
remains valid when the exact and approximate problems use that same
force. The observable functionals and core identities must then include
its actual contributions, including required time derivatives. The
unforced formulas above cannot be reused unchanged. A proposal that
sets a full residual equal to a physical force is a different, explicitly
forced construction and must establish that force's required smoothness
and assembly bounds. No artificial ban on forcing is imposed here.

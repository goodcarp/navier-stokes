# One finite-amplitude three-dimensional family passes both initial tests

For the **same** compact outer chiral seed and exact retuning from pass 6, every
amplitude `3/4 <= lambda <= 1` now has favorable full central pressure feedback
and a growing positive azimuthal-mean angular-momentum maximum, at viscosity
`nu=1/1000`. All nonlinear pressure interactions are included. This is an
initial-data and local-time result. It does not prove a return or blowup.

## The actual datum

Keep the established profiles `S_Phi`, `W_psi`, and unit outer swirl `v`.
The cutoff is `chi=1` on `s<=1/4`, `chi=0` on `s>=1`, and on the transition,
with `t=(4s-1)/3`, it is
`chi(s)=1/[1+exp(-1/t+1/(1-t))]`.
Write cylindrical coordinates as `(r,theta,z)`. Set

\[
 e(r)=\chi((r-21/100)^2/(7/50)^2),\qquad
 q(z)=\chi((z-4)^2/(3/5)^2)+\chi((z+4)^2/(3/5)^2),
\]
\[
 w=\operatorname{curl}\{e(r)q(z)\cos(4\theta+20r)e_z\},\qquad
 C_w=-\Pi(w,w),\quad C_v=-\Pi(v,v),\quad H_0=204/35,
\]
\[
 \boxed{A_\lambda=\sqrt{(H_0-\lambda^2 C_w)/C_v},\qquad
 u_{0,\lambda}=S_\Phi+W_\psi+A_\lambda v+\lambda w.}
 \tag{1}
\]

The unchanged pass 6 bounds are `0<Cw<3/2`, `780<A_lambda<1020`, exact
`p_zz(0)=-8`, and initial `b=Omega=1`, `b'=Omega'=2`, `beta'=0` for
`beta=b/Omega`. Oddness and fourfold rotation symmetry preserve the required
central matrix. The seed avoids the affine core and the axis.

These are the outer packets, not the different small inner receiver seed.

## Full pressure-derivative certificate

Let `T_lambda=p_zz'(0)+32`. The exact pass 6 expansion gives

\[
 T_\lambda< -m_0+\lambda^2
 \left[\frac{23199}{1000}C_w+D_\nu+A_\lambda E\right],
 \qquad m_0=\frac{290649}{8750}.
 \tag{2}
\]

The new full-kernel and interval calculations prove

\[
 \boxed{D_\nu<4-\frac{96}{5}C_w,\qquad E<\frac1{60}.}
 \tag{3}
\]

The `D_nu` bound includes the entire meridional/core pressure response and
viscosity. The `E` bound includes the full nonlocal cross pressure of the
strong swirl and the seed; it is not set to zero by symmetry. Its proof uses
an exact horizontal Green trial with both radial tails, and bounds the
whole-space residual by an angular-mode energy estimate. The outward
2048-panel range computation gives an upper bound below `0.016616141`.
This number is an upper bound, not a two-sided enclosure or the value of E.

Since `A_lambda>0`, (3) implies `A_lambda E<17` whether or not E is positive.
Keeping the useful `Cw` term yields

\[
 T_\lambda< -m_0+\lambda^2
 \left[\frac{3999}{1000}C_w+21\right]
 <-m_0+\frac{3999}{1000}\frac32+21
 =-\frac{435297}{70000}.
 \tag{4}
\]

The second step uses the positivity of the bracket and `lambda^2<=1`.
The actual full pressure therefore satisfies, throughout the amplitude interval,

\[
 \boxed{p_{zz}'(0)+32<-\frac{435297}{70000},\qquad
 \beta''(0)>\frac{435297}{140000}>3.109.}
 \tag{5}
\]

The relation between these quantities uses the flat-core initial curvature
cancellations. It does not neglect viscosity at later times.

## The same datum grows the positive global mean maximum

Let `G_lambda(t,r,z)` be the azimuthal mean of `r u_theta`. Choose any radius
`r_*` maximizing `r^2 chi(r^2/(63/200)^2)`, and take `z=4` (or `z=-4`).
The pass 6 proof places these circles in the flat seed envelope and identifies
their value with the positive global initial mean maximum. It proves

\[
 \boxed{\partial_tG_\lambda(0,r_*,4)>\frac{626}{35}>0}
 \tag{6}
\]

for every amplitude in (1), including the initial viscous curvature. At
`lambda=1` the stronger lower bound is `27509/315`.
Equation (6) concerns the azimuthal mean. No corresponding claim about the
pointwise angle-dependent maximum is made.

## Actual simultaneous local-time consequences

The family (1) is compact in every fixed Sobolev norm and consists of smooth
compact solenoidal data. Classical local existence and smooth dependence give
a common positive existence time. The strict margins in (5)--(6), continuity
of the actual solution jets, and compactness of the amplitude interval give
a common `tau>0` on which both signs persist.

Put `B0=435297/140000` and `g0=626/35`. Shrink that common interval so that
`beta_lambda''(t)>=B0/2`, `Omega_lambda'(t)>=1`, and
`partial_t G_lambda(t,r_*,4)>=g0/2`. Then, for `0<t<=tau`,

\[
 \boxed{\beta_\lambda(t)\ge1+\frac{435297}{560000}t^2,
 \quad \Omega_\lambda(t)\ge1+t,
 \quad \sup_{r,z}G_\lambda(t,r,z)
 \ge\sup_{r,z}G_\lambda(0,r,z)+\frac{313}{35}t.}
 \tag{7}
\]

There is no need to differentiate a moving maximizer: evaluate the evolved
mean at the fixed initial circle. For `h=b'-2b^2`, the central equations and
initial flatness give `h(0)=0` and `h'(0)=beta''(0)>B0`. A further common
shortening gives `b'(t)-2b(t)^2>=B0 t/2>0`.

These statements refer to the inherited solution of (1), with no reset or
fresh force. They give no useful numerical value of tau, no prescribed
finite gain, and no control showing entry into a successor profile class.
One cannot integrate the short-time Riccati inequality to a singularity
without proving it persists. That is a remaining mathematical obligation.

The active full Navier--Stokes goal remains open. The next task is a useful
quantitative interval with retained gain and a controlled inherited endpoint,
followed by an actual return or another compatible continuation mechanism.

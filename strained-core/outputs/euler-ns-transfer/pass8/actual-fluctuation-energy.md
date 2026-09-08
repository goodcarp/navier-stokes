# The actual initial mean transfer drains the nonaxisymmetric seed

This is an exact identity for the actual unforced Navier–Stokes solution
starting from the pass7 family, at viscosity 1/1000. It is not a frozen-phase
model. For the nonaxisymmetric kinetic energy K, it gives

\[
 \boxed{\frac{K'(0)}{K(0)}
 <-\frac{2953517}{67640}<-43.66.}
 \tag{1}
\]

Only the initial derivative is bounded here. No exponential decay over a
prescribed interval follows from (1), and finite energy loss does not rule
out concentration, transfer to another receiver, or a later regeneration.

## Full mean/fluctuation energy identity

Let P0 be azimuthal averaging of the vector field by rotation-equivariant
components, U=P0u, v=u-U, and K=1/2 integral |v|². P0 is an orthogonal
projection in physical L² and commutes with derivatives in the rotation-
equivariant sense, the Leray projection and the Laplacian. The usual whole-
space integration by parts gives

\[
 K'=-\int v_i v_j\partial_j U_i\,dx-\nu\int|\nabla v|^2\,dx.
 \tag{2}
\]

The full fluctuation pressure and cubic self-transport cancel in this
integrated identity. They do not vanish from the local covariance equation.
Smooth compact initial data and local finite-energy evolution justify the
integrations. At insertion v=lambda w and U=M+A_lambda v_out.

On the seed support, M=c(r e_r-2z e_z), c=7/5, and write
v_out=r C(r)q_v(z)e_theta, C(r)=chi(r²/a²), a=63/200.
The seed is horizontal with amplitude e(r)q_s(z). Its angular covariance is

\[
 \overline{w_r w_\theta}=-\frac{mk}{2r}e^2q_s^2,\qquad m=4,\ k=20.
\]

The actual strain contraction is
c|w|²+r A_lambda C'(r)q_v(z)w_r w_theta.
Consequently

\[
 \boxed{K'(0)=-2cK(0)-\nu\lambda^2\|\nabla w\|_2^2
 +\pi\lambda^2mk A_\lambda
       \left(\int r C'e^2\,dr\right)
       \left(\int q_s^2q_v\,dz\right).}
 \tag{3}
\]

The last term is strictly negative. The decreasing swirl profile and the
chosen useful chirality transfer energy from the fluctuation into the mean.
The horizontal pump strain also removes fluctuation energy initially.

## A rational fractional drain bound

The radial seed envelope is one on [7/50,7/25]. On
[63/400,7/25], C decreases from one to a value below one half:
(7/25)²/a²>5/8 and chi(5/8)=1/2. Therefore

\[
 -\int rC'e^2\,dr>\frac{63}{800}.
\]

The common flat axial plateau of q_s and q_v has total length 9/10, so
integral q_s²q_v dz >=9/10. Also A_lambda>780. The unchanged pass6
support/variation bounds give

\[
 K(0)=\frac{\pi\lambda^2}{2}R_0Z_0,\quad
 R_0\le6764/75,\quad Z_0\le12/5.
\]

Here R0 includes radial-envelope, radial-phase and angular derivatives;
it is not just one velocity component. Dividing the lower magnitude of
the swirl work in (3) by this upper energy bound yields

\[
 -\frac{\text{swirl work}}{K(0)}
 >\frac{mk\,780(63/800)(9/10)}
         {(1/2)(6764/75)(12/5)}
 =\frac{552825}{13528}.
\]

Adding 2c=14/5 and dropping only the favorable viscous dissipation proves
(1), uniformly for 3/4<=lambda<=1. Lambda cancels in this ratio.

This identifies the supplying reservoir for the initial mean gain. It
does not justify freezing or replenishing that reservoir during a return.
The next actual stage must account for its depletion and for whatever
coupled mechanism would regenerate a usable seed.

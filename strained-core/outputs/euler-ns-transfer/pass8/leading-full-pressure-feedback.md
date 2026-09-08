# The replacement seed passes all three actual initial requirements

For the **same new leading-envelope datum**, mean angular momentum, fluctuation
energy, and central pressure feedback all have favorable initial signs with
viscosity included. This is a candidate redesign, not a regenerated endpoint
of the previous solution. A useful quantitative interval and a return remain open.

## The new compact family

Keep the pass7 core and radial pump M=S_Phi+W_psi and unit outer swirl v.
Replace only the seed profile and its amplitude range by

\[
 e_L(r)=\chi(((r-7/50)/(1/10))^2),\quad
 w_L=\operatorname{curl}\{e_L(r)q_s(z)\cos(4\theta-20r)e_z\},
 \qquad 1/8\le\lambda\le1/4.
\]

The axial q_s is unchanged, with centers plus/minus four and width 3/5.
Set C_L=-Pi(w_L,w_L)>0 and use the named exact amplitude

\[
 A_\lambda=\sqrt{\frac{204/35-\lambda^2 C_L}{C_v}},\qquad
 u_{0,\lambda}=M+A_\lambda v+\lambda w_L.
 \tag{1}
\]

The [leading-envelope certificate](leading-envelope-alternative.md) proves
0<C_L<3/2, 900<A_lambda<1020, exact p_zz(0)=-8, and

\[
 \boxed{G_t(0,r_*,4)>\frac{8871}{800},\qquad
       \frac{K'(0)}{K(0)}>\frac{12853}{2205}.}
 \tag{2}
\]

K is the actual nonaxisymmetric kinetic energy, and G is the azimuthal
mean of r u_theta. The receiver is an initial positive global mean maximizer.
The seed has positive covariance, while the spatial derivative of that
covariance produces positive receiving torque. This permits extraction of
energy from the decreasing mean swirl and growth of a local mean maximum
at the same time. Total kinetic energy still dissipates.

## Complete pressure interaction bound

The exact degree and rotation-selection calculation is unchanged for this
new m=4 seed. With T_lambda=p_zz'(0)+32,

\[
 T_\lambda<-m_0+\lambda^2
 \left[\frac{23199}{1000}C_L+D_L+A_\lambda E_L\right],
 \qquad m_0=290649/8750.
 \tag{3}
\]

Its support r in [1/25,6/25], |z| in [17/5,23/5] lies inside the same
affine pump plateau and a smaller axial cone. The seed is still horizontal,
so the complete horizontal D kernel from pass7 gives

\[
 D_L\le-\frac{96}{5}C_L+\nu d_L.
\]

The actual leading-envelope interval calculation bounds R1<200000 and
the analytic energy budget gives R0<=392/5. Including both axial derivative
terms and the full tensor Q,

\[
 \nu d_L\le
 \frac6{1000}(5/17)^5
 [200000(12/5)+(392/5)(160/3)]
 =\frac{9078400}{1419857}<7.
 \tag{4}
\]

No contribution from the core or pressure is omitted in this D bound.

For E_L retain the exact decomposition E_loc+E_press. Since k=-20, the
local term is now **positive** and cannot use the favorable sign of the old
trailing seed. Its absolute value satisfies

\[
 |E_{\rm loc}|
 \le135|mk|\frac{(6/25)^3-(1/25)^3}{3}
       \frac{12/5}{(17/5)^7}<0.022635.
 \tag{5}
\]

This uses the actual kernel bound
|partial_r Q_theta|<=90r/(4 pi |z|^7), V0<=r, and compact supports.

The full nonlocal mixed-pressure term is bounded using the exact gradient
projection identity

\[
 E_{\rm press}
 =-4\langle\mathbb P_{\rm grad}(Qw_L),(w_L\cdot\nabla)v\rangle,\qquad
 |E_{\rm press}|\le\frac4m
       \|r\operatorname{div}(Qw_L)\|_2\|(w_L\cdot\nabla)v\|_2.
 \tag{6}
\]

The angular energy estimate in (6) is a whole-space bound. It neither drops
the pressure nor replaces it with a numerical finite-box solution.
On this support V/r<=1 and
|V_r|<=max(1,8(r/a)^2)<5. The seed is horizontal, so V_z does not
enter (w_L dot grad)v, and

\[
 \|(w_L\cdot\nabla)v\|_2\le
 5\|w_L\|_2\le5\sqrt{\pi(392/5)(12/5)}.
\]

The same separated Q-kernel bound as in pass7, evaluated for the new
envelope by one-dimensional outward range integration, gives

\[
 \|r\operatorname{div}(Qw_L)\|_2<0.002311336.
\]

Combining it with (5)--(6) bounds the full coefficient by

\[
 \boxed{|E_L|<0.303612507<1/3.}
 \tag{7}
\]

Displayed decimals are outward summaries; the executable certificate stores
and compares exact dyadic endpoints.

## Three simultaneous initial gains

From (3)--(7), positivity of A_lambda, and lambda²<=1/16,

\[
 T_\lambda<-m_0+
 \frac1{16}\left[\frac{3999}{1000}\frac32+7+\frac{1020}{3}\right]
 =-\frac{12493177}{1120000}<-11.15.
\]

The same odd/fourfold symmetry and initial affine core give

\[
 \boxed{\beta''(0)>
 \frac{12493177}{2240000}>5.57,\qquad \beta=b/\Omega.}
 \tag{8}
\]

Equations (2) and (8) apply to the same family (1). Its compact parameter
interval and smooth local dependence supply a common positive interval
on which the mean receiver increases, K increases, and beta gains positive
curvature. No useful numerical duration is enclosed here.

The old trailing-seed deceleration certificate does not apply to this
different datum. Conversely, these new initial gains do not give a way
to recreate it from the old solution or from its own later endpoint.
The steep envelope is essential to the torque sign, so the old constant-
amplitude Kelvin calculation cannot be propagated as this field's evolution.

The next task is a quantitative estimate for the actual coupled leading-
envelope solution, including its spatial covariance gradient, vertical
pressure response, mean deformation, viscosity and inherited endpoint.

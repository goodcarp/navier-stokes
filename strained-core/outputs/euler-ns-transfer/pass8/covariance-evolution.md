# Actual covariance, mean torque, and the second mean-angular-momentum jet

These identities apply to the smooth ordinary Navier–Stokes solution from
the actual compact Pass7 datum. They do not hold a wave polarization fixed.
All pressure terms below use the whole-space pressure. This note derives
initial jets and an exact initial loss of fluctuation energy; it does not
infer an interval of torque persistence from those jets alone.

## 1. Exact covariance equation

Use cylindrical components and azimuthal average. Put
\[
 U_i=\overline{u_i},\qquad v_i=u_i-U_i,\qquad
 R_{ij}=\overline{v_iv_j},\quad R=R_{r\theta},\quad Z=R_{z\theta},
 \quad p'=p-\overline p .
\]
The averaged components of \(v\) vanish, and the mean and fluctuation
fields are separately divergence free. Define
\[
 D_U=\partial_t+U_r\partial_r+U_z\partial_z,\qquad
 \Delta_0=\partial_{rr}+r^{-1}\partial_r+\partial_{zz}.
\]
The exact equation for the off-diagonal covariance is
\[
\begin{aligned}
 D_U R={}&-\frac1r\partial_r\!\left(r\overline{v_r^2v_\theta}\right)
           -\partial_z\overline{v_rv_\theta v_z}\\
 &-R_{rr}(\partial_rU_\theta+U_\theta/r)
  -R_{rz}\partial_zU_\theta-R\partial_rU_r-Z\partial_zU_r\\
 &+\frac{2U_\theta R_{\theta\theta}-U_rR}{r}
  +\frac{\overline{v_\theta^3}-\overline{v_r^2v_\theta}}{r}\\
 &-\overline{v_\theta\partial_rp'}
  -\frac1r\overline{v_r\partial_\theta p'}\\
 &+\nu\left[(\Delta_0-2r^{-2})R
       -2\overline{\nabla_{\rm sc}v_r\cdot\nabla_{\rm sc}v_\theta}\right].
\end{aligned}                                                        \tag{1}
\]
Here \(\nabla_{\rm sc}f=(\partial_rf,r^{-1}\partial_\theta f,\partial_zf)\)
is the scalar-component gradient. The displayed geometric terms are
therefore necessary. The two vector-Laplacian cross terms cancel under
azimuthal averaging, because their sum is an angular derivative of
\(v_r^2-v_\theta^2\). The pressure line is the full local pressure
transport/strain contribution; it is generally nonzero.

This follows by subtracting the exact mean equations from
\[
 D u_r=u_\theta^2/r-p_r+
 \nu(\Delta u_r-u_r/r^2-2\partial_\theta u_\theta/r^2),
\]
\[
 D u_\theta=-u_ru_\theta/r-p_\theta/r+
 \nu(\Delta u_\theta-u_\theta/r^2+2\partial_\theta u_r/r^2),
\]
multiplying by \(v_\theta,v_r\), and averaging. Incompressibility
converts fluctuation scalar advection into the first line of (1).

The mean angular momentum \(G=rU_\theta\) satisfies exactly
\[
 \boxed{D_U G=\nu\Delta_*G+\mathcal T,\qquad
 \Delta_*=\partial_{rr}-r^{-1}\partial_r+\partial_{zz},\qquad
 \mathcal T=-\frac1r\partial_r(r^2R)-r\partial_z Z.}                    \tag{2}
\]
The axial covariance flux in (2) cannot be dropped during evolution.

## 2. Evaluate the actual initial datum on its flat seed plateau

Near each maximizing circle \((r_*,z=\pm4)\), the initial seed envelopes
and outer axial swirl envelope equal one. The mean there is
\[
 (U_r,U_\theta,U_z)=(cr,A V_0(r),-2cz),\qquad
 c=7/5,\quad A=A_\lambda,\quad V_0=rC_0(r).
\]
The fluctuation has exact components
\[
 v_r=-\lambda m\sin(m\theta+kr)/r,\qquad
 v_\theta=\lambda k\sin(m\theta+kr),\qquad v_z=0,
 \quad m=4,\quad k=20 .
\]
Consequently,
\[
 R_0=-\frac{\lambda^2mk}{2r},\quad
 (R_{rr})_0=\frac{\lambda^2m^2}{2r^2},\quad
 (R_{\theta\theta})_0=\frac{\lambda^2k^2}{2},\quad
 Z_0=(R_{rz})_0=0,\qquad
 \mathcal T_0=\frac{\lambda^2mk}{2r}.
\]
All initial third moments in (1) vanish: the seed has one nonzero angular
frequency, so no product of three such factors has zero frequency.
This is an initial cancellation, not a closure at later times.

Let \(p_m(r,z)\) denote the actual unit \(v\)-\(w\) mixed pressure mode
from Pass7, with
\[
 -\Delta_m p_m=s_m,\qquad
 \Delta_m=\partial_{rr}+r^{-1}\partial_r+\partial_{zz}-m^2r^{-2}.
\]
The full initial pressure component with angular frequency \(m\) is
\(\lambda A\operatorname{Re}(p_m e^{im\theta})\).
The meridional pump/seed mixed advection has zero divergence exactly,
and all other base/seed overlap lies on that affine plateau.
The seed self-pressure has only frequencies \(0,2m\); these pair to
zero with one seed factor in the initial covariance equation.
Thus this \(p_m\) retains every pressure component that contributes.

Equation (1) gives the actual partial time derivative
\[
\boxed{\begin{aligned}
 R_t(0)={}&-cR_0
 +\lambda^2A\left[\frac{k^2V_0}{r}
       -\frac{m^2}{2r^2}(V_0'+V_0/r)\right]\\
 &-\nu R_0\left[2k^2+\frac{2m^2+1}{r^2}\right]\\
 &-\frac{\lambda^2A}{2}\operatorname{Re}
       \left[e^{-ikr}\left(ik\,\partial_rp_m+\frac{m^2}{r^2}p_m\right)\right].
\end{aligned}}                                                       \tag{3}
\]
The scalar-gradient contraction in the viscous line is
\[
 \overline{\nabla_{\rm sc}v_r\cdot\nabla_{\rm sc}v_\theta}
 =-\frac{\lambda^2mk}{2r}\left(k^2+\frac{m^2}{r^2}\right).
\]
A direct vector-Laplacian calculation verifies the same coefficient.

The initial vertical fluctuation acceleration is \(-\partial_zp'\)
on this plateau. Therefore the actual axial covariance begins with
\[
 \boxed{Z_t(0)=-\frac{\lambda^2A}{2}
       \operatorname{Re}\left[ik\,e^{-ikr}\partial_zp_m\right].}        \tag{4}
\]
The two-packet datum is even in \(z\), but it is not reflection symmetric
about the individual plane \(z=4\). There is no valid symmetry argument
for dropping (4), or its axial derivative at \(z=4\).

## 3. A full-pressure formula for the initial torque derivative

Differentiating both covariance fluxes in (2) gives
\[
\begin{aligned}
 \mathcal T_t(0)={}&-c\mathcal T_0
 +\nu\mathcal T_0\left[\frac{2m^2+1}{r^2}-2k^2\right]\\
 &+\lambda^2A\left[
  \frac{m^2}{2r}(V_0''+V_0'/r-V_0/r^2)
  -k^2(V_0'+V_0/r)\right]+\mathcal P_T ,
\end{aligned}                                                        \tag{5}
\]
where an exact use of the full Poisson equation collapses all pressure
second derivatives to
\[
 \boxed{\mathcal P_T=
 \frac{\lambda^2A}{2}\operatorname{Re}
 \left[e^{-ikr}\left\{
  (k^2r+m^2/r+ik)\partial_rp_m-ik r s_m\right\}\right].}                \tag{6}
\]
Specifically, differentiating (3) and (4) produces
\(ik r(p_{rr}+p_{zz})\); substituting
\(p_{rr}+p_{zz}=-s_m-p_r/r+m^2p_m/r^2\) cancels the remaining
\(p_m\) term. This simplification uses the axial flux rather than
neglecting it.

On the plateau the source is explicit:
\[
 e^{-ikr}s_m=\frac2r
 \left[ik V_0'-k^2V_0-\frac{m^2}{r}V_0'\right].
 \tag{7}
\]
Hence the source part of (6) is simply \(\lambda^2A k^2V_0'\).
Only one complex radial derivative of the actual full mixed pressure
remains to be bounded.

## 4. Exact second derivative at an initially maximizing fixed circle

Set \(G_0(r)=A rV_0(r)\) locally and
\[
 L_0=-cr\partial_r+\nu(\partial_{rr}-r^{-1}\partial_r).
\]
At any initial maximizing circle, \(G_0'=0\) and all initial axial
derivatives vanish on a neighborhood. Differentiating (2) at the
fixed circle removes the terms \(U_{r,t}G_{0,r}+U_{z,t}G_{0,z}\)
exactly. Since \(G_t(0)=L_0G_0+\mathcal T_0\), equations (5)-(6) give
\[
 \boxed{
 G_{tt}(0,r_*,\pm4)
 =L_0^2G_0+\frac{\lambda^2m^2}{2r_*^2}G_0''
 +\nu\mathcal T_0\left[\frac{2m^2+4}{r_*^2}-2k^2\right]
 +\mathcal P_T.}                                                     \tag{8}
\]
Here, at the maximum,
\[
 L_0^2G_0=c^2r_*^2G_0''
 -2c\nu r_*G_0'''
 +\nu^2\left(G_0''''-\frac{2G_0'''}{r_*}
                      +\frac{3G_0''}{r_*^2}\right).
 \tag{9}
\]
The two explicit \(c\mathcal T_0\) contributions cancel.
The curvature term \(\lambda^2m^2G_0''/(2r_*^2)\) is nonpositive
and can be substantial. Its size cannot be assessed independently of
the pressure term (6).

Equations (6)-(9) reduce the second jet to cutoff derivatives through
order four, a maximizing-radius enclosure, and one whole-space
pressure-gradient enclosure. The companion slab residual estimate
provides such a gradient bound around the horizontal trial.
Even a negative certified second jet does not alone rule out a later
gain; a stated time interval needs an evolution/remainder estimate.

## 5. One exact actual sign: the initial fluctuation reservoir loses energy

Define the actual fluctuation energy
\(K(t)=\frac12\int_{\mathbb R^3}|u-\overline u|^2\,dx\), where the mean
is understood as the axisymmetric vector-field average.
Its pressure work integrates to zero by fluctuation incompressibility.
For the actual Pass7 initial datum,
\[
 \boxed{\begin{aligned}
 K'(0)={}&-2cK(0)-\nu\lambda^2\|\nabla w\|_2^2\\
 &+\pi\lambda^2mkA
   \left(\int_0^\infty r C_0'(r)e(r)^2\,dr\right)
   \left(\int_{\mathbb R}q_s(z)^2q_v(z)\,dz\right)<0 .
 \end{aligned}}                                                       \tag{10}
\]
The last term is the mean-swirl shear production
\(-\int R_{r\theta}(U_{\theta,r}-U_\theta/r)\).
Since \(C_0'\le0\), it transfers energy out of this fluctuation packet;
the affine pump and viscosity also decrease its energy. All gradients
in the dissipative norm are full Euclidean gradients.

This exact integrated pressure cancellation does not cancel the
pointwise pressure-strain term in (3). It shows that the initially
favorable mean torque draws from a decreasing seed reservoir; it
does not prove that local transfer immediately stops, nor that a
different packet could not receive energy.

The standalone checker verifies the local covariance production,
both versions of the viscous term, the full pressure/axial-flux
collapse, the second mean jet, and the shear-production coefficient.

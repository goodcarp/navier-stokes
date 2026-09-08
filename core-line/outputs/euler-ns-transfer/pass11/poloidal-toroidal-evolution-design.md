# Solenoidal, axis-regular representations for actual evolution

The immediate diagnostic can use the new compatible MAC projection. A whole-space reconstruction for validation should instead be expressed in smooth scalar potentials or integrable continuous spectral amplitudes. A finite sum of unwindowed Bessel waves, or a periodic axial Fourier series extended to all space, is not a finite-energy whole-space velocity.

The most direct localized reconstruction for this datum is

\[
u=\nabla\times(T e_z)+\nabla\times\nabla\times(P e_z). \tag{1}
\]

Both potentials of the unchanged initial datum can be taken compact and smooth. Later approximations may also use compact potentials, provided the resulting full nonlocal residual and exterior pressure tail are included in validation. This is not a claim that the actual NS velocity remains compact.

## Exact cylindrical formulas and initial potentials

For the angular convention \(P=\sum_m P_m(r,z)e^{im\theta}\), and likewise for \(T\), let
\(\Delta_{h,m}=\partial_{rr}+r^{-1}\partial_r-m^2/r^2\). Then

\[
u_{r,m}=\partial_{rz}P_m+\frac{im}{r}T_m,
\quad u_{\theta,m}=\frac{im}{r}\partial_zP_m-\partial_rT_m,
\quad u_{z,m}=-\Delta_{h,m}P_m. \tag{2}
\]

Their divergence is identically zero, before any sampling or linear solve. The vertical vorticity is \(\omega_{z,m}=-\Delta_{h,m}T_m\).

For \(s=r^2+z^2\), write \(F_f(s)=\int_s^\infty f(\sigma)\,d\sigma\). With the exact current profiles \(\Phi=\chi+(7/5)\eta\), \(\psi=\chi\), \(a=63/200\), and the unchanged axial/radial envelopes, use

\[
P_0=\frac z2F_\Phi(s),
\]
\[
T_{0,\mathrm{mean}}=\frac16F_\psi(s)-\frac s3\psi(s)
 +\frac{1020a^2}{2}q_v(z)F_\chi(r^2/a^2),
\]
\[
T_{0,\mathrm{seed}}=\frac14 e(r)q_s(z)\cos(4\theta-20r). \tag{3}
\]

Here the subscript zero denotes the initial time in (3). The mean potential is angular mode zero; the positive mode-four seed coefficient is \(e(r)q_s(z)e^{-20ir}/8\). Direct differentiation gives exactly

\[
(u_r,u_z)=(-r[\Phi+2z^2\Phi'],\;2z[\Phi+r^2\Phi']),
\quad u_{\theta,\mathrm{core}}=r[\psi+(2s/3)\psi'],
\]

the outer swirl \(1020r\chi(r^2/a^2)q_v(z)\), and the specified leading seed. The antiderivatives vanish outside the relevant supports. Thus (3) is the same initial velocity, including the explicit amplitude 1020; it introduces no pressure retuning or replacement datum. Certified one-dimensional cutoff integration can enclose any antiderivative values needed for a numerical representation.

## Axis regularity should be built into the potentials

For every signed mode let \(n=|m|\), and represent

\[
P_m=r^n p_m(r^2,z),\qquad T_m=r^n t_m(r^2,z), \tag{4}
\]

with functions smooth in their displayed arguments. In Cartesian coordinates the angular factor becomes \((x+i y)^m\) for \(m\ge0\), and \((x-i y)^{|m|}\) for \(m<0\). Consequently (4) is smooth on the axis and (1) remains smooth there.

For \(s_h=r^2\) and \(f_m=r^n f(s_h,z)\), use the nonsingular identities

\[
\Delta_{h,m}f_m=4r^n[(n+1)f_s+s_hf_{ss}],
\]
\[
(\partial_r\pm m/r)f_m=(n\pm m)r^{n-1}f+2r^{n+1}f_s. \tag{5}
\]

When the coefficient \(n\pm m\) vanishes, remove that term symbolically. These identities imply the required helical regularity

\[
u_{r,m}+iu_{\theta,m}=r^{|m+1|}F_+(r^2,z),
\quad u_{r,m}-iu_{\theta,m}=r^{|m-1|}F_-(r^2,z),
\quad u_{z,m}=r^{|m|}F_z(r^2,z). \tag{6}
\]

For mode zero, \(u_r,u_\theta\) are odd in radius and \(u_z\) is even. Mode four requires powers 5, 3, and 4 in (6), respectively. Zero radial flux at the axis does not alone enforce these conditions. Near-axis nonlinear expressions should be evaluated from regular factored formulas or the Cartesian polynomial representation, not by separately dividing small, approximately canceling values by \(r\).

## Continuous Fourier–Hankel crosswalk

Use a unitary axial Fourier transform with frequency \(k\), and radial Hankel synthesis \(\int_0^\infty A(\kappa,k)J_n(\kappa r)\kappa\,d\kappa\). Use the same continuous radial frequency \(\kappa\) for the signed Bessel orders \(m+1,m-1,m\). The corresponding amplitudes of \(u_r+iu_\theta,u_r-iu_\theta,u_z\) are

\[
A_+=i\kappa(\widehat T-k\widehat P),\quad
A_-=i\kappa(\widehat T+k\widehat P),\quad
A_z=\kappa^2\widehat P. \tag{7}
\]

For arbitrary velocity amplitudes define

\[
D=(\kappa/2,-\kappa/2,ik),\quad
G=(-\kappa,\kappa,ik)^T,\quad K^2=\kappa^2+k^2.
\]

The exact divergence, gradient, projection, and diffusion symbols are

\[
\widehat{\operatorname{div}u}=DA,\qquad
\widehat{\nabla p}=G\widehat p,\qquad
\Pi A=\left(I+\frac{GD}{K^2}\right)A,\qquad
\widehat{\Delta u}=-K^2A. \tag{8}
\]

Here \(DG=-K^2\), and \(G\) is minus the adjoint of \(D\) in the kinetic metric \(W=\operatorname{diag}(1/2,1/2,1)\). Thus \(\Pi\) is a kinetic-energy-orthogonal projection. Formula (7) lies in its range. Signed Bessel orders matter: \(J_{-n}=(-1)^nJ_n\); replacing them by absolute orders without changing amplitude signs breaks (8), especially at mode zero.

With the angular convention above, the energy is

\[
E=\pi\sum_m\int_{\mathbb R}\int_0^\infty
\left[\frac{|A_+|^2+|A_-|^2}{2}+|A_z|^2\right]\kappa\,d\kappa\,dk
\]
\[
=\pi\sum_m\int\!\int
\kappa^2\bigl(|\widehat T|^2+K^2|\widehat P|^2\bigr)
\kappa\,d\kappa\,dk. \tag{9}
\]

The dissipation norm has the first integrand multiplied by \(2K^2\). These are full signed-mode sums; a positive-mode-only implementation must include the conjugate-mode factor of two.

At \(\kappa=0\), the inversion from velocity to these potentials is degenerate. The correct state space is the weighted potential norm in (9); it is not an unweighted L2 norm for \(P,T\). Constant horizontal Fourier modes are measure-zero in the continuous transform and do not supply finite-energy horizontally constant velocities. Their treatment in a discrete zero-frequency bin still needs an explicit convention. Avoid division by tiny \(\kappa^2\) by evolving velocity amplitudes in the range of \(\Pi\), or by using a velocity Gram matrix for the potential basis.

Independent Bessel-zero grids for orders \(m+1,m-1,m\) generally do not share radial frequencies and cannot be substituted into (8). The companion shared-frequency design/checker by `independent_direction` addresses this implementation issue. A discrete sum over isolated \(\kappa\) values is a sum of non-L2 Bessel waves. A finite-energy whole-space spectral reconstruction requires integrable amplitude functions, for example explicitly represented spectral wave packets, together with their transform/truncation errors. The same warning applies to isolated axial Fourier frequencies.

## Nonlinearity and potential evolution

The full cylindrical nonlinear term is

\[
F_r=u_r\partial_ru_r+(u_\theta/r)\partial_\theta u_r+u_z\partial_zu_r-u_\theta^2/r,
\]
\[
F_\theta=u_r\partial_ru_\theta+(u_\theta/r)\partial_\theta u_\theta
 +u_z\partial_zu_\theta+u_ru_\theta/r,
\]
\[
F_z=u_r\partial_ru_z+(u_\theta/r)\partial_\theta u_z+u_z\partial_zu_z. \tag{10}
\]

Modes multiply by convolution. The angular-basis terms in (10) must be retained. With input modes \(|m|\le M\), the quadratic term reaches \(2M\). An angular pseudospectral calculation needs more than \(3M\) samples to obtain the projected \(|m|\le M\) products without aliasing, and more than \(4M\) samples to recover all product modes without aliasing. Discarded modes remain part of the full residual.

Let \(p\) solve \(-\Delta p=\operatorname{div}F\) on the whole space. Then \(\Pi F=F+\nabla p\). The exact potential equations, with appropriate horizontal inverse-Laplacian gauges, are

\[
P_t=\nu\Delta P+\Delta_h^{-1}(\Pi F)_z,
\qquad
T_t=\nu\Delta T+\Delta_h^{-1}(\nabla\times F)_z. \tag{11}
\]

The toroidal equation eliminates pressure by curl; the poloidal equation still contains the full pressure projection. These formulas do not license a pressure-free planar evolution. In common-frequency amplitudes they read
\(\widehat P_t=-\nu K^2\widehat P-(\Pi\widehat F)_z/\kappa^2\),
\(\widehat T_t=-\nu K^2\widehat T-\widehat{(\nabla\times F)_z}/\kappa^2\).

## A localized Galerkin route without a radial pressure boundary

Choose compact functions \(p_m(s_h,z),t_m(s_h,z)\) in (4), using local polynomial or B-spline bases, with regular first radial spans and sufficiently smooth zero extensions. The induced vector basis functions \(\phi_a\) from (1) are exactly divergence free and axis regular. For a finite combination \(v=\sum c_a(t)\phi_a\), exact whole-space weak projection gives

\[
M\dot c+\nu Kc+C(c,c)=0,
\quad M_{ab}=\int\phi_a\cdot\phi_b,
\quad K_{ab}=\int\nabla\phi_a:\nabla\phi_b,
\quad C_a=\int\phi_a\cdot(v\cdot\nabla v). \tag{12}
\]

The mass matrix is positive on an independent velocity basis, and stiffness is nonnegative. Pure poloidal and toroidal fields are orthogonal in both these forms by horizontal integration by parts. Compact support gives locality and sparse radial/axial matrices. Remove any potential gauge that represents zero velocity. The velocity Gram norm, not an arbitrary potential coefficient norm, is the energy metric.

With exact integration and the full skew transport identity,
\(d(c^*Mc/2)/dt=-\nu c^*Kc\). Quadrature, angular truncation, time integration and any mass/stiffness approximations require their own error budgets. A finite Galerkin solution has zero residual only against the selected test space; its omitted full residual does not vanish.

Smoothness matters for the H4 bridge. A sufficient spatial requirement is \(P\in H^8\), \(T\in H^7\), giving \(v\in H^6\), with time derivatives sufficient for \(v_t\in H^4\). Tensor-product B-splines of degree at least eight for \(P\) and seven for \(T\), simple interior knots and appropriately smooth compact zero extensions can meet this requirement; ordinary cubic splines do not. A C-infinity basis is another option. Higher-derivative observables require correspondingly more regularity. The exact initial potentials (3) should be retained explicitly or their approximation error enclosed.

The current MAC periodic-cylinder run is an exploratory source of approximate trajectories. Its pressure solve is a different operation from whole-space reconstruction. Converting its data into (4) must enclose the representation error. An axial or radial cutoff should be applied inside the potentials so that (1) remains divergence free. The resulting join terms change the velocity and the NS residual and must be computed. Merely multiplying a velocity by a spatial cutoff does not preserve incompressibility.

## Compact reconstruction still has a noncompact residual

For a compact potential reconstruction with velocity support in \(|x|\le R\), define its full physical residual

\[
\mathcal R=v_t-\nu\Delta v+(v\cdot\nabla)v+\nabla p,
\qquad p=\partial_i\partial_jN*(v_iv_j),\quad N=\frac1{4\pi|x|}.
\]

Outside the support, \(\mathcal R=\nabla p\), normally of order \(|x|^{-4}\). This exterior contribution must be integrated in the H4 residual even when the approximate velocity is identically zero there. For
\(C_j=\sup_{|x|=1}\|\nabla^{j+3}N(x)\|_F\), \(E_v=\|v\|_2^2/2\), and \(R_* >R\),

\[
\|\nabla^j\mathcal R(x)\|_F
\le\frac{2E_vC_j}{(|x|-R)^{j+4}},\qquad |x|\ge R_*,
\]

and its squared L2 tail is at most \(4\pi(2E_vC_j)^2\) times

\[
\frac{h^{-2j-5}}{2j+5}
 +\frac{2Rh^{-2j-6}}{2j+6}
 +\frac{R^2h^{-2j-7}}{2j+7},\qquad h=R_*-R. \tag{13}
\]

Each \(C_j\) can be given an explicit analytic upper bound from derivatives of the Newton kernel; (13) is a recipe, not an evaluated certificate. It explains how a compact exactly solenoidal approximate field can enter the whole-space error bridge without asserting compactness of the actual solution. Local residual quadrature, pressure inside \(R_*\), and these tails all remain necessary. The eventual endpoint is the inherited evolved field; the initial packets and potentials are not reset there.

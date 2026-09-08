# An explicit short-time gain for the compact active core

**Proved estimate, conditional only on the explicit data and positive pressure margin below.** The compact datum from `work/pass2/active-core-pressure.md` produces a definite positive central vorticity gain over an explicitly specified time interval. The bounds are uniform for viscosity in a prescribed bounded interval, including the Euler endpoint. They retain the full pressure and viscosity. The gain can be extremely small; this is not a return, large-amplification, or regeneration theorem.

## 1. Fixed profiles and the positive margin

Fix the two smooth, divergence-free, compact profiles from the earlier calculation:
\[
 c_1(x,y,z)=\psi(x^2+y^2+z^2)(-y,x,0),
 \qquad v=\operatorname{curl}(\chi e_z),
\]
where \(\psi=1\) near the origin, \(\psi=0\) outside a ball of radius \(R\), and
\[
 \chi=\eta\!\left(\frac{x^2+2y^2+(z-d)^2}{\epsilon^2}\right)
       +\eta\!\left(\frac{x^2+2y^2+(z+d)^2}{\epsilon^2}\right).
\]
Here \(\eta\) is a nonzero smooth bump flat at its support boundary, \(0<\epsilon<d/5\), and \(R<d-\epsilon\). The two velocity supports are disjoint. Both velocities are odd under full spatial inversion. The packet support lies in \(x^2+y^2\le z^2/16\).

For the decaying pressure of \(v\), set
\[
 C_v=-\partial_{zz}p[v](0)>0,
 \qquad -\Delta p[v]=\partial_i\partial_j(v_iv_j).
\]
The earlier exact pressure calculation gives the explicit lower bound
\[
 \tag{1} C_v\ge c_v^-:=\frac{6}{4\pi}\int_{\mathbb R^3}\frac{|v(x)|^2}{|x|^5}\,dx>0.
\]
Take \(\Omega>0\), \(A\in\mathbb R\), and
\[
 \tag{2} u_0=\Omega c_1+Av,\qquad
 K=A^2C_v-\frac25\Omega^2>0.
\]
Every formula below remains valid with any explicitly known lower margin \(0<k_0\le K\) in place of \(K\) in the time and gain bounds. In particular, one can use \(k_0=A^2c_v^--2\Omega^2/5\) whenever this is positive. No numerical evaluation of an unknown solution is needed to specify that margin.

For an integer \(s\ge0\), use the following definite Sobolev norm, including all velocity components:
\[
 \tag{3} \|w\|_s^2=\sum_{|\alpha|\le s}\|\partial^\alpha w\|_{L^2(\mathbb R^3)}^2.
\]
Set \(C_{10}=\|c_1\|_{10}\), \(V_{10}=\|v\|_{10}\). Disjoint support also holds for their derivatives, so the initial norm is exactly
\[
 \tag{4} X_0=\|u_0\|_{10}
      =\big(\Omega^2C_{10}^2+A^2V_{10}^2\big)^{1/2}.
\]
These are finite integrals of the specified profiles. Narrow outer packets can make \(V_{10}\) very large; this cost is retained.

## 2. The quantitative theorem

Fix a maximum viscosity \(\nu_*\ge0\). Define the universal integers
\[
 C_E=65536,\qquad C_8=16384,\quad C_6=2048,\quad C_4=512,
\]
and the following explicit positive quantities:
\[
 \begin{aligned}
 V&=2X_0, & T_E&=\frac{1}{4C_EX_0},\\
 M_1&=3\nu_*V+C_8V^2,\\
 M_2&=(3\nu_*+2C_6V)M_1,\\
 M_3&=(3\nu_*+2C_4V)M_2+2C_4M_1^2.
 \end{aligned}
 \tag{5}
\]
For any \(0<k_0\le K\), put
\[
 \tag{6}
 \boxed{\displaystyle
 \tau=\min\left\{T_E,\frac{k_0}{M_2},\frac{\Omega k_0}{M_3}\right\}>0,
 \qquad g=\frac{k_0\tau^2}{4}>0.}
\]

For every \(\nu\in[0,\nu_*]\), the unique smooth solution of
\[
 u_t+u\cdot\nabla u+\nabla p=\nu\Delta u,
 \qquad\nabla\cdot u=0,
 \qquad u(0)=u_0
\]
exists on \([0,T_E]\). At \(\nu=0\) this denotes the corresponding Euler solution. For every \(t\in[0,\tau]\),
\[
 \tag{7}
 u(t,0)=0,\qquad
 \partial_z u_z(t,0)\ge\frac{k_0t}{2},
\]
and, writing \(\omega=\operatorname{curl}u\),
\[
 \tag{8}
 \partial_t\omega_z(t,0)\ge\Omega k_0t,
 \qquad
 \omega_z(t,0)\ge 2\Omega\left(1+\frac{k_0t^2}{4}\right).
\]
In particular,
\[
 \tag{9} \boxed{\omega_z(\tau,0)\ge2\Omega(1+g).}
\]
The origin is an actual material trajectory. The datum has finite energy and the external force is identically zero. Only the initial velocity is asserted to have compact support; the evolved velocity includes its actual nonlocal pressure response and viscous spreading.

There is also a genuine localized velocity gain, uniform in the radius of a radial receiver inside the initially rigid core. Let \(Jx=(-y,x,0)\), choose a fixed nonzero smooth \(\chi_0\ge0\) supported in \([0,1)\), and put
\[
 \phi_r(x)=\chi_0(|x|^2/r^2)Jx,\qquad
 a_r(t)=\frac{\langle u(t),\phi_r\rangle_{L^2}}
                 {\Omega\int\chi_0(|x|^2/r^2)|Jx|^2\,dx}.
\]
For every \(r>0\) whose receiver support lies inside the initially rigid core, the same \(\tau\) and \(g\) give
\[
 \tag{9a} a_r(t)\ge1+\frac{k_0t^2}{4}\quad(0\le t\le\tau),
 \qquad a_r(\tau)\ge1+g.
\]
This is a normalized velocity projection onto the specified compact divergence-free test field. No claim about the full velocity profile follows from this scalar projection alone.

## 3. Explicit Sobolev constants and lifespan

Here are elementary estimates supporting the numerical constants, rather than an unspecified Taylor interval. For the norm (3), scalar Sobolev embedding gives
\[
 \tag{10} \|h\|_\infty\le\|h\|_2.
\]
Here \(\|h\|_2\) denotes the order-two derivative norm (3), not the Lebesgue norm. For example, with the unitary Fourier transform, this norm dominates \(2^{-1/2}\|(1-\Delta)h\|_{L^2}\). The Fourier Cauchy–Schwarz bound has constant \(1/(2\sqrt\pi)<1\) for the derivative norm. In particular, a derivative of order \(q\) is bounded pointwise by the norm of order \(q+2\).

Retaining that scalar embedding constant gives the useful vector estimate
\[
 \tag{10a}
 \|\nabla w\|_{L^\infty,\mathrm{op}}
 \le \|\nabla w\|_{L^\infty,\mathrm{Frob}}
 \le\frac{\sqrt3}{2\sqrt\pi}\|w\|_3\le\|w\|_3.
\]
In the sum of order-two derivative norms of the first derivatives of \(w\), each derivative of \(w\) is counted at most three times. This proves the displayed factor and will prevent a receiver-radius loss.

Let \(P\) be the whole-space Leray projection and \(B(a,b)=P(a\cdot\nabla b)\). It is an \(L^2\) contraction commuting with spatial derivatives. Leibniz, (10), and Cauchy–Schwarz give, at each \(r\in\{4,6,8\}\),
\[
 \tag{11}
 \|B(a,b)\|_r\le C_r\|a\|_{r+2}\|b\|_{r+2},
 \qquad \|\Delta a\|_r\le3\|a\|_{r+2}.
\]
One sufficient product constant is \(3\,2^r\sqrt{\binom{r+3}{3}}\), which is smaller than the corresponding displayed \(C_r\). Indeed, for each differentiated product put the lower derivative factor in \(L^\infty\); its derivative order is at most \(\lfloor(r+1)/2\rfloor\). The other factor has at most \(r+1\) derivatives. Both are controlled by the indicated \(r+2\) norms. Summing the three directional products costs at most three, the binomial weights sum to \(2^{|\alpha|}\), and summing the output derivatives costs at most the square root of their count.

For \(X(t)=\|u(t)\|_{10}\), differentiation of the actual PDE, integration by parts, and incompressibility eliminate the pressure and the highest transport term in the energy identity. The viscous term is nonpositive. Every remaining commutator has both factors differentiated, with the lower derivative order at most five. Consequently
\[
 \frac12\frac{d}{dt}X^2
 +\nu\sum_{|\alpha|\le10}\|\nabla\partial^\alpha u\|_{L^2}^2
 \le C_EX^3,\qquad
 X'(t)\le C_EX(t)^2.
 \tag{12}
\]
The same counting gives \(3(2^{10}-1)\sqrt{\binom{13}{3}}<C_E\), so this specified \(C_E\) suffices. Thus, throughout the smooth lifespan,
\[
 \tag{13} X(t)\le\frac{X_0}{1-C_EX_0t}.
\]
The standard smooth local existence and continuation construction, applied with this a priori bound, supplies a solution beyond \(T_E\) and gives \(X(t)\le4X_0/3<V\) on \([0,T_E]\), uniformly in \(\nu\ge0\). The energy argument does not divide by viscosity. One can obtain the Euler endpoint by the same energy construction or by the uniform vanishing-viscosity limit. The requisite local theory and uniform energy method are also stated in [Tao's primary lecture notes, Theorem 1 and Corollary 3](https://terrytao.wordpress.com/2018/10/09/254a-notes-3-local-well-posedness-for-the-euler-equations/). Here (12)–(13) specify a concrete, conservative lifespan constant in the chosen norm.

Higher smoothness persists on this interval by the usual high-order energy estimates with \(\|\nabla u\|_\infty\) controlled by (13). The finite time derivatives used next already close within the tracked order-ten bound.

## 4. A uniform third-time-derivative bound, with viscosity and pressure retained

The projected equation is exactly
\[
 u_t=\nu\Delta u-B(u,u).
\]
It contains the full pressure through \(P\); replacing \(B\) by unprojected advection would be a different evolution. Differentiating this exact identity twice more gives
\[
 \begin{aligned}
 u_{tt}&=\nu\Delta u_t-B(u_t,u)-B(u,u_t),\\
 u_{ttt}&=\nu\Delta u_{tt}-B(u_{tt},u)
                -2B(u_t,u_t)-B(u,u_{tt}).
 \end{aligned}
 \tag{14}
\]
By (11), (13), and \(\nu\le\nu_*\), these imply on the entire interval \([0,T_E]\)
\[
 \tag{15}
 \|u_t\|_8\le M_1,\qquad
 \|u_{tt}\|_6\le M_2,\qquad
 \|u_{ttt}\|_4\le M_3.
\]
No support separation at positive time has been used. In particular, these bounds include every core/packet interaction that develops after insertion. Pointwise evaluation using (10) now yields
\[
 \tag{16}
 \left|\partial_t^2\partial_z u_z(t,0)\right|\le M_2,
 \qquad
 \left|\partial_t^3\omega_z(t,0)\right|\le2M_3.
\]
The first inequality takes one derivative of \(u_{tt}\); the second takes the two scalar derivative terms defining \((\operatorname{curl}u_{ttt})_z\). Both have more than the required two Sobolev derivatives available after differentiation.

## 5. Initial derivatives and the explicit gain proof

For completeness, the initial pressure identities underlying the sign are
\[
 p[u_0]=\Omega^2p[c_1]+A^2p[v],\qquad
 \partial_{zz}p[c_1](0)=\frac25,
\]
where the first equality uses disjoint support only at time zero. The core calculation retains the distributional contact term in the Hessian of \(1/(4\pi|x|)\). Near the origin the initial field is exactly the affine rotation \(\Omega(-y,x,0)\). There,
\[
 \Delta u_0=0,\qquad \omega_0=2\Omega e_z,\qquad \omega_t(0,\cdot)=0.
\]
Thus the exact vorticity equation, including viscosity, gives the initial jets
\[
 \begin{aligned}
 q(0)&=0,&q'(0)&=K,
    &&q(t):=\partial_z u_z(t,0),\\
 h(0)&=2\Omega,&h'(0)&=0,&h''(0)&=2\Omega K,
    &&h(t):=\omega_z(t,0).
 \end{aligned}
 \tag{17}
\]
In particular \(\Delta\omega_t(0,0)=0\), because \(\omega_t(0,\cdot)\) vanishes on an open core neighborhood. These initial equalities are independent of viscosity; the subsequent remainder bounds (15) are not.

The integral remainder and (16) give
\[
 \tag{18}
 q(t)\ge Kt-\frac{M_2t^2}{2},\qquad
 h'(t)\ge2\Omega Kt-M_3t^2.
\]
For \(t\le\tau\), the definition (6) makes the first expression at least \(k_0t/2\), and the second at least \(\Omega k_0t\). Integrating the latter from zero gives (8)–(9). Equivalently, the actual third-order remainder satisfies
\[
 |h(t)-2\Omega-\Omega Kt^2|\le\frac{M_3t^3}{3}.
\]
This supplies a checkable remainder bound on a specified interval, rather than an appeal to unspecified sufficiently small times.

Inversion symmetry of the data is preserved by uniqueness: \(-u(t,-x)\) solves the same initial-value problem. Hence \(u(t,0)=0\), as used to identify the fixed central point with a material trajectory. The estimates themselves concern the full coupled velocity, not a prescribed affine strain.

### Uniform-radius localized velocity estimate

Here is a direct proof of (9a), using the local receiver identity developed in `work/pass3/rotational-receiver.md`. It also supplies its quantitative remainder.

In the initially rigid neighborhood,
\[
 u_t(0)=\nabla\Psi,\qquad
 \Psi=\frac{\Omega^2}{2}(x^2+y^2)-p[u_0],\qquad \Delta\Psi=0.
\]
The affine rotation and the differentiated full PDE give, in that neighborhood,
\[
 \tag{18a} u_{tt}(0)=-2\Omega J\nabla\Psi+\nabla Q
\]
for a scalar \(Q\). In fact \(\nu\Delta u_t(0)=0\) there, and
\((u_0\cdot\nabla)\nabla\Psi=\nabla(u_0\cdot\nabla\Psi)+\Omega J\nabla\Psi\).
The gradient \(\nabla Q\) includes the actual differentiated pressure.

Each \(\phi_r\) is divergence free. Define \(I_r=\int\chi_0(|x|^2/r^2)|x|^2\,dx>0\). Radial symmetry gives the denominator \(2\Omega I_r/3\). The harmonic first-moment identity
\[
 \int\chi_0(|x|^2/r^2)x_i h(x)\,dx=\frac{I_r}{3}\partial_i h(0)
\]
holds whenever \(h\) is harmonic in the receiver ball; it follows by taking the degree-one spherical harmonic component, or the harmonic mean-value identity. Applied to \(h=\partial_x\Psi,\partial_y\Psi\), with \(\Delta\Psi=0\), it yields the exact normalized jets
\[
 \tag{18b} a_r(0)=1,\qquad a_r'(0)=0,\qquad a_r''(0)=\Psi_{zz}(0)=K.
\]

All time derivatives of the solution remain odd in space, so \(u_{ttt}(t,0)=0\). Equations (10a) and (15) give
\[
 |u_{ttt}(t,x)|\le M_3|x|.
\]
Consequently
\[
 \tag{18c}
 |a_r'''(t)|\le
 \frac{M_3\int\chi_0(|x|^2/r^2)|x|\,|Jx|\,dx}{2\Omega I_r/3}
 \le\frac{3M_3}{2\Omega},
\]
uniformly in \(r\). There is no negative power of receiver radius. Taylor's integral remainder therefore gives
\[
 \tag{18d} a_r(t)\ge1+\frac{Kt^2}{2}-\frac{M_3t^3}{4\Omega}
            \ge1+\frac{k_0t^2}{4}
 \quad\text{for }t\le\tau.
\]
This estimate includes the full evolved packet/core coupling in its remainder.

For \(\nu>0\), define the **rotational-mode Reynolds observable** at this fixed Eulerian receiver by
\[
 \mathcal R_{\rm rot}(t;r)=\frac{\Omega a_r(t)r^2}{\nu}.
\]
This definition uses the angular-speed coefficient extracted by the displayed test function. Let \(r_0\) be any allowed initial receiver radius and choose
\[
 q^2=\frac{1+g/2}{1+g}\in(0,1).
\]
Then the same actual solution satisfies
\[
 \tag{18e}
 \mathcal R_{\rm rot}(\tau;qr_0)
 \ge\left(1+\frac g2\right)\mathcal R_{\rm rot}(0;r_0).
\]
Thus one specified local velocity observable increases even when measured at a slightly smaller, fixed receiving radius. The smaller ball is not a material core. This statement does not control the orthogonal velocity components, localize vorticity into the receiver, provide a prescribed profile, or restore the inducing packet arrangement.

The estimate also implies actual weighted kinetic-energy and RMS-velocity gain, rather than only a projection gain. For \(\chi_r=\chi_0(|x|^2/r^2)\), let
\[
 E_r(t)=\frac12\int\chi_r|u(t)|^2,\qquad
 U_r(t)=\left(\frac{\int\chi_r|u(t)|^2}{\int\chi_r}\right)^{1/2}.
\]
Weighted Cauchy–Schwarz gives \(E_r(t)\ge a_r(t)^2E_r(0)\), since rigid rotation saturates that inequality initially. Thus \(U_r(t)\ge a_r(t)U_r(0)\), and \(U_r(0)=c_\chi\Omega r\). The same choice of \(q\) proves the analogue of (18e) for the standard RMS-based receiving Reynolds number \(rU_r(t)/\nu\). The receiver remains fixed in space and may exchange energy with the surrounding fluid.

## 6. A fixed positive gain on a normalized parameter range

The theorem is uniform on any fixed compact parameter range with a positive rotation and pressure margin. A scale-explicit version is useful. Fix the spatial profiles above, an amplitude scale \(S>0\), and write
\[
 \Omega=S\theta,\qquad A=Sa,\qquad \nu=S\lambda.
\]
Assume
\[
 \tag{19}
 0<\theta_-\le\theta\le\theta_+,\qquad |a|\le a_+,
 \qquad a^2C_v-\frac25\theta^2\ge\kappa_->0,
 \qquad 0\le\lambda\le\lambda_*.
\]
Set
\[
 \bar X=(\theta_+^2C_{10}^2+a_+^2V_{10}^2)^{1/2},\qquad \bar V=2\bar X,
\]
and compute \(\bar M_1,\bar M_2,\bar M_3\) by (5) with \(\nu_*\) replaced by \(\lambda_*\) and \(V\) by \(\bar V\). Define
\[
 \tag{20}
 \tau_* =\min\left\{\frac1{4C_E\bar X},
                   \frac{\kappa_-}{\bar M_2},
                   \frac{\theta_-\kappa_-}{\bar M_3}\right\},
 \qquad g_* =\frac{\kappa_-\tau_*^2}{4}>0.
\]
Then every datum in this normalized class, with every allowed viscosity, satisfies
\[
 \tag{21} \omega_z(\tau_*/S,0)\ge2S\theta(1+g_*),
 \qquad a_r(\tau_*/S)\ge1+g_*.
\]
To verify the scaling, the bounds for \(u_t,u_{tt},u_{ttt}\) scale respectively as \(S^2,S^3,S^4\), while \(K\) scales as \(S^2\). All three restrictions on the physical time therefore scale as \(S^{-1}\). Taking \(S=1\) gives the ordinary bounded parameter box with \(\nu\in[0,\nu_*]\). A fixed physical viscosity and a range \(S\ge S_{\min}>0\) fit by choosing \(\lambda_*\ge\nu_*/S_{\min}\).

This also closes the quantitative inputs of the one-time radius-selection corollary in `work/pass3/rotational-receiver.md`. Fix a radius-scaling exponent \(0<\alpha<1/2\), distinct from the normalized packet amplitude \(a\) in (19). For any permitted initial radius \(r_0\), continuity of the actual RMS velocity in the final receiver radius supplies some \(q\in(0,1)\) such that
\[
 \frac{\mathrm{Re}_{qr_0}(\tau_*/S)}{\mathrm{Re}_{r_0}(0)}=q^{-\alpha},\qquad
 \frac{U_{qr_0}(\tau_*/S)}{U_{r_0}(0)}=q^{-1-\alpha},\qquad
 \frac{E_{qr_0}(\tau_*/S)}{E_{r_0}(0)}=q^{1-2\alpha}.
\]
Indeed the continuous function \(q^{\alpha+1}U_{qr_0}(\tau_*/S)/U_{r_0}(0)\) tends to zero as \(q\downarrow0\) and exceeds one at \(q=1\). Moreover (10a) gives \(\|\nabla u(\tau_*/S)\|_\infty\le S\bar V\), so every selected root satisfies the uniform compact bracket
\[
 \tag{21a}
 \left(\frac{\theta_-}{\sqrt{3/2}\,\bar V}\right)^{1/(\alpha+2)}
 \le q\le(1+g_*)^{-1/(\alpha+2)}<1.
\]
No monotonicity or unique root is asserted. This is a comparison of observables at two times of one solution, not a map returning the evolved field to the original profile family.

This is a same-datum, actual-PDE finite gain statement for each member of the class, with a common relative lower gain. It does not identify one member's output with another member's input. For fixed \(\Omega\) and \(A\to\infty\), the displayed global-norm bounds can give only \(\tau=O(\Omega A^{-2})\) and a certified \(g=O(\Omega^2A^{-2})\); making the initial pressure acceleration large does not make this guaranteed gain large. That is a limitation of the bound, not an upper bound on the actual amplification.

## 7. Scope of the progress

The earlier local sign calculation now has an explicit uniform interval, monotone central axial vorticity, a localized rotational velocity gain, and a finite increase of the specified rotational-mode Reynolds observable at a slightly smaller fixed receiver, for ordinary unforced Navier–Stokes. This does not control a full turnover of the inducing packet unless the specified bounds happen to permit it. It does not preserve packet placement, establish a circulation increase, control the full receiver profile, prove robustness to inherited older fields, or supply a return to the starting geometry. Those are separate next-stage obligations.

The accompanying `verify_quantitative_core_stage.py` checks the numerical Sobolev counting constants, the polynomial derivative bounds, the gain/remainder inequalities, and amplitude scaling. The PDE proof is the energy and differentiated-equation argument above; the script is not a numerical fluid simulation or a formal proof assistant verification.

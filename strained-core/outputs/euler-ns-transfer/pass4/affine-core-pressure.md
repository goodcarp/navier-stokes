# Exact initial pressure feedback for a compact strained, rotating core

**Result.** For the radial vector-potential extension specified below, the central pressure Hessian is independent of the radial cutoff profile:
\[
 \boxed{D^2p[u_L](0)=\operatorname{diag}\left(
 -\frac{12}{7}b^2+\frac45\Omega^2,
 -\frac{12}{7}b^2+\frac45\Omega^2,
 -\frac{18}{7}b^2+\frac25\Omega^2\right).}
\]
Adding disjoint outer packets of axial pressure strength \(A^2C_v\) gives the exact initial central projections
\[
 \boxed{b'(0)=\frac12A^2C_v-\frac57b^2-\frac15\Omega^2,
 \qquad \Omega'(0)=2b\Omega.}
\]
These are initial derivatives of the actual compact-data, unforced Navier–Stokes solution. They are not a closed evolution law for the same profile family. Later pressure, curvature, packet geometry, and core/packet interactions remain unknown.

## 1. An explicit compact divergence-free affine extension

Let \(\psi\in C_c^\infty([0,\infty))\) equal one in a neighborhood of zero and vanish for \(s\ge R^2\). The notation means a smooth function on the nonnegative half-line, extendible smoothly across zero. Positivity or monotonicity is not needed for the pressure calculation. Put
\[
 J=\begin{pmatrix}0&-1&0\\1&0&0\\0&0&0\end{pmatrix},\qquad
 L=\operatorname{diag}(-b,-b,2b)+\Omega J,
 \quad s=|x|^2.
\]
Define
\[
 \tag{1} u_L=\operatorname{curl}\left[-\frac{\psi(s)}3\,x\times(Lx)\right].
\]
Since \(\operatorname{tr}L=0\), \(\operatorname{curl}[x\times(Lx)]=-3Lx\). The product rule therefore gives the exact formula
\[
 \tag{2}
 u_L=\left(\psi+\frac{2s}3\psi'\right)Lx
       -\frac23\psi'(x\cdot Lx)x.
\]
This is smooth, compactly supported, divergence free, odd under spatial inversion, and exactly \(Lx\) wherever \(\psi=1\). The symmetric and azimuthal parts are respectively
\[
 \begin{aligned}
 u_S&=\left(\psi+\frac{2s}3\psi'\right)Sx
             -\frac23\psi'(x\cdot Sx)x,\\
 u_W&=\Omega\left(\psi+\frac{2s}3\psi'\right)Jx,
 \qquad S=\operatorname{diag}(-b,-b,2b).
 \end{aligned}
 \tag{3}
\]
Thus the extension is axisymmetric: \(u_S\) is poloidal and \(u_W\) is azimuthal. The swirl cutoff is \(\psi+2s\psi'/3\), not \(\psi\). At \(b=0\) it belongs to the same radial-swirl class as the earlier rotating core; its pressure coefficient is again \(2\Omega^2/5\).

## 2. Pressure source, angular average, and contact term

Use the decaying whole-space pressure
\[
 -\Delta p=g,\qquad g=\partial_i u_j\,\partial_j u_i
                  =\operatorname{tr}[(\nabla u)^2],
 \qquad N(x)=\frac1{4\pi|x|},\quad p=N*g.
\]
At the origin \(g(0)=\operatorname{tr}(L^2)=6b^2-2\Omega^2\). The distributional Hessian formula is
\[
 \partial_{zz}N=\operatorname{PV}\frac{3z^2-|x|^2}{4\pi|x|^5}
                        -\frac13\delta_0.
\]
Writing \(\mu=z/|x|\) and letting angle brackets denote the uniform spherical average, this yields
\[
 \tag{4}
 p_{zz}(0)=-\frac13(6b^2-2\Omega^2)
       +\frac12\int_0^\infty
          \left\langle(3\mu^2-1)g(\sqrt s\,\cdot)\right\rangle\frac{ds}{s}.
\]
The contact term cannot be omitted. The angular contribution is identically zero near zero, since the initial field is affine there, so the radial integral is ordinary and convergent after angular integration.

Direct differentiation of (2) gives no \(b\Omega\) term in \(g\). In particular, the pressure cross term between this core's poloidal and azimuthal velocities vanishes identically. The precise angular average is
\[
 \left\langle(3\mu^2-1)g\right\rangle=b^2\mathcal A_S(s)+\Omega^2\mathcal A_W(s),
 \tag{5}
\]
where all cutoff derivatives below are with respect to \(s\):
\[
 \begin{aligned}
 \mathcal A_S(s)&=\frac{16s}{105}\left[
     -4s^2\psi'\psi''+6s\psi\psi''
     +2s(\psi')^2+21\psi\psi'\right],\\
 \mathcal A_W(s)&=\frac{16s}{135}
                 (2s\psi''+5\psi')(2s\psi'+3\psi).
 \end{aligned}
 \tag{6}
\]
The checker derives (6) from the full Cartesian gradient and uses \(\langle\mu^{2n}\rangle=(2n+1)^{-1}\); it does not assume a pressure coefficient.

## 3. All radial-profile dependence cancels

Let
\[
 I=\int_0^\infty s(\psi')^2\,ds<\infty.
\]
The exact endpoint conditions give
\[
 \begin{aligned}
 \int\psi\psi'\,ds&=-\frac12,\\
 \int s\psi\psi''\,ds&=\frac12-I,\\
 \int s^2\psi'\psi''\,ds&=-I.
 \end{aligned}
 \tag{7}
\]
For the strain part of the principal-value contribution in (4), substitution gives
\[
 \begin{aligned}
 \frac12\int\mathcal A_S\frac{ds}{s}
 &=-\frac{32}{105}(-I)
    +\frac{16}{35}\left(\frac12-I\right)
    +\frac{16}{105}I+\frac85\left(-\frac12\right)\\
 &=-\frac47.
 \end{aligned}
 \tag{8}
\]
For the swirl part,
\[
 \begin{aligned}
 \frac12\int\mathcal A_W\frac{ds}{s}
 &=\frac{32}{135}(-I)
    +\frac{16}{45}\left(\frac12-I\right)
    +\frac{16}{27}I+\frac89\left(-\frac12\right)\\
 &=-\frac4{15}.
 \end{aligned}
 \tag{9}
\]
The coefficients of \(I\) cancel exactly in both cases. Adding the contact term in (4) proves
\[
 \tag{10} p_{zz}[u_L](0)=-\frac{18}{7}b^2+\frac25\Omega^2.
\]
Axisymmetry makes the two transverse Hessian entries equal and all off-diagonal entries zero. Using \(\Delta p(0)=-6b^2+2\Omega^2\) therefore gives
\[
 \tag{11}
 D^2p[u_L](0)=\operatorname{diag}\left(
 -\frac{12}{7}b^2+\frac45\Omega^2,
 -\frac{12}{7}b^2+\frac45\Omega^2,
 -\frac{18}{7}b^2+\frac25\Omega^2\right).
\]
The coefficients do not depend on the radial transition width or support radius. High Sobolev norms and quantitative persistence times still depend strongly on that profile; their cancellation from this particular initial pressure Hessian does not remove those later costs.

## 4. Disjoint outer packets and the full initial matrix derivative

Let \(v\) be a fixed smooth compact divergence-free field whose support is disjoint from the core support and avoids the origin. Define
\[
 u_0=u_L+Av,\qquad H_v=D^2p[v](0),\qquad C_v=-(H_v)_{zz}.
\]
Initial disjoint support gives \(u_0\otimes u_0=u_L\otimes u_L+A^2v\otimes v\), so
\[
 \tag{12} D^2p[u_0](0)=D^2p[u_L](0)+A^2H_v.
\]
The old horizontal packets supported in the axial cone have
\[
 C_v\ge\frac{6}{4\pi}\int\frac{|v(x)|^2}{|x|^5}\,dx>0.
\]
No initial core/outer pressure interaction has been dropped: its quadratic cross tensor is exactly zero because the supports are disjoint.

Let \(u\) be the unique local smooth solution of the actual unforced NS equation, with any fixed \(\nu\ge0\) (Euler at zero viscosity). At time zero the core is affine, so both \(\nu\Delta u_0\) and its gradient vanish near the origin. Also \(u_0(0)=0\). Differentiating the full equation in space gives the exact initial identity
\[
 \tag{13} \frac{d}{dt}\nabla u(0,0)=-L^2-D^2p[u_0](0).
\]
Define the central projections at nearby times by
\[
 b(t)=\frac12\partial_z u_z(t,0),\qquad
 \Omega(t)=\frac12\big(\partial_xu_y-\partial_yu_x\big)(t,0).
\]
At time zero, (11)–(13) yield
\[
 \tag{14}
 \boxed{b'=\frac12A^2C_v-\frac57b^2-\frac15\Omega^2,
 \qquad \Omega'=2b\Omega.}
\]

For a nonaxisymmetric outer field, these two projections do not describe the entire matrix derivative. Since \(\operatorname{tr}H_v=0\), set
\[
 E_v=H_v-\operatorname{diag}(C_v/2,C_v/2,-C_v).
\]
This symmetric matrix has \((E_v)_{zz}=0\) and zero trace. The full result is
\[
 \tag{15}
 L'(0)=\operatorname{diag}(-b',-b',2b')+2b\Omega J-A^2E_v.
\]
Thus the previous elliptic outer envelopes can generate transverse anisotropy immediately. The axial formulas (14) remain correct, but treating them as a closed axisymmetric matrix system would omit \(E_v\).

For a fully axisymmetric option, replace the elliptic envelope by a circular pair:
\[
 \chi_{\rm out}(x,y,z)=
 \eta\!\left(\frac{x^2+y^2+(z-d)^2}{\epsilon^2}\right)
 +\eta\!\left(\frac{x^2+y^2+(z+d)^2}{\epsilon^2}\right),
 \qquad v=\operatorname{curl}(\chi_{\rm out}e_z).
\]
With the same \(\epsilon<d/5\) and disjoint-support condition, \(v\) is purely azimuthal, odd, and supported in the same cone, so \(C_v>0\) and \(E_v=0\). Axisymmetry and oddness are preserved by uniqueness. The central gradient retains the two-parameter matrix form at later smooth times. This preserves the form of the central matrix, not the compact radial extension profile or the entire pressure decomposition.

In that symmetric case, the exact later-time identities are
\[
 \begin{aligned}
 b'(t)&=-2b(t)^2-\frac12p_{zz}(t,0)
                    +\frac\nu2\partial_z\Delta u_z(t,0),\\
 \Omega'(t)&=2b(t)\Omega(t)+\frac\nu2\Delta\omega_z(t,0).
 \end{aligned}
 \tag{16}
\]
The pressure and curvature terms must be evolved and estimated. Replacing them at positive time by their initial values in (11)–(14) would impose an unjustified closure.

## 5. A quantitative initial cone test

Assume \(b>0\), \(\Omega>0\), and write \(\beta=b/\Omega\). The normalized strain ratio initially obeys
\[
 \tag{17}
 \beta'=\Omega\left[\frac{A^2C_v}{2\Omega^2}
                      -\frac{19}{7}\beta^2-\frac15\right].
\]
Therefore keeping the strain-to-rotation ratio from decreasing at insertion requires precisely
\[
 \tag{18}
 \boxed{\frac12A^2C_v\ge\frac{19}{7}b^2+\frac15\Omega^2.}
\]
This is stronger than merely making \(b'\ge0\). At strict inequality, the initial ratio increases; at equality, its initial derivative is zero. Neither fact proves that the corresponding region is invariant for the later PDE.

For comparison with the external-reservoir calculation, choose a fixed smooth compact axisymmetric cutoff \(\zeta\) that is one on a neighborhood of the outer support and zero near the core. Write
\[
 Q_{ij}(x)=-\partial_{zzij}N(x),\qquad
 C_F(t)=\int\zeta(x)Q_{ij}(x)u_i(t,x)u_j(t,x)\,dx,
 \qquad \mathcal R(t)=\frac{C_F(t)}{\Omega(t)^2}.
\]
Then \(C_F(0)=A^2C_v\). For circular outer packets the nondiffusive initial acceleration on their support is poloidal: the azimuthal self-advection produces a radial acceleration and the full axisymmetric pressure gradient is poloidal. The tensor \(Q\) has no azimuthal–poloidal coupling. Consequently the nondiffusive part of \(C_F'(0)\) vanishes pointwise. Since \(\Delta Q=0\) away from the origin and all derivatives of \(\zeta\) vanish near the initial packet support, integration by parts gives
\[
 \tag{19}
 C_F'(0)=2\nu A^2\int Q_{ij}v_i\Delta v_j
       =-2\nu A^2\int Q_{ij}\partial_kv_i\partial_kv_j.
\]
Every Cartesian derivative \(\partial_kv\) is horizontal, and \(Q\) is positive definite on horizontal vectors in the support cone, by the same kernel estimate that gives \(C_v>0\). Thus \(C_F'(0)<0\) for \(\nu>0\) and nonzero packets, and it equals zero at \(\nu=0\). In either case, for \(b>0\),
\[
 \tag{20} \mathcal R'(0)=\frac{C_F'(0)}{\Omega^2}-4b\mathcal R(0)<0.
\]
These statements concern this specified evolving functional; at later times its sign, geometry and derivative must be recomputed.

Set the cone slack
\[
 \mathcal F=\mathcal R-\frac{38}{7}\beta^2-\frac25.
\]
Then \(\beta'(0)=\Omega\mathcal F/2\). Consequently, **on this particular initial boundary** \(\mathcal F=0\), one has \(\beta'=0\) and \(\mathcal F'=\mathcal R'<0\). This is an outward tangent at the boundary of the fixed lower-reservoir/normalized-strain cone. It does not rule out a stage begun with positive slack, a different invariant region, or an evolved reservoir geometry. The reservoir derivative requires its own evolving functional; \(A^2C_v\) cannot simply be treated as a conserved or prescribed constant.

The outward derivative of this observable cone does not, by itself, determine \(\beta''(0)\). The equality \(\beta'=\Omega\mathcal F/2\) was established only at insertion; differentiating it requires the evolving core-pressure/profile defect. The full neutral-boundary condition is instead \(\beta''(0)=-(p_{zz}'(0,0)+32b^3)/(2\Omega)\), as derived with all pressure and curvature terms in `full-feedback-gate-audit.md`.

## 6. What the computation changes

The core now carries both rotation and axial strain. Its full compact pressure imposes the explicit costs \(5b^2/7\) and \(\Omega^2/5\) in the initial axial-strain equation. These coefficients are not tunable by changing the radial cutoff of this extension. The outer pressure can overcome them at insertion, but its strength and geometry are dynamical quantities. A viable invariant family must account for those dynamics and for the viscosity-generated curvature terms in (16), rather than resetting the core to rigid rotation after each gain.

`derive_affine_core_pressure.py` constructs the extension, checks divergence and the pressure source, performs the spherical average exactly, verifies both radial cancellations, and checks the resulting full matrix derivative. It is an exact symbolic calculation, not a PDE time evolution or proof of cone invariance.

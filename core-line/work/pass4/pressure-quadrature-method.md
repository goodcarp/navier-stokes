# Whole-space axisymmetric quadrature for the initial pressure derivative

Status: analytic formulation and audit guidance for evaluating one compact smooth initial datum. This removes periodic Poisson-box effects; it does not certify a numerical sign, evolve NS, or establish a regenerative stage. The formula retains the contact term and the full nonlocal initial acceleration.

## 1. Source, acceleration and the contact term

Let \(U=\nabla u_0\), with matrix entries \(U_{ij}=\partial_j u_{0,i}\), and set
\[
g=\operatorname{tr}(U^2),\qquad -\Delta p=g,
\qquad a=u_t(0)=\nu\Delta u_0-(u_0\cdot\nabla)u_0-\nabla p.
\]
This is equivalently \(a=\nu\Delta u_0-P\operatorname{div}(u_0\otimes u_0)\). Its pressure tail must remain in the computation. Differentiation gives
\[
h:=g_t(0)=2\operatorname{tr}(U\nabla a),\qquad
\boxed{p_{zz}'(0)=\partial_{zz}(-\Delta)^{-1}h(0).}
\]
Although \(a\) is noncompact, \(h\) is compactly supported: it contains the compact factor \(U\). Thus a finite radial integration domain enclosing the initial velocity support is sufficient.

Write spherical radius \(\rho=|x|\), \(\mu=z/\rho\), and use coefficients
\[
g_\ell(\rho)=\frac{2\ell+1}{2}\int_{-1}^{1}g(\rho,\mu)P_\ell(\mu)\,d\mu
\]
and the same convention for \(h_\ell\). If the datum has the stated inversion symmetry, both sources are even in \(\mu\), so only even degrees occur. The distributional Hessian identity
\[
\partial_{zz}\frac1{4\pi|x|}
=\operatorname{PV}\frac{3z^2-|x|^2}{4\pi|x|^5}-\frac13\delta_0
\]
gives the exact scalar endpoint formula
\[
\boxed{p_{zz}'(0)=-\frac{h(0)}3+\frac25\int_0^R\frac{h_2(\rho)}\rho\,d\rho.}
\]
The same formula with \(g\) instead of \(h\) gives the initial \(p_{zz}(0)\). Smoothness implies \(h_2(\rho)=O(\rho^2)\); the radial integral is regular. Dropping \(-h(0)/3\) would change the answer. In bilinear pressure notation this formula is \(2\Pi(u_0,a)\), with \(\Pi\) understood distributionally or through the source representation above.

## 2. Whole-space Legendre Green functions

For every retained degree, the decaying whole-space solution is
\[
p_\ell(\rho)=\frac1{2\ell+1}
\left[\rho^{-\ell-1}\int_0^\rho s^{\ell+2}g_\ell(s)\,ds
+\rho^\ell\int_\rho^R s^{1-\ell}g_\ell(s)\,ds\right].
\]
There is no periodic image, artificial Dirichlet wall, or free exterior boundary constant. For stable implementation define
\[
I_\ell(\rho)=\int_0^\rho s(s/\rho)^{\ell+1}g_\ell(s)\,ds,
\qquad
J_\ell(\rho)=\int_\rho^R s(\rho/s)^\ell g_\ell(s)\,ds.
\]
The ratio powers are at most one. Then
\[
\begin{aligned}
p_\ell&=\frac{I_\ell+J_\ell}{2\ell+1},\\
p_\ell'&=\frac{-(\ell+1)I_\ell+\ell J_\ell}{(2\ell+1)\rho},\\
p_\ell''&=\frac{(\ell+1)(\ell+2)I_\ell+\ell(\ell-1)J_\ell}{(2\ell+1)\rho^2}-g_\ell.
\end{aligned}
\]
This avoids numerical differentiation of pressure and explicitly satisfies the radial Poisson ODE. Analytic integration of a linear source interpolant on each radial panel is valid, with second-order source-interpolation error; a small Poisson trace residual does not bound that error.

Only the final central functional uses degree two. Computing \(h\) still requires the **full pressure Hessian on the velocity support**, hence all pressure degrees needed for convergence. Keeping only \(p_0,p_2\) in that intermediate solve would discard genuine interactions.

## 3. Pressure Hessian in cylindrical components

Let \(s=\sqrt{1-\mu^2}\); here \(s\) is an angular factor, not an integration radius. Reconstruct \(p_\rho,p_{\rho\rho},p_\mu,p_{\rho\mu},p_{\mu\mu}\) from the radial coefficients and Legendre derivatives. In the cylindrical orthonormal basis, with cylindrical radius \(r=\rho s\),
\[
\begin{aligned}
p_{rr}={}&s^2p_{\rho\rho}-\frac{2\mu s^2}{\rho}p_{\rho\mu}
+\frac{\mu^2s^2}{\rho^2}p_{\mu\mu}
+\frac{\mu^2}{\rho}p_\rho
+\frac{\mu(2-3\mu^2)}{\rho^2}p_\mu,\\
p_{zz}={}&\mu^2p_{\rho\rho}+\frac{2\mu s^2}{\rho}p_{\rho\mu}
+\frac{s^4}{\rho^2}p_{\mu\mu}
+\frac{s^2}{\rho}p_\rho-\frac{3\mu s^2}{\rho^2}p_\mu,\\
p_{rz}={}&\mu s p_{\rho\rho}+\frac{s(1-2\mu^2)}{\rho}p_{\rho\mu}
-\frac{\mu s^3}{\rho^2}p_{\mu\mu}
-\frac{\mu s}{\rho}p_\rho+\frac{s(3\mu^2-1)}{\rho^2}p_\mu,\\
H_{\theta\theta}={}&p_r/r=p_\rho/\rho-\mu p_\mu/\rho^2.
\end{aligned}
\]
The other components involving \(\theta\) are zero. These formulas have regular axial limits; in particular \(p_{rr}=p_r/r\) on the axis. At the origin use the regular solid-harmonic expansion rather than evaluating expressions containing \(1/\rho\).

## 4. Axisymmetric source and differentiated source

Write the velocity components as \((v,w,q)=(u_r,u_\theta,u_z)\). Its Cartesian gradient represented in the cylindrical frame is
\[
U=\begin{pmatrix}
v_r&-w/r&v_z\\
w_r&v/r&w_z\\
q_r&0&q_z
\end{pmatrix}.
\]
Thus
\[
g=v_r^2+(v/r)^2+q_z^2+2v_zq_r-2(w/r)w_r.
\]
The convective acceleration components are
\[
n_r=vv_r+qv_z-w^2/r,\quad
n_\theta=vw_r+qw_z+vw/r,\quad
n_z=vq_r+qq_z.
\]
The vector Laplacian has components
\[
(\Delta u)_r=(\partial_{rr}+r^{-1}\partial_r+\partial_{zz}-r^{-2})v,
\quad
(\Delta u)_\theta=(\partial_{rr}+r^{-1}\partial_r+\partial_{zz}-r^{-2})w,
\quad
(\Delta u)_z=(\partial_{rr}+r^{-1}\partial_r+\partial_{zz})q.
\]
Compute \(a=\nu\Delta u-n-\nabla p\), and construct \(\nabla a\) with the same cylindrical gradient pattern. Equivalently, a useful source decomposition is
\[
\boxed{h=\nu h_\nu+h_0,\qquad
h_\nu=2\operatorname{tr}(U\nabla\Delta u),\qquad
h_0=-2\operatorname{tr}(U^3)-(v\partial_r+q\partial_z)g-2\operatorname{tr}(UH).}
\]
It follows from \(\nabla((u\cdot\nabla)u)=U^2+(u\cdot\nabla)U\) in Cartesian coordinates and trace invariance. This is often more economical than constructing every acceleration derivative. In the pressure contraction,
\[
\operatorname{tr}(UH)=v_rp_{rr}+(v/r)(p_r/r)+q_zp_{zz}+(v_z+q_r)p_{rz}.
\]
The actual nonlocal pressure is retained even though the final source is compact.

Axis regularity should be enforced analytically. Smooth axisymmetric data have \(v=rV(r^2,z)\), \(w=rW(r^2,z)\), \(q=Z(r^2,z)\). For example, with \(\xi=r^2\),
\[
(\Delta u)_r=r(8V_\xi+4\xi V_{\xi\xi}+V_{zz}),\qquad
(\Delta u)_z=4Z_\xi+4\xi Z_{\xi\xi}+Z_{zz},
\]
and the azimuthal formula is identical to the radial one with \(W\). These regular expressions avoid subtracting separate \(r^{-2}\) terms near the axis.

## 5. Exact affine-inner-ball corrections

If \(u_0=Mx\) in \(\rho<r_*\), where \(M=bD+\Omega J\), then \(g=6b^2-2\Omega^2\) is constant there. For every \(\ell\ge1\), the pressure in that ball is \(p_\ell(\rho)=C_\ell\rho^\ell\), with
\[
C_\ell=\frac1{2\ell+1}\int_{r_*}^R s^{1-\ell}g_\ell(s)\,ds.
\]
Its monopole is a constant minus \(g\rho^2/6\). The center contact can be evaluated directly as
\[
h(0)=-2\operatorname{tr}(M^3)-2\operatorname{tr}(M H(0))
=12b b'_0-8b\Omega^2.
\]
On the initially neutral strain-ratio boundary, \(b'_0=2b^2\), this is \(24b^3-8b\Omega^2\). The viscous contribution vanishes throughout the initial affine ball.

There is also an exact inner formula for the final degree-two source. Since
\[
\partial_{zz}(\rho^\ell P_\ell(\mu))=\ell(\ell-1)\rho^{\ell-2}P_{\ell-2}(\mu),
\]
and \(\operatorname{tr}(bD\nabla^2p)=3b p_{zz}\) for harmonic pressure, only pressure degree four contributes to \(h_2\) in the affine ball:
\[
\boxed{h_2(\rho)=-72b C_4\rho^2,\qquad
C_4=\frac19\int_{r_*}^R s^{-3}g_4(s)\,ds.}
\]
This identity includes every higher pressure mode: their derivatives project into other degrees. Therefore a midpoint radial mesh starting at \(\rho_0>0\) omits the known contribution
\[
\frac25\int_0^{\rho_0}\frac{h_2(\rho)}\rho\,d\rho
=-\frac{72}{5}bC_4\rho_0^2.
\]
Using the analytic formula over the entire inner affine interval also prevents small-radius angular cancellation noise from being amplified by \(1/\rho\).

## 6. Viscosity and amplitude economy

For fixed data the initial pressure does not depend on viscosity. Hence
\[
p_{zz}'(0;\nu)=T_0+\nu T_1
\]
is exactly affine in \(\nu\), with each coefficient evaluated by the same contact-plus-degree-two functional applied to \(h_0\) or \(h_\nu\). For the affine-core family \(h_\nu\) vanishes in the inner ball, so \(T_1\) has no origin contact term and can be integrated over the cutoff shells and outer supports. Direct analytic profile derivatives are preferable to differencing two nearby viscosity evaluations.

An integration identity makes \(T_1\) substantially cheaper: define
\[
S_2(x)=\sum_{k=1}^3\operatorname{tr}\big((\partial_kU(x))^2\big).
\]
These are matrix squares, not squared Frobenius norms; no positivity is implied. Since
\[
h_\nu=\Delta g-2S_2,
\qquad
\partial_{zz}(-\Delta)^{-1}\Delta g(0)=-\partial_{zz}g(0),
\]
and an affine core has \(\partial_{zz}g(0)=0\) and \(S_2=0\) on an inner ball, the viscosity coefficient is exactly
\[
\boxed{T_1=-\frac45\int_0^R\frac{(S_2)_2(\rho)}\rho\,d\rho.}
\]
Thus only **second velocity derivatives** are needed. For a general nonaffine datum the extra term \(-\partial_{zz}g(0)\) and the contact contribution \(2S_2(0)/3\) must be retained.

To evaluate this scalar in the cylindrical frame, let \(J\) be the generator of rotation in the \((e_r,e_\theta)\) plane. Then
\[
\boxed{S_2=\operatorname{tr}(U_r^2)+\operatorname{tr}(U_z^2)
+r^{-2}\operatorname{tr}([J,U]^2).}
\]
Here \(U_r,U_z\) are ordinary derivatives of the displayed cylindrical matrix entries. The commutator term accounts for the changing cylindrical basis: \(U_{\rm Cartesian}=R_\theta U R_\theta^T\), so its angular derivative is \(R_\theta[J,U]R_\theta^T\). Omitting this term would not compute the Cartesian invariant. The available second derivatives of \((v,w,q)\) suffice to form all terms. Use axis-regular limits and the known zero value in the affine core to avoid cancellation at small \(r\). This formula computes the actual viscosity contribution without another pressure solve or fourth cutoff derivatives.

If \(u_0\) is a linear combination of fixed profiles, \(T_0\) is homogeneous cubic in their amplitudes and \(T_1\) is homogeneous quadratic. Precomputing those multilinear coefficients can make parameter sweeps economical and expose cancellations. This algebraic fact does not justify omitting cross interactions: the pressure part couples separated velocity supports nonlocally. A favorable \(T_0\) alone does not give a fixed-positive-viscosity sign without controlling \(\nu T_1\).

## 7. Audit and sign-control requirements

The current exploratory `evaluate_pressure_gate.py` was checked against these formulas: its cylindrical gradients, convective derivatives, radial Green panel integrals, Hessian reconstruction, inviscid \(h=2\operatorname{tr}(U\nabla a)\), and contact sign agree. No modifications were made in this audit.

For reliability:

1. Enclose the actual support. A meridional cutoff in \(r\le\sigma\), \(|z|\le d+h\) requires \(R>\sqrt{\sigma^2+(d+h)^2}\), as well as enclosure of the core. The expression \(d+h+0.25\) is adequate for some parameters, not all. Enforce disjoint core/packet support and the intended swirl plateau conditions separately.
2. Compare initial \(p_{zz}\) from the source formula with the analytic core coefficient and the sign-definite outer stress formula. Anchoring the center contact to those values does not repair an underresolved pressure Hessian elsewhere.
3. Refine radial panels, angular quadrature and Legendre cutoff independently. Resolve narrow radial/axial cutoff layers and oversample products before taking their degree-two projection. The identity \(\operatorname{tr}H=-\sum_{\ell\le L}g_\ell P_\ell\) means the reported Poisson trace discrepancy measures angular source truncation; it does not test radial Green accuracy or all individual Hessian components.
4. Use the exact inner correction above. Check the radial Poisson equation, axis limits, reflection symmetry, and the compact-source identity \(\int g\,dx=0\). The exterior monopole moment must vanish for an exactly compact divergence-free datum.
5. Retain the separate contact, interior radial integral and viscosity contributions in output. A cancellation-dominated total requires errors substantially below the final margin, not merely below the individual terms.

Convergence under refinement is evidence, not rigorous sign certification. A certified sign needs a computable error bound for angular projection, radial integration and arithmetic, including the reconstructed pressure Hessian. One useful propagation estimate is \(|\delta h|\le2\|U\|_F\|\delta H\|_F\) for pressure-field error alone, followed by the exact contact-plus-degree-two functional. Near zero, use the analytic affine-ball formula or a bound preserving \(\delta h_2=O(\rho^2)\); a uniform constant bound inserted into \(\int d\rho/\rho\) would be useless. The final interval must lie strictly on one side of zero before asserting a certified sign.

### Independent viscosity check

`verify_pressure_viscosity.py` supplies a callable `viscosity_from_second_jets(ev, amplitudes)` that reuses an existing evaluator, without solving pressure. Its symbolic checks passed for a generic six-parameter axisymmetric polynomial, including the Cartesian/cylindrical identity and a nonzero angular commutator contribution. A moderate exploratory run at \(n_r=400,n_\mu=512\), using the broader geometry \(\sigma=1,h=0.7,d=2\) and swirl ratio \(0.4\), gave:

| Profile | Second-jet integral | Independent stress formula or exact value |
|---|---:|---:|
| Meridional outer field | −260.9202331 | −260.9156713 |
| Swirl outer field | 0.1254760143 | 0.1253326623 |
| Strain core | −0.0525885 | 0 |
| Rotation core | 0.0402479 | 0 |

The outer formulas agree to relative discrepancies about \(1.75\times10^{-5}\) and \(1.14\times10^{-3}\), respectively. The nonzero core residuals expose unresolved radial cancellation at this moderate resolution; they are not interpreted as contradicting the exact zero core coefficients. Mixed-profile quadratic additivity agreed to roundoff on the same grid. These are independent quadrature diagnostics, not error bounds or a certified sign for a total stage coefficient.

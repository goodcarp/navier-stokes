# Audit of a mixed source/stress evaluation for the pressure derivative

**Accepted: the split is exact.** It removes derivatives of the acceleration from the outer integral without discarding its nonlocal pressure contribution. Noncompactness of the acceleration creates no boundary term because its product with the outer datum is compactly supported.

Let \(u_0=u_c+u_o\), where both pieces are smooth, compact, and divergence free, and the support of \(u_o\) is separated from the origin. Use the **full** initial acceleration
\[
 u_1=\nu\Delta u_0-P_L\operatorname{div}(u_0\otimes u_0),
 \qquad\nabla\cdot u_1=0.
\]
It need not be compactly supported. With the pressure polarization used in the preceding notes,
\[
 p_{zz}'(0)=2\Pi(u_c,u_1)+2\Pi(u_o,u_1).
\]

For any two smooth divergence-free fields \(a,b\), the exact identity
\[
 \operatorname{tr}(\nabla a\nabla b)=\partial_i\partial_j(a_i b_j)
\]
holds. The tensor \((u_o)_i(u_1)_j\) is smooth and compactly supported away from the evaluation point. Two integrations by parts therefore give
\[
 \boxed{2\Pi(u_o,u_1)
       =-2\int_{\mathbb R^3}Q_{ij}(x)(u_o)_i(x)(u_1)_j(x)\,dx,
 \qquad Q_{ij}=-\partial_{zzij}N.}
\]
There is no contact term at the origin because this tensor vanishes on an open neighborhood of it. There is no boundary term at spatial infinity, regardless of the pressure tail of \(u_1\). The symmetry of the fourth-derivative kernel supplies the required symmetric contraction; an additional factor of two from symmetrizing the tensor would be an error.

For the core piece define
\[
 h_c(x)=2\operatorname{tr}[\nabla u_c(x)\nabla u_1(x)].
\]
The factor two is already included. Retaining the distributional contact term gives
\[
 \boxed{2\Pi(u_c,u_1)=-\frac13h_c(0)
       +\operatorname{PV}\int\frac{3z^2-|x|^2}{4\pi|x|^5}h_c(x)\,dx.}
\]
Because \(u_o\) is zero near the origin, \(h_c(0)\) equals the total source derivative there. At the neutral affine insertion it is \(24b^3-8b\Omega^2\), so this term contains the entire contact contribution \(-8b^3+(8/3)b\Omega^2\).

For axisymmetric data write the conventional Legendre expansion
\[
 h_c(r,\mu)=\sum_{\ell\ge0}h_{c,\ell}(r)P_\ell(\mu),\qquad
 h_{c,2}(r)=\frac52\int_{-1}^{1}h_c(r,\mu)P_2(\mu)\,d\mu.
\]
If the core support lies in \(r\le1\), the core formula becomes
\[
 \boxed{2\Pi(u_c,u_1)=-\frac13h_c(0)
                         +\frac25\int_0^1\frac{h_{c,2}(r)}r\,dr.}
\]
Smoothness implies \(h_{c,2}(r)=O(r^2)\) near zero, so the angularly reduced radial integral is convergent. The factor \(2/5\) comes from \(3\mu^2-1=2P_2(\mu)\) and Legendre orthogonality.

Only the **final Hessian extraction** uses the degree-two coefficient. The pressure constructing \(u_1\) must retain the full modes of the full datum: truncating that pressure to degree two before multiplying by the core or outer profiles generally changes the answer. On the outer support the displayed stress integral needs \(u_1\), hence only the first spatial derivative of that full pressure, rather than its Hessian. The core term still requires \(\nabla u_1\).

The continuous equality assumes exact solenoidality. If an approximate acceleration has \(d=\nabla\cdot u_1\ne0\), then instead
\[
 \operatorname{tr}(\nabla u_o\nabla u_1)
 =\partial_i\partial_j((u_o)_i(u_1)_j)-u_o\cdot\nabla d.
\]
Thus a numerical comparison of source and stress formulas also checks the pressure-solve/divergence residual; it is not legitimate to infer their equality for an unresolved approximate field without that error control. This does not affect the exact proposed evaluation identity.

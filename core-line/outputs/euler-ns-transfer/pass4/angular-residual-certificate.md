# An angular-tail error certificate for the mixed pressure evaluation

**Proved analytic bounds; no numerical certification is claimed.** A high-angular source residual supported in a known ball controls the pressure-gradient error with an explicit inverse-angular-frequency factor. For the selected radial affine core, pressure modes above degree four have exactly zero contribution to the core part of the mixed evaluation. Consequently their entire pressure-derivative error is bounded by the outer stress estimate below. The required residual norms and the other numerical errors have not yet been enclosed rigorously.

## 1. A global energy estimate for a compact high-angular source

Let \(P_{\le L}\) denote exact orthogonal spherical-harmonic projection onto degrees \(0,\ldots,L\), at each radius. In the axisymmetric setting this is ordinary Legendre projection. Let
\[
 e=(I-P_{\le L})g,\qquad \operatorname{supp}g\subset B_R,
 \qquad \phi=(-\Delta)^{-1}e=N*e.
\]
The source is compact but the Newton potential is generally not. Rotational invariance of the Newton operator implies that \(\phi\) has the same high-angular restriction. Let \(\ell_*\) be its smallest allowed angular degree, and put \(\lambda_* =\ell_*(\ell_*+1)\). In general one may take \(\ell_*=L+1\). If the source is even under inversion, as for the selected odd velocity datum, only even angular degrees occur and one may take
\[
 \tag{1} \ell_*=2\left(\left\lfloor\frac L2\right\rfloor+1\right).
\]

Angular Poincare on every sphere and the spherical decomposition of the gradient give
\[
 \begin{aligned}
 \|\phi\|_{L^2(B_R)}^2
 &\le\lambda_*^{-1}\int_0^R r^2
                       \|\nabla_{S^2}\phi(r,\cdot)\|_2^2\,dr\\
 &\le\frac{R^2}{\lambda_*}\|\nabla\phi\|_{L^2(\mathbb R^3)}^2.
 \end{aligned}
 \tag{2}
\]
The Newton solution belongs to the homogeneous energy space: compact \(L^2\) data are in \(L^{6/5}\), and the weak Poisson solution has \(\nabla\phi\in L^2\). Testing the global equation against \(\phi\), justified in that energy space or by cutoff approximation, yields
\[
 \|\nabla\phi\|_2^2=\int_{B_R}e\phi
 \le\|e\|_2\|\phi\|_{L^2(B_R)}.
\]
Combining with (2) proves the explicit estimate
\[
 \tag{3} \boxed{\|\nabla(-\Delta)^{-1}e\|_2
       \le\frac{R}{\sqrt{\ell_*(\ell_*+1)}}\|e\|_2.}
\]
The global energy integral includes the noncompact potential tail. It has not been truncated at \(r=R\), so no zero-boundary condition at the source radius has been imposed.

With conventional axisymmetric coefficients
\(g_\ell(r)=\frac{2\ell+1}{2}\int_{-1}^{1}g(r,\mu)P_\ell(\mu)\,d\mu\), the exact tail norm is
\[
 \tag{4} \|e\|_2^2=4\pi\int_0^R r^2
          \sum_{\ell>L}\frac{|g_\ell(r)|^2}{2\ell+1}\,dr.
\]
This identity specifies the norm to certify. A numerical finite angular sum does not, by itself, bound the unevaluated tail.

## 2. The bound applies to the mixed estimator, with its exact definition retained

Keep the datum \(u_0=u_c+u_o\) fixed. The outer datum is compact and separated from the origin. Define the linear mixed functional
\[
 \begin{aligned}
 \mathcal M(w)&=\mathcal C(w)-2\int Q_{ij}(u_o)_iw_j,\\
 \mathcal C(w)&=\partial_{zz}\left[N*
                \big(2\operatorname{tr}(\nabla u_c\nabla w)\big)\right](0),
 \qquad Q_{ij}=-\partial_{zzij}N.
 \end{aligned}
 \tag{5}
\]
The core functional includes its distributional contact term. For the exact solenoidal acceleration \(u_1\), the previously audited integration by parts gives \(p_{zz}'(0)=\mathcal M(u_1)\).

Let \(p=N*g\), and first consider the exact angularly truncated pressure \(p^L=N*P_{\le L}g\), with exact radial solution and decay conditions. The corresponding acceleration is
\[
 u_1^L=\nu\Delta u_0-(u_0\cdot\nabla)u_0-\nabla p^L,
 \qquad u_1-u_1^L=-\nabla\phi.
\]
The approximation being bounded is **\(\mathcal M(u_1^L)\)**. In general \(\nabla\cdot u_1^L=-e\), so its all-source and mixed evaluations are not interchangeable. The comparison here is valid by linearity of the mixed functional; it does not assume that the truncated acceleration is solenoidal.

If \(\mathcal C(-\nabla\phi)=0\), the difference is exactly the outer term. Therefore Cauchy–Schwarz and (3) imply
\[
 \tag{6}
 \boxed{|p_{zz}'(0)-\mathcal M(u_1^L)|
 \le \frac{2R}{\sqrt{\ell_*(\ell_*+1)}}
             \|Q u_o\|_2\,\|(I-P_{\le L})g\|_2.}
\]
The weight \(Q u_o\) is square-integrable because its support avoids the singular point. The following section proves the needed exact core cancellation for every \(L\ge4\) in the selected profile class.

## 3. Exact degree-four cutoff in the radial affine core

The core is the radial vector-potential extension
\[
 u_c=bS_\psi+\Omega W_\psi,
 \quad S_\psi=\operatorname{curl}\left[-\frac{\psi(|x|^2)}3x\times(S_0x)\right],
 \quad S_0=\operatorname{diag}(-1,-1,2).
\]
The rotating part is azimuthal and the pressure gradient is poloidal. Their bilinear pressure source has zero trace identically, so \(\mathcal C\) has no rotating-core contribution to a pressure error. Only the strain part needs examination.

Write a pressure mode as \(\phi(r,\mu)=f(r)P_\ell(\mu)\), with \(\mu=z/r\). Let \(p=\psi(r^2)\), \(p_1=\psi'(r^2)\), \(p_2=\psi''(r^2)\), where cutoff derivatives are with respect to their squared-radius argument. The unit strain extension in cylindrical coordinates \((\varrho,z)\) is exactly
\[
 (S_\psi)_\varrho=-\varrho(p+2z^2p_1),\qquad
 (S_\psi)_z=2z(p+\varrho^2p_1).
\]
Set \(T=\operatorname{tr}(\nabla S_\psi\,D^2\phi)\). Direct differentiation followed by integration by parts in \(\mu\) gives
\[
 \tag{7}
 \int_{-1}^{1}(3\mu^2-1)T\,d\mu
 =\int_{-1}^{1}P_\ell(\mu)
               [f E_0+f' E_1+f''E_2]\,d\mu,
\]
where the exact three angular weights are
\[
 \begin{aligned}
 E_0={}&r^2p_2(-60\mu^4+48\mu^2-4)
       +p_1(120\mu^4-108\mu^2+12)\\
      &+\frac p{r^2}(135\mu^4-126\mu^2+15),\\
 E_1={}&r^3p_2(60\mu^4-48\mu^2+4)
       +rp_1(132\mu^4-108\mu^2+8)\\
      &+\frac p r(81\mu^4-66\mu^2+5),\\
 E_2={}&(3\mu^2-1)^2(2r^2p_1+p).
 \end{aligned}
 \tag{8}
\]
All weights are even polynomials of degree at most four in \(\mu\). In obtaining (7), if the original coefficients of \(P_\ell'\) and \(P_\ell''\) are \(B(\mu)\) and \(C(\mu)\), the boundary term is
\([C P_\ell'+(B-C')P_\ell]_{-1}^{1}\). The exact coefficients obey \(C=0\) and \(B-C'=0\) at both endpoints, so that boundary term vanishes. It is not necessary that \(B\) and \(C'\) vanish separately.

Orthogonality now makes (7) zero for every \(\ell>4\), at each radius before radial integration. For a smooth pressure mode, degrees \(\ell>2\) have zero Hessian at the origin. Thus the contact part of \(\mathcal C(-\nabla\phi)\) vanishes too. We obtain the exact selection rule
\[
 \tag{9} \boxed{\mathcal C(-\nabla\phi_\ell)=0\quad\text{for every }\ell>4.}
\]
Only degrees \(0,2,4\) can couple to this core. No assumption that \(f(r)\) is harmonic was used in (7)–(9); smoothness of the pressure at the origin is required for the pointwise contact interpretation. In the actual compact datum the initial core pressure source itself has degrees at most four, so the high-angular source also vanishes near the core and its Newton potential is harmonic there.

For smooth axisymmetric data, applying (9) to the angular tail proves the hypothesis of (6) whenever \(L\ge4\). At \(L=4\) and even parity the first omitted mode is six and \(\lambda_*=42\); for higher cutoffs use (1).

## 4. Low-mode radial and quadrature errors require separate bounds

An angular quadrature error in a retained coefficient is a low-mode error. A radial Poisson-solve residual also generally contains retained modes. Such errors cannot be assigned the \(1/\sqrt{\lambda_*}\) factor in (3) unless their actual angular support warrants it.

For comparison, any compact \(L^2\) residual \(r_{\rm low}\) supported in \(B_R\) has the safe viscosity-independent estimate
\[
 \tag{10}
 \|\nabla(-\Delta)^{-1}r_{\rm low}\|_2\le2R\|r_{\rm low}\|_2.
\]
Indeed, the three-dimensional Hardy inequality \(\|\phi/|x|\|_2\le2\|\nabla\phi\|_2\), followed by the same global energy identity, proves (10). The constant follows directly by integrating \(\operatorname{div}(x/|x|^2)=|x|^{-2}\) against \(\phi^2\) and applying Cauchy–Schwarz. This argument includes the far-field potential tail even when a monopole is present.

The corresponding **outer-part** bound is
\[
 \tag{11} |\delta\mathcal M_{\rm outer,low}|
 \le4R\|Q u_o\|_2\|r_{\rm low}\|_2.
\]
It is not a bound for the low-mode core contact and source terms. Those require a separate enclosure of the pressure modes \(0,2,4\), their radial derivatives, and the core integral. Errors in the decay/outer boundary condition may not be representable by a compact \(L^2\) residual at all; their harmonic contribution must also be tracked separately.

Thus a valid eventual certificate must combine the high-angular bound (6) with rigorous low-mode radial/source bounds, core contact/integration bounds, outer quadrature error, and the uncertainty in the specified datum and neutral tuning. The present lemma does not supply measured values of \(\|e\|_2\), \(\|Q u_o\|_2\), or any of those other enclosures. It supports a certification plan; it does not certify the current numerical sign.

`verify_angular_residual_certificate.py` checks the explicit angular weights, their boundary cancellations and degree cutoff, the contact selection rule, and the high/low prefactors. The energy and functional estimates are proved above; the checker does not replace residual-norm enclosures.

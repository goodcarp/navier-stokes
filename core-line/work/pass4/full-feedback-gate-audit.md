# Independent audit of the full second-order feedback gate

**Conclusion: the proposed gate and its constants are correct.** At a neutral affine-core insertion with \(b'=2b^2\), bending the actual normalized strain \(\beta=b/\Omega\) upward at second order requires
\[
 \boxed{\partial_t p_{zz}(0,0)<-32b^3.}
\]
Growth of one outer pressure functional does not substitute for this condition. The contact term in the full pressure derivative is essential.

## Exact equation and vanishing curvature jets

For the axisymmetric odd smooth solution, the central matrix is
\[
 L(t)=\operatorname{diag}(-b,-b,2b)+\Omega J.
\]
Define the actual central curvature quantities
\[
 B(t)=\tfrac12\partial_z\Delta u_z(t,0),\qquad
 W(t)=\tfrac12\Delta\omega_z(t,0).
\]
The full unforced NS equations give
\[
 b'=-2b^2-\tfrac12p_{zz}+\nu B,
 \qquad \Omega'=2b\Omega+\nu W,
\]
and consequently
\[
 \tag{1}
 \beta'=-\frac{p_{zz}+8b^2}{2\Omega}
              +\nu\left(\frac B\Omega-\frac{bW}{\Omega^2}\right).
\]
These are exact smooth-time identities; pressure is the pressure of the entire solution.

At insertion, the datum is exactly affine on an open core neighborhood. Hence \(\Delta u_0=0\) there and \(B(0)=W(0)=0\). Moreover
\[
 u_t(0,x)=-L^2x-\nabla p_0(x)
\]
on that neighborhood, and \(\Delta p_0=-\operatorname{tr}(L^2)\) is spatially constant there. It follows that \(\Delta u_t(0,\cdot)=0\), so \(B'(0)=W'(0)=0\) as well. This argument does not omit viscosity: its contribution and first time derivative vanish at this particular affine insertion.

Tuning \(b'(0)=2b^2\) therefore gives \(p_{zz}(0,0)=-8b^2\), \(\Omega'(0)=2b\Omega\), and \(\beta'(0)=0\). Differentiating (1) now yields
\[
 \tag{2}
 \boxed{\beta''(0)=-\frac{p_{zz}'(0,0)+32b^3}{2\Omega}.}
\]
The derivative of the denominator contributes zero because the initial numerator vanishes. The factor \(32\) is \(16bb'\) evaluated at \(b'=2b^2\). For \(\Omega>0\), the strict pressure inequality in the conclusion is equivalent to \(\beta''(0)>0\).

## The pressure derivative is a full nonlocal quantity

Let \(G_0=\nabla u_0\) and
\[
 u_{t,0}=\nu\Delta u_0-P\operatorname{div}(u_0\otimes u_0).
\]
For the decaying pressure convention \(-\Delta p=g\),
\[
 \tag{3}
 g_t(0,x)=2\operatorname{tr}\big[G_0(x)\nabla u_{t,0}(x)\big],
 \qquad p_t(0)=N*g_t(0),\quad N(x)=\frac1{4\pi|x|}.
\]
The product is a matrix trace without a transpose. Thus the scalar gate can be evaluated from the whole initial datum and its pressure projection, without first evolving the PDE. It is cubic in datum amplitudes plus viscosity times a quadratic expression. Compact support of \(\nabla u_0\) makes \(g_t(0)\) compactly supported even though the projected acceleration is nonlocal.

The required Hessian distribution gives
\[
 \tag{4}
 p_{zz}'(0,0)=
 \operatorname{PV}\int\frac{3z^2-|x|^2}{4\pi|x|^5}g_t(0,x)\,dx
 -\frac13g_t(0,0).
\]
At the tuned boundary,
\[
 \tag{5}
 g_t(0,0)=12bb'-4\Omega\Omega'
          =24b^3-8b\Omega^2.
\]
Therefore the contact contribution is exactly
\[
 -8b^3+\frac83b\Omega^2.
\]
Equivalently, the principal-value term alone must be less than \(-24b^3-(8/3)b\Omega^2\). Leaving out the contact term changes the actual gate. Neither an outer-swirl contribution nor even the full pressure contribution from a fixed outer region contains all the terms in (4).

## Why the exterior cone derivative is a different test

Let \(C_F(t)\) be the actual localized outer pressure functional, including all outer components needed for the initial decomposition, and define
\[
 \mathcal F=\frac{C_F}{\Omega^2}-\frac{38}{7}\beta^2-\frac25.
\]
For the specific compact extension, the equality \(\beta'(0)=\Omega\mathcal F(0)/2\) holds at insertion. Its unproved persistence must not be used to differentiate that equality. Introduce the explicit pressure-profile defect
\[
 \tag{6}
 D(t)=p_{zz}(t,0)+C_F(t)+\frac{18}{7}b(t)^2-\frac25\Omega(t)^2.
\]
The initial pressure calculation gives \(D(0)=0\). At later times, (1) reads exactly
\[
 \tag{7}
 \beta'=\frac\Omega2\mathcal F-\frac{D}{2\Omega}
       +\nu\left(\frac B\Omega-\frac{bW}{\Omega^2}\right).
\]
At a neutral affine insertion,
\[
 \tag{8} \boxed{\beta''(0)=\frac\Omega2\mathcal F'(0)-\frac{D'(0)}{2\Omega}.}
\]
For an explicit consistency check,
\[
 D'(0)=p_{zz}'(0,0)+C_F'(0)+\frac{72}{7}b^3-\frac85b\Omega^2.
\]
Substitution into (8) cancels \(C_F'\) and recovers (2) exactly. Thus a negative \(\mathcal F'(0)\) proves outward crossing of that particular observable cone, but does not by itself determine the sign of \(\beta''(0)\). Conversely, a positive outer-swirl derivative does not establish positive \(\beta''\). The core/profile defect and all other pressure contributions are part of the full feedback problem.

The associated checker verifies the quotient derivatives, contact coefficient, and exact cancellation between the cone and defect formulas. This is a finite initial-jet audit, not a proof of subsequent cone invariance.

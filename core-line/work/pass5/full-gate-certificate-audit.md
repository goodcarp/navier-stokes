# Audit of the complete initial Navier–Stokes feedback certificate

**Reconstruction accepted.** The complete initial scalar gate for the specified smooth compact datum contains the terms in the proposed formula and no additional core–outer, strain–swirl, or viscosity terms. Combining the exact kernel inequalities with the saved validated component enclosures gives
\[
\boxed{T=p_{zz}'(0)+32b^3<
       -\frac{290649}{8750}=-33.2170285714\ldots<0.}
\tag{1}
\]
This is a computer-assisted initial-time inequality for unforced, ordinary Navier–Stokes at viscosity \(1/1000\), and rescales to viscosity one. It is not a proof of a repeated amplification stage, global profile control, or finite-time blowup. The numerical premise is the range-integration certificates using the stated outward-rounded interval library; this audit is not a Lean formalization of that library or of the analytic identities.

## Exact datum and neutral tuning

Use the radial extension \(S_\Phi\) from the earlier pressure calculation:
\[
u_0=m+w,\quad m=S_\Phi,\quad \Phi=\psi+c\eta,\quad
w=\Omega W_\psi+Av,\qquad
b=\Omega=1,\quad c=\frac75,\quad \nu=\frac1{1000}.
\tag{2}
\]
Here \(\psi\) is the unchanged smooth core cutoff, equal to one near zero and zero for \(s\ge1\), and
\[
\eta(s)=-[1-\psi(4s/25)]
       \psi\!\left(\frac14+\frac3{44}(s-25)\right).
\]
The annulus is zero for \(R\le5/4\), equals \(-1\) on \(5/2\le R\le5\), and vanishes for \(R\ge6\). The actual exterior swirl is supported strictly inside this plateau:
\[
r\le63/200,\qquad71/20\le|z|\le89/20.
\]
The core and exterior velocities have disjoint supports; the annular strain and exterior swirl overlap intentionally. All components are smooth, compact, solenoidal, odd and axisymmetric. Near zero, \(u_0=Lx\) with
\(L=\operatorname{diag}(-b,-b,2b)+\Omega J\).

Define the amplitude **exactly** by
\[
C_v=\int Q_{\theta\theta}v_\theta^2\,dx>0,\qquad
A=\sqrt{\frac{H}{C_v}},\qquad
H=\frac{38}{7}b^2+\frac25\Omega^2=\frac{204}{35}.
\tag{3}
\]
A floating-point approximation to \(A\) would not establish exact neutrality. Equation (3) defines a finite real amplitude directly; no exact numerical quadrature for it is required.

The universal radial-strain identity applies to the sign-changing profile \(\Phi\), giving \(p[m]_{zz}(0)=-18/7\). Poloidal–azimuthal pressure cross terms vanish pointwise, and the disjoint swirls give \(p[w]_{zz}(0)=2\Omega^2/5-A^2C_v\). Thus
\[
p_0{}_{zz}(0)=-\frac{18}{7}+\frac25-\frac{204}{35}=-8.
\tag{4}
\]
This keeps the full pressure of every component.

## Exhaustive separation of the derivative

Define
\[
\Pi(a,d)=\partial_{zz}
       [N*\operatorname{tr}(\nabla a\nabla d)](0),\qquad
B(a,d)=\tfrac12\mathbb P[(a\cdot\nabla)d+(d\cdot\nabla)a].
\]
The full exact initial derivative is
\[
p_{zz}'(0)=2\Pi(u_0,\nu\Delta u_0-B(u_0,u_0)).
\tag{5}
\]
For axisymmetric fields, \(m\) is poloidal and \(w\) azimuthal. The acceleration from \(m,m\) or \(w,w\) is poloidal. Mixed transport is azimuthal and already divergence free. Poloidal–azimuthal pressure contractions vanish. Consequently (5) consists precisely of
\[
\begin{aligned}
&-2\Pi(m,B(m,m)),\\
&-2\Pi(m,B(w,w))-4\Pi(w,B(m,w)),\\
&2\nu\Pi(m,\Delta m)+2\nu\Pi(w,\Delta w).
\end{aligned}
\tag{6}
\]
These lines are respectively the complete meridional cubic, complete swirl feedback, and viscosity.

The first line is \(C_{300}[\Phi]\). It includes the core cubic, annular cubic, and both nonlocal interactions between them. It cannot be replaced by independent self terms. The pass 4 radial identity requires smoothness, compact support and \(\Phi(0)=1\), not positivity or monotonicity, so it applies to the combined profile.

In the second line, disjoint core/outer swirl supports eliminate every \(\Omega A\) term. Also \(B(S_\psi,v)=B(S_\eta,W_\psi)=0\), because their transport products vanish before Leray projection. The remaining terms are:

- The core mixed coefficient \(b\Omega^2t_{Wb}[\psi]\), with \(b=\Omega=1\).
- The response \(c\Omega^2t_{Wc}\) of the annulus to core-swirl pressure, which is exactly zero by the exterior \(R^{-3}P_2\) cancellation.
- The core response to exterior-swirl pressure and the full annular response to exterior swirl, including both local transport/backreaction terms and nonlocal pressure.

No pressure tail was discarded by disjoint support. Introduce the positive probability measure
\[
d\mathfrak m(x)=\frac{Q_{\theta\theta}(x)v_\theta(x)^2}{C_v}\,dx.
\]
The last item is exactly
\[
A^2C_v\int[k_{\rm core}+c\,k_{\rm pump}]\,d\mathfrak m.
\tag{7}
\]
Here \(k_{\rm core}\) is normalized by \(b=1\), and \(k_{\rm pump}\) is the complete local-plus-pressure annular kernel. Counting only its local part or adding its pressure twice would be incorrect.

For viscosity, \(\Delta S_\Phi=S_{4s\Phi''+14\Phi'}\), and the latter profile vanishes near zero. Polarizing the universal central-pressure identity gives \(\Pi(m,\Delta m)=0\). The radial core swirl likewise has \(2\Pi(W_\psi,\Delta W_\psi)=0\). Differentiation preserves the two swirl supports, so there is no viscous cross term. The complete viscous contribution is therefore
\[
\nu A^2d_v=\nu H\,\frac{d_v}{C_v}.
\tag{8}
\]
Equations (6)–(8) prove the full identity
\[
\boxed{
T=32+C_{300}[\Phi]+t_{Wb}[\psi]
 +H\int[k_{\rm core}+c\,k_{\rm pump}]\,d\mathfrak m
 +\nu H\,\frac{d_v}{C_v}.}
\tag{9}
\]

## Component bounds and exact aggregation

| Term | Bound | Evidence |
|---|---|---|
| \(C_{300}[\Phi]\) | \(<66\) | radial-strain-interval-1024.json; radial identity and streaming-radial-audit.md |
| \(t_{Wb}[\psi]\) | \(<4\) | core-twb-interval.json; actual enclosure is between \(2.5300\ldots\) and \(2.5971\ldots\) |
| \(k_{\rm core}\) | \(<-12/5\) | local-swirl-feedback-kernel.md; stronger exact upper bound \(-87238/35287\) |
| \(k_{\rm pump}\) | \(<-31/2\) | swirl-adjoint-kernel.md, Section 6; exact support/join envelope |
| \(d_v/C_v\) | \(<901\) | exact dyadic endpoints in local-interval-results.json |

The aggregate checker reads exact saved dyadic endpoints, rather than printed decimal approximations. Its rational comparisons establish \(C_v>0\), \(d_v<901C_v\), \(C_{300}<66\), and \(t_{Wb}<4\). The kernel estimates use exact rational support and join bounds. Their inequality directions are preserved because \(H,c,\nu>0\):
\[
\begin{aligned}
T
&<32+66+4+\frac{204}{35}
       \left[-\frac{12}{5}
             -\frac75\frac{31}{2}+\frac{901}{1000}\right]\\
&=-\frac{290649}{8750}<0.
\end{aligned}
\tag{10}
\]
The exact decimal is \(-33.2170285714\ldots\); the suggested \(-33.216\) was an arithmetic approximation. This margin does not rely on convergence of a truncated pressure computation.

## Actual solution and initial derivative conditions

Smooth compact divergence-free data have a local smooth solution. Uniqueness preserves odd axisymmetry. Define
\[
b(t)=\tfrac12\partial_z u_z(t,0),\quad
\Omega(t)=\tfrac12\omega_z(t,0),\quad\beta(t)=b(t)/\Omega(t).
\]
The exact central equations are
\[
b'=-2b^2-\tfrac12p_{zz}+\nu B,\qquad
\Omega'=2b\Omega+\nu W,
\]
where \(B=\partial_z\Delta u_z(t,0)/2\) and \(W=\Delta\omega_z(t,0)/2\).
Initially \(\Delta u_0=0\) near zero. Also
\[
u_t(0)=-L^2x-\nabla p_0,\qquad
\Delta p_0=-\operatorname{tr}(L^2)
\]
on that neighborhood, so \(\Delta u_t(0)=0\) there too. Hence
\[
B(0)=W(0)=B'(0)=W'(0)=0.
\]
These are initial identities; viscosity is not assumed absent later. Neutrality gives \(b'(0)=\Omega'(0)=2\) and \(\beta'(0)=0\). Differentiating the exact ratio equation, with its viscous terms, yields
\[
\boxed{\beta''(0)=-\frac{T}{2\Omega}
  >\frac{290649}{17500}=16.6085142857\ldots>0.}
\tag{11}
\]
Thus the actual solution initially bends toward a larger central strain-to-rotation ratio. Smoothness implies \(\beta(t)>\beta(0)\) on some sufficiently short positive interval. This certificate supplies no quantitative duration, invariant profile class, or return map.

## Ordinary viscosity one

Let \(u,p\) be the solution at \(\nu=1/1000\), and set \(\lambda=1000\):
\[
U(t,x)=\lambda u(\lambda t,x),\qquad
P(t,x)=\lambda^2p(\lambda t,x).
\tag{12}
\]
Direct substitution gives
\[
U_t+(U\cdot\nabla)U=-\nabla P+(\lambda\nu)\Delta U
                  =-\nabla P+\Delta U.
\]
The transformed datum remains smooth, compact and solenoidal, and there is no force. This converts viscosity by amplitude/time scaling; it is not the usual fixed-viscosity spatial symmetry.

The central amplitudes scale by \(\lambda\), the pressure Hessian by \(\lambda^2\), and the full gate by \(\lambda^3\):
\[
T_{\nu=1}=\lambda^3T_{\nu=1/1000}<0.
\]
Neutrality is preserved after this scaling, and \(\beta''\) scales by \(\lambda^2\). The certificate therefore supplies a favorable actual initial gate at viscosity one too, without claiming a later singularity or enlarging the certified time interval.

The file verify_full_gate_bound.py checks exact endpoint comparisons, rational aggregation, neutral tuning and scaling factors. The analytic decomposition and component range proofs remain its explicit mathematical dependencies.

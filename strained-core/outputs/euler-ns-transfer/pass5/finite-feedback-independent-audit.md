# Independent audit of the finite-feedback corollary

**Accepted without a substantive correction.** The explicit time restrictions in finite-feedback-corollary.md imply all of its central, fixed-region reservoir, modal, weighted-energy, and RMS inequalities. The order-ten Sobolev estimates provide the required derivatives of the actual solution and retain pressure and viscosity. The initial certificate is unchanged.

This audit checked the note against the pass 3 Sobolev argument, the pass 4 strained-receiver identity, and the pass 5 full-gate certificate. The author of the subsequent return-map note confirmed that no concurrent edits to these formulas were underway.

## Initial jets and the positive-time observables

At the certified neutral datum \(b=\Omega=1\), \(b'=\Omega'=2\). The local affine identities imply \(W(0)=W'(0)=0\), where \(W=\Delta\omega_z/2\). Therefore
\[
\Omega''(0)=2b'(0)\Omega(0)+2b(0)\Omega'(0)=8.
\]
The quotient formula gives \(\beta''(0)=b''(0)-\Omega''(0)\). For \(S=b'-2b^2\), this proves
\[
S(0)=0,\qquad S'(0)=b''(0)-8=\beta''(0)\ge\delta.
\]
The corollary correctly treats \(S>0\) and \(\beta'>0\) as distinct inequalities at later positive viscosity; it does not discard the subsequent curvature terms.

The swirl-stress ratio uses a fixed cutoff \(\zeta\), zero near the core and identically one on a neighborhood of the initial outer swirl support. The initial value and first derivative are therefore unaffected by cutoff derivatives. Its logarithmic lower bound is arithmetically exact:
\[
\frac{2381}{357}\frac75-4-\frac{901}{1000}
=4+\frac{22249}{51000}.
\]
This fixed-region swirl channel is not the total exterior pressure. The corollary retains that distinction.

## Derivative order and quotient constants

The pass 3 energy and product estimates are general estimates for smooth divergence-free solutions. Their use with the new datum is valid after replacing the initial norm by this datum's \(X_0\). They give
\[
\|u\|_{10}\le V,\quad \|u_t\|_8\le M_1,\quad
\|u_{tt}\|_6\le M_2,\quad \|u_{ttt}\|_4\le M_3
\]
on the stated \(T_E\), with the exact constants recorded in the corollary. In particular, the matrix-gradient embedding from pass 3 gives
\[
\|\nabla u\|_\infty\le V,\qquad
\|\nabla\partial_t^j u\|_\infty\le M_j,\quad j=1,2,3.
\]
Here the matrix norm can be the Frobenius norm; it controls the operator norm. The vector/matrix constants were explicitly accounted for in pass 3, so no additional component factor is missing. Central \(b^{(j)}\) and \(\Omega^{(j)}\) are consequently bounded by \(M_j\).

For \(t\le1/(2M_1)\), \(|b|\le3/2\) and \(\Omega\ge1/2\). The full third quotient derivative is
\[
\begin{aligned}
\beta'''={}&\frac{b'''}{\Omega}
-\frac{3b''\Omega'+3b'\Omega''+b\Omega'''}{\Omega^2}\\
&+\frac{6b'(\Omega')^2+6b\Omega'\Omega''}{\Omega^3}
-\frac{6b(\Omega')^3}{\Omega^4}.
\end{aligned}
\]
Bounding each term gives precisely
\[
|\beta'''|\le8M_3+96M_1M_2+192M_1^3.
\]
Likewise \(S''=b'''-4(b')^2-4bb''\) gives
\[
|S''|\le M_3+4M_1^2+6M_2.
\]
Thus the displayed \(B_\beta\) and \(B_S\) require no fourth time derivative and no Sobolev order beyond the supplied ladder.

For the reservoir, differentiating the fixed quadratic functional only in time gives
\[
|C'|\le2K_\zeta VM_1=N_1,\qquad
|C''|\le2K_\zeta(M_1^2+VM_2)=N_2.
\]
The azimuthal projection is only used as an \(L^2\) contraction, so no false assertion about its spatial derivatives near the axis enters. On \(C\ge H_0/2\), \(\Omega\ge1/2\),
\[
\left|\left(\frac{\mathcal R'}{\mathcal R}\right)'\right|
\le\frac{2N_2}{H_0}+\frac{4N_1^2}{H_0^2}+4M_2+8M_1^2.
\]
This is exactly the stated \(B_{\mathcal R}\).

## Every time restriction has the required direction

All bounds below hold on the common Sobolev lifespan and after the preceding denominator restrictions:

| Restriction | Consequence |
|---|---|
| \(t\le\delta/(2B_\beta)\) | \(\beta''(t)\ge\delta/2\), hence \(\beta'(t)\ge\delta t/2\) and \(\beta(t)\ge1+\delta t^2/4\) |
| \(t\le\delta/(2B_S)\) | \(S(t)\ge\delta t-B_St^2/2\ge3\delta t/4\), stronger than the claimed \(\delta t/2\) |
| \(t\le2/(3M_2)\) | \(\Omega'(t)\ge2-M_2t\ge4/3>1\); it also supplies the modal bound below |
| \(t\le H_0/(2N_1)\) | \(C(t)\ge H_0/2\) |
| \(t\le\varepsilon/(2B_{\mathcal R})\) | \(\mathcal R'/\mathcal R\ge4+\varepsilon/2>4\) |
| \(t\le8/[3(M_1^2+VM_2)]\) | \(E_r'(t)\ge8I_r>0\), uniformly in receiver radius |

The note's \(t\le1/(2M_1)\) restriction supplies the needed lower bound on \(\Omega\) throughout. No restriction is circular: its constants depend on the fixed datum and fixed cutoff, not on a norm evaluated at an unknown future solution.

No numerical enclosure of the new \(X_0\) has been supplied. Therefore the norm formula specifies a positive mathematical interval but does not certify a decimal duration or a useful-size gain. The corollary expressly states this limitation. Its qualitative positive interval requires only smoothness and the strictly favorable initial margins.

## Uniform receiver jets, energy and Reynolds comparisons

The condition \(\Delta u_t(0)=0\) holds throughout the initially affine ball, not merely at the center. Each component of \(u_t(0,\cdot)\) is harmonic there. The radial harmonic first-moment identity consequently gives the exact receiving derivatives
\[
\Omega_r'(0)=b_r'(0)=2
\]
from the actual central matrix derivative \(2(D+J)\), even though \(u_t\) may contain higher harmonic components.

Oddness is preserved for every time derivative. Hence
\[
|u_{tt}(t,x)|\le M_2|x|.
\]
Using \(\int\chi_r|x|^2=3I_r\), \(|Jx|\le|x|\) and \(|Dx|\le2|x|\) gives
\[
|\Omega_r''|\le\frac32M_2,\qquad |b_r''|\le M_2.
\]
There is no inverse receiver radius in these estimates. The restriction \(t\le2/(3M_2)\) yields both modal derivatives at least one.

The two weighted modes \(Jx,Dx\) are orthogonal, with squared norms \(2I_r,6I_r\). Weighted orthogonal projection therefore proves
\[
E_r(t)\ge I_r[\Omega_r(t)^2+3b_r(t)^2],
\qquad E_r(0)=4I_r.
\]
This establishes the stated energy and RMS lower gains. They concern the full evolved velocity, not merely its two projections.

Also \(E_r'(0)=16I_r\). Oddness and the bounds on \(u,u_t,u_{tt}\) give
\[
|E_r''|\le3I_r(M_1^2+VM_2),
\]
so the common interval gives \(E_r'>0\), and \(U_r'>0\) since the energy and receiver mass are positive.

Initially the affine RMS scales as \(U_r(0)=c_\chi r\). Thus for the chosen
\(q^2=(1+t/2)/(1+t)\), the uniform receiver estimate gives exactly
\[
\frac{qr\,U_{qr}(t)}{r\,U_r(0)}
\ge q^2(1+t)=1+t/2.
\]
The analogous \(r^2\Omega_r/\nu\) comparison follows identically. After choosing the endpoint time, the smaller receiver is fixed in space. These statements neither follow a material loop nor assert that vorticity is confined to the smaller ball.

## Scope and viscosity conversion

The finite corollary is a consequence for one actual solution. It does not supply a family of inherited states, an invariant velocity profile, a quantitative return time, or repeatability. The strict receiving gains are compatible with the new higher spatial jets discussed in the separate return-map note.

If this finite episode is transferred from viscosity \(1/1000\) to viscosity one by \(U(t,x)=1000u(1000t,x)\), its physical time interval becomes \([0,\tau/1000]\). Relative gains at corresponding times remain the same; one must not retain the unscaled physical duration or write \(1+t\) with the new physical time without replacing \(t\) by \(1000t\).

No change to finite-feedback-corollary.md or the initial certificate was necessary.

## Independent check of the next return calculation

The proposed cubic-jet defect is consistent with this finite corollary and has the stated coefficient and sign. Near the affine core, the radial meridional pressure has constant solid-harmonic degree-four coefficient
\[
h_4(0)=\frac{32}{21}\int_0^\infty(\Phi')^2\,ds.
\]
Since \(\partial_z^4H_4(0)=24\), its fourth axial derivative is
\[
p[m]_{zzzz}(0)=\frac{256}{7}\int_0^\infty(\Phi')^2\,ds.
\]
The degree-zero and degree-two radial coefficients are constant/quadratic there and add no fourth derivative. The radial core swirl also has only constant/quadratic pressure in this neighborhood; its fourth axial derivative is exactly zero.

The core transition alone obeys
\[
\int_0^\infty(\Phi')^2\,ds
\ge\int_{1/4}^1(\psi')^2\,ds
\ge\frac{[\psi(1)-\psi(1/4)]^2}{1-1/4}=\frac43.
\]
The derivative-free outer-swirl kernel gives
\[
A^2p[v]_{zzzz}(0)
=-H_0\int \frac{G(r^2/z^2)}{z^2}\,d\mathfrak m
\ge-H_0\frac{30}{(71/20)^2}.
\]
Combining these terms yields exactly
\[
\boxed{p_0{}_{zzzz}(0)\ge
\frac{256}{7}\frac43-\frac{204}{35}\frac{30}{(71/20)^2}
=\frac{3693184}{105861}>34.}
\]
The positive contribution from the annular joins is already included in \(\int(\Phi')^2\); discarding it here is a valid lower bound.

Initially \(u_0=Lx\), and its local convective acceleration is linear while its Laplacian vanishes. Therefore
\[
\partial_z^3u_{t,z}(0,0)=-p_0{}_{zzzz}(0),\qquad
\partial_z^3u_z(t,0)=-p_0{}_{zzzz}(0)t+O(t^2).
\]
The sign is negative. The \(O(t^2)\) remainder is legitimate from the same derivative ladder: \(u_{tt}\in H^6\) controls three spatial derivatives pointwise.

For the standard RMS-selected normalization
\(V_t(y)=q(t)^{1+\alpha}u(t,q(t)y)\), the scalar selection equation has a nonzero radius derivative \(\alpha+2\) at \((t,q)=(0,1)\). The implicit-function branch consequently satisfies \(q'(0)=-2/(\alpha+2)\). Its prefactor cancels the first-order change in the central matrix:
\[
\nabla V_t(0)=L+O(t^2).
\]
It cannot cancel the initially zero third spatial jet:
\[
\partial_{y_z}^3(V_t)_z(0)
=q(t)^{\alpha+4}\partial_z^3u_z(t,0)
=-p_0{}_{zzzz}(0)t+O(t^2).
\]
Thus a centered amplitude/radius normalization does not return the velocity to an exactly affine core even at first order in the elapsed time. This is a concrete defect that a stronger state class or return construction must accommodate. It does not refute returns allowing that defect, nor does it contradict any fixed-receiver gain in the finite corollary.

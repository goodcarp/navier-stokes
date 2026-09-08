# Independent audit of the full MAC diagnostic

Reviewed the current parity-corrected, axial-padding version of `evolve_full_mac.py`. `check_full_mac.py` completed with PASS and saved `full-mac-structural-checks.json`. The JSON pins the exact solver SHA256. Tests used 16×32, 32×64, and 64×128 radial/axial grids with modes 0,4,8 and the new `padding` setting; a separate 16×32 test retains the legacy `filter` setting. No production solver file was edited by this audit.

**Verdict:** the kinetic adjoints, orthogonal pressure projection, positive viscous form, constrained Crank–Nicolson signs, ±k factor reuse, and padded Lamb nonlinearity are consistent. The manufactured curl errors decrease under refinement. These results validate discrete identities and limited consistency tests; they do not provide a smooth whole-space reconstruction, a rigorous residual enclosure, or an endpoint theorem.

## Projection and kinetic metric

Let cell weights be \(V_i\), interior radial-face weights \(W_i=r_i\Delta r_{\mathrm{centers},i}\), and velocity mass \(M=\operatorname{diag}(W,V,V)\). With the exterior/axis radial-face degrees of freedom removed, summation by parts gives

\[
G_rp=\frac{p_{\rm right}-p_{\rm left}}{\Delta r_{\mathrm{centers}}},
\quad G_\theta p=im p/r,\quad G_zp=ikp,
\qquad G=-M^{-1}D^*V.
\]

Thus \(-DG\) is positive semidefinite in the cell-pressure metric. Solving \(-DGq=Da\) and returning \(a+Gq\) is the kinetic-orthogonal projection onto \(\ker D\). Axial frequency truncation commutes with this modewise projection and is itself orthogonal.

The zero pressure mode \(m=k=0\) has a constant-pressure gauge. The omitted final equation is redundant because \(\sum_iV_i(Du)_i=0\). Every divergence-free radial-face velocity in that mode is zero; explicitly zeroing it after projection is consistent with the projection rather than an additional physical constraint on another mode. Tests confirmed idempotence, orthogonality, contraction, and weighted self-adjointness. Relative defects were below 1.6e-15; divergence after projection was at relative roundoff compared with the input divergence.

## The viscous angular cross term has the correct sign

For a radial vector \(a=(a_r,a_\theta,a_z)\), set \(b=Ia_r\). The independently expanded radial/angular quadratic form is

\[
\sum_iV_i|(D_fa_r)_i|^2
+\sum_fW_{g,f}(|(Ea_\theta)_f|^2+|(Ea_z)_f|^2)
\]
\[
+\sum_i\frac{V_i}{r_i^2}
\left(|im b_i-a_{\theta,i}|^2+|im a_{\theta,i}+b_i|^2
 +m^2|a_{z,i}|^2\right). \tag{1}
\]

The off-diagonal block is therefore \(K_{r\theta}=2imI^TH\), as implemented, with its conjugate transpose in the other block. The angular two-component eigenvalues before gathering are \((m-1)^2\) and \((m+1)^2\); the form is nonnegative. Axial differentiation adds \(k^2M\). The independent sum-of-squares evaluation agreed with the assembled matrix to below 2.4e-16, and all tested generalized eigenvalues were positive.

The outer wall term in \(E\) enforces the stated zero tangential/axial boundary value in this gradient form. It is a no-slip finite-cylinder model, not an exterior harmonic boundary operator. The radial component is fixed to zero at the outer boundary separately. These choices must not be described as a whole-space pressure or viscous boundary condition.

## Constrained Crank–Nicolson and reuse at negative axial frequency

With \(\alpha=h\nu/2\), the saddle system has velocity block
\(A=M+\alpha(K+k^2M)\) and pressure block \(-D^*V=MG\). Its constraint is \(Du_{\rm new}=0\). For divergence-free input and exact solution, the discrete kinetic energy satisfies

\[
E(u_{\rm new})-E(u_{\rm old})
=-h\nu\,\mathcal D((u_{\rm new}+u_{\rm old})/2). \tag{2}
\]

This is the identity checked in the program; no endpoint trapezoidal dissipation approximation was substituted. Its numerical defect was below 6e-17 on the padding runs, and energy decreased in every tested viscous step.

For fixed angular mode, let \(S_z\) reverse only the axial velocity component. Then \(D_{-k}=D_kS_z\), the viscous/mass block commutes with \(S_z\), and the full negative-frequency saddle matrix is obtained by conjugation with \(\operatorname{diag}(S_z,I_p)\). Flipping the axial right-hand side before the positive-k solve and flipping the resulting axial velocity afterward is therefore correct; neither pressure nor the angular component is conjugated by this operation. Independent dense solves using the actual signed k, including the highest retained axial frequency on the small grid, agreed within 5.2e-15. The special zero-mode gauge solve also agreed.

The CN map is contractive, but this does not give strong damping or accurate dynamics for arbitrarily unresolved stiff modes. The complete CN/RK4/CN composition still needs time-step convergence and an error estimate. The nonlinear RK4 substep is not exactly energy preserving merely because its differential equation has zero kinetic work.

## Lamb nonlinearity, padding, and symmetry

The implemented scatter is the exact kinetic adjoint of the face-to-cell gather:
\(\langle a,SF\rangle_M=\langle Ia_r,F_r\rangle_V+
\langle a_\theta,F_\theta\rangle_V+\langle a_z,F_z\rangle_V\).
Therefore the cellwise identity \(u\cdot(u\times\omega)=0\) gives zero semidiscrete nonlinear work after angular Galerkin projection. Its sign is the correct NS Lamb sign, since \(-\Pi(u\cdot\nabla u)=\Pi(u\times\omega)\). The compatible pressure projection preserves this zero work for a discrete divergence-free state. The tested relative raw and projected work was below 4.3e-18 for the padding runs.

The new padding setting retains all axial integers \(|k|<N_z/2\), omitting the even-grid Nyquist coefficient, and evaluates products on at least \(3N_z/2\) points. This satisfies the strict non-aliasing inequality for the retained quadratic Galerkin product. The angular grid similarly satisfies \(N_\phi>3J\) for angular indices \(-J,\ldots,J\). The padded curl correctly recomputes axial frequencies from the padded array length. Fourier resampling down to the original grid is the adjoint of the upsampling on this retained state space in the physical \(L/N_z\) norm; the possible destination Nyquist merge cannot affect kinetic work because that state coefficient is zero and is removed by projection.

The independent reference comparison used a 12×32 grid with a 48-point axial product grid, against a 12×96 grid retaining twice as many angular modes. About 30.95% of the test input energy lies above the old \(N_z/3\) filter cutoff, so this genuinely exercises the newly retained frequencies. The largest retained axial integer is 15. The resampling-adjoint defect was 5.2e-18 and the retained nonlinear output agreed with the larger reference to 9.3e-16. The legacy filtering path was tested separately and continues to satisfy its own discrete identities.

Arbitrary angular phase rotation, axial grid translation, real mode-zero fields, and the Cartesian odd symmetry of the C4 class were preserved to roundoff by the tested projection/nonlinearity/viscous maps. A single small split step also retained discrete divergence and decreased energy in each tested case. This is a bounded path check, not a stability theorem for the proposed research interval.

## Axis consistency and remaining scope

The corrected curl now distinguishes odd radial dependence of \(u_\theta\) from even radial dependence of \(u_z\), and includes the outer wall derivative. The exact rigid-core field \(u_\theta=r\) has computed vertical vorticity two in the first four cells to within 2.3e-16. Its incompatible outer wall is deliberately excluded from that local check.

Smooth, wall-compatible toroidal manufactured fields produced the following relative weighted L2 curl errors:

| Radial cells | Mode 0 | Mode 4 |
|---:|---:|---:|
| 16 | 0.03651 | 0.09390 |
| 32 | 0.00936 | 0.04298 |
| 64 | 0.00237 | 0.01286 |

These are observed decreases, not a proved convergence rate. The fields were sampled analytically before projection, so the test isolates curl consistency instead of conflating it with an initial projection error.

Parity-aware derivatives and zero discrete divergence still do not enforce the full Cartesian regularity of every evolved mode. Mode four requires \(u_r+iu_\theta=O(r^5)\), \(u_r-iu_\theta=O(r^3)\), and \(u_z=O(r^4)\), with the corresponding smooth even radial factors. A potential-based reconstruction must impose and bound these conditions. The finite radial wall, periodic axial copies, initial representation changes, high angular modes, radial/axial truncation, and full-space residual/tails remain separate obligations. The audit found no remaining structural correction needed in the pinned version; it does not remove those continuum obligations.

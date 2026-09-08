# Independent solver strategy and shared-frequency projection reference

The compatible MAC calculation is a practical next diagnostic if its kinetic projection, viscous stability and axis convergence checks pass. The independent reference supplied here is a common-radial-frequency Fourier–Hankel projector. It checks the continuous operator and catches a concrete failure of mixing different Bessel-zero grids. For a later smooth whole-space approximate trajectory, compact, axis-factored poloidal/toroidal Galerkin functions give a more direct finite-energy representation with sparse radial linear algebra. Neither representation by itself supplies the nonlinear residual certificate or a return stage.

Companion files are `shared_hankel_projection.py`, `check_shared_hankel_projection.py`, and the generated `shared-hankel-projection-check.json`. The independently coordinated potential formulas, unchanged initial datum and full-space residual requirements are in `poloidal-toroidal-evolution-design.md`.

## The common-frequency operator

Use angular mode \(e^{im\theta}\), axial mode \(e^{i\kappa z}\), and helical components \(U_\pm=u_r\pm i u_\theta\). Transform them at the **same radial frequency** \(k\) with signed Bessel orders \(m+1,m-1,m\), respectively. Our convention is

\[
\widehat f_n(k)=\int_0^\infty rJ_n(kr)f(r)\,dr,
\qquad f(r)=\int_0^\infty kJ_n(kr)\widehat f_n(k)\,dk.
\]

This is the weighted version of the [DLMF Hankel transform and inverse](https://dlmf.nist.gov/10.22.E76). Applying the [Bessel ladder identities](https://dlmf.nist.gov/10.6.E6) gives

\[
G=(-k,k,i\kappa)^T,\quad D=(k/2,-k/2,i\kappa),\quad
DG=-(k^2+\kappa^2),\qquad
\Pi=I+\frac{GD}{k^2+\kappa^2}.
\]

In the velocity metric \(H=\operatorname{diag}(1/2,1/2,1)\), \(D=-G^*H\). Consequently \(D\Pi=0\), \(\Pi^2=\Pi\), and \(\Pi^*H=H\Pi\). Diffusion is the scalar multiplier \(-(k^2+\kappa^2)\), so it commutes with this continuous projector. These are exact algebraic identities, also checked symbolically.

If absolute transform orders are used, mode zero has orders \((1,1,0)\) and instead needs
\(G=(-k,-k,i\kappa)^T\), \(D=(k/2,k/2,i\kappa)\). The checker tests both conventions. Physical reality swaps the two helical components: \(U_{+,-m,-\kappa}=\overline{U_{-,m,\kappa}}\). For signed Hankel amplitudes this also contributes the Bessel sign: \(\widehat U_{+,-m,-\kappa}=(-1)^{m-1}\overline{\widehat U_{-,m,\kappa}}\), and \(\widehat U_{z,-m,-\kappa}=(-1)^m\overline{\widehat U_{z,m,\kappa}}\). Positive-mode storage must respect these relations.

Sampling each component at the zeros of its own Bessel order and applying the displayed projector indexwise is incorrect: equal indices are different \(k\)'s. For one exact Gaussian gradient at \(m=4\), the checker measures **9.1986% spurious surviving gradient** under that deliberately incompatible operation. This does not exclude compatible non-diagonal transforms between grids; it excludes pretending their frequency indices coincide.

## What the manufactured checks establish

The tests use \(W=\nabla\times(\psi e_z)+\nabla\phi\), with \(\phi=r^m e^{-0.8r^2}e^{i\kappa z}\), \(\psi=r^m e^{-1.1r^2}e^{i\kappa z}\), and \(\kappa=1.3\). The exact Gaussian Hankel pair follows from [DLMF 10.22.51](https://dlmf.nist.gov/10.22.E51). At modes 0, 4 and 8, 256-point Gauss quadrature in each of \(r\in[0,12]\), \(k\in[0,20]\) gives:

- exact-transform and gradient-removal relative errors below \(8\cdot10^{-15}\);
- physical-space recovery and projected roundtrip errors below \(4\cdot10^{-14}\);
- directly differentiated physical divergence below \(5\cdot10^{-15}\);
- Parseval error below \(2\cdot10^{-14}\).

An independent random-spectrum test checks kinetic orthogonality and energy partition. The checker also verifies the potential mass/stiffness symbols below. These are floating-point manufactured checks plus exact symbolic identities, not interval estimates. They do not establish resolution of the much sharper actual cutoff datum. Algebraic projection at the nodes cannot substitute for a physical transform roundtrip and physical derivative checks on that datum.

A finite quadrature sum of isolated Bessel waves is not a whole-space L2 velocity. For mathematical reconstruction, one can represent integrable spectral amplitude functions over bins and account for their integration errors, or reconstruct compact potentials. A finite periodic axial Fourier series likewise has infinite energy when repeated throughout the whole axial line. These qualifications do not prevent using either representation on a diagnostic domain, but they determine which problem it approximates.

## A sparse potential implementation

Write \(u=\nabla\times(Te_z)+\nabla\times\nabla\times(Pe_z)\), with scalar mode bases \(B_i(r)=r^{|m|}b_i(r^2)\). The regular factored basis builds smoothness at the axis into the velocity. Near the axis, evaluate the factored derivative identities in the companion design; avoid cancellation of separately singular terms.

Let \(L_m=\partial_{rr}+r^{-1}\partial_r-m^2/r^2\). For real compact radial basis functions, set

\[
A_{ij}=\int_0^\infty r\left(B_i'B_j'+\frac{m^2}{r^2}B_iB_j\right)dr,
\qquad B_{ij}=\int_0^\infty r(L_mB_i)(L_mB_j)dr,
\]
\[
C_{ij}=\int_0^\infty r\left[(L_mB_i)'(L_mB_j)'+\frac{m^2}{r^2}(L_mB_i)(L_mB_j)\right]dr.
\]

Up to the common angular/axial normalization, the exact blocks at axial frequency \(\kappa\) are

\[
M_T=A,\quad K_T=B+\kappa^2A,
\qquad M_P=B+\kappa^2A,\quad K_P=C+2\kappa^2B+\kappa^4A.
\]

The toroidal and poloidal sectors are orthogonal in both mass and stiffness. One quick verification uses their common-frequency velocity columns
\(t=(ik,ik,0)^T\), \(p=(-i\kappa k,i\kappa k,k^2)^T\):
\(t^*Hp=0\), \(t^*Ht=k^2\), \(p^*Hp=k^2(k^2+\kappa^2)\). Multiplication by \(k^2+\kappa^2\), Parseval and integration by parts give the displayed blocks. Compact outer support and axis regularity eliminate boundary terms, including for \(m=0\).

Local B-spline supports make these matrices banded. Their real factors can be shared between \(\kappa\) and \(-\kappa\). Remove any potential gauge representing zero velocity and use the velocity Gram norm; unweighted potential coefficients can be badly conditioned near zero frequency. In the axisymmetric zero axial mode, represent the correct zero-net-axial-flux compact class rather than silently introducing an incompatible constant axial velocity.

For an H4 residual involving \(\Delta u\), a sufficient choice is \(P\in H^8\), \(T\in H^7\), with the corresponding time regularity. Degree-nine simple-knot splines with sufficiently smooth zero extensions are one conservative common choice. Ordinary cubic splines do not meet this strong target. Knot spacing can resolve the narrow inner/packet transitions and remain coarser on broad pump annuli. A Galerkin evolution then solves the banded weak projected equations and computes the quadratic term on a dealiased angular/axial grid. Exact integration gives the weak energy identity; numerical quadrature and omitted angular modes must still be charged.

This compact-potential choice imposes an approximation space, not a physical impermeable wall. The true projected acceleration extends beyond compact velocity support through pressure. The resulting full residual includes that exterior pressure gradient. Any conversion from periodic MAC samples should put radial/axial windows **inside** the potentials and explicitly include the changed velocity, join residual, and pressure tail. The companion note gives a tail integral formula.

## Cost and practical choice

Dense common-grid Hankel quadrature needs roughly \(O(N_mN_zN_rN_k)\) transform work and a time-bandwidth budget growing with \(k_{\max}R\). A narrow radial transition and the distant pump annulus therefore make a uniformly fine global transform expensive even when the projection itself is diagonal. The following are memory illustrations, not resolution claims:

- Nine nonnegative C4 angular modes, \(N_k=512\), \(N_z=1024\): one three-component complex128 state is 216 MiB; three such buffers use 648 MiB.
- With \(N_r=768\), the 26 distinct absolute radial orders for those modes require about 78 MiB for one real transform-matrix cache. Derivative caches add comparable storage.
- A three-component real product grid with \(N_r=768\), padded \(N_z=1536\), and 32 C4 angular samples occupies 864 MiB. Streaming radial slabs avoids retaining the full product grid.
- Two potential coefficients on 256 local radial functions at the same nine angular modes and 1024 axial modes occupy 72 MiB per state. Banded linear factors and derivative workspaces add storage, but locality avoids a dense global radial transform.

Thus keep the common-frequency module as the free-space algebra and reconstruction benchmark, and use a compatible MAC run for the immediate trajectory diagnostic. Develop the local potential discretization if the trajectory survives convergence checks and merits whole-space reconstruction. No new basis choice can make an unresolved narrow initial profile reliable.

## Checks required of the immediate MAC evolution

For kinetic mass matrix \(M_h\), check \(D_h=-G_h^*M_h\) with the appropriate scalar pressure weights, and \(\Pi_h\)'s kinetic orthogonality and idempotence. Check the implemented vector diffusion has \(\langle u,L_hu\rangle_{M_h}\le0\). Then for \(u\) in the projected space,
\(\langle u,\Pi_hL_hu\rangle_{M_h}=\langle u,L_hu\rangle_{M_h}\); exact divergence–diffusion commutation is not required for this discrete dissipation statement.

Explicit RK4 must resolve the largest actual viscous eigenvalue, including \(m/r\) terms near the axis. Its negative-real stability endpoint is approximately 2.785; a bulk \(\nu\Delta t/h^2\) rule that ignores the angular/axis operator can miss the restriction. Transport and nonnormal effects require their own stability/step checks. A kinetic-adjoint pressure pair does not automatically make the nonlinear transport energy-conserving; monitor its discrete energy work separately from pressure and viscosity, with a compatible skew formulation if needed.

At mode four, smooth Cartesian velocity requires \(U_+=O(r^5)\), \(U_-=O(r^3)\), \(u_z=O(r^4)\), with the even smooth radial factors. Small discrete divergence or a single axis parity condition alone does not imply these orders. Check the normalized near-axis coefficients as resolution changes, using manufactured regular helical modes before interpreting a central derivative.

The impermeable outer wall changes pressure immediately. Even when the initial velocity and nonlinearity vanish near it, the wall projection imposes zero normal acceleration and hence the relevant homogeneous pressure Neumann condition there. The free-space or open periodic-cylinder pressure generally has a nonzero radial derivative at that radius. Therefore compare the initial pressure against a **matched wall** reference, or separate and estimate the wall correction. The independently computed axial image correction in pass10 isolates axial periodicity only; it is not a wall correction.

The first meaningful actual-datum comparisons are initial field interpolation, kinetic energy, fluctuation energy derivative, and central pressure curvature, followed by timestep and space refinement over a short interval. The Gaussian projector test is an independent operator control for these comparisons. A converged finite-cylinder trajectory still needs a continuous solenoidal reconstruction and a full-space residual bound before it can enter the proposed mathematical error bridge.

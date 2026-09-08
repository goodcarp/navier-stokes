# Implementable MAC blocks with a harmonic exterior

The finite-radius upgrade can retain the staggered interior and add one outer radial velocity degree of freedom per nonzero angular/axial mode. Its mass, boundary traces, pressure projection, viscous Gram and nonlinear weak scatter must change together. This note specifies that discrete model. `check_exterior_mac_blocks.py` independently assembles the blocks and checks their identities; it does not edit or call the production time integrator.

The exact continuum exterior mass and stiffness used below were derived in `exterior-harmonic-mass-audit.md`. The interior quadrature and finite differences are discrete choices. Positive energy and divergence identities for those choices do not establish their continuum consistency order, actual axis resolution, or an NS error certificate.

## 1. Coordinates, mode and state

Use a C4 angular mode \(m=0,4,8,\ldots\), axial frequency \(k\), radial faces \(0=\rho_0<\cdots<\rho_N=R\), and centers

\[
r_i=(\rho_i+\rho_{i+1})/2,\quad
h_i=\rho_{i+1}-\rho_i,\quad
V_i=(\rho_{i+1}^2-\rho_i^2)/2,\qquad 0\le i<N.
\]

For a nonzero pair \((m,k)\), store
\(a=(q_1,\ldots,q_N,t_0,\ldots,t_{N-1},z_0,\ldots,z_{N-1})\), with \(q_j=u_r(\rho_j)\), \(t_i=u_\theta(r_i)\), \(z_i=u_z(r_i)\). The axis value is \(q_0=0\), appropriate for the retained C4 modes. The additional variable is \(q_N=u_r(R)\). Let \(Q,T,Z\) be the linear maps selecting the full face vector \((q_1,\ldots,q_N)\) and the two center vectors, and let \(e_R a=q_N\).

Set \(\delta=R-r_{N-1}\), and use the exterior logarithmic derivative \(\kappa\) from the companion audit. Define

\[
\ell_\theta=-\frac{im}{R\kappa},\qquad
\ell_z=-\frac{ik}{\kappa},\qquad
t_R=\ell_\theta q_N,\quad z_R=\ell_z q_N.                \tag{1}
\]

These are boundary values defined by the model, not extra independent unknowns. The angular/axial Fourier normalization is common to all inner products and can be factored out.

## 2. Kinetic mass and compatible pressure pair

Interior radial face masses are \(w_j=\rho_j(r_j-r_{j-1})\), \(1\le j<N\), with boundary lump \(w_N=R\delta\). Define

\[
M_{\rm int}=\operatorname{diag}(w_1,\ldots,w_N,V_0,\ldots,V_{N-1},V_0,\ldots,V_{N-1}),
\]
\[
M=M_{\rm int}+\frac R\kappa e_R^*e_R.                  \tag{2}
\]

The exterior term already includes all three exterior velocity components. It must not be supplemented by another exterior tangent/axial kinetic mass.

For pressure weights \(V=\operatorname{diag}(V_i)\), set

\[
(Da)_i=\frac{\rho_{i+1}q_{i+1}-\rho_iq_i}{V_i}
       +\frac{im}{r_i}t_i+ikz_i,
\qquad G=-M^{-1}D^*V.                                 \tag{3}
\]

For the diagonal mass (2), this gives the interior radial pressure differences, the usual \(im p_i/r_i,ikp_i\), and the outer row

\[
(Gp)_N=-\frac{\kappa}{1+\kappa\delta}p_{N-1}.           \tag{4}
\]

Thus the existing tridiagonal \(-DG\) acquires the additional last diagonal entry
\(R\kappa/[V_{N-1}(1+\kappa\delta)]\). The pressure source divergence must also include the outer flux; changing the diagonal without retaining that degree of freedom is not this compatible model. Projection is

\[
\Pi=I-G(DG)^{-1}D.                                    \tag{5}
\]

It is orthogonal in \(M\), idempotent and divergence-free in the discrete sense. Equation (4) implies \(p_R=p_{N-1}/(1+\kappa\delta)\), and (1) applied to its radial gradient gives exactly
\(t_R=im p_R/R\), \(z_R=ikp_R\). This scalar/trace consistency is included in the checker.

## 3. A positive total viscosity Gram

Define the center gather and radial face derivative

\[
(Ia)_i=(q_i+q_{i+1})/2,\qquad
(B_ra)_i=(q_{i+1}-q_i)/h_i.                            \tag{6}
\]

Here \(I\) gathers only radial velocity; \(T,Z\) select center values. Extend the old gather and face derivative to include the new \(q_N\) in their last rows.

The scalar center derivative maps have interior rows

\[
(B_\theta a)_i=(t_{i+1}-t_i)/(r_{i+1}-r_i),\quad
(B_z a)_i=(z_{i+1}-z_i)/(r_{i+1}-r_i),\quad 0\le i<N-1,
\]

and outer rows

\[
(B_\theta a)_{N-1}=(\ell_\theta q_N-t_{N-1})/\delta,
\qquad
(B_z a)_{N-1}=(\ell_z q_N-z_{N-1})/\delta.              \tag{7}
\]

Their derivative quadrature weights are
\(W_g=\operatorname{diag}(w_1,\ldots,w_{N-1},R\delta)\). The old outer no-slip rows \(-t_{N-1}/\delta\), \(-z_{N-1}/\delta\) are **replaced** by (7). They must not remain as additional penalties.

Let \(H=\operatorname{diag}(V_i/r_i^2)\). The cylindrical angular derivative square maps are

\[
A_1=imI-T,\qquad A_2=I+imT,\qquad A_3=imZ.
\]

The interior horizontal gradient Gram is

\[
K_{h,\rm int}=B_r^*VB_r+B_\theta^*W_gB_\theta+B_z^*W_gB_z
               +\sum_{j=1}^3A_j^*HA_j.                \tag{8}
\]

This reproduces the existing angular coupling sign \(+2imI^*H\) in the radial–azimuthal block. It is manifestly Hermitian nonnegative. The new trace terms add only outer local couplings between \(q_N,t_{N-1},z_{N-1}\). The changed dependence on \(k\) is physical to this exterior closure.

The full exterior gradient energy coefficient is

\[
d_{\rm ext}=1+\frac{2R(m^2/R^2+k^2)}\kappa
                 +\frac{m^2}{R^2\kappa^2}.
\]

The complete stiffness is

\[
K=K_{h,\rm int}+k^2M_{\rm int}+d_{\rm ext}e_R^*e_R.     \tag{9}
\]

Equivalently, if the existing implementation adds \(k^2M\) separately, its radial/horizontal matrix must be
\(K_h=K_{h,\rm int}+(d_{\rm ext}-k^2R/\kappa)e_R^*e_R\).
Do not add \(k^2M\) and the full \(d_{\rm ext}\) on top of one another.

The interior stencil (8) follows the current low-order MAC quadrature pattern. It does not supply exact integral viscosity for every interpolated field, including the first small axis interval. The stronger axis treatment below makes the representation regular while preserving the algebra; its consistency still needs manufactured refinement tests.

## 4. Viscous advance and factor reuse

For a viscous substep of length \(\tau\), put \(\alpha=\nu\tau/2\). Solve the constrained Crank–Nicolson system

\[
\begin{pmatrix}M+\alpha K&D^*V\\VD&0\end{pmatrix}
\binom{a^{+}}{\pi}
=\binom{(M-\alpha K)a^{-}}0.                           \tag{10}
\]

Either consistent sign of the two pressure blocks is equivalent after changing the multiplier sign. For \(Da^-=Da^+=0\), taking the real pairing with \(a^++a^-\) gives

\[
\|a^+\|_M^2-\|a^-\|_M^2
 +\alpha\|a^++a^-\|_K^2=0.                            \tag{11}
\]

The nonnegative stiffness and this identity are exact properties of the discrete model. They do not require diffusion to commute with the pressure projection. Both \(M\) and \(K\) now depend on \((m,k)\). Factor caching must include those frequencies and the step coefficient.

Let \(S_z\) negate all axial velocity coefficients and leave radial/azimuthal coefficients fixed. Then
\(M(m,-k)=M(m,k)\), \(K(m,-k)=S_zK(m,k)S_z\), and \(D(m,-k)S_z=D(m,k)\). The current positive/negative axial factor reuse can therefore continue with the corresponding axial sign change. The checker verifies this identity including the exterior trace terms.

## 5. Nonlinear weak projection and observable energy

Use the gather \(\mathcal A a=(Ia,Ta,Za)\) on the common physical quadrature grid. After proper angular/axial dealiasing, let \(f=u_{\rm cell}\times\omega_{\rm cell}\) be the Lamb force. For each retained Fourier pair, scatter its coefficient by

\[
F_{m,k}=M_{m,k}^{-1}\mathcal A^*\operatorname{diag}(V,V,V)f_{m,k},
\qquad \dot a_{\rm nl}=\Pi_{m,k}F_{m,k}.               \tag{12}
\]

For a discretely divergence-free state, the full modal sum of its nonlinear energy work equals the cell-quadrature sum \(\operatorname{Re}\langle u_{\rm cell},u_{\rm cell}\times\omega_{\rm cell}\rangle=0\). This identity requires adjoint transform normalization and sufficiently dealiased products. It does not require the discrete curl to be an exact continuum curl, which remains a separate accuracy issue.

Because the outer inverse mass depends on \(k\), perform the frequency-dependent division after taking the force's axial Fourier coefficients. Retaining the old radial scatter with a fixed physical-space \(1/W\) at the outer face loses the weak identity. The outer radial gather receives half the last cell force, weighted by \(V_{N-1}\), then divided by \(W_R(m,k)\). The exterior itself has \(\omega=0\) in the model, so it contributes no Lamb-force integral; its changing energy is already included by (2) and the projection.

Compute kinetic energy with the Fourier-dependent \(M_{m,k}\), using the existing forward coefficient convention and positive-mode factor of two. A physical-space interior sum alone omits the exterior. A dissipation observable must use (9), including its exterior term. The scalar from Lamb projection is a Bernoulli pressure; a physical-pressure observable must subtract the appropriate \(|u|^2/2\) term, as explained in the companion audit.

## 6. The \((m,k)=(0,0)\) pair

Remove \(q_N\), set both exterior tangential/axial traces to zero, and omit all exterior mass/stiffness terms. The remaining radial variables are \(q_1,\ldots,q_{N-1}\), with \(q_0=q_N=0\). The divergence constraint then forces every radial variable to zero, cell by cell. Center \(t,z\) remain free inside, with their zero exterior trace realized by (7).

The scalar pressure operator has the ordinary constant gauge. Remove one multiplier or solve on the weighted mean-zero pressure space; do not divide by a fabricated \(\kappa\). This handling is tested explicitly, including the zero pressure-gradient gauge and disappearance of radial velocity after projection. The finite-energy harmonic exterior has no nonzero mean radial, azimuthal-circulation, or uniform axial velocity tail.

## 7. What it means to preserve the axis

The outer update leaves all inner flux rows unchanged and keeps \(q_0=0\). This preserves the existing axis boundary treatment. It does not, by itself, enforce the full smooth C4 orders
\(U_+=r^{m+1}F_+(r^2)\), \(U_-=r^{m-1}F_-(r^2)\), \(u_z=r^mF_z(r^2)\) for \(m\ge4\), or the odd horizontal/even axial mode-zero parity.

If those orders are to be structural rather than only checked under refinement, introduce a local full-column-rank injection \(a=C_m b\) from regular axis coefficients and unconstrained outer degrees. For example, on the first three face/center locations at mode four, use two even Taylor coefficients for each helical component and for \(u_z\):

\[
U_+=a_0r^5+a_1r^7,\quad U_-=b_0r^3+b_1r^5,
\quad u_z=c_0r^4+c_1r^6,
\]

with \(u_r=(U_++U_-)/2\), \(u_\theta=(U_+-U_-)/(2i)\), sampled at their own staggered locations. Other nodes remain independent. Scale the powers by the axis-patch radius for conditioning. Smooth interpolation from this local polynomial patch is part of a later continuous reconstruction.

Apply a congruence reduction everywhere:

\[
M_c=C_m^*MC_m,\quad K_c=C_m^*KC_m,\quad
D_c=DC_m,\quad \mathcal A_c=\mathcal A C_m,
\quad G_c=-M_c^{-1}D_c^*V.                             \tag{13}
\]

Equations (5), (10) and (12) then use these reduced operators. Positivity, weighted adjointness, nonlinear weak work and the axis representation are preserved simultaneously. Applying a separate axis filter **after** the old projection is not equivalent. A reduced pressure solve need no longer be the simple unmodified tridiagonal one; a banded/local KKT implementation can retain the sparse coefficient locality. The checker includes this mode-four reduction and verifies its full-rank positive mass, constrained projection and positive viscosity.

This axis prescription is for the retained C4 symmetry. A general Fourier solver containing \(|m|=1\) permits a nonzero radial component at the geometric axis as part of a smooth Cartesian horizontal translation; its trace representation needs the general signed-order helical treatment instead of blindly imposing the present \(q_0=0\).

## 8. Checks and scope

The independent small reference passes six mode/frequency cases, including the zero pair and both signs of a nonzero axial frequency. It checks the Gram against the explicitly summed squares; Hermitian positivity; scalar/velocity trace consistency; projection divergence, idempotence and kinetic orthogonality; nonlinear weak adjointness; and the constrained Crank–Nicolson energy identity. Projection defects are below \(2\cdot10^{-14}\), nonlinear work-pairing defects below \(7\cdot10^{-15}\), and relative CN energy defects below \(2\cdot10^{-15}\). The separate strict-axis congruence check has divergence defect below \(4\cdot10^{-14}\). These are floating-point consistency checks of exact algebraic formulas, not interval enclosures or continuum convergence tests.

The model removes the impermeable radial pressure wall and includes the exact energy of its chosen harmonic exterior. Actual viscous NS produces exterior vorticity, and its heat/vorticity tail is not represented by this constraint. Matching boundary values gives a potential H1 reconstruction but does not by itself match the higher derivatives needed for an H4 residual. A smooth whole-space reconstruction must account for the interior sampling error, interface smoothing and derivative jumps, exterior-vorticity/heat residual, noncompact pressure contribution, angular/radial truncation and axial periodicity. No packet or pump is replanted by this upgrade, and no repeated return or singularity claim follows from its discrete energy law.

# An exact degree-four pressure cutoff with a radial annular pump

Replacing the exterior meridional pump by a radial annular strain makes every pressure mode above degree four vanish from the **full initial mixed pressure-derivative evaluation**. The pressure itself generally has higher modes. This identity removes their contribution to this one scalar functional, provided the retained modes are computed exactly or enclosed rigorously.

Let \(P=S_\eta\) be the same radial strain extension as the core, with \(\eta(|x|^2)\) smooth, compactly supported, and zero near the origin. The exterior datum is
\[
u_o=cP+Av,
\]
where \(v\) is smooth, axisymmetric and purely azimuthal. Overlap between \(P\) and \(v\) is permitted. Keep the radial core \(u_c=bS_\psi+\Omega W_\psi\) and define
\[
\mathcal M(w)=\mathcal C(w)-2\int Q_{ij}(u_o)_i w_j,
\qquad Q_{ij}=-\partial_{zzij}N.
\tag{1}
\]
The core and exterior datum remain disjoint. For the full solenoidal initial acceleration
\[
u_1=\nu\Delta u_0-(u_0\cdot\nabla)u_0-\nabla p,
\]
the audited mixed formula gives \(p_{zz}'(0)=\mathcal M(u_1)\). The following statements concern this precisely defined mixed functional. Splitting its acceleration into separate summands does not require those summands individually to be solenoidal.

## Explicit annular weight

The exterior support avoids the origin, and \(\partial_jQ_{ij}=0\) there. Compact support therefore gives, for any smooth axisymmetric \(\phi\),
\[
2c\int QP\cdot\nabla\phi=-2c\int(Q:\nabla P)\phi.
\tag{2}
\]
There is no contact term or boundary term in (2). Put \(r=|x|\), \(\mu=z/r\), and take all derivatives of \(\eta\) with respect to \(s=r^2\). Direct Cartesian differentiation yields
\[
4\pi r^5 Q:\nabla S_\eta
=\eta(-315\mu^4+270\mu^2-27)
+r^2\eta'(-570\mu^4+468\mu^2-42)
+96r^4\eta''(\mu^2-\mu^4).
\tag{3}
\]
Equivalently,
\[
Q:\nabla S_\eta=\frac1{4\pi r^5}
                  \big[q_0(r)P_0(\mu)+q_2(r)P_2(\mu)+q_4(r)P_4(\mu)\big],
\]
where
\[
\begin{aligned}
q_0&=\frac{64}{5}r^4\eta'',\\
q_2&=-\frac{32}{7}r^2(3\eta'-2r^2\eta''),\\
q_4&=-\frac{24}{35}(105\eta+190r^2\eta'+32r^4\eta'').
\end{aligned}
\tag{4}
\]
For \(\phi=\sum_\ell\phi_\ell(r)P_\ell(\mu)\), equations (2)–(4) give the exact radial formula
\[
2c\int QP\cdot\nabla\phi
=-2c\sum_{\ell=0,2,4}\frac1{2\ell+1}
          \int_0^\infty r^{-3}q_\ell(r)\phi_\ell(r)\,dr.
\tag{5}
\]
The denominator and sign include the full three-dimensional volume integral and the normalization \(N=1/(4\pi r)\).

At every point away from the axis, the axisymmetric tensor \(Q\) preserves the azimuthal and poloidal subspaces. Thus \(Qv\) is azimuthal whereas \(\nabla\phi\) is poloidal, giving
\[
Qv\cdot\nabla\phi=0.
\tag{6}
\]
Continuity extends the identity across the axis. Therefore (5) is the entire exterior pressure-gradient contribution to \(\mathcal M(-\nabla\phi)\).

## Consequence for the full initial gate

The pass 4 core calculation already proves
\[
\mathcal C(-\nabla\phi_\ell)=0\qquad(\ell>4)
\]
for arbitrary smooth radial coefficients, without a harmonicity assumption. Combining this identity with (5) and (6) gives
\[
\boxed{\mathcal M(-\nabla p)
       =\mathcal M(-\nabla P_{\le4}p)
       =\mathcal M\!\left(-\nabla(-\Delta)^{-1}P_{\le4}g\right),}
\tag{7}
\]
where \(g=\operatorname{tr}[(\nabla u_0)^2]\). Odd pressure modes also vanish by parity in this functional, so only \(0,2,4\) are needed. Consequently,
\[
p_{zz}'(0)=
\mathcal M\!\left(\nu\Delta u_0-(u_0\cdot\nabla)u_0\right)
+\mathcal M\!\left(-\nabla(-\Delta)^{-1}P_{\le4}g\right).
\tag{8}
\]
The first term uses the actual datum directly and has no elliptic angular truncation. Equation (8) is an exact scalar reduction, not an approximation asserting that the velocity or pressure contains only three modes.

The truncated acceleration in (8) need not be divergence free; no generic source-to-stress conversion is being applied to it. Linearity of the already-defined mixed functional and the exact pressure selection identity are sufficient.

There is a further special cancellation here. Define the all-source functional
\[
\mathcal A(w)=\partial_{zz}\!\left[
       N*\big(2\operatorname{tr}(\nabla u_0\nabla w)\big)\right](0).
\]
The entire meridional datum is \(S_{b\psi+c\eta}\), one radial strain extension. The all-source pressure contribution is therefore its radial-strain contraction; the contraction of the arbitrary axisymmetric azimuthal datum with \(D^2p\) vanishes pointwise. The previously proved degree-four selection applies to the combined radial profile, including its correct contact at zero, and implies
\[
\mathcal A(-\nabla p)=\mathcal A(-\nabla P_{\le4}p).
\tag{9}
\]
Let \(F=\nu\Delta u_0-(u_0\cdot\nabla)u_0\). Since the exact acceleration is solenoidal,
\[
\boxed{\mathcal A(F-\nabla P_{\le4}p)
       =\mathcal A(F-\nabla p)
       =\mathcal M(F-\nabla p)
       =\mathcal M(F-\nabla P_{\le4}p).}
\tag{10}
\]
Thus all-source and mixed computations also agree **at exact angular cutoff four for this profile class**, although \(F-\nabla P_{\le4}p\) is generally not solenoidal. This equality follows from two separately proved scalar selection rules. It is not a general identity for nonsolenoidal accelerations. A discrepancy between the two numerical evaluations at cutoff four cannot be attributed to a true omitted high-angular pressure contribution; retained-mode projection, radial, differentiation, quadrature, or implementation errors remain possible.

## What this changes for certification

A certificate for this radial-annular datum does not need to bound the high-angular Newton tail for this initial scalar gate. It must still enclose the exact source projections \(g_0,g_2,g_4\), their radial Poisson solutions and endpoint conditions, the retained-mode core contact and radial integrals, the direct local-acceleration integrals, and the parameter/neutral-tuning uncertainty. An angular quadrature error in a retained coefficient cannot be discarded as a higher-mode term.

As an additional exact consistency check, the earlier universal central-pressure formula extends to every compact radial strain profile:
\[
p_{zz}[S_\eta](0)=-\frac{18}{7}\eta(0)^2=0.
\]
Also \(\Delta S_\eta=S_{\widetilde\eta}\), where
\(\widetilde\eta(s)=4s\eta''(s)+14\eta'(s)\).
Since both profiles vanish near zero, polarizing the universal identity gives zero annular self-viscous coefficient. These two cancellations concern the radial strain's self terms; they do not remove its nonlinear interaction with the rotating packet.

No claim of a certified feedback sign, positive-time profile preservation, or repeated amplification follows from (7)–(8). The companion file verify_radial_annular_pressure_cutoff.py checks the Newton kernel sign, the exact angular polynomial and its Legendre coefficients, the azimuthal cancellation, and the radial prefactor using exact symbolic arithmetic.

# Independent audit of full initial temporal coefficients

Files reviewed: `full_initial_evolution.py`, its finite-volume pilot, and
`high_order_pressure.py`. The continuous formulas and full Fourier
interactions are accepted. The finite-volume source repair creates a
specific inconsistency in the computed second coefficient; the later
unrepaired collocation option avoids that source alteration but is still
not a compatible divergence-free discretization. Neither option currently
certifies an actual positive-time whole-space NS value.

## 1. Sound continuous formulas and signs

Let \(u=u_0\), \(g=\operatorname{tr}((\nabla u)^2)\), and
\(F=\nu\Delta u-(u\cdot\nabla)u\). The correct first derivative is

\[
 -\Delta p_0=g,\qquad N_0=u_t(0)=F-\nabla p_0.
\]

For this divergence-free continuous coefficient,

\[
 \begin{split}
 C_1&=(N_0\cdot\nabla)u+(u\cdot\nabla)N_0,\\
 g_1&=2\operatorname{tr}(\nabla u\,\nabla N_0),\\
 -\Delta p_1&=g_1,\\
 N_1=u_{tt}(0)&=\nu\Delta N_0-C_1-\nabla p_1,\\
 \Delta N_0&=\Delta F+\nabla g.                               \tag{1}
 \end{split}
\]

The code has the correct pressure signs, factor two in \(g_1\), matrix
order in the trace, and both placements in \(C_1\). Its signed-mode loops
retain every interaction of initial modes \(0,\pm4\) with first-derivative
modes \(0,\pm4,\pm8\). Thus modes through \(\pm12\) are correctly present
in the second coefficient. The zero and positive Fourier coefficients have
the correct multiplicities. Keeping the original modes alone would have
been wrong; this implementation does not make that omission.

Although \(N_0\) has a nonlocal tail, \(g_1\) and \(C_1\) are supported
where the compact initial field or its gradient is supported. Both initial
pressure solves can therefore use an exterior condition outside those
supports in the exact periodic-axial problem. This does not justify keeping
that same source-support condition for arbitrary positive times.

## 2. Off-neutral core formula and the direct quotient check

For the actual affine initial core, \(b(0)=\Omega(0)=1\),
\(\Delta u_0=\Delta u_t(0)=0\) locally, and

\[
 b_1=-2-\tfrac12p_{0,zz},\quad \Omega_1=2,\quad
 b_2=-4b_1-\tfrac12p_{1,zz},\quad \Omega_2=2b_1+4.
\]

Consequently the exact off-neutral identities are

\[
 \beta_1=-4-\tfrac12p_{0,zz},\qquad
 \beta_2=-\tfrac12(p_{1,zz}+32)-10\beta_1,
 \qquad\beta=b/\Omega.                                      \tag{2}
\]

No pressure-neutral tuning is needed for (2). The selected rational
amplitude remains \(A=1020\).

However, the pressure-based formulas are not automatically the derivatives
of the numerically reconstructed polynomial. Extract from that same field

\[
 b_j=\tfrac12\partial_z(N_{j-1})_{z,0}(0),\qquad
 \Omega_j=\lim_{r\to0}(N_{j-1})_{\theta,0}(r,0)/r
 \quad(j=1,2),                                               \tag{3}
\]

and similarly extract \(b_0,\Omega_0\) from \(U\). Its direct quotient jets
are exactly

\[
 \begin{split}
 \beta'_h(0)&=b_1/\Omega_0-b_0\Omega_1/\Omega_0^2,\\
 \beta''_h(0)&=b_2/\Omega_0-b_0\Omega_2/\Omega_0^2
      -2b_1\Omega_1/\Omega_0^2+2b_0\Omega_1^2/\Omega_0^3.
 \end{split}                                                  \tag{4}
\]

Using the extracted \(b_0\), rather than replacing it by one silently, matters
when the axial Fourier derivative has visible error. Comparing (2) and (4)
is a strong local consistency test. The current `beta_second` field should
be understood as the pressure-based continuous-formula diagnostic until
that comparison passes. It does not determine the derivative of the
displayed polynomial by definition.

## 3. Concrete finite-volume source-repair inconsistency

The original finite-volume solve replaces its zero-angular, zero-axial
source by

\[
 \widetilde g=g-aB(r),\qquad
 B(r)=\chi((r/(R/2))^2).
\]

Even with exact continuous differentiation of this modified Poisson
problem, the reconstructed coefficient would satisfy

\[
 d:=\nabla\cdot N_0=-aB,
 \qquad
 \Delta N_0=\Delta F+\nabla\widetilde g
            =\Delta F+\nabla g-a\nabla B.                   \tag{5}
\]

The current `lapN0` uses the unrepaired \(g\). It therefore differs from the
Laplacian of its own repaired coefficient by \(+a\nabla B\), before adding
any finite-difference error. This changes the computed \(N_1\) by a
gradient contribution \(\nu a\nabla B\) unless a consistent projection
cancels it.

There is another exact missing term when applying (1) to a nonzero-divergence
coefficient:

\[
 \nabla\cdot C_1
   =2\operatorname{tr}(\nabla u\,\nabla N_0)+u\cdot\nabla d.
                                                               \tag{6}
\]

To project this modified first coefficient's differentiated right-hand side
to divergence zero, its scalar equation would instead contain
\(g_1+u\cdot\nabla d-\nu\Delta d\). That would be a different numerical
construction, not a justification for silently treating the modified
\(N_0\) as the true NS derivative. In addition, the second pressure solve
applies its own source repair.

The repair is axially constant, so it does not directly change
\(p_{0,zz}\). In the core where \(B=1\), it changes the horizontal
acceleration by \(-ar e_r/2\), leaves the axial acceleration unchanged, and
makes \(\nabla\cdot N_0=-a\). Thus agreement of one pressure-Hessian entry
cannot establish compatibility. If one otherwise applies the local
continuous differentiated equations to this modified coefficient, the
rotation jet becomes \(\Omega_2=2b_1+4-d\); the direct ratio jet differs
from (2) by \(+d\). This local observation is separate from, and generally
smaller than, the large derivative errors in the coarse pilot.

The pilot reports repair amplitudes approximately \(0.30246\) for \(p_0\)
and \(-845.43\) for \(p_1\), so the alteration is measurable. Recording the
repair is useful, but does not make (5) disappear.

## 4. Unrepaired high-order collocation option

The new `HighOrderCylinder` has the correct signs in the scalar Poisson
operator, pressure gradient and vector Laplacian. Its even scalar and odd
transverse radial extensions are consistent with the parity of these
\(m\in4\mathbb Z\) modes. The band layout contains the one-sided and
reflected stencil entries. Its exterior condition is

\[
 p'(r_b)+\kappa_m p(r_b)=0,\qquad
 \kappa_m=k\frac{K_{m-1}(kr_b)}{K_m(kr_b)}+m/r_b
\]

for \(k>0\), and \(\kappa_m=m/r_b\) for \(k=0,m>0\). These are the
correct decaying Bessel logarithmic derivatives; `kve` cancels the common
exponential factor. The last-node radius is used consistently.

For \(m=k=0\), imposing a pressure gauge at the outer node without repairing
the source avoids the inconsistency (5). A nonzero reported outer flux
still signals incompatibility with a decaying whole-radial zero mode:
outside the source it corresponds to a logarithmic component. Reporting
that flux is honest; it is not a whole-space boundary-error bound.

Further qualifications remain:

- Even/odd reflection alone does not impose every higher vanishing order
  in the smooth helical fields. Those orders and the singular cancellations
  should be tested near the axis as resolution changes.
- The collocation Laplacian is not the discrete divergence of its pressure
  gradient, and convection and differentiation need not obey exact discrete
  product rules. The code correctly does not assert a compatible projection.
- A seven-point stencil is exact on polynomials through degree six.
  Its generic one-sided second derivative is fifth-order accurate, although
  symmetry can give sixth order in the interior. The module's broad
  “sixth-order” description is not a proved global convergence order.
- Periodic axial products are not dealiased. The formal Fourier interactions
  are complete in angle at these two derivative orders, but axial sampling
  can alias products and their subsequent derivatives.

The new option is a useful independent numerical diagnostic. It removes one
specific source modification; it does not establish \(\nabla\cdot N_0=0\)
or make the formulas in (2) automatic identities of its grid polynomial.

## 5. Receiver and time-polynomial scope

The reconstructed polynomial \(U+TN_0+T^2N_1/2\) contains all modes through
twelve and has correctly weighted physical norms. Its finite-time norm is
not an actual NS norm without a temporal remainder or residual enclosure.
In particular, a quadratic approximation to rapid phase rotation can show a
large artificial norm gain even when the exact rotating signal has constant
amplitude. Here \(4AT\) is about four at \(T=10^{-3}\); small absolute time
alone does not justify truncating that phase evolution.

The subsequently added corotating diagnostic has the correct coefficient
identities for a rigid angular coordinate pullback:
\(\widetilde U_m(t)=e^{im\omega t}U_m(t)\), hence
\(\widetilde U'_m=N_{0,m}+im\omega U_m\) and
\(\widetilde U''_m=N_{1,m}+2im\omega N_{0,m}-(m\omega)^2U_m\).
Its least-squares first-jet rotation rate is also correctly minimized by
\(-\langle\partial_\theta u_0,N_0\rangle/\|\partial_\theta u_0\|_2^2\)
in the computed finite-cylinder norm. These identities do not bound the
remaining temporal Taylor error or turn the rotating-frame polynomial into
an actual NS endpoint.

The main first- and second-derivative receiver evaluations Fourier-interpolate
the axial coordinate to exactly \(z=4\). The polynomial receiver instead
uses a rounded axial array index. These coincide for the current even,
commensurate grids, including \(L=16,N_z=1024\). For arbitrary exposed
parameters they need not. Likewise the direct array access at `nz//2`
represents \(z=0\) only for even \(N_z\). A general implementation should
either enforce those restrictions or Fourier-interpolate all such values.

The companion read-only audit scripts separate the core jet terms and
compare (2) with (4). Their floating-point results are diagnostics, not
replacement certified values.

## 6. Quantified fd384 discrepancy

The read-only return-hook run at \(N_r=384,N_z=1024,L=16\) reproduced the
original polynomial and saved `full-initial-jets-audit384.json`. It finds

| Quantity | Reconstructed grid jet | Continuous pressure-based target |
|---|---:|---:|
| \(b_0\) | 0.9999626962 | 1 |
| \(b_1\) | 2.1157264001 | 2.1139663426 |
| \(b_2\) | −545.1629443 | 65.7440400 |
| \(\Omega_1\) | 2.0000000000 | 2 |
| \(\Omega_2\) | 8.1238098923 | 8.2279326852 |
| \(\beta''(0)\) | −553.7496552 | 57.0602419 |

The separate contributions to the reconstructed \(b_2\) are

| Term | Contribution |
|---|---:|
| \(\frac\nu2\partial_z(\Delta F)_z\) | −0.07294685 |
| \(\frac\nu2\partial_{zz}g\) | −0.04496550 |
| \(-\frac12\partial_z(C_1)_z\) | −619.24493722 |
| \(-\frac12p_{1,zz}\) | 74.19990537 |

The viscous terms, which vanish in the exact affine core, are much too small
to explain the discrepancy. The dominant error is in extracting the axial
derivative of the large outer convection product:

\[
 \partial_z(C_1)_z\big|_{\text{axial Fourier}}=1238.4898744,
 \qquad
 \partial_z(C_1)_z\big|_{\text{local centered stencil}}=16.8913381.
\]

The affine target is approximately \(16.91\). This isolates severe spectral
differentiation leakage/underresolution of the convection product; the
comparison alone does not distinguish every source of sampling and product
aliasing. The local finite difference is a diagnostic, not a certified
replacement. At this grid the axis divergence defect of \(N_0\) is only
about \(0.00352085\), so that local defect alone also does not explain a
six-hundred-unit error in \(b_2\).

The quotient formed from the directly extracted jets equals
\(0.9998021756829721\) at \(T=10^{-3}\), agreeing with the original
polynomial's reported \(0.9998021756829528\) to \(2\times10^{-14}\).
Thus the disagreement is real within the numerical reconstruction; it is
not a quotient-algebra mistake. It should be resolved by spatial/axial
resolution and consistent product differentiation before any direct
polynomial core gain is interpreted.

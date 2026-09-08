# Error-bound viability and a more selective commutator estimate

The uniform-in-time pass9 bound `560 V5` is numerically impractical for the
selected datum at `T=0.001`, even if a much more accurate evolution is
computed. This is a limitation of that bound, not evidence of NS
instability. A five-component derivative hierarchy gives a concrete next
validation target while retaining the higher derivatives needed for H4.

## 1. A rigorous initial norm lower bound without quadrature

Set \(a=63/200\), \(r_0=a/2=63/400\), \(d=9/10\), and \(A=1020\).
On \(0\le r\le r_0\), the axisymmetric outer azimuthal component is exactly
\(Arq_v(z)\). The two axial supports each have length \(d\), and their
plateaus have combined length \(d\). The compact core swirl is disjoint
from these axial supports. The poloidal field is orthogonal to the
azimuthal component, and the nonaxisymmetric seed is orthogonal under
angular averaging. Thus, for \(j=4,5\),

\[
 \|u_0\|_{H^j_{\rm der}}^2
 \ge A^2\frac{\pi r_0^4}{2}\|q_v^{(j)}\|_{L^2_z}^2
 \ge A^2\frac{\pi r_0^4}{2}(\pi/d)^{2j}d
 =\frac{A^2r_0^4\pi^{2j+1}}{2d^{2j-1}}.                    \tag{1}
\]

The second inequality applies the interval Dirichlet Poincare inequality
\(\|f'\|_2\ge(\pi/d)\|f\|_2\) repeatedly to each smooth, flat-ended
bump and its derivatives; its sine-series proof applies at every step.
Only the pure axial derivative is used, so no rotational invariance of the
derivative-sum H4/H5 norm is being assumed.

Even the elementary bound \(\pi>3\) gives the exact rational comparisons

\[
 \|u_0\|_{H^4_{\rm der}}^2>
 \frac{843075135}{64}>3600^2,\qquad
 \|u_0\|_{H^5_{\rm der}}^2>
 \frac{2341875375}{16}>12000^2.                              \tag{2}
\]

If the approximation is initialized exactly, any **uniform** bound
\(V_5\ge\sup_{[0,T]}\|v(t)\|_{H^5_{\rm der}}\) is therefore greater than
12000. The baseline linear exponent at \(T=10^{-3}\) exceeds 6720. Even
if initialization is known only in H4 with \(\|v(0)-u_0\|_{H^4}\le e_0\),
the same uniform coefficient obeys

\[
 560V_5T>0.56(3600-e_0).                                    \tag{3}
\]

For \(e_0\le1\), this exponent still exceeds 2015.

The linear part alone propagates \(e_0\) to \(e_0e^{560V_5T}\), and a
uniform residual bound \(\delta\) contributes
\(\delta(e^{560V_5T}-1)/(560V_5)\). The positive Riccati term only worsens
these bounds. At exact initialization, an exponent 6720 demands thousands
of decimal digits of residual control to certify an ordinary-sized endpoint
margin. Exact zero initial error and zero residual would avoid the issue,
but such an exact evolution representation has not been supplied.

This argument concerns the uniform `V5` route. A lower bound at time zero
alone does **not** lower-bound \(\int_0^T\|v(t)\|_{H^5}dt\). A proposed
time-dependent version must enclose that integral rather than infer it
from (2).

## 2. Retain the derivative hierarchy and cancel rigid rotation

Let \(u,v\) be smooth divergence-free whole-space velocities, with \(u\)
solving NS and the full projected residual of \(v\) denoted by

\[
 \mathcal R=v_t-\nu\Delta v+\mathbb P(v\cdot\nabla v).
\]

Put \(e=u-v\). Use the **ordered Cartesian derivative tensors**, with
their full Frobenius norms, and choose a fixed verification length
\(\ell>0\):

\[
 E_j=\ell^j\|\nabla^j e\|_{L^2},\quad j=0,\ldots,4,
 \qquad E_\ell=(\sum_{j=0}^4E_j^2)^{1/2},
 \qquad \delta_j=\ell^j\|\nabla^j\mathcal R\|_{L^2}.          \tag{4}
\]

Ordered derivatives give the mixed-derivative multiplicities needed for
rotational invariance. In particular the antisymmetric part of a velocity
gradient acts skew-symmetrically on each tensor index.
Explicitly,
\(E_\ell^2=\sum_{|\alpha|\le4}\ell^{2|\alpha|}
(|\alpha|!/\alpha!)\|\partial^\alpha e\|_2^2\); these multiplicities
range from one to twelve. At \(\ell=1\), the ordered norm lies between
the derivative-sum H4 norm and \(\sqrt{12}\) times that norm.

The quantities to enclose for the **entire approximating field** are

\[
 S=\|\operatorname{sym}\nabla v\|_{L^\infty,op},\qquad
 M_j=\|\nabla^j v\|_{L^\infty,Frob},\quad j=2,3,4,5.          \tag{5}
\]

They include the cutoff joins, nonaxisymmetric field and exterior tail.
Bounds only for the rotating plateau do not suffice.

With rows and columns numbered zero through four, define

\[
 B_\ell=
 \begin{pmatrix}
 S&0&0&0&0\\
 \ell M_2&2S&0&0&0\\
 \ell^2M_3&3\ell M_2&3S&0&0\\
 \ell^3M_4&4\ell^2M_3&6\ell M_2&4S&0\\
 \ell^4M_5&5\ell^3M_4&10\ell^2M_3&10\ell M_2&5S
 \end{pmatrix},\qquad q=(0,1,3,7,15)^T.                    \tag{6}
\]

The following componentwise estimate is sufficient:

\[
 D^+E_j\le(B_\ell E)_j+
        q_j\ell^{-5/2}E_\ell^2+\delta_j.                    \tag{7}
\]

The viscous term was dropped with its favorable sign. Equation (7) is a
derived error estimate, not an evaluated validation claim.

### Proof of the stated constants

Differentiate
\(e_t-\nu\Delta e+
 \mathbb P(v\cdot\nabla e+e\cdot\nabla v+e\cdot\nabla e)
 =-\mathcal R\)
with every ordered derivative of order \(j\), pair with that differentiated
error, and sum. Pressure drops out because each derivative of \(e\) is
divergence-free. The top transport terms cancel by integration by parts.

The first-derivative commutator from \(v\cdot\nabla e\) acts on the
\(j\) derivative indices. Its symmetric part costs at most \(jS E_j\).
The term with all derivatives on \(e\) in \(e\cdot\nabla v\) acts on
the velocity-component index and costs at most \(S E_j\). Antisymmetric
parts cancel in each case. This proves the diagonal \((j+1)S\), including
exact cancellation of a rigid-rotation contribution.

The remaining transport commutators have \(k\ge2\) derivatives on \(v\):

\[
 \binom jk\ell^{k-1}M_k E_{j-k+1}.
\]

The remaining \(e\cdot\nabla v\) terms, with \(k\ge1\) derivatives
on its gradient, contribute

\[
 \binom jk\ell^k M_{k+1}E_{j-k}.
\]

Tensor contraction is bounded by the product of Frobenius norms, without
an extra directional counting factor. Combining terms by Pascal's identity
gives \(\binom{j+1}{h}\ell^{j-h}M_{j-h+1}E_h\) for \(0\le h<j\),
which is precisely (6).

For self-advection, each commutator product has derivative orders adding
to \(j+1\le5\), with its smaller order at most two and its larger order at
most four. The scalar H2 embedding constant used in pass3 is less than one.
Summing componentwise over ordered tensors therefore gives
\(\|\nabla^a e\|_\infty\le\ell^{-a-3/2}E_\ell\) for \(a\le2\),
by rescaling \(x=\ell y\). Putting the lower derivative factor in
\(L^\infty\) bounds every weighted commutator product by
\(\ell^{-5/2}E_\ell^2\). The sum of its binomial weights is
\(2^j-1=q_j\). This proves (7), including its nonlinear constant.

## 3. The actual next computation

Prefer a validated integration of the five nonnegative majorants

\[
 z_j'=(B_\ell z)_j+q_j\ell^{-5/2}\sum_i z_i^2+\delta_j,
 \qquad z_j(0)\ge E_j(0).                                  \tag{8}
\]

All right-hand sides are monotone in the nonnegative coordinates when
using nonnegative upper bounds for (5), so the usual first-contact
comparison proves \(E_j\le z_j\). Higher derivatives of \(v\) enter only
**lower-triangular couplings**. In the linear part they feed higher error
derivatives from lower ones; they need not all be put into a common
exponential growth coefficient. Do not replace this system by merely
\(E'\lesssim S E\): that would omit the displayed \(M_2,\ldots,M_5\)
terms.

If a scalar radius is operationally simpler, the valid alternative is

\[
 D^+E_\ell\le a_\ell E_\ell+17\ell^{-5/2}E_\ell^2+
                  (\sum_j\delta_j^2)^{1/2},\qquad
 a_\ell=\lambda_{\max}\bigl((B_\ell+B_\ell^T)/2\bigr).       \tag{9}
\]

Here \(\sqrt{1+9+49+225}=\sqrt{284}<17\). The five-by-five matrix
eigenvalue needs an enclosure. This scalar reduction can lose the advantage
of the triangular hierarchy, so (8) is the preferred first test.

Choose \(\ell\) before the error propagation and optimize it against the
endpoint tests, rather than selecting it only to minimize one coefficient.
Smaller \(\ell\) reduces the high-derivative entries in (6), but increases
the nonlinear and pointwise conversion factors. No favorable value of
\(\ell\), interval, or residual budget has yet been certified for this
candidate. Anisotropic lengths may be considered later, but require their
own transformed-index strain matrices; the cancellation above should not
be assumed for an arbitrary anisotropic norm.

## 4. Endpoint conversion and remaining obligations

If (8) gives \(\rho_\ell=(\sum z_j(T)^2)^{1/2}\), then

\[
 \begin{split}
 \|e(T)\|_{L^2}&\le\rho_\ell,\\
 \|e(T)\|_\infty&\le\ell^{-3/2}\rho_\ell,\\
 \|\nabla e(T)\|_{\infty,Frob}&\le\ell^{-5/2}\rho_\ell,\\
 \|e(T)\|_{H^4_{\rm der}}&\le\max(1,\ell^{-4})\rho_\ell.
 \end{split}                                                  \tag{10}
\]

These constants are conservative; the scalar embedding constant below one
could be retained. Consequently:

- The mean receiver loses at most \(r_T\ell^{-3/2}\rho_\ell\).
- The fluctuation energy obeys
  \(K(u(T))\ge\frac12(\sqrt{2K(v(T))}-\rho_\ell)_+^2\).
- Each central \(b\) and \(\Omega\) error is bounded by
  \(\ell^{-5/2}\rho_\ell\), sufficient for the pass9 linear core-cone
  test with this substituted error.
- A successor class stated in the original H4 norm must use the last line
  of (10), followed by its already specified spatial/amplitude rescaling.

The inherited whole field, periodization correction, exterior reconstruction,
initial error and all numerical residuals remain necessary. Pressure has
been retained through the full Leray projection. The new estimate reduces
a demonstrable overestimate in the original error coefficient; it does not
validate the current time polynomial or eliminate the endpoint-return
problem.

`verify_sharp_error_hierarchy.py` checks the rational norm lower bounds,
all Pascal coefficients, tensor multiplicities and embedding bookkeeping,
nonlinear scaling and constants, and exact cancellation of rigid rotation
on tensor orders zero through four. These are algebraic checks; neither
(5) nor the error system (8) has been numerically enclosed along a proposed
evolution.

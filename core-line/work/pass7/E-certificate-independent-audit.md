# Independent audit of the E certificate and the common initial family

Accepted: the whole-space Green/residual construction in
certify_outer_E.py retains all pressure, and its range bounds certify
\(E<1/60\) from the stored 2048-panel outward-rounded computation.
This is an upper bound for one initial cubic coefficient. It is not a
pressure approximation used in the actual Navier–Stokes datum.

The derivations and symbolic checks are in
E-interaction-reduction.md and verify_E_interaction_reduction.py.
The independent checks below cover the additional separated inequalities,
range-integration logic, and exact dyadic comparison.

## 1. Kernel and measure constants

On the seed support \(0\le t=r^2/z^2\le(7/68)^2\),
\[
 |Q_{rr}-Q_{\theta\theta}|
 \le\frac{45r^2}{2\pi|z|^7},\qquad
 |Q_{rz}|\le\frac{15r}{\pi|z|^6}.
\]
They follow respectively from
\[
 Q_{rr}-Q_{\theta\theta}
 =\frac{15r^2(r^2-6z^2)}{4\pi(r^2+z^2)^{9/2}},\qquad
 Q_{rz}=-\frac{15rz(4z^2-3r^2)}{4\pi(r^2+z^2)^{9/2}}.
\]
Set \(H=(e'-e/r)^2+k^2e^2\). Equation (3) in the reduction note gives
\[
 |d_m|\le m\left[
 \frac{45r}{2\pi|z|^7}\sqrt H\,q_s+
 \frac{15}{\pi|z|^6}e\,|q_s'|\right].
\]
Minkowski's inequality in the physical measure \(\pi r^3\,dr\,dz\)
therefore gives precisely the script's nd expression, with radial
integrands \(r^5H\), \(r^3e^2\) and axial integrands
\(q_s^2/z^{14}\), \((q_s')^2/z^{12}\).

Similarly \(2\pi\int r|d_m||PF|\) is bounded by
\[
 m\left[
 45\left(\int r^2|P|\sqrt H\,dr\right)
     \left(\int q_s^2q_v/|z|^7\,dz\right)
 +30\left(\int re|P|\,dr\right)
     \left(\int|q_s q_s'|q_v/|z|^6\,dz\right)\right],
\]
which is exactly app. All axial quantities are even; doubling the
positive packet integral is correct for the norms, absolute pairing,
and local coefficient.

The local coefficient is exactly
\[
 E_{\rm loc}=-\frac{45}{2}mk
 \int_{\mathbb R}\int_0^\infty
 r^2 e^2 C_0\,\frac{q_s^2q_v}{|z|^7}
             \frac{6-t}{(1+t)^{9/2}}\,dr\,dz,\qquad t=r^2/z^2.
\]
Since \((6-t)/(1+t)^{9/2}\ge(6-t)/(1+t)^5\), and the latter decreases
throughout the cone, its lower bound at \(t_{\max}=(7/68)^2\) gives
the script's local_gain. Its lower interval endpoint is used when
subtracting it from an upper bound.

## 2. Whole-space range integration

The radial cells cover \([.07,.35]\) and the derivative-reduced Green
moments in the reduction note. Each exact prefix/suffix integral lies
in its outward range sum. Within a cell, the as-yet partial integral
lies in \([0,\mathrm{cell\ width}]\) times that cell's integrand range.
The same enclosing length interval may be used independently in the
left and right expressions: this discards correlation and enlarges
the range, without excluding the true pressure.

The radial norm includes the exact \(r^m\) tail below .07 and the exact
\(r^{-m}\) tail above .35. Although the outer swirl already ends at
.315, integrating the middle interval to .35 and starting the tail
there is equally exact. The tail constants include the entire
oscillatory moment; they have not been replaced by zero.

The axial chain rules for \(q_s''\), \(q_v''\), and
\(F''=q_s''q_v+2q_s'q_v'+q_sq_v''\) are correct.
Thus the residual norm is the complete
\[
 \|re\|_2=
 \left[\pi\int_0^\infty r^3|P|^2\,dr
             \int_{\mathbb R}|F''|^2\,dz\right]^{1/2}.
\]
The pressure-error multiplier \(2/m^2\) is correct. No finite-box
boundary or discarded pressure tail occurs.

The interval engine encloses the unchanged smooth cutoff and its
derivatives; rational panel endpoints and outward arithmetic are used
in the integrals. The result continues to rely on that elementary
range engine and the analytic identities, as do the earlier
certificates. Its printed intervals enclose upper-bound expressions,
not two-sided intervals for \(E\) or for the pressure pairing.

The saved upper dyadic endpoint implies
\[
 E<0.016616140729730<1/60.
\]
The decimal is explanatory; the independent checker compares the exact
dyadic endpoint to \(1/60\) using rational arithmetic.

## 3. The same finite family passes both initial tests

Use the separately derived bound
\[
 D_\nu<4-\frac{96}{5}C_w,
\]
and the earlier exact retuning and family bounds
\[
 T_\lambda<-m_0+\lambda^2
  \left[\frac{23199}{1000}C_w+D_\nu+A_\lambda E\right],\qquad
 m_0=\frac{290649}{8750},\quad
 0<C_w<\frac32,\quad 0<A_\lambda<1020.
\]
Since \(E<1/60\), positivity of \(A_\lambda\) gives
\(A_\lambda E<17\), regardless of the unknown actual sign of \(E\).
Keeping the favorable \(C_w\) term in \(D_\nu\) yields
\[
 \frac{23199}{1000}C_w+D_\nu+A_\lambda E
 <\frac{3999}{1000}\frac32+4+17
 =\frac{53997}{2000}.
\]
The last constant is positive. For every \(3/4\le\lambda\le1\),
\[
 \boxed{T_\lambda<
 -\frac{290649}{8750}+\frac{53997}{2000}
 =-\frac{435297}{70000}<-6.2185.}
\]
This joins the already proved initial mean-maximum increase and the
strict full initial pressure-feedback sign on that actual compact
three-dimensional family. It does not show that either sign persists
for a prescribed nonzero interval, that the profile regenerates, or
that repeated shrinking returns exist.

verify_E_certificate_audit.py replays the exact stored dyadic
comparison, the separated constants, and this rational common-family
budget. The E interval integrator was also rerun independently without
overwriting the stored certificate.

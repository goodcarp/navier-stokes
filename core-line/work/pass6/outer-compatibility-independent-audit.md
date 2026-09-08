# Independent audit of the outer-seed compatibility reduction

**Accepted.** The mean-maximum result and the full pressure-derivative test are correctly separated. The former proves a positive initial derivative of the actual azimuthal-mean maximum for \(3/4\le\lambda\le1\), with exact central pressure neutrality. It does not claim that the full gate \(T_\lambda<0\) has been proved at those amplitudes. No evaluation of the two new interaction coefficients was performed in this audit.

## Mean maximum and retuning

The seed envelope and its phase are compatible with the stated odd, fourfold symmetry and zero cylindrical angular mean. Its support avoids the axis and core, lies on the annular pump plateau, and its flat envelope covers the indicated exterior maximizing circles.

For this horizontal single-mode seed,
\[
\overline{w_rw_\theta}=-40E^2/r.
\]
The azimuthal-mean circulation equation therefore has the exact seed source
\[
-\frac1r\partial_r(r^2\overline{w_rw_\theta})
=\frac{40}{r}\partial_r(rE^2),
\]
so its value is \(40/r>0\) on the seed plateau. At the chosen mean-maximizing circles the mean advective terms vanish, the axial base cutoff is flat, and the stated radial diffusion term is exact. The nonzero Fourier mode has no direct mean viscous contribution; the angular pressure derivative averages to zero.

The positive mean maximum is indeed exterior: the core positive swirl satisfies \(r^2[\psi+(2/3)s\psi']\le1\) because \(\psi'\le0\), whereas the displayed outer lower bound exceeds nineteen. Evaluating at an initially maximizing circle with a strictly positive time derivative proves the asserted short-time increase of the supremum without differentiating a moving maximizer.

Here \(C_w=-\Pi(w,w)>0\), so the correct neutral retuning is
\[
A_\lambda^2C_v=H_0-\lambda^2C_w.
\]
This agrees with the earlier inner-seed convention \(J_w=\Pi(w,w)=-C_w\); there is no sign reversal between the two constructions. The positivity and amplitude bounds used in the mean calculation ensure \(A_\lambda>0\) throughout its interval. The mean-maximum inequality alone supplies neither the sign of \(p_{zz}'\) nor the second central strain-ratio derivative.

## All three cubic placements are present

At a fixed axisymmetric \(U\), expand
\[
F_\nu(U+\lambda w)
=2\Pi(U+\lambda w,\nu\Delta(U+\lambda w)
                  -B(U+\lambda w,U+\lambda w)).
\]
The quadratic coefficient in \(\lambda\) is exactly
\[
2\nu\Pi(w,\Delta w)
-2\{\Pi(U,B(w,w))+\Pi(w,B(U,w))+\Pi(w,B(w,U))\}.
\]
The unsymmetrized \(B\) requires both last placements; the definition of \(\mathcal L(U;w)\) includes them with the correct factor and sign.

A rotation by \(\pi/4\) fixes \(U\) and negates \(w\). The scalar functional is invariant under that rotation, so its odd powers of the seed amplitude vanish. Its viscous part has degree two and its inviscid part degree three. Taking \(U=M+A v\), linearity of \(\mathcal L\) in \(U\) leaves exactly
\[
\boxed{T_\lambda=T_{\rm base}(A_\lambda)
                         +\lambda^2(D_\nu+A_\lambda E).}
\]
There is no missing cubic seed term or term proportional to \(A_\lambda^2\lambda\). Conversely, angular mean zero does **not** eliminate \(A_\lambda E\): it contains two seed factors. Since \(A_\lambda\) is a parameter of the initial datum, evaluating the initial time derivative does not introduce a derivative of \(A_\lambda\) with respect to \(\lambda\).

## Direction and value of the sufficient inequality

The previous radial bound applies to every positive \(A\), even when that base field alone is not neutrally tuned:
\[
T_{\rm base}(A)<102-\frac{23199}{1000}A^2C_v.
\]
Direct substitution, rather than subtraction of upper bounds, gives
\[
T_\lambda<-m_0+\lambda^2
 \left[\frac{23199}{1000}C_w+D_\nu+A_\lambda E\right].
\]
At \(\lambda=3/4\), the upper bound \(C_w<3/2\) proves that
\[
\boxed{D_\nu+A_{3/4}E<
\frac{16m_0}{9}-\frac{23199}{1000}\frac32
=\frac{5093339}{210000}=24.2539952381\ldots}
\]
is sufficient for \(T_{3/4}<0\). This same amplitude has the separate mean-maximum lower derivative \(626/35>0\).

For independently enclosed coefficients, the proposed upper estimate
\[
D_+ +\max\{780E_+,1020E_+\}
\]
is valid for either sign of \(E_+\), because \(A_\lambda>0\) and lies between the two endpoints. It avoids an incorrect endpoint choice when \(E_+\) is negative.

The two notes expressly leave \(D_\nu\) and \(E\) unevaluated. Their sufficient inequality is a concrete remaining calculation, not a certificate that its premises already hold. Sustained joint feedback, sufficient accumulated transfer and return on inherited data remain additional tasks. The exact compatibility checker was independently read and rerun successfully.

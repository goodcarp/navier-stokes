# The remaining outer-seed pressure test reduces to two interaction coefficients

The outer family in [the mean-maximum calculation](outer-mean-maximum-target.md) has a verified initial mean-maximum gain for \(3/4\le\lambda\le1\), with exact central pressure neutrality. The remaining initial test is its full \(T_\lambda=p_{zz}'(0)+32\). This note reduces that test algebraically; it does not evaluate the two new coefficients.

Let \(M=S_\Phi+W_\psi\), retain the unit axisymmetric outer swirl \(v\), and let \(w\) be the outer fourfold seed. Write

\[
 u_\lambda=M+A_\lambda v+\lambda w,\qquad
 A_\lambda^2 C_v=H_0-\lambda^2C_w,\quad H_0=204/35.
\]

Use the unsymmetrized \(B(a,b)=\mathbb P(a\cdot\nabla b)\) and the pressure pairing \(\Pi\) from the preceding notes. Define the linear-in-\(U\) interaction

\[
 \mathcal L(U;w)=-2\{\Pi(U,B(w,w))
                  +\Pi(w,B(U,w))+\Pi(w,B(w,U))\}.
 \tag{1}
\]

Set

\[
 D_\nu=\mathcal L(M;w)+2\nu\Pi(w,\Delta w),\qquad
 E=\mathcal L(v;w).
 \tag{2}
\]

These definitions retain every local and nonlocal pressure acceleration. In particular, the term \(A_\lambda E\) may not be discarded merely because the seed has zero angular mean.

Rotation through \(\pi/4\) leaves \(M,v\) fixed and sends \(w\) to \(-w\). The axial pressure functional at the origin is rotation invariant. Its viscous part is quadratic and its Euler part cubic in velocity. Therefore all odd powers of \(\lambda\) vanish exactly. Expanding the full functional gives

\[
 \boxed{T_\lambda=T_{\rm base}(A_\lambda)
                   +\lambda^2(D_\nu+A_\lambda E).}
 \tag{3}
\]

There are no additional \(\lambda^3\), \(A_\lambda^2\lambda\), or linear seed terms. The coefficient \(E\) is the potentially important mean interaction of the strong axisymmetric swirl with two seed factors.

The old radial base calculation applies at every \(A>0\), with the same profiles:

\[
 T_{\rm base}(A)
 =32+C_{300}(\Phi)+t_{Wb}(\psi)
  +A^2 C_v\left\langle k_{\rm core}+\tfrac75k_{\rm pump}\right\rangle
  +\nu A^2d_v.
\]

The existing component bounds at \(\nu=1/1000\) imply

\[
 T_{\rm base}(A)<102-\frac{23199}{1000}A^2C_v.
\]

Substituting the exact retuning, and setting \(m_0=290649/8750\), yields the useful sufficient bound

\[
 \boxed{T_\lambda<
  -m_0+\lambda^2\left[\frac{23199}{1000}C_w+D_\nu+A_\lambda E\right].}
 \tag{4}
\]

This step uses the old coefficient inequality in its valid direction while multiplying by the positive retuned \(A_\lambda^2C_v\). It does not subtract two separately upper-bounded base values.

At the already torque-admissible amplitude \(\lambda=3/4\), \(C_w<3/2\), so one concrete sufficient target is

\[
 \boxed{D_\nu+A_{3/4}E<
  \frac{16m_0}{9}-\frac{23199}{1000}\frac32.}
 \tag{5}
\]

The right side is the explicit rational number \(5093339/210000\), approximately \(24.254\). The same amplitude already has an initial mean-maximum derivative exceeding \(626/35\). Proving (5) would therefore join those two initial requirements on one actual datum.

For example, certified upper bounds \(D_\nu\le D_+\), \(E\le E_+\), together with \(780<A_{3/4}<1020\), would suffice if the right side of

\[
 D_+ + \max\{780E_+,1020E_+\}
\]

is less than the threshold in (5). Sharper direct enclosure of \(D_\nu+A_{3/4}E\) may retain useful cancellation. Ordinary quadrature agreement would not certify either bound.

The broad \(H^{10}\) continuity estimate for an infinitesimal seed is not a substitute for this finite-amplitude test. Neither \(D_\nu\) nor \(E\) has been evaluated or enclosed in this pass. After (5), sustained joint feedback and mean transfer, the inherited endpoint, and repeated compatibility would still remain to be proved.

The exact checker verifies polynomial degree/rotation selection and the rational sufficient threshold. The pressure interactions in (1)–(2) are the actual next quantities to compute.

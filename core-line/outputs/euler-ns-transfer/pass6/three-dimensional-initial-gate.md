# A quantitative neighborhood of genuinely three-dimensional favorable data

The pass5 initial pressure certificate can be embedded in a compact, smooth, nonaxisymmetric family that also has a strict positive Reynolds-stress transfer into a rotational receiver. This removes the axisymmetric maximum-principle obstruction at the level of admissible initial data. It does not prove a repeatable three-dimensional stage.

## Family and exact neutral tuning

Let \(u_*=S_\Phi+W_\psi+A_0v\) be the certified neutral datum of pass5, with
\(\nu=1/1000,\ H_0=204/35,\ A_0^2C_v=H_0\).
Let \(w\) be the explicitly specified fourfold seed in
[the receiver calculation](rotational-receiver-torque.md), (6).
It is supported inside the original affine ball but vanishes near the origin.

Let \(p[w]\) be its actual decaying whole-space pressure and put
\(J_w=p[w]_{zz}(0)\). Define

\[
 A_\varepsilon=
   \sqrt{\frac{H_0+\varepsilon^2J_w}{C_v}},
 \qquad
 u_\varepsilon=S_\Phi+W_\psi+A_\varepsilon v+\varepsilon w.
 \tag{1}
\]

For sufficiently small \(\varepsilon\), the square root is positive. This is an exact definition, not a floating approximation to a tuned amplitude.

The functional \(p_{zz}(0)\) is invariant under rotations about the axis. The base is axisymmetric, while \(w\) has only angular modes \(m=\pm4\) in its cylindrical components. All bilinear base/seed contributions to this functional therefore vanish. Since \(p[v]_{zz}(0)=-C_v\),

\[
 \boxed{p[u_\varepsilon]_{zz}(0)=-8.}
 \tag{2}
\]

The actual initial field remains \(Dx+Jx\) near the origin. Oddness and fourfold rotational symmetry are preserved by uniqueness. They force the central velocity to be zero and its gradient to have the form \(bD+\Omega J\), even though the full flow is not axisymmetric.

The same central equations and affine-neighborhood curvature cancellations used in pass5 therefore apply at insertion:

\[
 b=\Omega=1,\quad b'=\Omega'=2,\quad
 \beta'=0,\qquad
 \beta''=-\tfrac12(p_{zz}'+32).
 \tag{3}
\]

Here the full three-dimensional pressure and viscous acceleration are used. Axisymmetry of the entire field is not required for these central identities once the stated discrete symmetry and local affine condition hold.

## Explicit norm-form bound for the pressure derivative

Use the pass3 derivative Sobolev norm
\(\|a\|_s^2=\sum_{|\gamma|\le s}\|\partial^\gamma a\|_2^2\).
For smooth solenoidal fields define

\[
 \Pi(a,b)=
 \partial_{zz}\left[N*\operatorname{tr}(\nabla a\nabla b)\right](0).
\]

The multiplier of \(\partial_{zz}N*\) has absolute value at most one. Scalar \(H^2\) embedding and Leibniz give

\[
 \boxed{|\Pi(a,b)|\le128\|a\|_4\|b\|_4.}
 \tag{4}
\]

Indeed there are nine gradient-product pairs, the sum of differentiated-product coefficients is at most four, and there are ten derivatives of order at most two. The lower derivative factor has order at most two and can be placed in \(L^\infty\) using its \(H^4\) norm. Thus the sufficient count is \(9\cdot4\sqrt{10}<128\).

Write \(B(a,b)=\mathbb P(a\cdot\nabla b)\). The previously proved bound is
\(\|B(a,b)\|_4\le512\|a\|_6\|b\|_6\).
The complete initial pressure derivative is the polynomial functional

\[
 F_\nu(u)=2\Pi(u,\nu\Delta u-B(u,u)).
 \tag{5}
\]

For \(X=\|u\|_{10}\), \(Y=\|v_0\|_{10}\) and
\(d=\|u-v_0\|_{10}\), bilinearity and (4) give

\[
 \boxed{|F_\nu(u)-F_\nu(v_0)|
 \le256d\left[3\nu(X+Y)+512(X^2+XY+Y^2)\right].}
 \tag{6}
\]

The name \(v_0\) in this estimate denotes an arbitrary comparison field, not the outer swirl profile in (1).

Let

\[
\begin{gathered}
 X_*=\|u_*\|_{10},\quad Z=X_*+1,\quad
 C_F=256(6\nu Z+1536Z^2),\\
 m_0=\frac{290649}{8750},\qquad
 d_*=\min\{1,m_0/(2C_F)\}.
\end{gathered}
\tag{7}
\]

If \(\|u_\varepsilon-u_*\|_{10}\le d_*\), then \(X,Y\le Z\), so the pass5 bound \(F_\nu(u_*)+32\le-m_0\) implies

\[
 \boxed{F_\nu(u_\varepsilon)+32\le-m_0/2,\qquad
 \beta_\varepsilon''(0)\ge\frac{290649}{35000}>0.}
 \tag{8}
\]

These are full initial pressure estimates, with no angular truncation of the new field.

For a completely specified sufficient range, set
\(W=\|w\|_{10}>0,\ V=\|v\|_{10}>0\) and
\(D=128A_0VW^2/H_0>0\). Define

\[
 \boxed{\varepsilon_0=
  \min\left\{
     1,\sqrt{\frac{H_0}{256W^2}},
     \frac{d_*}{2W},\sqrt{\frac{d_*}{2D}}
      \right\}>0.}
 \tag{9}
\]

By (4), \(|J_w|\le128W^2\). Therefore the second restriction gives
\(|\varepsilon^2J_w|\le H_0/2\), making (1) legitimate. Moreover,

\[
 |A_\varepsilon-A_0|
 =\frac{\varepsilon^2|J_w|}
      {C_v(A_\varepsilon+A_0)}
 \le \frac{128A_0W^2}{H_0}\varepsilon^2.
\]

Hence \(\|u_\varepsilon-u_*\|_{10}\le W|\varepsilon|+D\varepsilon^2\le d_*\)
for \(|\varepsilon|\le\varepsilon_0\).
Every nonzero parameter in this interval gives genuinely three-dimensional data satisfying (8), as well as the exact positive receiver transfer (7) in the companion note.

The quantities in (9) are finite integrals of explicitly specified smooth compact profiles. No decimal lower bound for \(\varepsilon_0\) has been computed. Its norm-based range can be extremely small. It does not show that the certified perturbation is large enough to pay the repeated-stage torque budget.

## What this does and does not provide

Smooth local existence and the strict margin in (8) give a positive interval with favorable central strain-ratio feedback for each chosen datum. The initial nonlinear receiver gain has the selected sign and is quadratic in the actual nonaxisymmetric seed amplitude. The construction and estimates survive conversion to viscosity one by the same amplitude/time scaling as pass5.

The symmetry obstruction has been removed from this initial-data family, but sustained torque, a controlled full terminal state, and repeated amplification remain unproved. The necessary energy budget for the torque does not shrink merely because (9) guarantees initial compatibility. The next decisive test must show a common nonempty amplitude/time window that both preserves full feedback and supplies enough **evolved** Reynolds stress to enter a usable successor class.

The companion checker covers the counting constants, polynomial Lipschitz estimate, tuning and parameter inequalities, central symmetry algebra, and retained half-margin. The analytic argument and the pass5 interval certificate remain explicit dependencies; this is not a proof-kernel replay.

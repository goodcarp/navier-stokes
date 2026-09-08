# A compact three-dimensional seed transfers angular momentum into a receiver

The swirl maximum principle rules out the former indefinitely rotating axisymmetric return class. This note supplies a concrete initial transfer mechanism outside that symmetry class. It also bounds its energetic cost. It does not prove that the transfer persists or regenerates the seed.

## Exact receiver identity

Let \(Jx=(-y,x,0)\), and choose a smooth nonincreasing radial weight
\(\chi_L(x)=\chi(|x|^2/L^2)\), where \(\chi\) is the unchanged pass5 cutoff.
The test field \(\phi_L=\chi_LJx\) is compact and divergence free. For any smooth unforced three-dimensional NS solution,

\[
 \frac d{dt}\int u\cdot\phi_L
 =\int u_i u_j\partial_j(\phi_L)_i+\nu\int u\cdot\Delta\phi_L .
 \tag{1}
\]

Pressure drops out by integration by parts. The antisymmetric constant matrix \(J\) also has zero quadratic contraction. Consequently the nonlinear part is exactly

\[
 Q_L(u)=\frac2{L^2}\int \chi'(|x|^2/L^2)
              (u\cdot x)(u\cdot Jx)\,dx .
 \tag{2}
\]

Decompose the cylindrical components into their angular means \(U\) and mean-zero fluctuations \(v\). Angular orthogonality gives

\[
 Q_L(U+v)=Q_L(U)+Q_L(v).
 \tag{3}
\]

The fluctuation part of the viscous receiver pairing is zero. This is a statement about the azimuthal mean of one actual field, not an ensemble of independently evolved solutions.

Write \(I_L=\int\chi_L x_1^2=\int\chi_L x_2^2\). The rotational coefficient is
\(\Omega_L=\int u\cdot\phi_L/(2I_L)\).
Its Reynolds-stress contribution to the derivative is therefore \(Q_L(v)/(2I_L)\).

## Explicit seed and strict sign

Use the compact seed of [the angular-momentum note](nonaxisymmetric-torque-budget.md):

\[
 w=\nabla\times[a(r,z)\cos(m\theta+kr)e_z],
 \quad
 w_r=-\frac{ma}{r}\sin(m\theta+kr),\quad
 w_\theta=-a_r\cos(m\theta+kr)+ka\sin(m\theta+kr),\quad w_z=0.
 \tag{4}
\]

Here \(a\) is smooth, supported away from the axis, and \(m\) is a nonzero integer. Its cylindrical component means vanish. The covariance is exactly
\(\overline{w_rw_\theta}=-mk a^2/(2r)\). Thus

\[
 \boxed{Q_L(w)=-\frac{2\pi mk}{L^2}
   \int_{\mathbb R}\int_0^\infty
      \chi'((r^2+z^2)/L^2)r^2a(r,z)^2\,dr\,dz.}
 \tag{5}
\]

This is strictly positive whenever \(mk>0\) and the nonzero envelope meets the decreasing collar of the receiver. No pressure solve or high-frequency approximation is involved.

A definite choice, using the same smooth cutoff everywhere, is

\[
\begin{gathered}
 L=\tfrac14,\quad m=4,\quad k=\frac{16}{3L},\\
 a(r,z)=L^2\chi\left(\frac{(r-3L/4)^2}{(L/16)^2}\right)
                 \chi\left(\frac{z^2}{(L/16)^2}\right).
\end{gathered}
\tag{6}
\]

The envelope is supported in
\(11L/16\le r\le13L/16,\ |z|\le L/16\), hence in the strict receiver transition
\(1/4<(r^2+z^2)/L^2<1\).
It vanishes in a neighborhood of the axis and the origin, so the Cartesian curl is globally smooth after extending it by zero.

The vector potential is even under \(x\mapsto-x\), because \(m=4\) and the axial envelope is even. Its curl is odd. It is invariant under rotations through \(\pi/2\), with the correct vector covariance, but is not axisymmetric.

For an axisymmetric base datum which equals \(Dx+Jx\) on this receiver, the exact initial coefficient for \(u_0=U_0+\varepsilon w\) is

\[
 \boxed{\Omega_L(0)=1,\qquad
  \Omega_L'(0)=2+\frac{\varepsilon^2}{2I_L}Q_L(w)>2
       \quad(\varepsilon\ne0).}
 \tag{7}
\]

The quadratic term is genuine transport into this fixed receiver. It does not increase total angular momentum of the entire isolated fluid. Adjusting a remote axisymmetric swirl amplitude, as in the [three-dimensional embedding](three-dimensional-initial-gate.md), leaves (7) unchanged.

Equation (7) does not state that the initial total RMS equals the affine base RMS: the fluctuation itself adds initial kinetic energy. Nor does the torque immediately improve every smaller receiver; a receiver whose entire test support lies inside the seed's central hole has zero initial fluctuation contribution.

## A frequency-independent energy budget

Since \(x\perp Jx\) and \(|Jx|\le |x|\),

\[
 2|(v\cdot x)(v\cdot Jx)|\le |x|^2|v|^2.
\]

Writing \(C_\chi=\sup_{s\ge0}s|\chi'(s)|<\infty\), (2) yields the exact universal bound

\[
 \boxed{|Q_L(v)|\le C_\chi\|v\|_2^2.}
 \tag{8}
\]

There is no frequency or derivative factor on the right. The same bound applies at every smooth time to the actual angular fluctuation.

Consequently, if a proposed stage needs a net amount \(D>0\) in the rotational coefficient specifically from this Reynolds-stress term at a fixed receiver, it necessarily pays

\[
 \boxed{\int_0^\tau\|v(t)\|_2^2\,dt
       \ge \frac{2I_LD}{C_\chi}.}
 \tag{9}
\]

The actual mean-flow advection and viscous terms still belong in the full receiver budget. They cannot be silently assigned to (9). The upper bound applies to the fluctuation contribution itself and does not assert a required value of \(D\) for every possible cascade.

For a normalized unit-size receiver with \(I_L\) bounded below, durations bounded above and perturbation energy tending to zero, (9) prevents a fixed positive transfer from that perturbation. Large \(m\) or \(k\) cannot circumvent the bound at fixed energy. A source whose favorable sign appears at insertion still needs a proof that its phase correlation, localization and energy persist through the required interval.

The symbolic checker verifies (1)–(5), the stated support and symmetry parameters, and the geometric quadratic bound behind (8). The integrated estimate is its direct time integral. No dynamical return is asserted.

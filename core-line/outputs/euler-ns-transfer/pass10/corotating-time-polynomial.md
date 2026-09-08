# A globally co-rotating full-NS Taylor field

## Exact coordinate transform

Fix the unchanged datum \(a_0=M+1020v+\tfrac14w_L\), viscosity \(\nu=10^{-3}\), and its whole-space NS time jets \(a_j\) defined in `time-polynomial-residual.md`. The frame rate below is a new constant \(\omega_f\); it is not the central swirl coefficient, which initially equals one.

Let \(R_\theta\) be rotation about the vertical axis and \(Jx=(-y,x,0)\). Define

\[
 (Q_t f)(x)=R_{-\omega_f t}f(R_{\omega_f t}x),\qquad
 Lf=(Jx)\cdot\nabla f-Jf.
\]

Then \(\partial_tQ_t=\omega_f LQ_t\). Cylindrical vector components of \(Lf\) are simply \(\partial_\theta f\): for the convention \(f=\sum_m f_m(r,z)e^{im\theta}\),

\[
 (Q_tf)_m=e^{im\omega_f t}f_m.
\]

The basis rotation \(-Jf\) is essential in Cartesian components. This transform rotates coordinates and vectors; it does not subtract the infinite-energy rigid-frame velocity \(\omega_f Jx\).

The whole-space Laplacian and Leray operator commute with rotations, and
\(Q_tB(a,b)=B(Q_ta,Q_tb)\). Therefore, with
\(N(u)=\nu\Delta u-B(u,u)\), the exact transformed solution \(w=Q_tu\) obeys

\[
 w_t=N(w)+\omega_f Lw. \tag{1}
\]

There is no omitted Coriolis or pressure term in (1). Introducing a relative-frame velocity by additionally subtracting a rigid rotation would be a different transformation and would require different equations.

## Actual transformed jets and exact residual

Write \(D_f=\nu\Delta+\omega_f L\). The exact transformed jets are

\[
 b_n=\sum_{j=0}^n\binom nj(\omega_f L)^{n-j}a_j.
 \tag{2}
\]

In particular, mode by mode,

\[
 b_{1,m}=a_{1,m}+im\omega_fa_{0,m},
\quad b_{2,m}=a_{2,m}+2im\omega_fa_{1,m}-m^2\omega_f^2a_{0,m}.
\]

One may compute \(b_3\) either by (2) or from

\[
 b_1=D_fb_0-B(b_0,b_0),
\]
\[
 b_2=D_fb_1-B(b_0,b_1)-B(b_1,b_0),
\]
\[
 b_3=D_fb_2-B(b_0,b_2)-B(b_2,b_0)-2B(b_1,b_1). \tag{3}
\]

These identities follow also from
\(LB(a,b)=B(La,b)+B(a,Lb)\). They use the full nonlocal NS jets, not jets obtained by keeping only angular advection.

For \(q_1=b_0+tb_1\) and \(q_2=b_0+tb_1+\tfrac12t^2b_2\), the residuals for (1) are exactly

\[
 \mathcal R_f(q_1):=q_{1,t}-D_fq_1+B(q_1,q_1)
 =-tb_2+t^2B(b_1,b_1), \tag{4}
\]
\[
 \boxed{\mathcal R_f(q_2)=
 -\tfrac12t^2b_3
 +\tfrac12t^3[B(b_1,b_2)+B(b_2,b_1)]
 +\tfrac14t^4B(b_2,b_2).} \tag{5}
\]

The physical approximations are

\[
 v_j(t)=Q_{-t}q_j(t),\qquad
 v_{j,m}(t)=e^{-im\omega_ft}q_{j,m}(t).
\]

They are smooth, exactly divergence free, full-space fields with exactly the selected initial datum. Their physical NS residual satisfies the exact identity

\[
 \boxed{v_{j,t}-N(v_j)=Q_{-t}\mathcal R_f(q_j).} \tag{6}
\]

Thus they resum a rigid angular phase without discarding any nonlinear interaction. All generated angular modes must be retained: starting with modes 0 and \(\pm4\), the first jet includes \(\pm8\) and the second includes \(\pm12\). Pressure and harmonic exterior tails belong to every coefficient in (2)–(6).

The finite angular support of each coefficient makes \(L\) a finite mode multiplier. For a general field, the coordinate factor in \((Jx)\cdot\nabla\) would require an additional domain/weighted-decay argument. No such unproved general boundedness is needed for these finite-mode jets.

## What the transform can and cannot improve

Rotations preserve angular means, axisymmetric weighted energies, and the total energy of all nonzero angular modes. For a smooth C4-equivariant field the central velocity-gradient matrix commutes with quarter-turn rotation, so its central swirl and axisymmetric strain coefficients are unchanged as well. Oddness and C4 equivariance are preserved. These are statements comparing a field to its exact rotated reconstruction.

Moreover, the mode-zero coefficients of \(q_2\) are exactly \(a_{0,0}+ta_{1,0}+\tfrac12t^2a_{2,0}\). Hence this construction has the same approximate angular mean and central core quantities as the ordinary quadratic physical Taylor field. It does not remove their genuine second-order change. The nonzero-mode amplitude evolution, and therefore its approximate fluctuation energy, can differ substantially from the ordinary polynomial.

For a manufactured pure phase \(u_m(t)=e^{-im\omega_ft}a_{0,m}\), ordinary linear and quadratic Taylor approximations multiply its energy respectively by

\[
 1+x^2,\qquad 1+x^4/4,\qquad x=m\omega_ft.
\]

The co-rotated approximation captures this exact phase with constant \(q\) and energy factor one. At \(m=4,\omega_f=1020,t=10^{-3}\), \(x=4.08\); the ordinary quadratic factor is 70.27565824. This is an exact manufactured illustration of a possible Taylor energy artifact, not a claim that the actual NS mode follows a pure phase or has this error factor.

The proposed \(\omega_f=1020\) can remove a large nearly rigid angular transport on part of the exterior packet. The base angular speed varies radially outside its flat region. Axial localization, pressure, meridional motion, shear, generated harmonics, and viscosity remain. No bound here proves that this frame rate gives a small enough residual for a useful stage.

For any chosen rotation-invariant Hilbert norm, the rate minimizing the first transformed jet alone is formally

\[
 \omega_f^*=-\frac{\operatorname{Re}\langle a_1,La_0\rangle}
 {\|La_0\|^2},\qquad \|La_0\|>0.
\]

This does not minimize the whole quadratic residual automatically; validated optimization would require bounds on the actual coefficients. The frame is a representation choice, not a new physical datum or reset.

## Sobolev norms and a safe use of the existing H4 bridge

Let the archive norm be
\(\|f\|_{D_s}^2=\sum_{|\alpha|\le s}\|\partial^\alpha f\|_2^2\), with Euclidean vector components. This derivative-sum norm is not invariant under arbitrary rotations. It must not be silently identified with a rotation-invariant Fourier Sobolev norm.

A convenient norm invariant under these vertical-axis rotations is

\[
 \|f\|_{T_s}^2=
 \sum_{h+k\le s}\sum_{j=0}^h\binom hj
 \|\partial_x^j\partial_y^{h-j}\partial_z^k f\|_2^2.
\]

Its Fourier weight is
\(\sum_{h+k\le s}(\xi_x^2+\xi_y^2)^h\xi_z^{2k}\). Consequently

\[
 \|f\|_{D_s}\le\|f\|_{T_s}
 \le\sqrt{\binom{s}{\lfloor s/2\rfloor}}\,\|f\|_{D_s},
\]

and uniformly in time and frame speed,

\[
 \|Q_tf\|_{D_4}\le\sqrt6\|f\|_{D_4},\qquad
 \|Q_tf\|_{D_5}\le\sqrt{10}\|f\|_{D_5}. \tag{7}
\]

Alternatively the standard Bessel Fourier norm with weight \((1+|\xi|^2)^s\) is invariant under all rotations. Its coefficient comparison gives \(D_4\le H^4_{\rm Fourier}\le\sqrt{24}D_4\) and \(D_5\le H^5_{\rm Fourier}\le\sqrt{60}D_5\). These are safe norm equivalences, not permission to reuse the numerical constants of a previously proved energy estimate in a different norm.

There is no need to rederive the pass9 physical \(D_4\) energy estimate just to use this approximation. If certified computations provide

\[
 \|q_2(t)\|_{D_5}\le W_5(t),\qquad
 \|\mathcal R_f(q_2)(t)\|_{D_4}\le\delta_f(t),
\]

then (6)–(7) give valid physical bounds
\(\|v_2\|_{D_5}\le\sqrt{10}W_5\),
\(\|\mathcal R(v_2)\|_{D_4}\le\sqrt6\delta_f\). The already derived bridge therefore allows

\[
 D^+E_4\le560\sqrt{10}\,W_5E_4+270E_4^2+\sqrt6\delta_f.
 \tag{8}
\]

No spurious factor proportional to \(|\omega_f|\) is introduced into this physical error estimate. One could instead derive a new estimate directly in an invariant norm, where the rotation generator is skew-adjoint, but its constants would need a separate proof.

Using the earlier product bound \(\|B(a,b)\|_{D_4}\le288\|a\|_{D_5}\|b\|_{D_5}\), coefficient norm bounds \(\|b_j\|_{D_s}\le\widetilde M_{j,s}\) give

\[
 \delta_f(t)=\tfrac12t^2\widetilde M_{3,4}
 +288t^3\widetilde M_{1,5}\widetilde M_{2,5}
 +72t^4\widetilde M_{2,5}^2,
\]

\[
 W_5(t)=\widetilde M_{0,5}+t\widetilde M_{1,5}
 +\tfrac12t^2\widetilde M_{2,5}.
\]

Large cancellations in the transformed jets should be enclosed directly if possible; estimating each summand in (2) by absolute values could discard the point of the transformation. Numerical coefficient, pressure, tail, and representation errors still need explicit budgets. The endpoint successor-class test applies to \(Q_{-T}q_2(T)\), unless that full class is itself rotation invariant. A norm bound or a pure-phase example alone certifies no positive duration and no regenerated endpoint.

# Full cylindrical Fourier representation for the selected compact datum

This note derives a representation of the full three-dimensional, unforced,
ordinary Navier–Stokes equations. It keeps the evolving mean, every retained
angular mode, axial transport and pressure. It does not prescribe an affine
background or replace the compact field by its planar slice. A finite Fourier
or radial discretization still needs an error estimate before it can certify
the corresponding whole-space solution.

## 1. Convention and full initial field

Write cylindrical components in the orthonormal basis
\((e_r,e_\theta,e_z)\). For a computational axial period \(L_z\), use

\[
 u_j(r,\theta,z)=\sum_{m\in4\mathbb Z}\sum_{\ell\in\mathbb Z}
 U_{j,m\ell}(r)e^{im\theta+i\kappa_\ell z},\qquad
 \kappa_\ell=2\pi\ell/L_z.                                      \tag{1}
\]

Coefficients are genuine Fourier coefficients: a positive angular mode and
its conjugate contribute \(2\operatorname{Re}(U_m e^{im\theta})\). The
whole-space axial transform replaces the sum over \(\ell\) by an integral,
with a consistently chosen transform normalization. All operator identities
below hold for either choice. Equation (1), when used with a finite period,
represents \(\mathbb R^2\times\mathbb T_{L_z}\), not \(\mathbb R^3\).

The selected datum is

\[
 u_0=S_\Phi+W_\psi+A v+\lambda w_L,\quad
 A=1020,\quad\lambda\in[1/8,1/4],\quad c=7/5,\quad\nu=1/1000,
 \qquad \Phi=\psi+c\eta.                                      \tag{2}
\]

Here \(s=r^2+z^2\), \(\psi(s)=\chi(s)\), and

\[
 \eta(s)=\chi(4s/25)-\chi((3s-64)/44).
\]

The unchanged smooth cutoff equals one for its argument at most \(1/4\),
zero for its argument at least one, and on the transition equals
\(\operatorname{logistic}(1/t-1/(1-t))\), \(t=(4s-1)/3\). Primes in
the following two formulas differentiate the argument \(s\):

\[
 (S_\Phi)_r=-r(\Phi+2z^2\Phi'),\qquad
 (S_\Phi)_z=2z(\Phi+r^2\Phi'),\qquad
 (W_\psi)_\theta=r(\psi+\tfrac23s\psi').                       \tag{3}
\]

Thus the compact swirl extension includes its cutoff-derivative term.
The outer profiles are

\[
 \begin{split}
 C(r)&=\chi\bigl(r^2/(63/200)^2\bigr),\qquad v=rC(r)q_v(z)e_\theta,\\
 q_v(z)&=\chi\bigl((z-4)^2/(9/20)^2\bigr)
          +\chi\bigl((z+4)^2/(9/20)^2\bigr),\\
 e(r)&=\chi\bigl((r-7/50)^2/(1/10)^2\bigr),\\
 q_s(z)&=\chi\bigl((z-4)^2/(3/5)^2\bigr)
          +\chi\bigl((z+4)^2/(3/5)^2\bigr),\\
 w_L&=\nabla\times[e(r)q_s(z)\cos(4\theta-20r)e_z].
 \end{split}                                                  \tag{4}
\]

Before taking the axial transform, the initial mode zero is (3) plus
\(A v\). The positive mode four is exactly

\[
 \Psi_4=\frac{\lambda}{2}e(r)q_s(z)e^{-20ir},\qquad
 U_{r,4}=\frac{4i}{r}\Psi_4,\quad
 U_{\theta,4}=-\partial_r\Psi_4,\quad U_{z,4}=0.                \tag{5}
\]

The negative coefficient is its componentwise conjugate. The initial field
has only modes \(0,\pm4\); its initial pressure source and convection have
only \(0,\pm4,\pm8\). Later times generate further multiples of four.
Keeping only the initial modes is not the full NS evolution.

The datum is smooth, divergence-free and supported inside the ball of radius
six. The nonaxisymmetric seed vanishes on an axis neighborhood. The complete
datum, including the compact core and both axial packets, is represented by
`initial_full_field.py`. Its coefficients and the normalization in (5) were
independently checked against the formulas here.

## 2. Exact cylindrical equations and mode convolution

For a single angular mode, with rows denoting output components and columns
denoting physical differentiation directions, the velocity gradient is

\[
 \mathcal G_m(U)=
 \begin{pmatrix}
 \partial_r U_r&(imU_r-U_\theta)/r&\partial_z U_r\\
 \partial_r U_\theta&(imU_\theta+U_r)/r&\partial_z U_\theta\\
 \partial_r U_z&imU_z/r&\partial_z U_z
 \end{pmatrix}.                                               \tag{6}
\]

Writing axial Fourier transforms after products, or including their
convolution indices explicitly, gives the exact identities

\[
 N_m=((u\cdot\nabla)u)_m=\sum_{a+b=m}\mathcal G_a(U_a)U_b,
 \qquad
 g_m=(-\Delta p)_m=\sum_{a+b=m}
        \operatorname{tr}(\mathcal G_a(U_a)\mathcal G_b(U_b)).  \tag{7}
\]

The pressure source is the trace of the matrix square, not the squared
Frobenius norm. For a divergence-free field it is also
\(g=\nabla\cdot N\). Two smooth fields with disjoint supports do have zero
bilinear initial pressure source and hence zero decaying bilinear initial
pressure. Their individual pressure fields still act across the supports
and create later interactions. Support separation therefore does not allow
one to omit those individual pressure fields or to discard cross terms
after the initial supports cease to be disjoint.

Equivalently, let

\[
 T_{j,m}=\sum_{a+b=m}
 [U_{r,a}\partial_r U_{j,b}
       +ib\,U_{\theta,a}U_{j,b}/r
       +U_{z,a}\partial_z U_{j,b}].
\]

Then \(N_r=T_r-\sum U_{\theta,a}U_{\theta,b}/r\),
\(N_\theta=T_\theta+\sum U_{r,a}U_{\theta,b}/r\), and
\(N_z=T_z\). In (7) the roles of the two summed indices differ from this
last display; exchanging the dummy indices gives exactly the same result.

For the initial real mode-zero field and positive coefficient \(U_4\),

\[
 \begin{split}
 g_0&=\operatorname{tr}(\mathcal G_0^2)
             +2\operatorname{tr}(\mathcal G_4\overline{\mathcal G_4}),\\
 g_4&=2\operatorname{tr}(\mathcal G_0\mathcal G_4),\qquad
 g_8=\operatorname{tr}(\mathcal G_4^2).
 \end{split}                                                  \tag{8}
\]

The zero-mode expression is real. Its second trace need not be positive.

## 3. Helical variables, diffusion and projection

Suppress the two mode indices and set

\[
 U_+=U_r+iU_\theta,\qquad U_-=U_r-iU_\theta,\qquad
 L_{n,\kappa}=\partial_{rr}+r^{-1}\partial_r-n^2/r^2-\kappa^2.
\]

Since \(u_x+i u_y=e^{i\theta}(u_r+i u_\theta)\), the Cartesian angular
orders of \(U_+\) and \(U_-\) are respectively \(m+1\) and \(m-1\).
The vector Laplacian is therefore

\[
 \mathcal L U=(L_{m+1,\kappa}U_+,L_{m-1,\kappa}U_-,L_{m,\kappa}U_z).
                                                               \tag{9}
\]

In the original components this says

\[
 (\Delta U)_r=L_{m,\kappa}U_r-U_r/r^2-2imU_\theta/r^2,
 \qquad
 (\Delta U)_\theta=L_{m,\kappa}U_\theta-U_\theta/r^2+2imU_r/r^2.
\]

Define the pressure gradient and velocity divergence in helical variables:

\[
 Gp=((\partial_r-m/r)p,(\partial_r+m/r)p,i\kappa p),             \tag{10}
\]

\[
 DU=\tfrac12[(\partial_r+(m+1)/r)U_+
                 +(\partial_r+(1-m)/r)U_-]+i\kappa U_z.         \tag{11}
\]

Their exact identities are

\[
 DG=L_{m,\kappa},\qquad
 \mathcal L G=G L_{m,\kappa},\qquad
 D\mathcal L=L_{m,\kappa}D.                                  \tag{12}
\]

The correct inner product on a radial mode is

\[
 \langle U,V\rangle_H=\int_0^\infty r
 [\tfrac12\overline{U_+}V_++\tfrac12\overline{U_-}V_-
                    +\overline{U_z}V_z]dr.                    \tag{13}
\]

With regularity at the axis and suitable decay, integration by parts gives
\(G^*=-D\). For finite radial endpoints it produces the boundary term
\([r\overline p\,U_r]\); it cannot simply be discarded.
Consequently the Leray projection is

\[
 \mathbb P=I-G L_{m,\kappa}^{-1}D
           =I+G(-L_{m,\kappa})^{-1}D.                         \tag{14}
\]

For a tentative field \(W\), solve \(L\phi=DW\) and set
\(\mathbb PW=W-G\phi\). For the physical pressure, solve
\(-Lp=DN\); the actual NS equation is

\[
 \partial_t U=\nu\mathcal L U-N-Gp
              =\nu\mathcal L U-\mathbb P N,\qquad DU=0.       \tag{15}
\]

The sign of a time-step correction potential and the sign of the physical
pressure Poisson equation must be kept distinct.

One may compute nonlinear helical components directly:

\[
 N_{\pm,m}=\sum_{a+b=m}
 [U_{r,a}\partial_r U_{\pm,b}
      +i(b\pm1)U_{\theta,a}U_{\pm,b}/r
      +U_{z,a}\partial_z U_{\pm,b}].                          \tag{16}
\]

For \(N_z\), the angular factor is \(ib\). The extra \(\pm1\) is the
derivative of the rotating vector basis, not an additional physical mode.

## 4. Axis, reality, symmetry and mapped radius

A smooth Cartesian vector field has the following regular form near the
axis, with smooth coefficient functions of \(r^2,z\):

\[
 U_+=r^{|m+1|}a_+(r^2,z),\quad
 U_-=r^{|m-1|}a_-(r^2,z),\quad
 U_z=r^{|m|}a_z(r^2,z),\quad p_m=r^{|m|}a_p(r^2,z).             \tag{17}
\]

For \(m=0\), \(U_\pm=O(r)\) and \(U_z,p\) are even. For \(m=4\),
\(U_+=O(r^5)\), \(U_-=O(r^3)\), and \(U_z,p=O(r^4)\).
Parity alone does not impose the required higher vanishing orders.
When \(f=r^n a\), \(n\ge0\),

\[
 L_{n,\kappa}f=r^n[a''+(2n+1)a'/r-\kappa^2a],                \tag{18}
\]

with \(a\) even and \(a'(0)=0\). This is one way to represent the removable
axis singularities. The ladder operators in (10) preserve the right orders;
for example \((\partial_r-m/r)r^m a(r^2)=O(r^{m+1})\) for \(m\ge0\).

Reality in cylindrical variables is ordinary conjugate symmetry. In helical
variables it instead reads

\[
 U_{-,-m,-\ell}=\overline{U_{+,m,\ell}},\qquad
 U_{z,-m,-\ell}=\overline{U_{z,m,\ell}}.                      \tag{19}
\]

Imposing separate scalar reality on each helical component is incorrect.
Fourfold rotational equivariance closes the cylindrical modes under
convolution. If the reduced angular coordinate is \(\varphi=4\theta\),
replace \(\partial_\theta\) by \(4\partial_\varphi\), including the basis
terms in (16). The present odd datum also has \(u_r,u_\theta\) even and
\(u_z\) odd in \(z\); the smooth evolution preserves this symmetry.

For a mapped grid \(r=R(\xi)\), \(J=R'(\xi)>0\),

\[
 L_{n,\kappa}=
 \frac1{rJ}\partial_\xi\left(\frac rJ\partial_\xi\right)
       -\frac{n^2}{r^2}-\kappa^2,
 \qquad r\,dr=rJ\,d\xi.                                    \tag{20}
\]

Ladders replace \(\partial_r\) by \(J^{-1}\partial_\xi\).
A weighted discretization should enforce a compatible discrete adjoint pair
\(G_h^*=-D_h\), and use \(L_h=D_hG_h\) for projection. Independently
discretizing each differential expression need not preserve (12). If the
discrete viscous step does not preserve \(D_hU=0\), a coupled Stokes step or
an appropriate projection is needed. A subsequent projection alone does not
prove the desired temporal accuracy; the splitting error remains to be
measured or bounded.

## 5. Free-space pressure and the exterior

For \((-L_{m,\kappa})p=f\), \(n=|m|\), and \(k=|\kappa|>0\), the
regular-at-zero, decaying radial Green solution is

\[
 p(r)=K_n(kr)\int_0^r I_n(ks)f(s)s\,ds
       +I_n(kr)\int_r^\infty K_n(ks)f(s)s\,ds.                \tag{21}
\]

There is no extra factor \(k\): the radial Wronskian, in the order
\((I_n(kr),K_n(kr))\), is \(-1/r\). The scalar equation and normalization
agree with the [NIST modified Bessel equation](https://dlmf.nist.gov/10.25.E1)
and [Wronskian identity](https://dlmf.nist.gov/10.28.E2), with the latter's
order reversed. When \(k=0\), \(n>0\), the corresponding formula is

\[
 p(r)=\frac1{2n}\left[
 r^{-n}\int_0^r s^{n+1}f(s)\,ds
       +r^n\int_r^\infty s^{1-n}f(s)\,ds\right].             \tag{22}
\]

For \(k=n=0\),

\[
 p'(r)=-\frac1r\int_0^r sf(s)\,ds.                          \tag{23}
\]

A decaying zero-mode pressure requires \(\int_0^\infty rf(r)dr=0\),
with a sufficiently decaying source, and its additive gauge can then be
fixed. This compatibility follows from the divergence source and vanishing
boundary flux in the exact problem; it is a useful discrete check. In this
same periodic-axial, angular-zero, axial-zero mode, regular decaying
incompressibility forces \(U_{r,00}=0\).

A homogeneous Bessel exterior condition at a finite radial boundary is
exact only if the omitted source truly vanishes outside that boundary.
Although (2) is compact initially, the actual velocity immediately develops
a nonlocal tail. Initially outside its support,

\[
 u_t(0,x)=-\nabla p_0(x),\qquad
 p_0=\partial_i\partial_j(4\pi|\cdot|)^{-1}*(u_{0i}u_{0j}).   \tag{24}
\]

Generically \(p_0=O(|x|^{-3})\), \(\nabla p_0=O(|x|^{-4})\);
additional moment cancellation can improve these rates. Thus a pressure
boundary justified by the initial support does not remain justified at
positive time. Mapping the radial domain to infinity avoids imposing an
artificial finite radius at the continuum level, but still requires
approximation and tail control near the mapped endpoint.

Taking \(L_z>12\) makes the initial datum fit inside one axial cell. It does
not make the periodic pressure equal to the whole-space pressure: periodic
images affect the interior immediately. Enlarging the period and testing
convergence is a diagnostic, not a rigorous whole-space error bound. A
whole-space certificate needs a free-space pressure treatment or a global
smooth divergence-free reconstruction with controlled exterior residual.
Multiplying a periodic velocity by an axial cutoff usually breaks
incompressibility; correcting that defect and its pressure effect is part of
the reconstruction, not an optional afterthought.

## 6. Energy, divergence and spectral checks

For normalization (1), Parseval gives

\[
 E=\tfrac12\|u\|_2^2
 =\pi L_z\sum_{m,\ell}\int_0^\infty r
 [\tfrac12|U_+|^2+\tfrac12|U_-|^2+|U_z|^2]dr.                \tag{25}
\]

The full Cartesian gradient dissipation is

\[
 \begin{split}
 \mathcal D=2\pi L_z\sum_{m,\ell}\int_0^\infty r\{&
 \tfrac12[|U_+'|^2+((m+1)^2/r^2+\kappa^2)|U_+|^2]\\
 &+\tfrac12[|U_-'|^2+((m-1)^2/r^2+\kappa^2)|U_-|^2]\\
 &+|U_z'|^2+(m^2/r^2+\kappa^2)|U_z|^2\}\,dr.
 \end{split}                                                  \tag{26}
\]

For smooth solutions with the stated boundary behavior,
\(E'=-\nu\mathcal D\). Both pressure work and the total convective work
vanish. Individual mean or fluctuation energies need not decrease. With
positive modes only, double their contribution to these sums; retain the
zero mode once. If integrating a quarter-angle wedge, multiply that physical
integral by four to recover the full energy.

The curl can provide an additional independent diagnostic:

\[
 \begin{split}
 \omega_+&=-\kappa U_+-i(\partial_r-m/r)U_z,\\
 \omega_-&=\kappa U_-+i(\partial_r+m/r)U_z,\\
 \omega_z&=\frac1{2i}[(\partial_r+(m+1)/r)U_+
                         -(\partial_r+(1-m)/r)U_-].
 \end{split}                                                  \tag{27}
\]

It satisfies \(D\operatorname{curl}=0\),
\(\operatorname{curl}G=0\), and
\(\operatorname{curl}\operatorname{curl}=GD-\mathcal L\).

An implementation should report these distinct checks:

1. Initial reconstruction, divergence, Fourier reality and axis vanishing
   orders, including the factor \(1/2\) in (5).
2. Discrete \(DG=L\), weighted adjointness and diffusion/divergence
   intertwining, measured independently on regular test modes.
3. Projection divergence, idempotence, self-adjointness and the orthogonal
   identity \(\|W\|^2=\|\mathbb PW\|^2+\|G\phi\|^2\).
4. Full energy versus integrated viscous loss, together with nonlinear and
   pressure work, unresolved-mode content and radial/axial boundary errors.
5. Pressure from \(\nabla\cdot N\) versus pressure from
   \(\operatorname{tr}((\nabla u)^2)\); their disagreement detects a
   divergence or differentiation error, not a discretionary pressure choice.

Exact Fourier Galerkin convolution preserves the continuous nonlinear energy
cancellation before radial and temporal discretization. A pseudospectral
implementation should use sufficient padding for quadratic products in both
angular and axial directions (for example the usual consistent 3/2 padding),
or another verified dealiasing rule. Radial products and quadrature need
their own accuracy check. Observed energy stability alone cannot rule out a
wrong pressure solver, omitted axial terms or unresolved modes.

`verify_cylindrical_NS_representation.py` independently checks the differential
identities, helical signs, nonlinear basis terms, mapped operator, initial
streamfunction divergence and selected exact regular polynomial modes. It is
an algebra check, not a convergence result or a validated NS time integrator.

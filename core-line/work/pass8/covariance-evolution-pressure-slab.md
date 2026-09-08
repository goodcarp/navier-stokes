# A whole-space radial pressure-gradient error smaller than five

For the exact unit mixed pressure \(p_m\) in the covariance evolution
note and the Pass7 horizontal trial \(p_h=P(r)F(z)\), put
\(q=p_m-p_h\). The stored Pass7 residual certificate implies the
uniform pointwise bound
\[
 \boxed{|\partial_rq(r,\pm4)|<5
 \quad\text{for }0\le r\le\frac{63}{200}\sqrt{\frac58}.}              \tag{1}
\]
This controls the complex \(m=4\) amplitude, including its phase.
It is a whole-space estimate with both packets and radial tails;
no finite computational boundary is involved.

## 1. The harmonic slab and the global energy already certified

The source and trial are exactly
\[
 -\Delta_m p_m=S(r)F(z),\qquad
 -\Delta_m(PF)=SF-PF'',\qquad
 -\Delta_m q=P(r)F''(z).
\]
The original axial cutoffs are flat for \(|z-4|\le9/40\); the negative
packet is disjoint from this interval. Hence \(q(r,z)e^{im\theta}\)
is harmonic in the entire infinite slab \(|z-4|<h=9/40\).
The same statement holds about \(z=-4\).

Use the physical real field \(\operatorname{Re}(q e^{im\theta})\)
for all three-dimensional norms. The whole-space angular energy
estimate from Pass7 gives
\[
 \|\nabla q_{\rm real}\|_2\le\|r(PF'')_{\rm real}\|_2/m
 =n_e/m.
\]
The saved outward upper endpoint has \(n_e<21.230615<85/4\).
The exact source is smooth, the trial is regular at the axis and
decays as \(r^{-4}\) radially, and the whole-space mixed pressure
has the corresponding decaying nonzero angular mode. The Hankel
and energy operations below therefore apply by smooth approximation
(in particular these fields are in \(L^2\) and have finite energy).

## 2. A harmonic-slab evaluation inequality with the exact normalization

For a pure complex angular mode \(q_m(r,z)\), use the Hankel transform
\[
 \widehat q(k,z)=\int_0^\infty r J_m(kr)q_m(r,z)\,dr,\qquad
 q_m(r,z)=\int_0^\infty kJ_m(kr)\widehat q(k,z)\,dk.
\]
The physical energy in a slab, with \(s=z-4\), is
\[
 \|\nabla q_{\rm real}\|_{L^2(|s|<h)}^2
 =\pi\int_0^\infty k\,dk\int_{-h}^h
      (|\widehat q_s|^2+k^2|\widehat q|^2)\,ds .
\]
Harmonicity gives
\(\widehat q(k,s)=a(k)\cosh(ks)+b(k)\sinh(ks)\).
Integration in the symmetric slab makes the cross term odd, so
\[
 \int_{-h}^h(|\widehat q_s|^2+k^2|\widehat q|^2)\,ds
 =k\sinh(2kh)(|a|^2+|b|^2).
\]
In particular, the global energy is at least
\(\pi\int k^2\sinh(2kh)|a(k)|^2\,dk\).
Hankel inversion and Cauchy–Schwarz now give the pointwise inequality
\[
 \boxed{|\partial_rq_m(r,4)|^2
 \le\frac{\|\nabla q_{\rm real}\|_2^2}{\pi}
 \int_0^\infty\frac{k^2|J_m'(kr)|^2}{\sinh(2hk)}\,dk.}               \tag{2}
\]
The factor is \(1/\pi\), not \(1/(2\pi)\); integrating the square of
the physical real angular mode contributes \(\pi\). Complex
Cauchy–Schwarz bounds the full complex amplitude and requires no
additional factor of two.

## 3. A global positive-weight bound on the Bessel derivative

For integer \(m\ge1\) and \(x\ge0\),
\[
 \boxed{|J_m'(x)|\le
       \frac{x^{m-1}}{2^m(m-1)!}.}                                  \tag{3}
\]
This sharper bound does not follow from a triangle estimate on
\((J_{m-1}-J_{m+1})/2\); it has a positive-weight proof.
In the Poisson integral write
\[
 J_m(x)=c_m x^m\int_{-1}^1 e^{ixt}w(t)\,dt,\qquad
 w(t)=(1-t^2)^{m-1/2}.
\]
Differentiate and integrate the term containing \(itxe^{ixt}\)
by parts. Since \(w(\pm1)=0\),
\[
 J_m'(x)=c_m x^{m-1}\int_{-1}^1 e^{ixt}
       [(m-1)w(t)-t w'(t)]\,dt .
\]
The new weight is
\[
 (m-1)(1-t^2)^{m-1/2}
 +(2m-1)t^2(1-t^2)^{m-3/2}\ge0 .
\]
It is integrable also when \(m=1\), and its integral equals
\(m\int w\). Taking its absolute integral and using the normalization
of the Poisson representation proves (3).

Substitution into (2) yields
\[
 |\partial_rq_m|^2\le
 \frac{n_e^2}{m^2\pi}
 \frac{r^{2m-2}}{2^{2m}((m-1)!)^2}
 \frac{2(2m)!}{(2h)^{2m+1}}
 \sum_{\substack{j\ge1\\j\ {\rm odd}}}j^{-(2m+1)} .                  \tag{4}
\]
The integral identity follows directly from the positive expansion
\(1/\sinh x=2\sum_{n\ge0}e^{-(2n+1)x}\); termwise integration is
justified by nonnegativity.

## 4. A rational bound for the actual datum

For \(m=4\), the scalar coefficient in (4) simplifies to
\[
 |\partial_rq_4|^2
 \le\frac{n_e^2}{16\pi}\frac{35}{4}
       \frac{r^6}{(2h)^9}
       \sum_{j\ {\rm odd}\ge1}j^{-9}.
\]
A rational integral comparison gives
\[
 \sum_{j\ {\rm odd}\ge1}j^{-9}
 \le1+3^{-9}+\frac12\int_3^\infty x^{-9}\,dx
 =1+\frac{19}{16\cdot3^9}<\frac{10001}{10000}.
\]
Use \(n_e<85/4\), \(\pi>25/8\),
\(2h=9/20\), and \(r^2\le(63/200)^2(5/8)\). Then
\[
 |\partial_rq_4|^2<
 \frac{(85/4)^2}{16}\frac8{25}\frac{35}{4}
 \frac{[(63/200)^2(5/8)]^3}{(9/20)^9}
 \frac{10001}{10000}
 =\frac{2380277273927}{95551488000}<25,
\]
which proves (1).

The horizontal trial has \(F=1\) at the target circles. Therefore
\[
 \partial_rp_m(r_*,\pm4)=P'(r_*)+\varepsilon,\qquad
 |\varepsilon|<5 .
\]
The derivative-reduced Green formula from Pass7 differentiates to
\[
 P'=-2C_0 b_s'
   +m(m-1)r^{-m-1}I_+
   -m(m+1)r^{m-1}I_-,
\]
\[
 I_+(r)=\int_0^r t^m C_0'(t)b_s(t)\,dt,\qquad
 I_-(r)=\int_r^\infty t^{-m}C_0'(t)b_s(t)\,dt.
\]
Every local \(C_0'b_s\) term cancels. On the actual seed plateau this is
\(P'=-40iC_0e^{20ir}+12r^{-5}I_+-20r^3I_-\).
Thus the existing radial moment ranges can directly enclose the trial
gradient; differentiating quadrature values is unnecessary.

The error contributed to the exact torque derivative (equation (6)
in the covariance note) is at most
\[
 \frac{5\lambda^2A}{2}
 \sqrt{(k^2r_*+m^2/r_*)^2+k^2}.
\]
This can be combined with a rigorous radial Green-derivative enclosure
and cutoff-jet bounds. It does not itself assign a sign to the second
mean jet or justify a finite time interval.

The independent checker replays the saved residual endpoint, the
positive-weight identity at \(m=4\), the exact slab-energy coefficient,
and the final rational inequality. A separate agent independently
checked the general positive-weight argument and the \(\pi\)
normalization in (2).

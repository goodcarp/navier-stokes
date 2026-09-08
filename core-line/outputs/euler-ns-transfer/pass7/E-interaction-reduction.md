# Exact reduction and a whole-space pressure certificate for the outer E coefficient

This note concerns only the coefficient \(E=\mathcal L(v;w)\) defined in the
Pass6 outer-pressure compatibility test. It retains the complete Leray
projection. The formulas below do not on their own certify a numerical value
or sign for \(E\), or any finite-time return.

## 1. Conventions and the complete pressure interaction

Let \(N(x)=1/(4\pi|x|)\), \(Q_{ij}=-\partial_{zzij}N\), and
\[
 \Pi(a,b)=\partial_{zz}\{N*\operatorname{tr}(\nabla a\nabla b)\}(0).
\]
For the compact outer fields used here, integration by parts gives
\(\Pi(a,b)=-\int Q_{ij}a_i b_j\). The Leray-projected bilinear notation in
Pass6 is \(B(a,b)=\mathbb P(a\cdot\nabla b)\). In particular,
\[
 E=-2\{\Pi(v,B(w,w))+\Pi(w,B(v,w))+\Pi(w,B(w,v))\}.             \tag{1}
\]

Write \(v=V(r,z)e_\theta\) and
\[
 w=\operatorname{Re}\{W(r,z)e^{im\theta}\},\qquad
 W_r=\frac{im}{r}\mathfrak b,\quad W_\theta=-\partial_r\mathfrak b,\quad
 W_z=0,\qquad \mathfrak b=a_s(r,z)e^{ikr}.
\]
The actual profiles have \(m=4,\ k=20\), \(a_s=e(r)q_s(z)\), and
\(V=V_0(r)q_v(z)\). The field \(Qv=Q_{\theta\theta}V e_\theta\) is
divergence free. Consequently the pressure gradient in \(B(w,w)\) pairs
to zero with \(Qv\). This removes exactly that pressure placement.

The two other placements combine into
\[
 C=(v\cdot\nabla)w+(w\cdot\nabla)v,\qquad
 \mathbb PC=C+\nabla p_{\rm mix},\qquad -\Delta p_{\rm mix}=\operatorname{div}C.
\]
Their complex amplitudes are
\[
 C_r=\frac{imV}{r}W_r-\frac{2V}{r}W_\theta,\qquad
 C_\theta=\frac{imV}{r}W_\theta+(V_r+V/r)W_r,\qquad C_z=0.
\]
Thus the full mixed pressure mode solves
\[
 -\Delta_m p_{\rm mix}=s_m,\qquad
 \Delta_m=\partial_{rr}+r^{-1}\partial_r+\partial_{zz}-m^2r^{-2},
\]
\[
 \boxed{s_m=\frac2r\left[\partial_r(V\partial_r\mathfrak b)
                          -\frac{m^2}{r}V_r\mathfrak b\right].}       \tag{2}
\]
The absence of \(z\)-derivatives in this source is a property of the two
horizontal velocities; it does not remove the \(z\)-direction from the
Poisson operator.

Since \(\operatorname{div}Q=0\) away from the origin and
\(\operatorname{div}w=0\), the amplitude of \(\operatorname{div}(Qw)\) is
\[
 \boxed{d_m=(Q_{rr}-Q_{\theta\theta})\partial_rW_r
                         +Q_{rz}\partial_zW_r.}                    \tag{3}
\]
The second term retains the vertical component of \(Qw\).

## 2. Local cancellation and the pressure term that remains

Angular averaging followed by one radial integration by parts gives
\[
 \boxed{E=E_{\rm loc}+E_{\rm press},}
\]
\[
 \boxed{E_{\rm loc}=6\pi mk
       \int_{\mathbb R}\int_0^\infty
                a_s(r,z)^2V(r,z)\partial_rQ_{\theta\theta}\,dr\,dz,}  \tag{4}
\]
\[
 \boxed{E_{\rm press}=-2\pi\operatorname{Re}
        \int_{\mathbb R}\int_0^\infty
                r\,\overline{d_m}\,p_{\rm mix}\,dr\,dz.}             \tag{5}
\]
For reproducibility, the three local averages used here are
\[
 \langle[(w\cdot\nabla)w]_\theta\rangle
   =-\frac{mk}{r}a_s(a_s)_r-\frac{mk}{2r^2}a_s^2,
\]
\[
 \langle w_r C_r\rangle=\frac{mkV}{r^2}a_s^2,\qquad
 \langle w_\theta C_\theta\rangle
   =-\frac{mk}{2r}a_s^2(V_r+V/r).
\]
The kernel identities are
\[
 Q_{\theta\theta}=\frac{3(4z^2-r^2)}
                        {4\pi(r^2+z^2)^{7/2}},\quad
 \partial_rQ_{\theta\theta}
       =\frac{15r(r^2-6z^2)}{4\pi(r^2+z^2)^{9/2}},
 \quad Q_{rr}-Q_{\theta\theta}=r\partial_rQ_{\theta\theta}.           \tag{6}
\]
For the actual nonnegative profiles, \(mk>0\) and \(r^2<6z^2\) throughout
their overlap. Thus \(E_{\rm loc}<0\). The sign of the remaining pressure
pairing cannot be dropped or inferred from this.

A meridional reflection maps the polar fields \(v\mapsto-v\) and
\(w_k\mapsto-w_{-k}\). Invariance of the axial pressure functional therefore
gives \(E(k)=-E(-k)\). This forces \(E(0)=0\), but does not force zero for
the chosen \(k=20\). Reflection in \(z\) also leaves an even integrand.

## 3. A rigorous residual estimate with the full pressure tail

All norms in this paragraph are physical three-dimensional norms. For a
real \(m\)-mode \(h=\operatorname{Re}(h_m e^{im\theta})\),
\[
 \|h\|_2^2=\pi\int r|h_m|^2\,dr\,dz.
\]
Let \(p_h\) be any regular whole-space finite-energy \(m\)-mode trial with
residual \(e=s+\Delta p_h\), so that the pressure error \(p-p_h\) solves
\(-\Delta(p-p_h)=e\). Assume \(r e\in L^2\).
The angular derivative and the energy identity imply
\[
 m\|(p-p_h)/r\|_2\le\|\nabla(p-p_h)\|_2,\qquad
 \|\nabla(p-p_h)\|_2\le\|r e\|_2/m.
\]
These identities follow first by approximation with finite-energy smooth
fields and then by the energy solution; compact support of the pressure
itself is not required. Consequently,
\[
 \boxed{|E_{\rm press}-E_{\rm press}[p_h]|
          \le\frac{2}{m^2}\|r d\|_2\|r e\|_2.}                     \tag{7}
\]
In complex-amplitude form the right side is
\[
 \frac{2\pi}{m^2}
   \left(\int r^3|d_m|^2\right)^{1/2}
   \left(\int r^3|e_m|^2\right)^{1/2}.
\]
A finite-box Dirichlet solution extended by zero does not automatically
qualify: its boundary distribution residual would also have to be bounded.
The following trial has no such boundary.

One alternative exact identity is
\(\operatorname{div}C=2\operatorname{div}[(w\cdot\nabla)v]\). If
\(\mathbb P_{\rm grad}\) is the orthogonal gradient projection, then
\[
 E_{\rm press}
 =-4\langle\mathbb P_{\rm grad}(Qw),(w\cdot\nabla)v\rangle,\qquad
 |E_{\rm press}|\le\frac4m\|rd\|_2\|(w\cdot\nabla)v\|_2.
\]
This is another complete-pressure bound, not an asserted favorable sign.

## 4. Exact separable source and a derivative-reduced radial Green formula

Set \(F(z)=q_s(z)q_v(z)\) and \(b_s(r)=e(r)e^{ikr}\). Equation (2) is
exactly \(s_m=S(r)F(z)\), where
\[
 S(r)=\frac2r\left[(V_0 b_s')'-\frac{m^2}{r}V_0'b_s\right].
 \tag{8}
\]
Let \(C_0(r)=V_0(r)/r=\chi(r^2/a_v^2)\), with \(a_v=63/200\), and
\(L_m=\partial_{rr}+r^{-1}\partial_r-m^2r^{-2}\).
The regular, decaying solution to \(-L_mP=S\) is
\[
 P(r)=\frac1{2m}\left[
 r^{-m}\int_0^r t^{m+1}S(t)\,dt+
 r^m\int_r^\infty t^{1-m}S(t)\,dt\right].                            \tag{9}
\]
Twice integrating each integral by parts, with flat compact profile
endpoints, gives the useful formula
\[
 \boxed{\begin{aligned}
 P(r)={}&-2C_0(r)b_s(r)
 -(m-1)r^{-m}\int_0^r t^m C_0'(t)b_s(t)\,dt\\
 &-(m+1)r^m\int_r^\infty t^{-m} C_0'(t)b_s(t)\,dt .
 \end{aligned}}                                                     \tag{10}
\]
Only one cutoff derivative remains. In detail, the integrals in (9) equal
\[
 2r^mV_0 b_s'-2m r^{m-1}V_0b_s
             -2m(m-1)\int_0^r t^m C_0'b_s\,dt
\]
and
\[
 -2r^{-m}V_0 b_s'-2m r^{-m-1}V_0b_s
             -2m(m+1)\int_r^\infty t^{-m}C_0'b_s\,dt,
\]
respectively. The \(b_s'\) boundary terms cancel in (9).
Direct symbolic substitution also verifies \(-L_mP=S\).

The globally regular trial \(p_h(r,z)=P(r)F(z)\) obeys
\[
 -\Delta_m p_h=SF-PF'',\qquad
 \boxed{e_m=s_m+\Delta_m p_h=P(r)F''(z).}                            \tag{11}
\]
Therefore its entire residual norm separates:
\[
 \boxed{\|r e\|_2^2=
 \pi\left(\int_0^\infty r^3|P(r)|^2\,dr\right)
       \left(\int_{\mathbb R}|F''(z)|^2\,dz\right).}                 \tag{12}
\]
Here \(F\) is smooth and compactly supported. The radial field has exact
tails, which must be included in (12). With \(r_*=7/100\),
\[
 P(r)=c_{\rm in}r^m\quad(0<r\le r_*),\qquad
 c_{\rm in}=-(m+1)\int_0^\infty t^{-m}C_0'(t)b_s(t)\,dt,
\]
\[
 P(r)=c_{\rm out}r^{-m}\quad(r\ge a_v),\qquad
 c_{\rm out}=-(m-1)\int_0^\infty t^mC_0'(t)b_s(t)\,dt .
\]
In particular,
\[
 \int_0^{r_*}r^3|P|^2\,dr
 =\frac{|c_{\rm in}|^2r_*^{2m+4}}{2m+4},\qquad
 \int_{a_v}^\infty r^3|P|^2\,dr
 =\frac{|c_{\rm out}|^2a_v^{4-2m}}{2m-4}.                           \tag{13}
\]
The exterior weighted norm is finite for \(m>2\), including \(m=4\).
Smooth flat transitions and these regular/decaying tails make the trial
admissible in the whole-space estimate (7).

For a certified upper bound, one may combine \(E_{\rm loc}\) with
\[
 2\pi\int r|d_m||P F|\,dr\,dz
 +\frac2{m^2}\|rd\|_2\|re\|_2,
\]
or rigorously enclose the signed trial pairing for a sharper result.
The finite radial integrals, cutoff integrals, and norm inputs must still
receive actual outward enclosures. Floating-point quadrature by itself
does not provide those enclosures.

The standalone checker verifies the angular local algebra, complete
mixed source, kernel/divergence identities, exact radial Green formula,
residual sign, and both tail antiderivatives. It does not certify any
profile quadrature or the remaining \(D_\nu\) coefficient.

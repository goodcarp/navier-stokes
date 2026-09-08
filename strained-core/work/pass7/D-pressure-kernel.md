# A favorable full pressure kernel and a certified bound for D_nu

For the unit outer chiral seed `w` from pass 6 and the fixed background `M=S_Phi+W_psi`, the coefficient in the pressure compatibility polynomial satisfies

\[
\boxed{D_\nu=\mathcal L(M;w)+2\nu\Pi(w,\Delta w)
<4-\frac{96}{5}C_w<4,\qquad \nu=\frac1{1000}.}
\tag{1}
\]

Here `C_w=-Pi(w,w)>0` is the exact seed pressure coefficient. This bound is below the earlier sufficient target `5093339/210000`. It does not assert that the separate coefficient `E=L(v;w)` vanishes; that is a distinct calculation. The proof retains the entire pressure acceleration through a whole-space adjoint and then bounds viscosity by a one-dimensional interval integral.

Use the unsymmetrized Leray nonlinearity `B(a,b)=P(a dot grad b)`, and

\[
\mathcal L(M;w)=-2\{\Pi(M,B(w,w))+
\Pi(w,B(M,w)+B(w,M))\}.
\]

## 1. Local mixed acceleration and the full quadratic kernel

Write `P(x)=(x,y,-2z)` for the **unit affine** pump on the seed support; the actual annular field is `S_eta`, which equals `P` there. The amplitude is `c=7/5`. The seed is horizontal, `w_z=0`, and has zero horizontal divergence. The compact mixed advection therefore obeys

\[
\operatorname{div}[(cP\cdot\nabla)w+cw]=0.
\]

Consequently its Leray projection is unchanged:

\[
B(M,w)+B(w,M)=c[(P\cdot\nabla)w+w].
\tag{2}
\]

The core part of `M` is disjoint from `w`. Formula (2) does not discard `B(w,w)`'s nonlocal pressure component.

Set `Q_ij=-partial_zzij N`, `N=1/(4 pi R)`, `R=|x|`. For solenoidal fields supported away from the evaluation point, `Pi(a,b)=-integral Q_ij a_i b_j`. Apply this identity to the mixed term and to the annular part of `Pi(M,B(w,w))`. In the latter write `B(w,w)=(w dot grad)w+grad p_w`. Integration by parts of the two local advection contributions gives, on horizontal vectors,

\[
K_{{\rm loc},ij}=-3(P\cdot\nabla)Q_{ij}.
\tag{3}
\]

For clarity, before simplifying this tensor is
`-2 sym grad(QP) - P dot grad Q + 2Q`. Since `partial_i Q_jk=partial_k Q_ij` and `partial_i P_k=delta_ik` for horizontal `i`, the two undifferentiated `Q` terms cancel and (3) follows. This cancellation uses the horizontal seed assumption.

The complete annular pressure term is also explicit. Let

\[
J_\eta=\frac1{4\pi}\sum_{\ell=0,2,4}\psi_\ell(R)P_\ell(z/R),
\qquad -\Delta J_\eta=\operatorname{div}(Q S_\eta).
\]

The globally regular, decaying adjoint from pass 5 retains both pump joins and all Newton tails. On its `eta=-1` plateau,

\[
\psi_0=\text{constant},\quad
\psi_2=\frac{16}{7R^3},\quad
\psi_4=\frac{36}{7R^3}+a_\eta R^4+\frac{b_\eta}{R^5},
\]

where the actual join moments satisfy `a_eta<0` and `0<b_eta<=5`. No freely chosen harmonic field is inserted. Green self-adjointness and two integrations by parts give

\[
2\int Q S_\eta\cdot\nabla p_w
=-2\int J_\eta\,\partial_{ij}(w_iw_j)
=-2\int (\partial_{ij}J_\eta)w_iw_j.
\]

Thus the **complete** pump kernel is

\[
\boxed{K_{{\rm pump},ij}=-3(P\cdot\nabla)Q_{ij}
-2\partial_{ij}J_\eta.}
\tag{4}
\]

The core sees the harmonic pressure of the separated seed. The already derived core response, including its contact contribution, is `18 p_w,zz/7-(M_psi/15)p_w,zzzz`, with `M_psi=5/8`. Since `p_w,zz=-integral Q_ij w_iw_j` and `p_w,zzzz=-integral partial_zz Q_ij w_iw_j`, this supplies

\[
\boxed{K_{{\rm core},ij}=-\frac{18}{7}Q_{ij}
+\frac1{24}\partial_{zz}Q_{ij}.}
\tag{5}
\]

Only the azimuthal mean of `p_w` contributes to this axisymmetric core functional. The nonzero angular pressure modes do not affect the displayed axial derivatives or the radial-core pairing. The core swirl `W_psi` contributes zero to the pressure-gradient pairing by the same angular selection and azimuthal/Hessian cancellation. Therefore

\[
\boxed{\mathcal L(M;w)=\int
(K_{{\rm core},ij}+cK_{{\rm pump},ij})w_iw_j\,dx.}
\tag{6}
\]

## 2. Explicit transverse kernels and conservative signs

Use the orthonormal cylindrical basis, put `t=z^2/R^2`, and define

\[
q(t)=35t^2-35t+4,\quad
p(t)=231t^3-336t^2+119t-6,
\]
\[
h(t)=\frac{21t^2-14t+1}{5t-1},\quad
j(t)=52920t^3-77840t^2+28280t-1504.
\]

The transverse entries are diagonal. Their exact ratios are

\[
\begin{aligned}
\frac{K_{{\rm core},rr}}{Q_{rr}}
&=-\frac{18}{7}+\frac5{8R^2}\frac{p(t)}{q(t)},\\
\frac{K_{{\rm core},\theta\theta}}{Q_{\theta\theta}}
&=-\frac{18}{7}+\frac5{8R^2}h(t),\\
\frac{K_{{\rm pump},rr}}{Q_{rr}}
&=-\frac{j(t)}{28q(t)}
+a_\eta R^7\frac{7t-3}{q(t)}
+\frac{5b_\eta}{4R^2}\frac{p(t)}{q(t)},\\
\frac{K_{{\rm pump},\theta\theta}}{Q_{\theta\theta}}
&=\frac{24}{7}+\left(\frac{5b_\eta}{4R^2}-10\right)h(t)
+a_\eta R^7.
\end{aligned}
\tag{7}
\]

For the actual outer seed, `R>=17/5` and `1-t<=49/4673<1/90`. Put `x=1-t` and `x_+=1/90`. Then

\[
q=4-35x+35x^2\in[4-35x_+,4],
\]
\[
p=8-140x+357x^2-231x^3\le8+357x_+^2,
\]
\[
j=1856-31360x+80920x^2-52920x^3
\ge1856-31360x_+-52920x_+^3>0.
\]

The `a_eta` terms in (7) are negative. Using `b_eta<=5`, the displayed rational bounds give

\[
\frac{K_{{\rm pump},rr}}{Q_{rr}}<-12,
\qquad \frac{K_{{\rm core},rr}}{Q_{rr}}<-\frac{12}{5}.
\]

Also `h` is increasing in this range and `19/10<h(t)<=2`. Its coefficient in the pump expression is negative. The same substitutions give

\[
\frac{K_{{\rm pump},\theta\theta}}{Q_{\theta\theta}}<-12,
\qquad \frac{K_{{\rm core},\theta\theta}}{Q_{\theta\theta}}<-\frac{12}{5}.
\]

Both `Q_rr,Q_thetatheta` are strictly positive on the seed cone. Since `w_z=0`, (6) proves

\[
\boxed{\mathcal L(M;w)\le-
\left(\frac{12}{5}+\frac75\,12\right)\int Q_{ij}w_iw_j
=-\frac{96}{5}C_w.}
\tag{8}
\]

No derivatives of the envelope enter (8). In particular, the pressure source's high angular modes have not been approximated by a pressure truncation; the adjoint pairing evaluates their contribution to this scalar functional exactly.

## 3. Viscosity reduces to a one-dimensional range certificate

The seed and its derivatives avoid the origin, and `Delta Q=0` there. Integration by parts gives

\[
d_w:=2\Pi(w,\Delta w)
=2\int Q_{ij}\partial_kw_i\partial_kw_j\,dx.
\tag{9}
\]

Every Cartesian derivative of `w` remains horizontal, so this is nonnegative on the support. The transverse kernel bound from pass 6 gives

\[
0\le d_w\le\frac6{\pi(17/5)^5}\|\nabla w\|_{L^2}^2.
\tag{10}
\]

Write the seed potential as `e(r)q(z)cos(m theta+kr)`, with `m=4,k=20` and the exact envelopes from pass 6. The integrated horizontal Hessian identity
`integral |D_h^2 potential|^2=integral |Delta_h potential|^2`
and exact angular integration show

\[
\|\nabla w\|_{L^2}^2=\pi(R_1Z_0+R_0Z_1),
\tag{11}
\]

where

\[
\begin{aligned}
R_1&=\int r\left[
\left(e''+\frac{e'}r-(400+16/r^2)e\right)^2
+400(2e'+e/r)^2\right]dr,\\
R_0&=\int r\left[(e')^2+(400+16/r^2)e^2\right]dr,\\
Z_0&=\int q^2dz,\qquad Z_1=\int(q')^2dz.
\end{aligned}
\]

All cylindrical angular-basis derivatives are included in (11). The already certified `|chi'|<4`, compact supports, and total variations give

\[
Z_0\le\frac{12}{5},\quad Z_1\le\frac{160}{3},\quad
R_0\le\frac{6764}{75}.
\tag{12}
\]

Only `R1` requires a new integral. Outward-rounded interval ranges of the unchanged cutoff and its first two derivatives, on 1024 adaptive rational panels at 35 decimal digits, enclose it by

\[
100952.1793303089676983<R_1<113430.747104983141991<120000.
\tag{13}
\]

The exact upper dyadic endpoint is
`1150324437215746995662639719374396981 * 2^(-103)`.
There is no floating-point quadrature-error estimate or numerical pressure solve in this enclosure.

Combining (10)--(13) gives the coarse fully rational certificate

\[
\nu d_w\le\frac6{1000}(\tfrac5{17})^5
\left[\frac{12}{5}\,120000
+\frac{6764}{75}\frac{160}{3}\right]<4.
\tag{14}
\]

Using the stored tighter endpoint gives `nu d_w<3.658516517079472128`. The interval printed by the script around this last upper-bound expression is **not** a two-sided interval for `nu d_w`; only its upper endpoint bounds that coefficient. Equations (8) and (14) prove (1).

The files `verify_D_pressure_kernels.py`, `certify_D_viscosity.py`, and `D-viscosity-certificate.json` contain the exact kernel checks, reproducible interval integral, and exact endpoints. The trust base is the analytical identities and the `mpmath.iv` range engine. There is no assertion here about later pressure evolution, persistence of torque, or a return map.

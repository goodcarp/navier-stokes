# A compact rotational receiving mode with exact positive second derivative

Status: an exact local lemma for the ordinary unforced Navier–Stokes solution constructed in `work/pass2/active-core-pressure.md`, followed by a conditional quantitative remainder estimate. No novelty, long-time regeneration, material-circulation amplification or blowup claim is made.

## Assumptions and observable

Let \(\nu>0\), and let \(u\) be the local smooth finite-energy NS solution on \(\mathbb R^3\) with smooth compactly supported datum from the active-core pressure lemma. In some ball \(B_{r_*}\),
\[
u_0(x)=\Omega Jx,\qquad \Omega>0,\qquad
J=\begin{pmatrix}0&-1&0\\1&0&0\\0&0&0\end{pmatrix}.
\]
The full datum is odd under spatial inversion, so \(u(t,0)=0\) throughout its smooth lifespan. Its decaying initial pressure satisfies
\[
K=-\partial_{33}p_0(0)=A^2 C_v-\frac25\Omega^2>0.
\]
All pressure here is the pressure of the **full** compact datum, including the outer packets and rotating core.

Fix a nonnegative nonzero \(\chi\in C_c^\infty([0,1))\). For \(0<r\le r_*\) with the support strictly inside the affine-core neighborhood, put
\[
\chi_r(x)=\chi(|x|^2/r^2),\qquad
\varphi_r(x)=\chi_r(x)Jx,\qquad
N_r=\int\chi_r|Jx|^2\,dx>0.
\]
Since \(x\cdot Jx=0\) and \(\operatorname{tr}J=0\), \(\nabla\cdot\varphi_r=0\). Define the signed rotational coefficient and normalized amplitude
\[
\Omega_r(t)=\frac{\langle u(t),\varphi_r\rangle}{N_r},
\qquad a_r(t)=\frac{\Omega_r(t)}\Omega.
\]
The integration region and weight are fixed in physical space. They do not move with the fluid.

**Exact conclusion, at every such radius:**
\[
\boxed{a_r(0)=1,\qquad a_r'(0)=0,\qquad a_r''(0)=K.}
\]
Consequently \(a_r(t)=1+Kt^2/2+o(t^2)\) at each radius. Below, a uniform third-derivative bound makes the remainder uniform over the entire radius family.

## Proof, retaining pressure and viscosity

The initial affine field has \(\Delta u_0=0\) in the core, and
\[
(u_0\cdot\nabla)u_0=(-\Omega^2x_1,-\Omega^2x_2,0).
\]
Therefore NS gives exactly
\[
u_t(0,x)=\nabla\psi(x),\qquad
\psi(x)=\frac{\Omega^2}{2}(x_1^2+x_2^2)-p_0(x).
\]
The pressure Poisson equation there is
\[
-\Delta p_0=\sum_{i,j}\partial_i u_{0,j}\,\partial_j u_{0,i}=-2\Omega^2.
\]
It follows that \(\Delta\psi=0\), \(\Delta u_t(0)=0\), and \(\psi_{33}(0)=K\). These identities hold on an open neighborhood of every receiver support, so differentiating and applying the Laplacian is justified. The pressure is not discarded; it is the source of \(K\).

Testing \(u_t(0)=\nabla\psi\) against \(\varphi_r\) gives \(a_r'(0)=0\). In fact, **every** compactly supported divergence-free test field in the initial affine core has zero initial velocity-projection derivative. The first-order change in pointwise strain is a harmonic-gradient effect, consistent with this cancellation.

Differentiate NS in time at zero:
\[
u_{tt}(0)=-(u_t(0)\cdot\nabla)u_0-(u_0\cdot\nabla)u_t(0)
-\nabla p_t(0)+\nu\Delta u_t(0).
\]
The viscosity term is exactly zero in the receiver neighborhood. Substituting the affine field and harmonic gradient, and using \(J^T=-J\), gives
\[
\begin{aligned}
u_{tt}(0)
&=-\Omega J\nabla\psi-\nabla^2\psi\,(\Omega Jx)-\nabla p_t(0)\\
&=-2\Omega J\nabla\psi
-\nabla\big((\Omega Jx)\cdot\nabla\psi+p_t(0)\big).
\end{aligned}
\]
Both pressure and the additional nonlinear gradient vanish after compact divergence-free testing. Hence
\[
\langle u_{tt}(0),\varphi_r\rangle
=-2\Omega\int\chi_r(x)(x_1\psi_1+x_2\psi_2)\,dx.
\]

For completeness, the radial harmonic first-moment identity used here is
\[
\int\chi_r(x)x_i h(x)\,dx=I_r\partial_i h(0),
\qquad I_r=\int\chi_r(x)x_i^2\,dx,
\]
for any harmonic \(h\) on a neighborhood of the receiver ball; radial symmetry makes \(I_r\) independent of \(i\). To prove it, \(\partial_i h\) is harmonic, so its volume mean on \(B_s\), followed by the divergence theorem, yields
\[
s^2\int_{S^2}h(s\theta)\theta_i\,d\theta
=\int_{B_s}\partial_i h\,dx
=\frac{4\pi s^3}{3}\partial_i h(0).
\]
Multiply by the radial weight and integrate in \(s\). The result is the identity above, since \(I_r=(4\pi/3)\int_0^r\chi(s^2/r^2)s^4\,ds\).

Apply this identity to \(h=\psi_i\), which is harmonic. Because \(\Delta\psi=0\) and \(N_r=2I_r\),
\[
\langle u_{tt}(0),\varphi_r\rangle
=-2\Omega I_r(\psi_{11}(0)+\psi_{22}(0))
=2\Omega I_rK=\Omega N_r K.
\]
Dividing by \(\Omega N_r\) proves \(a_r''(0)=K\). The local smooth solution justifies the time derivatives and Taylor expansion. No persistent separation of core and packets has been assumed.

## A uniform remainder sufficient for finite modal gain

Assume on a smooth time interval \([0,T_0]\) that
\[
M_3:=\sup_{0\le s\le T_0}\|\nabla\partial_t^3u(s)\|_{L^\infty(B_{r_*})}<\infty.
\]
Such a finite bound exists on a sufficiently short compact interval for the fixed smooth datum; extracting its dependence on the stage parameters is separate quantitative work. Oddness implies \(\partial_t^3u(s,0)=0\), whence
\[
|\partial_t^3u(s,x)|\le M_3|x|.
\]
As \(|Jx|\le|x|\) and \(N_r=(2/3)\int\chi_r|x|^2\,dx\),
\[
|a_r'''(s)|\le\frac{3M_3}{2\Omega}
\]
uniformly in \(r\). Taylor's integral remainder gives
\[
\boxed{\left|a_r(t)-1-\frac12Kt^2\right|
\le\frac{M_3}{4\Omega}t^3.}
\]
In particular, for \(0<t\le T_0\) and \(t\le\Omega K/M_3\), interpreting the latter restriction as vacuous if \(M_3=0\),
\[
a_r(t)\ge1+\frac14Kt^2>1
\]
at every receiver radius. This is a uniform small finite gain; it is not a lower bound valid for a full strain time unless the parameter dependence of \(M_3\) warrants that longer interval.

## Actual weighted velocity energy and a smaller fixed receiver

Let \(E_r(t)=\frac12\int\chi_r|u(t)|^2\,dx\). Weighted Cauchy–Schwarz gives the exact inequality
\[
E_r(t)\ge\frac12\Omega_r(t)^2N_r=a_r(t)^2E_r(0),
\]
because the initial rigid rotation saturates this inequality. Thus positive modal gain implies actual weighted kinetic-energy gain in that fixed receiver. It also bounds its weighted RMS velocity
\[
U_r(t)=\left(\frac{\int\chi_r|u(t)|^2}{\int\chi_r}\right)^{1/2}
\ge |a_r(t)|U_r(0),\qquad U_r(0)=c_\chi\Omega r,
\]
where \(c_\chi>0\) depends only on the fixed weight.

Fix \(0<q<1\). Comparing the initial radius \(r\) to the final radius \(qr\), and defining the RMS-based receiver Reynolds number \(\mathrm{Re}_r(t)=rU_r(t)/\nu\), yields
\[
\frac{U_{qr}(t)}{U_r(0)}\ge q\,a_{qr}(t),\qquad
\frac{\mathrm{Re}_{qr}(t)}{\mathrm{Re}_r(0)}\ge q^2a_{qr}(t).
\]
Accordingly, a positive lower bound \(a_{qr}(t)\ge G\) gives a strict receiver Reynolds gain if
\[
\boxed{G>q^{-2}.}
\]
The uniform small-time estimate supplies some such contraction: choose any
\[
(1+Kt^2/4)^{-1/2}<q<1.
\]
This proves a small finite gain of a precisely defined velocity/Re observable in a smaller fixed receiver, for the one actual NS solution. It does not prove gain by a prescribed substantial contraction, nor a scale-independent iteratable transfer map.

## Exact scale matching for one receiving time

The following consequence was proposed by the root agent and independently checked here. Fix one admissible time \(\tau>0\) with a uniform modal lower bound \(a_s(\tau)\ge G>1\) for every \(0<s\le r\). Choose any exponent \(0<a<1/2\). Define
\[
F(q)=q^a\frac{\mathrm{Re}_{qr}(\tau)}{\mathrm{Re}_r(0)}
=q^{a+1}\frac{U_{qr}(\tau)}{U_r(0)},\qquad 0<q\le1.
\]
This is continuous. Smooth bounded velocity implies \(F(q)\to0\) as \(q\downarrow0\), whereas \(F(1)\ge G>1\). The intermediate value theorem supplies at least one \(q\in(0,1)\) satisfying \(F(q)=1\). No monotonicity, uniqueness, or continuity of a chosen root in the datum is asserted or needed. For every such root,
\[
\boxed{\frac{\mathrm{Re}_{qr}(\tau)}{\mathrm{Re}_r(0)}=q^{-a},
\qquad\frac{U_{qr}(\tau)}{U_r(0)}=q^{-(1+a)}.}
\]
The weight mass is exactly \(\int\chi_{qr}=q^3\int\chi_r\). Therefore the weighted kinetic energies also obey the exact relation
\[
\boxed{\frac{E_{qr}(\tau)}{E_r(0)}=q^{1-2a}<1.}
\]
These are equalities between actual observables of one ordinary NS solution, comparing the initial radius \(r\) and the later fixed receiver radius \(qr\). The exponent \(a\) and output radius are selected for this episode; this calculation does not manufacture a new solution or alter its viscosity.

There are also explicit bounds on the selected ratio. The modal lower bound gives
\[
F(q)\ge Gq^{a+2},\qquad q\le G^{-1/(a+2)}<1
\]
at any root. Let
\[
M_1=\sup_{x\in B_r}\|\nabla u(\tau,x)\|_{\mathrm{op}}.
\]
Oddness yields \(|u(\tau,x)|\le M_1|x|\). The radial identity
\(\int\chi_s|Jx|^2=(2/3)\int\chi_s|x|^2\) then implies
\[
\frac{U_{qr}(\tau)}{U_r(0)}
\le\sqrt{\frac32}\frac{M_1}\Omega q.
\]
Thus every root lies in the compact interval
\[
\boxed{\left(\frac\Omega{\sqrt{3/2}\,M_1}\right)^{1/(a+2)}
\le q\le G^{-1/(a+2)}.}
\]
The existence argument guarantees consistency of these bounds. Quantitative parameter dependence of \(G\) and \(M_1\) is needed for a uniform family of stages; the fixed profile class below supplies it.

### A uniform contraction bracket on the normalized profile class

Section 6 of `work/pass3/quantitative-core-stage.md` supplies explicit \(\tau_*,g_*>0\) and \(\bar V<\infty\) for a fixed normalized parameter box with a positive pressure margin. In that notation \(\Omega=S\theta\), \(S>0\), \(\theta\ge\theta_->0\), and at the physical time \(\tau_*/S\),
\[
a_s(\tau_*/S)\ge1+g_*\quad(0<s\le r),\qquad
\|\nabla u(\tau_*/S)\|_{L^\infty,\mathrm{op}}\le S\bar V.
\]
Substitution into the radius bounds gives, for every member of that class and every selected scale-matching root,
\[
\boxed{0<q_-:=\left(\frac{\theta_-}{\sqrt{3/2}\,\bar V}\right)^{1/(a+2)}
\le q\le
q_+:=(1+g_*)^{-1/(a+2)}<1.}
\]
Thus one-time scale matching has a uniform positive contraction bracket for that fixed profile class, independently of the allowed initial receiver radius and amplitude scale \(S\). This uses the quantitative theorem's actual solution and norm bound, not a frozen pressure field. The bracket does not imply that the evolved velocity or outer packets lie in the same profile class, so it establishes no profile closure or return map.

## Material and regeneration limits

The receiver is **Eulerian and nonmaterial**. It can receive kinetic energy through its boundary and can sample a changing portion of the initially rotating material reservoir. This conclusion does not contradict the material-circulation cancellations in `work/pass3/material-core-transfer.md`. It neither establishes enhanced circulation of a transported loop nor shows that the final velocity retains a rigid rotational profile. The exact scale matching selects a radius after observing a single smooth episode; it is not a return map.

A repeatable stage additionally requires controlled terminal core shape, favorable outer packets, inherited history and total energy/force budgets. The reservoir or changing material footprint must be accounted for on restart. Smooth forcing remains available in broader C/D research, although this specific lemma uses an unforced solution.

## Finite algebra verification

`verify_rotational_receiver.py` passed 15 exact SymPy checks, including the full generic harmonic polynomial spaces of degrees two, four and six (dimensions five, nine and thirteen), their radial first moments, the differentiated nonlinear identity, the vanishing local viscous term, the weighted projection normalization and the scale exponents. The higher-degree tests check that nonconstant pressure Hessians contribute no hidden radial moment to the asserted second derivative. The script does not certify the pressure kernel calculation, local PDE existence, Sobolev estimates, or regeneration; those require the stated analytic arguments.

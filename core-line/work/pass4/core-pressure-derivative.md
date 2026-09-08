# Exact radial reduction of the core-only pressure derivative

For the actual profiles `S,W` in `affine-core-pressure.md`, set `c=A=0` in Section 4 of `full-pressure-feedback-gate.md`. The full whole-space initial derivative has the exact form

\[
 p_{zz}'(0)=t_{300}b^3+t_{Wb}b\Omega^2,
 \qquad \boxed{d_S=d_W=0}.
\]

The two cubic coefficients depend on the radial cutoff shape. They can both be evaluated by one-dimensional radial integrals, with the nonlocal pressure retained exactly. No time evolution, closure of the radial profile family, or sign certification for the full gate is asserted here.

## 1. Exact coefficient formulas

All primes below differentiate `s=|x|²`. Let

\[
 p(s)=\psi(s),\qquad f(s)=p(s)+\frac23sp'(s),
\]
\[
 I=\int_0^\infty sp(p')^2\,ds,\qquad
 J=\int_0^\infty s^2(p')^3\,ds,
\]
\[
 Q(s)=\int_0^s p(\tau)\,d\tau,\qquad
 K(s)=\int_0^s\tau^{7/2}
       \left[p(\tau)p'(\tau)+\frac79\tau(p'(\tau))^2\right]d\tau.
\]

Then

\[
 \boxed{t_{Wb}=\frac{44}{35}-\frac{16}{5}\int_0^\infty p'f^2\,ds
 =\frac{244}{105}-\frac{64}{15}I-\frac{64}{45}J,}
 \tag{1}
\]

and

\[
\boxed{
\begin{aligned}
 t_{300}={}&\frac{340}{49}-\frac{18368}{735}I-\frac{448}{105}J
       -\frac{256}{105}\int_0^\infty Q(s)(p'(s))^2\,ds\\
       &+\frac{1536}{49}\int_0^\infty s^{-7/2}p'(s)K(s)\,ds.
\end{aligned}}
\tag{2}
\]

These formulas require only `p,p'` and cumulative radial integrals. All integrands vanish near zero except harmless compact-profile terms, because `p=1` there. No principal-value quadrature remains. Radially dilating the same cutoff shape leaves these dimensionless cubic coefficients unchanged; changing its shape need not.

For a nonincreasing cutoff, (1) also gives the exact lower bound `tWb >= 44/35`, since `p' f² <= 0`. This is only a statement about the mixed core coefficient.

## 2. The full pressure that enters this reduction

Write the solid axisymmetric harmonics

\[
 H_2(x)=\frac{3z^2-s}{2},\qquad
 H_4(x)=\frac{35z^4-30sz^2+3s^2}{8}.
\]

The decaying pressures of the unit profiles have exactly the finite expansions

\[
 p[S]=h_0(s)+h_2(s)H_2(x)+h_4(s)H_4(x),\qquad
 p[W]=a_0(s)+a_2(s)H_2(x).
\]

Their radial functions are

\[
 h_0'=-p^2-\frac45spp'+\frac4{15}s^2(p')^2,
\]
\[
 h_2=-\frac27p^2+\frac4{21}s^{-5/2}
                  \int_0^s\tau^{7/2}(p'(\tau))^2d\tau,
\]
\[
\begin{aligned}
 h_4={}&\frac{48}{35}s^{-9/2}\int_0^s\tau^{7/2}p(\tau)p'(\tau)d\tau\\
       &+\frac{16}{15}s^{-9/2}\int_0^s\tau^{9/2}(p'(\tau))^2d\tau
        +\frac{32}{21}\int_s^\infty(p'(\tau))^2d\tau,
\end{aligned}
\]
\[
 a_0'=\frac13f^2,\qquad
 a_2=-\frac13s^{-5/2}\int_0^s\tau^{3/2}f(\tau)^2d\tau.
\tag{3}
\]

The constants in `h0,a0` are fixed by decay at infinity; their derivatives suffice for this calculation. The pressure tails in (3) are retained even though the initial velocities are compactly supported. In the affine neighborhood, `h4=(32/21)∫(p')²`, a profile-dependent harmonic pressure term invisible to the initial central Hessian but relevant to its next derivative.

To check (3), expand the source of `S` in Legendre polynomials:

\[
 g_S=G_0(s)+sG_2(s)P_2(z/\sqrt{s})+s^2G_4(s)P_4(z/\sqrt{s}),
\]
\[
\begin{aligned}
 G_0={}&6p^2+16spp'+\frac{16}{5}s^2pp''
              -\frac8{15}s^2(p')^2-\frac{32}{15}s^3p'p'',\\
 G_2={}&\frac8{21}(21pp'+6spp''+2s(p')^2-4s^2p'p''),\\
 G_4={}&-\frac{64}{35}(3pp''-13(p')^2-2sp'p'').
\end{aligned}
\]

For `W`, the corresponding coefficients are
`G0W=-2f²-(8/3)sff'` and `G2W=(8/3)ff'`.
For every degree `l`, the exact scalar Poisson equation is

\[
 4s h_l''+(4l+6)h_l'=-G_l.
\]

The displayed profiles satisfy this ODE, regularity at zero, and the decaying whole-space boundary condition. This supplies the nonlocal Leray contribution without projecting the evolved velocity onto the `S,W` profile span.

## 3. Cubic source and contact term

In cylindrical coordinates with transverse radius `r` and axial coordinate `z`,

\[
 S_r=-r(p+2z^2p'),\quad S_z=2z(p+r^2p'),\quad W_\theta=rf.
\]

Define the full inviscid initial accelerations by their amplitude coefficients:

\[
 V_{SS}=-(S\cdot\nabla)S-\nabla p[S],\quad
 V_{WW}=rf^2e_r-\nabla p[W],
\]
\[
 V_{SW}=-[(S\cdot\nabla)W+(W\cdot\nabla)S].
\]

The last vector is azimuthal and has no pressure term. The exact pressure-derivative sources are

\[
 g_{300}=2\operatorname{tr}(\nabla S\nabla V_{SS}),\qquad
 g_{Wb}=2\operatorname{tr}(\nabla S\nabla V_{WW})
          +2\operatorname{tr}(\nabla W\nabla V_{SW}).
\]

For this core-only datum the actual central derivatives are
`b'=−5b²/7−Ω²/5`, `Ω'=2bΩ`. Consequently

\[
 g_{300}(0)=-\frac{60}{7},\qquad g_{Wb}(0)=-\frac{52}{5}.
\]

The distributional Hessian contact contributions are therefore `20/7` and `52/15`. They are included in (1)–(2). In particular, these are different from the contact values at the full externally tuned neutral boundary.

The checker forms these full cylindrical sources, integrates their angular polynomial exactly using `⟨μ^(2k)⟩=1/(2k+1)`, and substitutes (3). For the mixed term the radial pressure integrals cancel by integration by parts and yield (1). For the strain cubic, the four remaining pressure moments reduce by integration by parts and Fubini to `Q,K`, yielding (2). Boundary terms use the flat central value `p(0)=1` and the compact outer support.

## 4. Why the two viscosity coefficients vanish exactly

Let `T_S[p]` denote the radial strain extension. Its vector potential is a radial function times a harmonic quadratic vector polynomial, so

\[
 \Delta T_S[p]=T_S[4sp''+14p'].
\]

The second profile vanishes near zero. The already established quadratic identity

\[
 \Pi(T_S[p],T_S[p])=-\frac{18}{7}p(0)^2
\]

therefore has zero derivative in that direction. This proves
`dS=2Π(S,ΔS)=0`.

Likewise `W=f(s)Jx` and

\[
 \Delta W=(4sf''+10f')Jx.
\]

The universal radial-swirl pressure identity is `Π(fJx,fJx)=2f(0)²/5`. Since the new amplitude vanishes near zero, its directional derivative gives `dW=0`.

These are exact initial scalar cancellations. The viscous velocity is generally nonzero in the cutoff transition, and no claim is made that the later curvature terms or viscosity contribution remain zero. Cross strain/swirl pressure terms vanish by axisymmetry and poloidal/azimuthal parity.

## 5. Independent numerical control for the root's flat cutoff

For the cutoff used in `evaluate_pressure_gate.py`, namely `p=1` for `s<=1/4`, `p=0` for `s>=1`, and

\[
 p(s)=\operatorname{logistic}\!\left(\frac1a-\frac1{1-a}\right),
 \qquad a=(s-1/4)/(3/4),
\]

independent one-dimensional floating-point quadrature gives

| Quantity | Value |
|---|---:|
| `t300` | `−6.7344633112` |
| `tWb` | `2.563122131694` |
| `dS,dW` | exactly `0,0` |
| core-only `pzz'`, `b=Ω=1` | `−4.1713411795` |

The radial grids with 2001, 4001, 8001 and 16001 points agree at this precision. The derivative-reduced formula gives `t300≈−6.73446331117`; the original full-source radial formula differs by less than `5e−11`, consistent with floating-point cancellation in its higher derivatives. This is a convergence check, not an interval bound. The script `evaluate_core_pressure_radial.py` uses the exact angular reduction and radial Green-function moments; it does not call the root's whole-space pressure solver.

As exact algebra checks on cutoff dependence, take polynomial proxies `p=(1−s)^n` on `[0,1]`, zero outside. These are **not** the admissible flat smooth cutoffs, and no viscosity conclusion is transferred from them. For their inviscid cubic functionals, all radial integrals are rational and the checker finds:

| Proxy | `t300` | `tWb` |
|---|---:|---:|
| `n=4` | `51414292252/11712375675` | `9188/4725` |
| `n=5` | `16529009732/3720401685` | `23908/12285` |

The formulas (1)–(2) agree with direct exact angular/source integration for these proxies. Their different values also rule out universal cubic coefficients for the flat smooth class: the proxies can be approximated in compactly supported `C^(2,α)`, `0<α<1`, by smooth radial cutoffs flat near zero, and these inviscid cubic functionals are continuous in that topology. The pressure Hessian continuity follows from the Hölder singular-integral estimate; the acceleration source uses at most two velocity derivatives.

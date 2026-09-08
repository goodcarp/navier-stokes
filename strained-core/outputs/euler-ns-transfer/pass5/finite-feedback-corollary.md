# A finite-time consequence of the certified initial feedback jets

For the fixed compact datum in the pass5 certificate, the favorable initial jets imply an actual positive interval of simultaneous strain-ratio improvement, rotation and receiving-velocity gain, and growth of the specified swirl-pressure ratio. This follows for the unforced smooth Navier–Stokes solution itself. The interval depends on this datum and its derivative bounds; no radial-profile preservation or repeated-stage conclusion is used.

## 1. Certified inputs and exact later-time observables

Take `b(0)=Omega(0)=1`, `c=7/5`, `nu=1/1000`, and the exactly tuned amplitude

\[
 A^2C_v=H_0:=\frac{204}{35}.
\]

Let `u` be the unique local smooth solution on `R³` from the specified smooth compact divergence-free datum. Axisymmetry and oddness are preserved by uniqueness, so the origin is a material trajectory and

\[
 \nabla u(t,0)=b(t)D+\Omega(t)J,
 \quad D=\operatorname{diag}(-1,-1,2),
 \quad Jx=(-y,x,0).
\]

Define the actual central quantities

\[
 \beta=b/\Omega,\qquad S=b'-2b^2,
 \qquad W=\frac12\Delta\omega_z(t,0).
\]

These do not assume that the velocity remains affine away from the origin. The certified gate and exact insertion identities give

\[
 \begin{gathered}
 T_{.001}\le-\frac{290649}{8750},\qquad
 \delta:=\frac{290649}{17500}>0,\\
 b'=\Omega'=2,\quad \beta'=0,\quad
 \beta''=-T_{.001}/2\ge\delta
 \qquad\text{at }t=0.
 \end{gathered}
 \tag{1}
\]

The exact vorticity equation at the center is

\[
 \Omega'=2b\Omega+\nu W.
\]

Initial affine-core identities give `W(0)=W'(0)=0`. Consequently

\[
 \Omega''(0)=8,\qquad
 b''(0)=8+\beta''(0),\qquad
 \boxed{S(0)=0,\quad S'(0)=\beta''(0)\ge\delta.}
 \tag{2}
\]

This is why the two favorable inequalities below both follow. At positive viscosity, `S>0` and `beta'>0` are not simply equivalent at later times; the curvature term has not been set to zero there.

Fix a nonnegative axisymmetric cutoff `zeta`, zero near the core and one on a neighborhood of the original swirl support. It may be supported in the exterior region where `Qtheta>0`. Define the actual fixed-region swirl stress and its normalized ratio by

\[
 C(t)=\int\zeta Q_{\theta\theta}u_\theta(t)^2\,dx,
 \qquad \mathcal R(t)=C(t)/\Omega(t)^2.
\]

Then `C(0)=R(0)=H0`. The certified diffusion quotient and the proved plateau multiplier imply

\[
 \left.\frac{\mathcal R'}{\mathcal R}\right|_{t=0}
 \ge\frac{2381}{357}\frac75-4-\frac{901}{1000}
 =4+\varepsilon,
 \qquad\varepsilon:=\frac{22249}{51000}>0.
 \tag{3}
\]

This concerns the defined swirl channel. It does not identify `C(t)` with the full exterior pressure contribution.

## 2. An unconditional local corollary for this fixed datum

By smoothness and (1)–(3), there exists `tau>0` such that, for `0<t<=tau`,

\[
 \boxed{\beta'(t)\ge\frac\delta2t>0,
 \qquad\beta(t)\ge1+\frac\delta4t^2,}
 \tag{4}
\]

\[
 \boxed{b'(t)-2b(t)^2\ge\frac\delta2t>0,
 \qquad\Omega'(t)\ge1,\quad\Omega(t)\ge1+t,}
 \tag{5}
\]

and

\[
 \boxed{\frac{\mathcal R'(t)}{\mathcal R(t)}\ge4,
 \qquad\mathcal R(t)\ge H_0e^{4t}.}
 \tag{6}
\]

In particular, the neutral initial strain-ratio boundary is followed by an actual interval of favorable feedback. These are inequalities along one ordinary unforced solution, with its full pressure and viscosity. Nothing in this argument guarantees that the interval reaches a full strain time or any prescribed positive duration.

The next section makes the dependence of such an interval explicit using the already proved pass3 derivative bounds.

## 3. A quantitative interval in terms of the new datum's Sobolev norm

Let

\[
 X_0=\|u_0\|_{H^{10}_{\rm der}},\qquad
 \|v\|_{H^s_{\rm der}}^2=\sum_{|\alpha|\le s}\|\partial^\alpha v\|_2^2.
\]

The energy/product estimates in [the pass3 derivative estimates](../pass3/quantitative-core-stage.md) apply to arbitrary smooth divergence-free data, not only its earlier rotating-core example. Reuse their explicit constants:

\[
\begin{aligned}
 V&=2X_0,&T_E&=(4\cdot65536X_0)^{-1},\\
 M_1&=3\nu V+16384V^2,\\
 M_2&=(3\nu+4096V)M_1,\\
 M_3&=(3\nu+1024V)M_2+1024M_1^2.
\end{aligned}
 \tag{7}
\]

On `[0,T_E]` these give

\[
 \|u\|_{10}\le V,\quad
 \|u_t\|_8\le M_1,\quad
 \|u_{tt}\|_6\le M_2,\quad
 \|u_{ttt}\|_4\le M_3.
\]

In particular `||grad ∂t^j u||infinity<=Mj` for `j=1,2,3`, with the conservative embedding constants already proved there. Hence `|b^(j)|,|Omega^(j)|<=Mj`.

Restrict first to `t<=1/(2M1)`. Then `|b|<=3/2` and `1/2<=Omega<=3/2`. Differentiating the quotient `beta=b/Omega` gives the valid uniform bound

\[
 |\beta'''|\le B_\beta:=8M_3+96M_1M_2+192M_1^3.
\]

Also

\[
 |S''|=|b'''-4(b')^2-4bb''|
 \le B_S:=M_3+4M_1^2+6M_2.
 \tag{8}
\]

For the reservoir, let

\[
 K_\zeta=\operatorname*{ess\,sup}|\zeta Q_{\theta\theta}|<\infty,
 \quad N_1=2K_\zeta VM_1,
 \quad N_2=2K_\zeta(M_1^2+VM_2).
\]

The azimuthal projection is bounded on `L²`, so direct differentiation of the fixed quadratic functional gives `|C'|<=N1`, `|C''|<=N2`. No spatial derivative of that projection is used. On `t<=H0/(2N1)`, `C>=H0/2`; therefore

\[
 \left|\frac{d}{dt}\left(\frac{\mathcal R'}{\mathcal R}\right)\right|
 \le B_{\mathcal R}:=
 \frac{2N_2}{H_0}+\frac{4N_1^2}{H_0^2}+4M_2+8M_1^2.
 \tag{9}
\]

One completely explicit sufficient interval is

\[
\boxed{
\begin{aligned}
 \tau=\min\bigg\{&T_E,\frac1{2M_1},
      \frac\delta{2B_\beta},\frac\delta{2B_S},
      \frac2{3M_2},\frac8{3(M_1^2+VM_2)},\\
      &\frac{H_0}{2N_1},\frac\varepsilon{2B_{\mathcal R}}\bigg\}>0.
\end{aligned}}
 \tag{10}
\]

A restriction with zero denominator is interpreted as absent. For the nonzero present datum all displayed ingredients are finite. Taylor's integral remainder gives (4)–(5). Equations (3), (9), (10) actually give the stronger logarithmic bound `R'/R>=4+epsilon/2`, hence (6).

**What is and is not quantitative here.** Formula (10) is a numerical-bound recipe using the fixed initial profile integrals and a simple exterior weight supremum. The initial sign certificate by itself does not enclose `X0`, and no new enclosure of this datum's `H10` norm has been performed in this pass. Consequently no decimal duration or useful-size gain is certified here. One can use an analytic upper bound for the cutoff derivatives and support volume instead of a large Sobolev quadrature, but a valid upper bound is still required to evaluate (10). The qualitative positive interval in Section 2 needs no such additional evaluation.

## 4. Uniform-radius receiving velocity and RMS gain

Choose a fixed nonzero smooth nonnegative radial weight `chi`, and let `chi_r(x)=chi(|x|²/r²)` have support in the original affine ball. Put

\[
 I_r=\int\chi_rx_1^2=\int\chi_rx_2^2=\int\chi_rx_3^2,
\]

\[
 \Omega_r(t)=\frac{\int\chi_ru\cdot Jx}{2I_r},
 \qquad b_r(t)=\frac{\int\chi_ru\cdot Dx}{6I_r}.
\]

The exact initial receiving jets are

\[
 \boxed{\Omega_r(0)=b_r(0)=1,
 \qquad\Omega_r'(0)=b_r'(0)=2,}
 \tag{11}
\]

at every allowed radius. Indeed `Delta u_t(0)=0` throughout the original affine ball, and the radial harmonic first-moment identity projects `u_t(0)` onto its actual central derivative `2(D+J)`. It does not remove the higher harmonic pressure terms from the velocity itself.

Oddness and the bound on `grad u_tt` give

\[
 |\Omega_r''|\le\frac32M_2,\qquad |b_r''|\le M_2,
\]

uniformly in `r`. Thus (10) implies

\[
 \boxed{\Omega_r'(t),b_r'(t)\ge1,
 \qquad\Omega_r(t),b_r(t)\ge1+t.}
 \tag{12}
\]

For the actual weighted energy and RMS velocity,

\[
 E_r(t)=\frac12\int\chi_r|u|^2,\qquad
 U_r(t)=\left(\frac{\int\chi_r|u|^2}{\int\chi_r}\right)^{1/2},
\]

the two initial modes are orthogonal, so

\[
 E_r(0)=4I_r,
 \qquad E_r(t)\ge I_r[\Omega_r(t)^2+3b_r(t)^2].
\]

Consequently

\[
 \boxed{E_r(t)\ge(1+t)^2E_r(0),
 \qquad U_r(t)\ge(1+t)U_r(0).}
 \tag{13}
\]

There is also uniform strict monotonicity over this interval. The exact initial energy derivative is `E_r'(0)=16I_r`, and

\[
 |E_r''(t)|\le3I_r(M_1^2+VM_2).
\]

The corresponding restriction in (10) gives `E_r'(t)>=8I_r>0`, and hence `U_r'(t)>0`. This estimate includes the energy in all components of the evolved velocity, not only its rotational projection.

At any chosen `0<t<=tau`, set

\[
 q^2=\frac{1+t/2}{1+t}<1.
\]

Since the initial affine RMS scales linearly with receiver radius, the same actual solution satisfies

\[
 \boxed{\frac{\operatorname{Re}_{qr}(t)}
                   {\operatorname{Re}_{r}(0)}
       =\frac{qr\,U_{qr}(t)}{r\,U_r(0)}
       \ge1+\frac t2>1,}
 \tag{14}
\]

where `Re_r=r U_r/nu`. The same bound holds for the rotational receiving observable `r² Omega_r/nu`.

These are comparisons at smaller **fixed Eulerian** receivers. They do not establish growing circulation of a material loop, confinement of the new vorticity to the smaller receiver, or a return to the starting exterior geometry.

## 5. Scope

The certificate therefore supports a real finite episode with favorable central feedback, increasing normalized swirl stress, and increasing local velocity observables, at positive viscosity and zero external force. Its persistence follows from strict derivative margins and smooth local theory. No inference about a uniform family of inherited states, a repeatable transfer map, indefinite shrinkage, or singularity formation follows from this corollary. In particular, the initial curvature cancellations are used only at the justified endpoint; subsequent curvature is controlled by the actual-solution remainder bounds.

# Independent audit of the outer E interval certificate and common margin

**Accepted for the saved 2048-panel certificate.** The range arithmetic in `certify_outer_E.py` validly encloses the expressions used to prove

\[
\boxed{E<\frac1{60}.}
\]

The exact stored upper endpoint is

\[
353386231037742657024752032401770567\,2^{-124}
=0.01661614072972984853428\ldots<\frac1{60}.
\]

This audit checked the implementation's domains, cutoff calls, Green moment inclusions, radial tails, norm factors, kernel bounds, and final inequality directions. The underlying horizontal Green identity was independently checked by the other agent's exact symbolic calculation in `verify_E_interaction_reduction.py`; it was not replaced by a numerical pressure approximation here. No repeated 2048-panel integration was needed for this audit. The stored dyadic comparison and common-interval arithmetic were checked directly as exact fractions.

## 1. Radial Green moments and both endpoint tails

The seed radial envelope is zero for `r<=7/100` and `r>=7/20`. The carrier cutoff is `C(r)=chi((r/a)^2)`, and `C'` in the script correctly includes `2r/a^2`. The two real components are those of `b(r)=e(r)exp(ikr)`. The trial pressure coefficient is

\[
P=-2Cb-(m-1)r^{-m}I_+(r)-(m+1)r^mI_-(r),
\]
\[
I_+(r)=\int_0^r t^mC'(t)b(t)dt,\qquad
I_-(r)=\int_r^\infty t^{-m}C'(t)b(t)dt.
\]

On each rational panel, the prefix contains every completed lower panel and the suffix contains every completed upper panel. Adding `[0,h]` times the respective integrand range encloses every partial-cell integral. The two partial lengths need not be treated as independent true variables: interval independence only enlarges the enclosure. Both left and right formulas are in the correct direction, and the full endpoint uncertainty is retained. There is no midpoint reset or sampled cumulative value.

For `r<r_min=7/100`, the exact form is `P=c_in r^m`, with `c_in=-(m+1)I_-(0)`. For `r>r_max=7/20`, it is `P=c_out r^(-m)`, with `c_out=-(m-1)I_+(infinity)`. Hence the two contributions to `integral r^3 |P|^2 dr` are exactly

\[
\frac{|c_{\rm in}|^2r_{\min}^{2m+4}}{2m+4},
\qquad
\frac{|c_{\rm out}|^2r_{\max}^{4-2m}}{2m-4}.
\]

These agree with the code, and the exterior integral converges for `m=4`. The trial and its first derivatives match the interior Green solution at the flat cutoff endpoints; no artificial boundary distribution is introduced. The residual norm therefore includes the entire radial tail rather than only the seed support.

## 2. Axial ranges and full-space residual estimate

The axial interval `[3.4,4.6]` contains the complete upper seed bump and the narrower unit-swirl bump. The script multiplies each integral by two for the disjoint reflected packet. All relevant squared or absolute-value expressions are even under reflection. Powers of `z` in this positive-axis integration represent powers of `|z|` in the full-space bounds.

For `F(z)=q(z)v_zprofile(z)`, the computed derivative
`F''=q'' v+2q' v'+q v''`
has the correct chain factors for widths `3/5` and `9/20`. Both cutoff tails are evaluated with the existing outward cutoff-jet engine, including its flat endpoints. No zero extension truncates a nonzero axial residual.

The exact mixed source is horizontally separated as `S_r(r)F(z)`. If `-Delta_h P=S_r`, then the full pressure trial `p_app=P F` has residual `P F''`, with an immaterial sign for its norm. A real nonzero angular mode contributes the factor `pi`, so

\[
\|r\,\mathrm{residual}\|_2^2
=\pi\left(\int_0^\infty r^3|P|^2dr\right)
\left(\int_{\mathbb R}|F''|^2dz\right).
\]

This is precisely the `ne` formula, with both radial tails included. For the whole-space pressure error `e_p`, the angular-mode energy inequality gives

\[
\|e_p/r\|_2\le\frac1{m^2}\|r\,\mathrm{residual}\|_2.
\]

Pairing with `d=div(Qw)` consequently costs at most
`2 ||r d||_2 ||r residual||_2/m^2`, matching the code's `error`. The factor is `2/m^2`, with no missing factor of two from the real mode or the reflected packet. The trial is regular at the axis and has finite energy; the estimate imposes no artificial finite-cylinder boundary condition.

## 3. Bounds for the pairing and the favorable local term

On the support, the exact kernel expressions imply

\[
|Q_{rr}-Q_{\theta\theta}|\le\frac{45}{2\pi}\frac{r^2}{|z|^7},
\qquad
|Q_{rz}|\le\frac{15}{\pi}\frac r{|z|^6}.
\]

Writing `H=(e_r-e/r)^2+k^2 e^2`, the triangle inequality yields

\[
\|rd\|_2\le m\sqrt\pi\left[
\frac{45}{2\pi}\sqrt{\left(\int r^5H\right)\left(\int q^2/|z|^{14}\right)}
+\frac{15}{\pi}\sqrt{\left(\int r^3e^2\right)\left(\int q_z^2/|z|^{12}\right)}\right].
\]

These are exactly `d_A,d_B,Z1,Z2` and `nd`. The bound for the trial pairing is

\[
2\left|\int d\,p_{\rm app}\right|
\le m[45\,\mathrm{trial_A}\,Y_1+30\,\mathrm{trial_B}\,Y_2],
\]

also matching the code. The `2*pi` angular-pairing factor is already included in the constants `45,30`.

The exact local negative contribution has the factor
`(6-t)/(1+t)^(9/2)`, where `t=r^2/z^2`. The implemented

\[
K=\frac{6-t_{\max}}{(1+t_{\max})^5},\qquad t_{\max}=(7/68)^2,
\]

is a valid **lower** bound: replacing exponent `9/2` by `5` decreases the expression for `t>=0`, and `(6-t)/(1+t)^5` decreases on this support. Therefore `local_gain` lower-bounds the magnitude of the favorable local term. Subtracting its interval uses its lower endpoint in the resulting upper endpoint, as required. The shared radial/axial quantities may be correlated, but the interval operations remain inclusion preserving.

`nd`, the absolute trial-pairing bound, and the resulting `E_upper_expression` are auxiliary bound expressions. Their displayed lower endpoints are not lower bounds on the actual coefficient `E`. In particular, the positive lower endpoint printed for the last expression does **not** prove `E>0`.

## 4. Interval arithmetic and assertion scope

The partition endpoints and cell widths are exact fractions. `mesh` uses rational arithmetic for its subdivision counts and introduces no gaps. Cutoff, trigonometric, square-root, prefix/suffix, and integral calculations remain in `mpmath.iv`. Squaring and summing the real/imaginary pressure ranges gives a nonnegative enclosing modulus. As with the earlier certificates, the trusted arithmetic library and audited analytic cutoff bounds form the trust base; this is not a formal proof-kernel replay.

At the time of the source review, the script's `PASS` test asserted the older compatibility threshold rather than the stronger explicit `E<1/60` statement. The exact stored endpoint above independently establishes the stronger claim. Adding a direct `E_upper_expression.upper<1/60` check makes that claim part of the executable certificate.

The source's expression `4+1020*E_upper_expression` is valid for the present saved result, whose upper bound is positive. For reusable runs with a negative certified upper bound, the amplitude range should instead be applied as `max(780*E_upper,1020*E_upper)`, or positivity of the selected upper bound should be asserted. This is a robustness point for other runs; it does not invalidate the current positive upper-bound calculation.

## 5. One uniform bound for the common amplitude interval

For the same exactly retuned family, the established inequalities are

\[
T_\lambda<-m_0+\lambda^2\left[\frac{23199}{1000}C_w+D_\nu+A_\lambda E\right],
\quad m_0=\frac{290649}{8750},
\]
\[
D_\nu<4-\frac{96}{5}C_w,\quad E<\frac1{60},\quad
0<A_\lambda<1020,\quad 0<C_w<\frac32.
\]

Multiplication by the positive exact amplitude gives `A_lambda E<A_lambda/60<17`, regardless of the actual sign of `E`. Hence

\[
\boxed{T_\lambda<-m_0+\lambda^2
\left[\frac{3999}{1000}C_w+21\right].}
\]

The bracket is positive, so `lambda^2<=1` and `C_w<3/2` imply, uniformly on `3/4<=lambda<=1`,

\[
\boxed{T_\lambda<-
\frac{290649}{8750}+\frac{3999}{1000}\frac32+21
=-\frac{435297}{70000}
=-6.218528571428571\ldots<-6.2185.}
\]

This joins the earlier initial global mean-maximum gain and the favorable central pressure-feedback sign on the same exact data. The lower amplitude restriction comes from the mean-maximum estimate; this pressure upper-bound step uses only `lambda<=1`. No conclusion about a uniform normalized duration, inherited-state return, or singularity follows from this common initial margin.

Finalization note: both executable hardening points were implemented. The E integrator now requires E<1/60 for PASS and uses the correct maximum over amplitude endpoints if a future upper bound is negative. The final packaged run rechecks these predicates.

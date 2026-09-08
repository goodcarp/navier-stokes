# Validated local compact integrals with the actual smooth cutoff

Status: range-enclosure quadrature implemented using `mpmath.iv` elementary interval operations, with a proof of the endpoint cutoff bounds. This retains the actual smooth datum. It supplies concrete enclosures for the original outer swirl coefficients and a positive local swirl-channel inequality. It does not certify the total nonlocal pressure derivative or a return map.

## Actual cutoff and endpoint derivative bounds

The cutoff is exactly
\[
f(x)=\begin{cases}1,&x\le1/4,\\
\displaystyle\frac1{1+\exp(-1/t+1/(1-t))},&1/4<x<1,\quad t=(4x-1)/3,\\
0,&x\ge1.
\end{cases}
\]
No polynomial replacement, finite smoothing layer, or exponential clipping changes it. Let \(p(t)=f((1+3t)/4)\) on the transition and \(z(t)=1/t-1/(1-t)\). Then
\[
p'=p(1-p)z',\qquad p''=p(1-p)((1-2p)(z')^2+z''),
\]
\[
z'=-t^{-2}-(1-t)^{-2},\qquad z''=2t^{-3}-2(1-t)^{-3}.
\]
Interior intervals use these identities and outward-rounded arithmetic. Monotonicity gives the tight value range \(p([a,b])\subseteq[p(b),p(a)]\); \(0\le p(1-p)\le1/4\) tightens derivative evaluation.

At the flat endpoints a direct interval substitution would encounter infinite reciprocal powers. The implementation instead covers \(0\le t\le\delta\), \(\delta=1/32\), by analytic bounds. For any \(0<T\le\delta\), define
\[
E_T=\exp(-1/T+1/(1-T)),\quad
M_1(T)=\frac43 E_T(T^{-2}+(1-T)^{-2}),
\]
\[
M_2(T)=\frac{16}{9}E_T\big[(T^{-2}+(1-T)^{-2})^2+2T^{-3}+2(1-T)^{-3}\big].
\]
Throughout \(0\le t\le T\), derivatives in the original \(x\) argument obey
\[
-M_1(T)\le f'(x)\le0,\qquad |f''(x)|\le M_2(T).
\]
To justify them, \(p(1-p)\le\exp(-1/t+1/(1-t))\), while \(t^{-m}e^{-1/t}\) is increasing for \(t<1/m\). Only \(m\le4\) is needed, so \(T\le1/32\) suffices. The remaining factors \((1-t)^{-j}\) also increase. Reflection \(p(1-t)=1-p(t)\) supplies the other endpoint bounds. The factors \(4/3\) and \(16/9\) account for the change of variable. Plateau derivatives are exactly zero.

`cutoff_jets` splits any argument interval across these regimes and takes an interval hull. Thus a box crossing a cutoff transition remains covered. The function provides values and the first two derivatives, enough for the local stress and convective formulas below.

## A genuinely bounded quadrature rule

`range_quadrature` accepts one-dimensional or tensor-product rational boxes. For each leaf box \(B\), interval arithmetic produces \(I_B\) containing the integrand at **every point** in that box, and therefore
\[
\int_B F\in |B|I_B.
\]
Summing the lower and upper endpoints with directed interval rounding encloses the whole integral. Adaptive bisection only decides where to reduce interval width; it is not used as an error estimator. No omitted tail or interpolation remainder is inferred from sample agreement.

The partition is exact rational arithmetic. Area factors are enclosed independently, and interval summation retains arithmetic rounding. On refinement the implementation subtracts the old lower/upper endpoints separately before adding new endpoint sums, avoiding the artificial uncertainty that subtracting an entire interval would introduce. Approximate floating priorities affect efficiency only.

The trusted numerical implementation is `mpmath.iv` for elementary arithmetic, exponentials, powers and interval endpoints. This is software-backed interval computation, not a Lean/kernel proof of that library. Every result also stores exact dyadic endpoints as `sign`, integer `mantissa`, and integer `exponent`, representing \((-1)^{\rm sign}\,\mathrm{mantissa}\,2^{\rm exponent}\). Those exact bounds, rather than a decimal display, are the reproducible enclosure record.

## The outer swirl integrals and their geometry

The unchanged outer swirl has the form
\[
w=r f(r^2/s_r^2)\left[f((z-d)^2/s_z^2)+f((z+d)^2/s_z^2)\right],
\quad s_r=63/200,\quad s_z=9/20,\quad d=4.
\]
The packets are disjoint. For the positive packet use \(\xi=r^2/s_r^2\in[0,1]\), \(\eta=(z-d)/s_z\in[-1,1]\). Including the azimuthal integral and the reflected packet gives the exact Jacobian \(2\pi s_r^2s_z\,d\xi\,d\eta\). This avoids axis singularities and fixes the integration domain.

Write \(q=r^2\), \(R^2=q+z^2\), \(W=w/r=f(\xi)f(\eta^2)\). The relevant horizontal pressure-tensor components are
\[
Q_{\theta\theta}=\frac{3(4z^2-q)}{4\pi R^7},\qquad
Q_{rr}=\frac{12R^4-105qz^2}{4\pi R^9}.
\]
The implemented exact scalar integrands are
\[
C_v=\int Q_{\theta\theta}qW^2\,dx,
\]
\[
d_v=2\int\left[Q_{rr}W^2+Q_{\theta\theta}\big((W+2qW_q)^2+qW_z^2\big)\right]dx.
\]
The \(Q_{rr}W^2\) term retains the azimuthal basis derivative; it cannot be dropped. All pressure kernels are smooth on these separated supports, and there is no contact term in these outer stress integrals.

At 35 decimal digits and 2048 adaptive rectangles for each integral, the actual obtained enclosures imply the safely rounded bounds
\[
\boxed{C_v\in[5.65032825024,\ 6.97463896085]\times10^{-6},}
\]
\[
\boxed{d_v\in[0.00200534946996,\ 0.00508692054732].}
\]
The exact bounds are in `local-interval-results.json`. These runs took about 14 and 20 seconds respectively, with memory proportional to the panel count. The bounds are intentionally coarse. They already establish, for the original cylindrical pump parameters \(c=2,b=1,\nu=1/100\),
\[
\boxed{\frac{2381}{357}c-4b-\nu\frac{d_v}{C_v}
\in[0.3360594190,\ 6.4637337450]\subset(0,\infty).}
\]
Together with the independently derived multiplier inequality, this certifies positivity of the initial **outer swirl pressure ratio** derivative for any nonzero swirl amplitude. It is not a bound for the full outer pressure derivative; meridional backreaction remains a separate term. The result continues to be useful if a new annular pump uses the same swirl datum, but its pump multiplier must be checked for that datum independently.

## Other local integrals and radial moments

The same file provides `pump_integrand('Cp')`, `pump_integrand('dP')`, and `pump_integrand('local_pzz_convection')` for the original cylindrical meridional field. They use regular variables \(u_r=rH(r^2,z)\), \(u_\theta=rW(r^2,z)\), \(u_z=Z(r^2,z)\), so all apparent axis divisions cancel analytically. The last callback encloses exactly
\[
2\int Q_{ij}u_i\big((u\cdot\nabla)u\big)_j\,dx,
\]
the outer contribution to \(p_{zz}'\) from acceleration \(-(u\cdot\nabla)u\). It does not include the separate nonlocal pressure-acceleration integral. These callbacks have been implemented but no expensive or tight \(C_P,d_P\) certification was run, following the decision to investigate an annular pump with possible exact cancellations.

For annular radial moments, `range_quadrature` also accepts a single list of rational panel boundaries and an interval-valued one-variable callback. Products of this unchanged cutoff and its derivatives can therefore be bounded directly. Split at known annular boundaries and cancel any removable singularity algebraically before interval evaluation. This gives a direct route to enclosing finite-degree angular reductions and radial moments without replacing the smooth profile.

The simple range method has first-order box-width convergence and can become expensive for very narrow relative bounds in two dimensions. Finite angular reduction to one dimension is preferable when exact. Higher-order validated Taylor-panel rules would require explicit bounds for their remainder derivatives; ordinary Simpson/Gauss convergence alone is not substituted for those bounds here.

`verify_interval_cutoff.py` passed independent checks: exact plateau/midpoint jets, enclosure of the known mass \(\int_0^1f=5/8\), and enclosure of \(\int_0^1f'=-1\). These tests check implementation consistency; the enclosure argument is the analytic cutoff bounds plus interval inclusion on every panel.

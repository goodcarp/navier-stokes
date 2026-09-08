# Independent audit of the initial mean-angular-momentum acceleration certificate

**Result: accepted within its stated initial-jet scope.** The code in
`certify_fixed_circle_acceleration.py`, the stored exact endpoints in
`fixed-circle-acceleration-certificate.json`, and equations (3)–(9) of
`covariance-evolution.md` are consistent. This audit uses the independently
established whole-space pressure premise
`|p_m,r(r_*,±4)−P′(r_*)|<5` from
`covariance-evolution-pressure-slab.md`. It does not replace that premise
by a numerical pressure approximation.

For the actual ordinary NS datum and every `3/4≤lambda≤1`, with its
exact retuned amplitude lying in `780<A_lambda<1020`, the certified
conclusion is

\[
 G_{tt}(0,r_*,\pm4)<-293333<-250000.
\]

The stored interval is an enclosure of an **upper-bound expression**.
Its lower endpoint is not asserted to be a lower bound on the actual
second derivative.

## Algebra and the actual pressure term

An independent symbolic recurrence, using `p′=p(1−p)z′`, verifies all four
logistic composition derivatives, including the fourth-order polynomial
`1−14p+36p²−24p³`. The derivatives of
`z=1/t−1/(1−t)` and the factors `(4/3)^j` are correct. Independent chain
rule differentiation also verifies every displayed derivative of
`g(r)=r² chi(r²/a²)` through order four.

Direct expansion of the transport–diffusion operator gives

\[
 L_0^2g-\left[c^2r^2g''-2c\nu rg'''
 +\nu^2(g''''-2g'''/r+3g''/r^2)\right]
 =(c^2r-3\nu^2/r^3)g',
\]

which vanishes at the exact maximum. Combining `L_0 T_0` with the
initial covariance derivative cancels the two `c T_0` terms and gives
the viscosity factor `(2m²+4)/r²−2k²`. The mean-shear contribution
reduces to `lambda² A m² g''/(2r²)` at `g′=0`.

The mixed source can also be checked directly by taking the cylindrical
divergence of the mixed base/seed convective acceleration. On the flat
envelopes it is exactly

\[
 e^{-ikr}s_m=\frac2r[ikV_0'-k^2V_0-(m^2/r)V_0'].
\]

Its contribution through `Re[−ikr e^(−ikr)s_m]/2` is `k²V_0′`.
At the exact root `V_0′=−C`, yielding the code's `−k² C` term, with
the indicated sign and no missing factor two. The pressure reduction
retains the axial covariance flux; dropping that flux would invalidate
the Poisson substitution. The real/imaginary multiplication implementing
`Re[e^(−ikr)(k²r+m²/r+ik)P′]` is correct.

## Root isolation and Green quadrature

The prior analytic argument puts every positive global radial maximum
in `(a/2,a sqrt(5/8))`, a subset of `(a/2,4a/5)`. The new interval
assertions cover that entire range: positive `g′` up to `0.2`, negative
`g′` from `0.23`, and strictly negative `g″` throughout `[0.2,0.23]`.
Endpoint signs isolate the unique root in

\[
 0.2154755753334<r_*<0.2154755753336.
\]

I reran those interval assertions successfully. The plateau condition
`e=1,e′=0` holds on the entire root enclosure. Fourth-order cutoff
formulas are used only on this enclosure, strictly inside the transition;
there is no untreated singular endpoint in that calculation.

The two Green moments integrate all of the support of `C′`. The left
moment ends at the lower root bound, the right moment starts at the upper
bound, and each missing partial cell is enclosed by its integrand range
times `[0,RHI−RLO]`. These contain the true variable-limit integrals,
without assuming independence or a sampled integrand. Outside
`[a/2,a]`, `C′=0`. The exact derivative of the Green formula cancels the
local `C′` boundary terms and gives

\[
 P'=-40iCe^{20ir}+12r^{-5}I_+-20r^3I_-.
\]

The pressure-error contribution is at most
`5 sqrt((k²r+m²/r)²+k²)/2` in the normalized seed coefficient. Its
interval upper endpoint, the transport–diffusion coefficient, and the
viscous torque factor are all strictly negative. Thus substituting the
*smallest* amplitudes `A=780` and `lambda²=9/16` gives an upper bound;
discarding the negative viscous torque also gives an upper bound.
The direction of every one of these endpoint choices is correct.

Exact rational decoding of the stored dyadic upper endpoint verifies
`Gtt_upper < −293333`. No new quadrature run was needed for this audit.
The trust base remains the displayed analytic identities and the
outward-rounded `mpmath.iv` implementation, rather than a formal proof
assistant or a convergence comparison.

## Additional fixed-height radial critical-branch calculation

The code now also records a valid implicit-function calculation. Since
`G_rr=A g″<0`, the initial radial critical point continues locally at
fixed `z=4`. For this branch,

\[
 r_c'=-G_{tr}/G_{rr},\qquad
 \frac{d^2}{dt^2}G(t,r_c(t),4)\big|_{t=0}
 =G_{tt}-G_{tr}^2/G_{rr}.
\]

The code's correction is therefore positive. Its formulas for `G_tr`
and `r_c′` are correct, and the exact stored upper endpoint after that
correction is below `−291459`. The interval bounds `85<G_t<186` and
`r_c′>1/5` are consistent with the same initial calculation.

This is a radial critical branch at a fixed height. It is not the
unrestricted global maximum over both radius and height, and it is not a
material circle. The code's final scope string is conservative but should
explicitly mention this additional fixed-height branch if these fields
are presented to a reader. Neither negative initial second derivative
proves a duration of gain, a sign reversal before `10^(−3)`, finite
cessation of the gain, or failure of a later regenerative mechanism.
Those conclusions need positive-time control or a justified remainder.

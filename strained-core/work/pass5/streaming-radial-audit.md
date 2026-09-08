# Independent audit of the cumulative radial interval certificate

The streaming construction in `certify_radial_strain.py` gives a valid interval enclosure of the radial cubic functional in `../pass4/core-pressure-derivative.md`. This audit checks the cumulative quadrature and the stated profile substitution. It assumes the separately derived analytic reduction and the outward-rounding behavior of `mpmath.iv`; it is not a formal proof of the full fluid argument.

## Profile and endpoints

For the unchanged cutoff `chi`, the profile is

\[
\Phi(s)=\chi(s)-c[1-\chi(4s/25)]\chi(v(s)),
\qquad v(s)=\tfrac14+\tfrac3{44}(s-25),\quad c=\tfrac75.
\]

Its three transitions are exactly `[1/4,1]`, `[25/16,25/4]`, and `[25,36]`. On these intervals respectively,

\[
(\Phi,\Phi')=(\chi,\chi'),\quad
(-c(1-\chi),\tfrac{4c}{25}\chi'),\quad
(-c\chi,-\tfrac{3c}{44}\chi'),
\]

where the last two cutoff arguments are `4s/25` and `v(s)`. Each affine argument maps its transition to `[1/4,1]`; reusing the same rationally partitioned cutoff jets is exact. The profile is zero on `[1,25/16]`, equals `-c` on `[25/4,25]`, and vanishes for `s>=36`. Therefore only the constant plateau changes `Q` between transitions, by `-c(25-25/4)`. It changes neither `K` nor the cubic integral.

For

\[
Q'=\Phi,\qquad
K'=s^{7/2}\big(\Phi\Phi'+\tfrac79s(\Phi')^2\big),
\]

the starting values `Q(1/4)=1/4` and `K(1/4)=0` are exact because `Phi=1` near zero. Every term in the remaining local cubic integrand contains `Phi'`, so there is no omitted integral before `1/4` or after `36`. The constant `340/49` already contains the boundary/contact contribution in the analytic identity; no new contact adjustment is introduced by the annulus. The identity does not require that `Phi` be nonnegative or monotone.

## Inclusion argument for each panel

Suppose intervals `Q_a,K_a` contain the exact cumulative values at the left endpoint `a` of `[a,a+h]`, and `F,H` enclose the ranges of `Phi,K'` on that panel. For every `s` in the panel,

\[
Q(s)\in Q_a+[0,h]F,\qquad K(s)\in K_a+[0,h]H.
\]

Indeed, the integral of a function with range `[L,U]` on a subinterval of length `ell` belongs to `ell[L,U]`, and `ell` lies in `[0,h]`. This argument permits either sign and does not assume independence of `Q`, `K`, `Phi`, and `Phi'`. Evaluating the local expression

\[
-\frac{18368}{735}s\Phi(\Phi')^2
-\frac{448}{105}s^2(\Phi')^3
-\frac{256}{105}Q(\Phi')^2
+\frac{1536}{49}s^{-7/2}\Phi'K
\]

with these interval ranges and multiplying by `h` encloses its panel integral. The next endpoint satisfies `Q_b in Q_a+hF` and `K_b in K_a+hH`. The implementation uses these fixed-length updates and retains the complete cumulative interval. Induction proves the final enclosure. Loss of correlations widens the result but cannot invalidate it. Since every panel begins at positive `s>=1/4`, no singular endpoint estimate for `s^(-7/2)` is needed.

The saved 1024-panel-per-transition result encloses the radial meridional cubic coefficient in `[36.077859259196359...,65.987726316817685...]`. The authoritative endpoints are the exact dyadic tuples in `radial-strain-interval-1024.json`. This statement alone says nothing about the full pressure gate or its persistence under evolution.

An optional exact tightening is `Q(1)=5/8`, following the cutoff reflection symmetry. The existing wider cumulative interval remains valid without that improvement.

## Separate core mixed coefficient

`certify_core_twb.py` applies the same actual cutoff range engine to

\[
t_{Wb}=\frac{44}{35}-\frac{16}{5}
\int_{1/4}^{1}\chi'(s)
\big[\chi(s)+\tfrac23s\chi'(s)\big]^2\,ds.
\]

At 35 decimal digits, 856 adaptive panels give

\[
2.53008839249903430705<t_{Wb}<2.59708098105375099748<4.
\]

The exact enclosing endpoints saved in `core-twb-interval.json` are

\[
\big[420383040390021368815411471639916255\,2^{-117},
863028186834300127712537570847562449\,2^{-118}\big].
\]

These certify the scalar coefficient for the unchanged smooth profile, conditional on the analytic identity and trusted interval arithmetic. They do not certify a full initial feedback sign by themselves, a positive-time invariant class, or repeated scale transfer.

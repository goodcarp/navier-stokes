# PASS6: the axisymmetric escape requires a second scale or loss of stage uniformity

2026-09-08. This note concerns the **same actual smooth unforced NS solution**, initially the compact odd axisymmetric PASS5 datum, at fixed physical viscosity `nu>0`. It does not assert regularity for arbitrary axisymmetric solutions. The new result below is a conditional obstruction to a bounded, smooth, normalized return class, including a meridional-dominant class.

## 1. The normalization and the missing smallness

Let `L_0=1`, `L_(m+1)=q_m L_m`, and set

`V_m(y,s)=L_m^(1+alpha) u(t_m+L_m^(2+alpha)s,L_m y)`.

The actual duration is `t_(m+1)-t_m=L_m^(2+alpha) delta_m`. These fixed stage changes of variables give viscosity `nu_m=L_m^alpha nu` and the exact successor relation

`V_(m+1)(y,0)=q_m^(1+alpha) V_m(q_m y,delta_m)`.

Write `r=sqrt(y_1^2+y_2^2)`, `Gamma=r V_theta`, `a=V_theta/r`, and `xi=omega_theta/r`, where `omega_theta=partial_z V_r-partial_r V_z`. Smoothness defines the quotient values on the axis. Their successor factors are respectively `q^alpha`, `q^(2+alpha)`, and `q^(3+alpha)`.

The axisymmetric angular-momentum maximum principle gives

`G_m(s):=||Gamma_m(s)||_infinity <= L_m^alpha G_0`.

Thus `Gamma_m` tends to zero for `alpha>0`, but **`G_m/nu_m <= G_0/nu` does not tend to zero by scaling**. Small normalized swirl alone cannot invoke a small-swirl regularity theorem at fixed viscosity. Lei–Zhang explicitly allow data-dependent smallness in their small-swirl result. Their equations (1-3), (1-4), and (1-8) also give the angular-momentum and `xi` equations used here. [Lei–Zhang, *Criticality of the axially symmetric Navier–Stokes equations*, PJM 2017](https://msp.org/pjm/2017/289-1/pjm-v289-n1-p06-s.pdf).

The separate `circulation-return-obstruction.md` already excludes a fixed positive rotational receiver. The following result tests removing that lower bound and allowing the receiver energy to become meridional.

## 2. A bounded meridional-return theorem

**Proposition.** Assume `alpha>0`, `0<q_m<=q_+<1`, `0<delta_m<=delta_+<infinity`, and that the entire normalized stages satisfy a common global Cartesian bound

`sup_m sup_(0<=s<=delta_m) max_(j, |beta|<=2) ||partial^beta V_(m,j)(s)||_infinity <= M < infinity`.

Assume all profiles are odd under `y -> -y`. Then

`V_m(.,0) -> 0 locally uniformly as m -> infinity`.

In particular, these stages cannot keep a nonzero fixed-radius receiver RMS. This applies even if the central rotational coefficient is allowed to tend to zero. It needs the **global bounds throughout each stage**; a C2 bound just inside the receiving ball is insufficient.

### Proof: a small angular momentum gives a small swirl source

At azimuth zero, axisymmetry gives `V_y(r,0,z)=V_theta(r,z)` and `V_y(0,0,z)=0`. Consequently

`a(r,z)=integral_0^1 partial_x V_y(tr,0,z) dt`,

`|partial_r a|<=M/2`, and `|partial_z a|<=M`.

For any `h>0`, `r>=h` gives `|a(r,z)|<=G/h^2`; for `0<=r<=h`, comparison with `a(h,z)` gives `|a(r,z)|<=G/h^2+Mh/2`. Choosing `h=(G/M)^(1/3)` yields

`||a||_infinity <= (3/2) M^(2/3) G^(1/3)`,

`S:=||partial_z(a^2)||_infinity <= 3 M^(5/3) G^(1/3)`.       (1)

The zero cases are immediate. All bounds are global in the axial variable. No a priori derivative bound for `Gamma` is being silently substituted for this interpolation step.

The exact meridional-vorticity equation is

`(partial_s+V_r partial_r+V_z partial_z) xi`

`   =nu_m(partial_rr+3/r partial_r+partial_zz)xi+partial_z(a^2)`.

Its scalar maximum principle therefore gives

`||xi_m(delta_m)||_infinity <= X_m+integral_0^delta_m S_m(s) ds`,

where `X_m=||xi_m(0)||_infinity`. The smooth axis is the regular radial origin for the displayed diffusion operator. Bounded C2 also makes `xi` bounded: `omega_theta` vanishes on the axis and its radial derivative is bounded by `2M`.

After the actual successor change of variables,

`X_(m+1) <= q_m^(3+alpha) [X_m+3 delta_+ M^(5/3) G_0^(1/3) L_m^(alpha/3)]`.       (2)

Put `Q=q_+^(3+alpha)`, `R=q_+^(alpha/3)`, and `D=3 delta_+ M^(5/3) G_0^(1/3)`. Then `0<Q<R<1` and

`X_m <= Q^m X_0+Q D (R^m-Q^m)/(R-Q) -> 0`.       (3)

Finally, the common global C2 bound gives subsequential C1 convergence on compact sets and a globally bounded limit. Equation (1) removes the swirl component. Equation (3) removes the remaining curl, since `omega_theta=r xi`; divergence also remains zero. Each Cartesian component of the limit is therefore an entire bounded harmonic function, hence constant. Oddness sets that constant to zero. Every convergent subsequence has this same limit, proving the stated local convergence and contradicting a fixed nonzero receiver RMS. This last step uses global boundedness to exclude nonconstant affine or harmonic limiting flows.

This is a self-contained finite-profile obstruction. It does not use an unproved assertion that small swirl always implies global regularity, or an assertion about all bounded ancient axisymmetric solutions.

## 3. What an axisymmetric meridional mechanism must actually change

Deleting the lower bound on normalized rotation is consequently not enough. To remain axisymmetric with `alpha>0`, a proposed return construction must lose at least one hypothesis of the proposition: a global normalized C2 bound, a global velocity bound, bounded normalized durations, a uniform contraction, or a nonzero fixed-radius receiver. The latter three changes also alter the intended cascade and must be tracked explicitly.

There is a directly testable source budget. For any stage, define

`I_m=integral_0^delta_m ||partial_z[(V_theta/r)^2]||_infinity ds`.

The exact inequality is

`X_(m+1) <= q_m^(3+alpha)(X_m+I_m)`.       (4)

For a candidate class that maintains `X_m>=x_*>0`, (4) implies `I_m` cannot tend to zero. At an individual endpoint the necessary bound is `I_m>=q_m^(-3-alpha)X_(m+1)-X_m`; its right-hand side can be negative for one stage, so only a positive right-hand side creates a positive one-stage requirement. Iterating the inequality gives the stronger cumulative test without assuming identical endpoint norms.

If `M_m` bounds the full-stage C2 norm, (1) gives

`I_m <= 3 delta_m M_m^(5/3) G_0^(1/3) L_m^(alpha/3)`.

Under bounded durations and a persistent positive `xi` return, this forces a subsequence with

`M_m >= c L_m^(-alpha/5)`       (5)

for some `c>0`. This deliberately weak lower bound is already incompatible with the proposed bounded smooth profile class. It is a necessary condition, not a sufficient instability criterion.

A concrete next eligible axisymmetric test is therefore a **second radial scale**, not a fresh affine core: track a shrinking normalized swirl layer on the actual evolved solution and verify (4) while retaining nonzero meridional receiver energy. For example the diagnostic scaling

`h_m^2=nu_m/lambda`, `a_m(r,z)=A_m(r/h_m,z)`

with fixed `lambda>0` gives `Gamma_m=h_m^2 R^2 A_m`, `V_theta=h_m R A_m`, and radial diffusion rate `nu_m/h_m^2=lambda`. It is consistent with `Gamma_m=O(nu_m)`, but viscosity remains leading order in the layer. For axially varying `A_m=O(1)`, `partial_z(A_m^2)` can remain order one while normalized swirl velocity vanishes. This algebra avoids the fixed-profile obstruction by allowing unbounded radial derivatives.

It does **not** prove the required return: an order-one source confined to a shrinking tube may transfer too little meridional velocity outside that tube. The decisive next lemma would have to prove, for the actual inherited layer and actual viscosity, both a positive lower bound on the source budget in (4) and a nonvanishing fixed-radius meridional receiver after the step. If only the source supremum survives, that is insufficient. The existing PASS5 initial jets provide neither statement. Nonaxisymmetric torque production is a separate eligible route analyzed by the other PASS6 note.

## 4. Setting alpha to zero also changes the target

At `alpha=0`, normalized viscosity and `Gamma` are unchanged by the successor scaling; the angular-momentum contraction obstruction disappears. Exact receiver matching then keeps `Re=L U/nu` constant. Physical amplitude could still grow like `L^(-1)`, so constant Reynolds number alone is not a regularity argument.

But a bounded full-profile critical return has an independent obstruction. Suppose `q_m<=q_+<1`, `delta_m<=delta_+`, and the **full normalized velocity**, at all points of every stage, is bounded by `M`. The accumulating time satisfies, for `t` in stage `m`,

`T-t <= delta_+ L_m^2/(1-q_+^2)`,

and hence

`||u(t)||_infinity <= M/L_m <= M sqrt[delta_+/(1-q_+^2)] (T-t)^(-1/2)`.       (6)

No lower duration bound is needed for this implication. A bound only on the receiver RMS, or only at discrete endpoint times, would not establish (6).

For strong axisymmetric solutions on `R^3`, Chen–Strain–Tsai–Yau Theorem 1.1 excludes a singularity under this full Type I velocity bound. Their hypotheses include initial `H^(1/2)` data, bounded initial `r u_theta`, and `p in L^(5/3)` in spacetime; the smooth compact initial datum and the finite-energy estimates supply these on a finite interval. Fixed `nu>0` reduces to their unit viscosity by a constant time/amplitude change. [Chen–Strain–Tsai–Yau, *Lower bounds on the blow-up rate of the axisymmetric Navier–Stokes equations II*, CPDE 2009, Theorem 1.1](https://arxiv.org/pdf/0709.4230).

Thus an axisymmetric critical mechanism must be Type II: a real violation of the full velocity bound (6), through unbounded normalized stage amplitudes or incompatible stage timing. Merely relabeling the PASS5 finite gain with `alpha=0` does not supply that mechanism.

## Scope and verification

The fixed rotational return is excluded by the separate angular-momentum memo. Here the stronger meridional proposition excludes a globally bounded C2 class for `alpha>0`; the cited Type I theorem excludes the globally bounded full-velocity class for `alpha=0`. Neither excludes all axisymmetric blowup scenarios. The remaining concrete axisymmetric target is a dynamically inherited second-scale layer satisfying the meridional source and transfer budget, with viscosity retained at that scale.

`python3 checks/verify_escape_route_audit.py` checks the exact quotient equations, normalization powers, interpolation constants, and recurrence sum. The compactness/Liouville argument and the application of the cited theorem are analytical steps, not computer-certified existence statements.

# PASS4: an iterated radius rule is a level set, not yet a dynamical cascade

2026-09-08. This analyzes repeated application of the PASS3 receiver rule to **one and the same** smooth NS solution. The center remains the origin, fixed by inversion symmetry; the receiver weight and exponent `a` are fixed throughout. No new initial datum is inserted.

## 1. Exact telescoping on the same solution

Let `chi` be the fixed nonnegative radial receiver weight, and define

`U(r,t)^2 = [integral chi(|y|^2)|u(t,r y)|^2 dy]/[integral chi(|y|^2)dy]`,

`Re(r,t)=r U(r,t)/nu`.

Suppose successive receiving times `t_m` and radii `L_m` satisfy `L_(m+1)=q_m L_m`, `0<q_m<1`, and the exact matching rule

`q_m^a Re(L_(m+1),t_(m+1))/Re(L_m,t_m)=1`.

Define `F(r,t)=r^a Re(r,t)=r^(a+1)U(r,t)/nu`. Then every selected pair lies on the **same level set**:

`F(L_m,t_m)=F(L_0,t_0)=C_0>0`.

In particular,

`Re_m=Re_0 (L_0/L_m)^a`,

`U_m=U_0 (L_0/L_m)^(1+a)`.

Thus repeating the scalar selection does not independently construct new amplification dynamics. It labels a constant level of an observable of the already-evolving field. If this level set really reaches arbitrarily small radii in finite time, the actual solution must supply the following diverging amplitudes and gradients.

## 2. A rigorous lower bound on radius while the solution stays regular

Put

`c_1^2=[integral chi(|y|^2)|y|^2dy]/[integral chi(|y|^2)dy]`.

Since `u(t,0)=0`, the mean-value theorem gives

`U(r,t)<=c_1 r M(r,t)`,

where `M(r,t)=sup_(|x|<=r)||grad u(t,x)||_op`. Therefore every selected radius must obey

`M(L_m,t_m) >= U_0 L_0^(1+a)/(c_1 L_m^(2+a))`.

On any compact time interval `[t_0,T_1]` where the same solution is smooth, set

`M_*=sup_{t_0<=t<=T_1, |x|<=L_0}||grad u(t,x)||_op < infinity`.

Then

`L_m >= [U_0 L_0^(1+a)/(c_1 M_*)]^(1/(2+a)) > 0`.

This rules out an artificial cascade produced solely by measuring successively smaller balls during a regular compact time interval. Equivalently, the total selected logarithmic contraction is bounded:

`sum_(j<m) log(1/q_j) <= [1/(2+a)] log(c_1 M_* L_0/U_0)`.

For the initial rigid core, `U_0=c_chi Omega L_0` and `c_chi=sqrt(2/3)c_1`, so the lower radius fraction is

`L_m/L_0 >= [sqrt(2/3) Omega/M_*]^(1/(2+a))`.

Even without a fixed zero-velocity center, bounded velocity would imply `Re(r,t)<=r||u(t)||_infty/nu ->0` as `r->0`, contradicting the required increasing Reynolds number. The symmetry gives the stronger gradient bound above and removes irrelevant uniform translation from this measurement.

For a smooth field with nonzero central gradient, radial averaging gives the sharper local asymptotic

`U(r,t)=r sqrt(c_1^2/3) ||grad u(t,0)||_F + O(r^3)`.

The odd symmetry removes the quadratic velocity Taylor term; the same RMS order follows from cancellation of odd radial moments under appropriate smoothness. Thus, near a regular point, the level set essentially obeys `r^(a+2)||grad u(t,0)||_F=constant`. It tracks actual gradient amplification.

## 3. What time and amplitude behavior a real iteration would require

If `L_m->0`, exact matching forces

`U_m ~ L_m^(-1-a) -> infinity`,

`||grad u(t_m)||_infty >= constant L_m^(-2-a) -> infinity`.

Because `U_m` is an RMS, the full velocity supremum is at least `U_m`. This is a genuine amplitude-growth requirement, unlike the bounded-velocity Euler construction.

Finite-time blowup additionally requires that the receiving times stay bounded. The matching rule alone provides no bound on `t_(m+1)-t_m`. If a **separate proved stage theorem** supplied

`delta_- L_m/U_m <= t_(m+1)-t_m <= delta_+ L_m/U_m`,

with `0<delta_-<=delta_+<infinity`, and a uniform contraction `q_m<=q_+<1`, then

`T-t_m ~ L_m^(2+a)`,

`U_m ~ (T-t_m)^[-(1+a)/(2+a)]`,

`||grad u(t_m)||_infty >= constant/(T-t_m)`.

Here the comparison constants depend on the stated duration and contraction bounds. The resulting velocity exponent exceeds `1/2` for `a>0`. These conclusions are conditional on the actual dynamical duration bounds; assigning turnover times to selected radii does not prove them. An infinite sequence with `prod q_m>0` has no vanishing length scale, and a sequence with `t_m->infinity` is not a finite-time blowup proof.

The PASS3 local gain and implicit radius branch do not themselves supply a uniform terminal class. In fact their small-time branch has `q(t)=1-O(t^2)` and a leading strain change of order `t`. On a regular finite interval any repeated choices must obey the positive radius lower bound in Section 2.

## 4. A stronger diagnostic: full RMS normalization reduces the central rotation fraction

For the first PASS3 episode, use its exact energy expansion at fixed radius `r`:

`beta_r=(1/2)[K+Q_r/(Omega^2 D_r)]`,

`D_r=integral chi_r|Jx|^2`, `Q_r=integral chi_r|grad H|^2`.

The matched branch satisfies

`q(t)=1-[beta_r/(a+2)]t^2+O(t^3)`.

The comparison state at the original length/amplitude is

`v(t,y)=q(t)^(1+a)u(t,q(t)y)`.

Its central axial vorticity is exactly multiplied by `q^(a+2)`. Since

`omega_z(t,0)=2Omega[1+(K/2)t^2]+O(t^3)`,

we obtain

`omega_z[v](t,0)/(2Omega)`

`   =1-[Q_r/(2Omega^2D_r)]t^2+O(t^3)`.

For `K>0`, `H` is a nonconstant harmonic function and `Q_r>0` for every nonzero nonnegative receiver weight. Thus the central rotation of the RMS-normalized successor **decreases at leading order**, even though the physical central vorticity and the selected receiver Reynolds number increase. The additional irrotational strain energy uses part of the RMS normalization.

Meanwhile its symmetric gradient remains `t Hess H(0)+O(t^2)`. A class requiring the original rotation fraction together with zero strain is therefore not invariant. Admitting strain can avoid the zero-strain restriction, but RMS normalization alone does not stop rotational strength from leaking into other components of the velocity.

This statement is a local expansion for the same actual solution, not a sign theorem for all later stages or all definitions of core quality.

## 5. Exact normalized matrix dynamics for a class that admits strain

Define the normalized stage field with fixed stage scales

`V_m(y,s)=u(t_m+(L_m/U_m)s,L_m y)/U_m`.

It solves the full NS equation with viscosity `Re_m^(-1)`. Let

`A(s)=grad_y V_m(s,0)`,

`Pi(s)=Hess_y P_m(s,0)`,

`C(s)=Delta_y grad_y V_m(s,0)`.

At the fixed material origin the exact matrix evolution is

`A'=-A^2-Pi+Re_m^(-1) C`,

`tr A=0`, `tr Pi=-tr(A^2)`.

The pressure Hessian and the higher spatial jet are outputs of the full normalized solution. They are not determined by `A` alone.

The reflection symmetry across `z=0` makes the central gradient block diagonal:

`A=diag_block(B,sigma)`, `tr B=-sigma`.

Write `w=(B_21-B_12)/2=omega_z[V_m](0)/2`, assume `w>0`, and denote the viscous jet terms by

`eta_w=Re_m^(-1) Delta_y w(0)`,

`eta_sigma=Re_m^(-1) Delta_y A_zz(0)`.

Without assuming axisymmetry or a circular transverse strain, the exact scalar equations are

`w'=sigma w+eta_w`,

`sigma'=-sigma^2-Pi_zz+eta_sigma`.

The first equation follows either from the normal vorticity equation on the invariant plane or from the two-by-two identity `skew(B^2)=(tr B)skew(B)`. For the dimensionless strain-to-rotation ratio `R=sigma/w`,

`R'=[-2sigma^2-Pi_zz+eta_sigma-R eta_w]/w`.

At the receiving time, the normalized successor matrix is

`A_next=q_m^(a+2) A_end`.

Therefore `sigma_next/w_next=sigma_end/w_end`: isotropic state rescaling cannot repair an unfavorable strain-to-rotation ratio. The transverse symmetric deviator divided by `w` likewise retains its ratio under this rescaling.

## 6. An actionable profile-class closure test

Here are exact tests a strain-admitting class would have to satisfy. They are more restrictive than the scalar RMS gain and can be checked on a candidate finite-stage construction.

First choose explicit constants `w_->0`, `w_+`, and `R_+>0`, plus bounds on transverse anisotropy and the required normalized spatial jets. Admit normalized initial profiles with receiver RMS equal to one, central `w in [w_-,w_+]` and `0<=sigma/w<=R_+`.

For the ratio interval to be invariant during the stage, the exact inward-pointing conditions at its two boundaries are

`at R=0:       -Pi_zz+eta_sigma >=0`,

`at R=R_+:     -Pi_zz+eta_sigma-R_+ eta_w <=2R_+^2 w^2`.

With continuous coefficients and `w>0`, these are the standard boundary viability conditions for the displayed scalar ratio equation; strict margins give a robust version. They must be proved from the actual pressure and viscous jets uniformly over the proposed class. Arbitrarily choosing `Pi_zz` to enforce them would replace the PDE by a prescribed-strain model.

For central rotation to lie in its allowed interval **after** the matched rescaling, the exact endpoint condition is

`w_- <= q_m^(a+2) w_end <= w_+`.

In particular, to return with at least the starting central rotation, one needs

`integral_stage [sigma+eta_w/w] ds >= (a+2)log(1/q_m)`.

The first short PASS3 branch fails that last requirement at leading order by the strictly positive term `Q_r/(2Omega^2D_r)t^2`. It may still lie in a wider class for one episode, but a repeated loss cannot be ignored or replenished by renaming the RMS amplitude.

A genuine uniform return lemma must then assert that, for **every** admitted normalized incoming profile and every allowed Reynolds number, an actual evolution over a duration `s in [delta_-,delta_+]` reaches an allowed radius `q in [q_-,q_+] subset (0,1)` and returns the complete normalized profile

`V_next(y)=q^(1+a)V_end(q y)`

to the same class, with the outer geometry, pressure control and jet bounds included. The ratio inequalities and endpoint rotation budget above are explicit necessary checks (and sufficient for the stated scalar bounds when their hypotheses hold). They do not replace control of the transverse anisotropy or of the full profile.

This supplies a concrete rejection rule: if a candidate stage gains RMS but violates the central rotation budget or drives its normalized strain ratio outside the proposed class, it has not regenerated that class. The present first episode supplies positive receiver gain but does not pass a fixed-rotation return test.

## Verification

`python3 work/pass4/verify_iteration_audit.py` checks telescoping exponents, the normalized rotation-loss coefficient, exact two-by-two matrix identities, strain-ratio evolution, and state-rescaling powers. The analytic gradient lower bound is the mean-value theorem and the material-center symmetry. No regularity failure or profile-class return is asserted without the additional dynamical hypotheses stated above.

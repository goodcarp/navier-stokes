# Independent whole-space reconstruction audit

8 September 2026. This audits the **representation and fitting operators**, not a Navier–Stokes trajectory or an endpoint certificate. The fixed initial field remains `u0=M+1020v+(1/4)w_L`, `nu=1/1000`. No earlier pass or published source is modified.

The present prototype can define a compact, exactly divergence-free, axis-regular **H7/C6 velocity approximation** when its coefficients and conjugate symmetries are interpreted as specified. Its degree-nine radial splines are not C-infinity. Neither the sampled nor continuous fit supplies the whole-space H4 residual, a certified approximation error, useful endpoint gains, or a successor class.

## 1. Operator signs and axis orders

Use genuine scalar Fourier coefficients `P_m(r,z)e^(im theta)` and `T_m(r,z)e^(im theta)`. For

`v=curl(T e_z)+curl curl(P e_z)`,

exact Cartesian differentiation gives

```
v_r = P_rz + (im/r) T,
v_theta = (im/r) P_z - T_r,
v_z = -L_m P,
L_m = partial_rr + r^-1 partial_r - m²/r².
```

The divergence vanishes identically, and `omega_z=-L_m T`. With `s=r²`, scalar factors `P_m=r^n p_m(s,z)`, `T_m=r^n t_m(s,z)`, `n=|m|`, obey

```
L_m[r^n f(s,z)] = 4 r^n [(n+1) f_s+s f_ss],
(partial_r ± m/r)[r^n f] = (n±m)r^(n-1)f+2r^(n+1)f_s.
```

Terms with `n±m=0` should be removed symbolically. The resulting helical velocities have orders

```
U_plus = v_r+i v_theta = r^|m+1| F_plus(r²,z),
U_minus= v_r-i v_theta = r^|m-1| F_minus(r²,z),
v_z = r^|m| F_z(r²,z).
```

Thus positive modes 0,4,16 require orders `(1,1,0)`, `(5,3,4)`, `(17,15,16)` respectively. These are stronger than zero radial flux or odd/even component parity. An instructive failure is `T=r^5 cos(4 theta)`: its velocity is homogeneous of degree four but not a Cartesian polynomial. A nonzero fifth velocity derivative is homogeneous of degree −1, so its square has a logarithmically divergent integral across the axis. The checker obtains `partial_x^5 v_y(1,1)=945 sqrt(2)/16`. This field looks small at the axis but cannot enter an H4 viscous-residual calculation requiring velocity H6.

The unchanged initial formulas do meet the stronger orders. Its nonzero seed mode even vanishes on an open axis neighborhood. The antiderivative potentials in [the pass11 design](../pass11/poloidal-toroidal-evolution-design.md) recover the exact compact datum, including the scalar Fourier factor `lambda/2`; amplitudes of real cosine waves must not be stored as positive Fourier coefficients without this half factor.

## 2. Exact axial localization terms

For a fixed scalar axial window q(z), the full identity is

```
v(qP,qT) = q v(P,T) + w,
w=q'(z) grad_h P,
```

or componentwise

```
v_r,new=q v_r+q' P_r,
v_theta,new=q v_theta+(im/r)q'P,
v_z,new=q v_z.
```

There is no toroidal q' term. These poloidal terms restore the divergence that multiplying velocity by q would create. The prototype's `evaluate_increment` and staggered joining expression have the correct signs and component locations.

The join must also enter derivatives and convection. For an unwindowed increment a=v(P,T), a fixed q gives

```
Delta(q a+w)=q Delta a+2q' partial_z a+q'' a
              +q' grad_h Delta P+2q'' grad_h P_z+q''' grad_h P,
partial_t(q a+w)=q a_t+q' grad_h P_t,
((q a+w).grad)(q a+w)
 =q²(a.grad)a+q q' a_z a
   +q(a.grad)w+q(w.grad)a+(w.grad)w.
```

For the selected reconstruction the complete velocity is `u0+q a+w`; all base/increment cross terms must be added when computing NS convection and pressure. Applying these formulas only to the increment does not replace u0 by a windowed field. A field retained through angular mode 16 generates quadratic residual/source modes through 32. Discarding those generated modes would not evaluate its full NS residual.

The Fourier synthesis is `exp(i k(z-z_origin))`, with the actual stored origin. This matters on odd grids. Negative axial frequencies of a positive angular mode are independent complex data; reality instead relates both signed indices together. Mode zero must remain real. Odd/C4 averaging uses signed coordinate permutations and can enforce the exact datum's symmetry, but its effect must be included in the reconstruction if not already preserved.

## 3. Spline support, gauge and Sobolev requirements

The prototype uses `r^m B(r²)` with degree-nine simple-knot splines, tenfold clamping at the endpoints, and the last nine basis functions removed. At the outer endpoint the removed functions are exactly those with vanishing orders zero through eight. The last retained function has order nine. Therefore the zero extension has eight continuous derivatives and weak derivatives through order nine in L2. Multiplication by r^m is harmless on bounded support, and clamping at s=0 imposes no additional vanishing requirement on the smooth factor p(s,z).

With a smooth compact axial window and finitely many Fourier modes, both potentials lie in H9 and C8 in Cartesian space. The velocity is therefore H7 and C6. This is more than the sufficient spatial condition `P in H8, T in H7`, which gives velocity H6 and hence its viscous term in H4. To obtain a full H4 time residual one also needs, for example, `P_t in H6, T_t in H5`. A saved endpoint has no time derivative and is not a space-time approximation. Ordinary cubic splines do not provide these residual derivatives, and C8 splines cannot be silently promoted to an all-order smooth-force construction.

For m=0, potentials constant in radius represent zero velocity on an unrestricted domain. Compact outer conditions remove the exact constant gauge; they do not prevent nearly constant, weakly resolved directions from being poorly conditioned in a sampled norm. Do not remove the m=k=0 poloidal sector: it represents nontrivial vertical profiles `v_z=-L_0P`. The sampled normal matrix should be assessed as an approximation problem as well as a linear solve.

Every compact solenoidal whole-space field has zero net axial flux through each entire horizontal plane: the flux is z-independent by divergence and zero outside the axial support. A periodic-cylinder state with a nonzero throughflow cannot match this constraint exactly on the same radial support. The selected odd state should have zero throughflow, but that is a check, not an automatic property of arbitrary MAC arrays.

A separate inverse-potential trap is worth retaining even though the **coupled fit avoids it**. For m>0, solving `-L_m P=f` with axis regularity and decay gives

```
P(r)= [r^-m integral_0^r rho^(m+1) f(rho) d rho
        +r^m integral_r^infinity rho^(1-m) f(rho) d rho]/(2m).
```

Compact f generally gives a harmonic `C(z)r^-m` tail, not compact P. For a compact velocity such tails of P and T can cancel through `T=-i C'(z)r^-m` outside. Cutting the potentials there produces new collar velocity, despite the original velocity being zero. For m=0 compact P requires `integral rho f(rho) d rho=0`; a nonzero moment produces a logarithmic potential. Independent inverse solves followed by truncation must pay these terms or enforce the appropriate moments.

## 4. Fitting operators: what was checked, and what failed

In the sampled fit, the P and T columns at axial frequency k are

```
A_P=(ik D_face, -k mB_cell/r, -L_cell),
A_T=(im B_face/r, -D_cell, 0).
```

The kinetic weights are `(W,V,V)`. Conjugating these columns yields the implemented RHS and the **positive** k-dependent cross block

`k [D_face^T W (mB_face/r)+(mB_cell/r)^T V D_cell]`.

Interleaving P/T coefficients gives half-bandwidth 19: degree-nine basis supports overlap only when scalar indices differ by at most nine. Independent full matrix assembly, dense QR, manufactured analytic polynomials and arbitrary off-grid evaluation found no sign/index defect. Core formulas in the continuous observer are also exact: `Omega=1-2 T_s(0,0)` and `b=1-2 P_sz(0,0)`. Consequently the reported extreme central jets in the finer sampled fit are not caused by an omitted energy factor or by the q' sign. They expose an invalid approximation-selection criterion.

There is no uniform implication from small sample error to axis accuracy or continuous Hs accuracy. For example choose epsilon below the first cell radius and

`T=(epsilon²/20)(1-r²/epsilon²)_+^10 q(z)`.

Then `v_theta=r(1-r²/epsilon²)_+^9 q(z)` vanishes at every MAC velocity sample, but its central rotation at q=1 is exactly one. Its horizontal squared L2 norm is `pi epsilon^4/380`; multiplying the potential by any amplitude leaves all samples zero while scaling actual energy quadratically. Its H4 norm is not controlled by the samples. A C-infinity cutoff gives the same counterexample principle. This example is not claimed to lie in a particular fixed finite spline trial space; it proves why a mesh-independent norm assertion needs additional information.

The replacement `continuous_potential_fit.py` uses a specified odd/even cubic interpolant as an L2 target. In true horizontal integration, the P/T cross term is a boundary derivative, so it vanishes. With

```
H_ij=integral r (B_i' B_j'+m² B_i B_j/r²) dr,
Q_ij=integral r (L_m B_i)(L_m B_j) dr,
```

the independent sectors have Gram matrices H and `Q+k²H`. On each union span, their polynomial degrees are at most `2m+35` and `2m+33`; the RHS degrees are at most `m+21`. Thus Gauss order `m+18` integrates them exactly in exact arithmetic, provided every potential and target knot is included. Floating arithmetic is still not an enclosure.

Small dense-QR checks of this replacement pass for m=0,4,16 on a different radial mapping, with maximum velocity discrepancy below `8e-14` and the expected L2 projection contraction. **The cubic target is not a smooth Cartesian field at higher modes merely because its components have parity.** It is an L2 interpolation choice. Continuous L2 projection prevents spurious growth of its unwindowed L2 norm relative to that target, but it does not uniformly control H4 or central jets as the basis grows. The axial q'P collar can also add energy after projection. Continuous derivative norms, or suitable explicitly priced regularization, remain necessary.

## 5. Continuous energy and actual Cartesian derivative norms

For the full signed angular sum,

`||v||_2²=2pi sum_m integral r [ (|U_plus|²+|U_minus|²)/2+|v_z|² ] dr dz`.

Positive angular modes instead have a conjugate factor two. Poloidal and toroidal sectors are orthogonal after horizontal integration. The fluctuation energy of a positive mode is consequently

`2pi integral [T* H T + P_z* H P_z + P* Q P] dz`.

The observer correctly uses this factor and, for the original m=4 field plus its increment, the cross factor `4pi Re integral r conjugate(u_initial,4).delta_u_4`. The full radial Gram is polynomial; the base cross integrals and axial q-window quadrature are not automatically exact. Its selected radial optimizer is not a global-maximum proof.

An exact way to compute Cartesian ordered derivative norms avoids a sum of unrelated cylindrical partial derivatives. Define `D_plus=partial_x+i partial_y` and `D_minus=partial_x-i partial_y`. On a scalar angular coefficient f_m,

```
D_plus:  (m,f_m) -> (m+1, (partial_r-m/r) f_m),
D_minus: (m,f_m) -> (m-1, (partial_r+m/r) f_m).
```

The order changes **at each application**. For a scalar f and integer j,

```
|nabla^j f|_F²
 =sum_(h=0)^j binom(j,h) 2^-h
    sum_(a=0)^h binom(h,a)
      |D_plus^a D_minus^(h-a) partial_z^(j-h) f|².
```

For velocity apply this to the Cartesian helices `v_x+i v_y=e^(i theta)U_plus`, `v_x-i v_y=e^(-i theta)U_minus`, with half weight each, plus v_z. Angular Parseval then gives the exact whole-space integral. The checker independently verifies the identity through j=4 against direct Cartesian polynomial derivatives, including their multinomial counts.

The pass9 derivative-sum H4 norm and the ordered-tensor norm are not identical. They satisfy

`||v||_H4,multi² <= sum_(j=0)^4 ||nabla^j v||_2² <=12 ||v||_H4,multi²`.

The factor 12 is the largest fourth-order Cartesian multinomial coefficient. Missing angular basis derivatives or missing multinomial weights would invalidate a residual estimate even if sampled divergence is tiny.

## 6. Finite energy and a usable exterior-pressure bound

Windowed compact potentials give finite-energy velocity. The actual NS solution need not remain compact. Unwindowed axial Fourier series and finite sums of point-frequency Bessel modes are not full-space L2 fields.

For a reconstruction supported in `r<=R_rad`, `|z|<=Z`, a containing spherical radius is `R_s=sqrt(R_rad²+Z²)` (also include the exact initial support). Using R_rad alone in a Newton tail bound is invalid. For N=1/(4pi|x|), whole-space pressure is `p=partial_i partial_j N*(v_i v_j)`. Outside the velocity support, the complete projected residual equals grad p and is normally nonzero, of order |x|^-4.

The ordered Newton derivative constants can be made explicit:

```
||nabla^n N(x)||_F² = (2n)!/[2^n (4pi)²] |x|^(-2n-2).
```

Indeed rotational invariance and homogeneity give a radial power; since every derivative of N is harmonic away from zero, `Delta ||nabla^nN||²=2||nabla^(n+1)N||²`, giving the recurrence `c_(n+1)=(n+1)(2n+1)c_n`. This proves the formula from n=0.

Set `C_j=sqrt((2j+6)!/2^(j+3))/(4pi)` and let E=||v||²/2 have a rigorous upper bound. For `|x|>R_s`,

`||nabla^j grad p(x)||_F <= 2E C_j/(|x|-R_s)^(j+4)`.

The squared L2 tail outside radius R_*>R_s is at most `4pi(2E C_j)²` times

`h^(-2j-5)/(2j+5)+2R_s h^(-2j-6)/(2j+6)+R_s² h^(-2j-7)/(2j+7)`,

where h=R_*−R_s. Sum j=0,...,4, with the correct verification-length weights if used. This is an explicit conditional tail enclosure, not a certificate from an unenclosed sampled energy. The pressure inside the integration radius, including distributional contact contributions when expressed with singular Hessian kernels, still needs full control.

Anchoring `v=u0+reconstructed[U_MAC(t)-U_MAC(0)]` can make v(0)=u0 exactly. It does not remove the initial discretization change from the subsequent PDE residual: the entire operator N(v), pressure, joins and time derivative must be recomputed on this anchored field. Later endpoints and errors remain inherited.

## 7. Checks and coordinated reference audit

[check_compact_potential_reconstruction.py](check_compact_potential_reconstruction.py) checks independent analytic P/T recovery for m=0,4,16, both axial grid parities, signed phases, axis/outer values, off-grid/collar evaluation, the sparse matrix against dense QR, Cartesian axis orders, the q' commutator and ordered derivative identities. [Its JSON](compact-potential-independent-checks.json) pins the checked revision. [check_continuous_potential_fit.py](check_continuous_potential_fit.py) separately checks the new continuous projection, exact vanishing P/T cross Gram, dense QR, and the invisible-axis example; [its JSON](continuous-potential-independent-checks.json) records the results. No production fit or time trajectory was rerun by these tests.

The [C4-adapted reference](c4-adapted-reference.md) was also read independently. Its central matrix symmetry, time-ordered off-axis comparison, cubic deformation estimate, and exact transformed equations check: the `2Bv` coupling, anisotropic diffusion metric and `-nu K y` linear viscous subtraction all survive. Its new H4-compatible half-Hölder remainder uses the correct `C_H=4sqrt(3)/pi`, Taylor factors 4/15 and 2/3, mixed-power comparison exponents, and central F/phase enclosures. These are conditional analytic results; finite-cylinder histories do not supply their actual-solution hypotheses.

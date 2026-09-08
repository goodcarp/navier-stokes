# An actual-solution-adapted reference for the fixed C4 datum

This note concerns the unchanged whole-space datum `u0=M+1020 v+(1/4)w_L`, viscosity `nu=1/1000`, and its actual unforced smooth Navier–Stokes solution on an interval on which it exists. No amplitude, field, phase, or outer packet is reset. The result is an exact comparison construction and a conditional local deformation bound; no useful stage, endpoint gain, or return is asserted.

Inputs are the fixed-amplitude identities in `work/pass9/fixed-amplitude-initialization.md`, the full Fourier source in `work/pass10/initial-full-field.md`, the exact coordinate caution in `work/pass10/corotating-time-polynomial.md`, and the current `NEXT_TARGET.md`. The structural inspiration is Claude U1's use of a reference driven by the actual solution. The comparison here is different: it retains actual matrix gradients and time ordering, and does not use an axisymmetric sign kernel.

## 1. Assumptions and conventions

Write `A_ij=partial_j u_i`. Assume `u,p` are a classical finite-energy whole-space solution, with all derivatives used below continuous on the stated tubes. Sufficient smooth Sobolev regularity and decay may be imposed instead. All derivative bounds in this note are hypotheses about the **actual** solution; the finite-cylinder diagnostics do not supply them. The pressure is the whole-space pressure, up to a function of time,

`-Delta p = tr[(grad u)^2]`.

The initial datum is equivariant under the quarter-turn Q about the z axis and odd under spatial inversion. Uniqueness preserves

`u(Qx,t)=Q u(x,t)`, `u(-x,t)=-u(x,t)`.

Thus `u(0,t)=0`; pressure can be chosen even and Q-invariant. For `M_j(t)` below, `||D^j u||` means the operator norm of the j-linear derivative tensor, with Euclidean input/output norms. A sup is always over an explicitly prescribed physical neighborhood or tube.

## 2. The actual central matrix and its coupling

The quarter-turn fixes the origin, so `A(t)=grad u(0,t)` commutes with Q. Incompressibility forces

```
A(t) = b(t) S + Omega(t) J,
S=diag(-1,-1,2),             J=[[0,-1,0],[1,0,0],[0,0,0]].                 (2.1)
```

This is an exact statement about the central first jet, including the nonaxisymmetric solution. It does not make the surrounding field axisymmetric. All central matrices commute because `[S,J]=0`; this special simplification is proved by symmetry rather than assumed for a general trajectory.

Differentiating the **full** NS equation at the stationary origin gives

```
A' + A^2 = -H + nu D,
H=Hess p(0,t),                     D=Delta grad u(0,t).                  (2.2)
```

Here H is symmetric, Q-invariant, and hence `diag(h,h,h_z)`, with `2h+h_z=-tr A^2`. In particular

```
b'     = -2 b^2 - p_zz(0,t)/2 + (nu/2) partial_z Delta u_z(0,t),
Omega' =  2 b Omega + (nu/2)[partial_x Delta u_y-partial_y Delta u_x](0,t). (2.3)
```

No scalar pressure closure is introduced. The quadratic stresses of all generated nonzero angular modes contribute to the pressure mean and thus to `p_zz(0,t)`. Initially `b=Omega=1`, but the chosen amplitude is not the old pressure-neutral value: its certified `b'(0)>2.02484`, `Omega'(0)=2` must be retained. The viscous terms vanish initially; they are not assumed to vanish later.

To display the inherited mean coupling exactly, let `P0` be the rotation average of vector fields, `U=P0 u`, `w=u-U`, `P=P0 p`, and `pi=p-P`. Then

```
U_t+(U.grad)U+grad P = nu Delta U - div P0(w tensor w),
w_t+(U.grad)w+(w.grad)U+grad pi
   = nu Delta w - [(w.grad)w-P0((w.grad)w)].                              (2.4)
```

Every field is divergence free. Consequently

```
-Delta P  = tr[(grad U)^2] + P0 tr[(grad w)^2],
-Delta pi = 2 tr[(grad U)(grad w)] + (I-P0)tr[(grad w)^2].                (2.5)
```

In cylindrical Fourier components, the second term in the first line is
`2 sum_(m>0) Re tr(G_m conjugate(G_m))`; this is a matrix trace, not a squared Frobenius norm. Initially the source has modes 0, ±4, ±8. Later modes are not truncated in (2.4)–(2.5). The mean is solution-adapted and receives the fluctuation stress; it is not prescribed as its initial pump/swirl profile.

For comparison, the exact fluctuation-energy identity is

```
(1/2)||w||_2^2 ' = -integral w^T sym(grad U) w - nu ||grad w||_2^2.       (2.6)
```

The opposite transfer appears in the mean-energy equation. There is no sign restriction on it. These identities identify where an apparently favorable local reference must still pay for the supplying field.

## 3. Exact reference and cancellation

Let F solve the actual matrix equation

```
F'=A(t) F,                         F(0)=I.                              (3.1)
```

For a general A this means a time-ordered propagator. Jacobi's formula gives `det F=1` from `tr A=0`. At the C4-fixed origin only, (2.1) yields the explicit formula

```
B0(t)=integral_0^t b(s)ds,         Theta(t)=integral_0^t Omega(s)ds,
F(t)=diag(exp(-B0),exp(-B0),exp(2B0)) R_Theta.                            (3.2)
```

Define the **actual** remainder and actual flow map by

`R(x,t)=u(x,t)-A(t)x`, `Phi_t(a)'=u(Phi_t(a),t)`, `Phi_0(a)=a`.

Then for `Y(a,t)=F(t)^(-1) Phi_t(a)`, direct differentiation gives

```
Y' = F^(-1) R(FY,t),                         Y(a,0)=a.                 (3.3)
```

The term `A Phi` cancels identically. There is no error between a predicted strain rate and the actual strain rate. In particular there is no sensitivity constant for a fitted pressure/strain functional multiplying `Y-a` merely because of a rate mismatch. This is the exact U1-type cancellation that transfers without an axisymmetric positivity argument.

The cancellation is not a bound on R. Also, `A x` and R separately need not be finite-energy global fields. They are only a local comparison decomposition of the unchanged finite-energy u; this note does not use `||R||_2` on all of space or insert an affine field into the physical datum.

For any time-dependent matrix A, its propagator satisfies the useful estimate

```
||F(t,s)|| <= exp integral_s^t lambda_max(sym A(tau)) d tau.              (3.4)
```

The backward propagator is bounded by the analogous integral of `-lambda_min(sym A)`. This follows by differentiating `|Fv|^2`; no commutation is needed. Thus rigid rotation is not itself an exponential norm cost. It can still rotate strain eigendirections and affect a noncommuting propagator. At the central C4 point, the exact norms in (3.2) depend on B0 and not Theta.

## 4. A cubic-remainder deformation bound

Oddness implies `D^2u(0,t)=0`. Hence on a physical ball `|x|<=R_*`, Taylor's integral remainder gives

```
|R(x,t)| <= M3(t)|x|^3/6,
||grad R(x,t)|| <= M3(t)|x|^2/2,
M3(t)=sup_(|x|<=R_*) ||D^3u(x,t)||.                                    (4.1)
```

Set

```
k3(t)=||F^(-1)(t)|| ||F(t)||^3 M3(t)/6,
I3(t)=integral_0^t k3(s)ds,                  D_r(t)=1-2r^2 I3(t).         (4.2)
```

**Local comparison theorem.** Fix `r>0` and T. Suppose `D_r(t)>0` for `0<=t<=T` and

`||F(t)|| r/sqrt(D_r(t)) < R_*` for `0<=t<=T`.                            (4.3)

Then for every `|a|<=r` and `0<=t<=T`,

```
|Y(a,t)|       <= r D_r(t)^(-1/2),
|Y(a,t)-a|     <= r [D_r(t)^(-1/2)-1],
||D_aY(a,t)-I||<= D_r(t)^(-3/2)-1,
det D_aY(a,t)   =1.                                                     (4.4)
```

**Proof.** Up to the first exit from the physical ball, (3.3) and (4.1) imply `d^+|Y|/dt <= k3 |Y|^3`. The solution of `z'=k3 z^3`, `z(0)=r` is `z=r/sqrt(D_r)`. Integrating `|Y'|<=k3 z^3=z'` gives the displacement bound. Condition (4.3), via a first-exit argument, makes these bounds valid through T.

For `Z=D_aY`, differentiating (3.3) gives

`Z'=F^(-1) (grad R)(FY) F Z`.

Its coefficient has norm at most `3k3 |Y|^2 <= 3k3 r^2/D_r`. Its integral is `-(3/2)log D_r`, so variation of constants gives both `||Z||<=D_r^(-3/2)` and the displayed `||Z-I||` bound. Finally `det Z=det(F^-1)det D_aPhi=1`. This last identity is exact incompressibility, not an estimated Jacobian defect. QED.

This is a new usable inequality, with explicit hypotheses: a material core remains close to the actual linear deformation if the accumulated cubic remainder is small. It uses no growth sign, no pressure ansatz, and no scalar closure. The difference from `exp integral ||grad u||` is substantive for a small core: the remaining exponent is governed by `r^2 I3`, while the actual linear deformation is retained exactly.

**Additional initial cancellation.** The unchanged datum is exactly affine on `|x|<1/2`. If `R_*<1/2` and `sup_(B_R*)||partial_t D^3u||<=H3(t)`, then

`M3(t)<=integral_0^t H3(s) ds`.

Thus, if `||F^-1||||F||^3<=K_F` and `H3<=H_*` on the slab,

`I3(t)<=K_F H_* t^2/12`.                                                (4.5)

The adapted material displacement begins at order `t^2 r^3`, not `t r^3`, on that initial core. More precisely, put

`q0(x)=p(x,0)-p(0,0)-(1/2)x^T H(0)x`.

Because the initial core is affine, its pressure source is constant there; q0 is harmonic, even, Q-invariant, and starts at quartic order. From NS and (2.2),

`R_t(x,0)=-grad q0(x)`, `Y(a,t)=a-(t^2/2)grad q0(a)+O(t^3)`               (4.6)

on compact subsets of the initial affine ball, under the corresponding temporal regularity. This identifies the actual nonlocal pressure distortion left after the central gradient is matched. It does not set it to zero.

### 4a. A version that the existing H4 validation radius can actually supply

An H4 error radius does **not** bound `||D^3u||_infinity` or `||partial_t D^3u||_infinity`. Thus it does not, by itself, price M3 or H3 in (4.1)–(4.5). The following mixed-power estimate avoids that missing derivative.

Let `u_app` be a smooth divergence-free, odd C4-equivariant whole-space reconstruction, and suppose the existing bridge gives

`||e(t)||_D4 <= E4(t)`, `e=u-u_app`,

where `||e||_D4^2=sum_(|alpha|<=4)||partial^alpha e||_2^2`. The reconstruction can be projected onto odd/C4 symmetry if necessary: those orthogonal coordinate/sign actions preserve this derivative norm and the exact solution's symmetry. Every other reconstruction requirement, including divergence, joins and tails, remains.

**Explicit endpoint modulus.** With unitary Fourier convention,

```
[D^2e]_(C^(0,1/2),operator) <= C_H ||e||_D4,
C_H=4 sqrt(3)/pi.                                                       (4.7)
```

For a vector-valued e and unit derivative directions, Fourier inversion and Cauchy–Schwarz bound the square of the difference by

`(2pi)^(-3) ||e||_H4,Fourier^2 integral |xi|^4 min(4,|h|^2|xi|^2)/(1+|xi|^2)^4 dxi`.

The radial integral is at most

`4pi integral_0^infinity r^(-2) min(4,|h|^2 r^2)dr = 16pi |h|`,

splitting at `r=2/|h|`. Thus the constant for the Fourier norm is `sqrt(2)/pi`. The coefficientwise polynomial bound `(1+|xi|^2)^4 <=24 sum_(|alpha|<=4)xi^(2alpha)` gives (4.7). Vector components cause no additional factor: Cauchy–Schwarz uses the Euclidean norm of the vector Fourier transform. This proof retains the endpoint exponent 1/2; replacing `min(2,|h||xi|)` by a single power before integration would lose it.

Since e is odd, `e(0)=0` and `D^2e(0)=0`. Taylor's integral formula consequently gives

```
|e(x)-grad e(0)x| <= (4 C_H/15) E4 |x|^(5/2),
||grad e(x)-grad e(0)|| <= (2 C_H/3) E4 |x|^(3/2).                       (4.8)
```

The factors are `integral_0^1(1-s)s^(1/2)ds=4/15` and `integral_0^1s^(1/2)ds=2/3`. For the smooth reconstruction define the computable quantity

`M3_app(t)=sup_(|x|<=R_*) ||D^3u_app(x,t)||`.

Crucially, decompose the actual remainder as

`R=[u_app-grad u_app(0)x]+[e-grad e(0)x]`.

The actual A in Section 3 is still matched exactly. Therefore

```
|Y'| <= k3_app |Y|^3 + k52 |Y|^(5/2),
k3_app=||F^-1||||F||^3 M3_app/6,
k52=(4 C_H/15) E4 ||F^-1||||F||^(5/2).                                 (4.9)
```

**Explicit mixed-error test.** Choose a label-space buffer `rho_b>r` and put

```
k_mix=k52+sqrt(rho_b) k3_app,
I_mix(t)=integral_0^t k_mix(s)ds,
D_mix(t)=1-(3/2)r^(3/2)I_mix(t),
z_mix(t)=r D_mix(t)^(-2/3).                                              (4.10)
```

If `D_mix>0`, `z_mix<rho_b`, and `||F||z_mix<R_*` throughout the interval, then

```
|Y|<=z_mix,
|Y-a|<=r[D_mix^(-2/3)-1],
||D_aY-I||<=D_mix^(-2)-1,                   det D_aY=1.                 (4.11)
```

Indeed, before exit `|Y|^3<=sqrt(rho_b)|Y|^(5/2)`, and the comparison equation is `z'=k_mix z^(5/2)`. The variational coefficient is bounded by

`3k3_app|Y|^2+(5/2)k52|Y|^(3/2) <=3k_mix z_mix^(3/2)`.

Its integral is `-2log D_mix`, yielding (4.11); the two strict buffer conditions justify continuation by first exit. This connects the **existing H4 error norm** to an actual local core deformation estimate without silently upgrading it to H5. It supplies an inequality to evaluate, not a successful numerical margin.

**Enclosing F without knowing the true central histories exactly.** Let `b_app,Omega_app` be extracted from the odd C4 reconstruction and let certified central errors give `|b-b_app|<=epsilon_b`, `|Omega-Omega_app|<=epsilon_Omega`. Put `epsilon_B=integral epsilon_b`, `epsilon_Theta=integral epsilon_Omega`. The two-parameter central matrices commute even between the exact and approximate fields. Hence

```
F=F_app diag(exp(-Delta B),exp(-Delta B),exp(2Delta B)) R_(Delta Theta),
||F||<=exp(2epsilon_B)||F_app||,
||F^-1||<=exp(2epsilon_B)||F_app^-1||.                                  (4.12)
```

For phase-sensitive comparison to the reconstructed linear reference,

```
||F-F_app|| <= ||F_app||[exp(2epsilon_B)-1
                       +exp(2epsilon_B) min(2,epsilon_Theta)].          (4.13)
```

Thus both the actual deformation coefficients and the residual phase error can be charged using the H4 bridge's central-gradient enclosure. Rotation errors must not be aligned away. A physical comparison follows from
`|Phi-F_app a|<=||F|| |Y-a|+||F-F_app||r`.

## 5. The transformed velocity does not obey unchanged NS

The trajectory cancellation must not be mistaken for cancellation of polarization, pressure or viscosity in the field equation. Define

```
v(y,t)=F^(-1)u(Fy,t)-B(t)y,       B=F^(-1)AF=F^(-1)F',
C=F^(-1)F^(-T),                  L_C=sum C_ij partial_i partial_j,
q(y,t)=p(Fy,t)-p(0,t)-(1/2)y^T F^T H F y.                               (5.1)
```

Here `v=F^-1 R(Fy)` is the relative velocity; trajectories Y obey `Y'=v(Y,t)`. `div_y v=0`, `v(0)=grad_y v(0)=0`, and v is odd. Put

`K(t)=grad_y(L_C v)(0,t)=F^(-1) D(t) F`.

The **exact** local transformed equations are

```
v_t+(v.grad_y)v+2B v = -C grad_y q + nu [L_C v-K(t)y],
-L_C q = 2 tr[B grad_y v] + tr[(grad_y v)^2].                            (5.2)
```

To check the first line, a trajectory `X=FY` has acceleration
`X''=F''Y+2F'v+F(v_t+v.grad v)`. Also `Delta_x u=F L_Cv`, and the quadratic pressure contributes `-F^-1 H F y`. Use `F''=(A'+A^2)F=(-H+nu D)F` to cancel the linear pressure and retain the viscous linear jet `-nu K y`. The second line follows by subtracting the constant pressure source `tr A^2` and using similarity invariance of trace.

In particular:

- The `2Bv` coupling survives. Matching material deformation does not freeze velocity polarization.
- Diffusion becomes anisotropic with metric C; it is not `nu Delta_y` except for an orthogonal F.
- The whole nonlocal q and its inherited harmonic content survive. Local Poisson data alone do not determine it.
- Removing the linear viscous jet is essential. Applying the Laplacian to a cubic remainder produces a linear term, which (2.2) has already assigned to A.

This transformed system is an identity for the actual field, not a separate finite-energy reference solution. Boundary/tail data for q and the actual matrices A,D must still come from the full solution.

## 6. Why the outer receiver needs full matrices and time ordering

C4 symmetry constrains a jet only at points fixed by Q. At an off-axis point it relates the four rotated jets by conjugation; it does not make an individual gradient a scalar strain plus rotation.

At the initial mean-angular-momentum maximum `(r_*,z=4)`, the actual mean field on the relevant pump plateau is

`U_r=c r`, `U_z=-2c z`, `U_theta=1020 r chi(r^2/a^2)`, `c=7/5`, `a=63/200`.

The maximum condition is `partial_r(r U_theta)=0`, so `partial_r U_theta=-U_theta/r`. Writing `omega_s=U_theta/r`, the initial **mean** gradient in its cylindrical orthonormal basis is

```
grad U = [[c,-omega_s,0],[-omega_s,c,0],[0,0,-2c]].                       (6.1)
```

Its horizontal eigenvalues are `c±omega_s`. This is a genuine large symmetric strain. Removing a rigid angular phase does not remove it. Along the mean trajectory, the basis itself rotates at `omega_s`; at that instant its frame generator is

```
(grad U)_cyl-omega_s J = [[c,0,0],[-2omega_s,c,0],[0,0,-2c]].             (6.2)
```

This explains the shear form in the earlier Kelvin diagnostic. Freezing (6.2) later would be a model assumption: the actual mean evolves by (2.4), including seed stress, axial pressure and viscosity. The actual full gradient also contains the phase-dependent seed gradient. Both the mean and full off-axis gradients must retain their time ordering.

There is an exact off-axis version of Section 3. Choose one actual material center `X_c'=u(X_c,t)`, define `A_c=grad u(X_c,t)`, `F_c'=A_c F_c`, and

`R_c(x,t)=u(X_c+x,t)-u(X_c,t)-A_c x`.

For particles relative to X_c, `Y'=F_c^-1 R_c(F_cY,t)`. With `M2=sup_tube ||D^2u||`,

```
k2=||F_c^-1||||F_c||^2 M2/2,             I2=integral_0^t k2,
|Y|<=r/(1-rI2),
|Y-a|<=r[(1-rI2)^(-1)-1],
||D_aY-I||<=(1-rI2)^(-2)-1,             det D_aY=1,                     (6.3)
```

provided `rI2<1` and the predicted physical tube remains inside the region where M2 is bounded. The proof is the same scalar comparison with a quadratic remainder. Oddness does not eliminate `D^2u` at an off-axis receiver. Nor is the radial maximizer of the evolving angular mean automatically a material center: a chosen Eulerian receiver has to be compared separately.

## 7. A solution-adapted mean reference, if the whole packet is needed

An alternative retains the full actual axisymmetric mean U rather than only a central affine jet. Let `Psi_t` be its divergence-free flow and define `Y=Psi_t^-1(Phi_t(a))`. Then exactly

`Y'=[D Psi_t(Y)]^-1 w(Psi_t(Y),t)`.                                    (7.1)

The mean advection cancels, but U still obeys (2.4). The frame deformation is bounded by symmetric strain, as in (3.4); its metric and spatial derivatives are not optional. This identity captures actual radial/axial drift and differential rotation, unlike a single rigid phase factor.

For an Eulerian formulation that keeps pressure explicit, define the time-dependent divergence-free linear operator

`L_U(t)z=nu Delta z-P_Leray[(U.grad)z+(z.grad)U]`.

On the zero-angular-mean subspace its actual propagator S_U retains full pressure/polarization, and

```
w(t)=S_U(t,0)w(0)-integral_0^t S_U(t,s) P_Leray(I-P0)[(w.grad)w](s) ds.   (7.2)
```

Its L2 propagator bound is

`||S_U(t,s)||_(2->2)<=exp integral_s^t sup_x lambda_max[-sym grad U]`.     (7.3)

The proof is the same integration by parts as (2.6). The time-ordered operator and the actual stress-driven U are essential. Formula (7.2) is not an independent approximate solution until the nonlinear remainder and the dependence of U on w are controlled. Its advantage is precise: the dominant actual mean evolution is carried in the propagator, rather than bounded again as a model-rate mismatch. Its disadvantage is also precise: obtaining S_U and its input mean still requires the full coupled trajectory.

## 8. Required bounds and an eligible next lemma

For the core comparison, a meaningful finite-slab certificate would supply:

1. Enclosures of the whole-space central history A(t), for example through (4.12) and the H4 error radius. Its dynamics remain the full identity (2.2), with pressure and viscous terms inherited from the validated NS field. Separate pointwise enclosures of H(t) and D(t) are **not** needed for the geometric bound (4.11); if (2.3) or the transformed field equation is used quantitatively, those terms need their own justified bounds. In particular H4 does not by itself give a pointwise bound on D=Delta grad u.
2. A physical ball/tube and either certified actual `M3(t)`/`H3(t)`, or the H4-compatible pair `M3_app(t),E4(t)` from Section 4a, plus bounds on F and F^-1. An H4 radius alone does not enclose actual M3/H3.
3. Strict inequalities (4.3) or the mixed-error buffer test (4.10)–(4.11), and a displacement/gradient tolerance that is actually useful for a named core observable.
4. Separate bounds on the inherited outer field, full receiving RMS, Fourier phases, and stresses in (2.4). The local theorem supplies none of these automatically.

The bounded next lemma is therefore concrete: **enclose the core remainder and establish (4.4), or use the existing H4 radius with the reconstruction's third derivatives to establish (4.11), with a prescribed deformation tolerance and actual central histories coupled to the whole-space approximation.** The mixed-error version states precisely how that reconstruction/error data can price the test; a list of favorable central jets cannot.

For an outer-packet reference, use (6.3) or (7.1)–(7.3), with full matrix/mean evolution and actual second spatial derivatives. A scalar replacement `A(t)=a(t)S` is excluded by (6.1) even at the initial mean maximum. The central two-coefficient matrix symmetry is valid, but its dynamics are not closed in b and Omega because of pressure, viscosity and generated higher jets.

No inequality here supplies retained endpoint gain, a stable successor class, or a renewed reservoir. It gives an exact cancellation, explicit local deformation tests, and the correct transformed equations to use when checking those further claims.

## 9. Verification

`check_c4_adapted_reference.py` passed 37 exact checks, recorded in `c4-adapted-reference-check.json`. These check the C4 commutant, allowed quartic pressure modes, full matrix acceleration with the factor two, anisotropic viscosity and pressure transformations for a manufactured noncommuting F, the cubic/quadratic/mixed comparison exponents, and the explicit Fourier/Taylor constants. The manufactured polynomial field is a local algebra test, not a finite-energy NS solution. No fluid norms, trajectory duration, or endpoint margin were supplied by this checker.

# GAP T — the scale-invariant Lipschitz stability of the axis-strain functional

Seat `gaps/gap-T-lipschitz`, DTC-2026-09-06.  Pre-registration in `PREREG.md`, written before any
script ran.  Every number below came out of a script in this folder that I wrote and ran
(`s1`…`s5`, re-asserted by `check_constants.py`); `SHA256SUMS` is computed, never typed.
Nothing outside this folder was written.  The algebra is exact (sympy); the numerics are refuters —
they falsify, they never prove.

---

## 0. Headline

The estimate GAP T asks for is **PROVED**, with explicit constants, in the form

> `|a[eta_0 o Phi^-1](0) - a[eta_0 o T_lambda^-1](0)|  <=  pi lambda M L [ 3 mu (1+mu_J)/(2(1-mu)^5) + 3 mu_J/8 ]`

and its scale-invariant corollary, for the bang-bang cap, is

> `|Delta| / a[eta_0 o T_lambda^-1](0)  <=  (15 pi/4) mu (1 + O(mu))  =  11.7809725 mu` ,

free of `lambda`, of `L`, and of `rho_0/R` — exactly the form `lower/prove-lagrangian` §4(2) needs
(`mu = O(c/L)` ⟹ relative error `O(c/L)`).

Three things the route as briefed did **not** anticipate, all established here:

1. **The displacement must be measured against the IMAGE point `|T_lambda x|`, not against `|x|`.**
   The brief's normalisation `|Phi-T_lambda|/|y| <= mu` costs a factor `max(lambda^-1, lambda^2)`
   because `T_lambda` is anisotropic (`|T_lambda x|/|x|` ranges over `[lambda^-2, lambda]`).
   Corollary 2 records the exact conversion; it is harmless at `lambda <= 3/2` (factor `2.25`).
2. **A Jacobian hypothesis is NOT optional and NOT a technicality.**  The registered control fires:
   at fixed relative displacement `mu` there are genuine `SO(4)`-equivariant `C^1` diffeomorphisms
   for which `|Delta|` is **unbounded**, growing like `phi_c^-1` as the axis-layer width `phi_c -> 0`.
   No bound of the form `C mu lambda M L` can hold on hypothesis (H2) alone.
3. **The whole radial direction is a null direction.**  For `Psi(u) = s(|u|) u` the `s^4` in the
   Jacobian cancels the `s^-4` in the kernel identically, and if `s` is log-periodic with period
   `L/k` the strain is *exactly unchanged* — for every `lambda`, every datum, and every amplitude.
   Radial rearrangement of the shells cannot move `a(0)` at all.

The integral the brief asked me to CHECK does converge, and better than asked: it is **exact**,
`int_S |eta_0| |T_lambda x|^-4 dx_5 = pi^3 M L / lambda` for the plateau (`<=` in general).

---

## 1. Setting and notation

5-D lift of axisymmetric no-swirl NS.  `x = (y,z)`, `y in R^4`, `r = |y|`, `rho = |x|`,
`phi` the polar angle from `+z`.  `eta = omega^theta/r`, `-Delta_5 psi_1 = eta`, `a = u^r/r`, and
(`rebuild/far-near-kernel-lemma` §0)

```
a[f](0) = int_{R^5} K(w) f(w) dw ,      K(w) = -(3/(8 pi^2)) w_z/|w|^5 ,   homogeneous of degree -4.
S = {rho_0 < |x| < R} ,   L = log(R/rho_0) ,   T_lambda(y,z) = (lambda y, lambda^-2 z) ,  J_{T_lambda} = lambda^2 .
dx_5 = 2 pi^2 rho^4 sin^3(phi) drho dphi .
```

**Datum (minimal hypothesis).**  `eta_0` measurable, supported in `S`, with

```
(D)   |eta_0(x)| <= M / r(x)      i.e.   |omega_0^theta| <= M  on S.
```

That is **all** that is used.  Oddness in `z` is *not* needed; bounded variation in the angle is
*not* needed; no regularity of `eta_0` beyond measurability is needed.  (Oddness and the taper do
enter elsewhere — see §6 — but not in Lemma T.)

**Map hypotheses.**  `Phi : S -> R^5` with

```
(H1)  Phi is a C^1 diffeomorphism onto its image, J_Phi := det D Phi > 0 ;
(H2)  |Phi(x) - T_lambda(x)| <= mu |T_lambda(x)|  for all x in S,   0 <= mu < 1 ;
(H3)  |J_Phi(x) - lambda^2| <= mu_J lambda^2      for all x in S,   mu_J >= 0 .
```

(H2) is the scale-invariant metric `|delta X|/|X|` of `prove-lagrangian` §4(2), with `X` the
**deformed** position.  (H3) is that seat's own bullet "the 5-D Jacobian and the log-radial measure
are preserved to `1 + O(c/L)`".

---

## 2. LEMMA T (PROVED)

> **LEMMA T.**  Under (D), (H1)–(H3), both `a[eta_0 o Phi^-1](0)` and `a[eta_0 o T_lambda^-1](0)`
> converge absolutely, and
>
> ```
> | a[eta_0 o Phi^-1](0) - a[eta_0 o T_lambda^-1](0) |
>          <=  pi lambda M L [ 3 mu (1 + mu_J) / (2 (1-mu)^5)  +  3 mu_J / 8 ] .        (T)
> ```
>
> In particular, if `mu_J <= mu <= mu_*` then `|Delta| <= C(mu_*) mu lambda M L` with
> `C(mu_*) = pi[ 3(1+mu_*)/(2(1-mu_*)^5) + 3/8 ]`, and `C(0+) = 15 pi/8 = 5.8904862`.

| `mu_*` | 0 | 0.01 | 0.05 | 0.10 | 0.20 | 0.25 | 0.50 |
|---|---|---|---|---|---|---|---|
| `C(mu_*)` | **5.8905** | 6.1829 | 7.5727 | 9.9566 | 18.4354 | 26.0006 | 227.3728 |

### Proof, step by step

**Step 1 (change of variables). PROVED.**  By (H1) and `eta_t = eta_0 o Phi^-1`,
`a[eta_0 o Phi^-1](0) = int_{Phi(S)} K(w) eta_0(Phi^-1 w) dw = int_S K(Phi(x)) eta_0(x) J_Phi(x) dx`,
and likewise for `T_lambda`.  Hence

```
Delta = int_S eta_0(x) [ (K(Phi(x)) - K(T_lambda x)) J_Phi(x) + K(T_lambda x)(J_Phi(x) - lambda^2) ] dx.
```

**Step 2 (kernel bounds, sharp). PROVED** (`s1_kernel_constants.py`, sympy, residual `0`).
`|K(w)| = (3/(8 pi^2)) |cos phi_w| |w|^-4 <= (3/(8 pi^2)) |w|^-4`, sharp on the axis.
`grad K(w) = -(3/(8 pi^2)) |w|^-5 ( e_z - 5 (w_z/|w|) w/|w| )` (verified componentwise in `R^5`),
and `|e_z - 5 c what|^2 = 1 + 15 c^2 <= 16`, so

```
sup_{|w|=1} |grad K| = 4 * 3/(8 pi^2) = 3/(2 pi^2) = 0.15198177546350666 ,
|grad K(w)| <= (3/(2 pi^2)) |w|^-5   (sharp, attained on the axis).
```

**Step 3 (the shell integral — EXACT, this is the step the brief asked me to check). PROVED.**
Put `I_lambda := int_S |eta_0(x)| |T_lambda x|^-4 dx_5`.  Substituting `u = T_lambda x`
(`dx = lambda^-2 du`, `r_x = lambda^-1 r_u`) and using (D),

```
I_lambda <= (M/lambda) int_{T_lambda S} du_5 / ( r_u |u|^4 ) .
```
In `(rho_u, phi)` the integrand is `2 pi^2 sin^2(phi) / rho_u`, and `T_lambda S` is
`{ rho_0/g(phi) < rho_u < R/g(phi) }` with `g(phi) = |T_lambda^-1 u|/|u|`, so the **radial log-length
is exactly `L` for every angle** — the anisotropy cancels.  Hence, exactly,

```
int_{T_lambda S} du_5/(r_u |u|^4) = 2 pi^2 L int_0^pi sin^2(phi) dphi = pi^3 L ,
I_lambda <= pi^3 M L / lambda ,   with EQUALITY iff |omega_0^theta| = M a.e. on S.     (3.1)
```
So the answer to "CHECK whether this converges to `C M log(R/rho_0)`" is **yes, with `C = pi^3` and
no `lambda`-loss**: the `r^-1` weight is integrable against the `r^3 dr dz` measure
(`int_0^pi sin^2 = pi/2`), and `|x|^-4` against `|x|^4 d|x|` is exactly the log.
Independent 2-D Gauss–Legendre recomputation of (3.1) at `lambda = 0.5, 1, 1.25, 1.5, 3`:
relative error `<= 1.54e-13`.  **KILL K1 did not fire.**

**Step 4 (the kernel difference). PROVED.**  By (H2) every point `w` of the segment
`[T_lambda x, Phi(x)]` obeys `|w| >= (1-mu)|T_lambda x| > 0`, so the segment misses the singularity
and the fundamental theorem of calculus along it gives

```
|K(Phi(x)) - K(T_lambda x)| <= mu |T_lambda x| * (3/(2 pi^2)) ((1-mu)|T_lambda x|)^-5
                            = (3 mu / (2 pi^2 (1-mu)^5)) |T_lambda x|^-4 .
```
With `J_Phi <= (1+mu_J) lambda^2` from (H3) and (3.1),
`|term I| <= (3 mu (1+mu_J) lambda^2 / (2 pi^2 (1-mu)^5)) I_lambda <= (3 pi mu (1+mu_J)/(2(1-mu)^5)) lambda M L`.

**Step 5 (the Jacobian term). PROVED.**  `|K(T_lambda x)| <= (3/(8 pi^2))|T_lambda x|^-4` and
`|J_Phi - lambda^2| <= mu_J lambda^2`, so `|term II| <= (3 mu_J lambda^2/(8 pi^2)) I_lambda
<= (3 pi mu_J/8) lambda M L`.

**Step 6 (assembly). PROVED.**  Add.  Absolute convergence of both functionals is the same
computation with `|K| <= (3/(8 pi^2))|w|^-4`: `|a[eta_0 o T_lambda^-1](0)| <= (3 pi/8) lambda M L`
and `|a[eta_0 o Phi^-1](0)| <= (3 pi (1+mu_J)/(8(1-mu)^4)) lambda M L`.  ∎

### Corollary 1 (the scale-invariant form the Lagrangian route consumes). PROVED.
For the bang-bang cap `omega_0^theta = -M sgn(z)`, Lemma 1 of `prove-lagrangian` gives
`a[eta_0 o T_lambda^-1](0) = (M/2) lambda L` exactly (re-verified here to `3.2e-14` at five `lambda`,
`s3`), so

```
|Delta| / a[eta_0 o T_lambda^-1](0)  <=  2 pi [ 3 mu(1+mu_J)/(2(1-mu)^5) + 3 mu_J/8 ]
                                     ->  (15 pi/4) mu = 11.7809725 mu     as mu, mu_J -> 0,
```
**independent of `lambda`, of `L` and of the geometry.**  For a general datum divide (T) by the
actual `a[eta_0 o T_lambda^-1](0)`.

### Corollary 2 (the brief's literal normalisation). PROVED.
If instead `|Phi(x) - T_lambda(x)| <= mu |x|`, then since
`|T_lambda x|/|x| in [min(lambda, lambda^-2), max(lambda, lambda^-2)]`, hypothesis (H2) holds with

```
mu~ = mu * max(lambda^-1, lambda^2) ,   and (T) applies with mu -> mu~ (needs mu~ < 1).
```
At `lambda <= 3/2` (the whole range `prove-lagrangian` uses) the loss is at most `lambda^2 = 2.25`.
This is a **repair of the task statement**, not a stylistic choice: the anisotropy of `T_lambda`
is exactly `lambda^3` between axis and equator, and it is real.

---

## 3. Two exact identities found on the way (ADDITIVE)

**(N) The radial null direction. PROVED (sympy, `s5_null_direction.py`).**
For `Psi(u) = s(|u|) u` composed after `T_lambda`, the integrand factor is exactly

```
(Z_Phi/|Phi|^5) J_Phi = z_u lambda^2 ( s + varrho s'(varrho) ) / varrho^5 ,     varrho = |u| ,
```
— the `s^4` of the Jacobian cancels the `s^-4` of the degree `-4` kernel *identically* (sympy
residual `0`).  Writing `h(varrho) = varrho (s-1)`, `s + varrho s' - 1 = h'`, and integrating by
parts in `d varrho/varrho`,

```
a_Phi(0) - a_{T_lambda}(0) = -(3/4) int_0^pi (cos phi sin^2 phi / g(phi)^5) omega_0^theta(phi)
                              [ s(R g) - s(rho_0 g) + int_{rho_0 g}^{R g} (s-1) d varrho/varrho ] dphi .
```
Hence **if `s` is log-periodic with period `L/k` and zero log-mean over a period, `Delta = 0`
exactly** — for every `lambda`, every datum, every amplitude, however large.  Verified to
`3.2e-15` relative at `k = 2, 4`, `mu = 0.05 … 0.5`, `lambda = 1, 1.5`.
*Consequence for the campaign:* the `O(c/L)` **shear between shells** in `prove-lagrangian` §4(2)
is a radial rearrangement — the least dangerous of the three error terms listed there.  All of the
Lipschitz sensitivity of `a(0)` lives in the angular / Jacobian part of the perturbation.

**(A) The exact half-period benchmark. PROVED (closed form + quadrature).**
For `s(varrho) = 1 + mu sin(pi log(varrho/rho_0)/L)` at `lambda = 1`, the bracket above is
`2 mu L/pi`, so with `kappa := -(3/4) int cos phi sin^2 phi omega_0^theta dphi / M`,

```
Delta = (2 kappa / pi) mu M L      exactly, for every mu (while the map is a diffeomorphism).
```
Measured `|Delta|/(mu M L) = 0.3183098861837907 = 1/pi` at `mu = 0.05, 0.2, 0.5, 0.9, 0.93` (and,
as the formula predicts, unchanged past the diffeomorphism limit `mu = 0.9342`), and
`0.1323076` vs the predicted `0.1323076` for the `7.5°` taper.  This is a **sharpness benchmark**:
in the radial direction the true Lipschitz constant is exactly `2 kappa/pi` (`= 1/pi` for the cap),
against the proved `C(0+) = 15 pi/8`: **Lemma T is conservative by a factor `18.5` there.**

---

## 4. Numerical verification of (T) — two map families x two data x three `(lambda,mu)`

`rho_0 = 1`, `R = 4096`, `L = 8.317766166719343`, `M = 1`.
Data: **D1** bang-bang cap `omega^theta = -M sgn(z)`; **D2** admissible `7.5°` axis taper.
Maps (both composed after `T_lambda`, both genuine `SO(4)`-equivariant `C^1` diffeomorphisms):
**A** radial half-period ripple `u -> u(1 + mu sin(pi log(|u|/rho_0)/L))`;
**B** half-plane angular shear `phi_u -> phi_u + mu_b sin^2(phi_u)`, `mu_b` set by `2 sin(mu_b/2) = mu`.
`mu` and `mu_J` are **measured** by a dense `2001 x 257` scan, not assumed.
Own Gauss–Legendre quadrature, panels split at `phi = pi/2` and at the taper corners; `n` vs `2n`.

| datum | map | `lambda` | `mu` | `mu_J` | `a_T` | `a_Phi` | `|Delta|` | bound (T) | slack | doubling |
|---|---|---|---|---|---|---|---|---|---|---|
| D1 | A | 1.00 | 0.05000 | 0.27715 | 4.158883 | 4.291264 | 0.132381 | 5.9505 | 44.95 | 3.3e-15 |
| D1 | B | 1.00 | 0.05000 | 0.12972 | 4.158883 | 4.143762 | 0.015121 | 4.1325 | 273.30 | 4.1e-15 |
| D1 | A | 1.25 | 0.10000 | 0.61259 | 5.198604 | 5.534224 | 0.335620 | 20.8840 | 62.23 | 3.9e-15 |
| D1 | B | 1.25 | 0.10000 | 0.26903 | 5.198604 | 5.123438 | 0.075166 | 13.8250 | 183.93 | 4.0e-15 |
| D1 | A | 1.50 | 0.20000 | 1.49421 | 6.238325 | 7.058219 | 0.819895 | 111.4690 | 135.96 | 4.8e-15 |
| D1 | B | 1.50 | 0.20000 | 0.57725 | 6.238325 | 5.885843 | 0.352481 | 65.0854 | 184.65 | 4.5e-15 |
| D2 | A | 1.00 | 0.05000 | 0.27715 | 4.156564 | 4.288872 | 0.132308 | 5.9505 | 44.98 | 3.5e-15 |
| D2 | B | 1.00 | 0.05000 | 0.12972 | 4.156564 | 4.141443 | 0.015121 | 4.1325 | 273.29 | 3.9e-15 |
| D2 | A | 1.25 | 0.10000 | 0.61259 | 5.178001 | 5.512245 | 0.334245 | 20.8840 | 62.48 | 2.7e-15 |
| D2 | B | 1.25 | 0.10000 | 0.26903 | 5.178001 | 5.102791 | 0.075209 | 13.8250 | 183.82 | 3.0e-15 |
| D2 | A | 1.50 | 0.20000 | 1.49421 | 6.124926 | 6.929422 | 0.804496 | 111.4690 | 138.56 | 1.3e-16 |
| D2 | B | 1.50 | 0.20000 | 0.57725 | 6.124926 | 5.770126 | 0.354801 | 65.0854 | 183.44 | 4.6e-15 |

**KILL K3 did not fire** (the bound holds in all 12 cases; worst slack `44.95`).
**KILL K4 did not fire** (worst doubling residual `4.8e-15`, threshold `1e-6`).
**KILL K5 did not fire**: `a[eta_0 o T_lambda^-1](0) = (M/2) lambda L` to `3.2e-14` at
`lambda = 0.7, 1, 1.25, 1.5, 2`.

Cross-check of the instrument against the two neighbouring seats: the `lambda = 1` value for D2
gives `kappa_delta = a_T/(ML) = 0.4997212305`, reproducing `prove-lagrangian` §2 (`0.4997212`) and
`far-near-kernel-lemma` §6 (`0.499721`) — a third, independent evaluation of the same constant.

---

## 5. THE CONTROL — where the estimate must fail, and does (KILL K6 FIRED)

Registered control C1: the angular shear with an **axis ramp of width `phi_c`**,

```
beta(phi) = mu_b * m(phi/phi_c) * m((pi-phi)/phi_c) * cos(phi) ,   m(t) = t^2(3-2t) on [0,1], 1 after.
```
This is a genuine `C^1` `SO(4)`-equivariant diffeomorphism (`1 + beta' >= 1 - mu_b > 0`), and its
image-relative displacement is `sup 2|sin(beta/2)| = 2 sin(mu_b/2) = mu`, **independent of `phi_c`**.
Only `mu_J` moves: `mu_J ~ (sin mu_b/phi_c)^3 -> infinity`.  On the bang-bang cap D1, `lambda = 1`:

| `mu_b` | `mu` | `phi_c` | `|Delta|` | `C_0 mu M L` | ratio `|Delta|/(mu M L)` | `mu_J` |
|---|---|---|---|---|---|---|
| 0.9 | 0.8699 | 1e-1 | 24.94 | 42.62 | 3.45 | 1.4e+04 |
| 0.9 | 0.8699 | 1e-2 | 191.95 | 42.62 | **26.53** | 1.0e+08 |
| 0.9 | 0.8699 | 1e-4 | 1.765e+04 | 42.62 | **2439.1** | 1.0e+16 |
| 0.9 | 0.8699 | 1e-8 | 1.762e+08 | 42.62 | **2.43e+07** | 1.0e+32 |
| 0.9 | 0.8699 | 1e-12 | 1.762e+12 | 42.62 | **2.43e+11** | 1.0e+48 |
| 0.5 | 0.4948 | 1e-2 | 31.08 | 24.24 | **7.55** | 1.2e+07 |
| 0.5 | 0.4948 | 1e-12 | 2.383e+11 | 24.24 | **5.79e+10** | 1.1e+47 |
| 0.3 | 0.2989 | 1e-2 | 6.16 | 14.64 | 2.48 | 1.6e+06 |
| 0.3 | 0.2989 | 1e-3 | 37.56 | 14.64 | **15.11** | 1.5e+10 |
| 0.3 | 0.2989 | 1e-12 | 3.409e+10 | 14.64 | **1.37e+10** | 1.5e+46 |

`C_0 = 15 pi/8 = 5.890486`.  Bold rows violate the linear corollary.  Fitted growth law over
`phi_c in [1e-12, 1e-2]` (own least squares on `log10`):

```
|Delta| = c(mu_b) phi_c^-p ,   p = 0.9979 (mu_b=0.9), 0.9935 (0.5), 0.9855 (0.3) ;
c = 1.83896, 0.27202, 0.04587 ;   c ~ mu_b^3.357 .
```
So `|Delta|` is **unbounded at fixed `mu`**: on hypotheses (H1)+(H2) alone no bound
`|Delta| <= C(mu) lambda M L` exists.  **(H3) is necessary.**  Lemma T itself is *not* violated —
its bound, evaluated at the measured `mu_J`, exceeds `|Delta|` in every row by `>= 5 orders`
(`9.18e+21` vs `1.76e+04` at `phi_c = 1e-4`) — which is exactly the point: the `mu_J` in (T) is
load-bearing, not decorative.

**Erratum vs PREREG.**  PREREG predicted the control's rate as `log(1/phi_c)` (from the
`beta` plateau `int d phi/phi`).  That is **wrong**: the measured rate is `phi_c^-1`, and it comes
from the *transition layer* `phi < phi_c`, where `beta' ~ mu_b/phi_c` and
`int_0^{phi_c} (beta^3/phi) beta' dphi ~ mu_b^4/phi_c`.  The registered *direction* (K6 must fire)
was right; the registered *rate* was wrong and is corrected here.

**Second erratum vs PREREG.**  PREREG's control C2 used the *full-log-period* radial ripple at
`mu >= 1`.  §3(N) shows that family is an exact **null direction** — `Delta = 0` identically — so
C2 as registered could never fire.  It was re-run on the non-null half-period ripple: there the
ratio is exactly `1/pi` for every `mu` up to and past the diffeomorphism limit `mu = 0.9342`, so
**pure radial perturbation never breaks the estimate**, at any amplitude.  Taken with C1, this
localises the entire failure mode of GAP T in the angular/Jacobian direction.

**The taper is immune (additive find).**  The same control on the admissible `7.5°` axis taper D2
**saturates** instead of diverging: `|Delta| = 13.19837, 13.19882, 13.19939, 13.19902, 13.20003,
13.20045` at `phi_c = 1e-2 … 1e-12` (spread `2.1e-3` over ten decades).  Admissibility
(`omega^theta -> 0` at the axis like `phi/delta`) kills the `1/phi` near-axis weight that drives the
divergence — the same mechanism as the `1/sin(phi)` collar constant of `far-near-kernel-lemma` §3(b)
and its `c_edge ~ (1/4) log(1/delta)` ceiling.  Beyond `phi_c ~ 1e-12` the evaluation is not
trustworthy in double precision (`J_5` reaches `1e+48`) and those rows are not reported.

---

## 6. Status per step, and what this does and does not close

| step | statement | status |
|---|---|---|
| S1 | change of variables; both functionals absolutely convergent | **PROVED** |
| S2 | `|K| <= 3/(8pi^2)|w|^-4`, `|grad K| <= 3/(2pi^2)|w|^-5`, both sharp | **PROVED** (sympy, residual 0) |
| S3 | `int_S |eta_0||T_lambda x|^-4 dx_5 = pi^3 M L/lambda` (`<=` in general) | **PROVED** (exact; quadrature `1.5e-13`) |
| S4 | segment misses the origin; mean-value bound | **PROVED** |
| S5 | Jacobian term | **PROVED** |
| **T** | **the estimate (T), with explicit constants** | **PROVED** |
| C1 | scale-invariant relative form, `11.7809725 mu` | **PROVED** (given Lemma 1 of `prove-lagrangian`) |
| C2 | conversion from the `|x|`-normalisation, factor `max(lambda^-1, lambda^2)` | **PROVED** |
| N | radial null direction | **PROVED** (sympy, residual 0) |
| A | `Delta = (2 kappa/pi) mu M L` for the half-period ripple | **PROVED** (closed form; quadrature `3e-6`) |
| — | necessity of (H3) | **ESTABLISHED BY COUNTEREXAMPLE** (numerical, `phi_c^-1`, 10 decades) |

**What this closes in GAP T.**  `lower/prove-lagrangian` §4(2) named two things: (a) that the flow
map on `[0,tau]` is `T_{lambda(s)}` up to relative `O(c/L)`, and (b) that the functional (1.1) is
stable under such a perturbation.  **(b) is now proved**, in the exact form needed, with the
constant `11.78` on the relative error, and with the hypothesis list made explicit.

**What it does not close.**  (a) is a PDE statement and is untouched here.  Lemma T is purely a
statement about the kernel and the two maps; it assumes (H2) and (H3) and proves nothing about the
true Navier–Stokes flow map.  §5 shows that assuming (H2) alone would be *false*, so the seat that
proves (a) must deliver **both** bullets — the displacement bound *and* the Jacobian bound — and,
by Corollary 2, must deliver the displacement bound **relative to the deformed position**, not to
the label.  With `prove-lagrangian`'s own §4(2) figures (`|theta| <= C_3 M rho`, `C_3 <= 0.71`,
`tau = c/(ML)`, `rho partial_rho lambda/lambda = -c/(2L)`) both hypotheses hold with
`mu, mu_J = O(c/L)`, giving relative error `<= 11.78 * O(c/L)` — i.e. `O(1/L)`, which is the order
`prove-lagrangian` already carries.  So `T` is downgraded from **GAP** to **PROVED-MODULO L3v**.

**Sharpness.**  The proved constant is conservative: measured `bound/|Delta|` is `45 … 273` on the
twelve verification rows, and the exact radial constant is `2 kappa/pi` against the proved
`15 pi/8` — a factor `18.5`.  Sharpening is not worth doing unless a downstream argument prices the
constant numerically; the structure (`mu` linearly, `mu_J` linearly, `lambda M L` scaling, and the
`(1-mu)^-5` blow-up as the map approaches the origin) is what the Lagrangian route consumes.

---

## 7. Files

| file | what it establishes |
|---|---|
| `PREREG.md` | pre-registration: claim, kill criteria K1–K6, controls, fixed settings |
| `lib_quad.py` | the shared 2-D quadrature for `a[eta_0 o Phi^-1](0)`, the map classes, the bound (T) |
| `s1_kernel_constants.py` / `s1_results.json` | sympy: `grad K`, the sharp constants `3/(8pi^2)`, `3/(2pi^2)`, `pi^3 L`, `P(lambda) = lambda` |
| `s2_lemma_constants.py` / `s2_results.json` | `C_0 = 15pi/8`, the `C(mu_*)` table, the term-I/term-II split |
| `s3_numeric_check.py` / `s3_results.json` | the 12-row verification of (T); K1, K3, K4, K5; `kappa_delta` cross-check |
| `s4_controls.py` / `s4_results.json` | control C1 (K6 FIRES), the `phi_c^-1` growth law, the taper immunity, control C2 |
| `s5_null_direction.py` / `s5_results.json` | the exact radial null direction (sympy residual 0 + machine-zero numerics) |
| `check_constants.py` | re-asserts every number displayed above against the JSONs |
| `SHA256SUMS` | computed, never typed |

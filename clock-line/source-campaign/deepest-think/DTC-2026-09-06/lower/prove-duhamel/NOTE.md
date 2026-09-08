# NOTE — lower seat `prove-duhamel` — DTC-2026-09-06

**Task.** Prove (ii) — `T(Re) <= c2/log Re` for the mollified tapered plateau — by the
Eulerian/Duhamel route: first-order growth at the innermost shell, plus a remainder bound
`|d_t^2 omega| <~ M^3 log^2(R/rho0) + nu(...)` showing the `O(tau^2)` remainder is a small
fraction of the `O(tau)` gain.

**Answer in one line.** (ii) is **NOT proved**, and the route as specified **cannot** prove it:
the `O(tau^2)` remainder is *not* a small fraction of the `O(tau)` gain — with the stated
`M^3 log^2` bound at an `O(1)` constant the Taylor polynomial never reaches `3/2` at all (the
remainder is already `1.5x` the gain at `x = 1/2`); the route closes only below the constant
**1/4**, i.e. only if `sup|beta| <= (2/3) a0^2`. What does close is a **Riccati** route that uses only the *sign*
of the second-order driver, not its size; that route is a theorem (Theorem D below), it reduces
(ii) to **one scalar inequality** along one trajectory, and it gives `c2 = 2` — not the frame's
`4 log(3/2) = 1.62186`, which is the value of an idealisation (frozen strain) and is not a bound.
Every step of the reduction is proved; the residual hypothesis is not, and is tested numerically
here (16 runs / 39 tracer histories; never violated inside the theorem's own window; the
reversed control fires).

Tier: receipts, single seat. Nothing is promoted. Every number came out of a script in this
folder that I wrote and ran; `SHA256SUMS` is recomputed, never typed. Numerics falsify, never
prove. Nothing outside this folder was written; the estate solver `nsring.py` is imported
read-only and its sha256 is recorded in every result file.

---

## 0. Notation

`nu > 0`; axisymmetric, no swirl; `omega = omega^theta e_theta`; `eta := omega^theta / r`;
`M0 := ||omega_0||_inf`; `a := u^r/r`; `L5 := d_rr + (3/r) d_r + d_zz` (the Laplacian of `R^5`);
`X(t)` a Lagrangian trajectory, `r(t)` its cylindrical radius, `D_t = d_t + u.grad`;

    beta := -(1/r) d_r p ,
    V    := (1/r)(d_rr + (1/r) d_r - 1/r^2 + d_zz) u^r ,
    W    := (L5 eta)/eta .

`L := log(R/rho0)`, `theta := a0 T_d`, `c2 := M0 T_d log Re_E`.

---

## 1. The identities (PROVED — sympy, `d1_clock_algebra.py`, all residuals identically 0)

Re-derived from 3D Navier–Stokes in **Cartesian** coordinates with the axisymmetric no-swirl
ansatz `u = u^r(r,z,t) e_r + u^z(r,z,t) e_z` — nothing quoted:

| id | statement | residual |
|---|---|---|
| I1 | `D_t u^r = -d_r p + nu (d_rr + (1/r)d_r - 1/r^2 + d_zz) u^r` | `0` |
| I2a | `D_t omega^theta = a omega^theta + nu(...)omega^theta` (from the Cartesian vorticity equation) | `0` |
| I2b | `omega^theta = r eta  =>  D_t eta = nu L5 eta` | `0` |
| I3a | `D_t a = D_t(u^r)/r - a^2` | `0` |
| I3b | **`D_t a = -a^2 + beta + nu V`** | `0` |

**I3b is the whole point.** The `a^2` in `D_t a = -a^2 + ...` is *not* a free compounding gain:
it is exactly the `-(u^r)^2/r^2` produced by differentiating `1/r` along the trajectory, and the
genuine second-order driver is the **radial pressure gradient** `beta = -(1/r) d_r p`. (The same
identity was found independently by `sharp/exact-first-order` §4; this is a second derivation from
a different starting point.)

Integrating I2b along `X`:

    d/dt log|omega^theta(X(t),t)| = a + nu W                                      (*)

so the entire question is a lower bound on `int_0^tau a ds` and an upper bound on the `eta`-loss.

---

## 2. What the brief's Duhamel/Taylor route actually needs (PROVED, `d1`)

Duhamel with the second derivative:
`omega(tau) = omega_0 + tau a0 omega_0 + int_0^tau (tau-s) beta omega ds`, and `|omega| <= (3/2)M`
up to the first crossing. With the hypothesis `sup|beta| <= C a0^2` on the window this gives
`omega(tau)/M >= 1 + x - (3/4) C x^2`, `x = a0 tau`, so `(3/2)` is reached iff

    max_x [ x - (3/4) C x^2 ] = 1/(3C) >= 1/2   <=>   **C <= 2/3** .

| `C = sup|beta|/a0^2` | 0 | 0.1 | 1/3 | 0.5 | 2/3 | 0.7 | 1 | 4 |
|---|---|---|---|---|---|---|---|---|
| `theta = a0 tau` | 0.5000 | 0.5203 | 0.5858 | 0.6667 | 1.0000 | FAILS | FAILS | FAILS |
| `c2 = 4 theta` | 2.000 | 2.081 | 2.343 | 2.667 | 4.000 | — | — | — |

**Translation into the brief's units.** `a0 = (M/2)L`, so `a0^2 = M^2L^2/4` and
`|d_t^2 omega| = |beta||omega| <= (3/8) C M^3 L^2`. Hence

> `C <= 2/3`  **⟺**  `|d_t^2 omega| <= (1/4) M^3 log^2(R/rho0)` on the window.

A bound `|d_t^2 omega| <~ M^3 log^2(R/rho0)` **with an O(1) constant does not close the argument**:
at constant 1 (`C = 4`) the `O(tau^2)` remainder is `1.5x` the `O(tau)` gain at `x = 1/2` and there
is **no** `x` at which the Taylor polynomial reaches 3/2. Remainder/gain `= (3/4) C x`:

| `C` | 0.25 | 0.5 | 2/3 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|
| remainder/gain at `x = 1/2` | 0.094 | 0.188 | 0.250 | 0.375 | 0.750 | 1.500 |

**This is the seat's first finding: the premise of the task ("show the `O(tau^2)` remainder is a
small fraction of the `O(tau)` gain") is false as a matter of arithmetic unless the constant in the
second-derivative bound is below 1/4.** And it should be: for the plateau the second-order driver
really is of size `a0^2` (§4 below measures `beta/a^2 = 1.53` at the origin, i.e. `C = 1.53`,
which is a factor 2.3 *outside* the Taylor route's window). The size of `beta` is not the enemy — its *sign* is the friend.

---

## 3. THEOREM D — the Riccati/Duhamel reduction (PROVED)

> **Theorem D.** Let `nu > 0` and let `u` be the smooth solution on `[0,T]` of Navier–Stokes with
> smooth finite-energy axisymmetric no-swirl data, `M0 = ||omega_0||_inf`. Let `X(t)` be the
> trajectory from `x0`, and put `a0 = a(x0,0) > 0`. Suppose that on `[0,tau]`, along `X`,
>
>   **(H2)**  `beta + nu V >= 0`      and      **(H3)**  `nu W >= -lambda` .
>
> Then `a(X(t),t) >= a0/(1 + a0 t)` and hence, by (*),
>
>   `|omega^theta(X(t),t)| >= |omega_0(x0)| (1 + a0 t) e^{-lambda t}` .
>
> Consequently `||omega(tau)||_inf >= (3/2) M0`, and therefore `T_d <= tau`, as soon as
> `(1 + a0 tau) e^{-lambda tau} >= (3/2) M0/|omega_0(x0)|`.
> In particular if `|omega_0(x0)| = M0` and `lambda = 0`, then **`M0 T_d <= 1/(2 a0)`**.

*Proof.* I3b + (H2) give the differential inequality `D_t a >= -a^2` along `X`; comparison with
`y' = -y^2`, `y(0) = a0` gives `a >= a0/(1+a0 t)`. Integrate (*) and use (H3):
`log|omega^theta(X(t))/omega_0(x0)| >= int_0^t a0/(1+a0 s) ds - lambda t = log(1+a0 t) - lambda t`.
Then `||omega(tau)||_inf >= |omega^theta(X(tau),tau)|`. ∎

**Weakened hypotheses, all computed in `d1`.**

* **(H2_q)** `beta + nu V >= -q^2 a0^2` (a *deficit* is allowed): then
  `a >= sqrt(B) tan(arctan(1/q) - sqrt(B) t)`, `B = q^2 a0^2`, and the conclusion survives iff
  `q <= sqrt(4/5) = 0.894427`. `theta(q) =` 0.5000, 0.5060, 0.5433, 0.6742 at `q =` 0, 0.2, 0.5,
  0.8 (closed form cross-checked against RK4 integration of the ODE to `2e-12`).
* **(H2\*)** `beta + nu V >= a^2` (strain non-decreasing): then `theta = log(3/2) = 0.405465`.
* viscosity: `(1+theta) e^{-eps theta} = 3/2`, `eps = lambda/a0`; `d theta/d eps|_0 = 3/4`
  (measured 0.750001); `theta =` 0.5000, 0.5415, 0.5914, 0.7389 at `eps =` 0, 0.05, 0.1, 0.2.

**Corollary (the constant).** For the plateau family, `a0 = kappa M L + O(M)` and
`log Re_E = 2L + O(1)`, so `c2 = M0 T_d log Re_E -> 2 theta / kappa`. With `kappa = 1/2`:

| hypothesis | `theta` | `c2` |
|---|---|---|
| (H2) `beta + nu V >= 0` | 1/2 | **2** |
| (H2\*) `beta + nu V >= a^2` | `log(3/2)` | **`4 log(3/2) = 1.62186`** |
| frozen strain (an idealisation, **not** a bound) | `log(3/2)` | 1.62186 |

With the measured offsets (`a0 = (M/2)L + 0.216773 M`, `log Re_E = 2L - 0.7031660`) the (H2)
constant is approached **from below**: `c2 =` 1.711, 1.849, 1.923, 1.961, 1.980 at `L =` 5, 10,
20, 40, 80 — so `c2 = 2` is an honest ceiling for the whole family, not an asymptote one has to
argue up to.

**Second finding: the frame's `c2 = 4 log(3/2) = 1.62186` is not a bound.** Freezing `a` at `a0`
is legitimate only if `a` does not decay, i.e. only under (H2\*), which is strictly stronger than
(H2). Under (H2) alone the correct constant is exactly **2**. (Numerically (H2\*) does hold for
this datum — §5 — so 1.62186 is likely to be right; but it needs the stronger hypothesis.)

---

## 4. Unconditional structure — why (H2) is the right thing to try (PROVED, `d5`, `d6`)

**L1 (exact, any axisymmetric no-swirl field, any `t`).** At the origin the elliptic factor of the
5D Biot–Savart strain kernel collapses (`k^2 = 0`, `2F1 = 1`, `J5 = (pi/2)/rho'^5`), giving the
closed form

    a(0,0,t) = int K omega^theta d^3x ,     K(r,z) = -(3/(8 pi)) r z / rho^5 ,
             = (3/4) int int omega^theta (-sin^2 phi cos phi) dlog(rho) dphi .

Checked against the grid Poisson solver at `N = 3..6`: rel `6.0e-3, 4.5e-3, 3.6e-3, 3.1e-3`
(the residual is the grid's treatment of the mollified inner edge; it is flat in `h` from `1/4` to
`1/32`, i.e. it is the *datum discretisation*, not the representation).

**L2 (exact).** `eta` odd in `z` ⟹ `psi1` odd ⟹ `u^z = 2 psi1 + r d_r psi1 = 0` on `z = 0`: the
equatorial plane is a **material surface** and the hemispheres never mix. `D_t eta = nu L5 eta`
is a maximum-principle equation, and on `{z > 0}` it carries the boundary value `eta = 0` (by
oddness), so `omega^theta <= 0` in `{z>0}` for **all** `t` if it is at `t = 0`. Since `K <= 0`
there too, `K omega^theta >= 0` **pointwise**, hence

    a(0,0,t) = (3/4) int int |omega^theta| sin^2 phi |cos phi| dlog(rho) dphi  >=  0   for all t,

a *positive* functional of the vorticity: every material element contributes with the same sign,
for all time. Measured: `max omega^theta` over the `z>0` half-domain is `0.00e+00` at every step
of every run; the reversed control flips both signs.

**L3 (exact, sympy).** Write L1 in material form: with `f = r^2 |z| / rho^5`,
`a(0,0,t) = (3/(8 pi)) int |eta_0(x_0)| f(X(t,x_0)) d^3x_0` (inviscid, `eta` conserved, `d^3x`
preserved). Under the pure axisymmetric strain `u = (a r, -2 a z)` that fills the empty hole,

    d log f / dt = 5 a (3 cos^2 phi - 1) ,

which **vanishes exactly at `cos^2 phi = 1/3`** — precisely the angle that maximises the weight
`sin^2 phi |cos phi|`. The dominant part of the strain functional is *stationary under its own
straining field*. This is the structural reason the strain does not decay, and it is exact.

**L4 (exact split of `beta` at the origin, `d6`).** At the origin `u = 0` (axis ∩ material
equatorial plane), so `D_t = d_t` there and `beta(0) = d_t a(0) + a(0)^2`. Using
`d_t omega^theta = a omega^theta - u.grad omega^theta` in L1 and integrating the transport term by
parts (`div u = 0`):

    d_t a(0) = P1 + P2 ,   P1 = int a K omega^theta d^3x  (stretching),
                           P2 = int omega^theta (u.grad K) d^3x  (transport),
    beta(0)  = a^2 + P1 + P2 .

By L2, `K omega^theta >= 0` pointwise, so **`P1 >= 0` wherever `a >= 0` on `supp omega`** — `P1` is
sign-definite. **`P2` is the only sign-indefinite term in `beta(0)`**, and the reduction is

> `beta(0) >= 0`  ⟸  `P2 >= -(a^2 + P1)` .

**L5 (log-range accounting model).** If the shell at log-radius `l` moves as `dl/dt = a(l) =
(M/2)(L-l)` and `omega^theta` intensifies as `e^{Delta}` (material conservation of `eta`), then L1
gives *exactly* `a(t) = (M/2)((1-s)/s)(e^{sL} - 1)`, `s = Mt/2`, whence `P1 -> a^2/2` and
`beta(0)/a^2 -> 3/2`. The plateau's log-range contracts by only `theta/L` over the window while the
vorticity on it intensifies — the intensification wins.

**Measured (`d6`, at `t = 0`, `h = 1/8`):**

| `N` | `L` | `a(0)` | `a^2` | `P1` | `P1/a^2` | `P2` | `P2/a^2` | `P1+P2` | `d_t a` (differenced) | `beta/a^2` |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 2.7726 | 1.36368 | 1.85962 | 0.78105 | 0.4200 | 0.23046 | 0.1239 | 1.01151 | 0.99832 | 1.5439 |
| 5 | 3.4657 | 1.70300 | 2.90023 | 1.25936 | 0.4342 | 0.31370 | 0.1082 | 1.57306 | 1.55933 | 1.5424 |

`P1/a^2` rises toward the model's `1/2`; the sign-indefinite `P2` is only `11–12%` of `a^2` **and
positive**; the split reproduces the independently differenced `d_t a(0)` to `1.3%`.
So `beta(0) = 1.53 a^2 > a^2 > 0` (grid-converged: `1.5797 / 1.5380 / 1.5299` at
`h = 1/4, 1/8, 1/16`): not merely (H2) but (H2\*), with a `53%` margin over (H2\*) and an
unbounded one over (H2).

---

## 5. Numerics — the residual hypothesis, tested (`d2`–`d4`)

**Datum** (admissible, mollified, smooth; `d2`/`dcommon.py`):

    omega^theta = -M tanh( sin(phi)/sin(delta) ) tanh( cos(phi)/w ) Theta(|x|) ,
    Theta(s) = (1/2)[tanh((s-rho0)/w0) - tanh((s-R)/w1)] ,  w0 = .25 rho0, w1 = .10 R,

with `eta = omega^theta/r` written in the axis-regular form `-(1/(s sin delta)) (tanh x)/x ...`,
so `eta` is bounded and smooth on the axis. `delta -> 0, w -> 0` is the bang-bang plateau.

**Solver.** The estate's validated `sharp/viscous-numerics/nsring.py`
(sha256 `3c8e93b2…`), imported read-only. Its Hill-vortex control is reproduced here exactly
(`rel psi1 = 6.116e-3 / 1.677e-3 / 4.405e-4`, Richardson `a_nose -> 0.19918` vs exact `1/5`, axis
peak `2.49832` vs exact `5/2`), and its operator control (`3/r -> 1/r`) **fires** at 183–185%.

**Datum validation against an independent instrument.** Fitting `a(s) = kappa log(R/s) + c` on the
datum's own support:

| `delta` | 30° | 15° | 7.5° | `w`-sweep at 7.5°: `w =` 0.30 / 0.20 / 0.12 / 0.08 |
|---|---|---|---|---|
| `kappa` here | 0.45317 | 0.47481 | 0.48137 | 0.45717 / 0.48137 / 0.49609 / **0.50091** |
| elliptic instrument (`sharp/refuter-correctness` §6) | 0.48363 | 0.49781 | 0.49972 | — |

so the datum reaches `kappa = 1/2` as the mollification sharpens, agreeing with the independent
elliptic-integral instrument; the residual at `w = 0.2` is entirely the equatorial smoothing.
`E/(M^2 R^5) = 0.1418961` constant to 6 digits over `N = 3..6` (sharp shell: `0.1724040`).

**The test.** By I3b, (H2) is *equivalent* to `a(X(t),t)(1 + a0 t)/a0 >= 1`, and it implies
`(r(t)/r0)/(1 + a0 t) >= 1`. Both are read straight off the run — no differentiation. Pre-declared
falsification rule (before reading any number): (H2) is refuted if either ratio dips below `1` by
more than `1e-2` on a forward run, and the test is non-diagnostic unless the reversed control
violates the floor.

**Coverage.** 16 runs, 39 tracer histories: `N = 3,4,5` octaves, `Re0 = M rho0^2/nu` from 1 to 400 (so `rho0/sqrt(nu/M)` from 1 to 20), two grid resolutions, and one sign-reversed control.

**TEST 1 — the Riccati floor.** Minimum over `t > 0` of the two ratios, per run (innermost tracer `s0 = 1.5 rho0`; `*` marks the reversed control):

| run | `log Re_E` | `a0` | `min_{t>0} a(t)(1+a0 t)/a0` | `min_{t>0} (r/r0)/(1+a0 t)` | window `a0 t_end` | verdict |
|---|---|---|---|---|---|---|
| `N3` | 8.021 | +0.65478 | 1.018324 | 1.000062 | 0.384 | HOLDS |
| `N4` | 9.407 | +0.99303 | 1.012035 | 1.000031 | 0.398 | HOLDS |
| `N5` | 10.794 | +1.33209 | 1.007479 | 1.000013 | 0.402 | HOLDS |
| `N4-fine` | 9.408 | +0.99100 | 1.006022 | 1.000009 | 0.397 | HOLDS |
| `N4-Re400` | 10.794 | +0.99303 | 1.012042 | 1.000031 | 0.392 | HOLDS |
| `N4-Re25` | 8.021 | +0.99303 | 1.012008 | 1.000031 | 0.435 | HOLDS |
| `N4-REVERSED(control)*` | 9.407 | -0.99303 | 0.112991 | 1.000038 | -0.746 | floor VIOLATED (control fires) |
| `N3-Re0=1` | 3.416 | +0.65478 | 1.001054 | 1.000000 | 0.873 | HOLDS |
| `N3-Re0=4` | 4.802 | +0.65478 | 1.004763 | 1.000003 | 0.874 | HOLDS |
| `N3-Re0=16` | 6.188 | +0.65478 | 1.018178 | 1.000062 | 0.495 | HOLDS |
| `N4-Re0=1` | 4.802 | +0.99303 | 1.001508 | 1.000000 | 0.994 | HOLDS |
| `N4-Re0=4` | 6.188 | +0.99303 | 1.006395 | 1.000008 | 0.633 | HOLDS |
| `N4-Re0=16` | 7.575 | +0.99303 | 1.011987 | 1.000031 | 0.447 | HOLDS |
| `N5-Re0=1` | 6.188 | +1.33209 | 1.001936 | 1.000001 | 0.832 | HOLDS |
| `N5-Re0=4` | 7.575 | +1.33209 | 1.007399 | 1.000012 | 0.531 | HOLDS |
| `N5-Re0=16` | 8.961 | +1.33209 | 1.007462 | 1.000012 | 0.438 | HOLDS |

Worst forward minimum over all 39 tracer histories and all times: **1.000000** (the floor is `1`; the solver's own max-principle violation is `1.5e-3..4.3e-3`). The reversed control reaches **0.1130**, so the test is diagnostic. **(H2) is not refuted by any of these runs.**

**TEST 2 — `beta/a^2` along the trajectory.** Theorem D needs (H2) only on `[0, tau]`, and `tau = theta/a0` with `theta <= 1/2` under (H2). Both the windowed and the unrestricted minima are reported; neither is withdrawn.

* worst `beta/a^2` inside `a0 t <= 1/2` (the theorem's window): **+0.0140** — **(H2) holds on every forward tracer inside its own window**;
* worst inside `a0 t <= 3/4`: **-0.2366**;
* worst with no restriction at all: **-0.2366**.

All negative values occur at `a0 t > 0.5` and only in the two runs at the deepest viscosity (`N = 3`, `Re0 = 1`) that **never reach `(3/2)M` at all** — by then the datum has lost 87% of its `eta` and the mechanism has failed viscously, so the window is fictitious there.

Per-run, innermost tracer (`s0 = 1.5 rho0`):

| run | `beta/a^2` early | min inside `a0 t <= 1/2` | min over the whole run | `eta_end/eta_0` | `lambda_eff/a0` |
|---|---|---|---|---|---|
| `N3` | 2.373 | +1.2819 | +1.282 | 0.9637 | 0.0962 |
| `N4` | 2.079 | +1.3241 | +1.324 | 0.9699 | 0.0767 |
| `N5` | 1.937 | +1.3505 | +1.351 | 0.9730 | 0.0683 |
| `N4-fine` | 2.103 | +1.3283 | +1.328 | 0.9849 | 0.0384 |
| `N4-Re400` | 2.080 | +1.3402 | +1.340 | 0.9796 | 0.0526 |
| `N4-Re25` | 2.074 | +1.2414 | +1.241 | 0.9170 | 0.1992 |
| `N3-Re0=1` | 2.055 | +0.2300 | -0.163 | 0.0683 | 3.0737 |
| `N3-Re0=4` | 2.315 | +0.7897 | +0.349 | 0.1669 | 2.0494 |
| `N3-Re0=16` | 2.354 | +0.9810 | +0.981 | 0.7509 | 0.5791 |
| `N4-Re0=1` | 1.939 | +0.8334 | +0.616 | 0.0681 | 2.7039 |
| `N4-Re0=4` | 2.049 | +1.0322 | +0.858 | 0.4222 | 1.3622 |
| `N4-Re0=16` | 2.071 | +1.2045 | +1.204 | 0.8702 | 0.3108 |
| `N5-Re0=1` | 1.857 | +1.0154 | +0.849 | 0.1280 | 2.4726 |
| `N5-Re0=4` | 1.916 | +1.1380 | +1.111 | 0.6114 | 0.9271 |
| `N5-Re0=16` | 1.932 | +1.2810 | +1.281 | 0.9131 | 0.2079 |

The `N`-trend at the viscous floor is the informative one: inside the window, `min beta/a^2 = 0.230` (`N=3`), `0.833` (`N=4`), `1.015` (`N=5`) — rising toward the origin value `1.53` as octaves are added, which is what the structure of §4 predicts.

`beta >= a^2` at every time on 25 of 36 forward tracers — i.e. on those the strain is non-decreasing and even (H2\*) holds, which is what licenses the frozen-strain constant `4 log(3/2)`.

**TEST 3 — the clock actually delivered.**

| run | `log Re_E` | `log(R/rho0)` | `M T_d` | `c2 = M T_d log Re_E` | `theta = a0 T_d(tracer)` |
|---|---|---|---|---|---|
| `N3` | 8.021 | 2.079 | 0.50609 | 4.0593 | 0.3688 |
| `N4` | 9.407 | 2.773 | 0.34491 | 3.2446 | 0.3772 |
| `N5` | 10.794 | 3.466 | 0.26181 | 2.8259 | 0.3815 |
| `N4-fine` | 9.408 | 2.773 | 0.34761 | 3.2701 | 0.3639 |
| `N4-Re400` | 10.794 | 2.773 | 0.33884 | 3.6573 | 0.3703 |
| `N4-Re25` | 8.021 | 2.773 | 0.37613 | 3.0170 | 0.4093 |
| `N3-Re0=1` | 3.416 | 2.079 | never | -- | -- |
| `N3-Re0=4` | 4.802 | 2.079 | never | -- | -- |
| `N3-Re0=16` | 6.188 | 2.079 | 0.64989 | 4.0216 | -- |
| `N4-Re0=1` | 4.802 | 2.773 | never | -- | -- |
| `N4-Re0=4` | 6.188 | 2.773 | 0.55402 | 3.4285 | -- |
| `N4-Re0=16` | 7.575 | 2.773 | 0.39065 | 2.9590 | -- |
| `N5-Re0=1` | 6.188 | 3.466 | 0.54240 | 3.3566 | -- |
| `N5-Re0=4` | 7.575 | 3.466 | 0.34558 | 2.6177 | -- |
| `N5-Re0=16` | 8.961 | 3.466 | 0.28371 | 2.5423 | 0.4245 |

`theta = a0 T_d` measured on the innermost tracer: mean **0.3851** (range 0.3639–0.4245) against Theorem D's (H2) bound `1/2` and the (H2\*) value `log(3/2) = 0.4055`. Every measured `theta` is **below `1/2`**, and most are below `log(3/2)` — the strain grows over the window, exactly as `beta > a^2` predicts, so the frozen-strain clock is conservative for this datum rather than optimistic.

**Reading `c2` honestly.** The frame's `c2` is defined at the viscous floor `rho0 = sqrt(nu/M)`, i.e. `Re0 = 1`. Runs at `Re0 >> 1` carry an inflated `log Re_E` for their `log(R/rho0)` and therefore **overstate** `c2`. At the floor itself (`Re0 = 1`, `N = 5`) the datum still reaches `(3/2)M` and gives `c2 = 3.357` at `log Re_E = 6.188`, falling with `N`; at `N = 3,4` with `Re0 = 1` the inner shells diffuse away before `(3/2)M` is reached (`T_d = never`) — the mechanism needs enough octaves before it beats its own viscous floor, which is the frame's `Re >= Re_*`.

**Where the `eta`-loss actually comes from.** On a plateau interior `nu L5 eta/eta = -nu/r^2` exactly, an `O(M)` rate at `r ~ sqrt(nu/M)`. The measured loss is several times larger, and the sweep shows why: at `Re0 = 100` the tracer at `s0 = 1.5 rho0` loses 3.0% while the one at `s0 = 3 rho0` loses 0.1% — the cost is the **mollification layer at the inner edge** (`w0 = rho0/4`), not the plateau. Tracking one or two `rho0` further out buys it back for an `O(M)` loss in `a0`, i.e. `O(1/L)` in `c2`. This is a real bookkeeping item for (H3) that the frame does not price.

---

## 6. Per-step status

| # | step | status | evidence |
|---|---|---|---|
| 1 | identities I1, I2a, I2b, I3a, I3b (in particular `D_t a = -a^2 + beta + nu V`) | **PROVED** | sympy, all residuals identically `0` (`d1`) |
| 2 | the Taylor/Duhamel route closes **iff** `sup\|beta\| <= (2/3) a0^2`, i.e. iff `\|d_t^2 omega\| <= (1/4) M^3 log^2(R/rho0)` on the window | **PROVED** (elementary) | `d1`; and the plateau's own `beta = 1.53 a^2` is OUTSIDE that window, so the brief's route is not merely unproved, it is **unavailable** for this datum |
| 3 | **Theorem D**: (H2)+(H3) ⟹ `\|omega(X(t))\| >= \|omega_0(x0)\|(1+a0 t)e^{-lambda t}`, so `M0 T_d <= 1/(2a0)` | **PROVED** | comparison for `y' = -y^2` + integration of (*); the (H2_q) closed form cross-checked against RK4 to `2e-12` |
| 4 | `c2 = 2` under (H2); `c2 = 4 log(3/2) = 1.62186` under (H2\*); with the measured offsets `c2 <= 2` at every finite `L` | **PROVED given** the two measured constants `a0 = kappa M L + c_a M` and `log Re_E = 2L + c_E` | `d1`, using `kappa = 1/2`, `c_a = 0.216773` (refuter's elliptic instrument) and `c_E = -0.7031660` (`sharp/exact-first-order` c6) |
| 5 | the frame's `c2 = 4 log(3/2)` is **not a bound** — freezing `a` needs (H2\*), strictly stronger than (H2) | **PROVED** | `d1` R0 note |
| 6 | L1 exact origin representation `a(0,0,t) = int K omega^theta`, `K = -(3/8pi) r z/rho^5`, valid at every `t` | **PROVED** | kernel collapse `k^2 = 0`; checked against an independent grid Poisson solve to `3.1e-3` at `N = 6` (`d5`) |
| 7 | L2 equatorial plane is material; `omega^theta <= 0` in `{z>0}` for all `t`; hence `a(0,0,t) >= 0` for all `t`, as a positive functional | **PROVED** | `u^z = 2 psi1 + r d_r psi1 = 0` on `z=0` for odd `eta`; maximum principle for `D_t eta = nu L5 eta` on `{z>0}` with boundary value 0. Measured: `max omega^theta` in `z>0` is `0.00e+00` at every step; reversed control flips it (`d5`, `d6`) |
| 8 | L3 the material weight `f = r^2\|z\|/rho^5` satisfies `d log f/dt = 5a(3cos^2 phi - 1)`, vanishing exactly at the weight-maximising angle `cos^2 phi = 1/3` | **PROVED** | sympy (`d5`) |
| 9 | L4 exact split `beta(0) = a^2 + P1 + P2`, with `P1 = int a K omega^theta >= 0` wherever `a >= 0` on `supp omega`; **`P2` is the only sign-indefinite term** | **PROVED** | integration by parts under `div u = 0`; split reproduces the independently differenced `d_t a(0)` to `0.7–2.3%` (`d6`) |
| 10 | `beta(0) >= 0` ⟸ `P2 >= -(a^2 + P1)` — the reduction of (H2) at the origin to one inequality on the transport term | **PROVED** (restatement of 9) | — |
| 11 | `P1 -> a^2/2` (log-range accounting) | **MODEL**, leading constant confirmed | `P1/a^2 = 0.4200, 0.4342` at `L = 2.77, 3.47`, grid-converged (`0.4351/0.4342/0.4341` at `h = 1/4,1/8,1/16`) (`d6`) |
| 12 | `\|P2\|` is small (`~0.10 a^2`, positive) | **NUMERICS ONLY** | `P2/a^2 = 0.1401 / 0.1082 / 0.1000` at `h = 1/4,1/8,1/16` (`d6`) — no a-priori bound attempted |
| 13 | (H2) **at the origin**, `t = 0`: `beta(0)/a^2 = 1.53` | **NUMERICS ONLY**, margin `1.53` over the required `0`, and `0.53` over (H2\*) | `1.5390/1.5383/1.5380` at `N = 3,4,5`; grid `1.5797/1.5380/1.5299` at `h = 1/4,1/8,1/16` (`d5`, `d6`) |
| 14 | transfer of the L1/L2/L4 structure from the origin to the material innermost shell (relative error `O(1/L)`) | **NOT PROVED** | plausibility only: the sources are at `rho >> \|x0\|`, and the refuter's independent instrument puts the difference at a *constant* `+0.216773 M`, i.e. relative `O(1/L)` |
| 15 | **(H2) along the trajectory, over the whole window** — the residual gap | **NOT PROVED** | 16 runs / 39 tracer histories: the equivalent floor `a(t)(1+a0t)/a0 >= 1` is never violated (worst forward minimum `1.000000`), and `beta/a^2 >= +0.0140` inside the theorem's own window `a0 t <= 1/2`. Outside that window `beta` does go negative (`-0.2366`) — in the two runs that never reach `(3/2)M` at all. Reversed control reaches `0.1130`, so the test is diagnostic. §5 |
| 16 | **(H3)** the viscous `eta`-loss along the trajectory | **NOT PROVED** | exact at `t=0` on a plateau: `nu L5 eta/eta = -nu/r^2` (an `O(M)` rate at `r ~ sqrt(nu/M)`, so `eps = O(1/L)` and the cost in `c2` is `O(1/log Re)`). Measured `eta`-loss is dominated instead by the **mollification layer** at the inner edge, not by the plateau interior — see §5 |
| 17 | **(ii)** `T(Re) <= c2/log Re` | **NOT PROVED** | — |
| 18 | the conflict with Bradshaw–Farhat–Grujić ARMA 2019 Thm 10 | **not engaged** | nothing here refutes Thm 10; this seat proves no upper bound on `T(Re)` |


---

## 6b. What would close the gap — and where it meets the BFG conflict

Everything except items 12–16 is proved. The single remaining *mathematical* target is an
a-priori bound on **one** integral:

    P2 = int omega^theta (u . grad K) d^3x ,   K = -(3/(8 pi)) r z / rho^5 ,  |grad K| <~ 1/rho^4 ,

together with its transfer from the origin to the material innermost shell (item 14, a relative
`O(1/L)` step). `beta(0) >= 0` follows from `P2 >= -(a^2 + P1)`, and `P1 >= 0` is already proved.

**This is the same integral the published conflict turns on.** Bradshaw–Farhat–Grujić's
Theorem 10 rests (p.10) on `int |f(x)|/(|x|+1)^4 dx <= c ||f||_BMO`, and the two refuter seats
showed that this fails on precisely this plateau family, the left side growing like
`log(R/rho0)` while the BMO norm stays bounded. `P2` is an integral of the *same* shape — a
bounded vorticity against a `rho^{-4}` kernel — and it too is `O(M^2 L^2)`, i.e. it carries a
factor `L` more than a naive `L^1`-against-`L^inf` count would give (`|u| ~ a rho ~ M L rho`
supplies one `L`, `int d^3x/rho^3 = 4 pi L` supplies the other), which is exactly the `log` that
BFG's step drops. Measured here it is `+0.10 a^2`: the *right* size, the *favourable* sign, and
`15x` smaller than the margin it has to fit inside. Whoever bounds `P2` bounds BFG's step too.

**Second target.** `nu V` and the `eta`-loss are priced here only at `t = 0` and only on a plateau
interior; the runs show the real cost is carried by the *mollification layer* at the inner edge,
not by the plateau, and that tracking `2–3 rho0` out instead of `1.5 rho0` removes almost all of
it at an `O(M)` cost in `a0`. A clean (H3) should be stated at a fixed distance from the layer,
not at the layer.

**What this seat does NOT claim.** No upper bound on `T(Re)` is proved here, so nothing here
refutes BFG Theorem 10, and "sharp up to constants" remains unavailable.

---

## 7. Files

| file | what it establishes |
|---|---|
| `d1_clock_algebra.py` / `d1_log.txt` / `d1_results.json` | sympy proof of I1, I2a, I2b, I3a, I3b (all residuals 0); the R0/R1/R2 clock constants; `C <= 2/3` ⟺ `|d_t^2 omega| <= (1/4) M^3 log^2`; the viscous correction |
| `dcommon.py` | the `delta`-tapered mollified plateau datum; read-only import of the estate solver |
| `d2_datum_probe.py` / `d2_log.txt` / `d2_results.json` | Hill control reproduced, operator control fires; `kappa(delta, w) -> 1/2` vs the independent elliptic instrument; `E/(M^2R^5)` |
| `d3_track.py` / `d3_log.txt` / `d3_results.json` | tracer runs `N = 3..5`, `Re0 = 25..400`, grid refinement, reversed control |
| `d3b_lowre.py` / `d3b_log.txt` / `d3b_results.json` | the same at the brief's viscous floor `Re0 = 1, 4, 16` |
| `d4_analyse.py` / `d4_log.txt` / `d4_results.json` | the pre-declared floor test, `beta/a^2` along trajectories, the delivered clock |
| `d5_structure.py` / `d5_log.txt` / `d5_results.json` | L1 representation vs grid solver; L2 sign law along a run; L3 sympy stationarity; `beta(0)/a^2` at the origin |
| `d6_split.py` / `d6_log.txt` / `d6_results.json` | the exact `beta(0) = a^2 + P1 + P2` split and its internal consistency check |

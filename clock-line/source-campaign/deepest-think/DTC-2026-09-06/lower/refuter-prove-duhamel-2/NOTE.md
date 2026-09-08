# NOTE — refuter seat `refuter-prove-duhamel-2` — DTC-2026-09-06

**Target.** `lower/prove-duhamel/` (Theorem D; the L1–L4 "unconditional structure"; the
16-run / 39-tracer (H2) test; the corollary `c2 = 2`).

**Verdict: REFUTED (MAJOR).** The seat's central theorem is correct and its lead negative
finding is correct, but **two of its three headline "KEY FINDS" fail**, **two steps carrying the
status `PROVED` are false as written**, and **the numeric gate does not test the hypothesis it
says it tests — its own data contain the counterexample.** Nothing here touches the seat's
correct refusals ((ii) NOT PROVED; BFG Thm 10 not engaged), which stand.

Everything below came out of a script in this folder that I wrote and ran. `SHA256SUMS` is
recomputed. Numerics falsify, never prove.

---

## 1. What survives (checked, not conceded)

| claim | how checked | result |
|---|---|---|
| identities I1, I2a, I2b, I3a, I3b | re-ran `d1_clock_algebra.py` from a copy | all residuals `0`, exit 0, RK4 cross-check `7e-13 … 2e-12` — **reproduces** |
| Taylor route closes iff `C = sup|beta|/a0^2 <= 2/3`, i.e. `|d_t^2 omega| <= (1/4)M^3 log^2` | re-derived: `omega'' = (beta+nuV)omega` (the `a^2` cancels), `max_x[x-(3/4)Cx^2]=1/(3C)>=1/2` | **correct.** This is the seat's real contribution: the brief's route is unavailable, not merely unproved |
| `K = -(3/(8 pi)) r z/rho^5` (L1) | **independent** direct 3D Biot–Savart on a circular filament, no stream function, no grid solve, no elliptic table (`x4`) | agrees to **1.6e-11** over 6 source positions — **CONFIRMED** |
| `a(0)` values | independent elliptic-ring instrument (`x6`) | `1.37117 / 1.71069` at `N=4/5` vs the seat's `1.36368 / 1.70300` (0.5%) — **CONFIRMED** |
| `kappa = 1/2` at the origin for the sharp shell | analytic: `a(0) = (3/4) M L * int_0^pi sgn(cos phi) sin^2 phi cos phi dphi`, angular integral `= 2/3` exactly (`x4b`) | `a(0) = (M/2)L` **EXACTLY**; so the origin-to-material-shell offset is exactly the frame's `+0.216773 M`. Confirms the seat's GAP-2 *size* claim |
| Theorem D's proof | comparison for `y'=-y^2` + integration of `(*)` | **valid** |
| `c2 = 2` under (H2) vs `4log(3/2)` under (H2\*); the frame's `1.62186` needs the stronger hypothesis | arithmetic | **correct** |
| the L2 sign law transfers to the equatorial plane | **attempted refutation, 28000-point kernel scan** (`x5`,`x5b`) | my refutation **FAILED** — the sign law holds there. Credit to the seat |
| `d6`'s numbers | recomputed inside `x3` | `P1/a^2 = 0.4200 (N=4), 0.4342 (N=5)`, `P2/a^2 = 0.1239 / 0.1082` — **reproduces exactly** |

---

## 2. REFUTATION 1 — "(H2) is *equivalent* to `a(t)(1+a0 t)/a0 >= 1`" is FALSE, and the seat's own runs are the counterexample

§5 of the target NOTE: *"By I3b, (H2) is **equivalent** to `a(X(t),t)(1 + a0 t)/a0 >= 1`"*, and TEST 1
is then run as the test of (H2).

The implication is one-way. `(H2)` on `[0,t]` ⟹ `F_a >= 1` by comparison; the converse is false,
because `a` can be *above* the Riccati floor from an earlier surplus while `beta + nu V < 0`.

**Counterexample from `d4_results.json` itself** (`x1`, part A):

| tracer | `min (beta+nuV)/a^2` | `min F_a` | `min F_r` | TEST 1 verdict |
|---|---|---|---|---|
| `N3-Re0=1`, `s0=1.5` | **−0.1632** | 1.001054 | 1.000000 | **HOLDS** |
| `N3-Re0=1`, `s0=2.0` | **−0.2366** | 1.001062 | 1.000000 | **HOLDS** |

(H2) is violated on those two tracers — the seat's own TEST 2 says so — and TEST 1 certifies HOLDS.

**How much violation the gate tolerates** (`x1b`, 112 ODE runs of `a' = -a^2 + B`, `B = b_p` on
`[0,t1]` then `b_n`, gate `= (min F_a >= 0.99 and min F_r >= 0.99)`): with `b_p = 1.5` — *less*
than the seat's own measured early `beta/a^2 ~ 1.9–2.4` — the gate returns **HOLDS at
`b_n = -3.0 a0^2`**, i.e. an arbitrarily large (H2) violation over the last 18% of the window.
The gate tests the **conclusion** of Theorem D, not its **hypothesis**.

---

## 3. REFUTATION 2 — FL-043: half the gate cannot fire, and it is the half that is reported

The pre-declared rule is *"(H2) is refuted if **either** ratio dips below 1 by more than 1e-2"*,
with `F_r = (r/r0)/(1+a0 t)` as the second ratio. Measured (`x1`, part B):

| | forward (36 tracers) | reversed control (3 tracers) | diagnostic? |
|---|---|---|---|
| `min F_a` | 1.001054 | **0.112991** | yes |
| `min F_r` | 0.999999776 | **1.000017** | **no — the control passes** |

`F_r` returns HOLDS on every run in the study *including the control built to make it fire*.
And the NOTE's headline — *"Worst forward minimum over all 39 tracer histories and all times:
**1.000000** (the floor is 1)"* — is `min(min F_a, min F_r) = 0.999999776`, i.e. it is the
`F_r` number: the non-diagnostic half, **rounded up through the floor it is being compared to**.
The diagnostic half's worst forward value is `1.001054`.

Third, `F_a >= 1` is automatic whenever `a` is non-decreasing (`a >= a0` and `1+a0 t >= 1`), and
`a` is non-decreasing on **25 of 36** forward tracer histories — the same 25/36 the NOTE reports
separately as "`beta >= a^2` at every time" (by I3b those are the *same statement*, `D_t a >= 0`,
not two corroborations). On those 25 the gate is vacuous.

---

## 4. REFUTATION 3 — the "SECOND KEY FIND" fails twice

> *"…gives `beta(0) = a^2 + P1 + P2` with `P1 >= 0` **proved**, so (H2) at the origin reduces to
> **ONE inequality on ONE integral**: `P2 >= -(a^2+P1)`."*

### 4a. The premise of `P1 >= 0` is false for the seat's own datum
`P1 = int a K omega^theta`. `K omega^theta >= 0` is proved (L2). `P1 >= 0` therefore needs
`a >= 0` on `supp omega` — stated in the status table (item 9) as the qualifier and nowhere
established. Measured (`x2`, the seat's datum, `delta = 7.5 deg`, `w = 0.12`, `h = 1/8`):

| `N` | `L` | `min a on supp` | `max a on supp` | fraction of supp with `a<0` | `P1` | `P1` restricted to `{a<0}` |
|---|---|---|---|---|---|---|
| 3 | 2.0794 | **−0.44743** | 1.03594 | **56.2 %** | 0.41783 | −0.01121 |
| 4 | 2.7726 | **−0.44981** | 1.37550 | **56.0 %** | 0.78105 | −0.01121 |
| 5 | 3.4657 | **−0.45043** | 1.71488 | **55.9 %** | 1.25936 | −0.01121 |
| 6 | 4.1589 | **−0.45063** | 2.05422 | **55.9 %** | 1.85278 | −0.01121 |

`a < 0` on **more than half the support** (the outer shells, `rho` from `~R/2` out to `~2.7R`,
where `|omega|` is still `0.999`), at every octave count and at every support threshold from
`1e-1` to `1e-6`. `P1` is positive **as a measured number**, not by a proved sign law: it is a
difference, whose negative piece (`−0.01121`, *`L`-independent*) happens to be small against a
piece that grows like `a^2 ~ L^2`. That is an asymptotic argument, not the claimed proof.

### 4b. The "exact split" is an EULER identity — a term is missing
`d6_split.py` derives the split from `d_t omega^theta = a omega^theta - u.grad omega^theta`
(its own docstring): **no viscous term.** But the quantity (H2) constrains is `beta + nu V`, and
at the origin I3b gives `beta + nu V = d_t a + a^2` with `d_t a` along the *Navier–Stokes*
evolution. The NS vorticity equation adds `nu(d_rr + (1/r)d_r - 1/r^2 + d_zz) omega^theta`, so
the exact split is

    beta(0) + nu V(0) = a^2 + P1 + P2 + P3 ,   P3 = nu * int K r (Delta_5 eta) d^3x ,

and `P3` is **negative**. Measured (`x3`, same recipe and same solver as `d6`):

| `N` | `Re0` | `P3/a^2` | `|P1+P2 − d_t a|/|d_t a|` | with `P3` |
|---|---|---|---|---|
| 4 | 100 | −0.0016 | 1.32e-2 | 1.03e-2 |
| 5 | 100 | −0.0010 | 8.81e-3 | 6.93e-3 |
| 4 | 4 | −0.0394 | 9.52e-2 | 1.58e-2 |
| 5 | 4 | −0.0253 | 6.25e-2 | 1.30e-2 |
| 4 | **1** | **−0.1576** | **4.89e-1** | 5.73e-2 |
| 5 | **1** | **−0.1011** | **2.70e-1** | 3.34e-2 |

At `Re0 = 100` — the only viscosity `d6` was run at — `P3` accounts for most of the `0.9–1.3 %`
residual the seat left unattributed. At `Re0 = 1`, i.e. `rho0 = sqrt(nu/M)`, **which is the only
regime in which the frame's `c2` is defined**, `P3/a^2 = −0.10 … −0.16`: `beta + nu V` at the
origin is `1.44 a^2` (N=5) / `1.39 a^2` (N=4), not the `1.54 a^2` the inviscid split reports, and
the split's residual falls from 27–49 % to 3–6 % once `P3` is put back.

So the reduction is **`P2 + P3 >= -(a^2 + P1)`** with `P1` not sign-definite and `P3` known
negative — three terms, not "ONE inequality on ONE integral".

---

## 5. REFUTATION 4 — GAP 2 is structural, not the `O(1/L)` size step the NOTE calls it

The NOTE prices the origin→material-shell transfer as *"a relative `O(1/L)` step"*, citing the
constant `+0.216773 M` offset. The **size** does transfer (§1 confirms it exactly). The **sign
structure** — which is the entire content of L2 and L4 — does not, because it lives precisely
where the vorticity vanishes.

Scan of the exact axisymmetric strain kernel (vortex-ring stream function via complete elliptic
integrals, validated against direct filament Biot–Savart to `2e-10`; `x5b`, 28 000 sources per
field point). "Violation" = kernel has the sign opposite to the origin law `sgn K = -sgn z`:

| field point | wrong-sign sources | worst magnitude |
|---|---|---|
| origin | 0 (proved) | — |
| equatorial plane `z=0`, `r0=1` | 128/28000 (all `~3e-8`, round-off) | 2.9e-8 |
| equatorial plane `z=0`, `r0=4` | 54/28000 (`~2e-9`) | 1.8e-9 |
| **`z0 = 0.5`** | **15053/28000 (54 %)** | **+1.14** |
| **`z0 = 2.0`** | **18239/28000 (65 %)** | **+0.25** |
| **near the axis, `(0.3, 0.95)`** | **15964/28000 (57 %)** | **+3.04** |

And on the equator the datum vanishes identically: `omega^theta = -M tanh(sin phi/sin delta)
tanh(cos phi/w) Theta(s)` gives `|omega| = 5.0e-16` at `phi = 90 deg` (`x6`, part B) — as it must
for *any* continuous field odd in `z`.

So: the origin carries no vorticity and never moves; the equatorial plane carries no vorticity;
and every point with `|omega_0| = M0` is strictly off the equator, where the kernel is
sign-indefinite with `O(1)` wrong-sign contributions. At the seat's **own tracer seed**
(`s = 1.5 rho0`, `phi = 45 deg`, from `d3_track.py`) the same strain functional is a difference of
comparable pieces (`x6`, part C; my instrument reproduces the seat's `a0 = 0.99303` as `1.00180`,
0.9 %):

| field point | `a = int K omega` | positive part | negative part | `|neg|/a` |
|---|---|---|---|---|
| origin (`N=4`) | 1.37117 | 1.37117 | 0.00000 | 0.000 |
| equator `s=1.5 rho0` (`omega=0` there) | 1.33421 | 1.33421 | −0.00000 | 0.000 |
| **seed `phi=45 deg` (`N=4`)** | 1.00180 | 1.25029 | **−0.24849** | **0.248** |
| **seed `phi=45 deg` (`N=5`)** | 1.33738 | 1.58156 | **−0.24417** | **0.183** |
| `phi=54.74 deg` (L3's angle, `N=5`) | 1.38373 | 1.57930 | −0.19557 | 0.141 |

There is no version of "`P1 >= 0`" available at any point that carries vorticity.

---

## 6. REFUTATION 5 — the "THIRD KEY FIND" (the BFG link) is an analogy sold as an identification

> *"`P2 = int omega^theta (u.grad K)`, `|grad K| ~ rho^-4`, is the **SAME shape** as the defective
> BFG step `int |f|/(|x|+1)^4 <= c||f||_BMO` … **same `rho^-4` weight** … Whoever bounds `P2`
> bounds BFG's step."*

Three defects.
1. `P2`'s integrand carries a **velocity factor**, absent from BFG's. With `|u| ~ |a| rho` the
   effective weight is `rho^{-3}`, not `rho^{-4}`. The NOTE's own §6b says exactly this
   (*"`|u| ~ a rho ~ M L rho` supplies one `L`, `int d^3x/rho^3 = 4 pi L` supplies the other"*),
   contradicting the `rho^-4` claim in the same document.
2. BFG's kernel `(|x|+1)^{-4}` carries a **fixed length scale** (the `+1`), and that broken
   scale-covariance is the refuters' located defect. `grad K` is homogeneous of degree `-4` with
   no scale. Different structural failure.
3. *"Whoever bounds `P2` bounds BFG's step"* is backwards: bounding one specific integral for one
   specific family does not bound an inequality asserted for all `f in BMO`.

The observation that both integrals carry an extra `log` is fair. The identification is not.

---

## 7. Smaller finds

* **Dimensional slip in the theorem statement.** *"`M0 T_d <= 1/(2 a0)`"* — LHS dimensionless,
  RHS a time. Correct: `T_d <= 1/(2 a0)`, i.e. `M0 T_d <= M0/(2 a0)`. It appears in the theorem,
  the abstract, and the `statement` field; the corollary uses the right quantity, so it is a typo,
  but it is the theorem's headline conclusion.
* **The corollary mixes two data.** It uses `c_E = -0.7031660`, which is `(2/5)log(0.1724040)`, the
  **sharp-shell** energy. The seat's own `d2` measures `E/(M^2R^5) = 0.1418961` for the mollified
  datum, giving `c_E = -0.7810641` (`x4c`). Mismatch `0.078` in `log Re_E`, in the direction that
  inflates `c2`.
* **`c2 = 2` is an unattained limit, not "an honest ceiling for the whole family".** `c2 = 2 theta/kappa`
  and `kappa = 1/2` only in the singular limit `delta, w -> 0` (where `||eta||_inf` diverges).
  The seat's own `kappa` measurements give ceilings `1/kappa =` **2.0774** for the datum the runs
  actually used (`w = 0.20`), 2.2067 at `delta = 30 deg`, 1.9964 at `w = 0.08`. Every admissible
  member built in the seat has `kappa < 1/2`, hence `c2 > 2`.
* **No measurement supports `c2 = 2`.** TEST 3 delivers `c2 = 2.5423 … 4.0593`, and at the frame's
  own viscous floor `Re0 = 1` only one run out of three reaches `(3/2)M` at all (`N=5`,
  `c2 = 3.3566`). The NOTE reports this honestly; the `key`'s *"which gives `c2 = 2`"* does not.

---

## 8. Corrected statement that survives

> **Theorem D (unchanged, correct).** `nu > 0`, smooth finite-energy axisymmetric no-swirl data,
> `a = u^r/r`, `beta = -(1/r)d_r p`, `V`, `W` as in the target. Along the trajectory `X` from `x0`
> with `a0 = a(x0,0) > 0`, if `beta + nu V >= 0` and `nu W >= -lambda` on `[0,tau]`, then
> `a(X(t),t) >= a0/(1+a0 t)` and `|omega^theta(X(t),t)| >= |omega_0(x0)|(1+a0 t)e^{-lambda t}`;
> hence `T_d <= tau` as soon as `(1+a0 tau)e^{-lambda tau} >= (3/2)M0/|omega_0(x0)|`, and with
> `|omega_0(x0)| = M0`, `lambda = 0`: **`T_d <= 1/(2 a0)`**, i.e. `M0 T_d <= M0/(2 a0)`.
>
> **Corollary (conditional).** For the tapered plateau with `a0 = kappa M L + c_a M` and
> `log Re_E = 2L + c_E`, `c2 = M0 T_d log Re_E -> 2 theta/kappa`; `theta = 1/2` under (H2),
> `log(3/2)` under (H2\*). `2` is the **infimum over the family as `delta, w -> 0`**, not attained;
> the ceiling for the datum actually run (`kappa = 0.48137`) is `2.0774`.
>
> **Unconditional (survives).** L1: `a(0,0,t) = int K omega^theta`, `K = -(3/8pi) r z/rho^5`, exact
> at every `t` (independently confirmed to `1.6e-11`); for the sharp shell `a(0) = (M/2)log(R/rho0)`
> **exactly**. L2: the equatorial plane is material and `omega^theta <= 0` in `{z>0}` for all `t`,
> so `a(0,0,t) >= 0` for all `t` — **and this sign structure extends to the equatorial plane, but to
> no point carrying vorticity.** L3 (sympy): `d log f/dt = 5a(3cos^2 phi - 1)` vanishes at
> `cos^2 phi = 1/3`.
>
> **Withdrawn.** (a) "(H2) is equivalent to the floor `F_a >= 1`" — one-way only, refuted by the
> seat's own `N3-Re0=1` tracers. (b) "`P1 >= 0` proved" — its premise `a >= 0` on `supp omega`
> fails on 56 % of the support. (c) "`beta(0) = a^2 + P1 + P2` exact" — Euler only; the NS identity
> is `beta(0) + nu V(0) = a^2 + P1 + P2 + P3` with `P3 < 0`, `P3/a^2 = -0.10 … -0.16` at `Re0 = 1`.
> (d) "(H2) at the origin reduces to ONE inequality on ONE integral" — three terms, two of unknown
> sign. (e) "`P2` is the same shape as BFG's defective step" — different homogeneity, different
> scale structure, and the implication runs the wrong way.
>
> **Unchanged.** The Taylor/Duhamel route is unavailable for this datum (`C <= 2/3` is required and
> `C ~ 2`). (ii) is NOT proved. Nothing here refutes Bradshaw–Farhat–Grujić Thm 10.

---

## 9. Remaining gap

`(H2)` — `beta + nu V >= 0` along a trajectory carrying `|omega_0| = M0` — is still the whole
problem, and after this audit it is **not** reduced to one scalar inequality: at the only points
where a reduction exists (the origin, the equatorial plane) the vorticity is zero, and at every
point where the vorticity is `M0` the strain functional is already a difference with a 14–25 %
negative part. What is genuinely established and worth carrying forward: I1–I3b; the `C <= 2/3`
obstruction to the Taylor route; Theorem D; the exact kernel `K` and the exact sharp-shell value
`a(0) = (M/2)log(R/rho0)`; and the equatorial-plane sign law, which is new and was not claimed.

## 10. Files

| file | what it does |
|---|---|
| `x1_gate_audit.py` / `x1_log.txt` / `x1_results.json` | the equivalence counterexample from the seat's own `d4_results.json`; the `F_r` FL-043 finding; the 25/36 vacuity count |
| `x1b_sweep.py` / `x1b_results.json` | 112 ODE runs: how much (H2) violation the gate certifies as HOLDS |
| `x2_P1_premise.py` / `x2_log.txt` / `x2_results.json` | `a < 0` on 56 % of `supp omega`; the negative part of `P1` |
| `x3_viscous_term.py` / `x3_log.txt` / `x3_results.json` | the dropped viscous term `P3`, at `Re0 = 400 … 1` |
| `x4_kernel_and_constants.py` / `x4_log.txt` / `x4_results.json` | independent filament Biot–Savart confirmation of `K`; sharp-shell `a(0) = (M/2)L`; the `c_E` / `kappa` audit |
| `x5_transfer.py`, `x5b_dense.py` (+ logs, json) | 28 000-source sign scan of the exact strain kernel; my failed refutation of the equatorial transfer |
| `x6_where_the_structure_lives.py` / `x6_log.txt` / `x6_results.json` | `omega == 0` on the equator; the positive/negative split of the strain functional at the seat's own tracer seed |
| `rcommon.py` | copy of the target's `dcommon.py` (datum + read-only import of `nsring.py`, sha256 `3c8e93b2…`) |
| `copy/` | copies of the target's `d1/d2/d5/d6` + `dcommon.py` for re-running; `copy/d1_rerun.txt` is the re-run log |

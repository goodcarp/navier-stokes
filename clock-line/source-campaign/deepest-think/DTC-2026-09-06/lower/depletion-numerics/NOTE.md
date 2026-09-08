# NOTE — lower seat `depletion-numerics` — DTC-2026-09-06

**Question (GENERATION_2_BRAINSTORM §B, "THE LOG IS A STOCK, NOT A FLOW").** Does the
multi-octave stack spend its `log(R/rho0)` of coherent strain on the FIRST doubling only, so
that `T_4 * M0` becomes independent of the octave count `N` while `T_{3/2} * M0` keeps falling
like `1/log`, and the coherence count `C(t) = max_x a / ||omega(t)||_inf` collapses to `O(1)`
within about one turnover?

**Answer at the strength the gate reached: the pre-registered KILL fires. DEPLETION IS KILLED
in its consequential half and CONFIRMED in its diagnostic half — and the two halves turn out
to be about two different strains.**

- `T_4 * M0` does **not** stall. It falls with `N` at least as fast as `T_{3/2} * M0` does:
  fitted `beta` (from `log T = alpha - beta log N`) is **1.452 for `T_4`** (`N = 4..7`) against
  **1.269 for `T_{3/2}`** (`N = 3..7`), and on the common `N`-set (`4,5,6,7`) the
  pre-registered variation factors are `V_{T4} = 2.252 > V_{T32} = 1.971`, i.e.
  `|V_{T4}/V_{T32} - 1| = 0.143 <= 0.30` — the KILL clause, verbatim. It fires on the own-set
  reading too (`0.233`), so the choice of `N`-set does not change it.
- `C(t)` **does** collapse, and its numerator `max_x a` does end up `N`-independent — but that
  is the wrong strain. `max_x a` sits **on the axis at the origin** at every recorded time and
  every `N`, while by the time `omega` has quadrupled its maximum has migrated to
  `s* = 9.4 ... 10.9`, deep inside the stack. The strain *at the maximum*, `a*`, is **not**
  depleted: at `T_4` it is 0.644 / 0.960 / 1.229 / 1.475 for `N` = 4 / 5 / 6 / 7 — still rising
  with the octave count, and very nearly linearly in it.

**The one-line physics.** The log is a stock *of the axis strain* and a flow *of the clock*:
the stack is spent from the inside out, the vorticity maximum is handed outward to successively
larger shells, and each of those still has the outer octaves above it to stretch it. The clock
never reads the axis.

Tier: receipts, single seat. Nothing is promoted. **Numerics falsify; they never prove.**
Every number below came out of a script in this folder that I wrote and ran (`analysis_log.txt`
is the raw capture); `SHA256SUMS` is recomputed by `shasum`, never typed. Gates are as written
in `PREREG.md` before the first run.

---

## 1. What was run

The solver is `nsring.py`, copied **byte-for-byte** from `sharp/viscous-numerics/nsring.py`
(sha256 `3c8e93b2da7c8f5437f74c3c39e67fb7b5ef6b29d92e01204e3176d817969f03`); their files were
not touched. The driver `run_dep.py` reproduces their `run_main.py` and changes only what
`PREREG.md §2` fixed in advance: stop level raised `2 M0 -> 8 M0`, horizon raised to
`T_max = 4.0/M0`, all four crossings recorded, `C(t)` recorded every step, and periodic
checkpointing so a long run can be harvested mid-flight. The `dt` cap deliberately keeps their
ORIGINAL horizon (`dt <= max(0.6, 4/N)/60`) so the early time-step sequence is identical to
theirs — which is what makes control **D1** real:

> **`T_{3/2}` and `T_2` reproduce the frozen sharp-seat values to `0.000e+00` relative
> (bit-for-bit) at every `N = 3, 4, 5, 6, 7`.** The driver is theirs; only the horizon is mine.

Exact system, no model and no closure: `eta = omega^theta/r`, `L5 psi1 = -eta`,
`Psi = r^2 psi1`, `u^r = r a` with `a = -d_z psi1`, `D_t eta = nu L5 eta`. Datum A (the smoothed
bang-bang shell over `N` octaves), `rho0 = M0 = 1`, `Re0 = M rho0^2/nu = 100`, `h = rho0/8`,
box `L = 1.5 R`, `z >= 0` half-domain. All times below are `t * M0`.

## 2. PRIMARY TABLE — datum A, `Re0 = 100`, `h = 1/8`, box `1.5 R`

| N | grid | `log Re_E` | `T_{3/2} M0` | `T_2 M0` | `T_4 M0` | `T_8 M0` | `T_4/T_{3/2}` | `C(0)` | peak `omega` | stop |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| 3 | 96² | 7.852 | 0.6667 | 1.1306 | 3.9070 † | — | 5.860 | 0.835 | 4.039 | horizon |
| 4 | 192² | 9.239 | 0.4477 | 0.7428 | **2.0197** | — | 4.511 | 1.112 | 5.675 | horizon |
| 5 | 384² | 10.625 | 0.3379 | 0.5548 | **1.4322** | — | 4.238 | 1.389 | 7.305 | horizon |
| 6 | 768² | 12.011 | 0.2716 | 0.4433 | **1.1043** | 2.7517 ‡ | 4.065 | 1.666 | 8.002 | reached 8 M0 |
| 7 | 1536² | 13.398 | 0.2272 | 0.3693 | **0.8969** | — ¶ | 3.948 | 1.944 | >= 4.111 | in flight |

† `N = 3`'s `T_4` is **excluded** on two independent pre-registered grounds: the argmax radius
at that crossing is `s* = 8.581`, `s*/L = 0.715 > 0.5` (**D5** box flag — the maximum has left
the datum's own outer edge, `log(R/s*) = -0.070`), and the crossing time `t = 3.907` is past
the point where the **D3** sign control has itself grown to `1.53 M0` (§5). It is recorded and
never used.

‡ `N = 6`'s `T_8` is likewise flagged: `s* = 57.19`, `s*/L = 0.596 > 0.5`. It is the only `T_8`
crossing in the whole sweep, so no `T_8` fit is reported at all (`PREREG §3` needs two).

¶ The `N = 7` run is a 1536² grid on a machine whose load average ran between 16 and **613**
during the session (other seats). In 900 steps and 1.6 h of wall clock it reached `t = 0.929`,
`||omega||_inf = 4.111` — enough for `T_{3/2}`, `T_2` and `T_4`. It was still running, short of
`8 M0`, when this note was closed, so its `T_8` is not measured and its row is the only one in
the table read from a mid-flight checkpoint. Its `T_4 = 0.8969` was **predicted as 0.8975 in
advance from `N = 4,5,6` alone** — see §8, written before the crossing landed.

**Crossing diagnostics** (`s*` = spherical radius carrying `sup|omega^theta|`, `a*` = the
strain there):

| N | level | `t` | `s*` | `s*/L` | `a*` | `C` there |
|--:|:--|--:|--:|--:|--:|--:|
| 4 | `T_{3/2}` | 0.4477 | 1.743 | 0.073 | 1.1181 | 0.9585 |
| 4 | `T_2` | 0.7428 | 2.581 | 0.108 | 1.0948 | 0.7075 |
| 4 | `T_4` | 2.0197 | 9.439 | 0.393 | 0.6439 | 0.2156 |
| 5 | `T_4` | 1.4322 | 10.314 | 0.215 | 0.9596 | 0.2758 |
| 6 | `T_4` | 1.1043 | 10.689 | 0.111 | 1.2293 | 0.3354 |
| 6 | `T_8` ‡ | 2.7517 | 57.188 | 0.596 | 0.4221 | 0.0917 |
| 7 | `T_2` | 0.3693 | 2.581 | 0.013 | 2.0803 | 1.1851 |
| 7 | `T_4` | 0.8969 | 10.939 | 0.057 | 1.4754 | 0.3947 |

**`a*` at `T_4` is 0.6439 → 0.9596 → 1.2293 → 1.4754 as `N` goes 4 → 5 → 6 → 7.** Whatever has been
depleted by the time the vorticity has quadrupled, it is not the strain that is doing the
stretching.

## 3. `C(t)` — the coherence count, and its numerator

`C(t) = max_x a(x,t) / ||omega(t)||_inf` at `t M0 = 0, 0.25, 0.5, 1, 1.5, 2`:

| N | 0 | 0.25 | 0.5 | 1.0 | 1.5 | 2.0 | first `t` with `C < 1.5` |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 3 | 0.8346 | 0.8476 | 0.8122 | 0.6035 | 0.4038 | 0.2873 | 0 (already below) |
| 4 | 1.1118 | 1.0780 | 0.9241 | 0.5175 | 0.3130 | 0.2190 | 0 (already below) |
| 5 | 1.3891 | 1.2731 | 0.9456 | 0.4400 | 0.2600 | 0.1770 | 0 (already below) |
| 6 | 1.6663 | 1.4256 | 0.9339 | 0.3851 | 0.2188 | 0.1458 | **0.2061** |
| 7 | 1.9436 | 1.5286 | 0.8672 | — | — | — | **0.2610** |

(`N = 7`'s row stops at `t = 0.929`, where that run stood when this note was closed.)

`C(0)` is exactly `a_max(0)` in these units, and it is linear in the octave count:
> `a_max(0) = 0.4000 log(R/rho0) + 0.0028`
against the first-order continuum value `0.5 log(R/rho0) + 0.2168` — the smoothed, axis-regular
discrete datum realises **80%** of the ideal bang-bang coefficient. Note that `C(0) < 1.5`
already at `N <= 5`, so the (S3) clause is informative only at `N = 6, 7`; both pass.

**The numerator is the informative half.** `a_max(t) = max_x a`:

| N | 0 | 0.25 | 0.5 | 1.0 | 1.5 | 2.0 |
|--:|--:|--:|--:|--:|--:|--:|
| 3 | 0.8346 | 0.9648 | 1.0861 | 1.1286 | 0.9743 | 0.8332 |
| 4 | 1.1118 | 1.3283 | 1.4607 | 1.2887 | 1.0369 | 0.8706 |
| 5 | 1.3891 | 1.6987 | 1.7703 | 1.3762 | 1.0740 | 0.8937 |
| 6 | 1.6663 | 2.0647 | 2.0066 | 1.4319 | 1.0988 | 0.9090 |
| 7 | 1.9436 | 2.4158 | 2.1822 | — | — | — |

Two facts, both against the naive reading of §B and both in favour of its diagnostic half:
1. `a_max` first **RISES** by 16–26% (peak near `t = 0.25`–`0.5`); it does not start spending.
2. By `t = 1.5` it has converged to `0.974 / 1.037 / 1.074 / 1.099` for `N = 3/4/5/6` — i.e.
   **`a_max(1.5)` is `N`-independent to 6% over a factor 8 in `R`**, and by `t = 2.0` to 4%.
   The axis strain *is* a stock and it *is* spent within about 1.5 turnovers, exactly as
   §B(ii) says.

And yet the clock does not care. `max_x a` is attained at the **origin on the axis**
(`s_amax = 0.000`) at every recorded time and every `N`, while `sup|omega^theta|` sits at
`s* ≈ 1.7` early and at `s* = 3.27 / 6.07 / 11.19 / 20.06` (for `N = 3/4/5/6`) by `t = 1.5`.
The two maxima are never in the same place, and `C(t)` falls mostly because its denominator
grows.

## 4. What the clock actually obeys (supplementary, NOT gated)

The sharp seat's invariant `T_{3/2} M log(R/s*) = 0.988 +- 0.018` is reproduced here
(0.9710, 0.9924, 0.9833, 0.9787, 0.9760 at `N = 3..7` — a 1.1% spread, an independent
reproduction of their headline) and **does not extend** to the higher levels:
`Q_{T4} = T_4 M log(R/s*)` is 1.066 / 1.622 / 1.976 / 2.206 at `N = 4/5/6/7` — it grows, it is
not an invariant.

What *does* extend is the `1/log Re_E` form itself. Fitting `1/(T_k M) = A_k log Re_E + B_k`:

| level | `A` | `1/A` | offset `b = -B/A` | max rel. residual | `N` used |
|:--|--:|--:|--:|--:|:--|
| `T_{3/2}` | 0.52306 | 1.9118 | 4.975 | 0.0033 | 3,4,5,6,7 |
| `T_2` | 0.32866 | 3.0427 | 5.150 | 0.0040 | 3,4,5,6,7 |
| `T_4` | 0.14908 | 6.7077 | 5.929 | 0.0034 | 4,5,6,7 |

so, to better than 0.4% at every point,
> `T_{3/2} = 1.912 / (M (log Re_E - 4.975))`,  `T_2 = 3.043 / (M (log Re_E - 5.150))`,
> `T_4 = 6.708 / (M (log Re_E - 5.929))`.

Read as a mean stretching rate `<a>_k := log(lambda_k)/T_k = alpha_k (log Re_E - b_k)`:
> `alpha = 0.21208` (`T_{3/2}`), `0.22781` (`T_2`), `0.20667` (`T_4`).

**The mean stretching rate that carries the vorticity from `M0` to `4 M0` is the same function
of `log Re_E` as the one that carries it from `M0` to `1.5 M0`: the two `alpha`s differ by
2.6%**, and the widest spread over the three levels is 10% (`T_2` is the high outlier).
Depletion, as a statement about the clock, needs `alpha` to fall sharply with the level. It
does not fall at all. (The `T_4` fit has four points and is not claimed as more than a fit.)

The same statement as an overhead over the frozen-initial-strain clock
`T_k^fr = log(lambda_k)/a_max(0)`:

| N | `T_{3/2}/T^fr` | `T_2/T^fr` | `T_4/T^fr` | `T_8/T^fr` |
|--:|--:|--:|--:|--:|
| 3 | 1.372 | 1.361 | (2.352 †) | — |
| 4 | 1.228 | 1.192 | 1.620 | — |
| 5 | 1.158 | 1.112 | 1.435 | — |
| 6 | 1.116 | 1.066 | 1.327 | (2.205 ‡) |
| 7 | 1.089 | 1.036 | 1.258 | — |

The quadrupling costs 26–62% more than holding the initial axis strain fixed, against 9–37% for
the first half-doubling — so there **is** a real, level-dependent depletion in the *constant*.
But every column **shrinks** as octaves are added. Nothing here becomes `N`-independent.

## 5. Controls (`PREREG.md §5`) — what fired and what did not

- **D1 continuity — PASSES exactly.** `T_{3/2}` and `T_2` reproduce the frozen sharp-seat
  values to `0.000e+00` relative at `N = 3,4,5,6,7`.
- **D2 maximum principle — PASSES.** `sup|eta|` overshoot above its `t = 0` value is
  `0.0e+00` to `2.2e-03` over every run in this folder (threshold `5e-3`); it is the
  non-monotonicity of the 3rd-order reconstruction, the same size the sharp seat measured.
- **D4 resolution — PASSES at `N = 4`.** `T_4 M0` = 2.0875 (`h = 1/6`), **2.0197** (`h = 1/8`),
  2.0485 (`h = 1/12`): a **3.4%** spread against a 10% threshold, and 1.4% between the two
  finest. `T_{3/2}` on the same three grids: 0.4552 / 0.4477 / 0.4476.
- **D5 box — PASSES.** `N = 4` with the box doubled (`lam = 3.0`, `L = 48`) gives
  `T_4 M0 = 2.0057` against 2.0197 at `lam = 1.5`: **0.7%**. The `s*/L > 0.5` flag fires only
  on `N = 3`'s `T_4` and `N = 6`'s `T_8`, both of which are excluded.
- **D3 sign control — FAILS at the extended horizon. This is the seat's real caveat and it is
  a correction to the sharp seat's C1 as well.** The sign-reversed `N = 5` datum, run to the
  full `T_max = 4.0`, dips to `0.9961` at `t = 0.135` and then creeps **up**: 1.005 (`t = 0.5`),
  1.026 (0.9), 1.065 (1.43), 1.120 (2.0), 1.256 (3.0), **1.499 (3.85)**, peak 1.583 at
  `t = 4.004`; it reaches `1.5 M0` at `t = 3.853`. The anti-datum grows too, slowly, at long
  times. The sharp seat's C1 ("the reversed sign must never reach `1.5 M`") held only because
  their horizon was `<= 0.8`; **it is not true on the horizon this seat needed**, and any
  future seat extending these runs must re-derive its own trust horizon.
  Applied here — the reversed-sign anti-signal at each crossing I use:

  | crossing | `t` | reversed-sign `\|\|omega\|\|` |
  |:--|--:|--:|
  | `N=4 T_{3/2}` | 0.4477 | 1.0030 |
  | `N=4 T_2` | 0.7428 | 1.0167 |
  | `N=4 T_4` | 2.0197 | **1.1226** |
  | `N=5 T_4` | 1.4322 | 1.0653 |
  | `N=6 T_4` | 1.1043 | 1.0364 |
  | `N=7 T_4` | 0.8969 | 1.0257 |
  | `N=3 T_4` (excluded) | 3.9070 | **1.5277** |

  Every `T_4` used is at least a factor 3.5 above the sign-blind background, but the `N = 4`
  point carries a 12% one. `N = 3`'s carries 53%, which is the second reason it is excluded.

## 6. Deviations from PREREG (L-09 — logged, thresholds not re-cut)

- **DEV-1 (`N = 7` horizon).** Machine load ran between 16 and 613 during the session; the
  1536² run took 1.6 h to reach `t = 0.929` and was **still running, short of `8 M0`**, when
  this note was closed. Its `T_{3/2}`, `T_2` and `T_4` are measured and used; its `T_8` is not.
  No INCONCLUSIVE-BY-REACH clause therefore attaches to the gate, which has all four
  `N = 4..7` `T_4` points. `N = 7`'s row is the only one read from a mid-flight checkpoint
  rather than a completed run.
- **DEV-2 (`T_8`).** Only `N = 6` reached `8 M0` (`T_8 = 2.7517`), and that crossing is box-
  flagged (`s*/L = 0.596`). `N = 3,4,5` peak at 4.04 / 5.68 / 7.31 and then viscosity turns the
  growth over, so `T_8` is **not measurable in this family at `Re0 = 100`** below `N = 6`. No
  `T_8` fit is reported. This is a reach limitation of the datum, not a result about depletion.
- **DEV-3 (`D4` at `N = 6`).** The pre-registered `N = 6, h = 1/6` resolution check was
  **cancelled** to leave CPU for the `N = 7` run on the loaded machine. `D4` therefore rests on
  the three-grid check at `N = 4` (3.4% spread) plus the `D5` box control. Recorded as not run.
- **DEV-4 (KILL clause read on a common `N`-set).** `PREREG §4` defines `V_k` "across
  `N = 3..7`" for both levels, but `T_{3/2}` is clean at more `N` than `T_4` is, so the own-set
  comparison is not like-for-like: `V_{T32}` would be computed over `N = 3..7` and `V_{T4}`
  over `N = 4,5,6,7`. Both readings are reported —
  own-set `|V_{T4}/V_{T32} - 1| = |2.2518/2.9346 - 1| = 0.233` and common-set (`N = 4..7`)
  `= |2.2518/1.9706 - 1| = 0.143`. **Both fire**, so the ambiguity this deviation was raised
  about turned out not to matter. (Before `N = 7`'s `T_4` landed the two readings disagreed —
  own-set 0.377, common-set 0.110 — and this seat had committed in writing to the common-set,
  like-for-like reading; the fourth point removed the disagreement.) On any common range
  `V_{T4} > V_{T32}`, which is *further* from depletion, not nearer, and the `beta` fits in §7
  do not depend on this choice at all.

## 7. VERDICT, applying `PREREG.md §4` verbatim

- **(S1)** `V_{T4} < 1.5`: `V_{T4} = 2.2518` (over `N = 4..7`) → **FALSE**.
- **(S2)** `V_{T32} > 2`: `V_{T32} = 2.9346` → TRUE.
- **(S3)** `C(t) < 1.5` by `t = 1.5/M0` for every `N`: first crossings at
  `t = 0, 0, 0, 0.2061, 0.2610` for `N = 3..7` → TRUE.
- **KILL** `|V_{T4}/V_{T32} - 1| <= 0.30`: `0.143` common-set (`N = 4..7`), `0.233` own-set
  → **TRUE both ways**.

> ## RECORDED VERDICT: **DEPLETION KILLED.**
> `T_4 M0` scales with the octave count at least as steeply as `T_{3/2} M0` does —
> `beta_{T4} = 1.452` (`N = 4..7`) against `beta_{T32} = 1.269` (`N = 3..7`) and
> `beta_{T2} = 1.319` (`N = 3..7`); the mean stretching rate to `4 M0` is the same function of
> `log Re_E` as the rate to `1.5 M0` to 2.6%; and the strain at the vorticity maximum at `T_4`
> still grows with `N` (0.644 → 0.960 → 1.229 → 1.475 at `N = 4..7`).
> §B's prediction (i) — "`T_{2^k} M0` → a constant independent of `N` for `k >= 2`" — is
> **refuted at this Reynolds number on this family**, over `N = 4,5,6,7` and a factor 64 in
> `Re_E`, with resolution (3.4%), box (0.7%) and sign controls all reported, and with the
> `N = 7` value predicted in advance to 0.07% from the `N = 4,5,6` fit (§8).
> §B's prediction (ii) — the coherence collapse — is **confirmed**, but of the **axis** strain,
> which is not the strain the clock reads.

**What this does NOT establish.** These are numerics: they falsify, they never prove. The
refutation is of §B's stated consequence, on this family, at `Re0 = 100`, over the octave range
that survived the box and sign controls, with four clean `T_4` points. It says nothing
about `Re -> infinity`, nothing about data outside the mollified plateau family, and nothing
against the frame's conjecture (ii): that seat's `T_{3/2}` clock is untouched here — it is
reproduced bit-for-bit and its `1/(log Re_E - b)` form is confirmed to 0.3% at every point.

**Consequence for the brainstorm's payoff claim.** §B proposed that BFG Theorem 10 might be
"FALSE for growth factors near 1 and TRUE (energy-free) for growth factor `>= M* ~ 4`". This
seat finds no numerical support for that resolution: the `4 M0` clock is as energy-dependent as
the `1.5 M0` clock, with the same `log Re_E` slope to 2.6%. On this evidence, if Theorem 10 is
wrong at `M = 1.5` it is wrong at `M = 4` too, and the two-clock reconciliation §B hoped for is
not available by this route.

## 8. The prediction, and its test — written before the crossing, settled after

`predict.py` (two independent routes, fitted on `N = 4,5,6` only, `predict_log.txt`):

| route | `T_4(N=7) M0` |
|:--|--:|
| A: `1/T_4 = A log Re_E + B` evaluated at `log Re_E = 13.3977` | 0.9009 |
| B: `T_4/T_{3/2}` linear in `1/N`, ratio(7) = 3.9356, times the measured `T_{3/2} = 0.2272` | 0.8941 |
| **prediction** | **0.8975** (routes differ by 0.007) |

Depletion instead requires `T_4(N=7) >= T_4(N=6) = 1.1043`: the two hypotheses were 23% apart
and the test was a single number. The prediction above was computed and written down while the
`N = 7` run stood at `t = 0.50`, `||omega|| = 2.52` — before it could cross.

> **MEASURED: `T_4(N=7) M0 = 0.8969.`** Predicted 0.8975 — **0.07% high**, and 19% below what
> depletion required. `predict.py` recomputes both from the stored runs, so the comparison is
> reproducible rather than asserted.

The `N = 7` process was left running past this point; it checkpoints every 100 steps to
`dep_n7_partial.json` and writes `dep_n7.json` on completion. The one number it can still add
is `T_8(N=7)`; together with `N = 6`'s box-flagged `T_8 = 2.7517` that would make two `T_8`
points and allow the same test one octave of growth higher. Re-running `python3 analyse.py` in
this folder picks it up automatically.

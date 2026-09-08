# NOTE — sharp seat `viscous-numerics` — DTC-2026-09-06

**Question.** Does the ACTUAL viscous axisymmetric no-swirl NS evolution of a dyadic N-ring
datum reach `(3/2) M` in a time `~ C/(M log Re_E)` — making Astra's continuation clock
`H0 = c/(M(1 + log_+ Re))` sharp up to constants — or does viscous decay of the inner rings /
deformation of the stack / the `O(M)` corrections prevent it?

**Answer, at the strength the gate reached.** The supposition **SURVIVES**. On the primary
datum the doubling time falls like `1/N` over seven octaves, and the product
`T32 * M * log Re_E` is monotone decreasing (8.69 -> 3.04) across a factor **4300 in Re_E**.
The clean invariant is not `T32 * log(R/rho0)` but
> **`T32 * M * log(R/s*) = 0.988 +- 0.018`**  (mean, sd over N = 2..7)
where `s*` is the radius actually carrying `sup|omega^theta|` at `T32`. That is a 1.9% spread
over a factor **64 in R** and **4300 in Re_E**. **Numerics falsify; they never prove.**

Tier: receipts, single seat (L-08). Nothing is promoted. Every number below came out of a
script in this folder that I wrote and ran; `SHA256SUMS` is recomputed, never typed.

---

## 1. What was actually solved

The exact no-swirl axisymmetric NS system in the estate's 5D lift, no model and no closure:
`eta = omega^theta/r`, `L5 = d_rr + (3/r)d_r + d_zz`, `L5 psi1 = -eta`, `Psi = r^2 psi1`,
`u^r = -(1/r) d_z Psi = r a` with `a = -d_z psi1`, `u^z = (1/r) d_r Psi`, and
`D_t eta = nu L5 eta`. Viscosity is carried, never modelled.

Discretisation follows the estate's validated persistence-toy (`psi1` on nodes, `eta` on cell
centres, face volume fluxes = differences of `Psi` so the discrete transport velocity is
EXACTLY divergence free, 3rd-order upwind-biased advection). The sparse-LU Poisson was
replaced by a separable **DST-I(z) + tridiagonal(r) direct solver**, which is what makes seven
octaves affordable: the `N = 7` run is a `1536 x 1536` grid, 373 steps, 1385 s.
`eta_0` is odd in `z` and the dynamics preserves that exactly, so the run is on `z >= 0`;
the half-domain solver was checked to agree with the full-domain solver to **3.6e-16**.

The mechanism is an identity, not an approximation: inviscidly `eta` is a material invariant,
so `omega^theta = r eta` grows **exactly** like the material radius, `d ln omega^theta/dt = a`.
Everything measured below is that identity plus what viscosity does to it.

## 2. Validations — each with a control that fires (`validate.py`, `validate.json`)

| check | result | control |
|---|---|---|
| **C3** exact Hill `psi1`, full domain | rel err 6.1e-3 / 1.7e-3 / **4.4e-4** at `h` = .08/.04/.02 | — |
| `a_nose` | .17061 / .18030 / .19446 -> Richardson **0.19918** | exact `1/5` |
| axis peak `/U` | 2.48320 / 2.49410 / **2.49832** | exact `5/2` |
| **C2** operator control (`3/r -> 1/r`) | **1.83 / 1.84 / 1.85 (i.e. 185% error)** at every `h` | fires |
| **C4** 5D diffusion law `d<|X|^2>/dt` | 0.20015 (`nu`=.02), 0.50036 (`nu`=.05) | exact `10 nu`, rel 7.6e-4 |
| half-domain vs full domain on odd data | **3.6e-16** | — |
| **max principle** `sup\|eta\|` non-increasing | exact to 1.5e-3 on the whole primary sweep (worst 4.3e-3, on the least viscous / coarsest run) | the 3rd-order reconstruction is not monotone; this is its size |

## 3. Datum (fixed in `PREREG.md` before any run)

**A (primary)** the brief's smoothed bang-bang shell, made axis-regular:
`eta_0 = -A (2z/s^2) Theta(s)`, so `omega^theta = -A sin(2 phi)` on `rho0 < s < R = 2^N rho0` —
`|omega^theta| <= A`, vanishing on the axis, maximal on the 45-degree cones at **every** octave.
`A` rescaled so the DISCRETE `sup|omega^theta|` at `t=0` is exactly `M = 1` (verified: 1.0 in
all 43 runs). `rho0 = M = 1`, so every time below is `t*M`.
**B (variant)** `N` discrete Gaussian dyadic ring-pairs on the same cone, self-similar cores.

Strain law at `t = 0` (`probe_cone.py`): for datum A, `a(s) = kappa log(R/s) + c` along the cone
with `kappa` = 0.437 / 0.415 / 0.411 at `N` = 3 / 5 / 6 and `c = -0.15`; for datum B,
`kappa` = 0.178 / 0.172 / 0.172. Both are linear in the octave count — the log law holds not
only at the axis (where G2's Biot-Savart got `kappa_0 = 1/2` for the pure bang-bang field) but
**on the field's own support**, which is what the mechanism needs.

## 4. PRIMARY RESULT — datum A, `Re0 = M rho0^2/nu = 100`, `h = rho0/8`, box `1.5R`

| N | R | grid | `Re_E` | `log Re_E` | `T32*M` | `T2*M` | `P = T32 M logRe_E` | `Q' = T32 M log(R/s*)` | octave of `s*` | `a*` at T32 | `sup\|eta\|` end |
|--:|--:|--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 2 | 32² | 152.5 | 5.027 | — (1.057 at t=4) | — | — | — | — | — | 0.604 |
| 2 | 4 | 48² | 640.1 | 6.462 | 1.3446 | — | 8.688 | 1.026 | 0.90 | 0.4058 | 0.691 |
| 3 | 8 | 96² | 2572 | 7.852 | 0.6667 | 1.1306 | 5.235 | 0.971 | 0.90 | 0.7623 | 0.741 |
| 4 | 16 | 192² | 1.029e4 | 9.239 | 0.4477 | 0.7428 | 4.136 | 0.992 | 0.80 | 1.1181 | 0.795 |
| 5 | 32 | 384² | 4.116e4 | 10.625 | 0.3379 | 0.5548 | 3.590 | 0.983 | 0.80 | 1.4504 | 0.839 |
| 6 | 64 | 768² | 1.646e5 | 12.011 | 0.2716 | 0.4433 | 3.263 | 0.979 | 0.80 | 1.7798 | 0.859 |
| 7 | 128 | 1536² | 6.585e5 | 13.398 | 0.2272 | 0.3693 | **3.044** | **0.976** | 0.80 | 2.1079 | 0.866 |

`E` is computed from the `t=0` velocity field, never assumed; `ell/R = 0.6337 +- 0.0015` across
the sweep, so `Re_E` is a measured quantity.

- **Fitted exponent** `log(T32*M) = 1.2004 - 1.4059 log N` over `N = 2..7`, i.e. **`beta = 1.406`**.
- `T2` obeys the same law: `T2*M*log(R/s*)` = 1.647, 1.647, 1.614, 1.597, 1.587 (its frozen-strain
  value is `(ln2/ln1.5) Q' = 1.69`).
- `T32*M*a*` = 0.546, 0.508, 0.501, 0.490, 0.483, **0.479** against the frozen-strain value
  `ln 1.5 = 0.4055` — i.e. the real evolution costs a **~18% overhead** over holding the strain
  fixed, and that overhead is *converging*, not growing, as octaves are added.
- **Which rings carry the maximum**: `s*` = 1.74–1.87, i.e. **octave 0.80–0.90 — the innermost
  ring, in every run**, unchanged from `N = 2` to `N = 7`. The stack does not hand the maximum
  outward at `Re0 = 100`.
- **Strain per octave measured DURING the run at `t = T32`** (not at `t=0`): increments
  0.3565, 0.3558, 0.3323, 0.3293, **0.3281** -> per e-fold 0.514, 0.513, 0.479, 0.475, **0.473**.
  Linear in `N`, flattening by `N = 6`. The stack's deformation costs ~8% of the log
  coefficient over `N = 2..7` and then stops costing.

### PRE-REGISTERED GATE, applied verbatim
`beta = 1.406 >= 0.7` **and** `P(N)` monotone decreasing with `max P = 8.688` — which exceeds
the stated bound 8 at the **single smallest** point `N = 2`, where `log(R/rho0) = 1.39` and
there is essentially no log to gain. So:
- the **KILL does not fire** (`beta` is 1.41, not `<= 0.3`; the largest `N` reached `1.5M` and
  `2M` while `N = 1` did not);
- the literal SURVIVE conjunction fails only on the `P <= 8` clause at `N = 2`.
**Recorded verdict: SURVIVE-WITH-CAVEAT.** `P` is bounded above, monotone decreasing, and its
asymptote is `2 Q' = 1.98`. I do not re-cut the threshold (L-05); see DEV-1.

## 5. The viscosity sweep — the sharpest thing this seat found

Fixed `N = 5` (`R = 32`), `nu` varied over a factor 1600 (`resweep`). This is the brief's own
spoiler — "viscous decay of the inner rings" — run directly.

| `Re0` | `nu` | `T32*M` | `s*` | `s*/sqrt(nu/M)` | `log(R/s*)` | `Q'` | `sup\|eta\|` end |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 1.0 | 0.9567 | 9.119 | 9.12 | 1.2553 | 1.201 | 0.070 |
| 2 | 0.5 | 0.6859 | 5.810 | 8.22 | 1.7061 | 1.170 | 0.147 |
| 4 | 0.25 | 0.5405 | 4.169 | 8.34 | 2.0380 | 1.102 | 0.260 |
| 16 | 0.0625 | 0.3911 | 2.259 | 9.03 | 2.6510 | 1.037 | 0.580 |
| 100 | 0.01 | 0.3379 | 1.743 | 17.4 | 2.9100 | 0.983 | 0.839 |
| 400 | 0.0025 | 0.3289 | 1.743 | 34.9 | 2.9100 | 0.957 | 0.895 |
| 1600 | 0.000625 | 0.3268 | 1.743 | 69.7 | 2.9100 | 0.951 | 0.911 |

Viscosity **does** destroy the inner rings — `sup|eta|` falls to 7% at `Re0 = 1` — and the
maximum migrates outward to `s* = 9.1`. But it migrates to exactly the scale it must:

> **`s* = 8.68 sqrt(nu/M)` (spread 8.2–9.1) whenever that exceeds the datum's inner cut**,
> and it pins to the datum's inner edge once `sqrt(nu/M)` drops below it.

The construction **self-selects** the brief's `rho0 >= sqrt(nu/M)` with a constant of ~8.7,
and the log it can still spend is (`closing.py`, fits the four viscosity-limited points to
`+-0.05` in the log)
> `log(R/s*) = (1/2) log Re_E - 1.705`,
so in the viscosity-limited regime the measured law is
> **`T32 * M * [ (1/2) log Re_E - 1.705 ] = 1.127`, i.e. `T32 = 2.255 / ( M ( log Re_E - 3.409 ) )`.**

That is the saturating form the brief asked about. `T32` changes by only **3.3%** between
`Re0 = 100` and `Re0 = 1600` (0.3379 -> 0.3268), so at the primary `Re0` the doubling time is
**not** viscosity-limited: it is set by the strain, and the strain is set by the octave count.

## 6. Variant B (discrete dyadic rings) — same law, different constants, and a real spoiler visible

`T32*M` = 1.7948 / 1.2039 / 0.9180 at `N` = 5 / 6 / 7 (`beta = 1.998`); `N <= 4` did not reach
`1.5M` even at `T = 3`. `Q' = T32*M*log(R/s*)` = 1.993, 2.200, 2.336 against the frozen-strain
prediction `ln(1.5)/kappa_B = 0.4055/0.172 = 2.36`. `P` = 16.4, 12.7, 11.0 — decreasing.
Here the max sits at **octave 2.33–2.40, NOT the innermost ring**: datum B's cores are
`sigma_k = 0.2 s_k`, so its inner rings are thin and viscously die (`sup|eta|` ends at 0.17–0.25).
**This is the brief's "viscous decay of the inner rings" spoiler actually firing** — and it
does not stop the log speed-up: the winner sits at a *fixed* octave index independent of `N`,
so it still sees `a ~ kappa log(R/s*) ~ kappa N ln 2`. Lowering `nu` (Re0 100 -> 400 at `N=5`)
moves `T32` 1.795 -> 1.554 and `s*` 5.27 -> 4.81, in the direction the mechanism predicts.

## 7. Controls, convergence, sensitivity

- **C1 sign control — FIRES.** The sign-reversed datum at `N = 5` (where `a < 0` at the initial
  argmax) reaches `sup|omega^theta| = 1.020` at `t = 0.80` and never `1.5`; the same datum with
  the correct sign reaches `2.0` at `t = 0.555`. Growth is not a property of the geometry, it is
  a property of the sign of the strain — as it must be, since `omega -> -omega` flips `a`.
- **C5 Lagrangian tracers** (N=4, to `t = 0.75`): material `eta` conserved to 8.0e-2 / 2.7e-2 /
  1.3e-2 / 3.9e-3 for seeds at `s = 1.5 / 2 / 3 / 6`, matching the `sup|eta|` decay of 0.795.
  The innermost tracer's `omega^theta` grew by 2.047 = (radial expansion 2.225) x (`eta` decay 0.920).
  The growth is radial expansion; nothing else.
- **C6 single octave.** `N = 1` (`R = 2`) never reaches `1.5M` in `t = 4` — it ends at 1.057.
  With no octaves there is no gain, exactly as `a(s) = 0.41 log(R/s) - 0.15 < 0` at `s = 1.5`, `R = 2`.
- **Grid convergence** (`T32*M`): `N=3`: .6736 (h=1/6) / **.6667** (1/8) / .6721 (1/12);
  `N=4`: .4552 / **.4477** / .4476; `N=5`: **.3379** (1/8) / .3362 (1/12);
  `N=6`: .2769 (1/6) / **.2716** (1/8). Converged to **0.02–2%**.
- **Domain size**: `N=3` at box `2.25R` gives `T32 = 0.6610` vs 0.6667 at `1.5R` — **0.9%**.
- **Viscosity at fixed N=4**: `Re0` 25 / 100 / 400 -> `T32` 0.5149 / 0.4477 / 0.4312.

## 8. What this establishes and what it does not

**Does.** On the actual viscous axisymmetric no-swirl NS evolution of this family, the time to
reach `(3/2)M` (and `2M`) is `Theta(1/(M log Re_E))` up to a constant, with the constant measured:
`T32 M log(R/s*) = 0.988 +- 0.018` across `N = 2..7`, and
`T32 M [(1/2) log Re_E - 1.705] = 1.127` in the viscosity-limited regime. None of the three
named spoilers stops it: viscous decay of the inner rings only *moves the winner outward to
`8.7 sqrt(nu/M)`*, deformation of the stack costs ~8% of the log coefficient and then stops,
and the `O(M)` corrections are a bounded additive constant in the log (`c = -0.15` for datum A).
Read against Astra's `H0 = c/(M(1+log_+ Re0))`: this family exhibits `sup|omega|` reaching
`(3/2)M0` on a time of the same shape, so **the logarithm in the record clock is not slack that
better estimates would remove — at least not on this family.**

**Does not.** (i) Numerics falsify only; nothing here is a proof, and a single seat's toy is not
a peer-confirmed object. (ii) This is the axisymmetric no-swirl sector, where `eta` obeys a
maximum principle and there is no blow-up — the family exhibits *fast doubling*, not *escape*;
it says the clock's log is sharp, not that the ladder closes. (iii) The measurement is of a
DOUBLING time, not of a continuation-time failure; the upper bound and this lower bound match
in shape and are both `O(1)`-in-the-constant, so "sharp up to constants" is the strongest
reading available, and the constants are not pinned against each other here.
(iv) `Re_E` reaches 6.6e5; the asymptotic regime is entered but not deeply.
(v) Prior art was NOT swept in this seat (no literature calls were made) — the Kim-Jeong /
Bang-Cheskidov attribution in the brief is carried unverified, and a literature seat is owed
before any of §4-§5 is cited anywhere.

## 9. Deviations (L-09)

- **DEV-1 (gate exhaustiveness).** `PREREG.md` §4 as written is not exhaustive: `beta >= 0.7`
  with `P > 8` falls through all three verdicts. I found this after a 3-point pilot
  (`results_quick.json`, `N = 2,3,4`) showed `P(2) = 8.69`. I added the label
  "SURVIVE-WITH-CAVEAT" for that cell **without moving any threshold**, and I record here that
  I saw `P(2) = 8.69` before adding it. The literal PREREG SURVIVE conjunction **fails**; the
  literal KILL **does not fire**. Both facts are reported above rather than resolved by
  re-cutting (L-05).
- **DEV-2 (horizon mis-sized for variant B).** `T_max = max(0.6, 4/N)` was set from datum A's
  strain constant. Datum B's `kappa` is 2.4x smaller, so at `T_max` none of `N = 2..7` had
  reached `1.5M` (`results_rings.json`, kept). Datum B was re-run at `T = 3` in a **separate**
  driver `run_extra.py` so that `run_main.py` stays frozen at the hash that produced the
  primary results. No threshold changed. Datum B is a variant, not the primary.
- **DEV-3 (cosmetic).** The tracer record stores `s0 * cos45` in the field labelled `s0`; the
  seeds were `s = 1.5, 2, 3, 6`. Affects a label in `results_controls.json`, no measurement.
- **DEV-4.** `N = 1` and (for datum B) `N <= 4` produce no `T32`, so the fit uses `N >= 2`
  (datum A) and `N >= 5` (datum B). Datum B's fit has only 3 points and is reported as such.
- **DEV-5.** The pre-registered `h = rho0/12` convergence check was run at `N = 3, 4, 5`;
  at `N = 6` only `h = rho0/6` vs `rho0/8` was affordable, and `N = 7` was run at `h = rho0/8`
  only. The `N = 7` row therefore carries the convergence estimate of its neighbours (~2%),
  not its own.

## 10. Files

`PREREG.md` (thresholds fixed before the first run) · `nsring.py` (grid, DST Poisson,
operators, data) · `validate.py` / `validate.json` (C2, C3, C4, half-domain) ·
`probe_kappa.py`, `probe_cone.py` (+ `.json`) (t=0 strain law) · `run_main.py` (frozen driver:
`main`, `n7`, `rings`, `controls`, `quick`) · `run_extra.py` (`ringsT`, `conv`, `resweep`) ·
`results_*.json` (full time series for every run) · `analyse.py` / `analysis_log.txt` /
`gates.json` (PREREG gates applied verbatim) · `closing.py` / `closing_log.txt` (the
`s* ~ sqrt(nu/M)` arithmetic) · `n7_log.txt` · `SHA256SUMS`.

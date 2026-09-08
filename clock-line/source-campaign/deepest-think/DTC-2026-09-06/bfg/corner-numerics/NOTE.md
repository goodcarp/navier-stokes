# NOTE — `bfg/corner-numerics` — DTC-2026-09-06

**Task.** Run the missing corner. The parent seat (`sharp/viscous-numerics`) swept the octave
count `N` at fixed `Re0 = M rho0^2/nu = 100`, along which `log Re_E - 2N log 2 = 3.69` is
constant to four digits — so `log Re_E` and `log(R/rho0)` are affinely degenerate and the sweep
cannot tell a logarithmic clock in `Re_E` from a clock in the octave count. This seat fixes `N`
and moves `nu`.

**Answer in one paragraph.** The pure-log law **is killed**, on all three pre-registered
clauses. `M T_d` is **not a function of `log Re_E` alone**: three runs at *identical*
`log Re_E = 13.3977` give `M T_d = 0.3268, 0.2659, 0.2272` (spread **36.5 %**) as
`log(R/rho0)` goes 3.466 → 4.159 → 4.852, against a grid uncertainty of **2.25 %**. Holding
`N = 6` and raising `Re0` from 1 to 400 (a factor 400 in `Re_E`, `Delta log Re_E = 5.99`) moves
`M T_d` only 0.5385 → 0.2659, a factor **2.03**, where the frozen law demands **3.57**. The
clock's real abscissa is `log(R/s*)`, `s*` the radius carrying the maximum, and `s*` sits at
`8.54 +- 1.01` viscous lengths `sqrt(nu/M)` until it pins to the datum's own inner cut. **The
log-free floor `A0 = 0.02` is untouched** — the smallest `M T_d` anywhere is 0.2232 — and I show
below that the floor clause as pre-registered **could not have fired at any reachable
resolution**, so this seat carries no evidence against BFG either way.
**Numerics falsify; they never prove.** Tier: receipts, single seat (L-08). Nothing promoted.

Provenance: I imported the parent seat's **frozen** `run()` (`run_main.py` sha256 `c131972a…`
over `nsring.py` sha256 `3c8e93b2…`) and changed nothing in it; only the job list and the output
directory are mine. Four of my runs duplicate runs they had already done and reproduce them
**bit-for-bit** (rel diff `0.00e+00` on `N=5,6,7` at `Re0=100`). Every number below came out of
a script in this folder that I wrote and ran.

---

## 1. What was run

Datum **A** (theirs): the smoothed bang-bang shell `omega^theta = -M sin 2 phi` on
`rho0 < |x| < R = 2^N rho0`, axis-regular, `rho0 = M = 1`, half-domain (`eta` odd in `z`),
`h = rho0/8`, box `1.5R`, `Tmax = 3.0`, exact axisymmetric no-swirl NS in the 5D lift
(`D_t eta = nu L5 eta`, no model, no closure). `T_d := T32` = first time `sup|omega^theta| = 1.5 M`.

`E` is computed from the `t = 0` velocity field and **does not depend on `nu`** (`build()` and
`solve()` never see it), so `log Re_E(N, Re0) = log Re_E(N, 100) + ln(Re0/100)` exactly. Checked:
max deviation over every production run **`1.78e-15`**. That is what makes the corner clean —
I can place runs at *identical* `Re_E` with *different* octave counts.

## 2. THE TABLE (all `h = rho0/8`; `M = rho0 = 1`, so `t` is `M t`)

| leg | N | `Re0` | `nu` | `sqrt(nu/M)` | `log Re_E` | `M T_d` | `M T_2` | `s*` | oct(`s*`) | `M T_d log(R/s*)` | `M T_d log Re_E` | `sup\|eta\|` end |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| B | 5 | 100 | 0.01 | 0.1000 | 10.6252 | 0.3379 | 0.5548 | 1.743 | 0.80 | 0.9833 | 3.5905 | 0.8395 |
| **A** | 6 | **1** | 1 | 1.0000 | 7.4063 | **0.5385** | 1.0346 | 7.871 | 2.98 | 1.1287 | 3.9886 | 0.1324 |
| **A** | 6 | **4** | 0.25 | 0.5000 | 8.7926 | **0.3812** | 0.6928 | 3.654 | 1.87 | 1.0913 | 3.3513 | 0.3448 |
| **A** | 6 | **16** | 0.0625 | 0.2500 | 10.1789 | **0.3011** | 0.5302 | 2.259 | 1.18 | 1.0071 | 3.0653 | 0.6404 |
| B | 6 | 25 | 0.04 | 0.2000 | 10.6252 | 0.2896 | 0.4924 | 2.138 | 1.10 | 0.9844 | 3.0770 | 0.7209 |
| **A** | 6 | **100** | 0.01 | 0.1000 | 12.0114 | **0.2716** | 0.4433 | 1.743 | 0.80 | 0.9787 | 3.2627 | 0.8586 |
| **A** | 6 | **400** | 0.0025 | 0.0500 | 13.3977 | **0.2659** | 0.4309 | 1.743 | 0.80 | 0.9581 | 3.5627 | 0.8988 |
| B | 7 | 6.25 | 0.16 | 0.4000 | 10.6252 | 0.2753 | 0.4855 | 2.895 | 1.53 | 1.0429 | 2.9247 | 0.5016 |
| C | 7 | 100 | 0.01 | 0.1000 | 13.3977 | 0.2272 | 0.3693 | 1.743 | 0.80 | 0.9760 | 3.0437 | 0.8656 |
| E | 7 | 400 | 0.0025 | 0.0500 | 14.7840 | 0.2232 | 0.3608 | 1.743 | 0.80 | 0.9591 | 3.3004 | 0.8962 |

Maximum principle (`D_t eta = nu L5 eta`): `sup|eta|` non-increasing in **every** run
(max end/start ratio 0.9373); no run voided.

## 3. Resolution (leg D, pre-registered)

| N | `Re0` | `h = rho0/6` | `h = rho0/8` | `h = rho0/12` | rel diff vs `1/8` |
|--:|--:|--:|--:|--:|--:|
| 6 | 400 | 0.27189 | **0.26592** | 0.26452 | **2.249 %** (1/6), 0.525 % (1/12) |
| 6 | 1 | 0.53855 | **0.53855** | — | **0.000 %** (4e-7 relative) |

The `Re0 = 1` pair agreeing to `4e-7` is fortuitous — its `s*` moves 7.945 → 7.871 and its `a*`
1.0777 → 1.0892 between the two grids, and its `T2` differs by `6.6e-4` — so I do **not** claim
that precision. The pre-registered rule takes the **maximum**, giving

>   `sigma_rel = 0.0225`,  tolerance `3 sigma_rel = 6.75 %`.

(Using the `Re0 = 1` pair instead would make the kill easier; the max is the conservative choice.)

## 4. PRE-REGISTERED VERDICT 1 — KILL-LOG: **FIRES** (all three clauses)

`H_log : M T_d = 1.8659/(log Re_E - 5.0730)`, constants frozen in `PREREG.md`, never refitted.

**K1 — leg-A residuals** (tolerance `+-6.75 %`):

| `Re0` | 1 | 4 | 16 | 100 | 400 |
|---|--:|--:|--:|--:|--:|
| measured | 0.5385 | 0.3812 | 0.3011 | 0.2716 | 0.2659 |
| `H_log` | 0.7997 | 0.5016 | 0.3654 | 0.2689 | 0.2241 |
| residual | **−32.66 %** | **−24.02 %** | **−17.59 %** | +1.01 % | **+18.64 %** |

Four of five exceed. The residual is monotone through zero at `Re0 = 100` — the signature of a
law fitted along a line and extrapolated off it, not of noise.

**K2 — the assigned "fails to decrease as the law predicts" clause.** Measured
`M T_d(Re0=1)/M T_d(Re0=400) = 2.0253`; `H_log` demands `3.5678`. **Off by 43.24 %**, i.e. 6.4×
the tolerance. `M T_d` *does* decrease — it decreases far too slowly.

**K3 — iso-`Re_E` families.** `H_log` (any law in `log Re_E` alone) demands these be equal:

| `log Re_E` | members `(N, Re0)` → `M T_d` | spread |
|---|---|--:|
| 10.6252 | (5,100) 0.3379 · (6,25) 0.2896 · (7,6.25) 0.2753 | **20.4 %** |
| 13.3977 | (5,1600)\* 0.3268 · (6,400) 0.2659 · (7,100) 0.2272 | **36.5 %** |

\*from the parent seat's `results_resweep.json`, read by my script. Every family is **monotone
decreasing in the octave count** at fixed `Re_E`. This clause needs no fitted constant at all:
a function of `log Re_E` cannot take three values at one `log Re_E`.

>   **KILL-LOG: FIRES.** `M T_d` is not a function of `log Re_E`, and the parent seat's fitted
>   `(c, d) = (1.8659, 5.0730)` is an artifact of its degenerate abscissa.

## 5. PRE-REGISTERED VERDICT 2 — KILL-FLOOR (`A0 = 0.02`): **does NOT fire**

Smallest `M T_d` anywhere in this seat: **0.2232** (`N = 7`, `Re0 = 400`), an order of magnitude
above the floor. Nothing here touches BFG's shape.

**And it never could have.** Under the measured `M T_d = 1.0538/log(R/s*)`, seeing `M T_d = 0.02`
needs `log(R/s*) = 52.7`; with `s* >= 8.54 sqrt(nu/M)` that needs
`log Re_E >= 2 x 52.7 + O(1) ≈ 105`, i.e. `Re_E >= 10^46`. The largest `Re_E` reached by this
seat or the parent is `6.6e5` (`log Re_E = 13.40`). **The KILL-FLOOR clause as pre-registered has
no power at any reachable resolution** — I record that as a finding about the test, not about the
world. Numerics cannot adjudicate the BFG conflict from this direction; only the proof step can.

## 6. What the corner says the clock actually is (descriptive; no threshold attached)

`H_oct` was *also* pre-registered with frozen constants (`s* = max(1.743, 8.68 sqrt(nu/M))`,
`Q' = 0.98` pinned / `1.10` viscous). Same discipline, same runs:

| law | max \|resid\| | rms resid |
|---|--:|--:|
| `H_log` (frozen) | **32.7 %** | **17.6 %** |
| `H_oct` (frozen) | 9.7 % | 5.0 % |

Free refits over the combined 19 runs (mine + the parent's, `h = rho0/8`):

| law | fit | rms | % of mean `T_d` |
|---|---|--:|--:|
| `c/(log Re_E - d)` | `c=2.1272, d=+4.1080` | 0.1256 | 27.1 % |
| `c/(log(R/rho0) - d)` | `c=1.3653, d=+0.3500` | 0.1584 | 34.2 % |
| `c/(w log Re_E + (1-w) log(R/rho0) - d)` | `w=0.62, c=1.4620, d=+3.4460` | 0.0358 | 7.7 % |
| **`c/log(R/s*)`** (one parameter) | **`c = 1.0538`** | **0.0406** | **8.8 %** |

`s*` is the measured location of the maximum, and it obeys

>   **`s* = 8.54 +- 1.01 sqrt(nu/M)`** on the 9 viscosity-limited runs (`Re0 <= 25`), and pins to
>   the datum's own inner cut (`s* = 1.74–1.87 rho0`) on the 10 runs with `Re0 > 25`.

`Q' = M T_d log(R/s*) = 1.0272 +- 0.0742` over all 19 (viscous branch `1.0849 +- 0.0693`,
pinned branch `0.9753 +- 0.0212`).

**The reframing.** Because `s* >= 8.54 sqrt(nu/M)` and `M R^2/nu = Re_E/(ell/R)^2` with the
measured `ell/R = 0.63389`,

>   `log(R/s*) <= (1/2) log Re_E - 1.6888`,

so the octave count is **capped** by `Re_E`, and the surviving statement is a **lower envelope**:

>   **`M T_d >= 2.1076 / ( log Re_E - 3.3775 )`**,  no run below it (0 of 19 by more than
>   `3 sigma_rel`), and **saturated to within ~5 %** by exactly the runs whose inner scale sits
>   at the viscous floor (`N=6 Re0=4/16/25`, `N=7 Re0=6.25`).

That is the same *shape* the parent seat conjectured and close to their own viscosity-limited
`closing.py` number `2.255/(log Re_E - 3.409)` — but the mechanism is different, and the
difference matters: **a logarithm in `Re_E` bounds the clock from below only through the octave
count it can pay for.** A datum whose inner cut sits far above `sqrt(nu/M)` (the parent's whole
primary sweep, `Re0 = 100`, `rho0 = 10 sqrt(nu/M)`) leaves Reynolds number on the table: at
`N = 6`, raising `Re0` 100 → 400 buys 2.1 %, and at `N = 7`, 100 → 400 buys 1.7 %. Speed is
bought with octaves, not with viscosity. Quantitatively: **octaves buy 2.18× more speed per unit
`log Re_E` than viscosity does.**

## 7. What this does and does not establish

**Does.** (i) The parent seat's `M T_d = 1.8659/(log Re_E - 5.0730)` is killed as a law: it fails
by up to 33 % off its own fitting line, and iso-`Re_E` families spread by 36.5 % against a 2.25 %
grid uncertainty. (ii) The clock's abscissa is `log(R/s*)`, with `s*` self-selecting the viscous
scale `8.54 sqrt(nu/M)`. (iii) A log-`Re_E` clock survives as a *lower envelope*
`M T_d >~ 2.11/(log Re_E - 3.38)`, tight to ~5 % on the critical line `rho0 ≍ sqrt(nu/M)` and
slack by up to 60 % off it. (iv) At fixed geometry `M T_d` **saturates** as `nu -> 0` — it does
not go to zero — so no version of this family exhibits fast doubling from small viscosity alone.

**Does not.** (i) Numerics falsify only; nothing here is a proof, and one seat's toy is not a
peer-confirmed object. (ii) This is the axisymmetric no-swirl sector — a maximum principle for
`eta`, no blow-up; a *doubling* time, not a continuation-time failure. (iii) **Nothing here bears
on BFG.** The floor test could not fire (§5); `M T_d = 0.02` is `Re_E ~ 10^46` away. The conflict
between Astra pass 8 and BFG Thm 8/10 is decided in the proof, not on a grid. (iv) The constant
`8.54` in `s* = 8.54 sqrt(nu/M)` has a `+-12 %` spread over 9 runs and is not a law. (v) Single
datum family (A); datum B was not re-run in the corner. (vi) `Re_E <= 6.6e5`: the asymptotic
regime is entered, not deeply. (vii) The refuted object is the parent seat's *fit*, not the
brief's mechanism — SYNTHESIS §1 (S3) already carried this run as "conjecture", and this seat
sharpens what the conjecture must say, it does not create or remove one.

## 8. Deviations (L-09)

- **DEV-1.** None of the thresholds moved. `sigma_rel` was defined before the runs as
  `max(0.02, measured)`; the measured value 0.0225 exceeded the floor, so 0.0225 was used, as
  written.
- **DEV-2.** `PREREG.md` §1 declared `N = 7` at the low-`Re0` end and `h = rho0/12` at `Re0 = 1`
  unaffordable. `N = 7` at `Re0 = 6.25` **was** affordable and was run (it is a leg-B member,
  declared in advance); `h = rho0/12` at `Re0 = 1` was not run, as declared.
- **DEV-3.** Leg E (`N = 7`, `Re0 = 400`), declared "if affordable", completed and is included.
- **DEV-4.** The K3 comparison at `log Re_E = 13.3977` uses the parent seat's `N=5, Re0=1600`
  run as a third member. That number is theirs, read by my script from their
  `results_resweep.json` (sha256 `f0cac209…`); it is not one of mine. The clause fires on my own
  two members alone (15.71 %).
- **DEV-5.** The machine was heavily loaded by other work during the runs (load average 37–678),
  which changes wall times only, not results — the solver is deterministic, as the four
  bit-for-bit reproductions of the parent seat's runs demonstrate.

## 9. Files

`PREREG.md` + `PREREG_SHA256.txt` (thresholds fixed and hashed before the first run) ·
`prereg_predict.py` / `prereg_predict_log.txt` / `prereg_predictions.json` (both laws'
predictions, computed before the runs) · `run_corner.py` (driver; imports the parent's frozen
`run()`) · `results_{a_re1,a_rest,b_n7,c_n7,d_conv}.json` (full time series for every run) ·
`logs/*.log` · `analyse_corner.py` / `analysis_log.txt` / `verdicts.json` (PREREG verdicts applied
verbatim) · `analyse2.py` / `analysis2_log.txt` (reproduction check, iso-`Re_E` test, free fits) ·
`analyse3.py` / `analysis3_log.txt` (classification, the envelope) · `inviscid_limit_log.txt` (§5b) ·
`inherited_SHA256SUMS.txt` (the parent seat's files I read) · `SHA256SUMS`.

## 5b. The `nu -> 0` limit at fixed datum (recorded because it cuts the other way)

| family | `Re0` sweep | `M T_d` | change over the last factor in `nu` |
|---|---|---|---|
| `N=5` (parent `resweep`) | 100 → 400 → 1600 | 0.3379 → 0.3289 → 0.3268 | **3.3 %** over ×16 |
| `N=6` (mine) | 100 → 400 | 0.2716 → 0.2659 | **2.1 %** over ×4 |
| `N=7` (mine) | 100 → 400 | 0.2272 → 0.2232 | **1.7 %** over ×4 |

For a **fixed datum shape**, `M T_d` has a **positive `nu -> 0` limit** — the inviscid value set
by the octave count. That is what a viscosity-uniform bound like BFG's Thm 8 would predict for
each fixed geometry, and this seat confirms it. Refuting BFG requires the *infimum over data* to
fall to zero, which needs `N -> infinity`, not `nu -> 0`. This seat's runs reach `M T_d = 0.2232`
at `N = 7`. **The two hypotheses are not separated by anything on this grid.**

# PREREG — `bfg/corner-numerics` — DTC-2026-09-06
### written and hashed BEFORE any production run (L-09: deviations logged, thresholds never re-cut)

Tier: receipts, single seat (L-08). Nothing here is promoted. **Numerics falsify, never prove.**

## 0. The gap this seat exists to close

The parent seat (`sharp/viscous-numerics`) swept the octave count `N` at **fixed**
`Re0 = M rho0^2 / nu = 100`. Along that sweep
`log Re_E - 2 N log 2 = 3.689, 3.693, 3.694, 3.694, 3.694, 3.694` — constant to four digits.
So `log Re_E` and `log(R/rho0)` are **affinely degenerate** on the entire primary sweep, and it
cannot tell a logarithmic clock in `Re_E` from a clock in the octave count. SYNTHESIS §1 (S3),
note 2: *"No run on disk has many octaves and an active viscous floor at once."*

This seat runs that corner: **fix `N`, move `nu`.**

## 1. Runs (fixed now; nothing added after the first launch)

Identical system, discretisation, datum and driver as the parent seat — I import their frozen
`run()` from `run_main.py` (sha256 `c131972a…`) over `nsring.py` (sha256 `3c8e93b2…`) and
change nothing in either file. Datum **A** = the smoothed bang-bang shell
`omega^theta = -M sin 2 phi` on `rho0 < |x| < R = 2^N rho0`, `rho0 = M = 1`, half-domain,
`h = rho0/8`, box `1.5R`, `Tmax = 3.0` (so nothing is horizon-truncated),
`T_d := T32` = first time `sup|omega^theta| = (3/2) M`.

| leg | runs | why |
|---|---|---|
| **A** (assigned) | `N=6`, `Re0 in {1, 4, 16, 100, 400}` | the viscous floor `sqrt(nu/M)` moves from `1.0 rho0` down to `0.05 rho0` at fixed `log(R/rho0) = 4.1589` |
| **B** (iso-`Re_E` ladder) | `N=6 Re0=25`, `N=7 Re0=6.25`, `N=5 Re0=100` | `E` (hence `ell`, hence `Re_E/nu`) does not depend on `nu` — `build()`/`solve()` never see it — so `log Re_E(N,Re0) = log Re_E(N,100) + ln(Re0/100)` EXACTLY, and these three runs sit at the SAME `log Re_E = 10.6252` with `log(R/rho0)` = 3.466 / 4.159 / 4.852 |
| **C** (iso-`Re_E` pair) | `N=7 Re0=100` (re-run) vs `N=6 Re0=400` (in leg A) | both at `log Re_E = 13.3977`, `log(R/rho0)` = 4.852 vs 4.159 |
| **D** (resolution) | `N=6 Re0=1` at `h = rho0/6`; `N=6 Re0=400` at `h = rho0/6` and `rho0/12` | the pre-registered convergence check, taken at BOTH ends of the viscosity sweep |
| **E** (range, if affordable) | `N=7 Re0=400` | pushes `log Re_E` to 14.78 |

`N = 7` at the low-`Re0` end and `h = rho0/12` at `Re0 = 1` are declared **unaffordable now**
(they are diffusion-`dt`-limited on a 1536² / 1152² grid); if they are not run that is not a
deviation, it is this line.

## 2. The two laws, with their constants FIXED BEFORE THE RUNS

**`H_log` — the pure logarithmic clock.** `M T_d = c / (log Re_E - d)` with
`(c, d) = (1.8659, 5.0730)`, the parent seat's fit to its six primary points. **Not refitted.**
Reproduced against those six points by `prereg_predict.py`: residuals
`+0.06, -0.69, -0.05, +0.55, +1.01, +1.36 %`.

**`H_oct` — the octave clock** (the alternative, indicative only, no KILL attached):
`M T_d = Q'/log(R/s*)`, `s* = max(1.743, 8.68 sqrt(nu/M))`, `Q' = 0.98` pinned / `1.10` on the
viscous branch — the parent seat's own measured `Q'` and `s* = 8.68 sqrt(nu/M)`.

**`H_floor` — a log-free floor.** `M T_d >= A0` with `A0 = 0.02`, regardless of `Re_E`.
This is the shape BFG's Thm 8/10 would force (`T >= 1/(c ||omega_0||_inf)`, `c` absolute).

Predicted `M T_d` (from `prereg_predict.py`, output in `prereg_predict_log.txt`):

| leg | N | Re0 | `log(R/rho0)` | `log Re_E` | `H_log` | `H_oct` | `H_floor` |
|---|--:|--:|--:|--:|--:|--:|--:|
| A | 6 | 1 | 4.1589 | 7.4063 | **0.7997** | 0.5506 | >= 0.02 |
| A | 6 | 4 | 4.1589 | 8.7926 | **0.5016** | 0.4088 | >= 0.02 |
| A | 6 | 16 | 4.1589 | 10.1789 | **0.3654** | 0.3250 | >= 0.02 |
| A | 6 | 100 | 4.1589 | 12.0114 | **0.2689** | 0.2720 | >= 0.02 |
| A | 6 | 400 | 4.1589 | 13.3977 | **0.2241** | 0.2720 | >= 0.02 |
| B | 6 | 25 | 4.1589 | 10.6252 | **0.3361** | 0.2720 | >= 0.02 |
| B | 7 | 6.25 | 4.8520 | 10.6252 | **0.3361** | 0.3049 | >= 0.02 |
| B | 5 | 100 | 3.4657 | 10.6252 | **0.3361** | 0.3368 | >= 0.02 |
| C | 7 | 100 | 4.8520 | 13.3977 | **0.2241** | 0.2281 | >= 0.02 |
| E | 7 | 400 | 4.8520 | 14.7840 | **0.1921** | 0.2281 | >= 0.02 |

## 3. Grid uncertainty — the denominator of every threshold

`sigma_rel :=` max over the leg-D convergence pairs of `|T32(h1) - T32(h2)| / T32(h = rho0/8)`,
**floored at `0.02`** (the parent seat measured 1.93 % at `N = 6`, `h = 1/6` vs `1/8`).
The pre-registered tolerance is `3 sigma_rel`, hence **>= 6 %**. If leg D returns a larger
`sigma_rel` the larger value is used (which is conservative *for* `H_log`, i.e. makes the
kill harder). The threshold is never cut downward.

## 4. PRE-REGISTERED VERDICTS

**KILL-LOG** — the pure-log law `H_log` is killed if ANY of:
- **K1** some leg-A run has `|MT_d(meas) - H_log| / H_log > 3 sigma_rel`;
- **K2** the measured ratio `MT_d(Re0=1)/MT_d(Re0=400)` at `N=6` differs from `H_log`'s
  ratio `0.7997/0.2241 = 3.568` by more than `3 sigma_rel` in relative terms
  (this is the assigned wording: `M T_d` "fails to decrease as the law predicts");
- **K3** two runs whose measured `log Re_E` agree to within `0.01` differ in `M T_d` by more
  than `3 sigma_rel` in relative terms (legs B and C).

If none fires, `H_log` **survives this corner** — which is a falsification test passed, not a proof.

**KILL-FLOOR** — the floor `A0 = 0.02` is killed if any run with converged resolution gives
`M T_d < 0.015`.

Both verdicts are reported verbatim whatever they say, including "the two kills disagree" or
"neither fires". No threshold moves after this file is hashed.

## 5. Also recorded for every run (fixed now)

`T32`, `T2`, `s*` (radius carrying `sup|omega^theta|` at `T32`) and its octave index,
`a*` at `T32`, the invariant `M T_d log(R/s*)`, `M T_d log Re_E`, `sup|eta|` end/start ratio
(**must be <= 1**: exact maximum principle for `D_t eta = nu L5 eta` — if it exceeds 1 by more
than the parent seat's worst measured overshoot 4.3e-3, that run is void), the measured
`log Re_E` against `log Re_E(N,100) + ln(Re0/100)` (must agree to 1e-9 — a check on the claim
that `E` is `nu`-independent), step count and wall time.

## 6. Discipline

Every number in `NOTE.md` comes from a script in this folder that I wrote and ran. I edit no
file outside this folder. `SHA256SUMS` is recomputed, never typed. Numerics falsify only.

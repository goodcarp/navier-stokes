# PREREG — lower seat `depletion-numerics` — DTC-2026-09-06
### written BEFORE any production run of this seat (L-09: deviations are logged, thresholds are never re-cut)

Tier: receipts, single seat. Nothing here is promoted. **Numerics falsify; they never prove.**

## 0. Hypothesis under test
GENERATION_2_BRAINSTORM §B, "THE LOG IS A STOCK, NOT A FLOW" (log-count depletion).
Claim: the multi-octave stack spends its `log(R/rho0)` of coherent strain on the FIRST
doubling only; after roughly one turnover the coherence is gone and the clock reverts to
`~1/M_new`. Two falsifiable consequences:
 (i) `T_4 * M0` becomes essentially independent of the octave count `N` while `T_{3/2} * M0`
     keeps falling like `1/log`;
 (ii) the **coherence count** `C(t) := max_x a(x,t) / ||omega(t)||_inf` (`a = u^r/r`) collapses
     from `~(1/2) log(R/rho0)` to `O(1)` within about one turnover, for every `N`.

## 1. System, code and datum — inherited verbatim, not re-derived
Solver: `nsring.py`, copied byte-for-byte from
`sharp/viscous-numerics/nsring.py` (sha256 `3c8e93b2da7c8f5437f74c3c39e67fb7b5ef6b29d92e01204e3176d817969f03`).
That seat's validations (C2 operator control, C3 exact Hill, C4 5D diffusion law, half-domain
vs full-domain 3.6e-16, `sup|eta|` maximum principle) are inherited and are NOT re-run here
except for the two cheap ones listed in §5. I do not edit their files.

Datum A (`datum_shell`), `rho0 = M = 1`, `R = 2^N`, `Re0 = M rho0^2/nu = 100` (`nu = 0.01`),
`h = rho0/8`, box `L = max(1.5 R, 4 rho0)`, `z >= 0` half-domain with `eta_0` odd in `z`,
amplitude rescaled so the DISCRETE `sup|omega^theta|(0) = M = 1` exactly. Same as the sharp
seat's primary sweep. `N in {3,4,5,6,7}`.

## 2. The ONLY changes to the run driver
(a) the stop level is raised from `2 M0` to `8 M0`, and the horizon from `max(0.6, 4/N)`
    to `T_max = 4.0 / M0` for every `N`;
(b) first-crossing times `T_{3/2}, T_2, T_4, T_8` are all recorded by linear interpolation
    between recorded steps, exactly as the sharp seat records `T_{3/2}` and `T_2`;
(c) `C(t) = max over the node grid of a(x,t)` divided by `||omega(t)||_inf` is recorded every
    step, together with `argmax` location of `a` and of `|omega|`;
(d) three early-stop rules, fixed now: stop when `8 M0` is crossed; stop when
    `||omega||_inf` has fallen to `<= 0.85` of its running maximum AND that running maximum
    is below `8 M0` (the run has peaked and will not reach the remaining levels); stop at a
    **wall-clock cap of 7200 s per run** (a cap hit is logged as a deviation, and the run is
    reported as "horizon-limited", never as a crossing).

## 3. Reported quantities
Table over `N = 3..7`: `Re_E`, `log Re_E`, `T_{3/2} M0`, `T_2 M0`, `T_4 M0`, `T_8 M0`,
the ratios `T_4/T_{3/2}` and `T_8/T_4`, `C(0)`, and `C(t)` sampled at `t M0 = 0, 0.25, 0.5,
1, 1.5, 2` (linear interpolation in `t`), plus `s*`/octave of the argmax at each crossing and
the box-proximity diagnostic `s*/L`.
Fitted exponents: least squares of `log(T_k M0) = alpha_k - beta_k log N` over those `N` in
3..7 that actually reached level `k`, for `k = 3/2, 2, 4, 8`.
Variation factor `V_k := max_N (T_k M0) / min_N (T_k M0)` over `N = 3..7`.

## 4. PRE-REGISTERED GATES (verbatim from the task brief; thresholds fixed now)
- **DEPLETION SURVIVES** iff all three hold:
  (S1) `V_4 < 1.5` (i.e. `T_4 M0` varies by less than a factor 1.5 across `N = 3..7`);
  (S2) `V_{3/2} > 2`;
  (S3) `C(t) < 1.5` for some `t <= 1.5/M0`, for every `N` in 3..7.
- **DEPLETION IS KILLED** iff `T_4 M0` scales like `T_{3/2} M0` across `N`, read as
  `|V_4 / V_{3/2} - 1| <= 0.30`.
- Otherwise **INCONCLUSIVE**. If any `N` in 3..7 fails to reach `4 M0` inside the horizon,
  (S1) cannot be evaluated on the full range and the verdict is **INCONCLUSIVE-BY-REACH**,
  reported with the sub-range that did reach it and with the reason (viscous peak vs horizon
  vs wall cap) named per run.
- Note recorded in advance: `C(0) ~ (1/2) N log 2 + O(1)`, so `C(0) < 1.5` already at `N = 3`
  and possibly `N = 4`; for those `N`, (S3) is passed trivially and is reported as such.
  The informative content of (S3) is at `N = 6, 7`, where `C(0)` should exceed 2.

## 5. Controls named in advance — each must fire
- **D1 continuity with the frozen sharp-seat run.** The `T_{3/2}` and `T_2` produced here at
  `N = 3..7` must reproduce the sharp seat's frozen values (`results_main.json`,
  `results_n7.json`) to better than 1% for `N = 3..6` (`N = 7` likewise). Any larger deviation
  means my driver differs from theirs in something that matters and every number is void.
- **D2 maximum principle.** `sup|eta|` must be non-increasing to <= 5e-3 relative over the
  whole extended run. It is an exact property of `D_t eta = nu L5 eta`.
- **D3 sign control.** `N = 5` with the sign-reversed datum must never reach `1.5 M0` over the
  full 4.0 horizon and its `||omega||_inf` must decrease.
- **D4 resolution.** `N = 4` at `h = 1/6` and `h = 1/12`, and `N = 6` at `h = 1/6`: `T_4 M0`
  must agree with the `h = 1/8` value to within 10%, or the `T_4` numbers are reported as
  resolution-limited.
- **D5 box.** `N = 4` at `lam = 3.0` (box doubled): `T_4 M0` within 10% of the `lam = 1.5`
  value. Additionally, any crossing whose argmax radius satisfies `s* > 0.5 L` is flagged
  CONTAMINATED and excluded from the fits (recorded, not deleted).

## 6. Discipline
Every number in NOTE.md comes from a script in this folder that I wrote and ran. `SHA256SUMS`
is recomputed by `shasum`, never typed. No threshold above is changed after the first
production run; deviations are logged in a DEVIATIONS section of NOTE.md.

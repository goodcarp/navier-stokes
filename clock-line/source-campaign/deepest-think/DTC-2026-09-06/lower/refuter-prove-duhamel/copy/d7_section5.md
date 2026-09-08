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

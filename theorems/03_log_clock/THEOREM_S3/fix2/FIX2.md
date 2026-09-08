# FIX2 — the five "beyond exactly eight items" issues, and the theorems/02 localisation slip

Seat `s3close/round2/THEOREM_S3/fix2`, sitting of 2026-09-08.
Laws: `TORMENT NEXUS/LAWS.md` (first 120 lines read).

Every number below came out of a script in **this** folder that I wrote and ran
(`f1_window.py`, `f2_datum_s3.py`, `f2_ramp_k2.py`, `f2b_measured.py`, `f4_scans.py`,
`f5_gamma_exist.py`); `SHA256SUMS` is computed with `shasum -a 256`, never typed.
Nothing outside `fix2/` was written: `THEOREM_S3.md`, its scripts, `hk2/`, `refute-u2/`,
`pmax-h3v/` and the external `alpoge-buckmaster-2026-09-08/` tree are all read-only here.

**Read** (L-14 declaration): `round2/THEOREM_S3/THEOREM_S3.md` in full and its four scripts;
`round2/THEOREM_S3/ADDENDUM_1_2026-09-08.md`; `round2/pmax-h3v/PROOF.md` §Task B (ϑ, f ≥ 4);
`hk2/PROOF.md` §0–§8 and `hk2/hk2lib.py`, `k2_datum.py`, `k3_bound.py`, `k4_direct.py`;
`round2/refute-u2/NOTE.md` §3–§4 and `r3_fixedpoint.py`; the other session's review at
`~/Desktop/forced-route-2026-09/docs/claude-scope-review.md` §1 and §2;
`campaign/external/alpoge-buckmaster-2026-09-08/THEOREM_forced_axisymmetric_on_axis_typeII.md`
§3.1–§3.2.

**Code imported** (byte copies in `imported/`, hashes in `SHA256SUMS`, each verified equal to
the hash in its source seat's own `SHA256SUMS`): `t1_results.json`, `t2_gamma_CR.py`,
`t3_budget.py`, `t2_results.json` from `THEOREM_S3/`; `hk2lib.py`, `k1_algebra.py`,
`k2_datum.py`, `k3_bound.py`, `k4_direct.py` from `hk2/`; `r3_fixedpoint.py` from `refute-u2/`.
L-14's inverse convention, declared: the objects under test (`t3_budget`'s viscous half and
`bootstrap`, `hk2`'s `k3` assembly) are **imported unchanged**, because the point is to test the
same object with the window and the datum corrected; everything I am testing a *claim* about
(`r_h`, `P_h`'s monotonicity, `E₀`, `𝔊₀`, `(Γ-off)`'s closure, the ramp slope, `K₂`'s datum
dependence) was re-implemented here from the mathematics.

---

## 0. STATUS TABLE

| # | issue | status | what changed |
|---|---|---|---|
| 1 | time-window coverage | **FIXED**, and it costs | the window is made self-consistent, `c ≥ c_*(1+ε(c,L))`; **the `ε ≤ 1/2` column is withdrawn as inadmissible** (§4 below); `L_*` (proved) moves `311511.9 → 1887775.3` |
| 2 | the ramp slope | **FIXED** | `ρΘ′` computed exactly; the claimed slope match is false; `K₂` recomputed for THEOREM_S3's own datum |
| 3 | the `λ` range | **FIXED** (with 2) | `K₂` re-run on `[1, λ_max]`, `λ_max` up to `2.0741`; the imported `161.7735` is a valid majorant, by `3.85×` |
| 4 | discrete scans labelled PROVED | **FIXED**: (a) ENCLOSED with certificate, (b) PROVED, (c) PROVED | and (b) yields the new structural cap `ε ≤ 0.19649` |
| 5 | `L_Γ* = 14260.5` | **FIXED** | renamed; the existence boundary is `14259.9`, a `0.0039 %` over-statement |
| 6 | Schwartz × indicator (theorems/02) | **FIXED** (text only) | corrected two-line argument in §6; the source file is **not** edited |

Direction of the whole sheet: **five of six corrections are conservative-neutral or better;
item 1 is a real loss.** The reach `L_*` gets worse by a factor `6.060` against THEOREM_S3's
headline, and the `ε ≤ 1/2` claim disappears entirely. Nothing here rescues (S3): the modulo
list of THEOREM_S3 §3 is untouched. **FL-000 stands.**

---

## 1. TIME-WINDOW COVERAGE — FIXED, and it is the expensive one

### 1.1 The defect, restated exactly

THEOREM_S3 §1.3 assumes every hypothesis and every bootstrap estimate **on `[0, τ]`**,
`τ := c_*/(ML)`, i.e. `θ = MLt ∈ [0, c_*]`, and concludes

```
    T_d(u_0)  <=  t_*  :=  c_*(1 + eps(L))/(M L)   =   tau (1 + eps)   >   tau .
```

`t3_budget.py` line 24 states the window as `theta in [0, c]` with `c = c_*`, and `eps_at()`
passes `CSTAR` at every call, so `λ_max = e^{3c_*/4}`, `r_h` on `[1, λ_max]`, `ell_loss`'s
`3 log λ_max`, the majorant ODE's integration range, `J_pow(c, ·)`, `σ_y`, `σ_z`, `c_G` and the
whole viscous budget are all evaluated on the **smaller** interval. The estimates are not
established on the interval on which the conclusion is drawn. `refute-u1` and `u1` inherit the
same slip; the review is right that it is a repair, not a rounding.

### 1.2 The repair

Run the budget on the window actually used, `θ ∈ [0, c]`, and impose

```
    (SC)      c  >=  c_* ( 1 + eps(c, L) )
```

with `eps(c,L)` the budget's own output **at that `c`**. Under the contradiction hypothesis
`T_d > c/(ML)` the range bound `M_s ≤ (3/2)M` holds on all of `[0, c/(ML))`, so every constant
in the chain is legitimate there; the clock then reaches `3/2` at `θ = c_*(1+ε(c,L))`, and (SC)
is exactly the statement that this happens inside the window where the estimates hold.
`ε(c,L)` is increasing in `c` (through `λ_max`, `μ`, `μ_J`, `c_G`), so the smallest admissible
window is the fixed point of (SC), and its `ε` is the honest `ε(L)`. `f1_window.py` solves it
by a scan-plus-bisection in `x := c/c_*`, re-reading `λ_max = e^{3c/4}`, `r_h` (f4's certified
lower bound at that `λ_max`), `ell_loss`, the ODE window and the viscous budget at each `c`.

`f1_window.py` also reports the cruder repair the review allows — majorise by a **fixed** larger
`c = c_*(1+ε_target)` and check the budget returns `ε ≤ ε_target` there.

**Control (frozen window).** With `c` frozen at `c_*` the instrument reproduces THEOREM_S3
§4.5's own column: `ε(1e6, proved) = 0.10621646` against `0.10622`;
`ε(1e5, measured) = 5.168244e-03` against `5.1682e-03`; `ε(1e6, measured) = 4.453336e-03`
against `4.4533e-03`; `ε(1e3, measured) = 0.1157152` against `0.11572`. So the only difference
below is the window, not the instrument.

### 1.3 A hard cap the enlarged window runs into (new)

A2/C6 need `P_h` increasing on `[1, λ_max]` for the quasimonotone comparison. §4(b) proves that
`P_h′` vanishes at

```
     lam_mono = 2.0769162        (P_h'(lam_mono) = -6.7e-16 by bisection)
```

and is negative beyond it. Since `λ_max = e^{3c/4}`, the window may not exceed

```
     c / c_*  <=  (4/3) log(lam_mono) / c_*  =  1.1964878        i.e.   eps <= 0.1964878 ,
     certified (by §4(b)'s certificate rather than by the exact root):  eps <= 0.1943662 .
```

**`ε ≤ 1/2` needs `λ_max = e^{3(1.5c_*)/4} = 2.499991 > λ_mono`.** On that window `P_h` is not
increasing, A2's monotonicity claim is false (not merely uncertified: the exact `P_h′` at
`λ = 2.4999908` is `-0.474818`), and the quasimonotone comparison of C6 has no hypothesis.
**THEOREM_S3's headline column `L_* = 311511.9 (eps <= 1/2)` is therefore withdrawn, not
merely re-priced.** The largest admissible target is `ε ≤ 0.194366`.

### 1.4 The corrected `ε(L)` and `L_*`

**`ε(L)`, self-consistent window, `f` minimised over `{0.05 … 32}`** (`∞` = no `c ≤ c_cap`
satisfies (SC) at that `L`, so the theorem says nothing there):

| column | `ε(1e3)` | `ε(1e4)` | `ε(1e5)` | `ε(1e6)` | `ε(1e7)` |
|---|---|---|---|---|---|
| all proved | `inf` | `inf` | `inf` | `inf` | `0.014141` |
| proved, `C_a` measured | `inf` | `inf` | `inf` | `inf` | `0.013835` |
| all measured | `inf` | `0.013023` | `5.1856e-03` | `4.4548e-03` | `4.3822e-03` |
| *(frozen window, THEOREM_S3 §4.5)* | *`0.11572` (meas.)* | *`0.012553` (meas.)* | *`5.1682e-03` (meas.)* | *`0.10622` (proved)* | |

The floor moves too: `ε → 4.3822e-03` as `L → ∞` in the measured column, against the frozen
window's `(1−2κ_δ)/(2κ_δ) = 4.3741737e-03`, because the self-consistent window settles at
`c = (1 + ε_floor)c_*` and carries the slightly larger `λ_max = 1.846949`.

**`L_*` and `log Λ_* = 2L_* + 2.9234859`, three columns, corrected.** The admissible targets
are `ε ≤ 0.194366` (the cap of §1.3, replacing the withdrawn `ε ≤ 1/2`) and `ε ≤ 0.1`.

| column | `L_*` (`ε ≤ 0.194366`) | `log Λ_*` | `L_*` (`ε ≤ 0.1`) | `log Λ_*` | THEOREM_S3, `ε ≤ 1/2` / `ε ≤ 0.1` |
|---|---|---|---|---|---|
| **all proved** | **`1887775.3`** | **`3775553.4`** | `1977293.7` | `3954590.4` | `311511.9` / `1056242.9` |
| proved, `C_a` measured | `1848675.9` | `3697354.6` | `1932164.9` | `3864332.8` | `302704.3` / `1025043.0` |
| all measured | `1532.97` | `3068.87` | `1691.84` | `3386.61` | `422.4` / `1119.3` |

Factors against THEOREM_S3: `6.0600` on the proved headline (`311511.9 → 1887775.3`, and the
same `6.0600` on `log Λ_*`, `623026.8 → 3775553.4`), `1.8720` on the proved `ε ≤ 0.1` column,
`3.6292` and `1.5115` on the two measured columns.

**With ADDENDUM_1's `f ≥ 4`** (which the theorem must carry, since `ϑ(f = 1) = 6.87` is vacuous):

| column | `L_*` (`ε ≤ 0.194366`) | `log Λ_*` | `L_*` (`ε ≤ 0.1`) | `log Λ_*` |
|---|---|---|---|---|
| **all proved** | **`1887786.7`** | **`3775576.3`** | `1977309.3` | `3954621.5` |
| proved, `C_a` measured | `1848687.4` | `3697377.6` | `1932180.6` | `3864364.0` |
| all measured | `1538.37` | `3079.67` | `1701.58` | `3406.09` |

The `f ≥ 4` constraint costs `6.07e-06` relative in the proved column (ADDENDUM_1 reported
`1.2e-05` on the frozen window) and `0.35 %` in the measured one.

**The cruder repair, for comparison.** Majorising by a *fixed* enlarged window
`c = c_*(1+ε_target)` and checking `ε(c,L) ≤ ε_target` gives, proved column,
`L_* = 1977297.9` at `ε ≤ 0.1` (`c = 1.1 c_*`, `λ_max = 1.958040`) and `2030134.3` at
`ε ≤ 0.194366` (`c = 1.194366 c_*`, `λ_max = 2.074226`). It is never better than the
self-consistent window and at the larger target it is `7.5 %` worse, because the fixed window
is larger than the one (SC) actually needs. Both are valid repairs; the self-consistent one is
the sharp one and is what the table above reports.


### 1.5 What the loss is made of

At `L = L_* = 1887775.3` the self-consistent window settles at `c = 1.1360550 c_*`
(`λ_max = 2.0016429`), and the whole loss is one term:

| quantity at `L = 1887775.3`, proved column, best `f = 0.1` | frozen `c = c_*` | self-consistent `c = 1.136 c_*` | ratio |
|---|---|---|---|
| `λ_max = e^{3c/4}` | `1.8420112` | `2.0016429` | `1.0867` |
| `r_h` on `[1, λ_max]` | `0.9186353` | `0.8661941` | `0.9429` |
| `C''` (`(Γ-off)`'s root) | `1464.5592` | `1999.7254` | `1.3654` |
| `C_R` (sup over all `φ`) | `853.0114` | `1182.6747` | `1.3865` |
| `μ` | `4.5168e-03` | `0.0100944` | `2.2349` |
| **`ε_T'`** | **`0.0476387`** | **`0.1158950`** | **`2.4328`** |
| `ε` | `0.0546318` | `0.1360550` | `2.4904` |

`ε_v = 3.04e-08` and `ε_ell = 1.66e-06` are invisible; `ε` is `ε_T'` plus the floor, and
`ε_T' = (2π/r_h)[3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8]` with
`μ(c)L = λ_max³(C_R + 2λ_om log λ_max)(λ_max² − 1)/(3/2)`. Enlarging the window by `13.6 %`
raises `λ_max` by `8.7 %`, and `λ_max³(λ_max²−1)` by `1.61`; `c_G = (λ + C''/L)c` rises with `c`,
so `Ĝ = (3π²/16)e^{pc_G}𝔊₀` and through it `C''` and `C_R` rise by another `1.37`; and `r_h`
falls by `5.7 %`. The product is `2.43` on `ε_T'`, and `L_*` must rise by about that factor to
bring `ε` back to target. **That is the whole of the item-1 loss**: it is a drive-side factor,
not an exponent — the same structural point THEOREM_S3 §4.6 makes about `C_R` — so it costs a
factor `6.06` and not `e^{...}`.

---

## 2 and 3. THE RAMP SLOPE, AND WHAT `K₂` ACTUALLY DEPENDS ON — FIXED

### 2.1 `ρΘ′`, exactly

For `Θ(ρ) = ½[tanh((u−ε_r)/ε_r) − tanh((u−L+ε_r)/ε_r)]`, `u = log(ρ/ρ₀)`, the chain rule gives
`ρ Θ′(ρ) = dΘ/du` **exactly** (sympy residual `0`, `f2_ramp_k2.py` §A), and

```
    dTheta/du = [ sech^2((u-eps_r)/eps_r) - sech^2((u-L+eps_r)/eps_r) ] / (2 eps_r) .
```

Its **maximum over `ρ`** is `1/(2ε_r) = 2`, attained at `u = ε_r`, i.e. at
`ρ = ρ₀ e^{ε_r} = 1.2840254 ρ₀` (the `L → ∞` derivative `d/du` of `dΘ/du` vanishes at `u = ε_r`;
sympy residual `0`). **At `ρ₀` itself**

```
    rho_0 Theta'(rho_0) = [ sech^2(1) - sech^2((L-eps_r)/eps_r) ] / (2 eps_r)
                        -> 2 sech^2(1) = 0.83994868          (0.8399486832280522 at L = 40)
```

which is `42.0 %` of the `2` claimed at `THEOREM_S3.md:103`. The reviewer's formula and its
`0.83995` are confirmed to eight digits.

`hk2`'s (D-B) ramp is `tanh` in `ρ`, not in `log ρ`:
`ρ₀Θ′_hk2(ρ₀) = ρ₀/(2w₀) − (exponentially small) = 1.99999999999813` at `w₀ = 0.25ρ₀`, and its
maximum of `ρΘ′_hk2` is `2.0306259` at `ρ = 1.030475 ρ₀`. So the two profiles agree in the
**maximum log-slope** (`2`) and disagree at `ρ₀` (`0.83995` against `2.0`). THEOREM_S3 §1.1's
sentence *"the inner slope agrees with `hk2` (D-B) to the digit"* compares the maximum of one
with the point value of the other. **The stated transfer justification is void.** (Nothing else
in the note rests on that sentence; what the transfer actually needs is §2.2.)

### 2.2 Which quantity `hk2`'s `K₂` depends on

Not a point slope. `hk2` §4–§5 make `K₂` depend on exactly three things:

1. **the shell total-variation density `𝒥`** of hypothesis (2.1),
   `|Dη|({ρ' < |x'| < ρ' + dρ'}) ≤ M ρ'² 𝒥 dρ'`, through Lemma 5.2's far/near integrals
   `Λ₄ = O(M𝒥/ρ)`, `Λ₅ = O(M𝒥/ρ²)`;
2. **the local suprema on `B(x,d)`**: `sup_B|∇η|`, `sup_B‖∇²η‖_F`, `sup_{∂B}|∇η|` (Cor. 4.3);
3. **the geometry**: `d`, the distance from the tracked point to the kink set (1.1).

The radial ramp enters (1) through the factor
`sup_ρ √((Θ_u − Θ)² W² + Θ² W_φ²)` and (2) only when the tracked point sits inside the ramp.
`hk2`'s reason for demanding a ramp at all is structural, not quantitative: a **sharp** radial
edge makes `η` jump, `a` a single-layer potential and `‖∇²a‖ ~ 1/dist` (`hk2` §8(c),
`K̂₂ f → 8`). Any mollification kills that; the *variable* the `tanh` is written in is
irrelevant to it. And under ADDENDUM_1's `f ≥ 4` the tracked point sits at `ρ_* = 5ρ₀`, where
both ramps are flat to `10^{-5}`, so (2) is the plateau's and the whole radial-ramp question
reduces to (1).

**`𝒥` for THEOREM_S3's own datum, computed here** (`f2_ramp_k2.py` §B, `hk2` Lemma 5.1 with the
ramp kept and the supremum over `ρ'` taken, as (2.1) requires):

| datum | `𝒥` | source |
|---|---|---|
| bare plateau, a.c. part | `39.478418 = 4π²` | this folder, `f6` — **not** `hk2` §8(a)'s `45.313074`; see the note below |
| (D-A) `δ = 7.5°`, sharp `sgn(z)` | `78.641147` | `hk2` §8(a) |
| (D-B) `δ = 7.5°`, `w = 0.20` | `65.625910` | `hk2` §8(a) |
| **THEOREM_S3's datum, `Θ ≡ 1` (plateau part)** | **`75.251104`** | this folder, `f2` §B |
| **THEOREM_S3's datum, uniform over `ρ'` (what (2.1) asks)** | **`96.249844`** | this folder, `f2` §B, attained at `u = L − 0.332`, in the **outer** ramp |

So THEOREM_S3's datum has a shell density `1.47×` `hk2`'s (D-B) and `1.22×` its (D-A). That
number was nowhere in the chain. It is not fatal — `Λ₄`, `Λ₅` are linear in `𝒥` — but it is the
quantity through which the datum swap actually propagates, and it moves the wrong way.

**Instrument control** (`f6_controls.py`): the same quadrature returns `39.162747` for `hk2`'s
(D-A) taper-only a.c. part against `hk2` §8(a)'s `39.162729` (relative `4.5e-07`) and
`65.625135` for its (D-B) against `65.625910` (relative `1.2e-05`, the finite-difference `W′`).
Those are the two rows that carry `hk2`'s bounds, and they reproduce.

*A small find in `hk2` while doing this.* Its §8(a) "bare plateau, no taper" reference row reads
`𝒥_ac = 2π²(√2 + arcsinh 1) = 45.313074`, computed in `k2_datum.py` from the integrand
`√(1 + cos²φ)/sin²φ`. Lemma 5.1's integrand — and `k2_datum.py`'s own `J_of_A`, which produces
the rows that are used — is `√(𝔥² + (1−t²)𝔥′²)`, which for the bare plateau `𝔥 = 1/sin φ` is
`1/sin²φ` exactly, because `|∇η_P| = M/r²` (`hk2` §8's own identity). That gives
`𝒥_ac = 2π²∫₀^π sin φ dφ = 4π² = 39.478418`. The two integrands are not the same function.
The row is a reference, used in no bound, so nothing downstream moves; recorded because the
number it displays (`45.313074`) is not the quantity its own Lemma 5.1 defines. (The numerical
coincidence that `4π²` is also the value of §8(a)'s equatorial-jump column is exactly that: the
jump term is `2|S³|𝔥(0⁺) = 4π²` too.)

### 2.3 `K₂` recomputed for THEOREM_S3's datum, over the full `λ` range (items 2 and 3)

`f2_datum_s3.py` implements THEOREM_S3 §1.1's field in the interface `hk2`'s instrument expects
(derivatives taken in `(ρ, φ)` and pushed to `(r,z)` by the exact chain rule; the self-test in
`f2_datum_s3.py` checks them against fourth-order central differences at six points under
`hk2/k2_datum.py` §1's own scaling convention, worst scaled relative error `1.24e-05`, and
`3.2e-09` or better at every point except the deep-plateau one where the second derivatives are
`O(10^{-5})` and the finite difference is itself at its noise floor -> `f2_datum_check.json`). `hk2`'s `k3_bound.lambda_bulk()`, `assemble()` and
`sups_on_ball_B()` are then called **unchanged**; only the datum object is swapped, and the
exclusion set is corrected: THEOREM_S3's datum has no radial edge and no `sgn(z)` jump, and its
kinks are the four cones `φ = δ`, `π/2 − δ_m`, `π/2 + δ_m`, `π − δ`, each transported by
`tan ↦ λ³ tan`.

```
 K2hat = K_2 rho_0/M  on the tube N_tau, THEOREM_S3's datum, phi_0 = 30 deg, L = 10
 ------------------------------------------------------------------------------------
   lambda      1.00     1.25     1.50    1.84201    2.07410       <- lam_max at c_*, and at c_cap
   f = 1     43.589   52.918   71.511   110.143    146.156
   f = 4     19.169   21.490   25.418    33.710     42.008        <- ADDENDUM_1's f >= 4
 ------------------------------------------------------------------------------------
 hk2, its own datum, f = 0:  66.662 (lam=1), 107.674 (1.25), 161.774 (1.5), no lam > 3/2
```

**Item 3 is settled two ways.** (i) The `λ` range: `t3_budget.py:44` imports `161.7735` as a
bound *for `λ ∈ [1, 3/2]`* while the clock needs `[1, λ_max]` with `λ_max = 1.8420112` (and, on
the corrected window of §1, up to `2.0741`). The table above supplies the missing range.
(ii) The datum: `hk2`'s constant is for `hk2`'s datum at `f = 0`; the table above is for
THEOREM_S3's datum at its own `f`. In both directions **the imported `161.7735` is a valid
majorant** for the theorem's datum at `f ≥ 4` over the whole corrected range, by a factor
`161.7735/42.008 = 3.851`. The budget's `K₂` input was therefore conservative, and the two
defects the review names are real but cost nothing: `(H-K2)` is worth `1.4e-06` of `L_*`
(THEOREM_S3 §4.6, unchanged here).

**Controls.**
`R5 (no `log`)`: quadrupling `L` from `10` to `40` moves `K̂₂` by `≤ 2.28 %` (worst cell
`f = 4, λ = 1.84201`; most cells `< 1 %`). A `log(R/ρ₀)` dependence would move it by a factor
`4`. The residual `2 %` is the discrete `d`-optimisation grid moving one node, not an `L`
dependence.
`R2 (quadrature)`: refining `(NS,NT,NC)` from `(250,100,100)` to `(450,180,180)` moves `K̂₂` by
`3.35e-04` relative.

**The measured value, and the slack** (`f2b_measured.py`). `hk2`'s `k4` zonal-harmonic
instrument is used with two declared adaptations: the Gegenbauer coefficients are integrated
with the profile's **kinks on panel boundaries** (a single Gauss rule, right for `hk2`'s
real-analytic (D-B), mis-resolves the high-`l` coefficients of a kinked profile, and those are
exactly what `∇²a` depends on), and `k4`'s PDE control is renormalised — at ADDENDUM_1's
`f ≥ 4` the tracked point is in the deep plateau where `∂_zη = O(10^{-5})` while `‖∇²a‖`
is `O(10^{-2})`, so `k4`'s `|Δ₅a − ∂_zη|/|∂_zη|` reports a `4e-4` absolute residual as `39`.

```
   lambda = 1, L = 10          K2hat measured   |grad a|   ||Hess a||   PDE resid / Hess scale
   f = 1, LMAX = 641              0.704969      0.346757    0.324368          1.06e-02
   f = 4, LMAX = 641              0.292697      0.108509    0.023130          1.39e-02
   control: a(series) vs a(independent 5-D kernel quadrature), relative:  3.6e-05 / 2.7e-05
   SLACK at lambda = 1:   f = 1 :  43.589 / 0.704969  =  61.83 x
                          f = 4 :  19.169 / 0.292697  =  65.49 x
   (hk2's own slack, its datum, f = 0, lambda = 1:  66.662 / 2.1903 = 30.4 x)
```

`K₂` measured is stable to five digits between `LMAX = 321` and `641`, and `a` agrees with a
completely independent 5-D kernel quadrature to `3e-05`; the PDE residual is still `1.1–1.4 %`
of the Hessian scale at `LMAX = 641`, so the measured `K₂` is reported to three digits and no
further. **`λ > 1` is not measurable with this instrument**: `η_λ = η₀ ∘ T_λ^{-1}` is not of the
separable form `G(ρ)H(t)` the zonal series needs. `hk2`'s own measured column has the same
limitation in a sharper form — `k4` evaluates the *untransported* `η₀` at the moved point, which
`hk2` §14 item 7 declares wrong for the bound — so `hk2`'s quoted "slack `30.0×` over the
window" compares a transported bound with an untransported measurement. Recorded, not repaired.

---

## 4. THE DISCRETE SCANS LABELLED PROVED — FIXED

The exact structure, derived here from scratch (`f4_scans.py`; six sympy residuals `0`).
Writing `D(v,λ) := 1 + (λ⁶−1)v²`, so that `Av² + B = λ^{-4}D`,

```
    P_h(lam) = 3 lam^10 int_0^1 H(v) v^2 D^{-5/2} dv ,        H(v) := h(arcsin v) ,
    d/dv [ v^3 D^{-3/2} ] = 3 v^2 D^{-5/2}      (EXACT)
```

so with `F_λ(v) := λ⁹v³D^{-3/2}` one has `F_λ(0) = 0`, `F_λ(1) = λ⁹(λ⁶)^{-3/2} = 1`, and

> **`dμ_λ := dF_λ` is a probability measure on `[0,1]`, and `Q(λ) := P_h(λ)/λ = E_{μ_λ}[H]`.**

Three things follow immediately.

* **`P_1(λ) = λ` exactly** (`H ≡ 1`): A2's control is now a theorem, not a coincidence.
  Numerically `|3λ⁹∫v²D^{-5/2}dv − 1| ≤ 2.4e-15` over `λ ∈ [1, 2.5]`.
* `∂_λ log F_λ = 9/λ − 9λ⁵v²/D = (9/λ)(1−v²)/D ≥ 0`, so `∂_λF_λ = 9λ⁸v³(1−v²)D^{-5/2} ≥ 0`:
  **`μ_λ` decreases stochastically as `λ` grows** — the strain moves the mass toward the axis,
  where `h` is small. That is the mechanism behind `r_h < 1`, stated for the first time.
* integrating by parts (both boundary terms vanish, `∂_λF_λ = 0` at `v = 0` and `v = 1`),

```
    Q'(lam) = - 9 lam^8 [ (1/delta) Aint(lam) - (1/delta_m) Bint(lam) ] ,
    Aint = int_0^{sin delta}   v^3 sqrt(1-v^2) D^{-5/2} dv ,
    Bint = int_{cos delta_m}^1 v^3 sqrt(1-v^2) D^{-5/2} dv ,
```

because `H′ = +1/(δ√(1−v²))` on `(0, sin δ)`, `0` on `(sin δ, cos δ_m)` and
`−1/(δ_m√(1−v²))` on `(cos δ_m, 1)`. Both integrals are majorised in **closed form** by
dropping `√(1−v²) ≤ 1`:

```
    int_0^V v^3 (1+K v^2)^{-5/2} dv = K^{-2}[ 2/3 + (1/3)S^{-3/2} - S^{-1/2} ] , S = 1+KV^2, K = lam^6-1
```

(sympy residual `0`; the `K → 0` limit is `V⁴/4`). Call them `Abar`, `Bbar`. `D` is increasing
in `λ` at every `v > 0`, so `Abar`, `Bbar` are **decreasing** in `λ` while `λ⁸`, `λ⁹` increase:
on any cell `[λ_i, λ_{i+1}]` a rigorous upper bound is `λ_{i+1}^{8}·(Abar(λ_i)/δ + Bbar(λ_i)/δ_m)`.
That monotone majorant is the whole certificate — no interval library is needed.

### (a) `r_h` — **ENCLOSED**, with a certificate

`|Q′| ≤ Lip(λ) := 9λ⁸[Abar/δ + Bbar/δ_m]`, a **proved** bound, evaluated cell by cell with the
monotone rule above. On `[1, λ_max]`, `λ_max = 1.8420112`, with `50000` cells
(`step = 1.6840e-05`, `max Lip = 0.3952580`):

```
     r_h  in  [ 0.9186364112 , 0.9186365334 ] ,      argmin at lam = lam_max
     THEOREM_S3's 61-point sample:  0.91863653      (inside the enclosure, at the upper end)
```

The sampled value was an upper bound on an infimum used as a favourable constant — the unsafe
direction — but it is unsafe by `1.2e-07`. Separately, `Q′ < 0` is **PROVED** on
`[1.3509082, λ_max]` by exhibiting `(1/δ)Aint > (1/δ_m)Bbar` there (`Bint ≤ Bbar`), so on that
stretch the infimum is at the right endpoint; and `Q(λ_max) = 0.91863653 < Q(1) = 0.99564488`.
Label: **ENCLOSED (certificate: monotone-majorant Lipschitz plus grid), argmin PROVED at
`λ_max`.**

Quadrature: Gauss–Legendre, `160` nodes per panel, panels split at the kinks `sin δ`, `cos δ_m`;
against `mpmath` at 40 digits the float `Q` agrees to `2.8e-10` at `λ = 1` and better above;
`Q′`'s closed form agrees with a central difference of `Q` to `2.4e-09`.

### (b) `P_h` increasing on `[1, λ_max]` — **PROVED**

```
    P_h'(lam) = Q + lam Q' = Q(lam) - 9 lam^9 [ (1/delta)Aint - (1/delta_m)Bint ]
              >= r_h - 9 lam^9 Abar(lam)/delta                (Bint >= 0, Aint <= Abar, Q >= r_h)
```

and `max_{[1,λ_max]} 9λ⁹Abar(λ)/δ = 0.5248602` (attained at `λ_max`), certified by the same
monotone subdivision, against `r_h ≥ 0.9186364112`. Hence

```
     min_{[1, lam_max]} P_h'  >=  0.3937762  >  0 .     PROVED.
```

THEOREM_S3's `0.01`-grid increment `0.0040871` corresponds to `P_h′ ≈ 0.409`; the exact
`P_h′(λ_max) = 0.397053`. **This is now a theorem, and it comes with the cap of §1.3**: the same
formula, solved for its zero, gives `λ_mono = 2.0769162`, and the certificate itself fails at
`c/c_* = 1.1943662`.

### (c) the profile suprema — **PROVED**

With `F(φ) := sgn(cos φ) g(φ)/sin φ`, `ρ|η₀|/M = |F|Θ` and
`ρ²|∇η₀|/M = √(F²(Θ − Θ_u)² + F′²Θ²)`. Write `Θ = σ₁ − σ₂` with `σ_i = (1 + tanh(·))/2 ∈ (0,1)`;
then, **exactly** (sympy residual `0`, `ε_r = 1/4`),

```
     Theta_u = 8[ sig1(1-sig1) - sig2(1-sig2) ] = 8 b (1-s) ,   b := sig1-sig2 = Theta,  s := sig1+sig2
     Theta - Theta_u = b (8s - 7) ,
     (rho^2|grad eta_0|/M)^2 = b^2 [ P (8s-7)^2 + S ] ,   P := F^2 ,  S := F'^2 ,
     0 <= b <= 1 ,   b <= s <= 2-b .
```

(`P` and `S` here are the two angular coefficients, not `P_h` and not `Q(λ)` of (a)–(b).)
`max_s (8s−7)² = (9−8b)²`, and `f(b) := b²[P(9−8b)² + S]` has
`f′(b) = 2b(2Py² − 9Py + S)`, `y := 9−8b ∈ [1,9]`. If `81P ≤ 8S` the quadratic has no real root
and `f` increases, so `f ≤ f(1) = P + S`; otherwise `f ≤ P·max_y(9−y)²y²/64 + S = 6.4072266 P + S`.
Taking the maximum over `φ` of each branch:

* on the bulk `δ < φ < π/2 − δ_m`, `F = 1/sin φ`, `F′ = −cos φ/sin²φ`, so `P + S = 1/sin⁴φ`,
  maximal at `φ = δ`: `1/sin⁴δ`; and `81P ≤ 8S` holds there (`4754.334 ≤ 27091.712`);
* everywhere the other branch applies, `6.4072P + S ≤ 183.914` (bulk, the `81P > 8S` stretch,
  which begins at `φ = 17.45°`), `≤ 376.187` (taper), `≤ 135.353` (equator) — all far below
  `1/sin⁴δ = 3445.159436`, and the `81P ≤ 8S` stretch of the bulk attains exactly
  `1/sin⁴δ` at `φ = δ` (`f6_controls.py`; `all_below_1_over_sin4_delta: true`).

Hence, **exactly**,

```
     E_0     = sup |F| Theta   =  tanh(2L-2)/sin delta   <=  1/sin delta   = 7.66129757554039
     Gfrak_0 = sup sqrt(...)   <=  1/sin^2 delta                           = 58.69548054098106
```

both attained (as suprema, `𝔊₀` essentially) at the taper kink `φ = δ` and its mirror, where the
one-sided `F′` jumps from `+0.3346697` to `−58.19333257 = −cos δ/sin²δ` and
`√(F² + F′²) = (1/sin²δ)√(sin²δ + cos²δ) = 1/sin²δ`. A `200000`-point scan of the certified
`φ`-bound returns `58.69548054` at `φ = 7.5°`, matching the closed form to `9e-10`.
THEOREM_S3's finite-difference grid values `7.6612976` and `58.695476` are the same numbers;
the second is `4.54e-06` **below** the true essential supremum
(relative `7.7e-08`), i.e. unsafe, and is replaced here by the closed form (this is ADDENDUM_1's `𝔊₀ = 1/sin²δ`, now with its proof). Every number
in §1 and §5 below uses the closed forms.

---

## 5. `L_Γ*` IS THE PICARD THRESHOLD, NOT THE EXISTENCE BOUNDARY — FIXED

### 5.1 What the two objects are

`t2_gamma_CR.gamma_off()` runs a damped simultaneous iteration in `(c_G, p)` from
`c_G = λc`, `p = 2`, and returns `None` — commented "no fixed point at that `L`" — when
`p c_G > 500`, `c_G > 60` or the iteration budget runs out. `t2.L_gamma_star()` bisects on that
`None`. So `14260.5` is the smallest `L` at which **that iteration converges**. The review and
`refute-u2` MINOR-1 are both right that this is not a nonexistence proof.

### 5.2 The direct analysis

In the variable `Γ̄ := λL + C″` (so `c_G = Γ̄c/L`), the closure is the scalar equation
`F(Γ̄) = Γ̄` with

```
   Chat_a  = (C_far + C_inner) lam + min( 2R_A lam/sigma_* + pi lam/8 , 2 R_A e^{c_G} E_0 )
   Ghat    = (3 pi^2/16) e^{p c_G} Gfrak_0 ,   A = (lam/2)L + Chat_a
   Gam_rad = max( A , A + 2Ghat + lam/2 , 2Chat_a + 2Ghat + lam/2 )
   p       = min( 3 , 1 + 2 Gam_rad/Gbar )        (an inner monotone fixed point, iterated up from p = 1)
   C''     = 2 Chat_a + Ghat + lam ,   F(Gbar) = lam L + C'' .
```

`F` is continuous and strictly increasing in `Γ̄`; `F(λL) = λL + C″ > λL`, so `F − id` starts
**positive**; and `F` grows like `exp(p c Γ̄/L)`, so `F − id → +∞`. Therefore a fixed point
exists **iff** `min_{Γ̄}(F − id) ≤ 0`, and a sign change is an intermediate-value-theorem
**proof** that a root exists. A supersolution is all the bootstrap needs — any `Γ̄` with
`F(Γ̄) ≤ Γ̄` gives `Γ(s) ≤ 2a(0,s) + C″M` with `C″ = Γ̄ − λL` — and the smallest root is the best
constant. `f5_gamma_exist.py` computes `min_{Γ̄}(F − id)` on a log grid and bisects on `L`.

**Control (L-14).** Driven with `u2`'s own datum norms (`E₀ = 7.66060196`, `𝔊₀ = 20.9197070`,
`c = 2log(3/2)`, `λ = 3/2`), an instrument written here from the algebra reproduces
`refute-u2` MINOR-1's directly-minimised boundaries and `u2` §8.2's constants:

| | this folder | `refute-u2` / `u2` | rel |
|---|---|---|---|
| existence boundary, `σ_* = 1/2` | `4904.7466` | `4904.7` | `9.5e-06` |
| existence boundary, `σ_* = 0` | `5273.0362` | `5273.0` | `6.9e-06` |
| `C″(10⁴)`, `σ_* = 1/2` | `582.2399924` | `582.24` | `1.3e-08` |
| `C″(10⁴)`, `σ_* = 0` | `801.0451149` | `801.05` | `6.1e-06` |

### 5.3 The answer for THEOREM_S3's datum

```
   c = c_*      L_Gamma^exist = 14259.948      L_Gamma^picard = 14260.507     over-statement 0.0039 %
   c = 1.05 c_*                16889.070                       16889.730
   c = 1.10 c_*                19959.907                       19960.391
   c = 1.194366 c_* (the cap)  27214.343                       27215.373
```

and the sign-change certificate at `c = c_*`:

```
   L = 0.9990 L_exist :  min (F - id) = +17.5229      no root
   L = 1.0000 L_exist :  min (F - id) = + 0.0200      tangency
   L = 1.0010 L_exist :  min (F - id) = -17.5005      root exists, by IVT
   L = 1.0100 L_exist :  root Gbar = 24893.144,  C'' = 3289.322,  c_G = 1.40773,  p = 2.36594
```

**Rename**: `L_Γ*` in THEOREM_S3 §3 M8, §4.3 and the theorem statement of §1.3 is
`L_Γ^picard = 14260.5`; the **existence boundary** is `L_Γ^exist = 14259.9`. The theorem's
hypothesis `L ≥ L_Γ*` should read `L ≥ L_Γ^exist = 14259.9`, and the over-statement it removes
is `0.0039 %` — a factor `485` smaller than `refute-u2`'s `1.89 %` for `u2`'s own datum, because THEOREM_S3's larger `𝔊₀` steepens `F` and pushes the tangency and the
Picard-divergence point together. **Direction: conservative.** The label was wrong; the number
was very nearly right, and the restriction is not binding anyway (`L_* = 1.89e+06 ≫ 1.43e+04`;
at the window `c = 1.1360550 c_*` that `L_*` actually sits at, `L_Γ^exist = 22487.19`).

For `L` above the boundary the direct root and the Picard limit agree to twelve digits
(`C″(10⁵) = 1537.3117911775` direct against `1537.3117911764` Picard;
`C″(10⁶) = 1467.9116771906` against `1467.9116771896`), so **no number in THEOREM_S3 §4.3
changes** — only the label and the threshold's meaning.

*A methodological note that belongs with this item.* My own `smallest_root()` — a coarse
geometric scan for the first `Γ̄` with `F(Γ̄) ≤ Γ̄` — reports "no root" at `1.001 L_exist`, where
`min(F − id) = −17.5` and a root demonstrably exists. That is the *same* failure mode as the
Picard threshold, in a different search. The reliable criterion is the sign of the minimum, and
the reported boundaries above use it.

---

## 6. THE LOCALISATION SLIP IN `theorems/02` — corrected text (the file is NOT edited)

`campaign/external/alpoge-buckmaster-2026-09-08/THEOREM_forced_axisymmetric_on_axis_typeII.md`,
Step S1's proof, lines ~175–178, reads: *"cut out `|y − x| <= |x|/4`, on which
`f 1_{|y-x|<=|x|/4}` is Schwartz with every seminorm `O_N(|x|^{-N})` by (5), and use
`||R_iR_j h||_{L^inf} <= C ||hhat||_{L^1} <= C ||h||_{W^{4,1}}`"*.

**The defect.** For a force that does not vanish on `∂B(x,|x|/4)`, the product `f·1_B` is
discontinuous across that sphere, hence not Schwartz and not even `C⁰`; `‖f 1_B‖_{W^{4,1}} = ∞`,
so the very norm the next inequality uses is infinite, and `\hat{f 1_B}` decays like `|ξ|^{-1}`
at best, so it is not in `L¹(R³)`. The step as written is unavailable. (The finding is
cosmetic in the sense that the conclusion (S1a) is true; it is the argument that is broken.)

**The corrected two lines.** Fix `|x| ≥ 1` and put `R := |x|/4`. Choose once and for all
`χ ∈ C_c^∞(R³)` with `0 ≤ χ ≤ 1`, `χ ≡ 1` on `B(0,1/2)`, `supp χ ⊂ B(0,1)`, and set
`χ_R(y) := χ((y − x)/R)`; split `f_j = h + (1 − χ_R)f_j` with `h := χ_R f_j(·,t)`.

> **Line 1 (the near piece).** `h ∈ C_c^∞(B(x,R)) ⊂ S`, and since `∂^α χ_R = O(R^{-|α|}) = O(1)`
> for `R ≥ 1/4` while every derivative of `f_j(·,t)` is `O_N((1+|y|)^{-N})(1+t)^{-K}` uniformly
> on `supp χ_R ⊂ {|y| ≥ 3|x|/4}` by (5),
> ```
>    ||h||_{W^{4,1}} <= C sum_{|a|+|b|<=4} R^{-|a|} ||d^b f_j(.,t)||_{L^1(B(x,R))}
>                    <= C_N (1+t)^{-K} |x|^{-N}     for every N
> ```
> (the volume factor `|B(x,R)| = O(|x|³)` is absorbed into `N`); hence, the multiplier
> `ξ_iξ_j/|ξ|²` being bounded by `1` and `(1+|ξ|)^{-4} ∈ L¹(R³)`,
> `||R_iR_j h||_{L^inf} <= ||hhat||_{L^1} <= C ||h||_{W^{4,1}} <= C_N (1+t)^{-K}|x|^{-N}`.

> **Line 2 (the far piece).** `(1 − χ_R)f_j` vanishes on `B(x, R/2)`, so on its support
> `|x − y| >= R/2 = |x|/8`: the kernel is non-singular there, **no principal value is needed**,
> and with `|k(w)| <= C_k|w|^{-3}` (`k` homogeneous of degree `−3`, smooth away from `0`)
> ```
>    | int k(x-y)(1-chi_R(y)) f_j(y,t) dy |  <=  C_k (|x|/8)^{-3} ||f_j(.,t)||_{L^1}
>                                            <=  512 C_k C (1+t)^{-K} |x|^{-3} .
> ```

Adding the two, and the `cδ` part of the kernel which is `O_N(|x|^{-N})` directly from (5),
gives `|R_iR_jf_j(x,t)| <= C(1+t)^{-K}(1+|x|)^{-3}` for `|x| >= 1`, which is (S1a).

Two remarks worth carrying into the file if it is ever revised. First, this replaces the whole
three-region split: the `χ_R` cut alone handles the singularity, and the far bound uses only
`f(·,t) ∈ L¹` with its `t`-decay, which (5) supplies — the `|y| ≤ |x|/2` versus `|y| > |x|/2`
split becomes unnecessary. Second, the note's own parenthetical (*"the earlier version said `k`
is integrable against the remaining decay, which is false"*) diagnosed the right problem and
then reached for an object — a Schwartz function cut by a sharp ball — that does not exist; the
smooth cutoff is what the diagnosis was asking for.

---

## GATE, FILES, AND WHAT THIS SHEET CANNOT CATCH

```
171 CHECKS, 171 PASS, 0 FAIL          (check_fix2.py)
```

`check_fix2.py` re-reads every displayed constant from the stored JSONs, checks that the string
actually appears in `FIX2.md` or `ADDENDUM_2_2026-09-08.md`, re-derives the ratios and factors
quoted between them rather than copying, and re-runs the `L_*` decomposition of §1.5 from
`f1_window`'s own `assemble_c` at the two windows.


| file | what it establishes |
|---|---|
| `f1_window.py` / `f1_results.json` | item 1: the self-consistent window; the corrected `ε(L)`, `L_*`, `log Λ_*`, three columns, with and without ADDENDUM_1's `f ≥ 4`; the frozen-window control against THEOREM_S3 §4.5 |
| `f2_datum_s3.py` / `f2_datum_check.json` | THEOREM_S3 §1.1's field in `hk2`'s interface: `Θ` in `log ρ`, the kinked angular profile, exact `(ρ,φ)` chain rule, the strained field, the four kink cones |
| `f2_ramp_k2.py` / `f2_results.json` | items 2 and 3: `ρΘ′` exactly; `𝒥` for the theorem's datum; `K₂` over `λ ∈ [1, 2.0741]` at `f = 1` and `f = 4`, `L = 10` and `40`; the `no-log` and quadrature controls |
| `f2b_measured.py` / `f2b_results.json` | the measured `K₂` at `λ = 1` with kink-aware Gegenbauer coefficients, two independent controls, and the slack |
| `f4_scans.py` / `f4_results.json` | item 4: the probability-measure identity, `P_1(λ) = λ`, `Q′`'s closed form, the certified `r_h`, the proof that `P_h` increases, `λ_mono`, and the closed forms for `E₀` and `𝔊₀` |
| `f6_controls.py` / `f6_results.json` | controls for the `𝒥` instrument against `hk2` §8(a); the three branch bounds of §4(c); `L_Γ^exist` at the window `L_*` sits at |
| `f5_gamma_exist.py` / `f5_results.json` | item 5: `(Γ-off)`'s closure analysed directly, the existence boundary, the sign-change certificate, the `u2`/`refute-u2` control |
| `imported/` | byte copies, with hashes, of every file imported from another seat |
| `check_fix2.py` | the gate: 171 checks over the two documents |
| `SHA256SUMS` | computed with `shasum -a 256`, never typed |

What this sheet **cannot** catch, stated so it is not mistaken for coverage. It repairs five
consistency and certification defects and one text slip. It does not touch THEOREM_S3 §3's
modulo list: (P-max), (H3-V), (V-strain), (λ-cap), (H5-Lip), (H\*), (Θ-tail) are exactly where
they were, and ADDENDUM_1's `f ≥ 4` constraint and refute-u2's MAJOR-1 and MAJOR-2 stand. It
does not check that `C_R` is the right constant, that `u2`'s P1 and P2 are the right
instruments, or that the modulo list is complete. The `K₂` measurement is a measurement:
numerics falsify, they never prove.

**FL-000 stands. Nothing here touches the headline problem.**

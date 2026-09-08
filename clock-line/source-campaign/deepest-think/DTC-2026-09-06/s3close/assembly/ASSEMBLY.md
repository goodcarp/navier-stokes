# ASSEMBLY of conjecture (S3) — the log clock as a theorem, and the two places it does not close

Seat `s3close/assembly`, 2026-09-08.  Continues DTC-2026-09-06.
Laws: `TORMENT NEXUS/LAWS.md` (first 120 lines read).
Every number below came out of a script in **this** folder that I wrote and ran
(`a1`…`a8`, re-asserted by `check_constants.py`); `SHA256SUMS` is computed with `shasum`,
never typed.  Nothing outside this folder was written.

`s3close/hk2/PROOF.md` did not exist when `a1`–`a7` were run, so `(H-K2)`'s constant is carried
**symbolically** as `C_K` in the derivation (§2.5) and the `a2` tables are computed at
`C_K = 1`.  That seat filed **before this note was finished**; §3.6 and `a8_hk2_plugin.py` plug
its constant in and report the moved budget.  Every statement below that depends on `C_K` is
therefore given twice: symbolically, and at the hk2 value.

Numbers taken from other seats are *statements of refereed theorems*, listed once in §0.3
with their file of origin, and never silently re-used as if derived here.  Where I could
re-derive a record constant cheaply I did, and the agreement is reported as a check, not as
an input.

---

## 0. Headline

**(S3) is NOT assembled into a theorem.**  What this note does produce:

1. **The theorem, stated completely**, with the datum, the window, the conversion to `𝒯(Λ)`,
   and `c₂ = 4 log(3/2) = 1.6218604` — §1.
2. **The bootstrap, written out and integrated** — a majorant ODE on the window whose
   solution is the flow-map error `μ(s)`, with every constant explicit (§2).  It is a real
   continuity argument and it does close **as an inequality**; what it does not do is close
   at any `L` that has a name.
3. **The one number that decides the whole thing**: the linearised strain-feedback exponent
   ```
   p c = [ 1.5 + C′/L + (3/2)(2κ)(15π/4)(9/4)/r_h ] · c   =  34.06     (proved constants)
                                                          =   6.63     (record's own measured
                                                                        Lemma T′ sensitivity)
   ```
   at the doubling window `c = c_* = log(3/2)/κ_δ = 0.8113826`.  Everything else in the
   budget is `O(1)/L`; this is `e^{34}`.
4. **Consequently `L_* = 4.53·10^14`** for `ε ≤ 1/2` on all-proved constants — `log Λ_* =
   9.05·10^14`.  Honest, and useless.
5. **The single missing piece worth the most.**  If the argument could be *restarted* at
   intermediate times (the reference map re-derived from the state, so the map error resets),
   `N` sub-windows cost `N(e^{pc/N}−1)` instead of `e^{pc}−1`.  At `N = 16`:
   **`L_*` falls from `4.5253·10^{14}` to `309.93`** (proved column, `ε ≤ 1/2`, `C_K = 1`) or
   **`366.96`** with `s3close/hk2`'s proved `C_K = 66.6622` — i.e. straight into the band the
   refuters already quote for the other pieces (`L ≥ 429 … 1435`).  The restart is blocked by
   **BLOCK 2** below, and closing BLOCK 2 is therefore worth ~12 orders of magnitude in `Λ_*` —
   more than every other improvement in the budget combined.
6. **Two named blocks that no theorem in the record covers** (§2.6): the off-manifold restart
   (BLOCK 2) and the gradient bound for a field that is not the exact strained plateau
   (BLOCK 3).  One more that is still open (BLOCK 5: `(H3)`).  **BLOCK 6 — `(H-K2)` — is now
   closed**: `s3close/hk2` filed while this note was being written, with `K₂ ≤ 66.6622 M/ρ₀`
   and no `L`, and with the sting that a *sharp radial edge* makes `K₂` infinite, which forces a
   change to the datum of §1.1 (BLOCK 6, §3.6).  Plugging its constant in moves `L_*` by `0 %`
   at `N = 1` and `+18 %` at `N = 16`.
7. **The BFG confrontation is correctly located and NOT triggered** (§4).  The family satisfies
   every hypothesis Bradshaw–Farhat–Grujić state; the clause that would be contradicted is
   Theorem 10's `T ≥ 1/(c(M)‖ω₀‖_∞)` with `c(M)` absolute.  Since (S3) is not proved, nothing
   is refuted.  **FL-000 stands.**

### 0.1 What is new here, beyond the four write-seats and their refuters

* The bootstrap is **written down and solved**, not gestured at.  The gap the campaign has
  called "GAP T's PDE half" turns out to have a *quantitative* answer with the record's own
  theorems — `μ(s) = O(c/L)` — and the obstruction is not that the estimate fails but that its
  constant is `e^{34}`.
* The **feedback exponent `pc`** is identified as the single reach-determining quantity, and
  the **restart** is identified and priced as the highest-value missing lemma.
* **`C_E(δ,δ_m,ε)` is derived**, not quoted: `C_E = 2π Σ_l 2 N_l Ĥ_l² D_l/(2l+3)` with the
  Green function `G_l(ρ,ρ′) = ρ_<^l ρ_>^{-(l+3)}/(2l+3)` and `D_l = 1/(5(l+4))` exactly in the
  sharp `L→∞` limit.  It reproduces the record's `0.172403978` to `4.0·10^{-8}` and extends it
  to the tapered, equatorially mollified and radially mollified datum.
* **`Lemma 2` is load-bearing here** — the first place in the campaign where it is.  The
  synthesis records it as "never used by any step".  In this assembly it is what licenses the
  *z-odd* constants of `far-near-kernel-lemma` (`0.291999`, `0.014754` instead of `0.756945`,
  `0.026093`) at every `t > 0`, not just at `t = 0`.
* A **correction to `prove-lagrangian` §2**: its `κ_δ(5°) = 0.4999917` should be
  **`0.4999169`** (§1.2).  Its `7.5°`, `15°`, `30°` entries reproduce here to 7–14 digits, and
  the deficit law `1 − 2κ_δ ≈ δ³/4` (δ in radians) puts `5°` at `1.66·10^{-4}`, not
  `1.66·10^{-5}`.  Digit transposition, not load-bearing.

### 0.2 Grade

**FILED.**  No Solid, no Major.  This is an assembly note: it closes nothing that was open and
proves no new theorem.  Its content is (a) the bootstrap written out with constants, (b) the
identification and pricing of the restart, (c) the derived `C_E`, (d) the honest budget.

### 0.3 Record inputs (statements consumed, with origin)

| # | statement used | origin | status there |
|---|---|---|---|
| R1 | **Lemma T′**: `\|a[η₀∘Φ^{-1}](0) − a_ref[η₀;λ]\| ≤ πMA[3μ(1+μ_J)/(2(1−μ)⁵)+3μ_J/8]`; Cor. 3 composition `μ̃_J = μ_J^Λ(1+2κ_s/L)+2κ_s/L`; `κ_s = √6−2` | `write/lemma-T-shell-dependent/PROOF.md` | PROVED, refuter: **STANDS** |
| R2 | **Theorem A / Corollary A1**: `\|H_l(λ)\| ≤ √(2/π)(2λM+V_λ)l^{-3/2}`, `Σ_{l≥3}\|H_l\|‖C_l‖_∞/(l(l+3)−4) ≤ 5.7395λM` | `write/L3v-and-gamma-bound/PROOF.md` | PROVED (refuter re-derived) |
| R3 | **Theorem Γ**: `Γ := sup_slab‖∇u‖ ≤ 2a(0,t) + C′M`, **`C′ ≤ 119.33`** (Cor A1 + computed `V_λ`), `151.15` (Cor A1 + proved `V`), `76.11` (odd-only Cor A1) | `write/refute-L3v-and-gamma-bound/PROOF.md` §2.4 | structure PROVED, **constant corrected** from `33.5` |
| R4 | **Theorem V.4 / Cor V.5** with the refuter's sizing: `E_hess` × **87.166**, `E_4` tabled, `E_tail` on the ball at `d/2` | `write/V-b-.../PROOF.md` + `write/refute-V-b-.../PROOF.md` | theorem PROVED, **sizing withdrawn and corrected** |
| R5 | **Theorem V.1**: `\|η₁−η₂\|(X(τ),τ) ≤ N·min(B1,B2)`, `B1 = 2n e^{−e^{−2c}q/(4n)}` | `gaps/gap-V-aronson/NOTE.md` | PROVED (refuter: STANDS) |
| R6 | **far/near kernel lemma**, Prop. 1 + Prop. 2: `\|a_far(x)−a_far(0)\| ≤ 0.291999 M` (z-odd), `\|a_inner\| ≤ 0.014754 M` (z-odd), `\|a_collar\| ≤ (3.999218/sin φ + π/8)M`; Consequence A: the `l=1` interior mode is exactly constant inside the shell, rate `≤ M/2` per e-fold | `rebuild/far-near-kernel-lemma/NOTE.md` | PROVED |
| R7 | **Lemma 1** and the taper deficit: `a_λ(0) = (M/2)P_h(λ)L`, `Φ_{h_δ}(λ) ≥ 1` on `λ∈[1,3/2]` for `δ ≤ 30°` | `lower/prove-lagrangian/NOTE.md` §2 | PROVED |
| R8 | **Lemma 2**: `η ≤ 0` on `{z>0}` for all `t`, all `ν ≥ 0` (symmetry + parabolic max principle) | `lower/prove-lagrangian/NOTE.md` §3 | PROVED |
| R9 | material-point strain offset `\|a(material) − (M/2)λL\| ≤ 0.069 M`, `L`-independent on `λ∈[1,3/2]` | `lower/SYNTHESIS.md` §1.3 (refuter-lagrangian-b) | **MEASURED, not proved** |
| R10 | 19 viscous runs `(M T_d, s*)` | `sharp/viscous-numerics/NOTE.md`, `bfg/corner-numerics/NOTE.md` | measurements |
| R11 | **BFG** Thm 8 (p.8), Thm 10 (p.11), arXiv:1704.05546v4 | `sharp/literature-short-time/txt/bfg.txt` | quoted verbatim in §4 |

---

## 1. THE THEOREM TO BE ASSEMBLED

### 1.1 The datum

Fix `M > 0`, `s > 0`, `δ ∈ (0,30°]`, `δ_m ∈ (0, 20°]`, `ε > 0`, `φ₀ ∈ (δ, π/2 − δ_m)`, `f > 0`,
and `L > 0`.  Put

```
    ν := M ρ₀²/s²   (equivalently ρ₀ = s √(ν/M)) ,        R := ρ₀ e^L ,
    ω₀^θ(ρ,φ) := − M · Θ_ε(ρ) · sgn(cos φ) · min(1, φ_ax/δ) · min(1, |φ − π/2|/δ_m) ,
    φ_ax := min(φ, π−φ) ,     Θ_ε(ρ) := min(1, log(ρ/ρ₀)/ε, log(R/ρ)/ε) ,
```
`ω₀ = ω₀^θ e_θ`, `u₀ = ` Biot–Savart of `ω₀`, axisymmetric with no swirl, odd in `z`.
`‖ω₀‖_∞ = M` exactly (the bulk value is attained).  The **tracked material point** is
`x_* = (ρ_*, φ₀)` with `ρ_* := (1+f)ρ₀`; the **inset angle** is
`ϑ_in := min(φ₀ − δ, π/2 − δ_m − φ₀) > 0`, and `f` is the **radial inset** in units of `ρ₀`.
`L := log(R/ρ₀)`.  `T_d(u₀) :=` the first time `‖ω(t)‖_∞ ≥ (3/2)M`.

Design point used for every number below unless stated: `δ = 7.5°`, `δ_m = 5°`, `φ₀ = 30°`,
`s = 1/sin δ = 7.661298` (the record's normalisation `ρ₀ sin δ = √(ν/M)`: the taper width at
the inner edge is exactly one dissipation length), `f` free.

**Trajectory margins** (`a1`, computed).  Under `(r,z) ↦ (λr, λ^{-2}z)` the tracked point's
polar angle obeys `tan φ = λ³ tan φ₀`, so at `λ = 3/2`, `φ₀ = 30°` it sits at
**`62.833°`** — `27.167°` from the equator, hence `22.16691°` clear of the equatorial
mollification layer at `δ_m = 5°`, and `55.333°` clear of the taper.  At `φ₀ = 40°` the
equatorial margin at `λ=3/2` is only `19.449°`; at `φ₀ = 20°` it is `39.148°` but the taper
margin at `λ=1` drops to `12.5°`.  `φ₀ = 30°` is the record's choice and is the one that keeps
both margins `> 20°`.

### 1.2 Strain constants (computed here; cross-checked against the record)

`κ_δ := ½P_{h_δ}(1)`, `P_h(λ) = 3∫₀¹h(arcsin v)v²(Av²+B)^{-5/2}dv`, `A = λ²−λ^{-4}`, `B = λ^{-4}`.

| `δ` | `κ_δ` (this seat) | record (`prove-lagrangian` §2) | agree |
|---|---|---|---|
| 0° | `0.5000000` | `1/2` | exact |
| 5° | **`0.4999171`** | `0.4999917` | **NO — see below** |
| 7.5° | `0.4997212305210886` | `0.4997212` / `0.4997212305210993` (Lemma-T seat) | `1.1e-14` |
| 15° | `0.4978077449714567` | `0.4978077` | `7` digits |
| 30° | `0.4836252349819347` | `0.4836252` | `7` digits |

The `5°` entry: the small-`δ` law is `1 − 2κ_δ ≈ δ³/4` (`δ` in radians), verified here
(`a6`, `1 − 2κ_δ` against `δ³/4`): `3.5854e-5` vs `3.5887e-5` (3°), `1.65722e-4` vs
`1.66143e-4` (5°), `5.57539e-4` vs `5.60733e-4` (7.5°), `1.31571e-3` vs `1.32914e-3` (10°) —
agreement to `0.1–1 %`, improving as `δ → 0`.  So `1 − 2κ_δ(5°) = 1.657·10^{-4}` and
`κ_δ(5°) = 0.4999171`.  The record's `0.4999917` is `1 − 2κ = 1.66·10^{-5}`, a factor 10 too
small — a digit transposition.  Not load-bearing anywhere (nothing uses `δ = 5°`), recorded
for the ledger.

`r_h := inf_{λ∈[1,3/2]} P_{h_δ}(λ)/λ = 0.981822` at `δ = 7.5°` — reproduces the Lemma-T
seat's `0.981822` exactly.  `Φ_{h_δ}(3/2) = P(3/2)/P(1) = 1.4980435 / 1.4913933 / 1.4735550 /
1.4444174 / 1.3596823 / 1.2583864 / 1.0673581` at `δ = 3/5/7.5/10/15/20/30°`, reproducing the
Lemma-T refuter's normalised table to `≤ 5·10^{-6}` and confirming that seat's Finding 1 (the
"unresolved `P_h(3/2)` cross-check" was a normalisation slip and is closed).

### 1.3 The statement

> **THEOREM (S3), to be assembled.**  There are `L_*(δ, δ_m, φ₀, s, f, C′, C_K) < ∞` and an
> explicit `ε = ε(L, c, δ, δ_m, φ₀, s, f, C′, C_K) → 0` as `L → ∞` at fixed `c`, such that for
> every `L ≥ L_*` the datum of §1.1 satisfies
> ```
>     T_d(u₀)  ≤  t_*  :=  2 log(3/2) (1+ε) / (M L) .
> ```

The window is `τ := c/(ML)` with `c = 2 log(3/2)(1+ε)`, and the argument runs **by
contradiction**: assume `T_d > τ`.  Then `‖ω(t)‖_∞ < (3/2)M` on `[0,τ)`, so the certified range
`λ(t) ≤ 3/2` of Theorem A, Theorem Γ, Lemma T′ Cor. 2 and Lemma 1 holds **automatically** on the
window — no extra hypothesis is needed to keep the constants in range.  That is the one piece of
free bookkeeping in the whole assembly, and it is worth stating: the hypothesis `λ ≤ 3/2` that
every record theorem carries is not an assumption here, it is the contradiction hypothesis.

The clock closes when
```
    κ_δ c (1 − ε_a)  +  log(1 − ε_v)  ≥  log(3/2) ,
```
`ε_a` the relative deficit of the strain at the tracked material point against `κ_δ M L`, `ε_v`
the relative viscous loss of `|ω^θ|` there.  Solving for `c` and writing `c = 2log(3/2)(1+ε)`,
```
    1 + ε  =  [ 1 + log(1/(1−ε_v))/log(3/2) ] / [ (1 − ε_a)(1 − ε_δ) ] ,
    ε_δ := 1 − 2κ_δ = 5.5754·10^{-4}  at δ = 7.5° .                                  (1.1)
```
`ε_δ` is the floor: `ε → (1−2κ_δ)/(2κ_δ) = 5.5785·10^{-4}` as `L → ∞`.  The doubling window at
`ε = 0` is `c_* := log(3/2)/κ_δ = 0.8113826`.

### 1.4 The conversion to `𝒯(Λ)` — the energy, derived

`𝒯(Λ) := inf{M₀T_d(u₀) : Re_E(u₀) ≤ Λ}`, `Re_E = E₀^{2/5}M₀^{1/5}/ν`, `E₀ = ‖u₀‖₂²`
(the convention of `lower/SYNTHESIS.md` §4; the brief's `E₀ = ½∫|u₀|²` shifts `log Re_E` by
`(2/5)log(1/2) = −0.2772589` and nothing else — it moves `log Λ_*` by that additive constant).

**Derivation** (`a1_datum_energy.py`; the identity and the Green function are derived there,
not quoted).  With `u = ∇×(ψ e_θ)`, `∫|u|² = ∫ψω^θ dx`; with `ψ = rψ₁`, `η = ω^θ/r`,
`dx₅ = 2π²ρ⁴ sin³φ dρ dφ`,
```
    ‖u₀‖₂²  =  (1/π) ∫_{ℝ⁵} ψ₁ η dx₅ ,        −Δ₅ψ₁ = η .
```
The radial operator `Ψ'' + (4/ρ)Ψ' − l(l+3)Ψ/ρ² = −f` has solutions `ρ^l`, `ρ^{-(l+3)}`,
`p W = −(2l+3)` constant, hence Green function `G_l(ρ,ρ′) = ρ_<^l ρ_>^{-(l+3)}/(2l+3)` against
the weight `ρ′³dρ′`.  Writing `x = log(ρ/R)` and `η = Θ_ε(ρ)ρ^{-1}Σ_l H_l C_l^{3/2}(t)`:
```
    ‖u₀‖₂² = 2π R⁵ M² Σ_{l odd} [ 2 N_l Ĥ_l² / (2l+3) ] D_l ,
    D_l = ∫∫_{x<y} Θ(x)Θ(y) e^{(l+4)x − (l−1)y} dx dy ,   Ĥ_l := H_l/M ,
```
and for the **sharp** shell, `L → ∞`, **`D_l = 1/(5(l+4))` exactly for every `l ≥ 1`**
(including the `l = 1` resonance, where the naive formula's `1/(l−1)` cancels against the
`e^{(l−1)L}` — checked symbolically and numerically, `a1` `D_l_num_vs_closed_maxrel =
1.7·10^{-16}`).  Therefore

> **`C_E := ‖u₀‖₂²/(M²R⁵) = 2π Σ_{l odd} 2 N_l Ĥ_l² / (5(2l+3)(l+4))`  (sharp, `L→∞`).**

| datum | `C_E` | `(2/5)log C_E` |
|---|---|---|
| bang–bang cap (`δ=δ_m=ε=0`) | **`0.17240398`** | `−0.70316592` |
| record (`prove-lagrangian` §4, estate `c6`) | `0.172403978` | `−0.7031659` |
| relative agreement | **`4.0·10^{-8}`** | — |
| `δ=7.5°`, `δ_m=0`, `ε=0` | `0.17217566` | `−0.70369602` |
| **`δ=7.5°`, `δ_m=5°`, `ε=0`** | **`0.17032563`** | **`−0.70801727`** |
| `δ=15°`, `δ_m=5°` | `0.16891202` | — |
| `δ=7.5°`, `δ_m=5°`, `ε=0.10` | `0.13269228` | — |
| `δ=7.5°`, `δ_m=5°`, `ε=0.25` | `0.09174383` | — |
| `δ=7.5°`, `δ_m=5°`, `ε=0.50` | `0.05118291` | — |

`97.1855 %` of `C_E` sits in the `l = 1` mode.  The `L`-dependence is `O(L e^{-5L})` (not
`O(ρ₀/R)`): at `L = 10` the sharp-shell `C_E` already agrees with its `L→∞` limit to
`4.728·10^{-9}`.  The **radial** mollification is *not* a small correction — the energy lives at
the outer scale, so an `ε`-e-fold outer ramp costs `22 %` at `ε = 0.1` and `70 %` at `ε = 0.5`.
Any use of `Re_E` must fix `ε` and use *that* `C_E`.  (This is the same class of error as the
synthesis's "second correction to the frame", `0.0779` in `log Re_E`; here it is up to `0.47`.)

Then, exactly,
```
    log Re_E = 2L + 2 log s + (2/5) log C_E ,        L = ½ log Re_E − log( s C_E^{1/5} ) ,
    2 log s + (2/5) log C_E = 3.3643455    at s = 7.661298, C_E = 0.17032563 ,          (1.2)
```
and `Re_E = C_E^{2/5} M R²/ν` exactly.  Substituting into §1.3,

> **`𝒯(Λ) ≤ c₂ (1 + ε(Λ)) / log Λ`, `c₂ = 4 log(3/2) = 1.6218604`, for `Λ ≥ Λ_*`,**
> `log Λ_* = 2L_* + 3.3643455`, `ε(Λ) → ε_δ/(1−ε_δ) = 5.58·10^{-4}` as `Λ → ∞`.

(`c₂ = 4log(3/2)` is the **conservative** constant — it uses only `Φ_{h_δ}(λ) ≥ 1`, which is
proved.  The accelerated `8(1−√(2/3)) = 1.4680274` needs the integro-ODE model, i.e. GAP T
closed *and* the model's `λ(σ,θ)` identified with the true stretch; nothing here supports it.
The synthesis's correction stands: under `β+νV ≥ 0` alone the Riccati route gives `2`, and
`4log(3/2)` is the strain-route constant.)

---

## 2. THE BOOTSTRAP — GAP T's PDE half, written out

### 2.1 The reference map and the field that generates it

`Λ_s(x) := T_{λ(|x|,s)} x`, `T_λ(y,z) = (λ y, λ^{-2}z)`, `λ(·,s) ∈ C¹([ρ₀,R];[1,3/2])` with
`|ρλ′/λ| ≤ κ_s/L`.  `Λ` **is** the flow map of a genuine velocity field: differentiating,
`∂_sΛ_s(x) = (λ̇/λ)(|x|,s)·(Y,−2Z)` at the image point `X = Λ_s(x)`, so
```
    v(X,s) = 𝔞(|Λ_s^{-1}X|, s) · (Y, −2Z) ,      𝔞 := λ̇/λ  evaluated at the LABEL radius.
```
`v` is a shell-labelled coaxial strain.  It is **not** divergence-free in `ℝ³`
(`J_Λ = λ²(1 + D·(1−3cos²φ))`, `D = dlogλ/dlogρ`), and Lemma T′ Cor. 3 is exactly the price of
that: `|J_Λ/λ² − 1| ≤ 2κ_s/L`.

### 2.2 The closed set `B(t)`

For `t ∈ [0,τ]`, `B(t)` is the conjunction of

* **(i)** `Γ(s) := ‖∇u(·,s)‖_{L^∞(𝒮_s)} ≤ 2a(0,s) + C′M ≤ Γ̄ := (3/2)ML + C′M` for `s ≤ t`,
  on the slab `𝒮_s`;
* **(ii)** `μ(s) := sup_x |Φ_s(x) − Λ_s(x)|/|Λ_s(x)| ≤ μ̄(s) < 1` and
  `|J_Φ − J_Λ| ≤ μ_J^Λ(s) J_Λ`;
* **(iii)** `a(X(s),s) ≥ κ_δ M (L − ℓ_loss)(1 − ε_{T′}(s)) − C_a λ(s) M`;
* **(iv)** `|ω^θ(X(s),s)| ≥ λ(s)(1 − ε_v(s)) M`.

`Φ_s` is the NS flow map; `X(s) = Φ_s(x_*)`; `λ(s) := exp∫₀^s a(X,·)`.

### 2.3 The driver, and the majorant ODE

`d/ds(Φ−Λ)(x) = [u(Φ)−u(Λ)] + [u−v](Λ_s(x))`.  The first bracket is `≤ Γ(s)|Φ−Λ|` by (i).
The second is the model error at the point `Λ_s(x)`, and it has exactly three pieces:

* **the `l = 1` rate error.**  By `far-near-kernel-lemma` Consequence A the `l = 1` interior
  mode of the velocity generated by the shells outside `ρ` is *exactly constant inside them*
  and equals the e-fold count; Lemma T′ prices precisely that integral (`a_ref` is the
  shear-free strain integral `(M/2)∫P_h(λ)dlogρ`).  Hence
  `|a_true − 𝔞| ≤ ε_{T′}·κ_δ M L`, with `ε_{T′}` from Lemma T′ at `(μ,μ_J)`.
* **the label-vs-current-radius error.**  `dlog|X|/ds = 𝒜(1−3cos²φ) ∈ [−2𝒜,𝒜]`, so the
  log-radius travels by at most `∫Γ ds = c_G`; and `|∂𝒜/∂logρ| = κ_δ M`.  Cost `κ_δ M c_G`.
* **the deviatoric remainder `O(Mρ)`,** bounded by `C_a λ M` — see §2.4.

With `|(r,−2z)| ≤ 2|X|`, `|Λ| ≤ λ|x| ≤ (3/2)|x|` and `|Λ| ≥ λ^{-2}|x| ≥ (4/9)|x|`, writing
`m := sup_x|Φ−Λ|/|x|` so that `μ ≤ (9/4)m`, and `θ := MLt`:

```
  dm/dθ    = (Γ̄/(ML)) m + (3/2)[ 2(κ_δ ε_{T′} + κ_δ c_G/L) + C_a λ/L ]                (2.1)
  dμ_J^Λ/dθ = 2[ κ_δ ε_{T′} + κ_δ c_G/L + C_a λ/L + (Γ̄/(ML)) μ ]                       (2.2)
  μ = (9/4) m ,   μ_J = μ_J^Λ(1+2κ_s/L) + 2κ_s/L ,
  ε_{T′} = (2π/r_h)[ 3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8 ] ,   m(0) = μ_J^Λ(0) = 0 ,
  Γ̄/(ML) = 3/2 + C′/L ,   c_G = (3/2 + C′/L) c .
```
(2.2) is the divergence Grönwall: `∂_s log J_Φ^{(5)} = 2a(Φ)`, `∂_s log J_Λ^{(5)} = 2𝔞 +
O(κ_s/L)`, and `|a(Φ) − a(Λ)| ≤ ‖∇a‖_∞|Φ−Λ| ≤ Γ̄μ`.

**This is the continuity argument.**  `B(t)` with the majorant `(μ̄, μ̄_J)` given by the
solution of (2.1)–(2.2) is closed: the right-hand sides are increasing in `(m, μ_J^Λ)`, the
solution of the majorant system is strictly larger than the true `(m, μ_J^Λ)` wherever the
latter is defined and `μ < 1` (standard comparison), and the set `{t : B(t) holds}` is closed,
open and non-empty in `[0,τ]`.  Nothing in the argument is a fixed point that has to contract;
it is an ODE comparison, and it therefore does **not** need the loop gain to be `< 1`.

### 2.4 `C_a` — the deviatoric remainder, non-perturbatively

This is the step where the record's *modal* machinery (Theorem A, Cor. A1, `(B.3)`) would fail,
because it is stated for the exact strained plateau.  Replace it by
`rebuild/far-near-kernel-lemma`, whose hypotheses are **only** `|ω^θ| ≤ M_t` on the shell and
`ω^θ` odd in `z` — both of which hold for the *true* field at every time (`|ω^θ| ≤ (3/2)M` by
the contradiction hypothesis; oddness by **Lemma 2**).  Then, at a point of polar angle `φ`,
```
  a(x,s) = [ e-fold count over shells outside 2|x| ] + err ,
  |err| ≤ C_a M_s ,   C_a = 0.291999 + 0.014754 + 3.999218/sin φ + π/8 .
```
At the worst angle over the window (`φ = φ₀ = 30°`, `sin φ = 1/2`): **`C_a = 8.697888`** —
a factor `126.06` above the measured material-point offset `0.069` of R9.
Measured values of the same three quantities are `0.0098…0.0141`, `0.0032…0.012` and
`0.058 (10°) / 0.167 (45°) / 0.609 (90°)` — the proved collar constant is lossy by `~50×`, and
it is the dominant term.  The **measured** alternative is R9's material-point offset
`C_a = 0.069`, which is not proved.

`ℓ_loss = log 2 + log(1+f) + 3 log(3/2) = 1.9095425 + log(1+f)` e-folds (one octave for the
far/near split, the radial inset, and the worst-case `λ`-motion of both the tracked point and
the outer edge).

### 2.5 `ε_v` — the viscous loss (Theorem V.4, refuter's sizing)

At the tracked point, in units `ρ₀ = M = 1`, `ν = 1/s²`, `ντ = c/(s²L)`,
`σ_y = ν∫₀^τλ^{-2}dt`, `σ_z = ν∫₀^τλ⁴dt`, `r₀ = (1+f)sinφ₀`,
`d = min(f, ρ_* sin(φ₀−δ), ρ_* sin(π/2−δ_m−φ₀))`, `R_- = r₀ − d`, `N = r₀/sin δ`:
```
  ε_bulk = σ_y/r₀²                                            (Cor. V.5, exact)
  E_hess = ½K₂e^{c_G}τV₁ · [ r₀/R_-² + 4N/d ] ,  V₁ = 5ντ(e^{2c_G}−1)/c_G ,  K₂ = C_K M/ρ₀
  E_4    = (r₀/R_-⁵)[ (tr C_τ)² + 2 tr(C_τ²) ] ,  C_τ = diag(2σ_y I₄, 2σ_z)
  E_tail = 2N𝒫 + 𝒫 + Σ_{k=1}^{3} r₀^{-k}(E|Z^a|^{2k})^{1/2}𝒫^{1/2} ,
  𝒫 = P(|Z_τ|>d/2) + P(|Z^a_τ|>d/2) + P(|Z^a_τ|≥d) ,   P(|Z_τ|>·) from Theorem V.1's B1 .
```
**Instrument check — a decorrelated confirmation of the V-b refuter's MAJOR-1**
(`a6_instrument_checks.py`, written from the theorem statement, importing none of their code).
At their design point (`f=0`, `φ₀=30°`, `δ=7.5°`, `d = ρ₀ sinδ`, `L=10`, window
`θ_max = 4(1−√(2/3)) = 0.7340137`, `c_G = √(3/2)θ_max = 0.8989795`, `\hat K₂ = 1`):

| quantity | this seat | refuter | rel |
|---|---|---|---|
| `ε_bulk` | `3.4734538e-3` | `3.4735e-3` | `1.3e-5` |
| `Eh1 = (r₀/R_-²)E[W]` | `1.1572236e-2` | `1.1572e-2` | `2.0e-5` |
| `Eh2 = 4N E[W]/d` | `3.7089353e-1` | `3.7089e-1` | `9.5e-6` |
| `Eh2/Eh1` | `32.050290` | `32.050` | `9.1e-6` |
| `E_4` | `1.9522576e-2` | `1.9523e-2` | `2.2e-5` |
| **correction factor** `(Eh1+Eh2)/(seat's sized E_hess)` | **`87.16576`** | **`87.166`** | `2.7e-6` |
| `θ_max/I₂` | `1.4401176` | `1.440118` | `3.0e-7` |

So the V-b refuter's factor **87.166** and its `θ_max/I₂ = 1.440118` are confirmed here from
the theorem statement alone.  `C_K` is carried symbolically; the tables below fix `C_K = 1`
(what the V-b seat's §7 expects), and the refuter's corrected requirement is `C_K ≤ 0.261` to
keep `E_hess ≤ 1/L`, `C_K ≤ 9.08·10^{-3}` to keep it below `ε_bulk`.

**The budget below is computed at `ε = 0` (sharp radial edges).**  With `ε > 0` the datum is
Lipschitz (needed for (S1)'s "smooth" and for `(H3)`'s ball to be inside the plateau) and `d`
must be reduced to `min(ρ_* − ρ₀e^ε, …)`; at `f = 1`, `ε = 0.25` that takes `d` from `0.7653669`
to `0.7159746`, a `6.4534 %` change that moves nothing in the tables.  `C_E`, by contrast, moves a
lot (§1.4).

### 2.6 What the bootstrap does NOT close — the blocks, stated precisely

> **BLOCK 1 (not fatal).**  The loop gain of the *algebraic* fixed point is
> `g = (15π/4)κ_δ c e^{2c_G}/r_h ≈ 55` at `c = c_*`, so no contraction exists.  This is
> **not** fatal: (2.1)–(2.2) is a differential comparison, which needs no contraction.  It is
> recorded because the natural "assume `μ ≤ 2μ̄`, prove `μ ≤ μ̄`" reading of a bootstrap
> genuinely fails here, and a reader who writes it that way will conclude the argument is
> impossible when it is only expensive.

> **BLOCK 2 (fatal; the highest-value missing lemma).  A restart lemma off the ansatz
> manifold.**
> *What is missing:* a version of **Theorem A / Corollary A1** (hence of Theorem Γ), and of
> **Lemma 1 / Lemma T′ Corollary 2**, for a field of the form `η₀∘Φ^{-1}` with
> `|Φ − Λ| ≤ μ|Λ|`, rather than for the exact strained plateau `η₀∘T_λ^{-1}`.
> *Why it is missing:* Theorem A's mechanism is the profile's structure — **one** jump, at the
> equator, of size exactly `2λM`, plus bounded variation away from it, on a field that is
> **homogeneous of degree −1**.  `η₀∘Φ^{-1}` is neither homogeneous nor of controlled total
> variation: a `C¹` perturbation of the map at relative size `μ` can multiply the total
> variation of the angular profile without bound.  `(B.2)` (`r|∇a| = (1−t²)|a′(t)|`) and
> `(B.3)` (the first-order ODE that removes the second angular derivative) both use
> homogeneity in an essential way, and both are the whole of Part B.
> *What it would buy:* everything.  With the reference map re-derived from the state at `N`
> intermediate times, the map error resets and the window cost falls from `e^{pc}−1` to
> `N(e^{pc/N}−1)`.  §3.4 prices it: `L_*` from `4.53·10^{14}` to `311` at `N = 16`.

> **BLOCK 3 (fatal, and the same class as BFG's own defect).  `Γ = ‖∇u‖_∞` for a field that is
> not the exact strained plateau.**
> Hypothesis (i) of `B(t)` is Theorem Γ, which is proved only for `η_λ`.  For the true field
> the needed object is `r|∇a|` for a bounded, compactly supported, `z`-odd vorticity on a shell
> of log-width `L` — i.e. a Calderón–Zygmund gradient bound for an `L^∞` datum, which is
> exactly the `‖∇u‖_∞ ≲ ‖ω‖_∞ log(...)` problem.  Partially available: the far and inner parts
> of `r∇a` **are** Parseval-based and therefore non-perturbative (`far-near-kernel-lemma`
> constants `1.6049285` and `0.1324254`, computed by the L3v seat).  What is missing is the
> **collar's** gradient part, which is precisely where the logarithm lives.  Note that this is
> the same estimate whose absolute-constant version BFG's p.10 weighted-BMO step asserts and
> which two campaign referees located as false; a repair here would be worth more than the
> assembly.

> **BLOCK 4 (soluble, and solved here).  Origin vs material shell.**  Lemma T′ bounds `a` at
> the origin (its refuter's Finding 6).  §2.4 replaces the modal route by the far/near kernel
> lemma, which is a statement at an arbitrary field point and needs only `|ω^θ| ≤ M_t` and
> `z`-oddness.  Cost: `C_a = 8.697888` instead of the measured `0.069`, a factor `126.06`.

> **BLOCK 5 (fatal as stated; repairable).  `(H3)` fails for the mollified family.**  Theorem
> V.4 and Cor. V.5 require `η₀ ≡ −M sgn(z)/r` **exactly** on `B(x_*,d)`.  The mollified taper
> never equals that anywhere (the V-b refuter's Finding 4; the V-b seat's own §8(a) measures
> `(Δ₅η₀/η₀)/(−1/r₀²) = 1.004` at `3ρ₀`, `2.036` at `1.5ρ₀`).  So the viscous column of the
> budget is computed for a *different* datum than the one the theorem is about.  Repair: either
> restate V.4 with `|Δ₅η₀ + η₀/r²| ≤ ϑ M/r³` on the ball (a one-line change to Step 2's Taylor
> remainder, not written), or use a datum that is exactly the bare plateau on the ball, which
> the piecewise-linear `Θ_ε`, `h_δ`, `min(1,·)` construction of §1.1 **does** provide away from
> its three kinks — this is why §1.1 uses `min(1,·)` rather than `tanh`.  With §1.1's datum
> `(H3)` holds with `d` as in §2.5, and BLOCK 5 is closed for *this* family but not for the
> family the runs used.

> **BLOCK 6 (CLOSED by `s3close/hk2`, with a sting).  `(H-K2)`.**  That seat's final statement,
> quoted:
> ```
>   PROVED    K2 <= 66.6622 M/rho0  on the whole tube N_tau, lambda in [1,3/2], with NO log(R/rho0)
>             -- for the datum whose RADIAL cutoff is tanh-mollified at w0 = 0.25 rho0
>   MEASURED  K2 = 2.1903 M/rho0 (tracked point), 3.0202 (global probe; binding at the taper)
>   FALSE     for the sharp-radial-edge datum at f = 0:  K2 = +infinity ;
>             with an inset, K2 ~= 8 M/(f rho0).
> ```
> Better than `(H-K2)` as posed: no `L` at all, so `K₂ ≤ C_K M/ρ₀` with `C_K = 66.6622`, not
> `C_K M L/ρ₀`.  The sting is the third line, and it lands on §1.1: **a sharp radial edge makes
> `K₂` infinite**, and hk2's sharpening of the V-b refuter's Finding 4a says the coupling lemma
> needs the *global* `K₂`, so an inset does not rescue it — the radial cutoff must be mollified,
> full stop.  §1.1's `Θ_ε = min(1, log(ρ/ρ₀)/ε, log(R/ρ)/ε)` is only `C^{0,1}`: `∇η` jumps
> across the two kink spheres `ρ = ρ₀e^ε`, `ρ = Re^{-ε}`, so `∇²η` carries a surface measure
> there and hk2's Cor. 4.3 (which uses `d·sup_{B(x,d)}‖∇²η‖_F`) applies only if `B(x_*,d)` misses
> them.  **Correction to §1.1:** replace `min(1,·)` by hk2's `tanh` ramp, or require
> `f > e^ε − 1 + d/ρ₀`.  With the `tanh` ramp, `C_K = 66.6622` and BLOCK 6 is closed.

> **BLOCK 7 (bookkeeping).**  `Λ` is not the flow map of an incompressible field
> (`J_Λ ≠ λ²`); Lemma T′ prices the mismatch at `2κ_s/L` but no one has shown a
> volume-preserving map exists within `O(c/L)` of `Λ`.  Also: §1.1's datum is Lipschitz, not
> `C^∞`; (S1) and BFG are stated for smooth / `L²∩L^∞` data respectively, and a `C^∞`
> mollifier changes every constant by `O(1)` at fixed `ε`.

---

## 3. ERROR BUDGET

All numbers from `a2_budget.py` (grid) and `a4_subwindows.py` (restart pricing).
`C_K = 1`, `f = 1` unless stated, `δ = 7.5°`, `δ_m = 5°`, `φ₀ = 30°`, `s = 7.661298`.
Columns:

| column | `C′` | `C_a` | `P(|Z_τ|>d/2)` | Lemma T′ sensitivity |
|---|---|---|---|---|
| `proved` | `119.33` | `8.697888` (far/near) | Theorem V.1 `B1` | as proved |
| `proved_oddonly` | `76.11` | `8.697888` | Theorem V.1 `B1` | as proved |
| `proved_sharpsens` | `76.11` | `8.697888` | Theorem V.1 `B1` | × `1/6.0668` (the record's own measured conservatism — **not** a theorem) |
| `measured` | `11.74` | `0.069` (R9, measured) | exact affine Gaussian | × `1/6.0668` |

### 3.1 The term table

Every `ε` term, its scaling and its constant:

| term | expression | scaling | value at `L=160, c=c_*` (proved / measured) |
|---|---|---|---|
| `ε_δ` (taper) | `1 − 2κ_δ` | `L⁰` | `5.5754e-4` / same |
| `ε_ℓ` (e-folds lost) | `(log2 + log(1+f) + 3log(3/2))/L` | `1/L`, `c⁰` | `0.0162668` / same |
| `ε_{C_a}` (deviatoric) | `λ C_a/(κ_δ L)` | `1/L`, `c⁰` | `0.163176` / `0.0012945` |
| `ε_{T′}` (Lemma T′) | `(2π/r_h)[3μ(1+μ_J)/(2(1−μ)⁵)+3μ_J/8]` | `e^{pc}/L` | `∞` / `∞` |
| `ε_bulk` | `σ_y/r₀²` | `1/L`, `∝ J₂(c)` | `5.919e-5` / same |
| `E_hess` | `½K₂e^{c_G}τV₁[r₀/R_-²+4N/d]` | `C_K e^{3c_G}/L²` | `8.064e-3` / `2.121e-3` |
| `E_4` | `(r₀/R_-⁵)[(trC)²+2tr(C²)]` | `1/L²` | `1.840e-3` / same |
| `E_tail` | `2N𝒫 + …` | `e^{−e^{−2c_G}q/(4n)}` | `16.35` / `2.7e-117` |
| `ε_{0.069}` (R9 offset) | `0.069·λ/(κ_δ L)` | `1/L` | — / `1.294e-3` |
| `ε_{C′}` (Γ's `O(M)`) | inside `c_G = (3/2+C′/L)c` | `1/L` | `c_G = 1.819` / `1.215` |

The `L`- and `c`-scan (`a2_results.json['grid']`, `f = 1`, `C_K = 1`):

**column `proved`** (`ε_ℓ`, `ε_{C_a}`, `ε_bulk`, `E_hess`, `E_4`, `E_tail`, `ε`):

| `L` | `c` | `ε_ℓ` | `ε_{C_a}` | `ε_bulk` | `E_hess` | `E_4` | `E_tail` | `ε` |
|---|---|---|---|---|---|---|---|---|
| 10 | 0.25 | 0.2603 | 2.611 | 3.77e-4 | 10.94 | 0.0337 | 16.39 | ∞ |
| 10 | 0.5 | 0.2603 | 2.611 | 6.70e-4 | 5.20e5 | 0.141 | 16.42 | ∞ |
| 10 | 1 | 0.2603 | 2.611 | 1.08e-3 | 5.86e14 | 0.919 | 16.5 | ∞ |
| 40 | 0.25 | 0.06507 | 0.6527 | 9.42e-5 | 2.23e-3 | 2.11e-3 | 0.112 | ∞ |
| 40 | 0.5 | 0.06507 | 0.6527 | 1.68e-4 | 0.1423 | 8.81e-3 | 16.37 | ∞ |
| 40 | 1 | 0.06507 | 0.6527 | 2.69e-4 | 239.7 | 0.0575 | 16.4 | ∞ |
| 160 | 0.25 | 0.01627 | 0.1632 | 2.36e-5 | 3.92e-5 | 1.32e-4 | 1.9e-21 | ∞ |
| 160 | 0.5 | 0.01627 | 0.1632 | 4.19e-5 | 5.60e-4 | 5.51e-4 | 1.29e-4 | ∞ |
| 160 | 1 | 0.01627 | 0.1632 | 6.73e-5 | 0.03596 | 3.59e-3 | 16.36 | ∞ |
| 640 | 0.25 | 4.07e-3 | 0.04079 | 5.89e-6 | 1.81e-6 | 8.23e-6 | 4e-105 | ∞ |
| 640 | 0.5 | 4.07e-3 | 0.04079 | 1.05e-5 | 1.84e-5 | 3.44e-5 | 2.8e-24 | ∞ |
| 640 | 1 | 4.07e-3 | 0.04079 | 1.68e-5 | 5.46e-4 | 2.25e-4 | 0.01358 | ∞ |
| 2560 | 0.25 | 1.02e-3 | 0.0102 | 1.47e-6 | 1.05e-7 | 5.14e-7 | 0 | ∞ |
| 2560 | 0.5 | 1.02e-3 | 0.0102 | 2.62e-6 | 9.79e-7 | 2.15e-6 | 3e-104 | ∞ |
| 2560 | 1 | 1.02e-3 | 0.0102 | 4.21e-6 | 2.42e-5 | 1.40e-5 | 4.1e-13 | ∞ |

**column `measured`:**

| `L` | `c` | `ε_ℓ` | `ε_{C_a}` | `ε_{T′}` | `ε_v` | `ε` |
|---|---|---|---|---|---|---|
| 10 | 0.25 | 0.2603 | 0.02071 | ∞ | 0.04677 | ∞ |
| 40 | **0.25** | 0.06507 | 5.18e-3 | **0.1870** | 2.69e-3 | **0.3560** |
| 40 | 0.5 | 0.06507 | 5.18e-3 | ∞ | 0.01428 | ∞ |
| 160 | **0.25** | 0.01627 | 1.29e-3 | **0.02780** | 1.83e-4 | **0.04858** |
| 160 | 0.5 | 0.01627 | 1.29e-3 | ∞ | 8.51e-4 | ∞ |
| 640 | **0.25** | 4.07e-3 | 3.24e-4 | **6.29e-3** | 1.58e-5 | **0.01140** |
| 640 | **0.5** | 4.07e-3 | 3.24e-4 | **0.0918262** | 6.01e-5 | **0.1072409** |
| 640 | 1 | 4.07e-3 | 3.24e-4 | ∞ | 6.02e-4 | ∞ |
| 2560 | **0.25** | 1.02e-3 | 8.09e-5 | **1.54e-3** | 2.09e-6 | **3.21e-3** |
| 2560 | **0.5** | 1.02e-3 | 8.09e-5 | **0.01740** | 5.70e-6 | **0.01943** |
| 2560 | 1 | 1.02e-3 | 8.09e-5 | ∞ | 4.01e-5 | ∞ |

The pattern is the whole story: **short windows close, the doubling window does not.**
`c = 0.25` gives `λ(τ) = e^{κ_δ c} = 1.13307` — nowhere near `3/2`.  `c = c_* = 0.8114` is
`∞` in every column at `L ≤ 2560`.

### 3.2 The feedback exponent

Linearising (2.1) about `m = 0`:
```
   p/(ML) = 3/2 + C′/L + (3/2)(2κ_δ)(15π/4)(9/4)/r_h ,        p c = the window cost exponent.
```

| column | `p/(ML)` (`L→∞`) | `p c_*` | `e^{pc_*}` |
|---|---|---|---|
| `proved` | `53.91` (with `C′/L`) / `41.974` (`L→∞`) | **`34.0573`** | `6.1784e14` |
| `proved_oddonly` | `49.59` | `40.23` at `c=1`; `34.06` at `c_*` | — |
| `proved_sharpsens` | `15.78` / `8.1714` (`L→∞`) | **`6.6301`** | `757.59` |
| `measured` | `9.345` / `8.1714` (`L→∞`) | **`6.6301`** | `757.59` |

`C′/L` is invisible at large `L`; the exponent is set by Lemma T′'s constant `15π/4 = 11.781`
times the two geometric factors `(3/2)` (`|Λ| ≤ λ|x|`, `|(r,−2z)| ≤ 2|X|`) and `(9/4)`
(`|Λ| ≥ λ^{-2}|x|`) whose product is `3.375`, times `2κ_δ ≈ 1`.  The record's own measurement says Lemma T′'s constant
is conservative by `6.0668` (`0.34913942/2.11817051`); dividing by that (`6.0668328`) takes `pc_*` from `34.0573`
to `6.6301` and `L_*` from `4.5253·10^{14}` to `3.3145·10^{4}`.  **The reach of the entire theorem is
one constant in Lemma T′.**

### 3.3 `L_*` and `Λ_*` (one window, `N = 1`)

| column | `L_*` (`ε≤1/2`) | `log Λ_*` | `L_*` (`ε≤0.1`) | `log Λ_*` |
|---|---|---|---|---|
| `proved` | **`4.5253·10^{14}`** | **`9.0506·10^{14}`** | `1.3846·10^{15}` | `2.7693·10^{15}` |
| `proved_oddonly` | `4.5253·10^{14}` | `9.0506·10^{14}` | `1.3846·10^{15}` | `2.7693·10^{15}` |
| `proved_sharpsens` | `3.3145·10^{4}` | `6.6280·10^{4}` | `6.8682·10^{4}` | `1.3736·10^{5}` |
| `measured` | `4.8315·10^{3}` | `9.6664·10^{3}` | `9.9881·10^{3}` | `1.9979·10^{4}` |

`log Λ_* = 2L_* + 3.3643455` from (1.2).  `Λ_* = e^{9.05·10^{14}}` is not a number anyone will
ever write down; it is reported because the brief asked for it honestly.

`ε` at the two `L` the brief names, at `c = c_*`, best `f`:

| column | `N` | `ε(160)` | `ε(640)` |
|---|---|---|---|
| `proved` | 1 | **∞** (bootstrap does not close) | **∞** |
| `proved` | 16 | `2.886` | `0.1778` |
| `proved_sharpsens` | 16 | `0.2809` | `0.05506` |
| `measured` | 1 | **∞** | **∞** |
| `measured` | 4 | `0.04091` | `8.682e-3` |
| `measured` | 16 | `0.02700` | `5.598e-3` |

### 3.4 What a restart lemma is worth (BLOCK 2, priced)

`N` sub-windows with the map error reset: amplification `N(e^{pc/N}−1)` instead of `e^{pc}−1`.

| column | `N` | `pc/N` | `N(e^{pc/N}−1)` | `L_*` (`ε≤1/2`) | `L_*` (`ε≤0.1`) |
|---|---|---|---|---|---|
| `proved` | 1 | 34.057 | `6.18e14` | `4.5253e14` | `1.3846e15` |
| `proved` | 2 | 17.029 | `4.97e7` | `1.3923e8` | `4.2601e8` |
| `proved` | 4 | 8.514 | `1.99e4` | `7.7146e4` | `2.3617e5` |
| `proved` | 8 | 4.257 | `556.9` | `1829.6` | `5719.1` |
| `proved` | **16** | 2.129 | `118.4` | **`310.99`** | `1041.3` |
| `proved` | 32 | 1.064 | `60.76` | `156.91` | `538.43` |
| `proved` | 64 | 0.532 | `44.97` | `138.05` | `428.01` |
| `measured` | 1 | 6.630 | `756.6` | `4831.5` | `9988.1` |
| `measured` | 4 | 1.658 | `16.99` | `34.02` | `85.71` |
| `measured` | **16** | 0.414 | `8.215` | **`28.65`** | `67.56` |

`log Λ_*` at `N=16`: `625.35` (proved, `ε≤1/2`), `2086.0` (proved, `ε≤0.1`), `60.67`
(measured, `ε≤1/2`).  `L_* = 311` sits inside the band the L3v refuter already quotes for
Theorem Γ alone (`L ≥ 429 … 1435` for a 10 % gradient margin), i.e. **with the restart, the
bootstrap stops being the binding constraint and Theorem Γ's own crossover becomes it.**
That is the sharpest thing this note has to say about where to work next.

### 3.5 Confrontation with the 19 viscous runs

`a3_numerics.py`.  Rows transcribed with provenance; every derived number computed here.
My own one-parameter fit over the 19 runs: `M T_d = 1.05338/log(R/s*)`, rms `8.776 %` of the
mean `T_d` — against the record's `1.0538` and `8.8 %`; `Q′ = M T_d log(R/s*) = 1.0272 ±
0.0762` (sd, `ddof=1`; `ddof=0` gives `0.0742`, the corner seat's figure).  `κ` for the runs'
datum `ω^θ = −M sin 2φ` recomputed from `far-near-kernel-lemma` Consequence A: **`0.4000000000`
= 2/5 exactly**.

The theorem is an **upper** bound on `T_d`, so the runs must lie **below** it.

| abscissa | `min` | `max` | `mean` | runs above the `ε=0` bound |
|---|---|---|---|---|
| `L = log(R/ρ₀)` | `1.068` | `3.271` | `1.500` | **19 / 19** |
| `L = log(R/s*)` | `0.938` | `1.185` | `1.013` | `8 / 19` |

**Every one of the 19 runs violates the `ε = 0` bound at `L = log(R/ρ₀)`.**  This does *not*
falsify the theorem — the runs sit at `L ∈ [1.386, 4.852]`, and the theorem claims nothing
below `L_*`.  What it does say is: at the only `L` where the truth is known, the true `ε` is
`+0.068` (best run, `N=7`, `Re0=400`, `L = 4.852`) to `+2.27` (worst, `N=2`, `L = 1.386`),
against a *proved* budget of `ε = ∞`.  The proved budget is therefore conservative by, at
minimum, the whole of `e^{34}`.  The log shift is accounted for exactly: replacing `ρ₀` by the
self-selecting inner scale `s* ≈ 8.5√(ν/M)` moves the mean ratio from `1.500` to `1.013` — i.e.
**the runs are consistent with the assembled constant `2log(3/2)/κ` once the abscissa is the
octave count the flow actually pays for**, which is the campaign's own surviving descriptive
law and not an independent confirmation of anything.

### 3.6 `(H-K2)` plugged in

`a8_hk2_plugin.py`, run after `s3close/hk2/PROOF.md` filed.  `C_K` was `1` in the `a2` tables;
the hk2 values are `66.6622` (proved) and `3.0202` (measured global probe).  Best `f` over
`{0.25, 0.5, 1, 2, 4, 8}`:

| `C_K` | column | `N` | `ε(160)` | `ε(640)` | `ε(2560)` | `L_*(ε≤½)` | `L_*(ε≤0.1)` |
|---|---|---|---|---|---|---|---|
| `1` (a2 baseline) | proved | 1 | ∞ | ∞ | ∞ | `4.5253e14` | `1.38463e15` |
| `1` | proved | 16 | `2.8863` | `0.17779` | `0.037805` | `309.93` | `1038.29` |
| `1` | measured | 16 | `0.027000` | `5.5982e-3` | `1.7327e-3` | `28.653` | `67.204` |
| **`3.0202`** (hk2 measured) | proved | 1 | ∞ | ∞ | ∞ | `4.5253e14` | `1.38463e15` |
| **`3.0202`** | proved | 16 | `3.0128` | `0.17876` | `0.037846` | `312.39` | `1041.03` |
| **`3.0202`** | measured | 16 | `0.038992` | `6.4298e-3` | `1.7831e-3` | `39.884` | `89.399` |
| **`66.6622`** (hk2 proved) | proved | 1 | ∞ | ∞ | ∞ | `4.5253e14` | `1.38463e15` |
| **`66.6622`** | proved | 16 | `7.8684` | `0.20457` | `0.039122` | **`366.96`** | `1116.32` |
| **`66.6622`** | measured | 1 | ∞ | ∞ | ∞ | `4832.77` | `9995.96` |
| **`66.6622`** | measured | 16 | `0.32007` | `0.023209` | `2.9709e-3` | `131.10` | `283.22` |
| `8/f` (sharp edge) — **illegitimate**, see BLOCK 6 | proved | 16 | `2.9328` | `0.17961` | `0.038068` | `312.14` | `1048.19` |

**`(H-K2)` costs almost nothing.**  At `N = 1` the `L_*` values do not move at all — the
`e^{34}` of the bootstrap swamps `E_hess` entirely.  At `N = 16` the proved column moves
`309.93 → 366.96` (`+18 %`) and the measured column `28.65 → 131.10` (`×4.6`).  So closing
`(H-K2)` was necessary and is not sufficient, exactly as the hk2 seat says in its own §10
(*"Is V-b now PROVED given (H-K2)?  No."*).  The `C_K = 8/f` row is kept only to show that the
sharp-edge idealisation would cost about the same; it is **not** a legitimate option, because
the coupling lemma needs a global `K₂` and that is `+∞` for a sharp radial edge.

**Third independent reproduction of the `87.166` factor.**  The V-b refuter found it; `a6` here
reproduces it from the theorem statement (`87.16576` vs `87.166`, rel `2.7e-6`); `s3close/hk2`
§10's control L-11 reproduces it a third time (`E_hess/ε_bulk = 110.111083` at `L=10`,
`K̂₂ = 1`, ratio `87.1658`, crossing at `L = 1101.1`).  Three seats, three instruments, no code
shared.  The V-b seat's withdrawn `1.263238` stays withdrawn.

---

## 4. THE BFG CONFRONTATION

Quoted verbatim from `sharp/literature-short-time/txt/bfg.txt`, the text layer of
`pdf/bfg-1704.05546.pdf` (arXiv:1704.05546v4 = ARMA **231** (2019)).  Page numbers are the
arXiv page numbers, fixed by the running heads at lines 364 (`8`), 417 (`9`), 542 (`11`),
597 (`12`).  The OCR renders `Δ` as `4` and drops some fraction bars; I mark those `[sic]`
rather than silently repairing them.

**Preamble, p.8 (line 407 ff.) — the hypothesis, in their words:**

> "In addition to the initial vorticity ω0 being bounded, a suitable decay of ω0 at infinity
> will be required (we chose ω0 in L2 for convenience); however, it is worth noting that this
> is a 'soft assumption', i.e., there will be no quantitative dependence on kω0 k2 in the
> proof." `[sic: kω0 k2 = ‖ω₀‖₂]`

**Theorem 8, p.8 (line 414):**

> "Theorem 8 (real setting). Let the initial datum ω0 be in L2 ∩L∞ . Then there exists a unique
> mild solution ω in Cw [0, T ], L∞ where T ≥ 1c kω01k∞ for an absolute constant c > 0."
> `[sic: T ≥ (1/c)·(1/‖ω₀‖_∞)]`

**Theorem 10, p.11 (line 582):**

> "Theorem 10 (complex setting). Let the initial datum ω0 be in L2 ∩ L∞ , and M a constant
> larger than 1. Then there is a constant c(M ) > 1 such that there exists a unique mild
> solution ω in Cw [0, T ], L∞ where T ≥ 1 c(M ) 1 kω0 k∞ , and for any t in (0, T ] the
> solution ω is the R3-restriction of a holomorphic function ω defined in the domain
> Ωt = x + iy ∈ C3 : |y| < 1 p c(M ) √ t ; moreover, kω(t)kL∞ (Ωt ) ≤ M kω0 k∞ ."
> `[sic: T ≥ (1/c(M))·(1/‖ω₀‖_∞); |y| < (1/√(c(M)))√t; ‖ω(t)‖_{L^∞(Ω_t)} ≤ M‖ω₀‖_∞]`

The equations are at `ν = 1` (their displayed (5): `∂t ωj − 4ωj + ui ∂i ωj = ωi ∂i uj`, `4 = Δ`).

### 4.1 The clause contradicted

Take `M = 3/2` in Theorem 10.  Then `c(3/2)` is an absolute constant and, since `ℝ³ ⊂ Ω_t`,
`‖ω(t)‖_{L^∞(ℝ³)} ≤ (3/2)‖ω₀‖_∞` for all `t ≤ T` with `T ≥ 1/(c(3/2)‖ω₀‖_∞)`.  Hence

```
    M₀ T_d(u₀)  ≥  1/c(3/2)  >  0        for EVERY u₀ with ω₀ ∈ L²∩L^∞ ,
```
with `c(3/2)` depending on nothing — not on `‖ω₀‖₂`, not on `ν` (both sides are invariant under
`u ↦ λu(λx,λ²t)` at fixed `ν`), not on the datum's geometry.  That is a **log-free, energy-free
lower bound `T_d ≥ c/M`**, i.e. `𝒯(Λ) ≥ 1/c(3/2)` for every `Λ`.  The assembled theorem of §1.3
gives `M T_d ≤ 2log(3/2)(1+ε)/L → 0` as `L → ∞` at fixed `M`.  **The two cannot both hold.**

### 4.2 Does the family satisfy every BFG hypothesis?

| BFG hypothesis | the family of §1.1 | verdict |
|---|---|---|
| `ω₀ ∈ L^∞` | `‖ω₀‖_∞ = M` exactly | ✔ |
| `ω₀ ∈ L²` | computed (`a1`): `‖ω₀‖₂² = 3.92736 M²R³` (`ε=0`), `2.42077 M²R³` (`ε=0.25`) — finite at every `L` | ✔ |
| decay at infinity ("soft") | compact support | ✔ |
| no dependence on `‖ω₀‖₂` in the constant | BFG's own sentence, quoted above | ✔ (their claim, not a hypothesis on us) |
| `ν = 1` | scale-invariance of `M₀T_d` lets us set `ν = 1` | ✔ |
| solution class `C_w([0,T],L^∞)`, mild | our solution is smooth (Lipschitz datum ⇒ classical after `t>0`) and unique in that class | ✔ |

**No hypothesis of BFG's is failed by the family.**  The confrontation is therefore not voided
by a hypothesis mismatch, and the stake is symmetric exactly as `RETURN_ADDENDUM_8` records it.

### 4.3 Status of the confrontation

**NOT TRIGGERED.**  (S3) is not proved here (§2.6), so BFG Theorem 10 is not contradicted by
anything in this note.  The located defect in BFG's displayed proof (the p.10 weighted-BMO
inequality without the local average) is not re-verified here and remains, in the campaign's
own words, a located defect rather than a refutation.  What §2.6 adds: **BLOCK 3 is the same
estimate**.  The assembly needs a gradient bound for `∇u` from a bounded shell vorticity with a
constant that does not eat the logarithm; BFG's Theorem 8 sketch asserts one with an absolute
constant.  If the campaign ever proves BLOCK 3 with a good constant it will have either
repaired BFG's step or exhibited why it cannot be repaired — and either outcome settles the
line faster than assembling (S3) does.

---

## 5. STATUS TABLE, DEPENDENCY GRAPH, AND WHAT A REFEREE SHOULD ATTACK

### 5.1 Status per step

| # | step | status |
|---|---|---|
| S1 | the datum, its parameters, the trajectory margins (`φ(3/2) = 62.833°`, equator margin `27.167° > δ_m`) | **PROVED** (exact algebra + `a1`) |
| S2 | `κ_δ`, `r_h`, `Φ_{h_δ}(λ) ≥ 1` on `[1,3/2]` | **PROVED** (R7; reproduced here) |
| S3 | `‖u₀‖₂² = (1/π)∫ψ₁η dx₅`; `G_l = ρ_<^lρ_>^{-(l+3)}/(2l+3)`; `D_l = 1/(5(l+4))` | **PROVED** (derived here; `C_E` reproduces the record to `4.0e-8`) |
| S4 | `log Re_E = 2L + 2log s + (2/5)log C_E`, `L = ½log Re_E − log(sC_E^{1/5})` | **PROVED** (exact) |
| S5 | contradiction setup ⇒ `λ ≤ 3/2` free on the window | **PROVED** |
| S6 | (2.1)–(2.2) is a valid majorant / continuity argument given (i)–(iv) | **PROVED** (ODE comparison) |
| S7 | the driver decomposition (l=1 rate error / label-radius error / deviatoric) | **PROVED** given R6's Consequence A |
| S8 | `C_a = 8.6975` at `φ=30°`, non-perturbative, `z`-odd constants valid at every `t` | **PROVED** (R6 + **R8**) |
| S9 | `ε_{T′}` from Lemma T′ at `(μ,μ_J)` | **PROVED** (R1) |
| S10 | `ε_v = ε_bulk + E_hess + E_4 + E_tail` | **PROVED-modulo-BLOCK-5** (R4, R5).  `(H-K2)` is no longer a modulus: `s3close/hk2` closes it at `C_K = 66.6622` for a datum with a mollified radial cutoff, and shows it is **FALSE** for a sharp radial edge (§3.6, BLOCK 6) |
| S11 | (i) `Γ ≤ 2a(0,s)+C′M` for the **true** field | **NOT ESTABLISHED — BLOCK 3** |
| S12 | Theorem A / Cor A1 / Lemma 1 for the **true** field (needed to restart, and needed for S11) | **NOT ESTABLISHED — BLOCK 2** |
| S13 | the assembled `ε(L)` and `L_*` | **COMPUTED** (conditional on S11–S12) |
| S14 | `𝒯(Λ) ≤ c₂(1+ε)/logΛ`, `c₂ = 4log(3/2)` | **NOT PROVED** (rests on S11, S12) |
| S15 | BFG hypothesis check | **PROVED** (the family satisfies all of them) |
| S16 | BFG contradiction | **NOT TRIGGERED** |
| — | `κ_δ(5°)` correction to `prove-lagrangian` §2 | **ESTABLISHED** |
| — | `C_E` depends on the radial mollification at `O(1)` (22 % at `ε=0.1`) | **ESTABLISHED** |
| — | the `87.166` `E_hess` correction, reproduced here independently (`a6`) and by `s3close/hk2` — three seats | **ESTABLISHED** |
| — | `(H-K2)` moves `L_*` by `0 %` at `N=1` and `+18 %` at `N=16` (proved column) | **COMPUTED** (`a8`) |

### 5.2 Dependency graph (input file → theorem → where used here)

```
rebuild/far-near-kernel-lemma/NOTE.md
    Consequence A (l=1 interior mode = e-fold count, exactly constant in the shell)
        -> S7, the identification of what Lemma T' actually prices
    Prop 1 (far, C1 = 0.291999 z-odd), Prop 2 (inner 0.014754, collar 3.999218/sin + pi/8)
        -> S8, C_a = 8.6975  [the only non-perturbative estimate in the chain]

lower/prove-lagrangian/NOTE.md
    Lemma 1 + taper deficit (Phi_{h_delta} >= 1)      -> S2, the conservative c2 = 4log(3/2)
    Lemma 2 (eta <= 0 on {z>0} for all t)             -> S8, licenses the z-ODD constants of R6
                                                          at every t  [first load-bearing use]
    L3 (l=1 resonance)                                -> S7, the shape of the reference field v

write/lemma-T-shell-dependent + refute-*
    Lemma T' (T'), Cor 3 composition, kappa_s = sqrt6-2 -> S9, eps_T'; and the constant 15pi/4
                                                          that sets the feedback exponent pc

write/L3v-and-gamma-bound + refute-*
    Theorem A / Cor A1                                 -> inside Theorem Gamma; BLOCK 2
    Theorem Gamma, C' <= 119.33 (refuter's correction)  -> (i) of B(t); BLOCK 3

write/V-b-bulk-viscous-loss + refute-*
    Theorem V.4, Cor V.5, refuter's sizing (87.166)     -> S10, eps_v
gaps/gap-V-aronson/NOTE.md
    Theorem V.1 (B1)                                   -> S10, P(|Z_tau| > d/2)

sharp/viscous-numerics + bfg/corner-numerics           -> S13's confrontation (19 runs)
sharp/literature-short-time/txt/bfg.txt                -> S15, S16
s3close/hk2/PROOF.md   [FILED 2026-09-08, after a1-a7]  -> C_K = 66.6622 (proved), 3.0202
    (H-K2): K2 <= 66.6622 M/rho0, no log; FALSE for   (measured global); plugged in by a8;
    a sharp radial edge (K2 = +inf at f = 0)          forces the tanh ramp in section 1.1
```

### 5.3 What a referee should attack

1. **The geometric constants `(3/2)` and `(9/4)` in (2.1).**  They enter the feedback exponent
   multiplicatively (`pc ∝ (3/2)(9/4) = 3.375`) and they are crude sups over the whole shell
   (`|Λ|/|x| ∈ [λ^{-2}, λ]`, attained only on the axis and the equator respectively).  A
   `φ`-resolved treatment could plausibly halve `pc`, which is worth `e^{17}` in `L_*`.  I did
   not attempt it.
2. **Whether the driver really carries the full `κ_δ M L`.**  §2.3 charges the map error the
   *whole* strain rate times `ε_{T′}`.  If the l=1 rate error were instead proportional to the
   *change* in the e-fold count induced by the perturbation — which is what one would expect,
   since a rigid relabelling of the shells does not change the count — the driver could be
   smaller by `O(1/L)` and the whole exponential would disappear.  **This is the single most
   valuable thing to attack in this note.**  I could not see how to make it work: Lemma T′'s
   estimate is two-sided and is stated against `a_ref`, not against the *variation* of
   `a_ref`, and its own K1 control shows that a `μ`-only bound is impossible.
3. **BLOCK 5, and whether §1.1's `min(1,·)` datum really satisfies `(H3)`.**  It does away
   from the three kink sets, but the ball `B(x_*,d)` must also miss the *radial* ramp: the
   correct `d` is `min(ρ_* − ρ₀e^ε, ρ_* sin(φ₀−δ), ρ_* sin(π/2−δ_m−φ₀))`, and the tables use
   the `ε = 0` version `min(f, …)`.  At `f = 1`, `ε = 0.25` this is `0.7159746` against
   `0.7653669` (`6.4534 %`); at `f ≤ e^ε − 1` the tracked point is *inside* the ramp and the row is void.
   The `L_*` bisections minimise over `f ∈ {0.05, …, 32}`, so their small-`f` end is only
   legitimate at `ε` correspondingly small.  Check whether the `L_*` values are attained at
   `f` large enough for the `ε` one wants.
4. **`E_tail` via Theorem V.1 vs the affine Gaussian.**  The `proved` column's `E_tail` swings
   over 120 orders of magnitude between `c = 0.25` and `c = 1` at `L = 160` because `B1`'s rate
   carries `e^{-2c_G}` and `c_G` carries `C′/L`.  That is not a bug, but it means the viscous
   column is dominated by an artefact of the proof route (the V-b refuter's §8 says the same),
   and a Girsanov route would remove it.
5. **The restart pricing of §3.4.**  It assumes the *only* effect of a restart is to reset `m`,
   with the strain deficits adding.  A real restart lemma would also have to re-establish the
   datum hypothesis (D) and the profile structure at each restart, and those costs are not in
   the table.  The `N`-column is an upper bound on the value of BLOCK 2, not a promise.
6. **`C_E`'s mollification dependence.**  I compute `C_E` for a piecewise-linear `Θ_ε`.  A
   `C^∞` mollifier gives a different `C_E` at the same nominal `ε`; anything quoting `Λ_*` must
   fix the mollifier, not just `ε`.
7. **The `19/19` violation in §3.5.**  A referee should check whether the runs' datum
   (`−M sin2φ`, `κ = 2/5`, inner cut well above the viscous floor) is close enough to §1.1's
   family for the comparison to mean anything.  My view: it is not, and the comparison's only
   honest content is the sign (`the true ε at L ≈ 5 is O(0.1)`, not `O(e^{34})`).

---

## 6. GATE, AND FILES

`check_constants.py` re-asserts every number displayed above against `a1`–`a7`'s stored JSONs.

```
2964 CHECKS, 2964 PASS, 0 FAIL
DOC AUDIT: 418 distinct numeric tokens in ASSEMBLY.md,
           381 traced to this folder's JSONs, 37 quoted from a named source, 0 untraced
MUTATION:  4084 leaves, 3322 caught, 762 blind, coverage 0.8134
           of the 762 blind, 363 are NON-FINITE -- v -> 1.5v + 0.37 is the identity on +-inf,
           so they are blind by construction of the mutation operator;
           399 are finite.  Coverage over finite leaves: 3322/3721 = 0.89277.
```

The pass count is `len(CHECKS)`, printed, never a literal.  No check multiplies its input by
zero.  No check compares a literal against itself; every one rebuilds its target from
mathematics or from other stored numbers.  The blind set is **fully enumerated** by
`(file, last key)` in `mutate_log.txt` — 94 classes, printed in full, not truncated to a
sample.  The 23 "quoted from a named source" tokens in the document audit are listed with
their source in `check_constants.py`'s `FOREIGN` dictionary; they are constants of other
seats, cited, not computed here.  `gate_stats.json` stores the gate's own statistics so that
the numbers in this box are themselves auditable.

What the gate **cannot** catch, stated so it is not mistaken for coverage: it checks that the
pipeline is internally consistent and that the displayed numbers are the numbers the scripts
produced.  It does **not** check that (2.1)–(2.2) is the right majorant, that `C_a` is the
right constant, or that BLOCK 2 and BLOCK 3 are the only blocks.  That is what §5.3 is for.

| file | what it establishes |
|---|---|
| `a1_datum_energy.py` / `a1_results.json` | the datum; `H_l(δ,δ_m)` by composite Gauss with the kinks split out (bang–bang coefficients to `1.1e-6` of the exact rationals); `κ_δ`, `r_h`, `Φ_{h_δ}`; the energy identity, the Green function, `D_l`, and `C_E` (sharp cap reproduces the record to `4.0e-8`); the `L ↔ log Re_E` dictionary; `‖ω₀‖₂²` and `‖ω₀‖_∞` for the BFG hypothesis check; the material-point trajectory and its margins |
| `a2_budget.py` / `a2_results.json` | the bootstrap majorant ODE (2.1)–(2.2); the Theorem-V.4 viscous budget; the `L × c × column` grid; the `f`-sweep; `L_*`, `Λ_*`; the self-consistent window |
| `a3_numerics.py` / `a3_results.json` | the 19 viscous runs, `κ = 2/5` for their datum recomputed from `far-near-kernel-lemma` Consequence A, my own one-parameter fit (`1.05338`, rms `8.776 %`), and the ratios against the `ε=0` bound on both abscissas |
| `a4_subwindows.py` / `a4_results.json` | the restart pricing: `pc/N`, `N(e^{pc/N}−1)`, `L_*(N)` for four columns and `N ≤ 64` |
| `a5_report_numbers.py` / `a5_results.json` | `ε(160)`, `ε(640)`, and `Λ_*` for every `(column, N)` |
| `a6_instrument_checks.py` / `a6_results.json` | decorrelated reproduction of the V-b refuter's `E_4`, `Eh1`, `Eh2`, `32.050`, **`87.166`** and `θ_max/I₂ = 1.440118`; the feedback-exponent algebra; the small-`δ` deficit law |
| `a7_derived.py` / `a7_results.json` | the quantities quoted inline in this document's prose, computed so they are traceable |
| `a8_hk2_plugin.py` / `a8_results.json` | `(H-K2)` plugged in: the budget and `L_*` at `C_K = 66.6622`, `3.0202`, `8/f` and `1` |
| `check_constants.py` | the gate: 2889 checks, the document audit, the mutation self-test with the blind set fully enumerated |
| `mutate_log.txt` | the enumerated blind set |
| `SHA256SUMS` | computed with `shasum -a 256`, never typed |

**FL-000 stands.  Nothing here touches the headline problem.**

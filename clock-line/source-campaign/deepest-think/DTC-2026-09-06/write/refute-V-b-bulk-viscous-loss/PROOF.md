# REFUTATION of `V-b-bulk-viscous-loss`

Seat `write/refute-V-b-bulk-viscous-loss`, DTC-2026-09-06.  Laws `TORMENT NEXUS/LAWS.md`.
Target: `write/V-b-bulk-viscous-loss/PROOF.md` (sha256 `707ade2a…`, verified) and its
scripts `b1`…`b5`, `check_constants.py`.  Every number below came out of a script in **this**
folder (`x1`…`x7`), run by me; the target's scripts were re-run from a copy in a scratch
directory, never in place.  Nothing outside this folder was written.

**VERDICT: REFUTED (MAJOR).**  The mathematics of Theorem V.4 is sound and every exact
identity reproduces independently.  What does not survive is the **sizing**: the two tables
that carry the seat's three headline claims (§5's `E_hess/ε_bulk` table with its `\hat K₂`
thresholds, and §5's tail/inset table) do **not** size the terms of Theorem V.4.  Corrected
against the theorem's own definitions the `∇²u` price is **87.17× larger at every `L`**, and
the tail price is **larger, not smaller, than `gap-V`'s** at the `L` where the comparison was
sold.  The status line "PROVED-MODULO-(H-K₂)" is therefore true of the theorem and false of
the campaign statement the seat draws from it.

---

## 0. What reproduces, and what I re-derived by my own route

| item | verdict |
|---|---|
| `SHA256SUMS` (12 files) | `shasum -a 256 -c` → **12/12 OK** |
| `b1` (exact identities), `b2` (exact Gaussian instrument), `b4` (`∇²u` necessity), `b5` (sizing) re-run from a copy | **bit-reproduce** the stored JSON and the displayed tables |
| `check_constants.py` | baseline **passes**; it is a real gate (see §6) |
| `r³Δ₅(1/r) = −1` (Cartesian ℝ⁵ **and** the `(r,z)` form) | **CONFIRMED** by my own sympy (`x4` §a, residual `-1` / `-1`) |
| `Δ₅η_P = −η_P/r²` | **CONFIRMED** (residual `0`) |
| `Hess(1/r)` spectrum `{0, −1/r³ ×3, 2/r³}`, trace `−1/r³` | **CONFIRMED** |
| `sup_{|v|=1}|∂_v^k(1/r)| = k!/r^{k+1}` via Legendre | **CONFIRMED** (residuals `0`, `k = 1…6`) |
| `Δ_y^k(1/r)` coefficients `−1, −3, −45, −1575` | **CONFIRMED** |
| Lemma V.3 propagator `Ψ_s = diag(e^{−μ(s)}×4, e^{2μ(s)})` | **CONFIRMED** (I substituted it back into `Ψ̇ = −S(τ−s)Ψ`; both block residuals `0`) |
| Lemma V.3′ `½C_τ:∇²η_P(x₀) = Mσ_y/r₀³`, `∂/∂σ_z ≡ 0` | **CONFIRMED** (residuals `0`) |
| Prop. 3.1 `V ≡ 0` for affine `X` | **CONFIRMED** in `n = 3` by my own construction (`['0','0','0']`), and `V = −2ε + O(ε²)` on a curved map |
| `I₂ = 4/5 − 16√6/135`, `I₄ = −4/7 + 27√6/28`, `θ_max = 4(1−√(2/3))`, `c = √(3/2)θ_max`, `C(τ) = 2log(3/2)` | **CONFIRMED** by my own sympy |
| Theorem V.1 Step 0 (the backward Feynman–Kac) | **CONFIRMED**: with `h(y,s) := η(y,τ−s)` one gets `∂_sh + L_sh = 0`, `L_s = −b(·,τ−s)·∇ + νΔ`, the same terminal-value problem the diffusion `dY = −b(Y,τ−s)ds + √(2ν)dW` solves; so `η(x,τ) = E[η₀(Y_τ)]` is right |
| Theorem V.4's proof steps 1–2 (coupling, Taylor, Isserlis, Cauchy–Schwarz) | **CONFIRMED** line by line; the inequalities are valid as written |

So this is not a refutation of the theorem.  It is a refutation of what the note says the
theorem costs.

---

## 1. FINDING 1 (MAJOR) — §5's `E_hess` is not Theorem V.4's `E_hess`.  Factor 87.17.

Theorem V.4 defines

```
   E_hess = (r₀/R_-²)·½K₂e^{c}τV₁  +  4N·E[W]/d ,        E[W] ≤ ½K₂e^{c}τV₁ .
```

§5 sizes it with **`E_hess/ε_bulk = ½e^{c}n(e^{2c}−1)/c · K₂τr₀`**, and `b5_sizing.py`'s
`hess_ratio()` implements exactly that.  Two independent errors:

**(1a) the denominator is the wrong `ε`.**  Dividing the first term (with `R_- = r₀`) by
`ε_bulk = ν∫₀^τ dt/r(t)² = ν τ I₂/(θ_max r₀²)` gives
`½e^{c}n(e^{2c}−1)/c · K₂τr₀ · (θ_max/I₂)`.  The published formula divides by `ντ/r₀²`
instead — i.e. by `prove-lagrangian`'s pricing, the very quantity §5 has just declared to be
`1.4401×` too big.  Measured (`x2` §a2): `Eh1(R_-=r₀)/(seat's value) = 1.440118` at
`L = 10, 40, 160, 640` — **exactly `θ_max/I₂`, at every `L`.**

**(1b) the second term is dropped, and it is the larger one.**  At the campaign design point
(`f = 0`, `φ₀ = 30°`, `d = ρ₀ sinδ`, `R_- = r₀ − d`):

| `L` | 10 | 40 | 160 | 640 |
|---|---|---|---|---|
| `Eh1 = (r₀/R_-²)E[W]` | 1.1572e-2 | 7.2326e-4 | 4.5204e-5 | 2.8253e-6 |
| `Eh2 = 4N·E[W]/d` | **3.7089e-1** | **2.3181e-2** | **1.4488e-3** | **9.0550e-5** |
| `Eh2/Eh1` | 32.050 | 32.050 | 32.050 | 32.050 |
| seat's sized `E_hess` | 4.3878e-3 | 2.7424e-4 | 1.7140e-5 | 1.0712e-6 |
| **(Eh1+Eh2)/(seat's)** | **87.166** | **87.166** | **87.166** | **87.166** |

`Eh2/Eh1 = 4N r₀ R_-²/(d r₀) `-type ratio is `L`-independent, so the error is a clean
constant factor `87.166` at every `L` (`x1` §1, `x2` §a2).  It is proportional to `K₂`, so it
is part of the `∇²u` price by the seat's own accounting; it is in `E_hess` in the seat's own
theorem statement.

**Consequence — the `\hat K₂` thresholds and §7's optimism.**  Largest `\hat K₂` (with
`K₂ = \hat K₂ M/ρ₀`) for which `E_hess ≤ ε_bulk`:

| `L` | 10 | 40 | 160 | 640 | 1200 |
|---|---|---|---|---|---|
| PROOF.md §5 | 0.792 | 3.166 | 12.666 | 50.663 | 94.99 |
| **Theorem V.4** | **0.00908** | **0.03633** | **0.14531** | **0.58123** | **1.0898** |

§7 remark 1 says: *"one expects `K₂ ≍ M/ρ₀`, i.e. `\hat K₂ = O(1)`, which by the table in §5
is already inside the admissible range for `L ≥ 13`."*  **The correct crossing is
`L = 1101.1`** (`x1` §2).  That is the single sentence in the note that makes (H-K₂) look
like a formality, and it is off by a factor `87`.

**Restated (H-K₂), honestly.**  With `K₂ = C M L/ρ₀` the ratio `E_hess/ε_bulk` is
`87.166 × 1.2632 × C = 110.1·C`, so to keep `E_hess ≤ ε_bulk` one needs `C ≤ 9.08e-3`, and to
keep `E_hess ≤ 1/L` one needs `C ≤ 9.08e-3/0.03473 = 0.261`.  (H-K₂) with an unspecified
`C` therefore does **not** deliver `O(1/L)` with a usable constant; the hypothesis the
campaign actually needs is `K₂ ≤ 0.26·M L/ρ₀`, which for `\hat K₂ = O(1)` means `L ≥ 38`
against the `1/L` budget and `L ≥ 1101` against the seat's own `ε_bulk` standard.

---

## 2. FINDING 2 (MAJOR) — §5's tail/inset table is not Theorem V.4's `E_tail`.

Theorem V.4's tail is built on
`𝒫 = P(|Z_τ| > d/2) + P(|Z^a_τ| > d/2) + P(|Z^a_τ| ≥ d)` — the **ball** event of (H3), at
threshold **`d/2`** (Step 1 needs `|Z_τ| ≤ d/2` *and* `W ≤ d/2` so that both `x₀+Z_τ` and
`x₀+Z^a_τ` sit inside `B(x₀,d)`), in `ℝ⁵`.  §5's table instead evaluates
`A·½exp(−d²/(2Var_n))` — a **one-dimensional half-space crossing at threshold `d`** along the
inner-edge normal.  These are different objects, and the difference is not small
(`x2` §b, exact affine Gaussian, `d = ρ₀ sinδ`, `f = 0`):

| `L` | seat's 1-D, `thr = d` | 1-D, `thr = d/2` | `P(|Z|>d)` in ℝ⁵ | `P(|Z|>d/2)` in ℝ⁵ |
|---|---|---|---|---|
| 10 | 9.1317e-2 | 3.2686e-1 | 2.4723e-1 | **8.6910e-1** |
| 20 | 1.6678e-2 | 2.1368e-1 | 4.1335e-2 | **6.1089e-1** |
| 40 | 5.5629e-4 | 9.1317e-2 | 1.7333e-3 | **2.4723e-1** |
| 80 | 6.1892e-7 | 1.6678e-2 | 4.5724e-6 | **4.1335e-2** |
| 160 | 7.6613e-13 | 5.5629e-4 | 4.5549e-11 | **1.7333e-3** |

(At `L = 10` the RMS of `|Z^a_τ|` is `0.875 ρ₀δ` while `d/2 = 0.5 ρ₀δ`: the escape event is
not rare, it is the typical event.)  The resulting `E_tail`, computed from the theorem's own
formula with the exact affine Gaussian (`x1` §3):

| `L` | 5 | 10 | 20 | 40 | 80 | 160 | 320 | 640 |
|---|---|---|---|---|---|---|---|---|
| §5's `A·P(edge)` | 1.139 | 0.4868 | 0.08890 | 2.965e-3 | 3.299e-6 | 4.084e-12 | 6.3e-24 | 1.5e-47 |
| **Theorem V.4's `E_tail`** | 22.88 | **17.64** | **11.17** | **4.392** | **0.7419** | **3.365e-2** | 2.079e-4 | 2.836e-7 |

**"One dissipation length suffices from `L = 24.603` (vs `1/L`)" becomes `L = 217.2`**
(`x1` §3).  The claim that this beats `gap-V`'s `L = 590` still holds — but by `2.7×`, not by
`24×`, and the number that beats it is not the one in the note.

**Where the damage comes from** (`x6`, four tail models, all with (H3)'s ball):

| `L` | as written | Cauchy–Schwarz repaired | affine-only (no coupling ⇒ no `d/2`) | §5's 1-D model |
|---|---|---|---|---|
| 10 | 1.7643e+1 | 1.7462e+1 | 2.2539e+0 | 4.8678e-1 |
| 20 | 1.1164e+1 | 1.1076e+1 | 3.7533e-1 | 8.8903e-2 |
| 40 | 4.3874e+0 | 4.3409e+0 | 1.5748e-2 | 2.9654e-3 |
| 80 | 7.4238e-1 | 7.2371e-1 | 4.2027e-5 | 3.2993e-6 |
| 160 | 3.3795e-2 | 3.0446e-2 | 4.1200e-10 | 4.0840e-12 |

So the `√𝒫` from Cauchy–Schwarz costs almost nothing; essentially **all** of the loss is the
`d/2` that the coupling step forces, i.e. it is the price of leaving the affine case — the
one thing §5's caveat says the table is exact for.  Even in the purely affine case
(third column, where the coupling is unnecessary and `𝒫 = P(|Z^a|≥d)` alone) the theorem's
tail at `L = 10` is `2.25`, not `0.487`: the remaining factor is the ball-versus-half-space
geometry of (H3), which the note never states it is replacing.

**The inset table.**  Running the seat's own search (`d = fρ₀`, its own targets) against the
theorem's own `E_tail + E_4`, with `φ₀` free and subject to (H3)'s `d < r₀` (`x4` §b):

| `L` | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|
| `f` — PROOF.md §5 (vs `ε_bulk`) | 0.25972 | 0.19228 | 0.14189 | 0.10441 | 0.07663 |
| `f` — Theorem V.4 (vs `ε_bulk`) | **not achievable** | 0.59059 | 0.39607 | 0.27875 | 0.19662 |
| `f` — Theorem V.4 (vs `1/L`) | 0.50770 | 0.37636 | 0.27528 | 0.20324 | 0.15006 |
| `f` — `gap-V` §5(ii) (vs `1/L`) | 0.9131 | 0.6590 | 0.4694 | 0.3382 | 0.2437 |

Headline claim 2 — *"the needed inset is `f = 0.2597 … 0.0766`, versus `gap-V`'s
`0.913 … 0.244`"* — compares a number that is not a bound on any term of Theorem V.4 with a
number that is a bound.  Against its own theorem and its own target the required inset at
`L = 20` is `0.591`, within `10 %` of `gap-V`'s `0.659`, and at `L = 10` the seat's own
standard cannot be met at all.  The corresponding `c₂` (`x7`):

| `L` | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|
| `c₂` — PROOF.md §5 | 1.50272 | 1.48105 | 1.47291 | 1.46985 | 1.46871 |
| `c₂` — Theorem V.4, `\hat K₂ = 0`, vs `1/L` | 1.53088 | 1.49186 | 1.47701 | 1.47143 | 1.46931 |
| `c₂` — Theorem V.4, `\hat K₂ = 1`, optimised geometry | **1.73196** | **1.58911** | 1.52617 | 1.49653 | 1.48214 |
| `c₂` — `gap-V` | 1.56987 | 1.50615 | 1.48229 | 1.47339 | 1.47003 |

*"the `c₂` inflation at `L = 10` falls from `1.5699/1.5843` to `1.5027`"* becomes
`1.5309` at `\hat K₂ = 0` and `1.7320` at the physically expected `\hat K₂ = 1` — i.e. with
the `∇²u` term switched on, **worse** than `gap-V`'s, because `E_hess ∝ 1/d` forces a huge
inset (`x5`: the optimiser goes to `f = 3.590`, `φ₀ = 75°`, at every `L`).
`8(1−√(2/3)) + O(1/L)` survives, as both notes say; the `O(1/L)` constant does not.

---

## 3. FINDING 3 (MAJOR) — `E_4` is never sized, and it too exceeds the leading term.

`E_4 = (r₀/R_-⁵)[(tr C_τ)² + 2tr(C_τ²)]` appears in the theorem and in no table.  At the
design point (`f = 0`, `φ₀ = 30°`, `d = ρ₀ sinδ`, so `R_- = r₀ − d = 0.3695ρ₀`):
`E_4 = 1.9523e-2` at `L = 10` against `ε_bulk = 3.4735e-3` — **5.62×** the quantity the whole
note is about — and `E_4/ε_bulk ∝ 1/L`, so it stays above `ε_bulk` until `L ≈ 56` (`x1`,
`breakdown_L10_f0`).  §10's status table lists every other term.

---

## 4. FINDING 4 (MINOR) — two hypotheses used but not stated

* **(H-K₂) is local; the proof uses `K₂` globally.**  §7 states (H-K₂) as a bound on
  `‖∇²₅b‖` *"on a `√(ντ)`-neighbourhood of the trajectory"*.  The coupling lemma (5.2) applies
  `‖∇b(x_s+θZ_s) − ∇b(x_s)‖ ≤ ½K₂|Z_s|` **pathwise, on every realisation**, and then takes
  `E[W]`; excursions of `Z` outside that neighbourhood carry weight in the expectation.  As
  written the theorem needs `K₂ = sup_{ℝ⁵×[0,τ]}‖∇²₅b‖` (which is what (H1) says).  The gap
  between (H1)'s global `K₂` and (H-K₂)'s local one needs a truncation argument that is not
  written.
* **(H3) is not satisfied by the campaign datum, at any `d > 0`.**  (H3) demands
  `η₀ ≡ η_P = −M sgn(z)/r` *exactly* on `B(x₀,d)`.  The mollified taper
  `ω^θ = −M tanh(sinφ/sinδ)tanh(cosφ/w)Θ(ρ)` never equals `−M sgn(z)` anywhere — the seat's
  own §8(a) measures `(Δ₅η₀/η₀)/(−1/r₀²) = 1.004` at `3ρ₀` and `2.036` at `1.5ρ₀`.  So
  Theorem V.4 and Corollary V.5 do not apply to the campaign's family as stated; §8's
  confrontation is run with a *different* instrument (an exact Gaussian average of the
  analytic datum), which is fine as a measurement but is not the theorem being tested.
  Neither point appears in the note's GAP list.

---

## 5. FINDING 5 — Proposition 6.1 is marked PROVED and its stated conclusion is false

Prop. 6.1's conclusion is *"the supremum of the loss over drifts of given `‖∇u‖_∞` is `+∞`"*,
argued by sending `κ → 0` in `b_κ(x) = βκ sin(x/κ)` through the leading-order formula (6.1).
Two problems:

1. **(6.1) is an expansion in `ν` whose parameter is `√(ντ)/κ`.**  It is validated in `b4`
   only for `κ ∈ [1,16]` at `√(2ντ) = 0.0447` — i.e. only in the direction where `‖b''‖`
   *decreases*.  The claim needs the opposite direction, `κ → 0`, precisely where the
   expansion has no footing.
2. **The family's own amplitude vanishes**: `‖b_κ‖_∞ = βκ → 0`, so as `κ → 0` the drift tends
   uniformly to `0` and the problem tends to pure diffusion with `η₀(x) = x`, whose loss is
   exactly `0`.  The loss must therefore turn over and vanish; the supremum over the family
   is **finite**, attained near `κ ≍ √(ντ)`.

MEASURED_X3_BLOCK

What survives, and is enough for the brief's question: **at fixed `ν` the loss is exactly
proportional to `‖b''‖` in the regime `κ ≫ √(ντ)`** (`b4`, `loss×κ` constant to `1.7e-4`
over `κ ∈ [1,16]`; reproduced here), so no bound on the loss can be a function of `‖∇u‖_∞`,
`ν`, `τ`, `‖∇η₀‖`, `‖∇²η₀‖` **alone that is valid uniformly in the fine structure of `b`**.
The corrected statement is a *non-existence of a `‖∇u‖`-only formula*, not a *`+∞`
supremum*.  The answer to the brief ("no `∇²u` needed?") is still **no**.

---

## 6. Controls, and the FL-043 question

* `check_constants.py` is a genuine gate: baseline passes, and dynamic mutation
  (`v ↦ 1.5v + 0.37`, one stored number at a time, 291 sampled across the five JSONs) is
  caught **107** times and missed **184** (`63.2 %`) — `x2` §d.  The misses are almost all
  stored intermediates that PROOF.md never displays (`T0_covariance/nu`, `cov_diag_mc/…`),
  which is the gate's stated scope.  **Not FL-043.**
* The seat's own R1–R4 can all fire, and R4 ("refinement must move the measurement toward the
  prediction") does fire.  **Not FL-043.**
* **But no gate in the folder can catch Findings 1–3**, because they are mismatches between
  PROOF.md's prose (`E_hess/ε_bulk = …`) and Theorem V.4's own definition of `E_hess`, and
  `check_constants.py` only checks that displayed literals equal what `b5` computed — never
  that what `b5` computed is what the theorem says.  That is the structural hole: the gate
  verifies transcription, not derivation.  My own instrument is separately controlled: the
  tail quadrature is Monte-Carlo checked at five parameter points, worst `|z| = 2.90` against
  a pre-declared `4σ` rule (`x2` §c) — a rule that can fire.

---

## 7. CORRECTED STATEMENT

> **V-b is PROVED-MODULO-(H-K₂′) as a theorem, and is NOT sized for the campaign.**
>
> Theorem V.4, Lemma V.3, Lemma V.3′, Prop. 3.1, Cor. V.5 and the exact identities (2.1)–(2.5)
> all stand as written and were independently re-derived.  In particular the leading term
> `s_C = ν∫₀^τ dt/r(t)²` with coefficient exactly `1` is EXACT, `σ_z` drops out identically,
> and `ε_bulk·L = 0.03473454` for the affine strain at `f = 0, φ₀ = 30°, δ = 7.5°`;
> `prove-lagrangian` §4(3) is conservative by `θ_max/I₂ = 1.440118`, and `gap-V-aronson`
> §5(iii)'s "optimistic by 6.2" is superseded.  Those four results are the note's real
> contribution and they survive.
>
> What must be withdrawn:
> * `E_hess/ε_bulk = ½e^{c}n(e^{2c}−1)/c·K₂τr₀` and the `\hat K₂` thresholds
>   `0.792 / 3.166 / 12.666 / 50.663`.  Correct: `× θ_max/I₂` and `+ 4N E[W]/d`, i.e.
>   `× 87.166`; thresholds `0.00908 / 0.03633 / 0.14531 / 0.58123`; and §7's "already inside
>   the admissible range for `L ≥ 13`" → **`L ≥ 1101.1`**.
> * the required hypothesis is not "(H-K₂) `K₂ ≤ C M L/ρ₀` for some `C`" but
>   **`K₂ ≤ 0.261·M L/ρ₀`** to keep `E_hess ≤ 1/L`, or `K₂ ≤ 9.08e-3·M L/ρ₀` to keep it below
>   `ε_bulk`; and (H1)'s `K₂` is global while §7's is local (Finding 4).
> * "`L_min = 24.603 / 48.337`" and the `f`-table `0.2597 … 0.0766` and the
>   `c₂ = 1.50272` at `L = 10`.  Correct, from the theorem: `L_min = 217.2` (vs `1/L`) at
>   `d = ρ₀ sinδ`; `f = 0.5077 / 0.3764 / 0.2753 / 0.2032 / 0.1501` (vs `1/L`, `\hat K₂ = 0`);
>   `c₂ = 1.5309` at `L = 10`, rising to `1.7320` once `\hat K₂ = 1` is switched on.
>   The comparison with `gap-V` should read: `gap-V`'s `min(B1,B2)` is loose by `~2.7×` in
>   `L_min` and `~1.8×` in `f`, not by `24×` and `3.5×`.
> * `E_4` must be tabled; it exceeds `ε_bulk` until `L ≈ 56`.
> * Prop. 6.1's "`sup = +∞`" (Finding 5); keep the `‖b''‖`-proportionality, drop the
>   unboundedness.
>
> Net effect on the campaign: **V-b does not yet buy `O(1/L)` with a usable constant.**  The
> `1/L` bulk term is now known exactly (`0.0347/L`), which is real progress; but the proved
> error terms around it exceed it until `L` in the hundreds (affine) or `L ≈ 1101` (general
> `b`, `\hat K₂ = 1`), so the "3.5 % of the `1/L` budget" reading is the leading term of an
> expansion whose remainder the note has not paid for.

## 8. Remaining gap

The single new open item this refutation creates is **not** (H-K₂) but the two-sided one
above it: Theorem V.4's `E_hess`, at `4N E[W]/d` with `N ≍ r₀/(ρ₀ sinδ)` and `d ≍ ρ₀ sinδ`,
carries an inverse dissipation length; and its `E_tail`, at `P(|Z_τ| > d/2)`, carries a
factor `4` in the Gaussian exponent that the coupling step forces.  Both are artefacts of the
*proof route* (truncate-and-couple on a ball), not of the physics.  A route that avoids them
— e.g. a Girsanov/Cameron–Martin change of measure from the affine reference instead of a
pathwise coupling, or a `C^{1,α}` Schauder estimate on the label-frame divergence-form
operator of Lemma V.0 (which has no drift and no zeroth-order term, so Krylov–Safonov applies
directly) — would keep the exact leading term and shrink the remainder from `O(K₂/d)` to
`O(K₂ r₀)`.  Until that is written, the honest per-step status of V-b is:

| step | status |
|---|---|
| exact identities (2.1)–(2.5), Lemma V.3, Lemma V.3′, Cor. V.5, Prop. 3.1 | **PROVED** (independently confirmed here) |
| Theorem V.4 as an inequality | **PROVED** (steps re-checked) |
| Theorem V.4 as a *useful* bound for the campaign | **NOT ESTABLISHED** — remainder exceeds the main term below `L ≈ 1101` (`\hat K₂ = 1`) |
| §5's `E_hess` sizing / `\hat K₂` thresholds | **WITHDRAWN** (Finding 1) |
| §5's tail/inset table and `L_min` | **WITHDRAWN as a sizing of Theorem V.4** (Finding 2) |
| `E_4` | **UNSIZED** (Finding 3) |
| Prop. 6.1 "`sup = +∞`" | **FALSE as stated**; the `‖b''‖`-proportionality stands (Finding 5) |
| (H-K₂) | still open, and now needs an explicit small constant |

**FL-000 stands.  Nothing here touches the headline problem.**

---

## 9. Files

| file | what it establishes |
|---|---|
| `x1_theoremV4_sizing.py` / `x1_results.json` | Theorem V.4's own `E_hess`, `E_4`, `E_tail` at the campaign design point, by my own route; the `87.166` factor; corrected `\hat K₂` thresholds; `L = 1101.1`; `L_min = 217.2` |
| `x2_breakdown_and_controls.py` / `x2_results.json` | term-by-term `E_tail`/`E_hess` breakdown; the `θ_max/I₂ = 1.440118` isolation; which convention change costs what; **Monte-Carlo control on my tail quadrature (worst `z = 2.90` against a `4σ` rule)**; mutation of the seat's gate (107 caught / 184 missed) |
| `x3_prop61_kappa_to_zero.py` / `x3_results.json` | Prop. 6.1 tested downward in `κ` on the seat's own PDE solver plus an independent backward-SDE Monte-Carlo with common random numbers |
| `x4_identities_and_inset.py` / `x4_results.json` | independent sympy re-derivation of (2.1)–(2.5), the Lemma V.3 propagator, Lemma V.3′ and Prop. 3.1 (my own `n = 3` construction); the corrected inset table |
| `x5_when_is_V4_nonvacuous.py` / `x5_results.json` | the theorem's total error minimised over `(f, φ₀)` subject to (H3); when it first drops below `1`, `1/L`, `ε_bulk` |
| `x6_repairable_tail.py` / `x6_results.json` | four tail models — as written / Cauchy–Schwarz repaired / affine-only / the seat's 1-D model — showing the damage is the coupling's `d/2`, not the Cauchy–Schwarz |
| `x7_c2_inflation.py` / `x7_results.json` | the `c₂` that follows from the corrected insets |
| `SHA256SUMS` | computed, never typed |

READ-ONLY inputs: `write/V-b-bulk-viscous-loss/*` (re-run from a copy in a scratch
directory), `lower/SYNTHESIS.md`, `gaps/gap-V-aronson/NOTE.md`,
`gaps/refute-gap-V-aronson/NOTE.md`.  Nothing outside this folder was written.

# REFUTER report — `write/L3v-and-gamma-bound`

Seat `write/refute-L3v-and-gamma-bound`, DTC-2026-09-06. Role: REFUTER.
Estate laws `TORMENT NEXUS/LAWS.md`. Every number below came out of a script in **this** folder
(`scripts/r1`…`r5`), written and run here; `SHA256SUMS` is computed, never typed. The target's
files were copied to `copy/` and re-run there; nothing outside this folder was written or edited.
Numerics falsify, they never prove; algebra marked EXACT is sympy-verified with residual `0`.

**Provenance.** `shasum -a 256 -c copy/SHA256SUMS` → all 19 files OK, so what I audited is what
the seat published. `python3 copy/check_constants.py` → `ALL 924 CHECKS PASS` (1.06 s), reproduced.

---

## 0. VERDICT

**REFUTED — MAJOR.** Part (a) survives in full and is independently re-derived here. Part (b)'s
*structure* survives — the reduction `‖∇u‖ ≤ 2|a| + r|∇a| + |ω^θ|` and the identity (B.3) that
removes the second angular derivative are correct and are the seat's real contribution. Part (b)'s
**headline constant is wrong**:

> **`C′ ≤ 33.5` is NOT obtained "using only proved inputs."** The column PROOF.md labels
> *"C′ (all inputs proved)"* substitutes the §A5 **COMPUTED** value of the L3v sum
> (`0.4280 … 0.7911`, whose tail contains a fitted residual-decay extrapolation) for `‖f‖_∞`,
> where the only thing Part A **proves** about `‖f‖_∞` is Corollary A1,
> `‖f‖_∞ ≤ 0.899159·√(2/π)(2λM+V_λ) = 3.62 … 6.46 λM`. That is a factor **6.1 … 10.1** larger.
> Substituting it and changing nothing else: **C′ ≤ 119.3** (Corollary A1 with the computed `V_λ`),
> **C′ ≤ 151.1** (Corollary A1 with the proved `V_λ ≤ 6λM`) — against the displayed **33.45**.
> The reach statement moves from `L ≥ 187 … 316` to **`L ≥ 679 … 899`** (up to **1435**), i.e.
> `log Re_E ≳ 1360 … 1800`, not `370 … 630`.

A second, independent arithmetic slip: **"43 … 80 (measured column)" is wrong** — the seat's own
`g8_results.json` holds `42.53 … 110.20`; the max is at `(δ=30°, λ=1.5)`, 38 % above the displayed
ceiling. The gate checks the crossover min/max for the other two columns and **not** for this one.

Both errors sit in the part of the pipeline the gate cannot see: they are among the 200 blind
mutation leaves, and PROOF.md §5's two claims about that set — that the 200 are "characterised"
and that "no leaf carrying a displayed constant is blind" — are **both false as written** (§4).

`survived` should be **false** for the (b) headline as stated, **true** for (a) and for (b)'s
structural reduction.

---

## 1. What I re-derived from scratch, and what survived

Instruments built here and used nowhere else: sympy with **symbolic exponents** (so the identity
checks are for a spanning family, not a sample), scipy QUADPACK adaptive quadrature, mpmath at 30
digits, and a 32-point composite Gauss rule at 800 panels per sub-interval (the seat's is 16-point).

| target statement | my route | result |
|---|---|---|
| `C_l^{3/2} = P_{l+1}′`, `N_l = (l+1)(l+2)/(l+3/2)`, `‖C_l^{3/2}‖_∞ = C_l(1) = (l+1)(l+2)/2` | sympy, `l ≤ 10` | residual **0** (`r1`) |
| `‖C_l‖_∞` attained at `t=1` | 2·10⁵ grid, `l < 60` | ratio max **1.0** |
| `H₁ = −5M/6, H₃ = 3M/40, H₅ = −247M/1680, H₇ = 1513M/40320, H₉ = −2773M/42240` | my own symbolic derivation `−2∫₀¹C_l√(1−t²)dt/N_l` | **exact match**, all five |
| the same by QUADPACK | adaptive, `epsabs=1e-13` | agrees to **≤ 2.8·10^{−16}** |
| `H_l(λ=1.5, δ=7.5°)` QUADPACK vs mpmath-30 | `l = 3,9,21,51` | **≤ 2.1·10^{−15}** rel. |
| **prefactor step** `l^{3/2}(l+3/2)/((l+1)(l+2)√(l+1)) ≤ 1` | **PROVED, not checked**: `(x+1)³(x+2)² − x³(x+3/2)² = 4x⁴ + (67/4)x³ + 25x² + 16x + 4`, every coefficient `≥ 0` ⟹ true for every `x > 0` | **stronger than the seat's** "checked over odd l < 10⁷" |
| `∫_{−1}^1 |t|(1−t²)^{−3/4}dt = 4` | exact | **4** |
| `V_λ ≤ 2λM(2 + √(sin θ_δ))` | QUADPACK, 20 (λ,δ) pairs | **holds everywhere**; equality (`V = 4λM`) exactly at δ = 0 |
| **THEOREM A** `|H_l| ≤ √(2/π)(2λM+V_λ)l^{−3/2}` | 20 (λ,δ) × 13 values of `l ∈ [3,301]`, my own `H_l` and my own `V_λ` | **holds everywhere**, worst ratio **0.4419** ⟹ slack **2.263** — matching the seat's reported 2.26–2.66 |
| Corollary A1 coefficient `(5/7)(ζ(3/2)−1−2^{−3/2}) = 0.8991585`, per-λM `5.7394` | own | **reproduced** |
| **B1** `∇₅b = a·diag(1,1,1,1,−2) + E`, `E_{zr} = E_{rz} − ω^θ` | sympy with `ψ₁ = q^i Z^j`, `q = |y|²`, **i, j symbolic** (spans every analytic axisymmetric ψ₁ by linearity) | residual = **zero matrix** |
| `tr ∇₅b = 2a` (so `b` is *not* 5-D divergence free — the `diag(1,1,1,1,−2)` is consistent) | same | residual **0** |
| `‖∇₅b‖_op = ‖∇₃u‖_op` | 400 random matrices of the proved shape | **1.3·10^{−16}** |
| `‖∇u‖_op ≤ 2|a| + r|∇a| + |ω^θ|` | 2000 random matrices | never violated; worst ratio **0.9985** — the bound is essentially **sharp**, not slack |
| `Δ₅(ρ^γ g) = ρ^{γ−2}[γ(γ+3)g + (1−t²)g″ − 4tg′]`, `∂_z(ρ^γ g)`, **B2**, **B3**, the integrating factor | sympy, `γ` and `n` symbolic | residual **0** in all four |
| kernel sups `K_{1/2}(0) = π/4`, `K_1(0) = 2/3`; `D′|_{p=1/2} = (3/2)[u√(1−u²)+arcsin u−π/2]`, `D(1) = 0` for both | sympy exact | **all confirmed**; hence `K_p′ ≤ 0` and (B.5) with `3π/4` and `2/3` is PROVED |

So Part A is sound and Theorem A / Corollary A1 do close `prove-lagrangian` §5's L3v row, exactly as
claimed. The structural half of Part B is sound. The refutation is entirely about **which** bound on
`‖f‖_∞` is fed into `C′`, and about two displayed numbers that the gate cannot see.

---

## 2. MAJOR-1 — `C′ ≤ 33.5` is not a proved constant

### 2.1 The substitution, at source

`copy/g8_assembly.py`:

```python
S_L3v = g5['L3v_sums'][key]['S']                       # the §A5 COMPUTED value
...
pr['adev_full']        = K_f * S_L3v + K_omega * pr['sup_omega_eff_env']
pr['r_grad_adev_full'] = 3 * S_L3v + 3 * pr['adev_full'] + pr['sup_s2_Wdev_env']
Cp_full = 2*(a1_offset + adev_full + edge_a + collar_a) + (r_grad_a1 + r_grad_adev_full + edge_grad + collar_grad) + lam
```

`Cp_full` is what PROOF.md prints in the column headed **"C′ (all inputs proved)"** and what §0
calls "**C′ ≤ 33.5** using only proved inputs" and what the `status` field calls "unconditional".

But `S` is built in `copy/g5_tail.py` as

```python
'S': p4 + jt + rt          # rt = g2's 'resid_tail_estimate'
```

with `rt` an extrapolation off a **fitted** decay exponent (`resid_p` = 1.998 at δ = 0, ≈ 2.50 at
δ > 0), and `S_lo/S_hi` a heuristic `[−|rt|, +3|rt|]` band. PROOF.md's own §A6 status line for A5 is
**COMPUTED**, not PROVED — and A5 is the input to the column labelled "all inputs proved."

The quantity Part A actually **proves** about the same object is Corollary A1:
`‖f‖_∞ ≤ Σ_{l≥3}|c_l|‖C_l‖_∞ ≤ 0.899159·√(2/π)(2λM + V_λ)`.

### 2.2 The corrected constant (`scripts/r3_honest_Cprime.py`, nothing else changed)

| case | seat's "fully proved" `C′` | `‖f‖` it used | **honest `C′`** (A1, computed `V_λ`) | `‖f‖` | **honest `C′`** (A1, proved `V ≤ 6λM`) |
|---|---|---|---|---|---|
| δ=0, λ=1 | 18.69 | 0.4280 | **75.99** | 4.3045 | **97.20** |
| δ=0, λ=1.5 | 33.39 | 0.6420 | **119.33** | 6.4568 | **151.15** |
| δ=7.5°, λ=1.5 | 33.42 | 0.6550 | **105.34** | 5.5206 | **150.99** |
| δ=15°, λ=1.5 | 33.45 | 0.7047 | **99.78** | 5.1919 | **150.29** |
| δ=30°, λ=1.5 | 32.61 | 0.7911 | **92.76** | 4.8611 | **148.16** |
| **max over the 12 cases** | **33.45** | — | **119.33** | — | **151.15** |

Inflation factor **3.57**.

**Reach.** `C′M ≤ 0.1·2a(0)` needs `L ≥ C′/(0.2κ)`:

| column | seat | honest |
|---|---|---|
| all-proved | `187 … 316` | **`679 … 899`** (A1 + computed `V`); up to **`1435`** (A1 + proved `V`) |

`log Re_E = 2L + O(1)` ⟹ `log Re_E ≳ 1360 … 1800`, not `370 … 630`.

### 2.3 The cheapest legitimate repair, and why it does not rescue 33.5

`H_l = 0` for even `l`, so Corollary A1 may sum over odd `l` only:
`Σ_{odd l≥3} l^{−3/2} = 0.6887612` in place of `Σ_{l≥3} l^{−3/2} = 1.2588220` — a free factor
**1.828**, which the seat did not take. The A1 coefficient becomes `0.4919723` (per-λM `3.1403`
instead of `5.7394`) and, re-running the assembly (`scripts/r3b_oddonly.json`):

**`C′ ≤ 76.11`, reach `L ≥ 429 … 583`.**

Still **2.3×** the displayed 33.5 and **2.3–3.1×** the displayed reach. The gap between the proved
envelope and the computed sum is ~6× and cannot be closed by bookkeeping: Corollary A1 bounds a sum
of absolute values by an `l^{−3/2}` envelope that is itself 2.26× lossy, while the true sum enjoys
the sign cancellation of `P_{l+1}(0)`.

### 2.4 What the corrected statement should say

> **THEOREM Γ (corrected).** On `Slab = {2λρ₀ ≤ |x| ≤ R/(2λ²)}`, for `λ ∈ [1,3/2]` and
> `δ ∈ {0°, 7.5°, 15°, 30°}`,
> `Γ := sup_slab ‖∇u‖_op ≤ 2a(0,t) + C′M` with
> * **`C′ ≤ 152` unconditionally** (Corollary A1 with `V_λ ≤ 6λM`), or **`C′ ≤ 120`** with the
>   computed `V_λ`, or **`C′ ≤ 77`** after the free odd-only refinement of Corollary A1;
> * `C′ ≤ 33.5` **only if** the §A5 computed value of the L3v sum is accepted as an input, i.e.
>   modulo a floating-point quadrature and a fitted residual-tail extrapolation — the same
>   epistemic status as A5, which is COMPUTED, not PROVED;
> * `C′ ≤ 18.5` with the measured `‖f‖_∞` and `C′ ≈ 4.3 … 11.8` fully measured (unchanged).
>
> Reach: `C′M ≤ 0.1·2a(0)` needs `L ≥ 429 … 1435` on the proved columns.

The *shape* of the result — `Γ ≤ 2a(0,t) + O(M)` uniformly in λ, with the `O(M)` explicit — is
untouched. Only the size of the `O(M)` moves, by a factor 2.3–4.5, and the crossover with it.

---

## 3. MAJOR-2 — the measured crossover range "43 … 80" does not track

PROOF.md §0: *"`C′M ≤ 0.1·(2a(0))` needs `L ≥ … 43 … 80` (measured column)."* The seat's own
`g8_results.json['crossover_L_for_10pct_measured']`:

```
min 42.533  (λ=1.00, δ=0°)          max 110.203  (λ=1.50, δ=30°)
```

`42.53, 72.74, 78.21 | 42.57, 73.06, 79.70 | 42.79, 74.99, 86.59 | 44.44, 85.11, 110.20`

"43 … 80" is the range over `δ ∈ {0°, 7.5°}` only. The theorem is asserted "for every λ ∈ [1,3/2]
and every δ ∈ {0°, 7.5°, 15°, 30°}", so the honest measured range is **`43 … 111`**.

`check_constants.py` asserts `crossover proved min/max` (75.7 / 178.5) and `crossover full min/max`
(187.1 / 315.9) but has **no** assertion on `crossover_L_for_10pct_measured`; all 12 of its leaves
are blind. That is exactly why a 924-check gate passed over it.

---

## 4. FL-043 — the gate's own claims about its blind set are false

`copy/check_constants.py` ends with `for m in misses[:40]: print(...)`. It prints **40** of the 200
blind leaves. PROOF.md §5 characterises all 200 from that printed sample. I enumerated all 200
(`scripts/r4_blind_leaves.py`, 1288 leaves total, 200 blind — reproducing the seat's count exactly):

| count | file | leaf key | PROOF.md §5 covers it? |
|---|---|---|---|
| 20 | g2 | `S_total_no_resid` | yes |
| 20 | g5 | `jump_tail_rem_bound` | yes |
| 20 | g5 | `S_lo` | yes |
| 20 | g5 | `resid_p` | yes |
| 12 | g4 | `partial_vs_cesaro_adev` | yes |
| 2 | g2,g4 | `secs` | yes |
| **20** | **g9** | **`l32_absH_at_4001`** | **no** |
| **20** | **g9** | **`analytic_V_upper`** | **no** |
| **12** | **g8** | **`proved_pieces.sup_s2_Wdev_env`** | **no** |
| **12** | **g8** | **`crossover_L_for_10pct_measured`** | **no** |
| **24** | **g8** | **`c_edge_taper_7.5deg[1..12]`, `c_edge_taper_15deg[1..12]`** | **no** |
| **11** | **g6** | **`supdC_check_ratio[0..10]`** | **no** |
| **4** | **g5** | **`sup_over_lambda[δ]`** | **no** |
| 3 | g6, g9 | `C_far_a_l1`, `C_in_grad_l3`, `szego_worst_ratio_n_plus_half` | no |

**106 of the 200 are outside the stated characterisation.** Two of them are load-bearing:

1. **`sup_s2_Wdev_env` is blind and it is a displayed constant** — PROOF.md (B.6) and the
   `constants` field both display `‖(1−t²)W_dev‖_∞ ≤ λM + (2/√3)|H₁|`, and its stored value feeds
   `pr['r_grad_adev_full']` and hence `C_prime_fully_proved`. It is blind because the gate
   recomputes `pr['r_grad_adev']` (the measured-norm variant) from its pieces but **never
   recomputes `pr['r_grad_adev_full']`**. So §5's "`C′ = the sum of its own pieces`" pipeline
   identity has an unchecked link precisely on the column the headline quotes. This directly
   contradicts §5's *"No leaf that carries a displayed constant is blind."*
2. **`crossover_L_for_10pct_measured` is blind** — MAJOR-2 above lives there.

Also worth recording: `g9['analytic_V_upper']` stores `4λ + (8/3)λ^{2.5}√δ`, which is **not** the
analytic envelope PROOF.md displays (`2λ(2+√(sin θ_δ))`). It is never used and never checked. Two
different `V` envelopes in one seat, one of them dead.

None of this is fatal — the gate is far better than no gate, and 84.5 % caught with a *computed*
pass count is real FL-043 discipline. But the §5 sentences must be corrected, the print truncated at
40 must be removed or the characterisation re-derived from the full list, and the two identities
above added.

---

## 5. MINOR-1 — (B.4) discards a homogeneous solution on an unstated hypothesis

The integrating-factor solution

`a_dev(t) = −(1−t²)^{−3/2}∫_{|t|}^1[3f(1−τ²)^{1/2} + (1−τ²)ω^θ_eff]dτ`

is the general solution **minus** `C(1−t²)^{−3/2}`. PROOF.md's entire justification for `C = 0` is
the parenthesis *"(a_dev being finite at one endpoint)"*. That is an assumption, and it is the
branch-selection step that makes (B.5) — the load-bearing bound of Part B — mean anything. Marked
PROVED; as written it is a sketch.

It is repairable in three lines, and I repaired it (`scripts/r2_partB_identities.py`):

> Put `s = 1−t`. The second-order equation `4f + (1−t²)f″ − 4tf′ = −W_dev` becomes
> `s(2−s)f″ + 4(1−s)f′ + 4f = −W_dev`, a **regular singular point** with indicial polynomial
> `2k(k−1) + 4k = 2k(k+1)` (sympy), roots **`k = 0, −1`**. Corollary A1 gives
> `Σ_{l≥3}|c_l|‖C_l‖_∞ < ∞`, so `f = Σ c_l C_l^{3/2}` is a uniform limit of continuous functions,
> hence **bounded**, hence on the `k = 0` branch (the `k = −1` branch is `~ s^{−1}`). Then
> `4f′ = −W_dev − 4f − s(2−s)f″ + 4sf′` gives `f′ = O(1)` when `W_dev` is bounded (δ > 0) and
> `f′ = O(s^{−1/2})` when `W_λ ~ (1−t²)^{−1/2}` (δ = 0); either way `(1−t²)f′ → 0`, so
> `a_dev(1) = −f(1)` is finite and `(1−t²)^{3/2}a_dev → 0`. ∎

Two consequences the seat should record: (i) **Part B is not independent of Part A even
structurally** — the branch selection uses Corollary A1; (ii) the seat's own falsification test
already caps any spurious `C`: their closed-form-vs-series agreement of `3.3·10^{−4}` at `|t| ≤ 0.99`
forces `|C| ≲ 10^{−6}`, since `(1−0.99²)^{−3/2} = 356`.

---

## 6. MINOR-2 — B7's `Γ = ML(1+O(1/L))` sentence overstates what Theorem Γ gives

PROOF.md §B7 and §3: *"`gap-V-aronson`'s standing assumption `Γ = ML` is now a theorem with a
priced `O(M)` correction"*, *"Theorem Γ supplies what that seat assumed: `Γ = ML(1 + O(1/L))`."*

Theorem Γ gives `Γ(t) ≤ 2κ(λ(t))ML + C′M`, and `2κ` runs from **1.000** (λ=1) to **1.500** (λ=3/2,
δ=0) / **1.4727** (λ=3/2, δ=7.5°). At the end of the window the bound is `1.47 ML`, i.e. **47 %
above** `ML` — not `ML(1+O(1/L))`.

The **integrated** statement is correct and is what `gap-V-aronson` (2.2) actually needs
(`c = ∫₀^τ‖∇₅b‖dt`, not the pointwise sup):
`c ≤ 2 log λ(τ) + C′c₀/L = 0.8109302 + C′c₀/L`. So the substance stands; only the sentence needs
rewriting to *"`∫₀^τ Γ dt ≤ 2log(3/2) + C′c₀/L`, which reproduces gap-V's `c = 2θ` with an explicit
`O(1/L)` correction; the pointwise `Γ` is up to `1.47 ML`, not `ML`."*

---

## 7. What I did NOT find

* No error in Theorem A, Corollary A1, or the exact coefficient work of A1/A2. Independently
  reproduced; the prefactor step is even upgradable from "checked" to PROVED (§1).
* No error in (B.1), (B.2), (B.3), the integrating factor, the kernel monotonicity, or the two
  kernel sups `3π/4` and `2/3`. All residual `0` on my own symbolic route with symbolic exponents.
* No error in the cross-seat agreements of §B8: I checked the cited neighbour values exist at
  source — `far-near-kernel-lemma` §§ carry `0.291999`, `0.014754`, `3.999218`, `π/8`, `+0.2167733`,
  `+0.2939`, `+0.1017`; `refuter-lagrangian-b` carries the seven partial sums `0.34990 … 0.41046`,
  the fitted `0.8254`, the estimate `0.4295` and `−0.01808 / +0.02838`. The seat's citations are
  faithful.
* The seat's own disclosures are honest as far as they go: gaps 1–6 correctly flag GAP T, the slab
  restriction, the un-reproduced λ>1 material-point row, the crudeness of the constants, and the
  edge sums' use of computed `H_l`. **What is missing from that list is the one that matters most:
  that `‖f‖_∞` in the "all-proved" column is also a computed input, not a proved one.**
* `V_bound_holds` reported `False` by my first pass at `(δ=0, λ=1.5)` is a floating-point artefact:
  `V = 6.0000000000000009` against the bound `2λ(2+√(sin 0)) = 4λ = 6`. The bound is attained with
  **equality** for the bang–bang cap, not violated.

### 7.1 Limits of this pass (what I did NOT independently re-derive)

* **B2's l = 1 constants** `sup_slab|a^{(1)} − κM log(ρ_out/ρ)| = 0.800000 κM` and
  `sup_slab r|∇a^{(1)}| = 1.146966 κM`: taken from the seat's `g6` and only checked for internal
  consistency by its own gate. They enter `C′` additively at `2(0.8κ) + 1.147κ ≈ 1.35 M`, i.e.
  ~4 % of the seat's `C′` and ~1 % of the corrected one, so they cannot change the verdict.
* **The edge sums** (`edge_a = 0.053 … 0.179 M`, `edge_grad = 0.42 … 1.06 M`) and the **ellipsoidal
  collar constants** (`0.2919985`, `0.0147543`, `1.6049285`, `0.1324254`): I confirmed the four
  displayed values exist at source in `far-near-kernel-lemma` and that the geometric factors
  `2^{−(l+4)}`, `2^{−(l−1)}` make the `l ≤ 200` truncation in `edge_sums()` harmless (the first
  dropped term is `O(2^{−204})`), but I did not rebuild the multipole/Parseval derivations.
* **`prove-lagrangian`'s `T_λ` model** and everything downstream of it (GAP T). Out of scope for
  both seats.
* The **λ > 1 material-point offsets** the seat records as NOT RUN: I did not build that instrument
  either, so that row stays PARTIAL.

---

## 8. Remaining gap after this refutation

1. **GAP T is untouched by both seats.** Everything about the time window remains conditional on
   `prove-lagrangian`'s `T_λ` model.
2. **The proved `C′` is 2.3–4.5× the advertised one**, so the "o(1) is genuine but the crossover is
   astronomically far" verdict is *worse* than the seat records: `L ≥ 429` at best, `L ≥ 1435` at
   the fully-proved end, against a campaign `L ≈ 8–10`. At `L = 10` the proved `O(M)` term is a
   **429 %–583 %** correction to `2a(0)` (odd-only A1), **679 %–899 %** (A1 as the seat states it),
   **972 %–1435 %** (A1 with the proved `V ≤ 6λM`) — against the seat's own **187 %–316 %**.
3. **Closing the 6× gap between Corollary A1 and the true sum** is the obvious next piece of work
   and is not attempted here: it needs an `l^{−3/2}` bound that keeps the sign of `P_{l+1}(0)`
   (alternating, so the sum telescopes) rather than an absolute-value envelope. Until that exists,
   `C′ ≈ 33.5` is a *measured* constant with a proved *shape*.
4. **The corner circles** `{z=0}∩{|x|=ρ₀}` and `{z=0}∩{|x|=R}` remain unpriced — the seat says so,
   and I confirm the octave factors `2^{−(l+4)}`, `2^{−(l−1)}` are the only thing making the edge
   sums converge.
5. The gate needs: the `r_grad_adev_full` pipeline identity, an assertion on
   `crossover_L_for_10pct_measured`, removal of the `misses[:40]` truncation (or a re-derived
   characterisation), and either a use or a deletion of `analytic_V_upper`.

---

## 9. Files

| file | what it establishes |
|---|---|
| `scripts/r1_foundations_and_thmA.py` / `r1_results.json` | independent A1 identities, the exact `H_l`, the prefactor step PROVED symbolically, `V_λ` and Theorem A at 20 (λ,δ) × 13 `l` |
| `scripts/r2_partB_identities.py` / `r2_results.json` | independent `∇₅b` structure (symbolic exponents), `tr = 2a`, `‖∇₅b‖=‖∇₃u‖`, the B1 bound's sharpness, `Δ₅(ρ^γ g)`, B2, B3, the integrating factor, the kernel sups, and the **Frobenius repair of (B.4)** |
| `scripts/r3_honest_Cprime.py` / `r3_results.json` | `C′` recomputed with Corollary A1 in place of the computed L3v sum |
| `scripts/r3b_oddonly.json` | `C′` with the free odd-only refinement of Corollary A1 |
| `scripts/r4_blind_leaves.py` / `r4_results.json` | all 200 blind mutation leaves enumerated and classified |
| `scripts/r5_l3v_sum_independent.py` / `r5_results.json` | independent two-sided evaluation of the L3v sum (32-point Gauss, 1600 panels, `l ≤ 12001`, analytic tail, QUADPACK cross-check) |
| `copy/` | the target seat's files, hash-verified, re-run in place |

# REFUTATION of `write/lemma-T-shell-dependent` — LEMMA T′ STANDS; four claims around it do not

Seat `write/refute-lemma-T-shell-dependent`, DTC-2026-09-06.
Adversarial re-check of `write/lemma-T-shell-dependent/PROOF.md`.
Every number below came out of `x1`–`x4` in **this** folder, written and run here; the target
seat's tree was re-run **from a copy** in a scratch directory, never in place, and its JSONs are
read read-only. Nothing outside this folder was written.

---

## 0. Verdict

| claim under test | verdict |
|---|---|
| **LEMMA T′** (the estimate, hypotheses (D),(H1),(H2),(H3), constants) | **STANDS.** Re-derived end to end by my own routes. |
| Steps S0–S6, Corollaries 1–4, variant (T′_Λ) — all marked PROVED | **STANDS.** None is a sketch; every constant tracks. |
| the kernel-sign correction to the task brief (`𝒦 = K(−·) = −3x_z/(8π²|x|⁵)`) | **STANDS** — re-derived from `−Δ₅ψ₁=η`, `a=−∂_zψ₁`, not quoted |
| the `κ_s ≈ 0.899c` correction to the task brief | **STANDS** (`κ_s=√6−2`, `μ_J L=2√6−4`, `=1.1086 c`) |
| the shear offset `0.21380335906432307 M`, `0.3491394233337655/L`, factor `6.0668` | **STANDS** (two of my own routes, one of which never touches `det DΛ`) |
| **"unresolved cross-check" on `P_{h_δ}(3/2)` (Corrections 4, gap 5)** | **REFUTED — it is a normalisation slip, not a discrepancy. Residual falls from 3.28e-2 to 3.93e-5.** |
| **"`p7` … an instrument sharing nothing with `p3`/`p4`" / "tests … `det DΛ`"** | **REFUTED — `p7` imports the same closed form for `det DΛ` that `lib5` uses.** |
| "`C(0+) = 15π/8`" quoted as the headline constant | **NOT IN FORCE in any of the 21 verification rows, nor in the campaign instance (0 of 21 have `μ_J ≤ μ`).** |
| "`a_ref` is the quantity step **P** of `prove-lagrangian` consumes" | **True at `σ = 0` only**; step P consumes it at every material shell `σ∈[0,1]`. Unstated reduction. |
| "12 non-trivial `Φ`" | 6 distinct maps × 2 values of `L`. |

**Severity: MINOR.** Nothing load-bearing is wrong. The lemma, its four corollaries, the variant
and every constant survive an independent re-derivation. What fails is one recorded "open
discrepancy" that is not open, one instrument-independence claim, and three scope/precision
overstatements.

---

## 1. Reproduction

`shasum -a 256 -c SHA256SUMS` on the target folder: **26/26 OK**. In a scratch copy,
`p1`,`p2`,`p3`,`p4`,`p5`,`p7` re-run clean; `check_constants.py` prints **553 checks, 553 PASS,
0 FAIL**; `p6_mutation.py` re-run gives **678 mutated, 678 caught, 0 blind, coverage 1.0000**;
`p8_doc_audit.py` re-run gives **295 distinct numbers, 0 untraced, 17 foreign**. Every printed
digit reproduces.

## 2. What I re-derived independently, and what it gave (`x1`, `x4`)

Routes chosen to be different from the seat's:

* **`det DΛ` (2.1) by finite differences of the actual 5×5 map** — no sympy, no hand-built
  matrix. On 400 random points with `λ = aρ^m` (which realises every pointwise `(λ, ρλ′/λ)`),
  `max rel = 5.99e-10` (the FD step's own error). **(2.1) confirmed.** I also confirm it by hand:
  choosing `y=(r,0,0,0)`, the `(y₁,z)` block has determinant
  `λ^{-1} − 2λ^{-2}λ′z²/ρ + λ^{-2}λ′r²/ρ`, and the three `S³` directions contribute `λ³`.
* **`𝒦` from the PDE, not from a quotation** (`x4`): `G₅ = 1/(8π²|w|³)` is harmonic in `ℝ⁵`
  (sympy residual `0`); `K = −∂_{w_z}G₅ = 3w_z/(8π²|w|⁵)` (residual `0`);
  `𝒦(x) = K(−x) = −3x_z/(8π²|x|⁵)` (residual `0`); and `a(0)/(ML) = +1/2` for the cap with `𝒦`,
  `−1/2` with the brief's `K`. **The seat's Correction 1 is right.**
* **The two shell identities by my own `mpmath` quadrature at five `λ` including `λ<1` and `λ=2.9`:**
  `∫₀^π sin²/g⁴ dφ = π/(2λ)` to `1.88e-31`; `∫₀^{π/2} cos sin²/g⁵ dφ = λ/3` to `0.0`.
* **`Q(λ)`**: the seat's closed form vs my quadrature, `3.59e-29`; `Q(1) = −1/15`;
  `Q(3/2) = −0.6812086838444756`.
* **The profile constants from scratch**: `λ̄ = 1.224744871391589 = √6/2`;
  `∫₀¹dσ/λ = 0.8277210825314643 = 5/9+√6/9`; `κ_s = 0.4494897427831781 = √6−2` (200 001-point
  grid max vs closed form, `2e-16`); `μ_J L = 0.8989794855663558`; `C_rel = 2.118170510687396`;
  `15π/8 = 5.890486225480862`; `15π/4 = 11.780972450961723`;
  `15√6π/16 = 7.214342794660484`; `μ_J L/c = 1.1085781089288433`.
* **`Step 0(b)`'s contraction constant, symbolically**: `∂log ĝ/∂log λ` is a convex combination of
  `+2` (at `ψ=0`) and `−1` (at `ψ=π/2`) — residual `0`; and `∂log g/∂log λ ∈ [−2,1]`. So the
  diffeomorphism condition is `L > 2κ_s = 0.8989794856`, met at every campaign `L`.
* **`a[η₀∘Λ⁻¹](0)` computed ENTIRELY IN THE `w`-VARIABLES** — the strongest check, because it never
  uses `det DΛ` at all. Inverting `Λ` (`|x| = |w| ĝ(ψ;λ(|x|))`) and integrating
  `∫_{Λ(S)}𝒦(w)η₀(Λ⁻¹w)dw = (3M/4)∫∫ λ(ρ)|cosψ|sin²ψ [1 − (∂log ĝ/∂log λ)D] dψ dlog ρ`
  reproduces the `x`-space value to `1.3e-6` (my `ψ`-quadrature straddles the kink of `|cosψ|`,
  which is the whole of the residual), and reproduces the offset `0.21381 M`, `L`-independent.
* **The closed-form offset** `(3M/2)∫₀¹(dlogλ/dσ)Q(λ)dσ = 0.21380335906432307`, matching the seat
  to every digit; relative `0.3491394233337655/L`; conservatism `2.118170510687396/0.349139… =
  6.066832815561182`.
* **All 21 verification rows re-checked with my own bound formula** and the seat's *measured*
  `(μ, μ_J, |Δ|)`: **21/21 hold**, my bound vs the seat's `max rel = 0.0`, worst slack `6.0668`
  for (T′) and `37.895` for (T′_Λ).
* **`P_{h_δ}` and `r_h`** reproduced: my `P_h(3/2)` differs from the seat's by `≤ 4.5e-7`;
  `κ_δ(7.5°) = 0.49972123052108863`.

**Conclusion of §2: the lemma is correct and every constant in it tracks.**

---

## 3. FINDING 1 (the only actual error) — the "unresolved cross-check" is a normalisation slip

The seat records, as Corrections item 4 and as **gap 5**:

> my `P_{h_δ}(3/2)` is `1.497990/1.491146/1.472733/1.442517/1.353721/1.245539/1.032403` against
> the synthesis's `1.4980/1.4914/1.4736/1.4444/1.3597/1.2584/1.0674` … the gap widens with `δ`
> (`3.4 %` at `30°`) … the discrepancy is **not** resolved in this note.

It resolves at source. `lower/prove-lagrangian/s2_taper.py` line 13 defines

```
Phi_h(lam) := P_h(lam)/P_h(1)
```

and `lower/SYNTHESIS.md` line 47 tabulates `Φ_{h_δ}(3/2)` — the **ratio**, not `P_{h_δ}(3/2)`.
The seat compared its unnormalised `P` against the other seat's normalised `Φ`. Dividing by
`P_h(1)` (`x2_taper_normalisation.py`, my own 30-dps quadrature):

| `δ` | `P_h(1)` | my `P_h(3/2)` | seat's | `Φ_h = P(3/2)/P(1)` | synthesis | rel, unnormalised | rel, normalised |
|---|---|---|---|---|---|---|---|
| 3° | 0.9999641 | 1.4979898 | 1.4979900 | 1.498044 | 1.4980 | 6.82e-06 | 2.90e-05 |
| 5° | 0.9998343 | 1.4911461 | 1.4911460 | 1.491393 | 1.4914 | 1.70e-04 | 4.49e-06 |
| 7.5° | 0.9994425 | 1.4727335 | 1.4727330 | 1.473555 | 1.4736 | 5.88e-04 | 3.05e-05 |
| 10° | 0.9986843 | 1.4425170 | 1.4425170 | 1.444417 | 1.4444 | 1.30e-03 | 1.21e-05 |
| 15° | 0.9956155 | 1.3537207 | 1.3537210 | 1.359682 | 1.3597 | 4.40e-03 | 1.30e-05 |
| 20° | 0.9897903 | 1.2455387 | 1.2455390 | 1.258386 | 1.2584 | 1.02e-02 | 1.08e-05 |
| 30° | 0.9672505 | 1.0324026 | 1.0324030 | 1.067358 | 1.0674 | 3.28e-02 | 3.93e-05 |

```
max relative residual, the seat's comparison : 3.2787e-02      ("3.4 % at 30 deg")
max relative residual, the correct comparison: 3.9253e-05      (= the synthesis's 4-decimal rounding)
```

So the two seats agree to the printed precision at **every** `δ`. The seat's own `P_h` numbers are
right (mine differ by `≤ 4.5e-7`); only the comparison was wrong. **Corrections item 4 must be
withdrawn and gap 5 deleted** — it is not a gap. `r_h` is unaffected: Corollary 2 legitimately uses
the *unnormalised* `P_h`, since `a_ref = (M/2)∫P_h(λ)dlog ρ ≥ (M/2)r_h A`.

Not load-bearing on (T′). But it is a "sources at source" miss: the definition is one line of the
other seat's own script.

---

## 4. FINDING 2 — `p7`'s independence is overstated in the one respect it is named for

PROOF.md §5.1 and the seat's file table say:

> `p3` and `p4` both use the meridian reduction … and the closed form (2.1) … `p7_montecarlo.py`
> avoids all of it … This tests the 5-D measure, the `SO(4)` reduction, `det DΛ` and the kernel in
> one shot. … an instrument sharing nothing with `p3`/`p4`.

`p7_montecarlo.py` line 38: `JL = lam**2*(1 + D*(sphi**2 - 2*cphi**2))`.
`lib5.py` `stage_Lambda`: `J = lam**2*(1 + D*(s**2 - 2*c**2))`.
Byte-for-byte the same closed form (2.1), on both sides of the comparison (`x4` checks the two
sources programmatically: `p7_independent_of_detDLambda = False`).

`p7` therefore does **not** test `det DΛ`; it assumes it. Nor is the meridian reduction avoided —
`𝒦`, `η₀` and `J_Λ` are all evaluated in the meridian variables `(λ, sinφ, cosφ)`; only the
*sampling* is 5-D. What `p7` genuinely checks is (a) the measure normalisation
`∫f dx₅ = L|S⁴|E[fρ⁵]`, `|S⁴| = 8π²/3` — a real and worthwhile check of `lib5`'s hardcoded
`2π²ρ⁴sin³φ`; (b) the `Λ`-specific reduced kernel `−(3/8π²)λ^{-2}cosφ/g⁵` against `lib5`'s generic
`cos(φ_f)e^{−4(u_f−u)}` form; (c) Monte-Carlo against Gauss–Legendre. That is worth stating, and
it is less than what is claimed.

No number is harmed: `det DΛ` is independently established by `p1`'s two symbolic routes, by the
refuter's `r1` finite differences at 400 points, and by mine at 400 points (`5.99e-10`).
The defect is in the advertised verification architecture, not in a result.

---

## 5. FINDING 3 — the headline constant is in force in none of the seat's own rows

`C(0+) = 15π/8` (absolute) and `15π/4 μ = 11.78097245 μ` (relative) are quoted in the Headline, in
the statement of Lemma T′, in Corollary 1 and in the status table. Both presume `μ_J ≤ μ`.
Measured over the seat's own stored rows (`x4`):

```
rows in which muJ <= mu : 0 of 21
```

and the campaign instance is `Φ = Λ`, where `μ = 0` and `μ_J = 2κ_s/L > 0` — the regime is not
merely untested, it is the opposite one. The operative constant everywhere the campaign uses the
lemma is `(3π/4)μ_J`, not `(15π/4)μ`.

This is **the same defect the refuter recorded as smaller finding #2 against `gap-T-lipschitz`**
("`μ_J > μ` in all twelve verification rows, so the clean corollary constant `C(0+) = 15π/8` …
is in force in none of them"). It is reproduced here without being noted. Not an error — (T′)
itself is two-parameter and is what is used — but the advertised constant is decorative.

Related, and smaller: `p4`'s "12 non-trivial `Φ`" are **6** distinct maps at 2 values of `L`
(`x4`: `n_distinct_nontrivial_maps = 6`). And the numerical verification of (T′) on maps other
than `Λ` has minimum slack `13.55`, so it could not detect a constant that is wrong by less than
`13×`; the constant's warrant is the symbolic derivation (which I reproduce), not those rows.

---

## 6. FINDING 4 — Lemma T′ is an origin statement; step **P** is a per-shell statement

The seat's statement calls `a_ref` "the quantity step **P** of `prove-lagrangian` consumes".
At `σ = 0` this is exactly right and I confirm it: step P's ODE is `∂_θF(σ,θ) = κH(σ,θ)`,
`H(σ,θ) = ∫_σ¹e^{F}dσ′`, and at `σ=0`, `H(0,θ) = ∫₀¹λ dσ = λ̄ = √6/2`, so
`a = ML·κH(0,θ) = (M/2)Lλ̄ = a_ref` for the cap (`κ = 1/2`). ✓

But step P consumes `H(σ,θ)` at **every** material shell `σ∈[0,1]`, and Lemma T′ bounds only
`a[·](0)` — a functional centred at the origin, whose whole proof (Steps 2–3) rests on `𝒦` being
homogeneous of degree `−4` about the point of evaluation and on `Λ` fixing that point. The
reduction of the shell-`σ` strain to an origin functional over the truncated shell
`{ρ<|x|<R}` is `prove-lagrangian`'s **L3/L3v** (`a = A(ρ) + O(M)`, `A(ρ) = (M/2)log(R/ρ)`) — which
is exactly the piece still SKETCH. Given L3v the extension is free (apply (T′) with `ρ₀ → ρ` and
`λ|_{[σ,1]}`; `μ_J = 2κ_s/L` and `C_rel` are unchanged), so this is an unstated reduction rather
than a hole. It belongs in the gap list; it is not there.

---

## 7. Gate quality (FL-043)

* `check_constants.py`: **553 checks, 553 PASS**, count is `len(CHECKS)`. No check multiplies its
  input by zero. `p6`'s measured mutation coverage **678/678 = 1.0000** reproduces. This is a
  genuine, well-built gate. Two cosmetic notes: `ck('p7[i] |z| all below 3', max|z|, 3.0, 3.0)`
  asserts `max|z| ≤ 6`, not `≤ 3`; and several one-sided clamp checks
  (`min(bound_single_rel,60)==60`, `min(P_at_1,1)==P_at_1`) are individually blind upward but are
  covered by companion recompute checks — which is why coverage is still 1.
* **`p8`'s document audit is 92 % powerful, not 100 %** (`x3`, my measurement, which the seat did
  not make): perturbing each of the 295 numbers in its last printed digit, one at a time,
  `271 caught, 24 missed → detection rate 0.9186`. The misses are the low-precision tokens, where a
  1822-element pool makes an accidental match likely — but they include
  `0.21097` (the tapered `|Δ|`, missed as `0.21090`) and the last digit of the headline
  `0.34913942333376546`. "0 untraced" is therefore weaker evidence than it reads; it should be
  quoted with its power.
* **Controls K1, K2, K3 all genuinely fire and none is vacuous.** K3's map is a diffeomorphism
  (`1+β′ ≥ 1−μ_b > 0`: both ramps have `β′ ≥ 0` because `cos φ` flips sign with the second ramp,
  and in the interior `β′ = −μ_b sin φ`), so the control is legitimate. K1 is close to tautological
  (with `a_ref` defined by the shear-free weight, a `μ`-only bound must fail at `Φ=Λ`) but it is
  the right tautology to record.
* Corollary 2 is labelled **PROVED** but `r_h` is an infimum over a 26-point `λ`-grid of step
  `0.02`, and "the infimum is at `λ = 3/2` in every case" is read off that grid. True on my
  quadrature too, but grid-established, not proved.

---

## 8. Corrected statement

> **LEMMA T′ — unchanged.** With `Λ(x) = T_{λ(|x|)}x`, `λ∈C¹([ρ₀,R];[1,3/2])` obeying
> `|ρλ′/λ| ≤ κ_s/L`, datum `|η₀| ≤ M/r` on `S`, `a_ref[η₀;λ] := ∫_S𝒦(Λx)η₀λ²dx₅`,
> `A := ∫λ dlog ρ`, and `Φ` obeying (H1),(H2),(H3):
> ```
> | a[η₀∘Φ⁻¹](0) − a_ref[η₀;λ] | ≤ π M A [ 3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8 ] .
> ```
> Corollaries 1–4 and the variant (T′_Λ) as written. `C_rel = 3π(√6−2)/2 = 2.1181705106873974`,
> derived. **All PROVED; all re-derived here.**
>
> **Amend the surrounding claims to:**
> 1. The `P_{h_δ}(3/2)` comparison is **resolved**: `prove-lagrangian` and the synthesis tabulate
>    `Φ_h := P_h(λ)/P_h(1)`; normalised, the two seats agree to `3.93e-5`, i.e. to the synthesis's
>    printed precision, at all seven `δ`. Corrections item 4 and gap 5 are withdrawn.
> 2. `p7` is an independent instrument **for the 5-D measure normalisation, the `Λ`-specific
>    reduced kernel, and the quadrature scheme** — not for `det DΛ`, whose closed form it shares
>    verbatim with `lib5`. `det DΛ` is decorrelated instead by `p1`'s two symbolic routes and by
>    two independent finite-difference checks (`r1`, and this seat's).
> 3. The operative constant in every campaign use is `(3π/4)μ_J`; `C(0+) = 15π/8` and `(15π/4)μ`
>    hold only for `μ_J ≤ μ`, a regime realised in **0 of 21** of the seat's own rows.
> 4. Lemma T′ bounds the strain **at the origin**. Step P consumes `H(σ,θ)` at every material shell;
>    the reduction of the shell-`σ` strain to an origin functional over `{ρ<|x|<R}` is `L3/L3v`,
>    still SKETCH. Add this to the gap list.
> 5. "12 non-trivial `Φ`" = 6 maps × 2 `L`. `p8`'s audit power is `0.9186`, measured.

## 9. Remaining gap

Unchanged and correctly stated by the seat: **the PDE half is the whole of what remains of GAP T** —
that the NS flow map on `[0,τ]`, `τ = c/(ML)`, satisfies (H2) and (H3_Λ) against `Λ` with
`μ, μ_J^Λ = O(c/L)` uniformly over the shell including both boundary layers (where
`prove-lagrangian` §4(2) itself records `sup|a−A(ρ)| = 0.7047 M` at the outer edge). It needs
**L3v** (SKETCH) and `Γ = ‖∇u‖_∞` (an unstated hypothesis, per the `gap-V` refuter). Nothing in
either note bears on it. Add to that list: (6) that a **volume-preserving** map exists within
`O(c/L)` of `Λ` (the seat's own gap 2 — `J_Λ = λ²(1+DW) ≠ λ²`, so the shell-dependent family is not
the flow map of an incompressible field); (7) the origin-vs-shell reduction of §6 above.
GAP T is still a gap. **FL-000 stands.** Nothing here touches conjecture (ii) or BFG ARMA 2019
Thm 10.

## 10. Files

| file | what it establishes |
|---|---|
| `x1_independent.py` / `x1_results.json` | `det DΛ` by finite differences of the real 5×5 map; the two shell identities and `Q` by my own `mpmath`; the profile constants from scratch; **`a[η₀∘Λ⁻¹](0)` computed in the `w`-variables, never touching `det DΛ`**; the offset and the conservatism factor |
| `x2_taper_normalisation.py` / `x2_results.json` | Finding 1: `Φ_h := P_h(λ)/P_h(1)`; residual `3.28e-2 → 3.93e-5` |
| `x3_p8_power.py` / `x3_results.json` | measured power of the seat's document audit: `271/295 = 0.9186` |
| `x4_kernel_and_rows.py` / `x4_results.json` | `𝒦` re-derived from the PDE (harmonicity, sign, `a(0)=+M L/2`); Step 0(b)'s `[−1,2]`; all 21 rows re-checked with my own formula; `p7`'s import of `det DΛ` |
| `SHA256SUMS` | computed, never typed |

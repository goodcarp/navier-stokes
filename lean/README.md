# `forced-route-2026-09/lean` — five formalised cores

Five small, self-contained Lean 4 / Mathlib developments backing the exponent, measure and
profile-integral steps of the notes in this repository. Each is sorry-free and axiom-clean.
**None of them formalises a PDE**; read "What is not formalised" below before quoting any of
them.

| File | Backs | Note |
|---|---|---|
| `Cores/CriticalityExponent.lean` | Lemma A, §1.2 and the exponent check in §1.3 | `theorems/01_two_fences.md` |
| `Cores/CircleMeasure.lean` | Lemma B, §2.2(iii) | `theorems/01_two_fences.md` |
| `Cores/CeilingExponent.lean` | the dissipation ceiling, §§2–4 | `analysis/01_dissipation_ceiling_of_the_layered_cascade.md` |
| `Cores/StrainedCore.lean` | the affine core generator, the two radial cancellations, the ratio test | `strained-core/work/pass4/affine-core-pressure.md`, `pass4/core-pressure-derivative.md`, `pass2/active-core-pressure.md` |
| `Cores/Receiver.lean` | stage time and gain, the scale-matching exponents, the Type II window | `strained-core/work/pass3/quantitative-core-stage.md`, `pass3/rotational-receiver.md`, `pass4/iteration-audit.md` |

## Toolchain and pin

* `lean-toolchain`: `leanprover/lean4:v4.33.1`
* Mathlib: `lake-manifest.json`, rev `0df444a360eaa60ab8c11dca51a86af692955474`, `inputRev` `v4.33.1`

Do not run `lake update`; the manifest is the pin.

## How to build

```sh
cd forced-route-2026-09/lean
export PATH="$HOME/.elan/bin:$PATH"
lake exe cache get      # fetch Mathlib oleans; never build Mathlib from source
lake build              # ~2 min from a warm cache
lake env lean Axioms.lean > AXIOMS.txt   # regenerate the axiom audit
```

`lake build` must not start compiling `Mathlib.*` modules. If it does, the olean cache is
missing: stop and run `lake exe cache get` first.

`Axioms.lean` is a standalone script, not part of the `Cores` library; `lake build` ignores
it.

## The exact statements

Everything lives under the namespace `ForcedRoute`.

### `Cores/CriticalityExponent.lean` — namespace `ForcedRoute.Criticality`

`theta α := (2 - α) / α`, and `Vel := EuclideanSpace ℝ (Fin 3)`,
`rescale α lam u p := lam ^ (α - 1) • u (lam • p.1, lam ^ α * p.2)`.

| Statement | Note |
|---|---|
| `theta_mem_Ico (h1 : 1 < α) (h2 : α ≤ 2) : theta α ∈ Ico (0:ℝ) 1` | §1.3 (a) |
| `theta_lt_one_iff (hα : 0 < α) : theta α < 1 ↔ 1 < α` | §1.3 (d) |
| `sobolev_index_interpolation (hα : α ≠ 0) (s : ℝ) : s + 1 - α/2 = (1 - theta α) * s + theta α * (s + α/2)` | §1.3 (b) |
| `two_div_one_sub_theta (h1 : 1 < α) : 2 / (1 - theta α) = α / (α - 1)` | §1.3 (c) |
| `iSup_norm_rescale (hlam : 0 < lam) (u) : (⨆ p, ‖rescale α lam u p‖) = lam ^ (α-1) * ⨆ p, ‖u p‖` | §1.2 (e) |
| `bddAbove_range_norm_rescale`, `iSup_norm_rescale_of_bddAbove` | (e), with the bounded hypothesis |
| `supNorm_scaling_invariant_iff (u) (hpos : 0 < ⨆ p, ‖u p‖) : (∀ lam > 0, ⨆ ‖u_lam‖ = ⨆ ‖u‖) ↔ α = 1` | §1.2, corollary |

Also `theta_nonneg`, `theta_lt_one`, `theta_two : theta 2 = 0`, `theta_one : theta 1 = 1`,
`one_sub_theta : 1 - theta α = 2*(α-1)/α`, `scalingMap_surjective`, `norm_rescale`.

Relation to the note. `theta_lt_one_iff` is the exact algebraic reason Young's inequality
with conjugate exponents `2/(1+θ)`, `2/(1-θ)` is available in Step 2 of §1.3 precisely when
`α > 1` — at `α = 1` one has `θ = 1` (`theta_one`) and the inequality is unavailable.
`two_div_one_sub_theta` is the Grönwall exponent `M^{α/(α-1)}`: `2` at `α = 2`, `3` at
`α = 3/2`, `→ ∞` as `α → 1+`. `supNorm_scaling_invariant_iff` is the sense in which `L^∞` is
critical exactly at `α = 1`.

The identity `iSup_norm_rescale` needs no boundedness hypothesis, because `sSup` of an
unbounded set of reals is the junk value `0` on both sides;
`iSup_norm_rescale_of_bddAbove` carries the boundedness hypothesis and additionally
certifies that both suprema are genuine suprema. It is proved through the bijection
`(x,t) ↦ (λ x, λ^α t)` (`scalingMap_surjective`), which identifies the two ranges of
pointwise norms before the scalar is pulled out.

### `Cores/CircleMeasure.lean` — namespace `ForcedRoute.CircleMeasure`

`Pt := EuclideanSpace ℝ (Fin 3)` with the axis taken to be the `x₂`-axis,
`coord i x := x i`, `circleAt r h := {x | x 2 = h ∧ (x 0)^2 + (x 1)^2 = r^2}`, and

```
rot θ x := !₂[cos θ * x 0 - sin θ * x 1, sin θ * x 0 + cos θ * x 1, x 2]
```

| Statement | Note |
|---|---|
| `lipschitzWith_coord (i) : LipschitzWith 1 (coord i)` | §2.2(iii), the projection |
| `hausdorffMeasure_image_coord_le (i) (s) : μH[1] (coord i '' s) ≤ μH[1] s` | via `LipschitzWith.hausdorffMeasure_image_le` |
| `Icc_subset_coord_image_circleAt (r h) : Icc (-r) r ⊆ coord 0 '' circleAt r h` | the explicit point `(a, √(r²-a²), h)` |
| `hausdorffMeasure_Icc (r) : μH[1] (Icc (-r) r) = ENNReal.ofReal (2*r)` | via `hausdorffMeasure_real` |
| `hausdorffMeasure_circleAt_pos (hr : 0 < r) (h) : 0 < μH[1] (circleAt r h)` | §2.2(iii), the measure step |
| `hausdorffMeasure_circle_pos (hr₀ : 0 < r₀) : 0 < μH[1] {x \| x 2 = 0 ∧ (x 0)^2 + (x 1)^2 = r₀^2}` | the shape used in §2.2 |
| `exists_rot_eq` : rotations act transitively on each circle of positive radius | §2.2(i), the geometric content |
| `circleAt_subset_of_rotInvariant` : a rotation-invariant `S` with a point off the axis contains the whole circle through it | §2.2(iii) |
| `hausdorffMeasure_pos_of_rotInvariant` : hence `0 < μH[1] S` | §2.2(iii), the contradiction |

Also `exists_cos_sin : c² + s² = 1 → ∃ θ, cos θ = c ∧ sin θ = s` (through
`Complex.norm_eq_one_iff`).

Relation to the note. §2.2(iii) argues: `P¹(Σ × {T}) = 0` from CKN; slicing parabolic
cylinders at `t = T` gives `H¹(Σ) = 0`; but a point of `Σ` at distance `r₀ > 0` from the
axis drags the whole circle into `Σ` by rotation invariance, and that circle has positive
`H¹`. The last two steps are what this file proves.

### `Cores/CeilingExponent.lean` — namespace `ForcedRoute.Ceiling`

| Statement | Note |
|---|---|
| `rpow_lt_rpow_exponent_iff (hlam : 1 < lam) (α β Q) : lam^(α*Q) < lam^(β/2) ↔ α*Q < β/2` | §3, the light-damping comparison |
| `damping_below_growth_iff (hlam : 1 < lam) (_hQ : 1 ≤ Q)` : the same, carrying `Q ≥ 1` | §3 |
| `mul_lt_half_iff (hQ : 0 < Q) : α*Q < β/2 ↔ α < β/(2*Q)` | §3, turning it into a ceiling |
| `isGreatest_ceiling (hβ : 0 ≤ β) : IsGreatest {x \| ∃ Q ≥ 1, x = β/(2*Q)} (β/2)` | §4, the optimisation over `Q` |
| `csSup_ceiling (hβ : 0 ≤ β) : sSup {x \| ∃ Q ≥ 1, x = β/(2*Q)} = β/2` | §4 |
| `ceiling_le_half (hβ : β ≤ 1) : β/2 ≤ 1/2` | §4, "supremum `1/2`" |
| `lt_half_of_admissible (0 ≤ β) (β ≤ 1) (1 ≤ Q) (α < β/(2*Q)) : α < 1/2` | §4, assembled |
| `lt_half_of_damping_below_growth` : the same, entered from `lam^(α*Q) < lam^(β/2)` | §§3–4 |

**One correction to the note's phrasing.** `isGreatest_ceiling` requires `0 ≤ β`, and this is
necessary: for `β < 0` the map `Q ↦ β/(2Q)` is increasing on `[1,∞)` and the supremum is `0`,
not `β/2`. The note's `β = 1/8` is positive, so nothing in the note is affected, but the
sentence "sup over `Q ≥ 1` of `β/(2Q)` is `β/2`" is false as stated for negative `β`. The
supremum is also a maximum, attained at `Q = 1`, which is why the file states `IsGreatest`
rather than only `sSup`.

### `Cores/StrainedCore.lean` — namespace `ForcedRoute.StrainedCore`

Backs `strained-core/work/pass4/affine-core-pressure.md` (with
`derive_affine_core_pressure.py`, `verify_affine_pressure_independent.py`),
`pass4/core-pressure-derivative.md`, and `pass2/active-core-pressure.md` (with
`verify_active_core_pressure.py`).

`Lmat b Ω := !![-b, -Ω, 0; Ω, -b, 0; 0, 0, 2*b]`, the affine generator
`L = diag(-b,-b,2b) + Ω J`.

| Statement | Note |
|---|---|
| `trace_Lmat : (Lmat b Ω).trace = 0` | §1, the hypothesis behind `curl[x × (Lx)] = -3Lx` |
| `Lmat_sq_zz : (Lmat b Ω * Lmat b Ω) 2 2 = 4*b^2` | §2 |
| `trace_Lmat_sq : (Lmat b Ω * Lmat b Ω).trace = 6*b^2 - 2*Ω^2` | §2, `g(0) = tr(L²)` |
| `affine_pressure_zz : -(Lmat b Ω * Lmat b Ω) 2 2 = -(4*b^2)` | the **homogeneous** affine `p_zz` |
| `affine_ne_core (b ≠ 0 ∨ Ω ≠ 0) : -(4*b^2) ≠ -(18/7)*b^2 + (2/5)*Ω^2` | the two values differ |

`Cutoff` bundles the radial profile `ψ` of §1 with `ψ'`, `ψ''`, a cutoff radius `T` in the
variable `s = |x|²`, `ψ(0) = 1` and `ψ(T) = ψ'(T) = 0`. Writing `I := ∫₀^T s ψ'(s)² ds`:

| Statement | Note |
|---|---|
| `Cutoff.moment_one : ∫₀^T ψψ' = -1/2` | (7), first line |
| `Cutoff.moment_two : ∫₀^T sψψ'' = 1/2 - I` | (7), second line |
| `Cutoff.moment_three : ∫₀^T s²ψ'ψ'' = -I` | (7), third line |
| `Cutoff.lincomb` : the four-term combination of the three moments | how (7) is used |
| `Cutoff.strain_radial : ½∫ 𝒜_S ds/s = -4/7` | (8) |
| `Cutoff.swirl_radial : ½∫ 𝒜_W ds/s = -4/15` | (9) |
| `Cutoff.core_pressure_zz : -(1/3)tr(L²) + b²·(8) + Ω²·(9) = -(18/7)b² + (2/5)Ω²` | (10) |
| `Cutoff.core_pressure_transverse : (-6b² + 2Ω² - p_zz)/2 = -(12/7)b² + (4/5)Ω²` | (11) |
| `Cutoff.core_pressure_zz_of_pure_rotation` : the `b = 0` case gives `(2/5)Ω²` | §1 of `pass2` |
| `Cutoff.radial_moment : ∫₀^{√T} r ψ(r²)ψ'(r²) dr = -1/4` | `rad` in `verify_active_core_pressure.py` |
| `angular_average : ½∫_{-1}^{1}(3μ²-1)(1-μ²)dμ = -4/15` | `ang` in the same checker |
| `rotating_core_coefficient : 2/3 - 4·(-4/15)·(-1/4) = 2/5` | `core = 2/3 - 4*ang*rad` |
| `ratio_second_derivative_algebraic` | the `deriv.subs(tuning)` line of `verify_affine_pressure_independent.py` |
| `ratio_second_derivative` | the same with real derivatives; `β''(0) = -(p_zz'(0) + 32b³)/(2Ω)` |

Also the plumbing `Cutoff.cont0`, `cont1`, `int1`–`int4`, `I_def`.

Relation to the notes. The three moments are **proved** from the fundamental theorem of
calculus with the antiderivatives `ψ²/2`, `sψψ'` and `s²(ψ')²/2` — they are not assumed — and
`strain_radial`, `swirl_radial` then hold for *every* profile in the class, the coefficient
of `I` cancelling exactly. That cancellation is the whole content of §3's "all radial-profile
dependence cancels". `affine_ne_core` records the point the note makes implicitly: the
compact extension's `-(18/7)b² + (2/5)Ω²` is not the homogeneous affine field's `-4b²`, so
the localisation is not a harmless truncation.

`ratio_second_derivative` takes the note's first-derivative relation

```
β'(t) = (-4 b(t)² - P(t)/2 + η_b(t) - (b(t)/Ω(t)) η_w(t)) / Ω(t)
```

as a **hypothesis** (it comes from the evolution equations (16) plus the profile-defect terms
of `full-feedback-gate-audit.md`, which need the PDE), and proves that at the neutral
insertion `b'(0) = 2b(0)²`, `Ω'(0) = 2b(0)Ω(0)`, `P(0) = -8b(0)²`, `η_b = η_w = η_b' = η_w' =
0`, the chain rule gives exactly `β''(0) = -(P'(0) + 32b(0)³)/(2Ω(0))`. The middle chain-rule
term `4b²/Ω² + p_zz/(2Ω²)` is annihilated precisely by `p_zz = -8b²`.

### `Cores/Receiver.lean` — namespace `ForcedRoute.Receiver`

Backs `strained-core/work/pass3/quantitative-core-stage.md` (with
`verify_quantitative_core_stage.py`), `pass3/rotational-receiver.md` (with
`verify_rotational_receiver.py`), and the exponent reading of `pass4/iteration-audit.md` §3.

`stageTime TE M2 M3 Ω k₀ := min TE (min (k₀/M2) (Ω k₀/M3))`, `gain k₀ t := k₀ t²/4`.

| Statement | Note |
|---|---|
| `stageTime_pos` : `τ > 0` | (6) |
| `gain_pos` : `g > 0`; `one_lt_one_add_gain` : `1 < 1 + g` | (6), (9) |
| `stageTime_mono` : `τ` is nondecreasing in the margin `k₀` | (6) |
| `gain_mono`, `gain_stageTime_mono` : so is `g` | (6) |
| `contraction_sq_mem_Ioo : (1+g/2)/(1+g) ∈ (0,1)`; `contraction_sq_mul` | §5 and the checker |
| `strict_Re_gain_iff (hq : 0 < q) : 1 < q²G ↔ 1/q² < G` | the boxed `G > q^{-2}` |
| `Re_ratio_gt_one : 1 < q^(-a)` | `Re_{qr}(τ)/Re_r(0) = q^{-a}` |
| `vel_ratio_gt_one : 1 < q^(-(1+a))` | `U_{qr}(τ)/U_r(0) = q^{-(1+a)}` |
| `energy_ratio_pos`, `energy_ratio_lt_one : 0 < q^(1-2a) < 1` | `E_{qr}(τ)/E_r(0) = q^{1-2a}` |
| `scale_match_Re : q·q^(-(1+a)) = q^(-a)` | `q*q**(-1-a) - q**(-a)` |
| `scale_match_energy : q^3·q^(-2(1+a)) = q^(1-2a)`, `vel_ratio_sq`, `scale_match_energy'` | `q**3*q**(-2*(1+a)) - q**(1-2*a)` |
| `quarter_exponents`, `quarter_Re_ratio`, `quarter_vel_ratio`, `quarter_energy_ratio` | the `a = 1/4` case: `q^{-1/4}`, `q^{-5/4}`, `q^{1/2}` |
| `bracket_upper : (1+g)^(-1/(a+2)) < 1`; `bracket_lower_pos` | (21a), the contraction bracket |
| `advective_exponent_gt_half_iff (-2 < a) : 1/2 < (1+a)/(2+a) ↔ 0 < a` | Type II end of the window |
| `advective_exponent_lt_one`, `finite_energy_iff : 0 < 1-2a ↔ a < 1/2`, `window`, `window_quarter` | the window |
| `sphere_average_rotational : ½∫_{-1}^{1}(1-μ²)dμ = 2/3` | `N_sphere` in the checker |
| `weighted_projection`, `rigid_rotation_saturates` | the weighted energy lower bound |
| `cubic_remainder_constant : (3M₃/(2Ω))t³/6 = M₃t³/(4Ω)` | the Taylor remainder constant |

Relation to the notes. All of it is arithmetic about the *formulas*. `q` is a hypothesis
here: the note obtains it from the intermediate value theorem applied to a continuous
function built from the actual solution, and that step is not formalised. Likewise `τ` and
`g` are the expressions of (6); that the solution satisfies (7)–(9) on `[0,τ]` is the note's
energy and differentiated-equation argument. `advective_exponent_gt_half_iff` says only that
the exponent `(1+a)/(2+a)` beats the type I rate exactly above `a = 0`; the note's
`T - t_m ∼ L_m^{2+a}` is itself conditional on duration bounds it does not prove.

## Axioms

See `AXIOMS.txt` for the full `#print axioms` printout of all 91 theorems. Every one of them reports

```
depends on axioms: [propext, Classical.choice, Quot.sound]
```

These are the three standard Lean 4 axioms, the axioms of Mathlib itself. `sorryAx` appears
nowhere. There is no `sorry`, no `admit` and no `native_decide` anywhere in `Cores/`, and
`lake build` completes with no warnings.

## What is **not** formalised

This is the important section. The five files are exponent algebra, one scaling identity,
one Hausdorff-measure lower bound, three one-dimensional profile integrals with their exact
linear combinations, and one chain rule. They are not proofs of the lemmas they support.

### The PDE content of Lemma A

Not formalised, at all. Specifically:

* No differential equation appears anywhere in `CriticalityExponent.lean`. `rescale` is
  applied to an arbitrary function `Vel × ℝ → Vel`. The fact that `u_λ` solves `(NS_α)`
  when `u` does — every term picking up `λ^{2α-1}` — is the direct check recorded in §1.2 of
  the note and is **not** verified in Lean.
* The energy estimate (Step 1), the `Ḣ^s` estimate (Step 2), the Kato–Ponce / fractional
  Leibniz product estimate (Grafakos–Oh on `ℝ³`, Bényi–Oh–Zhao on `𝕋³`), the
  Plancherel–Hölder interpolation `‖Λ^{s+θα/2}u‖₂ ≤ ‖Λ^s u‖₂^{1-θ} ‖Λ^{s+α/2}u‖₂^θ`,
  Young's inequality, Grönwall, and the continuation step (Step 3) are all outside Lean.
  What is verified is that the *indices* in those steps are what the note says they are.
* `⨆ p, ‖u p‖` is a supremum of pointwise norms of a bounded field on a bare product space
  with no measure. It is not an `L^∞` norm modulo null sets, and no `L^q` for finite `q`
  appears; `q_c(α) = 3/(α-1)` is not formalised.
* Lemma A's hypotheses — the strong-solution class, the force condition (A1), `s > 5/2`,
  `ν > 0`, divergence-freeness — do not appear.

### CKN and the parabolic measure

* **Mathlib has no parabolic Hausdorff measure `P^d`**, and none is defined here. What
  `CircleMeasure.lean` proves is the **Euclidean `H¹` step only**.
* The comparison used in §2.2(iii) — covering `Σ × {T}` by parabolic cylinders
  `B_{r_i} × (t_i - r_i², t_i)`, slicing at `t = T` to get a Euclidean cover of `Σ` by balls
  of the same radii, hence `P¹(Σ × {T}) = 0 ⇒ H¹(Σ) = 0` — is stated in the note and is
  **not formalised**.
* Caffarelli–Kohn–Nirenberg (1982) is cited, not proved: Main Theorem B (`P¹(S) = 0`), the
  ε-regularity proposition, the definition of a suitable weak solution, the local energy
  inequality with its `2(u·f)φ` term, the Vitali covering argument, and the force
  hypotheses `f ∈ L^q`, `q > 5/2`, `div f = 0`.
* No Navier–Stokes solution appears in `CircleMeasure.lean`. `S` is an arbitrary subset of
  `ℝ³`. That the blow-up set `Σ` of an axisymmetric solution *is* rotation invariant
  (step (i) of §2.2) is not formalised; neither is the non-emptiness of `Σ`, which on `ℝ³`
  needs the extra decay input the note flags rather than papers over.
* The circle's measure is shown only to be **positive**, not equal to `2π r₀`. Positivity is
  all §2.2(iii) uses.
* `rot θ` is defined by its coordinate formula. That it is an isometry, or generates `SO(2)`
  acting on `ℝ³`, is not proved; only transitivity on circles is.
* Tang–Yu, Ren–Wang–Wu, Katz–Pavlović, ESS and Seregin (§§1.5, 2.3) are untouched.

### The mechanism behind the ceiling

* None of the fluid mechanics of `analysis/01_…md`. The amplitude ODE
  `y' = [[-δD, Γ], [Γ, -D]] y` of §2, the eigenvalues `μ₊ = Γ - D` and
  `μ₊ = (-D + √(D² + 4Γ²))/2`, the growth rate `μ = √(-ζ₁(Jζ·G)/|ζ|²)` of §1, the retention
  condition, and the readings of the released schedules (`β = 1/8`, `Q_q = Q_* + q`,
  `q₀ ≥ 2²⁰`) are all outside this file.
* The note itself marks its §§1–4 as *derived* / *checked in a toy*, not as theorems, and
  they rest on Assumption 2 ("nothing else changes") — which §1 and §8 of the note already
  argue is where the released proof would break. Formalising the arithmetic does not upgrade
  the status of the model.
* Nothing here bounds `α` for Navier–Stokes. It bounds the ceiling *of the layered cascade
  as modelled in the note*.

### The strained core and the receiver

* **The pressure solve is assumed, not performed.** `p = N * g` with `N(x) = 1/(4π|x|)`, the
  distributional Hessian `∂_zz N = PV[(3z²-|x|²)/(4π|x|⁵)] - δ₀/3` with its contact term, the
  convergence of the principal value, and the decay condition fixing the constants are all
  outside Lean. The shape of `Cutoff.core_pressure_zz` — that `p_zz(0)` *is* `-(1/3)tr(L²)`
  plus `b²` times the strain radial integral plus `Ω²` times the swirl radial integral — is
  formula (4) of the note entered as the definition of what is being added up.
* **The angular average is assumed.** Equations (5)–(6), that
  `⟨(3μ²-1)g⟩ = b²𝒜_S + Ω²𝒜_W` with the displayed `𝒜_S`, `𝒜_W`, is the symbolic Cartesian
  differentiation in `derive_affine_core_pressure.py` (re-derived cylindrically in
  `verify_affine_pressure_independent.py`). It is **not** verified in Lean; what is verified
  is that *those integrands* integrate to `-4/7` and `-4/15`.
* **The localisation is absent.** `u_L = curl[-(ψ/3)x × (Lx)]` of (1), its
  divergence-freeness, formula (2), the poloidal/azimuthal split (3), the absence of a `bΩ`
  cross term in `g`, and the compact packet `v` with `C_v > 0` are not formalised. No vector
  field appears anywhere in `StrainedCore.lean`.
* **`∫₀^∞` is replaced by a compact endpoint, and `C_c^∞` by `C²`.** `Cutoff` carries a
  radius `T` with `ψ(T) = ψ'(T) = 0`, and integrates over `[0,T]`. For the note's compactly
  supported `ψ` this is the same number, but the improper integral is not formalised. Only
  two derivatives are asked for, because §3 uses only two.
* **The ratio test's `β'` relation is a hypothesis**, transcribed from
  `full-feedback-gate-audit.md`; deriving it needs the PDE. What is proved is that *given*
  it, the chain rule yields the note's `β''(0) = -(p_zz'(0) + 32b³)/(2Ω)`. **No sign is
  asserted**: `P'(0)` is neither evaluated nor bounded, so nothing about the cone follows.
* **The numerical lead is untouched.** `t₃₀₀ ≈ -6.7344633112` and
  `t_Wb ≈ 2.563122131694` of §5 of `core-pressure-derivative.md` are floating-point
  quadrature for one logistic cutoff. They need the cumulative integrals `Q(s) = ∫₀^s ψ` and
  `K(s) = ∫₀^s τ^{7/2}[ψψ' + (7/9)τ(ψ')²]`, half-integer powers, and a transcendental
  profile; they are **not** rational and **not** profile-independent — the note itself
  exhibits two polynomial proxies with different values. Nothing about them, nor about
  `d_S = d_W = 0`, the coefficient formulas (1)–(2), or the harmonic expansion (3) of that
  note, is formalised.
* **No solution, no receiver, no gain.** `Receiver.lean` contains no PDE, no `a_r(t)`, no
  `E_r`, `U_r`, `Re_r`. The "ratios" are the real numbers `q^{-a}`, `q^{-(1+a)}`, `q^{1-2a}`;
  that they equal observable ratios of one Navier–Stokes solution is the content of
  `rotational-receiver.md` and rests on local smooth existence, the initial pressure
  calculation, the harmonic first-moment identity and an IVT. **The existence of the root `q`
  is not formalised** — `q` is a hypothesis throughout.
* **`τ` and `g` are formulas only.** That the solution obeys (7), (8), (9) on `[0,τ]` needs
  the explicit Sobolev product constants `C₄, C₆, C₈, C_E`, the lifespan estimate and the
  third-derivative bound `M₃`; none of that is here, and `M₂`, `M₃` are just positive reals.
  The counting inequalities `9·2^{2r}·C(r+3,3) < C_r²` of the checker are omitted rather than
  formalised in isolation.
* **The Type II reading is conditional in the note and conditional here.** `T - t_m ∼
  L_m^{2+a}` holds in `iteration-audit.md` only if a separate stage theorem supplies the
  duration bounds and a uniform contraction; the note says so. `window` says two rational
  inequalities are simultaneously satisfiable on `0 < a < 1/2`, not that any solution
  realises them. Nothing here is a blow-up statement, and nothing here weakens FL-000.

### Also absent

* The Alpöge–Buckmaster preprints, which the note records as unverified ("as reported to
  me"), play no role here.
* No claim about Fefferman's (C) or (D) is formalised, and none follows from these files.

import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# The exponent inequality of the dissipation ceiling

This file formalises the **exponent arithmetic** of
`analysis/01_dissipation_ceiling_of_the_layered_cascade.md`, §§2-4: the comparison between a
layer's damping rate `D_q = ν λ_q^α |ζ_q|^α` and its growth rate `Γ_q ∼ λ^(β/2)`, and the
supremum of the resulting ceiling over the frequency-ratio schedule.

In the note, layer `q` sits at frequency `λ_q = λ_{q-1}^{Q_q}`, the growth-rate exponent is
`β/2` (with `β = 1/8` in both released schedules), and the light-damping condition
`D_q ≲ Γ_q` reads `λ^(α Q) < λ^(β/2)`, i.e. `α Q < β/2`, i.e. `α < β/(2Q)`.  Since `Q ≥ 1`
in every schedule, the largest ceiling is at `Q = 1`, giving `α < β/2 ≤ 1/2` for `β ≤ 1`.
The `1/2` is the amplitude-ODE supremum quoted in §4 for the both-fields-damped case, and it
sits a factor two below the bounded-velocity wall `α = 1` of Lemma A — the square root in the
mechanism.

## What is formalised

* `rpow_lt_rpow_exponent_iff` and `damping_below_growth_iff` : for `λ > 1` and reals
  `α, β, Q`, `λ^(α Q) < λ^(β/2) ↔ α Q < β/2`.  The second is the same statement carrying the
  `Q ≥ 1` side condition of the note, which the equivalence itself does not need.
* `mul_lt_half_iff` : for `Q > 0`, `α Q < β/2 ↔ α < β/(2Q)`.  This is the step that turns
  the frequency comparison into a ceiling on `α`.
* `isGreatest_ceiling` : for `β ≥ 0`, `β/2` is the **greatest** element of
  `{β/(2Q) : Q ≥ 1}` (attained at `Q = 1`), and `csSup_ceiling` : the supremum is `β/2`.
* `ceiling_le_half` : `β ≤ 1 ⇒ β/2 ≤ 1/2`.
* `lt_half_of_admissible` : combining these, if `0 ≤ β ≤ 1`, `Q ≥ 1` and `α < β/(2Q)`, then
  `α < 1/2`.

## What is NOT formalised

* **None of the fluid mechanics.**  The amplitude ODE `y' = [[-δD, Γ], [Γ, -D]] y` of §2,
  the eigenvalue computations `μ₊ = Γ - D` (both fields damped) and
  `μ₊ = (-D + √(D² + 4Γ²))/2` (momentum only), the identification of the growth rate
  `μ = √(-ζ₁ (Jζ·G)/|ζ|²)` in §1, the retention condition, and the reading of the released
  schedules (`β = 1/8`, `Q_q = Q_* + q`, `q₀ ≥ 2²⁰`) are all outside this file.  Only the
  inequality `λ^(αQ) < λ^(β/2) ↔ αQ < β/2` and the optimisation over `Q` are here.
* The note's own labelling applies: its §§1-4 are marked *derived* / *checked in a toy*, not
  theorems, and rest on the idealisation of Assumption 2 ("nothing else changes"), which §1
  and §8 of the note already argue is where the released proof would break.
* The hypothesis `0 ≤ β` in `isGreatest_ceiling` is **necessary and is an addition to the
  note's phrasing**: for `β < 0` the map `Q ↦ β/(2Q)` is increasing on `[1, ∞)` and its
  supremum is `0`, not `β/2`.  The note's `β = 1/8` is positive, so nothing is lost, but the
  claim "sup over `Q ≥ 1` of `β/(2Q)` is `β/2`" is false as stated for negative `β`.
* Nothing here bounds `α` for Navier-Stokes.  It bounds the ceiling *of the layered cascade
  as modelled in the note*; Lemma A of `theorems/01_two_fences.md` is the statement about the
  equation.
-/

namespace ForcedRoute
namespace Ceiling

open Set

/-- For a base `λ > 1`, comparing powers is comparing exponents. -/
theorem rpow_lt_rpow_exponent_iff {lam : ℝ} (hlam : 1 < lam) (α β Q : ℝ) :
    lam ^ (α * Q) < lam ^ (β / 2) ↔ α * Q < β / 2 :=
  Real.rpow_lt_rpow_left_iff hlam

/-- The light-damping comparison of §3, with the `Q ≥ 1` side condition of the note carried
explicitly.  The equivalence does not use `Q ≥ 1`; the optimisation in `isGreatest_ceiling`
does. -/
theorem damping_below_growth_iff {lam α β Q : ℝ} (hlam : 1 < lam) (_hQ : 1 ≤ Q) :
    lam ^ (α * Q) < lam ^ (β / 2) ↔ α * Q < β / 2 :=
  Real.rpow_lt_rpow_left_iff hlam

/-- `α Q < β/2` is the ceiling `α < β/(2Q)`. -/
theorem mul_lt_half_iff {α β Q : ℝ} (hQ : 0 < Q) : α * Q < β / 2 ↔ α < β / (2 * Q) := by
  have hQ' : Q ≠ 0 := ne_of_gt hQ
  have hkey : β / (2 * Q) - α = (β / 2 - α * Q) / Q := by field_simp
  constructor
  · intro h
    have hy : (0 : ℝ) < (β / 2 - α * Q) / Q := div_pos (by linarith) hQ
    rw [← hkey] at hy
    linarith
  · intro h
    have hy : (0 : ℝ) < β / (2 * Q) - α := by linarith
    rw [hkey] at hy
    have hz := mul_pos hy hQ
    rw [div_mul_cancel₀ _ hQ'] at hz
    linarith

/-- **The ceiling.**  For `β ≥ 0`, the family `{β/(2Q) : Q ≥ 1}` has greatest element `β/2`,
attained at `Q = 1`.

The hypothesis `0 ≤ β` cannot be dropped: for `β < 0` the supremum is `0`. -/
theorem isGreatest_ceiling {β : ℝ} (hβ : 0 ≤ β) :
    IsGreatest {x : ℝ | ∃ Q : ℝ, 1 ≤ Q ∧ x = β / (2 * Q)} (β / 2) := by
  constructor
  · exact ⟨1, le_refl 1, by norm_num⟩
  · rintro x ⟨Q, hQ, rfl⟩
    have hQ0 : (0 : ℝ) < Q := by linarith
    have hQ' : Q ≠ 0 := ne_of_gt hQ0
    have hkey : β / 2 - β / (2 * Q) = β * (Q - 1) / (2 * Q) := by field_simp
    have hnn : 0 ≤ β * (Q - 1) / (2 * Q) :=
      div_nonneg (mul_nonneg hβ (by linarith)) (by linarith)
    linarith [hkey ▸ hnn]

/-- The supremum form: `sup_{Q ≥ 1} β/(2Q) = β/2` for `β ≥ 0`. -/
theorem csSup_ceiling {β : ℝ} (hβ : 0 ≤ β) :
    sSup {x : ℝ | ∃ Q : ℝ, 1 ≤ Q ∧ x = β / (2 * Q)} = β / 2 :=
  (isGreatest_ceiling hβ).csSup_eq

/-- `β ≤ 1` puts the ceiling at or below `1/2`. -/
theorem ceiling_le_half {β : ℝ} (hβ : β ≤ 1) : β / 2 ≤ 1 / 2 := by linarith

/-- **The ceiling of the layered cascade, assembled.**  With `0 ≤ β ≤ 1` (the note has
`β = 1/8`) and any schedule exponent `Q ≥ 1`, every admissible dissipation order
`α < β/(2Q)` satisfies `α < 1/2`. -/
theorem lt_half_of_admissible {α β Q : ℝ} (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hQ : 1 ≤ Q)
    (h : α < β / (2 * Q)) : α < 1 / 2 := by
  have hub : β / (2 * Q) ≤ β / 2 := (isGreatest_ceiling hβ0).2 ⟨Q, hQ, rfl⟩
  have := ceiling_le_half hβ1
  linarith

/-- The same, phrased against the frequency comparison it comes from: if a layer at frequency
ratio `Q ≥ 1` is lightly damped, `λ^(α Q) < λ^(β/2)`, then `α < 1/2`. -/
theorem lt_half_of_damping_below_growth {lam α β Q : ℝ} (hlam : 1 < lam) (hβ0 : 0 ≤ β)
    (hβ1 : β ≤ 1) (hQ : 1 ≤ Q) (h : lam ^ (α * Q) < lam ^ (β / 2)) : α < 1 / 2 := by
  have hQ0 : (0 : ℝ) < Q := by linarith
  have h1 : α * Q < β / 2 := (rpow_lt_rpow_exponent_iff hlam α β Q).1 h
  exact lt_half_of_admissible hβ0 hβ1 hQ ((mul_lt_half_iff hQ0).1 h1)

end Ceiling
end ForcedRoute

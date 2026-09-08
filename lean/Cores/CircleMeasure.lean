import Mathlib.MeasureTheory.Measure.Hausdorff
import Mathlib.Analysis.Normed.Lp.MeasurableSpace
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.SpecialFunctions.Complex.Arg
import Mathlib.Analysis.Complex.Trigonometric

/-!
# A circle about the axis has positive one-dimensional Hausdorff measure

This file formalises the **measure step** of Lemma B, §2.2(iii) of
`theorems/01_two_fences.md`: the final contradiction in the proof that the blow-up set of an
axisymmetric solution of full Navier-Stokes lies on the axis.

The argument of §2.2 runs: CKN gives `P¹(Σ × {T}) = 0` for the parabolic Hausdorff measure;
slicing a parabolic cover at `t = T` gives a Euclidean cover of `Σ` by balls of the same
radii, so `H¹(Σ) = 0`; but if `Σ` contains a point at distance `r₀ > 0` from the axis then,
by rotation invariance, `Σ` contains the whole circle of radius `r₀` about the axis, and
`H¹(circle) = 2π r₀ > 0`.  Contradiction.

## What is formalised

Everything is in `EuclideanSpace ℝ (Fin 3)` with the axis of symmetry taken to be the
`x₂`-axis.

* `lipschitzWith_coord` : the coordinate map `x ↦ x i` is `1`-Lipschitz.
* `hausdorffMeasure_circleAt_pos` : for `r > 0` and any height `h`,
  `μH[1] {x | x 2 = h ∧ (x 0)² + (x 1)² = r²} > 0`.
  Route, exactly as instructed: `μH[1](proj '' C) ≤ μH[1](C)` by
  `LipschitzWith.hausdorffMeasure_image_le` with `K = 1`; `proj '' C ⊇ Icc (-r) r` by
  exhibiting the explicit point `(a, √(r² - a²), h)` on the circle; and `μH[1]` on `ℝ` is
  Lebesgue measure (`hausdorffMeasure_real`), so `μH[1](Icc (-r) r) = 2r > 0`.
* `hausdorffMeasure_circle_pos` : the statement in the exact shape used in §2.2, with
  `h = 0`.
* `rot` : the rotation by angle `θ` about the `x₂`-axis, **defined explicitly on
  coordinates**, and `exists_rot_eq` : the rotations act transitively on each circle
  `{x 2 = h, (x 0)² + (x 1)² = r² > 0}`.
* `circleAt_subset_of_rotInvariant` : if `S` is invariant under every `rot θ` and contains a
  point off the axis, then `S` contains the entire circle through that point.
* `hausdorffMeasure_pos_of_rotInvariant` : hence `μH[1] S > 0`.  This is the contradiction
  step of §2.2(iii).

## What is NOT formalised — read this before quoting the file

* **Mathlib has no parabolic Hausdorff measure `P^d`.**  What is proved here is the
  *Euclidean* `H¹` step only.  The comparison used in §2.2(iii) — that a cover of
  `Σ × {T}` by parabolic cylinders `B_{r_i} × (t_i - r_i², t_i)`, sliced at `t = T`, yields a
  Euclidean cover of `Σ` by balls of the same radii, so that `P¹(Σ × {T}) = 0 ⇒ H¹(Σ) = 0` —
  is stated in the note and **is not formalised here**.
* **CKN is not formalised.**  The partial-regularity input `P¹(S) = 0`
  (Caffarelli-Kohn-Nirenberg 1982, Main Theorem B), its `ε`-regularity Proposition, the
  suitable-weak-solution hypotheses and the local energy inequality, the Vitali covering
  argument, and the force hypotheses `f ∈ L^q`, `q > 5/2`, `div f = 0`, are all cited, not
  proved.
* **No Navier-Stokes solution appears in this file.**  `S` is an arbitrary subset of `ℝ³`;
  the fact that the blow-up set `Σ` of an axisymmetric solution *is* rotation invariant
  (step (i) of §2.2) is not formalised, nor is the non-emptiness of `Σ`, which on `ℝ³`
  needs an extra decay input flagged in the note.
* The measure computed for a circle here is only shown to be **positive**, not equal to
  `2π r₀`.  Positivity is all §2.2(iii) uses.
* `rot θ` is defined by its coordinate formula.  That it is an isometry, or that it
  generates `SO(2)` acting on `ℝ³`, is not proved; only transitivity on circles is.
-/

namespace ForcedRoute
namespace CircleMeasure

open MeasureTheory Set
open scoped ENNReal NNReal

noncomputable section

/-- Points of `ℝ³`; the axis of symmetry is the `x₂`-axis. -/
abbrev Pt : Type := EuclideanSpace ℝ (Fin 3)

/-- The `i`-th coordinate map. -/
def coord (i : Fin 3) (x : Pt) : ℝ := x i

/-- Coordinate projections of `ℝ³` are `1`-Lipschitz. -/
theorem lipschitzWith_coord (i : Fin 3) : LipschitzWith 1 (coord i) := by
  apply LipschitzWith.of_dist_le_mul
  intro x y
  simpa [coord] using PiLp.dist_apply_le x y i

/-- A `1`-Lipschitz image cannot increase `μH[1]`. -/
theorem hausdorffMeasure_image_coord_le (i : Fin 3) (s : Set Pt) :
    μH[(1 : ℝ)] (coord i '' s) ≤ μH[(1 : ℝ)] s := by
  simpa using (lipschitzWith_coord i).hausdorffMeasure_image_le zero_le_one s

/-- The circle of radius `r` about the `x₂`-axis, at height `h`. -/
def circleAt (r h : ℝ) : Set Pt := {x | x 2 = h ∧ (x 0) ^ 2 + (x 1) ^ 2 = r ^ 2}

/-- The first coordinate of the circle of radius `r` sweeps at least `[-r, r]`. -/
theorem Icc_subset_coord_image_circleAt (r h : ℝ) :
    Icc (-r) r ⊆ coord 0 '' circleAt r h := by
  rintro a ⟨ha₁, ha₂⟩
  have hnn : 0 ≤ r ^ 2 - a ^ 2 := by nlinarith
  refine ⟨!₂[a, Real.sqrt (r ^ 2 - a ^ 2), h], ⟨?_, ?_⟩, ?_⟩
  · simp
  · have e0 : (!₂[a, Real.sqrt (r ^ 2 - a ^ 2), h] : Pt) 0 = a := by simp
    have e1 : (!₂[a, Real.sqrt (r ^ 2 - a ^ 2), h] : Pt) 1 = Real.sqrt (r ^ 2 - a ^ 2) := by
      simp
    show (!₂[a, Real.sqrt (r ^ 2 - a ^ 2), h] : Pt) 0 ^ 2
        + (!₂[a, Real.sqrt (r ^ 2 - a ^ 2), h] : Pt) 1 ^ 2 = r ^ 2
    rw [e0, e1, Real.sq_sqrt hnn]
    ring
  · simp [coord]

/-- `μH[1]` on `ℝ` is Lebesgue measure, so an interval has the expected measure. -/
theorem hausdorffMeasure_Icc (r : ℝ) :
    μH[(1 : ℝ)] (Icc (-r) r) = ENNReal.ofReal (2 * r) := by
  rw [MeasureTheory.hausdorffMeasure_real, Real.volume_Icc]
  ring_nf

/-- **The measure step of Lemma B.**  A circle of positive radius about the axis has positive
one-dimensional Hausdorff measure in `ℝ³`. -/
theorem hausdorffMeasure_circleAt_pos {r : ℝ} (hr : 0 < r) (h : ℝ) :
    0 < μH[(1 : ℝ)] (circleAt r h) := by
  have h1 : μH[(1 : ℝ)] (Icc (-r) r) ≤ μH[(1 : ℝ)] (coord 0 '' circleAt r h) :=
    measure_mono (Icc_subset_coord_image_circleAt r h)
  have h2 : μH[(1 : ℝ)] (coord 0 '' circleAt r h) ≤ μH[(1 : ℝ)] (circleAt r h) :=
    hausdorffMeasure_image_coord_le 0 (circleAt r h)
  have h3 : (0 : ℝ≥0∞) < ENNReal.ofReal (2 * r) := ENNReal.ofReal_pos.2 (by linarith)
  calc (0 : ℝ≥0∞) < ENNReal.ofReal (2 * r) := h3
    _ = μH[(1 : ℝ)] (Icc (-r) r) := (hausdorffMeasure_Icc r).symm
    _ ≤ μH[(1 : ℝ)] (coord 0 '' circleAt r h) := h1
    _ ≤ μH[(1 : ℝ)] (circleAt r h) := h2

/-- The statement in the exact shape used in §2.2(iii). -/
theorem hausdorffMeasure_circle_pos {r₀ : ℝ} (hr₀ : 0 < r₀) :
    0 < μH[(1 : ℝ)] {x : Pt | x 2 = 0 ∧ (x 0) ^ 2 + (x 1) ^ 2 = r₀ ^ 2} :=
  hausdorffMeasure_circleAt_pos hr₀ 0

/-! ### Rotation invariance about the axis -/

/-- Rotation by the angle `θ` about the `x₂`-axis, written out on coordinates. -/
def rot (θ : ℝ) (x : Pt) : Pt :=
  !₂[Real.cos θ * x 0 - Real.sin θ * x 1, Real.sin θ * x 0 + Real.cos θ * x 1, x 2]

/-- Any point of the unit circle in `ℝ²` is `(cos θ, sin θ)` for some `θ`. -/
theorem exists_cos_sin {c s : ℝ} (h : c ^ 2 + s ^ 2 = 1) :
    ∃ θ : ℝ, Real.cos θ = c ∧ Real.sin θ = s := by
  have hz : ‖(⟨c, s⟩ : ℂ)‖ = 1 := by
    rw [Complex.norm_def, Complex.normSq_mk, show c * c + s * s = 1 by nlinarith]
    exact Real.sqrt_one
  obtain ⟨θ, hθ⟩ := (Complex.norm_eq_one_iff _).1 hz
  refine ⟨θ, ?_, ?_⟩
  · have := congrArg Complex.re hθ; rwa [Complex.exp_ofReal_mul_I_re] at this
  · have := congrArg Complex.im hθ; rwa [Complex.exp_ofReal_mul_I_im] at this

/-- **Transitivity.**  Rotations about the `x₂`-axis act transitively on each circle of
positive radius at a fixed height. -/
theorem exists_rot_eq {p q : Pt} (hp : 0 < (p 0) ^ 2 + (p 1) ^ 2)
    (hr : (q 0) ^ 2 + (q 1) ^ 2 = (p 0) ^ 2 + (p 1) ^ 2) (hz : q 2 = p 2) :
    ∃ θ : ℝ, rot θ p = q := by
  have hne : (p 0) ^ 2 + (p 1) ^ 2 ≠ 0 := ne_of_gt hp
  have key : ((p 0) * (q 0) + (p 1) * (q 1)) ^ 2 + ((p 0) * (q 1) - (p 1) * (q 0)) ^ 2
      = ((p 0) ^ 2 + (p 1) ^ 2) ^ 2 := by
    linear_combination ((p 0) ^ 2 + (p 1) ^ 2) * hr
  have hcs : (((p 0) * (q 0) + (p 1) * (q 1)) / ((p 0) ^ 2 + (p 1) ^ 2)) ^ 2
      + (((p 0) * (q 1) - (p 1) * (q 0)) / ((p 0) ^ 2 + (p 1) ^ 2)) ^ 2 = 1 := by
    rw [div_pow, div_pow, ← add_div, key, div_self (pow_ne_zero 2 hne)]
  obtain ⟨θ, hcos, hsin⟩ := exists_cos_sin hcs
  refine ⟨θ, ?_⟩
  ext i
  fin_cases i
  · show Real.cos θ * p 0 - Real.sin θ * p 1 = q 0
    rw [hcos, hsin]; field_simp; ring
  · show Real.sin θ * p 0 + Real.cos θ * p 1 = q 1
    rw [hcos, hsin]; field_simp; ring
  · show p 2 = q 2
    exact hz.symm

/-- **Rotation invariance forces a circle.**  If `S` is invariant under every rotation about
the `x₂`-axis and contains a point off the axis, then `S` contains the whole circle through
that point. -/
theorem circleAt_subset_of_rotInvariant {S : Set Pt}
    (hS : ∀ θ : ℝ, ∀ x ∈ S, rot θ x ∈ S) {p : Pt} (hpS : p ∈ S)
    (hp : 0 < (p 0) ^ 2 + (p 1) ^ 2) :
    circleAt (Real.sqrt ((p 0) ^ 2 + (p 1) ^ 2)) (p 2) ⊆ S := by
  rintro q ⟨hq2, hqr⟩
  have hr : (q 0) ^ 2 + (q 1) ^ 2 = (p 0) ^ 2 + (p 1) ^ 2 := by
    rw [hqr, Real.sq_sqrt hp.le]
  obtain ⟨θ, hθ⟩ := exists_rot_eq hp hr hq2
  exact hθ ▸ hS θ p hpS

/-- **The contradiction step of §2.2(iii).**  A rotation-invariant set containing a point off
the axis has positive one-dimensional Hausdorff measure. -/
theorem hausdorffMeasure_pos_of_rotInvariant {S : Set Pt}
    (hS : ∀ θ : ℝ, ∀ x ∈ S, rot θ x ∈ S) {p : Pt} (hpS : p ∈ S)
    (hp : 0 < (p 0) ^ 2 + (p 1) ^ 2) :
    0 < μH[(1 : ℝ)] S :=
  lt_of_lt_of_le (hausdorffMeasure_circleAt_pos (Real.sqrt_pos.2 hp) (p 2))
    (measure_mono (circleAt_subset_of_rotInvariant hS hpS hp))

end

end CircleMeasure
end ForcedRoute

import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Data.Real.Pointwise

/-!
# The criticality exponent `θ = (2 - α)/α`, and `L^∞` scaling

This file formalises the **exponent algebra** and the **scaling identity** that sit under
Lemma A of `theorems/01_two_fences.md` (§1.2 and the "exponent check" inside §1.3), for the
hypodissipative / hyperdissipative Navier-Stokes family

    ∂_t u + (u·∇)u + ∇p = -ν Λ^α u + f,   div u = 0,   Λ = |∇|.

Write `θ(α) := (2 - α)/α`.

## What is formalised

* `theta_mem_Ico` : for `1 < α ≤ 2`, `θ ∈ [0, 1)`.  This is item (a).
* `theta_lt_one_iff` : for `α > 0`, `θ < 1 ↔ α > 1`.  This is item (d), and it is the exact
  algebraic reason Young's inequality with conjugate exponents `2/(1+θ)`, `2/(1-θ)` is
  available in §1.3 precisely when `α > 1`.
* `sobolev_index_interpolation` : for every real `s`,
  `s + 1 - α/2 = (1-θ)·s + θ·(s + α/2)`.  This is item (b): the index `s + 1 - α/2` appearing
  in the commutator/Kato-Ponce term is the `θ`-interpolant between the energy index `s` and
  the dissipation index `s + α/2`.
* `two_div_one_sub_theta` : `2/(1-θ) = α/(α-1)`.  This is item (c): the Grönwall exponent on
  `‖u‖_{L^∞}` in Step 2 of §1.3.
* `iSup_norm_rescale`, `bddAbove_range_norm_rescale`, `iSup_norm_rescale_of_bddAbove` :
  the sup-norm scaling identity.  For `λ > 0` and `u_λ(x,t) := λ^(α-1) • u(λ·x, λ^α t)`,
  `⨆ ‖u_λ‖ = λ^(α-1) · ⨆ ‖u‖`.  This is item (e).
* `supNorm_scaling_invariant_iff` : if `⨆ ‖u‖ > 0`, then `⨆ ‖u_λ‖ = ⨆ ‖u‖` for every `λ > 0`
  iff `α = 1`.  This is the precise sense in which `L^∞` is the critical Lebesgue space
  exactly at `α = 1` (`q_c(α) = 3/(α-1) → ∞` as `α → 1+`).

## What is NOT formalised

Everything analytic.  In particular:

* No PDE appears in this file.  `rescale` is applied to an arbitrary function
  `u : E × ℝ → E`; the fact that `u_λ` solves the same equation when `u` does (every term of
  `(NS_α)` picking up `λ^{2α-1}`) is *not* proved here — it is the direct check recorded in
  §1.2 of the note.
* Lemma A itself — the energy/Grönwall argument of §1.3, the Kato-Ponce product estimate
  (Grafakos-Oh on `ℝ³`, Bényi-Oh-Zhao on `𝕋³`), the Plancherel-Hölder interpolation
  `‖Λ^{s+θα/2}u‖₂ ≤ ‖Λ^s u‖₂^{1-θ} ‖Λ^{s+α/2}u‖₂^θ`, Young's inequality, and the
  continuation step — is **not** formalised.  Only the exponent bookkeeping is.
* The domain here is a bare `EuclideanSpace ℝ (Fin 3) × ℝ` with no measure, no
  divergence-free condition, and no regularity: `⨆ ‖u‖` is the sup of the pointwise norm of a
  bounded field, not an `L^∞` norm modulo null sets.
* Nothing about `q_c(α) = 3/(α-1)` for finite `q` is formalised; only the `q = ∞` endpoint.
-/

namespace ForcedRoute
namespace Criticality

open Set

noncomputable section

/-- The interpolation exponent `θ(α) = (2 - α)/α` of §1.3. -/
def theta (α : ℝ) : ℝ := (2 - α) / α

/-! ### (a), (d): the range of `θ` -/

theorem theta_nonneg {α : ℝ} (hα : 0 < α) (h2 : α ≤ 2) : 0 ≤ theta α :=
  div_nonneg (by linarith) hα.le

/-- **(d)** For `α > 0`, `θ(α) < 1` exactly when `α > 1`.  This is the threshold that makes
Young's inequality available in Step 2 of §1.3. -/
theorem theta_lt_one_iff {α : ℝ} (hα : 0 < α) : theta α < 1 ↔ 1 < α := by
  rw [theta, div_lt_one hα]
  constructor <;> intro h <;> linarith

theorem theta_lt_one {α : ℝ} (h1 : 1 < α) : theta α < 1 :=
  (theta_lt_one_iff (by linarith)).2 h1

/-- **(a)** For `1 < α ≤ 2` one has `0 ≤ θ(α) < 1`. -/
theorem theta_mem_Ico {α : ℝ} (h1 : 1 < α) (h2 : α ≤ 2) : theta α ∈ Ico (0 : ℝ) 1 :=
  ⟨theta_nonneg (by linarith) h2, theta_lt_one h1⟩

/-- At the endpoint `α = 2` the exponent degenerates to `0`. -/
theorem theta_two : theta 2 = 0 := by norm_num [theta]

/-- At `α = 1` the exponent hits the forbidden value `1`. -/
theorem theta_one : theta 1 = 1 := by norm_num [theta]

/-! ### (b): the index `s + 1 - α/2` is the `θ`-interpolant -/

/-- **(b)** For every real `s`, the Kato-Ponce index `s + 1 - α/2` is the convex combination,
with weight `θ(α)`, of the energy index `s` and the dissipation index `s + α/2`. -/
theorem sobolev_index_interpolation {α : ℝ} (hα : α ≠ 0) (s : ℝ) :
    s + 1 - α / 2 = (1 - theta α) * s + theta α * (s + α / 2) := by
  rw [theta]
  field_simp
  ring

/-! ### (c): the Grönwall exponent -/

theorem one_sub_theta {α : ℝ} (hα : α ≠ 0) : 1 - theta α = 2 * (α - 1) / α := by
  rw [theta]
  field_simp
  ring

/-- **(c)** `2/(1 - θ) = α/(α - 1)`: the power of `‖u‖_{L^∞}` in the Grönwall bound of §1.3.
At `α = 2` this is `2`, at `α = 3/2` it is `3`, and it blows up as `α → 1+`. -/
theorem two_div_one_sub_theta {α : ℝ} (h1 : 1 < α) : 2 / (1 - theta α) = α / (α - 1) := by
  have hα : α ≠ 0 := by intro h; rw [h] at h1; linarith
  have hα1 : α - 1 ≠ 0 := by intro h; rw [sub_eq_zero] at h; rw [h] at h1; linarith
  rw [one_sub_theta hα]
  field_simp

example : 2 / (1 - theta 2) = 2 := by norm_num [theta]

/-! ### (e): the sup-norm scaling lemma -/

/-- The velocity values live in `ℝ³`. -/
abbrev Vel : Type := EuclideanSpace ℝ (Fin 3)

/-- The parabolic rescaling of §1.2: `u_λ(x,t) = λ^(α-1) u(λ x, λ^α t)`.

No PDE is involved: this is applied to an arbitrary field. -/
def rescale (α lam : ℝ) (u : Vel × ℝ → Vel) : Vel × ℝ → Vel :=
  fun p => (lam ^ (α - 1)) • u (lam • p.1, lam ^ α * p.2)

/-- The space-time dilation `(x,t) ↦ (λ x, λ^α t)` is a bijection for `λ > 0`. -/
theorem scalingMap_surjective {α lam : ℝ} (hlam : 0 < lam) :
    Function.Surjective (fun p : Vel × ℝ => (lam • p.1, lam ^ α * p.2)) := by
  rintro ⟨y, s⟩
  have h1 : lam ≠ 0 := ne_of_gt hlam
  have h2 : lam ^ α ≠ 0 := ne_of_gt (Real.rpow_pos_of_pos hlam α)
  exact ⟨(lam⁻¹ • y, (lam ^ α)⁻¹ * s), by
    simp [smul_smul, mul_inv_cancel₀ h1, mul_inv_cancel_left₀ h2]⟩

theorem norm_rescale {α lam : ℝ} (hlam : 0 < lam) (u : Vel × ℝ → Vel) (p : Vel × ℝ) :
    ‖rescale α lam u p‖ = lam ^ (α - 1) * ‖u (lam • p.1, lam ^ α * p.2)‖ := by
  rw [rescale, norm_smul, Real.norm_eq_abs, abs_of_pos (Real.rpow_pos_of_pos hlam _)]

/-- **(e)** The sup of `‖u_λ‖` is `λ^(α-1)` times the sup of `‖u‖`.

Proved through the bijection `(x,t) ↦ (λ x, λ^α t)` of `scalingMap_surjective`, so the two
ranges of pointwise norms coincide before the scalar is pulled out.

No boundedness hypothesis is needed for the *identity*, because `sSup` of an unbounded set of
reals is the junk value `0` on both sides.  For the identity to say what it is meant to say,
use `iSup_norm_rescale_of_bddAbove`, which also certifies that the two suprema are genuine
suprema. -/
theorem iSup_norm_rescale {α lam : ℝ} (hlam : 0 < lam) (u : Vel × ℝ → Vel) :
    (⨆ p, ‖rescale α lam u p‖) = lam ^ (α - 1) * ⨆ p, ‖u p‖ := by
  have hpos : (0 : ℝ) ≤ lam ^ (α - 1) := (Real.rpow_pos_of_pos hlam _).le
  have hrange : (Set.range fun p : Vel × ℝ => ‖u (lam • p.1, lam ^ α * p.2)‖)
      = Set.range fun q : Vel × ℝ => ‖u q‖ := by
    simpa [Function.comp_def] using
      (scalingMap_surjective (α := α) hlam).range_comp (fun q : Vel × ℝ => ‖u q‖)
  have hcomp : (⨆ p : Vel × ℝ, ‖u (lam • p.1, lam ^ α * p.2)‖) = ⨆ q, ‖u q‖ :=
    congrArg sSup hrange
  calc (⨆ p, ‖rescale α lam u p‖)
      = ⨆ p : Vel × ℝ, lam ^ (α - 1) * ‖u (lam • p.1, lam ^ α * p.2)‖ := by
        simp_rw [norm_rescale hlam]
    _ = lam ^ (α - 1) * ⨆ p : Vel × ℝ, ‖u (lam • p.1, lam ^ α * p.2)‖ :=
        (Real.mul_iSup_of_nonneg hpos _).symm
    _ = lam ^ (α - 1) * ⨆ q, ‖u q‖ := by rw [hcomp]

/-- A bounded field rescales to a bounded field. -/
theorem bddAbove_range_norm_rescale {α lam : ℝ} (hlam : 0 < lam) (u : Vel × ℝ → Vel)
    (hu : BddAbove (Set.range fun p => ‖u p‖)) :
    BddAbove (Set.range fun p => ‖rescale α lam u p‖) := by
  obtain ⟨C, hC⟩ := hu
  refine ⟨lam ^ (α - 1) * C, ?_⟩
  rintro _ ⟨p, rfl⟩
  show ‖rescale α lam u p‖ ≤ lam ^ (α - 1) * C
  rw [norm_rescale hlam]
  exact mul_le_mul_of_nonneg_left (hC (Set.mem_range_self _))
    (Real.rpow_pos_of_pos hlam _).le

/-- **(e), with the boundedness hypothesis.**  For a bounded `u` the rescaled field is
bounded and its supremum norm is `λ^(α-1)` times that of `u`. -/
theorem iSup_norm_rescale_of_bddAbove {α lam : ℝ} (hlam : 0 < lam) (u : Vel × ℝ → Vel)
    (hu : BddAbove (Set.range fun p => ‖u p‖)) :
    BddAbove (Set.range fun p => ‖rescale α lam u p‖) ∧
      (⨆ p, ‖rescale α lam u p‖) = lam ^ (α - 1) * ⨆ p, ‖u p‖ :=
  ⟨bddAbove_range_norm_rescale hlam u hu, iSup_norm_rescale hlam u⟩

/-- **Corollary of (e).**  For a field with positive supremum norm, the sup norm is invariant
under the parabolic rescaling for *every* `λ > 0` exactly when `α = 1`.

This is the formal content of "`L^∞` is the critical Lebesgue space exactly at `α = 1`"
(§1.2).  For `α > 1` the exponent `α - 1` is strictly positive, so zooming in (`λ → 0`) sends
`⨆ ‖u_λ‖ → 0`: a bounded velocity carries no scale-invariant information at a concentrating
singularity. -/
theorem supNorm_scaling_invariant_iff {α : ℝ} (u : Vel × ℝ → Vel)
    (hpos : 0 < ⨆ p, ‖u p‖) :
    (∀ lam : ℝ, 0 < lam → (⨆ p, ‖rescale α lam u p‖) = ⨆ p, ‖u p‖) ↔ α = 1 := by
  constructor
  · intro h
    have h2 := h 2 (by norm_num)
    rw [iSup_norm_rescale (by norm_num : (0:ℝ) < 2) u] at h2
    have hone : (2 : ℝ) ^ (α - 1) = 1 :=
      mul_right_cancel₀ hpos.ne' (by rw [one_mul]; exact h2)
    have hb : (1 : ℝ) < 2 := by norm_num
    have hle : α - 1 ≤ 0 := by
      have h' : (2 : ℝ) ^ (α - 1) ≤ (2 : ℝ) ^ (0 : ℝ) := by rw [hone, Real.rpow_zero]
      exact (Real.rpow_le_rpow_left_iff hb).1 h'
    have hge : (0 : ℝ) ≤ α - 1 := by
      have h' : (2 : ℝ) ^ (0 : ℝ) ≤ (2 : ℝ) ^ (α - 1) := by rw [hone, Real.rpow_zero]
      exact (Real.rpow_le_rpow_left_iff hb).1 h'
    linarith
  · intro hα lam hlam
    rw [iSup_norm_rescale hlam u, hα]
    simp

end

end Criticality
end ForcedRoute

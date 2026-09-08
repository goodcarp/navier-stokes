import Mathlib

noncomputable section
open Real Set MeasureTheory Filter Topology
open scoped ContDiff

namespace EulerBlowup

abbrev R2 := EuclideanSpace ℝ (Fin 2)

def e (i : Fin 2) : R2 := EuclideanSpace.single i (1 : ℝ)

def pd (i : Fin 2) (g : R2 → ℝ) (x : R2) : ℝ := fderiv ℝ g x (e i)

def pdt (g : ℝ → R2 → ℝ) (t : ℝ) (x : R2) : ℝ := derivWithin (fun τ => g τ x) (Ici 0) t

def gradSize (g : R2 → ℝ) (x : R2) : ℝ := ‖fderiv ℝ g x‖

def perpGrad (ψ : R2 → ℝ) (x : R2) : R2 := -(pd 1 ψ x) • e 0 + (pd 0 ψ x) • e 1

def curl2 (u : R2 → R2) (x : R2) : ℝ :=
  pd 0 (fun y => (u y) 1) x - pd 1 (fun y => (u y) 0) x

def div2 (u : R2 → R2) (x : R2) : ℝ :=
  pd 0 (fun y => (u y) 0) x + pd 1 (fun y => (u y) 1) x

def advect (u : R2 → R2) (g : R2 → ℝ) (x : R2) : ℝ := (u x) 0 * pd 0 g x + (u x) 1 * pd 1 g x

def InClassP (θ : R2 → ℝ) (u : R2 → R2) : Prop :=
  (∀ x, θ (-x) = -θ x) ∧ (∀ x, u (-x) = -u x)

structure SmoothForce {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (T : ℝ) (F : ℝ × E → ℝ) : Prop where
  smooth : ContDiffOn ℝ ∞ F (Ico 0 T ×ˢ univ)
  support : ∃ Rball : ℝ, ∀ t ∈ Ico 0 T, ∀ x : E, Rball < ‖x‖ → F (t, x) = 0
  bounded : ∀ k : ℕ, ∃ M : ℝ, ∀ z ∈ Ico 0 T ×ˢ (univ : Set E),
    ‖iteratedFDerivWithin ℝ k F (Ico 0 T ×ˢ univ) z‖ ≤ M

structure ClassicalBoussinesq (κ T : ℝ) (θ p : ℝ → R2 → ℝ) (u fu : ℝ → R2 → R2)
    (fθ : ℝ → R2 → ℝ) : Prop where
  smoothθ : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => θ z.1 z.2) (Ico 0 T ×ˢ univ)
  smoothu : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => u z.1 z.2) (Ico 0 T ×ˢ univ)
  smoothp : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => p z.1 z.2) (Ico 0 T ×ˢ univ)
  divFree : ∀ t ∈ Ico 0 T, ∀ x, div2 (u t) x = 0
  finiteEnergy : ∀ t ∈ Ico 0 T, MemLp (u t) 2 volume ∧ MemLp (θ t) 2 volume
  eqθ : ∀ t ∈ Ico 0 T, ∀ x, pdt θ t x + advect (u t) (θ t) x = fθ t x
  equ : ∀ t ∈ Ico 0 T, ∀ x, ∀ i : Fin 2,
    pdt (fun τ y => (u τ y) i) t x + advect (u t) (fun y => (u t y) i) x + pd i (p t) x
      = κ * θ t x * (if i = 1 then 1 else 0) + (fu t x) i

structure InLipschitzClass (κ T' : ℝ) (θ : ℝ → R2 → ℝ) (u : ℝ → R2 → R2)
    (fθ : ℝ → R2 → ℝ) (fu : ℝ → R2 → R2) : Prop where
  lip : LocallyLipschitzOn (Icc 0 T' ×ˢ univ) (fun z : ℝ × R2 => (θ z.1 z.2, u z.1 z.2))
  energy : ∃ M, ∀ t ∈ Icc 0 T', MemLp (u t) 2 volume ∧ MemLp (θ t) 2 volume ∧
    (eLpNorm (u t) 2 volume).toReal ≤ M ∧ (eLpNorm (θ t) 2 volume).toReal ≤ M
  gradBounded : ∃ M, ∀ t ∈ Icc 0 T', ∀ x, ‖fderiv ℝ (θ t) x‖ ≤ M ∧ ‖fderiv ℝ (u t) x‖ ≤ M
  divFree : ∀ t ∈ Icc 0 T', ∀ᵐ x ∂(volume : Measure R2), div2 (u t) x = 0
  solvesθAE : ∀ᵐ z : ℝ × R2 ∂(volume.restrict (Icc 0 T' ×ˢ univ)),
    pdt θ z.1 z.2 + advect (u z.1) (θ z.1) z.2 = fθ z.1 z.2
  solvesuWeak : ∀ φ : ℝ × R2 → R2, ContDiff ℝ ∞ φ → HasCompactSupport φ →
    tsupport φ ⊆ Ioo 0 T' ×ˢ univ → (∀ z : ℝ × R2, div2 (fun y => φ (z.1, y)) z.2 = 0) →
    ∫ z : ℝ × R2, ∑ i : Fin 2,
      ( (u z.1 z.2) i * deriv (fun τ => (φ (τ, z.2)) i) z.1
        + ∑ j : Fin 2, (u z.1 z.2) i * (u z.1 z.2) j * pd j (fun y => (φ (z.1, y)) i) z.2
        + (κ * θ z.1 z.2 * (if i = 1 then 1 else 0) + (fu z.1 z.2) i) * (φ z) i ) = 0

structure Staircase (T : ℝ) where
  t : ℕ → ℝ
  A : ℕ → ℝ
  t_zero : t 0 = 0
  t_strictMono : StrictMono t
  t_lt : ∀ n, t n < T
  t_tendsto : Tendsto t atTop (𝓝 T)
  A_pos : ∀ n, 0 < A n
  superexp : ∀ q : ℝ, ∀ᶠ n in atTop, A (n + 1) ≥ q ^ n * A n

namespace Staircase
def stage {T : ℝ} (S : Staircase T) (τ : ℝ) : ℕ := sSup {n | S.t n ≤ τ}
end Staircase

def Theorem01_Statement : Prop :=
  ∀ κ : ℝ, 0 < κ →
  ∃ (T : ℝ) (_ : 0 < T) (S : Staircase T)
    (θ p fθ : ℝ → R2 → ℝ) (u fu : ℝ → R2 → R2),
    ClassicalBoussinesq κ T θ p u fu fθ ∧ (∀ t ∈ Ico 0 T, InClassP (θ t) (u t)) ∧
    (∀ t ∈ Ico 0 T, InClassP (fθ t) (fu t)) ∧
    (∀ T' < T, 0 < T' → InLipschitzClass κ T' θ u fθ fu) ∧
    HasCompactSupport (θ 0) ∧ HasCompactSupport (u 0) ∧
    SmoothForce T (fun z => fθ z.1 z.2) ∧ (∀ i : Fin 2, SmoothForce T (fun z => (fu z.1 z.2) i)) ∧
    (∀ τ ∈ Ico 0 T, ∃ x, gradSize (θ τ) x ≥ S.A (S.stage τ) / 2) ∧
    Tendsto (fun τ => S.A (S.stage τ)) (𝓝[<] T) atTop ∧
    (∃ M, ∀ τ ∈ Ico 0 T, ∀ x, |θ τ x| ≤ M) ∧
    (∀ M : ℝ, ∃ᶠ τ in 𝓝[<] T, ∃ x, |curl2 (u τ) x| ≥ M) ∧
    (∀ T' < T, 0 < T' → ∀ θ' u', InLipschitzClass κ T' θ' u' fθ fu →
        θ' 0 = θ 0 → u' 0 = u 0 → ∀ t ∈ Icc 0 T', θ' t = θ t ∧ u' t = u t)

end EulerBlowup

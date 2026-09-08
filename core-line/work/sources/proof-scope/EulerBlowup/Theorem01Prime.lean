import Mathlib
import EulerBlowup.Theorem01
import EulerBlowup.Lit.PotentialTheory
import EulerBlowup.Lit.FlatConstants
import EulerBlowup.Lit.BogovskiiAnnulus

noncomputable section
open Real Set Filter Topology

namespace EulerBlowup
open S2

structure FlatWindowsUndressedODE : Prop where
  crossing : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ R : DressedPreCrossingRun ψ₀ c₀ Θ₀ w flatControl Dressing.trivial,
      R.sv ∈ Icc (1.90026 : ℝ) 1.90748 ∧ coastSigma ψ₀ (R.x R.sv) ∈ Icc (1.40024 : ℝ) 1.40752 ∧
      (R.x R.sv).V / (R.x R.sv).Θ ∈ Icc (0.29240 : ℝ) 0.29652
  W6 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ R : DressedPreCrossingRun ψ₀ c₀ Θ₀ w flatControl Dressing.trivial,
      ∀ s ∈ Icc 0 R.sv, (R.x s).V / (R.x s).Θ ≤ 0.32
  signs : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ R : DressedPreCrossingRun ψ₀ c₀ Θ₀ w flatControl Dressing.trivial,
      let σv := coastSigma ψ₀ (R.x R.sv)
      (∀ D ∈ Icc (0.0005 : ℝ) 1.51, ∀ r, FlatLandingSlope ψ₀ σv R.sv Dressing.trivial (R.x R.sv)
          D r → 0 < r) ∧
      (∀ D ∈ Icc (1.56 : ℝ) 1.75, ∀ r, FlatLandingSlope ψ₀ σv R.sv Dressing.trivial (R.x R.sv)
          D r → r < 0) ∧
      (∀ D₁ ∈ Icc (1.40 : ℝ) 1.66, ∀ D₂ ∈ Icc (1.40 : ℝ) 1.66, D₁ < D₂ → ∀ r₁ r₂,
        FlatLandingSlope ψ₀ σv R.sv Dressing.trivial (R.x R.sv) D₁ r₁ →
        FlatLandingSlope ψ₀ σv R.sv Dressing.trivial (R.x R.sv) D₂ r₂ →
          r₂ - r₁ ≤ -0.3097 * (D₂ - D₁))
  W7 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding → M.D ∈ Icc (1.47 : ℝ) 1.60
  W4 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding → M.R.sv + M.D ≤ 3.8
  W1 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding →
      (M.y (M.R.sv + M.D)).Φ ∈ Icc (2.12 : ℝ) 2.42
  W3 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding →
      (M.y (M.R.sv + M.D)).C ∈ Icc (3.2 : ℝ) 5.7
  W3_hull : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding →
      (M.y (M.R.sv + M.D)).C ∈ Icc (3.94 : ℝ) 4.06
  W2 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding →
      arctan (ψ₀ / (M.y (M.R.sv + M.D)).C) / ψ₀ ∈ Icc (0.175 : ℝ) 0.311
  W9 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding →
      ∀ s ∈ Icc M.R.sv (M.R.sv + M.D), flatBeta ψ₀ M.σv M.D (s - M.R.sv) (M.y s).C - 1 ≤ 0.35
  W8 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding →
      ∀ rT : ℝ → ℝ, (∀ᶠ D in 𝓝 M.D, FlatLandingSlope ψ₀ M.σv M.R.sv Dressing.trivial (M.R.x M.R.sv) D (rT D)) →
        ∀ d' : ℝ, HasDerivAt rT d' M.D → d' ≤ -0.29
  W11 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding →
      ∀ J : ℝ, IsFlatT3Sensitivity ψ₀ M.σv M.D M.R.sv M.R.x M.y J → |J| ≤ 2.5
  W5 : ∀ ψ₀ c₀ Θ₀ w, IsBoxPoint ψ₀ c₀ Θ₀ w →
    ∀ M : FlatManoeuvre ψ₀ c₀ Θ₀ w Dressing.trivial, M.IsLanding →
      ∀ s ∈ Icc 0 (M.R.sv + M.D), |iteratedDerivWithin 2 M.control (Ici 0) s| ≤ 15.5

structure LandingWindowsODE : Prop where
  undressed : FlatWindowsUndressedODE
  dressedBig : FlatWindowsDressedBig
  realizedBig : FlatCostRealizedBig
  dressedSmall : FlatWindowsDressedSmall

theorem flatWindowsUndressed_of_ODE (h : FlatWindowsUndressedODE) : FlatWindowsUndressed := by
  exact ⟨Lit.flatZ_mem, Lit.flatStep_rampCost, h.crossing, h.W6, h.signs, h.W7, h.W4, h.W1, h.W3, h.W3_hull, h.W2, h.W9, h.W8, h.W11, h.W5⟩

theorem landingWindowsW_of_dressed (hB : FlatWindowsDressedBig) (hR : FlatCostRealizedBig)
    (hS : FlatWindowsDressedSmall) : LandingWindowsW := by
  exact ⟨⟨Lit.flatZ_mem⟩, hB, hR, hS⟩

theorem landingWindowsW_of_ODE (h : LandingWindowsODE) : LandingWindowsW := by
  exact landingWindowsW_of_dressed h.dressedBig h.realizedBig h.dressedSmall

theorem theorem01' (hW : LandingWindowsW) (hU : FlatWindowsUndressed) : Theorem01_Statement := by
  exact theorem01 Lit.bogovskii_annulus_R2 Lit.potential_theory_R2 hW hU

theorem theorem01'' (hW : LandingWindowsODE) : Theorem01_Statement := by
  exact theorem01' (landingWindowsW_of_ODE hW) (flatWindowsUndressed_of_ODE hW.undressed)

end EulerBlowup

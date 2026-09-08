import EulerBlowup.Num.Cert.FinalPrime

noncomputable section
open Real Set MeasureTheory Filter Topology
open scoped ContDiff

open EulerBlowup in
theorem euler_smooth_force_blowup :
    ∃ (T : ℝ) (_ : 0 < T) (u f : ℝ → R3 → R3) (p : ℝ → R3 → ℝ),
      ClassicalEuler T u f p ∧
      HasCompactSupport (u 0) ∧
      (∀ i : Fin 3, SmoothForce T (fun z => (f z.1 z.2) i)) ∧
      (∃ g : ℝ → ℝ, Tendsto g (𝓝[<] T) atTop ∧ ∀ τ ∈ Ico 0 T, ∃ x, ‖curl3 (u τ) x‖ ≥ g τ) ∧
      (∀ M : ℝ, ∃ T' ∈ Ico 0 T, ∃ g : ℝ → ℝ, (∀ τ ∈ Icc 0 T', 0 ≤ g τ ∧ ∃ x, ‖curl3 (u τ) x‖ ≥ g τ) ∧
          IntervalIntegrable g volume 0 T' ∧ ∫ τ in (0 : ℝ)..T', g τ ≥ M) ∧
      (∀ T' < T, 0 < T' → InLipschitzClass3 T' u f) ∧
      (∀ T' < T, 0 < T' → ∀ u', InLipschitzClass3 T' u' f → u' 0 = u 0 →
          ∀ t ∈ Icc 0 T', u' t = u t) := by
  obtain ⟨T, hT, S, u, f, p, h1, _h2, h3, _h4, _h5, _h6, _h7, h8,
      _h9, _h10, _h11, _h12, h13, h14, h15, _h16, h17, h18⟩ :=
    EulerBlowup.Num.Cert.theorem01prime 1 one_pos
  refine ⟨T, hT, u, f, p, h1, h3, h8, ?_, h15, h17, h18⟩
  refine ⟨fun τ => S.A (S.stage τ) / (2 * 1), h14.atTop_div_const (by norm_num), ?_⟩
  intro τ hτ
  obtain ⟨x, _, hx⟩ := h13 τ hτ
  exact ⟨x, hx⟩

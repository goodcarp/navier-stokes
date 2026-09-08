import Mathlib
import EulerBlowup.LandingWindows
import EulerBlowup.S7.Readout
import EulerBlowup.S7.Uniqueness
import EulerBlowup.S6.Landing
import EulerBlowup.D.PDE
import EulerBlowup.S4.ProfilesExist
import EulerBlowup.S7.TimeRegularity
import EulerBlowup.S6.BaseRun

noncomputable section
open Real Set MeasureTheory Filter Topology
open scoped ContDiff

namespace EulerBlowup

open S4 S6 S7

def Theorem01_Body (κ : ℝ) : Prop :=
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

theorem theorem01_body_iff : Theorem01_Statement ↔ ∀ κ : ℝ, 0 < κ → Theorem01_Body κ := by
  exact Iff.rfl

theorem boussinesq_rescale_aux_pd (c : ℝ) (g : R2 → ℝ) (i : Fin 2) (x : R2) :
    pd i (fun y => c * g y) x = c * pd i g x := by
  unfold pd
  rw [show (fun y => c * g y) = c • g from rfl, fderiv_const_smul_field]
  rfl

theorem boussinesq_rescale_aux_pdt (c : ℝ) (θ : ℝ → R2 → ℝ) (t : ℝ) (x : R2) :
    pdt (fun τ y => c * θ τ y) t x = c * pdt θ t x := by
  unfold pdt
  exact derivWithin_const_mul_field c

theorem boussinesq_rescale_aux_advect (c : ℝ) (u : R2 → R2) (g : R2 → ℝ) (x : R2) :
    advect u (fun y => c * g y) x = c * advect u g x := by
  unfold advect; rw [boussinesq_rescale_aux_pd, boussinesq_rescale_aux_pd]; ring

theorem boussinesq_rescale_aux_gradSize (c : ℝ) (g : R2 → ℝ) (x : R2) :
    gradSize (fun y => c * g y) x = |c| * gradSize g x := by
  unfold gradSize
  rw [show (fun y => c * g y) = c • g from rfl, fderiv_const_smul_field, Pi.smul_apply, norm_smul, Real.norm_eq_abs]

theorem boussinesq_rescale_aux_lip {κ κ' c T' : ℝ} (hκc : κ * c = κ') {θ : ℝ → R2 → ℝ} {u : ℝ → R2 → R2}
    {fθ : ℝ → R2 → ℝ} {fu : ℝ → R2 → R2} (h : InLipschitzClass κ' T' θ u fθ fu) :
    InLipschitzClass κ T' (fun t x => c * θ t x) u (fun t x => c * fθ t x) fu := by
  obtain ⟨hlip, ⟨M, hM⟩, ⟨G, hG⟩, hdiv, hθae, hweak⟩ := h
  refine ⟨?_, ⟨max M (|c| * M), fun t ht => ?_⟩, ⟨max G (|c| * G), fun t ht x => ?_⟩, hdiv, ?_, ?_⟩
  ·
    set L : ℝ × R2 →L[ℝ] ℝ × R2 :=
      ContinuousLinearMap.prodMap (c • ContinuousLinearMap.id ℝ ℝ) (ContinuousLinearMap.id ℝ R2) with hL
    have hfun : (fun z : ℝ × R2 => (c * θ z.1 z.2, u z.1 z.2)) = L ∘ (fun z : ℝ × R2 => (θ z.1 z.2, u z.1 z.2)) := by
      funext z; simp [hL, smul_eq_mul]
    rw [hfun]
    intro x hx
    obtain ⟨K, t, ht, hKt⟩ := hlip hx
    exact ⟨_, t, ht, L.lipschitz.comp_lipschitzOnWith hKt⟩
  · obtain ⟨h1, h2, h3, h4⟩ := hM t ht
    refine ⟨h1, h2.const_mul c, h3.trans (le_max_left _ _), ?_⟩
    have : (fun x => c * θ t x) = c • θ t := rfl
    rw [this, eLpNorm_const_smul, ENNReal.toReal_mul, Real.enorm_eq_ofReal_abs, ENNReal.toReal_ofReal (abs_nonneg c)]
    exact (mul_le_mul_of_nonneg_left h4 (abs_nonneg c)).trans (le_max_right _ _)
  · obtain ⟨h1, h2⟩ := hG t ht x
    refine ⟨?_, h2.trans (le_max_left _ _)⟩
    have : (fun y => c * θ t y) = c • θ t := rfl
    rw [this, fderiv_const_smul_field, Pi.smul_apply, norm_smul, Real.norm_eq_abs]
    exact (mul_le_mul_of_nonneg_left h1 (abs_nonneg c)).trans (le_max_right _ _)
  · filter_upwards [hθae] with z hz
    rw [boussinesq_rescale_aux_pdt, boussinesq_rescale_aux_advect, ← mul_add, hz]
  · intro φ h1 h2 h3 h4
    have := hweak φ h1 h2 h3 h4
    have hκ : ∀ a : ℝ, κ * (c * a) = κ' * a := fun a => by rw [← mul_assoc, hκc]
    simp_rw [hκ]
    exact this

theorem boussinesq_rescale_aux_force (c T : ℝ) {F : ℝ × R2 → ℝ} (h : SmoothForce T F) :
    SmoothForce T (fun z => c * F z) := by
  obtain ⟨hs, ⟨Rb, hRb⟩, hb⟩ := h
  refine ⟨contDiffOn_const.mul hs, ⟨Rb, fun t ht x hx => by rw [hRb t ht x hx, mul_zero]⟩, fun k => ?_⟩
  obtain ⟨M, hM⟩ := hb k
  refine ⟨|c| * M, fun z hz => ?_⟩
  have hu : UniqueDiffOn ℝ (Ico (0:ℝ) T ×ˢ (univ : Set R2)) := (uniqueDiffOn_Ico 0 T).prod uniqueDiffOn_univ
  have hf : ContDiffWithinAt ℝ k F (Ico 0 T ×ˢ univ) z := (hs z hz).of_le (by exact_mod_cast le_top)
  have : (fun z => c * F z) = c • F := rfl
  rw [this, iteratedFDerivWithin_const_smul_apply hf hu hz, norm_smul, Real.norm_eq_abs]
  exact mul_le_mul_of_nonneg_left (hM z hz) (abs_nonneg c)

theorem boussinesq_rescale_aux_classical {κ κ' c T : ℝ} (hκc : κ * c = κ') {θ p : ℝ → R2 → ℝ} {u fu : ℝ → R2 → R2}
    {fθ : ℝ → R2 → ℝ} (h : ClassicalBoussinesq κ' T θ p u fu fθ) :
    ClassicalBoussinesq κ T (fun t x => c * θ t x) p u fu (fun t x => c * fθ t x) := by
  obtain ⟨hsθ, hsu, hsp, hdiv, hen, heqθ, hequ⟩ := h
  refine ⟨contDiffOn_const.mul hsθ, hsu, hsp, hdiv, fun t ht => ⟨(hen t ht).1, (hen t ht).2.const_mul c⟩,
    fun t ht x => ?_, fun t ht x i => ?_⟩
  · rw [boussinesq_rescale_aux_pdt, boussinesq_rescale_aux_advect, ← mul_add, heqθ t ht x]
  · rw [hequ t ht x i, ← hκc]; ring

theorem boussinesq_rescale_aux_staircase {T : ℝ} (c : ℝ) (hc : 0 < c) (S : Staircase T) :
    ∃ S' : Staircase T, (∀ τ, S'.stage τ = S.stage τ) ∧ ∀ n, S'.A n = c * S.A n := by
  refine ⟨⟨S.t, (fun n => c * S.A n), S.t_zero, S.t_strictMono, S.t_lt, S.t_tendsto, (fun n => mul_pos hc (S.A_pos n)),
    (fun q => ?_)⟩, fun τ => rfl, fun n => rfl⟩
  filter_upwards [S.superexp q] with n hn
  have := mul_le_mul_of_nonneg_left hn hc.le
  nlinarith

theorem boussinesq_rescale (κ κ' : ℝ) (hκ : 0 < κ) (hκ' : 0 < κ') (h : Theorem01_Body κ') : Theorem01_Body κ := by
  obtain ⟨T, hT, S, θ', p, fθ', u, fu, hcl, hP, hPf, hL, hcθ, hcu, hfθ, hfu, hread, htend, hbd, hvort, huniq⟩ := h
  set c : ℝ := κ' / κ with hc
  have hcpos : 0 < c := div_pos hκ' hκ
  have hκc : κ * c = κ' := by rw [hc]; field_simp
  have hc' : κ' * c⁻¹ = κ := by rw [hc, inv_div]; field_simp
  obtain ⟨S', hstage, hA⟩ := boussinesq_rescale_aux_staircase c hcpos S
  refine ⟨T, hT, S', fun t x => c * θ' t x, p, fun t x => c * fθ' t x, u, fu,
    boussinesq_rescale_aux_classical hκc hcl, ?_, ?_, fun T' hT' hT'0 => boussinesq_rescale_aux_lip hκc (hL T' hT' hT'0),
    ?_, hcu, ?_, hfu, ?_, ?_, ?_, hvort, ?_⟩
  · intro t ht; obtain ⟨h1, h2⟩ := hP t ht
    exact ⟨fun x => by show c * θ' t (-x) = -(c * θ' t x); rw [h1]; ring, h2⟩
  · intro t ht; obtain ⟨h1, h2⟩ := hPf t ht
    exact ⟨fun x => by show c * fθ' t (-x) = -(c * fθ' t x); rw [h1]; ring, h2⟩
  · show HasCompactSupport (fun x => c * θ' 0 x)
    exact hcθ.mul_left
  · exact boussinesq_rescale_aux_force c T hfθ
  · intro τ hτ
    obtain ⟨x, hx⟩ := hread τ hτ
    refine ⟨x, ?_⟩
    rw [hA, hstage, ge_iff_le]
    have h1 : gradSize (fun y => c * θ' τ y) x = |c| * gradSize (θ' τ) x := boussinesq_rescale_aux_gradSize c (θ' τ) x
    show _ ≤ gradSize (fun y => c * θ' τ y) x
    rw [h1, abs_of_pos hcpos]
    have := hx.le
    nlinarith
  · have : (fun τ => S'.A (S'.stage τ)) = fun τ => c * S.A (S.stage τ) := by funext τ; rw [hA, hstage]
    rw [this]; exact htend.const_mul_atTop hcpos
  · obtain ⟨M, hM⟩ := hbd
    exact ⟨|c| * M, fun τ hτ x => by show |c * θ' τ x| ≤ _; rw [abs_mul]; exact mul_le_mul_of_nonneg_left (hM τ hτ x) (abs_nonneg c)⟩
  · intro T' hT' hT'0 θ'' u'' hL'' h0θ h0u t ht
    have hL2 : InLipschitzClass κ' T' (fun t x => c⁻¹ * θ'' t x) u'' (fun t x => c⁻¹ * (c * fθ' t x)) fu :=
      boussinesq_rescale_aux_lip hc' hL''
    have hff : (fun t x => c⁻¹ * (c * fθ' t x)) = fθ' := by
      funext t x; rw [← mul_assoc, inv_mul_cancel₀ hcpos.ne', one_mul]
    rw [hff] at hL2
    have h0' : (fun x => c⁻¹ * θ'' 0 x) = θ' 0 := by
      funext x; rw [h0θ]; show c⁻¹ * (c * θ' 0 x) = θ' 0 x; rw [← mul_assoc, inv_mul_cancel₀ hcpos.ne', one_mul]
    obtain ⟨h1, h2⟩ := huniq T' hT' hT'0 (fun t x => c⁻¹ * θ'' t x) u'' hL2 h0' h0u t ht
    refine ⟨?_, h2⟩
    funext x
    have := congrFun h1 x
    show θ'' t x = c * θ' t x
    have h3 : c⁻¹ * θ'' t x = θ' t x := this
    rw [← h3, ← mul_assoc, mul_inv_cancel₀ hcpos.ne', one_mul]

theorem designed_classical_aux_pd_comm (f : R2 → ℝ) (hf : ContDiff ℝ ∞ f) (x : R2) (i j : Fin 2) :
    pd i (pd j f) x = pd j (pd i f) x := by
  have h2 : ContDiffAt ℝ 2 f x := (hf.of_le (by norm_cast)).contDiffAt
  have hsymm := h2.isSymmSndFDerivAt (by simp)
  have hdf : Differentiable ℝ (fderiv ℝ f) :=
    ((hf.fderiv_right (m := ∞) (by norm_cast)).differentiable (by simp))
  have key : ∀ a b : Fin 2, pd a (pd b f) x = fderiv ℝ (fderiv ℝ f) x (e a) (e b) := by
    intro a b
    unfold pd
    rw [fderiv_clm_apply (hdf x) (differentiableAt_const _)]
    simp
  rw [key, key]
  exact hsymm (e i) (e j)

theorem designed_classical_aux_pd_mul (f g : R2 → ℝ) (x : R2) (hf : DifferentiableAt ℝ f x)
    (hg : DifferentiableAt ℝ g x) (i : Fin 2) :
    pd i (fun y => f y * g y) x = pd i f x * g x + f x * pd i g x := by
  unfold pd
  rw [fderiv_fun_mul hf hg]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.coe_smul', Pi.smul_apply, smul_eq_mul]
  ring

theorem designed_classical_aux_curl_advect (u : R2 → R2) (hu : ContDiff ℝ ∞ u) (hdiv : ∀ x, div2 u x = 0) (x : R2) :
    curl2 (fun y => WithLp.toLp 2 (fun i : Fin 2 => advect u (fun z => (u z) i) y)) x = advect u (curl2 u) x := by
  set u0 : R2 → ℝ := fun y => (u y) 0 with hu0
  set u1 : R2 → ℝ := fun y => (u y) 1 with hu1
  have hc0 : ContDiff ℝ ∞ u0 := (EuclideanSpace.proj (𝕜 := ℝ) (ι := Fin 2) 0).contDiff.comp hu
  have hc1 : ContDiff ℝ ∞ u1 := (EuclideanSpace.proj (𝕜 := ℝ) (ι := Fin 2) 1).contDiff.comp hu
  have hd0 : Differentiable ℝ u0 := hc0.differentiable (by simp)
  have hd1 : Differentiable ℝ u1 := hc1.differentiable (by simp)
  have hpdc : ∀ (f : R2 → ℝ), ContDiff ℝ ∞ f → ∀ i, ContDiff ℝ ∞ (pd i f) := fun f hf i =>
    (hf.fderiv_right (m := ∞) (by norm_cast)).clm_apply contDiff_const
  have hpdd : ∀ (f : R2 → ℝ), ContDiff ℝ ∞ f → ∀ i, Differentiable ℝ (pd i f) := fun f hf i =>
    (hpdc f hf i).differentiable (by simp)
  have hcomp : ∀ j : Fin 2, (fun y => (WithLp.toLp 2 (fun i : Fin 2 => advect u (fun z => (u z) i) y) : R2) j) =
      fun y => u0 y * pd 0 (fun z => (u z) j) y + u1 y * pd 1 (fun z => (u z) j) y := by
    intro j; funext y; simp [advect, hu0, hu1]
  unfold curl2
  rw [hcomp 1, hcomp 0]
  show pd 0 (fun y => u0 y * pd 0 u1 y + u1 y * pd 1 u1 y) x - pd 1 (fun y => u0 y * pd 0 u0 y + u1 y * pd 1 u0 y) x =
    advect u (curl2 u) x
  have hadd : ∀ (f g : R2 → ℝ) (i : Fin 2), Differentiable ℝ f → Differentiable ℝ g →
      pd i (fun y => f y + g y) x = pd i f x + pd i g x := by
    intro f g i hf hg; unfold pd; rw [fderiv_fun_add (hf x) (hg x)]; rfl
  rw [hadd (fun y => u0 y * pd 0 u1 y) (fun y => u1 y * pd 1 u1 y) 0 (hd0.mul (hpdd u1 hc1 0)) (hd1.mul (hpdd u1 hc1 1)),
    hadd (fun y => u0 y * pd 0 u0 y) (fun y => u1 y * pd 1 u0 y) 1 (hd0.mul (hpdd u0 hc0 0)) (hd1.mul (hpdd u0 hc0 1))]
  rw [designed_classical_aux_pd_mul _ _ x (hd0 x) (hpdd u1 hc1 0 x),
    designed_classical_aux_pd_mul _ _ x (hd1 x) (hpdd u1 hc1 1 x),
    designed_classical_aux_pd_mul _ _ x (hd0 x) (hpdd u0 hc0 0 x),
    designed_classical_aux_pd_mul _ _ x (hd1 x) (hpdd u0 hc0 1 x)]
  have hcurl : curl2 u = fun y => pd 0 u1 y - pd 1 u0 y := by funext y; rfl
  have hadv : advect u (curl2 u) x = u0 x * pd 0 (curl2 u) x + u1 x * pd 1 (curl2 u) x := rfl
  have hsub : ∀ i : Fin 2, pd i (curl2 u) x = pd i (pd 0 u1) x - pd i (pd 1 u0) x := by
    intro i; rw [hcurl]
    show fderiv ℝ (fun y => pd 0 u1 y - pd 1 u0 y) x (e i) = fderiv ℝ (pd 0 u1) x (e i) - fderiv ℝ (pd 1 u0) x (e i)
    rw [fderiv_fun_sub ((hpdd u1 hc1 0) x) ((hpdd u0 hc0 1) x)]; rfl
  rw [hadv, hsub 0, hsub 1]
  have hdivx : pd 0 u0 x + pd 1 u1 x = 0 := hdiv x
  have hs1 : pd 0 (pd 1 u1) x = pd 1 (pd 0 u1) x := designed_classical_aux_pd_comm u1 hc1 x 0 1
  have hs0 : pd 1 (pd 0 u0) x = pd 0 (pd 1 u0) x := designed_classical_aux_pd_comm u0 hc0 x 1 0
  have e1 : pd 0 u1 x * pd 1 u1 x = -(pd 0 u0 x * pd 0 u1 x) := by
    have : pd 1 u1 x = -pd 0 u0 x := by linarith
    rw [this]; ring
  have e2 : pd 1 u0 x * pd 0 u0 x = -(pd 1 u1 x * pd 1 u0 x) := by
    have : pd 0 u0 x = -pd 1 u1 x := by linarith
    rw [this]; ring
  rw [hs1, hs0]
  linarith [e1, e2]

theorem designed_classical_aux_red {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ) (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T)
    (k : ℕ) (hk : 1 ≤ k) : T.StateRed k :=
  (hT.state k hk).1

theorem designed_classical_aux_mono {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ) (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) :
    Monotone T.t := by
  intro a b hab
  rcases eq_or_lt_of_le hab with h | h
  · rw [h]
  · exact ((designed_classical_aux_red NStar Cg T hT (b + 1) (by omega)).t_mono
      (show a ≤ b + 1 by omega) (show b ≤ b + 1 by omega) h).le

theorem designed_classical_aux_div_perpGrad (c : ℝ) (f : R2 → ℝ) (hf : ContDiff ℝ ∞ f) (x : R2) :
    div2 (fun y => c • perpGrad f y) x = 0 := by
  have hpdc : ∀ i : Fin 2, ContDiff ℝ ∞ (pd i f) := fun i =>
    (hf.fderiv_right (m := ∞) (by norm_cast)).clm_apply contDiff_const
  have hpd : ∀ i : Fin 2, Differentiable ℝ (pd i f) := fun i => (hpdc i).differentiable (by simp)
  have h0 : (fun y => (c • perpGrad f y) 0) = fun y => (-c) * pd 1 f y := by funext y; simp [perpGrad, e]
  have h1 : (fun y => (c • perpGrad f y) 1) = fun y => c * pd 0 f y := by funext y; simp [perpGrad, e]
  unfold div2
  rw [h0, h1]
  have e0 : pd 0 (fun y => (-c) * pd 1 f y) x = (-c) * pd 0 (pd 1 f) x := by
    rw [designed_classical_aux_pd_mul _ _ x (differentiableAt_const _) (hpd 1 x)]
    simp [pd]
  have e1 : pd 1 (fun y => c * pd 0 f y) x = c * pd 1 (pd 0 f) x := by
    rw [designed_classical_aux_pd_mul _ _ x (differentiableAt_const _) (hpd 0 x)]
    simp [pd]
  rw [e0, e1, designed_classical_aux_pd_comm f hf x 0 1]; ring

theorem designed_classical_aux_layer {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ)
    (hNuni : ∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) (hP : PotentialTheoryR2)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) (k : ℕ) (τ : ℝ) (hτ : τ ∈ Ico (T.t k) (T.t (k + 1)))
    (m : ℕ) (hm1 : 1 ≤ m) (hmk : m ≤ k + 1) :
    ContDiff ℝ ∞ (T.uL m τ) ∧ ContDiff ℝ ∞ (T.θL m τ) ∧ HasCompactSupport (T.θL m τ) ∧
    (∀ x, curl2 (T.uL m τ) x = T.omL m τ x) ∧ (∀ x, div2 (T.uL m τ) x = 0) ∧ MemLp (T.uL m τ) 2 volume := by
  have hmono := designed_classical_aux_mono NStar Cg T hT
  have hm0 : m ≠ 0 := by omega
  have hτ2 : τ ≤ T.t (k + 2) := hτ.2.le.trans (hmono (by omega))
  have ht' : T.t (k + 2) ∈ T.I (k + 2) := ⟨le_rfl, hmono (by omega)⟩
  obtain ⟨-, -, ⟨R, hKin, hAmp, K, hS4, hG, -⟩, -⟩ := hT.older (k + 2) (by omega) (T.t (k + 2)) ht' m hm1 (by omega)
  obtain ⟨tv, hSO⟩ := hT.sizes (k + 2) (by omega) (T.t (k + 2)) ht' m hm1 (by omega) R hKin
  obtain ⟨hDsel, hrun⟩ := hT.run (k + 1) (by omega)
  have htin : T.tIn m ≤ τ := (hmono (show m - 1 ≤ k by omega)).trans hτ.1
  have hera : τ ∈ Icc (T.bg m (T.t (k + 2))).tIn (T.bg m (T.t (k + 2))).τ := ⟨htin, hτ2⟩
  have hreal := hrun.realised m hm1 (by omega) R hKin τ ⟨htin, hτ2⟩
  set ζ := S4.omegaLayer Pr (T.P m) R (T.N (m - 1)) τ with hζdef
  have hueq : T.uL m τ = biotSavart ζ := funext fun y => (hreal y).2.2
  have hωeq : T.omL m τ = ζ := funext fun y => (hreal y).2.1
  have hθeq : T.θL m τ = S4.thetaLayer Pr (T.P m) R (T.N (m - 1)) τ := funext fun y => (hreal y).1
  have hN : T.N m = (T.P m).N := by unfold TowerData.N; rw [if_neg hm0]
  have hthr : S4.NthrUniform Pr (T.P m).β (Cg (T.P m).β) ≤ (T.P m).N := by
    rw [← hN]; exact (hNuni _ (T.P m).β_pos).trans (hT.thresholds m hm1)
  have hcutall := S4.nthrUniform_spec Pr (T.P m).β (Cg (T.P m).β) (T.P m) rfl hthr R K (T.N (m - 1)) hG hS4 hAmp
  have hslice : ∀ (F : ℝ → R2 → ℝ), ContDiffOn ℝ ∞ (fun z : ℝ × R2 => if (T.bg m (T.t (k + 2))).tIn ≤ z.1 then
      F z.1 z.2 else 0) (Iic (T.bg m (T.t (k + 2))).τ ×ˢ univ) → ContDiff ℝ ∞ (F τ) := by
    intro F hcut
    have hcomp : ContDiffOn ℝ ∞ ((fun z : ℝ × R2 => if (T.bg m (T.t (k + 2))).tIn ≤ z.1 then F z.1 z.2 else 0) ∘
        fun y : R2 => (τ, y)) univ :=
      hcut.comp (contDiffOn_const.prodMk contDiffOn_id) (fun y _ => ⟨show τ ≤ T.t (k + 2) from hτ2, mem_univ _⟩)
    have h2 : ContDiffOn ℝ ∞ (F τ) univ := hcomp.congr fun y _ => by
      simp only [Function.comp_apply, if_pos (show (T.bg m (T.t (k + 2))).tIn ≤ τ from htin)]
    exact contDiffOn_univ.1 h2
  have hsmω : ContDiff ℝ ∞ ζ := hslice _ hcutall.2.2.1
  have hsmθ : ContDiff ℝ ∞ (T.θL m τ) := by rw [hθeq]; exact hslice _ hcutall.2.1
  obtain ⟨Rχ, hRχ⟩ : ∃ r, tsupport Pr.χc ⊆ Metric.closedBall (0 : R2) r :=
    (Pr.hχc.compact : IsCompact (tsupport Pr.χc)).isBounded.subset_closedBall 0
  set ρ := max (max ((T.P m).gammaStar * (T.P m).ell) Rχ) 0 + 1 with hρ
  have hsuppω : tsupport ζ ⊆ Metric.ball (0 : R2) ρ := by
    refine (hSO.supp_om τ hera).trans ?_
    intro x hx
    rw [Metric.mem_ball, dist_zero_right]
    rcases hx with h | h
    · have := Metric.mem_closedBall.1 h; rw [dist_zero_right] at this; rw [hρ]
      linarith [le_max_left ((T.P m).gammaStar * (T.P m).ell) Rχ, le_max_left (max ((T.P m).gammaStar * (T.P m).ell) Rχ) 0]
    · have := Metric.mem_closedBall.1 (hRχ h); rw [dist_zero_right] at this; rw [hρ]
      linarith [le_max_right ((T.P m).gammaStar * (T.P m).ell) Rχ, le_max_left (max ((T.P m).gammaStar * (T.P m).ell) Rχ) 0]
  have hcsω : HasCompactSupport ζ :=
    HasCompactSupport.of_support_subset_isCompact (isCompact_closedBall (0 : R2) ρ)
      ((subset_tsupport _).trans (hsuppω.trans Metric.ball_subset_closedBall))
  have hcsθ : HasCompactSupport (T.θL m τ) := by
    rw [hθeq]
    exact HasCompactSupport.of_support_subset_isCompact (isCompact_closedBall (0 : R2) _)
      ((subset_tsupport _).trans (hSO.supp_th τ hera))
  obtain ⟨hNPsm, -, -⟩ := hP.newtonian ζ hsmω hcsω
  have husm : ContDiff ℝ ∞ (biotSavart ζ) := by
    have hpd : ∀ i : Fin 2, ContDiff ℝ ∞ (pd i (newtonPotential ζ)) := fun i =>
      (hNPsm.fderiv_right (m := ∞) (by norm_cast)).clm_apply contDiff_const
    show ContDiff ℝ ∞ (perpGrad (newtonPotential ζ))
    unfold perpGrad
    exact ((hpd 1).neg.smul contDiff_const).add ((hpd 0).smul contDiff_const)
  have hmean : ∫ y, ζ y = 0 := hSO.meanZero τ hera
  have hρpos : 0 < ρ := by
    rw [hρ]; linarith [le_max_right (max ((T.P m).gammaStar * (T.P m).ell) Rχ) 0]
  have hL2 : MemLp (biotSavart ζ) 2 volume := by
    obtain ⟨C, hC⟩ := biotSavart_memLp hP ρ hρpos
    obtain ⟨M, hM⟩ : ∃ M, ∀ y, |ζ y| ≤ M := by
      obtain ⟨M, hM⟩ := hsmω.continuous.norm.bddAbove_range_of_hasCompactSupport hcsω.norm
      exact ⟨M, fun y => by simpa [Real.norm_eq_abs] using hM (Set.mem_range_self y)⟩
    exact (hC ζ hsmω hsuppω hmean M hM).1
  refine ⟨by rw [hueq]; exact husm, hsmθ, hcsθ, fun x => ?_, fun x => ?_, by rw [hueq]; exact hL2⟩
  · rw [hueq, hωeq]; exact biotSavart_curl hP ζ hsmω hcsω x
  · rw [hueq]; exact biotSavart_div hP ζ hsmω hcsω x
theorem designed_classical_aux_curl_finsum (G : ℕ → R2 → R2) (s : Finset ℕ) (x : R2)
    (hG : ∀ m ∈ s, Differentiable ℝ (G m)) :
    curl2 (fun y => ∑ m ∈ s, G m y) x = ∑ m ∈ s, curl2 (G m) x ∧
    div2 (fun y => ∑ m ∈ s, G m y) x = ∑ m ∈ s, div2 (G m) x := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [curl2, div2, pd]
  | insert a s ha ih =>
    have hGa : Differentiable ℝ (G a) := hG a (Finset.mem_insert_self a s)
    have hGs : ∀ m ∈ s, Differentiable ℝ (G m) := fun m hm => hG m (Finset.mem_insert_of_mem hm)
    obtain ⟨ih1, ih2⟩ := ih hGs
    have hS : Differentiable ℝ (fun y => ∑ m ∈ s, G m y) := by
      have : (fun y => ∑ m ∈ s, G m y) = ∑ m ∈ s, G m := by funext y; simp
      rw [this]; exact Differentiable.sum fun m hm => hGs m hm
    simp only [Finset.sum_insert ha]
    have hc : ∀ j : Fin 2, ∀ u : R2 → R2, Differentiable ℝ u → Differentiable ℝ (fun y => (u y) j) :=
      fun j u hu => (EuclideanSpace.proj (𝕜 := ℝ) (ι := Fin 2) j).differentiable.comp hu
    have hsplit : ∀ j : Fin 2, (fun y => (G a y + ∑ m ∈ s, G m y) j) = fun y => (G a y) j + (∑ m ∈ s, G m y) j := by
      intro j; funext y; simp
    constructor
    · rw [← ih1]
      unfold curl2 pd
      rw [hsplit 1, hsplit 0, fderiv_fun_add ((hc 1 _ hGa) x) ((hc 1 _ hS) x), fderiv_fun_add ((hc 0 _ hGa) x) ((hc 0 _ hS) x)]
      simp only [ContinuousLinearMap.add_apply]; ring
    · rw [← ih2]
      unfold div2 pd
      rw [hsplit 0, hsplit 1, fderiv_fun_add ((hc 0 _ hGa) x) ((hc 0 _ hS) x), fderiv_fun_add ((hc 1 _ hGa) x) ((hc 1 _ hS) x)]
      simp only [ContinuousLinearMap.add_apply]; ring

theorem designed_classical_aux_curl_add2 (f g : R2 → R2) (x : R2) (hf : Differentiable ℝ f) (hg : Differentiable ℝ g) :
    curl2 (fun y => f y + g y) x = curl2 f x + curl2 g x ∧ div2 (fun y => f y + g y) x = div2 f x + div2 g x := by
  have hc : ∀ j : Fin 2, ∀ u : R2 → R2, Differentiable ℝ u → Differentiable ℝ (fun y => (u y) j) :=
    fun j u hu => (EuclideanSpace.proj (𝕜 := ℝ) (ι := Fin 2) j).differentiable.comp hu
  have hsplit : ∀ j : Fin 2, (fun y => (f y + g y) j) = fun y => (f y) j + (g y) j := by intro j; funext y; simp
  constructor
  · unfold curl2 pd
    rw [hsplit 1, hsplit 0, fderiv_fun_add ((hc 1 _ hf) x) ((hc 1 _ hg) x), fderiv_fun_add ((hc 0 _ hf) x) ((hc 0 _ hg) x)]
    simp only [ContinuousLinearMap.add_apply]; ring
  · unfold div2 pd
    rw [hsplit 0, hsplit 1, fderiv_fun_add ((hc 0 _ hf) x) ((hc 0 _ hg) x), fderiv_fun_add ((hc 1 _ hf) x) ((hc 1 _ hg) x)]
    simp only [ContinuousLinearMap.add_apply]; ring

theorem designed_classical_aux_stage {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ) (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T)
    (τ : ℝ) (hτ : τ ∈ Ico 0 T.Tstar) : ∃ k : ℕ, τ ∈ Ico (T.t k) (T.t (k + 1)) :=
  S7.readout_grad_aux_stage_exists NStar Cg T hT τ hτ

theorem designed_classical_aux_finsum {F' : Type*} [AddCommMonoid F'] [TopologicalSpace F']
    {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ) (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) (k : ℕ) (τ : ℝ)
    (hτ : τ ∈ Ico (T.t k) (T.t (k + 1))) (f : ℕ → F') :
    (∑' m, if T.tIn m ≤ τ then f m else 0) = ∑ m ∈ Finset.range (k + 2), f m := by
  have hmono := designed_classical_aux_mono NStar Cg T hT
  rw [tsum_eq_sum (s := Finset.range (k + 2)) (fun m hm => ?_)]
  · refine Finset.sum_congr rfl fun m hm => ?_
    rw [if_pos]
    unfold TowerData.tIn
    exact (hmono (show m - 1 ≤ k by have := Finset.mem_range.1 hm; omega)).trans hτ.1
  · rw [if_neg]
    unfold TowerData.tIn
    exact not_le.2 (lt_of_lt_of_le hτ.2 (hmono (by have := Finset.mem_range.not.1 hm; omega)))

set_option maxHeartbeats 800000 in
theorem designed_classical_aux_fields {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ)
    (hNuni : ∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) (hP : PotentialTheoryR2)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) (τ : ℝ) (hτ : τ ∈ Ico 0 T.Tstar) :
    ContDiff ℝ ∞ (S7.Designed.vel T τ) ∧ (∀ x, curl2 (S7.Designed.vel T τ) x = S7.Designed.omega T τ x) ∧
    (∀ x, div2 (S7.Designed.vel T τ) x = 0) ∧ MemLp (S7.Designed.vel T τ) 2 volume ∧
    ContDiff ℝ ∞ (S7.Designed.theta T τ) ∧ HasCompactSupport (S7.Designed.theta T τ) ∧
    MemLp (S7.Designed.theta T τ) 2 volume := by
  classical
  obtain ⟨k, hk⟩ := designed_classical_aux_stage NStar Cg T hT τ hτ
  have hτ0 : 0 ≤ τ := hτ.1
  have hL := fun m (hm1 : 1 ≤ m) (hmk : m ≤ k + 1) => designed_classical_aux_layer NStar Cg hNuni hP T hT k τ hk m hm1 hmk
  have hvel : S7.Designed.vel T τ = fun x => T.uSt τ x + ∑ m ∈ Finset.range (k + 2), T.uL m τ x := by
    funext x; unfold S7.Designed.vel; rw [designed_classical_aux_finsum NStar Cg T hT k τ hk]
  have hom : S7.Designed.omega T τ = fun x => T.omSt τ x + ∑ m ∈ Finset.range (k + 2), T.omL m τ x := by
    funext x; unfold S7.Designed.omega; rw [designed_classical_aux_finsum NStar Cg T hT k τ hk]
  have hth : S7.Designed.theta T τ = fun x => ∑ m ∈ Finset.range (k + 2), T.θL m τ x := by
    funext x; unfold S7.Designed.theta; rw [designed_classical_aux_finsum NStar Cg T hT k τ hk]
  have hH : ContDiff ℝ ∞ Pr.H := Pr.hH.smooth
  have hpdH : ∀ i : Fin 2, ContDiff ℝ ∞ (pd i Pr.H) := fun i =>
    (hH.fderiv_right (m := ∞) (by norm_cast)).clm_apply contDiff_const
  have hpg_sm : ContDiff ℝ ∞ (perpGrad Pr.H) := by
    unfold perpGrad
    exact ((hpdH 1).neg.smul contDiff_const).add ((hpdH 0).smul contDiff_const)
  have hst_sm : ContDiff ℝ ∞ (T.uSt τ) := by
    show ContDiff ℝ ∞ fun x => deriv T.ψ τ • perpGrad Pr.H x
    exact (contDiff_const (c := deriv T.ψ τ)).smul hpg_sm
  have hst_curl : ∀ x, curl2 (T.uSt τ) x = T.omSt τ x := fun x =>
    S7.designed_vorticity_at_parking_aux_curl_smul_perpGrad _ _ x hH
  have hst_div : ∀ x, div2 (T.uSt τ) x = 0 := fun x => designed_classical_aux_div_perpGrad _ _ hH x
  have hst_cs : HasCompactSupport (T.uSt τ) := by
    have hHc : HasCompactSupport Pr.H := Pr.hH.compact
    have hpg : HasCompactSupport (perpGrad Pr.H) := by
      refine HasCompactSupport.of_support_subset_isCompact hHc ?_
      intro x hx
      rw [Function.mem_support] at hx
      by_contra hnot
      apply hx
      have h0 : ∀ i : Fin 2, pd i Pr.H x = 0 := by
        intro i
        have : fderiv ℝ Pr.H x = 0 := fderiv_of_notMem_tsupport ℝ hnot
        simp [pd, this]
      simp [perpGrad, h0]
    refine HasCompactSupport.of_support_subset_isCompact hpg ?_
    intro x hx
    rw [Function.mem_support] at hx
    apply subset_tsupport
    rw [Function.mem_support]
    intro h0
    apply hx
    show deriv T.ψ τ • perpGrad Pr.H x = 0
    rw [h0, smul_zero]
  have hb_u : T.uL 0 τ = perpGrad (T.basePsi τ) := funext fun x => hT.base.u0 τ hτ0 x
  have hψb : ContDiff ℝ ∞ (T.basePsi τ) := by
    have h1 : ContDiff ℝ ∞ fun x : R2 => Real.cos (inner ℝ (T.xi 0 τ) x) :=
      Real.contDiff_cos.comp (innerSL ℝ (T.xi 0 τ)).contDiff
    have h2 : ContDiff ℝ ∞ fun x : R2 => Pr.χ0 (rot (-(T.ψ τ)) x) :=
      Pr.hχ.smooth.comp (S7.designed_vorticity_at_parking_aux_rot_contDiff _)
    have : T.basePsi τ = fun x => -(T.Om0 τ / T.base.N0 ^ 2) * Real.cos (inner ℝ (T.xi 0 τ) x) *
        Pr.χ0 (rot (-(T.ψ τ)) x) := rfl
    rw [this]
    exact (contDiff_const.mul h1).mul h2
  have hpdb : ∀ i : Fin 2, ContDiff ℝ ∞ (pd i (T.basePsi τ)) := fun i =>
    (hψb.fderiv_right (m := ∞) (by norm_cast)).clm_apply contDiff_const
  have hb_sm : ContDiff ℝ ∞ (T.uL 0 τ) := by
    rw [hb_u]; unfold perpGrad
    exact ((hpdb 1).neg.smul contDiff_const).add ((hpdb 0).smul contDiff_const)
  have hb_curl : ∀ x, curl2 (T.uL 0 τ) x = T.omL 0 τ x := fun x => (hT.base.om0 τ hτ0 x).symm
  have hb_div : ∀ x, div2 (T.uL 0 τ) x = 0 := by
    intro x; rw [hb_u]
    have := designed_classical_aux_div_perpGrad 1 (T.basePsi τ) hψb x
    simpa using this
  obtain ⟨Rχ0, hRχ0⟩ : ∃ r, tsupport Pr.χ0 ⊆ Metric.closedBall (0 : R2) r :=
    (Pr.hχ.compact : IsCompact (tsupport Pr.χ0)).isBounded.subset_closedBall 0
  have hcut0 : ∀ (G : R2 → ℝ), (∀ x, Pr.χ0 (rot (-(T.ψ τ)) x) = 0 → G x = 0) → HasCompactSupport G := by
    intro G hG
    refine HasCompactSupport.of_support_subset_isCompact (isCompact_closedBall (0 : R2) Rχ0) ?_
    intro x hx
    rw [Function.mem_support] at hx
    have hχx : Pr.χ0 (rot (-(T.ψ τ)) x) ≠ 0 := fun h => hx (hG x h)
    have := hRχ0 (subset_tsupport _ (Function.mem_support.2 hχx))
    rw [Metric.mem_closedBall, dist_zero_right, S7.parked_node_vorticity_aux_rot_norm] at this
    rwa [Metric.mem_closedBall, dist_zero_right]
  have hb_θcs : HasCompactSupport (T.θL 0 τ) := by
    apply hcut0; intro x hx; rw [hT.base.theta0 τ hτ0 x, hx, mul_zero]
  have hb_θsm : ContDiff ℝ ∞ (T.θL 0 τ) := by
    have : T.θL 0 τ = fun x => T.base.Theta0 * Real.sin (inner ℝ (T.xi 0 τ) x) * Pr.χ0 (rot (-(T.ψ τ)) x) :=
      funext fun x => hT.base.theta0 τ hτ0 x
    rw [this]
    exact (contDiff_const.mul (Real.contDiff_sin.comp (innerSL ℝ (T.xi 0 τ)).contDiff)).mul
      (Pr.hχ.smooth.comp (S7.designed_vorticity_at_parking_aux_rot_contDiff _))
  have hb_ψcs : HasCompactSupport (T.basePsi τ) := by
    apply hcut0; intro x hx
    show -(T.Om0 τ / T.base.N0 ^ 2) * Real.cos (inner ℝ (T.xi 0 τ) x) * Pr.χ0 (rot (-(T.ψ τ)) x) = 0
    rw [hx, mul_zero]
  have hb_ucs : HasCompactSupport (T.uL 0 τ) := by
    rw [hb_u]
    refine HasCompactSupport.of_support_subset_isCompact hb_ψcs ?_
    intro x hx
    rw [Function.mem_support] at hx
    by_contra hnot
    apply hx
    have h0 : ∀ i : Fin 2, pd i (T.basePsi τ) x = 0 := by
      intro i
      have : fderiv ℝ (T.basePsi τ) x = 0 := fderiv_of_notMem_tsupport ℝ hnot
      simp [pd, this]
    simp [perpGrad, h0]
  have hall : ∀ m ∈ Finset.range (k + 2), ContDiff ℝ ∞ (T.uL m τ) ∧ (∀ x, curl2 (T.uL m τ) x = T.omL m τ x) ∧
      (∀ x, div2 (T.uL m τ) x = 0) ∧ MemLp (T.uL m τ) 2 volume ∧ ContDiff ℝ ∞ (T.θL m τ) ∧
      HasCompactSupport (T.θL m τ) := by
    intro m hm
    have hmk : m ≤ k + 1 := by have := Finset.mem_range.1 hm; omega
    rcases Nat.eq_zero_or_pos m with h0 | hm1
    · subst h0
      exact ⟨hb_sm, hb_curl, hb_div, hb_sm.continuous.memLp_of_hasCompactSupport hb_ucs, hb_θsm, hb_θcs⟩
    · obtain ⟨h1, h2, h3, h4, h5, h6⟩ := hL m hm1 hmk
      exact ⟨h1, h4, h5, h6, h2, h3⟩
  have hdiff : ∀ m ∈ Finset.range (k + 2), Differentiable ℝ (T.uL m τ) := fun m hm =>
    (hall m hm).1.differentiable (by simp)
  have hSdiff : Differentiable ℝ (fun y => ∑ m ∈ Finset.range (k + 2), T.uL m τ y) := by
    have : (fun y => ∑ m ∈ Finset.range (k + 2), T.uL m τ y) = ∑ m ∈ Finset.range (k + 2), T.uL m τ := by funext y; simp
    rw [this]; exact Differentiable.sum fun m hm => hdiff m hm
  have hstd : Differentiable ℝ (T.uSt τ) := hst_sm.differentiable (by simp)
  have hθcs : HasCompactSupport (S7.Designed.theta T τ) := by
    rw [hth]
    have key : ∀ s : Finset ℕ, (∀ m ∈ s, HasCompactSupport (T.θL m τ)) →
        HasCompactSupport (fun x => ∑ m ∈ s, T.θL m τ x) := by
      intro s hs
      induction s using Finset.induction_on with
      | empty => simp [HasCompactSupport, tsupport]
      | insert a s ha ih =>
        have h1 : HasCompactSupport (T.θL a τ) := hs a (Finset.mem_insert_self a s)
        have h2 := ih fun m hm => hs m (Finset.mem_insert_of_mem hm)
        have : (fun x => ∑ m ∈ insert a s, T.θL m τ x) = fun x => T.θL a τ x + ∑ m ∈ s, T.θL m τ x := by
          funext x; rw [Finset.sum_insert ha]
        rw [this]; exact h1.add h2
    exact key _ fun m hm => (hall m hm).2.2.2.2.2
  have hθsm : ContDiff ℝ ∞ (S7.Designed.theta T τ) := by rw [hth]; exact ContDiff.sum fun m hm => (hall m hm).2.2.2.2.1
  refine ⟨?_, fun x => ?_, fun x => ?_, ?_, hθsm, hθcs, hθsm.continuous.memLp_of_hasCompactSupport hθcs⟩
  · rw [hvel]; exact hst_sm.add (ContDiff.sum fun m hm => (hall m hm).1)
  · rw [hvel, hom, (designed_classical_aux_curl_add2 _ _ x hstd hSdiff).1,
      (designed_classical_aux_curl_finsum _ _ x hdiff).1, hst_curl x]
    simp only
    congr 1
    exact Finset.sum_congr rfl fun m hm => (hall m hm).2.1 x
  · rw [hvel, (designed_classical_aux_curl_add2 _ _ x hstd hSdiff).2, (designed_classical_aux_curl_finsum _ _ x hdiff).2,
      hst_div x, zero_add]
    exact Finset.sum_eq_zero fun m hm => (hall m hm).2.2.1 x
  · rw [hvel]
    refine (hst_sm.continuous.memLp_of_hasCompactSupport hst_cs).add ?_
    have : (fun x => ∑ m ∈ Finset.range (k + 2), T.uL m τ x) = ∑ m ∈ Finset.range (k + 2), T.uL m τ := by funext y; simp
    rw [this]
    exact memLp_finsetSum' _ fun m hm => (hall m hm).2.2.2.1
theorem designed_classical {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ)
    (hNuni : ∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) (hκ : 0 < κ) (hP : PotentialTheoryR2)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T)
    (fu : ℝ → R2 → R2) (hfu : ∀ t ∈ Ico 0 T.Tstar, ∀ x, curl2 (fu t) x = S7.Designed.fOmega T t x)
    (hfu_sm : ∀ i : Fin 2, ContDiffOn ℝ ∞ (fun z : ℝ × R2 => fu z.1 z.2 i) (Ico 0 T.Tstar ×ˢ univ))
    (hθ : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => S7.Designed.theta T z.1 z.2) (Ico 0 T.Tstar ×ˢ univ))
    (hu : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => S7.Designed.vel T z.1 z.2) (Ico 0 T.Tstar ×ˢ univ))
    (hfθ : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => S7.Designed.fTheta T z.1 z.2) (Ico 0 T.Tstar ×ˢ univ)) :
    ∃ p : ℝ → R2 → ℝ,
      ClassicalBoussinesq κ T.Tstar (S7.Designed.theta T) p (S7.Designed.vel T) fu (S7.Designed.fTheta T) := by
  classical
  have hred1 : T.StateRed 1 := (hT.state 1 le_rfl).1
  have hT0 : 0 < T.Tstar := by
    obtain ⟨hbdd, -, -, -⟩ := S7.Tstar_le NStar Cg T hT
    have h1 : T.t 1 ≤ T.Tstar := le_ciSup hbdd 1
    have : T.t 0 < T.t 1 := hred1.t_mono (show (0:ℕ) ≤ 1 by norm_num) (show (1:ℕ) ≤ 1 from le_rfl) (by norm_num)
    rw [hred1.t_zero] at this
    linarith
  set θ := S7.Designed.theta T with hθdef
  set u := S7.Designed.vel T with hudef
  set Om := S7.Designed.omega T with hOmdef
  have hF := fun t (ht : t ∈ Ico 0 T.Tstar) => designed_classical_aux_fields NStar Cg hNuni hP T hT t ht
  have slice : ∀ {F' : Type} [NormedAddCommGroup F'] [NormedSpace ℝ F'] (G : ℝ × R2 → F'),
      ContDiffOn ℝ ∞ G (Ico 0 T.Tstar ×ˢ univ) → ∀ t ∈ Ico 0 T.Tstar, ContDiff ℝ ∞ (fun x => G (t, x)) := by
    intro F' _ _ G hG t ht
    have := hG.comp (contDiffOn_const.prodMk contDiffOn_id) (fun (y : R2) (_ : y ∈ univ) => (⟨ht, mem_univ _⟩ : (t, y) ∈ Ico 0 T.Tstar ×ˢ univ))
    exact contDiffOn_univ.1 this
  have huc : ∀ i : Fin 2, ContDiffOn ℝ ∞ (fun z : ℝ × R2 => u z.1 z.2 i) (Ico 0 T.Tstar ×ˢ univ) :=
    fun i => contDiffOn_euclidean.1 hu i
  have slab : ∀ (G : ℝ × R2 → ℝ), ContDiffOn ℝ ∞ G (Ico 0 T.Tstar ×ˢ univ) → ∀ T₂ < T.Tstar,
      ContDiffOn ℝ ∞ G (Icc 0 T₂ ×ˢ univ) := fun G hG T₂ hT₂ =>
    hG.mono (prod_mono (Icc_subset_Ico_right hT₂) subset_rfl)
  set A : ℝ → R2 → R2 := fun t x => (κ * θ t x) • (e 1 : R2) with hA
  set Pv : ℝ → R2 → R2 := fun t x => WithLp.toLp 2 fun i : Fin 2 => pdt (fun τ y => (u τ y) i) t x with hPv
  set Av : ℝ → R2 → R2 := fun t x => WithLp.toLp 2 fun i : Fin 2 => advect (u t) (fun y => (u t y) i) x with hAv
  set g : ℝ → R2 → R2 := fun t x => A t x + fu t x - (Pv t x + Av t x) with hg
  have hA_sm : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => A z.1 z.2) (Ico 0 T.Tstar ×ˢ univ) := by
    show ContDiffOn ℝ ∞ (fun z : ℝ × R2 => (κ * θ z.1 z.2) • (e 1 : R2)) _
    exact (contDiffOn_const.mul hθ).smul contDiffOn_const
  have hfu_sm : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => fu z.1 z.2) (Ico 0 T.Tstar ×ˢ univ) :=
    contDiffOn_euclidean.2 fun i => hfu_sm i
  have hPv_sm : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => Pv z.1 z.2) (Ico 0 T.Tstar ×ˢ univ) := by
    refine contDiffOn_euclidean.2 fun i => ?_
    show ContDiffOn ℝ ∞ (fun z : ℝ × R2 => pdt (fun τ y => (u τ y) i) z.1 z.2) _
    refine S7.lemma751_fields_aux_local T.Tstar (fun z : ℝ × R2 => pdt (fun τ y => (u τ y) i) z.1 z.2) fun τ' hτ'0 hτ' => ?_
    have hT₂ : (τ' + T.Tstar) / 2 < T.Tstar := by linarith
    exact S7.lemma751_fields_aux_pdt (fun τ y => (u τ y) i) ((τ' + T.Tstar) / 2) τ' hτ'0 (by linarith)
      (slab _ (huc i) _ hT₂)
  have hAv_sm : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => Av z.1 z.2) (Ico 0 T.Tstar ×ˢ univ) := by
    refine contDiffOn_euclidean.2 fun i => ?_
    show ContDiffOn ℝ ∞ (fun z : ℝ × R2 => advect (u z.1) (fun y => (u z.1 y) i) z.2) _
    refine S7.lemma751_fields_aux_local T.Tstar (fun z : ℝ × R2 => advect (u z.1) (fun y => (u z.1 y) i) z.2)
      fun τ' hτ'0 hτ' => ?_
    exact S7.lemma751_fields_aux_advect (fun τ y => (u τ y) i) u τ' (slab _ (huc i) _ hτ')
      (hu.mono (prod_mono (Icc_subset_Ico_right hτ') subset_rfl))
  have hg_sm : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => g z.1 z.2) (Ico 0 T.Tstar ×ˢ univ) :=
    (hA_sm.add hfu_sm).sub (hPv_sm.add hAv_sm)
  have hcurl : ∀ t ∈ Ico 0 T.Tstar, ∀ x, curl2 (g t) x = 0 := by
    intro t ht x
    obtain ⟨hu_t, hcurl_u, hdiv_u, -, hθ_t, -, -⟩ := hF t ht
    have hA_t : ContDiff ℝ ∞ (A t) := slice (fun z => A z.1 z.2) hA_sm t ht
    have hfu_t : ContDiff ℝ ∞ (fu t) := slice (fun z => fu z.1 z.2) hfu_sm t ht
    have hPv_t : ContDiff ℝ ∞ (Pv t) := slice (fun z => Pv z.1 z.2) hPv_sm t ht
    have hAv_t : ContDiff ℝ ∞ (Av t) := slice (fun z => Av z.1 z.2) hAv_sm t ht
    have d := fun {f : R2 → R2} (h : ContDiff ℝ ∞ f) => h.differentiable (by simp)
    have cA : curl2 (A t) x = κ * pd 0 (θ t) x := by
      have h1 : (fun y => (A t y) 1) = fun y => κ * θ t y := by funext y; simp [hA, e]
      have h0 : (fun y => (A t y) 0) = fun y => (0:ℝ) := by funext y; simp [hA, e]
      unfold curl2
      rw [h1, h0]
      have : pd 1 (fun _ : R2 => (0:ℝ)) x = 0 := by simp [pd]
      rw [this, sub_zero]
      unfold pd
      rw [fderiv_const_mul ((hθ_t.differentiable (by simp)) x)]; rfl
    have cfu : curl2 (fu t) x = S7.Designed.fOmega T t x := hfu t ht x
    have cAv : curl2 (Av t) x = advect (u t) (Om t) x := by
      rw [show Av t = fun y => WithLp.toLp 2 fun i : Fin 2 => advect (u t) (fun z => (u t z) i) y from rfl,
        designed_classical_aux_curl_advect (u t) hu_t hdiv_u x]
      congr 1
      exact funext hcurl_u
    have cPv : curl2 (Pv t) x = pdt Om t x := by
      set T₂ := (t + T.Tstar) / 2 with hT₂
      have htT₂ : t < T₂ := by rw [hT₂]; linarith [ht.2]
      have hT₂T : T₂ < T.Tstar := by rw [hT₂]; linarith [ht.2]
      have hT₂0 : 0 < T₂ := by rw [hT₂]; linarith [ht.1]
      have htI : t ∈ Icc 0 T₂ := ⟨ht.1, htT₂.le⟩
      have comm : ∀ i j : Fin 2, pd i (fun y => pdt (fun τ z => (u τ z) j) t y) x =
          derivWithin (fun s => pd i (fun z => (u s z) j) x) (Icc 0 T₂) t := by
        intro i j
        have hK := slab _ (huc j) T₂ hT₂T
        have h := Kit.iteratedFDerivWithin_eq_zero_of_flat_aux_comm1 0 T₂ hT₂0 (fun z : ℝ × R2 => u z.1 z.2 j) hK t htI x i
        have hset : (fun y => pdt (fun τ z => (u τ z) j) t y) = fun y => derivWithin (fun s => u s y j) (Icc 0 T₂) t := by
          funext y
          have := S7.timeReg_bgTimeData_aux_set2 (fun s => u s y j) T₂ t htT₂ 1
          simp only [iteratedDerivWithin_one] at this
          unfold pdt; exact this.symm
        rw [hset]
        exact h
      unfold curl2
      rw [comm 0 1, comm 1 0]
      have hdiffble : ∀ i j : Fin 2, DifferentiableWithinAt ℝ (fun s => pd i (fun z => (u s z) j) x) (Icc 0 T₂) t := by
        intro i j
        have h := (S7.lemma751_fields_aux_pd (fun s z => (u s z) j) T₂ (slab _ (huc j) T₂ hT₂T) i)
        have hsl : ContDiffOn ℝ ∞ (fun s => pd i (fun z => (u s z) j) x) (Icc 0 T₂) :=
          h.comp (contDiffOn_id.prodMk contDiffOn_const) (fun s hs => ⟨hs, mem_univ _⟩)
        exact (hsl.differentiableOn (by simp)) t htI
      rw [← derivWithin_fun_sub (hdiffble 0 1) (hdiffble 1 0)]
      have hset2 : derivWithin (fun s => pd 0 (fun z => (u s z) 1) x - pd 1 (fun z => (u s z) 0) x) (Icc 0 T₂) t =
          derivWithin (fun s => Om s x) (Icc 0 T₂) t := by
        apply derivWithin_congr
        · intro s hs
          have := (hF s ⟨hs.1, lt_of_le_of_lt hs.2 hT₂T⟩).2.1 x
          exact this
        · exact (hF t ht).2.1 x
      rw [hset2]
      have := S7.timeReg_bgTimeData_aux_set2 (fun s => Om s x) T₂ t htT₂ 1
      simp only [iteratedDerivWithin_one] at this
      unfold pdt; exact this
    have hsum : curl2 (g t) x = curl2 (A t) x + curl2 (fu t) x - (curl2 (Pv t) x + curl2 (Av t) x) := by
      have h1 := (designed_classical_aux_curl_add2 (A t) (fu t) x (d hA_t) (d hfu_t)).1
      have h2 := (designed_classical_aux_curl_add2 (Pv t) (Av t) x (d hPv_t) (d hAv_t)).1
      have hsub : ∀ (f k : R2 → R2), Differentiable ℝ f → Differentiable ℝ k →
          curl2 (fun y => f y - k y) x = curl2 f x - curl2 k x := by
        intro f k hf hk
        have := (designed_classical_aux_curl_add2 (fun y => f y - k y) k x (hf.sub hk) hk).1
        simp only [sub_add_cancel] at this
        linarith
      have h3 := hsub (fun y => A t y + fu t y) (fun y => Pv t y + Av t y) ((d hA_t).add (d hfu_t)) ((d hPv_t).add (d hAv_t))
      have hgt : g t = fun y => (A t y + fu t y) - (Pv t y + Av t y) := rfl
      rw [hgt, h3, h1, h2]
    rw [hsum, cA, cfu, cPv, cAv]
    unfold S7.Designed.fOmega
    ring
  obtain ⟨p, hp_sm, hp_grad⟩ := Kit.gradient_of_curl_free_param (Ico 0 T.Tstar) (convex_Ico 0 T.Tstar)
    ⟨T.Tstar / 2, by rw [interior_Ico]; exact ⟨by linarith, by linarith⟩⟩ g hg_sm hcurl
  refine ⟨p, ⟨hθ, hu, hp_sm, fun t ht => (hF t ht).2.2.1, fun t ht => ⟨(hF t ht).2.2.2.1, (hF t ht).2.2.2.2.2.2⟩,
    fun t ht x => rfl, fun t ht x i => ?_⟩⟩
  have hgi := hp_grad t ht x i
  have : g t x i = κ * θ t x * (if i = 1 then 1 else 0) + fu t x i -
      (pdt (fun τ y => (u τ y) i) t x + advect (u t) (fun y => (u t y) i) x) := by
    fin_cases i <;> simp [hg, hA, hPv, hAv, e]
  rw [this] at hgi
  show pdt (fun τ y => (S7.Designed.vel T τ y) i) t x + advect (S7.Designed.vel T t) (fun y => (S7.Designed.vel T t y) i) x +
      pd i (p t) x = κ * S7.Designed.theta T t x * (if i = 1 then 1 else 0) + fu t x i
  rw [hgi]; ring

theorem designed_parity_aux_pd_odd (g : R2 → ℝ) (hg : ∀ x, g (-x) = -g x) (i : Fin 2) (x : R2) :
    pd i g (-x) = pd i g x := by
  have h1 := S0.lemma02_step4_pd_reflect_aux g i x
  have h2 : (fun y => g (-y)) = -g := by funext y; rw [hg]; rfl
  have hneg : pd i (-g) x = -pd i g x := by unfold pd; rw [fderiv_neg]; rfl
  rw [h2, hneg] at h1
  have : pd i g (- -x) = pd i g x := by rw [neg_neg]
  linarith

theorem designed_parity_aux_rot_neg (a : ℝ) (x : R2) : rot a (-x) = -rot a x := by
  unfold rot; simp only [PiLp.neg_apply]
  rw [show (Real.cos a * -x 0 - Real.sin a * -x 1) = -(Real.cos a * x 0 - Real.sin a * x 1) by ring,
    show (Real.sin a * -x 0 + Real.cos a * -x 1) = -(Real.sin a * x 0 + Real.cos a * x 1) by ring,
    neg_smul, neg_smul, neg_add]

theorem designed_parity_aux_base {κ : ℝ} {Pr : Profiles} (T : TowerData κ Pr) (hB : T.BaseOK) (s : ℝ) (hs : 0 ≤ s) :
    (∀ x, T.θL 0 s (-x) = -T.θL 0 s x) ∧ (∀ x, T.uL 0 s (-x) = -T.uL 0 s x) := by
  have hΨeven : ∀ x : R2, T.basePsi s (-x) = T.basePsi s x := by
    intro x; unfold TowerData.basePsi; rw [inner_neg_right, Real.cos_neg, designed_parity_aux_rot_neg, Pr.hχ.even]
  refine ⟨fun x => ?_, fun x => ?_⟩
  · rw [hB.theta0 s hs, hB.theta0 s hs, inner_neg_right, Real.sin_neg, designed_parity_aux_rot_neg, Pr.hχ.even]; ring
  · rw [hB.u0 s hs (-x), hB.u0 s hs x]; exact S6.wf_base_aux_perpGrad_odd _ hΨeven x

theorem designed_parity_aux_layer {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) (k : ℕ) (τ : ℝ) (hτ : τ ∈ Ico (T.t k) (T.t (k + 1)))
    (m : ℕ) (hm1 : 1 ≤ m) (hmk : m ≤ k + 1) :
    (∀ x, T.θL m τ (-x) = -T.θL m τ x) ∧ (∀ x, T.uL m τ (-x) = -T.uL m τ x) := by
  have hmono := designed_classical_aux_mono NStar Cg T hT
  have hτ2 : τ ≤ T.t (k + 2) := hτ.2.le.trans (hmono (by omega))
  have ht' : T.t (k + 2) ∈ T.I (k + 2) := ⟨le_rfl, hmono (by omega)⟩
  obtain ⟨-, -, ⟨R, hKin, -, K, hS4, -, -⟩, -⟩ := hT.older (k + 2) (by omega) (T.t (k + 2)) ht' m hm1 (by omega)
  obtain ⟨hDsel, hrun⟩ := hT.run (k + 1) (by omega)
  have htin : T.tIn m ≤ τ := (hmono (show m - 1 ≤ k by omega)).trans hτ.1
  have hera : τ ∈ Icc (T.bg m (T.t (k + 2))).tIn (T.bg m (T.t (k + 2))).τ := ⟨htin, hτ2⟩
  have hreal := hrun.realised m hm1 (by omega) R hKin τ ⟨htin, hτ2⟩
  have hpar := L.layer_parity Pr (T.P m) R K (T.N (m - 1)) hS4 τ hera
  refine ⟨fun x => ?_, fun x => ?_⟩
  · rw [(hreal (-x)).1, (hreal x).1]; exact hpar.1 x
  · rw [(hreal (-x)).2.2, (hreal x).2.2]; exact hpar.2.2 x

theorem designed_parity_aux_fields {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) (τ : ℝ) (hτ : τ ∈ Ico 0 T.Tstar) :
    InClassP (S7.Designed.theta T τ) (S7.Designed.vel T τ) := by
  classical
  obtain ⟨k, hk⟩ := designed_classical_aux_stage NStar Cg T hT τ hτ
  have hth : ∀ x, S7.Designed.theta T τ x = ∑ m ∈ Finset.range (k + 2), T.θL m τ x := by
    intro x; unfold S7.Designed.theta; rw [designed_classical_aux_finsum NStar Cg T hT k τ hk]
  have hvel : ∀ x, S7.Designed.vel T τ x = T.uSt τ x + ∑ m ∈ Finset.range (k + 2), T.uL m τ x := by
    intro x; unfold S7.Designed.vel; rw [designed_classical_aux_finsum NStar Cg T hT k τ hk]
  have hlay : ∀ m ∈ Finset.range (k + 2),
      (∀ x, T.θL m τ (-x) = -T.θL m τ x) ∧ (∀ x, T.uL m τ (-x) = -T.uL m τ x) := by
    intro m hm
    have hmk : m ≤ k + 1 := by have := Finset.mem_range.1 hm; omega
    rcases Nat.eq_zero_or_pos m with h0 | hpos
    · rw [h0]; exact designed_parity_aux_base T hT.base τ hτ.1
    · exact designed_parity_aux_layer NStar Cg T hT k τ hk m hpos hmk
  have hHeven : ∀ x : R2, Pr.H (-x) = Pr.H x := fun x => Pr.hH.radial _ _ (norm_neg x)
  have hst : ∀ x, T.uSt τ (-x) = -T.uSt τ x := by
    intro x; unfold TowerData.uSt; rw [S6.wf_base_aux_perpGrad_odd _ hHeven, smul_neg]
  refine ⟨fun x => ?_, fun x => ?_⟩
  · rw [hth, hth, ← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl fun m hm => (hlay m hm).1 x
  · rw [hvel, hvel, neg_add, ← Finset.sum_neg_distrib, hst]
    congr 1
    exact Finset.sum_congr rfl fun m hm => (hlay m hm).2 x

theorem designed_parity {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ) (Kexp : ℕ)
    (hNdom : ∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg Kexp) Cg β ≤ NStar β)
    (hNuni : ∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T)
    (fu : ℝ → R2 → R2) (hfu_odd : ∀ t ∈ Ico 0 T.Tstar, ∀ x, fu t (-x) = -fu t x) :
    (∀ t ∈ Ico 0 T.Tstar, InClassP (S7.Designed.theta T t) (S7.Designed.vel T t)) ∧
    (∀ t ∈ Ico 0 T.Tstar, InClassP (S7.Designed.fTheta T t) (fu t)) := by
  have h1 : ∀ t ∈ Ico 0 T.Tstar, InClassP (S7.Designed.theta T t) (S7.Designed.vel T t) :=
    fun t ht => designed_parity_aux_fields NStar Cg T hT t ht
  refine ⟨h1, fun t ht => ⟨fun x => ?_, hfu_odd t ht⟩⟩
  have hθ := (h1 t ht).1
  have hu := (h1 t ht).2
  have hpdt : pdt (S7.Designed.theta T) t (-x) = -pdt (S7.Designed.theta T) t x := by
    unfold pdt
    have hev : (fun τ => S7.Designed.theta T τ (-x)) =ᶠ[𝓝[Ici 0] t] (-fun τ => S7.Designed.theta T τ x) := by
      filter_upwards [inter_mem_nhdsWithin (Ici (0:ℝ)) (Iio_mem_nhds ht.2)] with τ hτ
      exact (h1 τ ⟨hτ.1, hτ.2⟩).1 x
    rw [hev.derivWithin_eq (hθ x), derivWithin.neg]
  have hadv : advect (S7.Designed.vel T t) (S7.Designed.theta T t) (-x) =
      -advect (S7.Designed.vel T t) (S7.Designed.theta T t) x := by
    unfold advect
    rw [hu x, designed_parity_aux_pd_odd _ hθ, designed_parity_aux_pd_odd _ hθ]
    simp only [PiLp.neg_apply]; ring
  unfold S7.Designed.fTheta
  rw [hpdt, hadv]; ring

theorem designed_lipschitz_aux_cmul {U : Set (ℝ × R2)} (hU : IsOpen U) (f g : ℝ × R2 → ℝ)
    (hf : ContinuousOn f U) (hg : Continuous g) (hgs : tsupport g ⊆ U) :
    Continuous (fun z => f z * g z) := by
  rw [continuous_iff_continuousAt]
  intro z
  by_cases hz : z ∈ U
  · exact ((hf z hz).continuousAt (hU.mem_nhds hz)).mul hg.continuousAt
  · have hz' : z ∉ tsupport g := fun h => hz (hgs h)
    have h0 : g =ᶠ[𝓝 z] 0 := notMem_tsupport_iff_eventuallyEq.1 hz'
    have h1 : (fun _ : ℝ × R2 => (0 : ℝ)) =ᶠ[𝓝 z] fun w => f w * g w := by
      filter_upwards [h0] with w hw
      simp [hw]
    exact continuousAt_const.congr h1

theorem designed_lipschitz_aux_integrable {U : Set (ℝ × R2)} (hU : IsOpen U) (f g : ℝ × R2 → ℝ)
    (hf : ContinuousOn f U) (hg : Continuous g) (hgc : HasCompactSupport g) (hgs : tsupport g ⊆ U) :
    Integrable (fun z => f z * g z) := by
  refine (designed_lipschitz_aux_cmul hU f g hf hg hgs).integrable_of_hasCompactSupport ?_
  exact (hgc.mul_left : HasCompactSupport (f * g))

theorem designed_lipschitz_aux_pdt (g : ℝ → R2 → ℝ) (t : ℝ) (x : R2) (ht : 0 < t)
    (hg : DifferentiableAt ℝ (fun z : ℝ × R2 => g z.1 z.2) (t, x)) :
    pdt g t x = fderiv ℝ (fun z : ℝ × R2 => g z.1 z.2) (t, x) ((1 : ℝ), (0 : R2)) := by
  have hc : HasDerivAt (fun τ : ℝ => ((τ, x) : ℝ × R2)) ((1 : ℝ), (0 : R2)) t :=
    (hasDerivAt_id t).prodMk (hasDerivAt_const t x)
  have h := hg.hasFDerivAt.comp_hasDerivAt t hc
  have h' : HasDerivAt (fun τ => g τ x) (fderiv ℝ (fun z : ℝ × R2 => g z.1 z.2) (t, x) ((1 : ℝ), (0 : R2))) t := h
  unfold pdt
  rw [derivWithin_of_mem_nhds (Ici_mem_nhds ht), h'.deriv]

theorem designed_lipschitz_aux_pdx (g : ℝ → R2 → ℝ) (t : ℝ) (x : R2) (j : Fin 2)
    (hg : DifferentiableAt ℝ (fun z : ℝ × R2 => g z.1 z.2) (t, x)) :
    pd j (g t) x = fderiv ℝ (fun z : ℝ × R2 => g z.1 z.2) (t, x) ((0 : ℝ), e j) := by
  have hc : HasFDerivAt (fun y : R2 => ((t, y) : ℝ × R2)) (ContinuousLinearMap.inr ℝ ℝ R2) x :=
    ((hasFDerivAt_const t x).prodMk (hasFDerivAt_id x))
  have h := hg.hasFDerivAt.comp x hc
  have h' : HasFDerivAt (g t) ((fderiv ℝ (fun z : ℝ × R2 => g z.1 z.2) (t, x)).comp
      (ContinuousLinearMap.inr ℝ ℝ R2)) x := h
  unfold pd
  rw [h'.fderiv]
  rfl

theorem designed_lipschitz_aux_test (φ : ℝ × R2 → ℝ) (hφ : ContDiff ℝ ∞ φ) (t : ℝ) (x : R2) (j : Fin 2) :
    deriv (fun τ => φ (τ, x)) t = fderiv ℝ φ (t, x) ((1 : ℝ), (0 : R2)) ∧
    pd j (fun y => φ (t, y)) x = fderiv ℝ φ (t, x) ((0 : ℝ), e j) := by
  have hd : DifferentiableAt ℝ φ (t, x) := (hφ.differentiable (by norm_cast)).differentiableAt
  constructor
  · have hc : HasDerivAt (fun τ : ℝ => ((τ, x) : ℝ × R2)) ((1 : ℝ), (0 : R2)) t :=
      (hasDerivAt_id t).prodMk (hasDerivAt_const t x)
    exact (hd.hasFDerivAt.comp_hasDerivAt t hc).deriv
  · have hc : HasFDerivAt (fun y : R2 => ((t, y) : ℝ × R2)) (ContinuousLinearMap.inr ℝ ℝ R2) x :=
      ((hasFDerivAt_const t x).prodMk (hasFDerivAt_id x))
    have h' : HasFDerivAt (fun y => φ (t, y)) ((fderiv ℝ φ (t, x)).comp (ContinuousLinearMap.inr ℝ ℝ R2)) x :=
      hd.hasFDerivAt.comp x hc
    unfold pd; rw [h'.fderiv]; rfl

theorem designed_lipschitz_aux_supp {E' F' : Type*} [Zero E'] [Zero F'] [TopologicalSpace E']
    (g : ℝ × R2 → E') (h : ℝ × R2 → F') (h0 : ∀ z, g z = 0 → h z = 0) :
    tsupport h ⊆ tsupport g ∧ (HasCompactSupport g → HasCompactSupport h) := by
  have hs : Function.support h ⊆ Function.support g := fun z hz hg => hz (h0 z hg)
  exact ⟨closure_mono hs, fun hc => hc.mono hs⟩

theorem designed_lipschitz_aux_prodDiff {U : Set (ℝ × R2)} (hU : IsOpen U) (f g : ℝ × R2 → ℝ)
    (hf : DifferentiableOn ℝ f U) (hg : Differentiable ℝ g) (hgs : tsupport g ⊆ U) (z v : ℝ × R2) :
    DifferentiableAt ℝ (fun w => f w * g w) z ∧
      fderiv ℝ (fun w => f w * g w) z v = fderiv ℝ f z v * g z + f z * fderiv ℝ g z v := by
  by_cases hz : z ∈ U
  · have hfd : DifferentiableAt ℝ f z := hf.differentiableAt (hU.mem_nhds hz)
    refine ⟨hfd.mul (hg z), ?_⟩
    rw [fderiv_fun_mul hfd (hg z)]
    simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.smul_apply, smul_eq_mul]
    ring
  · have hz' : z ∉ tsupport g := fun h => hz (hgs h)
    have h0 : g =ᶠ[𝓝 z] 0 := notMem_tsupport_iff_eventuallyEq.1 hz'
    have h1 : (fun w => f w * g w) =ᶠ[𝓝 z] fun _ => (0 : ℝ) := by
      filter_upwards [h0] with w hw
      simp [hw]
    refine ⟨(h1.differentiableAt_iff).2 (differentiableAt_const _), ?_⟩
    rw [h1.fderiv_eq, fderiv_of_notMem_tsupport ℝ hz', image_eq_zero_of_notMem_tsupport hz']
    simp

set_option maxHeartbeats 400000 in
theorem designed_lipschitz_aux_weak (κ Tst T' : ℝ) (hT' : T' ≤ Tst)
    (θ p : ℝ → R2 → ℝ) (u fu : ℝ → R2 → R2)
    (hu : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => u z.1 z.2) (Ico 0 Tst ×ˢ univ))
    (hp : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => p z.1 z.2) (Ico 0 Tst ×ˢ univ))
    (hθ : ContinuousOn (fun z : ℝ × R2 => θ z.1 z.2) (Ico 0 Tst ×ˢ univ))
    (hfu : ∀ i : Fin 2, ContinuousOn (fun z : ℝ × R2 => fu z.1 z.2 i) (Ico 0 Tst ×ˢ univ))
    (hdiv : ∀ t ∈ Ico 0 Tst, ∀ x, div2 (u t) x = 0)
    (hequ : ∀ t ∈ Ico 0 Tst, ∀ x, ∀ i : Fin 2,
      pdt (fun τ y => (u τ y) i) t x + advect (u t) (fun y => (u t y) i) x + pd i (p t) x
        = κ * θ t x * (if i = 1 then 1 else 0) + (fu t x) i)
    (φ : ℝ × R2 → R2) (hφ : ContDiff ℝ ∞ φ) (hφc : HasCompactSupport φ)
    (hφs : tsupport φ ⊆ Ioo 0 T' ×ˢ univ) (hφdiv : ∀ z : ℝ × R2, div2 (fun y => φ (z.1, y)) z.2 = 0) :
    ∫ z : ℝ × R2, ∑ i : Fin 2,
      ( (u z.1 z.2) i * deriv (fun τ => (φ (τ, z.2)) i) z.1
        + ∑ j : Fin 2, (u z.1 z.2) i * (u z.1 z.2) j * pd j (fun y => (φ (z.1, y)) i) z.2
        + (κ * θ z.1 z.2 * (if i = 1 then 1 else 0) + (fu z.1 z.2) i) * (φ z) i ) = 0 := by
  haveI hHaar : (volume : Measure (ℝ × R2)).IsAddHaarMeasure := Measure.prod.instIsAddHaarMeasure _ _
  set U : Set (ℝ × R2) := Ioo 0 Tst ×ˢ univ with hUdef
  have hUo : IsOpen U := isOpen_Ioo.prod isOpen_univ
  have hUS : U ⊆ Ico 0 Tst ×ˢ univ := prod_mono Ioo_subset_Ico_self subset_rfl
  set ui : Fin 2 → ℝ × R2 → ℝ := fun i z => u z.1 z.2 i with huidef
  set φi : Fin 2 → ℝ × R2 → ℝ := fun i z => φ z i with hφidef
  set P : ℝ × R2 → ℝ := fun z => p z.1 z.2 with hPdef
  set v0 : ℝ × R2 := ((1 : ℝ), (0 : R2)) with hv0
  set vx : Fin 2 → ℝ × R2 := fun j => ((0 : ℝ), e j) with hvx
  have hφi_sm : ∀ i, ContDiff ℝ ∞ (φi i) := fun i => contDiff_euclidean.1 hφ i
  have hφi_d : ∀ i, Differentiable ℝ (φi i) := fun i => (hφi_sm i).differentiable (by norm_cast)
  have hui_sm : ∀ i, ContDiffOn ℝ ∞ (ui i) U := fun i => (contDiffOn_euclidean.1 hu i).mono hUS
  have hui_d : ∀ i, DifferentiableOn ℝ (ui i) U := fun i => (hui_sm i).differentiableOn (by norm_cast)
  have hui_c : ∀ i, ContinuousOn (ui i) U := fun i => (hui_sm i).continuousOn
  have hP_sm : ContDiffOn ℝ ∞ P U := hp.mono hUS
  have hP_d : DifferentiableOn ℝ P U := hP_sm.differentiableOn (by norm_cast)
  have hDui_c : ∀ i (w : ℝ × R2), ContinuousOn (fun z => fderiv ℝ (ui i) z w) U := fun i w =>
    ((hui_sm i).continuousOn_fderiv_of_isOpen hUo (by norm_cast)).clm_apply continuousOn_const
  have hDP_c : ∀ (w : ℝ × R2), ContinuousOn (fun z => fderiv ℝ P z w) U := fun w =>
    (hP_sm.continuousOn_fderiv_of_isOpen hUo (by norm_cast)).clm_apply continuousOn_const
  have hφU : tsupport φ ⊆ U := hφs.trans (prod_mono (Ioo_subset_Ioo le_rfl hT') subset_rfl)
  have hφi_s : ∀ i, tsupport (φi i) ⊆ U ∧ HasCompactSupport (φi i) := fun i => by
    have h := designed_lipschitz_aux_supp φ (φi i) (fun z hz => by simp [hφidef, hz])
    exact ⟨h.1.trans hφU, h.2 hφc⟩
  have hDφi : ∀ i (w : ℝ × R2), Continuous (fun z => fderiv ℝ (φi i) z w) ∧
      tsupport (fun z => fderiv ℝ (φi i) z w) ⊆ U ∧ HasCompactSupport (fun z => fderiv ℝ (φi i) z w) := by
    intro i w
    have hc : Continuous (fun z => fderiv ℝ (φi i) z w) :=
      ((hφi_sm i).continuous_fderiv (by norm_cast)).clm_apply continuous_const
    have h := designed_lipschitz_aux_supp (fderiv ℝ (φi i)) (fun z => fderiv ℝ (φi i) z w)
      (fun z hz => by simp [hz])
    have h2 : tsupport (fun z => fderiv ℝ (φi i) z w) ⊆ tsupport (φi i) := tsupport_fderiv_apply_subset ℝ w
    exact ⟨hc, h2.trans (hφi_s i).1, h.2 ((hφi_s i).2.fderiv ℝ)⟩
  have hin : ∀ i z, φi i z ≠ 0 → z ∈ U := fun i z hz => (hφi_s i).1 (subset_tsupport _ hz)
  have hA : ∀ i, ∫ z, pdt (fun τ y => u τ y i) z.1 z.2 * φi i z = - ∫ z, ui i z * fderiv ℝ (φi i) z v0 := by
    intro i
    have h1 : (fun z : ℝ × R2 => pdt (fun τ y => u τ y i) z.1 z.2 * φi i z) =
        fun z => fderiv ℝ (ui i) z v0 * φi i z := by
      funext z
      by_cases hz : φi i z = 0
      · simp [hz]
      · have hzU := hin i z hz
        have hzt : 0 < z.1 := hzU.1.1
        have hd : DifferentiableAt ℝ (ui i) z := (hui_d i).differentiableAt (hUo.mem_nhds hzU)
        rw [designed_lipschitz_aux_pdt (fun τ y => u τ y i) z.1 z.2 hzt hd]
    rw [h1]
    have hibp := integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable (μ := volume) (v := v0)
      (f := ui i) (g := φi i)
      (designed_lipschitz_aux_integrable hUo _ _ (hDui_c i v0) (hφi_sm i).continuous (hφi_s i).2 (hφi_s i).1)
      (designed_lipschitz_aux_integrable hUo _ _ (hui_c i) (hDφi i v0).1 (hDφi i v0).2.2 (hDφi i v0).2.1)
      (designed_lipschitz_aux_integrable hUo _ _ (hui_c i) (hφi_sm i).continuous (hφi_s i).2 (hφi_s i).1)
      (fun z hz => (hui_d i).differentiableAt (hUo.mem_nhds ((hφi_s i).1 hz)))
      (fun z _ => hφi_d i z)
    linarith
  have hB : ∀ i, ∫ z, advect (u z.1) (fun y => u z.1 y i) z.2 * φi i z =
      - ∫ z, ∑ j : Fin 2, ui i z * ui j z * fderiv ℝ (φi i) z (vx j) := by
    intro i
    have hj : ∀ j : Fin 2, ∫ z, ui j z * pd j (fun y => u z.1 y i) z.2 * φi i z =
        - ∫ z, ui i z * (fderiv ℝ (ui j) z (vx j) * φi i z + ui j z * fderiv ℝ (φi i) z (vx j)) := by
      intro j
      set g : ℝ × R2 → ℝ := fun z => ui j z * φi i z with hgdef
      have hg := fun z w => designed_lipschitz_aux_prodDiff hUo (ui j) (φi i) (hui_d j) (hφi_d i) (hφi_s i).1 z w
      have h1 : (fun z : ℝ × R2 => ui j z * pd j (fun y => u z.1 y i) z.2 * φi i z) =
          fun z => fderiv ℝ (ui i) z (vx j) * g z := by
        funext z
        by_cases hz : φi i z = 0
        · simp [hgdef, hz]
        · have hzU := hin i z hz
          have hd : DifferentiableAt ℝ (ui i) z := (hui_d i).differentiableAt (hUo.mem_nhds hzU)
          rw [designed_lipschitz_aux_pdx (fun τ y => u τ y i) z.1 z.2 j hd]
          simp only [hgdef]; ring
      have h2 : (fun z => ui i z * (fderiv ℝ (ui j) z (vx j) * φi i z + ui j z * fderiv ℝ (φi i) z (vx j))) =
          fun z => ui i z * fderiv ℝ g z (vx j) := by
        funext z; rw [(hg z (vx j)).2]
      rw [h1, h2]
      have hgs : tsupport g ⊆ U ∧ HasCompactSupport g := by
        have h := designed_lipschitz_aux_supp (φi i) g (fun z hz => by simp [hgdef, hz])
        exact ⟨h.1.trans (hφi_s i).1, h.2 (hφi_s i).2⟩
      have hI1 : Integrable (fun z => fderiv ℝ (ui i) z (vx j) * g z) := by
        have : (fun z => fderiv ℝ (ui i) z (vx j) * g z) = fun z => (fderiv ℝ (ui i) z (vx j) * ui j z) * φi i z := by
          funext z; simp only [hgdef]; ring
        rw [this]
        exact designed_lipschitz_aux_integrable hUo _ _ ((hDui_c i (vx j)).mul (hui_c j)) (hφi_sm i).continuous
          (hφi_s i).2 (hφi_s i).1
      have hI2 : Integrable (fun z => ui i z * fderiv ℝ g z (vx j)) := by
        have : (fun z => ui i z * fderiv ℝ g z (vx j)) =
            fun z => (ui i z * fderiv ℝ (ui j) z (vx j)) * φi i z + (ui i z * ui j z) * fderiv ℝ (φi i) z (vx j) := by
          funext z; rw [(hg z (vx j)).2]; ring
        rw [this]
        exact (designed_lipschitz_aux_integrable hUo _ _ ((hui_c i).mul (hDui_c j (vx j))) (hφi_sm i).continuous
          (hφi_s i).2 (hφi_s i).1).add
          (designed_lipschitz_aux_integrable hUo _ _ ((hui_c i).mul (hui_c j)) (hDφi i (vx j)).1
          (hDφi i (vx j)).2.2 (hDφi i (vx j)).2.1)
      have hI3 : Integrable (fun z => ui i z * g z) := by
        have : (fun z => ui i z * g z) = fun z => (ui i z * ui j z) * φi i z := by
          funext z; simp only [hgdef]; ring
        rw [this]
        exact designed_lipschitz_aux_integrable hUo _ _ ((hui_c i).mul (hui_c j)) (hφi_sm i).continuous
          (hφi_s i).2 (hφi_s i).1
      have hibp := integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable (μ := volume) (v := vx j)
        (f := ui i) (g := g) hI1 hI2 hI3
        (fun z hz => (hui_d i).differentiableAt (hUo.mem_nhds (hgs.1 hz)))
        (fun z _ => (hg z (vx j)).1)
      linarith
    have hsplit : (fun z : ℝ × R2 => advect (u z.1) (fun y => u z.1 y i) z.2 * φi i z) =
        fun z => ui 0 z * pd 0 (fun y => u z.1 y i) z.2 * φi i z + ui 1 z * pd 1 (fun y => u z.1 y i) z.2 * φi i z := by
      funext z; simp only [advect, huidef]; ring
    have hI : ∀ j : Fin 2, Integrable (fun z => ui j z * pd j (fun y => u z.1 y i) z.2 * φi i z) := by
      intro j
      have h1 : (fun z : ℝ × R2 => ui j z * pd j (fun y => u z.1 y i) z.2 * φi i z) =
          fun z => (fderiv ℝ (ui i) z (vx j) * ui j z) * φi i z := by
        funext z
        by_cases hz : φi i z = 0
        · simp [hz]
        · have hzU := hin i z hz
          have hd : DifferentiableAt ℝ (ui i) z := (hui_d i).differentiableAt (hUo.mem_nhds hzU)
          rw [designed_lipschitz_aux_pdx (fun τ y => u τ y i) z.1 z.2 j hd]; ring
      rw [h1]
      exact designed_lipschitz_aux_integrable hUo _ _ ((hDui_c i (vx j)).mul (hui_c j)) (hφi_sm i).continuous
        (hφi_s i).2 (hφi_s i).1
    rw [hsplit, integral_add (hI 0) (hI 1), hj 0, hj 1]
    have hdivz : ∀ z, ui i z * (fderiv ℝ (ui 0) z (vx 0) * φi i z) + ui i z * (fderiv ℝ (ui 1) z (vx 1) * φi i z) = 0 := by
      intro z
      by_cases hz : φi i z = 0
      · simp [hz]
      · have hzU := hin i z hz
        have hd : ∀ j, DifferentiableAt ℝ (ui j) z := fun j => (hui_d j).differentiableAt (hUo.mem_nhds hzU)
        have h0 := hdiv z.1 (hUS hzU).1 z.2
        simp only [div2] at h0
        rw [designed_lipschitz_aux_pdx (fun τ y => u τ y 0) z.1 z.2 0 (hd 0),
          designed_lipschitz_aux_pdx (fun τ y => u τ y 1) z.1 z.2 1 (hd 1)] at h0
        have : ui i z * (fderiv ℝ (ui 0) z (vx 0) * φi i z) + ui i z * (fderiv ℝ (ui 1) z (vx 1) * φi i z) =
            ui i z * φi i z * (fderiv ℝ (fun z : ℝ × R2 => u z.1 z.2 0) z (vx 0) +
              fderiv ℝ (fun z : ℝ × R2 => u z.1 z.2 1) z (vx 1)) := by
          simp only [huidef]; ring
        rw [this, h0, mul_zero]
    have hIuu : ∀ j : Fin 2, Integrable (fun z => ui i z * (ui j z * fderiv ℝ (φi i) z (vx j))) := by
      intro j
      have : (fun z => ui i z * (ui j z * fderiv ℝ (φi i) z (vx j))) = fun z => (ui i z * ui j z) * fderiv ℝ (φi i) z (vx j) := by
        funext z; ring
      rw [this]
      exact designed_lipschitz_aux_integrable hUo _ _ ((hui_c i).mul (hui_c j)) (hDφi i (vx j)).1
        (hDφi i (vx j)).2.2 (hDφi i (vx j)).2.1
    have hIdu : ∀ j : Fin 2, Integrable (fun z => ui i z * (fderiv ℝ (ui j) z (vx j) * φi i z)) := by
      intro j
      have : (fun z => ui i z * (fderiv ℝ (ui j) z (vx j) * φi i z)) = fun z => (ui i z * fderiv ℝ (ui j) z (vx j)) * φi i z := by
        funext z; ring
      rw [this]
      exact designed_lipschitz_aux_integrable hUo _ _ ((hui_c i).mul (hDui_c j (vx j))) (hφi_sm i).continuous
        (hφi_s i).2 (hφi_s i).1
    have e1 : ∀ j : Fin 2, ∫ z, ui i z * (fderiv ℝ (ui j) z (vx j) * φi i z + ui j z * fderiv ℝ (φi i) z (vx j)) =
        (∫ z, ui i z * (fderiv ℝ (ui j) z (vx j) * φi i z)) + ∫ z, ui i z * (ui j z * fderiv ℝ (φi i) z (vx j)) := by
      intro j
      rw [← integral_add (hIdu j) (hIuu j)]
      congr 1; funext z; ring
    have e2 : (∫ z, ui i z * (fderiv ℝ (ui 0) z (vx 0) * φi i z)) + ∫ z, ui i z * (fderiv ℝ (ui 1) z (vx 1) * φi i z) = 0 := by
      rw [← integral_add (hIdu 0) (hIdu 1)]
      simp only [hdivz, integral_zero]
    have e3 : ∫ z, ∑ j : Fin 2, ui i z * ui j z * fderiv ℝ (φi i) z (vx j) =
        (∫ z, ui i z * (ui 0 z * fderiv ℝ (φi i) z (vx 0))) + ∫ z, ui i z * (ui 1 z * fderiv ℝ (φi i) z (vx 1)) := by
      rw [← integral_add (hIuu 0) (hIuu 1)]
      congr 1; funext z; simp only [Fin.sum_univ_two]; ring
    rw [e1 0, e1 1, e3]
    linarith
  have hC : (∫ z, pd 0 (p z.1) z.2 * φi 0 z) + ∫ z, pd 1 (p z.1) z.2 * φi 1 z = 0 := by
    have hi : ∀ i : Fin 2, ∫ z, pd i (p z.1) z.2 * φi i z = - ∫ z, P z * fderiv ℝ (φi i) z (vx i) := by
      intro i
      have h1 : (fun z : ℝ × R2 => pd i (p z.1) z.2 * φi i z) = fun z => fderiv ℝ P z (vx i) * φi i z := by
        funext z
        by_cases hz : φi i z = 0
        · simp [hz]
        · have hzU := hin i z hz
          have hd : DifferentiableAt ℝ P z := hP_d.differentiableAt (hUo.mem_nhds hzU)
          rw [designed_lipschitz_aux_pdx p z.1 z.2 i hd]
      rw [h1]
      have hibp := integral_mul_fderiv_eq_neg_fderiv_mul_of_integrable (μ := volume) (v := vx i)
        (f := P) (g := φi i)
        (designed_lipschitz_aux_integrable hUo _ _ (hDP_c (vx i)) (hφi_sm i).continuous (hφi_s i).2 (hφi_s i).1)
        (designed_lipschitz_aux_integrable hUo _ _ hP_sm.continuousOn (hDφi i (vx i)).1 (hDφi i (vx i)).2.2
          (hDφi i (vx i)).2.1)
        (designed_lipschitz_aux_integrable hUo _ _ hP_sm.continuousOn (hφi_sm i).continuous (hφi_s i).2 (hφi_s i).1)
        (fun z hz => hP_d.differentiableAt (hUo.mem_nhds ((hφi_s i).1 hz)))
        (fun z _ => hφi_d i z)
      linarith
    have hI : ∀ i : Fin 2, Integrable (fun z => P z * fderiv ℝ (φi i) z (vx i)) := fun i =>
      designed_lipschitz_aux_integrable hUo _ _ hP_sm.continuousOn (hDφi i (vx i)).1 (hDφi i (vx i)).2.2
        (hDφi i (vx i)).2.1
    have hdz : ∀ z, P z * fderiv ℝ (φi 0) z (vx 0) + P z * fderiv ℝ (φi 1) z (vx 1) = 0 := by
      intro z
      have h0 := hφdiv z
      simp only [div2] at h0
      rw [(designed_lipschitz_aux_test (φi 0) (hφi_sm 0) z.1 z.2 0).2,
        (designed_lipschitz_aux_test (φi 1) (hφi_sm 1) z.1 z.2 1).2] at h0
      rw [← mul_add, h0, mul_zero]
    rw [hi 0, hi 1, ← neg_add, ← integral_add (hI 0) (hI 1)]
    simp only [hdz, integral_zero, neg_zero]
  set fi : Fin 2 → ℝ × R2 → ℝ := fun i z => κ * θ z.1 z.2 * (if i = 1 then 1 else 0) + fu z.1 z.2 i with hfidef
  have hfi_c : ∀ i, ContinuousOn (fi i) U := fun i =>
    ((continuousOn_const.mul (hθ.mono hUS)).mul continuousOn_const).add ((hfu i).mono hUS)
  have hpdt_eq : ∀ i, (fun z : ℝ × R2 => pdt (fun τ y => u τ y i) z.1 z.2 * φi i z) =
      fun z => fderiv ℝ (ui i) z v0 * φi i z := by
    intro i; funext z
    by_cases hz : φi i z = 0
    · simp [hz]
    · have hzU := hin i z hz
      have hd : DifferentiableAt ℝ (ui i) z := (hui_d i).differentiableAt (hUo.mem_nhds hzU)
      rw [designed_lipschitz_aux_pdt (fun τ y => u τ y i) z.1 z.2 hzU.1.1 hd]
  have hadv_eq : ∀ i, (fun z : ℝ × R2 => advect (u z.1) (fun y => u z.1 y i) z.2 * φi i z) =
      fun z => (fderiv ℝ (ui i) z (vx 0) * ui 0 z + fderiv ℝ (ui i) z (vx 1) * ui 1 z) * φi i z := by
    intro i; funext z
    by_cases hz : φi i z = 0
    · simp [hz]
    · have hzU := hin i z hz
      have hd : DifferentiableAt ℝ (ui i) z := (hui_d i).differentiableAt (hUo.mem_nhds hzU)
      simp only [advect, huidef]
      rw [designed_lipschitz_aux_pdx (fun τ y => u τ y i) z.1 z.2 0 hd,
        designed_lipschitz_aux_pdx (fun τ y => u τ y i) z.1 z.2 1 hd]
      ring
  have hp_eq : ∀ i, (fun z : ℝ × R2 => pd i (p z.1) z.2 * φi i z) = fun z => fderiv ℝ P z (vx i) * φi i z := by
    intro i; funext z
    by_cases hz : φi i z = 0
    · simp [hz]
    · have hzU := hin i z hz
      have hd : DifferentiableAt ℝ P z := hP_d.differentiableAt (hUo.mem_nhds hzU)
      rw [designed_lipschitz_aux_pdx p z.1 z.2 i hd]
  have hIpdt : ∀ i, Integrable (fun z : ℝ × R2 => pdt (fun τ y => u τ y i) z.1 z.2 * φi i z) := fun i => by
    rw [hpdt_eq i]
    exact designed_lipschitz_aux_integrable hUo _ _ (hDui_c i v0) (hφi_sm i).continuous (hφi_s i).2 (hφi_s i).1
  have hIadv : ∀ i, Integrable (fun z : ℝ × R2 => advect (u z.1) (fun y => u z.1 y i) z.2 * φi i z) := fun i => by
    rw [hadv_eq i]
    exact designed_lipschitz_aux_integrable hUo _ _ (((hDui_c i (vx 0)).mul (hui_c 0)).add
      ((hDui_c i (vx 1)).mul (hui_c 1))) (hφi_sm i).continuous (hφi_s i).2 (hφi_s i).1
  have hIp : ∀ i, Integrable (fun z : ℝ × R2 => pd i (p z.1) z.2 * φi i z) := fun i => by
    rw [hp_eq i]
    exact designed_lipschitz_aux_integrable hUo _ _ (hDP_c (vx i)) (hφi_sm i).continuous (hφi_s i).2 (hφi_s i).1
  have hIf : ∀ i, Integrable (fun z => fi i z * φi i z) := fun i =>
    designed_lipschitz_aux_integrable hUo _ _ (hfi_c i) (hφi_sm i).continuous (hφi_s i).2 (hφi_s i).1
  have hIt : ∀ i, Integrable (fun z => ui i z * fderiv ℝ (φi i) z v0) := fun i =>
    designed_lipschitz_aux_integrable hUo _ _ (hui_c i) (hDφi i v0).1 (hDφi i v0).2.2 (hDφi i v0).2.1
  have hIuu : ∀ i j, Integrable (fun z => ui i z * ui j z * fderiv ℝ (φi i) z (vx j)) := fun i j =>
    designed_lipschitz_aux_integrable hUo _ _ ((hui_c i).mul (hui_c j)) (hDφi i (vx j)).1 (hDφi i (vx j)).2.2
      (hDφi i (vx j)).2.1
  have hsumI : ∀ i, (∫ z, pdt (fun τ y => u τ y i) z.1 z.2 * φi i z) +
      (∫ z, advect (u z.1) (fun y => u z.1 y i) z.2 * φi i z) + (∫ z, pd i (p z.1) z.2 * φi i z) =
      ∫ z, fi i z * φi i z := by
    intro i
    have h12 : Integrable (fun z : ℝ × R2 => pdt (fun τ y => u τ y i) z.1 z.2 * φi i z +
        advect (u z.1) (fun y => u z.1 y i) z.2 * φi i z) := (hIpdt i).add (hIadv i)
    have hpt : (fun z : ℝ × R2 => pdt (fun τ y => u τ y i) z.1 z.2 * φi i z +
        advect (u z.1) (fun y => u z.1 y i) z.2 * φi i z + pd i (p z.1) z.2 * φi i z) = fun z => fi i z * φi i z := by
      funext z
      by_cases hz : φi i z = 0
      · simp [hz]
      · have hzU := hin i z hz
        have heq := hequ z.1 (hUS hzU).1 z.2 i
        simp only [hfidef]
        rw [← heq]; ring
    rw [← hpt, integral_add h12 (hIp i), integral_add (hIpdt i) (hIadv i)]
  have hF : ∀ i, ∫ z, (ui i z * fderiv ℝ (φi i) z v0 + ∑ j : Fin 2, ui i z * ui j z * fderiv ℝ (φi i) z (vx j) +
      fi i z * φi i z) = ∫ z, pd i (p z.1) z.2 * φi i z := by
    intro i
    have hIs : Integrable (fun z => ∑ j : Fin 2, ui i z * ui j z * fderiv ℝ (φi i) z (vx j)) := by
      simp only [Fin.sum_univ_two]; exact (hIuu i 0).add (hIuu i 1)
    have h12 : Integrable (fun z => ui i z * fderiv ℝ (φi i) z v0 +
        ∑ j : Fin 2, ui i z * ui j z * fderiv ℝ (φi i) z (vx j)) := (hIt i).add hIs
    rw [integral_add h12 (hIf i), integral_add (hIt i) hIs]
    linarith [hA i, hB i, hsumI i]
  have ht1 : ∀ i (z : ℝ × R2), deriv (fun τ => (φ (τ, z.2)) i) z.1 = fderiv ℝ (φi i) z v0 := fun i z =>
    (designed_lipschitz_aux_test (φi i) (hφi_sm i) z.1 z.2 0).1
  have ht2 : ∀ i j (z : ℝ × R2), pd j (fun y => (φ (z.1, y)) i) z.2 = fderiv ℝ (φi i) z (vx j) := fun i j z =>
    (designed_lipschitz_aux_test (φi i) (hφi_sm i) z.1 z.2 j).2
  have hint_eq : (fun z : ℝ × R2 => ∑ i : Fin 2,
      ( (u z.1 z.2) i * deriv (fun τ => (φ (τ, z.2)) i) z.1
        + ∑ j : Fin 2, (u z.1 z.2) i * (u z.1 z.2) j * pd j (fun y => (φ (z.1, y)) i) z.2
        + (κ * θ z.1 z.2 * (if i = 1 then 1 else 0) + (fu z.1 z.2) i) * (φ z) i )) =
      fun z => ∑ i : Fin 2, (ui i z * fderiv ℝ (φi i) z v0 +
        ∑ j : Fin 2, ui i z * ui j z * fderiv ℝ (φi i) z (vx j) + fi i z * φi i z) := by
    funext z
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [ht1, ht2, huidef, hφidef, hfidef]
  have hFint : ∀ i, Integrable (fun z => ui i z * fderiv ℝ (φi i) z v0 +
      ∑ j : Fin 2, ui i z * ui j z * fderiv ℝ (φi i) z (vx j) + fi i z * φi i z) := by
    intro i
    have hIs : Integrable (fun z => ∑ j : Fin 2, ui i z * ui j z * fderiv ℝ (φi i) z (vx j)) := by
      simp only [Fin.sum_univ_two]; exact (hIuu i 0).add (hIuu i 1)
    exact ((hIt i).add hIs).add (hIf i)
  rw [hint_eq, integral_finsetSum _ (fun i _ => hFint i)]
  rw [Fin.sum_univ_two, hF 0, hF 1]
  exact hC

theorem designed_lipschitz_aux_gsum {κ : ℝ} {Pr : Profiles} {F' : Type*} [AddCommMonoid F'] [TopologicalSpace F']
    (NStar Cg : ℝ → ℝ) (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) (k : ℕ) (τ : ℝ)
    (hτ : τ < T.t (k + 1)) (f : ℕ → F') :
    (∑' m, if T.tIn m ≤ τ then f m else 0) = ∑ m ∈ Finset.range (k + 2), if T.tIn m ≤ τ then f m else 0 := by
  have hmono := designed_classical_aux_mono NStar Cg T hT
  refine tsum_eq_sum (s := Finset.range (k + 2)) (fun m hm => ?_)
  rw [if_neg]
  unfold TowerData.tIn
  exact not_le.2 (lt_of_lt_of_le hτ (hmono (by have := Finset.mem_range.not.1 hm; omega)))

theorem designed_lipschitz_aux_layerU {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ)
    (hNuni : ∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) (hP : PotentialTheoryR2)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) (k : ℕ) (m : ℕ) (hm1 : 1 ≤ m) (hmk : m ≤ k + 1) :
    ∃ ρ : ℝ, 0 < ρ ∧ ∃ G : ℝ → R2 → ℝ,
      ContDiffOn ℝ ∞ (fun z : ℝ × R2 => G z.1 z.2) (Icc 0 (T.t (k + 2)) ×ˢ univ) ∧
      (∀ i : Fin 2, ContDiffOn ℝ ∞ (fun z : ℝ × R2 => if T.tIn m ≤ z.1 then T.uL m z.1 z.2 i else 0)
        (Icc 0 (T.t (k + 2)) ×ˢ univ)) ∧
      ContDiffOn ℝ ∞ (fun z : ℝ × R2 => if T.tIn m ≤ z.1 then T.θL m z.1 z.2 else 0) (Icc 0 (T.t (k + 2)) ×ˢ univ) ∧
      (∀ τ ∈ Icc 0 (T.t (k + 2)), ∀ x : R2, ρ ≤ ‖x‖ →
        G τ x = 0 ∧ (if T.tIn m ≤ τ then T.θL m τ x else 0) = 0) ∧
      (∀ τ ∈ Icc 0 (T.t (k + 2)), ContDiff ℝ ∞ (G τ) ∧ tsupport (G τ) ⊆ Metric.ball 0 ρ ∧ ∫ y, G τ y = 0 ∧
        (fun x => if T.tIn m ≤ τ then T.uL m τ x else 0) = biotSavart (G τ)) := by
  have hmono := designed_classical_aux_mono NStar Cg T hT
  have hm0 : m ≠ 0 := by omega
  have ht' : T.t (k + 2) ∈ T.I (k + 2) := ⟨le_rfl, hmono (by omega)⟩
  obtain ⟨-, -, ⟨R, hKin, hAmp, K, hS4, hG, -⟩, -⟩ := hT.older (k + 2) (by omega) (T.t (k + 2)) ht' m hm1 (by omega)
  obtain ⟨tv, hSO⟩ := hT.sizes (k + 2) (by omega) (T.t (k + 2)) ht' m hm1 (by omega) R hKin
  obtain ⟨hDsel, hrun⟩ := hT.run (k + 1) (by omega)
  have hN : T.N m = (T.P m).N := by unfold TowerData.N; rw [if_neg hm0]
  have hthr : S4.NthrUniform Pr (T.P m).β (Cg (T.P m).β) ≤ (T.P m).N := by
    rw [← hN]; exact (hNuni _ (T.P m).β_pos).trans (hT.thresholds m hm1)
  have hcutall := S4.nthrUniform_spec Pr (T.P m).β (Cg (T.P m).β) (T.P m) rfl hthr R K (T.N (m - 1)) hG hS4 hAmp
  have htin_eq : (T.bg m (T.t (k + 2))).tIn = T.tIn m := rfl
  have hτ_eq : (T.bg m (T.t (k + 2))).τ = T.t (k + 2) := rfl
  obtain ⟨Rχ, hRχ⟩ : ∃ r, tsupport Pr.χc ⊆ Metric.closedBall (0 : R2) r :=
    (Pr.hχc.compact : IsCompact (tsupport Pr.χc)).isBounded.subset_closedBall 0
  set ρ := max (max ((T.P m).gammaStar * (T.P m).ell) Rχ) 0 + 1 with hρ
  have hρpos : 0 < ρ := by
    rw [hρ]; linarith [le_max_right (max ((T.P m).gammaStar * (T.P m).ell) Rχ) 0]
  have hρ1 : (T.P m).gammaStar * (T.P m).ell < ρ := by
    rw [hρ]
    linarith [le_max_left ((T.P m).gammaStar * (T.P m).ell) Rχ, le_max_left (max ((T.P m).gammaStar * (T.P m).ell) Rχ) 0]
  have hρ2 : Rχ < ρ := by
    rw [hρ]
    linarith [le_max_right ((T.P m).gammaStar * (T.P m).ell) Rχ, le_max_left (max ((T.P m).gammaStar * (T.P m).ell) Rχ) 0]
  set G : ℝ → R2 → ℝ := fun τ x => if T.tIn m ≤ τ then S4.omegaLayer Pr (T.P m) R (T.N (m - 1)) τ x else 0 with hGdef
  have hsub : Icc 0 (T.t (k + 2)) ×ˢ (univ : Set R2) ⊆ Iic (T.bg m (T.t (k + 2))).τ ×ˢ univ :=
    prod_mono (fun τ hτ => hτ.2) subset_rfl
  have hreal : ∀ τ, T.tIn m ≤ τ → τ ≤ T.t (k + 2) → ∀ y,
      T.θL m τ y = S4.thetaLayer Pr (T.P m) R (T.N (m - 1)) τ y ∧
      T.omL m τ y = S4.omegaLayer Pr (T.P m) R (T.N (m - 1)) τ y ∧
      T.uL m τ y = biotSavart (S4.omegaLayer Pr (T.P m) R (T.N (m - 1)) τ) y :=
    fun τ h1 h2 y => hrun.realised m hm1 (by omega) R hKin τ ⟨h1, h2⟩ y
  refine ⟨ρ, hρpos, G, ?_, ?_, ?_, ?_, ?_⟩
  · exact hcutall.2.2.1.mono hsub
  · intro i
    refine ((hcutall.2.2.2 i).mono hsub).congr fun z hz => ?_
    by_cases h : T.tIn m ≤ z.1
    · simp only [htin_eq, if_pos h]
      rw [(hreal z.1 h hz.1.2 z.2).2.2]; rfl
    · simp only [htin_eq, if_neg h]
  · refine (hcutall.2.1.mono hsub).congr fun z hz => ?_
    by_cases h : T.tIn m ≤ z.1
    · simp only [htin_eq, if_pos h]; rw [(hreal z.1 h hz.1.2 z.2).1]
    · simp only [htin_eq, if_neg h]
  · intro τ hτ x hx
    by_cases h : T.tIn m ≤ τ
    · have hera : τ ∈ Icc (T.bg m (T.t (k + 2))).tIn (T.bg m (T.t (k + 2))).τ := ⟨h, hτ.2⟩
      constructor
      · simp only [hGdef, if_pos h]
        apply image_eq_zero_of_notMem_tsupport
        intro hmem
        rcases hSO.supp_om τ hera hmem with h' | h'
        · have := Metric.mem_closedBall.1 h'; rw [dist_zero_right] at this; linarith
        · have := Metric.mem_closedBall.1 (hRχ h'); rw [dist_zero_right] at this; linarith
      · simp only [if_pos h]
        rw [(hreal τ h hτ.2 x).1]
        apply image_eq_zero_of_notMem_tsupport
        intro hmem
        have := Metric.mem_closedBall.1 (hSO.supp_th τ hera hmem); rw [dist_zero_right] at this; linarith
    · simp [hGdef, if_neg h]
  · intro τ hτ
    by_cases h : T.tIn m ≤ τ
    · have hera : τ ∈ Icc (T.bg m (T.t (k + 2))).tIn (T.bg m (T.t (k + 2))).τ := ⟨h, hτ.2⟩
      have hGτ : G τ = S4.omegaLayer Pr (T.P m) R (T.N (m - 1)) τ := by
        funext x; simp only [hGdef, if_pos h]
      have hsl : ContDiff ℝ ∞ (G τ) :=
        Kit.paramLinearOp_contDiffOn_aux_slice 0 (T.t (k + 2)) G (hcutall.2.2.1.mono hsub) τ hτ
      rw [hGτ] at hsl
      rw [hGτ]
      refine ⟨hsl, ?_, hSO.meanZero τ hera, ?_⟩
      · refine (hSO.supp_om τ hera).trans ?_
        intro x hx'
        rw [Metric.mem_ball, dist_zero_right]
        rcases hx' with h' | h'
        · have := Metric.mem_closedBall.1 h'; rw [dist_zero_right] at this; linarith
        · have := Metric.mem_closedBall.1 (hRχ h'); rw [dist_zero_right] at this; linarith
      · funext x; simp only [if_pos h]; exact (hreal τ h hτ.2 x).2.2
    · have hGτ : G τ = fun _ => 0 := by funext x; simp only [hGdef, if_neg h]
      rw [hGτ]
      refine ⟨contDiff_const, by simp [tsupport], by simp, ?_⟩
      rw [lemma48_g_aux_biotSavart_zero]
      funext x; simp only [if_neg h]

theorem designed_lipschitz_aux_opn1 (g : R2 → ℝ) (x : R2) (M : ℝ) (hM : 0 ≤ M)
    (h : ∀ j : Fin 2, |pd j g x| ≤ M) : ‖fderiv ℝ g x‖ ≤ 2 * M := by
  refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun v => ?_
  rw [lemma42_size_ak_clm_aux (fderiv ℝ g x) v, Real.norm_eq_abs]
  have hv : ∀ j : Fin 2, |v j| ≤ ‖v‖ := fun j => by
    rw [EuclideanSpace.norm_eq, Fin.sum_univ_two, Real.norm_eq_abs, Real.norm_eq_abs]
    rcases j with ⟨j, hj⟩
    interval_cases j
    · calc |v 0| = Real.sqrt (|v 0| ^ 2) := (Real.sqrt_sq (abs_nonneg _)).symm
        _ ≤ Real.sqrt (|v 0| ^ 2 + |v 1| ^ 2) := Real.sqrt_le_sqrt (by nlinarith [sq_nonneg (v 1)])
    · calc |v 1| = Real.sqrt (|v 1| ^ 2) := (Real.sqrt_sq (abs_nonneg _)).symm
        _ ≤ Real.sqrt (|v 0| ^ 2 + |v 1| ^ 2) := Real.sqrt_le_sqrt (by nlinarith [sq_nonneg (v 0)])
  have h0 := h 0; have h1 := h 1
  unfold pd at h0 h1
  calc |v 0 * fderiv ℝ g x (e 0) + v 1 * fderiv ℝ g x (e 1)|
      ≤ |v 0| * |fderiv ℝ g x (e 0)| + |v 1| * |fderiv ℝ g x (e 1)| := by
        refine (abs_add_le _ _).trans ?_; rw [abs_mul, abs_mul]
    _ ≤ ‖v‖ * M + ‖v‖ * M := by
        gcongr
        · exact hv 0
        · exact hv 1
    _ = 2 * M * ‖v‖ := by ring

theorem designed_lipschitz_aux_opn2 (F : R2 → R2) (x : R2) (hF : DifferentiableAt ℝ F x) (M : ℝ) (hM : 0 ≤ M)
    (h : ∀ i j : Fin 2, |pd j (fun y => F y i) x| ≤ M) : ‖fderiv ℝ F x‖ ≤ 4 * M := by
  have hcomp : ∀ i : Fin 2, ∀ v, (fderiv ℝ F x v) i = fderiv ℝ (fun y => F y i) x v := by
    intro i v
    have h1 : HasFDerivAt (⇑(EuclideanSpace.proj (𝕜 := ℝ) (ι := Fin 2) i) ∘ F)
        ((EuclideanSpace.proj (𝕜 := ℝ) (ι := Fin 2) i).comp (fderiv ℝ F x)) x :=
      (EuclideanSpace.proj (𝕜 := ℝ) (ι := Fin 2) i).hasFDerivAt.comp x hF.hasFDerivAt
    have h2 : (fun y => F y i) = ⇑(EuclideanSpace.proj (𝕜 := ℝ) (ι := Fin 2) i) ∘ F := rfl
    rw [h2, h1.fderiv]; rfl
  have hci : ∀ i : Fin 2, ‖fderiv ℝ (fun y => F y i) x‖ ≤ 2 * M := fun i =>
    designed_lipschitz_aux_opn1 (fun y => F y i) x M hM (h i)
  refine ContinuousLinearMap.opNorm_le_bound _ (by positivity) fun v => ?_
  have hn : ‖fderiv ℝ F x v‖ ≤ |(fderiv ℝ F x v) 0| + |(fderiv ℝ F x v) 1| := by
    rw [EuclideanSpace.norm_eq, Fin.sum_univ_two, Real.norm_eq_abs, Real.norm_eq_abs]
    rw [Real.sqrt_le_left (by positivity)]
    nlinarith [abs_nonneg ((fderiv ℝ F x v) 0), abs_nonneg ((fderiv ℝ F x v) 1)]
  calc ‖fderiv ℝ F x v‖ ≤ |(fderiv ℝ F x v) 0| + |(fderiv ℝ F x v) 1| := hn
    _ = |fderiv ℝ (fun y => F y 0) x v| + |fderiv ℝ (fun y => F y 1) x v| := by rw [hcomp 0, hcomp 1]
    _ ≤ 2 * M * ‖v‖ + 2 * M * ‖v‖ := by
        gcongr
        · rw [← Real.norm_eq_abs]; exact (fderiv ℝ (fun y => F y 0) x).le_of_opNorm_le (hci 0) v
        · rw [← Real.norm_eq_abs]; exact (fderiv ℝ (fun y => F y 1) x).le_of_opNorm_le (hci 1) v
    _ = 4 * M * ‖v‖ := by ring

theorem designed_lipschitz_aux_L2 {E' : Type*} [NormedAddCommGroup E'] (g : R2 → E') (M ρ : ℝ) (hM : 0 ≤ M)
    (hb : ∀ x, ‖g x‖ ≤ M) (hρ : ∀ x : R2, ρ ≤ ‖x‖ → g x = 0) :
    eLpNorm g 2 volume ≤ ENNReal.ofReal (M * (volume (Metric.closedBall (0 : R2) ρ)).toReal ^ (1 / (2 : ℝ))) := by
  have hpt : ∀ x, ‖g x‖ ≤ ‖(Metric.closedBall (0 : R2) ρ).indicator (fun _ => M) x‖ := by
    intro x
    by_cases hx : x ∈ Metric.closedBall (0 : R2) ρ
    · rw [indicator_of_mem hx, Real.norm_eq_abs, abs_of_nonneg hM]; exact hb x
    · rw [Metric.mem_closedBall, dist_zero_right, not_le] at hx
      rw [hρ x hx.le, norm_zero]; exact norm_nonneg _
  refine (eLpNorm_mono hpt).trans ?_
  rw [eLpNorm_indicator_const measurableSet_closedBall (by norm_num) (by norm_num)]
  have hv : volume (Metric.closedBall (0 : R2) ρ) = ENNReal.ofReal (volume (Metric.closedBall (0 : R2) ρ)).toReal :=
    (ENNReal.ofReal_toReal (measure_closedBall_lt_top).ne).symm
  rw [hv, ENNReal.toReal_ofReal ENNReal.toReal_nonneg, ENNReal.toReal_ofNat,
    ENNReal.ofReal_rpow_of_nonneg ENNReal.toReal_nonneg (by norm_num : (0 : ℝ) ≤ 1 / 2), Real.enorm_eq_ofReal_abs,
    abs_of_nonneg hM, ← ENNReal.ofReal_mul hM]

theorem designed_lipschitz_aux_n2 (y : R2) : ‖y‖ ≤ |y 0| + |y 1| := by
  rw [EuclideanSpace.norm_eq, Fin.sum_univ_two, Real.norm_eq_abs, Real.norm_eq_abs]
  rw [Real.sqrt_le_left (by positivity)]
  nlinarith [abs_nonneg (y 0), abs_nonneg (y 1)]

theorem designed_lipschitz_aux_base {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ) (T : TowerData κ Pr)
    (hT : TowerSpec NStar Cg T) :
    ∃ Rb : ℝ, ∀ τ, 0 ≤ τ → ∀ x : R2, Rb ≤ ‖x‖ → T.uL 0 τ x = 0 ∧ T.θL 0 τ x = 0 ∧ T.uSt τ x = 0 := by
  obtain ⟨Rχ0, hRχ0⟩ : ∃ r, tsupport Pr.χ0 ⊆ Metric.closedBall (0 : R2) r :=
    (Pr.hχ.compact : IsCompact (tsupport Pr.χ0)).isBounded.subset_closedBall 0
  obtain ⟨RH, hRH⟩ : ∃ r, tsupport Pr.H ⊆ Metric.closedBall (0 : R2) r :=
    (Pr.hH.compact : IsCompact (tsupport Pr.H)).isBounded.subset_closedBall 0
  refine ⟨max Rχ0 RH + 1, fun τ hτ x hx => ?_⟩
  have hxχ : Rχ0 < ‖x‖ := by linarith [le_max_left Rχ0 RH]
  have hxH : RH < ‖x‖ := by linarith [le_max_right Rχ0 RH]
  have hχ0 : ∀ y : R2, Rχ0 < ‖y‖ → Pr.χ0 (rot (-(T.ψ τ)) y) = 0 := by
    intro y hy
    apply image_eq_zero_of_notMem_tsupport
    intro hmem
    have := Metric.mem_closedBall.1 (hRχ0 hmem)
    rw [dist_zero_right, S7.parked_node_vorticity_aux_rot_norm] at this
    linarith
  have hψsupp : x ∉ tsupport (T.basePsi τ) := by
    intro hmem
    have hsub : tsupport (T.basePsi τ) ⊆ Metric.closedBall (0 : R2) Rχ0 := by
      refine closure_minimal (fun y hy => ?_) Metric.isClosed_closedBall
      rw [Metric.mem_closedBall, dist_zero_right]
      by_contra h
      apply hy
      show -(T.Om0 τ / T.base.N0 ^ 2) * Real.cos (inner ℝ (T.xi 0 τ) y) * Pr.χ0 (rot (-(T.ψ τ)) y) = 0
      rw [hχ0 y (not_le.1 h), mul_zero]
    have := Metric.mem_closedBall.1 (hsub hmem)
    rw [dist_zero_right] at this; linarith
  have hHsupp : x ∉ tsupport Pr.H := by
    intro hmem
    have := Metric.mem_closedBall.1 (hRH hmem)
    rw [dist_zero_right] at this; linarith
  have hpg : ∀ (F : R2 → ℝ), x ∉ tsupport F → perpGrad F x = 0 := by
    intro F hF
    have h0 : ∀ i : Fin 2, pd i F x = 0 := by
      intro i
      have : fderiv ℝ F x = 0 := fderiv_of_notMem_tsupport ℝ hF
      simp [pd, this]
    simp [perpGrad, h0]
  refine ⟨?_, ?_, ?_⟩
  · rw [hT.base.u0 τ hτ x]; exact hpg _ hψsupp
  · rw [hT.base.theta0 τ hτ x, hχ0 x hxχ, mul_zero]
  · show deriv T.ψ τ • perpGrad Pr.H x = 0
    rw [hpg _ hHsupp, smul_zero]

theorem designed_lipschitz_aux_bs (hP : PotentialTheoryR2) (ζ : R2 → ℝ) (hζ : ContDiff ℝ ∞ ζ)
    (hc : HasCompactSupport ζ) : ContDiff ℝ ∞ (biotSavart ζ) := by
  obtain ⟨hNPsm, -, -⟩ := hP.newtonian ζ hζ hc
  have hpd : ∀ i : Fin 2, ContDiff ℝ ∞ (pd i (newtonPotential ζ)) := fun i =>
    (hNPsm.fderiv_right (m := ∞) (by norm_cast)).clm_apply contDiff_const
  show ContDiff ℝ ∞ (perpGrad (newtonPotential ζ))
  unfold perpGrad
  exact ((hpd 1).neg.smul contDiff_const).add ((hpd 0).smul contDiff_const)

theorem designed_lipschitz {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ) (Kexp : ℕ)
    (hNdom : ∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg Kexp) Cg β ≤ NStar β)
    (hNuni : ∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) (hκ : 0 < κ) (hP : PotentialTheoryR2)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T) (fu : ℝ → R2 → R2) (p : ℝ → R2 → ℝ)
    (hsol : ClassicalBoussinesq κ T.Tstar (S7.Designed.theta T) p (S7.Designed.vel T) fu (S7.Designed.fTheta T))
    (hfu_sm : ∀ i : Fin 2, ContDiffOn ℝ ∞ (fun z : ℝ × R2 => fu z.1 z.2 i) (Ico 0 T.Tstar ×ˢ univ))
    (hfu_supp : ∃ Rball : ℝ, ∀ t ∈ Ico 0 T.Tstar, ∀ x : R2, Rball < ‖x‖ → fu t x = 0) :
    ∀ T' < T.Tstar, 0 < T' →
      InLipschitzClass κ T' (S7.Designed.theta T) (S7.Designed.vel T) (S7.Designed.fTheta T) fu := by
  intro T' hT'lt hT'pos
  classical
  haveI hHaar : (volume : Measure (ℝ × R2)).IsAddHaarMeasure := Measure.prod.instIsAddHaarMeasure _ _
  have hmono := designed_classical_aux_mono NStar Cg T hT
  obtain ⟨k, hk⟩ := designed_classical_aux_stage NStar Cg T hT T' ⟨hT'pos.le, hT'lt⟩
  have hkt : T' < T.t (k + 1) := hk.2
  have hk2 : T.t (k + 1) ≤ T.t (k + 2) := hmono (by omega)
  have ht0 : T.t 0 = 0 := (hT.state 1 le_rfl).1.t_zero
  have hIco : Icc 0 T' ⊆ Ico 0 T.Tstar := Icc_subset_Ico_right hT'lt
  have hS : Icc 0 T' ×ˢ (univ : Set R2) ⊆ Ico 0 T.Tstar ×ˢ univ := prod_mono hIco subset_rfl
  have hS2 : Icc 0 T' ×ˢ (univ : Set R2) ⊆ Icc 0 (T.t (k + 2)) ×ˢ univ :=
    prod_mono (Icc_subset_Icc_right (hkt.le.trans hk2)) subset_rfl
  have hτH : ∀ τ ∈ Icc 0 T', τ ∈ Icc 0 (T.t (k + 2)) := fun τ hτ => ⟨hτ.1, hτ.2.trans (hkt.le.trans hk2)⟩
  have hLay0 : ∀ m : ℕ, ∃ ρ : ℝ, ∃ G : ℝ → R2 → ℝ, (1 ≤ m ∧ m ≤ k + 1) → (0 < ρ ∧
      ContDiffOn ℝ ∞ (fun z : ℝ × R2 => G z.1 z.2) (Icc 0 (T.t (k + 2)) ×ˢ univ) ∧
      (∀ i : Fin 2, ContDiffOn ℝ ∞ (fun z : ℝ × R2 => if T.tIn m ≤ z.1 then T.uL m z.1 z.2 i else 0)
        (Icc 0 (T.t (k + 2)) ×ˢ univ)) ∧
      ContDiffOn ℝ ∞ (fun z : ℝ × R2 => if T.tIn m ≤ z.1 then T.θL m z.1 z.2 else 0) (Icc 0 (T.t (k + 2)) ×ˢ univ) ∧
      (∀ τ ∈ Icc 0 (T.t (k + 2)), ∀ x : R2, ρ ≤ ‖x‖ →
        G τ x = 0 ∧ (if T.tIn m ≤ τ then T.θL m τ x else 0) = 0) ∧
      (∀ τ ∈ Icc 0 (T.t (k + 2)), ContDiff ℝ ∞ (G τ) ∧ tsupport (G τ) ⊆ Metric.ball 0 ρ ∧ ∫ y, G τ y = 0 ∧
        (fun x => if T.tIn m ≤ τ then T.uL m τ x else 0) = biotSavart (G τ))) := by
    intro m
    by_cases hm : 1 ≤ m ∧ m ≤ k + 1
    · obtain ⟨ρ, hρ, G, h⟩ := designed_lipschitz_aux_layerU NStar Cg hNuni hP T hT k m hm.1 hm.2
      exact ⟨ρ, G, fun _ => ⟨hρ, h⟩⟩
    · exact ⟨0, fun _ _ => 0, fun h => absurd h hm⟩
  choose ρ G hLay using hLay0
  have hL := fun m (hm : m ∈ Finset.range (k + 1)) =>
    hLay (m + 1) ⟨by omega, by have := Finset.mem_range.1 hm; omega⟩
  obtain ⟨Rb, hRb⟩ := designed_lipschitz_aux_base NStar Cg T hT
  set Rall : ℝ := (∑ m ∈ Finset.range (k + 1), ρ (m + 1)) + |Rb| + 1 with hRall
  have hρle : ∀ m ∈ Finset.range (k + 1), ρ (m + 1) ≤ Rall := by
    intro m hm
    have h1 : ρ (m + 1) ≤ ∑ m ∈ Finset.range (k + 1), ρ (m + 1) :=
      Finset.single_le_sum (f := fun m => ρ (m + 1)) (fun m' hm' => (hL m' hm').1.le) hm
    rw [hRall]; linarith [abs_nonneg Rb]
  have hRb_le : Rb ≤ Rall := by
    rw [hRall]
    have : 0 ≤ ∑ m ∈ Finset.range (k + 1), ρ (m + 1) := Finset.sum_nonneg fun m' hm' => (hL m' hm').1.le
    linarith [le_abs_self Rb]
  set θd := S7.Designed.theta T with hθd
  set ud := S7.Designed.vel T with hud
  set gv : ℕ → ℝ → R2 → R2 := fun m τ x => if T.tIn m ≤ τ then T.uL m τ x else 0 with hgv
  set w : ℝ → R2 → R2 := fun τ x => T.uSt τ x + T.uL 0 τ x with hw
  have hguard0 : ∀ τ, 0 ≤ τ → T.tIn 0 ≤ τ := fun τ hτ => by unfold TowerData.tIn; simp [ht0, hτ]
  have hθsum : ∀ τ ∈ Icc 0 T', ∀ x, θd τ x = T.θL 0 τ x +
      ∑ m ∈ Finset.range (k + 1), (if T.tIn (m + 1) ≤ τ then T.θL (m + 1) τ x else 0) := by
    intro τ hτ x
    simp only [hθd, S7.Designed.theta]
    rw [designed_lipschitz_aux_gsum NStar Cg T hT k τ (lt_of_le_of_lt hτ.2 hkt), Finset.sum_range_succ',
      if_pos (hguard0 τ hτ.1), add_comm]
  have husum : ∀ τ ∈ Icc 0 T', ∀ x, ud τ x = w τ x + ∑ m ∈ Finset.range (k + 1), gv (m + 1) τ x := by
    intro τ hτ x
    simp only [hud, S7.Designed.vel, hw, hgv]
    rw [designed_lipschitz_aux_gsum NStar Cg T hT k τ (lt_of_le_of_lt hτ.2 hkt), Finset.sum_range_succ',
      if_pos (hguard0 τ hτ.1)]
    abel
  have hθsm : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => θd z.1 z.2) (Icc 0 T' ×ˢ univ) := hsol.smoothθ.mono hS
  have husm : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => ud z.1 z.2) (Icc 0 T' ×ˢ univ) := hsol.smoothu.mono hS
  have hwsm : ∀ i : Fin 2, ContDiffOn ℝ ∞ (fun z : ℝ × R2 => (w z.1 z.2) i) (Icc 0 T' ×ˢ univ) := by
    intro i
    have h1 : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => (ud z.1 z.2) i -
        ∑ m ∈ Finset.range (k + 1), (if T.tIn (m + 1) ≤ z.1 then T.uL (m + 1) z.1 z.2 i else 0)) (Icc 0 T' ×ˢ univ) :=
      (contDiffOn_euclidean.1 husm i).sub (ContDiffOn.sum fun m hm => ((hL m hm).2.2.1 i).mono hS2)
    refine h1.congr fun z hz => ?_
    have := husum z.1 hz.1 z.2
    have hcomp : (ud z.1 z.2) i = (w z.1 z.2) i + ∑ m ∈ Finset.range (k + 1), (gv (m + 1) z.1 z.2) i := by
      rw [this]; simp
    have hgi : ∀ m, (gv (m + 1) z.1 z.2) i = (if T.tIn (m + 1) ≤ z.1 then T.uL (m + 1) z.1 z.2 i else 0) := by
      intro m; simp only [hgv]; split_ifs <;> simp
    rw [hcomp]; simp only [hgi]; ring
  have hθvan : ∀ τ ∈ Icc 0 T', ∀ x : R2, Rall ≤ ‖x‖ → θd τ x = 0 := by
    intro τ hτ x hx
    rw [hθsum τ hτ x, (hRb τ hτ.1 x (hRb_le.trans hx)).2.1, zero_add]
    refine Finset.sum_eq_zero fun m hm => ?_
    exact ((hL m hm).2.2.2.2.1 τ (hτH τ hτ) x ((hρle m hm).trans hx)).2
  have hwvan : ∀ i : Fin 2, ∀ τ ∈ Icc 0 T', ∀ x : R2, Rall ≤ ‖x‖ → (w τ x) i = 0 := by
    intro i τ hτ x hx
    obtain ⟨h1, -, h3⟩ := hRb τ hτ.1 x (hRb_le.trans hx)
    simp only [hw, h1, h3, add_zero]; simp
  have hGvan : ∀ m ∈ Finset.range (k + 1), ∀ τ ∈ Icc 0 T', ∀ x : R2, ρ (m + 1) ≤ ‖x‖ → G (m + 1) τ x = 0 :=
    fun m hm τ hτ x hx => ((hL m hm).2.2.2.2.1 τ (hτH τ hτ) x hx).1
  obtain ⟨Mθ0, hMθ0, hθb0⟩ := Kit.paramLinearOp_contDiffOn_aux_bound 0 T' Rall θd hθsm hθvan []
  have hθb1 : ∀ j : Fin 2, ∃ M : ℝ, 0 ≤ M ∧ ∀ τ ∈ Icc 0 T', ∀ x : R2, |pd j (θd τ) x| ≤ M := fun j =>
    Kit.paramLinearOp_contDiffOn_aux_bound 0 T' Rall θd hθsm hθvan [j]
  choose Mθ1 hMθ1 hθb1 using hθb1
  have hwb0 : ∀ i : Fin 2, ∃ M : ℝ, 0 ≤ M ∧ ∀ τ ∈ Icc 0 T', ∀ x : R2, |(w τ x) i| ≤ M := fun i =>
    Kit.paramLinearOp_contDiffOn_aux_bound 0 T' Rall (fun τ x => (w τ x) i) (hwsm i) (hwvan i) []
  choose Mw0 hMw0 hwb0 using hwb0
  have hwb1 : ∀ i j : Fin 2, ∃ M : ℝ, 0 ≤ M ∧ ∀ τ ∈ Icc 0 T', ∀ x : R2, |pd j (fun y => (w τ y) i) x| ≤ M := fun i j =>
    Kit.paramLinearOp_contDiffOn_aux_bound 0 T' Rall (fun τ x => (w τ x) i) (hwsm i) (hwvan i) [j]
  choose Mw1 hMw1 hwb1 using hwb1
  have hGb0 : ∀ m : ℕ, ∃ M : ℝ, 0 ≤ M ∧ (m ∈ Finset.range (k + 1) → ∀ τ ∈ Icc 0 T', ∀ x : R2, |G (m + 1) τ x| ≤ M) := by
    intro m
    by_cases hm : m ∈ Finset.range (k + 1)
    · obtain ⟨M, h0, h⟩ := Kit.paramLinearOp_contDiffOn_aux_bound 0 T' (ρ (m + 1)) (G (m + 1))
        ((hL m hm).2.1.mono hS2) (hGvan m hm) []
      exact ⟨M, h0, fun _ => h⟩
    · exact ⟨0, le_rfl, fun h => absurd h hm⟩
  choose MG0 hMG0 hGb0 using hGb0
  have hGb1 : ∀ m : ℕ, ∀ j : Fin 2, ∃ M : ℝ, 0 ≤ M ∧ (m ∈ Finset.range (k + 1) → ∀ τ ∈ Icc 0 T', ∀ x : R2,
      |pd j (G (m + 1) τ) x| ≤ M) := by
    intro m j
    by_cases hm : m ∈ Finset.range (k + 1)
    · obtain ⟨M, h0, h⟩ := Kit.paramLinearOp_contDiffOn_aux_bound 0 T' (ρ (m + 1)) (G (m + 1))
        ((hL m hm).2.1.mono hS2) (hGvan m hm) [j]
      exact ⟨M, h0, fun _ => h⟩
    · exact ⟨0, le_rfl, fun h => absurd h hm⟩
  choose MG1 hMG1 hGb1 using hGb1
  have hBS : ∀ m : ℕ, ∃ CL CD : ℝ, 0 ≤ CL ∧ 0 ≤ CD ∧ (m ∈ Finset.range (k + 1) → ∀ τ ∈ Icc 0 T',
      ContDiff ℝ ∞ (gv (m + 1) τ) ∧
      eLpNorm (gv (m + 1) τ) 2 volume ≤ ENNReal.ofReal CL ∧ ∀ x, ‖fderiv ℝ (gv (m + 1) τ) x‖ ≤ CD) := by
    intro m
    by_cases hm : m ∈ Finset.range (k + 1)
    · have hρp : 0 < ρ (m + 1) := (hL m hm).1
      obtain ⟨C2, hC2⟩ := biotSavart_memLp hP (ρ (m + 1)) hρp
      obtain ⟨C1, hC1⟩ := biotSavart_fderiv_bounded hP (ρ (m + 1)) hρp
      refine ⟨C2 * MG0 m, max (C1 * (MG1 m 0 + MG1 m 1)) 0, mul_nonneg C2.coe_nonneg (hMG0 m), le_max_right _ _,
        fun _ τ hτ => ?_⟩
      obtain ⟨hGs, hGsupp, hGmean, hGbs⟩ := (hL m hm).2.2.2.2.2 τ (hτH τ hτ)
      have hgvτ : gv (m + 1) τ = biotSavart (G (m + 1) τ) := hGbs
      have hGc : HasCompactSupport (G (m + 1) τ) :=
        HasCompactSupport.of_support_subset_isCompact (isCompact_closedBall (0 : R2) (ρ (m + 1)))
          ((subset_tsupport _).trans (hGsupp.trans Metric.ball_subset_closedBall))
      have hb0 : ∀ y, |G (m + 1) τ y| ≤ MG0 m := hGb0 m hm τ hτ
      have hb1 : ∀ y, ∀ i : Fin 2, |pd i (G (m + 1) τ) y| ≤ MG1 m 0 + MG1 m 1 := by
        intro y i
        have h0 := hGb1 m 0 hm τ hτ y; have h1 := hGb1 m 1 hm τ hτ y
        fin_cases i
        · exact h0.trans (by linarith [hMG1 m 1])
        · exact h1.trans (by linarith [hMG1 m 0])
      rw [hgvτ]
      refine ⟨designed_lipschitz_aux_bs hP _ hGs hGc, ?_, fun x => (hC1 _ hGs hGsupp _ hb1 x).trans (le_max_left _ _)⟩
      refine (hC2 _ hGs hGsupp hGmean _ hb0).2.trans (le_of_eq ?_)
      rw [ENNReal.ofReal_mul C2.coe_nonneg, ENNReal.ofReal_coe_nnreal]
    · exact ⟨0, 0, le_rfl, le_rfl, fun h => absurd h hm⟩
  choose CL CD hCL hCD hBS using hBS
  set VC : ℝ := (volume (Metric.closedBall (0 : R2) Rall)).toReal ^ (1 / 2 : ℝ) with hVC
  have hVC0 : 0 ≤ VC := by positivity
  set Mw1s : ℝ := Mw1 0 0 + Mw1 0 1 + Mw1 1 0 + Mw1 1 1 with hMw1s
  have hMw1s0 : 0 ≤ Mw1s := by have := hMw1 0 0; have := hMw1 0 1; have := hMw1 1 0; have := hMw1 1 1; positivity
  have hMw1le : ∀ i j : Fin 2, Mw1 i j ≤ Mw1s := by
    have := hMw1 0 0; have := hMw1 0 1; have := hMw1 1 0; have := hMw1 1 1
    refine Fin.forall_fin_two.2 ⟨Fin.forall_fin_two.2 ⟨?_, ?_⟩, Fin.forall_fin_two.2 ⟨?_, ?_⟩⟩ <;>
      simp only [hMw1s] <;> linarith
  set ME : ℝ := max (Mθ0 * VC) ((Mw0 0 + Mw0 1) * VC + ∑ m ∈ Finset.range (k + 1), CL m) with hME
  set MD : ℝ := max (2 * (Mθ1 0 + Mθ1 1)) (4 * Mw1s + ∑ m ∈ Finset.range (k + 1), CD m) with hMD
  have hwcont : ∀ τ ∈ Icc 0 T', Differentiable ℝ (w τ) ∧ Continuous (w τ) := by
    intro τ hτ
    have hu := (designed_classical_aux_fields NStar Cg hNuni hP T hT τ (hIco hτ)).1
    have heq : w τ = fun x => ud τ x - ∑ m ∈ Finset.range (k + 1), gv (m + 1) τ x := by
      funext x; rw [husum τ hτ x]; abel
    have hd : Differentiable ℝ (w τ) := by
      rw [heq]
      exact (hu.differentiable (by norm_cast)).sub
        (Differentiable.fun_sum fun m hm => ((hBS m hm τ hτ).1.differentiable (by norm_cast)))
    exact ⟨hd, hd.continuous⟩
  refine ⟨?_, ⟨ME, fun t ht => ?_⟩, ⟨MD, fun t ht x => ?_⟩, fun t ht => ?_, ?_, ?_⟩
  ·
    have h := (hsol.smoothθ.prodMk hsol.smoothu).mono hS
    exact (h.of_le (by norm_cast)).locallyLipschitzOn (convex_Icc 0 T' |>.prod convex_univ)
  ·
    obtain ⟨hmu, hmθ⟩ := hsol.finiteEnergy t (hIco ht)
    refine ⟨hmu, hmθ, ?_, ?_⟩
    · have hud_eq : ud t = fun x => w t x + ∑ m ∈ Finset.range (k + 1), gv (m + 1) t x := funext (husum t ht)
      have hw2 : eLpNorm (w t) 2 volume ≤ ENNReal.ofReal ((Mw0 0 + Mw0 1) * VC) := by
        refine designed_lipschitz_aux_L2 (w t) (Mw0 0 + Mw0 1) Rall (by linarith [hMw0 0, hMw0 1]) (fun x => ?_) (fun x hx => ?_)
        · exact (designed_lipschitz_aux_n2 _).trans (add_le_add (hwb0 0 t ht x) (hwb0 1 t ht x))
        · obtain ⟨h1, -, h3⟩ := hRb t ht.1 x (hRb_le.trans hx)
          simp only [hw, h1, h3, add_zero]
      have hs2 : eLpNorm (fun x => ∑ m ∈ Finset.range (k + 1), gv (m + 1) t x) 2 volume ≤
          ENNReal.ofReal (∑ m ∈ Finset.range (k + 1), CL m) := by
        have hfe : (fun x => ∑ m ∈ Finset.range (k + 1), gv (m + 1) t x) = ∑ m ∈ Finset.range (k + 1), gv (m + 1) t := by
          funext x; simp [Finset.sum_apply]
        rw [hfe]
        refine (eLpNorm_sum_le (fun m hm => ((hBS m hm t ht).1.continuous.aestronglyMeasurable)) (by norm_num)).trans ?_
        rw [ENNReal.ofReal_sum_of_nonneg (fun m _ => hCL m)]
        exact Finset.sum_le_sum fun m hm => (hBS m hm t ht).2.1
      have htot : eLpNorm (ud t) 2 volume ≤ ENNReal.ofReal ((Mw0 0 + Mw0 1) * VC + ∑ m ∈ Finset.range (k + 1), CL m) := by
        rw [hud_eq, ENNReal.ofReal_add (by have := hMw0 0; have := hMw0 1; positivity) (Finset.sum_nonneg fun m _ => hCL m)]
        refine (eLpNorm_add_le (hwcont t ht).2.aestronglyMeasurable ?_ (by norm_num)).trans (add_le_add hw2 hs2)
        exact (continuous_finsetSum _ fun m hm => (hBS m hm t ht).1.continuous).aestronglyMeasurable
      refine ENNReal.toReal_le_of_le_ofReal (by positivity) (htot.trans ?_)
      exact ENNReal.ofReal_le_ofReal (le_max_right _ _)
    · have hθ2 : eLpNorm (θd t) 2 volume ≤ ENNReal.ofReal (Mθ0 * VC) :=
        designed_lipschitz_aux_L2 (θd t) Mθ0 Rall hMθ0 (fun x => by
          have := hθb0 t ht x; rwa [← Real.norm_eq_abs] at this) (hθvan t ht)
      refine ENNReal.toReal_le_of_le_ofReal (by positivity) (hθ2.trans ?_)
      exact ENNReal.ofReal_le_ofReal (le_max_left _ _)
  ·
    constructor
    · refine (designed_lipschitz_aux_opn1 (θd t) x (Mθ1 0 + Mθ1 1) (by linarith [hMθ1 0, hMθ1 1]) fun j => ?_).trans
        (le_max_left _ _)
      fin_cases j
      · exact (hθb1 0 t ht x).trans (by linarith [hMθ1 1])
      · exact (hθb1 1 t ht x).trans (by linarith [hMθ1 0])
    · have hwd := (hwcont t ht).1
      have hgd : ∀ m ∈ Finset.range (k + 1), HasFDerivAt (gv (m + 1) t) (fderiv ℝ (gv (m + 1) t) x) x :=
        fun m hm => (((hBS m hm t ht).1.differentiable (by norm_cast)) x).hasFDerivAt
      have hsumd : HasFDerivAt (fun y => ∑ m ∈ Finset.range (k + 1), gv (m + 1) t y)
          (∑ m ∈ Finset.range (k + 1), fderiv ℝ (gv (m + 1) t) x) x := HasFDerivAt.fun_sum hgd
      have hall : HasFDerivAt (ud t) (fderiv ℝ (w t) x + ∑ m ∈ Finset.range (k + 1), fderiv ℝ (gv (m + 1) t) x) x := by
        have hud_eq : ud t = fun y => w t y + ∑ m ∈ Finset.range (k + 1), gv (m + 1) t y := funext (husum t ht)
        rw [hud_eq]
        exact (hwd x).hasFDerivAt.add hsumd
      rw [hall.fderiv]
      have hw4 : ‖fderiv ℝ (w t) x‖ ≤ 4 * Mw1s :=
        designed_lipschitz_aux_opn2 (w t) x (hwd x) Mw1s hMw1s0 fun i j => (hwb1 i j t ht x).trans (hMw1le i j)
      refine (norm_add_le _ _).trans ((add_le_add hw4 ((norm_sum_le _ _).trans
        (Finset.sum_le_sum fun m hm => (hBS m hm t ht).2.2 x))).trans (le_max_right _ _))
  ·
    exact ae_of_all _ (hsol.divFree t (hIco ht))
  ·
    have hmeas : MeasurableSet (Icc 0 T' ×ˢ (univ : Set R2)) := measurableSet_Icc.prod MeasurableSet.univ
    exact (ae_restrict_iff' hmeas).2 (ae_of_all _ fun z hz => hsol.eqθ z.1 (hIco hz.1) z.2)
  ·
    intro φ hφ hφc hφs hφdiv
    exact designed_lipschitz_aux_weak κ T.Tstar T' hT'lt.le θd p ud fu hsol.smoothu hsol.smoothp
      hsol.smoothθ.continuousOn (fun i => (hfu_sm i).continuousOn) hsol.divFree hsol.equ φ hφ hφc hφs hφdiv

theorem designed_data {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ)
    (hNuni : ∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) (T : TowerData κ Pr)
    (hT : TowerSpec NStar Cg T) :
    HasCompactSupport (S7.Designed.theta T 0) ∧ HasCompactSupport (S7.Designed.vel T 0) ∧
    S7.Designed.vel T 0 = 0 := by
  classical
  have hred1 : T.StateRed 1 := (hT.state 1 le_rfl).1
  have ht0 : T.t 0 = 0 := hred1.t_zero
  have hmono : Monotone T.t := by
    intro a b hab
    rcases eq_or_lt_of_le hab with h | h
    · rw [h]
    · exact ((hT.state (b + 1) (by omega)).1.t_mono (show a ≤ b + 1 by omega) (show b ≤ b + 1 by omega) h).le
  have ht1 : 0 < T.t 1 := by
    have := hred1.t_mono (show 0 ≤ 1 from Nat.zero_le _) (show 1 ≤ 1 from le_rfl) (by omega)
    rwa [ht0] at this
  have hins : ∀ m, T.tIn m ≤ 0 ↔ m ≤ 1 := by
    intro m
    unfold TowerData.tIn
    constructor
    · intro h
      by_contra hm
      push_neg at hm
      have : T.t 1 ≤ T.t (m - 1) := hmono (by omega)
      linarith
    · intro hm
      have : T.t (m - 1) = T.t 0 := by congr 1; omega
      rw [this, ht0]
  have hvan := S7.layer_vanishes_at_insertion NStar Cg hNuni T hT 1 le_rfl
  have htIn1 : T.tIn 1 = 0 := by unfold TowerData.tIn; simpa using ht0
  have hθ : S7.Designed.theta T 0 = T.θL 0 0 := by
    funext x
    unfold S7.Designed.theta
    rw [tsum_eq_sum (s := Finset.range 2) (fun m hm => ?_)]
    · rw [Finset.sum_range_succ, Finset.sum_range_one, if_pos ((hins 0).2 (by omega)), if_pos ((hins 1).2 le_rfl)]
      have : T.θL 1 0 x = 0 := by have := (hvan x).1; rwa [htIn1] at this
      rw [this, add_zero]
    · rw [if_neg ((hins m).not.2 (by have := Finset.mem_range.not.1 hm; omega))]
  have hderiv0 : deriv T.ψ 0 = 0 := by
    by_cases hd : DifferentiableAt ℝ T.ψ 0
    ·
      have hfree := hT.base.free_window
      have hL : 0 < T.base.L1free / Real.sqrt (κ * T.base.A0) := by
        rw [← hT.base.manStart1]
        have := (hred1.man_mem 1 le_rfl le_rfl).1
        rwa [ht0] at this
      have hconst : T.ψ =ᶠ[𝓝[Ici (0 : ℝ)] 0] fun _ : ℝ => T.ψ 0 := by
        have hmem : Icc 0 (T.base.L1free / Real.sqrt (κ * T.base.A0)) ∈ 𝓝[Ici (0 : ℝ)] (0 : ℝ) := Icc_mem_nhdsGE hL
        filter_upwards [hmem] with s hs
        rw [hfree s hs, hfree 0 ⟨le_rfl, hL.le⟩]
      have hwithin : derivWithin T.ψ (Ici 0) 0 = 0 := by
        have : derivWithin T.ψ (Ici 0) 0 = derivWithin (fun _ : ℝ => T.ψ 0) (Ici 0) 0 :=
          Filter.EventuallyEq.derivWithin_eq hconst rfl
        rw [this]; simp
      rw [← hd.hasDerivAt.hasDerivWithinAt.derivWithin (uniqueDiffWithinAt_Ici 0), hwithin]
    · exact deriv_zero_of_not_differentiableAt hd
  have hOm0 : T.Om0 0 = 0 := by unfold TowerData.Om0; simp
  have hbase : T.uL 0 0 = 0 := by
    funext x
    rw [hT.base.u0 0 le_rfl x]
    have : T.basePsi 0 = fun _ => 0 := by
      funext y; unfold TowerData.basePsi; rw [hOm0]; simp
    rw [this]
    simp [perpGrad, pd]
  have hu : S7.Designed.vel T 0 = 0 := by
    funext x
    unfold S7.Designed.vel
    have hst : T.uSt 0 x = 0 := by
      show deriv T.ψ 0 • perpGrad Pr.H x = 0
      rw [hderiv0, zero_smul]
    rw [tsum_eq_sum (s := Finset.range 2) (fun m hm => ?_)]
    · rw [Finset.sum_range_succ, Finset.sum_range_one, if_pos ((hins 0).2 (by omega)), if_pos ((hins 1).2 le_rfl),
        hst, hbase]
      have : T.uL 1 0 x = 0 := by have := (hvan x).2.2; rwa [htIn1] at this
      rw [this]; simp
    · rw [if_neg ((hins m).not.2 (by have := Finset.mem_range.not.1 hm; omega))]
  have hθsupp : HasCompactSupport (T.θL 0 0) := by
    have hfun : T.θL 0 0 = fun x => T.base.Theta0 * Real.sin (inner ℝ (T.xi 0 0) x) * Pr.χ0 (S6.rot (-(T.ψ 0)) x) :=
      funext fun x => hT.base.theta0 0 le_rfl x
    rw [hfun]
    have hχ : HasCompactSupport Pr.χ0 := Pr.hχ.compact
    obtain ⟨Rχ, hRχ⟩ : ∃ r, tsupport Pr.χ0 ⊆ Metric.closedBall (0 : R2) r := hχ.isBounded.subset_closedBall 0
    refine HasCompactSupport.of_support_subset_isCompact (isCompact_closedBall (0 : R2) Rχ) ?_
    intro x hx
    rw [Function.mem_support] at hx
    have hχx : Pr.χ0 (S6.rot (-(T.ψ 0)) x) ≠ 0 := fun h => hx (by rw [h, mul_zero])
    have hmem : S6.rot (-(T.ψ 0)) x ∈ tsupport Pr.χ0 := subset_tsupport _ (Function.mem_support.2 hχx)
    have := hRχ hmem
    rw [Metric.mem_closedBall, dist_zero_right, S7.parked_node_vorticity_aux_rot_norm] at this
    rwa [Metric.mem_closedBall, dist_zero_right]
  refine ⟨by rw [hθ]; exact hθsupp, by rw [hu]; exact HasCompactSupport.zero, hu⟩

theorem designed_pair_solves {κ : ℝ} {Pr : Profiles} (NStar Cg : ℝ → ℝ) (Kexp : ℕ)
    (hNdom : ∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg Kexp) Cg β ≤ NStar β)
    (hNuni : ∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β)
    (hκ : 0 < κ) (hP : PotentialTheoryR2) (hB : BogovskiiAnnulusR2)
    (T : TowerData κ Pr) (hT : TowerSpec NStar Cg T)
    (fu : ℝ → R2 → R2) (hfu : ∀ t ∈ Ico 0 T.Tstar, ∀ x, curl2 (fu t) x = S7.Designed.fOmega T t x)
    (hfu_odd : ∀ t ∈ Ico 0 T.Tstar, ∀ x, fu t (-x) = -fu t x)
    (hfu_supp : ∃ Rball : ℝ, ∀ t ∈ Ico 0 T.Tstar, ∀ x : R2, Rball < ‖x‖ → fu t x = 0)
    (hfu_sm : ∀ i : Fin 2, ContDiffOn ℝ ∞ (fun z : ℝ × R2 => fu z.1 z.2 i) (Ico 0 T.Tstar ×ˢ univ)) :
    ∃ p : ℝ → R2 → ℝ,
      ClassicalBoussinesq κ T.Tstar (S7.Designed.theta T) p (S7.Designed.vel T) fu (S7.Designed.fTheta T) ∧
      (∀ t ∈ Ico 0 T.Tstar, InClassP (S7.Designed.theta T t) (S7.Designed.vel T t)) ∧
      (∀ t ∈ Ico 0 T.Tstar, InClassP (S7.Designed.fTheta T t) (fu t)) ∧
      (∀ T' < T.Tstar, 0 < T' →
        InLipschitzClass κ T' (S7.Designed.theta T) (S7.Designed.vel T) (S7.Designed.fTheta T) fu) ∧
      HasCompactSupport (S7.Designed.theta T 0) ∧ HasCompactSupport (S7.Designed.vel T 0) := by
  have hfields := fun τ (hτ0 : 0 < τ) (hτ : τ < T.Tstar) => S7.lemma751_fields_aux_fields NStar Cg hNuni T hT τ hτ0 hτ
  have hθ : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => S7.Designed.theta T z.1 z.2) (Ico 0 T.Tstar ×ˢ univ) :=
    S7.lemma751_fields_aux_local T.Tstar (fun z : ℝ × R2 => S7.Designed.theta T z.1 z.2) fun τ h0 h => (hfields τ h0 h).1
  have hu : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => S7.Designed.vel T z.1 z.2) (Ico 0 T.Tstar ×ˢ univ) :=
    S7.lemma751_fields_aux_local T.Tstar (fun z : ℝ × R2 => S7.Designed.vel T z.1 z.2) fun τ h0 h => (hfields τ h0 h).2.2
  have hfθ : ContDiffOn ℝ ∞ (fun z : ℝ × R2 => S7.Designed.fTheta T z.1 z.2) (Ico 0 T.Tstar ×ˢ univ) := by
    refine S7.lemma751_fields_aux_local T.Tstar (fun z : ℝ × R2 => S7.Designed.fTheta T z.1 z.2) fun τ h0 h => ?_
    have hT₂ : (τ + T.Tstar) / 2 < T.Tstar := by linarith
    have hτ2 : τ < (τ + T.Tstar) / 2 := by linarith
    obtain ⟨hθ2, -, hu2⟩ := hfields ((τ + T.Tstar) / 2) (by linarith) hT₂
    have h1 := S7.lemma751_fields_aux_pdt (S7.Designed.theta T) ((τ + T.Tstar) / 2) τ h0 hτ2 hθ2
    have h2 := S7.lemma751_fields_aux_advect (S7.Designed.theta T) (S7.Designed.vel T) ((τ + T.Tstar) / 2) hθ2 hu2
    have hsub : Icc 0 τ ×ˢ (univ : Set R2) ⊆ Icc 0 ((τ + T.Tstar) / 2) ×ˢ univ :=
      prod_mono (Icc_subset_Icc le_rfl hτ2.le) subset_rfl
    show ContDiffOn ℝ ∞ (fun z : ℝ × R2 => pdt (S7.Designed.theta T) z.1 z.2 +
      advect (S7.Designed.vel T z.1) (S7.Designed.theta T z.1) z.2) (Icc 0 τ ×ˢ univ)
    exact h1.add (h2.mono hsub)
  obtain ⟨p, hsol⟩ := designed_classical NStar Cg hNuni hκ hP T hT fu hfu hfu_sm hθ hu hfθ
  obtain ⟨hclass, hclassF⟩ := designed_parity NStar Cg Kexp hNdom hNuni T hT fu hfu_odd
  have hlip := designed_lipschitz NStar Cg Kexp hNdom hNuni hκ hP T hT fu p hsol hfu_sm hfu_supp
  obtain ⟨hθ0, hu0, -⟩ := designed_data NStar Cg hNuni T hT
  exact ⟨p, hsol, hclass, hclassF, hlip, hθ0, hu0⟩

theorem towerExists_discharged_aux (hB : BogovskiiAnnulusR2) (hP : PotentialTheoryR2) (hW : LandingWindowsW) (hU : FlatWindowsUndressed) :
    ∃ Pr : Profiles, ∀ κ : ℝ, 0 < κ → κ ≤ 1e5 → |Real.log κ| ≤ 30 →
      ∃ Cg : ℝ → ℝ, ∀ Nextra : ℝ → ℝ, ∃ NStar : ℝ → ℝ,
        (∀ β : ℝ, 0 < β → Nextra β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg 1) Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrRates Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrForce Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrShear Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrReg Pr Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrTail Pr Cg β ≤ NStar β) ∧
        ∃ T : S6.TowerData κ Pr, S6.TowerSpec NStar Cg T := by
  obtain ⟨Pr, hBump⟩ := S4.profiles_exist
  refine ⟨Pr, fun κ hκ hκ' hlog => ?_⟩
  obtain ⟨C₁, hbase⟩ := S6.base_exists κ hκ hκ' Pr hU hP
  obtain ⟨Cg, hE, Cε, hCε, Nrow, hbg, hR1, Ct, hCt, hTay⟩ :=
    D.pde_rows_discharge (κ := κ) (Pr := Pr) hκ hP (fun β => max C₁ (8 / β + 4))
  have hCg : ∀ β : ℝ, 0 < β → β ≤ 1 / 80 → 8 / β + 4 ≤ Cg β := fun β _ _ => (le_max_right _ _).trans (hE β)
  have hC1 : C₁ ≤ Cg (1 / 80) := (le_max_left _ _).trans (hE _)
  refine ⟨Cg, fun Nextra => ?_⟩
  obtain ⟨NStar, hNdom, hNuni, hNrow, hNabs, hNrate, hNforce, hNshear, hNreg, hNtail, hNextra, hNcross⟩ :
      ∃ NStar : ℝ → ℝ,
        (∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg 1) Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → Nrow β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S6.NthrAbs β (Cε β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrRates Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrForce Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrShear Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrReg Pr Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrTail Pr Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → Nextra β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrCross hP Pr β (Cg β) ≤ NStar β) := by
    refine ⟨fun β => max (max (max (max (S6.NStar (S6.thr4Of Pr Cg 1) Cg β) (S4.NthrUniform Pr β (Cg β)))
        (max (Nrow β) (S6.NthrAbs β (Cε β)))) (max (max (S4.NthrRates Pr β (Cg β)) (S4.NthrForce Pr β (Cg β)))
        (max (S7.NthrShear Pr β (Cg β)) (S7.NthrReg Pr Cg β)))) (max (max (S7.NthrTail Pr Cg β) (Nextra β)) (S4.NthrCross hP Pr β (Cg β))),
      ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> intro β _ <;>
      simp only [le_max_iff, le_refl, true_or, or_true]
  have hCross : S6.CrossRow (Pr := Pr) Cg NStar := by
    intro p hβ hN κ' B R K Nprev hG hS hA t ht
    exact S4.nthrCross_spec hP Pr p.β (Cg p.β) p rfl ((hNcross p.β hβ).trans hN) R K Nprev hG hS hA t ht
  obtain ⟨T₁, hB1, hS1, hCB1, hSR1, hAbs1, hSJ1, hN1, -, -, hL1, hlogs1, -, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ :=
    hbase Cg Cε hC1 hCε Nrow ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩ NStar hNrow
  obtain ⟨T, hT⟩ := S6.towerExists NStar Cg Cε Nrow 1 hNdom hNuni hNrow hNabs hCross hCg hCε ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩
    hP hB hBump hκ hκ' hlog ⟨T₁, hB1, hS1, hN1, hCB1, hL1, hlogs1, hSR1, hAbs1, hSJ1, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ hW
  exact ⟨NStar, hNextra, hNdom, hNuni, hNrate, hNforce, hNshear, hNreg, hNtail, T, hT⟩

theorem towerExists_discharged (hB : BogovskiiAnnulusR2) (hP : PotentialTheoryR2) (hW : LandingWindowsW) (hU : FlatWindowsUndressed) :
    ∃ Pr : Profiles, ∀ κ : ℝ, 0 < κ → κ ≤ 1e5 → |Real.log κ| ≤ 30 →
      ∃ (NStar Cg : ℝ → ℝ) (T : S6.TowerData κ Pr), S6.TowerSpec NStar Cg T := by

  obtain ⟨Pr, hBump⟩ := S4.profiles_exist
  refine ⟨Pr, fun κ hκ hκ' hlog => ?_⟩
  obtain ⟨C₁, hbase⟩ := S6.base_exists κ hκ hκ' Pr hU hP
  obtain ⟨Cg, hE, Cε, hCε, Nrow, hbg, hR1, Ct, hCt, hTay⟩ :=
    D.pde_rows_discharge (κ := κ) (Pr := Pr) hκ hP (fun β => max C₁ (8 / β + 4))
  have hCg : ∀ β : ℝ, 0 < β → β ≤ 1 / 80 → 8 / β + 4 ≤ Cg β := fun β _ _ => (le_max_right _ _).trans (hE β)
  have hC1 : C₁ ≤ Cg (1 / 80) := (le_max_left _ _).trans (hE _)
  obtain ⟨NStar, hNdom, hNuni, hNrow, hNabs, hNcross⟩ :
      ∃ NStar : ℝ → ℝ,
        (∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg 1) Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → Nrow β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S6.NthrAbs β (Cε β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrCross hP Pr β (Cg β) ≤ NStar β) := by
    refine ⟨fun β => max (max (max (S6.NStar (S6.thr4Of Pr Cg 1) Cg β) (S4.NthrUniform Pr β (Cg β)))
        (max (Nrow β) (S6.NthrAbs β (Cε β)))) (S4.NthrCross hP Pr β (Cg β)), ?_, ?_, ?_, ?_, ?_⟩ <;> intro β _ <;>
      simp only [le_max_iff, le_refl, true_or, or_true]
  have hCross : S6.CrossRow (Pr := Pr) Cg NStar := by
    intro p hβ hN κ' B R K Nprev hG hS hA t ht
    exact S4.nthrCross_spec hP Pr p.β (Cg p.β) p rfl ((hNcross p.β hβ).trans hN) R K Nprev hG hS hA t ht
  obtain ⟨T₁, hB1, hS1, hCB1, hSR1, hAbs1, hSJ1, hN1, -, -, hL1, hlogs1, -, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ :=
    hbase Cg Cε hC1 hCε Nrow ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩ NStar hNrow
  obtain ⟨T, hT⟩ := S6.towerExists NStar Cg Cε Nrow 1 hNdom hNuni hNrow hNabs hCross hCg hCε ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩
    hP hB hBump hκ hκ' hlog ⟨T₁, hB1, hS1, hN1, hCB1, hL1, hlogs1, hSR1, hAbs1, hSJ1, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ hW
  exact ⟨NStar, Cg, T, hT⟩

theorem tower_readout (hB : BogovskiiAnnulusR2) (hP : PotentialTheoryR2) (hW : LandingWindowsW) (hU : FlatWindowsUndressed) :
    ∃ Pr : Profiles, ∀ κ : ℝ, 0 < κ → κ ≤ 1e5 → |Real.log κ| ≤ 30 →
      ∃ (NStar Cg : ℝ → ℝ) (T : S6.TowerData κ Pr) (h : S7.StaircaseFacts T), S6.TowerSpec NStar Cg T ∧
        (∀ τ ∈ Ico 0 T.Tstar, ∃ x,
          gradSize (S7.Designed.theta T τ) x ≥ (S7.staircaseOf T h).A ((S7.staircaseOf T h).stage τ) / 2) ∧
        Tendsto (fun τ => (S7.staircaseOf T h).A ((S7.staircaseOf T h).stage τ)) (𝓝[<] T.Tstar) atTop ∧
        (∃ Mθ : ℝ, ∀ t ∈ Ico 0 T.Tstar, ∀ x, |S7.Designed.theta T t x| ≤ Mθ) ∧
        (∀ M : ℝ, ∃ᶠ τ in 𝓝[<] T.Tstar, ∃ x, |curl2 (S7.Designed.vel T τ) x| ≥ M) := by
  obtain ⟨Pr, hBump⟩ := S4.profiles_exist
  refine ⟨Pr, fun κ hκ hκ' hlog => ?_⟩
  obtain ⟨C₁, hbase⟩ := S6.base_exists κ hκ hκ' Pr hU hP
  obtain ⟨Cg, hE, Cε, hCε, Nrow, hbg, hR1, Ct, hCt, hTay⟩ :=
    D.pde_rows_discharge (κ := κ) (Pr := Pr) hκ hP (fun β => max C₁ (8 / β + 4))
  have hCg : ∀ β : ℝ, 0 < β → β ≤ 1 / 80 → 8 / β + 4 ≤ Cg β := fun β _ _ => (le_max_right _ _).trans (hE β)
  have hC1 : C₁ ≤ Cg (1 / 80) := (le_max_left _ _).trans (hE _)
  obtain ⟨NStar, hNdom, hNuni, hNrow, hNabs, hNcross, hNshear⟩ :
      ∃ NStar : ℝ → ℝ,
        (∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg 1) Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → Nrow β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S6.NthrAbs β (Cε β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrCross hP Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrShear Pr β (Cg β) ≤ NStar β) := by
    refine ⟨fun β => max (max (max (max (S6.NStar (S6.thr4Of Pr Cg 1) Cg β) (S4.NthrUniform Pr β (Cg β)))
        (max (Nrow β) (S6.NthrAbs β (Cε β)))) (S4.NthrCross hP Pr β (Cg β))) (S7.NthrShear Pr β (Cg β)),
        ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> intro β _ <;>
      simp only [le_max_iff, le_refl, true_or, or_true]
  have hCross : S6.CrossRow (Pr := Pr) Cg NStar := by
    intro p hβ hN κ' B R K Nprev hG hS hA t ht
    exact S4.nthrCross_spec hP Pr p.β (Cg p.β) p rfl ((hNcross p.β hβ).trans hN) R K Nprev hG hS hA t ht
  obtain ⟨T₁, hB1, hS1, hCB1, hSR1, hAbs1, hSJ1, hN1, -, -, hL1, hlogs1, -, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ :=
    hbase Cg Cε hC1 hCε Nrow ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩ NStar hNrow
  obtain ⟨T, hT⟩ := S6.towerExists NStar Cg Cε Nrow 1 hNdom hNuni hNrow hNabs hCross hCg hCε ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩
    hP hB hBump hκ hκ' hlog ⟨T₁, hB1, hS1, hN1, hCB1, hL1, hlogs1, hSR1, hAbs1, hSJ1, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ hW
  have hsf : S7.StaircaseFacts T := S7.staircase_facts NStar Cg T hT
  obtain ⟨hg1, hg2⟩ := S7.readout_grad NStar Cg 1 hNdom hNuni hP T hT hsf
  exact ⟨NStar, Cg, T, hsf, hT, hg1, hg2, S7.theta_bounded NStar Cg 1 hNdom hNuni T hT,
    S7.readout_vorticity NStar Cg 1 hNdom hNuni hNshear hP T hT⟩

theorem theorem01_at (hB : BogovskiiAnnulusR2) (hP : PotentialTheoryR2) (hW : LandingWindowsW) (hU : FlatWindowsUndressed) :
    Theorem01_Body 1e5 := by
  have hκ5 : (0 : ℝ) < 1e5 := by norm_num
  have hlog : |Real.log (1e5 : ℝ)| ≤ 30 := by
    have h0 : 0 ≤ Real.log (1e5 : ℝ) := Real.log_nonneg (by norm_num)
    rw [abs_of_nonneg h0]
    have h1 : Real.log (1e5 : ℝ) < Real.log ((2 : ℝ) ^ 17) := Real.log_lt_log (by norm_num) (by norm_num)
    rw [Real.log_pow] at h1
    have h2 := Real.log_two_lt_d9
    push_cast at h1
    nlinarith
  obtain ⟨Pr, hBump⟩ := S4.profiles_exist
  obtain ⟨C₁, hbase⟩ := S6.base_exists (1e5 : ℝ) hκ5 le_rfl Pr hU hP
  obtain ⟨Cg, hE, Cε, hCε, Nrow, hbg, hR1, Ct, hCt, hTay⟩ :=
    D.pde_rows_discharge (κ := (1e5 : ℝ)) (Pr := Pr) hκ5 hP (fun β => max C₁ (8 / β + 4))
  have hCg : ∀ β : ℝ, 0 < β → β ≤ 1 / 80 → 8 / β + 4 ≤ Cg β := fun β _ _ => (le_max_right _ _).trans (hE β)
  have hC1 : C₁ ≤ Cg (1 / 80) := (le_max_left _ _).trans (hE _)
  have hSS : S4.StructuredShearTimeRows Pr := S4.structuredShearTimeRows Pr
  obtain ⟨NStar, hNdom, hNuni, hNrow, hNabs, hNrate, hNforce, hNshear, hNreg, hNtail, hNss, hNcross⟩ :
      ∃ NStar : ℝ → ℝ,
        (∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg 1) Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → Nrow β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S6.NthrAbs β (Cε β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrRates Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrForce Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrShear Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrReg Pr Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrTail Pr Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrSS Pr hSS β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrCross hP Pr β (Cg β) ≤ NStar β) := by
    refine ⟨fun β => max (max (max (max (max (S6.NStar (S6.thr4Of Pr Cg 1) Cg β) (S4.NthrUniform Pr β (Cg β)))
        (max (Nrow β) (S6.NthrAbs β (Cε β)))) (max (max (S4.NthrRates Pr β (Cg β)) (S4.NthrForce Pr β (Cg β)))
        (max (S7.NthrShear Pr β (Cg β)) (S7.NthrReg Pr Cg β)))) (max (S7.NthrTail Pr Cg β) (S4.NthrSS Pr hSS β (Cg β)))) (S4.NthrCross hP Pr β (Cg β)),
      ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> intro β _ <;>
      simp only [le_max_iff, le_refl, true_or, or_true]
  have hCross : S6.CrossRow (Pr := Pr) Cg NStar := by
    intro p hβ hN κ' B R K Nprev hG hS hA t ht
    exact S4.nthrCross_spec hP Pr p.β (Cg p.β) p rfl ((hNcross p.β hβ).trans hN) R K Nprev hG hS hA t ht
  obtain ⟨T₁, hB1, hS1, hCB1, hSR1, hAbs1, hSJ1, hN1, -, -, hL1, hlogs1, -, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ :=
    hbase Cg Cε hC1 hCε Nrow ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩ NStar hNrow
  obtain ⟨T, hT⟩ := S6.towerExists NStar Cg Cε Nrow 1 hNdom hNuni hNrow hNabs hCross hCg hCε ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩
    hP hB hBump hκ5 le_rfl hlog
    ⟨T₁, hB1, hS1, hN1, hCB1, hL1, hlogs1, hSR1, hAbs1, hSJ1, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ hW
  obtain ⟨fu, Rball, hcurl, hodd, hsupp, hSFθ, hSFu⟩ :=
    S7.lemma751 NStar Cg 1 hNdom hNuni hNrate hNforce hNshear hNreg hNtail hSS hNss hP hB T hT
  obtain ⟨p, hsol, hclass, hclassF, hlip, hθ0, hu0⟩ :=
    designed_pair_solves NStar Cg 1 hNdom hNuni hκ5 hP hB T hT fu hcurl hodd
      ⟨Rball, fun t ht x hx => (hsupp t ht x hx).1⟩ (fun i => (hSFu i).smooth)
  have hsf := S7.staircase_facts NStar Cg T hT
  obtain ⟨hgrad, htend⟩ := S7.readout_grad NStar Cg 1 hNdom hNuni hP T hT hsf
  obtain ⟨Mθ, hMθ⟩ := S7.theta_bounded NStar Cg 1 hNdom hNuni T hT
  have hvort := S7.readout_vorticity NStar Cg 1 hNdom hNuni hNshear hP T hT
  have hTpos : 0 < T.Tstar := by rw [← hsf.t_zero]; exact hsf.t_lt 0
  have hfθc : ∀ T' < T.Tstar, ContinuousOn (fun z : ℝ × R2 => S7.Designed.fTheta T z.1 z.2) (Icc 0 T' ×ˢ univ) :=
    fun T' hT' => hSFθ.smooth.continuousOn.mono (prod_mono (Icc_subset_Ico_right hT') subset_rfl)
  have hsubT : ∀ T' < T.Tstar, Icc 0 T' ×ˢ (univ : Set R2) ⊆ Ico 0 T.Tstar ×ˢ univ :=
    fun T' hT' => prod_mono (Icc_subset_Ico_right hT') subset_rfl
  have hfuc : ∀ T' < T.Tstar, ContinuousOn (fun z : ℝ × R2 => fu z.1 z.2) (Icc 0 T' ×ˢ univ) := by
    intro T' hT'
    have h0 := (hSFu 0).smooth.continuousOn.mono (hsubT T' hT')
    have h1 := (hSFu 1).smooth.continuousOn.mono (hsubT T' hT')
    have h := (h0.smul (continuousOn_const (c := (e 0 : R2)))).add (h1.smul (continuousOn_const (c := (e 1 : R2))))
    refine h.congr fun z _ => ?_
    show fu z.1 z.2 = fu z.1 z.2 0 • e 0 + fu z.1 z.2 1 • e 1
    ext j; fin_cases j <;> simp [e]
  have hL2 : ∀ T' < T.Tstar, ∃ C : ℝ, ∀ t ∈ Icc 0 T', MemLp (S7.Designed.fTheta T t) 2 volume ∧ MemLp (fu t) 2 volume ∧
      (eLpNorm (S7.Designed.fTheta T t) 2 volume).toReal ≤ C ∧ (eLpNorm (fu t) 2 volume).toReal ≤ C := by
    intro T' hT'
    set R : ℝ := max Rball 1 with hRdef
    have hR1 : 1 ≤ R := le_max_right _ _
    have hK : IsCompact (Icc 0 T' ×ˢ Metric.closedBall (0 : R2) R) := isCompact_Icc.prod (isCompact_closedBall _ _)
    have hKsub : Icc 0 T' ×ˢ Metric.closedBall (0 : R2) R ⊆ Icc 0 T' ×ˢ (univ : Set R2) := prod_mono subset_rfl (subset_univ _)
    obtain ⟨Bθ, hBθ⟩ := hK.exists_bound_of_continuousOn ((hfθc T' hT').mono hKsub)
    obtain ⟨Bu, hBu⟩ := hK.exists_bound_of_continuousOn ((hfuc T' hT').mono hKsub)
    have hzero : ∀ t ∈ Icc 0 T', ∀ x : R2, x ∉ Metric.closedBall (0 : R2) R → fu t x = 0 ∧ S7.Designed.fTheta T t x = 0 := by
      intro t ht x hx
      rw [Metric.mem_closedBall, dist_zero_right, not_le] at hx
      exact hsupp t ⟨ht.1, ht.2.trans_lt hT'⟩ x ((le_max_left _ _).trans_lt hx)
    have hmeas : MeasurableSet (Metric.closedBall (0 : R2) R) := measurableSet_closedBall
    have hμ : volume (Metric.closedBall (0 : R2) R) ≠ ⊤ := (isCompact_closedBall _ _).measure_lt_top.ne
    set B : ℝ := max (max Bθ Bu) 0 with hBdef
    set g : R2 → ℝ := (Metric.closedBall (0 : R2) R).indicator (fun _ => B) with hgdef
    have hgLp : MemLp g 2 volume := memLp_indicator_const 2 hmeas B (Or.inr hμ)
    have hdomθ : ∀ t ∈ Icc 0 T', ∀ x, ‖S7.Designed.fTheta T t x‖ ≤ ‖g x‖ := by
      intro t ht x
      by_cases hx : x ∈ Metric.closedBall (0 : R2) R
      · rw [hgdef, indicator_of_mem hx, Real.norm_of_nonneg (le_max_right _ _)]
        exact (hBθ (t, x) ⟨ht, hx⟩).trans ((le_max_left _ _).trans (le_max_left _ _))
      · rw [(hzero t ht x hx).2, norm_zero]; exact norm_nonneg _
    have hdomu : ∀ t ∈ Icc 0 T', ∀ x, ‖fu t x‖ ≤ ‖g x‖ := by
      intro t ht x
      by_cases hx : x ∈ Metric.closedBall (0 : R2) R
      · rw [hgdef, indicator_of_mem hx, Real.norm_of_nonneg (le_max_right _ _)]
        exact (hBu (t, x) ⟨ht, hx⟩).trans ((le_max_right _ _).trans (le_max_left _ _))
      · rw [(hzero t ht x hx).1, norm_zero]; exact norm_nonneg _
    have hsliceθ : ∀ t ∈ Icc 0 T', Continuous (S7.Designed.fTheta T t) := fun t ht =>
      (hfθc T' hT').comp_continuous (Continuous.prodMk_right t) fun x => ⟨ht, mem_univ x⟩
    have hsliceu : ∀ t ∈ Icc 0 T', Continuous (fu t) := fun t ht =>
      (hfuc T' hT').comp_continuous (Continuous.prodMk_right t) fun x => ⟨ht, mem_univ x⟩
    refine ⟨(eLpNorm g 2 volume).toReal, fun t ht => ⟨?_, ?_, ?_, ?_⟩⟩
    · exact hgLp.of_le (hsliceθ t ht).aestronglyMeasurable (Filter.Eventually.of_forall (hdomθ t ht))
    · exact hgLp.of_le (hsliceu t ht).aestronglyMeasurable (Filter.Eventually.of_forall (hdomu t ht))
    · exact ENNReal.toReal_mono hgLp.eLpNorm_ne_top (eLpNorm_mono fun x => hdomθ t ht x)
    · exact ENNReal.toReal_mono hgLp.eLpNorm_ne_top (eLpNorm_mono fun x => hdomu t ht x)
  refine ⟨T.Tstar, hTpos, S7.staircaseOf T hsf, S7.Designed.theta T, p, S7.Designed.fTheta T, S7.Designed.vel T, fu,
    hsol, hclass, hclassF, hlip, hθ0, hu0, hSFθ, hSFu, hgrad, htend, ⟨Mθ, hMθ⟩, hvort, ?_⟩
  intro T' hT' hT'0 θ' u' h' h0 h0'
  exact S7.uniqueness_lipschitz hP 1e5 T' hT'0 _ θ' _ u' _ fu (hfθc T' hT') (hfuc T' hT') (hL2 T' hT')
    (hlip T' hT' hT'0) h' h0 h0'

theorem theorem01 (hB : BogovskiiAnnulusR2) (hP : PotentialTheoryR2) (hW : LandingWindowsW) (hU : FlatWindowsUndressed) :
    Theorem01_Statement := by
  rw [theorem01_body_iff]
  intro κ hκ
  exact boussinesq_rescale κ 1e5 hκ (by norm_num) (theorem01_at hB hP hW hU)

def Lemma751Statement : Prop :=
  ∀ {Pr : Profiles} (Cg : ℝ → ℝ), ∃ Nss : ℝ → ℝ, ∀ {κ : ℝ} (NStar : ℝ → ℝ) (Kexp : ℕ),
    (∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg Kexp) Cg β ≤ NStar β) →
    (∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) →
    (∀ β : ℝ, 0 < β → S4.NthrRates Pr β (Cg β) ≤ NStar β) →
    (∀ β : ℝ, 0 < β → S4.NthrForce Pr β (Cg β) ≤ NStar β) →
    (∀ β : ℝ, 0 < β → S7.NthrShear Pr β (Cg β) ≤ NStar β) →
    (∀ β : ℝ, 0 < β → S7.NthrReg Pr Cg β ≤ NStar β) →
    (∀ β : ℝ, 0 < β → S7.NthrTail Pr Cg β ≤ NStar β) →
    (∀ β : ℝ, 0 < β → Nss β ≤ NStar β) →
    PotentialTheoryR2 → BogovskiiAnnulusR2 →
    ∀ (T : TowerData κ Pr), TowerSpec NStar Cg T →
      ∃ fu : ℝ → R2 → R2, ∃ Rball : ℝ,
        (∀ t ∈ Ico 0 T.Tstar, ∀ x, curl2 (fu t) x = S7.Designed.fOmega T t x) ∧
        (∀ t ∈ Ico 0 T.Tstar, ∀ x, fu t (-x) = -fu t x) ∧
        (∀ t ∈ Ico 0 T.Tstar, ∀ x : R2, Rball < ‖x‖ → fu t x = 0 ∧ S7.Designed.fTheta T t x = 0) ∧
        SmoothForce T.Tstar (fun z : ℝ × R2 => S7.Designed.fTheta T z.1 z.2) ∧
        (∀ i : Fin 2, SmoothForce T.Tstar (fun z : ℝ × R2 => fu z.1 z.2 i))

theorem theorem01_modulo_751 (hB : BogovskiiAnnulusR2) (hP : PotentialTheoryR2) (hW : LandingWindowsW) (hU : FlatWindowsUndressed)
    (h751 : Lemma751Statement) :
    Theorem01_Statement := by
  rw [theorem01_body_iff]
  intro κ hκ
  suffices h5 : Theorem01_Body 1e5 from boussinesq_rescale κ 1e5 hκ (by norm_num) h5
  have hκ5 : (0 : ℝ) < 1e5 := by norm_num
  have hlog : |Real.log (1e5 : ℝ)| ≤ 30 := by
    have h0 : 0 ≤ Real.log (1e5 : ℝ) := Real.log_nonneg (by norm_num)
    rw [abs_of_nonneg h0]
    have h1 : Real.log (1e5 : ℝ) < Real.log ((2 : ℝ) ^ 17) := Real.log_lt_log (by norm_num) (by norm_num)
    rw [Real.log_pow] at h1
    have h2 := Real.log_two_lt_d9
    push_cast at h1
    nlinarith
  obtain ⟨Pr, htower⟩ := towerExists_discharged_aux hB hP hW hU
  obtain ⟨Cg, hCg⟩ := htower 1e5 hκ5 le_rfl hlog
  obtain ⟨Nss, h751'⟩ := @h751 Pr Cg
  obtain ⟨NStar, hNss, hNdom, hNuni, hNrate, hNforce, hNshear, hNreg, hNtail, T, hT⟩ := hCg Nss
  obtain ⟨fu, Rball, hcurl, hodd, hsupp, hSFθ, hSFu⟩ :=
    h751' NStar 1 hNdom hNuni hNrate hNforce hNshear hNreg hNtail hNss hP hB T hT
  obtain ⟨p, hsol, hclass, hclassF, hlip, hθ0, hu0⟩ :=
    designed_pair_solves NStar Cg 1 hNdom hNuni hκ5 hP hB T hT fu hcurl hodd
      ⟨Rball, fun t ht x hx => (hsupp t ht x hx).1⟩ (fun i => (hSFu i).smooth)
  have hsf := S7.staircase_facts NStar Cg T hT
  obtain ⟨hgrad, htend⟩ := S7.readout_grad NStar Cg 1 hNdom hNuni hP T hT hsf
  obtain ⟨Mθ, hMθ⟩ := S7.theta_bounded NStar Cg 1 hNdom hNuni T hT
  have hvort := S7.readout_vorticity NStar Cg 1 hNdom hNuni hNshear hP T hT
  have hTpos : 0 < T.Tstar := by rw [← hsf.t_zero]; exact hsf.t_lt 0
  have hfθc : ∀ T' < T.Tstar, ContinuousOn (fun z : ℝ × R2 => S7.Designed.fTheta T z.1 z.2) (Icc 0 T' ×ˢ univ) :=
    fun T' hT' => hSFθ.smooth.continuousOn.mono (prod_mono (Icc_subset_Ico_right hT') subset_rfl)
  have hsubT : ∀ T' < T.Tstar, Icc 0 T' ×ˢ (univ : Set R2) ⊆ Ico 0 T.Tstar ×ˢ univ :=
    fun T' hT' => prod_mono (Icc_subset_Ico_right hT') subset_rfl
  have hfuc : ∀ T' < T.Tstar, ContinuousOn (fun z : ℝ × R2 => fu z.1 z.2) (Icc 0 T' ×ˢ univ) := by
    intro T' hT'
    have h0 := (hSFu 0).smooth.continuousOn.mono (hsubT T' hT')
    have h1 := (hSFu 1).smooth.continuousOn.mono (hsubT T' hT')
    have h := (h0.smul (continuousOn_const (c := (e 0 : R2)))).add (h1.smul (continuousOn_const (c := (e 1 : R2))))
    refine h.congr fun z _ => ?_
    show fu z.1 z.2 = fu z.1 z.2 0 • e 0 + fu z.1 z.2 1 • e 1
    ext j; fin_cases j <;> simp [e]
  have hL2 : ∀ T' < T.Tstar, ∃ C : ℝ, ∀ t ∈ Icc 0 T', MemLp (S7.Designed.fTheta T t) 2 volume ∧ MemLp (fu t) 2 volume ∧
      (eLpNorm (S7.Designed.fTheta T t) 2 volume).toReal ≤ C ∧ (eLpNorm (fu t) 2 volume).toReal ≤ C := by
    intro T' hT'
    set R : ℝ := max Rball 1 with hRdef
    have hR1 : 1 ≤ R := le_max_right _ _
    have hK : IsCompact (Icc 0 T' ×ˢ Metric.closedBall (0 : R2) R) := isCompact_Icc.prod (isCompact_closedBall _ _)
    have hKsub : Icc 0 T' ×ˢ Metric.closedBall (0 : R2) R ⊆ Icc 0 T' ×ˢ (univ : Set R2) := prod_mono subset_rfl (subset_univ _)
    obtain ⟨Bθ, hBθ⟩ := hK.exists_bound_of_continuousOn ((hfθc T' hT').mono hKsub)
    obtain ⟨Bu, hBu⟩ := hK.exists_bound_of_continuousOn ((hfuc T' hT').mono hKsub)
    have hzero : ∀ t ∈ Icc 0 T', ∀ x : R2, x ∉ Metric.closedBall (0 : R2) R → fu t x = 0 ∧ S7.Designed.fTheta T t x = 0 := by
      intro t ht x hx
      rw [Metric.mem_closedBall, dist_zero_right, not_le] at hx
      exact hsupp t ⟨ht.1, ht.2.trans_lt hT'⟩ x ((le_max_left _ _).trans_lt hx)
    have hmeas : MeasurableSet (Metric.closedBall (0 : R2) R) := measurableSet_closedBall
    have hμ : volume (Metric.closedBall (0 : R2) R) ≠ ⊤ := (isCompact_closedBall _ _).measure_lt_top.ne
    set B : ℝ := max (max Bθ Bu) 0 with hBdef
    set g : R2 → ℝ := (Metric.closedBall (0 : R2) R).indicator (fun _ => B) with hgdef
    have hgLp : MemLp g 2 volume := memLp_indicator_const 2 hmeas B (Or.inr hμ)
    have hdomθ : ∀ t ∈ Icc 0 T', ∀ x, ‖S7.Designed.fTheta T t x‖ ≤ ‖g x‖ := by
      intro t ht x
      by_cases hx : x ∈ Metric.closedBall (0 : R2) R
      · rw [hgdef, indicator_of_mem hx, Real.norm_of_nonneg (le_max_right _ _)]
        exact (hBθ (t, x) ⟨ht, hx⟩).trans ((le_max_left _ _).trans (le_max_left _ _))
      · rw [(hzero t ht x hx).2, norm_zero]; exact norm_nonneg _
    have hdomu : ∀ t ∈ Icc 0 T', ∀ x, ‖fu t x‖ ≤ ‖g x‖ := by
      intro t ht x
      by_cases hx : x ∈ Metric.closedBall (0 : R2) R
      · rw [hgdef, indicator_of_mem hx, Real.norm_of_nonneg (le_max_right _ _)]
        exact (hBu (t, x) ⟨ht, hx⟩).trans ((le_max_right _ _).trans (le_max_left _ _))
      · rw [(hzero t ht x hx).1, norm_zero]; exact norm_nonneg _
    have hsliceθ : ∀ t ∈ Icc 0 T', Continuous (S7.Designed.fTheta T t) := fun t ht =>
      (hfθc T' hT').comp_continuous (Continuous.prodMk_right t) fun x => ⟨ht, mem_univ x⟩
    have hsliceu : ∀ t ∈ Icc 0 T', Continuous (fu t) := fun t ht =>
      (hfuc T' hT').comp_continuous (Continuous.prodMk_right t) fun x => ⟨ht, mem_univ x⟩
    refine ⟨(eLpNorm g 2 volume).toReal, fun t ht => ⟨?_, ?_, ?_, ?_⟩⟩
    · exact hgLp.of_le (hsliceθ t ht).aestronglyMeasurable (Filter.Eventually.of_forall (hdomθ t ht))
    · exact hgLp.of_le (hsliceu t ht).aestronglyMeasurable (Filter.Eventually.of_forall (hdomu t ht))
    · exact ENNReal.toReal_mono hgLp.eLpNorm_ne_top (eLpNorm_mono fun x => hdomθ t ht x)
    · exact ENNReal.toReal_mono hgLp.eLpNorm_ne_top (eLpNorm_mono fun x => hdomu t ht x)
  refine ⟨T.Tstar, hTpos, S7.staircaseOf T hsf, S7.Designed.theta T, p, S7.Designed.fTheta T, S7.Designed.vel T, fu,
    hsol, hclass, hclassF, hlip, hθ0, hu0, hSFθ, hSFu, hgrad, htend, ⟨Mθ, hMθ⟩, hvort, ?_⟩
  intro T' hT' hT'0 θ' u' h' h0 h0'
  exact S7.uniqueness_lipschitz hP 1e5 T' hT'0 _ θ' _ u' _ fu (hfθc T' hT') (hfuc T' hT') (hL2 T' hT')
    (hlip T' hT' hT'0) h' h0 h0'

theorem theorem01_conditional (hB : BogovskiiAnnulusR2) (hP : PotentialTheoryR2) (hW : LandingWindowsW) (hU : FlatWindowsUndressed)
    (hL1111 : ∀ Pr : Profiles, S4.StructuredShearTimeRows Pr) :
    Theorem01_Statement := by
  rw [theorem01_body_iff]
  intro κ hκ
  suffices h5 : Theorem01_Body 1e5 from boussinesq_rescale κ 1e5 hκ (by norm_num) h5
  have hκ5 : (0 : ℝ) < 1e5 := by norm_num
  have hlog : |Real.log (1e5 : ℝ)| ≤ 30 := by
    have h0 : 0 ≤ Real.log (1e5 : ℝ) := Real.log_nonneg (by norm_num)
    rw [abs_of_nonneg h0]
    have h1 : Real.log (1e5 : ℝ) < Real.log ((2 : ℝ) ^ 17) := Real.log_lt_log (by norm_num) (by norm_num)
    rw [Real.log_pow] at h1
    have h2 := Real.log_two_lt_d9
    push_cast at h1
    nlinarith
  obtain ⟨Pr, hBump⟩ := S4.profiles_exist
  obtain ⟨C₁, hbase⟩ := S6.base_exists (1e5 : ℝ) hκ5 le_rfl Pr hU hP
  obtain ⟨Cg, hE, Cε, hCε, Nrow, hbg, hR1, Ct, hCt, hTay⟩ :=
    D.pde_rows_discharge (κ := (1e5 : ℝ)) (Pr := Pr) hκ5 hP (fun β => max C₁ (8 / β + 4))
  have hCg : ∀ β : ℝ, 0 < β → β ≤ 1 / 80 → 8 / β + 4 ≤ Cg β := fun β _ _ => (le_max_right _ _).trans (hE β)
  have hC1 : C₁ ≤ Cg (1 / 80) := (le_max_left _ _).trans (hE _)
  have hSS : S4.StructuredShearTimeRows Pr := hL1111 Pr
  obtain ⟨NStar, hNdom, hNuni, hNrow, hNabs, hNrate, hNforce, hNshear, hNreg, hNtail, hNss, hNcross⟩ :
      ∃ NStar : ℝ → ℝ,
        (∀ β : ℝ, 0 < β → S6.NStar (S6.thr4Of Pr Cg 1) Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrUniform Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → Nrow β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S6.NthrAbs β (Cε β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrRates Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrForce Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrShear Pr β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrReg Pr Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S7.NthrTail Pr Cg β ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrSS Pr hSS β (Cg β) ≤ NStar β) ∧
        (∀ β : ℝ, 0 < β → S4.NthrCross hP Pr β (Cg β) ≤ NStar β) := by
    refine ⟨fun β => max (max (max (max (max (S6.NStar (S6.thr4Of Pr Cg 1) Cg β) (S4.NthrUniform Pr β (Cg β)))
        (max (Nrow β) (S6.NthrAbs β (Cε β)))) (max (max (S4.NthrRates Pr β (Cg β)) (S4.NthrForce Pr β (Cg β)))
        (max (S7.NthrShear Pr β (Cg β)) (S7.NthrReg Pr Cg β)))) (max (S7.NthrTail Pr Cg β) (S4.NthrSS Pr hSS β (Cg β)))) (S4.NthrCross hP Pr β (Cg β)),
      ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> intro β _ <;>
      simp only [le_max_iff, le_refl, true_or, or_true]
  have hCross : S6.CrossRow (Pr := Pr) Cg NStar := by
    intro p hβ hN κ' B R K Nprev hG hS hA t ht
    exact S4.nthrCross_spec hP Pr p.β (Cg p.β) p rfl ((hNcross p.β hβ).trans hN) R K Nprev hG hS hA t ht
  obtain ⟨T₁, hB1, hS1, hCB1, hSR1, hAbs1, hSJ1, hN1, -, -, hL1, hlogs1, -, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ :=
    hbase Cg Cε hC1 hCε Nrow ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩ NStar hNrow
  obtain ⟨T, hT⟩ := S6.towerExists NStar Cg Cε Nrow 1 hNdom hNuni hNrow hNabs hCross hCg hCε ⟨hbg, hR1⟩ ⟨Ct, hCt, hTay⟩
    hP hB hBump hκ5 le_rfl hlog
    ⟨T₁, hB1, hS1, hN1, hCB1, hL1, hlogs1, hSR1, hAbs1, hSJ1, hEsc1, hrun0, hdrive0, hslot0, hdress0⟩ hW
  obtain ⟨fu, Rball, hcurl, hodd, hsupp, hSFθ, hSFu⟩ :=
    S7.lemma751 NStar Cg 1 hNdom hNuni hNrate hNforce hNshear hNreg hNtail hSS hNss hP hB T hT
  obtain ⟨p, hsol, hclass, hclassF, hlip, hθ0, hu0⟩ :=
    designed_pair_solves NStar Cg 1 hNdom hNuni hκ5 hP hB T hT fu hcurl hodd
      ⟨Rball, fun t ht x hx => (hsupp t ht x hx).1⟩ (fun i => (hSFu i).smooth)
  have hsf := S7.staircase_facts NStar Cg T hT
  obtain ⟨hgrad, htend⟩ := S7.readout_grad NStar Cg 1 hNdom hNuni hP T hT hsf
  obtain ⟨Mθ, hMθ⟩ := S7.theta_bounded NStar Cg 1 hNdom hNuni T hT
  have hvort := S7.readout_vorticity NStar Cg 1 hNdom hNuni hNshear hP T hT
  have hTpos : 0 < T.Tstar := by rw [← hsf.t_zero]; exact hsf.t_lt 0
  have hfθc : ∀ T' < T.Tstar, ContinuousOn (fun z : ℝ × R2 => S7.Designed.fTheta T z.1 z.2) (Icc 0 T' ×ˢ univ) :=
    fun T' hT' => hSFθ.smooth.continuousOn.mono (prod_mono (Icc_subset_Ico_right hT') subset_rfl)
  have hsubT : ∀ T' < T.Tstar, Icc 0 T' ×ˢ (univ : Set R2) ⊆ Ico 0 T.Tstar ×ˢ univ :=
    fun T' hT' => prod_mono (Icc_subset_Ico_right hT') subset_rfl
  have hfuc : ∀ T' < T.Tstar, ContinuousOn (fun z : ℝ × R2 => fu z.1 z.2) (Icc 0 T' ×ˢ univ) := by
    intro T' hT'
    have h0 := (hSFu 0).smooth.continuousOn.mono (hsubT T' hT')
    have h1 := (hSFu 1).smooth.continuousOn.mono (hsubT T' hT')
    have h := (h0.smul (continuousOn_const (c := (e 0 : R2)))).add (h1.smul (continuousOn_const (c := (e 1 : R2))))
    refine h.congr fun z _ => ?_
    show fu z.1 z.2 = fu z.1 z.2 0 • e 0 + fu z.1 z.2 1 • e 1
    ext j; fin_cases j <;> simp [e]
  have hL2 : ∀ T' < T.Tstar, ∃ C : ℝ, ∀ t ∈ Icc 0 T', MemLp (S7.Designed.fTheta T t) 2 volume ∧ MemLp (fu t) 2 volume ∧
      (eLpNorm (S7.Designed.fTheta T t) 2 volume).toReal ≤ C ∧ (eLpNorm (fu t) 2 volume).toReal ≤ C := by
    intro T' hT'
    set R : ℝ := max Rball 1 with hRdef
    have hR1 : 1 ≤ R := le_max_right _ _
    have hK : IsCompact (Icc 0 T' ×ˢ Metric.closedBall (0 : R2) R) := isCompact_Icc.prod (isCompact_closedBall _ _)
    have hKsub : Icc 0 T' ×ˢ Metric.closedBall (0 : R2) R ⊆ Icc 0 T' ×ˢ (univ : Set R2) := prod_mono subset_rfl (subset_univ _)
    obtain ⟨Bθ, hBθ⟩ := hK.exists_bound_of_continuousOn ((hfθc T' hT').mono hKsub)
    obtain ⟨Bu, hBu⟩ := hK.exists_bound_of_continuousOn ((hfuc T' hT').mono hKsub)
    have hzero : ∀ t ∈ Icc 0 T', ∀ x : R2, x ∉ Metric.closedBall (0 : R2) R → fu t x = 0 ∧ S7.Designed.fTheta T t x = 0 := by
      intro t ht x hx
      rw [Metric.mem_closedBall, dist_zero_right, not_le] at hx
      exact hsupp t ⟨ht.1, ht.2.trans_lt hT'⟩ x ((le_max_left _ _).trans_lt hx)
    have hmeas : MeasurableSet (Metric.closedBall (0 : R2) R) := measurableSet_closedBall
    have hμ : volume (Metric.closedBall (0 : R2) R) ≠ ⊤ := (isCompact_closedBall _ _).measure_lt_top.ne
    set B : ℝ := max (max Bθ Bu) 0 with hBdef
    set g : R2 → ℝ := (Metric.closedBall (0 : R2) R).indicator (fun _ => B) with hgdef
    have hgLp : MemLp g 2 volume := memLp_indicator_const 2 hmeas B (Or.inr hμ)
    have hdomθ : ∀ t ∈ Icc 0 T', ∀ x, ‖S7.Designed.fTheta T t x‖ ≤ ‖g x‖ := by
      intro t ht x
      by_cases hx : x ∈ Metric.closedBall (0 : R2) R
      · rw [hgdef, indicator_of_mem hx, Real.norm_of_nonneg (le_max_right _ _)]
        exact (hBθ (t, x) ⟨ht, hx⟩).trans ((le_max_left _ _).trans (le_max_left _ _))
      · rw [(hzero t ht x hx).2, norm_zero]; exact norm_nonneg _
    have hdomu : ∀ t ∈ Icc 0 T', ∀ x, ‖fu t x‖ ≤ ‖g x‖ := by
      intro t ht x
      by_cases hx : x ∈ Metric.closedBall (0 : R2) R
      · rw [hgdef, indicator_of_mem hx, Real.norm_of_nonneg (le_max_right _ _)]
        exact (hBu (t, x) ⟨ht, hx⟩).trans ((le_max_right _ _).trans (le_max_left _ _))
      · rw [(hzero t ht x hx).1, norm_zero]; exact norm_nonneg _
    have hsliceθ : ∀ t ∈ Icc 0 T', Continuous (S7.Designed.fTheta T t) := fun t ht =>
      (hfθc T' hT').comp_continuous (Continuous.prodMk_right t) fun x => ⟨ht, mem_univ x⟩
    have hsliceu : ∀ t ∈ Icc 0 T', Continuous (fu t) := fun t ht =>
      (hfuc T' hT').comp_continuous (Continuous.prodMk_right t) fun x => ⟨ht, mem_univ x⟩
    refine ⟨(eLpNorm g 2 volume).toReal, fun t ht => ⟨?_, ?_, ?_, ?_⟩⟩
    · exact hgLp.of_le (hsliceθ t ht).aestronglyMeasurable (Filter.Eventually.of_forall (hdomθ t ht))
    · exact hgLp.of_le (hsliceu t ht).aestronglyMeasurable (Filter.Eventually.of_forall (hdomu t ht))
    · exact ENNReal.toReal_mono hgLp.eLpNorm_ne_top (eLpNorm_mono fun x => hdomθ t ht x)
    · exact ENNReal.toReal_mono hgLp.eLpNorm_ne_top (eLpNorm_mono fun x => hdomu t ht x)
  refine ⟨T.Tstar, hTpos, S7.staircaseOf T hsf, S7.Designed.theta T, p, S7.Designed.fTheta T, S7.Designed.vel T, fu,
    hsol, hclass, hclassF, hlip, hθ0, hu0, hSFθ, hSFu, hgrad, htend, ⟨Mθ, hMθ⟩, hvort, ?_⟩
  intro T' hT' hT'0 θ' u' h' h0 h0'
  exact S7.uniqueness_lipschitz hP 1e5 T' hT'0 _ θ' _ u' _ fu (hfθc T' hT') (hfuc T' hT') (hL2 T' hT')
    (hlip T' hT' hT'0) h' h0 h0'

end EulerBlowup

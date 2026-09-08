import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# The rotational receiver: stage time, gain, and the scale-matching exponents

This file formalises the **exponent and inequality arithmetic** of

* `strained-core/work/pass3/quantitative-core-stage.md` (with
  `verify_quantitative_core_stage.py`), §2: the stage time
  `τ = min{T_E, k₀/M₂, Ω k₀/M₃}` and the gain `g = k₀τ²/4`;
* `strained-core/work/pass3/rotational-receiver.md` (with
  `verify_rotational_receiver.py`), the scale-matching exponents `q^{-a}`, `q^{-(1+a)}`,
  `q^{1-2a}` and the identities among them;
* the advective-stage-time window recorded in `strained-core/work/pass4/iteration-audit.md`
  §3 (`T - t_m ∼ L_m^{2+a}`, `U_m ∼ (T-t_m)^{-(1+a)/(2+a)}`), whose Type II reading is
  "`(1+a)/(2+a) > 1/2` exactly when `a > 0`".

Nothing here is a PDE statement.  Read "What is NOT formalised" before quoting any of it.

## What is formalised

### Stage time and gain (§2 of `quantitative-core-stage.md`, equation (6))

* `stageTime TE M₂ M₃ Ω k₀ = min TE (min (k₀/M₂) (Ω k₀/M₃))` and
  `gain k₀ t = k₀ t²/4`.
* `stageTime_pos` : positive when `T_E, M₂, M₃, Ω, k₀` are.  This is the `τ > 0` of the
  boxed formula (6).
* `gain_pos` : `0 < g` when `0 < k₀` and `t ≠ 0`, and `one_lt_one_add_gain` : `1 < 1 + g`.
  This is the sense in which (9), `ω_z(τ,0) ≥ 2Ω(1+g)`, is a gain.
* `stageTime_mono` : `τ` is monotone nondecreasing in the pressure margin `k₀` (the other
  parameters fixed), `gain_mono` : `g` is monotone in `k₀` and in `t` on nonnegative data,
  and `gain_stageTime_mono` : hence `g(k₀, τ(k₀))` is monotone in `k₀`.  A larger certified
  margin never certifies a smaller gain.

### The contraction ratio (§5 of `quantitative-core-stage.md`, and the checker)

* `contraction_sq_mem_Ioo` : for `g > 0`, `q² = (1+g/2)/(1+g) ∈ (0,1)`, and
  `contraction_sq_mul` : `q²(1+g) = 1+g/2`.  These are the two assertions
  `0 < q_squared < 1` and `q_squared*(1+gain) == 1+gain/2` of
  `verify_quantitative_core_stage.py`.
* `strict_Re_gain_iff` : for `q > 0`, `1 < q² G ↔ 1/q² < G`, which is the boxed condition
  `G > q^{-2}` of `rotational-receiver.md`.

### The scale-matching exponents (`rotational-receiver.md`, "Exact scale matching")

For `0 < q < 1` and `0 < a < 1/2`:

| Statement | Note |
|---|---|
| `Re_ratio_gt_one : 1 < q ^ (-a)` | `Re_{qr}(τ)/Re_r(0) = q^{-a} > 1` |
| `vel_ratio_gt_one : 1 < q ^ (-(1+a))` | `U_{qr}(τ)/U_r(0) = q^{-(1+a)} > 1` |
| `energy_ratio_pos`, `energy_ratio_lt_one` | `E_{qr}(τ)/E_r(0) = q^{1-2a} ∈ (0,1)` |
| `scale_match_Re : q * q^(-(1+a)) = q^(-a)` | `q*q**(-1-a) - q**(-a)` in the checker |
| `scale_match_energy : q^(3:ℝ) * q^(-2*(1+a)) = q^(1-2a)` | `q**3*q**(-2*(1+a)) - q**(1-2*a)` |
| `vel_ratio_sq : (q^(-(1+a)))^2 = q^(-2*(1+a))` | the weight-mass step `∫χ_{qr} = q³∫χ_r` |

with the `a = 1/4` specialisations `quarter_exponents`, `quarter_Re_ratio`,
`quarter_vel_ratio`, `quarter_energy_ratio` : `q^{-1/4} > 1`, `q^{-5/4} > 1`,
`q^{1/2} ∈ (0,1)`.

`bracket_upper` : for `g > 0` and `a > 0`, the note's upper bracket endpoint
`q₊ = (1+g)^{-1/(a+2)}` is `< 1`; `bracket_lower_pos` : the lower endpoint is positive.

### The Type II / finite-energy window

With advective stage time `∼ r^{2+a}`, so that `U ∼ (T-t)^{-(1+a)/(2+a)}`:

* `advective_exponent_gt_half_iff (h : -2 < a) : 1/2 < (1+a)/(2+a) ↔ 0 < a`.  Above `a = 0`
  the velocity blows up faster than the type I rate `(T-t)^{-1/2}`, so such a scenario is
  Type II in the KNSS sense used in `theorems/02_forced_axisymmetric_on_axis_typeII.md`.
* `advective_exponent_lt_one (h : -2 < a) : (1+a)/(2+a) < 1`.
* `finite_energy_iff : 0 < 1 - 2*a ↔ a < 1/2`, the other end: the stage energy exponent
  `q^{1-2a}` is contracting exactly below `a = 1/2`.
* `window` : both at once on `0 < a < 1/2`.

### Small exact identities from `verify_rotational_receiver.py`

* `sphere_average_rotational : ½∫_{-1}^{1}(1-μ²)dμ = 2/3`, the normalisation
  `⟨|Jx|²⟩ = (2/3)⟨|x|²⟩` on the unit sphere.
* `weighted_projection` and `rigid_rotation_saturates`, the two lines behind
  `E_r(t) ≥ ½Ω_r(t)²N_r = a_r(t)²E_r(0)`.
* `cubic_remainder_constant : (3M₃/(2Ω))·t³/6 = M₃t³/(4Ω)`, the Taylor remainder constant.

## What is NOT formalised

* **No PDE, no solution, no receiver.**  `a_r(t)`, `Ω_r(t)`, `E_r`, `U_r`, `Re_r` do not
  appear.  The three "ratios" above are the real numbers `q^{-a}`, `q^{-(1+a)}`, `q^{1-2a}`;
  that they *equal* those observable ratios is the content of `rotational-receiver.md` and
  rests on local smooth existence, the full initial pressure calculation, the harmonic
  first-moment identity, and the intermediate value theorem applied to `F(q)`.  None of that
  is here.
* **The existence of the root `q` is not formalised.**  The note gets `q ∈ (0,1)` with
  `F(q) = 1` from the IVT applied to a continuous function built from the actual solution.
  Here `q` is a hypothesis.  Nothing below asserts that any such `q` exists.
* **`τ` and `g` are not shown to do anything.**  `stageTime` and `gain` are the *formulas*
  (6); that the solution satisfies (7), (8), (9) on `[0,τ]` is the note's energy and
  differentiated-equation argument, which needs the explicit Sobolev product constants
  `C₄, C₆, C₈, C_E`, the lifespan estimate, and the third-derivative bound `M₃`.  None of
  those is formalised, and neither `M₂` nor `M₃` is here anything but a positive real.
* **The Sobolev counting constants are not checked.**  The inequalities
  `9·2^{2r}·C(r+3,3) < C_r²` of `verify_quantitative_core_stage.py` are decidable rational
  arithmetic, but they are inputs to an estimate that is not formalised, so they are omitted
  rather than formalised in isolation.
* **The Type II reading is conditional in the note and conditional here.**  `T - t_m ∼
  L_m^{2+a}` holds in `iteration-audit.md` only *if* a separate stage theorem supplies the
  duration bounds `δ₋ L_m/U_m ≤ t_{m+1}-t_m ≤ δ₊ L_m/U_m` and a uniform contraction.  The
  note says so explicitly; `advective_exponent_gt_half_iff` is arithmetic about the exponent,
  not evidence for the hypothesis.
* **Nothing here is a blow-up statement.**  In particular `window` says two rational
  inequalities are simultaneously satisfiable, not that any solution realises them.
-/

open intervalIntegral MeasureTheory Set

namespace ForcedRoute
namespace Receiver

/-! ## 1. Stage time and gain -/

/-- The stage time `τ = min{T_E, k₀/M₂, Ω k₀/M₃}` of equation (6) of
`quantitative-core-stage.md`. -/
noncomputable def stageTime (TE M2 M3 Om k0 : ℝ) : ℝ :=
  min TE (min (k0 / M2) (Om * k0 / M3))

/-- The gain `g = k₀ τ²/4` of equation (6). -/
noncomputable def gain (k0 t : ℝ) : ℝ := k0 * t ^ 2 / 4

/-- `τ > 0`, the positivity asserted in the boxed formula (6). -/
theorem stageTime_pos {TE M2 M3 Om k0 : ℝ} (hTE : 0 < TE) (hM2 : 0 < M2) (hM3 : 0 < M3)
    (hOm : 0 < Om) (hk0 : 0 < k0) : 0 < stageTime TE M2 M3 Om k0 :=
  lt_min hTE (lt_min (div_pos hk0 hM2) (div_pos (mul_pos hOm hk0) hM3))

/-- `g > 0`. -/
theorem gain_pos {k0 t : ℝ} (hk0 : 0 < k0) (ht : t ≠ 0) : 0 < gain k0 t := by
  have : 0 < t ^ 2 := pow_pos (abs_pos.mpr ht) 2 |>.trans_le (le_of_eq (sq_abs t))
  exact div_pos (mul_pos hk0 this) (by norm_num)

/-- `1 + g > 1`: the growth factor of (9), `ω_z(τ,0) ≥ 2Ω(1+g)`, is genuinely bigger than
one. -/
theorem one_lt_one_add_gain {k0 t : ℝ} (hk0 : 0 < k0) (ht : t ≠ 0) : 1 < 1 + gain k0 t := by
  have := gain_pos hk0 ht; linarith

/-- **`τ` is monotone in the pressure margin `k₀`.**  Every branch of the minimum is
nondecreasing in `k₀`, so a larger certified margin never shortens the certified stage. -/
theorem stageTime_mono {TE M2 M3 Om k0 k1 : ℝ} (hM2 : 0 < M2) (hM3 : 0 < M3) (hOm : 0 ≤ Om)
    (h : k0 ≤ k1) : stageTime TE M2 M3 Om k0 ≤ stageTime TE M2 M3 Om k1 := by
  refine min_le_min (le_refl TE) (min_le_min ?_ ?_)
  · gcongr
  · gcongr

/-- `g` is monotone in `k₀` and in `t`, on nonnegative data. -/
theorem gain_mono {k0 k1 t u : ℝ} (hk0 : 0 ≤ k0) (hk : k0 ≤ k1) (ht : 0 ≤ t) (htu : t ≤ u) :
    gain k0 t ≤ gain k1 u := by
  have h1 : t ^ 2 ≤ u ^ 2 := by nlinarith
  have h2 : (0:ℝ) ≤ t ^ 2 := sq_nonneg t
  have : k0 * t ^ 2 ≤ k1 * u ^ 2 := by nlinarith
  unfold gain; linarith

/-- Composing the two: the certified gain `g(k₀, τ(k₀))` is monotone in the margin `k₀`. -/
theorem gain_stageTime_mono {TE M2 M3 Om k0 k1 : ℝ} (hTE : 0 < TE) (hM2 : 0 < M2)
    (hM3 : 0 < M3) (hOm : 0 < Om) (hk0 : 0 < k0) (h : k0 ≤ k1) :
    gain k0 (stageTime TE M2 M3 Om k0) ≤ gain k1 (stageTime TE M2 M3 Om k1) :=
  gain_mono hk0.le h (stageTime_pos hTE hM2 hM3 hOm hk0).le
    (stageTime_mono hM2 hM3 hOm.le h)

/-! ## 2. The contraction ratio -/

/-- `q² = (1+g/2)/(1+g) ∈ (0,1)` for `g > 0`.  This is the assertion
`0 < q_squared < 1` of `verify_quantitative_core_stage.py`. -/
theorem contraction_sq_mem_Ioo {g : ℝ} (hg : 0 < g) : (1 + g/2) / (1 + g) ∈ Ioo (0:ℝ) 1 := by
  have h1 : (0:ℝ) < 1 + g := by linarith
  constructor
  · exact div_pos (by linarith) h1
  · rw [div_lt_one h1]; linarith

/-- `q²(1+g) = 1+g/2`, the checker's `q_squared*(1+gain) == 1+gain/2`. -/
theorem contraction_sq_mul {g : ℝ} (hg : 0 < g) : ((1 + g/2) / (1 + g)) * (1 + g) = 1 + g/2 := by
  have h1 : (1:ℝ) + g ≠ 0 := by positivity
  field_simp

/-- **The strict Reynolds-gain condition.**  Since the note's lower bound is
`Re_{qr}(t)/Re_r(0) ≥ q² a_{qr}(t) ≥ q² G`, a lower bound `G` gives a strict gain exactly
when `q² G > 1`, i.e. when `G > q^{-2}` — the boxed condition of `rotational-receiver.md`. -/
theorem strict_Re_gain_iff {q G : ℝ} (hq : 0 < q) : 1 < q ^ 2 * G ↔ 1 / q ^ 2 < G := by
  have h : (0:ℝ) < q ^ 2 := pow_pos hq 2
  rw [div_lt_iff₀ h]
  constructor <;> intro h' <;> nlinarith

/-! ## 3. The scale-matching exponents -/

/-- `Re_{qr}(τ)/Re_r(0) = q^{-a} > 1` for `0 < q < 1` and `a > 0`. -/
theorem Re_ratio_gt_one {q a : ℝ} (hq0 : 0 < q) (hq1 : q < 1) (ha : 0 < a) :
    1 < q ^ (-a) := by
  rw [Real.one_lt_rpow_iff_of_pos hq0]
  exact Or.inr ⟨hq1, by linarith⟩

/-- `U_{qr}(τ)/U_r(0) = q^{-(1+a)} > 1` for `0 < q < 1` and `a > 0`. -/
theorem vel_ratio_gt_one {q a : ℝ} (hq0 : 0 < q) (hq1 : q < 1) (ha : 0 < a) :
    1 < q ^ (-(1 + a)) := by
  rw [Real.one_lt_rpow_iff_of_pos hq0]
  exact Or.inr ⟨hq1, by linarith⟩

/-- `E_{qr}(τ)/E_r(0) = q^{1-2a} > 0`. -/
theorem energy_ratio_pos {q a : ℝ} (hq0 : 0 < q) : 0 < q ^ (1 - 2*a) :=
  Real.rpow_pos_of_pos hq0 _

/-- `E_{qr}(τ)/E_r(0) = q^{1-2a} < 1` for `0 < q < 1` and `a < 1/2`: the energy in the
smaller receiver is a strictly smaller multiple of the original.  This is why the note calls
the window `a < 1/2` the finite-energy end. -/
theorem energy_ratio_lt_one {q a : ℝ} (hq0 : 0 < q) (hq1 : q < 1) (ha : a < 1/2) :
    q ^ (1 - 2*a) < 1 := by
  rw [Real.rpow_lt_one_iff_of_pos hq0]
  exact Or.inr ⟨hq1, by linarith⟩

/-- **The RMS/Reynolds scale match.**  `q · (velocity ratio) = Reynolds ratio`; this is the
checker's `q*q**(-1-a) - q**(-a)`. -/
theorem scale_match_Re {q : ℝ} (hq : 0 < q) (a : ℝ) : q * q ^ (-(1 + a)) = q ^ (-a) := by
  have h : q ^ (1:ℝ) * q ^ (-(1 + a)) = q ^ ((1:ℝ) + -(1 + a)) := (Real.rpow_add hq _ _).symm
  rw [Real.rpow_one] at h
  rw [h]; congr 1; ring

/-- **The energy scale match.**  Weight mass `q³` times the square of the velocity ratio is
the energy ratio; this is the checker's `q**3*q**(-2*(1+a)) - q**(1-2*a)`. -/
theorem scale_match_energy {q : ℝ} (hq : 0 < q) (a : ℝ) :
    q ^ (3:ℝ) * q ^ (-2*(1 + a)) = q ^ (1 - 2*a) := by
  rw [← Real.rpow_add hq]; congr 1; ring

/-- The square of the velocity ratio, as an ordinary square. -/
theorem vel_ratio_sq {q : ℝ} (hq : 0 < q) (a : ℝ) :
    (q ^ (-(1 + a))) ^ (2:ℕ) = q ^ (-2*(1 + a)) := by
  rw [← Real.rpow_natCast (q ^ (-(1 + a))) 2, ← Real.rpow_mul hq.le]
  congr 1; push_cast; ring

/-- The energy scale match in the note's own shape:
`E ratio = q³ · (U ratio)²`. -/
theorem scale_match_energy' {q : ℝ} (hq : 0 < q) (a : ℝ) :
    q ^ (3:ℝ) * (q ^ (-(1 + a))) ^ (2:ℕ) = q ^ (1 - 2*a) := by
  rw [vel_ratio_sq hq, scale_match_energy hq]

/-! ### The `a = 1/4` specialisations -/

/-- At `a = 1/4` the velocity and energy exponents are `-5/4` and `1/2` (the Reynolds
exponent `-a` is `-1/4` by inspection). -/
theorem quarter_exponents : -(1 + (1/4 : ℝ)) = -(5/4) ∧ 1 - 2*(1/4 : ℝ) = 1/2 := by norm_num

/-- `q^{-1/4} > 1`. -/
theorem quarter_Re_ratio {q : ℝ} (hq0 : 0 < q) (hq1 : q < 1) : 1 < q ^ (-(1/4) : ℝ) :=
  Re_ratio_gt_one hq0 hq1 (by norm_num)

/-- `q^{-5/4} > 1`. -/
theorem quarter_vel_ratio {q : ℝ} (hq0 : 0 < q) (hq1 : q < 1) : 1 < q ^ (-(5/4) : ℝ) := by
  have := vel_ratio_gt_one (a := 1/4) hq0 hq1 (by norm_num)
  norm_num at this ⊢
  exact this

/-- `q^{1/2} ∈ (0,1)`. -/
theorem quarter_energy_ratio {q : ℝ} (hq0 : 0 < q) (hq1 : q < 1) :
    0 < q ^ ((1/2) : ℝ) ∧ q ^ ((1/2) : ℝ) < 1 := by
  refine ⟨Real.rpow_pos_of_pos hq0 _, ?_⟩
  have := energy_ratio_lt_one (a := 1/4) hq0 hq1 (by norm_num)
  norm_num at this ⊢
  exact this

/-! ### The contraction bracket -/

/-- The upper bracket endpoint `q₊ = (1+g)^{-1/(a+2)}` is `< 1` for `g > 0`, `a > 0`. -/
theorem bracket_upper {g a : ℝ} (hg : 0 < g) (ha : 0 < a) :
    (1 + g) ^ (-(1/(a + 2))) < 1 := by
  have h1 : (0:ℝ) < 1 + g := by linarith
  have ha2 : (0:ℝ) < a + 2 := by linarith
  rw [Real.rpow_lt_one_iff_of_pos h1]
  exact Or.inl ⟨by linarith, neg_neg_iff_pos.mpr (by positivity)⟩

/-- The lower bracket endpoint is positive. -/
theorem bracket_lower_pos {x a : ℝ} (hx : 0 < x) : 0 < x ^ (1/(a + 2)) :=
  Real.rpow_pos_of_pos hx _

/-! ## 4. The Type II / finite-energy window -/

/-- **The Type II end of the window.**  With advective stage time `∼ r^{2+a}`, the note's
`U_m ∼ (T-t_m)^{-(1+a)/(2+a)}` beats the type I rate `(T-t)^{-1/2}` exactly when `a > 0`. -/
theorem advective_exponent_gt_half_iff {a : ℝ} (ha : -2 < a) :
    1/2 < (1 + a)/(2 + a) ↔ 0 < a := by
  have h2 : (0:ℝ) < 2 + a := by linarith
  rw [lt_div_iff₀ h2]
  constructor <;> intro h <;> linarith

/-- The same exponent is always below `1`; the note's `‖∇u(t_m)‖_∞ ≳ 1/(T-t_m)` is a
separate statement. -/
theorem advective_exponent_lt_one {a : ℝ} (ha : -2 < a) : (1 + a)/(2 + a) < 1 := by
  have h2 : (0:ℝ) < 2 + a := by linarith
  rw [div_lt_one h2]; linarith

/-- **The finite-energy end of the window.**  The stage energy exponent `1 - 2a` is positive
exactly below `a = 1/2`. -/
theorem finite_energy_iff {a : ℝ} : 0 < 1 - 2*a ↔ a < 1/2 := by
  constructor <;> intro h <;> linarith

/-- **The window.**  On `0 < a < 1/2` both ends hold: the velocity exponent exceeds `1/2`
(Type II) and the energy exponent is positive (contracting energy). -/
theorem window {a : ℝ} (h0 : 0 < a) (h1 : a < 1/2) :
    1/2 < (1 + a)/(2 + a) ∧ 0 < 1 - 2*a :=
  ⟨(advective_exponent_gt_half_iff (by linarith)).mpr h0, finite_energy_iff.mpr h1⟩

/-- At `a = 1/4` the window is non-empty, with velocity exponent `5/9`. -/
theorem window_quarter : (1 + (1/4:ℝ))/(2 + 1/4) = 5/9 ∧ 1/2 < (5/9:ℝ) ∧ 0 < 1 - 2*(1/4:ℝ) := by
  norm_num

/-! ## 5. Small exact identities from `verify_rotational_receiver.py` -/

/-- `⟨|Jx|²⟩ = (2/3)⟨|x|²⟩` on the unit sphere, in the form `½∫_{-1}^{1}(1-μ²)dμ = 2/3`.
This is the checker's `N_sphere - 2/3`. -/
theorem sphere_average_rotational :
    (1/2) * (∫ mu in (-1:ℝ)..1, (1 - mu ^ 2)) = 2/3 := by
  have key : ∀ x : ℝ, HasDerivAt (fun t : ℝ => t - t ^ 3 / 3) (1 - x ^ 2) x := fun x =>
    ((hasDerivAt_id x).sub ((hasDerivAt_pow 3 x).div_const 3)).congr_deriv (by norm_num)
  rw [integral_eq_sub_of_hasDerivAt (fun x _ => key x)
    ((by fun_prop : Continuous fun x : ℝ => 1 - x ^ 2).intervalIntegrable _ _)]
  norm_num

/-- The weighted orthogonal decomposition behind `E_r(t) ≥ ½Ω_r(t)²N_r`: substituting
`P = ω_r N` into `2E - 2ω_r P + ω_r² N` gives `2E - ω_r² N`. -/
theorem weighted_projection (E N om : ℝ) :
    2*E - 2*om*(om*N) + om ^ 2 * N = 2*E - om ^ 2 * N := by ring

/-- The initial rigid rotation saturates the weighted energy normalisation:
with `ω_r = Ω a` and `E₀ = Ω²N/2`, one has `ω_r²N/2 = a²E₀`. -/
theorem rigid_rotation_saturates (Om N a : ℝ) :
    (Om * a) ^ 2 * N / 2 = a ^ 2 * (Om ^ 2 * N / 2) := by ring

/-- The Taylor remainder constant: a uniform third derivative `|a_r'''| ≤ 3M₃/(2Ω)` gives the
recorded cubic remainder `M₃t³/(4Ω)`.  This is the checker's last `zero(...)` line. -/
theorem cubic_remainder_constant {M3 Om t : ℝ} (hOm : Om ≠ 0) :
    (3/2 * M3 / Om) * t ^ 3 / 6 = M3 * t ^ 3 / (4 * Om) := by
  field_simp; ring

end Receiver
end ForcedRoute

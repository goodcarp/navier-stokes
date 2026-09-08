import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.Notation

/-!
# The strained rotating core: affine algebra, the radial cancellations, and the ratio test

This file formalises the **elementary, exact** parts of the strained-core sitting:

* `strained-core/work/pass4/affine-core-pressure.md` (with
  `derive_affine_core_pressure.py` and `verify_affine_pressure_independent.py`), the affine
  generator `L = diag(-b,-b,2b) + Ω J` and the profile-independence of the central pressure
  Hessian of its compact divergence-free extension;
* `strained-core/work/pass2/active-core-pressure.md` (with
  `verify_active_core_pressure.py`), the rotating-core coefficient `2Ω²/5`;
* `strained-core/work/pass4/core-pressure-derivative.md` (with
  `verify_affine_pressure_independent.py`), the ratio test at the neutral insertion.

Nothing here is a PDE statement.  Read "What is NOT formalised" before quoting any of it.

## What is formalised

### The affine generator (§1 of `affine-core-pressure.md`)

* `Lmat b Ω = !![-b, -Ω, 0; Ω, -b, 0; 0, 0, 2*b]`, the matrix `diag(-b,-b,2b) + Ω J`.
* `trace_Lmat : (Lmat b Ω).trace = 0`, the hypothesis under which
  `curl [x × (Lx)] = -3 L x` in (1)-(2) of the note.
* `Lmat_sq_zz : (Lmat b Ω * Lmat b Ω) 2 2 = 4 * b ^ 2`, and
  `trace_Lmat_sq : (Lmat b Ω * Lmat b Ω).trace = 6 * b ^ 2 - 2 * Ω ^ 2`, which is the value
  `g(0) = tr (L²)` of the pressure source at the origin in §2.
* `affine_pressure_zz : -(Lmat b Ω * Lmat b Ω) 2 2 = -(4 * b ^ 2)`, the axial pressure
  Hessian entry of the **homogeneous** affine field `u = L x` on all of `ℝ³` (for which the
  pressure is the quadratic `p(x) = -x ⬝ L² x / 2`, so `D²p = -sym(L²)`).  It is independent
  of `Ω`.
* `affine_ne_core : (b, Ω) ≠ (0,0) → -(4*b^2) ≠ -(18/7)*b^2 + (2/5)*Ω^2`.  The homogeneous
  and the localised values are *different numbers*: the compact extension of §1 is not a
  harmless truncation of `L x`, and the note's `-(18/7) b² + (2/5) Ω²` is a statement about
  the extension, not about `L`.

### The profile moments and the two radial cancellations (§3 of `affine-core-pressure.md`)

`Cutoff` bundles the note's radial profile `ψ` with its first two derivatives, a radius `T`
past which `ψ` and `ψ'` vanish, and `ψ(0) = 1`.  With `I = ∫₀^T s ψ'(s)² ds`:

| Statement | Note |
|---|---|
| `Cutoff.moment_one : ∫₀^T ψ ψ' = -1/2` | (7), first line |
| `Cutoff.moment_two : ∫₀^T s ψ ψ'' = 1/2 - I` | (7), second line |
| `Cutoff.moment_three : ∫₀^T s² ψ' ψ'' = -I` | (7), third line |
| `Cutoff.strain_radial : ½∫ 𝒜_S ds/s = -4/7` | (8) |
| `Cutoff.swirl_radial : ½∫ 𝒜_W ds/s = -4/15` | (9) |

The three moments are proved from the fundamental theorem of calculus (with the
antiderivatives `ψ²/2`, `s ψ ψ'` and `s²(ψ')²/2`); they are *not* assumed.  The two radial
identities are then genuine theorems for **every** profile in the class: the coefficient of
`I` cancels, which is the content the note calls "all radial-profile dependence cancels".

`Cutoff.core_pressure_zz` and `Cutoff.core_pressure_transverse` assemble (4), (8), (9) into
(10) and (11):

    p_zz  = -(1/3)·tr(L²) + b²·(½∫𝒜_S ds/s) + Ω²·(½∫𝒜_W ds/s) = -(18/7) b² + (2/5) Ω²,
    p_xx  = p_yy = (-tr(L²) - p_zz)/2                          = -(12/7) b² + (4/5) Ω².

### The rotating core of `active-core-pressure.md` §1

* `angular_average : ½∫_{-1}^{1} (3μ²-1)(1-μ²) dμ = -4/15`, the note's `ang`, an exact
  polynomial integral over an interval.
* `Cutoff.radial_moment : ∫₀^{√T} r ψ(r²) ψ'(r²) dr = -1/4`, the note's `rad`, obtained from
  `moment_one` by the substitution `s = r²`.
* `rotating_core_coefficient : (2:ℝ)/3 - 4 * (-4/15) * (-1/4) = 2/5`, the note's line
  `∂_zz p[c](0) = 2/3 Ω² - 4/15 Ω² = 2/5 Ω²`, and
  `Cutoff.core_pressure_zz_of_pure_rotation`, the same number reached from the `b = 0` case
  of the strained calculation.

### The ratio test (`core-pressure-derivative.md`, `verify_affine_pressure_independent.py`)

* `ratio_second_derivative_algebraic` : the chain-rule sum of the checker's `deriv.subs`
  line, as an identity in `b, Ω, p_zz'`.
* `ratio_second_derivative` : the same with real derivatives.  Given differentiable
  `b, Ω, P, η_b, η_w` at `0` and the note's first-derivative relation

      β'(t) = (-4 b(t)² - P(t)/2 + η_b(t) - (b(t)/Ω(t)) η_w(t)) / Ω(t)

  holding as an identity of functions, together with the neutral insertion
  `b'(0) = 2b(0)²`, `Ω'(0) = 2b(0)Ω(0)`, `P(0) = -8b(0)²` and `η_b = η_w = η_b' = η_w' = 0`
  at `0`, the second derivative of `β` at `0` is `-(P'(0) + 32 b(0)³)/(2 Ω(0))`.

  The relation for `β'` is a **hypothesis**, transcribed from the note; see below.

## What is NOT formalised

This is the important section.  What is below is one matrix computation, three
one-dimensional integrals with their exact linear combinations, and one chain rule.  None of
it is a theorem about Navier-Stokes.

* **The pressure solve.**  `p = N * g` with `N(x) = 1/(4π|x|)`, the distributional Hessian
  identity `∂_zz N = PV[(3z² - |x|²)/(4π|x|⁵)] - δ₀/3` and its contact term, the convergence
  of the principal value, and the decay condition fixing the constants are all **assumed**.
  The factor `-(1/3)·tr(L²)` and the shape `½∫⟨(3μ²-1)g⟩ ds/s` of formula (4) enter
  `core_pressure_zz` as the *definition of what is being added up*, not as derived facts.
* **The angular average.**  Equations (5)-(6), i.e. that `⟨(3μ²-1) g⟩ = b² 𝒜_S + Ω² 𝒜_W`
  with `𝒜_S`, `𝒜_W` the displayed expressions in `ψ, ψ', ψ''`, is the symbolic Cartesian
  differentiation performed by `derive_affine_core_pressure.py` (and re-derived
  cylindrically by `verify_affine_pressure_independent.py`).  It is **not** verified here.
  What is verified is that *those* integrands integrate to `-4/7` and `-4/15`.
* **The localisation itself.**  The extension `u_L = curl[-(ψ/3) x × (Lx)]` of (1), its
  divergence-freeness, formula (2), the poloidal/azimuthal split (3), and the absence of a
  `bΩ` cross term in `g` are not formalised.  No vector field appears in this file at all.
* **The vanishing of the `I`-coefficient is proved; its interpretation is not.**
  `strain_radial` and `swirl_radial` hold for every `Cutoff`, but "the coefficients do not
  depend on the radial transition width or support radius" is a statement about a family of
  physical fields, and no such family is constructed here.
* **Convergence at infinity is replaced by a compact endpoint.**  The note writes
  `∫₀^∞`; `Cutoff` carries a radius `T` with `ψ T = ψ' T = 0`, and every integral is over
  `[0, T]`.  Since `ψ` is compactly supported in the note this is the same number, but the
  improper integral is not formalised.
* **`ψ ∈ C_c^∞` is weakened to `C²` with compact endpoint data.**  `Cutoff` asks only for
  two derivatives and continuity of the second.  Nothing in §3 needs more; nothing in §3 is
  claimed to need more.
* **The ratio test's `β'` relation is a hypothesis.**  `ratio_second_derivative` assumes the
  displayed formula for `β'(t)`.  That formula comes from the evolution equations (16) of
  `affine-core-pressure.md` together with the profile-defect terms `η_b, η_w` of
  `full-feedback-gate-audit.md`; deriving it needs the PDE, which is absent here.  What is
  proved is that *given* it, the chain rule produces exactly the note's
  `β''(0) = -(p_zz'(0) + 32b³)/(2Ω)` — including that the `4b²/Ω² + p_zz/(2Ω²)` term is
  annihilated by `p_zz = -8b²` at the neutral point.
* **The numerical lead is untouched.**  `t₃₀₀ ≈ -6.7344633112` and `t_Wb ≈ 2.563122131694`
  of §5 of `core-pressure-derivative.md` are floating-point quadrature for one specific
  logistic cutoff.  They need `Q(s) = ∫₀^s ψ`, `K(s) = ∫₀^s τ^{7/2}[ψψ' + (7/9)τ(ψ')²]`,
  half-integer powers and cumulative integrals of a transcendental profile.  They are **not**
  rational, they are not profile-independent (the note itself exhibits two polynomial proxies
  with different values), and nothing about them is formalised here.  Neither is the claim
  `d_S = d_W = 0`, nor the exact coefficient formulas (1)-(2), nor the harmonic expansion (3).
* **No sign is asserted.**  `ratio_second_derivative` computes `β''(0)` as a formula in
  `P'(0)`; it does not evaluate or bound `P'(0)`, so it decides nothing about the cone.
-/

open intervalIntegral MeasureTheory

namespace ForcedRoute
namespace StrainedCore

/-! ## 1. The affine generator `L = diag(-b,-b,2b) + Ω J` -/

/-- The affine generator `L = diag(-b, -b, 2b) + Ω J` of §1 of `affine-core-pressure.md`,
`J` being the generator of rotation in the `x`-`y` plane. -/
def Lmat (b Om : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  !![-b, -Om, 0; Om, -b, 0; 0, 0, 2*b]

/-- `tr L = 0`.  This is the hypothesis under which `curl[x × (Lx)] = -3 L x`, i.e. under
which the vector potential of (1) reproduces `L x` near the origin. -/
theorem trace_Lmat (b Om : ℝ) : (Lmat b Om).trace = 0 := by
  simp [Lmat, Matrix.trace_fin_three]
  ring

/-- `(L²)_zz = 4b²`, independent of `Ω`: the swirl does not reach the axial diagonal entry of
`L²`. -/
theorem Lmat_sq_zz (b Om : ℝ) : (Lmat b Om * Lmat b Om) 2 2 = 4 * b ^ 2 := by
  simp [Lmat]
  ring

/-- `tr(L²) = 6b² - 2Ω²`.  This is `g(0)`, the value at the origin of the pressure source
`g = ∂_i u_j ∂_j u_i = tr[(∇u)²]` of §2 of `affine-core-pressure.md`. -/
theorem trace_Lmat_sq (b Om : ℝ) : (Lmat b Om * Lmat b Om).trace = 6 * b ^ 2 - 2 * Om ^ 2 := by
  simp [Lmat, Matrix.trace_fin_three]
  ring

/-- **The pure affine pressure Hessian entry.**  For the homogeneous field `u = L x` on all of
`ℝ³` with `tr L = 0`, the decaying-modulo-affine pressure is the quadratic
`p(x) = -x ⬝ L² x / 2`, so `D²p = -sym(L²)` and `(D²p)_zz = -(L²)_zz = -4b²`.

This is **not** the note's core value; see `affine_ne_core`. -/
theorem affine_pressure_zz (b Om : ℝ) : -(Lmat b Om * Lmat b Om) 2 2 = -(4 * b ^ 2) := by
  rw [Lmat_sq_zz]

/-- **The homogeneous and localised axial pressures are different numbers.**  For any
`(b, Ω) ≠ (0,0)`, the pure affine value `-4b²` differs from the note's compact-extension
value `-(18/7)b² + (2/5)Ω²`: the first is `≤ 0` and the second `≥ 0` exactly when the other
is not, and they meet only at the origin.  The localisation of §1 is therefore not a
harmless truncation of `L x`, and `-(18/7)b² + (2/5)Ω²` is a **profile-class constant of the
extension**, established in `core_pressure_zz` below from the extension's own angular
average, not a property of `L`. -/
theorem affine_ne_core {b Om : ℝ} (h : b ≠ 0 ∨ Om ≠ 0) :
    -(4 * b ^ 2) ≠ -(18/7) * b ^ 2 + (2/5) * Om ^ 2 := by
  intro hEq
  have hb2 : (0:ℝ) ≤ b ^ 2 := sq_nonneg b
  have hO2 : (0:ℝ) ≤ Om ^ 2 := sq_nonneg Om
  have hb : b ^ 2 = 0 := by nlinarith
  have hO : Om ^ 2 = 0 := by nlinarith
  rcases h with h | h
  · exact h (by nlinarith [sq_abs b])
  · exact h (by nlinarith [sq_abs Om])

/-! ## 2. The radial profile and its three moments -/

/-- A compact radial cutoff profile, as in §1 of `affine-core-pressure.md`: `ψ` is a function
of `s = |x|²`, equal to one at the origin and vanishing together with its first derivative
from the radius `T = R²` outwards.

The note's class is `C_c^∞([0,∞))`; only `C²` and the two endpoint values are used in §3, so
that is all that is asked for here.  `psi1`, `psi2` are `ψ'` and `ψ''`. -/
structure Cutoff where
  /-- The radial profile `ψ`, a function of `s = |x|²`. -/
  psi : ℝ → ℝ
  /-- `ψ'`. -/
  psi1 : ℝ → ℝ
  /-- `ψ''`. -/
  psi2 : ℝ → ℝ
  /-- The cutoff radius `T = R²` in the variable `s`. -/
  T : ℝ
  /-- The cutoff radius is positive. -/
  hT : 0 < T
  /-- `psi1` is the derivative of `psi`. -/
  hd1 : ∀ s, HasDerivAt psi (psi1 s) s
  /-- `psi2` is the derivative of `psi1`. -/
  hd2 : ∀ s, HasDerivAt psi1 (psi2 s) s
  /-- `ψ''` is continuous, so every integrand below is. -/
  hcont2 : Continuous psi2
  /-- `ψ` is one near the origin; only the value at `0` is used. -/
  hzero : psi 0 = 1
  /-- `ψ` vanishes from the cutoff radius outwards. -/
  hendA : psi T = 0
  /-- so does `ψ'`. -/
  hendB : psi1 T = 0

namespace Cutoff

variable (c : Cutoff)

theorem cont0 : Continuous c.psi :=
  Differentiable.continuous fun s => (c.hd1 s).differentiableAt

theorem cont1 : Continuous c.psi1 :=
  Differentiable.continuous fun s => (c.hd2 s).differentiableAt

theorem int1 : IntervalIntegrable (fun s : ℝ => s ^ 2 * c.psi1 s * c.psi2 s) volume 0 c.T :=
  (((continuous_pow 2).mul c.cont1).mul c.hcont2).intervalIntegrable _ _

theorem int2 : IntervalIntegrable (fun s : ℝ => s * c.psi s * c.psi2 s) volume 0 c.T :=
  ((continuous_id.mul c.cont0).mul c.hcont2).intervalIntegrable _ _

theorem int3 : IntervalIntegrable (fun s : ℝ => s * c.psi1 s ^ 2) volume 0 c.T :=
  (continuous_id.mul (c.cont1.pow 2)).intervalIntegrable _ _

theorem int4 : IntervalIntegrable (fun s : ℝ => c.psi s * c.psi1 s) volume 0 c.T :=
  (c.cont0.mul c.cont1).intervalIntegrable _ _

/-- `I = ∫₀^T s ψ'(s)² ds`, the one profile integral that survives the substitutions in §3 of
`affine-core-pressure.md` — and then cancels. -/
noncomputable def I : ℝ := ∫ s in (0:ℝ)..c.T, s * c.psi1 s ^ 2

theorem I_def : c.I = ∫ s in (0:ℝ)..c.T, s * c.psi1 s ^ 2 := rfl

/-- `∫₀^T ψ ψ' ds = -1/2`, the first endpoint identity of (7).  Antiderivative `ψ²/2`, with
`ψ(0) = 1` and `ψ(T) = 0`. -/
theorem moment_one : ∫ s in (0:ℝ)..c.T, c.psi s * c.psi1 s = -(1/2) := by
  have key : ∀ s : ℝ, HasDerivAt (fun t => c.psi t * c.psi t / 2) (c.psi s * c.psi1 s) s :=
    fun s => (((c.hd1 s).mul (c.hd1 s)).div_const 2).congr_deriv (by ring)
  rw [integral_eq_sub_of_hasDerivAt (fun s _ => key s) c.int4, c.hendA, c.hzero]
  norm_num

/-- `∫₀^T s² ψ' ψ'' ds = -I`, the third endpoint identity of (7).  Antiderivative
`s²(ψ')²/2`, which vanishes at both endpoints. -/
theorem moment_three : ∫ s in (0:ℝ)..c.T, s ^ 2 * c.psi1 s * c.psi2 s = -c.I := by
  have key : ∀ s : ℝ, HasDerivAt (fun t : ℝ => t ^ 2 * (c.psi1 t * c.psi1 t) / 2)
      (s * c.psi1 s ^ 2 + s ^ 2 * c.psi1 s * c.psi2 s) s := fun s =>
    (((hasDerivAt_pow 2 s).mul ((c.hd2 s).mul (c.hd2 s))).div_const 2).congr_deriv
      (by simp only [Pi.mul_apply]; ring)
  have hz : (∫ s in (0:ℝ)..c.T, (s * c.psi1 s ^ 2 + s ^ 2 * c.psi1 s * c.psi2 s)) = 0 := by
    rw [integral_eq_sub_of_hasDerivAt (fun s _ => key s) (c.int3.add c.int1), c.hendB]
    norm_num
  rw [integral_add c.int3 c.int1] at hz
  rw [c.I_def]
  linarith

/-- `∫₀^T s ψ ψ'' ds = 1/2 - I`, the second endpoint identity of (7).  Antiderivative
`s ψ ψ'`, which vanishes at both endpoints; combined with `moment_one`. -/
theorem moment_two : ∫ s in (0:ℝ)..c.T, s * c.psi s * c.psi2 s = 1/2 - c.I := by
  have key : ∀ s : ℝ, HasDerivAt (fun t : ℝ => t * (c.psi t * c.psi1 t))
      (c.psi s * c.psi1 s + s * c.psi1 s ^ 2 + s * c.psi s * c.psi2 s) s := fun s =>
    ((hasDerivAt_id s).mul ((c.hd1 s).mul (c.hd2 s))).congr_deriv
      (by simp only [Pi.mul_apply, id_eq]; ring)
  have hz : (∫ s in (0:ℝ)..c.T,
      (c.psi s * c.psi1 s + s * c.psi1 s ^ 2 + s * c.psi s * c.psi2 s)) = 0 := by
    rw [integral_eq_sub_of_hasDerivAt (fun s _ => key s) ((c.int4.add c.int3).add c.int2),
      c.hendA, c.hendB]
    norm_num
  rw [integral_add (c.int4.add c.int3) c.int2, integral_add c.int4 c.int3, c.moment_one] at hz
  rw [c.I_def]
  linarith

/-- The four-term linear combination of the moments, which is the only way (7) is used in
(8) and (9). -/
theorem lincomb (a b d e : ℝ) :
    (∫ s in (0:ℝ)..c.T, (a * (s ^ 2 * c.psi1 s * c.psi2 s) + b * (s * c.psi s * c.psi2 s)
        + d * (s * c.psi1 s ^ 2) + e * (c.psi s * c.psi1 s)))
      = a * (-c.I) + b * (1/2 - c.I) + d * c.I + e * (-(1/2)) := by
  have j1 := c.int1.const_mul a
  have j2 := c.int2.const_mul b
  have j3 := c.int3.const_mul d
  have j4 := c.int4.const_mul e
  rw [integral_add ((j1.add j2).add j3) j4, integral_add (j1.add j2) j3,
    integral_add j1 j2, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul, c.moment_one, c.moment_two, c.moment_three,
    ← c.I_def]

/-! ## 3. The two radial cancellations -/

/-- **The strain radial cancellation**, equation (8) of `affine-core-pressure.md`:

    ½ ∫₀^∞ 𝒜_S(s) ds/s = -4/7,   𝒜_S(s) = (16 s/105)[-4s²ψ'ψ'' + 6sψψ'' + 2s(ψ')² + 21ψψ'].

The coefficient of the profile integral `I` cancels exactly; the value is the same for every
`Cutoff`. -/
theorem strain_radial :
    (1/2) * (∫ s in (0:ℝ)..c.T, (16/105) * (-4 * (s ^ 2 * c.psi1 s * c.psi2 s)
        + 6 * (s * c.psi s * c.psi2 s) + 2 * (s * c.psi1 s ^ 2)
        + 21 * (c.psi s * c.psi1 s))) = -(4/7) := by
  have h : (∫ s in (0:ℝ)..c.T, (16/105) * (-4 * (s ^ 2 * c.psi1 s * c.psi2 s)
        + 6 * (s * c.psi s * c.psi2 s) + 2 * (s * c.psi1 s ^ 2)
        + 21 * (c.psi s * c.psi1 s)))
      = ∫ s in (0:ℝ)..c.T, ((-64/105) * (s ^ 2 * c.psi1 s * c.psi2 s)
        + (96/105) * (s * c.psi s * c.psi2 s) + (32/105) * (s * c.psi1 s ^ 2)
        + (336/105) * (c.psi s * c.psi1 s)) := by
    apply integral_congr; intro s _; ring
  rw [h, c.lincomb]; ring

/-- **The swirl radial cancellation**, equation (9) of `affine-core-pressure.md`:

    ½ ∫₀^∞ 𝒜_W(s) ds/s = -4/15,   𝒜_W(s) = (16 s/135)(2sψ'' + 5ψ')(2sψ' + 3ψ).

Again the coefficient of `I` cancels exactly. -/
theorem swirl_radial :
    (1/2) * (∫ s in (0:ℝ)..c.T, (16/135) * (2 * s * c.psi2 s + 5 * c.psi1 s)
        * (2 * s * c.psi1 s + 3 * c.psi s)) = -(4/15) := by
  have h : (∫ s in (0:ℝ)..c.T, (16/135) * (2 * s * c.psi2 s + 5 * c.psi1 s)
        * (2 * s * c.psi1 s + 3 * c.psi s))
      = ∫ s in (0:ℝ)..c.T, ((64/135) * (s ^ 2 * c.psi1 s * c.psi2 s)
        + (96/135) * (s * c.psi s * c.psi2 s) + (160/135) * (s * c.psi1 s ^ 2)
        + (240/135) * (c.psi s * c.psi1 s)) := by
    apply integral_congr; intro s _; ring
  rw [h, c.lincomb]; ring

/-! ## 4. The central pressure Hessian of the compact extension -/

/-- **Equation (10) of `affine-core-pressure.md`.**  Adding the contact term `-tr(L²)/3` of
(4) to the two radial cancellations gives

    p_zz[u_L](0) = -(18/7) b² + (2/5) Ω²

for every radial profile in the class.

The shape of the sum — that `p_zz(0)` *is* `-(1/3) tr(L²)` plus `b²` times the strain radial
integral plus `Ω²` times the swirl radial integral — is formula (4) together with the angular
average (5)-(6) of the note, and is **assumed**, not derived; see the file header.  What is
proved is that this combination evaluates to the stated rational coefficients, for every
`Cutoff`. -/
theorem core_pressure_zz (b Om : ℝ) :
    -(1/3) * (6 * b ^ 2 - 2 * Om ^ 2)
      + b ^ 2 * ((1/2) * (∫ s in (0:ℝ)..c.T, (16/105) * (-4 * (s ^ 2 * c.psi1 s * c.psi2 s)
          + 6 * (s * c.psi s * c.psi2 s) + 2 * (s * c.psi1 s ^ 2)
          + 21 * (c.psi s * c.psi1 s))))
      + Om ^ 2 * ((1/2) * (∫ s in (0:ℝ)..c.T, (16/135) * (2 * s * c.psi2 s + 5 * c.psi1 s)
          * (2 * s * c.psi1 s + 3 * c.psi s)))
      = -(18/7) * b ^ 2 + (2/5) * Om ^ 2 := by
  rw [c.strain_radial, c.swirl_radial]; ring

/-- **Equation (11).**  Axisymmetry makes the two transverse entries equal, and the trace is
`Δp(0) = -tr(L²) = -6b² + 2Ω²`; hence `p_xx = p_yy = (-6b² + 2Ω² - p_zz)/2`. -/
theorem core_pressure_transverse (b Om : ℝ) :
    (-6 * b ^ 2 + 2 * Om ^ 2 - (-(18/7) * b ^ 2 + (2/5) * Om ^ 2)) / 2
      = -(12/7) * b ^ 2 + (4/5) * Om ^ 2 := by ring

/-- The `b = 0` case: the compact strained core reduces to the rotating core of
`active-core-pressure.md`, with the same coefficient `2Ω²/5`. -/
theorem core_pressure_zz_of_pure_rotation (Om : ℝ) :
    -(1/3) * (6 * (0:ℝ) ^ 2 - 2 * Om ^ 2)
      + (0:ℝ) ^ 2 * ((1/2) * (∫ s in (0:ℝ)..c.T, (16/105) * (-4 * (s ^ 2 * c.psi1 s * c.psi2 s)
          + 6 * (s * c.psi s * c.psi2 s) + 2 * (s * c.psi1 s ^ 2)
          + 21 * (c.psi s * c.psi1 s))))
      + Om ^ 2 * ((1/2) * (∫ s in (0:ℝ)..c.T, (16/135) * (2 * s * c.psi2 s + 5 * c.psi1 s)
          * (2 * s * c.psi1 s + 3 * c.psi s)))
      = (2/5) * Om ^ 2 := by
  have := c.core_pressure_zz 0 Om
  linarith [this]

/-! ## 5. The rotating core of `active-core-pressure.md` §1 -/

/-- The note's radial factor `rad = -1/4`: `∫₀^R r ψ(r²) ψ'(r²) dr = -1/4`, with `R = √T`.
This is `moment_one` after the substitution `s = r²`. -/
theorem radial_moment :
    ∫ r in (0:ℝ)..Real.sqrt c.T, r * (c.psi (r ^ 2) * c.psi1 (r ^ 2)) = -(1/4) := by
  have hsub : (∫ r in (0:ℝ)..Real.sqrt c.T, (2 * r) • ((fun s => c.psi s * c.psi1 s) ∘
      (fun r : ℝ => r ^ 2)) r) = ∫ s in ((0:ℝ) ^ 2)..(Real.sqrt c.T ^ 2), c.psi s * c.psi1 s :=
    intervalIntegral.integral_deriv_smul_comp
      (fun x _ => (hasDerivAt_pow 2 x).congr_deriv (by norm_num))
      (by fun_prop) (c.cont0.mul c.cont1)
  rw [Real.sq_sqrt c.hT.le] at hsub
  simp only [smul_eq_mul, Function.comp_apply] at hsub
  have h2 : (∫ r in (0:ℝ)..Real.sqrt c.T, 2 * (r * (c.psi (r ^ 2) * c.psi1 (r ^ 2))))
      = ∫ s in ((0:ℝ) ^ 2)..c.T, c.psi s * c.psi1 s := by
    rw [← hsub]; apply integral_congr; intro s _; ring
  rw [intervalIntegral.integral_const_mul] at h2
  norm_num at h2
  rw [c.moment_one] at h2
  linarith

end Cutoff

/-- The note's angular factor `ang = -4/15`: the uniform spherical average of
`(3μ² - 1)(1 - μ²)` is `½∫_{-1}^{1} (3μ²-1)(1-μ²) dμ = -4/15`.  This is
`verify_active_core_pressure.py`'s `ang`. -/
theorem angular_average :
    (1/2) * (∫ mu in (-1:ℝ)..1, (3 * mu ^ 2 - 1) * (1 - mu ^ 2)) = -(4/15) := by
  have key : ∀ x : ℝ, HasDerivAt (fun t : ℝ => -3 * t ^ 5 / 5 + 4 * t ^ 3 / 3 - t)
      ((3 * x ^ 2 - 1) * (1 - x ^ 2)) x := fun x =>
    (((((hasDerivAt_pow 5 x).const_mul (-3)).div_const 5).add
      (((hasDerivAt_pow 3 x).const_mul 4).div_const 3)).sub (hasDerivAt_id x)).congr_deriv
      (by norm_num; ring)
  rw [integral_eq_sub_of_hasDerivAt (fun x _ => key x)
    ((by fun_prop : Continuous fun x : ℝ => (3 * x ^ 2 - 1) * (1 - x ^ 2)).intervalIntegrable
      _ _)]
  norm_num

/-- **The rotating-core coefficient of `active-core-pressure.md` §1**:
`∂_zz p[c](0) = (2/3)Ω² - 4·ang·rad·Ω² = (2/5)Ω²`, with `ang = -4/15` (`angular_average`)
and `rad = -1/4` (`Cutoff.radial_moment`).  The line is `core = 2/3 - 4*ang*rad` of
`verify_active_core_pressure.py`. -/
theorem rotating_core_coefficient :
    (2:ℝ)/3 - 4 * (-(4/15)) * (-(1/4)) = 2/5 := by norm_num

/-! ## 6. The ratio test at the neutral insertion -/

/-- **The ratio test, as arithmetic.**  This is the `deriv.subs(tuning)` line of
`verify_affine_pressure_independent.py`: the five chain-rule terms of

    β̇ = (-4b² - p_zz/2 + η_b - (b/Ω) η_w)/Ω

at `b' = 2b²`, `Ω' = 2bΩ`, `p_zz = -8b²`, `η_b = η_w = η_b' = η_w' = 0` sum to
`-(p_zz' + 32b³)/(2Ω)`.  The middle term vanishes: `4b²/Ω² + p_zz/(2Ω²)` is `0` exactly
because `p_zz = -8b²` at the neutral point. -/
theorem ratio_second_derivative_algebraic {b Om pzz pzz' : ℝ} (hOm : Om ≠ 0)
    (hpzz : pzz = -8 * b ^ 2) :
    (-8 * b / Om) * (2 * b ^ 2) + (4 * b ^ 2 / Om ^ 2 + pzz / (2 * Om ^ 2)) * (2 * b * Om)
      + (-1 / (2 * Om)) * pzz' = -(pzz' + 32 * b ^ 3) / (2 * Om) := by
  subst hpzz; field_simp; ring

/-- **The ratio test.**  Let `b, Ω, P (= p_zz), η_b, η_w` be differentiable at `0`, with
`Ω(0) ≠ 0`, and suppose the note's first-derivative relation

    β'(t) = (-4 b(t)² - P(t)/2 + η_b(t) - (b(t)/Ω(t)) η_w(t)) / Ω(t)

holds for all `t` (this is the hypothesis `hbeta`; it is transcribed from
`full-feedback-gate-audit.md` and is **not** derived here).  At the neutral insertion

    b'(0) = 2 b(0)²,  Ω'(0) = 2 b(0) Ω(0),  P(0) = -8 b(0)²,  η_b = η_w = η_b' = η_w' = 0,

the second derivative of `β = b/Ω` at `0` is

    β''(0) = -(P'(0) + 32 b(0)³)/(2 Ω(0)),

which is the boxed formula at the end of §5 of `affine-core-pressure.md`. -/
theorem ratio_second_derivative {beta b Om P Eb Ew : ℝ → ℝ} {b1 Om1 P1 Eb1 Ew1 : ℝ}
    (hb : HasDerivAt b b1 0) (hOm : HasDerivAt Om Om1 0) (hP : HasDerivAt P P1 0)
    (hEb : HasDerivAt Eb Eb1 0) (hEw : HasDerivAt Ew Ew1 0) (hOm0 : Om 0 ≠ 0)
    (hbeta : ∀ t, HasDerivAt beta
      ((-4 * b t ^ 2 - P t / 2 + Eb t - (b t / Om t) * Ew t) / Om t) t)
    (hb1 : b1 = 2 * b 0 ^ 2) (hOm1 : Om1 = 2 * b 0 * Om 0) (hP0 : P 0 = -8 * b 0 ^ 2)
    (hEb0 : Eb 0 = 0) (hEw0 : Ew 0 = 0) (hEb1 : Eb1 = 0) (hEw1 : Ew1 = 0) :
    HasDerivAt (deriv beta) (-(P1 + 32 * b 0 ^ 3) / (2 * Om 0)) 0 := by
  have hderiv : deriv beta = fun t =>
      (-4 * b t ^ 2 - P t / 2 + Eb t - (b t / Om t) * Ew t) / Om t :=
    funext fun t => (hbeta t).deriv
  rw [hderiv]
  have hnum := ((((hb.pow 2).const_mul (-4)).sub (hP.div_const 2)).add hEb).sub
    ((hb.div hOm hOm0).mul hEw)
  refine (hnum.div hOm hOm0).congr_deriv ?_
  subst hb1; subst hOm1; subst hEb1; subst hEw1
  simp only [Pi.sub_apply, Pi.add_apply, Pi.mul_apply, Pi.div_apply, Pi.pow_apply,
    hP0, hEb0, hEw0]
  norm_num
  field_simp
  ring

end StrainedCore
end ForcedRoute

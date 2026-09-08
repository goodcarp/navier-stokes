# prove-lagrangian — the Lagrangian route to the upper half of the doubling clock

Sub-seat `lower/prove-lagrangian`, DTC-2026-09-06. Every number below came out of a script in this
folder that I wrote and ran (`s1`…`s7`, re-asserted by `check_constants.py`); `SHA256SUMS` is
computed, never typed. Nothing outside this folder was written. Numerics falsify, never prove;
where a statement is an identity it is marked EXACT and derived.

---

## 0. Headline

The obstacle the brief asks me to clear — *does the strain on the innermost material shell survive
the deformation of the datum over `[0,τ]`?* — turns out not to need a perturbative
`(1 + O(c))` estimate at all. For this datum the deformation is **exactly the axisymmetric strain
`T_λ`, and its effect on the strain is exactly a factor λ**:

> **LEMMA 1 (EXACT).** For the plateau `η₀ = −M sgn(z)/r` on `ρ₀<|x|<R` and the volume-preserving
> map `T_λ : (r,θ,z) ↦ (λr, θ, λ^{-2}z)` with `η` transported (`ω^θ ↦ λω^θ`),
>
> ```  a_λ(0) = (M/2) λ log(R/ρ₀)        exactly, for every λ > 0.  ```

The λ is the vortex-stretching factor and the geometric redistribution is **exactly neutral**. So
the strain never decays under its own leading-order motion; it grows. Closing the resulting
integro-ODE gives a closed form, a constant strictly better than the campaign's first-order one,
and a model blow-up time.

| quantity | value | status |
|---|---|---|
| `Φ(λ) = a_λ(0)/a_0(0)` for the plateau | `= λ` | **EXACT** (three independent routes) |
| conservative constant (uses only `Φ ≥ 1`) | `c₂ = 4 log(3/2) = 1.6218604` | PROVED-MODULO (T,V,R) |
| accelerated constant (uses `Φ = λ`) | **`c₂ = 8(1−√(2/3)) = 1.4680274`** | PROVED-MODULO (T,V,R) |
| admissible taper δ = 7.5°, model | `c₂ = 1.4709` | same |
| model blow-up of the closed ODE | `t = 4/(M log(R/ρ₀))` | model only — **not** an NS claim |

**The theorem is not proved.** Two steps (T and V of §4) are genuine gaps and I name the exact
estimates they need. What *is* proved is the part the brief singled out as the risk.

---

## 1. Set-up and notation

Axisymmetric, no swirl. `η = ω^θ/r`, `Ψ = r²ψ`, and the estate's 5D lift

```
Δ₅ψ = −η ,  Δ₅ = ∂_rr + (3/r)∂_r + ∂_zz   (Laplacian on ℝ⁵ = ℝ⁴_r × ℝ_z)
a := u^r/r = −∂_zψ ,  u^z = 2ψ + r∂_rψ ,  D_tη = νΔ₅η ,  D_t ω^θ = a ω^θ + ν r Δ₅η .
```
`G₅ = 1/(8π²|x|³)`, so with `y ∈ ℝ⁵`, `ρ = |y|`, `t = cos φ = z/ρ`,

```
a(0) = (3/8π²) ∫ (−z_y)/|y|⁵ η(y) dy₅ = (3/4) ∬ (−cos φ sin³φ) η dρ dφ .          (1.1)
```

**Datum.** `ω^θ₀ = −M sgn(z) h_δ(φ)` on `ρ₀<ρ<R`, `h_δ(φ) = min(1, φ_ax/δ)`, `φ_ax = min(φ, π−φ)`;
then `‖ω₀‖_∞ = M` and `‖η₀‖_∞ ≤ M/(ρ₀ sin δ)` (admissible — the bare bang-bang is not). Mollify at
scale `ερ` (scale-adapted). Put `L = log(R/ρ₀)`, and fix the viscous floor by
`ρ₀ δ = √(ν/M)`, so the *finest feature of the datum* — the taper width at the inner edge — is
exactly the dissipation length.

Tracked material point `x(t)`: starts at `ρ = ρ₀`, `φ₀ = 30°` (so `h_δ(φ₀) = 1`, `|ω^θ₀| = M`).
Inviscidly `ω^θ(x(t),t) = M exp(∫₀^t a)` exactly, since `η` is materially conserved and
`ṙ = a r`.

---

## 2. LEMMA 1 — the strain-transport identity (EXACT, PROVED)

`T_λ` commutes with isotropic dilation and is volume preserving in ℝ³; in the 5D lift its Jacobian
is `J₅ = λ²` and `η` is transported. Changing variables in (1.1) to the **material** label, the
image of `−z_y/|y|⁵` times `J₅` is `−z/(λ²r²+λ^{-4}z²)^{5/2}` in the *original* coordinates, so with
`Q(φ,λ) = λ²sin²φ + λ^{-4}cos²φ` the ρ- and φ-integrals **factorise exactly** (no endpoint error,
because the strained plateau still has the scale-invariant radial profile `1/ρ`):

```
a_λ(0) = (M/2) · P_h(λ) · L ,   P_h(λ) = 3 ∫₀¹ h(v) v² (A v² + B)^{-5/2} dv ,
v = sin φ ,  A = λ² − λ^{-4} ,  B = λ^{-4} ,  A + B = λ² .                        (2.1)
```

**For `h ≡ 1` the integral is elementary and gives `P₁(λ) = λ`.** The antiderivative is
`v³ / (3B (Av²+B)^{3/2})`; at `v = 1` it is `1/(3B(A+B)^{3/2}) = λ⁴/(3λ³) = λ/3`.
`s1_exact_lemma.py` verifies `d/dv[antiderivative] − integrand ≡ 0` **symbolically** (sympy residual
`0`), `P₁(λ) − λ ≡ 0` symbolically, and `|3∫ − λ|/λ = 0.0` in 30-digit mpmath at seven λ.

Three independent confirmations that this is not a bookkeeping artefact:

| route | what it computes | result |
|---|---|---|
| (i) closed form / mpmath | `P₁(λ)` | `= λ`, residual `0` (exact) |
| (ii) Lagrangian quadrature | `P₁(λ)` | rel. dev. from λ ≤ **1.4e-12** at λ = 1…3 |
| (iii) **Eulerian**, never uses `J₅` — builds `η_λ = η₀∘T_λ^{-1}` on a grid and integrates (1.1) | `a_λ(0)/(ML)` | `0.50005, 0.62500, 0.74999, 0.99999` vs `λ/2 = 0.5, 0.625, 0.75, 1.0`; worst rel **1.1e-4** (grid) |
| (iv) **Gegenbauer projection** of the strained profile (`s6`) — `κ(λ) = −3H₁(λ)/5` | `κ(λ)` | `0.500000, 0.550000, 0.625000, 0.750000` at λ = 1, 1.1, 1.25, 1.5; rel err **1.9e-8** |

**First-order coefficient.** `P'_h(1) = 3∫h v²(−5)(3v²−2)dv / (…)`. For the plateau this is exactly
`+1` (the weight `v²` and the sign change of `1−3cos²φ` at `cos²φ = 1/3` conspire); computed:
`1.0000000000`. For the campaign's simulated profile `ω^θ = −M sin2φ`, `Φ'(1) = 10/7 = 1.42857`.
For the 7.5° taper, `0.99504`. **`Φ' > 0` in every case: the strain accelerates.**

**Taper (EXACT deficit bound).** `P_{h_δ}(λ) = λ − D(λ,δ)` with
`0 ≤ D ≤ w³ / (B (A w² + B)^{3/2})`, `w = sin δ` — the same antiderivative. At δ = 7.5°, λ = 3/2 the
true deficit is `0.02727` and the bound `0.10042`, so `P_{h_δ}(3/2) ≥ 1.39958`. Verified
(`s2_taper.py`): `Φ_{h_δ}(λ) := P_{h_δ}(λ)/P_{h_δ}(1) ≥ 1` on `λ ∈ [1, 3/2]` for **every** taper
tested (δ = 3°, 5°, 7.5°, 10°, 15°, 20°, 30°), with `Φ_{h_δ}(3/2) = 1.4736 (7.5°)`, `1.0674 (30°)`.

`κ_δ := ½P_{h_δ}(1)` = **0.4999917 (5°), 0.4997212 (7.5°), 0.4978077 (15°), 0.4836252 (30°)** —
independent reproduction of `refuter-correctness` §6 (`0.49972, 0.49781, 0.48363`) to 5 digits,
and `κ = 2/5` exactly for `−M sin 2φ`, reproducing `refuter-scope-priorart` §2.

---

## 3. LEMMA 2 — sign preservation (PROVED, non-perturbative)

The datum is odd in `z` (`η₀` odd) and the axisymmetric no-swirl NS system preserves that symmetry
exactly, so `η(·,t)` is odd and `η = 0` on `{z = 0}`, which is a material plane. `η` solves the
drift-diffusion equation `∂_tη + u·∇₅η = νΔ₅η` on `{z>0}` with zero boundary data, so by the
parabolic maximum principle `η ≤ 0` on `{z>0}` for all `t`, for every `ν ≥ 0`.

Consequently the integrand of (1.1) is **pointwise non-negative for all time**:
`(−z_y)η(y,t) ≥ 0`. Hence `a(0,t) ≥ 0` unconditionally, and `a(0,t)` is bounded below by its
restriction to *any* sub-region of sources. This is what makes the argument robust: no cancellation
can occur, only omission.

---

## 4. The Lagrangian argument, step by step, with status

### (1) `t = 0`. **PROVED.**
`a(0,0) = (M/2)·2κ_δ·L` by (2.1). On the innermost *material* shell, the full four-region C¹-matched
Gegenbauer solution (`s5_velocity_structure.py`, an instrument independent of both the campaign's
elliptic-integral Biot–Savart and its 5D quadrature) gives, at `ρ = ρ₀`, `φ = 10°`:

| R/ρ₀ | `a(0,0)` | `(M/2)log(R/ρ₀)` | `a(material) − (M/2)log R` |
|---|---|---|---|
| 64 | 2.07944154 | 2.07944154 | **0.2168076** |
| 256 | 2.77258872 | 2.77258872 | **0.2167801** |
| 1024 | 3.46573590 | 3.46573590 | **0.2167784** |
| 4096 | 4.15888308 | 4.15888308 | **0.2167783** |

`a(0,0) − (M/2)log(R/ρ₀) = −1.24e-9` at every R (machine/truncation), and the material offset (R = 4096)
`0.2167783` reproduces `refuter-correctness` §5's `0.216773 ± 2e-6` to **5e-6** by a completely
different method. Coefficients: `H₁ = −5/6`, `H₃ = 3/40`, `H₅ = −247/1680` reproduced to 1e-11
(estate `c6`).

### (2) Persistence of the strain on `[0,τ]`. **This is where the brief expected the work.**

**L3 (the velocity is a per-shell pure strain). PROVED for the initial datum, modulo an absolute
mode-sum bound that is verified numerically.**
`η` homogeneous of degree −1 ⟹ the source `ρ^{-1}C_l^{3/2}` resonates with the radial operator
`Δ₅(ρ^γ C_l) = [γ(γ+3) − l(l+3)]ρ^{γ-2}C_l` **only at `l = 1`** (`γ = 1`: `4 − 4 = 0`). Since
`Δ₅(ρ log ρ · C₁) = 5ρ^{-1}C₁` and `H₁ = −5M/6`, `C₁ = 3t`,

```
ψ_{l=1} = −A(ρ) z ,   A(ρ) = (M/2) log(R/ρ)   ⟹   u = A(ρ)·(r, −2z) + O(Mρ) ,
```
i.e. **the entire `O(M log)` part of the velocity is the uniform axisymmetric strain**, and every
`l ≥ 3` mode is non-resonant, contributing `f_l = H_lρ/(l(l+3)−4)`, i.e. `ψ_dev = Mρ S(t)` — a
velocity of size `O(Mρ)`, a factor `1/L` below. Explicit constants, computed (`s5`, `s4`):

* `sup_t |a − A(ρ)| / M = 0.4197` from the `l ≥ 3` bulk modes + the `l = 1` remainder `−(M/2)cos²φ`;
* `sup_t |u^z − (−2A z)| / (Mρ) = 0.19143`;
* over the **whole** shell including both boundary layers, `sup |a − A(ρ)| = 0.7047 M` (attained at
  the outer edge), and `0.3187 M` one octave in from each end;
* `|H_l| ~ l^{-3/2}` (`H_801 = −7.33e-5`), so `Σ_{l≥3}|H_l|‖C_l‖_∞/(l(l+3)−4) = 0.3958` converges
  absolutely; partial sums of `S(t)` are stable to 1e-6 by `l = 201`.

**L3v (the same on the DEFORMED field). SKETCH; numerically verified, not proved.**
`T_λ` commutes with dilation, so the strained plateau is again homogeneous of degree −1 and the
whole L3 argument applies verbatim with the strained angular profile. `s6_deformed_field.py`
computes its Gegenbauer coefficients: `κ(λ) = −3H₁(λ)/5 = λ/2` to 1.9e-8 (route (iv) above), and the
L3 remainder bound `Σ_{l≥3}|H_l(λ)|‖C_l‖_∞/(l(l+3)−4)` = **0.3958, 0.4354, 0.4948, 0.5938** at
λ = 1, 1.1, 1.25, 1.5 (0.6069 for the 7.5° taper at λ = 3/2) — **bounded uniformly on the range the
theorem uses.** What is missing is a proof of that uniform bound rather than a computation of it.

**T (trajectory closure). GENUINE GAP — SKETCH.**
Given L3v, along a trajectory the velocity is `A(ρ)(r,−2z) + θ(y)` with `|θ| ≤ C₃Mρ`, `C₃ ≤ 0.71`.
Over `τ = c/(ML)` this gives, for each material shell:
* the pure-strain part integrates to exactly `T_{λ(s)}`, `λ(s) = exp(∫₀^τ A)` — coaxial strains
  commute, so the composite map *is* a single `T_λ`;
* `log ρ` moves by at most `2∫a dt = 2 log(3/2) = 0.811` e-folds over the whole run, so
  `A(ρ)` changes along the trajectory by `≤ (M/2)(0.811)`, i.e. **relative `O(1/L)`**;
* the remainder `θ` displaces a point by `≤ C₃Mρτ = C₃ρ c/L`, i.e. **relative `O(c/L)`**;
* the shear between shells: `ρ ∂_ρλ/λ = −τM/2 = −c/(2L)`, so the 5D Jacobian and the log-radial
  measure are preserved to `1 + O(c/L)`.

*The estimate that is missing* is the stability of the functional (1.1) under an `O(c/L)`
relative perturbation of the flow map in the scale-invariant metric `|δX|/|X|` — i.e.
`|a[η∘Φ⁻¹] − a[η∘T_λ⁻¹]| ≤ C ‖ |δX|/|X| ‖_∞ · (M L)`. It is a Lipschitz estimate for a
Calderón–Zygmund-type kernel in log coordinates; I did not prove it. Everything downstream is
conditional on it.

**P (the integro-ODE). PROVED given T and Lemma 1.**
Write `σ = s/L ∈ [0,1]` for the material shell, `θ = MLt`, `F(σ,θ) = ∫₀^t a` along its trajectory,
`λ = e^F`. The strain at shell `σ` is generated by the shells outside it (the inside contributes
`O(M)`: the 5D kernel falls like `(ρ'/ρ)⁴`, so the inner sum converges geometrically), and by
Lemma 1 each outer shell's contribution is multiplied by exactly `Φ(λ) = λ`. Hence

```
∂_θF(σ,θ) = κ H(σ,θ) + O(1/L),   H(σ,θ) = ∫_σ¹ e^{F(σ',θ)}dσ' ,  F(·,0) = 0 .
```
Differentiating, `∂_θ∂_σH = −e^F ∂_θF = κ H ∂_σH = ∂_σ(κH²/2)`, and `H(1,θ) ≡ 0` kills the
integration constant, so the system **closes**:

```
∂_θH = (κ/2)H² ,  H(σ,0) = 1−σ   ⟹   H = (1−σ)/(1 − (1−σ)κθ/2) ,
F(σ,θ) = −2 log(1 − (1−σ)κθ/2) ,   λ(σ,θ) = (1 − (1−σ)κθ/2)^{-2} .              (4.1)
```
`s3_ode.py` verifies (4.1) against the defining integral identity (max |H − ∫_σ¹e^F| ≤ 4e-14 at
θ = 0…2) and against a 4001-point method-of-lines solve of the integro-ODE (max |F_num − F_closed|
= 2.4e-10 at θ = 0.734, 2.5e-9 at θ = 1.5). `s7_tapered_ode.py` re-derives it by direct integration
with the tabulated `P_{h_δ}` and reproduces the bang-bang value to 7e-10.

### (3) Viscosity. **SKETCH — GENUINE GAP.**
Exactly, on a plateau `Δ₅(1/r) = 2/r³ − 3/r³ = −1/r³` (sympy: `r³Δ₅(1/r) = −1`), so
`D_t log|η| = −ν/r²` with **no** `1/σ²` core penalty. At the tracked point, with
`ρ₀δ = √(ν/M)` and `r = ρ sin φ`, the penalty in units of `M` is
`δ²/sin²φ = 0.0685 (φ=30°)`, `0.0216 (φ=62.9°)` at `ρ = ρ₀` and falls like `ρ^{-2}` outward;
integrated over `τ = c/(ML)` it costs `O(1/L)`. The diffusion length is
`√(ντ) = ρ₀δ√(c/L) ≪ ρ₀δ = ` the finest feature, for `L ≫ c/δ²`.

*The estimate that is missing:* a Gaussian (Aronson-type) upper bound for
`∂_tη + u·∇₅η = νΔ₅η` whose constants depend on the drift only through `‖∇u‖_∞ τ = O(c)`,
giving `|η(x(t),t) − η_transported| ≤ ‖η₀‖_∞ exp(−d²/(Cντ))` with `d` the material-frame distance
from the tracked point to the nearest datum discontinuity (`d ≳ ρ₀ δ`). I did not write it.

### (4) Conclusion.
`ω^θ(x(t),t) = M h_δ(φ₀) exp(∫₀^t a) · (1 − viscous loss)`, and `‖ω(t)‖_∞ ≥ |ω^θ(x(t),t)|`.

**Conservative form** (uses only `Φ_h(λ) ≥ 1`, i.e. Lemma 1 in the weak form the taper certainly
satisfies): `∂_θF ≥ κ_δ(1−σ) − C_⋆/L`, so `F(0,θ) ≥ (κ_δ − C_⋆/L)θ` and

```
M₀T_d ≤ log(3/2)/(κ_δ − C_⋆/L)/L        ⟹    T(Re) ≤ (4 log(3/2) + o(1))/log Re
                                                     = (1.6218604 + o(1))/log Re .
```

**Accelerated form** (uses `Φ = λ`, (4.1) with `κ = κ_δ`): `λ(0,θ) = 3/2` at `κ_δθ = 2(1−√(2/3))`,

```
M₀T_d ≤ 2(1−√(2/3))/(κ_δ L) · (1+O(1/L))  ⟹  T(Re) ≤ (8(1−√(2/3)) + o(1))/log Re
                                                     = (1.4680274 + o(1))/log Re ,
```
using `log Re_E = 2L + 2log(1/δ) − 0.7031659 + O(ε,δ)` (`E = 0.172403978 M²R⁵`, estate `c6`,
confirmed by `refuter-scope-priorart` to 2.5e-7), so `L = ½log Re_E + O(1)`.
Direct integration with the tapered profile (`s7`) gives `c₂ = 1.46889 (5°), 1.47086 (7.5°),
1.47451 (10°), 1.48818 (15°), 1.59019 (30°)` — all strictly below `4log(3/2)`.

Trajectory bookkeeping at `λ = 3/2`: the inner shell moves at most `0.811` e-folds in `log ρ`, and
its polar angle obeys `tan φ = tan φ₀ · λ³`, so `φ₀ = 30° ↦ 62.83°` — still `27°` from the equatorial
mollification layer and `55°` from the taper.

---

## 5. Status per step

| step | statement | status |
|---|---|---|
| **L1** | `a_λ(0) = (M/2)λ log(R/ρ₀)` exactly; taper deficit bound; `Φ_{h_δ} ≥ 1` on `[1,3/2]` | **PROVED** (symbolic + 3 independent numerical routes) |
| **L2** | sign preservation ⟹ `a(0,t) ≥ 0` for all `t`, all `ν ≥ 0` | **PROVED** (symmetry + parabolic max principle) |
| **L3** | `u = A(ρ)(r,−2z) + O(Mρ)` for the initial datum, `C₃ ≤ 0.705` | **PROVED-MODULO** absolute convergence of `Σ|H_l|‖C_l‖_∞/l²` (computed = 0.3958, not proved) |
| **L3v** | same, uniformly for the deformed field, `λ ∈ [1,3/2]` | **SKETCH** (computed ≤ 0.607; needs a uniform-in-λ bound) |
| **T** | flow map on `[0,τ]` = `T_{λ(s)}` up to relative `O(c/L)`; and (1.1) is stable under that | **GAP** — needs the scale-invariant Lipschitz stability of (1.1) stated in §4(2) |
| **P** | the integro-ODE closes; closed form (4.1) | **PROVED** given T, L1 |
| **V** | viscosity costs `O(1/L)` | **SKETCH** — needs the Aronson bound stated in §4(3) |
| **R** | mollification: `‖ω₀‖_∞ = M(1+O(ε))`, `E`, `Re_E` bookkeeping | **SKETCH** (routine, not written) |
| **THEOREM** | `T(Re) ≤ c₂/log Re`, `c₂ = 1.6219` (conservative) / `1.4680` (accelerated) | **PROVED-MODULO T, V, R** |

**Therefore: (ii) is not proved, and Bradshaw–Farhat–Grujić ARMA 2019 Thm 10 is not refuted here.**
What this seat closes is the specific worry the brief raised — that the strain degrades as the
configuration deforms. It does not: it grows, exactly.

---

## 6. Self-attack (the hour), and what survived

1. **"The `Φ = λ` identity is an artefact of the Lagrangian Jacobian `J₅ = λ²`."** — Refuted.
   Route (iii) builds `η_λ` on an Eulerian grid and never touches `J₅`; agrees to 1.1e-4 (grid
   resolution). Route (iv) projects onto Gegenbauer modes; agrees to 1.9e-8.
2. **"The identity needs the bang-bang profile, and the admissible taper breaks it."** — Survives.
   The exact deficit bound is `w³/(B(Aw²+B)^{3/2})`; the monotonicity `Φ_{h_δ} ≥ 1` on `[1,3/2]`
   holds for every δ up to 30°, and the tapered ODE moves `c₂` by 0.2% at 7.5°.
3. **"Material rotates to the equator, where the weight `|cos φ|sin²φ` vanishes — that must kill
   the strain."** — It does not: the exact φ-integral (2.1) *includes* that rotation, and the
   compensating inward motion near the axis (`ρ ↦ λ^{-2}ρ`) exactly balances it. This is the
   content of `∫₀¹v²(Av²+B)^{-5/2}dv = λ/3`.
4. **"Near-axis material collapses onto the origin and the `λ^{10}` kernel amplification blows the
   taper deficit up."** — Real effect, bounded: at λ = 3/2, δ = 7.5° the deficit is 0.0273 against
   λ = 1.5 (bound 0.1004). It is the reason δ must be a small constant (or `δ ≫ √(c/L)`), not why
   the argument fails.
5. **"The campaign measured `D_t a < 0`."** — That was the *ring stack*
   (`exact-first-order` §4), where rings self-propel and disperse. For the plateau the model gives
   `D_t a|₀ = M²L²/8 > 0` and `β = D_ta + a² = (3/8)M²L²`. The campaign's own **plateau** runs
   corroborate the sign: their measured per-e-fold strain rises `0.411 → 0.473` during the run
   (`viscous-numerics` §4), and `T₃₂·M·a* = 0.479 > log(3/2) = 0.4055`.
6. **"Then the model should match the six viscous runs."** — It over-predicts the speed.
   Model `T·a*/∫a = 1.1086`; measured `0.479/0.4055 = 1.1814` (model conservative there).
   Model `Q' = M T₃₂ log(R/s*) = 2(1−√(2/3))/κ = 0.8930` at their `κ = 0.411`; measured
   `0.976 ± 0.018` — the model is **9% optimistic** at their `L = log(R/s*) ≈ 2.9`. That is the
   size of the `O(1/L)` terms the model drops (and their `sup|η|` falls 14% to viscosity). So the
   numerics corroborate the *sign* of the acceleration and bound its *size* from above; they do not
   support `c₂ = 1.468` at finite `L`, and the conservative `4 log(3/2)` is the value I would
   defend in print.
7. **"The closed ODE blows up at `θ = 2/κ` — that is a blow-up claim."** — It is not, and must not
   be quoted as one. `λ → ∞` means the shells collapse onto `z = 0`; the `O(1/L)` error terms and
   L3v's uniformity are controlled only for `λ ≤ O(1)`, and `A(ρ)` stops being slowly varying along
   trajectories once `∫a dt ~ L`. The model is used only up to `λ = 3/2`, i.e. `κθ = 0.367`, less
   than a fifth of the way to its own singularity.
8. **"`a` at the origin is not `a` where the vorticity is."** — Handled: the C¹-matched field gives
   `a(material, ρ₀, 10°) − (M/2)log(R/ρ₀) = +0.21678 M`, a *favourable* bounded offset, confirming
   `refuter-correctness` §5 independently.

---

## 7. What is new here (beyond the two construct seats and the two refuters)

1. **Lemma 1**, `a_λ(0) = (M/2)λ log(R/ρ₀)` exactly — the geometric redistribution under the
   datum's own strain is exactly neutral, so the strain scales precisely like the stretched
   vorticity. This replaces the perturbative step (2) the brief asked for with an identity.
2. **The closed integro-ODE (4.1)** and `λ(σ,θ) = (1 − (1−σ)κθ/2)^{-2}`.
3. **`c₂ = 8(1−√(2/3)) = 1.4680274`**, strictly below the campaign's first-order `4log(3/2)`,
   with the acceleration ratio `θ_accel/θ_frozen = 2(1−√(2/3))/log(3/2) = 0.9051502` in closed form.
4. **A third instrument** (Gegenbauer + four-region C¹ matching) reproducing `a(0) = (M/2)log(R/ρ₀)`
   to 1.2e-9 and the material offset `0.2167783` against the refuter's `0.216773 ± 2e-6`.
5. **Lemma 2**: `a(0,t) ≥ 0` for all `t` and all `ν ≥ 0`, unconditionally — the first
   non-perturbative statement in this line about `t > 0`.
6. The exact identification of the two remaining estimates (§4 T and V) as the whole gap.

## 8. Files

| file | what it establishes |
|---|---|
| `s1_exact_lemma.py` / `s1_results.json` | Lemma 1 three ways; `P₁(λ) = λ`; `κ` for all profiles; `Φ'(1)` |
| `s2_taper.py` / `s2_results.json` | exact taper-deficit bound; `Φ_{h_δ} ≥ 1` on `[1,3/2]` for δ ≤ 30° |
| `s3_ode.py` / `s3_results.json` | the integro-ODE closes; closed form vs method-of-lines; `θ₃₂`, `c₂` |
| `s4_constants.py` / `s4_results.json` | every displayed constant; `Re_E` bookkeeping; `sup|a−A|`; viscous-run cross-check |
| `s5_velocity_structure.py` / `s5_results.json` | Lemma 3: `H_l`, the `l=1` resonance, `O(Mρ)` constants, full C¹ matching, `0.2167783` |
| `s6_deformed_field.py` / `s6_results.json` | self-attack: `κ(λ) = λ/2` from the deformed field's modes; L3 remainder uniform in λ |
| `s7_tapered_ode.py` / `s7_results.json` | the model run with the admissible tapered datum |
| `check_constants.py` | asserts every number displayed in this NOTE |

# BLOCK 3, the hypothesis (Γ-off): `‖∇u‖_∞` for the TRUE field on the window

Seat `s3close/round2/u2`, DTC-2026-09-06, sitting of 2026-09-08.
Laws: `TORMENT NEXUS/LAWS.md`. Every number below came out of a script in **this**
folder (`u1`…`u5`), written and run here; `SHA256SUMS` is computed with `shasum`, never
typed. Nothing outside this folder was written or read into a script. `round2/u1` was not
read (leak guard). Numerics falsify, they never prove; the algebra marked EXACT is
sympy-verified with residual `0`.

---

## 0. Headline

**BLOCK 3 closes, on the whole of `ℝ⁵`, without homogeneity, without the exact plateau,
without Theorem A, and without the flow map.**

> **THEOREM (Γ-off).**  Let `η` be a bounded classical solution of `D_tη = νΔ₅η` on
> `ℝ⁵×[0,τ]` for the 5-D lift of axisymmetric no-swirl NS, with `z`-odd datum `η₀`
> satisfying `E₀ := sup|x||η₀|/M < ∞`, `𝔊₀ := sup|x|²|∇η₀|/M < ∞`. Suppose on `[0,s]`
> `‖ω^θ(·,σ)‖_∞ ≤ λM` (`λ ≤ 3/2`) and `Γ(σ) := ‖∇u(·,σ)‖_{L^∞(ℝ⁵)} ≤ Γ̄`. Put
> `c_G := ∫₀^sΓ`, `c_R := ∫₀^sΓ_rad` (`Γ_rad` = §4), `p := 1 + 2c_R/c_G ≤ 3`. Then for
> **every** `x ∈ ℝ⁵`
> ```
>   ‖∇u(x,s)‖_op  ≤  2 a(0,s) · (1 + C₁/L)  +  C″ M ,        C₁ = 0 ,
>   C″ = 2 Ĉ_a  +  (3π²/16) e^{p c_G} 𝔊₀  +  λ ,
>   Ĉ_a = (0.29199853 + 0.01475427)λ
>         + min( 3.9992184 λ/σ_* + πλ/8 ,  3.9992184 e^{c_R} E₀ ) ,
> ```
> `σ_*` any angular floor with `sin φ(x) ≥ σ_*` (the second collar branch needs none, so
> `σ_* = 0` is admissible and the bound is global). `3π²/16 = 1.8505508252042546` is
> `C_K · π⁴/2` — the exact Riesz composition constant in `ℝ⁵`, computed in `u1`.

For the campaign datum (`δ = 7.5°`, equatorial mollifier `w = 0.20`, radial ramp
`ε_r = 0.25`, i.e. `w₀ = 0.25ρ₀` at the inner edge): `E₀ = 7.66060196`, `𝔊₀ = 20.9197070`
(`u2`). At `λ = 3/2`, `μ ≤ 0.05`, `σ_* = 1/2` (which contains the whole tube `N_τ`), the
window fixed point `c_G = (λ + C″/L)c`, `c = 2log(3/2)`, gives

| column | `L_*` | `C″` at `L = 10⁴` | `c_G` | `p` |
|---|---|---|---|---|
| proved, map-free, `σ_* = 1/2` | **`4997.6`** | **`582.24`** | `1.2636` | `2.1068` |
| proved, map-free, `σ_* = 0` (whole `ℝ⁵`) | **`5313.8`** | **`801.05`** | `1.2814` | `2.1097` |
| crude (`p ≡ 3`), `σ_* = 1/2` | `9908.9` | `3604.7` | `1.5087` | — |
| map-sharpened (`ν = 0`), `μ = 0.05`, `σ_* = 1/2` | `736.5` | `361.34` | `1.2457` | `2.0652` |

**`μ` is not priced because `μ` does not enter.**  The proved column never uses `Φ_s`; the
map is needed only to sharpen the radial factor, where `μ = 0 → 0.02 → 0.05` moves `L_*`
from `670.0` to `696.2` to `736.5` and `C″(10⁴)` from `329.53` to `342.05` to `361.34`.

**Five things this note establishes beyond the assembly's BLOCK 3 statement.**

1. The hypothesis is provable **for `‖∇u‖` globally on `ℝ⁵`**, not on a slab. That is what
   Theorem V.1 and the map Grönwall actually need (both take a global `sup‖∇b‖`); the
   assembly's `𝒮_s` was a hole, and the axis wedge — where the far/near collar bound
   `4/sin φ` diverges — is closed by a *radius-localised maximum principle* for `ρ|η|`.
2. **`C₁ = 0`.** No multiplicative `1 + C₁/L` correction is needed: `a_far(0,s) ≤ a(0,s)`
   holds exactly, because `κ(ρ',s) ≥ 0` at every radius and every time (Lemma 2).
3. **Theorem A is not used, anywhere.** Nor is `(B.2)`, `(B.3)`, `(B.4)`, the L3v sum, the
   Gegenbauer expansion, homogeneity, or the plateau's profile. §9.
4. The route the brief asks for — BV transport with a pushforward inequality — is **false
   for `ν > 0`** (§5.4), and the correct replacement is cheaper: two weighted maximum
   principles that are **exactly sharp on the reference strain** (`u5`: relative error
   `0.0`, both).
5. `Ĉ_a(λ=1, σ_*=1/2) = 8.697888775` reproduces `ASSEMBLY §2.4`'s `C_a = 8.697888` to
   seven digits from the far/near constants re-derived here, and the shell TV density
   `𝒥_ac = 65.62590996` reproduces `hk2 §8(a)`'s `65.625910` — two decorrelated
   confirmations of neighbouring seats.

**Grade.** The chain is PROVED modulo the two items named in §11: the contradiction
hypothesis `‖ω‖ ≤ λM` (which is (S3)'s own, not an extra assumption) and the standard
approximate-maximum argument for the two weighted principles. The cost is a large `C″`:
`L_*` moves from the refuted-Γ route's `679 … 1435` to `≈ 5.0·10³`, a factor `3.5 … 7.4`.
That is the price of dropping homogeneity, and it is not the assembly's binding constraint
(`L_* = 4.53·10¹⁴` from BLOCK 2's feedback).

---

## 1. Setting, datum, hypotheses, and the set

`ℝ⁵ = ℝ⁴_y × ℝ_z`, `r = |y|`, `ρ = |x|`, `t = cos φ = z/ρ`, `s_φ = sin φ = r/ρ`.
`η = ω^θ/r`, `−Δ₅ψ₁ = η`, `a = u^r/r = −∂_zψ₁`, `u^z = 2ψ₁ + r∂_rψ₁`,
`b(y,z) = (a y, u^z)`, `D_t = ∂_s + b·∇₅`, `D_tη = νΔ₅η`.

**Datum (D-C).**  The campaign datum of `hk2` (D-B) with the radial ramp written in the
log variable, so both radial edges carry the same *relative* width and the inner edge
agrees with `hk2`'s `tanh` ramp exactly (`ρΘ'(ρ₀) = 1/(2ε_r) = 2 = ρ₀/(2w₀)` at
`w₀ = 0.25ρ₀`):

```
  ω^θ₀ = −M · tanh(sin φ/sin δ) · tanh(cos φ/w) · Θ(ρ) ,      δ = 7.5°, w = 0.20 ,
  Θ(ρ) = ½[ tanh((log(ρ/ρ₀) − ε_r)/ε_r) − tanh((log(ρ/R) + ε_r)/ε_r) ] ,  ε_r = 0.25 ,
  η₀ = −(M/ρ) 𝔥(φ) Θ(ρ) ,   𝔥 = tanh(s_φ/sin δ) tanh(t/w)/s_φ ,   L = log(R/ρ₀) .
```
This is `BLOCK 6`'s required family (a sharp radial edge makes `K₂` infinite; the ramp is
mandatory). It is `C^∞`, `z`-odd, `‖ω^θ₀‖_∞ ≤ M`.

**Hypotheses.**

* **(H1)** `η` is a bounded classical solution on `ℝ⁵×[0,τ]`, `η(·,0) = η₀` as above.
* **(H2)** `η₀` is odd in `z`. (Then `η(·,s)` is odd for all `s`, and by **Lemma 2** of
  `lower/prove-lagrangian` §3, `η ≤ 0` on `{z>0}` for all `s` and all `ν ≥ 0`.)
* **(H3)** the contradiction hypothesis of (S3): `‖ω^θ(·,σ)‖_∞ ≤ λ(σ)M ≤ (3/2)M` for
  `σ ∈ [0,τ]`. This is not an extra assumption — §3.
* **(H4)** the bootstrap posture of `B(t)`: `Γ(σ) := ‖∇u(·,σ)‖_{L^∞(ℝ⁵)} ≤ Γ̄` for
  `σ ≤ s`; `c_G := ∫₀^sΓ ≤ Γ̄τ`, `τ = c/(ML)`, `c = 2log(3/2)(1+ε)`.
* **(H5)** `E₀ := sup_x|x||η₀(x)|/M` and `𝔊₀ := sup_x|x|²|∇η₀(x)|/M` are finite.

**The set the theorem holds on.**  All of `ℝ⁵`, at every `s ∈ [0,τ]`. The constant `C″`
depends on the evaluation angle only through `min(3.9992184λ/σ_* + πλ/8, 3.9992184e^{c_R}E₀)`
in `Ĉ_a`; the second branch has no angular dependence at all, so the theorem needs no
angular exclusion. The two branches cross at `sin φ = 0.05834932` (`u3`), i.e. `φ = 3.345°`
from the axis, so on `{sin φ ≥ 1/2}` — which contains the tube `N_τ` at every `λ ∈ [1,3/2]`
(`hk2 §8(b)`: the tracked point runs `φ = 30° → 62.833°`) and the great majority of the
shell — the smaller branch is in force.

This matters. **Theorem V.1 takes `Γ` as a global `sup_{ℝⁿ}‖∇b‖`**, and the map Grönwall
(2.1)–(2.2) of `ASSEMBLY §2.3` takes `sup_x` over the shell including the axis. A slab
statement would not have discharged either.

---

## 2. STEP 1 — the exact reduction (`u1`, EXACT)

For any axisymmetric no-swirl field, in the frame `y = (r,0,0,0)`:

```
  div₅ b = 2a                                                              (sympy residual 0)
  ∇₅b = a·diag(1,1,1,1,−2) + E ,  E = [[ r∂_r a , r∂_z a ],
                                       [ r∂_z a − ω^θ , −r∂_r a ]]        (sympy residual 0)
  ‖∇u‖_op = ‖∇₅b‖_op ≤ 2|a| + r|∇a| + |ω^θ|                                     (2.1)
```
`u1` re-derives all three for a generic non-polynomial `ψ₁ = F(r²,z)` (spanning family,
residual the zero matrix), and checks (2.1) over `4000` random `(a,Q,P,ω)`: worst ratio
`0.99893185`, no violation. The bound is essentially sharp, not slack — the same finding
the L3v refuter reports (`0.9985` over 2000 samples), reproduced here independently.

(2.1) holds for **any** axisymmetric no-swirl field. It is the whole of Part B of the L3v
seat that does not survive off the plateau; (2.1) itself does.

So `Γ(s) ≤ sup_x [ 2|a(x,s)| + r|∇a(x,s)| + |ω^θ(x,s)| ]`, and the three terms are §5–§7.

---

## 3. STEP 2 — `|ω^θ| ≤ λM` on the window

`T_d(u₀)` is defined as the first time `‖ω(t)‖_∞ ≥ (3/2)M`, and (S3) runs by contradiction
from `T_d > τ`. Hence on `[0,τ)`, `‖ω^θ(·,s)‖_∞ < (3/2)M` **by the definition of the
window**; write `M_s ≤ λ(s)M` with `λ(s) ≤ 3/2`. `ASSEMBLY §1.3` states this explicitly and
calls it the one piece of free bookkeeping in the assembly. Status: **PROVED given the
contradiction hypothesis**, which is (S3)'s own.

*The route the brief suggests does not work and is not needed.* "max principle for `η`
times `r ≤ λρ_*(1+…)`" gives `|ω^θ| = r|η| ≤ r‖η₀‖_∞`, and `‖η₀‖_∞ = M/(ρ₀ sin δ)`, not
`M/ρ₀` — the taper makes `η₀` largest **on the axis**, where `h_δ(φ)/sin φ → 1/sin δ`
(`u2`: `sup ρ|η₀|/M = tanh(1/w)/sin δ = 7.66060196`, against `1/sin δ = 7.66129758`).
That route would cost a factor `7.66` and is strictly worse. **Brief-layer correction
(L-17): `|η| ≤ M/ρ₀` is false for the tapered datum; the true constant is `M/(ρ₀ sin δ)`.**

---

## 4. STEP 3a — the radial rate `Γ_rad` (`u1`, `u5`)

**Lemma 4.1 (linear growth of the lift). PROVED.**  `b(0,s) = 0`, hence
`|b(x,s)| ≤ Γ(s)|x|`.
*Proof.* `η` odd in `z` ⟹ `ψ₁` odd in `z` ⟹ `ψ₁(0) = 0` ⟹ `u^z(0) = 2ψ₁(0) = 0`; and
`b_i = a y_i → 0`. `u1` item B returns `b(0) = (0,0)` symbolically for a generic `z`-odd
`ψ₁`. Then `b(x) = ∫₀¹(∇₅b)(tx)·x\,dt`. ∎

**Lemma 4.2 (the radial rate). PROVED.**  With
`Γ_rad(s) := sup_x Λ_max(sym ∇₅b(x,s))`,
```
  x·b(x,s) = |x|² ∫₀¹ x̂·sym(∇₅b)(tx,s) x̂ dt  ≤  Γ_rad(s)|x|² ,
  Λ_max(sym ∇₅b) = max( a , ½[ −a + √( (3a+2Q)² + (2P−ω)² ) ] ) ,   Q = r∂_ra, P = r∂_za, ω = ω^θ .
```
*Proof.* The transverse `ℝ³ ⊂ ℝ⁴` carries the eigenvalue `a`; the `(ê_r,ê_z)` block of the
symmetric part is `[[a+Q, P−ω/2],[P−ω/2, −2a−Q]]`, trace `−a`, so its top eigenvalue is
`½(−a + √((3a+2Q)²+(2P−ω)²))`. `u5` (a): `2·10⁵` random samples, worst relative error
against `eigvalsh` **`1.02·10^{-15}`**. ∎

**Corollary 4.3 (the envelope). PROVED.**  If `|a| ≤ A`, `a ≥ −Ĉ_aM` (Lemma 2 + §7),
`|Q|,|P| ≤ Ĝ`, `|ω| ≤ λM`, then
```
  Γ_rad ≤ max( A , A + 2Ĝ + λM/2 , 2Ĉ_aM + 2Ĝ + λM/2 ) .                       (4.1)
```
`u5` (b): `2·10⁵` random points in the admissible box, **`0` violations**, worst ratio
`0.99842` (so (4.1) is attained).

**Why this is the whole game.** `Γ̄ = 2A + Ĝ + λM` and `A = a(0,s) + Ĉ_aM ≍ (M/2)λL`, so
`Γ_rad/Γ̄ → 1/2` as `L → ∞` and the exponent `p = 1 + 2Γ_rad/Γ̄ → 2`. For the pure
axisymmetric strain `b = (ar,−2az)` the ratio is **exactly** `1/2` (`u5`:
`Γ_rad/Γ = 0.5`), which is why §5's two principles are exactly sharp there.

---

## 5. STEP 3b — the two weighted maximum principles

These replace the brief's BV-transport step. They are pointwise, need no flow map, and are
sharp on the reference strain.

**Lemma 5.0 (no zeroth-order term). EXACT.**  `η` is a genuine scalar on `ℝ⁵` obeying the
**advective** equation with the full Cartesian `Δ₅`, so `∂_j` commutes with `Δ₅` and
```
   D_t(∂_jη) = −(∂_jb)·∇₅η + νΔ₅(∂_jη) .                                       (5.1)
```
`div₅b = 2a ≠ 0` does **not** appear: it would appear only in the conservation form
`∂_sη + div₅(bη) = νΔ₅η + 2aη`, which is a different and here unnecessary bookkeeping —
the same two answers `gap-V-aronson §2.3` gives, confirmed here symbolically (`u1` items C:
both residuals `0`). This settles the brief's question: **no zero-order term appears in the
lift.**

**Lemma 5.1 (the weight identities). EXACT.**  In `ℝ⁵`,
`Δ₅(ρ^k g) = ρ^kΔ₅g + 2kρ^{k−2}x·∇g + k(k+3)ρ^{k−2}g`, hence with `V = ρg`, `W = ρ²g`,
```
   ρ  Δ₅g = Δ₅V − (2/ρ²) x·∇V − 2V/ρ² ,      ρ² Δ₅g = Δ₅W − (4/ρ²) x·∇W − 2W/ρ² .   (5.2)
```
(`u1` item D: four residuals `0`.) **Both rearrangements carry the same good sign
`−2/ρ²`**, and the first-order terms vanish at an interior extremum. That is the whole
mechanism.

> **THEOREM P1 (radius-localised sup of `η`).  PROVED.**
> ```
>    sup_x |x| |η(x,s)|  ≤  e^{c_R} sup_x |x| |η₀(x)|  =  e^{c_R} E₀ M ,   c_R := ∫₀^s Γ_rad .
> ```
*Proof.* `D_t|η| ≤ νΔ₅|η|` (Kato; no `Γ` term, because `D_tη = νΔ₅η` exactly).
`D_t|x| = (x·b)/|x| ≤ Γ_rad|x|` by Lemma 4.2. So with `V := |x||η|`,
`D_tV ≤ Γ_rad V + ν|x|Δ₅|η| = Γ_radV + ν[Δ₅V − (2/ρ²)x·∇V − 2V/ρ²]` by (5.2). `V` is
continuous on `ℝ⁵`, `V(0) = 0` (since `η` is bounded), and `V → 0` at infinity — datum
(D-C) has `Θ(ρ) = O(ρ^{−2/ε_r}) = O(ρ^{−8})` at large `ρ` and `O(ρ^{8})` at small `ρ`, and
the heat–transport semigroup preserves that decay up to Gaussian tails. At a positive
interior maximum
`∇V = 0`, `Δ₅V ≤ 0`, so `d/ds (sup V) ≤ Γ_rad sup V` in the sense of the upper Dini
derivative (Hamilton's trick). Grönwall. ∎

> **THEOREM P2 (radius-localised sup of `∇η`).  PROVED.**
> ```
>    sup_x |x|² |∇η(x,s)|  ≤  e^{c_G + 2c_R} sup_x |x|²|∇η₀(x)|  =  e^{p c_G} 𝔊₀ M ,
>    p := 1 + 2 c_R/c_G  ≤  3 .
> ```
*Proof.* `g := ∇₅η` obeys (5.1), so `D_t|g| ≤ ‖∇₅b‖_op|g| + νΔ₅|g| ≤ Γ|g| + νΔ₅|g|` (Kato).
`D_t(|x|²) = 2x·b ≤ 2Γ_rad|x|²`. With `W := |x|²|g|`,
`D_tW ≤ (Γ + 2Γ_rad)W + ν[Δ₅W − (4/ρ²)x·∇W − 2W/ρ²]` by (5.2). Same maximum-principle
step. ∎

**Sharpness (`u5` (c), pre-declared control).** For the pure strain, `η_λ = η₀∘T_λ^{-1}`
has exactly `sup ρ|η_λ| = λ·sup|α||η₀|` and `sup ρ²|∇η_λ| = λ⁴·sup|α|²|∇η₀|`
(`‖T_λ‖ = λ`, `‖T_λ^{-1}‖ = λ²`, `u5` (d)). The principles give `e^{c_R} = λ` and
`e^{c_G+2c_R} = λ⁴`: **relative error `0.0` for both, at `λ = 1.25` and `λ = 3/2`.**
The crude versions (`e^{c_G} = λ²`, `e^{3c_G} = λ⁶`) are lossy by `λ` and `λ²`. So all of
the loss in the theorem sits in `Γ_rad ≤ Γ`, none in the weights.

### 5.4 What the brief's BV route gives, and where it fails

The brief asks for `|Dη(s)| ≤ e^{c_G+Cc}·Φ_{s#}|Dη₀|` **as measures**. Three findings.

* **(a) For `ν = 0` the pushforward identity is exact and the factor is
  `‖(DΦ)^{−T}‖·J_Φ ≤ e^{c_G}·e^{2c_G}`, not `e^{c_G}.**  With `f = η₀∘Φ^{-1}`,
  `∫ϕ d|Df| = ∫ϕ(Φ(α))|(DΦ(α))^{−T}∇η₀(α)| J_Φ(α)dα`, so the pushforward carries the
  Jacobian as well as the cotangent factor, and `∂_s log J_Φ = 2a(Φ) ∈ [−2Γ,2Γ]`. The
  brief's `‖(DΦ)^{−T}‖` alone is not the density.
* **(b) For `ν > 0` the stated inequality is FALSE.**  `Φ_{s#}|Dη₀|` is supported on
  `Φ_s(supp|Dη₀|)`; the equatorial jump of `η₀` lives on the `4`-plane `{z=0}`, a
  Lebesgue-null set which the flow preserves. For `ν > 0` the solution is smooth and
  `|Dη(s)|` is absolutely continuous with mass spread over a layer of thickness `≍√(νs)`
  **off** that null set. No constant multiple of a measure supported on `{z=0}` dominates
  it. The corrected statement must be a mollified one, and the brief's clause cannot be
  used as written.
* **(c) The correct global statement is available and costs exactly `e^{3c_G}`, with the
  diffusion contributing nothing** (so the brief's extra `e^{Cc}` has `C = 0`). From (5.1)
  and Kato, `D_t|g| ≤ Γ|g| + νΔ₅|g|`; integrating over `ℝ⁵` and using
  `∫b·∇|g| = −∫(div₅b)|g| = −∫2a|g|` with `|a| ≤ Γ`,
  ```
      d/ds ∫|∇η(·,s)| dx  ≤  3Γ ∫|∇η|dx        ⟹    TV(η(s)) ≤ e^{3c_G} TV(η₀) .
  ```
  This matches the `ν = 0` transport factor `e^{c_G}·e^{2c_G}` of (a) exactly, so diffusion
  is free. What it does **not** give is the *shell* density `𝒥_s` the far-field integral
  needs; localising it in radius costs a weight `χ` with `χ(w)e^{3w}` integrable, i.e.
  `|χ'|/χ ≥ 4`, hence a factor `e^{7c_G}` — worse than P2's `e^{pc_G}` with `p ≤ 3`. **That
  is why this note uses the pointwise route instead of the BV route.**

The price of the pointwise route against the TV route is `|S⁴|𝔊₀/𝒥 = 26.319·20.920/65.626
= 8.39` at `c_G = 0` — the concentration of `|∇η₀|` in the taper transition is not
exploited. The price of the localised TV route is the extra exponent `e^{(7−p)c_G}`. At
`c_G = λ·2log(3/2) = 1.2164` (`λ=3/2`) that is `e^{4c_G} = 129.8` for `p = 3` and
`e^{5c_G} = 437.89` for `p = 2`, so the pointwise route wins by **`15.47×`** and
**`52.19×`** respectively.

---

## 6. STEP 4 — the kernel route for `r|∇a|` (`u1`, `u3`)

`a = K*η` with `K(w) = −∂_zG₅(w) = 3w_z/(8π²|w|⁵)`, `G₅ = 1/(8π²|x|³)`,
`|K(w)| ≤ C_K|w|^{-4}`, `C_K = 3/(8π²) = 0.0379954439`, `C_K|S⁴| = 1` exactly (`u1`).
Since `∇η ∈ L¹∩L^∞`, `∇a = K*∇₅η` (differentiation under the convolution), hence

```
   |∇a(x,s)|  ≤  C_K ∫_{ℝ⁵} |x−x'|^{-4} |∇η(x',s)| dx'
              ≤  C_K e^{p c_G} 𝔊₀ M ∫_{ℝ⁵} |x−x'|^{-4}|x'|^{-2} dx'      (P2)
              =  C_K e^{p c_G} 𝔊₀ M · (π⁴/2)/|x|
```
by the **Riesz composition** in `ℝ⁵`: `∫|x−x'|^{-α}|x'|^{-β}dx' = C|x|^{n−α−β}` with
`n=5, α=4, β=2` gives `C = π^{5/2}Γ(½)Γ(3/2)Γ(½)/(Γ(2)Γ(1)Γ(2)) = π⁴/2 = 48.704545517`.
`u1` verifies it two ways: the classical Gamma-function formula, and a direct
spherical-average quadrature (`48.70454551715526` vs `48.70454551700121`, relative
`3.16·10^{-12}`). Therefore, with `r = ρ sin φ`,

> **(6.1)   `r|∇a(x,s)| ≤ (3π²/16) e^{p c_G} 𝔊₀ M · sin φ ≤ 1.8505508 e^{p c_G} 𝔊₀ M`,
> `3π²/16 = C_K·π⁴/2`.**

**What this route does not need.**

* **No split at `d`.** `hk2`'s Cor. 4.3 splits `∇a` at a distance `d` because it works with
  the *measure* `|Dη|`, which can be singular. Under P2 the density is a bounded function
  with `|x'|^{-2}` decay and `|x−x'|^{-4}` is locally integrable in `ℝ⁵`
  (`∫_{|w|<d}|w|^{-4}dw = |S⁴|d`), so the whole integral is absolutely convergent. **The
  transported layer sets, their positions to relative `μ`, the collar quadrature and the
  shell TV density `𝒥_s` all drop out.**
* **No log.** The bound is `O(1/ρ)` with an absolute constant. Structurally this is
  `hk2`'s Remark 5.3 — the `ML` lives in `a`, never in `∇a` — obtained here without the
  modal decomposition and without Lemma 5.2's far/near split: the Riesz constant is finite
  because `4 + 2 > 5`, full stop.
* **No principal value, no `C^α` modulus, no mollification scale, no `1/δ_m²`.**

**Cross-check of magnitude.** At `c_G = 0`, (6.1) gives `r|∇a| ≤ 1.8506·20.920 = 38.72 M`
against `hk2 §8(c)`'s `|∇a| ≤ 4.6795 M/ρ₀` at `λ = 1`, i.e. `r|∇a| ≤ 2.340M` at
`sin φ = ½`, for the same datum by the TV route — a ratio `16.55`, against the predicted
`(|S⁴|𝔊₀/𝒥)·(1/½) = 16.78`; the residual `1.4 %` is the difference between the two radial
ramps. The measured value is `r|∇a| = 0.258 M` (`hk2 §9`), reproduced independently in §9
below.

---

## 7. STEP 5 — `2|a|` against `a(0,s)` (`u3`)

Write the far/near split of `rebuild/far-near-kernel-lemma` at `2ρ` and `ρ/2`:
`a = a_far + a_inner + a_collar`. Its hypotheses are **only** `|ω^θ| ≤ M_s` and `z`-oddness
— both true for the true field at every time, by (H3) and (H2)+Lemma 2. All three constants
are **re-derived here** from their own definitions (Cauchy–Schwarz against Parseval on the
interior/exterior multipole series, and the bathtub constant), not imported:

| constant | this seat (`u3`) | `far-near-kernel-lemma` |
|---|---|---|
| far Taylor remainder, `z`-odd | **`0.29199853254882`** | `0.291999` |
| inner multipole, `z`-odd | **`0.014754271582666`** | `0.014754` |
| collar coefficient `2R_A`, `R_A = (2⁵−2^{-5})^{1/5}` | **`3.9992184446453`** | `3.999218` |
| collar axis piece | `π/8 = 0.3926991` | `π/8` |

**(a) The far part and `C₁ = 0`.**  `a_far(0,s) = ∫_{2ρ}^{∞}κ(ρ',s)dρ'/ρ'` and
`a(0,s) = ∫_0^{∞}κ(ρ',s)dρ'/ρ'`, where `κ(ρ',s) = −(3/5)g₁(ρ',s) =
−(3/4)∫_0^π ω^θ(ρ',φ,s)cos φ sin²φ dφ` is the `l = 1` interior mode, which by
*Consequence A* of that lemma is **exactly constant inside its own shell**. By **Lemma 2**,
`ω^θ ≤ 0` on `{z>0}` and (by oddness) `≥ 0` on `{z<0}`, so `ω^θ cos φ ≤ 0` everywhere and
`κ(ρ',s) ≥ 0` at every radius and every time. Hence
```
     a_far(0,s) ≤ a(0,s)          for every x, every s.                       (7.1)
```
That is where `C₁ = 0` comes from: no multiplicative slack is needed because the omitted
shells contribute non-negatively. This is the **second** load-bearing use of Lemma 2 in the
campaign (the assembly's §2.4 is the first).

**(b) The collar, made axis-safe by P1.**  `hk2`/`far-near`'s collar bound
`|a_collar| ≤ (3.9992184/sin φ + π/8)λM` diverges on the axis, and that divergence is
**real** for a datum with no taper (`far-near §6`: `c_edge ~ ¼log(1/φ)`). P1 removes it:
on `A = {ρ/2 < |x'| < 2ρ}`, `|η(x',s)| ≤ e^{c_R}E₀M/|x'| ≤ 2e^{c_R}E₀M/ρ`, and by
rearrangement `∫_A|x−x'|^{-4}dx' ≤ |S⁴|R_A` with `R_A = (2⁵−2^{-5})^{1/5}ρ = 1.99960922ρ`
(the annulus has the volume of a ball of that radius), so with `C_K|S⁴| = 1`
```
    |a_collar(x,s)| ≤ 2 R_A · 2 e^{c_R} E₀ M /2 = 3.9992184 e^{c_R} E₀ M ,      (7.2)
```
**with no angular dependence at all.**  The two branches cross at
`sin φ = 0.05834932` (`u3`), `φ = 3.345°`.

**(c) Assembled.**
```
   |a(x,s)| ≤ a(0,s) + Ĉ_a M ,
   Ĉ_a = (0.29199853+0.01475427)λ + min( 3.9992184λ/σ_* + πλ/8 , 3.9992184 e^{c_R}E₀ ) .  (7.3)
```
Numbers (`u3`, at `c_R = c_G = λ·2log(3/2)`):

| `λ` | `Ĉ_a`, `σ_* = 1/2` | `Ĉ_a`, `σ_* = sin 7.5°` | `Ĉ_a`, `σ_* = 0` (whole `ℝ⁵`) |
|---|---|---|---|
| 1.00 | **`8.697888775`** | `31.338654` | `69.238699` |
| 1.25 | `10.872399` | `39.173318` | `84.807500` |
| 1.50 | `13.046833` | `47.007982` | `103.858160` |

The `λ = 1`, `σ_* = 1/2` entry is `ASSEMBLY §2.4`'s `C_a = 8.697888` — **reproduced to
seven digits from constants re-derived in this folder**, a decorrelated confirmation of
both the far/near seat and the assembly.

---

## 8. STEP 6 — assembly, and the window fixed point (`u3`)

From (2.1), (6.1), (7.3) and (H3):

> **`Γ(s) ≤ 2a(0,s) + C″M`,  `C″ = 2Ĉ_a + (3π²/16)e^{pc_G}𝔊₀ + λ`,  `C₁ = 0`.**

`C″` is not a free constant: hypothesis (i) of `B(t)` feeds back through `c_G`.
`a(0,s) ≤ (λM/2)L` (since `κ ≤ M_s/2`), so `Γ̄ = λML + C″M` and
`c_G = Γ̄τ = (λ + C″/L)c`, with `p = 1 + 2Γ_rad/Γ̄` from (4.1). `u3` solves the fixed point
in `(Ĝ, c_G, p)` and bisects for the smallest `L` at which it exists.

### 8.1 `C″` at `L → ∞` (feedback frozen at `c_G = λc`, `c = 2log(3/2)`)

| `λ` | `c_G` | `Ĉ_a` (`σ_*=1/2`) | `Ĝ` (`p=3`) | `Ĝ` (`p=2`) | `Ĝ` (map, `μ=0.05`) | `C″` (`p=2`) |
|---|---|---|---|---|---|---|
| 1.00 | `0.8109` | `8.6979` | `440.97` | `195.98` | `96.03` | **`214.38`** |
| 1.25 | `1.0137` | `10.8724` | `810.10` | `293.98` | `183.77` | **`316.97`** |
| 1.50 | `1.2164` | `13.0468` | `1488.26` | `440.97` | `324.11` | **`468.56`** |

Axis-inclusive (`σ_* = 0`), `p = 2`: `C″ = 335.46 / 464.84 / 650.18` at `λ = 1/1.25/1.5`.

### 8.2 The fixed point: `L_*` and `C″(L)`

`σ_* = 1/2` (contains `N_τ`):

| column | `λ` | `L_*` | `C″(10³)` | `C″(10⁴)` | `c_G(10⁴)` | `p(10⁴)` |
|---|---|---|---|---|---|---|
| proved (map-free) | 1.00 | `2189.7` | — | `232.82` | `0.8298` | `2.0629` |
| proved | 1.25 | `3277.2` | — | `361.48` | `1.0430` | `2.0790` |
| proved | 1.50 | **`4997.6`** | — | **`582.24`** | `1.2636` | `2.1068` |
| crude `p≡3` | 1.50 | `9908.9` | — | `3604.69` | `1.5087` | — |
| map, `μ=0` | 1.50 | `670.0` | — | `329.53` | `1.2431` | `2.0591` |
| map, `μ=0.02` | 1.50 | `696.2` | — | `342.05` | `1.2441` | `2.0615` |
| map, `μ=0.05` | 1.50 | `736.5` | — | `361.34` | `1.2457` | `2.0652` |

`σ_* = 0` (the whole of `ℝ⁵`):

| column | `λ` | `L_*` | `C″(10⁴)` | `c_G(10⁴)` | `p(10⁴)` |
|---|---|---|---|---|---|
| proved | 1.00 | `2435.2` | `362.85` | `0.8404` | `2.0635` |
| proved | 1.25 | `3577.2` | `526.76` | `1.0564` | `2.0803` |
| proved | 1.50 | **`5313.8`** | **`801.05`** | `1.2814` | `2.1097` |
| map, `μ=0.05` | 1.50 | `1172.3` | `557.90` | `1.2616` | `2.0654` |

`σ_* = sin 7.5°`: `L_* = 5088.2` (proved, `λ=3/2`), `C″(10⁴) = 658.12`; map `μ=0.05`,
`L_* = 788.2`, `C″(10⁴) = 431.16`.

### 8.3 Pricing `μ`

`μ` enters **nowhere in the proved column**. It enters only if one replaces P2's radial
factor `e^{2c_R}` by the map bound `((1+μ)λ)²` (valid at `ν = 0`: `ρ = |Φ(α)| ≤ (1+μ)λ|α|`
and `‖(DΦ)^{-1}‖ ≤ e^{c_G}`). Then:

| `μ` | 0 | 0.02 | 0.05 |
|---|---|---|---|
| `L_*` (`λ=3/2`, `σ_*=1/2`) | `670.0` | `696.2` | `736.5` |
| `C″(10⁴)` | `329.53` | `342.05` | `361.34` |
| `C″` at `L→∞` | `321.5704` | `333.4470` | `351.7030` |

So `μ ≤ 0.05` is worth `+9.37 %` on `C″` (`L→∞`) and `+9.92 %` on `L_*`. **The map is a
convenience here, not a hypothesis.**

---

## 9. STEP 7 — the numerical confrontation (`u4`)

**The instrument (built here; shares no code with `hk2`'s `k3/k4/k9` or L3v's `g`-series).**
Polar coordinates centred on the *evaluation* point remove the kernel singularity exactly,
with no regularisation and no principal value: writing `w = x − x'`,
`ŵ = (n sin α, cos α)`, `n = cos β ê₁ + sin β ê_⊥ ∈ S³`,
```
   dw = ξ⁴dξ · sin³α dα · 4π sin²β dβ ,      K(w) = 3cos α/(8π²ξ⁴)
   ⟹  K(w) dw = (3/2π) cos α sin³α sin²β dξ dα dβ        (the ξ⁴ cancels exactly)
   a(x)     = (3/2π) ∫∫∫ cos α sin³α sin²β · η(x') dβ dα dξ ,
   ∂_j a(x) = (3/2π) ∫∫∫ cos α sin³α sin²β · (∂_jη)(x') dβ dα dξ ,
```
`r' = √(r² − 2rξ sin α cos β + ξ²sin²α)`, `z' = z − ξ cos α`.

**Its controls (pre-declared before any number was read).**

| control | result |
|---|---|
| **R1** `a(0)` against the exact `l = 1` closed form `κ_prof·∫Θ dlog ρ` (a route that shares nothing with the quadrature) | relative **`1.09·10^{-8}`** (`L=10`), **`2.63·10^{-9}`** (`L=40`) |
| **R2** grid refinement (`α: 200→400`, `β: 32→64`, `ξ`-panels `0.75→0.5` decades) | `a` moves `1.9·10^{-6}` … `2.9·10^{-8}`; `∂_ra` `≤4.0·10^{-4}`; `∂_za` `≤9.3·10^{-5}` |
| **R3** `κ_prof` for the sharp-taper profile against `ASSEMBLY §1.2` | `0.49972123052108863` vs `0.4997212305210886` — **16 digits** |
| **R4** `L`-independence of `r|∇a|` (a `log(R/ρ₀)` would move it by `4×` from `L=10` to `L=40`) | `0.56703` → `0.56705` (`λ=1`), `0.96387` → `0.96387` (`λ=3/2`) — relative **`3.53·10^{-5}`** and **`1.25·10^{-8}`** |

**The perturbed map.**  `Φ = Λ∘(I+ψ)`, `Λ = T_λ`, `ψ(α) = (ψ_r ŷ, ψ_z)` with
`ψ_r = ε ρ (1−t²) sin(2π u)`, `ψ_z = ε ρ t(1−t²) cos(2π u)`, `u = log(ρ/ρ₀)/L` — smooth,
axisymmetric, and `z`-parity preserving (`ψ_r` even, `ψ_z` odd in `z`), so `η` stays `z`-odd
and Lemma 2 still applies. `ε` is normalised so that `max(sup|ψ|/|α|, ‖Dψ‖_op) = ‖ψ‖_{C¹}`:

| `‖ψ‖_{C¹}` | `sup|ψ|/|α|` | `‖Dψ‖_op` | `μ = sup|Φ−Λ|/|Λ| ≤ λ³sup|ψ|/|α|`, `λ=1` | `λ=3/2` |
|---|---|---|---|---|
| 0.02 | `0.0166542` (`L=10`) / `0.0172925` (`L=40`) | `0.0200000` | `0.01665` | `0.05621` |
| 0.05 | `0.0416355` / `0.0432312` | `0.0500000` | `0.04164` | `0.14052` |

So the `μ ≤ 0.05` of the brief corresponds to `‖ψ‖_{C¹} = 0.05` at `λ = 1` and to
`‖ψ‖_{C¹} ≈ 0.018` at `λ = 3/2`; the `‖ψ‖_{C¹} = 0.05`, `λ = 3/2` row is `μ = 0.14`, i.e.
**three times outside** the brief's map hypothesis, and is reported as a stress case.

**The measurement.**  Slab sample: `4` radii log-spaced in `[2ρ₀, R/2]` × `3` angles
(`φ = 30°, 60°, 90°`), all with `sin φ ≥ 1/2`; `Γ` from (2.1) as the exact operator norm
`max(|a|, σ_max([[a+Q,P],[P−ω,−2a−Q]]))`, not the envelope.

| `L` | `λ` | field | `a(0)` | `max Γ` | `max Γ − 2a(0)` | `max r|∇a|` |
|---|---|---|---|---|---|---|
| 10 | 1.00 | exact | `4.50998` | `8.49649` | **`−0.52347`** | `0.56703` |
| 10 | 1.00 | `‖ψ‖=0.02` | `4.50886` | `8.50075` | `−0.51697` | `0.57478` (`+1.37 %`) |
| 10 | 1.00 | `‖ψ‖=0.05` | `4.50634` | `8.50532` | `−0.50736` | `0.58633` (`+3.40 %`) |
| 10 | 1.50 | exact | `6.69585` | `12.57542` | **`−0.81629`** | `0.96387` |
| 10 | 1.50 | `‖ψ‖=0.02` | `6.69520` | `12.57948` | `−0.81092` | `0.97093` (`+0.73 %`) |
| 10 | 1.50 | `‖ψ‖=0.05` | `6.69356` | `12.58392` | `−0.80321` | `0.98142` (`+1.82 %`) |
| 40 | 1.00 | exact | `18.75201` | `36.98055` | **`−0.52347`** | `0.56705` |
| 40 | 1.00 | `‖ψ‖=0.05` | `18.74756` | `36.97815` | `−0.51696` | `0.59040` (`+4.12 %`) |
| 40 | 1.50 | exact | `27.84065` | `54.86501` | **`−0.81629`** | `0.96387` |
| 40 | 1.50 | `‖ψ‖=0.05` | `27.83646` | `54.86332` | `−0.80961` | `0.96977` (`+0.61 %`) |

**Three things the table says.**

* `max Γ − 2a(0)` is **`L`-independent to five digits** (`−0.52347` at `L = 10` and at
  `L = 40`; `−0.81629` likewise): the theorem's structure `Γ ≤ 2a(0,s) + C″M` with `C″`
  free of `L` is confirmed, and the true `C″` on the sampled slab is *negative*
  (`−0.523M` at `λ=1`, `−0.816M` at `λ=3/2`) — the slab points sit at radii where the
  outer log count is short of `a(0)` by more than the `O(M)` remainder.
* **A `5 %` `C¹` perturbation of the flow map moves `Γ` by `0.1 %` and `r|∇a|` by
  `0.6–4.1 %`.** The quantity BLOCK 3 is about is not sensitive to the map at the
  brief's `μ`, which is the numerical form of the analytic finding that `μ` does not enter
  the proof.
* **`max r|∇a| < λM` at every sampled point, for the exact map and for both perturbed
  maps** (`0.56703`, `0.58633` at `λ=1`; `0.96387`, `0.98142` at `λ=3/2`). `L3v §B3`'s
  measured consequence for the *exact strained plateau* is `sup_t(1−t²)|a'_dev| = λM`
  **exactly**; so the mollified datum and its `C¹`-perturbed transports stay inside the
  homogeneous field's own sharp value. This is the same check `hk2 §9` runs on its
  `r|∇a| = 0.258` at one point, extended here to twelve points, two `λ`, two `L` and two
  perturbations.

**The slack.**  At `L = 10` and `L = 40` the window fixed point does not exist
(`L_* ≈ 5·10³`), so the comparison is made against the **feedback-frozen** constants of
§8.1 (`c_G = λ·2log(3/2)`, `p = 2`), which is what the theorem gives in the limit:

| quantity | `λ = 1` | `λ = 3/2` |
|---|---|---|
| bound on `r|∇a|`, `(3π²/16)e^{2c_G}𝔊₀` | `195.98 M` | `440.97 M` |
| measured `max r|∇a|` (exact map) | `0.56703 M` | `0.96387 M` |
| **slack** | **`345.6×`** | **`457.5×`** |
| bound `C″` (`σ_*=1/2`) | `214.38 M` | `468.56 M` |
| measured `max Γ − 2a(0)` | `−0.523 M` | `−0.816 M` |
| bound on `Γ` at `L = 40` | `251.88 M` | `524.24 M` |
| measured `max Γ` at `L = 40` | `36.98 M` | `54.87 M` |
| **slack in `Γ` at `L = 40`** | **`6.81×`** | **`9.55×`** |

The slack in `Γ` falls like `1 + C″/(2a(0))`, i.e. `→ 1` as `L → ∞`; at the existence
boundary `L = L_* = 4997.6` (`λ=3/2`, `σ_*=1/2`, where the fixed point gives
`C″ = 1035.59`, `c_G = 1.3844`, `p = 2.3544`) it is **`1.149`**, and at
`L_* = 5313.8` (`σ_*=0`, `C″ = 1387.74`) it is **`1.188`** — using the measured
`a(0) = 0.696 ML` and `Γ = 2a(0) − 0.816M` of the table above. The `457×` in `r|∇a|` is the honest number for the *sub-estimate*,
and it decomposes as `e^{2c_G} = 11.4` (the window, sharp on the reference strain — `u5`),
`|S⁴|𝔊₀/𝒥 = 8.4` (pointwise versus total-variation density, §5.4), and `≈ 4.8`
(`|·|`-inside-the-kernel — the same step `hk2 §9` prices at `30×` for `K₂`).

### 9.1 A by-product finding: `ASSEMBLY §1.3`'s floor `ε_δ` omits the equatorial mollifier

The control R1 needs `κ_prof = −(3/5)g₁/M`, so it prices the `l = 1` strain rate of the
*actual* datum. `u6` isolates the features (`κ = 1/2` exactly for bang-bang):

| profile | `κ` | `1 − 2κ` | clock floor `ε = (1−2κ)/(2κ)` |
|---|---|---|---|
| bang–bang | `0.5000000000` | `0` | `0` |
| sharp taper `δ = 7.5°` only | `0.4997212305` | `5.5754·10^{-4}` | `5.5785·10^{-4}` |
| `tanh` taper `δ = 7.5°` only | `0.4984964209` | `3.0072·10^{-3}` | `3.0163·10^{-3}` |
| linear equatorial `δ_m = 5°` only | `0.4981012077` | `3.7976·10^{-3}` | — |
| `tanh` equatorial `w = 0.20` only | `0.4762377156` | `4.7525·10^{-2}` | — |
| **`ASSEMBLY §1.1`'s datum** (sharp taper `7.5°` + linear `δ_m = 5°`) | **`0.4978224382`** | **`4.3551·10^{-3}`** | **`4.3742·10^{-3}`** |
| campaign (D-B)/(D-C) (`tanh` taper + `tanh` `w=0.20`) | `0.4747343554` | `5.0531·10^{-2}` | `5.3221·10^{-2}` |

`ASSEMBLY §1.3` writes `ε_δ := 1 − 2κ_δ = 5.5754·10^{-4}` and calls it "the floor:
`ε → (1−2κ_δ)/(2κ_δ) = 5.5785·10^{-4}` as `L → ∞`". That is `κ` for the **taper alone**.
`§1.1`'s own datum also carries `min(1,|φ−π/2|/δ_m)` with `δ_m = 5°`, which costs
`3.7976·10^{-3}` — **6.8× more than the taper** — so the correct floor for the sealed datum
is `ε → 4.3742·10^{-3}`, a factor **`7.84`** above the displayed value. For the datum the
viscous runs and `hk2` actually use (`tanh(cos φ/0.20)`, i.e. `δ_m ≈ 11.5°` in angle) the
floor is `5.3221·10^{-2}`, a factor `95`. Nothing fails: `t_* = 2log(3/2)(1+ε)/(ML)` is
`0.44 %` (not `0.056 %`) above its `ε = 0` value for the sealed datum, and `5.3 %` for the
run datum. But the displayed floor is wrong, and `δ_m` is the term that sets it, not `δ`.
Filed for the ledger; not load-bearing for BLOCK 3.

---

## 10. What is closed

**Theorem A is no longer needed for BLOCK 3.**  `ASSEMBLY §5.2`'s dependency graph routes
`write/L3v-and-gamma-bound`'s Theorem A / Cor. A1 into two places: "inside Theorem Γ" and
"BLOCK 2". Theorem Γ was the only consumer of the `l^{-3/2}` jump estimate for hypothesis
(i), through `‖f‖_∞ = Σ_{l≥3}|c_l|‖C_l‖_∞` in `(B.5)`/`(B.6)`. This note proves hypothesis
(i) without any modal expansion: the chain is (2.1) → P1/P2 → Riesz → far/near. So

* **Theorem A, Corollary A1, the L3v sum, `(B.2)`, `(B.3)`, `(B.4)`, `(B.5)`, `(B.6)`, the
  slab, homogeneity and the Gegenbauer machinery are all off the critical path for
  BLOCK 3.** They remain the *sharp* route for the exact plateau (`C′ ≤ 33.5` claimed,
  `≤ 119.33 … 151.15` after the refuter's correction, against `C″ ≈ 582` here), and the
  L3v refuter's MAJOR-1 therefore no longer propagates into `B(t)`(i).
* **Theorem A is still needed for BLOCK 2** (the restart lemma off the ansatz manifold),
  which is where the assembly's `L_* = 4.53·10¹⁴` actually comes from.
* `rebuild/far-near-kernel-lemma` becomes **the** non-perturbative estimate in the chain,
  used twice (§7 and, through `C_a`, in `ASSEMBLY §2.4`). Its three `z`-odd constants are
  re-derived here.

**What BLOCK 3 items remain.**

1. **The constant.** `C″ ≈ 582` (`σ_*=1/2`) / `801` (global) at `L = 10⁴`, `λ = 3/2`,
   against a measured `Γ − 2a(0)` of order `M` (§9). The slack is in `Γ_rad ≤ Γ̄/2 + O(M)`
   (which is sharp), in `sup|x|²|∇η|` versus the true concentration (a factor `8.4`
   against the TV density, §5.4), and in `|·|`-inside-the-kernel (the same step `hk2 §9`
   prices at `30×`).
2. **`L_*` ≈ `5.0·10³`** against the exact-plateau route's `679 … 1435`. A factor
   `3.5 … 7.4`. Not binding against BLOCK 2.
3. **The map-sharpened column is proved only at `ν = 0`.** Extending
   `sup ρ²|∇η| ≤ ((1+μ)λ)²e^{c_G}𝔊₀M` to `ν > 0` needs a localised (Theorem V.1-style)
   version of P2, not written here.
4. **`Γ_rad` is bounded by (4.1), which is attained** (worst ratio `0.99842`). Sharpening
   it needs `a ≥ 0` pointwise, which is false in general (only `a(0,s) ≥ 0` is proved).
5. Everything downstream of hypothesis (i) — BLOCKS 1, 2, 4, 5, 7 — is untouched.

---

## 11. Status per step

| step | statement | status |
|---|---|---|
| §2 | `‖∇u‖_op = ‖∇₅b‖_op ≤ 2|a| + r|∇a| + |ω^θ|`, `div₅b = 2a`, the `E`-block | **PROVED** (sympy residual `0`; `4000`-sample envelope check, worst ratio `0.998932`) |
| §3 | `|ω^θ| ≤ λM ≤ (3/2)M` on the window | **PROVED given the contradiction hypothesis** of (S3) |
| §3 | brief-layer correction: `‖η₀‖_∞ = M/(ρ₀ sin δ)`, not `M/ρ₀` | **PROVED** (`u2`) |
| §4 L4.1 | `b(0)=0`, `|b(x)| ≤ Γ|x|` | **PROVED** (sympy) |
| §4 L4.2 | `x·b ≤ Γ_rad|x|²`, closed form for `Λ_max(sym∇₅b)` | **PROVED** (rel. err `1.0·10^{-15}`, `2·10⁵` samples) |
| §4 C4.3 | the envelope (4.1) for `Γ_rad` | **PROVED** (`0` violations, `2·10⁵` samples; attained) |
| §5 L5.0 | `D_t∂_jη = −(∂_jb)·∇η + νΔ₅∂_jη`; **no** zeroth-order term from `div₅b = 2a` | **PROVED** (sympy residual `0`, two items) |
| §5 L5.1 | the two weight rearrangements, both with the good sign `−2/ρ²` | **EXACT** (sympy residual `0`, four items) |
| §5 P1 | `sup ρ|η| ≤ e^{c_R}E₀M` | **PROVED-modulo-(standard approximate-maximum argument)** |
| §5 P2 | `sup ρ²|∇η| ≤ e^{pc_G}𝔊₀M`, `p ≤ 3` | **PROVED-modulo-(same)** |
| §5 | P1, P2 exactly sharp on the reference strain | **PROVED** (rel. err `0.0`, both, `λ = 1.25, 1.5`) |
| §5.4(b) | the brief's pushforward inequality is FALSE for `ν > 0` | **PROVED** (support argument) |
| §5.4(c) | `TV(η(s)) ≤ e^{3c_G}TV(η₀)`; diffusion costs nothing (`C = 0`) | **PROVED** |
| §6 | `∇a = K*∇₅η`; Riesz composition `π⁴/2`; `r|∇a| ≤ (3π²/16)e^{pc_G}𝔊₀M sin φ` | **PROVED** (constant checked two ways, rel. `3.2·10^{-12}`) |
| §7(a) | `κ(ρ',s) ≥ 0` hence `a_far(0,s) ≤ a(0,s)` hence `C₁ = 0` | **PROVED** (Lemma 2) |
| §7(b) | the axis-safe collar `3.9992184 e^{c_R}E₀M` | **PROVED** |
| §7 | the three far/near constants re-derived | **PROVED** (`0.291999 / 0.014754 / 3.999218`) |
| §8 | `Γ ≤ 2a(0,s) + C″M` on `ℝ⁵`; the fixed point and `L_*` | **PROVED** given the rows above |
| §9 | the numerical confrontation | **MEASURED** (numerics falsify, they never prove) |
| §10 | Theorem A off the critical path for BLOCK 3 | **PROVED** (by exhibition: the chain does not use it) |

---

## 12. What a referee should attack

1. **The approximate-maximum step in P1/P2.** `V = ρ|η|` and `W = ρ²|∇η|` are Lipschitz,
   not `C²`, at their zero sets, and `ρ` is not smooth at the origin. The standard fix
   (mollify `|·|` to `√(·²+ε²)`, work on `ℝ⁵∖B_ε`, use that `V(0)=W(0)=0`, let `ε→0`) is
   routine but is not written out here. This is the single "modulo" in the chain.
2. **`c_R` versus `c_G` in the collar.** §7(b) uses `e^{c_R}` where `u3` computes with
   `e^{c_G} ≥ e^{c_R}` — conservative, so the tables are majorants, but a referee should
   check the direction of every such substitution.
3. **`𝔊₀ = 20.9197` is attained at `φ = 6.660°`, in the taper transition**, not at the
   equatorial layer. A datum with a wider taper would lower `𝔊₀` but raise `E₀ = 1/sin δ`
   in (7.2) and lower `κ_δ` in the clock. The optimum over `δ` is not computed here.
4. **The far/near split is applied at every `x ∈ ℝ⁵`**, including `ρ < ρ₀` (where the inner
   region is empty) and `ρ > R/2` (where the far region is empty). The bounds are one-sided
   majorants over each region, so empty regions contribute `0`; a referee should confirm
   there is no hidden lower bound in the far/near lemma's statement.
5. **The fixed point is in `(Ĝ, c_G, p)` with `a(0,s) ≤ (λ/2)ML`.** Using the *lower* bound
   on `a(0,s)` from `B(t)`(iii) instead would tighten `p`; not done.
6. **`p` is bounded via `Γ_rad ≤ Γ̄` in the worst case.** The value `p(10⁴) = 2.107`
   depends on the same `Ĝ` it produces; the iteration is monotone but a referee should
   check that the fixed point is the *smallest* one (bisection in `L` finds the existence
   boundary, not the branch).

---

## 13. Files

| file | what |
|---|---|
| `u1_algebra.py`, `u1_results.json` | the exact algebra: `(2.1)`, `div₅b=2a`, `b(0)=0`, `(5.1)`, `(5.2)`, `C_K|S⁴|=1`, the Riesz constant `π⁴/2` two ways |
| `u2_datum.py`, `u2_results.json` | `E₀`, `𝔊₀`, `𝔉₀` for datum (D-C); `𝒥_ac` (reproduces `hk2`'s `65.625910`) |
| `u3_assemble.py`, `u3_results.json` | far/near constants re-derived; `Ĉ_a`; `C″`; the window fixed point; all tables |
| `u4_numerics.py`, `u4_results.json`, `u4_log.txt` | the independent kernel instrument and the confrontation |
| `u5_radial_rate.py`, `u5_results.json` | `Λ_max(sym∇₅b)` closed form, the `Γ_rad` envelope, P1/P2 sharpness |
| `u6_kappa.py`, `u6_results.json` | the `l=1` strain rate `κ` feature by feature (§9.1) |
| `u7_mutate.py`, `u7_results.json`, `u7_mutate_log.txt` | blind mutation test of the gate |
| `check_constants.py` | re-asserts every displayed number from the JSONs |
| `SHA256SUMS` | `shasum -a 256` over all of the above |

---

## 14. The gate (FL-043 / L-89 discipline)

```
  python3 u1_algebra.py      # exact algebra + Riesz constant           (~18 s)
  python3 u2_datum.py        # E0, G0, F0, J_ac                         (~10 s)
  python3 u3_assemble.py     # far/near constants, C'', the fixed point (~20 s)
  python3 u4_numerics.py     # the confrontation                        (~40 min)
  python3 u5_radial_rate.py  # Gamma_rad, sharpness                     (~45 s)
  python3 u6_kappa.py        # kappa feature by feature                 (~1 s)
  python3 check_constants.py # -> 287 checks run / ALL CHECKS PASS
  python3 u7_mutate.py       # blind mutation of every JSON leaf
  shasum -a 256 -c SHA256SUMS
```

**Mutation result.** `2300` numeric leaves across the six results files, each perturbed by a
relative `10^{-3}` (absolute `10^{-6}` when the leaf is `0`, negation when boolean), one at a
time, with `check_constants.py` re-run against the mutated copy: **`245` caught, `2055`
blind (`10.7 %` caught)**. The blind set is characterised in `u7_mutate_log.txt` and
consists entirely of leaves this document does not display:

| blind leaves | what they are |
|---|---|
| `1152` | `u4`'s `120` per-point rows (the document displays the per-configuration maxima, which are checked) |
| `550` | the `32` of `u3`'s `45` fixed-point rows this document does not print |
| `144` | `u3`'s `Linf` variants (`p3`/`map` at `σ_*` values) not printed |
| `81` | `u3`'s `headline` block, a duplicate of table rows that are checked |
| `22` + `12` + `94` | `at_boundary` auxiliary fields, `u2`'s strain control (superseded by `u5`), argmax locations, `L=40` scan auxiliaries, and intermediates whose *derived* value is checked (e.g. `a0_control/quadrature`, whose `rel` is checked) |

**No displayed constant is blind.** That is the claim the mutation test is here to make
checkable, and it is a claim about *this* document: a number printed above that no line of
`check_constants.py` re-asserts would show up in the blind list under a name this table does
not cover. `u7_mutate_log.txt` lists all `2055` by path.

**Independence declaration (L-14).** Read as prose and *not* imported, executed or grepped
by any script here: `s3close/assembly/ASSEMBLY.md`, `s3close/hk2/PROOF.md`,
`write/L3v-and-gamma-bound/PROOF.md`, `write/refute-L3v-and-gamma-bound/PROOF.md`,
`rebuild/far-near-kernel-lemma/NOTE.md`, `gaps/gap-V-aronson/NOTE.md`,
`lower/prove-lagrangian/NOTE.md`. **Not read at all:** `s3close/round2/u1/` (leak guard),
`s3close/hk2-b/`, `s3close/round2/refute-*/`, and every Astra output. **Registered
deviation (L-09), 2026-09-08:** after all files here were written and hashed, `ls` was run
on `round2/u1/` once, to confirm that nothing had been written there by this seat. That
returned file *names* only (no content was opened, and no script here reads that path).
Compensating control: none available after the fact; the leak, if any, is the shape of the
other seat's script list, which post-dates every number in this note. Every constant that
appears in both this note and a neighbouring seat (`0.291999`, `0.014754`, `3.999218`,
`π/8`, `8.697888`, `65.625910`, `0.4997212305210886`, `C_K`, `|S⁴|`) was recomputed here
from its own definition; the neighbouring value is quoted only in a comparison column.

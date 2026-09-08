# LEMMA T′ — the scale-invariant Lipschitz stability of the axis-strain functional
# against a SHELL-DEPENDENT strain

Seat `write/lemma-T-shell-dependent`, DTC-2026-09-06.
Task: write out repair **(R1)** of `gaps/refute-gap-T-lipschitz` §4 — state and prove Lemma T
against the shell-dependent reference `Λ(x) = T_{λ(|x|)}x`, derive the constant, and check every
number against the refuter's `r2` figures.

Every number below came out of a script in **this** folder that I wrote and ran
(`p1`…`p8`; `p1`–`p5`,`p7` re-asserted by `check_constants.py`, `p6` measuring that gate’s
mutation coverage and `p8` auditing this document against the JSONs); `SHA256SUMS` is computed,
never typed.
Nothing outside this folder was written; no other seat's files were modified. The refuter's
`r2_results.json` is **read** (read-only) by `p3_constants.py` so that the comparison is against
their stored numbers rather than against a transcription of them.

The algebra marked PROVED is sympy-exact (residual `0`) or a written argument reproduced here in
full; the numerics are refuters — they falsify, they never prove.

---

## 0. Headline, and what it is not

**Proved here.** With `Λ(x) = T_{λ(|x|)}x`, `λ : [ρ₀,R] → [1,3/2]` of class `C¹` with
`|ρλ′(ρ)/λ(ρ)| ≤ κ_s/L`, and with the **reference value**

```
a_ref[η₀;λ]  :=  ∫_S 𝒦(Λx) η₀(x) λ(|x|)² dx₅          ( = (M/2)∫_{ρ₀}^R P_h(λ(ρ)) dlog ρ
                                                          = (M/2) A   for the bang-bang cap )
A  :=  ∫_{ρ₀}^R λ(ρ) dlog ρ ,
```

for every `C¹` diffeomorphism `Φ` of the shell obeying `|Φ(x) − Λ(x)| ≤ μ|Λ(x)|` (`μ<1`) and
`|J_Φ(x) − λ(|x|)²| ≤ μ_J λ(|x|)²`, and every datum with `|η₀| ≤ M/r` on `S`:

> **(T′)**  `| a[η₀∘Φ⁻¹](0) − a_ref[η₀;λ] |  ≤  π M A [ 3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8 ]`

and, for the bang-bang cap (where `a_ref = (M/2)A` exactly),

> **(T′-rel)**  `|Δ| / a_ref  ≤  2π[ 3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8 ]  →  (15π/4)μ = 11.78097245 μ`
> as `μ, μ_J → 0` with `μ_J ≤ μ`, **free of `λ`, of `L`, of `ρ₀/R`, and of the profile `λ(·)`.**

The constants are **the same** as the single-`T_λ` Lemma T of `gaps/gap-T-lipschitz`
(`15π/8 = 5.890486225` absolute, `15π/4 = 11.78097245` relative). What changes is *what they are
constants of*: the anisotropy weight `λ M L` of the single-map lemma is replaced by the
**strain integral** `M A = M ∫λ dlog ρ`, and that replacement is exactly what makes the
shell-dependent map admissible with `μ = 0`.

**The instance the campaign needs (`Φ = Λ`, the integro-ODE profile).**
`μ = 0` identically; the *only* hypothesis constant that survives is the inter-shell shear, and it
is computed, not fitted:

```
κ_s   = sup_ρ L|ρλ′/λ|          = 2(√(3/2) − 1) = √6 − 2     = 0.4494897427831781
μ_J   = sup_x |J_Λ/λ² − 1|      = 2κ_s/L        = (2√6−4)/L  = 0.8989794855663562 / L
⇒  |Δ|/a_ref  ≤  (3π/4)μ_J  =  (3π/2)κ_s/L  =  3π(√6−2)/2 / L  =  2.1181705106873974 / L .
```

Measured true error: `0.34913942333376546/L` (own quadrature, two independent routes agreeing to
a relative `1.35e-13`; §4) — the proved constant is conservative by the factor **6.0668**.

**What this is NOT.** Lemma T′ is a statement about the **kernel and two maps**. It bounds the
*sensitivity of the strain functional* to a perturbation of the flow map. It says **nothing**
about the Navier–Stokes flow map. The identification of the true flow map with `Λ` to relative
`O(c/L)` — i.e. producing the `μ` and the `μ_J` that (T′) consumes — is the **separate** piece
`L3v` + `Γ` of `lower/prove-lagrangian` §4(2)–(3), untouched here and untouched by every seat of
this round. GAP T therefore remains a gap; what closes with this note is the half the refuter
labelled "cheap, verified numerically here, unwritten" (`refute-gap-T-lipschitz` §8, item 1).

---

## 1. Setting, notation, and the kernel convention

5-D lift of axisymmetric no-swirl Navier–Stokes (`lower/prove-lagrangian` §1,
`rebuild/far-near-kernel-lemma` §0). `x = (y,z)`, `y ∈ ℝ⁴`, `r = |y|`, `ρ = |x|`, `φ` the polar
angle from `+z` (`z = ρ cos φ`, `r = ρ sin φ`). `η = ω^θ/r`, `−Δ₅ψ₁ = η`, `a = u^r/r = −∂_zψ₁`,
`D_tη = νΔ₅η`.

```
dx₅ = 2π² ρ⁴ sin³φ  dρ dφ                            (|S³| = 2π², dy = r³dr dS³, dr dz = ρ dρ dφ)
a(x) = ∫ K(x−x′)η(x′)dx′ ,   K(w) = 3 w_z/(8π²|w|⁵)  (far-near-kernel-lemma §0)
⇒  a[f](0) = ∫_{ℝ⁵} 𝒦(x) f(x) dx ,   𝒦(x) := K(−x) = −(3/(8π²)) x_z/|x|⁵ ,   homogeneous of degree −4.
```

> **Correction to the task brief.** The brief writes the origin kernel as `K(w)=3w_z/(8π²|w|⁵)`.
> That is `K`, not `𝒦 = K(−·)`; used at the origin it reverses the sign of `a(0)`. `p1_kernel.py`
> §4 checks the convention by producing `a(0) = +(M/2)L` for the campaign datum
> `η₀ = −M sgn(z)/r` (sympy, exact). The gradient bound the brief quotes is unaffected
> (`|∇𝒦| = |∇K|`).

**Shell and datum.** `S = {ρ₀ < |x| < R}`, `L = log(R/ρ₀)`, `σ = log(ρ/ρ₀)/L ∈ [0,1]`.

```
(D)   η₀ measurable, supp η₀ ⊂ S,  |η₀(x)| ≤ M/r(x)     (i.e. |ω₀^θ| ≤ M on S).
```
That is all that is used. Oddness in `z` is not needed; no regularity beyond measurability is
needed. The campaign datum `ω₀^θ = −M sgn(z)h_δ(φ)` satisfies (D) with equality where `h_δ = 1`.

**Reference map.** `T_λ(y,z) = (λy, λ^{-2}z)` (volume-preserving in ℝ³; `J₅ = λ²`), and

```
Λ(x) := T_{λ(ρ)} x = (λ(ρ)y, λ(ρ)^{-2}z) ,      λ ∈ C¹([ρ₀,R]; [1,3/2]) ,
(S)   |ρ λ′(ρ)/λ(ρ)| ≤ κ_s/L      for all ρ ∈ [ρ₀,R]        (the inter-shell shear).
```

Two elementary quantities used throughout:

```
g(φ;λ)  := |T_λ x|/|x|      = ( λ² sin²φ + λ^{-4} cos²φ )^{1/2}
ĝ(ψ;λ)  := |T_λ^{-1} w|/|w| = ( λ^{-2} sin²ψ + λ⁴ cos²ψ )^{1/2}
A       := ∫_{ρ₀}^{R} λ(ρ) dlog ρ = L ∫₀¹ λ dσ            (the strain integral)
```

**Reference functional (this is the object step P of `prove-lagrangian` consumes).**

```
a_ref[η₀;λ] := ∫_S 𝒦(Λx) η₀(x) λ(ρ)² dx₅ .
```
It is `a[η₀∘Λ⁻¹](0)` with the Jacobian replaced by its **shear-free** value `λ(ρ)²`, i.e. the
value obtained by pretending each shell is strained by its own `T_{λ}` with no coupling between
shells. `prove-lagrangian`'s integro-ODE (4.1) is closed on exactly this quantity: its right-hand
side `κH(σ,θ) = κ∫_σ¹ e^{F}dσ′` is a per-shell sum of `Φ(λ)=λ` contributions with **no** shear
term. Lemma T′ is what prices the difference.

---

## 2. LEMMA T′

> **LEMMA T′.** Let `λ ∈ C¹([ρ₀,R];[1,3/2])` satisfy (S) with `2κ_s/L < 1`, let `η₀` satisfy (D),
> and let `Φ : S → ℝ⁵` satisfy
>
> ```
> (H1)  Φ is a C¹ diffeomorphism onto its image, J_Φ := det DΦ > 0 ;
> (H2)  |Φ(x) − Λ(x)| ≤ μ |Λ(x)|     for all x ∈ S,   0 ≤ μ < 1 ;
> (H3)  |J_Φ(x) − λ(ρ)²| ≤ μ_J λ(ρ)²  for all x ∈ S,   μ_J ≥ 0 .
> ```
>
> (No constraint on the *range* of `λ` is used in (T′) itself; `λ ∈ [1,3/2]` enters only in
> Step 0's condition `2κ_s/L < 1` and in Corollary 2's `r_h`.)
>
> Then `a[η₀∘Φ⁻¹](0)` and `a_ref[η₀;λ]` both converge absolutely, with
> `|a_ref| ≤ (3π/8) M A` and `|a[η₀∘Φ⁻¹](0)| ≤ (3π/8)(1+μ_J)(1−μ)^{-4} M A`, and
>
> ```
> | a[η₀∘Φ⁻¹](0) − a_ref[η₀;λ] |  ≤  π M A [ 3μ(1+μ_J)/(2(1−μ)⁵)  +  3μ_J/8 ] .        (T′)
> ```
>
> In particular, if `μ_J ≤ μ ≤ μ_*` then `|Δ| ≤ C(μ_*) μ M A` with
> `C(μ_*) = π[3(1+μ_*)/(2(1−μ_*)⁵) + 3/8]`, `C(0+) = 15π/8 = 5.890486225`.

| `μ_*` | 0 | 0.01 | 0.05 | 0.10 | 0.20 | 0.25 | 0.50 |
|---|---|---|---|---|---|---|---|
| `C(μ_*)` | **5.890486** | 6.182895 | 7.572683 | 9.956617 | 18.435381 | 26.000558 | 227.372768 |

(`p3`/`p4` recompute this table from `π[3(1+μ)/(2(1−μ)⁵)+3/8]`; it is the same function as
`gap-T-lipschitz` §2, and the two tables agree.)

> **VARIANT (T′_Λ) — the literal form of the task brief. PROVED.**
> If instead the hypotheses are stated **against `Λ` itself**, i.e.
> `|Φ − Λ| ≤ μ|Λ|` and `|J_Φ − J_Λ| ≤ μ_J J_Λ`, then, using
> `J_Λ ≤ (1+2κ_s/L)λ²` from (2.2) in Steps 4–5 in place of (H3),
>
> ```
> | a[η₀∘Φ⁻¹](0) − a[η₀∘Λ⁻¹](0) |  ≤  C(μ,μ_J) · M · L ,
> C(μ,μ_J) = π λ̄ (1 + 2κ_s/L) [ 3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8 ] ,   λ̄ := A/L = ∫₀¹λ dσ ≤ 3/2 .
> ```
> As `μ_J ≤ μ → 0` and `L → ∞`, `C → (15π/8)λ̄ μ`; for the integro-ODE profile
> `λ̄ = √(3/2)` and `(15π/8)λ̄ = 15√6 π/16 = 7.21434279466049`. This is the estimate GAP T names verbatim. It is
> **vacuous when `Φ = Λ`** (both sides are `0`), which is precisely why the campaign needs the
> `a_ref` form: what step P consumes is the shear-free strain integral, not `a[η₀∘Λ⁻¹](0)`.
> Verified numerically in §5 on 12 non-trivial maps (worst slack `37.90`).

### Proof

**Step 0 (Λ is a `C¹` diffeomorphism of `S` onto its image). PROVED.**
*(a) Positive Jacobian.* `Λ` is `SO(4)`-equivariant, so on the meridian half-plane it is
`(r,z) ↦ (λ(ρ)r, λ(ρ)^{-2}z)`, and the three `S³`-directions orthogonal to `y` are scaled by
`λ`. Hence `J_Λ = λ³ · det ∂(r′,z′)/∂(r,z)`, and (`p1_kernel.py` §2, sympy, residual `0`, by two
independent routes — autodifferentiation of the full 5×5 map on the family `λ = aρ^m`, which
realises every pointwise pair `(λ, ρλ′/λ)`, and a hand-built matrix with free `(λ,λ′)`):

```
J_Λ(x) = det DΛ(x) = λ(ρ)² [ 1 + (ρλ′(ρ)/λ(ρ)) (sin²φ − 2cos²φ) ] .                    (2.1)
```
Since `sin²φ − 2cos²φ = 1 − 3cos²φ ∈ [−2,1]` (sympy identity, residual `0`), (S) gives

```
λ² (1 − 2κ_s/L)  ≤  J_Λ  ≤  λ² (1 + 2κ_s/L) ,        so   J_Λ > 0  when  2κ_s/L < 1.       (2.2)
```
*(b) Injectivity.* Let `w = Λ(x)`, let `ψ` be the polar angle of `w` and `ρ = |x|`. Then
`x = T_{λ(ρ)}^{-1}w`, so `ρ = |w| ĝ(ψ;λ(ρ))`, i.e. with `u = log ρ`,
`u = log|w| + log ĝ(ψ; λ(e^u))`. Differentiating and using
`∂ log ĝ/∂ log λ = (−λ^{-2}sin²ψ + 2λ⁴cos²ψ)/ĝ²` (sympy, residual `0`), a convex combination of
`−1` and `+2`, hence `|∂ log ĝ/∂ log λ| ≤ 2`, the map `u ↦ log|w| + log ĝ(ψ;λ(e^u))` is
Lipschitz on `[log ρ₀, log R]` with constant `≤ 2κ_s/L < 1`: a contraction, so it has **at most
one** fixed point. Thus `ρ` — and then `x = T_{λ(ρ)}^{-1}w` — is determined by `w`. `Λ` is
injective; with (a) and the inverse function theorem it is a `C¹` diffeomorphism onto its image. ∎
*Condition:* `L > 2κ_s`. For the integro-ODE profile `2κ_s = 0.8989794856`, so this is satisfied
at every `L` the campaign uses (`L ≥ 8.3`).

**Step 1 (change of variables). PROVED.**
`η_Φ := η₀∘Φ⁻¹` on `Φ(S)`, `0` outside (this is the definition of the transported field the
statement uses; nothing is asserted about how a fluid would transport it). By (H1),
`a[η_Φ](0) = ∫_{Φ(S)}𝒦(w)η₀(Φ^{-1}w)dw = ∫_S 𝒦(Φ(x))η₀(x)J_Φ(x)dx`. Subtracting the definition
of `a_ref` and adding and subtracting `𝒦(Λx)J_Φ(x)`,

```
Δ := a[η₀∘Φ⁻¹](0) − a_ref
   = ∫_S η₀(x) [ (𝒦(Φx) − 𝒦(Λx)) J_Φ(x)  +  𝒦(Λx)(J_Φ(x) − λ(ρ)²) ] dx .              (2.3)
```
Absolute convergence of every integral written is Step 6.

**Step 2 (kernel bounds, sharp). PROVED** (`p1_kernel.py` §1, sympy, residual `0`).
`|𝒦(w)| = (3/(8π²))|cos φ_w| |w|^{-4} ≤ (3/(8π²))|w|^{-4}`, sharp on the axis; and
`∇𝒦(w) = −(3/(8π²))|w|^{-5}(e_z − 5(w_z/|w|)ŵ)` componentwise in `ℝ⁵`, with
`|e_z − 5c ŵ|² = 1 + 15c² ≤ 16` (`c = w_z/|w|`), so

```
sup|w|⁴|𝒦(w)|   = 3/(8π²) = 0.037995443865876666 ,
sup|w|⁵|∇𝒦(w)| = 3/(2π²) = 0.15198177546350666 ,   both attained on the axis.           (2.4)
```

**Step 3 (the shell integral against the SHELL-DEPENDENT reference — this is the new step).
PROVED, exact.**
The identity that makes it work (`p2_identities.py`, sympy exact + 40-digit quadrature,
max relative error `0.0`):

```
J(λ) := ∫₀^π sin²φ / g(φ;λ)⁴ dφ  =  π/(2λ)      EXACTLY, for every λ > 0.               (2.5)
```
*Proof of (2.5).* `d/dφ[ arctan(√(a/b) tan φ)/√(ab) ] = 1/(b cos²φ + a sin²φ)` (sympy, residual
`0`), so `∫₀^π dφ/(b cos²φ + a sin²φ) = π/√(ab)` for `a,b>0` (two half-periods). Differentiating
in `a` under the integral, `∫₀^π sin²φ/(b cos²φ + a sin²φ)² dφ = (π/2)a^{-3/2}b^{-1/2}`. Put
`a = λ², b = λ^{-4}`; then `b cos²φ + a sin²φ = g(φ;λ)²` and the right side is
`(π/2)λ^{-3}λ² = π/(2λ)`. ∎

Now, **with no change of variables at all** — this is why the shell-dependence costs nothing —
write the reference weight `λ(ρ)²`, use (D), `|Λx| = ρ g(φ;λ(ρ))` and `dx₅ = 2π²ρ⁴sin³φ dρdφ`:

```
I := ∫_S |η₀(x)| λ(ρ)² |Λx|^{-4} dx₅
  ≤ ∫_{ρ₀}^{R}∫₀^π  (M/(ρ sinφ)) λ(ρ)² ρ^{-4} g(φ;λ(ρ))^{-4} · 2π² ρ⁴ sin³φ  dφ dρ
  =  2π² M ∫_{ρ₀}^{R} λ(ρ)² [ ∫₀^π sin²φ / g(φ;λ(ρ))⁴ dφ ] dlog ρ
  =  2π² M ∫ λ(ρ)² · π/(2λ(ρ)) dlog ρ
  =  π³ M ∫_{ρ₀}^{R} λ(ρ) dlog ρ  =  π³ M A ,      with EQUALITY iff |ω₀^θ| = M a.e. on S. (2.6)
```
The `φ`-integral is performed at each fixed `ρ` with `λ = λ(ρ)`, so a `ρ`-dependent `λ` is no
obstruction whatsoever: the anisotropy contributes exactly the factor `1/λ(ρ)`, which cancels one
power of the weight `λ(ρ)²` and leaves the **strain density `λ(ρ)` per e-fold** — the same
`λ`-per-e-fold that `far-near-kernel-lemma` §1 Consequence A and `prove-lagrangian` Lemma 1
produce for the value itself. Numerically re-checked in `p4_verify.py` (§0(iii)): max relative
error `1.32e-15` at `L = 8.318 … 400`.

*(Relation to the single-map lemma.* `gap-T-lipschitz` Step 3 bounds the **unweighted**
`I_λ = ∫_S|η₀||T_λx|^{-4}dx₅ ≤ π³ML/λ` and carries the factor `λ²` outside. Mine carries the
weight inside because the weight is now `ρ`-dependent. For constant `λ` the two agree exactly:
`λ² · (π³ML/λ) = π³MλL = π³MA`. So (2.6) contains the single-map step as the special case
`λ ≡ const`.*)

**Step 4 (the kernel difference). PROVED.**
By (H2), every point `w` of the segment `[Λx, Φ(x)]` satisfies
`|w| ≥ |Λx| − |Φ(x) − Λx| ≥ (1−μ)|Λx| > 0`, so the segment misses the singularity of `𝒦` and the
fundamental theorem of calculus along it, with (2.4), gives

```
|𝒦(Φx) − 𝒦(Λx)| ≤ |Φx − Λx| · sup_{[Λx,Φx]}|∇𝒦| ≤ μ|Λx| · (3/(2π²))((1−μ)|Λx|)^{-5}
                 = (3μ / (2π²(1−μ)⁵)) |Λx|^{-4} .                                        (2.7)
```
With `J_Φ ≤ (1+μ_J)λ(ρ)²` from (H3) and then (2.6),

```
|term I| ≤ (3μ(1+μ_J)/(2π²(1−μ)⁵)) ∫_S|η₀|λ²|Λx|^{-4}dx₅ ≤ (3π μ(1+μ_J)/(2(1−μ)⁵)) M A .  (2.8)
```

**Step 5 (the Jacobian term). PROVED.**
`|𝒦(Λx)| ≤ (3/(8π²))|Λx|^{-4}` and `|J_Φ − λ²| ≤ μ_J λ²`, so by (2.6) again

```
|term II| ≤ (3μ_J/(8π²)) ∫_S|η₀|λ²|Λx|^{-4}dx₅ ≤ (3π μ_J/8) M A .                        (2.9)
```

**Step 6 (assembly and absolute convergence). PROVED.**
Adding (2.8) and (2.9) in (2.3) gives (T′). Absolute convergence is the same computation with
`|𝒦| ≤ (3/(8π²))|·|^{-4}`: `|a_ref| ≤ (3/(8π²))·π³MA = (3π/8)MA`, and, using
`|Φx| ≥ (1−μ)|Λx|` and `J_Φ ≤ (1+μ_J)λ²`,
`|a[η₀∘Φ⁻¹](0)| ≤ (3π/8)(1+μ_J)(1−μ)^{-4} M A`. ∎

### Corollary 1 (relative form). PROVED.
For the bang-bang cap `ω₀^θ = −M sgn(z)` the inequality (D) is an equality, and
(`p2_identities.py`, exact) `∫₀^{π/2}cos φ sin²φ/g(φ;λ)⁵dφ = λ/3`, so

```
a_ref = (3M/4)∫dlog ρ ∫₀^π |cos φ| sin²φ / g(φ;λ(ρ))⁵ dφ = (M/2)∫λ(ρ)dlog ρ = (M/2) A ,   (2.10)
```
whence
```
|Δ|/a_ref ≤ 2π[ 3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8 ]  →  (15π/4)μ = 11.78097245 μ  (μ_J ≤ μ → 0),
```
**independent of `λ`, of `L`, of `ρ₀/R`, and of the profile `λ(·)`.**

### Corollary 2 (the δ-tapered datum). PROVED.
For `ω₀^θ = −M sgn(z)h_δ(φ)`, (2.10) reads `a_ref = (M/2)∫P_{h_δ}(λ(ρ))dlog ρ` with
`P_h(λ) = 3∫₀¹h(arcsin v)v²(Av²+B)^{-5/2}dv`, `A = λ²−λ^{-4}`, `B = λ^{-4}`. Since `h ≤ 1`,
(2.6) still holds with the same `A`, so with
`r_h := inf_{λ∈[1,3/2]} P_{h_δ}(λ)/λ`,

```
|Δ|/a_ref ≤ (2π/r_h)[ 3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8 ] .
```
Computed here (`p2`, `mpmath`, 25 dps): `r_h = 0.998660 / 0.994097 / 0.981822 / 0.961678 /
0.902480 / 0.830359 / 0.688268` at `δ = 3/5/7.5/10/15/20/30°` (the infimum is at `λ = 3/2` in
every case). So the taper costs **≤ 1.9 % on the constant at `δ ≤ 7.5°`**, the range
`prove-lagrangian` runs.
*Instrument cross-check:* `P_{h_δ}(1)/2 = κ_δ(7.5°) = 0.4997212305210993` here,
`0.4997212305210992` (refuter `r1`), `0.4997212305210886` (`gap-T-lipschitz` `s3`),
`0.4997212` (`prove-lagrangian` §2) — four independent evaluations.
*Unresolved cross-check (recorded, not load-bearing):* my `P_{h_δ}(3/2)` is
`1.497990 / 1.491146 / 1.472733 / 1.442517 / 1.353721 / 1.245539 / 1.032403` against the
synthesis's `1.4980 / 1.4914 / 1.4736 / 1.4444 / 1.3597 / 1.2584 / 1.0674`. The first two agree to
`4–5` digits and the gap widens with `δ` (`3.4 %` at `30°`); the two computations use different
routes (direct quadrature here, Gegenbauer projection of the deformed field there) and the
discrepancy is **not** resolved in this note. It does not touch Lemma T′, whose `A` uses `λ`, not
`P_h`, and it only makes `r_h` slightly pessimistic if theirs is right.

### Corollary 3 (composition — how the PDE half plugs in). PROVED.
Suppose the PDE half delivers the hypotheses **against `Λ`** rather than against the shear-free
weight, i.e.
```
(H2_Λ)  |Φ − Λ| ≤ μ|Λ| ,        (H3_Λ)  |J_Φ − J_Λ| ≤ μ_J^Λ J_Λ .
```
By (2.2), `|J_Λ − λ²| ≤ (2κ_s/L)λ²` and `J_Λ ≤ (1+2κ_s/L)λ²`, so
`|J_Φ − λ²| ≤ μ_J^Λ J_Λ + |J_Λ − λ²| ≤ [ μ_J^Λ(1 + 2κ_s/L) + 2κ_s/L ] λ²`, i.e. (H3) holds with

```
μ̃_J  =  μ_J^Λ (1 + 2κ_s/L)  +  2κ_s/L ,                                                (2.11)
```
and (T′) applies with `(μ, μ̃_J)`. Hence: *if* the PDE half gives `μ = γ c/L` and
`μ_J^Λ = γ′ c/L`, then

```
|a[η₀∘Φ⁻¹](0) − a_ref| / a_ref  ≤  (15π/4)(γ c/L) + (3π/4)(γ′ c/L + 2κ_s/L) + O(1/L²) = O(1/L),
```
which is the order `prove-lagrangian` already carries. **This is the only sentence that connects
Lemma T′ to GAP T, and its hypothesis — the PDE half — is not proved anywhere.**

### Corollary 4 (the instance: `Φ = Λ`, the integro-ODE profile). PROVED + measured.
`Φ = Λ` gives `μ = 0` exactly, and by (2.2), `μ̃_J = 2κ_s/L` **and this is attained** (at
`σ` where `|dlogλ/dσ|` is maximal and `φ ∈ {0,π}`). So

```
|a[η₀∘Λ⁻¹](0) − a_ref| ≤ (3π/8)(2κ_s/L) M A = (3π/4)(κ_s/L) M A ,
|a[η₀∘Λ⁻¹](0) − a_ref| / a_ref ≤ (3π/2) κ_s / L .                                        (2.12)
```

---

## 3. The integro-ODE instance: every constant, derived

`prove-lagrangian` (4.1): `λ(σ,θ) = (1 − (1−σ)κθ/2)^{-2}`, `σ = log(ρ/ρ₀)/L`; terminal time
`κθ = 2(1 − √(2/3))` (so that `λ(0) = 3/2`). All of the following are exact (sympy) and are
re-asserted by `check_constants.py` (`p3_constants.py`):

| quantity | exact | value |
|---|---|---|
| `κθ` (terminal) | `2 − 2√6/3` | `0.36700683814454793` |
| `λ(σ=0)` , `λ(σ=1)` | `3/2` , `1` | `1.5` , `1.0` |
| `λ̄ := ∫₀¹λ dσ` (so `A = λ̄ L`) | `√6/2 = √(3/2)` | `1.224744871391589` |
| `∫₀¹dσ/λ` | `5/9 + √6/9` | `0.8277210825314643` |
| `κ_s := sup_ρ L\|ρλ′/λ\| = sup_σ\|dlogλ/dσ\|` | `√6 − 2 = 2(√(3/2)−1)` | `0.4494897427831781` |
| `μ_J · L = 2κ_s` | `2√6 − 4` | `0.8989794855663562` |
| `C_rel := (3π/2)κ_s` | `3π(√6−2)/2` | `2.1181705106873974` |

*Derivation of `κ_s`.* `dlogλ/dσ = −κθ/(1 − (1−σ)κθ/2)`, monotone, maximal in modulus at `σ=0`
where the denominator is `1 − κθ/2 = √(2/3)`; hence
`κ_s = κθ/√(2/3) = 2(1−√(2/3))·√(3/2) = 2(√(3/2)−1) = √6−2`. (Verified against the maximum over a
2001-point `σ`-grid to `< 1e-13`.)
*Derivation of `μ_J·L`.* By (2.1)–(2.2), `sup|J_Λ/λ²−1| = sup|ρλ′/λ| · sup|sin²φ−2cos²φ| = 2κ_s/L`.
*Derivation of `C_rel`.* (2.12): `(3π/4)·(2κ_s/L) = (3π/2)κ_s/L`.

> **Correction to the task brief.** The brief states "the integro-ODE gives `κ_s ≈ 0.899c`".
> `0.899` is **`μ_J·L = 2κ_s`**, a pure number, not a multiple of the window `c`; and `κ_s` itself
> is `0.4494897428`. If one insists on writing it against `c := 2log(3/2) = 0.8109302162` (the
> total log-radius travel of `prove-lagrangian` §4(2), which is the `c` the refuter used), then
> `μ_J·L = 1.1085781089 c` — the refuter's "constant `1.109 c`". Neither `κ_s` nor `μ_J L` is
> `0.899 c`.

### 3.1 Check against `refute-gap-T-lipschitz` `r2_results.json` (their stored numbers)

`p3_constants.py` reads their JSON directly. `a_ref` here `= (M/2)λ̄L`; `a_true` is my
prediction `a_ref + 0.21380335906432307 M` (§4 below); `bound/a_true` is (2.12) divided by their
`a_true`.

| `L` | `μ_J·L` derived | `μ_J·L` (r2) | `a_true` predicted | `a_true` (r2) | bound/`a_true` derived | bound/`a_true` (r2) | `×L` |
|---|---|---|---|---|---|---|---|
| 8.318 | 0.8989795 | 0.8989795 | 5.3073741 | 5.3073741 | 0.2443976 | 0.2443976 | 2.03284 |
| 20 | 0.8989795 | 0.8989795 | 12.4612521 | 12.4612521 | 0.1040914 | 0.1040914 | 2.08183 |
| 50 | 0.8989795 | 0.8989795 | 30.8324251 | 30.8324251 | 0.0420696 | 0.0420696 | 2.10348 |
| 100 | 0.8989795 | 0.8989795 | 61.4510469 | 61.4510469 | 0.0211080 | 0.0211080 | 2.11080 |
| 400 | 0.8989795 | 0.8989795 | 245.1627776 | 245.1627776 | 0.0052908 | 0.0052908 | 2.11632 |

```
max |μ_J·L derived − r2|        = 8.44e-15
max relative |a_true pred − r2| = 5.67e-15
max relative |bound − r2|       = 8.91e-15
```

So the refuter's measured band **"relative error ≤ 2.03/L … 2.12/L at `L = 8.3 … 400`"** is
exactly `C_rel/L` divided by `a_true/a_ref = 1 + 0.34913942/L`; its limit is the derived constant
`C_rel = 2.1181705106873974`. **Derived, not fitted.** (Against `a_ref` itself the number is
`C_rel/L` at every `L`, to `1e-12`.)

---

## 4. What the shear actually costs (sharpness of Corollary 4)

`a[η₀∘Λ⁻¹](0)` differs from `a_ref` by exactly the Jacobian-shear term of (2.1). Writing
`W = sin²φ − 2cos²φ`, `D(ρ) = dlogλ/dlogρ`, and (`p2`, sympy closed form)

```
Q(λ) := ∫₀^{π/2} cos φ sin²φ (sin²φ − 2cos²φ) / g(φ;λ)⁵ dφ
      = λ(−2λ¹²√(λ⁶−1) + 9λ⁹ asinh√(λ⁶−1) − 8λ⁶√(λ⁶−1) + √(λ⁶−1)) / (3√(λ⁶−1)(λ¹²−2λ⁶+1)) ,
Q(1) = −1/15   (exact limit) ,   Q(3/2) = −0.68120868384 ,
```
one gets, for the bang-bang cap,

```
a[η₀∘Λ⁻¹](0) − a_ref = (3M/2) ∫₀¹ (dlogλ/dσ) Q(λ(σ)) dσ  =  +0.21380335906432307 M ,   (4.1)
```
**independent of `L`** — hence a relative error `0.34913942333376546/L`. Two independent
instruments agree: `p3` (1-D `mpmath` quadrature in `σ` against the closed-form `Q`, giving
`0.21380335906432307 M`) and `p4`/`p5` (2-D Gauss–Legendre over `(σ,φ)` of the full 5-D
integrand, giving `0.213803359064352 M` at `L = 8.318` and stable to a relative spread of
`7.14e-12` over `L = 8.3 … 400`) — relative agreement between the two routes `1.35e-13`. And **three** independent computations of
`a[η₀∘Λ⁻¹](0)` itself agree at `L = 8.317766166719343`: `5.307374086126348` (refuter `r2`),
`5.3073740861` (my closed-form route `a_ref + offset`), `5.3073740861` (my 2-D quadrature) —
`5.7e-15` relative. The refuter's own `r2` measured the offset as `0.3351/L … 0.3488/L`,
converging to `≈ 0.349/L`; the closed form (4.1) is what it converges to.

```
proved  C_rel = 2.1181705         actual 0.3491394        conservatism factor 6.0668 .
```
Note the sign: the shear **raises** the true strain above the shear-free reference, so the
integro-ODE's `a_ref` is a lower bound for `a[η₀∘Λ⁻¹](0)` on this profile. That direction is
favourable for the campaign's lower bound on the strain, but Lemma T′ is deliberately two-sided —
it is the two-sidedness the refuter's R-A1 identified as the real burden.

---

## 5. Numerical verification of (T′) on maps that are NOT `Λ`

Own quadrature (`lib5.py`, `p4_verify.py`): all maps are `SO(4)`-equivariant, so the 5-D Jacobian
of a meridian map is `J₅ = (r′/r)³ det ∂(r′,z′)/∂(r,z)` and stages compose multiplicatively; all
radii are carried as `log ρ` so that `L = 400` does not overflow. `μ` and `μ_J` are **measured** on
a `2001 × 1501` grid, never assumed. `ρ₀ = 1`, `M = 1`, bang-bang cap unless stated.

Instrument cross-checks (`p4` §0): `a[T_λ](0)/(ML) − λ/2` max error `2.07e-14` at
`λ = 0.7,1,1.25,1.5,2`; `det DΛ` against the closed form of `p1`, max relative `0.0` over 500
random points; Step-3 identity (2.6) max relative `1.32e-15`; `κ_δ(7.5°) = 0.4997212305210993`.
Gauss–Legendre `200` vs `400` nodes: `1.6e-14`.

**`Φ = Λ` (`μ = 0`; only the shear enters).**

| `L` | `μ` | `μ_J` | `μ_J·L` | `\|Δ\|/a_ref · L` | bound/`a_ref` · `L` | slack |
|---|---|---|---|---|---|---|
| 8.318 | 0 | 0.108079 | 0.8989795 | 0.349139 | 2.118171 | 6.067 |
| 20 | 0 | 0.044949 | 0.8989795 | 0.349139 | 2.118171 | 6.067 |
| 50 | 0 | 0.017980 | 0.8989795 | 0.349139 | 2.118171 | 6.067 |
| 100 | 0 | 0.008990 | 0.8989795 | 0.349139 | 2.118171 | 6.067 |
| 400 | 0 | 0.002247 | 0.8989795 | 0.349139 | 2.118171 | 6.067 |

**`Φ ≠ Λ`** — `Λ` composed with a radial log-ripple `u ↦ u(1+μ_a sin(kπ log|u|/L))` and/or an
angular shear `φ ↦ φ + β sin²φ`, and one pre-composed case (`shear` applied *before* `Λ`, so the
two do not commute):

| `L` | map | `μ` | `μ_J` | `\|Δ\|` | bound (T′) | slack |
|---|---|---|---|---|---|---|
| 8.32 | `Λ` + ripple(0.02,k=1) | 0.02000 | 0.21142 | 0.28230 | 3.8241 | 13.55 |
| 8.32 | `Λ` + ripple(0.05,k=3) | 0.05000 | 0.41624 | 0.27958 | 9.3886 | 33.58 |
| 8.32 | `Λ` + shear(0.02) | 0.02000 | 0.15410 | 0.21090 | 3.0753 | 14.58 |
| 8.32 | `Λ` + shear(0.08) | 0.07998 | 0.32921 | 0.16758 | 11.6934 | 69.78 |
| 8.32 | `Λ` + ripple + shear | 0.05048 | 0.38851 | 0.21581 | 9.0217 | 41.80 |
| 8.32 | shear + `Λ` (pre-composed) | 0.04999 | 0.18583 | 0.20661 | 5.9082 | 28.60 |
| 50.00 | `Λ` + ripple(0.02,k=1) | 0.02000 | 0.12193 | 0.60490 | 15.9597 | 26.38 |
| 50.00 | `Λ` + ripple(0.05,k=3) | 0.05000 | 0.29858 | 0.54818 | 45.7548 | 83.47 |
| 50.00 | `Λ` + shear(0.02) | 0.02000 | 0.06753 | 0.19603 | 11.6876 | 59.62 |
| 50.00 | `Λ` + shear(0.08) | 0.07998 | 0.23140 | 0.06916 | 59.8098 | 864.80 |
| 50.00 | `Λ` + ripple + shear | 0.05048 | 0.29609 | 0.20244 | 45.8202 | 226.34 |
| 50.00 | shear + `Λ` (pre-composed) | 0.04999 | 0.13518 | 0.16931 | 30.9175 | 182.61 |

**(T′_Λ), the literal brief form** — same 12 maps, hypotheses and difference both taken against
`Λ` (`μ_J` here is `sup|J_Φ/J_Λ − 1|`):

| `L` | map | `μ` | `μ_J` vs `J_Λ` | `\|Δ_Λ\|` | bound (T′_Λ) | slack |
|---|---|---|---|---|---|---|
| 8.32 | `Λ` + ripple(0.02,k=1) | 0.02000 | 0.10439 | 0.06850 | 2.6881 | 39.24 |
| 8.32 | `Λ` + ripple(0.05,k=3) | 0.05000 | 0.28393 | 0.06577 | 8.1892 | 124.51 |
| 8.32 | `Λ` + shear(0.02) | 0.02000 | 0.05075 | 0.00290 | 1.9116 | 658.24 |
| 8.32 | `Λ` + shear(0.08) | 0.07998 | 0.21208 | 0.04623 | 10.6435 | 230.24 |
| 8.32 | `Λ` + ripple + shear | 0.05048 | 0.28080 | 0.00200 | 8.1897 | 4090.12 |
| 8.32 | shear + `Λ` | 0.04999 | 0.12616 | 0.00720 | 5.5482 | 770.86 |
| 50.00 | `Λ` + ripple(0.02,k=1) | 0.02000 | 0.10409 | 0.39110 | 14.8207 | 37.90 |
| 50.00 | `Λ` + ripple(0.05,k=3) | 0.05000 | 0.27650 | 0.33438 | 44.5369 | 133.19 |
| 50.00 | `Λ` + shear(0.02) | 0.02000 | 0.05075 | 0.01778 | 10.5567 | 593.86 |
| 50.00 | `Λ` + shear(0.08) | 0.07998 | 0.21208 | 0.28296 | 58.7784 | 207.72 |
| 50.00 | `Λ` + ripple + shear | 0.05048 | 0.27875 | 0.01136 | 45.0374 | 3963.82 |
| 50.00 | shear + `Λ` | 0.04999 | 0.12910 | 0.04450 | 30.9112 | 694.70 |

**(T′_Λ) holds in every row; worst slack `37.90`.**

**Tapered datum (`7.5°`), with `a_ref` computed as `(M/2)∫P_{h_δ}(λ)dlog ρ` by the same
instrument:** `Φ=Λ` gives `|Δ| = 0.21097`, bound `1.2971`, slack `6.15` (at both `L=8.32` and
`L=50`, both `L`-independent as they must be); `Λ`+shear(0.05) gives slack `34.82` / `320.78`.

**(T′) holds in every row. Worst slack `6.067`** (the `Φ=Λ` rows — the ones the campaign uses).

### 5.1 An instrument that shares nothing with `p3`/`p4` — raw 5-D Monte-Carlo

`p3` and `p4` both use the meridian reduction `J₅ = (r′/r)³ det ∂(r′,z′)/∂(r,z)` and the closed
form (2.1) — the same structure the refuter's `r2` uses. `p7_montecarlo.py` avoids all of it:
`N = 8 000 000` points per `L`, direction uniform on `S⁴` (a normalised 5-D Gaussian), `log ρ`
uniform on `[0,L]`, and `∫_S f dx₅ = L·(8π²/3)·E[f ρ⁵]`. This tests the 5-D measure, the `SO(4)`
reduction, `det DΛ` and the kernel in one shot.

| `L` | quantity | Monte-Carlo | reference | `z` |
|---|---|---|---|---|
| 8.3178 | `a[η₀∘Λ⁻¹](0)` | `5.326003 ± 0.010264` | `5.307374` (Gauss–Legendre) | `+1.81` |
| 8.3178 | `a_ref` | `5.110606 ± 0.009357` | `5.093571` (exact `(M/2)L√(3/2)`) | `+1.82` |
| 8.3178 | Step-3 integral | `316.3973 ± 0.2833` | `315.8653` (exact `π³MA`) | `+1.88` |
| 50 | `a[η₀∘Λ⁻¹](0)` | `30.829647 ± 0.055639` | `30.832425` | `−0.05` |
| 50 | `a_ref` | `30.615847 ± 0.054760` | `30.618622` | `−0.05` |
| 50 | Step-3 integral | `1898.4750 ± 1.6668` | `1898.7389` | `−0.16` |

`max|z| = 1.878` over six comparisons. (The three `L = 8.3178` deviations are the *same* draw, so
they are one `+1.9σ` fluctuation, not three.) No disagreement.

---

## 6. Controls — a gate that cannot fail certifies nothing (FL-043)

Registered before running, in `p5_controls.py`:

**K1 — the `μ_J` term is load-bearing. FIRED.** Delete `3μ_J/8` from (T′). On `Φ = Λ`, `μ = 0`
exactly while `|Δ| = 0.213803 M`, so the surviving bound is `0.000000` and is violated at every
`L ∈ {8.318, 20, 50, 100, 400}`. No bound of the form `C(μ)·μ·M A` can hold. **(H3) is not
optional in the shell-dependent lemma either.**

**K2 — the single-`T_λ̄` reference really is vacuous. FIRED** (own instrument; re-derives the
refuter's central claim without importing their numbers). Optimising `λ̄` over `[1,1.5]` to
minimise the sup relative displacement:

| `L` | `λ̄*` | `μ` | `μ_J` | single-`λ` Lemma-T bound / `a_true` |
|---|---|---|---|---|
| 8.318 | 1.1760 | 0.385344 | 0.8028 | 70.514 |
| 20 | 1.1760 | 0.385344 | 0.7001 | 67.972 |
| 50 | 1.1760 | 0.385344 | 0.6562 | 66.848 |
| 100 | 1.1760 | 0.385344 | 0.6416 | 66.469 |
| 400 | 1.1760 | 0.385344 | 0.6306 | 66.182 |

`μ` is **identical across `L` to `0.0` spread** (my grid gives `0.385344`, the refuter's finer
`λ̄`-grid gives `0.385329`; the `1.50e-5` difference is the `λ̄`-grid spacing, `0.002` vs `0.001`).
So against any single `T_λ̄` the estimate certifies a `66–71×` relative error at every `L` —
vacuous, exactly as the refuter found. **This is why Lemma T′ is needed at all.**

**K3 — `μ_J` cannot be dropped even at small `μ`. FIRED.** `Λ` composed with the axis-ramp shear
`β(φ) = μ_b m(φ/φ_c)m((π−φ)/φ_c)cos φ`, `m(t)=t²(3−2t)`, whose relative displacement is
`2|sin(β/2)|` and is `φ_c`-independent while `μ_J ~ φ_c^{-3}`:

| `μ_b` | `φ_c` | `μ` | `μ_J` | `\|Δ\|` | `μ`-only bound | violated? |
|---|---|---|---|---|---|---|
| 0.3 | 1e-1 | 0.29740 | 4.34e+02 | 3.3909 | 108.18 | no |
| 0.3 | 1e-2 | 0.29886 | 1.81e+06 | 8.5746 | 109.98 | no |
| 0.3 | 1e-3 | 0.29888 | 1.64e+10 | 50.8019 | 110.00 | no |
| 0.3 | 1e-4 | 0.29888 | 1.62e+14 | 464.3991 | 110.00 | **yes** |
| 0.3 | 1e-5 | 0.29888 | 1.62e+18 | 4591.7832 | 110.00 | **yes** |
| 0.9 | 1e-5 | 0.86993 | 1.11e+20 | 236993.66 | 2097664.55 | no |

`μ` spread over the `μ_b = 0.3` rows: `1.5e-7`. Fitted growth on the asymptotic subrange
(`φ_c ≤ 1e-3`): `|Δ| ≈ 0.05835 φ_c^{-0.9781}` (full range `φ_c^{-0.7997}`; `gap-T-lipschitz` §5
measured `φ_c^{-0.9855}` for the same `μ_b` over ten decades, and this note reproduces that rate
from a different starting map). The **full** bound (T′), evaluated at the measured `μ_J`, holds in
every row (by 4–20 orders) — which is the point: the `μ_J` in (T′) is load-bearing, not decorative.

**Mutation coverage of `check_constants.py` (`p6_mutation.py`).** Every stored number in
`p1…p5_results.json` was mutated, one at a time, by `v ↦ 1.5v + 0.37` (booleans flipped), and the
gate re-run. Result: **678 stored numbers mutated, 678 caught, 0 blind — coverage 1.0000.**

`check_constants.py` executes **553** checks (count is `len(CHECKS)`, printed, never a
literal);
no check multiplies its input by zero; every check consumes a stored number and compares it against
a value rebuilt from other stored numbers or from mathematics.

**Audit of this document against the JSONs (`p8_doc_audit.py`).** Every decimal number printed in
PROOF.md outside code fences was matched, to the precision it is printed at, against a stored
number or an explicitly listed function of stored numbers: **295 distinct numbers, 0 untraced**.
A further **17** numbers are quoted verbatim from other seats (`prove-lagrangian`,
`gap-T-lipschitz`, `refute-gap-T-lipschitz`, the synthesis) and are listed in
`p8_results.json` as such — they are cited, not computed here.

---

## 7. Status per step

| step | statement | status |
|---|---|---|
| S0 | `Λ` is a `C¹` diffeomorphism onto its image when `2κ_s/L < 1` (positive Jacobian + contraction ⟹ injective) | **PROVED** |
| S1 | change of variables; the two-term split (2.3) | **PROVED** |
| S2 | `\|𝒦\| ≤ (3/8π²)\|w\|^{-4}`, `\|∇𝒦\| ≤ (3/2π²)\|w\|^{-5}`, both sharp | **PROVED** (sympy, residual 0) |
| S3 | `∫₀^π sin²φ/g(φ;λ)⁴dφ = π/(2λ)`; hence `∫_S\|η₀\|λ²\|Λx\|^{-4}dx₅ ≤ π³MA` (equality iff `\|ω₀^θ\|=M`) | **PROVED** (exact; quadrature `1.3e-15`) |
| S4 | segment misses the origin (needs `μ<1`); mean-value bound | **PROVED** |
| S5 | Jacobian term | **PROVED** |
| **T′** | the estimate (T′), explicit constants | **PROVED** |
| C1 | relative form for the cap, `11.78097245 μ` | **PROVED** |
| C2 | tapered datum, factor `1/r_h`, `r_h(7.5°) = 0.981822` | **PROVED** (`P_h` numeric) |
| C3 | composition `(H2_Λ),(H3_Λ) ⟹ (H3)` with `μ̃_J = μ_J^Λ(1+2κ_s/L) + 2κ_s/L` | **PROVED** |
| C4 | `Φ = Λ`: `μ = 0`, `μ̃_J = 2κ_s/L`, `\|Δ\|/a_ref ≤ (3π/2)κ_s/L = 2.1181705/L` | **PROVED** |
| §3 | `κ_s = √6−2`, `μ_J L = 2√6−4`, `λ̄ = √6/2`, `C_rel = 3π(√6−2)/2` | **PROVED** (sympy exact) |
| §4 | the true shear cost `= 0.2138034 M`, i.e. `0.3491394/L` relative | **MEASURED** (two instruments, `5e-15`) — not proved |
| K1,K2,K3 | necessity of `μ_J`; vacuity of the single-`λ` reference; `\|Δ\|` unbounded at fixed `μ` | **ESTABLISHED BY COUNTEREXAMPLE** (numerical) |
| — | that the NS flow map satisfies (H2),(H3) with `μ, μ_J^Λ = O(c/L)` | **NOT ADDRESSED** (`L3v` + `Γ`) |

---

## 8. What this closes, and what it does not

**Closes.** The refuter's item 1 of §8 — "the restatement/decomposition above — cheap, verified
numerically here, unwritten". It is now written, with the constant `2.1181705106873974` **derived**
from `κ_s = √6−2` rather than read off a table, and with the hypothesis list, the diffeomorphism
condition `L > 2κ_s`, the composition rule (2.11), and the tapered-datum factor `1/r_h` all stated.
Lemma T (single map) remains true and is now seen to be the special case `λ ≡ const` of (T′);
the two have identical constants.

**Does not close — and this is the whole of GAP T that remains.**
1. **The PDE half.** That the Euler/NS flow map of this datum on `[0,τ]`, `τ = c/(ML)`, satisfies
   (H2) and (H3_Λ) with `μ, μ_J^Λ = O(c/L)` **uniformly over the whole shell including both
   boundary layers**, where `prove-lagrangian` §4(2) itself records `sup|a − A(ρ)| = 0.7047 M` at
   the outer edge against `0.3187 M` one octave in. Nothing in this note bears on it. It needs
   `L3v` (still SKETCH: the uniform-in-`λ` mode-sum bound `Σ_{l≥3}|H_l(λ)|‖C_l‖_∞/(l(l+3)−4)`,
   whose displayed value `0.3958` the synthesis shows is truncation-dependent with true sum
   `≈ 0.4295`) **and** `Γ = ‖∇u‖_∞`, which the `gap-V` refuter records as an unproven hypothesis
   nowhere stated as one.
2. **`Λ` is not the flow map of an incompressible field.** `J_Λ = λ²(1 + D W) ≠ λ²`, so the
   shell-dependent family is not volume-preserving in `ℝ³` (the single `T_λ` is). The true flow is.
   (T′) prices that mismatch — it is exactly the `2κ_s/L` in (2.11) — but it does **not** show a
   volume-preserving map exists within `O(c/L)` of `Λ`.
3. Nothing here touches conjecture (ii), BFG Thm 10, or FL-000. **FL-000 stands.**

---

## 9. Corrections this note carries

1. **Kernel sign (task brief).** The origin kernel is `𝒦(x) = −3x_z/(8π²|x|⁵)`, i.e. `K(−x)`; the
   brief's `+3w_z/(8π²|w|⁵)` is `K` itself and reverses the sign of `a(0)`. Checked in `p1` §4.
2. **`κ_s ≈ 0.899c` (task brief).** `0.8989794856` is `μ_J·L = 2κ_s`, a pure number; `κ_s = √6−2 =
   0.4494897428`; and against `c = 2log(3/2)` the ratio is `1.1085781089`, not `0.899`.
3. **The refuter's `μ_J(local) = 0.8990/L` and band `2.03/L … 2.12/L` are exactly reproduced** and
   are now identified in closed form: `μ_J L = 2√6−4`, and the band is
   `C_rel/L ÷ (1 + 0.3491394/L)` with `C_rel = 3π(√6−2)/2 = 2.1181705106873974`. Their `2.12/L` at
   `L = 400` is the tail of that convergence, not a fitted constant.
4. **A `P_{h_δ}(3/2)` cross-check does not close** at `δ ≥ 10°` (up to `3.4 %` at `30°`) against the
   synthesis's tabulated `Φ_{h_δ}(3/2)`; the `λ = 1` values agree to `1e-14` across four seats.
   Recorded as unresolved; not load-bearing here.
5. **The refuter's `μ = 0.385329` against a single `T_λ̄` is confirmed independently** at
   `0.385344` (my coarser `λ̄`-grid), `L`-independent to `0.0` spread over `L = 8.3 … 400`.

---

## 10. Files

| file | what it establishes |
|---|---|
| `p1_kernel.py` / `p1_results.json` | `𝒦`, `∇𝒦`, the sharp constants `3/(8π²)`, `3/(2π²)`; `det DΛ` by two independent symbolic routes; the ranges `[−2,1]` and `[−1,2]`; the kernel-sign check |
| `p2_identities.py` / `p2_results.json` / `p2_log.txt` | `∫₀^π sin²/g⁴ = π/(2λ)` exactly (+40-dps quadrature, `0.0`); `∫₀^{π/2}cos sin²/g⁵ = λ/3`; the closed form of `Q(λ)` and `Q(1) = −1/15`; `P_{h_δ}` and `r_h` |
| `p3_constants.py` / `p3_results.json` | `κθ`, `λ̄ = √6/2`, `∫dσ/λ`, `κ_s = √6−2`, `μ_J L = 2√6−4`, `C_rel = 3π(√6−2)/2`, the shear offset, and the row-by-row comparison against the refuter's stored `r2_results.json` |
| `lib5.py` | own 2-D reduction of the 5-D functional, the map stages and their exact Jacobians, the (T′) bound |
| `p4_verify.py` / `p4_results.json` | instrument cross-checks; (T′) verified on `Φ = Λ` and on 12 non-trivial `Φ` and 4 tapered rows; convergence |
| `p5_controls.py` / `p5_results.json` | controls K1 (`μ_J` necessary), K2 (single-`λ` vacuous), K3 (`|Δ|` unbounded at fixed `μ`) |
| `p6_mutation.py` / `p6_results.json` / `p6_log.txt` | measured mutation coverage of the gate |
| `p7_montecarlo.py` / `p7_results.json` / `p7_log.txt` | raw 5-D Monte-Carlo, no meridian reduction: an instrument sharing nothing with `p3`/`p4` |
| `p8_doc_audit.py` / `p8_results.json` | audit of every decimal number in PROOF.md against the stored JSONs |
| `check_constants.py` | 553 checks re-asserting every number above against the JSONs |
| `SHA256SUMS` | computed, never typed |

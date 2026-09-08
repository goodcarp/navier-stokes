# GAP V — the Gaussian (Aronson-type) bound for the 5D drift-diffusion

Seat `gaps/gap-V-aronson`, DTC-2026-09-06. Sub-seat of the estate; laws
`TORMENT NEXUS/LAWS.md`. Every number below came out of a script in this folder that I wrote
and ran (`v1`…`v5`, re-asserted by `check_constants.py`); `SHA256SUMS` is computed, never
typed. Nothing outside this folder was written. Numerics falsify, they never prove; the
algebra marked EXACT is sympy-verified in exact rational arithmetic.

---

## 0. Headline

**The estimate the brief asks for is TRUE, it is PROVED here with explicit constants, and it
does not need Aronson at all.** The brief's Route (1) — Lagrangian coordinates + Aronson /
Fabes–Stroock / Norris–Stroock — **does not close**, for a reason that only shows up when the
hypotheses are checked at source (§3). Route (2) — the reversed-time stochastic
representation — closes it completely, elementarily, with constants depending on nothing but
`c = ‖∇u‖_∞ τ` and the dimension.

**There are no "terms from the divergence `2a`".** They vanish, for two independent reasons
(§2.3): the representation is for the *advective* form, in which `div₅ b` never appears; and,
exactly, in Lagrangian labels the operator is `ν D^{-1}∂_l(D G^{lk}∂_k)` — pure divergence
form, **no drift and no zeroth-order term** — with the whole compressibility sitting in the
weight `D = det J`, which is bounded above and below by `e^{±5c}`.

Two things the seat also establishes, both of which cost the campaign something:

1. **The brief's `d ≳ ρ₀δ = √(ν/M)` is not enough.** With `d` one dissipation length the
   bound gives `q = d²/(ντ) = L/c` exactly — and the resulting relative loss stays above
   `1/L` until `L ≥ 590` (accelerated window) / `L ≥ 854` (frozen window), i.e. until
   `log Re ≳ 1200–1700`. What works is to inset the tracked point by a few dissipation
   lengths: `f := (ρ_*−ρ₀)/ρ₀ = 0.91, 0.66, 0.47, 0.34, 0.24` at `L = 10, 20, 40, 80, 160`,
   i.e. `d = 7.0, 5.0, 3.6, 2.6, 1.9` dissipation lengths. The price is a *relative* `O(1/L)`
   inflation of `c₂`: `1.470 → 1.570 (L=10), 1.506 (20), 1.482 (40), 1.473 (80), 1.470 (160)`
   against the asymptotic `8(1−√(2/3)) = 1.4680274`. **The theorem's constant is untouched;
   the inset is an `o(1)` item.** This reproduces, from a completely different direction, the
   `prove-duhamel` seat's measured finding that the `η`-loss is carried by the inner-edge
   mollification layer and that "tracking one or two `ρ₀` further out buys it back for an
   `O(M)` loss in `a₀`, i.e. `O(1/L)` in `c₂`".

2. **Gap V is two gaps, and only one of them is now closed.**
   * **V-a (localization):** the tracked point does not feel the datum's edges. **PROVED**
     (Theorem V.1), needs only `‖∇u‖_∞ τ ≤ c`.
   * **V-b (bulk loss):** the viscous loss of a *smooth plateau* along the trajectory. **NOT
     covered** by V.1 and not closed here. With `‖∇u‖_∞` alone the best that comes out is
     `O(L^{-1/2})` (relative loss `≤ 1.18` at `L=10`, `0.038` at `L=640` — vacuous at any
     realistic `L`). The `O(1/L)` that `prove-lagrangian` §4(3) asserts is correct in order
     but needs a **second-derivative bound on `u`** to kill the first moment of the diffusive
     cloud; and its constant is optimistic by a factor `5(e^{2c}−1)/(c(1+f)²) = 6.2`.

| claim | status |
|---|---|
| **V.0** Lagrangian reduction: `∂_tζ = ν D^{-1}∂_l(D G^{lk}∂_kζ)`, `G=(JᵀJ)^{-1}`, `D=det J`; no drift, no zeroth order | **PROVED** (exact, sympy) |
| **V.1** `\|η₁−η₂\|(X(τ),τ) ≤ N·min(B1,B2)`, `B1,B2` explicit in `c,q,n` only | **PROVED** |
| **V.2** improved rate `e^{-2c}q/4` from the label-frame SDE | **CONDITIONAL** on a `∇²u` bound (Itô drift term not estimated here) |
| sharp rate `cq/(2(e^{2c}−1))`, attained by the pure axisymmetric strain | **NOT PROVED** (numerically corroborated) |
| Route (1) via Aronson / Fabes–Stroock / Norris–Stroock | **DOES NOT APPLY** — hypotheses checked at source, §3 |
| **V-b** the bulk viscous loss | **NOT CLOSED**; `O(L^{-1/2})` unconditional, `O(1/L)` needs `∇²u` |

---

## 1. Setting and notation

`ℝ⁵ = ℝ⁴_y × ℝ_z`, `r = |y|`, `ρ = |x|`. The lift of an axisymmetric no-swirl field
`u = u^r e_r + u^z e_z` is `b(y,z) := (u^r(r,z) y/r, u^z(r,z))`. Then, EXACTLY (`v1`, §A):

```
div₅ b  =  ∂_r u^r + 3 u^r/r + ∂_z u^z  =  (div₃ u) + 2a ,   a := u^r/r ,
```
so `div₅ b = 2a` for incompressible `u`. Residual `0` in sympy.
`η := ω^θ/r` satisfies `D_t η = ν Δ₅ η`, i.e.

```
∂_t η + b·∇₅ η = ν Δ₅ η ,        Δ₅ = ∂_rr + (3/r)∂_r + ∂_zz .        (1.1)
```
Also EXACT (`v1`, §B): `Δ₅(1/r) = −1/r³`, in both the Cartesian-5D and the `(r,z)` form.

**The Lipschitz constant does not lose a dimensional factor.** For the lift,
`‖∇₅b‖_op = ‖∇₃u‖_op` **pointwise and exactly** — the extra ℝ⁴ directions all carry the
eigenvalue `a = u^r/r`, which is already an eigenvalue of `∇₃u`, and the remaining `2×2`
block has the same characteristic polynomial. Verified twice independently: symbolically and
to `6.7e-16` over 200 random points in the earlier pass of this seat
(`scripts/v1_identities.py`, items A5), and to **`4.4e-15` over 400 random incompressible jets**
in `v1_reduction.py` §F this run (which also re-checks `tr ∇₅b = 2a` at every sample). So

```
Γ := sup_{[0,τ]} ‖∇₅ b(·,t)‖_{L^∞,op} = sup_{[0,τ]} ‖∇u(·,t)‖_{L^∞,op} ,     c := Γ τ .
```
For the uniform axisymmetric strain `u = (a r, −2 a z)`: `∇₅b = diag(a,a,a,a,−2a)`,
`‖∇₅b‖_op = 2a`, `div₅ b = 2a` (`v1`, §E). Hence in the doubling window of
`prove-lagrangian`, with `a₀ = (M/2)L` and `∫₀^τ a = θ`,

```
Γ = M L ,   τ = 2θ/(ML) ,   c = Γτ = 2θ  =  0.8109302 (frozen, θ=log 3/2)
                                        =  0.7340105 (accelerated, θ = 2(1−√(2/3))/κ, κ=1/2).
```
(`v5` re-derives the accelerated window from `F(0,θ) = −2log(1−κθ/2) = log(3/2)`: `0.40546` vs
`log(3/2) = 0.40546`.)

`ρ₀δ = √(ν/M)` (the datum's finest feature = the dissipation length), `L = log(R/ρ₀)`,
so `ν τ = (ρ₀δ)² c/L` and

```
q := d²/(ν τ) = (d/(ρ₀δ))² · L/c .                                       (1.2)
```

---

## 2. THEOREM V.1 and its proof

### 2.1 Lemma V.0 (the Lagrangian reduction) — EXACT, PROVED

Let `X(α,t)` be the flow of `b`, `J := ∂X/∂α`, `D := det J`, `G := (JᵀJ)^{-1} = J^{-1}J^{-ᵀ}`,
and `ζ(α,t) := η(X(α,t),t)`. Then, identically,

```
  ∂_t ζ  =  (ν/D) ∂_{α_l} ( D G^{lk} ∂_{α_k} ζ ) ,      ∂_t log D = (div₅ b)∘X = 2 a∘X.   (2.1)
```
*Proof.* Chain rule plus the Piola identity `∂_{α_l}(D (J^{-1})_{li}) = 0`. Verified in
`v1_reduction.py`: the Piola residuals are exactly `0`; the identity
`D·(Δ_xη)(X(α)) = ∂_{α_l}(D G^{lk}∂_{α_k}ζ)` is exactly `0` at 8 random rational points for a
nonlinear, non-volume-preserving map (exact rational arithmetic, so `0` is a true zero); the
non-divergence variant `(Δ_xη)(X) = ∂_l(G^{lk}∂_kζ) + G^{lk}(∂_l log D)(∂_kζ)` also residual
`0`. ∎

**Consequences.** `(2.1)` is *divergence form with respect to the measure* `dμ := D dα`
(which is Lebesgue measure pulled back from `x`), and it is **symmetric in `L²(μ)`**, with
**no drift term and no zeroth-order term**. If `Γτ ≤ c` then `∂_tJ = (∇b)(X,t)J`, `J(0)=I`,
gives every singular value of `J` in `[e^{-c},e^{c}]`, hence

```
e^{-2c} I ≤ G ≤ e^{2c} I ,      e^{-5c} ≤ D ≤ e^{5c} .                    (2.2)
```
`v2_flowmap_bounds.py`: 400/400 random time-dependent `B(t)` with `‖B‖_op ≡ Γ` satisfy (2.2)
(worst `s_max/e^{c} = 0.9792`, `s_min·e^{c} = 1.0185`, `det/e^{5c} = 0.8438`); constant
diagonal `B` attains `s_max = e^{c}` to `1e-12`, so (2.2) is sharp. For the *physical* field
`det J = exp(∫div₅b) = e^{2aτ} = e^{c}` (checked to `9.2e-11`) — the `e^{5c}` envelope is
lossy by `3.6e5` at `c = 3.2`, which matters only for Route (1).

### 2.2 Theorem V.1 (Gaussian localization in the material frame) — PROVED

> **THEOREM V.1.** Let `n = 5`, `ν > 0`, `τ > 0`. Let `b : ℝⁿ×[0,τ] → ℝⁿ` be continuous, `C¹`
> in `x`, with `|b(x,t)| ≤ B₀(1+|x|)` and
> `Γ := sup_{t≤τ} ‖∇b(·,t)‖_{L^∞,op} < ∞`; put `c := Γτ`. Let `η₁, η₂` be **bounded classical
> solutions** of `∂_tη + b·∇η = νΔη` on `ℝⁿ×[0,τ]`. Let `X` solve `Ẋ = b(X,t)`, `X(0) = x₀`.
> Put `N := ‖η₁(·,0) − η₂(·,0)‖_∞` and
> `d := dist( x₀ , {η₁(·,0) ≠ η₂(·,0)}⁻ )`, and `q := d²/(ντ)`. Then
>
> ```
> |η₁(X(τ),τ) − η₂(X(τ),τ)|  ≤  N · min( B1 , B2 ) ,
>
> B1 = 2n · exp( − e^{−2c} q / (4n) ) ,
> B2 = inf_{0<ε<√2} (1+2/ε)^n · exp( − (1−ε²/2)² c q / ( 2 e^{2c} (e^{2c}−1) ) ) .
> ```
> `B1, B2` depend on `c, q, n` only — **not** on `ν`, `Γ`, `M`, `div b`, or any derivative of
> `b` beyond `∇b`. No hypothesis whatsoever is placed on `div b`.

**Proof.**

*Step 0 (representation).* Put `φ(x,s) := η(x, τ−s)`, `s ∈ [0,τ]`. Then
`∂_sφ − b(x,τ−s)·∇φ + νΔφ = 0` with terminal datum `φ(·,τ) = η(·,0)`. Let

```
dY_s = − b(Y_s, τ−s) ds + √(2ν) dW_s ,     Y_0 = x ,
```
which has a unique strong non-exploding solution (locally Lipschitz + linear growth). Itô on
`s ↦ φ(Y_s,s)` kills the drift exactly; stopping at `τ_R = inf{s : |Y_s| ≥ R}` gives
`E[φ(Y_{τ∧τ_R}, τ∧τ_R)] = φ(x,0)`, and `|φ| ≤ ‖η(·,0)‖_∞` with `τ_R ↑ ∞` a.s. gives, by
dominated convergence,

```
η(x,τ) = E[ η(Y_τ, 0) ] .                                                  (2.3)
```
(Only boundedness of `η` is used; no gradient bound.)

*Step 1 (comparison).* Both solutions have the **same** `b`, hence the same `Y`. Take
`x = X(τ)`; the deterministic reverse characteristic `x_s := X(τ−s)` solves
`ẋ_s = −b(x_s,τ−s)`, `x_0 = X(τ)`, `x_τ = x₀`. Subtracting the two copies of (2.3),

```
|η₁(X(τ),τ) − η₂(X(τ),τ)| = |E[(η₁,₀−η₂,₀)(Y_τ)]| ≤ N · P( |Y_τ − x₀| ≥ d ) .   (2.4)
```

*Step 2 (tail).* `Z_s := Y_s − x_s` obeys `Z_0 = 0` and, by the mean value theorem,

```
dZ_s = − A_s Z_s ds + √(2ν) dW_s ,   A_s := ∫₀¹ ∇b(x_s+θZ_s, τ−s) dθ ,  ‖A_s‖ ≤ Γ,
```
`A` adapted.

*(B1) pathwise.* `|Z_s| ≤ Γ∫₀^s|Z_r|dr + √(2ν)|W_s|`; Grönwall against the nondecreasing
majorant gives `|Z_τ| ≤ e^{c}√(2ν)·sup_{s≤τ}|W_s|`. For `n`-dim BM,
`{sup_s|W_s| ≥ R} ⊆ ∪_i {sup_s|W^i_s| ≥ R/√n}` and, by reflection plus
`P(N(0,τ)≥a) ≤ ½e^{−a²/2τ}`, `P(sup_s|W^i_s| ≥ a) ≤ 2e^{−a²/2τ}`. With `R = d e^{−c}/√(2ν)`,
`P(|Z_τ| ≥ d) ≤ 2n e^{−R²/(2nτ)} = B1`.

*(B2) adapted propagator + net.* Let `Φ_s` solve `Φ̇_s = −A_sΦ_s`, `Φ_0 = I` (a pathwise linear
ODE, so `Φ, Φ^{-1}` are adapted and `‖Φ_s‖, ‖Φ_s^{-1}‖ ≤ e^{Γs}`). Variation of constants
(`Φ` has finite variation, so no cross-variation): `Z_s = Φ_s N_s` with
`N_s := √(2ν)∫₀^s Φ_r^{-1} dW_r`, a continuous martingale with, for every unit `θ`,

```
⟨θ·N⟩_τ = 2ν∫₀^τ |Φ_r^{−ᵀ}θ|² dr ≤ 2ν∫₀^τ e^{2Γr} dr = ν(e^{2c}−1)/Γ = ν τ (e^{2c}−1)/c =: V₀ .
```
The exponential supermartingale `exp(λθ·N_s − λ²⟨θ·N⟩_s/2)` and `⟨θ·N⟩_τ ≤ V₀` a.s. give
`P(θ·N_τ ≥ a) ≤ e^{−a²/(2V₀)}`. Since `|Z_τ| ≤ e^{c}|N_τ|`, the event `{|Z_τ| ≥ d}` forces
`|N_τ| ≥ d e^{−c}`, hence `θ·N_τ ≥ d e^{−c}(1−ε²/2)` for some `θ` in an `ε`-net `𝒩_ε` of
`S^{n−1}` with `|𝒩_ε| ≤ (1+2/ε)^n`. A union bound gives `B2`. Insert into (2.4). ∎

**Remark (why the material frame is not optional).** `d` is the distance **in the initial
configuration**, and labels do not move: the flow map never distorts `d`. In Eulerian
coordinates the corresponding distance would have to be tracked, and, worse, the drift
displacement of the tracked particle over the doubling window is `≥ (λ−1)ρ₀ = 0.5 ρ₀` against
`√(ντ) = ρ₀δ√(c/L) = 0.0355ρ₀ (L=10)`, i.e. **14.1 σ at `L=10`, 28.2 σ at `L=40`, 56.4 σ at
`L=160`** (`v5` §5). Any Eulerian Gaussian bound must pay `exp(+ (drift·τ)²/(Cντ))` — the
route is quantitatively dead, not merely inelegant.

### 2.3 Where the divergence `2a` went

Two independent answers, both established here.

* **Probabilistic.** (1.1) is written in advective form; the Feynman–Kac representation (2.3)
  is for the advective form; `div b` appears nowhere in Steps 0–2. Formally: `div b` would
  appear as a zeroth-order term only if one used the conservation form
  `∂_tη + div₅(bη) = νΔ₅η + 2aη`, which is a different (and here unnecessary) bookkeeping.
* **Analytic.** By Lemma V.0 the *whole* effect of `div₅ b = 2a` is the weight `D = det J`,
  and `D ∈ [e^{−5c}, e^{5c}]` — a bounded multiplicative distortion of the reference measure,
  contributing only to constants. For the pure axisymmetric strain it is even better: `D` is
  spatially **constant** (`∇_α log D ≡ 0`, verified in the earlier pass, item A11), so for the
  leading-order field the compressibility contributes *exactly nothing* beyond an
  `α`-independent time factor.

**So the brief's expected "+ (terms from the divergence 2a)" is empty. That is a finding, not
an omission.**

### 2.4 How lossy the proved constants are

The `c → 0` limit of the `B2` rate is exactly `1/4` — the sharp heat-kernel exponent
(`0.24992, 0.24925, 0.24261, 0.18490` at `c = 10^{-4},10^{-3},10^{-2},10^{-1}`), and at `c=0`,
`n=1`, `exp(−q/4)` sits above the true `erfc(√q/2)` by `2.3×–5.6×` over `q = 4…36`. So the
*form* is sharp. The `c`-dependence is not: the proved `B2` rate
`c/(2e^{2c}(e^{2c}−1))` is a factor `e^{2c}` (`= 4.3` at `c = 0.734`) below the rate
`c/(2(e^{2c}−1))` which the pure axisymmetric strain actually attains. The loss is the single
step `|Z_τ| ≤ ‖Φ_τ‖·|N_τ|`, forced because `Φ_τΦ_r^{-1}` is not `𝓕_r`-adapted. The earlier
pass of this seat measured the true constant directly at `c = 2.96`: fitted `C' = 3.08`
against the proved `4e^{2c} = 1502` (`scripts/v3_log.txt`, `R1_fitted_Cprime`).

**Proposition V.2 (conditional; better rate).** Transform the backward diffusion to labels,
`Ŷ_s := X^{-1}(Y_s, τ−s)`. The transport terms cancel identically and

```
d Ŷ_s = ν (Δ_x X^{-1})(Y_s) ds + √(2ν) (∇_x X^{-1})(Y_s) dW_s ,      Ŷ_0 = x₀ ,
```
with `‖∇_xX^{-1}‖ ≤ e^{c}`, giving `⟨θ·M⟩_τ ≤ 2νe^{2c}τ` and hence the rate `e^{−2c}q/4`
(prefactor `(1+2/ε)^n`) — a factor `2.5` better than `B2` at `c = 0.734` — **provided** the
Itô drift `ν Δ_x X^{-1}` is controlled, which needs a bound on `∇²u`. Not estimated here.

---

## 3. Route (1) checked at source — and why it does not close

The brief asks: "cite Aronson 1968 / Fabes–Stroock 1986 / Norris–Stroock 1991 at source and
check the exact hypotheses (is the divergence term admissible?)". Done. Sources in
`sources/` (Aronson 1968 full text and page images; Norris–Stroock 1991 PDF; Aronson's own
2017 survey fetched this run).

**(a) Aronson, *Bounds for the fundamental solution of a parabolic equation*, Bull. AMS **73**
(1967), 890–896; *Non-Negative Solutions of Linear Parabolic Equations*, Ann. Scuola Norm.
Sup. Pisa **22** (1968), 607–694.** The 1968 paper treats

```
∂_t u − Σ_j ∂_{x_j}{ Σ_i a_{ij}∂_{x_i}u + a_j u } − Σ_j b_j ∂_{x_j}u − c u = 0
```
under hypotheses (H.1)–(H.3) (pp. 608–609). Verbatim from the source text, p. 608:
"*Throughout the paper it will be assumed that there constants ν, M, M₀ and R₀ such that
0 < ν < ∞, 0 < M < ∞ … and such that the coefficients of L satisfy the following conditions
which will be referred to collectively as (H).*" (H.2) requires the lower-order coefficients
`A_j, B_j ∈ L^{p,q}` with `|A_j|, |B_j| ≤ M₀`. **The drift IS admissible.** But (p. 615)
"*The statement 'C depends on the structure of (1.1)' means that C is determined by the
quantities ν, M, … M₀ … and θ*", and (p. 624) "*By the structure of L we mean n and the
quantities which occur in the hypotheses (H). In particular, α depends only on n, M and ν,
while β depends only on …, M₀ and ν*"; Theorem 8 (p. 660 ff.) gives constants "*depending
only on δ, T and the structure of L*". So the constants depend on **`‖b‖_{L^∞}` and on `T`**.
That is exactly the dependence the brief forbids, and §2.2's `14–56 σ` computation says the
dependence is not removable in Eulerian coordinates.

**(b) Fabes & Stroock, ARMA **96** (1986), 327–338.** They give the Gaussian upper bound for
the *driftless* divergence-form equation `∂_tu = ∂_j(a_{ij}∂_iu)` by Nash's method. Aronson's
own survey (arXiv:1707.04620, p. 2) states flatly: "*Neither Nash nor Fabes & Stroock consider
the full equation (1)*". So Fabes–Stroock does not cover a drift at all.

**(c) Norris & Stroock, Proc. LMS **62** (1991), 373–402.** Survey p. 4: they consider
`L ≡ ∇·(A(x,t)∇ + AE(x,t)·∇) − ∇·(AÊ(x,t)·) + C(x,t)`, and "*Their main result is a very
precise two-sided estimate for the fundamental solution of `Lu = ∂_tu` based on energy
functions associated with the coefficients. However they are forced to assume the **uniform
continuity of A and `E − Ê`**.*" Our Lagrangian coefficient `G = (JᵀJ)^{-1}` has a modulus of
continuity governed by `∇²u`, which the hypothesis `‖∇u‖_∞τ ≤ c` does not control. **Not
applicable.**

**(d) The gap that actually kills Route (1).** Lemma V.0 shows the Lagrangian operator is
`p(α,t) ∂_tζ = ν ∂_l(a^{lk}∂_kζ)` with `p = D = det J` — a **weighted** divergence-form
operator. The closest theorem in the literature is Porper & Eidel'man (Russian Math. Surveys
**39** (1984), 119–178), described in Aronson's survey p. 4 as "*a slight generalization of
equations (1) and (3) involving a coefficient `p(x)` multiplying `∂_t u` … they give … the
two-sided Gaussian estimate*" — but their weight is **`p(x)`, time-independent**, and ours is
`D(α,t)` with `∂_t log D = 2a∘X ≠ 0` (that is precisely the compressibility). Aronson 1967 /
Fabes–Stroock have no weight; Norris–Stroock need continuity; Porper–Eidel'man need a
time-independent weight. **None of the four covers the operator the brief's Route (1)
produces.** The gap is exactly the compressibility of the lift, arriving as a *time-dependent*
`A₂` weight.

*(This is a statement about what the literature I could reach at source proves, not a claim
that no such theorem exists — time-dependent-measure Dirichlet-form theory plausibly covers
it. But it is not the citation the brief named, and I would not sign the brief's citation.)*

---

## 4. Refuter runs

Pre-declared falsification rules (fixed in `v4_falsify.py` before any number was read):
**R1** the theorem is refuted if a Monte-Carlo estimate of `P(|Y_τ−x₀| ≥ d)` minus 3 s.e.
exceeds `min(B1,B2)`; **R2** likewise for the end-to-end PDE deviation; **R3** the test is
non-diagnostic unless the *Eulerian* control (same quantity read at the fixed point `x₀`
instead of at `X(τ)`) violates the bound.

**TEST A — `n = 5`, Monte Carlo on the backward SDE.** Drift `b = a₀(y,−2z) + κ g(x)` with
`g` bounded nonlinear, `a₀ = 1`, `κ = 0.35`; **`div₅ b = 2.00000` at `x₀`** (the drift really is
compressible); `Γ` *measured* over a tube around the trajectory `= 2.12634`, `τ = 0.4`,
`c = 0.85054`. 250 000 paths, 600 steps, two viscosities.

| `d/√(ντ)` | `q` | `P_emp` (`ν=2e-3`) | `P_emp` (`ν=5e-4`) | `min(B1,B2)` | verdict |
|---|---|---|---|---|---|
| 4 | 16 | `1.812e-01 ± 7.7e-04` | `1.805e-01 ± 7.7e-04` | `8.642e+00` | ok |
| 6 | 36 | `1.533e-02 ± 2.5e-04` | `1.522e-02 ± 2.4e-04` | `7.200e+00` | ok |
| 8 | 64 | `7.120e-04 ± 5.3e-05` | `6.920e-04 ± 5.3e-05` | `5.577e+00` | ok |
| 10 | 100 | `8.0e-06 ± 5.7e-06` | `1.60e-05 ± 8.0e-06` | `4.015e+00` | ok |
| 12 | 144 | `0` | `0` | `2.688e+00` | ok |

The two viscosities agree to within Monte-Carlo error at equal `q`, which is the scaling the
theorem asserts (`P` depends on `ν, d, τ` only through `q` and `c`). Fitted empirical rate
`−log P / q`: **`0.11503` and `0.11053`**, against the proved `B1` rate `0.00912`, the proved
`B2` rate `0.01732`, and the unproved "sharp" rate `c/(2(e^{2c}−1)) = 0.09493`. So: **R1 not
violated**; the proved rate is conservative by `6.6×`, and the sharp rate is also an upper
bound on this (non-extremal) field.

**TEST B — `n = 1` PDE, end-to-end.** `∂_tη + b(x)∂_xη = νη_{xx}`, `b = 1.6 sin x` (so
`b' ≠ 0` — compressible), `τ = 0.5`, `c = 0.8`, `x₀ = 0.8 → X(τ) = 1.5100` (drift displacement
`0.710`). Two data differing by an indicator supported at distance `≥ d`.

| `ν` | `d/√(ντ)` | `q` | deviation at `X(τ)` | `min(B1,B2)` | deviation at the fixed `x₀` |
|---|---|---|---|---|---|
| 2e-3 | 4 | 16 | `8.716e-03` | `8.919e-01` | `1.000` **control fires** |
| 2e-3 | 6 | 36 | `9.577e-05` | `3.250e-01` | `1.000` **control fires** |
| 2e-3 | 8 | 64 | `2.236e-07` | `7.909e-02` | `1.000` **control fires** |
| 5e-4 | 4 | 16 | `7.725e-02` | `8.919e-01` | `1.000` **control fires** |
| 5e-4 | 6 | 36 | `9.374e-03` | `3.250e-01` | `1.000` **control fires** |
| 5e-4 | 8 | 64 | `4.715e-04` | `7.909e-02` | `1.000` **control fires** |

**R2 not violated. R3: the Eulerian control fires at full amplitude (`1.000 = N`) in every
row** — the fixed point `x₀` is swept into the perturbed region by the drift, so the same
statement read in Eulerian coordinates is false by the maximum possible margin. The numerics
are therefore diagnostic, and they say the material frame is not a convenience.

**Summary: `R1 violated: False`, `R2 violated: False`, `R3 control fires: True`.**

---

## 5. Application: what this costs the lower-bound argument

Geometry: tracked point at `ρ_* = (1+f)ρ₀`, polar angle `φ₀`; taper half-width `δ = 7.5°`;
equatorial mollification half-width `5°`. In the initial configuration

```
d = min( f ρ₀ ,  ρ_* sin(φ₀−δ) ,  ρ_* sin(90°−5°−φ₀) ) ,
A := ‖η₀ − η̃₀‖_∞ / |η₀(x_*)|  ≤  2 (1+f) sin φ₀ / sin δ ,
```
(using `‖η₀‖_∞ = M/(ρ₀ sin δ)` exactly and `|η₀(x_*)| = M/(ρ_* sin φ₀)`), and the relative
`η`-loss from the datum's edges is `≤ A · min(B1,B2)` at `q = (d/(ρ₀δ))² L/c` by (1.2).

**(i) The brief's `d = ρ₀δ = √(ν/M)`: `q = L/c` exactly, and it is not enough.**
`A = 2 sin30°/sin7.5° = 7.6613`.

| `L` | 5 | 10 | 20 | 40 | 80 | 160 | 320 | 640 |
|---|---|---|---|---|---|---|---|---|
| `q = L/c` (`c=0.734`) | 6.8 | 13.6 | 27.2 | 54.5 | 109.0 | 218.0 | 436.0 | 871.9 |
| `A·min(B1,B2)` | 70.8 | 65.5 | 56.0 | 40.9 | 21.8 | 6.22 | 0.505 | 3.6e-4 |
| target `1/L` | 0.200 | 0.100 | 0.050 | 0.025 | 0.0125 | 0.0063 | 0.0031 | 0.0016 |

Smallest `L` at which one dissipation length suffices: **`L = 590` (`c = 0.734`), `L = 854`
(`c = 0.811`)**, i.e. `log Re_E ≈ 1180 / 1708`. **The brief's `d ≳ ρ₀δ` does not close the
gap at any Reynolds number anyone will ever quote.**

**(ii) The inset that does work.** Smallest `f` with `A·min(B1,B2) ≤ 1/L` (accelerated
window `c = 0.73401`; the frozen window `c = 0.81093` needs ~16% more `f`):

| `L` | `φ₀` | `f = (ρ_*−ρ₀)/ρ₀` | `d/(ρ₀δ)` | `q` | `c₂` inflation | `c₂` |
|---|---|---|---|---|---|---|
| 10 | 45° | 0.9131 | 6.98 | 663 | 1.0694 | 1.5699 |
| 20 | 45° | 0.6590 | 5.03 | 691 | 1.0260 | 1.5062 |
| 40 | 30° | 0.4694 | 3.59 | 701 | 1.0097 | 1.4823 |
| 80 | 30° | 0.3382 | 2.58 | 728 | 1.0037 | 1.4734 |
| 160 | 30° | 0.2437 | 1.86 | 755 | 1.0014 | 1.4700 |

(`c₂` inflation `= L/(L − log(1+f))`, leading order: the datum and hence `log Re_E = 2L+O(1)`
is unchanged, only the tracked shell's outer log-range shrinks.) Asymptotic value
`8(1−√(2/3)) = 1.4680274` is recovered. Under the conditional Proposition V.2 the same `f`
values drop to `0.60, 0.43, 0.31, 0.22, 0.16`.

**(iii) V-b — the part V.1 does not do.** Write `η(X(τ),τ) = E[η₀(Y_τ)]` and split at `d`:

```
| η(X(τ),τ) − η₀(x₀) | / |η₀(x₀)|  ≤  (Λ_d/|η₀(x₀)|) · E|Z_τ|  +  A·min(B1,B2) ,
E|Z_τ| ≤ √( n ν τ (e^{2c}−1)/c ) ,     Λ_d = sup_{B(x₀,d)}|∇₅η₀| = M/r_min² ,
```
giving the **unconditional** first-moment bound `[r₀/r_min²]·δ√(5(e^{2c}−1)/L) = O(L^{-1/2})`
and, if the first moment cancels (which needs `∇²u`), the second-moment bound
`5δ²(e^{2c}−1)/(L r₀²) = O(L^{-1})`:

| `L` | 10 | 20 | 40 | 80 | 160 | 320 | 640 |
|---|---|---|---|---|---|---|---|
| first-moment (proved, `∇u` only) | 1.184 | 0.531 | 0.273 | 0.154 | 0.0924 | 0.0580 | 0.0376 |
| second-moment (needs `∇²u`) | 0.0156 | 0.0104 | 0.00658 | 0.00397 | 0.00230 | 0.00129 | 0.00070 |
| `prove-lagrangian`'s `ν τ/r²` | 0.00252 | 0.00126 | 0.00063 | 0.00031 | 0.00016 | 0.00008 | 0.00004 |

Three readings:
* the `O(L^{-1/2})` bound is still `o(1)`, so **`c₂ = 8(1−√(2/3)) + o(1)` survives even with
  `∇u` alone** — but the `o(1)` is `L^{-1/2}`, not `L^{-1}`, and is numerically vacuous
  (`1.18`) below `L ≈ 40`;
* the `O(1/L)` that `prove-lagrangian` §4(3) asserts is the right order but needs `∇²u`;
* and its constant is optimistic by `5(e^{2c}−1)/(c(1+f)²) = 6.2` — the note prices the
  instantaneous `ν/r²` and not the `e^{2c}` spreading of the diffusive cloud.
The earlier pass of this seat measured the first moment of `Z_τ` directly and found it equal
to `‖∇²b‖ ν τ²` to within `0.91–1.14` (`scripts/v3_log.txt`, `R2_first_moment`,
`m_over_predB`), confirming that the `∇²u` hypothesis is not an artefact of the proof: the
first moment really is of that size, and really does vanish only when `b` is linear.

---

## 6. Self-attack

1. **"The Feynman–Kac step needs bounded `b`; yours grows like `MLρ`."** — Handled: linear
   growth gives non-explosion, and the stopping argument in Step 0 needs only `|η| ≤ N`, not a
   gradient bound. The estimate itself never uses `‖b‖_∞`, only `‖∇b‖_∞`.
2. **"`d` should be measured in the deformed frame, and the flow map distorts it."** — No:
   `d` is the distance between *labels*, and labels do not move. This is the whole reason the
   material formulation is the right one, and it is why the estimate has no `e^{c}` in `d`.
3. **"The compressibility must appear."** — Refuted twice (§2.3), once probabilistically and
   once by the exact identity (2.1).
4. **"Aronson closes this; you are reproving a 1968 theorem."** — Checked at source (§3):
   Aronson 1968's constants depend on `‖b‖_∞` and `T`; Fabes–Stroock has no drift;
   Norris–Stroock needs uniform continuity of `A`; Porper–Eidel'man needs a time-independent
   weight. The Lagrangian reduction produces a **time-dependent** weight, which none of them
   covers, and the Eulerian route is `14–56 σ` off.
5. **"Then your constants must be terrible."** — They are, by a factor `e^{2c} ≈ 4.3` in the
   rate and by `~500` in the prefactor (measured `C' = 3.08` vs proved `4e^{2c} = 1502` at
   `c = 2.96`). That is why the required inset `f` is `O(1)` rather than `O(δ)`, and it is the
   single most improvable number in this note.
6. **"Then Gap V is closed."** — No. Only V-a. V-b is open, and §5(iii) is the honest price:
   `O(L^{-1/2})` with `∇u` alone (enough for the `+o(1)` statement, vacuous at `L ≤ 40`),
   `O(1/L)` only with a `∇²u` bound.
7. **"The window `c` is a free parameter; shrink it."** — `c = 2θ` and `θ` is the *doubling*
   requirement `∫a = log(3/2)`; it cannot be shrunk without giving up the `3/2` threshold.
   Subdividing `[0,τ]` does not help: after the first sub-window the difference `η₁−η₂` is no
   longer supported at distance `d`, so the lemma cannot be composed.
8. **"`‖∇₅b‖` should cost a dimensional factor over `‖∇u‖`."** — It does not: the lift's
   gradient has the same operator norm as `∇₃u`, exactly (§1).

---

## 7. Files

| file | what it establishes |
|---|---|
| `v1_reduction.py` / `v1_results.json` | EXACT: `div₅b = div₃u + 2a`; `Δ₅(1/r) = −1/r³`; Piola; **Lemma V.0** (residual exactly `0` in rational arithmetic); `∇₅b = diag(a,a,a,a,−2a)` |
| `v2_flowmap_bounds.py` / `v2_results.json` | the ellipticity/weight bounds (2.2), 400/400 random trials, sharpness, `det J = e^{∫div₅b}` to `9.2e-11` |
| `v3_tail_constants.py` / `v3_results.json` | the explicit constants `B1`, `B2`, the `ε`-optimisation, the `c→0` sharpness check against `erfc`, the crossover, and `q` needed per target |
| `v4_falsify.py` / `v4_log.txt` / `v4_results.json` | refuter runs R1 (n=5 SDE Monte Carlo), R2 (n=1 PDE end-to-end), R3 (Eulerian control) |
| `v5_application.py` / `v5_results.json` | §5: the brief's `d = ρ₀δ` table, the required inset, the `c₂` inflation, V-b, the Eulerian `σ` count |
| `check_constants.py` | asserts every number displayed above (`FAILURES: none`) |
| `SHA256SUMS` | computed, `shasum -a 256 -c` clean over all 36 entries |
| `scripts/` | an earlier, interrupted pass of this same seat (kept as independent corroboration: item A6 = Lemma V.0, A5 = the operator-norm identity, A11 = `∇log D ≡ 0` for the pure strain, `R1_fitted_Cprime = 3.08`, `R2` first moment `≈ ‖∇²b‖ντ²`) |
| `sources/` | Aronson 1968 (full text + page images), Norris–Stroock 1991, and the Aronson 2017 survey fetched this run |

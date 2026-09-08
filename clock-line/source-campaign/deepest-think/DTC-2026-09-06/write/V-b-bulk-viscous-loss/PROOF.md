# V-b — the bulk viscous loss of `η` on the innermost material shell over the window

Seat `write/V-b-bulk-viscous-loss`, DTC-2026-09-06.  Laws `TORMENT NEXUS/LAWS.md`.
Every number below came out of a script in this folder (`b1`…`b5`, re-asserted by
`check_constants.py`); `SHA256SUMS` is computed, never typed.  Nothing outside this folder
was written; the campaign's `d3_results.json` was read, and its sha256 is recorded in
`b3_results.json`.  Statements marked EXACT are sympy-verified in exact arithmetic.
Numerics falsify; they never prove.

---

## 0. Headline and status

**V-b closes at `O(1/L)` — with the exact constant, not an order — for the affine
(pure axisymmetric strain) motion, and it closes for the true motion under one further
hypothesis on `‖∇²u‖`.  That hypothesis cannot be removed: I exhibit a family in which
`‖∇u‖` is held fixed, `Δη₀ ≡ 0` (so every metric/Laplacian term vanishes identically), and
the loss is nonzero and exactly proportional to `‖∇²u‖`.**

Three results the campaign can use immediately:

1. **The first-order loss is not `ντ/r²(1+O(c))`; it is `ν∫₀^τ dt/r(t)²` EXACTLY**, with
   coefficient `1`, no dimensional factor `n`, and no `e^{2c}`.  For the accelerated window
   this equals `(ρ₀δ/r₀)²·I₂/L` with `I₂ = 4/5 − 16√6/135 = 0.5096901` EXACT, i.e.
   `ε_bulk·L = 0.03473454` at `f = 0`, `φ₀ = 30°`, `δ = 7.5°` — **3.5 % of the `1/L` budget.**
   `prove-lagrangian` §4(3)'s pricing `δ²/sin²φ · c/L` is **conservative by exactly
   `θ_max/I₂ = 1.4401`**, not "optimistic by 6.2" as `gap-V-aronson` §5(iii) states.
2. **The tail is far smaller than `gap-V` computed.**  In the material frame the affine
   displacement is exactly Gaussian with the *anisotropic* covariance
   `diag(2σ_y I₄, 2σ_z)`, `σ_z/σ_y = I₄/I₂ = 3.5130747`.  One dissipation length suffices
   from **`L = 24.60`** (against the `1/L` budget) or **`L = 48.34`** (against the much
   stricter budget `ε_bulk`) — not `gap-V`'s `L = 590 / 854`.  If one insets anyway, the
   needed inset is `f = 0.2597, 0.1923, 0.1419, 0.1044, 0.0766` at `L = 10…160`, versus
   `gap-V`'s `0.913…0.244` and its refuter's `1.083…0.303`; the `c₂` inflation at `L = 10`
   falls from `1.5699 / 1.5843` to **`1.5027`**.
3. **Where `∇²u` enters, exactly.**  In the label frame the operator of Lemma V.0 is
   `νG^{lk}∂_l∂_k + νV^k∂_k` with `V^k = D^{-1}∂_l(DG^{lk})`, and `V ≡ 0` whenever `X(·,t)` is
   affine in the label (for `n = 1`, `V = −X''/X'³`, so the converse holds too).  `V` is a
   contraction of `∂²_αX`, hence `O(K₂τ)` with `K₂ := ‖∇²₅b‖_∞`.  Acting on `∇η₀ = O(M/r²)`,
   against the metric term's `∇²η₀ = O(M/r³)`, the ratio is `K₂τr`.  So `∇²u` is a
   **relative** correction of size `K₂τr`, never a leading term — and the requirement is only

   > **(H-K₂)**  `K₂ ≤ C M L/ρ₀ = C‖∇u‖_∞/ρ₀` on a `√(ντ)`-neighbourhood of the trajectory:
   > *`u` has no structure finer than `ρ₀`.*  This is much weaker than the size one actually
   > expects (`K₂ ≍ M/ρ₀`, §7), and it is strictly stronger than `‖∇u‖_∞` alone, which is
   > **provably** insufficient (§6).

| claim | status |
|---|---|
| **Prop. 3.1** label operator `= νG:∇² + νV·∇`, `V ≡ 0` for affine `X` (`n=1`: iff); `‖V‖ ≤ C(n)e^{Cc}K₂τ` | **PROVED** (exact; `b1` §F) |
| **Lemma V.3** affine strain ⟹ `Z_τ` exactly Gaussian, `C_τ = 2ν∫₀^τ diag(λ^{-2}I₄, λ⁴)du` | **PROVED** (exact; MC-confirmed `b2` T0) |
| **Lemma V.3′** `½C_τ:∇²η_P(x₀) = −η_P(x₀)·ν∫₀^τ dt/r(t)²` — coefficient exactly 1, `σ_z` drops out | **PROVED** (EXACT, `b1` §E) |
| `sup_{|v|=1}|∂_v^k(1/r)| = k!/r^{k+1}` (sharp Taylor remainder) | **PROVED** (EXACT, `b1` §C) |
| **Theorem V.4** the full bound with `E_hess`, `E_4`, `E_tail` explicit | **PROVED** |
| **Corollary V.5** `ε_bulk = (ρ₀δ/r₀)² I₂/L` for the campaign window | **PROVED** |
| **`∇²u` is necessary** (no `‖∇u‖`-only bound exists) | **PROVED** (§6) + numerically confirmed to 1.7e-4 |
| **(H-K₂)** itself — a bound on `‖∇²u‖` for the actual mollified-shell field | **NOT PROVED — open, and `prove-lagrangian` L3's mode sum does not deliver it** (§7) |
| the `3.0 % / 0.1 %` split in the campaign runs | **explained, and the `3.0 %` is `~3×` discretisation** (§8) |

**FL-000 stands.  Nothing here touches the headline problem.**

---

## 1. Setting, notation, standing hypotheses

`ℝ⁵ = ℝ⁴_y × ℝ_z`, `r = |y|`, `ρ = |x|`, `φ` the polar angle from the `z`-axis, `n = 5`.
The lift of an axisymmetric no-swirl `u = u^r e_r + u^z e_z` is `b(y,z) = (u^r y/r, u^z)`;
`η = ω^θ/r` obeys

```
    ∂_t η + b·∇₅ η = ν Δ₅ η ,     Δ₅ = ∂_rr + (3/r)∂_r + ∂_zz .                    (1.1)
```

`X(α,t)` is the flow of `b`, `J = ∂X/∂α`, `D = det J`, `G = (JᵀJ)^{-1}`, `ζ(α,t) = η(X,t)`.

**Standing hypotheses.**

* **(H1)** `b` is `C^{1,1}` in `x`, continuous in `t`, `|b| ≤ B₀(1+|x|)`;
  `Γ := sup_{[0,τ]}‖∇b‖_{L^∞,op}`, `c := Γτ`;
  `K₂ := sup_{[0,τ]}‖∇²₅b‖_{L^∞}` (as the Lipschitz constant of `∇₅b` in operator norm).
  By `gap-V-aronson` §1, `‖∇₅b‖_op = ‖∇₃u‖_op` pointwise, so `Γ = ‖∇u‖_∞`.  **No such
  identity is claimed for `K₂`**: `K₂` is a property of the lift `b`, and translating it into
  a `3`-D quantity (it involves `∂²u` together with `∂(u^r/r)` terms from the `y/r` factor) is
  not done here.  Where I write `‖∇²u‖` informally I mean `K₂ = ‖∇²₅b‖_∞`.
* **(H2)** `η` is a bounded classical solution of (1.1) on `ℝ⁵×[0,τ]`, `‖η₀‖_∞ < ∞`.
* **(H3)** (the *bulk* hypothesis) there are `x₀` with `z₀ > 0`, `r₀ = r(x₀)`, and `0 < d < r₀`
  such that on the ball `B(x₀,d)` the datum coincides with the bare plateau
  `η_P(x) := −M sgn(z)/r`.  Put `R_- := r₀ − d = inf_{B(x₀,d)} r`.
* **(H4)** `ν > 0`, `τ > 0`.

**Campaign normalisation.**  `ρ₀δ = √(ν/M)` (the datum's finest feature = the dissipation
length), `L = log(R/ρ₀)`, `θ = MLt`, and by `prove-lagrangian` (4.1) at the tracked shell
`σ = 0`, `λ(θ) = (1 − κθ/2)^{-2}` with `κ = 1/2`; the doubling window ends at
`θ_max = 4(1−√(2/3)) = 0.7340137` where `λ = 3/2`, so `τ = θ_max/(ML)` and `ντ = (ρ₀δ)²θ_max/L`.
Along the tracked trajectory `r(t) = λ(t) r₀`.
The window constants (`b5`):

```
  Γ(t) = 2a(0,t) = M L λ(t)^{1/2}  ⟹  c := (sup_t Γ)·τ = √(3/2)·θ_max = 0.8989795 ,
                                       C(τ) := ∫₀^τ Γ dt = 2 log(3/2) = 0.8109302 ,
                                       Γ(0)·τ = θ_max      = 0.7340137 .
```
(the refuter of `gap-V` is right that `0.7340105` is `Γ(0)τ`, not `(sup Γ)τ`; I use the
proved value `c = 0.8989795` everywhere a Grönwall constant is needed).

---

## 2. Exact identities (all sympy-verified, `b1_identities.py`)

```
  r³ Δ₅ (1/r) = −1                       EXACT, in the (r,z) form and in Cartesian ℝ⁵     (2.1)
  Δ₅ η_P = − η_P / r²      for  η_P = −M sgn(z)/r ,  z ≠ 0                                (2.2)
  Hess(1/r) has eigenvalues { 2/r³ , (−1/r³)×3 , 0 } ;  trace = −1/r³ ; ‖·‖_op = 2/r³      (2.3)
  sup_{|v|=1} |∂_v^k (1/r)| = k!/r^{k+1} ,  attained at v = ±ŷ                             (2.4)
  Δ_y^k (1/r) = c_k r^{-(2k+1)} ,  c_1..c_4 = −1, −3, −45, −1575                           (2.5)
```

(2.4) is the sharp Taylor-remainder constant and comes from the Legendre generating
function: with `|y₀| = r`, `μ = v·ŷ`, `(r²+2sμr+s²)^{-1/2} = Σ_k s^k P_k(−μ)/r^{k+1}`, and
`|P_k| ≤ 1`.  Residuals `0` for `k = 1…6`.

**(2.2) is the whole content of "the plateau is a bulk eigenfunction".**  It says the
instantaneous logarithmic loss rate of `η` at a plateau point is exactly `−ν/r²` — no `1/σ²`
core penalty, no dimensional factor.  What §§3–5 do is to control what happens to this over
a *window*, when the field is no longer the plateau.

---

## 3. The label frame: exactly where `∇²u` sits

**Lemma V.0** (`gap-V-aronson`, re-verified there in exact arithmetic; used as an input):

```
    ∂_t ζ = (ν/D) ∂_{α_l} ( D G^{lk} ∂_{α_k} ζ )  —  divergence form in dμ = D dα,
    no drift, no zeroth-order term;  and  e^{-2c}I ≤ G ≤ e^{2c}I,  e^{-5c} ≤ D ≤ e^{5c}.
```

**Proposition 3.1 (the `∇²u` locus).  PROVED.**

> Expand Lemma V.0's operator:
> ```
>       (ν/D)∂_l(D G^{lk}∂_k ζ)  =  ν G^{lk} ∂_l∂_k ζ  +  ν V^k ∂_k ζ ,
>       V^k := D^{-1} ∂_l ( D G^{lk} ) .                                              (3.1)
> ```
> (i) `V ≡ 0` identically whenever `X(·,t)` is affine in `α`.  For `n = 1`,
> `V = −X''/X'³` exactly, so there `V ≡ 0` ⟺ `X` affine.
> (ii) `V` is a contraction of `∂²_αX` alone; and `‖∂²_α X(·,t)‖ ≤ K₂ t e^{3c}`, so
> `‖V‖ ≤ C(n) e^{C'c} K₂ τ`.

*Proof.* (i) If `X = A(t)α + a(t)` then `J = A(t)`, `D` and `G` are `α`-independent, so every
`∂_l` in (3.1) annihilates them.  Verified symbolically for a sheared dilation in `n = 2` and
for a generic rational `5×5` matrix in `n = 5`: `V = 0` exactly (`b1` §F).  For `n = 1`,
`J = D = X'`, `G = X'^{-2}`, so `V = (1/X')∂_α(1/X') = −X''/X'³` — residual `0` in exact
arithmetic — whence the converse in that case.  (For `n ≥ 2` I prove only the forward
implication; `V ≡ 0` is then a genuine PDE constraint that I have not shown forces affineness.)
(ii) `V^k` is built from `∂_lG` and `∂_l log D`, both contractions of `∂_αJ = ∂²_αX` with
`J^{-1}`.  Exhibited on `X = (α₁+εα₂², α₂+εα₁α₂)`: `V = (−2ε + O(ε²), O(ε²))`, i.e. exactly
first order in the curvature of the map (`b1` §F).
Finally `P := ∂²X/∂α²` obeys `∂_tP^m_{ij} = (∇²b)^m_{pq}J^q_jJ^p_i + (∇b)^m_pP^p_{ij}`, `P(0)=0`,
and `‖J‖ ≤ e^c`; Grönwall gives `‖P(t)‖ ≤ K₂ e^{2c}∫₀^t e^{Γ(t−s)}ds ≤ K₂ t e^{3c}`. ∎

**Consequence (the answer to the brief's question).**  Insert `ζ ≈ ζ₀ = η₀` into (3.1) at the
tracked point.  The second-order term is `νG:∇²η₀ = −ν η₀/r²·(1+O(c))` — **the metric, and
only the metric.**  The first-order term is `νV·∇η₀ = O(νK₂τ)·O(M/r²)`.  Their ratio is

```
        |ν V·∇η₀| / |ν G:∇²η₀|  ≍  K₂ τ r .                                          (3.2)
```

So the brief's hope — "the metric enters, not `∇²u`" — is **half right**: the metric carries
the leading term, `∇²u` carries a correction of relative size `K₂τr`, and that correction is
not removable (§6).  §5 proves the same statement probabilistically, with constants.

---

## 4. The affine case: the loss is EXACTLY `ν∫dt/r²`

**Lemma V.3 (affine displacement).  PROVED.**

> Let `b(x,t) = S(t)x` with `S(t) = a(t)·diag(1,1,1,1,−2)` (the uniform axisymmetric strain,
> which by `prove-lagrangian` L3 is exactly the `O(M log)` part of `u`).  Let `Y` be the
> backward diffusion of Theorem V.1 Step 0, `dY_s = −b(Y_s,τ−s)ds + √(2ν)dW_s`, `Y_0 = X(τ)`,
> and `x_s := X(τ−s)` the reverse characteristic, `Z_s := Y_s − x_s`.  Then `Z_τ` is
> **exactly Gaussian, mean zero**, with
> ```
>     C_τ = 2ν ∫₀^τ diag( λ(u)^{-2}·I₄ , λ(u)^{4} ) du  =: diag(2σ_y I₄ , 2σ_z) ,
>     λ(u) = exp ∫₀^u a .                                                             (4.1)
> ```

*Proof.*  `Z` obeys `dZ = −S(τ−s)Z ds + √(2ν)dW`, a linear SDE with **deterministic**
coefficients and `Z_0 = 0`; hence `Z_s = √(2ν)Ψ_s∫₀^sΨ_p^{-1}dW_p` with `Ψ̇ = −S(τ−s)Ψ`,
`Ψ_0 = I` — a Wiener integral of a deterministic integrand, so Gaussian with mean `0` and
`C_τ = 2ν∫₀^τ (Ψ_τΨ_p^{-1})(Ψ_τΨ_p^{-1})ᵀ dp`.  With `μ(s) = ∫_{τ−s}^{τ}a`,
`Ψ_s = diag(e^{-μ(s)}×4, e^{2μ(s)})`, so `Ψ_τΨ_p^{-1} = diag(λ(τ−p)^{-1}×4, λ(τ−p)²)`;
substituting `u = τ−p` gives (4.1). ∎
*(Cross-check `b2` T0: 200 000 paths, 800 steps, `ν = 2e-3`, `τ = 0.4`, `a₀ = 1` reproduce
the five diagonal entries of `C_τ` to `5.1e-5 … 3.6e-3` relative — within the Monte-Carlo
standard error `√(2/N) = 3.2e-3` — with `max|off-diagonal|/tr = 7.4e-4` and every component
mean within `2.02` standard errors of `0`.)*

**Lemma V.3′ (the identification).  EXACT.**

```
    ½ C_τ : ∇²η_P(x₀)  =  M σ_y / r₀³  =  − η_P(x₀) · σ_y/r₀²  =  − η_P(x₀) · ν∫₀^τ dt/r(t)² ,
```
because `η_P` has no `z`-dependence off `{z=0}`, so `∂_z²η_P = 0` and **`σ_z` drops out of the
bulk term identically**; and `Δ_yη_P = −η_P/r²` by (2.2); and `r(u) = λ(u)r₀` turns
`σ_y/r₀² = (ν/r₀²)∫λ^{-2}du` into `ν∫du/r(u)²`.  Residual `0` in exact arithmetic, and
`∂(½C:∇²η_P)/∂σ_z = 0` identically (`b1` §E).

**Numerical confirmation of both, on the exact instrument** (`b2`, quadrature; the 4 `y`
directions collapse to a noncentral-`χ₄` density, the `z` direction to a normal, so
`E[η₀(x₀+Z)]` is a 2-D quadrature with no PDE solve and no Monte-Carlo error):

| `s := σ_y/r₀²` | `1e-2` | `3e-3` | `1e-3` | `3e-4` | `1e-4` | `3e-5` | `1e-5` |
|---|---|---|---|---|---|---|---|
| `rel. loss / s` | 1.015826 | 1.004569 | 1.001508 | 1.000451 | 1.000150 | 1.000045 | 1.000015 |
| `(rel. loss − s)/s²` | 1.5826 | 1.5231 | 1.5076 | 1.5023 | 1.5007 | 1.5002 | 1.4996 |

The leading coefficient is `1` and the next is `3/2`, matching the formal series
`Φ(s) = 1 − s − (3/2)s² − …` generated by (2.5) (`c₂/2! = 3/2`).  And the `σ_z`-independence
is exact to 11 digits over a `1000×` range in `σ_z/σ_y` (`b2` T3: `1.001507566422e-03` at
`σ_z/σ_y = 0.1, 1, 5, 25`, `1.001507599e-03` at `100`, the last digit being quadrature drift).

---

## 5. THEOREM V.4 — the bound, with all constants

**Moment lemma.**  Under (H1), for the backward `Z` of §4 with a general `b`:
`d|Z|² = −2ZᵀA_sZ\,ds + 2√(2ν)Z·dW + 2nν\,ds` with `A_s = ∫₀¹∇b(x_s+θZ_s,τ−s)dθ`,
`‖A_s‖ ≤ Γ`, so

```
    v(s) := E|Z_s|² ≤ nν(e^{2Γs}−1)/Γ ≤ nντ(e^{2c}−1)/c =: V₁   for s ≤ τ .          (5.1)
```

**Coupling lemma (this is where `∇²u` is paid for, once).**  Let `Z^a` solve the *affine*
equation `dZ^a = −A^a_sZ^a ds + √(2ν)dW`, `A^a_s := ∇b(x_s,τ−s)`, driven by the **same** `W`,
`Z^a_0 = 0`.  Since `A_s − A^a_s = ∫₀¹[∇b(x_s+θZ_s)−∇b(x_s)]dθ` has norm `≤ ½K₂|Z_s|`,

```
    d(Z−Z^a) = −A^a_s(Z−Z^a)ds − (A_s−A^a_s)Z_s ds     ⟹  pathwise
    |Z_τ − Z^a_τ| ≤ e^{c}·(K₂/2)∫₀^τ |Z_s|² ds =: W ,   E[W] ≤ ½ K₂ e^{c} τ V₁ .     (5.2)
```
`Z^a` is exactly Gaussian, mean zero, covariance `C_τ = 2ν∫₀^τΦ_{τ,p}Φ_{τ,p}ᵀdp`
(`Φ` the propagator of `−A^a`); for the pure axisymmetric strain this is (4.1).

> **THEOREM V.4.**  Assume (H1)–(H4).  Write `s_C := −½C_τ:∇²η_P(x₀)/η₀(x₀)`,
> `𝒫 := P(|Z_τ| > d/2) + P(|Z^a_τ| > d/2) + P(|Z^a_τ| ≥ d)`, and
> `N := ‖η₀‖_∞ r₀/M ( = ‖η₀‖_∞/|η₀(x₀)| )`.  Then
>
> ```
>   | η(X(τ),τ)/η₀(x₀) − (1 − s_C) |  ≤  E_hess + E_4 + E_tail ,
>
>   E_hess = (r₀/R_-²)·½K₂e^{c}τV₁            +  2N·(2E[W]/d)          ,
>   E_4    = (r₀/R_-⁵)·[ (tr C_τ)² + 2 tr(C_τ²) ]                      ,
>   E_tail = 2N·𝒫 + 𝒫 + Σ_{k=1}^{3} r₀^{-k} (E|Z^a_τ|^{2k})^{1/2} 𝒫^{1/2}  .
> ```
> All three depend on `b` only through `c = Γτ` and `K₂τ`.

*Proof.*  By Theorem V.1 Step 0 (proved in `gap-V-aronson`, refuter-confirmed),
`η(X(τ),τ) = E[η₀(Y_τ)] = E[η₀(x₀+Z_τ)]`.

*Step 1 (replace `Z` by `Z^a`).*  On the event `𝒢 := {|Z_τ| ≤ d/2} ∩ {W ≤ d/2}` both
`x₀+Z_τ` and `x₀+Z^a_τ` lie in `B(x₀,d)`, where `η₀ = η_P` and `|∇η_P| ≤ M/R_-²`; so
`|η₀(x₀+Z_τ) − η₀(x₀+Z^a_τ)|1_𝒢 ≤ (M/R_-²)W`.  Off `𝒢` bound by `2‖η₀‖_∞`, and
`P(𝒢^c) ≤ P(|Z_τ|>d/2) + P(W>d/2) ≤ P(|Z_τ|>d/2) + 2E[W]/d` (Markov).  Dividing by
`|η₀(x₀)| = M/r₀` gives the two pieces of `E_hess` (using (5.2) for `E[W]`).

*Step 2 (the affine average).*  Let `P₃(ζ)` be the third-order Taylor polynomial of `η_P` at
`x₀`.  On `{|Z^a| < d}`, `η₀ = η_P` and, by (2.4) with `k = 4`, the remainder obeys
`|η_P(x₀+ζ) − P₃(ζ)| = |(1/4!)∂_v^4η_P(ξ)||ζ|^4 ≤ M|ζ|⁴/R_-⁵`.  Because `Z^a` is Gaussian and
mean zero, `E[P₃(Z^a)] = η_P(x₀) + ½C_τ:∇²η_P(x₀) = η₀(x₀)(1 − s_C)` exactly (the linear and
cubic terms have zero expectation).  Hence
```
  E[η₀(x₀+Z^a)] = η₀(x₀)(1−s_C) + E[(η_P−P₃)(Z^a)1_{<d}]
                   − E[P₃(Z^a)1_{≥d}] + E[η₀(x₀+Z^a)1_{≥d}] ,
```
and `|E[(η_P−P₃)1_{<d}]| ≤ (M/R_-⁵)E|Z^a|⁴ = (M/R_-⁵)[(tr C)²+2tr(C²)]` (Isserlis),
`|E[η₀1_{≥d}]| ≤ ‖η₀‖_∞ P(|Z^a|≥d)`, and by Cauchy–Schwarz and (2.4),
`|E[P₃1_{≥d}]| ≤ M P/r₀ + M Σ_{k=1}^{3}(E|Z^a|^{2k})^{1/2}P^{1/2}/r₀^{k+1}` (the `k = 0` term
needs no Cauchy–Schwarz).  Dividing by `M/r₀` gives `E_4` and `E_tail`. ∎

> **COROLLARY V.5 (the campaign window).**  For the pure axisymmetric strain with
> `λ(θ) = (1−κθ/2)^{-2}`, `κ = 1/2`, `θ = MLt`, `τ = θ_max/(ML)`, `ρ₀δ = √(ν/M)`,
> `r₀ = (1+f)ρ₀ sin φ₀`:
> ```
>    s_C  =  ν ∫₀^τ dt/r(t)²  =  (ρ₀δ/r₀)² · I₂ / L      EXACTLY,
>    I₂ = ∫₀^{θ_max} λ^{-2} dθ = 4/5 − 16√6/135 = 0.5096901045590 ,
>    I₄ = ∫₀^{θ_max} λ^{4}  dθ = −4/7 + 27√6/28 = 1.7905793948266 ,   σ_z/σ_y = 3.5130747 .
> ```
> At `f = 0`, `φ₀ = 30°`, `δ = 7.5°`: `s_C·L = 0.03473454` for every `L`.

| `L` | 10 | 20 | 40 | 80 | 160 | 320 | 640 |
|---|---|---|---|---|---|---|---|
| `ε_bulk = s_C` | 3.4735e-3 | 1.7367e-3 | 8.6836e-4 | 4.3418e-4 | 2.1709e-4 | 1.0855e-4 | 5.4273e-5 |
| `prove-lagrangian`'s `ντ/r²` | 5.0022e-3 | 2.5011e-3 | 1.2505e-3 | 6.2527e-4 | 3.1264e-4 | 1.5632e-4 | 7.8159e-5 |
| `target 1/L` | 0.1 | 0.05 | 0.025 | 0.0125 | 0.00625 | 0.003125 | 0.0015625 |

The ratio of rows 2 and 1 is `θ_max/I₂ = 1.4401` at every `L`: **`prove-lagrangian` §4(3) is
conservative by 44 %, not optimistic by 6.2×.**  `gap-V-aronson` §5(iii)'s "constant
optimistic by `5(e^{2c}−1)/(c(1+f)²) = 6.2`" comes from bounding
`E[Zᵀ∇²η_PZ] ≤ ‖∇²η_P‖_{op}E|Z|²` instead of contracting the trace: `‖∇²η_P‖_op = 2M/r³`
while the *contraction against the actual covariance* is `M σ_y/r³` — the anisotropy of
`C_τ` and the tracelessness of `Hess(1/r)` off the radial direction do the whole job.

**The tail, sized in the material frame** (`b5`; `σ_y = (ρ₀δ)²I₂/L`, `σ_z = (ρ₀δ)²I₄/L`; the
inner edge's outward normal at `φ₀ = 30°` has variance
`V_n = 2(ρ₀δ)²(I₂sin²φ₀ + I₄cos²φ₀)/L`; `A := ‖η₀‖_∞/|η₀(x_*)| + 1 = (1+f)sinφ₀/sinδ + 1`):

| `L` | 5 | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|---|
| `A·P(edge)` at `d = ρ₀δ` | 1.139 | 0.4868 | 0.08890 | 2.965e-3 | 3.299e-6 | 4.084e-12 |

**One dissipation length suffices from `L = 24.603` (vs `1/L`) or `L = 48.337` (vs the far
stricter `ε_bulk`)** — against `gap-V`'s `L = 590 / 854`.  If one insets anyway so that the
edge tail is below the *bulk* term:

| `L` | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|
| `f` (this seat) | 0.2597 | 0.1923 | 0.1419 | 0.1044 | 0.0766 |
| `f` (`gap-V` §5(ii)) | 0.9131 | 0.6590 | 0.4694 | 0.3382 | 0.2437 |
| `f` (its refuter, corrected `c`) | 1.0829 | 0.7922 | 0.5673 | 0.4149 | 0.3031 |
| `c₂ = 8(1−√(2/3))·L/(L−log(1+f))` | 1.50272 | 1.48105 | 1.47291 | 1.46985 | 1.46871 |

**Caveat, stated plainly.**  The tail table is computed with the *exact affine* Gaussian.
For the true motion the coupling (5.2) transfers it, but only up to a large-deviation control
of `W` that I have **not** written; the rigorous fallback for a general `b` is Theorem V.1's
`min(B1,B2)`, i.e. `gap-V`'s larger `f`.  What is unconditional is that both tables give
`c₂ = 8(1−√(2/3)) + O(1/L)`.

**The `∇²u` correction, sized** (`b5`; `E_hess/ε_bulk = ½e^{c}n(e^{2c}−1)/c · K₂τr₀`,
prefactor `= 34.42001` at `c = 0.8989795`, `n = 5`; `K₂ = \hat K₂ M/ρ₀`):

| `\hat K₂` | `L = 10` | `40` | `160` | `640` |
|---|---|---|---|---|
| 1 | 1.2632 | 0.3158 | 0.0790 | 0.0197 |
| 3 | 3.7897 | 0.9474 | 0.2369 | 0.0592 |
| 10 | 12.6324 | 3.1581 | 0.7895 | 0.1974 |

Largest `\hat K₂` for which `E_hess ≤ ε_bulk`: `0.792 (L=10), 3.166 (40), 12.666 (160),
50.663 (640)` — i.e. the admissible `K₂` grows like `L`, so **`K₂ ≤ C M L/ρ₀` keeps
`E_hess = O(ε_bulk) = O(1/L)`** for every `C`, and `K₂ ≤ C M/ρ₀` makes it `O(1/L²)`.

---

## 6. `∇²u` cannot be removed — a proof, and a measurement

**Proposition 6.1.  PROVED.**  There is no bound of the form
`|η(X(τ),τ) − η₀(x₀)| ≤ F(‖∇u‖_∞τ, ν, τ, ‖∇²η₀‖, …)` that is uniform over drifts with a given
`‖∇u‖_∞`, even for data with `Δη₀ ≡ 0`.

*Proof.*  Take `n = 1`, `η₀(x) = x` (so `Δη₀ ≡ 0`: every metric/Laplacian term in §§3–5
vanishes identically and `s_C = 0`).  Then `η(X(τ),τ) = E[Y_τ] = x₀ + E[Z_τ]`, so the loss is
exactly the first moment `m_τ := E[Z_τ]`.  From `dZ = −A_sZ ds + √(2ν)dW`,
`A_s = b'(x_s) + ½b''(x_s)Z_s + O(K₃Z_s²)`, one gets `m'_s = −b'(x_s)m_s − ½b''(x_s)E[Z_s²]`
to leading order in `ν`, `m_0 = 0`, hence

```
   m_τ = −½ ∫₀^τ exp(−∫_s^τ b'(x_u)du) · b''(x_s) · v_s ds  + O(ν²) ,
   v_s = 2ν ∫₀^s exp(−2∫_p^s b'(x_u)du) dp ,     x_s = X(τ−s) .                      (6.1)
```
Now the family `b_κ(x) = βκ sin(x/κ)` has `‖b'_κ‖_∞ = β` **independent of `κ`** and
`‖b''_κ‖_∞ = β/κ`, and (6.1) gives `m_τ ∝ 1/κ ≠ 0`.  Letting `κ → 0` at fixed `β` makes
`m_τ` arbitrarily large while `‖∇u‖_∞` is fixed, so the supremum of the loss over drifts of
given `‖∇u‖_∞` is `+∞`. ∎

**Measurement** (`b4`, a Fourier-spectral solve of the PDE itself — the probabilistic
representation is never used, so this is an independent instrument).  `β = 1.6`, `τ = 0.5`,
label `x₀ = 0.8κ`:

| `κ` | `‖b''‖ = β/κ` | measured loss | `loss × κ` | prediction (6.1) | meas/pred |
|---|---|---|---|---|---|
| 1 | 1.6 | 2.325352e-4 | 2.325352e-4 | 2.324953e-4 | 1.0002 |
| 2 | 0.8 | 1.162527e-4 | 2.325053e-4 | 1.162477e-4 | 1.0000 |
| 4 | 0.4 | 5.812446e-5 | 2.324978e-4 | 5.812384e-5 | 1.0000 |
| 8 | 0.2 | 2.906200e-5 | 2.324960e-4 | 2.906192e-5 | 1.0000 |
| 16 | 0.1 | 1.453097e-5 | 2.324955e-4 | 1.453096e-5 | 1.0000 |

`loss × κ` is constant to `1.7e-4` relative across a `16×` sweep at fixed `‖∇u‖_∞ = β`:
**the loss is exactly proportional to `‖∇²u‖` and vanishes with it.**  A `ν` sweep gives
`loss/ν = 0.1163270, 0.1162874, 0.1162676, 0.1162576, 0.1162527` at
`ν = 8e-3 … 5e-4` — proportional to `ν`, as (6.1) says.  Resolution control:
`N = 512…4096`, `nt = 2000…16000` change the answer by `1.13e-9` relative
(`2.3253515702e-4, 2.3253515703e-4, 2.3253515692e-4, 2.3253515718e-4`).

The agreement of an independently derived leading-order formula (6.1) with a
representation-free PDE solve to 4–5 digits is also a check on the Feynman–Kac step of
Theorem V.1 and on the whole first-moment analysis of §5.

---

## 7. What (H-K₂) would take — the honest open item

(H-K₂) is `K₂ = ‖∇²₅b‖_{L^∞}` on a `√(ντ)`-neighbourhood of the trajectory, `≤ C M L/ρ₀`.
Two remarks:

1. **The `O(M log)` part of `u` contributes nothing to `K₂`.**  By `prove-lagrangian` L3 the
   entire `O(ML)` velocity is the *uniform* strain `A(ρ)(r,−2z)` with `A = (M/2)log(R/ρ)`,
   which is affine to leading order; its second derivative is `O(M/ρ)`, not `O(ML/ρ)`,
   because only `A′ = −M/(2ρ)` survives.  So one expects `K₂ ≍ M/ρ₀`, i.e. `\hat K₂ = O(1)`,
   which by the table in §5 is already inside the admissible range for `L ≥ 13`.
2. **L3's mode sum does not prove it.**  L3 controls `Σ_{l≥3}|H_l|‖C_l^{3/2}‖_∞/(l(l+3)−4)`
   with `|H_l| ~ l^{-3/2}` and `‖C_l^{3/2}‖_∞ = C_l^{3/2}(1) ≍ l²`; the sum converges for the
   *velocity* but the corresponding sum for `∇²u` gains two more powers of `l` and
   **diverges** for the bang-bang profile.  So (H-K₂) genuinely needs the mollification:
   `∇u` is a Calderón–Zygmund operator on `ω`, so `‖∇²u‖_∞ ≲ ‖∇ω‖_{C^α}/α`, and the datum's
   mollification scale `ε` must be held fixed as `L → ∞` (it is, in the campaign's family)
   *and* `‖∇ω(t)‖_{C^α}` must be propagated over `[0,τ]` — a standard but unwritten estimate
   which would cost `e^{Cc}` and nothing in `L`.

**This is the one genuine gap left in V-b.**  It is a hypothesis about the *velocity field*,
not about the diffusion, and it is much weaker than what `gap-V-aronson` and its refuter
describe ("needs a `∇²u` bound to kill the first moment") — the first moment does not have to
be killed, only shown to be `O(ML/ρ₀)·ντ²/r`.

---

## 8. Confrontation with the campaign runs (`b3`)

The campaign's `d3_results.json` (sha256 recorded in `b3_results.json`) tracks three material
points at `s₀ = 1.5, 2, 3 ρ₀`, `φ = 45°`, on the mollified taper datum
`ω^θ = −M tanh(sinφ/sinδ) tanh(cosφ/w) Θ(ρ)`, `Θ = ½[tanh((ρ−ρ₀)/w₀) − tanh((ρ−R)/w₁)]`,
`w₀ = ρ₀/4`, `w = 0.12`, `δ = 7.5°`.  (The analytic datum used below reproduces the
campaign's own `build_taper` array — imported read-only — to `8.807e-4` maximum relative
difference over `11 376` bulk cells at `h = 0.0625`, the residual being their sub-cell
averaging; `b3` T5.)  `prove-duhamel` §5 reports "at `Re₀ = 100` the tracer at
`1.5ρ₀` loses 3.0 % while the one at `3ρ₀` loses 0.1 %".

**(a) The datum's own structure.**  The *exact* bulk rate `Δ₅η₀/η₀` at the tracer's own
starting point, against the plateau value `−1/r₀²` (`b3` T1):

| `s₀/ρ₀` | 1.25 | 1.5 | 2.0 | 2.5 | 3.0 | 4.0 | 6.0 |
|---|---|---|---|---|---|---|---|
| `(Δ₅η₀/η₀)/(−1/r₀²)` | 4.351 | 2.036 | 1.041 | 1.005 | 1.004 | 1.004 | 1.004 |

So the plateau identity (2.2) is exact to `0.4 %` from `2.5ρ₀` outward, and is off by `2.04×`
at `1.5ρ₀` — the point sits **two mollification widths** from the inner edge.

**(b) The exact affine prediction versus the measurement.**  Feeding the run's own `ν`, its
own measured stretch `λ(t) = r(t)/r₀`, and the analytic datum into the exact Gaussian
instrument of §4 (`b3` T3):

| run | `h` | `ν` | `s₀` | `s_C` | exact loss | measured | meas/exact |
|---|---|---|---|---|---|---|---|
| N4 | 0.125 | 0.01 | 1.5 | 2.348e-3 | 9.975e-3 | 3.008e-2 | **3.015** |
| N4 | 0.125 | 0.01 | 2.0 | 1.407e-3 | 1.564e-3 | 3.294e-3 | 2.106 |
| N4 | 0.125 | 0.01 | 3.0 | 6.805e-4 | 6.855e-4 | 7.389e-4 | **1.078** |
| N4-fine | 0.0625 | 0.01 | 1.5 | 2.348e-3 | 9.949e-3 | 1.511e-2 | **1.518** |
| N4-fine | 0.0625 | 0.01 | 2.0 | 1.407e-3 | 1.563e-3 | 1.675e-3 | 1.072 |
| N4-fine | 0.0625 | 0.01 | 3.0 | 6.800e-4 | 6.850e-4 | 5.765e-4 | 0.842 |
| N4-Re400 | 0.125 | 0.0025 | 3.0 | 1.682e-4 | 1.693e-4 | 1.946e-4 | 1.149 |
| N4-Re25 | 0.125 | 0.04 | 3.0 | 2.897e-3 | 2.934e-3 | 2.201e-3 | 0.750 |

**At `3ρ₀` the theorem reproduces the measured `η`-loss to within the solver's own accuracy,
across a 16× range in `ν`** (ratios `1.078, 0.842, 1.149, 0.750`).

**(c) The `3.0 %` is mostly discretisation.**  The absolute excess `measured − exact` at
`s₀ = 1.5` is `2.010e-2` (N4, `h = 0.125`) and `5.157e-3` (N4-fine, `h = 0.0625`): halving `h`
divides it by **3.90** (and by **15.39** at `s₀ = 2.0`).  Over a `16×` range in `ν` at fixed
`h` the same excess changes by only `1.82×` (`3.285e-2, 2.010e-2, 1.808e-2` at
`ν = 0.04, 0.01, 0.0025`) while the physical loss changes by `21.6×`.  A viscous effect scales
with `ν`; a discretisation error does not, and shrinks with `h`.  Both signatures say the
same thing.

**Corrected reading of the campaign's own number.**  The physically correct statement is:
*at `1.5ρ₀` the loss is `4.25×` the plateau rate because the point sits two mollification
widths from the inner edge (a datum-structure effect, exactly computable); at `3ρ₀` it is the
plateau rate to `0.4 %`.  The `3.0 %` figure quoted in `prove-duhamel` §5 and repeated in
`lower/SYNTHESIS.md` is inflated roughly `3×` by the `h = 0.125` discretisation; the converged
value is `≈ 1.0 %`.*  The qualitative conclusion drawn from it — track one or two `ρ₀`
further out, at `O(1/L)` cost in `c₂` — stands, and §5's inset table prices it much more
cheaply than `gap-V` did.

---

## 9. Self-attack (pre-declared falsification rules, fixed before any number was read)

* **R1** — *Lemma V.3 is refuted if a Monte-Carlo of the backward SDE with the affine strain
  disagrees with `C_τ` by more than `4×` the Monte-Carlo standard error on any entry.*
  Ran: worst relative error `3.63e-3` against s.e. `3.16e-3`, i.e. `1.15 σ`; every mean within
  `2.02 σ` of zero; `max|off-diagonal|/tr = 7.4e-4`.  **Not violated.**  (This rule *can*
  fire: the predicted `σ_z/σ_y = 3.59` for that field is a `3.6×` anisotropy that a wrong
  propagator would miss by orders of magnitude.)
* **R2** — *Lemma V.3′ is refuted if the exact quadrature's relative loss divided by
  `s = σ_y/r₀²` fails to converge to `1`, or if it depends on `σ_z`.*  Ran: `→ 1` with
  residual `1.5s²`, and `σ_z`-independent to 11 digits over `σ_z/σ_y ∈ [0.1,100]`.
  **Not violated.**  (Fires if the coefficient were `n`, `2`, `e^{2c}`, or anything but 1.)
* **R3** — *Proposition 6.1 is refuted if `loss × κ` is not constant, or if the loss does not
  vanish as `‖b''‖ → 0` at fixed `‖b'‖`.*  Ran: constant to `1.7e-4` over `κ ∈ [1,16]`.
  **Not violated** — and this is the rule that would have saved the "no `∇²u` needed" hope
  had it been true; it says the hope is false.
* **R4** — *the confrontation with the runs is non-diagnostic unless refinement moves the
  measurement toward the prediction.*  It does: `3.015 → 1.518` at `s₀ = 1.5`,
  `2.106 → 1.072` at `s₀ = 2.0`.  **Control fires.**

Attacks I made on my own statement and what happened:

1. *"Your affine reduction assumes the velocity IS the uniform strain, which it is not."* —
   Correct, and that is exactly what §5's coupling lemma pays for: `Z` is compared with the
   affine `Z^a` **pathwise**, and the entire price is one term of size `½K₂e^cτV₁`.  The
   affine reference in Theorem V.4 is `A^a_s = ∇b(x_s,τ−s)`, the true Jacobian along the
   reverse characteristic, not the idealised strain; Corollary V.5 then specialises.
2. *"`E_4` has `R_-⁵` in it; as `d → r₀` it blows up."* — Yes.  `R_- = r₀ − d` and the bound
   degrades if the ball reaches the axis.  For the campaign `d ≤ f ρ₀ ≪ r₀`, so
   `r₀/R_- ≤ 1 + O(f)`.  Stated as a hypothesis, not hidden.
3. *"Odd moments vanish only because `Z^a` is Gaussian; for the true `Z` the cubic term is
   `O(s^{3/2})`, larger than `O(s²)`."* — True, and that is why the proof does **not** Taylor
   the true `Z`: the coupling replaces `Z` by an exactly Gaussian `Z^a` first, and the whole
   non-Gaussianity is charged once, to `E_hess`.
4. *"`σ_z` is bigger than `σ_y` by `3.5×`; you cannot ignore the `z` spread."* — For the
   *bulk* term I do not ignore it, I prove it drops out identically (`∂_z²η_P ≡ 0` off the
   equator; `b1` §E residual `0`, `b2` T3 to 11 digits).  For the *tail* it is carried
   explicitly: the edge-normal variance in §5 uses `I₂sin²φ₀ + I₄cos²φ₀`, and the `z`
   direction is the binding one at `φ₀ = 30°`.
5. *"You claim `L_min = 24.6` where `gap-V` claims 590.  One of you is wrong."* — Both are
   right about their own object.  `gap-V` uses Theorem V.1's `min(B1,B2)`, whose rate is
   `0.0083–0.0253` per unit `q`; the exact affine Gaussian rate is `1/2` per unit
   `d²/Var`.  The factor `20–60` in the rate is the `e^{c}‖Φ_τ‖` step that V.1 is forced into
   because `Φ_τΦ_r^{-1}` is not adapted.  For the affine motion that step is unnecessary.
   My number is conditional on affineness (see the caveat in §5); theirs is unconditional.
6. *"Then V-b is closed."* — **No.** (H-K₂) is not proved (§7), and the tail table is affine.
   What is closed: the identification of the leading term with the exact constant, the
   structural location of `∇²u`, its necessity, and the whole bound modulo (H-K₂).
7. *"Your campaign comparison uses the campaign's own trajectories, so it is not
   independent."* — It uses their `r(t)` (a measured input) and their `ν`; the predicted loss
   is computed by an instrument (exact Gaussian quadrature on an analytic datum) that shares
   no code with their solver.  The `h`-refinement and `ν`-sweep controls are what make it
   diagnostic, and they are theirs, not mine.
8. *"`Γ = 2a(0,t)` is an unstated hypothesis (the refuter of `gap-V` §5)."* — It is inherited,
   and I flag it: every `c`-dependent constant here (`E_hess`, `E_tail`) rides on it.  The
   `ε_bulk` identity of Corollary V.5 does **not** — it uses only `λ(t)`, i.e. the measured
   stretch of the tracked shell, so the headline number `0.03473454/L` survives any
   `O(1)` error in `Γ`.

---

## 10. Status per step

| step | status |
|---|---|
| (2.1)–(2.5) exact identities | **PROVED** (sympy, exact) |
| Prop. 3.1 (affine `X` ⟹ `V ≡ 0`; `n=1` converse; `‖V‖ ≤ C e^{Cc}K₂τ`) | **PROVED** |
| (3.2) the ratio `K₂τr` | **PROVED** (from Prop. 3.1) |
| Lemma V.3 (exact Gaussian, exact `C_τ`) | **PROVED**; MC cross-check at `1.15σ` |
| Lemma V.3′ (`½C:∇²η_P = −η_P ν∫dt/r²`) | **PROVED, EXACT**; quadrature to 5 digits |
| Moment lemma (5.1) | **PROVED** |
| Coupling lemma (5.2) | **PROVED** (pathwise) |
| Theorem V.4 | **PROVED**, all constants explicit |
| Corollary V.5 (`ε_bulk = (ρ₀δ/r₀)²I₂/L`) | **PROVED, EXACT** |
| Prop. 6.1 (`∇²u` necessary) | **PROVED**; measured to `1.7e-4` |
| leading-order first-moment formula (6.1) | **DERIVED, leading order in ν**; matches a representation-free PDE solve to 4–5 digits |
| the tail table of §5 | **PROVED for affine `b`**; transfer to general `b` needs a large-deviation bound on `W` — **not written** |
| **(H-K₂)** `‖∇²u‖ ≤ C M L/ρ₀` | **NOT PROVED — the open item** |
| campaign confrontation (§8) | **measured**; numerics falsify only |

**What V-b now needs, in one line:** a bound `K₂ = ‖∇²₅b(·,t)‖_{L^∞(N_τ)} ≤ C M L/ρ₀` for
`t ∈ [0,τ]` on the mollified taper family — equivalently `‖∇ω(t)‖_{C^α}` propagated over the
window at fixed mollification scale.  Nothing else.

---

## 11. Files

| file | what it establishes |
|---|---|
| `b1_identities.py` / `b1_results.json` | (2.1)–(2.5); `Hess(1/r)` spectrum; `∂_v^k = k!P_k(−μ)/r^{k+1}`; `½C:∇²η_P` identification (residual `0`) and `σ_z`-independence; Prop. 3.1 (`V = 0` for affine `n=2,5`, `V = −2ε+O(ε²)` for a curved map) |
| `b2_gaussian_exact.py` / `b2_results.json` | T0 Monte-Carlo check of `C_τ`; T1 the exact quadrature `rel.loss/s → 1`, `(rel.loss−s)/s² → 3/2`; T3 `σ_z`-independence to 11 digits |
| `b3_datum_and_runs.py` / `b3_results.json` | T5 the analytic datum vs the campaign's array (`8.8e-4`); exact `Δ₅η₀/η₀` for the campaign datum; the exact affine prediction vs the measured `η`-loss on 6 runs; T4 the `h`- and `ν`-signature of the residual |
| `b4_hessian_necessity.py` / `b4_results.json` | the `κ`-sweep at fixed `‖∇u‖` (`loss×κ` constant to `1.7e-4`), the `ν`-sweep, the resolution control, and (6.1) verified against a spectral PDE solve |
| `b5_sizing.py` / `b5_results.json` | `I₂`, `I₄` exact; the window constants `c`, `C(τ)`, `Γ(0)τ`; the `ε_bulk` table; the tail and inset tables; `L_min`; the `E_hess/ε_bulk` table and the `\hat K₂` thresholds |
| `check_constants.py` | re-asserts every number displayed above against the results files |
| `SHA256SUMS` | computed over this folder |

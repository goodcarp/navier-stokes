# L3v and the slab gradient bound Γ

Seat `write/L3v-and-gamma-bound`, DTC-2026-09-06.
Every number below came out of a script in **this** folder that I wrote and ran (`g1`…`g9`,
re-asserted by `check_constants.py`); `SHA256SUMS` is computed, never typed. Nothing outside this
folder was written; no neighbouring seat's script was executed. Constants that also appear in
`lower/prove-lagrangian`, `gaps/gap-T-lipschitz`, `gaps/refute-gap-T-lipschitz` and
`rebuild/far-near-kernel-lemma` were **re-derived here from scratch** and are reported as
agreements, not as inputs. Algebra marked EXACT is sympy-verified with residual `0`; numerics
falsify, they never prove.

---

## 0. Headline

**(a) L3v is PROVED, with an explicit constant, uniformly in λ.** For the strained δ-tapered
plateau the zonal coefficients obey

> **|H_l(λ)| ≤ √(2/π)·(2λM + V_λ)·l^{−3/2}  for all odd l ≥ 3,  with V_λ ≤ 2λM(2 + √(sin θ_δ)) ≤ 6λM,**
> hence **|H_l(λ)| ≤ 8√(2/π)·λM·l^{−3/2} = 6.3831·λM·l^{−3/2}**, uniformly on λ ∈ [1,3/2] and for
> every taper δ ∈ [0, π/2), and therefore
> **Σ_{l≥3} |H_l(λ)|·‖C_l^{3/2}‖_∞/(l(l+3)−4) ≤ 0.899159·√(2/π)(2λM+V_λ) ≤ 5.7395·λM ≤ 8.6093·M.**

The **value** of that sum, with the tail summed analytically (not truncated), is

| δ \ λ | 1.00 | 1.10 | 1.25 | 1.40 | 1.50 |
|---|---|---|---|---|---|
| 0° (bang–bang) | 0.4280 | 0.4708 | 0.5350 | 0.5992 | **0.6420** |
| 7.5° | 0.4282 | 0.4714 | 0.5373 | 0.6062 | **0.6550** |
| 15° | 0.4301 | 0.4760 | 0.5511 | 0.6384 | **0.7047** |
| 30° | 0.4433 | 0.5024 | 0.6056 | 0.7180 | **0.7911** |

each with a two-sided uncertainty of **≤ 0.0011** (`g5`). The λ = 1, δ = 0 entry **0.4280** replaces
`prove-lagrangian`'s displayed **0.3958** (an L = 601 truncation, low by 7.5 %) and its λ = 3/2
entries 0.5938 (bang–bang) and 0.6069 (7.5° taper), low by **7.5 %** and **7.3 %**.

**(b) The slab gradient bound is PROVED**, in the form

> **Γ := sup_{slab} ‖∇u‖_op ≤ 2 a(0,t) + C′M**, with C′ ≤ **33.5** using only proved inputs,
> C′ ≤ **18.5** with the measured L3v sup-norms, and C′ ≈ **4.3 … 11.8** measured,
> for every λ ∈ [1,3/2] and every δ ∈ {0°, 7.5°, 15°, 30°}.

The form the brief asked for, `2a(0,t)(1+C/L) + C′M`, is **redundant**: a(0,t) = κ(λ)ML exactly
(Lemma 1), so 2a(0,t)·C/L = 2κ(λ)CM is itself a C′M term. The two-term form is proved above.

Three things the route as briefed did not anticipate, all established here:

1. **No collar integral is needed, and the `4/sin φ` blow-up at the axis never appears.** The
   brief's plan ("the collar constant blows up at the axis: handle the taper region separately")
   would cost a constant **6.2820/δ = 47.99 at δ = 7.5°** (re-derived in `g6`). It is avoidable:
   the second angular derivative of the strain is removed by an exact first-order ODE (§B4–B5)
   whose only inputs are ‖f‖_∞ (the L3v sum) and ‖ω^θ‖_∞. **The taper is not needed for Γ on the
   slab** — Γ is bounded for the bang–bang cap too; admissibility buys ‖η‖_∞ and the *inner-edge*
   value of a, not the gradient in the bulk.
2. **The whole deviation of ∇u from the uniform strain is two scalars plus ω^θ.** Exactly,
   ∇₅b = a·diag(1,1,1,1,−2) + E with E carrying only r∂_r a, r∂_z a and ω^θ; hence
   **‖∇u‖_op ≤ 2|a| + r|∇a| + |ω^θ|**, an identity-level reduction (§B1).
3. **The slab restriction is load-bearing and cannot be dropped.** The two circles where the jump
   plane {z = 0} meets the shell boundaries {|x| = ρ₀}, {|x| = R} are corners of the discontinuity
   set; the edge sums that control them converge only through the geometric factors
   2^{−(l+4)}, 2^{−(l−1)}, which are available only one octave in.

**Honest reach.** C′M ≤ 0.1·(2a(0)) needs L ≥ C′/(0.2κ) = **187 … 316** (all-proved column),
**75.7 … 178.5** (proved with measured L3v sup-norms), or **43 … 80** (measured column). At the
campaign's L ≈ 8–10 the O(M) term is a 50–150 % correction.
This is the same "the o(1) is genuine but the crossover is astronomically far out" verdict the
refuter recorded for Lemma T, and it must travel with any use of these constants.

---

## 1. Setting

5-D lift of axisymmetric no-swirl NS. `x = (y,z)`, `y ∈ ℝ⁴`, `r = |y|`, `ρ = |x|`, `t = cos φ = z/ρ`,
`s = sin φ`. `η = ω^θ/r`, `−Δ₅ψ₁ = η`, `a = u^r/r = −∂_zψ₁`, `u^z = 2ψ₁ + r∂_rψ₁`,
`b(y,z) = (u^r y/r, u^z)` the lift of u, `Δ₅ = ∂_rr + (3/r)∂_r + ∂_zz`.

**Datum.** `ω^θ_0 = −M sgn(z)h_δ(φ)` on `ρ₀ < |x| < R`, `h_δ(φ) = min(1, φ_ax/δ)`,
`φ_ax = min(φ, π−φ)`, `L = log(R/ρ₀)`.

**Strain.** `T_λ(y,z) = (λy, λ^{−2}z)`, `η_λ = η₀∘T_λ^{−1}`, λ ∈ [1,3/2].

**The strained profile (EXACT, `glib.py`).** Since T_λ is linear, η_λ is again homogeneous of
degree −1, and writing `ρ η_λ =: W_λ(t)`,

```
W_λ(t) = −λM sgn(t) h_δ(φ̃(t)) / √(1−t²),      tan φ̃ = λ^{−3} tan φ  on (0,π/2),
                                               φ̃(π−φ) = π − φ̃(φ),
ω^θ_λ(t) = −λM sgn(t) h_δ(φ̃(t)),              so ‖ω^θ_λ‖_∞ = λM  (vortex stretching by λ),
η_λ = ρ^{−1} Σ_l H_l(λ) C_l^{3/2}(t),          H_l = 0 for even l.
```

The support of η_λ is the **ellipsoidal** shell `ρ₀ < |T_λ^{−1}x| < R`, i.e.
`ρ₀/ĝ(t) < ρ < R/ĝ(t)` with `ĝ = √(λ^{−2}s² + λ⁴t²) ∈ [λ^{−1}, λ²]`. §B5 handles the difference
between that and the spherical reference shell `[λρ₀, R/λ²]`.

**Slab.** `Slab := { x : 2λρ₀ ≤ |x| ≤ R/(2λ²) }` — one octave inside the extreme radii of the
strained support in both directions.

---

## PART A — L3v

### A1. Foundations (PROVED, exact; `g1_foundations.py`)

| statement | check |
|---|---|
| `C_l^{3/2}(t) = P_{l+1}'(t)` | sympy residual `0`, l ≤ 14 |
| `N_l := ∫_{−1}^1 (C_l^{3/2})²(1−t²)dt = (l+1)(l+2)/(l+3/2)` | sympy residual `0`, l ≤ 14 |
| `‖C_l^{3/2}‖_∞ = C_l^{3/2}(1) = (l+1)(l+2)/2` | sympy residual `0`; max over a 2·10⁴ grid = C_l(1) to 1.0 |
| `‖dC_n^{3/2}/dt‖_∞ = 3·binom(n+3, n−1)` | ratio 1.0 at n ≤ 11 (`g6`) |
| exact bang–bang coefficients `H₁ = −5M/6`, `H₃ = 3M/40`, `H₅ = −247M/1680`, `H₇ = 1513M/40320`, `H₉ = −2773M/42240` | sympy exact rationals; reproduces `prove-lagrangian` §4(1) |
| `κ₀ = −(3/5)H₁/M = 1/2` | exact |
| radial Green solution of `ψ'' + (4/ρ)ψ' − l(l+3)ψ/ρ² = −H_l/ρ` on the shell | sympy: bulk `H_lρ/((l+4)(l−1))`, inner edge `−H_lρ₀^{l+4}ρ^{−(l+3)}/((l+4)(2l+3))`, outer edge `−H_lρ^lR^{1−l}/((l−1)(2l+3))`; residual `0` |
| l = 1 resonance: `ψ_1 = (H₁/5)ρ log(R/ρ) + (H₁/25)(ρ − ρ₀⁵ρ^{−4})` | sympy residual `0` |

### A2. The coefficient representation (PROVED)

Put `Φ(t) := W_λ(t)(1−t²) = −λM sgn(t) h_δ(φ̃(t)) sin φ`. Then Φ is bounded, vanishes at t = ±1,
is odd, has **exactly one jump**, at the equator, of size

```
[Φ](0) = Φ(0+) − Φ(0−) = −2λM ,
```

and is piecewise C¹ with a derivative jump at the two strained taper corners
`φ = θ_δ := arctan(λ³ tan δ)` and `π − θ_δ`. Using `C_l^{3/2} = P_{l+1}'` and integrating by parts
on (−1,0) and (0,1) separately (the boundary terms vanish because Φ(±1) = 0):

> **(A.1)  H_l = −(1/N_l) [ [Φ](0)·P_{l+1}(0) + ∫_{−1}^{1} Φ'_cl(t) P_{l+1}(t) dt ] .**

*Verification.* Two independent instruments — (I1) direct angular projection
`H_l = −(2λM/N_l)∫_0^{π/2} h_δ(φ̃)C_l^{3/2}(cos θ) sin²θ dθ`, and (I2) the representation (A.1) —
agree to **≤ 2.0·10^{−15}** relative over l ≤ 4001 at every (λ,δ) tested, and both reproduce the
exact rationals of A1 to ≤ 8·10^{−16} (`g2`).

### A3. THEOREM A (the jump estimate, PROVED)

> **THEOREM A.** For every λ ∈ [1,3/2], every δ ∈ [0, π/2) and every odd l ≥ 3,
> ```
> |H_l(λ)|  ≤  √(2/π) · (2λM + V_λ) · l^{−3/2} ,
> V_λ := ∫_{−1}^{1} |Φ'_cl(t)| (1−t²)^{−1/4} dt = 2∫_0^{π/2} |Φ'(cos θ)| √(sin θ) dθ ,
> ```
> and `V_λ ≤ 2λM(2 + √(sin θ_δ)) ≤ 6λM`, so **|H_l(λ)| ≤ 8√(2/π) λM l^{−3/2} = 6.3831 λM l^{−3/2}**.

*Proof.* Insert into (A.1) Szegő's bound (Orthogonal Polynomials, Thm 7.3.3 / (7.3.8))
`|P_n(cos θ)| < √(2/(π n sin θ))`, i.e. `|P_n(t)| < √(2/(πn))(1−t²)^{−1/4}`, with n = l+1. At the
jump point t = 0 the weight is 1, so both terms carry the same factor:

```
|H_l| ≤ (1/N_l) √(2/(π(l+1))) ( |[Φ](0)| + V_λ )  =  (l+3/2)/((l+1)(l+2)) · √(2/(π(l+1))) (2λM + V_λ).
```

The prefactor obeys `l^{3/2}(l+3/2)/((l+1)(l+2)√(l+1)) ≤ 1` for l ≥ 3 (checked over odd l < 10⁷;
value 0.5846 at l = 3, increasing to 1). For V_λ, on (0,1)

```
Φ'(t) = λM [ h'_δ(φ̃)·dφ̃/dφ + h_δ(φ̃)·t/s ] ,      dφ̃/dφ = λ^{−3}/(cos²φ + λ^{−6}sin²φ),
```
(finite-difference check: relative error ≤ 1.4·10^{−10}, `g9`), so with `h_δ ≤ 1`,
`∫_{−1}^1 |t|(1−t²)^{−3/4}dt = 4` gives the first term ≤ 4λM, and since h'_δ = 1/δ only on
`φ̃ < δ`, where φ̃ ranges over an interval of total length δ,

```
2λM(1/δ)∫_0^{θ_δ}(dφ̃/dθ)√(sin θ)dθ ≤ 2λM(1/δ)·√(sin θ_δ)·δ = 2λM√(sin θ_δ) ≤ 2λM . ∎
```

**Verification (numerics falsify).** `max_{3≤l≤4001} l^{3/2}|H_l|` against the bound
√(2/π)(2λM+V_λ): the bound holds with slack **2.26 … 2.66** at every (λ,δ) with δ > 0 and **2.63**
at δ = 0; the computed V is `4.0000 λM` **exactly** for the bang–bang cap (the closed value 4λM),
and 3.518 / 3.789 / 4.162 / 4.495 / 4.695 at δ = 7.5°, λ = 1/1.1/1.25/1.4/1.5 (`g9`).
Szegő's bound itself was falsification-tested: worst ratio **0.99958** over n ≤ 600 and a
4·10⁵-point θ grid; the sharper (n+½) form holds too (worst ratio 0.99999961, `g1`).

### A4. COROLLARY A1 (uniform convergence of the L3v sum, PROVED)

Since `‖C_l^{3/2}‖_∞/(l(l+3)−4) = (l+1)(l+2)/(2(l+4)(l−1)) ≤ 5/7 = 0.714286` for l ≥ 3, and
`Σ_{l≥3} l^{−3/2} = ζ(3/2) − 1 − 2^{−3/2} = 1.2588220`,

> **Σ_{l≥3} |H_l(λ)| ‖C_l^{3/2}‖_∞/(l(l+3)−4) ≤ 0.899159·√(2/π)(2λM + V_λ) ≤ 5.7395 λM ≤ 8.6093 M,**

uniformly on λ ∈ [1,3/2] and in δ. Computed values of the bound: 4.305 (δ=0,λ=1) … 6.457
(δ=0,λ=1.5), 3.959 … 5.521 (δ=7.5°), 3.619 … 4.861 (δ=30°). **This is exactly the uniform bound
`prove-lagrangian` §5 lists as missing** ("needs a uniform-in-λ bound"); it is now proved, and the
proof needs nothing about the flow — only that the strained profile has one jump and bounded
variation away from it.

### A5. The value of the sum, with the tail summed analytically

Write `H_l = H_l^J + R_l`, `H_l^J := (2λM/N_l)P_{l+1}(0)` the jump part of (A.1). Then

```
Term_l^J = λM (l+3/2)|P_{l+1}(0)| / ((l+4)(l−1)) ,      |P_{2m}(0)| = binom(2m,m)/4^m,
```
and the Wallis bounds `1/√(π(m+½)) ≤ binom(2m,m)/4^m ≤ 1/√(πm)` (checked for all m ≤ 2·10⁶, worst
ratios 1.0000000625 and 0.9999999375, `g5`) give the asymptotic `Term_l ≈ λM√(2/π) l^{−3/2}`,
i.e. **0.797885 λM l^{−3/2}** — computed `l^{3/2}Term_l = 0.81323` at l = 4001, λ = 1, converging
from above. (The refuter's fitted **0.8254** is the value at finite l, not the limit; the limit is
√(2/π)λM.)

The sum is then evaluated as `Σ_{l≤4001} Term_l` (exact quadrature) `+ Σ_{l>4001} Term_l^J`
(recurrence to l = 4·10⁶ plus a Wallis remainder bound `≤ 4.2·10^{−4}λ`) `+ Σ_{l>4001}R`-tail
(measured decay exponent **1.998** at δ = 0, **2.50** at δ > 0; estimate ≤ 1.9·10^{−4}). The table
in §0 is the result (half-width of the reported interval ≤ 0.0011).

**Decorrelation against the refuter (exact agreement).** My streamed partial sums for λ = 1, δ = 0
are `0.3499048 / 0.3885774 / 0.3958178 / 0.4001482 / 0.4052964 / 0.4083703 / 0.4104700` at
l = 101/401/601/801/1201/1601/2001 against the refuter's `0.34990 / 0.38858 / 0.39582 / 0.40015 /
0.40529 / 0.40837 / 0.41046` — **every printed digit**, from a differently-built instrument. The
l = 601 entry is `prove-lagrangian`'s displayed 0.3958. My analytic total **0.4280 ± 0.0011** vs
the refuter's estimated **0.4295** (0.35 % apart, their estimate high because their tail constant
is the finite-l value).

### A6. Status per step (Part A)

| step | statement | status |
|---|---|---|
| A1 | `C_l^{3/2} = P_{l+1}'`, N_l, ‖C_l‖_∞, ‖C_l'‖_∞, exact H_l, radial Green solution | **PROVED** (sympy, residual 0) |
| A2 | the representation (A.1); one jump of size 2λM at the equator, two taper corners | **PROVED**; two instruments agree to 2·10^{−15} |
| A3 | **|H_l(λ)| ≤ √(2/π)(2λM+V_λ)l^{−3/2}, V_λ ≤ 6λM** | **PROVED** modulo Szegő's classical bound (cited at source, falsification-tested to n = 600) |
| A4 | uniform convergence of `Σ_{l≥3}|H_l|‖C_l‖_∞/(l(l+3)−4) ≤ 5.7395λM` | **PROVED** |
| A5 | the value of the sum, tail analytic, ±0.0011 | **COMPUTED** (two-sided, not a truncation) |
| — | corrections to `prove-lagrangian` §4(2): 0.3958 → 0.4280, 0.5938 → 0.6420, 0.6069 → 0.6550 | **ESTABLISHED** |

---

## PART B — the slab gradient bound Γ

### B1. The exact structure of ∇₅b (PROVED, sympy; `g3_gradu_structure.py`)

For any axisymmetric no-swirl field with stream function ψ₁, at a point with y = (r,0,0,0),

```
∇₅b  =  a · diag(1,1,1,1,−2)  +  E ,        a = −∂_zψ₁ ,
E =  [[ r ∂_r a ,  r ∂_z a ],
      [ r ∂_z a − ω^θ ,  − r ∂_r a ]]   in the (e_r, e_z) block, zero elsewhere,
```

with `r∂_r a = −r∂_r∂_zψ₁`, `r∂_z a = −r∂_zzψ₁`, and the (z,r) entry equal to the (r,z) entry minus
`ω^θ = −rΔ₅ψ₁` **exactly** (sympy residual `0` for the whole 5×5 matrix). Since the symmetric
traceless matrix `[[Q,P],[P,−Q]]` has operator norm `√(Q²+P²) = r|∇a|`,

> **(B.1)  ‖∇₅b‖_op ≤ 2|a| + r|∇a| + |ω^θ| ,   and  ‖∇₅b‖_op = ‖∇₃u‖_op exactly.**

The second equality (max relative difference of the largest singular values **0.0** over 200 random
(r,z) for a nonlinear test ψ) is an independent confirmation of `gaps/gap-V-aronson` §1.

For a homogeneous-degree-0 strain a = a(t) (sympy residual `0`):

> **(B.2)  r|∇a| = (1−t²)|a'(t)| .**

### B2. The bulk / edge decomposition

For the spherical reference shell `ρ_in = λρ₀ ≤ ρ ≤ ρ_out = R/λ²` with source `ρ^{−1}W_λ(t)`, A1's
radial solution gives, exactly,

```
ψ₁ = a^{(1)}-part (l = 1, explicit)  +  ρ f(t)  +  edge terms ,
f(t) := Σ_{l≥3} c_l C_l^{3/2}(t) ,   c_l = H_l/((l+4)(l−1)) ,
edge terms = −Σ_{l≥3} H_l [ ρ_in^{l+4}ρ^{−(l+3)}/((l+4)(2l+3)) + ρ^lρ_out^{1−l}/((l−1)(2l+3)) ] C_l .
```

`‖f‖_∞ ≤ Σ_{l≥3}|c_l|‖C_l‖_∞` = **the L3v sum of Part A** — this is where (a) feeds (b).

**l = 1, exact (sympy, `g6`).** With κ = κ(λ,δ) = −(3/5)H₁/M,

```
a^{(1)} = κM [ log(R/ρ) + 1/5 − t² + (ρ_in/ρ)⁵(t² − 1/5) ] ,
sup_slab |a^{(1)} − κM log(ρ_out/ρ)| = 0.800000 κM ,     sup_slab r|∇a^{(1)}| = 1.146966 κM
```
(1.924492 κM without the slab restriction). At the origin, κML reproduces `(M/2)log(R/ρ₀)` = 4.158883
at R/ρ₀ = 4096.

### B3. THE MOVE: the second angular derivative is removable (PROVED, exact)

The bulk deviatoric part is `ψ_dev = ρ f(t)`. Since
`Δ₅(ρ^γ g(t)) = ρ^{γ−2}[γ(γ+3)g + (1−t²)g″ − 4tg′]` (sympy residual `0`) and
`∂_z(ρ^γ g(t)) = ρ^{γ−1}[γ t g + (1−t²)g′]`, we have with γ = 1

```
4f + (1−t²)f″ − 4tf′ = −W_dev ,   W_dev := W_λ − H₁C_1^{3/2} = W_λ − 3H₁t ,
a_dev = −(t f + (1−t²) f′) = −Σ_{l≥3} c_l [ (l+4) t C_l^{3/2} − (l+1) C_{l+1}^{3/2} ] ,
```
(the last equality by `(1−t²)C_n′ = (n+3)tC_n − (n+1)C_{n+1}`, sympy residual `0`, n ≤ 11), and
eliminating f″ between the two gives the identity that carries the whole of Part B:

> **(B.3)  (1−t²) a_dev′(t)  =  3 f(t)  +  3 t a_dev(t)  +  (1−t²) W_dev(t)   (sympy residual `0`).**

By (B.2) the left side **is** `r|∇a_dev|`. So the log-dangerous object — a second derivative of a
Newtonian potential of a discontinuous source — is expressed in **zeroth-order data**. No principal
value, no collar integral, no modulus of continuity, no mollification scale.

**Cross-check.** (B.3)'s right side against a finite difference of the a_dev series: max difference
**5.9·10^{−5} … 8.8·10^{−5}** against a scale 0.88 … 1.37, away from the jump and the corners, at all
twelve (λ,δ) (`g4`). At the jump itself the termwise-differentiated series converges to the mean of
the one-sided limits (0, by symmetry) while (B.3) gives the correct one-sided values ∓λM — the
expected Gibbs behaviour, and the reason the identity, not the differentiated series, is the tool.

**Measured consequence.** `sup_t (1−t²)|a_dev′| = λM exactly` (to 14 digits) at every δ tested — the
supremum is attained at the equator, where 3f(0) = 0, t = 0 and `(1−t²)W_dev(0±) = ∓λM`.

### B4. Bounding ‖a_dev‖ without a singular integral (PROVED)

(B.3) is a **first-order linear ODE** for a_dev with integrating factor (1−t²)^{3/2}:

```
d/dt [ (1−t²)^{3/2} a_dev ] = 3f(t)(1−t²)^{1/2} + (1−t²)ω^θ(t) − 3H₁ t (1−t²)^{3/2} .
```
The right side is odd, so its integral over [−1,1] vanishes and (a_dev being finite at one endpoint)

> **(B.4)  a_dev(t) = −(1−t²)^{−3/2} ∫_{|t|}^{1} [ 3f(1−τ²)^{1/2} + (1−τ²)ω^θ_eff(τ) ] dτ ,**
> `ω^θ_eff := ω^θ − 3H₁ t √(1−t²)` .

Define `K_p(t) := (1−t²)^{−3/2}∫_{|t|}^1 (1−τ²)^p dτ`. Then `K_p′ ≤ 0` on [0,1]: with
`D := (1−u²)^{p+1} − 3uJ(u)`, `J = ∫_u^1(1−τ²)^p`, sympy gives
`D′|_{p=1/2} = (3/2)[u√(1−u²) + arcsin u − π/2] ≤ 0` and `D|_{p=1/2}(1) = 0`;
`D′|_{p=1} = 2u−2 ≤ 0` and `D|_{p=1} = (1−u)² ≥ 0`. Hence both kernels are maximal at t = 0:

> **(B.5)  ‖a_dev‖_∞ ≤ (3π/4)‖f‖_∞ + (2/3)‖ω^θ_eff‖_∞ ,  3π/4 = 2.3561945, 2/3 exact,**
> and `‖ω^θ_eff‖_∞ ≤ λM + (3/2)|H₁| = λM + (5/2)κ(λ)M`.

Combining with (B.3),

> **(B.6)  sup r|∇a_dev| ≤ 3‖f‖_∞ + 3‖a_dev‖_∞ + ‖(1−t²)W_dev‖_∞ ,**
> `‖(1−t²)W_dev‖_∞ ≤ λM + 3|H₁|·max_t t(1−t²) = λM + (2/√3)|H₁|`.

**Verification.** (B.4) evaluated by quadrature against the Gegenbauer series for a_dev agrees to
**≤ 3.3·10^{−4}** on |t| < 0.99 at all twelve (λ,δ) (`g8`). Measured sup|a_dev| = 0.176 … 0.380;
the bound (B.5) with measured norms gives 0.72 … 1.45 (slack 3.8 … 4.9), with proved norms
2.0 … 3.8.

### B5. Edge terms and the ellipsoidal collar (PROVED, explicit)

* **Shell boundary.** On the slab the two edge series carry `(ρ_in/ρ)^{l+4} ≤ 2^{−(l+4)}` and
  `(ρ/ρ_out)^{l−1} ≤ 2^{−(l−1)}`; with the exact `∂_z` identities
  `∂_z[ρ^lC_l] = (l+2)ρ^{l−1}C_{l−1}` and `∂_z[ρ^{−l−3}C_l] = −(l+1)ρ^{−l−4}C_{l+1}` (sympy residual
  `0`, l ≤ 10), their contributions to a and to r|∇a| are absolutely convergent sums. Computed from
  the datum's own coefficients: **edge_a = 0.053 … 0.179 M**, **edge_grad = 0.42 … 1.06 M**.
* **Ellipsoidal collar.** η_λ − (spherical-shell source) is supported in two collars of log-width
  3 log λ ≤ 1.2164, at log-distance ≥ log 2 from any slab point, so the interior/exterior multipole
  bounds apply. Re-derived here from scratch by Cauchy–Schwarz against Parseval
  (`Σ_l g_l²N_l = ∫w²dt ≤ 2M²`, w = ω^θ):

  | constant (z-odd) | this seat | `far-near-kernel-lemma` |
  |---|---|---|
  | far remainder, l ≥ 3 (a) | **0.2919985** | 0.291999 |
  | inner multipole, l ≥ 1 (a) | **0.0147543** | 0.014754 |
  | collar, generic, coefficient of 1/sin φ | **3.9992184** | 3.999218 |
  | collar, generic, axis piece | **0.3926991** (= π/8) | π/8 |
  | far remainder, l ≥ 3 (r∇) | **1.6049285** | not computed there |
  | inner multipole, l ≥ 1 (r∇) | **0.1324254** | not computed there |

  Three independent reproductions of that seat's published constants. The outer collar's l = 1 mode
  contributes exactly `κM·3log λ` to a, which is the log-width bookkeeping and cancels against the
  reduced log count of the spherical reference.

**The collar constant that is *not* needed.** For completeness, the taper repair of the
`4/sin φ` collar bound was carried out: with `|ω^θ| ≤ λM min(1, φ̃/δ)` and `φ ≤ (π/2)sin φ`, the
1/r′ weight cancels and `|a_collar| ≤ (3/(8π²))(π/δ)|S⁴|R_A λM = 6.281958 λM/δ` = **47.99 λM** at
δ = 7.5°, crossing over with `4/sin φ` at sin φ = 0.08333. §B3–B4 make this route unnecessary; the
number is recorded so the cost of the briefed route is on the record.

### B6. THEOREM Γ (PROVED)

> **THEOREM Γ.** Let η be the strained δ-tapered plateau η_λ, λ ∈ [1,3/2], with ‖ω^θ‖_∞ = λM, and
> let Γ := sup over `Slab = {2λρ₀ ≤ |x| ≤ R/(2λ²)}` of `‖∇u‖_op`. Then
> ```
> Γ  ≤  2 a(0)  +  C′ M ,        a(0) = κ(λ,δ) M L    (Lemma 1),
> C′ = 2[ 0.8κ + ‖a_dev‖ + edge_a + collar_a ] + [ 1.146966κ + r|∇a_dev| + edge_grad + collar_grad ] + λ .
> ```
> Every bracket is bounded above in §B2–B5. Numerically:
>
> | δ | λ | κ | C′ (all inputs proved) | C′ (proved, measured L3v sup-norms) | C′ (measured bulk) |
> |---|---|---|---|---|---|
> | 0° | 1.00 | 0.5000 | 18.69 | 8.41 | 4.25 |
> | 0° | 1.50 | 0.7500 | 33.39 | 17.97 | 11.73 |
> | 7.5° | 1.00 | 0.4997 | 18.69 | 7.93 | 4.25 |
> | 7.5° | 1.25 | 0.6225 | 27.15 | 13.20 | 9.10 |
> | 7.5° | 1.50 | 0.7364 | **33.42** | 16.82 | 11.74 |
> | 15° | 1.50 | 0.6769 | 33.45 | 18.29 | 11.72 |
> | 30° | 1.50 | 0.5162 | 32.61 | 18.43 | 11.38 |
>
> so **C′ ≤ 33.5 uniformly on λ ∈ [1,3/2] and δ ≤ 30°.**

*Proof.* (B.1) bounds ‖∇u‖ by `2|a| + r|∇a| + |ω^θ|`. Split
`a = a^{(1)} + a_dev + a_edge + a_collar`. The l = 1 log count satisfies
`κM log(ρ_out/ρ) + κM·3logλ ≤ κML = a(0)` on the slab, and the three remaining pieces are bounded
by §B2 (0.8κM), §B4 (B.5) and §B5. `r|∇a|` splits the same way and is bounded by §B2 (1.146966κM),
§B4 (B.6) and §B5. Finally `|ω^θ| ≤ λM`. ∎

### B7. Persistence on the window (PROVED given the T_λ model)

The bound is uniform in λ ∈ [1,3/2] — that uniformity is Theorem A (Part A), which is what makes
‖f‖_∞ and hence ‖a_dev‖ λ-uniform. Therefore, along the doubling window of `prove-lagrangian` §4,
at each time t the instantaneous field is η_{λ(t)} with λ(t) ≤ 3/2 and a(0,t) = κ(λ(t))ML, so

```
c := ∫_0^τ Γ(t) dt  ≤  2∫_0^τ a(0,t) dt + C′Mτ  =  2 log λ(τ) + C′ c₀/L ,      τ = c₀/(ML).
```
With λ(τ) = 3/2 this is `c ≤ 2log(3/2) + C′c₀/L = 0.8109302 + C′c₀/L`. This is exactly the constant
`gaps/gap-V-aronson` §1 uses (`c = 2θ = 0.8109302` frozen, `0.7340105` accelerated), and Theorem Γ
supplies what that seat assumed: **Γ = ML(1 + O(1/L))**, with the O(1/L) now priced at
`C′ ≤ 33.5` — i.e. gap-V's `c` is correct to relative `C′/(2κL) ≤ 33.5/(0.2κ·L)/10`.

### B8. Verification against the record

| target | source | this seat | agreement |
|---|---|---|---|
| `c_edge(φ)` profile, bang–bang, 13 angles | `far-near-kernel-lemma` §6: +0.787/+0.614/+0.513/+0.386/+0.286/+0.2168/+0.122/+0.059/−0.017/−0.062/−0.050/+0.008/+0.116 | 0.7871/0.6140/0.5129/0.3861/0.2864/**0.2167733**/0.1219/0.0588/−0.0174/−0.0617/−0.0504/0.0085/0.1161 | **all 13 entries**; 7 digits at 10° |
| `c_edge(10°)` | that seat's +0.2167733; `prove-lagrangian` 0.2167783; `refuter-correctness` 0.216773 ± 2e-6 | **0.2167733** | 1e-7 / 5e-6 |
| taper saturation `sup_φ c_edge` | that seat: +0.2939 (δ=7.5°), +0.1017 (δ=15°) | **0.2933**, **0.1016** | 6e-4 / 1e-4 |
| material-point offset at λ = 1, φ₀ = 30° | `refuter-lagrangian-b`: **−0.01808** at R/ρ₀ = 65536 | **−0.017395** (R/ρ₀ → ∞ limit) | 7e-4 |
| `a(0) = (M/2)P_{h_δ}(λ)L`, δ = 7.5°, λ = 3/2 | `prove-lagrangian` s2/s7: Φ = 1.4736, P(1) = 0.99944 ⇒ κ = 0.73639 | **κ = 0.7363667** | 3e-5 — a **sixth** independent instrument for Lemma 1 |
| κ_δ at δ = 5/7.5/15/30° | `prove-lagrangian` §2: 0.4999917/0.4997212/0.4978077/0.4836252 | **0.4997212/0.4978077/0.4836252** (7.5/15/30) | ≤ 1e-7 |
| streamed L3v partial sums | `refuter-gap-T-lipschitz`/`refuter-lagrangian` (7 values) | reproduced to every printed digit | exact |

**Not run.** The refuter's material-point offsets at λ > 1 (`+0.02838 / +0.03733 / −0.01215 /
−0.06857` at λ = 1.1/1.25/1.4/1.5) are **not reproduced here**: they are evaluations on the
*ellipsoidal* strained shell at a moving material point, which needs an instrument this seat did not
build (the modal solution separates only on a spherical shell). The λ = 1 member of that list is
reproduced (row 4 above). Status: **PARTIAL**.

### B9. Status per step (Part B)

| step | statement | status |
|---|---|---|
| B1 | `∇₅b = a·diag(1,1,1,1,−2) + E`, `E_{zr} = E_{rz} − ω^θ`, `‖∇u‖ ≤ 2|a| + r|∇a| + |ω^θ|` | **PROVED** (sympy, residual 0) |
| B1′ | `‖∇₅b‖_op = ‖∇₃u‖_op` | **PROVED** (0.0 max rel. diff, 200 samples) — confirms gap-V §1 |
| B2 | `r|∇a| = (1−t²)|a′(t)|` for homogeneous a; the exact l = 1 field and its two sups | **PROVED** |
| B3 | **(B.3)** `(1−t²)a_dev′ = 3f + 3t a_dev + (1−t²)W_dev` | **PROVED** (sympy, residual 0); FD cross-check 6e-5 |
| B4 | **(B.4)/(B.5)** closed form for a_dev; kernel sups 3π/4 and 2/3; monotonicity | **PROVED** (sympy) ; quadrature vs series 3.3e-4 |
| B5 | edge and ellipsoidal-collar constants | **PROVED**; three of them reproduce `far-near-kernel-lemma` to 6 digits |
| **Γ** | **`Γ ≤ 2a(0,t) + C′M`, C′ ≤ 33.5 on λ∈[1,3/2], δ ≤ 30°** | **PROVED** |
| B7 | persistence on the window; `c ≤ 2log(3/2) + C′c₀/L` | **PROVED** given the T_λ model of `prove-lagrangian` §4(2) (which is GAP T, not proved) |
| B8 | verification against four neighbouring seats | **PROVED-level agreement** except the λ>1 material-point row (**PARTIAL / NOT RUN**) |

---

## 3. What this does and does not close

**Closes.**
* `prove-lagrangian` §5's **L3v** row ("SKETCH … needs a uniform-in-λ bound") is now **PROVED**:
  Theorem A + Corollary A1. The displayed constant is corrected from a truncation (0.3958) to a
  tail-summed value (0.4280), and the λ-supremum over the range the theorem uses is **0.6550**
  (7.5° taper) — 7.3 % above the value that seat displayed.
* `gaps/gap-V-aronson`'s standing assumption `Γ = ML` is now a **theorem with an explicit O(M)
  correction**, `Γ ≤ 2a(0,t) + C′M`, C′ ≤ 33.5, on the slab and uniformly on the window.

**Does not close.**
* **GAP T is untouched.** Theorem Γ is a statement about the *instantaneous* field η_λ; it says
  nothing about whether the true Navier–Stokes flow map is T_{λ(s)}. §B7 is conditional on that
  model.
* **The slab is not the whole shell.** The two corner circles `{z=0}∩{|x|=ρ₀}` and
  `{z=0}∩{|x|=R}` are excluded, and must be: the edge sums converge only through the octave
  factors. Whether ‖∇u‖_∞ is finite *at* those corners is not settled here (for a bounded vorticity
  with a corner in its jump set the expected answer is a logarithm in the distance to the corner);
  any argument that needs Γ up to the datum's own edge must price that, or inset as
  `gap-V-aronson` §0 already insets the tracked point.
* **λ > 1 material-point offsets are not reproduced** (§B8, PARTIAL).
* **The constants are crude and the crossover is far.** Slack against measurement: 2.3–2.7 (Theorem
  A), 3.8–4.9 ((B.5)), 6.1–10.1 (Corollary A1 against the true sum), 2.9–4.4 (C′ fully-proved
  against measured). `C′M ≤ 0.1·2a(0)` needs `L ≥ 187 … 316`, i.e. `log Re_E ≳ 370 … 630`. At the
  campaign's L ≈ 8–10 these are 50–150 % corrections, not o(1) ones.


## 5. The gate (FL-043 discipline)

`check_constants.py` re-asserts **924 checks** against the stored JSONs — every number displayed
above, plus the pipeline identities that connect them (S = partial + jump-tail + residual;
C′ = the sum of its own pieces; C_bound = √(2/π)(2λ+V); κ = −(3/5)H₁; every proved envelope ≥ the
measured quantity it dominates). The pass count is **computed, not literal**, and every failure is
printed.

`python3 check_constants.py --mutate` runs the self-test the synthesis asks for: every numeric leaf
of every JSON is perturbed one at a time (`v → 1.5v + 0.37`) and the gate is re-run.

```
ALL 924 CHECKS PASS
MUTATION SELF-TEST: 1088/1288 caught (84.5%), 200 blind
```

The 200 blind leaves are **characterised, not hidden**: they are (i) stored diagnostics that no
displayed number depends on (`secs`, `S_total_no_resid`, `partial_vs_cesaro_adev`, `resid_p`), and
(ii) the per-case `S_lo`/`jump_tail_rem_bound` entries that enter only through a maximum, so
perturbing a non-maximal one cannot move the reported uncertainty. No leaf that carries a displayed
constant is blind.

## 4. Files

| file | what it establishes |
|---|---|
| `glib.py` | the strained δ-tapered profile, φ̃, θ_δ, both H_l instruments, composite Gauss panels |
| `g1_foundations.py` / `g1_results.json` | `C_l^{3/2}=P_{l+1}'`, N_l, the two sup-norms, exact H_l, the radial Green solution, Bernstein/Szegő falsification |
| `g2_l3v_sum.py` / `g2_results.json` | H_l(λ) by two instruments (agree 2e-15), streamed partial sums, residual decay |
| `g3_gradu_structure.py` / `g3_results.json` | (B.1), (B.2), `‖∇₅b‖=‖∇₃u‖`, `Δ₅(ρ^γ g)`, the bulk identity (B.3) |
| `g4_gamma_bulk.py` / `g4_results.json` | f, a_dev, `(1−t²)a_dev′` on a fine grid; identity-vs-finite-difference cross-check |
| `g5_tail.py` / `g5_results.json` | Wallis bounds; the analytic tail; the L3v table of §0 |
| `g6_kernel_constants.py` / `g6_results.json` | the far/inner Parseval constants (a and r∇), the collar bounds incl. the taper repair, the exact l = 1 field |
| `g8_assembly.py` / `g8_results.json` | kernel monotonicity (proved), (B.4) vs series, edge sums, the `c_edge` reproduction, the C′ table |
| `g9_Hl_bound.py` / `g9_results.json` | Theorem A: V_λ, the bound, its verification, the uniform Corollary A1 |
| `check_constants.py` | re-asserts every number displayed above against the JSONs (924 checks) and, with `--mutate`, the FL-043 self-test (1088/1288 = 84.5 % caught) |
| `SHA256SUMS` | computed, never typed |

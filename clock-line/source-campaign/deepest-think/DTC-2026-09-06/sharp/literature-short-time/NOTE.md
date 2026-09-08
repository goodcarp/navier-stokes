# literature-short-time — DTC-2026-09-06 / sharp

Scope: literature at source for the short-time / doubling-time question. All PDFs
downloaded from arXiv into `pdf/`, text-extracted with `pdftotext -layout` into `txt/`.
Page numbers below are **PDF pages of the arXiv file in `pdf/`**, located by
`scripts/pagefind.sh`. Quotes are transcribed from the pdftotext output (math
re-typeset into ASCII; wording verbatim).

## Bottom line (four answers)

1. **Kim–Jeong 2111.14078 is Euler only.** Its "Key Lemma" (Lemma 2.1, p.4) is a
   *static* Biot–Savart identity, not an evolution lemma. Its short-time window is
   `A T <= c_0` with `A = sup_t ||omega(t)||_inf` — i.e. **`T ~ c/M`, no logarithm**.
   The logarithm appears only through the *number of rings*, not through Re.
2. **Bang–Cheskidov 2605.16502 is Euler only** — no viscous statement, no viscous
   remark, and no Navier–Stokes reference in its bibliography. Its inflation theorem
   is qualitative in time ("for any delta > 0"), for every Lorentz exponent q > 1.
3. **No published lower bound shows a logarithm is necessary** in a fixed-factor
   vorticity-doubling clock, for NS *or* for Euler. The only rigorous NS
   vorticity-growth-ratio result is Luo 2504.08288 Thm 1.3 (Apr 2025), and it lives at
   `M_0 t_* -> infinity` (long clock), so it says nothing about `c/M_0` scales.
4. **Yes — such a statement is published: Bradshaw–Farhat–Grujić Thm 10** gives
   `||omega(t)||_inf <= M ||omega_0||_inf` on `[0, 1/(c(M)||omega_0||_inf)]` with the
   `L^2` hypothesis explicitly declared non-quantitative. **If correct as stated it
   contradicts the necessity of the `log Re` factor.** I located the exact step in
   their (self-described) *sketch* where the energy/log enters and is dropped — see
   §4 below. Flag: my reading of a proof sketch, not a published erratum.

---

## 1. Kim–Jeong, arXiv:2111.14078v2 = J. Funct. Anal. 283 (2022) 109673

Title (arXiv abs page): "A simple proof of ill-posedness for incompressible Euler
equations in critical Sobolev spaces", Junha Kim and In-Jee Jeong; v1 28 Nov 2021,
v2 17 Jul 2022. Published title (as cited by Bang–Cheskidov, [KJ22], p.98):
"A simple ill-posedness proof for incompressible Euler equations in critical Sobolev
spaces, J. Funct. Anal. 283:109673, 2022."

**Equation studied** (p.2, eq. (1.4)): the axisymmetric-without-swirl reduction of the
**incompressible Euler** equations in R^d, `d >= 3`,
`∂_t ω + (u^r ∂_r + u^d ∂_d) ω = (d-2) (u^r / r) ω`. Viscosity appears nowhere. The
only Navier–Stokes mentions in the entire paper are in the introduction's citation list
(p.1, "enstrophy growth in the Navier–Stokes equations ([14, 13])") and refs [13],[14],
[15] — Jeong–Yoneda enhanced-dissipation papers and Kato. **Euler only.**

### 1a. The "outer-region" Key Lemma — verbatim (p.4, §2.2)

> **Lemma 2.1.** We impose the following assumptions on `ω ∈ H^{d/2} ∪ (L^∞ ∩ L^2)(R^d)`:
> - Odd with respect to the last coordinate: `ω(y_h, y_d) = -ω(y_h, -y_d)`;   (2.3)
> - For any `y ≠ 0` satisfying either `|y_h| = 0` or `y_d = 0`, there exists an open
>   neighborhood of `y` such that `ω` vanishes.
>
> Then, there exists a constant `C > 0` such that for any `x ∈ R^d` with
> `r = |x_h| >= x_d > 0`, we have
>
>     | u^r(x)/r  -  (1/((d-1)|B_d|)) ∫_{Q(x)} (|y_h| y_d / |y|^{d+2}) ω(y) dy |  <=  C B_1(x)     (2.4)
>
> and
>
>     | u^d(x)/x_d  +  (1/|B_d|) ∫_{Q(x)} (|y_h| y_d / |y|^{d+2}) ω(y) dy |  <=  C B_2(x),        (2.5)
>
> where `Q(x) := { y ∈ R^d ; |y_h| >= 4|x_h| }` and `B_1, B_2` are non-negative functions
> with the upper bound
>
>     B_1(x) <= min( ||∇ω||_{L^d(R^d)} , ||ω||_{L^∞(R^d)} ),
>     B_2(x) <= min{ (1 + log(r/x_d))^{(d-1)/d} ||∇ω||_{L^d} , (1 + log(r/x_d)) ||ω||_{L^∞} }.

Notes: `Q(x)` — the region `|y_h| >= 4|x_h|` — is the "outer region". The lemma is
**static in `t`** (no time variable, no time scale); it is applied to `ω(t,·)` for each
fixed `t`. The remainder for `u^r/r` is `O(||ω||_∞)` with **no log** — this is exactly
the `a = (M/2) log(R/rho_0) + O(M)` split in the brief, with the log carried entirely by
the outer integral. They credit the form to Kiselev–Šverák (2D) and Bourgain–Li, and
say (p.4): "Our version of Key Lemma is sharp in that the remainder is bounded using
only critical norms of the Euler equations."

### 1b. The short-time window (§3 "Short time dynamics", pp.6–9)

Standing hypothesis (p.6, (3.1)–(3.2)): a solution on `[0,T]` with
`sup_{t∈[0,T]} min{ ||Λ^{d/2} ω(t)||_{L^2}, ||ω(t)||_{L^∞} } <= A`, and then verbatim:

> "We shall let `T` be small so that `A T <= c_0` for some absolute constant `c_0 > 0`,
> which is determined at the end of this section."  (p.6)

**So the master time scale is `T <= c_0 / A`, `A` = sup of the vorticity max — no
logarithm anywhere.** The scale-dependent refinement is:

> **Lemma 3.3** (p.8). Assume that `ω` be a solution to (1.4) satisfying (3.2) with initial
> data (2.1). Then, there exists a sequence `{T_n}_{n>=n_0}` with
>
>     T_n := min{ T , c_1 (1-α) n^{-1+α} }
>
> for some absolute constant `c_1 > 0` such that
> `r/2 <= Φ^r(t,x) <= 2r` and `x_d/2 <= Φ^d(t,x) <= 2 x_d`   (3.8)
> for all `x ∈ Ω_n` with `n >= n_0` and `t ∈ [0, T_n]`.

with data (2.1) = the dyadic ring stack `ω_0 = Σ_{k=n_0}^{m} k^{-α} φ(r/(1/8)^{k-1}, z/(1/8)^k)`
and `I_n(t) ≃ I_n(0) = C n^{-α}` on `[0,T_n]`, `Σ_{k=√n}^{n-1} ∫_0^{T_k} I_k dt >= log n^{c_2}` (3.10, p.9, `c_2 < 1/4`).

**Quantitative payload — Proposition 4.1 (p.10), verbatim:**

> **Proposition 4.1.** We consider the sequence of `C^∞`–smooth initial data `ω_0^{(m)}` given
> by (2.1), where `n_0 < m < ∞` and `0 < α < c_2`, with `c_2 > 0` being the constant from (3.10).
> Given any `ε > 0` and `q > 1/α`, by taking `n_0` large, we have uniform bounds
> `||ω_0^{(m)}||_{L^1 ∩ L^∞(R^3)} <= ε`, `||r^{-1} ω_0^{(m)}||_{L^{3,q}(R^3)} <= ε` for all `m >= n_0`.
> Then, there exists `M = M(n_0) > 0` such that for all `m > M`, the unique global in time
> solution `ω^{(m)}` with initial data `ω_0^{(m)}` satisfies
>
>     sup_{t ∈ [0, T(m)]} ||ω^{(m)}(t,·)||_{L^∞(R^3)}  >=  (1/4) m^{c_2 - α}
>
> for `T(m) := c_1 (1-α) m^{(1/2)(-1+α)}`, where `c_1` is from (3.9).

(Prop 6.1, p.15, is the `Ḣ^{d/2}`-analogue: `sup_{[0,T(m)]} ||Λ^{d/2}ω||_{L^2}^2 >= (c_3/2) Σ_{n=m/2}^{m} n^{c_4-2α}`
at `T(m) = c_1(1-α)(m/2)^{(1/2)(-1+α)}`.)

### 1c. What Prop 4.1 does and does not say about a log clock (my arithmetic, `scripts/kj_scale_audit.py`)

Bookkeeping only on their own formulas; no new PDE claim. With `N ≈ m` rings and
`Λ := log(R_outer/R_inner) = (m-n_0) log 8`, `M_0 = ||ω_0||_∞ = n_0^{-α}`:

- **The stretching rate at ring `k` is the log-enhanced one.** `a_k = Σ_{n<k} n^{-α}` and
  `a_k / (M_0 Λ)` = 0.472, 0.453, 0.433, 0.414 for `α=0.02`, `k=10,10^2,10^3,10^4`
  (`n_0=2`) — i.e. `≈ 1/log 8 = 0.481`, matching the brief's `κ_0 = 1/2` per e-fold.
- **`T_k` *is* the reciprocal of that rate.** `a_k T_k / c_1` = 0.795, 0.979, 0.998, 0.9998
  over the same `k`. So Lemma 3.3's window is exactly `T_k ≈ c_1/(M_0 Λ)`: **KJ's own
  scale-dependent clock has the `1/(M log)` form.** But it is the window on which the
  *weak inner ring* (amplitude `k^{-α}`) is trackable, not a doubling time for `||ω||_∞`.
- **Prop 4.1's sup-norm clock is `Λ^{-(1-α)/2}`, not `Λ^{-1}`.** Fitting
  `T(m) = c/(M_0 Λ^β)` over `m ∈ [10^6, 10^12]` gives `β` = 0.4900, 0.4750, 0.4500, 0.3800
  for `α` = 0.02, 0.05, 0.10, 0.24. Since Prop 4.1 needs `α < c_2 < 1/4`,
  **`β ∈ (3/8, 1/2)` — at most half a power of the log.** The upper bound's clock has `β = 1`.
- **At a *fixed* growth factor KJ gives no `Λ`-dependence at all.** Growth is `m^{c_2}`, so a
  factor-2 gain needs only `m = 2^{1/c_2}` = 32 / 18 / 16.9 rings for `c_2` = 0.20 / 0.24 / 0.245.
  Equivalently `M_0 T ≲ G^{-(1-α)/(2 c_2)}` with `(1-α)/(2c_2) > 2`: a **power law in the growth
  factor**, not a log in Re.
- **KJ Prop 4.1 does refute the log-free clock for Euler**: `M_0 T(m)/c_1` = 3.4e-2, 1.3e-3,
  4.9e-5, 1.8e-6 for `m = 10^3 … 10^12` (`α=0.05`, `c_2=0.20`) while the growth factor
  `A(m)/M_0` = 0.73, 2.06, 5.79, 16.33. So `sup ||ω||` exceeds `(3/2) M_0` at times
  `≪ c/M_0`. At infinite Re, *some* correction is necessary; KJ does not say it is one log.

---

## 2. Bang–Cheskidov, arXiv:2605.16502v1 [math.AP] 15 May 2026

"Sharp Ill-Posedness of the Euler Equations in Lorentz Spaces", Jeaheang Bang and
Alexey Cheskidov (Westlake University), 98 pp.

**Euler only.** Equations (p.2, (1)–(2)) are `∂_t u + u·∇u + ∇p = 0`, `div u = 0`. Greps
for `navier|viscos|viscous` over the full text hit **only three bibliography lines**
(Bedrossian–Vicol GSM 225, Kato 1972, Tsai GSM 192) — **no viscous theorem, no viscous
remark, no Navier–Stokes discussion anywhere in the body.**

**Theorem 1.1 (Norm inflation), p.8, verbatim:**

> Fix any `q > 1`. For any `ε̃, δ, A > 0`, there exists an axisymmetric initial datum
> `ω_0 ∈ C_c^∞(R^3)` such that
> `||ω_0||_{L^∞ ∩ L^1(R^3)} + ||ω_0/r||_{L^{3,q}(R^3)} <= ε̃`
> and the unique global-in-time solution `ω(x,t)` to the vorticity equation (4) with the
> initial condition `ω|_{t=0} = ω_0` satisfies
> `sup_{0 <= t <= δ} ||ω(·,t)||_{L^∞(R^3)} >= A`.

**Theorem 1.2 (Instantaneous blow-up), p.11, verbatim:**

> For any `ε̃ > 0` and `q > 1`, there exists an axisymmetric initial datum
> `ω_0 ∈ C^∞(R^3 \ {0})` with compact support such that
> `||ω_0||_{L^∞ ∩ L^1(R^3)} + ||ω_0/r||_{L^{3,q}(R^3)} <= ε̃`, and the following properties hold:
> (1) For any `T > 0`, there does not exist a Yudovich-type weak solution `(u,ω)` on `(0,T)` of
> (18) with the initial data `ω_0` in the sense of Definition 1.1;
> (2) There exists a distributional solution `(u,ω)` in `(0,∞)` of (18) with the same initial
> data `ω_0` in the sense of Definition 1.1 such that
> `ess sup_{0 < t < T} ||ω(t)||_{L^∞(R^3)} = +∞` for any `T ∈ (0,∞)`.   (20)

**Answering the brief's question literally: `q` = every Lorentz exponent `q > 1`; time =
"any `δ`" (Thm 1.1) / "any `T`" (Thm 1.2); growth = from `||ω_0||_∞ <= ε̃` to `A`
arbitrarily large (Thm 1.1) or `+∞` (Thm 1.2). There is no quantified time-vs-amplitude
relation and no viscous statement.** Theorem 1.2 uses infinitely many rings, i.e.
`R/rho_0 = ∞`; the viscous cut `rho_0 >= sqrt(ν/M)` truncates that stack, so Thm 1.2 has
no NS analogue as stated.

**Their description of the KJ mechanism (§1.5, p.7, verbatim in part):**

> "Kim–Jeong's key lemma turns the outer-region mechanism in the Biot–Savart law into an
> estimate usable for short times after `t = 0`: it writes the stretching rate `u^r(x,t)/r`
> as a sign-definite outer-region integral up to controllable errors. … Kim–Jeong overcame
> this difficulty by combining the key lemma with a scale-dependent stability statement for
> the outer rings: on a scale-dependent time interval, each outer ring remains localized in
> a comparable annulus, so that the contributions by the outer rings stay comparable to its
> initial value."

And their own mechanism — the answer to the brief's "does something prevent it?" question,
in the *inviscid* setting (abstract, p.1):

> "vortex stretching weakens its own future forcing: as a ring amplifies, incompressibility
> flattens it, the aspect ratio collapses, and the induced stretching coefficient is
> geometrically depleted."

The explicit inviscid depletion scaling they derive (p.6): `u^r/r ≲ ||ω||_∞ H/R` — the
stretching rate is cut by the aspect ratio `H/R` as the ring flattens.

---

## 3. Published lower bounds for NS with bounded vorticity

**Searched:** "vorticity doubling time"; "lifespan lower bound bounded vorticity
Navier-Stokes logarithm"; "norm inflation Navier-Stokes bounded vorticity short time";
Kim–Jeong follow-ups 2022–2026; Choi–Jeong; citing work on BFG.

**Result: nothing establishes that a logarithm is necessary.** The relevant items found:

- **Luo, arXiv:2504.08288, "Sharp norm inflation for 3D Navier-Stokes equations in
  supercritical spaces" (11 Apr 2025).** Its own background statement (p.3): *"no rigorous
  lower bounds for (1.5) [`|ω(t)|_∞/|ω_0|_∞`] have been established for solutions of the
  Navier-Stokes equations to date. In other words, there remains no proof of growth—let
  alone sustained growth or blowup—for solutions of (1.1)."* Then **Theorem 1.3 (p.3),
  verbatim:** *"For any `M > 0`, there exist a time `t_* > 0` and a smooth solution `u` of
  (1.1) on `[0,t_*]` such that its vorticity `ω = ∇ × u` satisfies `|ω(t_*)|_{L^∞}/|ω_0|_{L^∞} >= M`."*
  This is a genuine NS statement, but at the **wrong clock**: with their parameters
  (p.9, (2.10)–(2.11): `ν = μ^{1-b}`, `t_* = ε^{-N-2} μ^{-1-2/p+s} ν^{-1/p}`), and the
  bounds (6.11)/(6.12) at `p=1`, one gets `|ω_0|_∞ t_* ≈ ε^{-N}` and
  `|ω(t_*)|_∞ t_* ≈ ε^{-N-4}`, both `→ ∞`. So the growth happens **long after** `1/M_0`,
  and Thm 1.3 constrains neither a `c/M_0` nor a `c/(M_0 log Re)` clock. Mechanism is
  mixing/un-mixing (forward and backward cascade), **not** an axisymmetric ring stack.
- Luo's Theorem 1.1 covers `Ḃ^s_{p,q}` with `-3 < s - 3/p < -1` (supercritical side of the
  critical line), `s ≠ 0`. Bounded vorticity (`u ∈ Ḃ^1_{∞,∞}`-type) is far on the
  *subcritical* side and is **not** covered.
- **Lim–Jeong, arXiv:2409.19497v2, ARMA 249:32 (2025), Thm 1.1 (p.1):** for **Euler**
  axisymmetric no-swirl with `ω_0^θ` compactly supported and `||r^{-1} ω_0^θ||_∞ < ∞`,
  `||ω(t,·)||_{L^∞} <= A(ω_0^θ)(1+|t|)^{4/3}` for all `t ∈ R` (Childress's conjectured rate;
  a `(1+|t|)^{3/2}` variant appears at p.6). This is a **long-time upper** bound and uses
  energy + conserved quantities — the opposite direction and the wrong time regime.
- **Egamberganov–Yao, arXiv:2512.13456 (15 Dec 2025)**, Euler, axisymmetric no-swirl:
  radial moment `≳ t/log t`, `||ω(·,t)||_{L^p} ≳ t^{1/4}` in limsup — again long-time,
  Euler, and the *first* power-law `L^p` growth for smooth compactly supported data.
- **Cao–Fan–Qin, arXiv:2511.03171v4 (3 Aug 2026)**, Euler, anti-parallel class: linear
  growth of the vorticity maximum on a large proportion of every dyadic interval. Long-time.
- **Grujić, arXiv:2607.08866 (2026), "Logarithmic Depletion of Vortex Stretching and
  Singularity Evasion in the 3D Navier-Stokes Equations"** — unrefereed preprint by the
  same author as [BFG]; a *conditional* regularity result (vorticity direction in a
  log-weighted BMO class), not a lifespan lower bound. Flagged, not used.

**Also checked and negative:** Jeong–Yoneda (Math. Ann. 380 (2021) 2041–2072;
arXiv:2001.02333, 2012.14621) is enhanced/anomalous dissipation, not lifespan. No
Kim–Jeong viscous sequel exists as of this search. Bang–Cheskidov's acknowledgements
(p.97) thank In-Jee Jeong for posing the problem "beyond Danchin's regime" — still Euler.

---

## 4. Bradshaw–Farhat–Grujić, arXiv:1704.05546v4 (7 Sep 2018) = ARMA (2019), DOI 10.1007/s00205-018-1314-5

Setting: `∂_t u + (u·∇)u = ν Δu - ∇p + f` (p.1, (1)); the vorticity work is at `ν = 1`
(p.8, (5): `∂_t ω_j - Δ ω_j + u_i ∂_i ω_j = ω_i ∂_i u_j`). Theorem 8/10 are stated
in §2 as an ingredient, not as the paper's main result.

**Theorem 8 (real setting), p.8, verbatim:**

> Let the initial datum `ω_0` be in `L^2 ∩ L^∞`. Then there exists a unique mild solution
> `ω` in `C_w([0,T], L^∞)` where `T >= (1/c) (1/||ω_0||_∞)` for an absolute constant `c > 0`.

**Theorem 10 (complex setting), p.11, verbatim:**

> Let the initial datum `ω_0` be in `L^2 ∩ L^∞`, and `M` a constant larger than 1. Then
> there is a constant `c(M) > 1` such that there exists a unique mild solution `ω` in
> `C_w([0,T], L^∞)` where `T >= (1/c(M)) (1/||ω_0||_∞)`, and for any `t` in `(0,T]` the
> solution `ω` is the `R^3`-restriction of a holomorphic function `ω` defined in the domain
> `Ω_t = { x + iy ∈ C^3 : |y| < (1/sqrt(c(M))) sqrt(t) }`; moreover,
> `||ω(t)||_{L^∞(Ω_t)} <= M ||ω_0||_∞`.

**Their explicit no-energy-dependence claim, p.8, verbatim:**

> "In addition to the initial vorticity `ω_0` being bounded, a suitable decay of `ω_0` at
> infinity will be required (we chose `ω_0` in `L^2` for convenience); however, it is worth
> noting that this is a 'soft assumption', i.e., there will be no quantitative dependence
> on `||ω_0||_2` in the proof."

### Answers to question (4)

- **Existence time vs growth bound.** Theorem 8 is *existence time only*: `T >= 1/(c||ω_0||_∞)`.
  **Theorem 10 is the one that also carries the growth bound**: `||ω(t)||_{L^∞(Ω_t)} <= M ||ω_0||_∞`
  on `[0, 1/(c(M)||ω_0||_∞)]` — sup over the complex strip, so a fortiori over `R^3`.
- **So yes: a published statement gives `||ω(t)||_∞ <= C ||ω_0||_∞` on `[0, c/||ω_0||_∞]`
  with no energy dependence.** It is BFG Theorem 10, ARMA 2019, p.11 of arXiv v4.
  Under NS scaling `ω ↦ λ^2 ω(λx, λ^2 t)` both sides are invariant, so the statement is
  `ν`-uniform and **carries no `log Re` whatsoever.** Taken at face value it makes the
  brief's `H_0 = c/(M_0(1 + log_+ Re_0))` **not sharp**, and answers the sharpness question
  from the *upper*-bound side rather than the lower-bound side.

### The step where I think the energy/log was dropped — FLAG, my reading of their sketch

They label §2's argument "Sketch of the proof" (p.8) and "we present a sketch here …
for completeness. What follows is a modification of the exposition given in [Ku2]"
(Kukavica, JDE 194 (2003) 39–50). The vortex-stretching estimate (pp.9–10) runs:

> "The main ingredients are: a pointwise estimate on the heat kernel `G`,
> `G(x,t) <= c sqrt(t)/(|x| + sqrt(t))^4` …, a property of a scalar-valued BMO function `f`
> **featuring a suitable decay at infinity** (e.g., being in the closure of the test functions
> in the uniformly-local `L^p` for some `p`, `1 <= p < ∞`), `∫ |f(x)|/(|x|+1)^4 dx <= c ||f||_{BMO}`
> (**in general, one has to subtract a local average**, for example, over a unit ball `B`, in
> which case the inequality takes the form `∫ |f(x) - (1/|B|)∫_B f|/(|x|+1)^4 dx <= c ||f||_{BMO}` [St])"

leading to (10), p.10: `∫_0^t∫ G(x-y,t-s) |ω| |∂u| <= c t (sup_{s} ||ω^{(n)}(s)||_∞)^2`,
hence `sup ||ω^{(n+1)}||_∞ <= c||ω_0||_∞ + c t (sup||ω^{(n)}||_∞)^2` and `T ≈ 1/(c||ω_0||_∞)`.

Two things about that chain:

1. `∫ |f|/(|x|+1)^4 dx <= c ||f||_{BMO}` is **false with an absolute constant** — as they
   themselves note, the correct inequality subtracts a local average. The parenthetical
   decay hypothesis ("closure of the test functions in the uniformly-local `L^p`") is
   *qualitative*: it makes the integral finite but does **not** supply a constant depending
   only on `||f||_{BMO}`. The honest form is
   `∫ |f|/(|x|+1)^4 dx <= c ( ||f||_{BMO} + |f_{B_1}| )`.
2. After the heat-kernel rescaling `y = x + sqrt(t-s) z` the "unit ball" is
   `B(x, sqrt(t-s))`, so the dropped term is `|avg_{B(x, sqrt(t-s))} ∇u|`. For
   `ω ∈ L^∞ ∩ L^2`, this average is **not** `≲ ||ω||_∞`; the standard Biot–Savart splitting
   gives `|avg_{B(x,ρ)} ∇u| ≲ ||ω||_∞ log(L/ρ) + (large-scale term controlled by the energy)`,
   with `L` the outer/energy scale. With `ρ ~ sqrt(t-s) ~ sqrt(1/M)` and `L` set by `E`,
   `log(L/ρ)` is exactly a `log Re`. **The dropped term is the log.** It is also exactly the
   term whose absence they advertise in the "no quantitative dependence on `||ω_0||_2`" remark.

I found **no published erratum and no citing paper** contesting Theorem 8/10 (searches on
the ARMA article, Grujić's page, and citation listings returned nothing). So the register-honest
statement is: **the log-free clock with a constant-factor growth bound is in the published
record (BFG Thm 10, ARMA 2019), and its only available proof is a sketch with a step I can
identify as the place the `log Re` would have to enter.** Resolving this is a decision point
for the campaign, not something this pass settles. The two directly relevant precedents they
cite for the `p = ∞` case — Giga–Inui–Matsui (Hokkaido preprint, 1998) and Kukavica (JDE 194
(2003) 39–50) — are **velocity-formulation** `L^∞` results (`T ~ 1/||u_0||_∞^2`), a different
scaling, and do not by themselves yield the vorticity statement.

---

## Table

| # | Statement | Hypotheses | Time scale | Euler / NS | Citation with page |
|---|---|---|---|---|---|
| 1 | **Key Lemma**: `u^r(x)/r = (1/((d-1)|B_d|)) ∫_{Q(x)} (\|y_h\|y_d/\|y\|^{d+2}) ω dy + O(B_1)`, `Q(x)={\|y_h\|>=4\|x_h\|}`, `B_1 <= min(\|\|∇ω\|\|_{L^d}, \|\|ω\|\|_∞)`; same for `u^d/x_d` with `B_2 <= min((1+log(r/x_d))^{(d-1)/d}\|\|∇ω\|\|_{L^d}, (1+log(r/x_d))\|\|ω\|\|_∞)` | `ω ∈ H^{d/2} ∪ (L^∞∩L^2)(R^d)`, odd in `y_d`, vanishing near the axis and the `y_d=0` plane; `r=\|x_h\| >= x_d > 0` | **none** (static Biot–Savart identity; no `t`) | **Euler** | Kim–Jeong, arXiv:2111.14078v2, **Lemma 2.1, p.4** (§2.2); = JFA 283 (2022) 109673 |
| 2 | Standing short-time hypothesis: "We shall let `T` be small so that `A T <= c_0` for some absolute constant `c_0`" | `A = sup_{[0,T]} min{\|\|Λ^{d/2}ω\|\|_2, \|\|ω\|\|_∞}` | **`T <= c_0/A`, no log** | **Euler** | Kim–Jeong, **p.6** (§3, after (3.2)) |
| 3 | **Lemma 3.3** (ring stability): `∃ T_n := min{T, c_1(1-α)n^{-1+α}}` with `r/2 <= Φ^r(t,x) <= 2r`, `x_d/2 <= Φ^d <= 2x_d` for `x ∈ Ω_n`, `t ∈ [0,T_n]` | dyadic ring data (2.1); `1/q < α < 1`; solution satisfying (3.2) | **`T_n = c_1(1-α) n^{α-1}`** — scale-dependent; numerically `a_n T_n → c_1`, i.e. `T_n ≈ c_1/(M_0 log(R/rho_0))` | **Euler** | Kim–Jeong, **Lemma 3.3, p.8**; ratio check `scripts/kj_scale_audit.py` |
| 4 | **Prop 4.1**: `sup_{[0,T(m)]} \|\|ω^{(m)}\|\|_∞ >= (1/4) m^{c_2-α}` from `\|\|ω_0\|\|_{L^1∩L^∞} <= ε` | `n_0 < m < ∞`, `0 < α < c_2 < 1/4`, `q > 1/α`, `m > M(n_0)` | **`T(m) = c_1(1-α) m^{(α-1)/2}`**; fitted `T ~ c/(M_0 Λ^β)`, `β ∈ (3/8, 1/2)`; `M_0 T → 0` | **Euler** | Kim–Jeong, **Prop 4.1, p.10** |
| 5 | **Prop 6.1** (`Ḣ^{d/2}` version): `sup_{[0,T(m)]} ∫\|Λ^{d/2}ω\|^2 >= (c_3/2)Σ_{n=m/2}^{m} n^{c_4-2α}` | (2.1) with `1/2 < α < (1+c_4)/2`, `m > M(n_0)` | `T(m) = c_1(1-α)(m/2)^{(α-1)/2}` | **Euler** | Kim–Jeong, **Prop 6.1, p.15** |
| 6 | **Thm 1.1 (Norm inflation)**: for every `q>1` and every `ε̃, δ, A > 0` there is `ω_0 ∈ C_c^∞` with `\|\|ω_0\|\|_{L^∞∩L^1} + \|\|ω_0/r\|\|_{L^{3,q}} <= ε̃` and `sup_{0<=t<=δ} \|\|ω(t)\|\|_∞ >= A` | axisymmetric, no swirl; multi-ring conical data (14); **every `q > 1`** | **"any `δ > 0`"** — qualitative, no rate | **Euler** | Bang–Cheskidov, arXiv:2605.16502v1, **Thm 1.1, p.8** |
| 7 | **Thm 1.2 (Instantaneous blow-up)**: same data class; (1) no Yudovich-type weak solution on any `(0,T)`; (2) a distributional solution with `ess sup_{0<t<T}\|\|ω(t)\|\|_∞ = +∞` for every `T` | as above with **infinitely many rings**; `ω_0 ∈ C^∞(R^3\{0})`, compact support | **`T → 0^+`** (instantaneous) | **Euler** | Bang–Cheskidov, **Thm 1.2, p.11** |
| 8 | Self-slowdown: `u^r/r ≲ \|\|ω\|\|_∞ H/R` — flattening depletes the stretching coefficient | axisymmetric ring of radial extent `R`, height `H` | — (mechanism, not a clock) | **Euler** | Bang–Cheskidov, **§1.4, p.6** |
| 9 | **Thm 8**: unique mild solution `ω ∈ C_w([0,T],L^∞)`, **`T >= 1/(c\|\|ω_0\|\|_∞)`**, `c` absolute. *Existence time only* — no growth bound | `ω_0 ∈ L^2 ∩ L^∞`; `L^2` declared a "soft assumption", "no quantitative dependence on `\|\|ω_0\|\|_2`" | **`T >= c/M_0`, no log, no energy** | **NS** (`ν=1`) | Bradshaw–Farhat–Grujić, arXiv:1704.05546v4, **Thm 8, p.8**; ARMA (2019) 10.1007/s00205-018-1314-5 |
| 10 | **Thm 10**: same, plus analyticity in `\|y\| < sqrt(t/c(M))` and **`\|\|ω(t)\|\|_{L^∞(Ω_t)} <= M \|\|ω_0\|\|_∞`** | `ω_0 ∈ L^2 ∩ L^∞`, `M > 1`, `c(M) > 1` | **`T >= 1/(c(M)\|\|ω_0\|\|_∞)`** with growth factor `M` | **NS** | Bradshaw–Farhat–Grujić, **Thm 10, p.11**. ⚠ proof is a self-described *sketch*; the BMO local-average term (pp.9–10) is dropped — that term is where `log Re` enters |
| 11 | **Thm 1.3**: `∃ t_* > 0` and smooth NS solution with `\|ω(t_*)\|_∞ / \|ω_0\|_∞ >= M`, any `M`; the first rigorous NS vorticity growth ratio | `u_0 ∈ C_c^∞`, small in a prescribed supercritical Besov/Sobolev norm; mixing construction | **`M_0 t_* ≈ ε^{-N} → ∞`** (long clock) — says nothing at `c/M_0` scales | **NS** | Luo, arXiv:2504.08288, **Thm 1.3, p.3**; params (2.10)–(2.11), **p.9**; proof **§6.2** |
| 12 | Background statement: *"no rigorous lower bounds for `\|ω(t)\|_∞/\|ω_0\|_∞` have been established for solutions of the Navier-Stokes equations to date"* | — | — | **NS** | Luo, arXiv:2504.08288, **p.3** |
| 13 | **Thm 1.1**: `\|\|ω(t,·)\|\|_{L^∞} <= A(ω_0^θ)(1+\|t\|)^{4/3}` for all `t ∈ R` (Childress rate) | axisymmetric no swirl, `ω_0^θ` compactly supported, `\|\|r^{-1}ω_0^θ\|\|_∞ < ∞` | **long time**; upper bound; uses energy + conserved quantities | **Euler** | Lim–Jeong, arXiv:2409.19497v2, **Thm 1.1, p.1**; ARMA 249:32 (2025) |
| 14 | Danchin's global well-posedness (the theorem the above ill-posedness results are sharp against) | `ω_0 ∈ L^∞ ∩ L^{3,1}`, `ω_0/r ∈ L^{3,1}(R^3)` | global | **Euler** | Danchin, Russ. Math. Surv. 62(3) (2007) 475–496, as stated in Bang–Cheskidov **(8), p.4** |

Legend: `M_0 = ||ω_0||_∞`, `Λ = log(R_outer/R_inner)`, `Re_E = M ℓ^2/ν`.

## Files

- `pdf/` — the five source PDFs as downloaded (`jeong-kim-2111.14078.pdf`,
  `bang-cheskidov-2605.16502.pdf`, `bfg-1704.05546.pdf`, `2504.08288.pdf`,
  `2409.19497.pdf`, `2511.03171.pdf`, `2512.13456.pdf`).
- `txt/` — `pdftotext -layout` extractions.
- `scripts/pagefind.sh` — locates the PDF page containing a given string.
- `scripts/kj_scale_audit.py` — the arithmetic behind §1c (run with `python3`; also
  `python3 -c "import sys;sys.path.insert(0,'scripts');import kj_scale_audit as k;k.q4()"`).

# REFUTER of `s3close/round2/u1` — the slaved reference map

Seat `s3close/round2/refute-u1`, sitting of 2026-09-08.  Single refuter.
Laws: `TORMENT NEXUS/LAWS.md` (first 120 lines read).
Target: `s3close/round2/u1/PROOF.md` (46,887 bytes) and `u1`…`u7`, `check_constants.py`.

**Nothing in `u1/` was edited.**  Every script of `u1` was re-run from a byte copy in a scratch
directory.  My own instruments sit beside this note — `r1_constants.py` (the angular strain
constants from scratch), `r3_shells2.py` (the integro-ODE with a non-quantised shift),
`r4_CR.py` (the velocity remainder with its geometry kept), `r6_budget_sharp.py` (the budget at
the corrected `C_R`).  They import nothing from `u1`…`u7` except, in `r6` only, `u5_budget`
used as the *object under test* — L-14's inverse convention: import the object when the point
is to test that same object, re-implement when the point is to test the same claim.
`r6_budget_sharp.py` expects `u5_budget.py` on `sys.path`; edit the one line at its head.

**Read (read-only):** `u1/PROOF.md` in full; `u1/u1…u7`, `check_constants.py`;
`s3close/assembly/ASSEMBLY.md` §0–3.5; `rebuild/far-near-kernel-lemma/NOTE.md` in full;
`write/lemma-T-shell-dependent/PROOF.md` §0–4; `lower/prove-lagrangian/NOTE.md` §3–4;
`s3close/hk2/PROOF.md` §0; `write/L3v-and-gamma-bound/PROOF.md` §B5–B6.
That last file is the one that matters most below and **`u1` did not read it** — its declared
list does not contain it, and it quotes that seat's `G_far = 1.6049285`, `G_inner = 0.1324254`
second-hand from `ASSEMBLY` §2.6 BLOCK 3.

---

## 0. VERDICT

**The headline survives.  `pc = O(1)`, and `e^{34}` is genuinely gone.**

The redefinition does what the note claims.  `𝔞(ρ,s) := ∫_{2ρ}^∞ (−(3/5)g₁[ω^θ(·,s)])(ρ′) dlogρ′`
is a functional of the *true* solution, so it is well defined without reference to the
bootstrap; the `l = 1` rate error is identically zero rather than estimated away; the majorant
(2.3) is a linear ODE with an `O(1/L)` forcing and **no `ε_{T′}` anywhere on its right-hand
side**; and the a priori cap `λ_max = e^{3c/4}` comes from the contradiction hypothesis alone,
so there is no hidden circularity.  I re-derived the driver, both bounds of the majorant, the
closed form, the Jacobian equation and the `l = 1` identity independently, and they hold.

**What does not survive is the size of `L_*`.**  The theorem in §7 assumes **(A-cone)** — that
the majorant is established on `{δ ≤ φ ≤ π−δ}` only — and in the same paragraph quotes
`L_* = 16161.986`, a number computed at `φ = φ₀ = 30°`, which is not the supremum over that
cone.  The note's own `u7` already contains the consistent figure, `52437.352`.  Everything in
§0's headline table, §6.3, §6.4, §6.5, §6.6 and §6.7 is at `φ₀` and is therefore not the value
the theorem's own hypotheses support.

**Severity roll-up:** no FATAL. **3 MAJOR**, **11 MINOR**.  `pc = O(1)` **survives**;
`L_* = 1.6·10⁴` **does not** — the defensible figures are `52437.4` (the note's own hypotheses,
honestly applied) and `53653.8` (my sharpened `C_R`, a genuine supremum over **every** angle
with no cone restriction at all, F3).  Both are still ten orders below `ASSEMBLY`'s
`4.53·10^{14}`, so the qualitative conclusion — the restart, BLOCK 2, is no longer the binding
constraint — stands.

**(S3) is not proved.  `FL-000` stands.**

### 0.1 Reproduction

All seven scripts re-run from a copy produce **bit-identical** JSONs
(`u1`…`u7_results.json`, `diff` clean), and `check_constants.py` returns
`95 CHECKS, 95 PASS, 0 FAIL`, `468 numeric tokens, 0 untraced`.  I independently reproduce
`κ_δ(7.5°) = 0.4997212`, `P₁(λ) = λ` (exact), `r_h[1,3/2] = 0.9818226`,
`r_h[1,λ_max] = 0.9199400` (u1: `0.9199394`), `λ_max = 1.8377407`, `C_a(30°) = 8.6978881`,
`C_a(90°) = 4.6986701`, `C_R(30°) = 29.024165`, `C_R(7.5°) = 98.289489`,
`ℓ_shell = 1.9102211`, `ℓ_loss(f=1,μ=0) = 3.2119052`, `μ(c)L = 446.2052034`, and the six
`L_*`/`log Λ_*` entries of §6.4.  The arithmetic of the note is sound.

---

## 1. FINDINGS

### F1 — **MAJOR**.  The theorem's `L_*` is inconsistent with the theorem's own hypothesis.

§7 states, as a hypothesis, **(A-cone)**: *"the majorant is established on `{δ ≤ φ ≤ π−δ}`
only (§6.8)"*, and then concludes *"with `L_* = 16161.986` on all-proved constants
(`G_collar = 0`, `φ = φ₀`)"*.  Those two clauses are incompatible.  `μ` and `μ_J` are suprema
over `x ∈ S` (Lemma T′ (H2), (H3)); if the region on which the majorant holds is the cone
`{δ ≤ φ ≤ π−δ}`, the drive constant is `sup` of `C_R(φ)` over that cone, i.e. `C_R(δ)`, and
§6.8's own table gives the answer.

> **CORRECTED STATEMENT.**  Under (A-cone) as written, and with everything else exactly as in
> §7, `L_* = 52437.352` (all-proved), `log Λ_* = 1.0487807·10⁵`, and `L_* = 9884.164`
> (all-measured).  The improvement over `ASSEMBLY`'s `4.5253·10^{14}` is `8.6301·10⁹`, not
> `2.7999653·10^{10}`; the `G_collar` and `λ_max` sensitivity tables of §6.5, the `C_K` sweep of
> §6.6 and the two variant tables of §6.7 all move by the same factor `3.2445`.

I reproduce both numbers by driving `u5_budget.assemble` through the same `G_collar` slot `u7`
uses (`L_*(proved) = 16162.0` at `φ₀`, `52437.4` on the cone; `L_*(measured) = 600.3` and
`9884.2`).  The note *has* the right number in `u7_results.json`; it does not headline it.

The note is not hiding this — §6.8 states the problem plainly and calls it *"a real hole … where
I would attack this note first"*.  The finding is that the *theorem statement and every headline
number* were nevertheless written at `φ₀`.

### F2 — **MAJOR**.  §7.1's status table contradicts itself: T12 is `PROVED`, T17 is `NOT ESTABLISHED`, and T12 depends on T17.

| row | printed | correct |
|---|---|---|
| T12 | Lemma T′ applied once … **PROVED (R1)** | **PROVED-modulo-(A-cone)** |
| T14 | (5.1) and the quasimonotone comparison … **PROVED** | **PROVED-modulo-(A-cone), and modulo T2** |
| T17 | the angular supremum … **NOT ESTABLISHED** | correct |
| T18 | the budget … **COMPUTED (conditional on T7, T8, T17)** | correct |

Lemma T′'s hypotheses (H2), (H3) are quantified *for all `x ∈ S`*, and Step 4 uses (H2)
pointwise inside `∫_S`.  A `μ` that is `+∞` on `{φ < δ}` is not a `μ`.  T12 cannot carry a
`PROVED` while T17 carries `NOT ESTABLISHED`; nor can T14, which is T12 applied at every shell.
Under the brief's default (`refuted = true` if uncertain) these two rows are refuted as printed.

### F3 — **MAJOR** (a repair, in the note's favour).  The `1/sin φ` in `|R_y|` is spurious, the `log(1/sin φ)` in `|R_z|` is exactly removable, and with a constant already on the record the angular supremum is **finite** — so (A-cone) is closeable and `C_R` is a factor `1.67` too large at `φ₀`.

Three separate points, all checked in `indep/r4_CR.py` (sympy + quadrature).

**(a) The `r`-component does not diverge.**  `R_y = (a(X) − 𝔞(|X|))·Y` and `|Y| = r = ρ sin φ`,
so
```
   |R_y| ≤ C_a(φ) M_s · ρ sin φ = M_s ρ [ (C1 + C2in + π/8) sin φ + 3.999218 ] ,
   sup_{φ∈(0,π)} C_a(φ) sin φ  =  C_a(90°)  =  4.6986701      (attained at the equator).
```
The collar's `3.999218/sin φ` cancels against `|Y| = ρ sin φ` **exactly**.  §3.1's
`|R_y| ≤ C_a M_s|X|` throws that away, and §6.8's whole `C_a(φ) → ∞` story is, for the
`r`-component, an artefact of dropping `sin φ`.  (`ASSEMBLY` §2.4 drops it too.)

**(b) The `log(1/sin φ)` in `C_R` is not there.**  §3.2 bounds the middle integral by
`|z| log(1/sin φ)`.  The integral is elementary and finite:
```
   ∫₀^{|z|} log(ρ/ρ_ζ) dζ  =  |z| − r arctan(|z|/r)  ≤  |z|        EXACT (sympy, residual 0)
```
so the term is `(1 − (r/|z|)arctan(|z|/r)) M_s|z| ≤ M_s|z|`; at `φ₀ = 30°` it is `0.3954`, not
`0.6931472`, and it is bounded uniformly in `φ` instead of diverging.

**(c) The collar's axis divergence is already repaired on the record, with a proved constant.**
`write/L3v-and-gamma-bound/PROOF.md` §B5, *"The collar constant that is not needed"*: for a
datum with `|ω^θ| ≤ λM min(1, φ̃/δ)`, using `φ ≤ (π/2)sin φ` so that the `1/r′` weight cancels,
```
   |a_collar| ≤ (3/(8π²))(π/δ)|S⁴| R_A λM = 6.281958 λM/δ = 47.99 λM   at δ = 7.5° ,
```
crossing over with `4/sin φ` at `sin φ = 0.08333`.  `u1` does not read that file and writes
(§6.8) that closing the taper cone *"needs `far-near-kernel-lemma` §6's observation … turned
into a proved constant"* — it already is one.

**Assembled.**  With `min(3.999218/sin φ, 6.281958/δ)` for the collar and the geometry of
(a),(b), the exact `ζ`-integrals give a **uniform** majorant:
```
   sup_{φ ∈ (0,π)}  C_R(φ)  =  100.6122          (no cone restriction at all)
   C_R(φ₀ = 30°)             =   17.4007         (u1: 29.024165)
   sup_{φ ≥ δ}  C_R(φ)       =   40.2293         (u1's cone value: 98.289489)
```
| `C_R` | `μ(c)L = λ_max³(C_R+2logλ_max)(λ_max²−1)` | `L_*` proved (`ε≤½`) |
|---|---|---|
| `29.0242` — u1 headline, `φ₀` | `446.205` | `16161.99` |
| `98.2895` — u1 cone, `φ = δ` | `1468.205` | `52437.35` |
| `17.4007` — sharpened, `φ₀` | `274.703` | `10075.1` |
| `40.2293` — sharpened, sup over the cone | `611.535` | `22030.1` |
| **`100.6122` — sharpened + taper-repaired collar, sup over ALL `φ`** | `1502.477` | **`53653.8`** |

What is still owed to close (A-cone) outright is a transport statement — that the true field
keeps a taper, `|ω^θ(·,s)| ≤ M_s min(1, Cφ/δ)` on the taper cone — which is far weaker than
BLOCK 3, and is plausible because the strain moves material *away* from the axis
(`tan φ_new = λ³ tan φ_old`).  Note also that the taper repair removes the last obstruction to
`(H2)` itself: *without* it, the sharpened `μ(φ) ≈ 59.0/(φL)` still exceeds `1` in a shrinking
cone `φ ≲ 59/L` (`0.21°` at `L = 16162`), where Lemma T′'s `(1−μ)^{-5}` is vacuous; *with* it
`C_R` is capped at `100.6122`, so `μ ≤ 1502.5/L < 1` at every `L ≥ 1503` and the bad cone does
not exist.

### F4 — **MINOR**.  `κ_δ` omits the equatorial mollification of the datum it is applied to.

`ASSEMBLY` §1.1's datum — which §7's theorem adopts verbatim — carries
`min(1, |φ−π/2|/δ_m)` with `δ_m = 5°`.  `κ_δ = ½P_{h_δ}(1)` is computed from the **axis taper
alone**.  Recomputed here from scratch (mpmath, 30 dps) with the full angular profile:
```
   κ (axis taper only, 7.5°)      = 0.4997212362      (u1: 0.4997212305)
   κ (full datum, 7.5° and 5°)    = 0.4978226672      (u1's own u2 computes 0.4978224382
                                                       as `kappa_full_7.5_5.0` and never uses it)
   relative difference            = 3.8137·10^{-3}
```
Every consequence is in the optimistic direction:
`c_* = log(3/2)/κ` becomes `0.8144770` (u1: `0.8113826`); `c₂ = 2log(3/2)/κ` becomes
`1.6289540` (u1: `1.6227652`); `λ_max = e^{3c/4}` becomes `1.8447`; `r_h[1,λ_max]` becomes
`0.9198409`.  Inherited from `ASSEMBLY` §1.2, which has the identical defect.  Not load-bearing
(0.4 % on the clock constant), but the theorem's `c₂` as printed is wrong for its own datum.

### F5 — **MINOR**.  `P_h` increasing is a finite-grid check labelled `PROVED`, and it is load-bearing.

Row T2 prints `PROVED`.  The evidence is `u2`'s `np.diff(P_TAB) > 0` on a grid and `u4`'s
`assert (np.diff(P_TAB) > 0).all()` on a 1051-point table, with `"P_h_argmax_above": 2.075` a
hard-coded literal.  Monotonicity is what makes §5.1's discretised system quasimonotone
(`∂f_i/∂x_j = (1−ε)½P_h′(λ_j)λ_jΔσ ≥ 0`), so the comparison principle — and with it T14, §5.2
and §5.3 — rests on it.  My own scan (mpmath, 25 dps, step `0.005`) agrees: the first decrease
is at `λ = 2.08`.  A proof is two lines and should replace the grid:

> `P₁(λ) = λ` exactly, so `dμ_λ(v) := 3v²(Av²+B)^{-5/2}dv` (`A = λ²−λ^{-4}`, `B = λ^{-4}`) is a
> positive measure on `[0,1]` of total mass `λ`; `h_δ(arcsin v)` is non-decreasing in `v`;
> `μ_λ` is stochastically increasing in `λ` on `[1,λ_max]`, hence `P_{h_δ}(λ) = ∫h dμ_λ` is
> increasing there.

### F6 — **MINOR**.  §5.2's 200-shell table is quoted to seven digits and is converged to three.

I integrated (5.1) independently (`indep/r3_shells2.py`), treating the shift `d` **continuously**
— cumulative trapezoid plus interpolation of the antiderivative at `σ+d` — rather than rounding
`σ_i + d` to a cell boundary, and Richardson-checked it at `N = 200…3200` (stable to `1·10^{-8}`):
```
   theta at lambda(0) = 3/2, system (A):
        L = 40 :  0.7796043   (u1, N=200 : 0.7816769,  +0.27 %)
        L = 160:  0.7460319   (u1, N=200 : 0.7479297,  +0.25 %)
        L = 640:  0.7380526   (u1, N=200 : 0.7399101,  +0.25 %)
```
`u1`'s cut is quantised to a cell: at `L = 640`, `d/Δσ = 0.597`, so the scheme's effective `d`
is `0` or one cell.  My `θ_C` reproduces `u1`'s to eleven digits (`0.8543167310`,
`0.8232696872`, `0.8158573431`) — the bias is in (A) alone.  §5.4's *percentages* survive
because comparator (D) shifts with (A) (I get `1.08 %`, `0.38 %`, `0.19 %` against `u1`'s
`1.0852265 %`, `0.369 %`, `0.1945463 %`), and §5.3's `c_*` — the number the theorem uses — is
untouched.  Recorded because the table reads as a converged computation and is quoted as one.

### F7 — **MINOR**.  §6.2's `ε_{T′} ≈ (15π/4)μ/r_h = 12.8062478 μ` assumes `μ_J = μ`, which is false here by a factor 20.

Lemma T′'s `(T′-rel)` limit `(15π/4)μ` is stated *for `μ_J ≤ μ → 0` with `μ_J ≤ μ`*, and is
attained at `μ_J = μ`.  Here the note's own numbers are `μ·L = 446.2052` and
`μ_J·L = 22.6532`.  The correct linearisation is
```
   ε_{T′} = (2π/r_h)[ 3μ/2 + 3μ_J/8 ] + O(μ²) = 10.2450 μ + 2.5613 μ_J .
```
`eps_Tprime()` in `u3`/`u5` uses the exact expression, so **no computed number moves**; the
prose's closure criterion *"the window closes when `μ ≲ 0.04`"* should read `μ ≲ 0.048`.

### F8 — **MINOR**.  T9's "agreeing with the solver to `1e-13`" is not true, and is not what the gate checks.

```
   closed form  λ_max³(C_R + 2 log λ_max)(λ_max² − 1) = 446.2052033630385
   u5 solver    feedback.proved.mu_times_L            = 446.2052364675315
   relative                                            = 7.4191·10^{-8}
```
`check_constants.py` asserts this at `rtol = 1e-6`.  The closed form itself I confirm and it is
exact: at `L → ∞`, `G = 3/2` and `e^{Gc} = e^{3c/2} = λ_max²`, so
`m(c)L = λ_max(C_R + 2logλ_max)(λ_max² − 1)` and `μ = λ_max² m`.

### F9 — **MINOR**.  The Lipschitz constant is that of the 5-D lift `b`; the hypothesis (Γ-off) is about `∇u` on an undefined "slab"; and the bound is needed on a segment, not on two trajectories.

§2.2 sets `Γ(s) := ‖∇b(·,s)‖_{L^∞(slab)}` and then invokes (Γ-off), stated for
`‖∇u(s)‖_{L^∞(slab)}`.  The identification is in fact **exact for the operator norm** and is
worth writing down: in the frame `(ê_r,ê_θ,ê_z)` the 3-D gradient is
`[[∂_r u^r,0,∂_z u^r],[0,u^r/r,0],[∂_r u^z,0,∂_z u^z]]`, and the 5-D `∇b` is the same `2×2`
block plus three diagonal copies of `a = u^r/r`, so `‖∇b‖_op = ‖∇u‖_op` identically.  Two
things are still missing: (i) that sentence; (ii) the domain — `|b(Φ)−b(Λ)| ≤ Γ|Φ−Λ|` needs the
gradient bound on the whole segment `[Λ_s x, Φ_s x]` for every `x ∈ S`, and "slab" is defined
nowhere in `u1`, in `ASSEMBLY` §2.2, or in R3.  (Theorem Γ's own slab is
`{2λρ₀ ≤ |x| ≤ R/(2λ²)}`, which excludes both ends of the shell.)

### F10 — **MINOR**.  (V-strain)'s "`≤ 0.32 %` in `L_*`" is an assumption, not a bound, and §6.7's "it is bookkeeping" over-states what was shown.

Theorem V.4 is pointwise, on a tube around the tracked point.  `a[·](0)` is a weighted integral
of `η` over the **whole** shell, and the viscous defect at the axis taper, at the equatorial
layer and at the two radial edges is not covered by it.  §6.7 charges the *tracked point's*
`ε_v` to the strain integral and reports the result as the cost of closing the gap.  The
scaling is favourable (`√(ντ) = ρ₀δ√(c/L) ≪ ρ₀δ`, and the relative defect at radius `ρ` falls
like `√(ντ)/(ρδ)`), so the conclusion is probably right; but §6.7 measures what it would cost
*if the answer is what we assume*, which is not an upper bound.  Keep it as the carried gap
§4.2 already names, and delete "it is not a load-bearing gap".

### F11 — **MINOR**.  The theorem's datum is `tanh`-mollified radially; the budget is at `ε = 0` and `ℓ_loss` carries no radial-ramp cost.

§7 adopts hk2's correction (`w₀ = 0.25ρ₀`), which removes `≈ log(1.25) = 0.2231` e-folds at the
inner edge from `a_ref = (M/2)∫Θ(ρ)P_h(λ)dlogρ`.  `ℓ_loss = 3.2119052` (at `f = 1`, `μ = 0`)
does not include it.  Size `0.2231/L`, i.e. `1.4·10^{-5}` at `L_*`.  Inherited: `ASSEMBLY` §2.5
says the same ("the budget below is computed at `ε = 0`").

### F12 — **MINOR**.  T1's trajectory margins are computed for `Λ`, not for `Φ`.

`φ(λ_max) = 74.407°` and the `15.593°` equator margin come from `tan φ = λ³tan φ₀`, the
*reference* map.  The true point sits inside a relative `μ` ball of it.  At `L = L_*`,
`μ = 0.0277`, so the angle can slip by `≈ 1.6°`: `sin φ ≥ 0.4756` instead of `0.5`, i.e.
`C_a` up by `4.7 %`, and `d` (the distance to the nearest datum discontinuity, which the whole
viscous column rests on) down by the same order.  It does not break the margins at the `L` that
matters, but T1's `PROVED` is for `Λ`.

### F13 — **MINOR** (gate hygiene).  One of the 95 checks is an algebraic tautology.

`check_constants.py`:
```
   DISPLAY["exp_ratio"] = math.exp(34.05725232185504 - 1.5*CST)
   chk("ratio of the two exponentials",
       math.exp(34.05725232185504)/math.exp(1.5*CST), DISPLAY["exp_ratio"], rtol=1e-9)
```
This asserts `e^{a−b} = e^a/e^b` and tests nothing about the mathematics.  §8's claim *"No check
compares a literal against itself"* is therefore not quite right; one of the five that were
caught on the first run came back in a different disguise.

### F14 — **MINOR**.  `u7`'s cone pricing raises `C_R` but leaves `C_a` at `φ₀` in the Jacobian equation.

`μ_J = sup_x|J_Φ/λ² − 1|` is also a shell supremum, so the `C_a` in
`dw/dθ = (2λ_ω/L)[C_a + ½Λ_rad]` must be `C_a(δ) = 31.3387` on the cone, not
`C_a(30°) = 8.6979`.  `u7` feeds the extra only through the `G_collar` slot, i.e. only into
`C_R`.  The effect is small — `μ_J·L` goes from `22.65` to `≈ 77.8`, moving `ε_{T′}` by about
`1.3 %` — so `52437.352` is understated but not materially.  (`ε_{C_a} = λ_ωC_a/(κ_δL)` at
`φ₀` **is** correct: that one is the deficit at the tracked point, whose trajectory does stay
in `[30°, 74.4°]`.)

---

## 2. WHAT I TRIED TO BREAK AND COULD NOT (severity NONE)

Taken in the order the brief set.

**(1) §1.1–1.2, §4.1 — is `𝔞` the exterior `l = 1` strain, is `Λ` the flow of `𝔞·DX`, and is the
origin identity exact?  YES, YES, and YES; nothing has been moved into `ℓ_loss`.**

* `F(ρ′,s) = −(3/5)g₁[ω^θ(·,s)](ρ′)` is *verbatim* the far/near lemma's Consequence A shell
  functional (`Φ[w](0) = −(3/5)g₁ = −(3/4)∫₀^π w cos φ sin²φ dφ`, `|·| ≤ M/2`), and the split at
  `2ρ` is that lemma's own split (`a_far(x) := ∫_{2ρ}^R Φ[w(ρ′·)](x/ρ′) dρ′/ρ′`).  I re-derived
  the normalisation from `N₁ = 12/5`, `C₁^{3/2}(t) = 3t`: `−(3/5)g₁ = −(3/4)∫w cos φ sin²φ dφ`.
* **(P1)** `𝔞 ≥ 0`: Lemma 2 gives `η ≤ 0` on `{z>0}` and `z`-oddness gives `ω^θ cos φ ≤ 0`
  everywhere, so the integrand of `Φ[w](0)` is pointwise `≤ 0` and `F ≥ 0`.  Correct — and note
  the positivity holds after restriction to **any** set (the indicator is `≥ 0`), so §4.2 step 2
  does not actually need `S^out` to be reflection-invariant.
* **(P2)** `∂𝔞/∂logρ = −F(2ρ,s)` is exact (differentiation in the lower limit), `|F| ≤ M_s/2` is
  Consequence A, and `κ_s := sup L|ρ∂_ρλ/λ| ≤ ∫₀^τ(M_σ/2)dσ·L = (3/4)c = 0.6085369` follows.
  Step 0's condition `2κ_s/L < 1` is `L > 1.2170739`.  All correct.
* **(P3)** `𝔞(ρ,s) ≤ (M_s/2)log(R/2ρ) ≤ (3/4)ML`, so `λ ≤ e^{3c/4} = 1.8377407`; `𝔞` is
  decreasing in `ρ` because `F ≥ 0`.  Correct, and — the load-bearing point — it uses **only**
  the contradiction hypothesis `M_s ≤ (3/2)M`, so `λ_max` is genuinely a priori and the loop
  really is broken.
* `∂_sΛ_s(x) = (∂_sλ/λ)DΛ_s(x) = 𝔞(|x|,s)DΛ_s(x)` is exact (`∂T_λ/∂λ = λ^{-1}DT_λ`).  So `Λ` is
  the flow of `v(X,s) = 𝔞(|Λ_s^{-1}X|,s)DX`.  §1.2's headline *"and that field is the `l = 1`
  mode"* is loose — `v` is a **shell-labelled** coaxial strain, not an `l = 1` solid harmonic,
  and it is not 3-D divergence-free (the note's own T4 gives `J_Λ = λ²[1+(ρλ′/λ)(1−3cos²φ)]`,
  and `ASSEMBLY` §2.1 says so explicitly).  Nothing downstream uses more than the flow identity,
  so this is wording, not a defect.
* **§4.1 is exact and nothing is hidden.**  At `x = 0` the interior expansion's `|ξ|^{l−1}` kills
  every `l ≥ 2` term, so `a(0,s) = ∫_{ρ₀}^{R}F dlogρ′` and
  `𝔞(ρ,s) = a(0,s) − ∫_{ρ₀}^{2ρ}F dlogρ′` are identities.  The `l ≥ 2` interior remainder of the
  far shells does **not** appear in them and has **not** been moved into `ℓ_loss` — it is in
  `R`, priced by `C1 = 0.291999` inside `C_a`.  And §4.2 does not use the identity at all: it
  uses the *restriction* route (step 2), `𝔞(ρ,s) = ∫_{2ρ}^∞F ≥ ∫_{2ρ}^∞F^{out} = ∫_0^∞F^{out} =
  a[η1_{Φ(S^out)}](0)`, which is exact given `F ≥ F^{out} ≥ 0` and the label cut.  The brief's
  worry does not land.

**(2) §2.1–2.3 — the driver, re-derived.  Three terms, and no `ε_{T′}`.**

I re-derived (2.2) from scratch:
`∂_s(Φ−Λ) = b(Φ) − 𝔞(|x|)DΛ = [b(Φ)−b(Λ)] + R(Λ) + [𝔞(|Λ|)−𝔞(|x|)]DΛ`, using
`b(X) = 𝔞(|X|)DX + R(X)` as the *definition* of `R`.  (§2.1's sentence *"by Consequence A the
shells outside `2|X|` contribute exactly `𝔞(|X|,s)·DX` at `X`"* is wrong as literally written —
they contribute that plus the `l ≥ 2` far Taylor remainder — but since `R` is defined as the
remainder, the split is trivially an identity and the mis-statement costs nothing; §3.1 puts the
remainder where it belongs.)  **No term proportional to Lemma T′'s sensitivity survives.**  The
three bounds check out:
```
   |b(Φ)−b(Λ)| ≤ Γ|Φ−Λ|                                             (subject to F9)
   |R(Λ)|      ≤ C_R M_s|Λ| ≤ C_R M_s λ_max|x|
   |[𝔞(|Λ|)−𝔞(|x|)]DΛ| ≤ (M_s/2)·2 log λ_max · |DΛ| ,  |DΛ| ≤ 2|Λ|
      (verified: |DΛ|² = λ²r² + 4λ^{-4}z² ≤ 4(λ²r² + λ^{-4}z²) = 4|Λ|²)
   ⇒ dm/dθ = G m + λ_max λ_ω (C_R + 2 log λ_max)/L ,  μ = λ_max² m  (|Λ| ≥ λ^{-2}|x|) .
```
`G = (3/2)(1+C₁/L) + C″/L` follows from `2a(0,s) ≤ 2(M_s/2)L ≤ (3/2)ML`.  The closed form is
exact (F8) and I reproduce `μ(c)L = 446.2052034` and `pc(L=640) = 1.4106013`.  The Jacobian
equation is right and **Corollary 3 is genuinely not used**: I re-derived
`∂_s log J_Φ^{(5)} = div₅b(Φ) = 2a(Φ)` from `div₅b = 4a + r∂_ra + ∂_zu^z` and 3-D
incompressibility `r∂_ra + ∂_zu^z = −2a`, and `∂_s log λ² = 2𝔞(|x|,s)` is the definition of `λ`,
so `(H3)` is obtained directly against `λ²` and the `2κ_s/L` of (2.11) is not paid.  This is the
cleanest part of the note.

**(3) §3.2 — the `z`-component, and the `r`-component pointwise.**

`∂_zu^z = −2a − r∂_ra` is correct (3-D incompressibility with `u^r = ar`); `u^z(r,0) = 0`
follows from Lemma 2's `z`-oddness; the three-term split of `R_z` and the monotonicity
`sinφ_ζ ≥ sinφ` are correct; `r∂_r𝔞(ρ) = sin²φ·∂𝔞/∂logρ` gives the `½`.  The `r`-component
**is** covered by the far/near lemma at every field point of the shell — its hypotheses are only
`|ω^θ| ≤ M_s` and `z`-oddness, both true for the true field at every `t` — including the
collar's angular dependence, which is real and not an artefact of the estimate
(`far-near-kernel-lemma` §6 exhibits the exact `c_edge` diverging on the axis for the bang-bang
datum).  What is *not* right is dropping the geometric factors: see F3.  §3.4's finding against
`ASSEMBLY` (that (2.1) charges the whole **vector** remainder to the **scalar** constant `C_a`)
**stands** — but the understatement factor is `C_R/C_a = 2.0006` with the sharpened constants,
not the printed `3.3369210`.

**(4) §4.2 — Lemma T′ applied once.**

The reference profile does satisfy Lemma T′'s `C¹`/shear hypothesis: `|ρλ′/λ| ≤ κ_s/L` with
`κ_s = (3/4)c` by (P2), and Step 0's `2κ_s/L < 1` holds at every `L` in play.  (Two small
technical points: `λ(·,s) ∈ C¹` needs `F(·,s)` continuous, which for `t > 0` is fine but at
`t = 0` the datum's `Θ_ε` kinks make `λ` only Lipschitz — Step 0's contraction argument needs
only Lipschitz, so this is harmless and should be said.)  `μ = λ_max²m` and
`μ_J = sup|J_Φ/λ²−1|` are the right objects for (H2),(H3), Corollary 2 applies to `η₀1_{S^out}`
with `A^out = ∫_{ρ_c}^Rλdlogρ` and `a_ref^out = (M/2)∫_{ρ_c}^RP_h(λ)dlogρ`, and the Jacobian is
handled correctly.  **Cor. 3 is not used** — confirmed.  The one place this is not proved is the
angular supremum (F1, F2, F3).

**(5) §5 — the comparison principle and the numerics.**

The two-line proof is correct as written.  The discretised system is not merely quasimonotone,
it is **cooperative with `f_i` independent of `x_i`** (because `d > 0` forces `j ≠ i`), so the
argument is even easier than claimed.  `P_h` increasing is what makes `∂f_i/∂x_j ≥ 0` — see F5.
The subtracted `−C_visc/L` is profile-independent and cancels between `x` and `y`, so it does
**not** reverse the direction; what it does require, silently, is `1 − ε_{T′} ≥ 0`, or else
`∂f_i/∂x_j ≤ 0` and quasimonotonicity fails — which is the closure condition anyway, so no row
of the budget is affected.  My independent integration confirms the comparison qualitatively
(A ≥ B, A ≥ C, A < D) and the percentages, and disagrees with the absolute `θ_A` at the `0.25 %`
level: F6.  §5.4's refusal of the brief's `t_*` formula is correct and I confirm the arithmetic
(`2∫₁^{3/2}dλ/(λP_h(λ)) = 0.6698262`; the brief's expression replaces `H(0,θ) = √λ` by
`P_h(λ(0))`, over-estimating the rate).  L-17 refusal properly registered.

**(6) §6 — the budget.**

`u5_budget.py` is `a2_budget.py` with `bootstrap()`, `COLUMNS`, `assemble()` replaced, as
declared; `viscous_budget()`, `gauss_tail_aniso()`, `J_pow()`, `C_a_proved()`, `eps_Tprime()`
and `RECORD` are carried over unchanged, so the viscous half is the assembly's own instrument.
I re-ran it from a copy: `u5_results.json` is bit-identical.  I re-assembled `ε` by hand at
`measured`, `L = 640`, `f = 1` from the §6.2 term table and get `0.431782` against the printed
`0.431788`, and the closure condition `ε_a ≤ 0.33296` at `L ≈ 16162` reproduces
`L_*(proved) = 16161.986` to three digits from an independent hand calculation.  `ε(640)`,
`ε(2560)`, all six `L_*` and all six `log Λ_* = 2L_* + 3.3643455` reproduce.  The pricing of
(A-cone) is **not** consistent — F1, F14.  §6.6's conclusion (that `(H-K2)` moves `L_*` by
`+0.015 %` proved and `+5.6 %` measured) is arithmetically right and unaffected by everything
above.

**(7) §7 — the `PROVED` rows.**  T3, T4, T5, T10, T11, T13, T15 are genuinely proved.  T1 is
proved for `Λ`, not `Φ` (F12).  T2's monotonicity clause is a grid check (F5).  T6 is proved
pointwise but its shell supremum is `+∞` as written and is finite when the geometry is kept
(F3).  T9's stated precision is wrong (F8).  T12 and T14 cannot be `PROVED` while T17 is `NOT
ESTABLISHED` (F2).  T7, T8, T16, T17, T18, T19 are honestly labelled.

---

## 3. REMAINING GAPS, restated

In descending order of what they cost.

1. **(Γ-off)** — `‖∇u(s)‖_{L^∞} ≤ 2a(0,s)(1+C₁/L) + C″M` for the *true* field.  This is
   `ASSEMBLY` BLOCK 3 verbatim, which that note calls fatal.  `u1` carries it as "being proved
   by a parallel seat".  It supplies the entire exponent `G`; without it there is no (2.3).
   Note that §2.4's conjugated form (2.6) does **not** use it — that is the structurally better
   estimate, and §6.7 shows it costs a factor `3.0` in `L_*` at present constants only because
   of the crude `|T_λ^{-1}R| ≤ λ²|R|`.
2. **(R-z)** — `|r∂_r[a − 𝔞]| ≤ G_collar M_s` on the collar.  The gradient half of BLOCK 3,
   restricted to the collar.  Enters `C_R` additively, so it is priced linearly (`×5` costs
   `×2.409` in `L_*`), which is the note's real structural gain and is correct.
3. **(A-cone)** — closeable; see F3.  What is owed is (i) the geometry `|R_y| ≤ C_a(φ)M_sρ sinφ`
   and the exact `ζ`-integral, both elementary; (ii) the L3v taper-repaired collar constant
   `6.281958/δ`, already proved; (iii) a transport statement that the true field keeps a taper.
   Then `C_R = 100.6122` uniformly, with no cone restriction.
4. **(V-strain)** — the viscous defect of the strain *functional*, not of `η` pointwise.  Not
   done in the record, not done in `u1`, and §6.7's pricing is an assumption (F10).
5. **(H3-V)** — `ASSEMBLY` BLOCK 5, inherited.
6. **The self-consistent cap.**  `λ ≤ (3/2)(1+O(1/L))` instead of `e^{3c/4} = 1.8377`.  §6.5
   prices it at `×1.94`.  It looks cheap: material on a shell where `h_δ = 1` has
   `|ω^θ| = Mλ_material ≤ (3/2)M`, and `λ_material ≥ λe^{-C_ac/L}`.  Worth writing.
7. **`ε_{T′}` under the time integral** (§7.2 item 5) — a further factor `(e^{Gc}−1)/(Gc)`.
8. **The restart (BLOCK 2)** — still open, and, as `u1` correctly says, no longer binding.

---

## 4. WHAT I WOULD DO NEXT

1. Rewrite §3 with the geometric factors kept (F3 (a),(b)) and the L3v collar repair (F3 (c)).
   That deletes (A-cone) as a hypothesis, replaces every `L_*` in the note by an *honest*
   supremum, and — because `C_R` is in the drive and not in the exponent — costs only `≈ ×3.4`
   against the current headline while removing an unproved restriction.
2. Recompute `κ` with the full angular profile (F4) and re-issue `c_*`, `c₂`, `λ_max`, `r_h`.
3. Prove `P_h′ > 0` (F5, two lines) and re-run §5.2 with a non-quantised `d` (F6).
4. Then the note's grade claim — FILED, no Solid, no Major — is right for the *content*, and the
   headline table becomes defensible.

---

`FL-000` stands.  Nothing here touches the headline problem.

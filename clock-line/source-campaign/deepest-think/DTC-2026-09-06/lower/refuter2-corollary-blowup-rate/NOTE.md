# NOTE — REFUTER 2 of `corollary-blowup-rate` — DTC-2026-09-06 / lower

Target: `../corollary-blowup-rate/`.
A first refuter note already exists at `../refuter-corollary-blowup-rate/`. I read it after
forming my own attack list; where we agree I say so, and I do not re-report its findings as
mine. My three independent contributions are §2 (the FL-043 blindness, with an explicit
witness the gate cannot see), §3 (the sharp explicit form the attempt's (c) is a weakening
of), and §4 (a strictly stronger published competitor than the one refuter 1 found, plus
closure of refuter 1's own open gap 3).

## 0. Verdict

**REFUTED — MAJOR.** Two independent grounds, neither of which touches the mathematics of
(a)/(b)/(c):

* **the novelty verdict is wrong at source** (§4): a published, peer-reviewed 3D
  Navier–Stokes lower bound at rate exactly `(T-t)^{-1}`, **log-free**, **for every `t<T`**,
  with an **absolute constant**, in a space with the **same scaling as `‖ω‖_∞`**, has existed
  since 2015 — McCormick–Olson–Robinson–Rodrigo–Vidal-López–Zhou, Theorem 4.1. It sat in the
  attempt's own search neighbourhood and in the file refuter 1 downloaded. Three of the four
  properties the attempt's §7 lists as row 0's advantages are matched by it;
* **the block-(B) gate is blind in the direction that matters** (§2): I exhibit an explicit
  point at which the mutant "`1/5 → 0.199`" is FALSE, and the attempt's sampler returns
  `0 violations` on **20 out of 20** seeds. FL-043: the gate returns the same verdict for a
  true statement and for a false one that differs in the load-bearing exponent.

**The proof itself survives.** I re-derived every step against the input theorem's own text
and could not break any of them; §3 shows the one place where (c) is *lossy*, not wrong, and
gives the sharp replacement.

## 1. Reproduction, and what survived my attack

`shasum -a 256 -c SHA256SUMS` in the attempt's folder: **14/14 OK**. Both scripts re-run from
my copy `rerun/` reproduce `s1_log.txt` and `s2_log.txt` **bit-for-bit** (`diff` empty;
`rerun/my_s1.txt`, `rerun/my_s2.txt`).

I read the input theorem myself
(`~/<core-line-working-tree>/2026-09-06/new-chat/outputs/navier-stokes-pass8/LOGARITHMIC_RECORD_CLOCK.md`).

* **Line 2 is licensed, and by more than the attempt claims.** §1 boxes
  `H_0 = c/(M_0(1+log_+ R_0))` and says the solution *"continues through the interval of
  length"* `H_0`. §5 discharges the maximal-endpoint case explicitly: *"Restarting
  sufficiently close to the endpoint and using uniqueness extends the same solution beyond
  it. This contradicts maximality … the mild restart retains finite `L^2` norm … the usual
  `H^m` energy estimate preserve[s] the original strong regularity across the former
  endpoint … Thus the actual solution reaches `H`."* So the class question I came to attack —
  that the Giga–Inui–Matsui restart delivers only a bounded-velocity mild solution, while
  `T_max` is defined in `C([0,T);H^m_σ)`, so Lines 3–4 would need an unstated
  persistence-of-regularity step — **is answered inside the input theorem**, not by the
  corollary. No gap. (It is the input theorem's step, so it inherits the input theorem's
  unverified status, not the corollary's.)
* **Consequence, and a calibration point.** Because §5 already concludes "the actual solution
  reaches `H`" *for the actual solution*, `T_max ≥ t + H(t)` is immediate; the five-line
  contradiction of §3.1 is an unwinding of a definition, not a step. The mathematical content
  of the corollary is entirely §3.3 (the bootstrap) and §4 (the integration). The attempt says
  as much ("exactly as strong as its input"), and the literature says it in general —
  McCormick et al., p.2: *"In fact this result, and all subsequent lower bounds, are a
  consequence of upper bounds on the local existence time"* — but the note presents §3.1 as
  "PROOF, five-line continuation argument", which reads as content it is not.
* **The claim "never uses the `(3/2)M_0` cap" misdescribes the input.** In §5 the continuation
  is not separable from the cap: the argument runs *inside* the first-exit interval of
  `‖ω‖_∞ < 2M_0`, and it is the cap bootstrap that excludes a first exit before `H`. What is
  true is that the corollary *invokes* only the continuation conclusion. Harmless, but the
  sentence claims a modularity the source does not have.
* Lines 0, 1, 5, form (b), form (c)'s bootstrap and §4's antiderivative: all checked by hand
  and re-derived independently (§3, §5 below). Correct.

## 2. FL-043 — the numerics gate is blind where it matters (MAJOR, new)

The attempt's block (B) is a 4,000,000-draw uniform search over `(E_0, ν, s, c_1, M)` reporting
`(b) true and (c) false = 0`, plus a control (`11` violations without the `min(c_1,1)` guard).
Refuter 1 ran eight structural mutants, 8/8 fired, and flagged that the exponent mutant M1
fired only 5 times. I pushed on exactly that.

**P1 — the control is NOT seed-fragile** (`scripts/r1_fl043_power.py`). 40 seeds, `N = 4e6`:
violations `min 6 / median 13 / max 22`, **zero seeds returned 0**. So refuter 1's implicit
worry is unfounded and the attempt's "the guard is load-bearing" stands. (Reported because it
is a check I ran that came out *for* the attempt.)

**P2 — power curve on the exponent.** Perturb `1/5 → 1/5+ε` inside
`Re_* = E_0^{2/5}/(ν s^{1/5})`, seed 20260906, `N = 4e6`:

| ε | −0.05 | −0.01 | −0.001 | −0.0001 | 0 | +0.01 | +0.03 |
|---|---|---|---|---|---|---|---|
| violations | 288 | **1** | **0** | **0** | 0 | 0 | 4 |

A 5 % error in the load-bearing exponent is caught by **one draw in four million**; a 0.5 %
error by none.

**P3 — the gate cannot see the constant at all.** Replacing `c_1` by `c_1/2` or by `c_1/10^6`
(true but weaker) returns `0 violations` — PASS. The gate certifies the *direction* of (c),
nothing about its constant.

**P4 — the ensemble is not near the boundary.** Of the 2,053,207 draws with (b) true, only
**148 (0.0072 %)** have `M/m ≤ 1.01`; median `M/m = 8.4e6`, geometric mean `4.7e7`.

**P5 — the decisive one: an explicit witness the gate provably cannot see**
(`scripts/r3_exact_gate_and_accumulation.py`). Because `f(M) = c_1/(s(1+log_+(K M^{1/5})))` is
**strictly decreasing** in `M` (block (1), symbolic), (b) is *equivalent* to `M ≥ M_fix`, the
root of `M = f(M)`. That converts the whole question into a 3-parameter decision problem with
no sampling. Under `M = μ/s`, `κ = K s^{-1/5}` it becomes 2-parameter. Then:

> take `K = s^p` exactly. Then `K s^{-p} = 1`, the mutant's `log_+` switches off and its
> threshold is `m_p = c_1/s`, its largest possible value, while `κ = s^{p-1/5} → ∞` makes
> `M_fix` fall strictly below `c_1/s`. Any `M ∈ [M_fix, c_1/s)` satisfies (b) and violates the
> mutant. **So every mutant `p < 1/5` is false**, for a reason the uniform sampler cannot reach.

Witness for `p = 0.199`, printed by the script and re-checked in the original variables:

    s = T-t = 1e-4000 ,  K = E_0^{2/5}/ν = 1e-796 ,  c_1 = 0.5 ,  M = 0.051981/s * (1+1e-9)
    (b):  M s (1+log_+(K M^{1/5})) = 0.500000001 >= c_1 = 0.5      TRUE
    mutant threshold  m_p = c_1/s = 0.5/s ,  M s = 0.051981 < 0.5   VIOLATED

The witness needs `log10(1/s) ≳ 4/(1/5-p)`, i.e. **4000 decades** at `p = 0.199`, and
`log10 K ≈ −796`. The attempt's sampler draws `log10(1/s) ≤ 12` and `log10 K ∈ [−7.2, 11.2]`.
Nothing in the box can find it. Empirically (block (3), `N = 2e6`, 20 seeds):

| mutant | zeros / 20 seeds |
|---|---|
| `p = 0.19` (false) | **8 / 20** |
| `p = 0.199` (false) | **20 / 20** |

**FL-043 answer: YES for the exponent family — the gate returns PASS for the truth and for
false near-neighbours alike, and its verdict on `p = 0.19` is decided by the seed.** The
correct gate is the exact one above; block (1) of `r3` re-verifies (c) itself on a
`κ × c_1` grid with no sampling and no seed.

This does not damage the corollary — §3.3 is proved analytically, and the exact procedure
confirms (c) — but the note's evidential claim ("4,000,000 samples … the §3.3 step survives")
is worth far less than it reads, and the same sampler would have certified a false statement.

## 3. Form (c) is not "the explicit form"; it is a bounded weakening of one (MINOR)

`scripts/r2_sharp_explicit_form.py`. Since (b) `⟺ M ≥ M_fix`, the **sharp** explicit
consequence of (b) is `M ≥ M_fix`, not the attempt's
`m = c_1/(s(1+log_+(K s^{-1/5})))`. On a 4×4×3 grid in `(K, s, c_1)` the loss is

    max M_fix/m = 2.086984          (never an order; (c) is correct and lossy)

One line of the same bootstrap recovers most of it — feed `M ≥ c_1/s` into the logarithm
instead of `M ≥ 1/s`:

> **(c′)**  `‖ω(t)‖_∞ ≥ c_1 / ( (T-t)(1 + [ log K + (1/5)log(1/(T-t)) + (1/5)log c_1 ]_+ ) )`

verified `≤ M_fix` on every grid point. Effective coefficient `κ(s) = M s log(1/s)/c_1`:

| `T-t` | 1e-6 | 1e-12 | 1e-30 | 1e-100 | 1e-10000 |
|---|---|---|---|---|---|
| (c) as written | 3.671 | 4.234 | 4.663 | 4.894 | 4.9989 |
| (c′) | 4.183 | 4.555 | 4.812 | 4.942 | 4.9994 |
| sharp `M_fix` | 4.489 | 4.833 | 4.995 | 5.025 | 5.0013 |

All three tend to 5, so the headline `(5c_1+o(1))/((T-t)log(1/(T-t)))` is right for each. The
attempt's own form is the worst of the three at every reachable scale.

## 4. The novelty verdict is wrong at source (MAJOR)

### N1. A published NS lower bound at rate `(T-t)^{-1}`, log-free, every `t`, absolute constant

**McCormick, Olson, Robinson, Rodrigo, Vidal-López & Zhou**, *Lower bounds on blowing-up
solutions of the 3D Navier–Stokes equations in `Ḣ^{3/2}`, `Ḣ^{5/2}` and `Ḃ^{5/2}_{2,1}`*,
arXiv:1503.04323 (14 Mar 2015). Verified at source in `pdf/`, `txt/` (page numbers by
per-page `pdftotext`), and against the arXiv abstract page:

* **Theorem 4.1, PDF p.8** — *"Suppose that `u` is a classical solution of the Navier–Stokes
  … with maximal existence time `T`. Then `‖u(t)‖_{Ḃ^{5/2}_{2,1}} ≥ c/(T − t)`."*
* **Theorem 3.1, PDF p.6** — `‖u(T−t)‖²_{Ḣ^{3/2}} ≥ c_{3/2}^{-2} t^{-1}`.
* **Theorem 3.2, PDF p.7** — `limsup (T−t)‖u(t)‖_{Ḣ^{5/2}} ≥ c` (limsup only; the paper is
  explicit that this is the *weak* form and that the strong `Ḣ^{5/2}` bound is open).
* Abstract, p.1: *"we prove the strong lower bounds `‖u(t)‖_{Ḣ^{3/2}} ≥ c(T−t)^{-1/2}` and
  `‖u(t)‖_{Ḃ^{5/2}_{2,1}} ≥ c(T−t)^{-1}`."*

`Ḃ^{5/2}_{2,1}` for the velocity is `Ḃ^{3/2}_{2,1}` for the vorticity, and **has exactly the
scaling of `‖ω‖_∞`** (both scale as `λ²` under `u ↦ λu(λx)`). So Theorem 4.1 is: **Navier–Stokes,
every `t < T`, rate `(T−t)^{-1}` with no logarithm, constant absolute (no `E_0`, no data — it
is the constant of the Besov nonlinear estimate), in the vorticity's own scaling class.**

Against the attempt's §7 list of row 0's four advantages over row 2 (Euler / limsup /
non-invariant constant), Theorem 4.1 matches or beats three: it is NS, it holds at every `t`,
its constant is absolute, and its rate is *better* than row 0's by the whole logarithm. This
is strictly stronger evidence than the Cortissoz–Montero–Pinilla row refuter 1 raised (whose
constant does depend on `‖u_0‖_{L^2}` and which carries the log) — and it was in the same file.

**It still does not subsume row 0, and I checked why.** `‖ω‖_{Ḃ^{3/2}_{2,1}}` is an `L²`-based
norm with 3/2 derivatives; it is *not* controlled by `(E, ‖ω‖_∞)` in either direction
(`Ḃ^{3/2}_{2,1} ↪ L^∞` runs the wrong way for transferring a lower bound). Same for
`Ḣ^{3/2}` and for Robinson–Sadowski–Silva's `Ḣ^s` family. So none of the published NS lower
bounds transfers to `‖ω(t)‖_{L^∞}`, and row 0 transfers to none of them.

### N2. Refuter 1's open gap 3, closed

Refuter 1 read Cortissoz–Montero–Pinilla (2014) and Robinson–Sadowski–Silva (2012) only
through McCormick et al.'s quotation and asked for them at source. Independent confirmation:

* **Robinson, Sadowski & Silva**, *Lower bounds on blow up solutions of the three-dimensional
  Navier–Stokes equations in homogeneous Sobolev spaces*, **J. Math. Phys. 53 (2012) 115618**
  (AIP; Warwick WRAP record): `‖u(T−t)‖_{Ḣ^s} ≥ c_s t^{-(2s-1)/4}` for `1/2 < s < 5/2`,
  `s ≠ 3/2`. Real paper, real authors, **velocity `Ḣ^s`** — which settles refuter 1's N3 in
  refuter 1's favour and confirms the attempt's "the names are a fabrication" is wrong.
* **Cortissoz, Montero & Pinilla**, *On lower bounds for possible blow-up solutions to the
  periodic Navier–Stokes equation*, **J. Math. Phys. 55 (2014) 033101**: optimal-rate `Ḣ^s`
  lower bounds, `1/2 < s < 5/2`, with logarithmic corrections, constant depending on
  `‖u_0‖_{L^2}`.

Both are velocity Sobolev bounds; neither is a vorticity bound. The genre is exactly the
attempt's — and, per McCormick et al. p.2, is derived by exactly the attempt's method.

### N3. The attempt did not read its own input theorem's prior-art audit

`.../navier-stokes-pass8/PRIOR_ART_SCOPE.md` §2 already names BFG Thm 8/10 with page numbers
and the same caveat, and §3 names Kozono–Ogawa–Taniuchi (Kyushu J. Math. 57 (2003) 303–324)
**Theorem 1, p.306, with an explicit lifespan `T_* = C_ε‖a‖_{B^0_{∞,∞}}^{-2/(1-ε)}`** — a
lifespan bound, i.e. a member of the genre, not merely the rate-free continuation criterion
the attempt's §6.2 calls it. Dualised through `‖u‖_{B^0_{∞,∞}} ≲ ‖u‖_∞ ≤ CE^{1/5}M^{3/5}` it
lands back on the `5/6` baseline, so nothing changes — but the attempt ran eight fresh
WebSearches without opening the audit that shipped with its own input.

### What actually survives as novel

> A pointwise-in-`t` lower bound on **`‖ω(t)‖_{L^∞}`** for 3D Navier–Stokes with `ν > 0`,
> universal constant, all data-dependence inside a scaling-invariant Reynolds number.
> The quantity is what is unmatched — not the rate, not the pointwise-in-`t` property, not
> the constant-quality, all three of which McCormick et al. Thm 4.1 already has in print with
> a *better* rate. The single published statement in the same quantity is BFG ARMA 2019
> Thm 8, which subsumes it and whose proof is a self-described sketch with a located defect:
> **conflict unresolved**, as the attempt says.

## 5. §4's accumulation bound: correct, and its constant is unreachable

`r3` block (4). By substitution `r = e^{-v}` the integral is elementary and I evaluated it by
`mpmath` quadrature at 60 digits against the closed form: agreement to 9 decimals at every
`V ∈ {10, 10², 10³, 10⁶, 10⁸}` (e.g. `V=1e8`: quadrature `8.405621441`, closed form
`8.405621441`). So `∫_0^t‖ω‖_∞ ≥ 5c_1 loglog(1/(T-t)) − C(E_0,ν,T)` is exactly right.

The effective coefficient `I/(c_1 log V)`:

| `T-t` | 1e-6 | 1e-12 | 1e-100 | 1e-10000 | 1e-10^7 |
|---|---|---|---|---|---|
| coefficient | 2.52 | 2.83 | 3.54 | 4.20 | 4.53 |

It first reaches **4.5 at `T-t ≈ 10^{-4.24e6}`**. Refuter 1 flagged this qualitatively; the
solve is here. The `− C(E_0,ν,T)` carries it, so the statement is true, but "`5c_1`" is not an
approximation to anything that can be written down.

## 6. Minor

* **m1.** `s1` block (E)'s ratio table sets both unknown constants to 1 and then tabulates a
  *dimensional* ratio across four `(E_0,ν)` pairs. Only the fitted slope is a claim; the note
  says so. Fine, but the table invites the reading it disclaims.
* **m2.** The attempt's "SOURCE-INTEGRITY WARNING, propagate this" is half wrong in exactly
  the way refuter 1 says, and §4/N2 here confirms the names at source. It should be reissued
  as *"a garbled citation is a lead, not a phantom"* before it propagates.
* **m3.** `c_1` is called *universal*; the input theorem §1 calls it *dimensional*. Equal under
  the verified scaling invariance, but it is the input's word to change, not the corollary's.

## 7. Corrected statement that survives

> **Mathematics — unchanged, conditional on theorem (i).** For `ν>0`, `u_0 ∈ H^m`, `m ≥ 4`,
> divergence-free finite-energy, blowing up at `T = T_max < ∞`, for every `t ∈ [0,T)`:
> (a), (b), (c) exactly as written, and
> `∫_0^t‖ω‖_∞ ds ≥ 5c_1 loglog(1/(T−t)) − C(E_0,ν,T)`.
> Replace (c) by the sharper **(c′)** of §3, or by the exact `M ≥ M_fix`, at no cost.
> The proof is one line, not five: the input theorem's §5 concludes that *the actual solution*
> reaches `t + H(t)`, i.e. `T − t ≥ H(t)`, which is (a) rearranged.
>
> **Novelty — corrected, NOT NOVEL AS FRAMED.** The rate shape, the pointwise-in-`t` property
> and the absolute constant are all already in print for 3D Navier–Stokes, with a *better*
> (log-free) rate, in a space of the same scaling: McCormick et al., arXiv:1503.04323,
> Thm 4.1, p.8. The method is the standard one, named as such in that paper's introduction.
> What is unlocated in print is the **quantity** — a pointwise `‖ω(t)‖_{L^∞}` lower bound for
> NS — with the standing exception of BFG ARMA 2019 Thm 8, which subsumes it and is disputed.
>
> **Numerics — downgraded.** Block (B) is not evidence about the exponent. Report it as a
> direction check only, or replace it with the exact procedure of `r3` block (1).

## 8. Remaining gap

1. Theorem (i) is still not re-verified — by the attempt, by refuter 1, or by me. Everything
   above is conditional on it, and its §5 continuation is the single step the corollary leans on.
2. **BFG Thm 8/10 is still unresolved**, and it is the whole novelty question. Nothing here
   moves it. Note the stake: if it stands, the corollary is subsumed *and* the campaign's
   conjecture (ii) is false; the two hang together.
3. I did not pull the CMP (2014) or RSS (2012) PDFs themselves — I confirmed them from the
   AIP/WRAP records and McCormick et al.'s verbatim quotation. Better than refuter 1's
   single-source position, still not the PDFs.
4. Refuter 1's N2 (an NS transcription of Ingimarson–Kukavica Thm 2.1) is checked
   ingredient-by-ingredient by refuter 1 and not written out by anyone. Unchanged.
5. Sharpness — whether the logarithm in (a) is necessary — is conjecture (ii) and is untouched.
6. My §2 P5 witness is exhibited at 60-digit precision at one point of a family; I did not
   prove `M_fix < c_1/s` for all `κ>1` symbolically, only evaluated it.

## 9. Files

    NOTE.md                                     this note
    SHA256SUMS                                  recomputed, never typed
    scripts/r1_fl043_power.py   r1_log.txt      FL-043: seed robustness, power curve, tightness
    scripts/r2_sharp_explicit_form.py r2_log.txt  M_fix, the loss in (c), the sharpening (c')
    scripts/r3_exact_gate_and_accumulation.py r3_log.txt
                                                exact decision procedure; the p=0.199 witness;
                                                accumulation quadrature at 60 digits
    rerun/                                      copy of the attempt's scripts + my_s1.txt, my_s2.txt
    pdf/mcc-1503.04323.pdf  txt/mcc-1503.04323.txt   the N1 source (copied from refuter 1's
                                                download; page numbers re-located here, and the
                                                statements re-checked against arxiv.org/abs/1503.04323)

# NOTE — REFUTER of `corollary-blowup-rate` — DTC-2026-09-06 / lower

Target: `../corollary-blowup-rate/` (NOTE.md, scripts, pdf/, txt/).
Method: verified the attempt's SHA256SUMS, copied `scripts/` to `rerun/scripts/` and re-ran
both scripts from the copy; then attacked (1) the five-line continuation argument against the
input theorem's own text, (2) the numerics under FL-043 with a mutant battery, (3) the
novelty sweep at source, including a fresh literature pass.

## 0. Verdict

**REFUTED — MAJOR, on the novelty verdict only. The proof survives intact.**

The corollary (a),(b),(c) and the §4 loglog accumulation bound are correct, conditional on
theorem (i), and I could not break any step. What fails is the second half of the task: the
literature comparison. Three defects, one of which the attempt's own downloaded files
already contained the refutation of.

## 1. What survives (I tried to break it and could not)

**Reproduction.** `shasum -a 256 -c SHA256SUMS` in the attempt's folder: 14/14 OK.
Re-running both scripts from my copy reproduces `s1_log.txt` and `s2_log.txt`
**bit-for-bit** (`diff` empty; `rerun/scripts/my_s1.txt`, `my_s2.txt`).

**Line 2 is licensed at source.** I read the input theorem
(`~/<core-line-working-tree>/2026-09-06/new-chat/outputs/navier-stokes-pass8/LOGARITHMIC_RECORD_CLOCK.md`).
Its §1 says verbatim: *"There is a dimensional constant c>0 such that the solution
**continues through the interval of length** H_0 = c/(M_0(1+log_+ R_0))"*, and §5 ("First exit
and continuation") supplies that continuation via the Giga–Inui–Matsui restart, ending
*"Thus the actual solution reaches H and (2)–(3) follow."* So the corollary's Line 2 is not
an inference the attempt smuggled in; it is the input theorem's own conclusion. §5 also
confirms `0<c<1` ("Let A=1+log_+R_0 and H=c/(M_0A) with 0<c<1 … Choose a universal c small
enough"), so the `min(c_1,1)` guard of §3.3 is legitimately dropped. Lines 0,1,3,4,5 are
standard and correct; the argument uses only the continuation half and never refers to T
inside H(t), as claimed.

**§3.3 bootstrap.** Checked by hand: `M<m ≤ c_1/s < 1/s ⟹ Re_E < c_1^{1/5}Re_* ≤ Re_*
⟹ log_+Re_E ≤ log_+Re_*`, fed back into (b). Valid.

**§4 accumulation.** I re-derived the constant by **numerical quadrature** in `v=log(1/(T-t))`
(`scripts/r3_accumulation_quadrature.py`), not by the attempt's symbolic antiderivative.
The exact value is `5c_1·log(1+log K+V/5) − 5c_1·log(1+log K)`, i.e.
`5c_1·loglog(1/(T-t)) − C(E_0,ν,T)` exactly as stated; for K=1, V=1e8 my quadrature gives
84.0562 against the closed form's 84.0562. The `log_+`-inactive branch is correctly excluded:
where `log_+` activates, the bracket is `1+(log K+v/5) > 1 > 0`, so `F` is defined there.
**Caveat, not an error:** the convergence is `1/log V`. At `(T-t)=1e-100` the effective
coefficient is **3.54**, not 5. The "− C(E_0,ν,T)" carries that, so the statement is right,
but "(5c_1+o(1))" reads stronger than it is at any reachable scale.

**FL-043 on block (B): answered NO — the gate is not blind.**
`scripts/r1_fl043_mutants.py` re-runs the attempt's exact sampler (same seed, same 4e6 draws,
`(b)` true on 2,053,207 = 51.33%) against eight mutants of the explicit form (c):

| mutant | violations | |
|---|---|---|
| (c) TRUE | 0 | PASS |
| M1 exponent 1/4 in Re_* | 5 | fires |
| M2 exponent 1/6 in Re_* | 91 | fires |
| M3 E_0 exponent 1/5 | 2,277 | fires |
| M4 no `+1` in the bracket | 581,986 | fires |
| M5 no `min(c_1,1)` (attempt's own control) | 11 | fires |
| M6 numerator `5cc` | 92,688 | fires |
| M7 `log` for `log_+` | 2,060 | fires |
| M8 ν multiplied not divided | 80,068 | fires |

8/8 fire, the truth passes. **But note M1: 5 hits in 4,000,000 (1.3e-6).** The gate's power is
carried by structural mutations; against the exponent class that actually matters it is a
hair from blind. The "4,000,000 samples, 0 violations" headline overstates the test.

**§5's classical baseline 5/6 is right, and better-founded than the attempt argued.**
The attempt asserted `q=∞` is the best Leray exponent without checking.
`scripts/r2_independent_checks.py` block (1): the induced exponent is `5(q−3)/(6(q−2))`,
derivative `5/(6(q−2)²) > 0`, so it is increasing and `sup = 5/6` at `q=∞`. Confirmed.
Block (2) reaches the same baseline by a **second, disjoint route with no Leray and no Giga**:
the naive Duhamel lifespan for the vorticity mild equation with `‖u‖_∞ ≤ CE^{1/5}M^{3/5}` gives
`t_* = ν/(E^{2/5}M^{6/5})`, hence `M·t_* − 1/Re_E = 0` **exactly**, and pushing `c/(M·Re_E)`
through the same continuation argument returns `M ≥ ν^{5/6}E^{-1/3}s^{-5/6}` — row 5' again.
So the number is confirmed twice; only the word **"the only classical route"** stays unsupported.

## 2. What is refuted — the novelty verdict (MAJOR)

### N1. The "closest true analogue" is not Euler. It is a 2014 Navier–Stokes theorem with the same rate shape.

**Cortissoz, Montero & Pinilla, J. Math. Phys. 55, 033101 (2014)** proved, for 3D
**Navier–Stokes**, for every `t < T`:

        ||u(t)||_{Hdot^{3/2}} >= c / sqrt((T-t)|log(T-t)|)
        ||u(t)||_{Hdot^{5/2}} >= c / ((T-t)|log(T-t)|)          <-- row 0's exact rate shape
        "where in both cases c depends on ||u_0||_{L^2}."

Quoted verbatim from **McCormick–Olson–Robinson–Rodrigo–Vidal-López–Zhou, arXiv:1503.04323**
(SIAM J. Math. Anal. 48 (2016)), lines 95–101 of my extraction `txt/mcc-1503.04323.txt`;
independently confirmed by a second search. That is: **Navier–Stokes, not Euler; pointwise in
`t`, not `limsup`; constant depending only on the energy** — three of the four properties the
attempt's §7 lists as row 0's advantages over row 2. It is a different, incomparable quantity
(`‖u‖_{Ḣ^{5/2}} = ‖ω‖_{Ḣ^{3/2}}`, and `Ḣ^{3/2} ⊄ L^∞ ⊄ Ḣ^{3/2}`), so it neither implies nor is
implied by row 0 — but it is unambiguously closer than an Euler `limsup`, and the sweep's
statement "the closest true analogue in the literature is **Euler**" is false.

The attempt had this in its hands. Its own `txt/rss-1503.03063.txt` cites both
**[5] Cortissoz–Montero–Pinilla (2014)** (line 458) and **[11] Robinson–Sadowski–Silva,
J. Math. Phys. 53 (2012) 115618** (line 474), and its **line 44** contains the string
`(T − t) |log (T − t)|` in a displayed Navier–Stokes lower bound. The sweep grepped that file
for `vortic|omega|curl`, got 0, and stopped.

### N2. "First viscous counterpart of row 3" is unsupported: the Ingimarson–Kukavica proof is viscosity-blind.

I read their §3 in full (`../corollary-blowup-rate/txt/ik-2603.17431.txt`, lines 149–260). The
proof of Thm 2.1 has exactly four ingredients, and each survives adding `νΔ` with a
favourable sign:

1. `(3.2)  d/dt‖ω‖_∞ ≲ ‖Du‖_∞‖ω‖_∞` "following Lagrangian trajectories" — for NS the parabolic
   maximum principle (equivalently the positivity-preserving unit-mass kernel `U(t,s)` of the
   input theorem's own eq. (15)) gives the same bound on the upper Dini derivative;
2. `(3.3)` Kozono–Taniuchi `‖Du‖_∞ ≤ C(1+‖ω‖_∞(1+log_+‖u‖_{H^3}))` — a functional inequality,
   equation-free;
3. `d/dt‖u‖_{H^3} ≲ ‖Du‖_∞‖u‖_{H^3}` — for NS the extra term is `−2ν‖∇u‖²_{H^3} ≤ 0`;
4. everything after (3.6) is ODE manipulation (reparametrise by `τ=∫‖ω‖_∞`, Gronwall, two logs).

BKM holds for NS. So an NS version of their (2.1) and (2.2) follows by transcription of a
published proof. Row 3' is therefore not "the first viscous counterpart of row 3", and §7's
"(i) row 2 is Euler, row 0 is NS, so neither implies the other" conceals that the competitor
row 0 must actually beat is the NS transcription, against which only two of the four listed
advantages survive: **every `t` instead of `limsup`**, and **a universal constant instead of
one depending on `‖u_0‖_{H^3}`**. (Advantage (iii), "a factor 5 in the constant", compares
`5c_1` with `1/C(‖u_0‖_{H^3})` — incomparable unknowns; it is not an advantage.)
Registered honestly: I checked each ingredient survives, I did not write the transcription out.

### N3. The "SOURCE-INTEGRITY WARNING, propagate this" is itself half wrong, and would propagate an error.

The attempt says of the machine summary of arXiv:1503.03063: *"Theorem 1.1 (Robinson,
Sadowski, Silva) … That is a fabrication."* The **arXiv number and the vorticity claim** are
indeed wrong. The **names are not invented**: Robinson, Sadowski & Silva are the authors of
*"Lower bounds on blow up solutions of the three-dimensional Navier–Stokes equations in
homogeneous Sobolev spaces"*, **J. Math. Phys. 53 (2012) 115618**, the founding paper of
exactly this genre — cited as reference **[11] at line 474 of the attempt's own extraction**,
and named in Cortissoz–Montero's opening paragraph ("based on ideas presented by Robinson,
Sadowski and Silva in [11]"). The summary **conflated a paper with its cited predecessor**.
The attempt's grep established only that `omega` is absent from the *wrong* paper; it never
checked whether the named authors exist. The correct lesson is *"a garbled citation is a lead,
not a phantom — chase it"*, and chasing it lands on N1. Propagating "the names are a
fabrication" would tell the estate to discard a real, directly on-point reference.

## 3. Minor findings

* **m1.** `s2` block 2's "three exact zeros" is **two** distinct checks: the dict entries for
  forms (a) and (b) are the *identical* sympy expression. All three use plain `log`, not
  `log_+`; and since each factor (`M(T-t)`, `Re_E`) is separately invariant, the test can only
  fail on a mistyped exponent. Non-vacuous, but much weaker than "verified symbolically,
  three exact zeros" suggests.
* **m2.** Transferring BFG Thm 8 to a lower bound needs an unstated step: their solution is
  mild in `C_w([0,T],L^∞)`, and identifying it with the maximal `H^m` strong solution (so that
  `T_max ≥ T`) is asserted, not argued. This slightly *weakens* row 1's threat — a point the
  attempt declined to take in its own favour.
* **m3.** IK's own proof text says the constant depends on `‖u_0‖_{H^3}` **and `T_*`** (their
  (3.5) and the line after); their theorem statement drops `T_*`. The note follows the
  statement. Faithful, incomplete.
* **m4.** The input theorem calls `c` a *"dimensional constant"* (§1); the corollary calls
  `c_1` *"universal"*. Given full scaling invariance these agree, but the note should say so
  rather than silently upgrade the word.
* **m5.** McCormick et al.'s introduction states that Leray's bound *"and all subsequent lower
  bounds, are a consequence of upper bounds on the local existence time"*. That is the
  corollary's own five-line argument, named in the literature as the standard route. The note
  claims no methodological novelty, but the novelty verdict would be better calibrated by it.

## 4. Corrected statement that survives

> **Proved (conditional on theorem (i), whose §1 explicitly supplies the continuation):**
> forms (a), (b), (c) and `∫_0^t‖ω‖_∞ ds ≥ 5c_1 loglog(1/(T-t)) − C(E_0,ν,T)`, exactly as
> written. Nothing in §2–§5 needs repair.
>
> **Novelty verdict, corrected — NOT NEW AS FRAMED; A NARROWER NOVELTY SURVIVES.**
> The corollary is the standard dual of a lifespan theorem, so its content is entirely its
> input's. Against the literature:
> (i) it is **implied by** BFG ARMA 2019 Thm 8 if that stands — conflict unresolved, as the
>     attempt says;
> (ii) it is **neither implied by nor implies** Cortissoz–Montero–Pinilla, J. Math. Phys. 55
>     (2014) 033101, which already gives, for **Navier–Stokes**, `‖u(t)‖_{Ḣ^{5/2}} ≥
>     c(‖u_0‖_2)/((T-t)|log(T-t)|)` for **every** `t<T` — the same rate shape, same equation,
>     energy-only constant, in an incomparable norm;
> (iii) an NS version of Ingimarson–Kukavica Thm 2.1/(2.1) follows from their **published
>     Euler proof by transcription**, so §4 here is not the first viscous counterpart.
>
> What remains genuinely unlocated in print is narrower and should be claimed as exactly this:
> *a pointwise-in-time lower bound on `‖ω(t)‖_{L^∞}` for 3D Navier–Stokes whose constant is
> universal and whose entire data-dependence sits inside a scaling-invariant Reynolds number* —
> with the standing exception of BFG Thm 8, which subsumes it.

## 5. Remaining gap

1. The whole novelty question still hangs on **BFG Thm 8/10**. Nothing here moves it.
2. The NS transcription of IK Thm 2.1 is **checked ingredient-by-ingredient, not written out**;
   somebody should write it and compare its constant to `5c_1`.
3. I read Cortissoz–Montero–Pinilla (2014) and Robinson–Sadowski–Silva (2012) **only through
   McCormick et al.'s verbatim quotation and Cortissoz–Montero's reference list**, not from
   their own PDFs. Both should be pulled at source before the corrected verdict is banked.
4. Theorem (i) itself is not re-verified here either (I read its statement and §5, not its
   §2–§4 kernel estimates).
5. Sharpness (conjecture (ii)) is untouched, as the attempt says.

## 6. Files

    NOTE.md                                   this note
    SHA256SUMS
    scripts/r1_fl043_mutants.py     r1_log.txt   FL-043 mutant battery on the attempt's block (B)
    scripts/r2_independent_checks.py r2_log.txt  Leray-family optimality + disjoint route to 5/6
    scripts/r3_accumulation_quadrature.py r3_log.txt  loglog constant by quadrature, branch check
    rerun/scripts/                            copy of the attempt's scripts + my_s1.txt, my_s2.txt
    pdf/mcc-1503.04323.pdf  txt/mcc-1503.04323.txt   McCormick et al., the source for N1

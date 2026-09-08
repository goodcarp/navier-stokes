# REFUTER — gap-V-aronson

Seat `gaps/refute-gap-V-aronson`, DTC-2026-09-06. Laws `TORMENT NEXUS/LAWS.md`.
Target: `gaps/gap-V-aronson/` (NOTE.md + v1…v5 + check_constants.py + SHA256SUMS + sources/).
Every number below came out of `r1`–`r3` in this folder, which I wrote and ran; scripts were
re-run from a copy at
`/private/tmp/.../scratchpad/work/`. Nothing outside this folder was written.

---

## 0. Verdict

**REFUTED — MAJOR.** The mathematics survives; the evidence apparatus and the application do not.

| the attempt's claim | my finding |
|---|---|
| **Theorem V.1** (`N·min(B1,B2)`, constants in `c,q,n` only) | **STANDS.** Proof checked step by step (Feynman–Kac; MVT propagator; Grönwall; exponential supermartingale; ε-net). `B1 = 2n e^{−e^{−2c}q/4n}` and `B2` re-derived independently and both algebraic reductions confirmed. |
| **Lemma V.0** (Lagrangian operator is exactly `νD^{-1}∂_l(DG^{lk}∂_k)`, no drift, no zeroth order) | **STANDS.** `v1_reduction.py` re-run: Piola residuals `[0,0,0]`, reduction residual exactly `0` at all 8 rational points, `‖∇₅b‖_op−‖∇₃u‖_op = 4.441e-15` over 400 jets. |
| **"the divergence `2a`" term is empty** | **STANDS.** |
| **Route (1) (Aronson / Fabes–Stroock / Norris–Stroock) does not close** | **Conclusion stands — but the seat did not establish it. §1 below.** |
| **"REFUTER RUNS … R1 not violated, R3 control fires"** | **FAILS. R1 has no failure mode at all; R3 is forced by construction. §2.** |
| **§5 application (`f`-table, `c₂` inflation, `L_min = 590/854`)** | **FAILS — computed at the wrong window constant. §3.** |
| **headline "The estimate the brief asks for is TRUE, it is PROVED here"** | **OVERCLAIM.** The estimate `prove-lagrangian` §4(3) asks for is V-b, which the same note concedes is open. §4. |
| **`Γ = ‖∇u‖_∞ = 2a(0,t)`** | **UNPROVEN HYPOTHESIS**, nowhere stated as one. §5. |

Nothing here overturns Theorem V.1. Everything here bears on whether the seat's deliverable can
be signed as it stands. It cannot.

---

## 1. The source chain for §3 is broken (`r3`, part a)

The note says (§3): *"Sources in `sources/` (Aronson 1968 full text and page images;
Norris–Stroock 1991 PDF; **Aronson's own 2017 survey fetched this run**)"*, lists the survey again
in §7, and quotes it **three times verbatim** — those three quotes are the entire basis for the
dismissals of Fabes–Stroock (§3b), Norris–Stroock (§3c) and Porper–Eidel'man (§3d), i.e. for the
verdict "Route (1) DOES NOT APPLY". What is actually on disk:

| file | bytes | what it really is |
|---|---|---|
| `sources/aronson1968.pdf` | 5 556 660 | real PDF ✓ |
| `sources/aronson1968.txt` | 131 825 | OCR text ✓ |
| `sources/norris-stroock-1991.pdf` | 8 349 962 | real PDF ✓ (never quoted — see below) |
| `sources/gh.pdf` | 7 715 | **HTML: arXiv "No document for '1707.04620v2'"** |
| `sources/gh2.pdf` | 126 | **HTML 404 page** |
| `sources/zhang1997.pdf` | 3 038 | **HTML "Client Challenge" interstitial** |

The survey was never retrieved. `SHA256SUMS` certifies all three of those HTML pages as sources and
`shasum -a 256 -c` reports 36/36 OK — **a manifest that hashes a 404 and reports clean is a control
that cannot fail.**

I retrieved the paper myself (`arXiv:1707.04620v1`, 5 pp., saved here as
`aronson_survey_1707.04620v1.pdf/.txt`) and checked the three quotes:

* **Norris–Stroock** — *"However they are forced to assume the uniform continuity of A and E − Ê."*
  **verbatim ✓** (survey, final paragraph).
* **Porper–Eidel'man** — *"a slight generalization of equations (1) and (3) involving a coefficient
  `p(x)` multiplying `∂_t u`"* **verbatim ✓**.
* **Fabes–Stroock** — the note prints *"Neither Nash **nor** Fabes & Stroock consider the full
  equation (1)"*. The source reads *"Neither Nash **or** Fabes & Stroock consider the full equation
  (1)."* **Not verbatim** (substance unaffected).

So the §3 verdict is **correct** — I have now checked it at source — but the seat asserted it from a
document it did not have, and one of its three "verbatim" quotes is altered. Two further gaps in
that section:

* **Norris–Stroock is never quoted from `norris-stroock-1991.pdf`**, which *is* on disk. The
  hypothesis "uniform continuity of `A`" is taken from Aronson's one-sentence description of them.
* **"their weight is `p(x)`, time-independent"** (§3d, the step the note calls "the gap that
  actually kills Route (1)") is an inference from the survey's *notation*. The survey does not say
  the weight is time-independent, Porper–Eidel'man 1984 is not on disk, and the same sentence says
  they treat "the analogue of equation (1)" — the full drift equation — "under conditions similar
  to (H)". The kill-step is unverified at its own source.

The Aronson 1968 quotes *are* at source in `aronson1968.txt` (OCR): p. 608 "Throughout the paper it
will be assumed…"; p. 624 "By the structure of L we mean n and the quantities which occur in the
hypotheses (H). In particular, α depends only on n, M and ν, while β depends only on …, M₀ and ν";
Theorem 8 at line 2833. The note silently cleans up the OCR without saying so, but the substance
matches. The core Aronson finding — constants depend on `M₀ = ‖b‖` and on `T` — is confirmed.

## 2. The falsification apparatus is vacuous (`r3`, part b)

The note's pre-declared rules are **R1** "refuted if `P_emp − 3·s.e. > min(B1,B2)`" (TEST A, the
`n=5` Monte Carlo) and **R3** "non-diagnostic unless the Eulerian control violates the bound".

* **R1 cannot fire.** `P_emp` is a probability, so `P_emp ≤ 1`. The five TEST A rows have
  `min(B1,B2) = 8.642, 7.200, 5.577, 4.015, 2.688` — **every bound exceeds 1**. R1 is satisfied
  identically, for every conceivable dataset, at every `q` tested, for any drift whatsoever.
  Reporting "R1 violated: False" over a vacuous bound is a pass with no information content. (The
  *scaling* observation in the same table — that the two viscosities agree at equal `q` — is real
  and does carry information; it is just not a test of the bound.)
* **R3 cannot fail.** The control reads the deviation at the fixed Eulerian point `x₀`, and the
  datum difference is `1{|x−x₀| ≥ d}` — a hole of radius `d ≤ 0.2530` around `x₀` — while the drift
  displaces `x₀` by `0.7100` over the window, **2.8× the hole's own radius**. The hole is advected
  clear of `x₀` by construction, so the control reads `1.000 = N` in every row necessarily. It
  demonstrates that the material frame matters; it does not test anything.
* **R2 is a genuine test** (bounds `0.892 / 0.325 / 0.079` against measurements `8.7e-3 … 4.7e-4`;
  tightest margin 11.5×) and it passed. It is the only one of the three that could have fired.

So the note's line "**Summary: R1 violated False, R2 violated False, R3 control fires True**"
overstates the evidence by two thirds.

## 3. §5 is computed at the wrong window constant (`r1`)

Theorem V.1 defines `Γ := sup_{t≤τ}‖∇b(·,t)‖_{L^∞,op}` and `c := Γτ`; every `e^{c}` in the proof
comes from Grönwall against that **constant** majorant. `v5_application.py` line 31 sets
`c_acc = theta_acc = 4(1−√(2/3)) = 0.7340137`, which is `Γτ` only if `Γ = ML`, i.e. only if the
strain never grows. But the accelerated window is *defined* by the strain growing. From
`prove-lagrangian`'s own closed form (4.1), at the tracked shell `σ = 0`,

```
a(θ) = M L κ H(0,θ) ,   H(0,θ) = 1/(1 − κθ/2) ,   κ = 1/2 ,   θ = M L t ,
‖∇₅b‖_op = 2a  (the note's own §1)   ⟹   Γ(θ) = M L · H(0,θ) ,
```
so `Γ` grows by `H(0,θ_acc) = √(3/2) = 1.2247449` across the window (the *stretch* `λ` reaches
`3/2`, but `a` is driven by the average over outer shells `H`, not by `λ`; `2∫a dt = 2log(3/2)`,
which `r1` asserts). Hence

| constant | value | what it is |
|---|---|---|
| `c` used by the note | `0.7340137` | `Γτ` with `Γ` frozen at its `t=0` value — **not** `sup Γ · τ` |
| `c = (sup Γ)·τ` | **`0.8989795`** | the theorem **as proved**, `= √(3/2)·θ_acc` |
| `C(τ) = ∫₀^τ Γ dt` | **`0.8109302` `= 2log(3/2)`** | the sharper constant available from a time-dependent Grönwall — **a repair the note does not make** |

(`ντ` is physical and unchanged, so `q = (d/ρ₀δ)² L/θ_acc` in all three columns; only the rate moves.)
Consequences (`r1`, verified against the note's own table, which my harness reproduces exactly in the
first column):

```
rate (−log bound / q)          B1        B2       best-rate loss vs the note
  as the note computes      0.011519  0.025310        —
  theorem as proved         0.008282  0.014780      1.712x
  time-dependent repair     0.009877  0.020243      1.250x

required inset f            L=10     L=20     L=40     L=80     L=160
  note                     0.9131   0.6590   0.4694   0.3382   0.2437
  theorem as proved        1.0829   0.7922   0.5673   0.4149   0.3031
  time-dependent repair    0.9886   0.7237   0.5185   0.3785   0.2727

c2 = 8(1-sqrt(2/3)) inflated by L/(L - log(1+f))
  note                     1.5699   1.5061   1.4823   1.4734   1.4700
  theorem as proved        1.5843   1.5121   1.4847   1.4744   1.4705
```
Also: the note's §1 sentence "`τ = 2θ/(ML)`, `c = Γτ = 2θ`" is internally inconsistent — it is
applied as `c = 2θ` for the frozen window (`θ = log 3/2`) and as `c = θ` for the accelerated one.
Two different `θ`'s under one symbol.

**What survives:** the inflation is `L/(L − log(1+f))` and `f` grows only like `1/√rate`, so the
inflation is `O(1/L)` for every constant in the sweep. **`c₂ = 8(1−√(2/3)) + o(1)` is untouched.**
The tables are wrong; the conclusion drawn from them is right.

## 4. The headline does not match the brief

`prove-lagrangian` §4(3) names the missing estimate as
`|η(x(t),t) − η_transported| ≤ ‖η₀‖_∞ exp(−d²/(Cντ))` — the viscous solution against its own
*inviscid transport*. Theorem V.1 bounds `η₁ − η₂`, two solutions of the **same** equation with
different data. That is V-a. It is not the brief's estimate. The note's own §0(2) and §5(iii) say so
("V-b … NOT CLOSED", `O(L^{-1/2})` unconditional), which directly contradicts its own opening
sentence, *"The estimate the brief asks for is TRUE, it is PROVED here."* Both cannot be filed.

A second, repairable, statement gap: §5 compares `η₀` with an unnamed `η̃₀`, and never says that
`η̃` must solve the **linear** equation with the *same* `u`. If `η̃` is read as a second Navier–Stokes
solution — the natural reading of "the datum with its edges removed" — then `b₁ ≠ b₂` (the velocity is
`u[η]` by Biot–Savart) and Theorem V.1 does not apply at all. The repair is one sentence; it is not
in the note.

## 5. `Γ` is an unstated hypothesis (`r2`)

Theorem V.1 needs the **global** spatial sup of `‖∇₅b‖_op` over `ℝ⁵ × [0,τ]`. §5 supplies the
**on-axis** value `2a(0,t)`, justified only by the identity `‖∇₅b‖_op = 2a` for a *uniform*
axisymmetric strain. Neither note bounds `‖∇u‖_{L^∞}` for the actual mollified-shell field:
`prove-lagrangian` L3 bounds the *velocity* deviation (`sup|a − A(ρ)| = 0.7047M`,
`sup|u^z + 2Az|/(Mρ) = 0.19143`) and a mode sum on `ψ` — no derivative of the deviation velocity.
`∇u` is a Calderón–Zygmund operator on `ω` and is not controlled by `‖ω‖_∞`.

Writing `Γ = g·2a(0,t)` and sweeping the unknown `g` (`r2`; `g = 0.8165` reproduces the note):

```
 g       c        L=10                 L=40                 L=160
 0.8165  0.73402  f= 0.907 c2=1.5694   f=0.469 c2=1.4823    f=0.241 c2=1.4700   <- the note
 1.0000  0.89898  f= 1.076 c2=1.5837   f=0.567 c2=1.4847    f=0.300 c2=1.4704
 1.5000  1.34847  f= 1.844 c2=1.6394   f=0.913 c2=1.4922    f=0.478 c2=1.4716
 2.0000  1.79796  f= 3.589 c2=1.7319   f=1.458 c2=1.5018    f=0.761 c2=1.4732
 4.0000  3.59592  f=29.639 c2=2.2318   f=14.85 c2=1.5770    f=7.176 c2=1.4876
```
An inset always exists and the `c₂` inflation stays `O(1/L)` even at `g = 4`, so **the asymptotic
claim is robust to `Γ`** — this is the one place where the note's headline conclusion is stronger
than its own argument warrants but still true. What is not robust is every displayed number: an
undetermined `O(1)` factor in `Γ` moves `f` at `L = 10` from `0.91` to `30`.

## 6. Smaller items

1. *"sharp rate `cq/(2(e^{2c}−1))`, **attained by the pure axisymmetric strain**"* — not shown. The
   TEST A field is not a pure strain, and the measured rate there (`0.11503`) **exceeds** the
   "sharp" value (`0.09493`), so that run is evidence the quantity is an upper bound, not that it is
   attained.
2. `scripts/` is described as "an **earlier, interrupted pass of this same seat** … kept as
   **independent** corroboration". Same seat, same session lineage — it is a re-run, not
   decorrelation, and two of §2.4/§5(iii)'s constants (`C' = 3.08`, first moment `≈ ‖∇²b‖ντ²`) rest
   on it alone.
3. `check_constants.py` prints `FAILURES: none` and does verify every displayed number against the
   results files — but `c_acc` is hard-coded in `v5_application.py`, so the check is a
   self-consistency test and cannot see §3's error. Confirmed by re-running it: it passes.
4. Everything in `v1`/`v2`/`v3` reproduces on re-run from a copy; the exact-arithmetic content of
   Lemma V.0 is sound.

## 7. Corrected statement

> **Theorem V.1 is correct as stated and proved.** `Lemma V.0` is exact. The Route (1) verdict is
> correct and is now checked at source *here* rather than there.
> For the accelerated doubling window the window constant is `c = (sup_t Γ)τ = √(3/2)·4(1−√(2/3)) =
> 0.8989795`, or `C(τ) = ∫₀^τ Γ dt = 2log(3/2) = 0.8109302` after the (unmade) time-dependent
> Grönwall repair — **not** `0.7340105`; the required insets are
> `f = 1.083, 0.792, 0.567, 0.415, 0.303` at `L = 10,20,40,80,160`, and `c₂ ≤ 1.5843` at `L = 10`,
> `1.4705` at `L = 160`.
> These rest on the unstated hypothesis `Γ := ‖∇u‖_{L^∞(ℝ⁵×[0,τ])} = 2a(0,t)`, which neither note
> proves; an undetermined `O(1)` factor there rescales every inset but leaves
> `c₂ = 8(1−√(2/3)) + o(1)` intact.
> **Gap V is not closed.** V-a is closed *for a prescribed common `b`*; V-b — the brief's actual
> estimate — is open at `O(L^{-1/2})`, and `O(1/L)` still needs a `∇²u` bound.
> The R1/R3 results carry no evidential weight and the §3 citations must be re-cited from the
> retrieved survey (`aronson_survey_1707.04620v1.pdf`, here) and from Porper–Eidel'man 1984 at
> source before the Route (1) verdict can be signed.

## 8. Files

| file | what it does |
|---|---|
| `r1_window_c.py` / `r1_results.json` | derives `Γ(θ)` from `prove-lagrangian` (4.1), computes the three window constants, the three rate laws, and re-derives §5(ii)'s inset table in each (first column reproduces the note exactly) |
| `r2_gamma_sensitivity.py` / `r2_results.json` | sweep in the unknown multiplier `g` of `Γ`; the `f` and `c₂` tables, and the largest `g` for which the `L=10` inset stays below `f = 5, 2, 1` |
| `r3_sources_and_controls.py` / `r3_results.json` | classifies every file in `gap-V-aronson/sources/`; checks the three §3 quotes against the retrieved survey; failure-mode audit of R1/R2/R3 |
| `aronson_survey_1707.04620v1.pdf` / `.txt` | the survey the target note quotes but does not have — retrieved here from arXiv |
| `SHA256SUMS` | computed over this folder |

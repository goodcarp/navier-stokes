# SYNTHESIS — DTC-2026-09-06 / lower

Seat: SYNTHESIZER. Inputs: four attempt seats (`prove-lagrangian`, `prove-duhamel`,
`depletion-numerics`, `corollary-blowup-rate`) and five refuter notes
(`refuter-lagrangian`, `refuter-lagrangian-b`, `refuter-prove-duhamel`,
`refuter-prove-duhamel-2`, `refuter-corollary-blowup-rate`, `refuter2-corollary-blowup-rate`).
Working files of this seat: `synth/s1_constants.py`, `synth/s2_gate_recheck.py` and their JSON.
Nothing outside `lower/` was written; nothing inside another seat's folder was modified.

This note states only what survives refutation. Where an attempt's headline and its own
status table disagree, the status table governs (three of the four headlines were written
stronger than the seat's own per-step status; the refuters caught all three).

FL-000 stands. Nothing in this round touches the headline problem.

---

## 1. STATUS OF CONJECTURE (ii): **SKETCH**

Not PROVED. Not PROVED-MODULO-a-named-lemma. **SKETCH** — two independent partial routes,
each terminating in an unproved hypothesis that is the analytic heart of the claim, plus
(on the Lagrangian route) two further unwritten steps.

### 1.1 What conjecture (ii) would say

> There exist c₂ < ∞ and Re_* such that T(Re) ≤ c₂/log Re for Re ≥ Re_*, realized by the
> mollified axisymmetric no-swirl δ-tapered plateau family: ω^θ_0 = −M sgn(z) h_δ(φ) on
> ρ₀ < |x| < R, mollified at scale ε·ρ, ρ₀ = δ^{-1}√(ν/M), R = ρ₀e^L.

**This is not established.** No seat proved an upper bound on T(Re) of any shape.

### 1.2 Route A (Lagrangian, `prove-lagrangian`) — what is PROVED

**LEMMA 1 (PROVED, exact; confirmed by five independent instruments).** For
η₀ = −M sgn(z)h(φ)/r on ρ₀<|x|<R and the coaxial volume-preserving strain
T_λ:(r,z)↦(λr, λ^{-2}z) with η transported,

    a_λ(0) = (M/2)·P_h(λ)·L ,  L = log(R/ρ₀),  P_h(λ) = 3∫₀¹ h(v)v²(Av²+B)^{-5/2}dv,
    A = λ²−λ^{-4}, B = λ^{-4};   and for h ≡ 1,  P_1(λ) = λ  EXACTLY.

Instruments agreeing: sympy antiderivative (residual identically 0); mpmath 30-dps (rel 0.0);
Lagrangian quadrature (1.4e-12); an independent Eulerian rebuild (1.1e-4 grid);
Gegenbauer projection of the deformed profile (1.9e-8); and the refuter's own 5D Gegenbauer
instrument (3.1e-7 to 6.9e-6 at λ = 1…1.5, R/ρ₀ = 65536).

**Taper deficit (PROVED).** D ≤ w³/(B(Aw²+B)^{3/2}), w = sin δ; hence Φ_{h_δ}(λ) ≥ 1 on
λ∈[1,3/2] for every δ ≤ 30°. Φ_{h_δ}(3/2) = 1.4980/1.4914/1.4736/1.4444/1.3597/1.2584/1.0674
at δ = 3/5/7.5/10/15/20/30°.

**LEMMA 2 (PROVED, unconditional, all ν ≥ 0).** The datum is odd in z; axisymmetric no-swirl
NS preserves that; {z=0} is a material plane on which η = 0; by the parabolic maximum principle
η ≤ 0 on {z>0} for all t. Hence the integrand of a(0) is pointwise non-negative and
a(0,t) ≥ 0 for all t. **Deflated by refuter:** this is the classical axisymmetric maximum
principle in new packaging, and it is never used by any step of the argument.

**Closed integro-ODE (PROVED *given* the model).** With F = ∫₀ᵗ a, H(σ,θ) = ∫_σ¹ e^F dσ′,
∂_θH = (κ/2)H², H(σ,0) = 1−σ, so λ(0,θ) = (1−κθ/2)^{-2} and λ = 3/2 at κθ = 2(1−√(2/3)).
Verified against a 4001-point method-of-lines solve to 2.4e-10 and against the defining
integral identity to 5e-14.

### 1.3 Route A — what is REFUTED (MAJOR)

**R-A1. The headline "step (2) is not needed; the perturbative step is replaced by an
identity" is FALSE**, and the seat's own status table says so. Lemma 1 is an identity about
the *hypothetical* field η₀∘T_λ^{-1}. Using it requires
 (i) **L3v** — the deformed field is still a per-shell pure strain + O(Mρ) (status: SKETCH); and
 (ii) **GAP T** — the true flow map equals T_{λ(σ)} to relative O(c/L) *and* a[η](0) is
 Lipschitz-stable under that perturbation with constant C·M·L (status: unproved).
The brief asked only for a **one-sided** lower bound a ≥ κML(1−O(c)). GAP T is a **two-sided
Lipschitz estimate for a Calderón–Zygmund-type kernel in log-radial coordinates**. The
perturbative burden was relocated and made strictly stronger, not removed.

**R-A2. Lemma 1 contains no cancellation** (deflation of the stated mechanism). Eulerian:
η₀ = −M sgn(z)h/r depends on r alone; T_λ scales r by one constant; the image shell is
{ρ₀ < |w|G(t) < R} with *both* radii scaled by 1/G(t), so the log-thickness is L in every
direction. a is linear in η with a scale-free kernel ⇒ a_λ(0) = λ·a₀(0). Deleting G(t) from
the integrand leaves the answer bit-identical (0.49999999955921703 / 0.6249999994490216 /
0.7499999993388258 / 2.4999999977960865 at λ = 1/1.25/1.5/5, with and without G). So the
advertised mechanism ("rotation toward the equator versus inward collapse near the axis is
exactly neutral") misdescribes the seat's own identity, and "the strain scales precisely like
the stretched vorticity" is circular inside the ansatz.

**R-A3. The accelerated constant c₂ = 8(1−√(2/3)) = 1.4680274 has no empirical support and is
contradicted at 5.3 σ by the only measurement of it.** Against viscous-numerics'
Q′ = T_{3/2}·M·log(R/s*) = 0.988 ± 0.018 with measured κ = 0.437/0.415/0.411 at N = 3/5/6
(this seat's own recomputation, `synth/s1_results.json`):

| model | Q′(κ=0.437) | z | Q′(κ=0.415) | z | Q′(κ=0.411) | z |
|---|---|---|---|---|---|---|
| frozen, θ=log(3/2) | 0.92784 | −3.34 | 0.97702 | −0.61 | **0.98653** | **−0.08** |
| accelerated, θ=2(1−√(2/3)) | 0.83983 | −8.23 | 0.88435 | −5.76 | 0.89296 | −5.28 |

Advertised gain 1 − 1.4680274/1.6218604 = **9.485 %**; the accelerated model's shortfall against
the data at the finest resolution = 1 − 0.892961/0.988 = **9.619 %**. The entire gain the
headline sells is the size of the model's disagreement with the only measurement of it.

**R-A4. FL-043 — the Route-A gate cannot fail on most of what it reports.** Verified directly
by this seat (`synth/s2_gate_recheck.py`):
 - one check line reads `ck('sum |H| ||C||/l^2 (lam=1)', s5[...]*0+0.3958, 0.3958, 1e-9)` — the
   stored value is multiplied by **zero**; the s5 key it names holds 0.011311110612413619,
   a 35× mismatch that passes;
 - the summary is `print(("ALL %d CHECKS PASS" % 0) …)` — it always prints **"ALL 0 CHECKS
   PASS"** (log shows 88 PASS lines, 0 FAIL, final line `ALL 0 CHECKS PASS`);
 - five further checks recompute the constant they are compared against, taking no stored input;
 - dynamic mutation (v ↦ 1.5v+0.37, one number at a time): 416 stored numbers, gate catches 95,
   **blind to 321 (77.2 %)** — including Φ′(1) = 1.0000000000192764, the sign of the whole
   acceleration claim.

**What Route A gained that is NEW and survives** (from the refuter, not the attempt):
at the tracked **material** point (ρ₀, φ₀ = 30°) moved by T_λ,
a(material) − (M/2)λL = −0.01808 / +0.02838 / +0.03733 / −0.01215 / −0.06857 at
λ = 1/1.1/1.25/1.4/1.5 (R/ρ₀ = 65536), stable to three digits at R/ρ₀ = 256 and 4096 — i.e.
the offset is **O(M) and L-independent, uniformly on λ∈[1,3/2]**. That closes, *numerically*,
one piece the attempt had folded into GAP T. It is measured, not proved (PIN 5).

### 1.4 Route B (Riccati, `prove-duhamel`) — what is PROVED

**Exact identities (PROVED, sympy, five residuals identically 0, re-derived from Cartesian 3D
NS with no quoting):**

    D_t a = −a² + β + νV ,  a = u^r/r ,  β = −(1/r)∂_r p ,
    V = (1/r)(∂_rr + (1/r)∂_r − 1/r² + ∂_zz)u^r ,  D_t η = ν Δ₅ η .

**THEOREM D (PROVED).** ν>0, smooth finite-energy axisymmetric no-swirl data,
M₀ = ‖ω₀‖_∞, X(t) the trajectory from x₀, a₀ = a(x₀,0) > 0. If on [0,τ] along X
**(H2)** β + νV ≥ 0 and **(H3)** νW ≥ −λ (W = Δ₅η/η), then

    a(X(t),t) ≥ a₀/(1+a₀t)  and  |ω^θ(X(t),t)| ≥ |ω₀(x₀)|(1+a₀t)e^{−λt}.

Hence T_d ≤ τ as soon as (1+a₀τ)e^{−λτ} ≥ (3/2)M₀/|ω₀(x₀)|; with |ω₀(x₀)| = M₀ and λ = 0,
**T_d ≤ 1/(2a₀)**, i.e. M₀T_d ≤ M₀/(2a₀).
*(Dimensional correction, refuter-confirmed: the attempt wrote "M₀T_d ≤ 1/(2a₀)", which is
dimensionally inconsistent. The corrected form is above.)*

**The Duhamel/Taylor route the brief suggested is UNAVAILABLE, not merely unproved (PROVED).**
Along the trajectory ω″ = (β+νV)ω (the a² cancels), so the Taylor route closes iff
C = sup|β+νV|/a₀² ≤ 2/3, i.e. iff |∂_t²ω| ≤ (1/4)M³log²(R/ρ₀). The datum's measured C ≈ 1.5–2.4.
At the brief's suggested C = 4 the Taylor polynomial never reaches 3/2 at all.

**L1 (PROVED, independently confirmed to 1.6e-11 by direct filament Biot–Savart).**
a(0,0,t) = ∫K ω^θ d³x with K = −(3/8π) r z/ρ⁵, exact for every axisymmetric no-swirl field at
every t; and for the sharp shell ∫₀^π sgn(cos φ)sin²φ cos φ dφ = 2/3 exactly, so
a(0) = (M/2)log(R/ρ₀) **exactly** — which fixes the origin-to-material-shell offset at exactly
the frame's +0.216773 M.

**Sign law, extended (PROVED + refuter extension).** K·ω^θ ≥ 0 pointwise at the origin for all
t; the refuter's attempted refutation of the equatorial case **failed** (credit to the seat):
the sign law also holds at every field point **on** the equatorial plane (128/28000 sources
violate at ~3e-8 = round-off). It fails at 54–65 % of sources with O(1) magnitude at every
**off**-equator point.

### 1.5 Route B — what is REFUTED (MAJOR)

**R-B1. "(H2) is EQUIVALENT to a(1+a₀t)/a₀ ≥ 1" is FALSE — one-way only.** The floor test
measures Theorem D's **conclusion**, not its hypothesis. The seat's own `d4_results.json`
contains the counterexample: tracers N3-Re0=1 at s₀=1.5 and 2.0 have
min(β+νV)/a² = −0.1632 and −0.2366 (H2 violated) while the gate certifies HOLDS at
min F_a = 1.001054 and 1.001062. A 112-run ODE sweep returns HOLDS at β = −3.0a₀² once β was
1.5a₀² early — less than the seat's own measured 1.9–2.4.

**R-B2. FL-043 again.** The F_r half of the gate passes on the sign-reversed control
(min F_r = 1.000017 ≥ 0.99) — it **cannot fire**. The NOTE's headline "worst forward minimum …
1.000000" **is** that non-diagnostic quantity, actual value 0.99999978, rounded up through the
floor it is compared to; the diagnostic half's worst is 1.001054. And F_a ≥ 1 is automatic
whenever a is non-decreasing, true on 25/36 forward tracers — the same 25/36 reported
separately as "β ≥ a²", which by the identity D_ta = −a²+β+νV is the **identical** statement,
not a second corroboration.

**R-B3. "P1 ≥ 0 proved" fails at its premise.** a < 0 on **56 %** of supp(ω) (min a = −0.4504
against max +1.7149 at N=5), at every N = 3…6 and every support threshold 1e-1…1e-6. P1 > 0 is
a measured number whose negative piece (−0.01121, L-independent) is small against a piece
growing like L² — an asymptotic argument, not a sign law.

**R-B4. The "exact split" β(0) = a² + P1 + P2 is the EULER identity.** The NS identity is
β(0) + νV(0) = a² + P1 + P2 + **P3**, P3 = ν∫K r Δ₅η < 0. At Re0 = 100, P3/a² = −0.0010…−0.0016
and accounts for most of the 0.9–1.3 % residual the seat left unattributed. **At the frame's
viscous floor Re0 = 1 — the only regime where c₂ is defined — P3/a² = −0.1011 (N=5) /
−0.1576 (N=4)**, and restoring it drops the split's residual from 27–49 % to 3–6 %.

**R-B5. Therefore (H2) at the origin does NOT reduce to one inequality on one integral.** It
reduces to P2 + P3 ≥ −(a² + P1): three terms, two of unknown sign, one known negative.

**R-B6. "Whoever bounds P2 bounds BFG's step" is an analogy sold as an identification.**
P2 carries a velocity factor, so its effective weight is ρ^{-3}, not ρ^{-4} (the seat's own
§6b says so, contradicting its key); BFG's (|x|+1)^{-4} carries a fixed length scale that
∇K does not; and the implication runs the wrong way.

**R-B7. GAP 2 (origin → material shell) is STRUCTURAL, not O(1/L).** The size transfers exactly
(+0.216773 M against a leading (M/2)L). The **sign structure** does not: it is supported
precisely on the set where ω ≡ 0 — the origin, and the equatorial plane, on which any
continuous field odd in z vanishes identically (measured 5.0e-16 for this datum). At the
seat's own tracer seed φ = 45°, where |ω| = 0.982M, the same strain functional already has
negative part/a = 0.248 (N=4) and 0.183 (N=5).

### 1.6 The strongest surviving statement about (ii)

> **CONDITIONAL THEOREM (survives).** For the mollified δ-tapered plateau family with
> a₀ = κML + c_a M and log Re_E = 2L + c_E, if along the trajectory carrying |ω₀| = M₀ over the
> window a₀t ≤ 1/2 one has (H2) β + νV ≥ 0 and (H3) νW ≥ 0, then
> T(Re) ≤ (c₂ + o(1))/log Re with **c₂ = 2θ/κ**, θ = 1/2, i.e. **c₂ = 1/κ**.
> Under the strictly stronger (H2\*) β + νV ≥ a², θ = log(3/2) and c₂ = 4log(3/2) = 1.6218604.
> Under the full T_λ model (Route A, conditional on GAP T, GAP V, L3v, R),
> c₂ = 8(1−√(2/3)) = 1.4680274.
> None of (H2), (H2\*), GAP T, GAP V is proved.

**CORRECTION TO THE FRAME (established by proof; unchallenged by any refuter).** The frame's
"c₂(first order) = 4 log(3/2) = 1.6219" is **not a bound**. It requires β + νV ≥ a² (strain
non-decreasing), not merely β + νV ≥ 0. Under β + νV ≥ 0 the correct constant is exactly
**2** (this seat's recomputation: 2θ/κ = 2·0.5/0.5 = 2.000000). c₂ = 2 is the **infimum** over
the family as δ,w → 0, where ‖η‖_∞ diverges — not attained: for the datum the 16 runs actually
used (κ = 0.48137) the ceiling is **2.0774**; at δ = 30° it is 2.2067.

**SECOND CORRECTION TO THE FRAME.** The frame's log Re_E = 2L − 0.703 uses the sharp shell's
E = 0.172403978 M²R⁵ ⇒ c_E = −0.7031659. The datum the runs actually used measures
E = 0.1418961 M²R⁵ ⇒ c_E = −0.7810641 (this seat: `synth/s1_results.json`). Mismatch
**0.0779 in log Re_E, in the direction that inflates c₂**.

---

## 2. STATUS OF BFG (ARMA 2019, arXiv:1704.05546v4) THEOREM 10: **UNTOUCHED**

Not "refuted as stated by a proved theorem". Not "contradicted by numerics only". **Untouched.**

- Both proof seats explicitly decline the refutation. `prove-lagrangian`: "CONSEQUENCE FOR (iii):
  NOT ASSERTED … this seat does NOT refute Bradshaw–Farhat–Grujić Thm 10." `prove-duhamel`:
  "NOT ENGAGED … This seat proves no upper bound on T(Re), so it does not refute Thm 10."
  Both refuters confirm the disclaimers are correct.
- The contradiction is **real and correctly located**: Theorem D + (H2) would give
  M₀T_d ≤ M₀/(2a₀) ~ 1/(κL) → 0 against Thm 10's energy-independent M₀T_d ≥ 1/c(3/2) > 0, same
  data class (smooth finite-energy 3D NS), same quantity (sup-norm growth time). It is not
  **triggered**, because (H2) is unproved.
- The located defect in BFG's proof (p.10: ∫|f|/(|x|+1)⁴dx ≤ c‖f‖_BMO without the local average;
  on the plateau family the left side grows like log(R/ρ₀) while the BMO norm stays bounded)
  was found by two campaign referees and is **not re-verified in this round**. A defect located
  in a published proof is not a refutation of the theorem. BFG's own p.8 introduces the section
  as "we present a sketch here"; no published erratum and no citing paper contesting it were
  found by either literature seat.
- **The stake is symmetric and should be recorded as such.** If BFG Thm 8/10 stands: conjecture
  (ii) is **false**, and the corollary of §4 is **subsumed**. If (ii) is proved: BFG Thm 10 is
  refuted as stated, and the corollary becomes half of a matched pair. The two hang together.
  Neither seat can settle it alone.

---

## 3. DEPLETION VERDICT: **KILLED**, with a pre-registered gate and an advance prediction

`GENERATION_2_BRAINSTORM` §B ("the log is a stock, not a flow") is **refuted in its
consequential half and confirmed in its diagnostic half** — and the two halves concern two
different strains, at two different places in the flow.

**Gate (pre-registered in `depletion-numerics/PREREG.md` before any production run, never re-cut):**

| clause | requirement | measured | fires? |
|---|---|---|---|
| S1 | V_{T4} < 1.5 | **2.2518** | FALSE |
| S2 | V_{T32} > 2 | **2.9346** | TRUE |
| S3 | C < 1.5 by t = 1.5 for every N | 0, 0, 0, 0.2061, 0.2610 | TRUE |
| KILL | \|V_{T4}/V_{T32} − 1\| ≤ 0.30 | **0.143** (common set) / **0.233** (own set) | TRUE both ways |

**VERDICT: DEPLETION KILLED.**

**The numbers.** Datum A, Re0 = Mρ₀²/ν = 100, h = ρ₀/8, box 1.5R:

| N | log Re_E | T_{3/2}M₀ | T₂M₀ | T₄M₀ | T₄/T_{3/2} | C(0) |
|---|---|---|---|---|---|---|
| 3 | 7.852 | 0.6667 | 1.1306 | 3.9070 *(excluded)* | 5.860 | 0.835 |
| 4 | 9.239 | 0.4477 | 0.7428 | 2.0197 | 4.511 | 1.112 |
| 5 | 10.625 | 0.3379 | 0.5548 | 1.4322 | 4.238 | 1.389 |
| 6 | 12.011 | 0.2716 | 0.4433 | 1.1043 | 4.065 | 1.666 |
| 7 | 13.398 | 0.2272 | 0.3693 | 0.8969 | 3.948 | 1.944 |

Fitted exponents log(T_kM₀) = α − β log N: β_{T32} = **1.2690** (N = 3…7),
β_{T2} = 1.3185, β_{T4} = **1.4517** (N = 4…7). The quadrupling clock falls with the octave
count **at least as steeply** as the half-doubling clock does — the opposite of §B's prediction.

**The diagnostic half is confirmed, and explains why it doesn't matter.** C(t) = max_x a/‖ω‖_∞
does collapse, and its numerator max_x a does become N-independent
(a_max(1.5) = 0.974/1.037/1.074/1.099 at N = 3/4/5/6, a 6 % spread over a factor 8 in R). But
max_x a sits **on the axis at the origin** at every recorded time and every N, while by T₄ the
vorticity maximum has migrated to s* = 9.4…10.9 deep inside the stack, where the local strain
is **not** depleted: a* = 0.6439/0.9596/1.2293/1.4754 at N = 4/5/6/7, still rising near-linearly
with N. The stack is spent from the inside out and the maximum is handed outward to
successively larger shells, each of which still has the outer octaves above it. **The clock
never reads the axis.**

**Advance prediction (written and logged while N=7 stood at t=0.50, ‖ω‖=2.52).** Two
independent routes gave T₄(N=7) = 0.9009 and 0.8941, mean **0.8975**; depletion required
≥ 1.1043. Measured **0.8969** — prediction 0.07 % high, 19 % below what depletion required.

**Supplementary (not gated).** 1/(T_kM₀) = A_k log Re_E + B_k fits to better than 0.4 % at every
point: T_{3/2} = 1.912/(M(log Re_E − 4.975)), T₂ = 3.043/(M(log Re_E − 5.150)),
T₄ = 6.708/(M(log Re_E − 5.929)). Mean stretching rates
α = 0.21208 / 0.22781 / 0.20667 — the T₄ and T_{3/2} coefficients differ by **2.6 %**.
The sharp seat's invariant reproduces: Q_{T32} = T_{3/2}M log(R/s*) = 0.9710/0.9924/0.9833/
0.9787/0.9760 at N = 3…7 (1.1 % spread against their 0.988 ± 0.018) — but it does **not**
extend: Q_{T4} = 1.066/1.622/1.976/2.206.

**Consequence for the payoff §B hoped for: none.** No numerical support for "BFG Thm 10 false
near growth factor 1, true at factor ≥ 4". The 4M₀ clock is as energy-dependent as the 1.5M₀
clock, same log Re_E slope to 2.6 %.

**Caveats that must travel with the verdict.**
1. **D3 sign control FAILS at the extended horizon — a correction to the sharp seat's C1.** The
   sign-reversed N=5 datum creeps up: 1.005 (t=0.5), 1.120 (2.0), 1.499 (3.85), peak 1.583 at
   t = 4.004; it reaches 1.5M₀ at t = 3.853. The anti-datum grows too at long times. The sharp
   seat's "reversed sign never reaches 1.5M" held only because its horizon was ≤ 0.8. Every T₄
   used sits ≥ 3.5× above this background; N=3's T₄ carries a 53 % sign-blind component and is
   **excluded** on two independent pre-registered grounds (s*/L = 0.715 > 0.5, and the control
   itself at 1.5277).
2. One Reynolds number (Re0 = 100), one datum family. Not tested against Re0 variation; the
   sharp seat's own sweep shows T_{3/2} saturating below Re0 ≈ 16.
3. Tested at k = 2 (quadrupling) only, never at k = 3: T₈ is essentially unmeasurable in this
   family at Re0 = 100 (N = 3/4/5 peak at 4.04/5.68/7.31 and viscosity turns them over; only
   N=6 crossed, box-flagged at s*/L = 0.596).
4. Four points in the T₄ fit; L = 2.08–3.47, far from asymptotic.
5. Controls passing: D1 (bit-exact reproduction of the sharp seat, 0.000e+00 relative at every
   N), D2 (sup|η| overshoot ≤ 2.2e-3), D4 (3.4 % over three grids; 1.4 % between the two
   finest), D5 (0.7 % on a doubled box). DEV-3: D4 at N=6 cancelled for CPU. DEV-4: PREREG's
   "across N = 3…7" was not like-for-like; **both readings reported**, and before N=7 landed
   they disagreed — the seat had committed in writing to the reading that fires.

**Numerics falsify, never prove.** This kills §B's mechanism. It does not establish (ii), and
the campaign's own six viscous runs still do not exclude a log-free clock with a large constant
(saturating fit with floor M T_d ≥ 0.011 fits 4× better).

---

## 4. THE COROLLARY: **PROVED (conditional on theorem (i)); NOVELTY REFUTED**

### 4.1 Statement (survives both refuters intact)

> **COROLLARY.** ν>0, u₀ divergence-free finite-energy in H^m(ℝ³), m ≥ 4 integer, u the maximal
> strong solution in C([0,T_max);H^m_σ(ℝ³)), blowing up at T := T_max < ∞. E(t) = ‖u(t)‖₂²,
> M(t) = ‖ω(t)‖_∞, Re_E(t) = E(t)^{2/5}M(t)^{1/5}/ν, E₀ = E(0), K = E₀^{2/5}/ν, c₁∈(0,1) the
> constant of theorem (i). Then for **every** t ∈ [0,T):
> (a) M(t) ≥ c₁/((T−t)(1+log₊Re_E(t)));
> (b) M(t) ≥ c₁/((T−t)(1+log₊(E₀^{2/5}M(t)^{1/5}/ν)));
> (c) M(t) ≥ c₁/((T−t)(1+log₊(K(T−t)^{-1/5}))) ⇒ M(t) ≥ (5c₁+o(1))/((T−t)log(1/(T−t)));
> and ∫₀ᵗ‖ω(s)‖_∞ ds ≥ 5c₁ loglog(1/(T−t)) − C(E₀,ν,T).
> All three forms are scaling-invariant under u ↦ λu(λx,λ²t) at fixed ν (three exact symbolic zeros).

**Sharpenings established by refuter 2, adopt them:** (b) is **equivalent** to M ≥ M_fix, the
unique root of M = c₁/((T−t)(1+log₊(KM^{1/5}))), since the right side is strictly decreasing in
M (symbolic). The cheap sharper explicit form is
(c′) M(t) ≥ c₁/((T−t)(1 + [log K + (1/5)log(1/(T−t)) + (1/5)log c₁]₊)), checked ≤ M_fix at every
grid point; max M_fix/m = 2.086984 over a 4×4×3 grid.

**The proof is one line, not five.** Theorem (i)'s own §5 concludes that *the actual solution*
reaches t+H(t) ("Thus the actual solution reaches H"), so T−t ≥ H(t) is (a) rearranged. The
attempt's Lines 0–5 are correct but Lines 2–4 unwind a definition. Content is §3.3 (bootstrap)
and §4 (accumulation) only.

**Both refuters tried to break the proof and could not.** The class attack (GIM restart gives a
bounded-velocity mild solution while T_max is defined in C([0,T);H^m_σ), so persistence of
regularity would be an unstated step) is **answered inside theorem (i)'s own §5**. §3.3's
bootstrap and §4's antiderivative were re-derived independently, the latter by numerical
quadrature rather than the symbolic route (V = 1e8: 84.0562 both ways, 9 decimals).

### 4.2 Novelty: **NOT NOVEL AS FRAMED** (refuted twice, independently, at MAJOR)

The attempt's verdict was "NOVEL-CONDITIONAL, closest true analogue is Euler
(Ingimarson–Kukavica 2026, limsup only)". That is false in three ways, and two of the three
refutations were sitting in files the attempt had already downloaded.

- **N1 (refuter 1).** *Cortissoz–Montero–Pinilla, J. Math. Phys. 55, 033101 (2014)*: for 3D
  **Navier–Stokes**, for every t < T, ‖u(t)‖_{Ḣ^{5/2}} ≥ c/((T−t)|log(T−t)|) — row 0's exact rate
  shape, NS not Euler, pointwise not limsup, constant depending only on ‖u₀‖_{L²}. Quoted
  verbatim from McCormick et al.'s extraction. The attempt's own `txt/rss-1503.03063.txt`
  cites it (line 458) and contains the string `(T − t) |log (T − t)|` in a displayed NS lower
  bound at line 44. The sweep grepped that file for `vortic|omega|curl`, got 0, and stopped.
- **N2 (refuter 2, stronger).** *McCormick–Olson–Robinson–Rodrigo–Vidal-López–Zhou,
  arXiv:1503.04323 (SIAM J. Math. Anal. 48 (2016)), **Theorem 4.1, PDF p.8**, verbatim:*
  "Suppose that u is a classical solution of the Navier-Stokes … with maximal existence time T.
  Then ‖u(t)‖_{Ḃ^{5/2}_{2,1}} ≥ c/(T−t)." Ḃ^{5/2}_{2,1} for velocity = Ḃ^{3/2}_{2,1} for
  vorticity and has **exactly the scaling of ‖ω‖_∞**. So a published NS lower bound at rate
  (T−t)^{-1}, **log-free**, for **every** t < T, with an **absolute** constant, in the
  vorticity's own scaling class, has existed since 2015. Three of the four advantages the
  attempt claimed are matched and the rate is beaten by the whole logarithm.
- **N3 (both).** The same paper's introduction, p.2: "In fact this result, and all subsequent
  lower bounds, are a consequence of upper bounds on the local existence time" — the
  corollary's own method, named in the literature as the standard route.
- **N4 (refuter 1).** "First viscous counterpart of row 3" is unsupported: the
  Ingimarson–Kukavica proof of Thm 2.1 is **viscosity-blind** — all four ingredients (Lagrangian
  Grönwall on ‖ω‖_∞, Kozono–Taniuchi, the H³ estimate, and ODE manipulation) survive adding νΔ
  with a favourable sign, so an NS version follows by transcription of a published proof.

**What survives as novel, narrowly.** No published pointwise-in-t lower bound on
‖ω(t)‖_{L^∞} for 3D Navier–Stokes with ν > 0, universal constant, and all data-dependence
inside a scaling-invariant Reynolds number, was located — **with the standing exception of BFG
Thm 8**, which is strictly stronger (log-free, energy-free, ν-uniform) and **implies** the
corollary. None of the located NS bounds transfers to ‖ω‖_∞ in either direction: they are all
L²-based Sobolev/Besov norms not controlled by (E, ‖ω‖_∞) either way (Ḃ^{3/2}_{2,1} embeds in
L^∞, which runs the wrong way for a lower bound).

### 4.3 Two further corrections that must travel

- **The "(5c₁+o(1))" coefficient is unreachable.** Independently recomputed by this seat
  (`synth/s1_results.json`, matching refuter 2 to all printed digits): effective coefficient
  M·s·log(1/s)/c₁ under form (c) is **3.6713 / 4.2339 / 4.6625 / 4.8937** at
  T−t = 1e-6 / 1e-12 / 1e-30 / 1e-100. Convergence is 1/log V. The loglog form's effective
  coefficient first reaches 4.5 at T−t = 10^{−4.24e6}.
- **FL-043 on block (B).** The 4-million-sample "0 violations" gate is a **direction check
  only**. It passes every weakening of the constant (c₁ → c₁/10⁶: 0 violations); an exponent
  error of 5 % in Re_* is caught by one draw in four million; and refuter 2 exhibited a
  **false** mutant, 1/5 → 0.199, that returns 0 violations on **20/20 seeds** — with an explicit
  60-digit witness (s = 1e-4000, K = 1e-796, c₁ = 0.5, M = 0.051981/s·(1+1e-9): (b) true,
  mutant threshold violated). Replace it with the exact decision procedure ((b) ⇔ M ≥ M_fix),
  which is seed-free and sampling-free.
- **A fabricated citation was caught and must be propagated as a warning.** A WebFetch summary
  of arXiv:1503.03063 reported "Theorem 1.1 (Robinson, Sadowski, Silva): ‖ω(t)‖_∞ ≥ C(T−t)^{-1}
  for 3D Navier–Stokes". The paper is Cortissoz–Montero, about Ḣ^s norms of the **velocity** on
  𝕋³ at ν=1; `grep -c -i -E "vortic|omega|curl"` over its full extraction returns **0**. Had it
  been believed, the corollary would have been wrongly declared subsumed. Refuter 2 corrected
  the attempt's over-correction: Robinson–Sadowski–Silva, J. Math. Phys. 53 (2012) 115618 **is**
  a real paper by those authors (velocity Ḣ^s), so "the names are a fabrication" is itself wrong
  — the *attribution to that arXiv number and to vorticity* was the fabrication.

---

## 5. WHAT WOULD CLOSE EACH REMAINING GAP

**GAP T (Route A; the whole theorem).** Two things, neither written anywhere in the campaign:
(a) a proof that the Euler/NS flow map of this datum on [0,τ], τ = c/(ML), equals the
shell-wise coaxial strain T_{λ(σ)} up to sup|Φ(y)−T_λ(y)|/|y| = O(c/L) — noting the shell-wise
family is **not** volume-preserving in ℝ³ and is not the image of a single T_λ, whereas the true
flow is incompressible; and (b) the scale-invariant Lipschitz estimate
|a[η∘Φ^{-1}](0) − a[η∘T_λ^{-1}](0)| ≤ Cμ M L for the Calderón–Zygmund-type functional
a[η](0) = (3/8π²)∫(−z_y)|y|^{-5}η dy₅.
*Cheap next test before attempting either:* rerun the refuter's `gegen5d.py` with q_l built
from a **σ-dependent** λ(σ) = (1−(1−σ)κθ/2)^{-2} and compare against κH(σ,θ). Both the attempt
and both refuters have only ever strained **every shell by the same λ**; the integro-ODE needs
the shell-dependent map. This is the cheapest unrun test in the round.

**GAP V (Route A; viscosity).** An Aronson-type Gaussian upper bound for
∂_tη + u·∇₅η = νΔ₅η whose constants depend on the drift only through ‖∇u‖_∞·τ = O(c), giving
|η(x(t),t) − η_transported| ≤ ‖η₀‖_∞ exp(−d²/(Cντ)) with d ≳ ρ₀δ = √(ν/M). Note the plateau
identity r³Δ₅(1/r) = −1 that makes the bulk penalty ν/r² is exact **only in the bulk**; near
the z=0 jump and the taper corners Δ₅η is not controlled by it, and the runs show the actual
η-loss is dominated by the **mollification layer**, not the plateau (3.0 % at 1.5ρ₀ vs 0.1 % at
3ρ₀ at Re0 = 100). A clean (H3) must be stated at a fixed distance from that layer, at an O(M)
cost in a₀ that the frame does not price.

**L3v uniformity.** sup_{λ∈[1,3/2]} Σ_{l≥3}|H_l(λ)|·‖C_l^{3/2}‖_∞/(l(l+3)−4) < ∞, plus
|H_l| ≲ l^{-3/2} (a jump-discontinuity estimate; easy). And the displayed constant must stop
being truncation-dependent: the attempt's 0.3958 truncates at LMAX = 601; streamed to l = 2001
the partial sums are 0.34990/0.38858/0.39582/0.40015/0.40529/0.40837/0.41046 at
l = 101/401/601/801/1201/1601/2001, terms ≈ 0.8254 l^{-3/2}, **true sum ≈ 0.4295** — 8.5 % above
the displayed value, and exactly the number the gate is blind to.

**(H2) (Route B; the whole theorem).** After the refutation this is **not** one scalar
inequality. What is needed: an a-priori bound on P2 + P3 ≥ −(a²+P1) **at a point that carries
vorticity**. The origin and the equatorial plane are exactly the sets where ω ≡ 0, so the whole
sign apparatus must be rebuilt from scratch at, e.g., L3's stationary angle cos²φ = 1/3. At
every such point the strain functional is already a difference with a 14–25 % negative part and
the kernel has the wrong sign on 54–65 % of source space with O(1) magnitude.

**Bookkeeping (R).** ‖ω₀‖_∞ = M(1+O(ε)), E = C_E M²R⁵, Re_E — admitted unwritten by the seat and
not re-derived by any refuter. Must use the **mollified** datum's C_E, not the sharp shell's
(0.0779 in log Re_E).

**BFG Thm 10.** Only two things settle it: (a) a proof of (ii), which refutes it as stated; or
(b) an independent audit of BFG's §pp.8–11 that either repairs the BMO step or confirms the
defect is load-bearing. Route (b) is cheap relative to (a) and has not been run in this round —
the two campaign referees' location of the defect is the only evidence, and it has never been
re-verified by a third seat.

**Theorem (i).** Still **not re-verified by anyone** — not the attempt, not either refuter, not
this seat. Everything in §4 is conditional on it. Its §5 continuation is the single step the
corollary leans on; §§2–4 (kernel/BMO/entropy estimates) have been read by no one on this side.
This is the largest unaudited dependency in the round and should be the next seat.

**Gate repair (FL-043 ×2, before any of these numbers is quoted onward).** Route A: remove the
`*0+` mask, point that check at the value it actually names, delete or re-source the five
literal-vs-literal checks, source the viscous cross-check from the viscous-numerics data rather
than a hardcoded 0.479, and fix `("ALL %d CHECKS PASS" % 0)`. Route B: drop the F_r half as
non-diagnostic and build a control that fires on F_a under a milder perturbation than sign
reversal. Corollary: replace the sampler with the exact decision procedure.

---

## 6. GRADE CANDIDATE: **FILED**

Named against the estate rubric quoted verbatim from
`Tricritical Exploration/ESTATE_GRADING_2026-08-04.md` §0 (L-53, L-54), with PIN 4 (novelty is
a precondition) and PIN 5 (grade what is closed) governing.

| item | status | grade | governing reason |
|---|---|---|---|
| Conjecture (ii) | SKETCH | — | not closed; no upper bound on T(Re) proved |
| Theorem D (Riccati comparison) | PROVED | **Filed** | correct and checkable, but elementary comparison over the classical identity D_ta = −a²+β+νV; its hypothesis (H2) is verified in no instance |
| Taylor route unavailable (C ≤ 2/3 vs measured C ≈ 2) | PROVED | **Filed** | a real, closed negative that redirects the campaign (L-57) |
| Lemma 1 (a_λ(0) = (M/2)λL) | PROVED, 5 instruments | **Filed** | correct; deflated to an algebraic consequence of "η₀ depends on r alone" |
| L1 kernel collapse + sign law (+ equatorial extension) | PROVED | **Filed** | classical maximum principle repackaged; the extension is new but never load-bearing |
| Material-point strain uniform in λ, L | measured | below Filed | PIN 5: measured-but-unproved |
| Depletion KILLED | pre-registered, controlled, advance-predicted | **Filed** | numerics falsify, never prove; a clean negative with a 0.07 % advance prediction |
| Corollary (a)(b)(c)+loglog | PROVED conditional on (i) | **Filed** | PIN 4: rate shape, pointwise-in-t property, absolute constant and method all already published for NS (McCormick et al. Thm 4.1; Cortissoz–Montero–Pinilla); BFG Thm 8 strictly stronger and unresolved |
| Two FL-043 gate failures + one fabricated citation caught | methodological | **Filed** | ledger content, not mathematics |

**No Solid. No Major.** A Major would require changing what the field can do; the round's one
candidate for that — refuting a published ARMA theorem by proving (ii) — did not happen, and
both proof seats say so in their own words. This is the **seventh consecutive** estate pass
returning zero Majors.

**The conditional path to a higher grade, stated so it can be priced.** If theorem (i) is
independently audited **and** (ii) is proved, the pair (matching upper and lower logarithmic
clocks, refuting BFG Thm 10 as stated) is a genuine Major candidate — it would settle a
published ARMA theorem and pin the log. Neither half is close: (ii) needs GAP T or (H2), and
theorem (i) has been audited by nobody. If instead BFG Thm 10 stands, (ii) is **false** and the
corollary is **subsumed** — in that branch this round's surviving deliverables are the depletion
kill, the Taylor-route-unavailable negative, and the two gate failures.

---

## 7. LEDGER ROWS THIS ROUND EARNS

- **FL-043 (two fresh instances).** A gate that multiplies its stored input by zero and prints a
  hardcoded pass count; a gate whose diagnostic half cannot fire on its own control and whose
  headline quotes the non-diagnostic half rounded up through its floor. Verified here directly.
- **New candidate law — "a gate must fail on the number it certifies."** 77.2 % of one seat's
  stored numbers do not affect its verdict; a second seat's sampler passes a mathematically
  false mutant on 20/20 seeds. Mutation coverage should be *measured and reported*, not asserted.
- **New candidate law — "a literature sweep reads the citations of the papers it downloads."**
  Both novelty refutations were inside files the attempt had already fetched; one was found by
  a grep that returned 0 and stopped.
- **Source-integrity warning (propagate).** A WebFetch summary invented a vorticity theorem and
  attributed it to real authors at a real arXiv number. Every load-bearing citation must be
  re-read from the downloaded PDF, and the correction to the correction must travel too (the
  authors are real; the attribution was not).
- **Frame corrections ×2.** c₂ = 4log(3/2) requires β+νV ≥ a², not β+νV ≥ 0 (under which the
  constant is exactly 2, infimum-not-attained, 2.0774 for the datum actually run); and the
  frame's c_E is the sharp shell's, 0.0779 in log Re_E away from the mollified datum's, in the
  direction that inflates c₂.
- **Correction to `sharp/viscous-numerics` C1.** The sign-reversed control does reach 1.5M₀, at
  t = 3.853; the original claim held only inside a horizon ≤ 0.8.

# refute-kukavica-at-source — adversarial pass over `gaps/kukavica-at-source`

Refuter seat, DTC-2026-09-06.  Task: try to break the `kukavica-at-source` attempt.  Every
quotation below was read by me from **my own** 220-dpi renders of the PDF (`png/r*-*.png`), not
from the target seat's renders.  Every number came out of a script I wrote and ran
(`scripts/rf1_independent.py`, `scripts/rf2_bmo_limits.py`; raw output `scripts/rf*_results.json`).
Scripts belonging to the target seat were copied to `copy/` and re-run there.  Nothing outside this
folder was written.

## Verdict

**NOT REFUTED.**  Every load-bearing claim survives.  Three defects, all MINOR, none touching a
load-bearing claim; one of them is a real hypothesis-tracking omission and is written up in §4.

---

## 1. Integrity and reproduction

* `shasum -a 256 -c SHA256SUMS` in the target folder: **114/114 OK**, exit 0.
* `SHA256SUMS.round1`: **55** entries (not 54, as the attempt says); exactly one FAILED — `./NOTE.md`,
  which the attempt itself discloses (it appended §7).  The other 54 are unchanged.  So the
  substance of the claim holds and the count is off by one.
* Both target scripts re-run from `copy/`: `v1_verify_numbers.py` and `v2_duhamel_endpoint.py`
  produce output **bit-identical** to the stored `v1_results.json` / `v2_results.json`.
  (This is a weak control — same code, same machine — so §3 recomputes everything independently.)
* The PDF: 492 472 bytes, sha256 `dd07f7fe9a797843d144504d621471b302d7961d60bb8fc54d30b64bc2e9bf52`,
  SHA-1 base32 `ZOMHLCABYDYCMM5HK2Q32XSJJNIGJH6B`, 14 pages, producer `pstill 1.11`,
  CreationDate 2003-07-15 20:55:14, `pdftotext` yields **14 bytes** (no text layer).  All as claimed.
  `pdf/kukavica-JDE194-2003-a28-authorcopy.pdf` and `pdf/wayback-a28.pdf` are byte-identical.

## 2. Every quotation re-read at source, from my own renders

I rendered pp. 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 at 220 dpi and read them as images.
**Every quotation in the attempt is verbatim correct.**  Specifically confirmed:

* **p.1** title *On local uniqueness of weak solutions of the Navier-Stokes system with bounded
  initial data*; Igor Kukavica, USC; keywords; MSC 35Q30, 76D05, 35K15; running head.
* **p.3** the FJR sentence, the integral form `u(·,t)=∫₀^t k*(u⊗u)ds+e^{tΔ}u₀`, and
  `|k(x,t)| ≤ C(n)/(|x|+t^{1/2})^{n+1}`.  (Note (1.1) on p.3 is written in **non**-divergence form
  `∂_tu−Δu+u·∇u+∇π=0`; the divergence form is (2.1) on p.4, as the attempt says.)
* **p.4** (2.1) `∂_tu_k − Δu_k + ∂_j(u_ju_k) + ∂_kπ = 0`, `∂_ju_j = 0`.
* **p.5** Theorem 2.2 with `T ≥ 1/C‖u₀‖²_{L^∞(ℝⁿ)}` and the transgalilean characterisation.
* **p.6** Lemma 3.1 with kernel `∇G` / `∂_jG`, conclusion `‖u(·,t)‖_{L^∞} ≤ Ct^{1/2}‖f‖_{L^∞(0,T,BMO)}`,
  `|∇G(x,t)| ≤ C/(|x|+t^{1/2})^{n+1}`, `‖π‖_BMO ≤ C‖f‖_{L^∞(ℝⁿ)}`, `π = K_{ij}f = ∂_i(−Δ)^{-1}∂_jf`.
* **p.7** the weighted-BMO estimate **with the local average subtracted**, `(cf. [S, Lemma 3.9] or
  [St2, p. 141])`; (3.4) with `∂_jG` on `u_ju_k`, `∂_kG` on `π`, plain `G` only on `u_{0k}`; (3.5);
  Proposition 3.2 `T ≥ 1/(CM₀²)`; the "completness sake" typo.
* **p.8** the iteration and (3.6).
* **pp.9, 10** the de Rham / Lemma 3.3–3.4 arguments — **no heat kernel appears at all** on these two
  pages, which the attempt did not say but which strengthens its exhaustiveness claim.
* **p.11** (3.10) with `ũ_kũ_j ∂_jG` and `π̃₀ ∂_kG`; undifferentiated `G` only in the two
  `(1/ε)ψ′` time-mollifier terms (linear in `ũ`) and against `u_{0k}`; the "heat instead of the
  Oseen kernel" sentence.
* **p.12** the displayed limit, exactly as the attempt transcribes it.
* **pp.13–14** `[FJR]` = Fabes–Jones–Rivière, Arch. Rational Mech. Anal. **45** (1972) 222–240;
  `[S]` = C. Sadosky, *Interpolation of Operators and Singular Integrals*, Marcel Dekker 1979;
  `[St2]` = E.M. Stein, *Harmonic Analysis: Real-Variable Methods, Orthogonality, and Oscillatory
  Integrals*, Princeton 1993.

So the headline is confirmed at source: **every nonlinear kernel in Kukavica is `∂G`**, and the
weighted-BMO inequality he states is **the subtracted one, and it is the only form in the paper**.

**The BFG side, also at source** (`pdf/bfg-1704.05546v4.pdf`, 24 pp., sha256
`c263fc6e…f693c`), all confirmed verbatim: the `[GIM, Ku2]` / "modification of the exposition given
in [Ku2]" sentences (p.8); (5) non-divergence `∂_tω_j − 4ω_j + u_i∂_iω_j = ω_i∂_iu_j` (p.9); (8) and
(9) carrying **undifferentiated `G`** against every right-hand term (pp.9–10); the main-line
`∫|f(x)|/(|x|+1)^4 dx ≤ c‖f‖_BMO` with the subtracted form only in a parenthesis carrying `[St]`,
the main line carrying **no citation** (p.10); Theorem 8 `T ≥ (1/c)‖ω₀‖_∞^{-1}` (p.8);
`G(x,t) ≤ c√t/(|x|+√t)^4` cited `[LR]` (p.10); bibliography `[Ku2]` with the local-uniqueness title
and `[St]` = Stein, *Harmonic Analysis* (p.24).  The typographic slip is real: (8) and (9) print
`ω_i ∂_j u_j`, the estimate two lines later uses `ω_i^{(n)} ∂_i u_j^{(n)}`.

**Provenance of the title correction, at source.**  The archived index
(`html/wayback-kukavica-pdf-index.html`) reads: item **29** "On local uniqueness of solutions of the
Navier-Stokes equations with bounded initial data, *J. Diff. Eq* **194** (2003), 39-50" →
`<a HREF="a28.pdf">`; item **28** "Length of vorticity nodal sets …, Comm. PDE 28 (2003), 771-793" →
`a29.pdf`; item **20** "On the dissipative scale for the Navier-Stokes equation, *Indiana Univ.
Math. J.* **48** (1999), 1057-1081" → `a20.pdf`; item **18** "Self-similar variables and the complex
Ginzburg-Landau equation" → `a19.pdf`.  The brief's proposed title is a different paper, and the
file numbering does not track the list.  Crossref, Unpaywall and OpenAlex all return the
local-uniqueness title, vol. 194, pp. 39-50, published-print 2003-10.

**The version-of-record failure log is honest.**  `logs/probe_vor.tsv` has 14 rows with the HTTP
codes the attempt reports (ScienceDirect 403 ×3 with block code `CPE9` present in the saved body,
Elsevier API 406, CORE 429, fatcat 000, ResearchGate 403, the rest 200 with no full text).  Crossref
does carry an Elsevier open-access licence for this DOI starting **2013-07-17** with
`content-version: vor`; Unpaywall `oa_status: bronze`, `has_repository_copy: false`.  All as stated.

## 3. Every number recomputed independently

`scripts/rf1_independent.py` (mpmath, 30 dps, adaptive quadrature — not numpy Gauss-Legendre or
Gauss-Hermite) plus hand-derived closed forms.

* **`∫G dz = 1` exactly; `∫∂₁G dz = 0` exactly** (odd integrand); **`∫|∇G|dz = 2/√(πτ)`** —
  derived by hand (`∫₀^∞ r³e^{-r²/4τ}dr = 8τ²` ⟹ `16πτ/(4πτ)^{3/2} = 2/√(πτ)`) and matched to
  30 dps at τ = 1, 10⁻², 10⁻⁴: `1.12837916709551257389615890312`, `11.28379167…`, `112.83791671…`.
  **Confirms [V1] / [K1].  Mean one against mean zero.**
* **`∫_{ℝ³}(|z|+a)^{-4}dz = 4π/(3a)`** derived by hand (`∫_a^∞(u−a)²u^{-4}du = 1/(3a)`) and matched
  to 30 dps at a = 0.1, 1, 3.  Hence exactly
  `∫₀^t∫(|z|+√(t−s))^{-4}dz ds = (8π/3)√t` and `∫₀^t∫√(t−s)(|z|+√(t−s))^{-4}dz ds = (4π/3)t`;
  my quadrature returns `8.3775804095727819101` and `4.18879020478639098462` at t = 0.25, 1, 4.
  **Confirms [V3].  The exponent really is the kernel's fingerprint.**
* **Weight equivalence.** `log f` stationary ⟺ `4/(r+1) = 4r³/(r⁴+1)` ⟺ `r³ = 1` ⟺ `r = 1`, so
  `r = 1` is the *unique* positive critical point and `sup (r+1)⁴/(r⁴+1) = f(1) = 8`, `f(0) = 1`.
  **Confirms [V4] and supplies the proof the attempt only asserted.**
* **`[K2]`** (in the target's §4, not in its constants block): grid search gives
  `sup|∇G|(|x|+1)⁴ = 0.9231098781` at r = 2.81082 and `sup G(|x|+1)⁴ = 0.7109819927` at r = 2.37228,
  against the target's `0.9231098779` / `0.7109819926`.  Confirmed.

### 3a. The weighted-BMO table — one correction to the target's precision, and two new exact numbers

`scripts/rf2_bmo_limits.py`.  My **first** attempt (rf1, breakpoints `[0,1,10,R−1,∞]`) *failed*: at
R = 10⁶⁴ it returned `J₀ = 12.6098`, which breaks the monotonicity of `J₀` in `R` — so the failure
was mine, not the target's.  Re-done on log-spaced panels, and cross-checked against the exact
`R → ∞` limits, which depend on no panel choice:

* `avg_{B₁}h_R = 1 − A/log R` with **`A = 3∫₀¹log(1+r)r²dr = 2log2 − 5/6 = 0.552961027786557285501…`**
  (closed form derived by hand, matched to 30 dps).  Reproduces every `avg_B1` entry in the target's
  table exactly.
* **`J₀(R) ↑ 4π·π/(2√2) = √2·π² = 13.9577283992777590683151218723…`**, monotonically.  The target's
  column 10.1492 → 12.0399 → 12.9987 → 13.4782 → 13.7180 → 13.8379 is monotone and converging to it.
  This is a control that could have failed and did not.
* **`J₁(R)·log R → 4π∫₀^∞|A − log(1+r)|r²/(r⁴+1)dr = 10.3196706278586567517…`**  The target's flat
  value `10.31967` is this limit.

My log-panel table: `J₀` = 10.1491837463, 12.0398800767, 12.9987360226, 13.4782322075,
**13.7179803034**, **13.8378543513**; `J₁·log R` = 10.1933806, 10.3184201, 10.3196767, 10.3196768,
10.3196768, 10.3196768.  Agreement with the target to **9–10 significant figures** at every R.

> **Defect (MINOR, precision).** The target prints `J₁·log R = 10.3196681` and
> `J₀ = 13.8378543607` — 9–11 digits — where only ~7 are converged (mine: 10.3196768,
> 13.8378543513; limit: 10.3196706).  Its own prose says "6–7 significant figures", so the prose is
> honest and only the constants block over-prints.  Nothing downstream depends on those digits.

### 3b. The counterexample's BMO norm — asserted by the target, closed here

The claim "the un-subtracted inequality is **false** on this family" hinges on `‖h_R‖_BMO ≤ C/log R`
with `C` independent of `R`.  **Neither the target's round 1 nor its round 2 computes anything for
this**; both assert it.  Without it, `J₀·log R` growing proves nothing.  Closed two ways:

* *Analytically.* `h_R = (1/log R)·(log R − min(log(1+|x|), log R))`.  Truncation against a constant
  is bounded on BMO, so `‖h_R‖_BMO ≤ 2‖log(1+|x|)‖_BMO / log R` — `C` manifestly `R`-free.
* *Numerically, as a control that can fail.*  A lower bound from the mean oscillation over concentric
  balls gives **`‖h_R‖_BMO · log R ≥ 0.245254288168531…`**, flat across R = 10², 10⁴, 10⁸, 10¹⁶,
  10³², 10⁶⁴ (`R5` in `rf1_results.json`).  Had this drifted, the family would have died.

Hence `J₀ / ‖h_R‖_BMO ≳ (13.96 / 2‖log(1+|x|)‖_BMO)·log R → ∞`: the un-subtracted inequality fails
by a full factor `log R`, and the subtracted one survives at `10.31967`.  **The target's conclusion
is right and is now actually established.**

## 4. The one real gap: a hypothesis of BFG's that the comparison never states

BFG p.10 reads, verbatim as printed:

> "a property of a scalar-valued `BMO` function `f` featuring a suitable decay at infinity (e.g.,
> being in the closure of the test functions in the uniformly-local `L^p` for some `p`, `1 ≤ p < ∞`),
> `∫ |f(x)|/(|x|+1)^4 dx ≤ c‖f‖_BMO`"

The target's comparison table (§4) and its returned statement record BFG's inequality **stripped of
that decay hypothesis**.  It is quoted nowhere in the target's NOTE.  This matters twice:

1. It is the hypothesis that excludes the one-line counterexample (`f ≡ const`: `‖f‖_BMO = 0`, left
   side positive).  A reader of the target's NOTE alone would think BFG had made a beginner's error.
   They did not; their inequality is false for a subtler reason.
2. It is what makes `h_R` the *right* counterexample: `h_R` is bounded and compactly supported,
   hence lies in the closure of the test functions in uniformly-local `L^p` for every `p`, so it
   **satisfies BFG's stated hypothesis** and still kills the inequality by a factor `log R`.

So the refutation survives — it is strengthened — but a refuter is required to check that every
hypothesis is stated, and this one was not.  The corrected statement in §5 restores it.

Two lesser slips: `SHA256SUMS.round1` has **55** entries, not 54; and `png/verify/` holds all
**14** renders (v1…v14) plus `a29-01.png`, not the ten the target's file list names.  Both are
understatements or miscounts, neither affects a claim.

## 5. Corrected statement

Unchanged from the target except: (a) BFG's inequality must be quoted **with** its decay hypothesis
("a scalar-valued BMO function f featuring a suitable decay at infinity, e.g. being in the closure
of the test functions in the uniformly-local L^p for some p, 1 ≤ p < ∞"), and the counterexample
noted to satisfy it; (b) the weighted-BMO constants are good to ~7 significant figures, with exact
limits `J₀ ↑ √2π² = 13.95772839927775906832` and `J₁·log R → 10.31967062785865675166`,
`avg_{B₁}h_R = 1 − (2log2 − 5/6)/log R`; (c) `SHA256SUMS.round1` has 55 entries.

## 6. Remaining gap

The Elsevier version of record is still not obtained, and both passes say so.  What we hold is the
author's accepted manuscript (created 2003-07-15, ~7 weeks before the 2003-09-03 online date); page
numbers are manuscript pages, not JDE 39–50.  I made no further acquisition attempt; the target's
14-route log is complete and its reading — a bot block, not a paywall — is supported by Crossref's
2013-07-17 `vor` open-access licence.  A human with a browser remains the cheapest route.

Separately, the target's own §7.5 open item stands exactly as it flags it: both of BFG's nonlinear
terms are divergence-form (`u_i∂_iω_j = ∂_i(u_iω_j)`, `ω_i∂_iu_j = ∂_i(ω_iu_j)`), so (8) *could* be
rewritten with a mean-zero kernel, at the price of needing `u_iω_j ∈ L^∞(0,T,BMO)`; BMO is not an
algebra and `u` is not controlled in `L^∞` by `ω₀ ∈ L²∩L^∞` without the logarithm at issue.  I did
not attempt to close it either.

## 7. Files

`NOTE.md` · `scripts/rf1_independent.py`, `scripts/rf1_results.json` · `scripts/rf2_bmo_limits.py`,
`scripts/rf2_results.json` · `png/r{1,3,4,5,6,7,8,9,10,11,12,13,14}-*.png` (my own 220-dpi renders,
the read-by-eye evidence) · `copy/` (the target's two scripts, re-run there, with their output) ·
`SHA256SUMS`.

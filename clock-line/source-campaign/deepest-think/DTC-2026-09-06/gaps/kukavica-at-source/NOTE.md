# kukavica-at-source — Kukavica, J. Differential Equations 194 (2003) 39–50, obtained and read

Sub-seat `gaps/kukavica-at-source`, DTC-2026-09-06.  Task: get the Kukavica paper **at source** —
the paper BFG say their §2 is "a modification of the exposition given in [Ku2]" — and settle
whether Kukavica's original carries the kernel `∇G` (integral zero, so the local average may be
subtracted for free) or the undifferentiated `G`, and with or without a subtracted average in the
weighted-BMO inequality.

**Sci-Hub was not used.**  Every route tried is logged in `log/round1.tsv`, `log/round2.tsv`
(URL, HTTP code, bytes, content-type, saved file).  Every number tagged `[K#]` came out of
`scripts/s1_kernel_numbers.py`, which I wrote and ran; raw output `scripts/s1_results.json`.
`SHA256SUMS` is computed, never typed.  Nothing outside this folder was written.

---

## 0. Headline

**OBTAINED.**  The author's own PDF of the paper, from his (retired) USC homepage
`www-bcf.usc.edu/~kukavica/pdf/a28.pdf`, recovered through the Internet Archive.
14 pages, produced 2003-07-15, cover sheet + abstract + 12 pages of body (matching the journal's
39–50).  It is the **author's accepted manuscript, not the Elsevier version of record** — see §2
for why I am confident it is the same paper, and §5 for the caveat this places on page numbers.

Three answers, all verbatim at source:

1. **The kernel in Kukavica's original is `∇G` — every single nonlinear occurrence.**
   Lemma 3.1, the mild equation (3.4) and the iteration all carry `∂_jG(x−y,t−s)` / `∂_kG(x−y,t−s)`.
   The undifferentiated `G` appears only against the *initial datum* `u_{0k}` (the heat semigroup)
   and in the §3 limiting arguments.  There is no undifferentiated `G` against a nonlinear term
   anywhere in the paper.
2. **The weighted-BMO inequality appears in Kukavica exactly once, and it is the SUBTRACTED form:**
   `∫ |f(x) − (1/|B₁|)∫_{B₁} f| / (|x|^{n+1}+1) dx ≤ C‖f‖_BMO` (manuscript p.7), cited to Sadosky
   and Stein.  The un-subtracted form BFG display in their main line **does not occur in Kukavica.**
   What BFG relegate to a parenthesis *is* Kukavica's actual inequality.
3. **Existence input and time:** the local existence for `L^∞` **velocity** data is
   **Fabes–Jones–Rivière**, ARMA 45 (1972) 222–240, taken as known and re-sketched; the time is
   `T ≥ 1/(C‖u₀‖²_{L^∞})` (Theorem 2.2 and Proposition 3.2).  Quadratic in the datum, log-free,
   energy-free — i.e. the GIM/velocity clock, not a vorticity clock.

And one new structural fact, verified numerically here:

> **The exponent in the existence time is the kernel's fingerprint.**  Kukavica's estimate carries
> `t^{1/2}` and yields `T ≳ ‖u₀‖_∞^{-2}`; BFG's carries `t` and yields `T ≳ ‖ω₀‖_∞^{-1}`.  The two
> powers are exactly the two kernels' own Duhamel scalings — `∫₀^t∫(|z|+√(t−s))^{-4}dz ds = (8π/3)√t`
> against `∫₀^t∫√(t−s)(|z|+√(t−s))^{-4}dz ds = (4π/3)t` **[K3]**.  Moving the derivative off the
> kernel is not a presentational change; it is what converts an inverse-square clock into an
> inverse-linear one.

---

## 1. Title — the brief's proposed title is wrong (settled at source)

The brief asks to verify "On the dissipative scale for the Navier–Stokes equation?".  It is **not**
that paper.  Five independent sources agree:

| source | says |
|---|---|
| the obtained PDF's own title page (`png/ku/p-01.png`) | *On local uniqueness of weak solutions of the Navier-Stokes system with bounded initial data* |
| Crossref for DOI 10.1016/S0022-0396(03)00153-0 (`raw/probe_crossref.json`) | same title; JDE **194**, pages **39-50**, published-print 2003-10 |
| Unpaywall / Semantic Scholar (`raw/s2.json`, prior seat's `unpaywall-kukavica.json`) | same title |
| zbMATH Zbl 1050.35070 (prior seat) | same title |
| **BFG's own bibliography**, arXiv:1704.05546v4 p.24 | "[Ku2] I. Kukavica, *On local uniqueness of weak solutions of the Navier-Stokes system with bounded initial data*, J. Differential Equations **194**, 39 (2003)." |

*On the dissipative scale for the Navier-Stokes equation* is a **different** Kukavica paper —
**Indiana Univ. Math. J. 48 (1999), 1057–1081**, item 20 of his own publication list
(`pdf/kukavica-cv.pdf` line 407, read by the prior seat; the same list, as archived at
`www-bcf.usc.edu/~kukavica/pdf/index.html`, is saved here as `raw/kukavica_pdf_index.html`).
BFG do **not** cite it.  So the paper BFG modify is the local-uniqueness paper, and that is the one
obtained below.

---

## 2. Acquisition — what was tried, what worked, and why the copy is the right one

**The route that worked:** author's homepage → Internet Archive.

* The live homepage (`https://dornsife.usc.edu/igor-kukavica/publications/`, HTTP 200, 326 700 B,
  `raw/dornsife_pubs.html`) carries **no PDF links at all** — confirming the prior seat.
* His old page `www-bcf.usc.edu/~kukavica/` is dead (connection error).  The Wayback CDX index for
  `www-bcf.usc.edu/~kukavica*` (`raw/wb_cdx_bcf.txt`) lists a directory of his own paper PDFs,
  `pdf/a01.pdf … pdf/a82.pdf`.
* `pdf/a28.pdf`, snapshot 2017-06-28, **is** the paper: its title page, keywords
  ("Navier-Stokes equations, weak solutions, mild solutions, nonuniqueness"), 2000 MSC
  "35Q30, 76D05, 35K15", and Theorem 2.2 all match, respectively, the Crossref title, and the
  zbMATH record's MSC list (35Q30, 76D05, 35K15) and its statement of the main result
  ("`ũ(x,t) = u(x − Φ(t),t) + φ(t)`").  The file's own PDF creation date is **2003-07-15**;
  the article was published online 2003-09-03 (Unpaywall).
* **The file numbering is idiosyncratic and was NOT trusted.**  `a29.pdf` is *Length of vorticity
  nodal sets…* (Comm. PDE 28) and `a19.pdf` is *Self-similar variables and the complex
  Ginzburg-Landau equation* — neither matches its list position.  Identification was made from
  the title page, not from the file name.
* **Integrity (`scripts/verify_provenance.sh`, re-runnable).**  Our bytes reproduce the Wayback
  CDX SHA-1 digest of the archived object exactly — `ZOMHLCABYDYCMM5HK2Q32XSJJNIGJH6B`,
  492 472 bytes — and a second, independently addressed snapshot (timestamp `20170628022520`,
  percent-encoded URL) is byte-identical: sha256 `dd07f7fe…c2e9bf52` both times.

**Routes that failed** (all logged with codes):

| route | result |
|---|---|
| ScienceDirect landing page and `/pdfft` PDF endpoint | **403** both (Cloudflare block page, 1.2 MB of HTML) — despite Unpaywall/S2 calling it `BRONZE` open access |
| Elsevier text-mining API (`api.elsevier.com`, the link Crossref advertises) | **400** (needs a key) |
| CORE search | **403** (needs a key) |
| fatcat / `api.fatcat.wiki` release lookup | connection error (unreachable from here — same as the prior seat) |
| scholar.archive.org search | 200 but a JS shell only, 4 428 B, no results in the HTML |
| ResearchGate publication search | **403** |
| USC digital library search | **404** |
| BASE search | 200, no full text |
| arXiv author query (`au:"Kukavica_I"`, 100 results; and `au:Kukavica AND abs:"bounded initial data"`) | **no preprint of this paper exists on arXiv** (0 entries for the targeted query) |
| author's `dreamhosters.com` site | empty root (0 B), `/publications.html` 404 |
| live `www-bcf.usc.edu` | dead host |
| Sci-Hub | **not attempted — disallowed** |

**Not obtained: the Elsevier version of record.**  See §5.

---

## 3. What the paper says — verbatim, with manuscript page numbers

Page numbers below are the **manuscript's own printed page numbers** = PDF page numbers
(p.1 title/cover, p.2 abstract, p.3 §1 Introduction … p.13–14 references).
Page images: `png/ku/p-01.png` … `png/ku/p-14.png`.  A machine OCR of all 14 pages is in
`txt/kukavica-JDE194-2003-OCR.txt` (Apple Vision, `scripts/ocr_vision.swift`); the PDF has **no
embedded text layer** (`pdffonts` returns nothing — the pages are rasterised), so the OCR is a
search index only and **every quotation below was read from the rendered page image by eye.**

### 3.1 The local-existence input for `L^∞` data (p.3)

> "Fabes, Jones, and Riviere proved in [FJR] existence and uniqueness of local mild solutions of the
> Navier-Stokes system if `n < p ≤ ∞` (`p = ∞` included).  They showed that the integral form of
> (1.1)-(1.2) is
>
>   `u(·,t) = ∫₀^t k * (u ⊗ u) ds + e^{tΔ}u₀`
>
> (where `e^{tΔ}u₀` denotes the solution of the heat equation with the initial datum `u₀`) with a
> kernel `k` which satisfies
>
>   `|k(x,t)| ≤ C(n)/(|x| + t^{1/2})^{n+1} ,   (x,t) ∈ ℝⁿ × (0,∞).`"

Note the kernel is a **derivative** kernel (Oseen), of homogeneity `−(n+1)`, hit against the
**quadratic** `u⊗u`, and that the reference is FJR = *The initial value problem for the
Navier-Stokes equations with data in `L^p`*, Arch. Rational Mech. Anal. **45** (1972), 222–240
(p.13).

Later (p.11), on his own construction: "This part follows [FJR] rather closely except that we work
with the heat instead of the Oseen kernel."

### 3.2 The existence time (pp.5, 7, 8)

**Theorem 2.2** (p.5), verbatim:

> "*There exists at least one solution `u` of the initial value problem on `[0,T)` where
> `T ≥ 1/C‖u₀‖²_{L^∞(ℝⁿ)}`.  Also, `ũ` is a solution of the initial value problem on `[0,T)` if and
> only if there exists a function `φ ∈ L^∞_loc([0,T),ℝⁿ)` with `lim_{t→0+}φ(t) = φ(0) = 0` such that
> `ũ(x,t) = u(x − Φ(t),t) + φ(t)`, a.e. `(x,t) ∈ S_T` where `Φ(t) = ∫₀^t φ(s)ds`.*"

**Proposition 3.2** (p.7), verbatim:

> "*There exists a mild solution `u` of the Navier-Stokes system on `(0,T)` where*
> `T ≥ 1/(C M₀²)`.
> *This solution is unique within the class of mild solutions and continuity with respect to initial
> data holds.*"
>
> with `M₀ = ‖u₀‖_{L^∞}` (p.7), and: "This proposition is proven in [FJR]; however, for completness
> sake, we include the proof of existence."

So the existence time is **quadratic** in the datum, `T ≳ ‖u₀‖_∞^{-2}`, log-free and energy-free,
for the **velocity**.  This is exactly the Giga–Inui–Matsui clock the prior seat established at
source, now confirmed in the second of BFG's two cited precedents.

### 3.3 The nonlinear estimate — the kernel is `∇G` and the average IS subtracted (pp.6–8)

**Lemma 3.1** (p.6), verbatim:

> "*Let `T > 0` and `f_j ∈ L^∞(0,T,BMO)` for `j = 1,…,n`.  Then the integral*
> `∬_{S_T} |∇G(x−y,t−s)| |f(y,s)| dy ds`
> *is finite for every `(x,t) ∈ S_T = ℝⁿ × (0,T)`.  The function*
> `u(x,t) = ∫₀^t ∫ ∂_jG(x−y,t−s) f_j(y,s) dy ds`
> *satisfies `u ∈ C([0,T),L^∞(ℝⁿ))` with*
> `‖u(·,t)‖_{L^∞(ℝⁿ)} ≤ C t^{1/2} ‖f‖_{L^∞(0,T,BMO)} .`
> *The function `u` solves (3.3) in `D′(S_T)`.*"

and immediately after:

> "This lemma is a simple consequence of the inequality
> `|∇G(x,t)| ≤ C/(|x| + t^{1/2})^{n+1} ,  (x,t) ∈ ℝⁿ × (0,∞)`"

continuing on p.7:

> "and the estimate
> `∫_{ℝⁿ} |f(x) − (1/|B₁|)∫_{B₁} f| / (|x|^{n+1} + 1) dx ≤ C‖f‖_{BMO}`
> where `B₁` is the ball centered at 0 with radius 1 (cf. [S, Lemma 3.9] or [St2, p. 141])."

`[S]` = C. Sadosky, *Interpolation of Operators and Singular Integrals*, Marcel Dekker 1979;
`[St2]` = E.M. Stein, *Harmonic Analysis*, Princeton 1993 (p.14).

**The mild equation (3.4)** (p.7), verbatim in structure:

```
u_k(x,t) = − ∫₀^t ∫ ∂_j G(x−y,t−s) u_j(y,s) u_k(y,s) dy ds
           − ∫₀^t ∫ ∂_k G(x−y,t−s) π(y,s)           dy ds
           +      ∫ G(x−y,t)      u_{0k}(y)          dy                     (3.4)
```
with `π = K_{ij}(u_iu_j)` (3.5), `‖π‖_BMO ≤ C‖f‖_{L^∞}` (p.6).

**The iteration** (p.8), verbatim:

> `‖u^{(m+1)}‖_{L^∞(S_t)} ≤ C t^{1/2} ‖u^{(m)}‖²_{L^∞(S_t)} + C₀M₀ ,   t > 0`
>
> … "Using induction, we get `‖u^{(m)}‖_{L^∞(S_T)} ≤ 2C₀M₀`, `m ∈ ℕ` provided
> `0 < T ≤ 1/(C M₀²)`   (3.6)  with `C` sufficiently large."

**Every nonlinear kernel is `∂G`.**  The undifferentiated `G` occurs in (3.4) only against `u_{0k}`,
and in the §3 technical limits (pp.11–13).  I checked all 14 pages: `BMO` occurs on pp.3, 4, 6, 7, 8
only (`grep` over the OCR), and the weighted inequality occurs **once**, on p.7, in the subtracted
form quoted above.

---

## 4. The comparison with BFG §2 (arXiv:1704.05546v4, pp.8–11)

| | **Kukavica 2003 (the original)** | **BFG §2 (the "modification")** |
|---|---|---|
| formulation | velocity–pressure, `u₀ ∈ L^∞` | vorticity–velocity, `ω₀ ∈ L²∩L^∞` |
| kernel on the nonlinear term | `∂_jG(x−y,t−s)` — **mean zero** | `G(x−y,t−s)` — **mean one** |
| pointwise kernel bound | `\|∇G(x,t)\| ≤ C/(\|x\|+t^{1/2})^{n+1}` | `G(x,t) ≤ c√t/(\|x\|+√t)^4` (cited to [LR]) |
| the BMO integrand | `f_j ∈ L^∞(0,T,BMO)` (the pressure / the quadratic), a **signed** BMO function | `\|∂_iu_j\|` — the **absolute value**, non-negative |
| weighted-BMO inequality used | `∫\|f − (1/\|B₁\|)∫_{B₁}f\|/(\|x\|^{n+1}+1) ≤ C‖f‖_BMO` — **average subtracted**, cited [S],[St2] | `∫\|f\|/(\|x\|+1)^4 ≤ c‖f‖_BMO` — **average not subtracted**; the subtracted form appears only in a parenthesis, "(in general, one has to subtract a local average …)" |
| power of `t` produced | `t^{1/2}` | `t` |
| resulting time | `T ≳ ‖u₀‖_∞^{-2}` (Thm 2.2, Prop 3.2) | `T ≥ (1/c)‖ω₀‖_∞^{-1}` (Thm 8) |

**Answer to the question the brief poses:** Kukavica's original has **`∇G`**, and the weighted-BMO
inequality is used there **with the local average subtracted**.  So the prior seat's structural
diagnosis (`bfg/source-adjudication` §1c, `[N2]`) is confirmed verbatim at source, and can be
sharpened in two ways:

1. **The subtraction is not merely "available for free" in Kukavica — it is the only form he
   states.**  BFG's parenthetical is not a caveat about an unusual case; it is a description of the
   inequality the source they are modifying actually uses.  Their displayed main-line inequality
   occurs nowhere in Kukavica.
2. **Two independent things block restoring it in BFG's chain**, not one.  (i) The kernel: `∂_jG`
   annihilates constants, so `∫∂_jG·f = ∫∂_jG·(f−c)` for any `c` and the subtraction is exact —
   `∫_{ℝ³}∂_1G dz = 0` (exactly, by oddness; `0.0` in 120-node Gauss–Hermite) while
   `∫_{ℝ³}G dz = 1` to 4e-16 **[K1]**.  BFG's `G` does not annihilate constants.  (ii) The
   integrand: BFG take absolute values *before* applying the inequality, so the object they feed it
   is `|∂_iu_j| ≥ 0`, whose local average is not controlled by `‖∂u‖_BMO` at all.  Even with a
   mean-zero kernel, `|f|` is the wrong argument.

**Numerics that pin the two claims [K5], [K3]** (`scripts/s1_kernel_numbers.py`):

* On the family `h_R(x) = (1 − log(1+|x|)/log R)₊` (`‖h_R‖_BMO ≤ C/log R`), with the weight
  `1/(|x|^{n+1}+1)` Kukavica actually writes:

  | `R` | `J₀ = ∫\|h_R\|·w` | `J₀·log R` | `J₁ = ∫\|h_R − avg_{B₁}h_R\|·w` | `J₁·log R` |
  |---|---|---|---|---|
  | 1e2 | 10.1492 | 46.739 | 2.21346 | 10.1934 |
  | 1e4 | 12.0399 | 110.891 | 1.12031 | 10.3184 |
  | 1e8 | 12.9987 | 239.446 | 0.560222 | 10.3197 |
  | 1e16 | 13.4782 | 496.556 | 0.280111 | 10.31967 |
  | 1e32 | 13.7180 | 1010.778 | 0.140056 | 10.31967 |
  | 1e64 | 13.8379 | 2039.222 | 0.0700278 | 10.31967 |

  The **un**subtracted integral times `log R` grows linearly (46.7 → 2039.2); the **subtracted** one
  is flat at `10.3197`.  The subtraction is not a formality — on this family it is the whole
  inequality.  (Third confirmation of the prior seat's `[N4]`, and the first to show the *repaired*
  form surviving on the same family.)
* Duhamel scalings, `n = 3`: `∫₀^t∫(|z|+√(t−s))^{-4}dz ds` divided by `√t` equals
  **8.377580409573** at `t = 0.25, 1, 4` (`= 8π/3 = 8.377580409573`), while
  `∫₀^t∫√(t−s)(|z|+√(t−s))^{-4}dz ds` divided by `t` equals **4.188790204786** at the same three
  times (`= 4π/3`).  Kukavica's `Ct^{1/2}` and BFG's `ct` are precisely these two.
* The two weights are equivalent, so nothing turns on the difference: `1 ≤ (|x|+1)^4/(|x|^4+1) ≤ 8`,
  the sup attained at `|x| = 1` **[K4]**.
* Sharp constants (parabolic scaling, evaluated at `t = 1`):
  `sup |∇G|·(|x|+√t)^4 = 0.9231098779` and `sup G·(|x|+√t)^4/√t = 0.7109819926` **[K2]**; and
  `∫|∇G(·,τ)|dz = 2/√(πτ)` to 2e-12 at `τ = 1, 1e-2, 1e-4` **[K1]** (independent reproduction of the
  prior seat's `[N2c]`).

**What this does not do.**  It does not decide BFG Theorem 8/10.  It removes the last "NOT OBTAINED"
from the prior seat's §3 chain and converts its (1c) from a structural inference into a
source-verified statement.  The adjudication stands where that seat left it: **(b)** — unsupported
by the displayed argument, uncontradicted by any published NS result.  The size of the repair
(one factor `Λ`, `[N5]`) is that seat's number, not re-derived here.

---

## 5. Caveats, stated plainly

1. **This is the author's accepted manuscript, not the Elsevier version of record.**  I could not
   get past the ScienceDirect Cloudflare block (403 on both the landing page and the PDF endpoint),
   and the Elsevier TDM API needs a key.  The manuscript is dated 2003-07-15, ~7 weeks before the
   article went online.  I have no way to exclude copy-editing or referee-driven changes between
   this file and pp.39–50 as printed.  The title, keywords, MSC, main theorem and the abstract all
   match the published metadata and the zbMATH review, so I regard the *mathematical content* as
   the published content — but I say so as a judgement, not as a verified fact.
2. **Therefore the page numbers in §3 are manuscript pages, not journal pages.**  Do not cite them
   as "JDE 194, p. N".  Cite as: Kukavica, JDE 194 (2003) 39–50; quotations from the author's
   accepted manuscript, `www-bcf.usc.edu/~kukavica/pdf/a28.pdf` (Internet Archive, 2017-06-28),
   manuscript p. N.
3. **The PDF carries no text layer.**  All quotations were transcribed by eye from 170-dpi page
   renders; the OCR file is an index, and its math is mangled.  Anyone re-checking should re-read
   the page images, which are in this folder.
4. **The published BFG §2 text is still unread** (paywalled) — unchanged from the prior seat.  The
   comparison in §4 is against arXiv:1704.05546v4, which is labelled the final version.

---

## 6. Files

| file | what it is |
|---|---|
| `pdf/kukavica-JDE194-2003-a28-authorcopy.pdf` | **the paper** — author's accepted manuscript, 14 pp., sha256 `dd07f7fe…` |
| `raw/a19.pdf`, `png/a19-01.png` | the numbering cross-check (a19 = *Self-similar variables and the complex Ginzburg-Landau equation*, i.e. the file numbers do not track the publication list) |
| `png/ku/p-01.png … p-14.png` | 170-dpi renders — the read-by-eye evidence for every quotation |
| `txt/kukavica-JDE194-2003-OCR.txt` | Apple Vision OCR of all 14 pages (search index only) |
| `raw/kukavica_pdf_index.html` | the author's own archived publication index (`pdf/index.html`, 2017-06-27) |
| `raw/dornsife_pubs.html`, `raw/wb_cdx_*.txt`, `raw/sd_*.{html,pdf}`, `raw/s2.json`, `raw/openalex.json`, `raw/arxiv_*.xml`, … | every route's response as saved |
| `log/round1.tsv`, `log/round2.tsv` | the acquisition log (route, URL, HTTP, bytes, content-type) |
| `scripts/fetch_round1.sh`, `scripts/fetch_round2.sh` | the probes, re-runnable |
| `scripts/verify_provenance.sh` | re-derives the CDX SHA-1 match and the second-snapshot byte identity |
| `scripts/ocr_vision.swift` | the OCR tool (Apple Vision, accurate mode) |
| `scripts/s1_kernel_numbers.py`, `scripts/s1_results.json` | every `[K#]` number |
| `SHA256SUMS` | computed over everything above |

---

# 7. Independent verification pass (second run of this seat, 2026-09-06)

Everything above is round 1.  This section is a second, deliberately decorrelated pass over the
same task: re-acquire the paper from scratch, re-read every load-bearing quotation by eye from
fresh renders, recompute every number with different quadrature, and make one more honest attempt
at the version of record.  New probe log: `logs/probe_vor.tsv`.  New numbers: `scripts/v1_results.json`
(from `scripts/v1_verify_numbers.py`) and `scripts/v2_results.json` (from
`scripts/v2_duhamel_endpoint.py`), both written and run in this pass.  Round 1's manifest is
preserved verbatim as `SHA256SUMS.round1`; `SHA256SUMS` is recomputed over the whole folder.

**Verdict of the pass: round 1 is confirmed in full.  No claim in §§0–6 was found to be wrong.**
Three things are added (§7.3) and one round-1 caveat is re-confirmed as still open (§7.5).

## 7.1 Re-acquisition, addressed differently

Run again from nothing, without reusing round 1's saved bytes:

| route | result this pass |
|---|---|
| `https://dornsife.usc.edu/profile/igor-kukavica/` | 200, 311 396 B. Full publication list in the page text; **zero PDF links** (the only `href`s are USC navigation, `mailto:kukavica@usc.edu`, and the WP oEmbed endpoints). Confirms round 1. |
| `https://sites.usc.edu/kukavica/` | 200, 16 348 B — but the body is a **USC Shibboleth SSO error page** ("Login failed … Stale Request"), not his site. Not a publication route. |
| `http://www-bcf.usc.edu/~kukavica/`, `https://kukavica.usc.edu/` | dead host (curl exit, HTTP 000) |
| `http://kukavica.dreamhosters.com/` | 200, **0 bytes** |
| Wayback CDX, `matchType=prefix`, `www-bcf.usc.edu/~kukavica/` | the same archived directory: `pdf/a01.pdf … pdf/a82.pdf`, `pdf/n47–n54.pdf`, `pdf/index.html`, `webcv.pdf` |

`pdf/index.html` (snapshot `20170627135430`, saved here as `html/wayback-kukavica-pdf-index.html`)
settles the file-numbering question at source.  Verbatim from that page:

> "29. On local uniqueness of solutions of the Navier-Stokes equations with bounded initial data,
> *J. Diff. Eq* **194** (2003), 39-50. … `<a HREF="…S0022039603001530">published version</a>`
> (`<a HREF="a28.pdf">pdf</a>`)"

so **list item 29 → file `a28.pdf`**: the numbering really is offset, exactly as round 1 warned, and
the correct file is the one round 1 used.  The offset is not even uniform: I downloaded `a29.pdf`
too (`pdf/wayback-a29.pdf`, sha256 `d95dced9…4a`) and its title page (`png/verify/a29-01.png`) reads
*Length of vorticity nodal sets for solutions of the 2D Navier-Stokes equations* — list item **28**
(Comm. PDE **28** (2003), 771–793).  So item 29 → `a28`, item 28 → `a29`; the file names do not track
the list, and identification has to be made from the title page, which is what both passes did.  (The same page, item 20, is
"*On the dissipative scale for the Navier-Stokes equation, Indiana Univ. Math. J.* **48** (1999),
1057-1081", `a20.pdf` — the brief's proposed title, a different paper.  Independently: Crossref for
DOI 10.1016/S0022-0396(03)00153-0 returns title "On local uniqueness of weak solutions of the
Navier–Stokes system with bounded initial data", volume 194, page "39-50", published-print 2003-10;
and BFG's own bibliography line, arXiv:1704.05546v4 p.24, reads
"[Ku2] I. Kukavica, *On local uniqueness of weak solutions of the Navier-Stokes system with bounded
initial data*, J. …".)

**Byte identity, three ways.**  I fetched `a28.pdf` through a *different* Wayback snapshot address
than round 1 used (`https://web.archive.org/web/20170628022413id_/…`, vs round 1's `…022417` and
the percent-encoded `…022520`).  Result: 492 472 bytes, sha256
`dd07f7fe9a797843d144504d621471b302d7961d60bb8fc54d30b64bc2e9bf52` — **identical to round 1's
`pdf/kukavica-JDE194-2003-a28-authorcopy.pdf`** (`shasum -a 256` over both files in one call), and
its SHA-1 in base32 is `ZOMHLCABYDYCMM5HK2Q32XSJJNIGJH6B`, matching the CDX digest round 1 recorded.
My copy is kept as `pdf/wayback-a28.pdf`.  `pdfinfo`: 14 pages, 612×792, producer "pstill 1.11",
CreationDate 2003-07-15 20:55:14.

*Timing note for anyone re-running:* the Internet Archive went to its "Temporarily Offline" page
partway through this pass (a CDX query at 15:44 returned that page, minutes after the 15:42 fetch
succeeded).  The route is not reliably available on demand.

## 7.2 Every load-bearing quotation re-read by eye

Fresh 200-dpi renders, `png/verify/v{1,3,4,5,6,7,8,9,11,12}-*.png`, made by me from
`pdf/wayback-a28.pdf` and read as images (the PDF has no text layer; `pdftotext` yields 14 bytes for
the whole file).  Confirmed, verbatim as printed:

* **p.1** — title *On local uniqueness of weak solutions of the Navier-Stokes system with bounded
  initial data*; Igor Kukavica, Department of Mathematics, USC; keywords "Navier-Stokes equations,
  weak solutions, mild solutions, nonuniqueness"; MSC "35Q30, 76D05, 35K15"; running head "Local
  uniqueness with bounded initial data".  (The manuscript's title carries "weak solutions", the same
  as Crossref and as BFG's bibliography; his own web list drops the word.)
* **p.3** — "Fabes, Jones, and Riviere proved in [FJR] existence and uniqueness of local mild
  solutions of the Navier-Stokes system if `n < p ≤ ∞` (`p = ∞` included)", the integral form
  `u(·,t) = ∫₀^t k*(u⊗u)ds + e^{tΔ}u₀`, and `|k(x,t)| ≤ C(n)/(|x|+t^{1/2})^{n+1}`.
* **p.4** — equation (2.1) in **divergence form**: `∂_tu_k − Δu_k + ∂_j(u_ju_k) + ∂_kπ = 0`,
  `∂_ju_j = 0`, with `u₀ ∈ L^∞(ℝⁿ)`.
* **p.5** — **Theorem 2.2**: "*There exists at least one solution `u` of the initial value problem on
  `[0,T)` where `T ≥ 1/C‖u₀‖²_{L^∞(ℝⁿ)}`*", plus the transgalilean characterisation
  `ũ(x,t) = u(x−Φ(t),t)+φ(t)`.  "Here and in the sequel, `C` denotes a generic positive constant
  depending only on `n`."
* **p.6** — `‖π‖_BMO ≤ C‖f‖_{L^∞(ℝⁿ)}`; **Lemma 3.1** exactly as round 1 quotes it, kernel
  `∇G` / `∂_jG`, conclusion `‖u(·,t)‖_{L^∞} ≤ Ct^{1/2}‖f‖_{L^∞(0,T,BMO)}`; and
  `|∇G(x,t)| ≤ C/(|x|+t^{1/2})^{n+1}`.
* **p.7** — the weighted-BMO estimate, **with the local average subtracted**:
  `∫_{ℝⁿ} |f(x) − (1/|B₁|)∫_{B₁}f| / (|x|^{n+1}+1) dx ≤ C‖f‖_BMO`, "where `B₁` is the ball centered
  at 0 with radius 1 (cf. [S, Lemma 3.9] or [St2, p. 141])."  Then "The existence of mild (integral)
  solutions is known.  For completeness sake, we state the result and sketch the iteration
  argument."; the mild equation (3.4) with `∂_jG` on `u_ju_k`, `∂_kG` on `π`, plain `G` on `u_{0k}`;
  (3.5) `π = K_{ij}(u_iu_j)`; **Proposition 3.2** `T ≥ 1/(CM₀²)` with `M₀ = ‖u₀‖_{L^∞}`, and "This
  proposition is proven in [FJR]; however, for completness sake, we include the proof of existence."
  (the typo "completness" is his).
* **p.8** — the iteration `‖u^{(m+1)}‖_{L^∞(S_t)} ≤ Ct^{1/2}‖u^{(m)}‖²_{L^∞(S_t)} + C₀M₀`,
  `‖u^{(1)}‖ ≤ C₀M₀`, "Using induction, we get `‖u^{(m)}‖_{L^∞(S_T)} ≤ 2C₀M₀`" provided
  `0 < T ≤ 1/(CM₀²)` (3.6).
* **p.9, p.11, p.12** — the limiting/de Rham arguments.  Every nonlinear occurrence still carries a
  derivative on the kernel: (3.10) p.11 has `ũ_kũ_j ∂_jG` and `π̃₀ ∂_kG`; the displayed limit on p.12
  is `ũ_k − ∫ũ_{0k}G + ∫₀^t∫ ũ_jũ_k ∂_jG + ∫₀^t∫ π̃₀ ∂_kG = 0`.  The undifferentiated `G` appears
  there only against the **initial datum** and inside the **time-mollifier** terms
  `(1/ε)ψ'(·)·ũ_k·G`, which are linear in `ũ`.  p.11 also carries the sentence round 1 quotes:
  "This part follows [FJR] rather closely except that we work with the heat instead of the Oseen
  kernel."

So round 1's two headline answers stand, now read twice by two passes off byte-identical bytes:
**the kernel on every nonlinear term in Kukavica is `∇G`**, and **the weighted-BMO inequality he uses
is the subtracted one, and it is the only form that appears in his paper.**

## 7.3 Three things at source that round 1 did not record

**(a) Why the kernels differ is a divergence-form fact, not a stylistic one.**  Kukavica's (2.1) is
written in divergence form, `∂_j(u_ju_k)`, so in passing to Duhamel the derivative is integrated by
parts **onto the kernel** — that is the entire provenance of `∂_jG` in (3.4).  BFG's (5),
arXiv:1704.05546v4 p.9, is the **non**-divergence vorticity system,
`∂_tω_j − Δω_j + u_i∂_iω_j = ω_i∂_iu_j`, and their (8) p.9 puts the **undifferentiated** `G` against
**all three** right-hand terms — advection and vortex-stretching alike.  The difference is therefore
not that a mean-zero kernel is unavailable in their setting; it is that they never integrate by
parts.  Anyone repairing §2 should say which of the two moves they are making.

**(b) BFG's own citation for the parenthetical is Kukavica's citation for the real inequality.**
BFG attach `[St]` to the *subtracted* form only, and `[St]` in their bibliography (v4 p.24) is
"E.M. Stein, *Harmonic Analysis: Real-Variable Methods, Orthogonality, and Oscillatory Integrals*" —
the same book Kukavica cites as `[St2, p. 141]` for the subtracted inequality he actually uses.
**The un-subtracted inequality BFG display in their main line carries no citation at all.**

**(c) A typographic caveat when quoting BFG's (8)/(9).**  As printed in v4 (pp.9–10), the third term
of (8) and of (9) reads `ω_i ∂_j u_j` — which is `ω·(∇·u) = 0`.  Two lines later, in the estimate,
the same term is written `ω_i^{(n)} ∂_i u_j^{(n)}`.  The estimate's form is the intended one.  Not
load-bearing, but quote it as printed and note the slip rather than silently correcting it.

## 7.4 The numbers, recomputed with different quadrature

`scripts/v1_verify_numbers.py` uses Gauss–Legendre on the compactified radius `r = s/(1−s)` (plus a
symmetric Cartesian Riemann sum for the signed gradient integral), where round 1's
`s1_kernel_numbers.py` used Gauss–Hermite.  Every rule is reported at two resolutions.
`scripts/v2_duhamel_endpoint.py` re-does the Duhamel integrals after removing the `1/√(t−s)`
endpoint singularity by `s = t − σ²`.

**[V1] kernel means, `n = 3`.**  `∫G(·,τ)dz = 1.0000000000000013 / 1.0000000000000049 /
0.9999999999999993` at `τ = 1, 10⁻², 10⁻⁴`; `∫∂₁G dz = −1.29e−16, −8.07e−16, −5.58e−15`, i.e. zero to
machine precision at all three; `∫|∇G|dz = 1.1283791670955152, 11.283791670955168, 112.83791670955407`
against `2/√(πτ) = 1.1283791670955126, 11.283791670955125, 112.83791670955125`.  **Mean one against
mean zero, reproduced.**

**[V2] the weighted-BMO family** `h_R(x) = (1 − log(1+|x|)/log R)₊` on `ℝ³`, weight `1/(|x|⁴+1)`
(Kukavica's weight), `J₀` un-subtracted, `J₁` with the `B₁`-average subtracted, at `n = 2400`:

| `R` | `avg_{B₁}h_R` | `J₀` | `J₀·log R` | `J₁` | `J₁·log R` |
|---|---|---|---|---|---|
| 1e2 | 0.8799260385 | 10.1491838721 | 46.7387190 | 2.2134627486 | 10.1933727 |
| 1e4 | 0.9399630192 | 12.0398801727 | 110.8913944 | 1.1203073373 | 10.3184119 |
| 1e8 | 0.9699815096 | 12.9987360905 | 239.4455676 | 0.5602218641 | 10.3196681 |
| 1e16 | 0.9849907548 | 13.4782322449 | 496.5564264 | 0.2801109321 | 10.3196681 |
| 1e32 | 0.9924953774 | 13.7179803221 | 1010.7781439 | 0.1400554660 | 10.3196681 |
| 1e64 | 0.9962476887 | 13.8378543607 | 2039.2215789 | 0.0700277330 | 10.3196681 |

`J₀·log R` grows linearly (46.7 → 2039.2, a factor 43.6 across a factor-32 range of `log R`); `J₁·log R`
is flat at `10.31967`.  Agreement with round 1's table to 6–7 significant figures at every `R`.
**The subtraction is the whole inequality on this family, and the un-subtracted form is false.**

**[V3] the two Duhamel scalings** (desingularised, `scripts/v2_results.json`):
`∫₀^t∫(|z|+√(t−s))^{-4}dz ds / √t = 8.377580569016, 8.377579999536, 8.377580367455`
at `t = 0.25, 1, 4`, against `8π/3 = 8.377580409572781`; and
`∫₀^t∫√(t−s)(|z|+√(t−s))^{-4}dz ds / t = 4.188790204786723, 4.188790204784745, 4.188790204786042`
against `4π/3 = 4.188790204786390`.  (The un-desingularised V3 block in `v1_results.json` gives
`8.3715069` for the first — a 0.07 % endpoint-quadrature error, not a discrepancy; the radial closed
form `∫_{ℝ³}(|z|+a)^{-4}dz = 4π/(3a)` is checked separately at `a = 0.1, 1, 3`, agreeing to 12-13 significant figures.)
**Kukavica's `Ct^{1/2}` and BFG's `ct` are exactly these two, and the exponent is the kernel's
fingerprint** — round 1's structural point, reproduced.

**[V4]** `sup_x (|x|+1)⁴/(|x|⁴+1) = 8.0` attained at `|x| = 1`, `inf = 1.0` at `|x| = 0`: the two
weights are equivalent, so nothing in the comparison turns on which is written.

## 7.5 One observation, flagged as NOT established

Since `∇·u = 0` and `∇·ω = 0`, both right-hand terms of BFG's (5) are themselves divergence-form:
`u_i∂_iω_j = ∂_i(u_iω_j)` and `ω_i∂_iu_j = ∂_i(ω_iu_j)`.  So (8) *could* be rewritten with `∂_iG`
against `u⊗ω`, a mean-zero kernel, and Kukavica's Lemma 3.1 would then hand back a `t^{1/2}` clock
with the subtracted inequality legitimately in play.  What that costs is a hypothesis on the
**product**: Lemma 3.1 needs its argument in `L^∞(0,T,BMO)`, and here the argument is `u_iω_j`.
`BMO` is not an algebra and `L^∞·BMO ⊄ BMO`, and `u` itself is not controlled in `L^∞` by
`ω₀ ∈ L²∩L^∞` without exactly the logarithm at issue.  **I did not attempt to close this and do not
claim it can be closed.**  I record it because it relocates the difficulty from the kernel to the
product, which is where a repair attempt should start.

## 7.6 Still not obtained: the Elsevier version of record

Fresh attempt this pass, logged in `logs/probe_vor.tsv` (`scripts/probe_vor.sh`, re-runnable),
responses in `raw2/`:

| route | HTTP | note |
|---|---|---|
| ScienceDirect landing / `/pdf` / `/pdfft` | **403** ×3 | 1.2 MB Cloudflare block page each; block code `CPE9` appears in the body |
| `https://doi.org/10.1016/S0022-0396(03)00153-0` | 200 | 2 712 B Elsevier redirect interstitial (`articleSelectSinglePerm`), not the article |
| Elsevier article API (`api.elsevier.com/content/article/doi/…`) | **406** | needs a key |
| CORE v3 search API | **429** | rate-limited without a key |
| OpenAIRE (`api.openaire.eu`, 1 hit) | 200 | **no repository copy exists**: every `<webresource>` is the DOI or `zbmath.org/2003922`; `<license>` = "Elsevier Non-Commercial" |
| Unpaywall / Semantic Scholar / OpenAlex | 200 | all three name the **same single** OA location, the ScienceDirect PDF; `oa_status: bronze`, `has_repository_copy: false` |
| BASE | 200 | no full text |
| ResearchGate | **403** | |
| scholar.archive.org | 200 | 4 428 B JS shell |
| `api.fatcat.wiki` | connection error | unreachable from here, third pass in a row |
| WebFetch on the ScienceDirect landing page | **403** | body not retrieved |
| Sci-Hub | **not attempted — disallowed** | |

So round 1's caveat 5.1/5.2 stands unchanged: what we hold is the **author's accepted manuscript**
(created 2003-07-15, ~7 weeks before the online publication date 2003-09-03), not pp. 39–50 as
printed.  Cite it as round 1 says, with manuscript page numbers.  Crossref does record an Elsevier
open-access user licence for this DOI starting 2013-07-17 with `content-version: vor`, i.e. the
version of record *is* meant to be freely readable on ScienceDirect; the obstacle is a bot block,
not a paywall, and a human with a browser would very likely get it in one click.  That is the
cheapest remaining route and it is an operator action, not one I can take from here.

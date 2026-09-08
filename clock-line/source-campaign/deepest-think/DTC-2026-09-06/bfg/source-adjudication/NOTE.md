# BFG source-adjudication — Kukavica 2003, Giga–Inui–Matsui 1998, the ARMA version, citing work, errata

Seat: `bfg/source-adjudication` (DTC-2026-09-06).  Task: adjudicate **at source** whether
Bradshaw–Farhat–Grujić Theorem 10 (the log-free vorticity clock) is (a) proved by a valid
published argument, (b) unsupported by its displayed proof but not contradicted, or
(c) contradicted by a published result.

**Register.** Nothing below is peer-confirmed outside the estate.  Every number tagged `[N#]`
comes from `scripts/adjudication_numbers.py`, which I wrote and ran; its raw output is
`scripts/adjudication_numbers.out`.  Sources were read at source (PDF or publisher page) except
where explicitly marked NOT OBTAINED.

---

## 0. Verdict

**(b), sharpened.**  Theorem 10 as stated is **not supported by the argument BFG display**, and it
is **not contradicted by any published Navier–Stokes result I could find**.  Three things are new
from this pass:

1. The prior art BFG cite for the log-free `p = ∞` case — Giga–Inui–Matsui 1998 and Kukavica 2003 —
   is, **by BFG's own sentence and by GIM's own theorem**, the *velocity–pressure* result
   `T₀ ≥ C‖u₀‖_∞^{-2}`.  It is log-free and energy-free, but it is **not** the vorticity statement,
   and it does not imply it: the best vorticity clock I can derive from it is
   `M₀T ≳ Re_ω^{-1/3}`, a **power-law** penalty, far weaker than either BFG's `Re⁰` or the
   campaign's `1/log Re` **[N3]**.
2. There is a **structural reason** the transplant from Kukavica/GIM to vorticity is not automatic,
   and it is exactly at the disputed step: in the velocity/divergence-form mild equation the kernel
   is `∇G`, which has **integral zero**, so the local average may be subtracted for free; in BFG's
   (8) the vortex-stretching term carries the **undifferentiated** kernel `G`, whose integral is
   **one**, against a **non-negative** integrand `|ω||∇u|` **[N2]**.  Subtracting a local average is
   precisely what is *not* available there — which is why their parenthetical "(in general, one has
   to subtract a local average)" is not a harmless aside in this setting even though it is harmless
   in the setting they are modifying.
3. Restoring the omitted term changes the fixed-point condition by exactly one factor `Λ`:
   `T ≍ 1/(4c₀c₁ M Λ)` **[N5]**.  `Λ = 1` is BFG's Theorem 8; `Λ = 1 + log₊(ℓ/√(ν/M))` is the
   campaign's logarithmic clock.  The two statements differ by nothing else.

**Nobody has re-proved Theorem 8/10, nobody has contested it, and there is no erratum.**  Two
published papers (Grujić–Xu, JMFM 2024 and Ann. PDE 2025) **use it verbatim as a black box** — the
`c(M)` time window in their escape-time machinery *is* BFG's Theorem 10 — and one published paper
(Albritton–Bradshaw, *Nonlinearity* 2022, with a BFG co-author) contests BFG's **headline**
scaling-gap claim while leaving §2 untouched.

**A live route to (c) exists and is not closed.**  BFG Thm 10 is exactly ν-uniform **[N1]**; the
Euler analogue is false along the Kim–Jeong / Bang–Cheskidov norm-inflation families, where
`M₀ × (doubling time) → 0` **[N6]**.  A justified inviscid limit on the Euler existence interval
would therefore contradict Thm 10.  **I have not justified that limit and do not assert it.**

---

## 1. The table

| # | Source | What it proves | Time scale | Energy / log dependence | Verbatim quote + page |
|---|---|---|---|---|---|
| 1 | **Giga–Inui–Matsui**, *On the Cauchy problem for the NSE with nondecaying initial data*, Hokkaido Univ. Preprint Series in Math. **410** (May 1998), DOI 10.14943/83556 (scanned; obtained at source, `pdf/gim-1998-re410.pdf`) | Local existence + uniqueness of a mild **velocity** solution for `u₀ ∈ BUC` (or `L^∞`), `R^n`, `n ≥ 2`, `ν = 1`. Data may be **nondecaying** (constants, periodic). | `T₀ ≥ C‖u₀‖_∞^{-2}` | **none of either.** No decay hypothesis at all, hence no energy; no log. Riesz-unboundedness on `L^∞` is bypassed via Carpio's `H¹`(Hardy) bound on `√t ∇e^{tΔ}ℙ`, *not* via a weighted BMO integral. | **Thm 1, p.2**: "*(i) For `u₀ ∈ BUC` there exists `T₀ > 0` and a unique solution `u = u(t)` of (INT) such that `u ∈ C([0,T₀]; BUC)`.*" — **Remark (i), p.3**: "*For a lower estimate for `T₀ > 0`, we get `T₀ ≥ C ‖u₀‖_∞^{-2}` with a numerical constant `C > 0`. This follows from the way of construction in §3.*" — **Remark (vi), p.4**: "*sup_{0≤t≤T₀} ‖u(t)‖_∞ ≤ 2‖u₀‖_∞, sup_{0≤t≤T₀} t^{1/2}‖∇u(t)‖_∞ ≤ C‖u₀‖_∞ with a constant `C` depending only on `n`.*" — **Remark (iii), p.3**: "*If we replace `BUC` by `L^∞` in Theorem 1, we obtain same results as Theorem 1.*" (preprint pp.1–4 = PDF pp.3–6) |
| 2 | **Kukavica**, *On local uniqueness of weak solutions of the Navier–Stokes system with bounded initial data*, **J. Differential Equations 194 (2003) 39–50**, DOI 10.1016/S0022-0396(03)00153-0 | **NOT OBTAINED** (see §2). Established at source about it: **velocity** formulation, `u₀ ∈ L^∞(R^n)`; local existence of a *mild* solution is taken as **known** input; the paper's own theorem is the characterisation of non-uniqueness of **weak** solutions up to a drift. | — (its own theorem carries no existence-time claim) | — | **zbMATH review Zbl 1050.35070** (Klaus Deckelnick), read at source via the zbMATH API: "*where the initial function `u₀` belongs to `L^∞(R^n, R^n)`. It is known that this problem has a unique local mild solution, while weak solutions are not unique. The main result of the paper states that `ũ` is a weak solution on `[0,T)` if and only if there exists a function `φ ∈ L^∞_loc([0,T), R^n)` with `lim_{t→0⁺} φ(t) = 0` such that `ũ(x,t) = u(x − Φ(t), t) + φ(t)`.*" |
| 3 | **BFG**, arXiv:1704.05546**v4** (7 Sep 2018; comments field: "*final version; to appear in ARMA*"), **Thm 8, p.8** | Unique mild **vorticity** solution in `C_w([0,T],L^∞)`, `ν = 1`. Existence time only. | `T ≥ (1/c)(1/‖ω₀‖_∞)`, `c` absolute | **claims none of either**, explicitly | "*Let the initial datum `ω₀` be in `L²∩L^∞`. Then there exists a unique mild solution `ω` in `C_w([0,T], L^∞)` where `T ≥ (1/c)(1/‖ω₀‖_∞)` for an absolute constant `c > 0`.*" — and p.8: "*this is a 'soft assumption', i.e., there will be no quantitative dependence on `‖ω₀‖₂` in the proof.*" |
| 4 | **BFG**, arXiv:1704.05546v4, **Thm 10, p.11** | Same **plus** complex-analytic extension **plus** the growth bound | `T ≥ (1/c(M))(1/‖ω₀‖_∞)`, growth factor `M` | claims none of either; **ν-uniform** by rescaling **[N1c]** | "*… `Ω_t = {x + iy ∈ C³ : |y| < (1/√(c(M)))√t}`; moreover, `‖ω(t)‖_{L^∞(Ω_t)} ≤ M‖ω₀‖_∞`.*" |
| 5 | **BFG published**: *Arch. Ration. Mech. Anal.* **231** (2019) **1983–2005**, DOI 10.1007/s00205-018-1314-5, online 21 Sep 2018 | **Full text NOT OBTAINED** (paywalled; Springer page read at source for metadata + full reference list). | — | — | Springer landing page: "*Published: 21 September 2018 · Volume 231, pages 1983–2005 (2019)*". **Note the campaign's cite "ARMA 233 (2019)" is wrong — it is volume 231.** No correction notice on the page; Crossref shows `update-to = None`, `updated-by = None`, `relation = {}`. arXiv v4 posted **14 days before** online publication and is labelled the final version, so the presumption is that §2 is identical — **but I could not verify the published §2 text at source.** |
| 6 | **Grujić–Xu**, *Asymptotic Criticality of the NS Regularity Problem*, arXiv:1911.00974 → **J. Math. Fluid Mech. 26 (2024)**, DOI 10.1007/s00021-024-00888-x (open access) | **Uses** Thm 10 as a quoted black box; does **not** re-prove it. Their own level-`k` analogue (Thm 2.5) is proved, and **does** carry an `L^p` norm in the existence time. | quotes `T ≥ 1/(c(M)‖ω₀‖_∞)`; their escape-time window is `s ∈ (t + 1/(4c(M)‖ω(t)‖_∞), t + 1/(c(M)‖ω(t)‖_∞))` | quoted statement: none. Their **own** Thm 2.5 time `T*` depends on `‖ω₀‖_p` **and** `‖D^kω₀‖_∞` | **Thm 2.1 (Guberović [12] and Bradshaw et al. [3])**, arXiv p.6: "*Let the initial datum `u₀ ∈ L^∞` (resp. `ω₀ ∈ L^∞ ∩ L¹`). Then, for any `M > 1`, there exists a constant `c(M)` such that there is a unique mild solution `u` (resp. `ω`) in `C_w([0,T],L^∞)` where `T ≥ 1/(c(M)²‖u₀‖²_∞)` (resp. `T ≥ 1/(c(M)‖ω₀‖_∞)`) …*" — note they write the hypothesis as `L^∞ ∩ L¹`, BFG write `L² ∩ L^∞`. **Thm 2.5**, arXiv p.12: "*Assume the initial value `ω₀ ∈ L^∞(R³) ∩ L^p(R³)` where `1 ≤ p < 3`. Fix `k ∈ N`, `M > 1` … `T* = C(M)·min{ 2^{-k}‖ω₀‖_p^{k/(k+d/p)}·(‖D^kω₀‖_∞^{(d/p)/(k+d/p)} + ‖ω₀‖_p)^{-1}, … }*" |
| 7 | **Grujić–Xu**, *Time-Global Regularity of the NS System with Hyper-Dissipation: Turbulent Scenario*, arXiv:2012.05692 → **Ann. PDE 11:9 (2025)**, DOI 10.1007/s40818-025-00199-y | Same use, same window, same attribution | same | same | **Thm 2.12 (Farhat et al. [14] and Bradshaw et al. [4])**, arXiv p.15: "*resp. `s = s(t) ∈ (t + 1/(4c(M)‖ω(t)‖_∞), t + 1/(c(M)‖ω(t)‖_∞))`*" |
| 8 | **Albritton–Bradshaw**, *Remarks on sparseness and regularity of Navier–Stokes solutions*, arXiv:2110.02187v2 → ***Nonlinearity* 35 (2022)**, DOI 10.1088/1361-6544/ac62de | **Contests BFG's headline claim** (the algebraic reduction of the scaling gap), **not** Thm 8/10. Gives an independent, analyticity-free proof of the sparseness⇒regularity implication at the **velocity** `L^p` level, `d < p ≤ ∞`. | their local theory: `‖u₀‖_{L^p} T̄^{(1/2)(1−d/p)} = c_p`; at `p = ∞`, `T̄ = (c_∞/‖u₀‖_∞)²` | log-free, energy-free — again the **velocity** clock, matching GIM | Abstract: "*Second, we analyze the claims in [BFG19, GX19] that a priori estimates on the sparseness of the vorticity and higher velocity derivatives reduce the 'scaling gap' in the regularity problem.*" — §1.2, p.5: "*We offer a different interpretation of [BFG19, GX19] which suggests that, in a certain reasonable sense, the 'scaling gap' is not improved beyond the energy class.*" — §1, p.2: "*we show for a broad class of concrete examples of vector fields that the framework introduced in [GX19] does not rule out singularities beyond those ruled out by membership in `L^∞_t L²_x`.*" — p.4, (1.11): "*`‖u₀‖_{L^p} T̄^{(1/2)(1−d/p)} = c_p`*". **No vorticity-`L^∞` analogue is given, and §2's local theory is never questioned.** |
| 9 | **Farhat–Grujić–Leitmeyer**, *Erratum to: The Space `B^{-1}_{∞,∞}`, Volumetric Sparseness, and 3D NSE*, **J. Math. Fluid Mech. 19 (2017) 525–527**, DOI 10.1007/s00021-016-0295-0 | An erratum **exists for the velocity companion paper [FGL]**, not for BFG. Full text NOT OBTAINED (paywalled; Springer landing page read at source). Its only two references are Bahouri–Chemin–Danchin and Iyer–Kiselev–Xu, i.e. it concerns the Besov/mix-norm sparseness side, **not** the `L^∞` local theory. | — | — | Springer landing page: "*Erratum · Published: 28 September 2016 · Volume 19, pages 525–527 (2017)*"; "*The Original Article was published on 06 September 2016*". BFG's reference list cites only the **original** [FGL], not the erratum. |
| 10 | **Euler side** (Kim–Jeong JFA 283 (2022) 109673 = arXiv:2111.14078; Bang–Cheskidov arXiv:2605.16502) — as established by the `literature-short-time` seat, not re-read here | Norm inflation for 3D **Euler** from small smooth axisymmetric no-swirl data | `M₀ × T_double → 0` along the family **[N6]** | — | See `../sharp/literature-short-time/NOTE.md` rows 4 and 6. Relevant because BFG Thm 10 is **ν-uniform [N1c]**. |

---

## 2. Acquisition record — what I got and what I could not

**Obtained at source:**

* `pdf/gim-1998-re410.pdf` — Hokkaido preprint #410, direct from `eprints.lib.hokudai.ac.jp`.
  It is a **scanned (CCITT bilevel) PDF**; `pdftotext` recovers only the HUSCAP cover sheet, so
  pages 1–4 of the preprint (PDF pages 3–6) were read as rendered images and transcribed by hand.
  There is no OCR tool on this machine (`tesseract`/`ocrmypdf` absent).
* Springer landing pages for the ARMA article, the JMFM 2024 article, and the FGL erratum
  (metadata, abstracts, full reference lists) — read at source in the browser pane.
* `pdf/cite-*.pdf` + `txt/cite-*.txt` — the open-access citing works (arXiv/Springer).
* `pdf/kukavica-vicol-KV01.pdf` — Kukavica–Vicol, *J. Dyn. Diff. Eq.* 20 (2008) 719–732, the direct
  follow-up by the same author, from `cims.nyu.edu/~vicol/`.
* `pdf/kukavica-cv.pdf` — Kukavica's CV (item 29 confirms the 2003 entry).
* zbMATH review Zbl 1050.35070 via the zbMATH Open API (`html/zb-kukavica.json`).
* Crossref, OpenAlex, Unpaywall, Semantic Scholar metadata (`html/*.json`).

**NOT OBTAINED — Kukavica, JDE 194 (2003) 39–50.**  Routes tried, all failed:
ScienceDirect article page and `/pdf` endpoint (HTTP 403 to WebFetch; Cloudflare bot block
`CPE00001` in the browser pane); Elsevier text-mining API (400, needs a key); Unpaywall
(`oa_status: bronze`, the **only** OA location is that same ScienceDirect PDF);
Semantic Scholar (same URL, abstract elided by publisher); OpenAlex (same);
CORE (needs a key); fatcat/scholar.archive.org (`api.fatcat.wiki` unreachable from here — connection
timed out); the author's USC publications page and CV (no PDF links; his `dreamhosters.com`
root is a 5-byte stub); Google-style searches (they repeatedly mis-resolve to the 2008
Kukavica–Vicol `BMO^{-1}` paper, which is a **different** paper).
**Consequence:** I cannot quote Kukavica's §-level text, so I cannot state verbatim whether his
exposition uses `∫|f|/(|x|+1)^{n+1} dx ≤ c‖f‖_BMO` with or without the local average.  What I
*can* establish at source is below.

**NOT OBTAINED — the ARMA published full text** (paywalled, $39.95).  See table row 5 for what is
established instead.

**NOT OBTAINED — the FGL erratum full text** (paywalled).  See row 9.

---

## 3. Question (1): what Kukavica actually provides, and why the transplant is not automatic

Three source-level facts, none of which require the 2003 full text:

**(1a) BFG themselves say the cited prior art is velocity–pressure.**  arXiv v4, p.8, verbatim:

> "Since the Riesz transforms are not bounded on `L^∞`, to obtain an estimate in the case `p = ∞`
> without a logarithmic correction requires a different argument (on the real level), see, e.g.,
> **[GIM, Ku2]** and [Gu] in the real and the complex setting, respectively, **in the realm of the
> velocity-pressure formulation of the 3D NSE.**
> **The `L^∞`-argument within the vorticity-velocity description requires a modification**, and we
> present a sketch here, including an estimate on the vortex-stretching term, for completeness.
> What follows is a modification of the exposition given in [Ku2]."

So the log-free precedent is, by their own statement, for **velocity**.  What is transplanted is an
*exposition*, and the transplant is explicitly a *modification*, sketched, "for completeness".

**(1b) The velocity precedent's time scale is `T ≳ ‖u₀‖_∞^{-2}`, and it is log- and energy-free.**
That is GIM Remark (i), quoted verbatim in row 1, obtained at source.  It is independently
re-stated in 2022 by Albritton–Bradshaw (row 8, eq. (1.11) at `p = ∞`).  Both are **velocity**.

**(1c) The velocity route is structurally immune to the disputed step; the vorticity route is not.**
In the velocity mild equation the nonlinear term carries the kernel `∇e^{(t−s)Δ}ℙ` — GIM's (INT),
Albritton–Bradshaw's (1.13), Kukavica's setting.  `∫_{R³} ∇G(z,τ) dz = 0` (odd), while
`∫_{R³} G(z,τ) dz = 1` and `∫_{R³}|∇G(z,τ)| dz = 2/√(πτ)` **[N2b, N2c]**:

```
tau=1        int G dz = 1.000000000000   int |grad G| dz = 1.12837917   = 2/sqrt(pi tau)
tau=1e-2     int G dz = 1.000000000000   int |grad G| dz = 11.28379167  = 2/sqrt(pi tau)
tau=1e-4     int G dz = 1.000000000000   int |grad G| dz = 112.83791671 = 2/sqrt(pi tau)
```

A kernel of mean zero annihilates constants, so in the velocity setting **the local average may be
subtracted for free** and the honest inequality `∫|f − f_B|/(|x|+1)^4 ≤ c‖f‖_BMO` is exactly what one
needs.  BFG's third term in (8) is `∫₀^t∫ G(x−y,t−s) ω_i ∂_i u_j` — the **undifferentiated** kernel,
mean one, applied to the **non-negative** integrand `|ω||∇u|`.  There is nothing to subtract.  This
is, structurally, the whole content of the phrase "requires a modification", and it is exactly where
the modification is not carried out.

**(1d) The inequality as displayed is false, third independent family [N4].**  I built my own
compactly supported family `h_R(x) = (1 − log(1+|x|)/log R)₊` (so `‖h_R‖_BMO ≤ ‖log(1+|x|)‖_BMO/log R
≤ C/log R`, and `h_R` is in the closure of test functions in every uniformly-local `L^p`, i.e. it
satisfies BFG's decay proviso):

```
R=1e2   int h_R/(1+|x|)^4 = 2.548370   1/log R = 2.171e-01   I*log R =  11.74
R=1e4   int h_R/(1+|x|)^4 = 3.355141   1/log R = 1.086e-01   I*log R =  30.90
R=1e8   int h_R/(1+|x|)^4 = 3.771897   1/log R = 5.429e-02   I*log R =  69.48
R=1e16  int h_R/(1+|x|)^4 = 3.980344   1/log R = 2.714e-02   I*log R = 146.64
R=1e32  int h_R/(1+|x|)^4 = 4.084567   1/log R = 1.357e-02   I*log R = 300.96
R=1e64  int h_R/(1+|x|)^4 = 4.136679   1/log R = 6.786e-03   I*log R = 609.60
```

The left side rises to `4π/3 = 4.188790` while the BMO seminorm falls like `1/log R`; the ratio grows
linearly in `log R`.  This is a **third** family, built independently of `refuter-correctness`'s
shell field and `refuter-scope-priorart`'s `f_L`, and it agrees with both.  BFG's decay proviso does
not rescue the inequality.

**(1e) The cost is exactly one factor [N5].**  With the local average restored the iteration reads
`sup‖ω^{(n+1)}‖_∞ ≤ c₀M + c₁ t Λ (sup‖ω^{(n)}‖_∞)²`; the ball `{‖ω‖ ≤ 2c₀M}` is invariant iff

```
t <= 1/(4*Lambda*M*c0*c1)
```

`Λ = 1` reproduces Theorem 8 verbatim.  `Λ = 1 + log₊(ℓ/√(ν/M))` is a logarithmic clock.  **Nothing
else in the chain changes.**

**(1f) What clock the cited prior art *does* buy for vorticity [N3].**  Biot–Savart splitting at
radius `ρ` gives `‖u₀‖_∞ ≲ ρM₀ + ρ^{-1/2}‖ω₀‖₂`, optimised at `ρ* = (‖ω₀‖₂/M₀)^{2/3}·2^{1/3}/2`, so
`‖u₀‖_∞ ≲ (3·2^{1/3}/2) M₀^{1/3}‖ω₀‖₂^{2/3}`.  Feeding that into GIM's `T ≥ C‖u₀‖_∞^{-2}`:

```
M_0 * T  >~  (2*2^{1/3}/9) * M_0^{1/3} * ||om_0||_2^{-4/3}  =  (2*2^{1/3}/9) * Re_om^{-1/3}
Re_om := ||om_0||_2^4 / ||om_0||_inf     (dimensionless: scaling ratio = 1 exactly)
```

So the honest clock inherited from the prior art BFG cite is `M₀T ≳ Re^{-1/3}` — a **power-law**
degradation.  BFG's Theorem 8 asserts `Re⁰`.  The campaign's S1 asserts `1/(1+log Re)`.  The gap
between the cited precedent and the asserted theorem is therefore not a technicality; it is the
entire content of §2.

---

## 4. Question (3): ARMA vs arXiv v4

* ARMA **231** (2019) **1983–2005**, online 21 Sep 2018 (Springer page, read at source).
  The campaign's "ARMA 233 (2019)" is **wrong**; correct it to 231.
* arXiv v4 was posted **7 Sep 2018**, comments "*final version; to appear in ARMA*" — 14 days before
  online publication.  v1/v2 (Apr 2017) are 19 KB; v3 (Jan 2018) and v4 are ~1 MB.
* Crossref: `update-to = None`, `updated-by = None`, `relation = {}` → **no registered correction**.
* Springer's landing page shows no correction notice.
* **I could not read the published §2.**  The strongest available proxy that the *statement* survived
  refereeing unchanged is Grujić–Xu's 2024 (JMFM) and 2025 (Ann. PDE) restatements by the same senior
  author, which reproduce `T ≥ 1/(c(M)‖ω₀‖_∞)`, the `c(M)` window and the `M‖ω₀‖_∞` bound exactly
  (row 6, row 7).  One discrepancy: they write the hypothesis as `ω₀ ∈ L^∞ ∩ L¹`, BFG write
  `ω₀ ∈ L² ∩ L^∞`.

---

## 5. Question (4): who cites Thm 8 / Thm 10, and what they take from it

OpenAlex lists 18 citing works for `W2608160632` (`html/oa-bfg-citedby.json`); Semantic Scholar lists
1.  I fetched every open-access one (`scripts/fetch_citing.sh`) and grepped each for the citation.

**Uses it, load-bearingly, without re-proving it — 2 papers, both Grujić–Xu:**

* *Asymptotic Criticality* (JMFM 2024): **Theorem 2.1** is BFG Thm 10 quoted, attributed
  "(Guberović [12] and Bradshaw et al. [3])"; **Theorem 3.4** then runs the escape-time argument on the
  window `s ∈ (t + 1/(4c(M)‖ω(t)‖_∞), t + 1/(c(M)‖ω(t)‖_∞))`.  That window **is** Thm 10's time scale.
  Their own proved level-`k` result (Thm 2.5) is stated only for `ω₀ ∈ L^∞ ∩ L^p` and its `T*`
  **does** contain `‖ω₀‖_p`.  *(My own observation, flagged as extrapolation, not their claim: the
  `k = 0` specialisation of their (2.18) reads `T* ≈ C(M)/(‖ω₀‖_∞ + ‖ω₀‖_p)` — i.e. the only
  fully-proved analogue in this literature carries an `L^p` norm exactly where BFG's does not.)*
* *Time-Global Regularity … Hyper-Dissipation* (Ann. PDE 11:9, 2025): **Theorem 2.12**, identical
  window, identical attribution.

**Cites the paper but not §2 — 3 papers examined:** Do–Farhat–Grujić–Xu, *Oscillations and
integrability of the vorticity* (arXiv:1801.09040 → IUMJ 69 (2020)); Bradshaw–Tsai, *On the local
pressure expansion* (arXiv:2001.11526); the axisymmetric NSE review (arXiv:2101.04905), which mentions
only the scaling-gap headline.

**Contests the paper — 1 published paper, but not §2:** Albritton–Bradshaw, *Nonlinearity* 35 (2022)
(row 8).  This is significant context — a BFG co-author publicly arguing the headline conclusion does
not do what it claims — but the disputed local-existence ingredient is not part of their critique, and
their own replacement local theory is again the **velocity** clock `T̄ = (c_∞/‖u₀‖_∞)²`.

**Re-proves Thm 8 or Thm 10 — none found.**
**Contests Thm 8 or Thm 10 — none found.**

Coverage caveat: arXiv has no full-text search API, so my citation net is OpenAlex + Semantic Scholar
+ targeted web search.  Three OA-listed items resolved only to publisher HTML behind IOP/Springer
paywalls and were not read: *Remarks on the separation of Navier–Stokes flows* (Nonlinearity 2024),
*Regularity, Uniqueness and the Relative Size of Small and Large Scales in SQG Flows* (JMFM 2025),
*A regularity criterion for 3D NSE in 'dynamically restricted' local Morrey spaces* (2021).

---

## 6. Question (5): errata

* **BFG (ARMA 231, 2019): no erratum, no corrigendum.**  Crossref `update-to`/`updated-by` both null;
  no notice on the Springer page; no citing work refers to one.
* **The companion paper does have one**: Farhat–Grujić–Leitmeyer, *Erratum to: The Space
  `B^{-1}_{∞,∞}`, Volumetric Sparseness, and 3D NSE*, JMFM **19** (2017) 525–527,
  DOI 10.1007/s00021-016-0295-0 (published 28 Sep 2016).  Full text paywalled; from its landing page
  its only two references are Bahouri–Chemin–Danchin and Iyer–Kiselev–Xu, so it concerns the
  Besov/mix-norm side rather than the `L^∞` local theory.  BFG's own reference list cites the
  original [FGL] and **not** the erratum, even though the erratum predates BFG v1 by seven months.
  Bibliographic, not load-bearing.

---

## 7. Adjudication (one paragraph)

**Theorem 10 as stated is (b): unsupported by its displayed proof, and not contradicted by any
published Navier–Stokes result.**  BFG present §2 as a *sketch* that *modifies* an exposition whose
original — Kukavica 2003, and behind it Giga–Inui–Matsui 1998 — is a **velocity–pressure** result with
the different, weaker-in-currency clock `T₀ ≥ C‖u₀‖_∞^{-2}`, which I verified verbatim at source in
GIM and which BFG themselves label as belonging to "the realm of the velocity-pressure formulation".
The single inequality that carries the modification, `∫|f|/(|x|+1)^4 dx ≤ c‖f‖_BMO`, is false with an
absolute constant even under their stated decay proviso — three independent families now say so, mine
included **[N4]** — and the reason it is harmless in the velocity setting and not in theirs is
structural and checkable: the velocity mild equation's kernel `∇G` has integral zero and therefore
tolerates subtracting the local average, while BFG's vortex-stretching term carries the
undifferentiated kernel `G`, integral one, against a non-negative integrand **[N2]**.  Restoring the
omitted term costs exactly one factor and turns Theorem 8 into a logarithmic clock **[N5]**, and the
best clock actually inherited from the cited prior art is the much weaker `M₀T ≳ Re^{-1/3}` **[N3]**.
Against that, the theorem has stood in a refereed ARMA paper for seven years with no erratum, and two
further refereed papers (JMFM 2024, Ann. PDE 2025) use it as a black box in a load-bearing way, while
the one published critique of BFG (Albritton–Bradshaw, *Nonlinearity* 2022, co-authored by BFG's own
first author) attacks the headline scaling-gap claim and leaves §2 alone.  **No published result
contradicts it**, and I have no counterexample at the level of an actual Navier–Stokes solution, so
(c) is not reached.  The nearest live route to (c) — and the thing I would put next on the docket —
is that Theorem 10 is exactly ν-uniform **[N1]** while its Euler analogue is false along the
Kim–Jeong / Bang–Cheskidov norm-inflation families **[N6]**; if the inviscid limit can be justified on
the Euler existence interval for their smooth compactly supported data, that is a contradiction at
the level of actual Navier–Stokes solutions.  **I have not justified that limit and do not assert
it.**  Operationally, for the campaign: the honest public statement remains that S1's logarithmic
clock is in open conflict with a published theorem whose only available proof has an identified gap,
and that neither side may be asserted as settled.

---

## Files

* `pdf/`, `txt/` — sources as downloaded and their `pdftotext -layout` extractions.
* `html/` — API responses (Crossref, OpenAlex, Unpaywall, Semantic Scholar, zbMATH) as saved.
* `scripts/adjudication_numbers.py` — every `[N#]` number; output in `scripts/adjudication_numbers.out`.
* `scripts/fetch_citing.sh` — the citing-works download.
* `scripts/find_oa.py` — OA-location lookup for the citing works.
* `SHA256SUMS` — hashes of everything above.

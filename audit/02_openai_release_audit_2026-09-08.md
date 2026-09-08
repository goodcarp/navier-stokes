# Audit of the 2026-09-08 release: forced Navier-Stokes blowup and unforced Euler blowup

Date: 2026-09-08. Seat: HOLD UP SHIPS / forced route.
Materials audited: `navier-stokes.pdf` / `.txt` (165 pp), `euler.pdf` / `.txt` (45 pp), and the Lean
repository `NavierStokesAndEuler` (Lean 4.34.0-rc2, 2486 `.lean` files), all as received.
Instruments: `../theorems/01_two_fences.md` (Lemma A, Lemma B, the section 4 checklist, Fefferman
(4)-(11) quoted verbatim) and `../theorems/02_forced_axisymmetric_on_axis_typeII.md`.
Line citations are to the `.txt` extractions; page numbers are the printed pages.

---

## 1. Verdict

The Navier-Stokes paper claims exactly Fefferman's forced alternatives (C) and (D) and nothing
weaker or stronger: for every `nu > 0` a force `f in C_c^inf(R^3 x (0,inf))`, a compact `K`, and
smooth `(u,p)` on `R^3 x [0,1)` from rest with `supp u(.,t) u supp p(.,t) subset K`,
`sup_{0<=t<1} ||u(t)||_{L^2} < inf` and `limsup_{t->1} ||u(t)||_{L^inf} = inf`, and "consequently"
no smooth bounded-energy solution on `R^3 x [0,inf)` with the same force and datum
(navier-stokes.txt:36-46, p.1); the periodic case is Corollary 10.6 (:6471-6478, p.125-126). The
Euler paper claims unforced finite-time breakdown from `u_0 in C_{c,sigma}^inf(R^3)`, with
`limsup ||grad u||_{L^inf} = inf` and `int_0^{T_*} ||curl u||_{L^inf} dt = inf`, and no claim
about `||u||_{L^inf}` (euler.txt:33-37, p.1). The Lean certificate certifies the two Navier-Stokes
nonexistence sentences (that is, (C) and (D)) and two Euler statements, one of which does carry the
constructive content. Against our fences: the Navier-Stokes claim passes Lemma A by construction,
since the velocity is unbounded with the explicit rate `sup_x|u| ~ tau^{-1/2-h}`, `A = 1/2 + h`,
`0 < h < 1/100` (:312, p.7; :779, p.15; :6413, p.124), so the fence that kills bounded-velocity
cascades does not bite. Read through the fence-checklist lens, the singular point is the spatial
origin, hence on the axis, and the singular set is a single space-time point, which is what Lemma B
and the `02` class would demand; but Lemma B does not apply, because the full flow is not
axisymmetric (only the leading field is; the pulses carry nonzero integer angular frequencies,
:3880-3881, p.75). The rate is Type II by the smallest margin the construction carries,
`sqrt(1-t) sup|u| >= c tau^{-h}` with `h < 1/100`, so no Type I exclusion bites, and swirl is the
growing component: on-axis, Type II, swirl-carrying, exactly the corridor `02` predicted for the
axisymmetric sub-class, reached by a non-axisymmetric construction that is therefore not bound by
that note's hypotheses. We found no fatal or major defect in either paper. Every major-grade
candidate raised against the mathematics was refuted or downgraded on adjudication. The one
surviving MAJOR is about the Lean certificate's coverage, not about any mathematical claim.

---

## 2. Fence checklist (`01_two_fences.md` section 4), item by item

| # | Item | Verdict | Basis |
|---|------|---------|-------|
| 1 | Is `sup_{t<T} ||u||_{L^inf}` finite? (decisive; Lemma A) | **PASS** (unbounded) | Theorem 1.1 asserts `limsup ||u(t)||_{L^inf} = inf` (:40-41, p.1); rate `tau^{-A}`, `A = 1/2+h` (:312, p.7), realised along `x_tau = (sqrt(2 X_in tau),0,0)` by (3.6) (:779, p.15) and (10.21) (:6413, p.124) |
| 2 | Does `||u||_{L^3}` diverge? (expected, unproved for `f != 0`) | **PASS**, by our computation; paper never does it | Core volume `tau^{3/2-h}` and speed `tau^{-1/2-h}` (:815-816, p.16; :355-364, p.8) give `||u||_{L^3} ~ tau^{-4h/3}`, exponent below 0.0134. Exterior swirl `~ r^{-1-2h}` from (4.29) (:1697-1699, p.33) gives the same exponent |
| 3 | Where is the singular point; is the flow axisymmetric? | **PASS, non-binding** | Full flow is not axisymmetric (:3880-3881, p.75), so Lemma B does not apply; the singular point is nonetheless the origin, on the axis, and the singular set is one space-time point, inside CKN |
| 4 | Rate (Type I exclusion) | **PASS** as Type II | `sqrt(1-t) sup_x|u| >= c tau^{-h} -> inf` from (3.6)/(10.21); margin over Type I is exactly `tau^{-h}`, `h < 1/100` |
| 5 | Energy bookkeeping (decisive) | **PASS**, stronger than requested | Lemma 10.4 proves `||u(t)||_2^2 + 2 int_0^t ||grad u||_2^2 <= F(t)^2` directly from the equation (:6255-6257, p.121); `E_core ~ tau^{1/2-3h} -> 0`, `int D_core < inf`, both under `h < 1/6` (:815-822, p.16) |
| 6a | Force in (5) / (9), quantifiers honoured jointly | **PASS** | Compact support in space and time makes (5) trivial with a single fixed smooth `f` (:6193-6195, :6242-6244, p.120); no layer-by-layer constant survives |
| 6b | Data in the Clay class | **PASS** | `u_0 = 0` satisfies (4) vacuously; (C) quantifies existentially over the data (:38, p.1) |
| 6c | Nonexistence in class, with the uniqueness step supplied | **PASS**; the step our note said was always missing is present | Lemma 10.5 (:6277-6279, p.121): smooth on `R^3 x [0,T]`, `v in L^inf_t L^2_x`, no enstrophy, no decay of the competitor, no energy inequality, no Leray-Hopf property. Pressure recovered by Riesz transforms with `P`'s spatial growth unrestricted (:6280-6281) |
| 6d | `nu > 0` fixed; `f` a single function of `(x,t)` | **PASS** | Rescaling `u_nu(x,t) = sqrt(nu) u(x/sqrt(nu), t)` (10.22) (:6431-6436, p.124), permitted by the quantifier order in (C)/(D) and by our own item 6(d) |
| — | Lemma A applicability | applies, not violated | `alpha = 2`, `nu > 0`, force in (5) |
| — | Lemma B applicability | does not apply | flow not axisymmetric |
| — | `02` axisymmetric class | does not apply, but profile matches | on-axis, Type II, swirl-carrying |

No item is UNDETERMINED at the level the checklist asks about. Every item that our own instrument
called decisive (1, 3, 5, 6) is passed with text we could locate and check.

---

## 3. Surviving findings

| id | lens | sev | Corrected claim | Evidence |
|----|------|-----|-----------------|----------|
| SF-05 | lean-faithfulness | **MAJOR** | The two certified Navier-Stokes statements assert only Clay (C) and (D). Three distinctive clauses of the printed Theorem 1.1 appear in no certified statement: the compact `K` with `supp u(.,t) u supp p(.,t) subset K`, `sup_{0<=t<1}||u(t)||_{L^2} < inf`, and `limsup_{t->1}||u(t)||_{L^inf} = inf`. They live only in project-local vocabulary (`SpeedUnboundedAtOne`, `velocity_support`, `pressure_support`), which no independent reference file states, so the project's rendering of them is refereed by nobody outside the project. Coverage defect, not unsoundness: those fields are load-bearing inside the proof of the refereed (C), so trivialising them would break the derivation rather than silently pass it. The Euler side does expose the analogous clauses to the checker, which shows the asymmetry is a choice | `ComparatorChallenges/NavierStokes.lean:272-284` (only entries in `NavierStokes.json` `theorem_names`, only NS entries in `formalization.yaml` `main_results`); `NavierStokes/ProblemStatement.lean:94-96`; `NavierStokes/R3CompactCandidate.lean:24-37`; contrast `ComparatorChallenges/Euler.lean:170-183`; navier-stokes.txt:36-41 vs :43-44 |
| SF-02 | lean-faithfulness | MINOR | `formalization.yaml` alignment rows name source statements strictly stronger than the Lean declarations they map to: "Theorem 1.1 (Navier-Stokes on R^3)" and "Corollary 10.6" are constructive, `navier_stokes_breakdown_R3` / `_periodic` are the (C)/(D) consequences. The same file's `main_results` and `README.md` describe them correctly, so this is one metadata row overstating, not a misrepresentation of content. The constructive NS content is proved but only in project-local definitions, so the independently-referenced surface certifies only (C) and (D) | `formalization.yaml:109-116` vs `:49-62`; `ComparatorChallenges/NavierStokes.lean:272-284`; `NavierStokes/R3ActualCandidate.lean:18-21` |
| P-05 | provenance | MINOR | Same defect seen from the repository side, with the sharpest instance: the certified candidate bundle carries no clause for `sup_{0<=t<1}||u(t)||_{L^2} < inf`; the module that does state the printed theorem calls itself "a target proposition, not an asserted theorem" and is never discharged; the lemma proving the energy clause is sorry-free and has zero call sites | `NavierStokes/R3CompactCandidate.lean:24-37`; `NavierStokes/R3/ProblemStatement.lean:146-150`, `:217`; `NavierStokes/R3/CompactEnergy.lean:341-343` |
| P-01 | provenance | MINOR | The repository and the two manuscripts are bound to each other by no stable identifier in either direction: the yaml names the papers by title with no URL, DOI, version, date or hash; the PDF sha256 values appear nowhere in the tree; neither paper mentions Lean, a formalization or a repository. The tree pins everything else by commit, which makes the omission conspicuous. Consequence is bounded: the challenge files are self-contained over Mathlib and anchored to Fefferman's numbered conditions, so the Lean side is checkable without the PDFs; what is lost is the paper-to-Lean audit and any detection of manuscript drift | `formalization.yaml:17-29`, `:32`; repo-wide grep for hash/doi/arxiv returns only three claymath.org links at `README.md:16-18`; single root commit, message "." |
| P-02 | provenance | MINOR | The shipped Comparator instructions omit the `systemd-run --property=RestrictAddressFamilies=~AF_UNIX` wrapper that upstream states its guarantee for, and which upstream says guards a landrun vulnerability. landrun itself is retained, so the sandbox is present; the hardening layer is not. Platform-specific, so the fix is a caveat plus the Linux command | `ComparatorChallenges/README.md:3-9` against the upstream Comparator README |
| P-03 | provenance | MINOR | The two READMEs, read in the order presented, lead a checker to `lake build` the solution libraries unsandboxed before running the check, which is what upstream assumption 2 forbids; the default targets glob in both named solution modules, and an existing `.lake` means the solution is not rebuilt under the sandbox. The checking recipe taken alone is clean; the missing item is the ordering caveat or a fresh-checkout instruction | `README.md:29-42`; `lakefile.toml:3`, `:17-24`; `ComparatorChallenges/NavierStokes.json`, `Euler.json` |
| FR-11 / F-02 | clay-statement | NOTE | Alternative (D) is proved for competitors whose pressure as well as velocity is 1-periodic. The paper states this openly and grounds it in the Clay erratum; the reference formalization encodes the same hypothesis, so the class proved is the accepted post-errata (D). The hypothesis is genuinely load-bearing: with `P = c(t).x + periodic` the displayed Gronwall does not close. Anyone quoting (D) from the uncorrected text should know the proof depends on the erratum. On `R^3` the same gauge freedom is eliminated with no extra hypothesis by Lemma 10.5 | navier-stokes.txt:6466-6467 (p.125), :6521-6526 (p.126); `ComparatorChallenges/NavierStokes.lean:266-269`; `01_two_fences.md:464-465` |
| SF-04 | lean-faithfulness | NOTE | The Lean force class is Fefferman (5) plus smoothness on `R^3 x [0,inf)`; compact spatial support, compact time support and vanishing near `t = 0` do not reach the certified statement. Direction is safe (the force is existentially quantified, so the certified statement is the weaker one) and the width is required for the statement to be Clay (C) as posed | `ComparatorChallenges/NavierStokes.lean:175-177`, `:185-191` vs navier-stokes.txt:36 |
| EUL-05 | euler-claim | NOTE | Euler Theorem 1.1 localizes nothing: `0 < T_*(u_0) <= T_infinity <= S_base`, and the paper states the lifespan may end strictly earlier. The diverging gradients sit at the origin only along the approximating sequence; each `U_j` is an exact smooth solution. A downstream summary saying "blowup at the origin at time `T_infinity`" would overstate on both counts | euler.txt:2252-2254 (p.44), :2189-2190, :2246-2250 (p.43-44) |
| FR-16 | fence-checklist | NOTE (scope) | Audit boundary: Theorem 4.6 with Appendix C, Lemma 7.4 with Proposition 7.5, Proposition 9.6 with Lemma 9.8, and Lemma 5.4 were not verified. The paper also never checks itself against Type I exclusion or Liouville results (zero occurrences of either phrase); CKN and the `L^3` endpoint appear only in the p.1 history paragraph, where it is noted that CKN allows isolated singularities and that the `L^3` endpoint is for the unforced problem | navier-stokes.txt:69-73 (p.1), :702-703 (p.14), :5849-5851 (p.113), :3017 (p.57) |

---

## 4. Findings raised and then refuted or downgraded

* **F-03** (Lemma 10.3's support `K x [0,2]` contradicts `C_c^inf(R^3 x (0,inf))`): **refuted outright**; "support contained in" is a superset containment, jointly satisfiable, and the proof discharges both ends (:6236-6238, p.120).
* **F-02** (alternative (D) proved in a strictly narrower class): downgraded to NOTE; `p` is a solution unknown, so the erratum can only bind solutions, and no global smooth competitor with nonzero linear pressure is exhibited.
* **F-04** ("Consequently" is false as a connective): downgraded; two of the three allegedly hidden Lemma 10.5 hypotheses are in Theorem 1.1 itself and the third is discharged by the scaling at :6454-6459.
* **F-13** (Corollary 10.6's uniqueness step under-specified): downgraded; (D) carries no energy condition, so "two smooth periodic solutions" names the whole competitor class, and the length asymmetry tracks the weaker hypotheses of Lemma 10.5.
* **F-14** (`H^3` maximal-lifespan asserted, not proved): stands as a NOTE only; true, standard, and load-free, with one quotation-hygiene caution at :6422.
* **FR-12** (no upper bound on `||u||_{L^inf}`): downgraded; true but the named consequence was false, since compact support makes the KNSS side hypothesis `|u| <= C_0/r` trivial, and KNSS is inapplicable anyway (non-axisymmetric, forced).
* **FR-13** (`h` not quantified): downgraded; `h` is pinned from above by an explicit chain (`h < min{1/100, lambda, e^{-T_d}}`, :1761, p.34), and the `r^{-1}` criteria are pass/fail hypotheses that any fixed `h > 0` fails outright.
* **S-01** (scope wording, "Millennium problem statement"): downgraded; the phrase is the title of reference [8] and the paper makes no prize claim; the finding also misattributed `u_0 = 0` to (C).
* **S-03** (thin bibliography): downgraded; the modern blowup-construction literature is cited in the companion Euler paper where it belongs; the weak-strong ask misread Lemma 10.5, and KNSS/CSTY exclude Type I, not Type II.
* **SF-01, SF-03, SF-06** (Theorem 1.1 not formalized at all; zero datum absent; `L^2` clause absent everywhere): all downgraded; the construction is proved in-tree, the certified target is Clay (C) which has no from-rest clause, and the energy clause is proved in a parallel structure. The residue is coverage, which is what SF-05 and P-05 now carry.
* **P-04** (Euler challenge self-authored): downgraded; disclosed in the challenge file's own header, and Formal Conjectures has no Euler statement to copy. Residue is one clarifying clause in `ComparatorChallenges/README.md`.
* **P-06** (upstream pin inconsistent): downgraded; miscounted files and links, and the two challenge/solution definition blocks are byte-identical, so the adapted text is not ambiguous. Residual: pin the Euler challenge file.
* **EUL-01** (the sign `m . Mv > 0` is a single point of failure): downgraded; the sign is proved with an explicit margin in section 4.5 (:1707-1718, p.33-34), uniformity over labels is inside (4.5), and the two regimes are exhaustive.
* **EUL-02** (global-in-space uniformity of the `O(k^{-1/4})` errors): downgraded; all shift bounds are global `H^6(R^3)` norms, the mean pressure Hessian is obtained pointwise without an inverse Laplacian, and every `k`-carrying coefficient is compactly supported.
* **EUL-03** (stage-independence of exponents is a convention): downgraded; Propositions 3.1 and 4.1 are single-step statements with no stage index anywhere in section 4, and (3.16) resets the child constant.
* **EUL-04** (packets are exact only through a vanishing-viscosity limit): downgraded; velocity converges strongly in every `H^j_loc`, the pressure is re-derived from the strong limit, and the viscous term is killed by a deliberate two-derivative margin.
* **EUL-06** (`delta` scale reversed): downgraded to a reading-pipeline note; the `.txt` inverts `delta_j^{-1}`, and (5.16) plus the seed line fix the intended reading.

---

## 5. Weakest joints, by the mechanism lens

**Navier-Stokes.** (1) *Theorem 4.6 with Appendix C, pp.32-33 and 157-165.* One profile must be
regular and `eta`-analytic at the axis, have exactly zero stress outside `[X_a, X_b]`, be the exact
radial heat solution for `X >= X_b`, match five cumulative radial moments, and keep a uniform strict
cone margin including `v_s > 2` on the closed annulus. The move that buys the cone is Proposition
C.2: an `N log X` shear modulation of amplitude `O(N^{-1})` that changes radial derivatives at order
one, followed by a localized correction restoring all five moments. Check that the restoring
correction preserves the cone margin, the edge flatness and the analytic collar. (2) *Proposition
9.6 with Lemma 9.8, pp.107-113.* The induction closes on margins as thin as `0.18 - 2 kappa_s >
0.17`; completeness of the residual inventory is decisive, and the stage-independence of the
derivative loss `K_m` in (9.18) is what licenses Lemma 5.4 (note that `C_{j,m}` and `P_{j,m}` do
depend on the stage, and Lemma 5.4's shrinking cutoffs are the second leg). (3) *Proposition B.2,
pp.144-146.* The profile equations are Cauchy-Kovalevskaya in `eta`, so analyticity is forced, not
convenient; every later smooth edit must sit strictly outside the analytic collar.

**Euler.** (1) *Proposition 3.1 via Lemma 3.3, pp.10-24*, where exactness of each packet rests on a
vanishing-viscosity limit; the load-bearing step is the uniform-in-`(m, nu)` a priori inequality
(3.54) at :1062-1069 and its bootstrap closure, not the limit passage. (2) *Proposition 4.1, Claims
1 and 3, pp.28-35*, the sign `m . Mv > 0` for `tau >= 1` at every label, together with the
exponential separation (4.12) covering the complementary regime; proved, but the whole
pressure-Hessian budget hangs on it. (3) *Sections 5.5 and 5.7, pp.38-43*, the exponent race in
which `log(l_j^{-1}) = x_{j-1}/j^{7/2}` must beat `c C_* log K_h = c C_* x_{j-1}/(j-1)^4`; the
margin grows like `j^{1/2}` and the base case is secured only by choosing `J` after `c` and `C_*`,
with no numeral attached to `c`.

---

## 6. Lean faithfulness: what is certified versus what is claimed

**Certified.** Four declarations. On the Navier-Stokes side, `navier_stokes_breakdown_R3` and
`navier_stokes_breakdown_periodic` (`ComparatorChallenges/NavierStokes.lean:272-284`), which are
Fefferman (C) and (D) verbatim, copied from an independent third-party statement file at a pinned
commit; the definition blocks are byte-identical between challenge and solution sides. On the Euler
side, `euler_breakdown_R3` and `exists_compact_smooth_euler_singularity`
(`ComparatorChallenges/Euler.lean:170-184`), the second of which does carry the constructive
content: compactly supported nonzero smooth divergence-free datum, `0 < Tstar <= 1`, a maximality
biconditional, `limsup` of the `C^1` norm equal to top, and a divergent vorticity integral.

**Not certified.** The construction half of Navier-Stokes Theorem 1.1: compact support of `u` and
`p` in a fixed `K`, the uniform `L^2` bound, and the `L^inf` blowup. Those exist in the tree, in
project-local vocabulary, and are proved; they are simply not what the independent reference file
states or what Comparator compares against.

**Loophole ledger, all closed at source level.** Viscous sign pinned to the analyst's Laplacian by
the repo's own bridge lemma; the Bochner-integral junk value in the energy clause closed by an
adjacent `MemLp` field, and in any case a hypothesis on the excluded competitor; `toL2`'s junk
fallback closed by an adjacent integrability field; `derivWithin` on `Set.Ici 0` is the honest
one-sided notion and `ContDiffOn` on `univ x Ici 0` is non-vacuous; the Euler witness's
extension-by-zero cannot manufacture an infinite limsup, because the filter localizes strictly
inside the lifespan; the solution class is inhabited, so no theorem is discharged by vacuity.
Feature census over 2486 files: zero `axiom`, `native_decide`, `unsafe`, `implemented_by`, `extern`,
`opaque`, `set_option`, `decide +kernel`, `ofReduceBool`; the only `sorry`s are the four intended
challenge placeholders. Permitted axioms are the standard three.

**Weaker renderings on the Euler side.** The `C^1` blowup is stated as the sum of the velocity and
gradient suprema rather than the gradient supremum alone, and maximality of the lifespan is asserted
within a class that also demands all-order Sobolev regularity of the strong time derivative, with no
uniqueness clause. **One fragility.** `Euler/Solution.lean:40-41` installs a local order instance so
the restatement elaborates like the challenge; that is the single hand-steered point of statement
identity. **Caveat.** Nothing here was built; the zero-sorry and axiom claims in
`formalization.yaml` are read off the file and the `#print axioms` lines, not off a compiler run.
`formalization.yaml` self-reports review status "self-assessed".

---

## 7. Provenance checklist, and what a replicator must run

| Item | State |
|------|-------|
| Paper-to-repository binding | **absent** in both directions (P-01) |
| Repository history | single unsigned root commit, message ".", no tags, no releases, issues disabled, Apache-2.0 |
| Toolchain pin | `leanprover/lean4:v4.34.0-rc2`, published; Mathlib pinned at the matching tag `85e3a25e`, its release job green, so `lake exe cache get` is honest for Mathlib only |
| Upstream reference pin | NS side pinned consistently at `8bf45ed7`; Euler solution-side copy pins `8323e878`; three unpinned branch links |
| Checker versions | landrun, lean4export, nanoda_bin: **unrecorded** |
| CI / committed check output | none |
| Sandbox instructions | landrun retained; upstream's `systemd-run` hardening omitted (P-02) |
| Build-then-check ordering | unwarned; default targets glob in both solution modules (P-03) |
| Dead weight | 75 modules and about 26k lines built but unreachable from either result, including a 41-module `NavierStokes.R3*` cluster; lock file still names the old package `fluidEquations` |
| Scale / cost | 2486 files, 616,276 lines on top of Mathlib, no cache of its own; no machine, memory or timing guidance published |

**A replicator must:** (a) start from a **fresh checkout** and do **not** run `lake build` first;
(b) `lake exe cache get`; (c) install pinned versions of landrun, lean4export and nanoda_bin,
recording them, since the repository does not; (d) run each Comparator config inside the upstream
`systemd-run --property=RestrictAddressFamilies=~AF_UNIX ... -- bash -c 'lake env <comparator>
<config.json>'` wrapper, on Linux; (e) read the two challenge files as the statements actually
certified, and treat `formalization.yaml`'s alignment rows as source pointers, not as claims about
what was proved; (f) separately audit `NavierStokes/R3CompactCandidate.lean` and
`NavierStokes/ProblemStatement.lean` by hand, since nothing outside the project refereed them; and
(g) hash both PDFs and record the hashes alongside the commit, since nobody else has.

---

## 8. What would change our assessment

1. **A defect in sections 4 through 9 of the Navier-Stokes paper.** Nothing in the (C)/(D)
   interface, the checklist, or the mechanism ledger caught this construction out. If it falls, it
   falls in the profile construction or the correction induction. Concretely: if (3.4) degrades from
   flatness to *every* order `N` to *some* `N`, Lemma 10.2 fails and there is no smooth force at
   all; if Theorem 3.1(ii) fails on the cutoff transition annulus, the localization produces a
   residual that does not extend through `t = 1`; if any interaction omitted from Proposition 9.3's
   inventory gains less than `1/10` in the exponent, the induction does not close; if `K_m` in
   (9.18) turns out to depend on the stage, Lemma 5.4's summation is unlicensed.
2. **A reading of the Clay erratum on which pressure periodicity is not a condition on admissible
   solutions.** Then (D) is proved in a narrower class than (D), and the paper's periodic claim
   would need either a competitor with nonzero linear pressure ruled out, or a referee's ruling. We
   judge this reading untenable, but it is a referee's call, not ours.
3. **A build.** If the tree fails to compile, or if `#print axioms` reports anything beyond the
   standard three, or if Comparator rejects either statement match, the certificate collapses and
   only the human papers remain.
4. **Either constructive clause being exposed to the checker.** If the NS constructive statement
   were stated in the independent reference file and certified, SF-05 and P-05 close and the
   Navier-Stokes certificate reaches parity with the Euler one.
5. **Independent confirmation of the Euler literature claims.** The novelty framing rests on two
   2026 preprints (Chen; Shkoller) we could not check offline. If those already cover smooth
   compactly supported data on boundaryless `R^3`, the Euler paper's positioning changes, though not
   its correctness.
6. **A paper revision.** With no hash binding, a later manuscript could diverge from what was
   audited here without anything in the tree registering the change.

One reading worth stating flatly for downstream use: the Euler paper is a **different mechanism**
and could not have been upgraded into the Navier-Stokes one. Its conclusion is unbounded `C^1` norm
with no claim of unbounded velocity (euler.txt:33-37, p.1), which places it on the bounded-velocity
side of Lemma A. The Navier-Stokes route is a concentrating self-similar swirling vortex with
genuinely unbounded velocity. Anyone describing the pair as "Euler blowup upgraded to
Navier-Stokes" has it backwards. That the two papers sit on opposite sides of Lemma A, exactly where
Lemma A says they must, is the strongest external corroboration our fences received from this
reading.

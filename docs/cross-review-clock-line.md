# Scientific scope review of the historical forced-route archive

This is a bounded read-only review of the existing `forced-route-2026-09` repository at the publication handoff identified as commit `a37d92e`. It reviews claim scope and selected concrete consistency issues; it is not a fresh proof audit of every cited paper, an independent Lean build, or a replication of every script.

All file and line references below refer to the historical repository layout. If the unified archive relocates those files, retain a path map. Historical documents should remain identifiable as historical; a corrected current summary should take precedence over their unrevised headlines.

## Recommended current summary

The archive contains analytic necessary-condition arguments under explicit solution and forcing hypotheses; elementary Lean formalizations and their stored axiom audit; a conditional logarithmic amplification-clock assembly with unresolved analytic and numerical-certification issues; statement-level audits of external forced-fluid formalizations; and replicated finite algebra, interval and exploratory numerical checks from the strained-core work. It does not establish Navier–Stokes blowup, an iterated return mechanism, a complete formalization of the PDE arguments, or a refutation of an external regularity theorem.

## 1. The log-clock assembly is not a proved sharp clock

The latest historical assembly already says **“(S3) is not proved”** at `theorems/03_log_clock/THEOREM_S3/THEOREM_S3.md:69`. The label “all proved” attached to its numerical budget at lines 34, 173, 220 and 469 must therefore not become a current theorem claim. Row G1 at line 265 gives the safer description: computed, conditional on unresolved items.

The concrete unresolved points in its own list include:

- **Unmet exact-plateau premise:** lines 305–314 say the V.4 hypothesis requires an exact local plateau, while the selected tanh radial profile never equals that plateau. A perturbative replacement is proposed but not proved there. As written, assuming that premise for the selected datum does not yield an applicable example.
- **Global viscous strain error:** lines 316–320 distinguish the existing pointwise estimate near a tracked point from the required shell-wide strain-functional error, and state the latter is not done.
- **Maximum-principle passage:** lines 292–303 leave the regularization/approximate-maximum argument unwritten. The collar and gradient estimates depend on it.
- **Regularity and tails:** lines 328–352 leave the angular Lipschitz-to-classical passage and noncompact tanh-tail corrections outstanding.

Consequently the abstract's suggestion that the approximate-maximum step “carries all the others” (`THEOREM_S3.md:540–543`) and the suggestion that it is the entire distance to a BFG refutation (lines 523–527) overstate what has been discharged. The actual list contains distinct hypotheses and an explicit incompatible datum premise.

The object called a doubling clock is the first time the vorticity reaches **3/2 of the nominal amplitude**, defined at line 108. Even a completed finite-amplification result would not be a singularity or a repeating cascade. The datum is axisymmetric **without swirl** (lines 93–98), and the proposed result concerns one finite amplification time.

### Additional issues beyond “exactly eight items”

These follow directly from the displayed definitions and code, so the unified summary should not repeat “nothing else is unproved” (`THEOREM_S3.md:290`). They are recorded as unresolved consistency/certification issues, not repaired here.

1. **Time-window coverage.** The theorem assumes estimates only on \([0,\tau]\), \(\tau=c_*/(ML)\), at lines 151–163, but concludes a time \(t_*=\tau(1+\varepsilon)>\tau\) at line 168. The cap is computed at \(c=c_*\) (line 142), and the budget still uses that window (lines 407 and 457; `t3_budget.py:24`). The notes do not show the needed estimates on the larger claimed interval. A repair must recompute or majorize the time-dependent budget and cap on the actual conclusion interval.

2. **An exact derivative mismatch.** For the shifted log-radius tanh profile at `THEOREM_S3.md:91`, writing \(u=\log(\rho/\rho_0)\), one obtains
   \[
   \rho_0\Theta'(\rho_0)=\frac{\operatorname{sech}^2(1)-\operatorname{sech}^2(1-L/\varepsilon_r)}{2\varepsilon_r}.
   \]
   This is not \(1/(2\varepsilon_r)=2\), as claimed at line 103. At \(\varepsilon_r=1/4\), its large-L limit is \(2\operatorname{sech}^2(1)\approx0.83995\). Thus the asserted exact slope match to the separately defined radial ramp is not available. This does not by itself falsify every proposed upper bound, but it invalidates that stated transfer justification.

3. **Range mismatch in the imported Hessian constant.** `t3_budget.py:44` explicitly imports `161.7735` as a bound for \(\lambda\in[1,3/2]\). `THEOREM_S3.md:242` repeats that range, while the clock's cap is about 1.842 at line 142. A bound on the larger range, or a smaller justified cap, is needed. The underlying H-K2 datum also uses a different radial profile (`01_HK2_second_derivative_bound.md:115–121` versus the assembly's log-radius formula), so its numerical constant is not automatically transferable unchanged.

4. **Discrete scans are not certified extrema.** `t1_datum.py:98–106` computes the purported infimum \(r_h\) by sampling 61 points, and lines 115–132 check monotonicity only through step-0.01 increments. Lines 156–183 compute profile suprema using finite differences and a finite grid. Yet A1/A2 are labelled PROVED at `THEOREM_S3.md:234–235`. A sampled minimum is generally an upper bound on a true infimum, and a sampled maximum is generally a lower bound on a true supremum: precisely the unsafe directions when used as favorable rigorous constants. Analytic derivative arguments or outward-certified global enclosures are required. The earlier referee already flags the monotonicity issue at `u1_REFEREE/NOTE.md:179–186`; it was not removed merely by another scan.

5. **A failed fixed-point iteration is not a nonexistence proof.** `t2_gamma_CR.py:72–76` says it returns failure when iteration runs away and parenthetically identifies this with no fixed point. The summary calls the resulting threshold an existence threshold (`THEOREM_S3.md:354–360`). The earlier `u2_REFEREE/NOTE.md:166–179` explicitly distinguishes iteration convergence from actual existence. A numerical failure threshold alone does not certify the absence of a supersolution/fixed point.

The “all measured” and “all proved” decimal columns should therefore be described as **exploratory evaluations using different constant choices**, conditional on the analytic scheme. They are not validated quantitative PDE bounds in this review. The title “sharp” should remain a historical title, not an editorial endorsement of optimality.

## 2. Necessary conditions: preserve the exact hypotheses

`theorems/01_two_fences.md:64–95` gives the bounded-velocity continuation argument for dissipation order \(1<\alpha\le2\), in an explicit strong solution class, with force integrability through an interval beyond the proposed singular time. This is an analytic argument using a cited product estimate; it is not a Lean PDE theorem. In particular, a smooth bounded force on the whole spatial space does not automatically have the required L2/Sobolev integrability, as the note itself explains.

`theorems/02_forced_axisymmetric_on_axis_typeII.md:72–83` explicitly carries the strong/energy hypothesis H*. Its theorem at lines 100–137 makes the distinctions that the unified summary must retain:

- The on-axis conclusion is conditional on that solution class; nonemptiness additionally uses uniform spatial decay (lines 109–112).
- The exclusion of the joint space/time CSTY-form bound is stated without the separate Dec hypothesis (lines 116–121).
- Exclusion of a velocity Type I bound depending only on \(T-t\) requires **Dec** as well (lines 123–129).
- The swirl-free conclusion assumes both swirl-free initial data and swirl-free forcing (lines 134–137).

The compressed reading at lines 139–144 obscures the Dec qualification when it says the global velocity must grow faster than the Type I rate. Quote the precise theorem parts instead. The paired referee report is campaign-internal, and its initial verdict concerns an earlier statement (`02_REFEREE_REPORT.md:8–16`); the later H* repairs should not be mistaken either for an unchanged refuted statement or for external journal peer review.

A small textual proof repair is still visible in `02_forced_axisymmetric_on_axis_typeII.md:175–178`: a function multiplied by the sharp ball indicator is called Schwartz. That is false in general. The localization argument should use a smooth cutoff with its derivative seminorms controlled. This is a local repair issue, not a new conclusion about the Navier–Stokes problem.

## 3. The dissipation ceiling is an architecture/model diagnostic

`analysis/01_dissipation_ceiling_of_the_layered_cascade.md:3–8` explicitly says it is not a theorem. Its frozen diagonal damping model and “nothing else changes” idealization are at lines 12–26. The reoptimized one-half ceiling is an amplitude-ODE supremum (lines 97–118), and the infinite-regularity structural incompatibility is marked conjectural outside the fixed released architecture (lines 156–165).

Safe summary: the analysis identifies severe damping and correction-budget obstructions to a direct transfer of the specified layered scheme. Unsafe summary: it proves that no positive-dissipation blowup, no alternative smooth forcing construction, or no full NS mechanism is possible. Formalizing the exponent arithmetic does not strengthen the model's PDE scope.

## 4. The archive's own Lean results are elementary cores

The detailed `lean/README.md` is unusually explicit and should control the unified description. Its opening says that none of the five files formalizes a PDE (lines 3–6), and lines 233–348 enumerate what is outside Lean.

I counted **91** theorem declarations in `lean/Cores/*.lean`, **91** `#print axioms` commands, and **91** entries in the stored `lean/AXIOMS.txt`. The root `README.md:12` still says 33 and is stale. The stored axiom report names only `propext`, `Classical.choice`, and `Quot.sound`; my source string scan found no `sorry`, `admit`, or `native_decide` in `Cores/`. This review did **not** rerun Lean, so a current publication build claim should be tied to a new build log or clearly attributed to the archived report. “No additional axioms/sorry in the audited core statements” is clearer than implying axiom-free mathematics.

Concrete limits:

- `CriticalityExponent.lean:119–123` rescales an arbitrary function. It does not prove PDE covariance or a continuation theorem. The real supremum identity has an explicit unbounded-set default-value caveat at lines 138–180; use its bounded version when describing a norm statement.
- `CircleMeasure.lean` supplies the Euclidean circle/rotation measure argument. CKN, the suitable-weak-solution theory and the parabolic-to-Euclidean measure step are not formalized (`lean/README.md:258–278`).
- `StrainedCore.lean:350–370` evaluates the specified radial-integral combination. Its comment expressly assumes the identification with pressure and the angular reduction. The chain rule at lines 467–474 assumes the derivative relation and gives no feedback sign or global solution.
- `Receiver.lean:326–337` contains arithmetic identities, including a ring identity named `weighted_projection`. This is not a formal proof of the underlying Hilbert-space projection inequality, the actual receiver's existence or a gain theorem. Those exclusions are made explicit at `lean/README.md:327–342`.

Safe current wording: **five Lean developments with 91 elementary algebraic/geometric/integral statements and a stored standard-axiom report, supporting selected steps of the analytic notes; no Navier–Stokes theorem formalized.**

## 5. External Lean audit is separate from compilation and statement matching

`audit/01_lean_statement_audit.md:3` states that its audit did not typecheck the external sources. Lines 13–16 and 53–65 distinguish statement-fidelity issues, provenance, build-target coverage and pending compilation/comparator checks.

Its reported external Euler statement is a weakening of the paper's geometric theorem: the trusted challenge omits axisymmetry, torus geometry and several circulation conclusions. The Boussinesq challenge uses an existential datum rather than the paper's specified datum; the separate affinecore project has a different statement. Three project directories should not be advertised as one-to-one Lean proofs of the three announced fluid models.

The `replication/` directory in the reviewed checkout contains only its README. That README labels the Euler build **in progress** at lines 1 and 12, says the comparator was not attempted at line 10, and lists the planned outputs at line 13. This archived evidence does not support “independently rebuilt and replayed all released proofs.” Challenge `sorry` placeholders are also different from proof dependence on `sorryAx`; the relevant Solution theorem and comparator result must be checked, as the audit explains.

The archive's own 91-core build, the external Euler build, and an external comparator/nanoda replay are three distinct verification claims. Do not merge them into one green badge.

## 6. The 94-run strained-core replication is reproducibility evidence

`strained-core/REPLICATION.md:5–8` accurately scopes the rerun as scripts-as-written, not a theorem audit. Lines 22–24 report 94 runs passing and 16 scripts not rerun; the detailed exclusions are at lines 165–184. Some entries are aggregate runners, resolution variants or exploratory evaluations rather than distinct theorems.

The important safeguards to retain are at lines 194–203: one regenerated result still says its requested target was not met, and a coarse pressure calculation deliberately reproduces the previously identified underresolved result. “All executed programs exited zero” is a reproducibility claim, not a declaration that every output is accurate or every research target is achieved. Its actual interval endpoint reproductions are stronger, separately identified evidence for those finite bounds.

The pass-nine planar slice is explicitly a reduced model at line 161. It must remain distinct from the later full three-dimensional finite-cylinder pilot and from a certified whole-space NS trajectory. The imported strained-core README also describes a passes-2–9 snapshot (`strained-core/README.md:3–8`), so it should not silently be relabelled as a replication of subsequently added passes 10–11.

## Publication handling

Preserve the original notes, scripts, critiques and stored verification records with their provenance. Add this scope review or an equivalent current errata notice near the unified README's status statement. Describe internal “refereed,” “proved,” “PASS,” and “all measured” labels in their specific contexts. Mark new mathematical repairs as new revisions with their own checks rather than editing the historical record into apparent prior completeness.

This review found meaningful, reusable mathematics and substantial reproducibility work. Its scope findings require a conservative current synopsis; they do not diminish exact results by association with unclosed parts, and they do not turn an unclosed proposal into a proof.

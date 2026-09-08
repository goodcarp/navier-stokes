# A3 collision scan — 2026-09-08 — Alpöge–Buckmaster release (Euler / Boussinesq / IPM with smooth forcing)

Source: cims.nyu.edu/~tristanb/{statement,euler,ipm,boussinesq}.pdf, read at source 2026-09-08 (PDFs + SHA256SUMS filed
under `sources/pdfs/alpoge-buckmaster-2026-09-08/`). Lean repo: github.com/tristanbuckmaster/fluid_lean, commit
d0124689 (2026-09-08 00:07 -0400), shallow clone read in the session scratchpad only.

## What was released
- Euler (112 pp): Thm 1.1 — forced axisymmetric Euler WITH swirl on R³, smooth compactly supported data in a torus
  away from the axis, force f ∈ C^∞(R³×[0,T*]) supported in the same torus; solution smooth on [0,T*), unique in the
  locally-Lipschitz finite-energy class; circulation and meridional velocity stay BOUNDED; ‖∇Γ‖∞, ‖ω‖∞ → ∞ and
  ∫‖ω‖∞ dt = ∞. Mechanism = Córdoba–Martínez-Zoroa layered cascade (bounded amplitude, growing derivatives).
- Boussinesq (76 pp): same mechanism on R², forces C^∞ through T*, temperature bounded, ‖∇θ‖∞ → ∞, limsup ‖ω‖∞ = ∞.
- IPM (57 pp, with Coiculescu): T², space-time smooth force, gradient blowup, ρ(t) converges in C^η, η<1.
- Statement (4 pp): claims hypo-dissipative NS blowup too (Lean unfinished, not released); reports OpenAI told him an
  internal model has forced full-NS blowup on R³ and T³ with smooth force = Clay options (C)/(D), ~100 pp, unseen,
  prompted "in the past few days"; no public artifact as of this scan (web search 2026-09-08: nothing posted).

## Lean certificates (read, not built)
- Three projects: affinecore, boussinesq-blowup, euler-blowup. Challenge.lean states each theorem over Mathlib only;
  Solution.lean proves it; comparator.json restricts axioms to propext/Classical.choice/Quot.sound. sorry count = 0
  outside the three deliberate Challenge.lean placeholders; axiom count = 0. Lean 4.32.2, Mathlib pinned 81a5d257.
- Euler Challenge statement certifies the CORE only: smooth solution, compactly supported data, smooth force with all
  mixed derivatives bounded and support in a fixed ball, ‖ω‖∞ → ∞ (as a minorant g), ∫‖ω‖∞ = ∞, Lipschitz-class
  uniqueness. NOT in the Lean statement: axisymmetry, swirl, torus support, bounded circulation/velocity, ‖∇Γ‖ → ∞,
  any rate. (formalization.yaml says so itself.)
- Audit items: (1) IPM paper's SHA-256 placeholders "[BOUSSINESQ SHA-256], [EULER SHA-256]" are UNFILLED in the
  posted PDF — the paper is not hash-bound to the repo; (2) affinecore README: Challenge.lean "has not yet been read
  against the theorem stated above by the responsible maintainer"; (3) build needs Mathlib from source for 4.32.2,
  ~100–150 GB RAM, hours on 24–36 cores — no laptop replication; no independent replication exists yet.

## Fences (walls already in FL-004, applied to this object)
Dissipation written |∇|^α as in CMZZ (ARMA 2026): CMZZ blowup for α < (22−8√7)/9 ≈ 0.093 with rough force.
- W6 (LPS/Serrin): bounded velocity is critical at α = 1 and subcritical for α > 1 ⇒ a bounded-velocity cascade
  (which this mechanism is) cannot blow up for α > 1. Navier–Stokes is α = 2.
- W3 (CKN): at α = 2 an axisymmetric singular set is rotation-invariant, so an off-axis singular point gives a
  circle of singular points, contradicting P¹(S) = 0 ⇒ axisymmetric NS singularities lie ON the axis. The released
  Euler blowup lives in a torus away from the axis.
⇒ Any forced full-NS (α = 2) blowup must have unbounded velocity (ESS: ‖u‖_{L³} → ∞) and, if axisymmetric, sit on
  the axis and be Type II (KNSS/CSTY). That is the (K)/clock-line territory, not the CMZ program's.

## Ledger consequences (operator/registrar decisions owed)
- P-F1 wording "certified smooth no-boundary Euler blowup inside the window": literal reading FIRES (Lean-certified,
  smooth data, R³, window to 2027-08-01); intended reading (P-Model's phrase, unforced) does NOT. Seal-design
  ambiguity on "forced" — adjudicate; recommend applying the pre-committed ACTION (consume within one session)
  regardless of Brier scoring.
- Q4/R5-B for THIS certified Euler structure: viscous continuation exists only for small α (CMZZ, AB's claimed
  hypo-dissipative result) and is fenced at α ≤ 1 by W6, at α = 2 off-axis by W3.
- FL-000 unchanged.

## Addendum 2026-09-08 (WIN sitting, refereed external results; files under campaign/external/alpoge-buckmaster-2026-09-08/)
- FENCES_forced_NS_blowup.md — Lemma A PROVED (bounded velocity ⇒ continuation for |∇|^α, 1<α≤2; CMZZ's u ∈ L^∞_{t,x} verbatim), Lemma B assembled from CKN (on-axis at α=2); 3-lens refutation, 20 corrections applied; checklist items 2 and 4 are unforced-only (ESS, KNSS).
- SCOUT_dissipation_ceiling_of_the_layered_cascade.md — mechanism ceiling α < β/(2Q) (amplitude ODE, sup 1/2); released schedules Q_q = Q_*+q ⇒ no fixed α>0 survives all stages under diagonal damping; C^∞ force needs k_q→∞ while the papers' budgets cap k_q by Q_q ⇒ C^∞ force and fixed α incompatible in the released correction scheme (scope: assumptions 1–2); CMZZ α₀ located in their §4 (s=0 root). 3-lens refutation, 16 corrections.
- LEAN_STATEMENT_AUDIT.md — 33 agents: no loopholes; Euler and Boussinesq certificates are undocumented weakenings (bridge discards proved axisymmetry/swirl/torus/∇Γ clauses; Boussinesq datum existential); no hash binding; affinecore default build omits Challenge/Solution; comparator runs not committed.
Prediction on record: AB's unreleased hypo-dissipative result cannot be the released architecture plus damping; expect a different scheme or a finite-regularity force.

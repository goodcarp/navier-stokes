# Finite-time blow-up for the three-dimensional incompressible Euler equations with smooth force — Lean 4 formalisation

This repository contains a Lean 4 / Mathlib proof of the theorem stated below. The proof is complete: no file other
than the trusted statement file [`Challenge.lean`](Challenge.lean) contains a `sorry` (that file states the theorem with a
placeholder proof, as the comparator / Palomar template requires), and the proof rests on no axiom beyond Lean's standard
three (`propext`, `Classical.choice`, `Quot.sound`). `Challenge.lean` is the only file a reader must trust;
[comparator](https://github.com/leanprover/comparator) checks the proof against it.

> Research artifact. Not maintained; issues are welcome, pull requests are not accepted.

## The theorem

**Theorem.** There exist T > 0, maps u, f : ℝ × ℝ³ → ℝ³ and p : ℝ × ℝ³ → ℝ with the following properties,
every time t below ranging over the half-open interval [0, T). u and p are C^∞ on [0, T) × ℝ³, u(t, ·) is
divergence-free and square-integrable, and ∂ₜu + (u·∇)u + ∇p = f holds at every point (the time derivative is
taken within [0, ∞), hence one-sided at t = 0); u(0, ·) has compact support; every component of f is C^∞ on
[0, T) × ℝ³ with every mixed space–time derivative bounded there, and f vanishes outside a fixed ball; there is a
function g with g(t) → +∞ as t ↑ T such that at every t the vorticity ω = curl u satisfies ‖ω(t, x)‖ ≥ g(t) at some
point x — in particular ‖ω(t, ·)‖_∞ → ∞; for every M there exist T′ ∈ [0, T) and an integrable function h ≥ 0 on
[0, T′] such that at every s ∈ [0, T′] one has h(s) ≤ ‖ω(s, x)‖ at some point x, and ∫₀^{T′} h ≥ M — in particular
∫₀^T ‖ω(t, ·)‖_∞ dt = +∞; and for every T′ with 0 < T′ < T, u belongs on [0, T′] to the finite-energy Lipschitz
class (velocity fields locally Lipschitz on [0, T′] × ℝ³ with sup_t (‖u(t)‖_{L²} + ‖∇u(t)‖_{L^∞}) < ∞, divergence-free
almost everywhere, solving the momentum equation with force f in weak, pressure-free form) and every member of that
class with the same initial velocity coincides with u on [0, T′].

In words: there exist a time T > 0, smooth compactly supported finite-energy initial data, and a smooth force vanishing
outside a fixed ball with every derivative bounded up to time T, which together produce a classical solution of the
incompressible Euler equations on [0, T) × ℝ³ that is unique in the natural strong class on every shorter interval and
whose vorticity becomes unbounded as t ↑ T.

## Building and checking

Toolchain: `leanprover/lean4:v4.32.2`. Mathlib is pinned in `lake-manifest.json` (commit
`81a5d257c8e410db227a6665ed08f64fea08e997`) together with its own dependencies. Do not run `lake update`.

```
lake build                               # Mathlib publishes no prebuilt objects for Lean 4.32.2 (= 4.32.0 plus two kernel
                                         # soundness fixes), so Mathlib is compiled from source, then this project
lake env lean scripts/PrintAxioms.lean   # prints the axioms of the main theorem
```

Expected: no errors; exactly one `declaration uses 'sorry'` warning, from `Challenge.lean` (the placeholder statement
file), and none from `EulerBlowup/`, `vendor/` or `Solution.lean`; every line of the axiom printout reads
`[propext, Classical.choice, Quot.sound]`. The build is large (about 1,100 modules); plan for a machine with on the
order of 100 GB of memory and a few hours at two dozen parallel jobs rather than a laptop.

`Challenge.lean` states the theorem using Mathlib only; `Solution.lean` proves that statement from the library's main
theorem; `comparator.json` is the configuration for [comparator](https://github.com/leanprover/comparator), which
type-checks the statement file independently, checks that the solution inhabits exactly that statement, restricts
axioms to the standard three, and replays the proof.

## Layout

- `EulerBlowup/` — the construction and its analysis, and under `EulerBlowup/Num/` interval-arithmetic certificates,
  each an inequality between explicit rationals decided by the kernel; the fixed-point arithmetic they use is itself
  proved sound in Lean.
- `vendor/cm24-r2/` — a self-contained supporting library (kernel estimates, transport, uniqueness) built by this
  project on top of Mathlib and included as a build arrangement.
- `Challenge.lean`, `Solution.lean`, `comparator.json` — statement, solution bridge, comparator configuration
  (Palomar template layout); `scripts/PrintAxioms.lean` prints the axioms.

Authorship: all Lean code in this directory, the trusted statement file `Challenge.lean` included, was written by
Claude (Anthropic) under the direction of Levent Alpöge, who wrote no Lean by hand and read `Challenge.lean` against the
theorem stated above. `formalization.yaml` (`automation`, `review`) records the same account.

## Licence

Apache License 2.0; see `LICENSE` and `NOTICE`.

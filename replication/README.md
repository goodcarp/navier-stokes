# Independent replication of the euler-blowup Lean certificate (in progress)

Target: github.com/tristanbuckmaster/fluid_lean, commit d0124689230b58b4f86e7b90ac59de06404b3b6b (2026-09-08 00:07 -0400), project `euler-blowup`
(theorem `euler_smooth_force_blowup` in Challenge.lean; comparator.json permits propext, Quot.sound, Classical.choice; enable_nanoda true).

Host: MacBook Pro (2019), Intel i7-9750H, 6 physical cores, 17 GB RAM, macOS 24.6.0. Toolchain: leanprover/lean4:v4.32.2 (installed by elan 2026-09-08).
Mathlib pinned by the project's lake-manifest.json (rev 81a5d257c8e410db227a6665ed08f64fea08e997, tag v4.32.0); no prebuilt objects exist for this toolchain, so Mathlib was compiled from source.
Command: `LEAN_NUM_THREADS=4 nice -n 19 lake build` (script run_build.sh, memory watchdog at 2 GB free), started 2026-09-08 06:34:47Z.
Then: `lake env lean scripts/PrintAxioms.lean` (expected: every line `[propext, Classical.choice, Quot.sound]`; exactly one `declaration uses 'sorry'` warning, from Challenge.lean).
Comparator replay (leanprover/comparator with nanoda) not yet attempted on this host.

Status log (UTC): [4085/4141] modules at 09:24Z; build alive.
Outputs to be added here when the build ends: build.log (tail), PrintAxioms output, sha256 of the checked-out tree, wall time.
What this replication checks: that the shipped Lean sources compile against the pinned Mathlib on an independent host and that the main theorem depends on the standard axioms only.
What it does not check: statement faithfulness (see ../audit/), the comparator's statement-match and nanoda kernel replay (needs the comparator tool), and anything about the paper's prose.

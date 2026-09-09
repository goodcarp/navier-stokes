# Independent replication of the euler-blowup Lean certificate (first attempt lost at 9,739 of 9,769 modules; moving to a rented host)

Target: github.com/tristanbuckmaster/fluid_lean, commit d0124689230b58b4f86e7b90ac59de06404b3b6b (2026-09-08 00:07 -0400), project `euler-blowup`
(theorem `euler_smooth_force_blowup` in Challenge.lean; comparator.json permits propext, Quot.sound, Classical.choice; enable_nanoda true).

Host: MacBook Pro (2019), Intel i7-9750H, 6 physical cores, 17 GB RAM, macOS 24.6.0. Toolchain: leanprover/lean4:v4.32.2 (installed by elan 2026-09-08).
Mathlib pinned by the project's lake-manifest.json (rev 81a5d257c8e410db227a6665ed08f64fea08e997, tag v4.32.0); no prebuilt objects exist for this toolchain, so Mathlib was compiled from source.
Command: `LEAN_NUM_THREADS=4 nice -n 19 lake build` (script run_build.sh, memory watchdog at 2 GB free), started 2026-09-08 06:34:47Z.
Then: `lake env lean scripts/PrintAxioms.lean` (expected: every line `[propext, Classical.choice, Quot.sound]`; exactly one `declaration uses 'sorry'` warning, from Challenge.lean).
Comparator replay (leanprover/comparator with nanoda) not yet attempted on this host.

Status log (UTC): [4085/4141] modules at 09:24Z (Mathlib phase). Project phase: the build was killed three times by the memory watchdog while other work ran on the same laptop, restarted from cache each time at two threads and then one, and reached [9739/9769] at 2026-09-09 02:12Z, inside the ring modules (`EulerBlowup/Ring/S6r/`), which are the slow kernel-checked interval certificates. At 11:52Z the host rebooted; the build tree lived in a temporary directory and was lost. No axiom printout was obtained. Wall time to that point: about 19.5 hours of compute.
What was learned: the sources compile module by module against the pinned Mathlib on an independent host through 9,739 of 9,769 modules with no error; this is not a replication, because the last thirty modules and the axiom printout were never reached.
Next: the same recipe on a rented 32-core, 128 GB CPU host (the authors' README asks for 100 to 150 GB), using `replicate-on-cloud.sh` in this folder, which builds all three fluid_lean projects and the other 2026-09-08 release, prints axioms, and clones the comparator for the kernel replay. Outputs to be added here: build.log (tail), PrintAxioms output, comparator verdicts, sha256 of the checked-out trees, wall time.
What this replication checks: that the shipped Lean sources compile against the pinned Mathlib on an independent host and that the main theorem depends on the standard axioms only.
What it does not check: statement faithfulness (see ../audit/), the comparator's statement-match and nanoda kernel replay (needs the comparator tool), and anything about the paper's prose.

# (working name: forced-route-2026-09) — DRAFT, not published

Theorems, analyses, and certificate reviews produced on 2026-09-08 around the Alpöge–Buckmaster release
(finite-time blowup for 3D Euler, Boussinesq and IPM with smooth forcing; Lean certificates at
github.com/tristanbuckmaster/fluid_lean, commit d0124689) and the forced route to the Clay problem.

Layout
- theorems/    proved statements with proofs (01: the two fences; 02: forced axisymmetric Clay (C) class is on-axis, Type II, needs swirl (refereed); 03: the log doubling clock, assembled modulo eight named items (THEOREM_S3), with (H-K2), U1, U2 and their referee reports)
- analysis/    the dissipation ceiling of the layered cascade (refereed); the event record
- audit/       statement-faithfulness audit of the three Lean certificates (33 seats, two refuters per finding); source hashes
- replication/ independent build of euler-blowup: toolchain, build log, PrintAxioms output, comparator output [in progress]
- lean/        formalizations of the elementary cores (33 theorems, standard axioms, Lean 4.33.1 / Mathlib pinned)
- strained-core/ the Codex half (passes 2-9 of the strained-core pressure-feedback line, verbatim, hashed) with a second-host replication of all check scripts (94/94)
- page/        the publication page (HTML, OpenHell identity; copy is a draft for the author)

Grades are campaign-internal. No claim of a Clay solution is made anywhere in this repository.
Every number in these notes came out of a script checked in beside it; every citation was read at source.

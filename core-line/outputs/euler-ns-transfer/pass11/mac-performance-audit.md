# Full MAC performance audit

These are timing diagnostics of the finite-cylinder numerical operator. They
provide no approximation-error or Navier–Stokes existence certificate. No
`evolve_full_mac.py` edits were made by this audit. Only one small full step was
run; the tests at radial size 384 constructed selected sparse matrices, without
initializing or evolving a large velocity field. Concurrent machine load was
not controlled, so relative structural findings are more useful than projected
production runtimes.

## 1. Warm-up and steady step are different costs

`profile_mac_step.py` profiles the requested case `nr=96,nz=256,J=4`, padding,
`R=8,L=16,mapping=4,nu=.001,dt=2e-5`. It explicitly builds the 640 distinct
Stokes factors before profiling exactly one step. The timed step created no new
factor and returned finite values.

| Operation | Wall time |
|---|---:|
| Constructor | 0.152 s |
| Initial-field evaluation/projection | 2.632 s |
| Cold 640-factor preparation | 24.315 s |
| One cached full step | 5.736 s |
| Four nonlinear evaluations, included in step | 4.061 s |
| Two diffusion half-steps, included in step | 1.590 s |

The cached nonlinear work accounts for about 71% of the profiled step. FFT
kernels alone consume 1.993 s of self-time; the 2550 individual SuperLU solves
consume 0.902 s. `diffuse` also spends 0.447 s in its Python body, with additional
allocation overhead. The cold factorization kernel consumes only 1.348 s of the
24.315 s preparation; repeated sparse block assembly, format conversion, and
index checks dominate that preparation. Nested profiler times are not additive.

Evidence: `mac-step-timing96.json`. The profiler adds overhead; subsequent
unprofiled nonlinear timings below were appreciably smaller.

## 2. The alternative ordering improves solves, without a fill collapse

`profile_stokes_permutation.py` reconstructs the same saddle matrix and compares
the existing cell ordering with `(theta_i,z_i,right_face_i,pressure_i)`.
Both omit the same final pressure unknown for the zero-zero gauge. All 40 tested
matrix pairs at radial size 96 have relative solution differences below
`4e-16`; relative linear residuals are below `3e-15`.

The current factors are already sparse. For dimension 383, the largest tested
`nnz(L)+nnz(U)` is 3426 with the current ordering, and 3519 with the alternative.
The maximum fill ratios are respectively 1.866 and 1.833, over matrices of
different base sparsity. Thus these data reject a large fill explosion as the
current bottleneck; they do not reject benefits from a better triangular-solve
layout. For example, at `m=4,kindex=1`, current/alternative factor storage is
3404/3475 entries at size 96 and 13843/14139 entries at size 384.

`profile_stokes_solve_layout.py` alternates the measured ordering over eight
repetitions and compares two separate right-hand sides with a two-column solve.
For selected radial-size-384 matrices:

| Axial index | Current, separate | Alternative, separate | Current, batched | Alternative, batched |
|---|---:|---:|---:|---:|
| 1 | 1.272 ms | 0.494 ms | 0.745 ms | 0.356 ms |
| 32 | 1.299 ms | 0.442 ms | 0.817 ms | 0.343 ms |

This is a repeatable solve-layout benefit on these sampled matrices, not a
guaranteed full-step speedup. The alternative factor has slightly more entries.
The same tests at size 96 show the same qualitative ordering. Matrix-only tests
used a tiny axial constructor cache; the physical tested wavenumbers were
assembled explicitly, so they are the same matrices as in the full grid.

Evidence: `stokes-permutation-timing96.json` and
`stokes-solve-layout-timing.json`. The sampled factor comparison does not by
itself test the negative-wavenumber sign transformation in a new batched
diffusion implementation; that implementation should retain the existing
vertical-component similarity transformation and run the operator checks.

## 3. Four FFT workers give a modest direct improvement

`profile_fft_workers.py` evaluates the unchanged nonlinear operator on the same
initialized small field, alternating warmed one-worker and four-worker runs.

| Workers | Three wall-time samples | Median |
|---|---|---:|
| 1 | 0.5411, 0.5057, 0.4747 s | 0.5057 s |
| 4 | 0.3606, 0.3809, 0.3819 s | 0.3809 s |

The measured speedup is 1.328 and the resulting arrays are exactly equal in
this test. The source hash remained unchanged during this comparison. Evidence:
`fft-worker-timing96.json` records the imported nonlinear method as well as its
source hash, allowing comparison if the solver is subsequently optimized.

## 4. Concrete next optimizations, preserving the chosen operator

1. Cache the sparse Stokes block structure and update its wavenumber-dependent
   values, reducing repeated cold assembly. Keep the same mass, constraints,
   gauge, and matrix entries.
2. Use the alternative permutation and pair the positive/negative axial
   right-hand sides for the shared factor. Reuse scratch buffers. Recheck
   diffusion, weighted divergence, and dissipation after implementing it.
3. Use four FFT workers where the machine is not simultaneously oversubscribed.
   The measured gain is 25% in elapsed nonlinear time on the small case.
4. Fuse repeated spectral transforms and reuse the gathered velocity in the
   curl evaluation. The present nonlinear function gathers each field twice.
   A real angular transform is another possible reduction, provided the same
   retained modes, padding, Fourier normalization, and conjugate symmetry are
   preserved and compared against the existing operator.

No reduction of resolution, retained modes, or padding is justified by this
performance audit. The full evolution and its validation remain separate work.

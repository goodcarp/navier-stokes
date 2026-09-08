# Independent fast MAC audit

Date: 2026-09-08. Result: **PASS for the projected finite-cylinder operators on their retained Fourier state space.** No continuum evolution, spatial convergence, useful duration, or whole-space endpoint has been certified by this check.

`check_fast_full_mac.py` compares `FastMAC` against `MAC` and independently assembles signed-frequency dense Stokes saddle systems. Results are in `fast-full-mac-checks.json`. The run took 11.63 seconds. Solver/checker hashes were unchanged before and after the run:

| File | SHA-256 |
|---|---|
| `evolve_full_mac.py` | `2e671a3599fc05ed3d0fdef657a75bb0e9dceeab05fb1b96c4c5fc7243a7fe95` |
| `fast_full_mac.py` | `3a7f76fcdaa98ccc39a3215380236f1154b09677bda7afca96c728ce1a5eca75` |
| `check_fast_full_mac.py` | `f91ea27b659aa40176dc05d5cef14dce84552dbe36f059d32b9f0a476184a895` |
| imported `check_full_mac.py` | `4acc133afb512c8dd3c477fd072b5be80ada3f1fde2eafcd4659222b5510cbb4` |

The six parameter sets have `(nr,nz,J,dealias)` equal to `(12,31,0,padding)`, `(12,32,1,padding)`, `(16,33,2,padding)`, `(28,64,3,padding)`, `(12,31,2,filter)`, and `(16,32,2,filter)`. They include odd/even axial grids, m=0 and m=4,8,12, independent complex positive/negative axial coefficients, and a radial block extending across the implementation's 24-row chunk boundary. Inputs are projected random fields, scaled for the small full-step test. They have real m=0 coefficients in axial physical space and zero excluded Fourier modes.

## Measured comparisons

| Check | Largest relative discrepancy |
|---|---:|
| Projected nonlinear output, fast versus reference | 7.59e-16 |
| Raw nonlinear output after applying the same axial restriction | 7.26e-16 |
| Constrained CN output, fast versus reference | 3.45e-16 |
| CN versus independently assembled **signed-k** dense saddle solve | 1.23e-14 |
| Independently assembled matrix versus modified sparse template | 6.22e-17 |
| Complete CN/RK4/CN step, fast versus reference | 5.64e-16 |
| CN at h=0 versus retained projected input | 3.38e-16 |
| Weighted CN midpoint energy identity | 7.45e-17 |
| Projected nonlinear kinetic work | 1.38e-17 |
| Arbitrary angular phase equivariance | 7.08e-16 |
| Fractional axial shift equivariance | 8.90e-16 |
| Cartesian odd/C4 symmetry after a full step | 7.16e-16 |

The largest normalized divergence of the projected nonlinearity is 1.20e-14; after CN or a complete step it is below 1.94e-15. All tested CN updates and small complete steps dissipate kinetic energy. This establishes those discrete checks, rather than a stability theorem for arbitrary step sizes.

## Fourier fusion

For positive angular index j, the real-field partner of coefficient `(j,k)` is `(-j,-k)` with conjugated value. The fast assignments to `negative_indices=(-ki)%nwork` and `-j` implement precisely this relation. Reversing only the angular index would be incorrect. The m=0 input already has axial Hermitian symmetry. With forward-normalized Fourier coefficients, the fused inverse/forward pair has the same product normalization as the reference's axial resampling followed by angular transforms.

Padding uses a product grid large enough for retained quadratic interactions; it keeps the input high axial frequencies that the older two-thirds state restriction discarded. Odd grids include their largest positive and negative integer frequencies, and even grids omit the Nyquist input. The tests compare both signs of the highest retained frequency and include arbitrary fractional axial shifts, which are sensitive to incorrect two-dimensional reality pairing.

**Raw-output API distinction:** `FastMAC.nonlinear(U,project=False)` already sets excluded axial output modes to zero; `MAC.nonlinear(U,project=False)` returns them until a later projection. Their unrestricted raw outputs therefore differ (up to 0.703 in these test states). After the same axial restriction they agree to roundoff. This distinction does not change the projected time-stepping operator or its work against a retained state. It would matter if a caller treated the fast raw array as the complete unfiltered product or used its norm to estimate unresolved residuals. Likewise, this audit does not claim arbitrary unfiltered even-N inputs are supported equivalently: their Nyquist treatment differs from SciPy resampling. The current initial projection and all time-stepping stages satisfy the retained-input requirement.

## Sparse template and signed axial solves

Writing `M=diag(mass)`, the constrained CN saddle matrix is

`[[M + alpha*(K_j+k^2 M), -D(k)^* V],[-V D(k), 0]]`.

Only the velocity diagonal has quadratic k dependence; only pressure/axial-velocity cross entries have linear k dependence. The template is assembled with unit k, then its identified diagonal entries receive `alpha*k^2*mass` and both axial pressure couplings are multiplied by k. The checker verifies these sparse index sets are unique and disjoint and compares the entire permuted matrix against an independent dense construction. It also tests alpha=0 and the m=k=0 pressure gauge omission. The altered right-face ordering is a permutation of the same unknowns.

The matrix K_j does not mix axial velocity with radial or azimuthal velocity. With `S_z` reversing axial velocity, `D(-k)=D(k)S_z`, while the velocity block commutes with `S_z`. Thus the negative-k saddle problem is obtained from the positive-k factor by reversing the axial RHS and solution. The paired solve does that on its second column. For positive FFT storage index k, its partner is `nz-k`, which is valid on both tested grid parities. Dense solves with the actual signed k independently check this reuse, including the highest retained pairs.

The inherited pressure weights, viscosity Gram form, radial wall and axis stencils are unchanged. Their earlier structural review remains relevant. In particular, this speedup does not turn the finite no-slip/periodic diagnostic into a full-space field, enforce all higher-order axis regularity, or bound the reconstruction and residual errors required by the endpoint bridge.

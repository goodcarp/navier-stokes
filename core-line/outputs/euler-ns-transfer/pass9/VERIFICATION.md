# Verification scope

All eight programs passed. The packaged output is saved in
[CHECK_RESULTS.txt](CHECK_RESULTS.txt).

The runner executes eight finite programs:

1. General-envelope covariance, pressure-flux and mean-jet algebra.
2. A newly recomputed full-pressure initial-jet interval enclosure.
3. Exact rational replay of its bounds and dependency paths.
4. H4/H7 residual majorants, core formulas and endpoint-error algebra.
5. Exact affine-envelope finite-time gain algebra and its new interval moment.
6. Fixed-amplitude off-neutral core and initial-gain rational inequalities.
7. Analytic-solution, incompressibility and circulation checks for the slice solver.
8. Sensitivity comparisons of the stored nonlinear slice runs.

Run `python3 checks/run_checks.py`. The interval engine and established
cutoff, pressure and root certificates are imported from unchanged sibling
pass4–pass8 packages. Dyadic endpoints are compared as exact rational values.
Analytic proof steps, whole-space pressure identities, Sobolev embeddings
and local existence remain part of the trust base. No Lean or other formal
proof-kernel replay is performed.

The actual-NS certificate concerns initial derivatives only. Intervals for
pressure-error upper-bound expressions do not bound the actual error from
below. Fixed-height radial extrema are distinct from a global spatial
maximum. The new off-neutral amplitude uses its correct core identity.

The affine finite-time result is an exact linear reduced-model theorem.
The nonlinear slice computations are finite-radius planar diagnostics with
prescribed strain. Their production runs use 768/1536 radial cells, paired
time refinement at A=950 and A=1020, plus angular and outer-boundary
comparisons at A=950. The solver checks include known Gaussian evolution and a known Poisson
solution with second-order convergence, plus discrete divergence and
circulation tests. These checks do not enclose continuum, boundary, tail,
nonlinear evolution or compact three-dimensional errors.

Stored slice energy integrals exclude exterior harmonic tails. Small
residual circulation is an initial numerical representation error, with a
nonintegrable mean energy tail if interpreted literally on the whole plane.
Raw finite-grid data are therefore not already an admissible whole-space
approximation for the H4 theorem. A smooth divergence-free reconstruction
and a full residual bound are still required.

The H4/H7 theorems are conditional validation tools. Their residuals,
majorants and endpoint-class margins have not been enclosed for the selected
full 3D solution. No full compact 3D time integration, useful actual stage
duration, inherited return or blowup proof is claimed.

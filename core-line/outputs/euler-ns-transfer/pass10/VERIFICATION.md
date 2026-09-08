# Verification scope

Run `python3 checks/run_checks.py`. The runner copies the code and reference
data into a temporary directory, then executes the listed checks without
modifying this archive. It reads the unchanged pass9 initial-jet certificate
for the linear-ansatz obstruction.

The checks cover exact time-polynomial identities and the H4 obstruction;
cylindrical/helical algebra; the complete initial field and independent
Cartesian differences; initial coefficient and quotient identities;
manufactured scalar pressure convergence; co-rotating reconstruction and
norm factors; the exact omitted-image tail; and stored full-field results
against independent pressure and energy integrals. A ninth program verifies
the rational initial norm lower bounds, the sharper derivative hierarchy's
coefficients, its scaling and endpoint bookkeeping, and rigid-rotation
cancellation through tensor order four.

The full production grids are recorded, not rerun by the short verification
runner. Reproduction commands, for example, are:

```
python3 experiments/full_initial_evolution.py --nr 768 --nz 2048 --L 16 --pressure fd
python3 experiments/full_initial_evolution.py --nr 768 --nz 4096 --L 32 --pressure fd
python3 experiments/initial_whole_space_integrals.py --orders 32 48 72
python3 experiments/initial_periodic_images.py --help
```

Production JSON retains the diagnostics as originally recorded; subsequent
code adds direct-jet and rotating-frame diagnostics. The legacy finite-volume
pilot also retains its discovered source-repair inconsistency and is not
an accepted NS coefficient computation. Numerical pressure convergence is
separate from exact discrete incompressibility, axis regularity, and a
smooth whole-space reconstruction.

Gaussian quadrature comparisons and Fourier/finite-difference grids are
not rigorous continuum enclosures. The image-tail certificate bounds only
the omitted images; it does not enclose the quadrature of the included
images. The time-polynomial obstruction concerns the specified validator,
not actual-flow failure. Rotation formulas are exact, but their practical
coefficient/residual norms have not been enclosed.

The new derivative hierarchy is a conditional whole-space energy estimate.
Its initial norm lower bound concerns a uniform/supremum coefficient; a
time-zero lower bound alone does not lower-bound a time integral. Neither
the hierarchy's time-dependent coefficients and residuals nor its comparison
ODE have been enclosed for a proposed evolution.

Analytic PDE identities, classical smooth local theory, Sobolev embeddings,
and the named earlier cutoff/initial-jet certificates remain in the trust
base. No Lean or other formal proof-kernel replay was performed. No actual
useful full 3D duration, inherited return or blowup is certified.

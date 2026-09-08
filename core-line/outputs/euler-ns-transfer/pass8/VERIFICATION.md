# Verification scope

Eleven finite programs cover the actual covariance/torque algebra, the
harmonic-slab full-pressure bound, a newly recomputed old-family initial-jet
certificate, its independent audit, the actual energy identity, the separate
Kelvin model and exact endcap estimate, the exact exterior-pressure identities, and the
two newly recomputed leading-envelope certificates with an independent
same-family audit.

All eleven programs passed; the saved run is in [CHECK_RESULTS.txt](CHECK_RESULTS.txt).
Run `python3 checks/run_checks.py`. It regenerates:

- The 2048-panel radial-root/Green-moment/fourth-jet certificate.
- The 1024-panel leading-envelope torque/extraction/viscosity certificate.
- The 1024-panel leading-envelope full-pressure norm certificate.

The interval engine and earlier residual/cutoff certificates are read from
the unchanged sibling pass5, pass6 and pass7 packages. Exact dyadic endpoints
are compared using rational arithmetic. Proofs also depend on the displayed
analytic integrations, the Hankel/slab energy argument, the established
local existence theory, and earlier profile/retuning certificates.
These are not Lean or other formal proof-kernel replays.

Intervals for upper-bound expressions are not two-sided enclosures of the
underlying pressure or acceleration. Fixed-height radial maxima are not the
unrestricted global maximum. The old trailing-family jets are not asserted
for the new leading family.

Model-only statements include the Kelvin covariance's integrated lifetime
budget and frozen mean-deformation estimates. The experiments directory
contains exploratory whole-space component quadratures; agreement between
grids is not an interval certificate. No actual Navier–Stokes time integration
was run. The new common positive-time interval follows qualitatively from
strict margins and smooth dependence, without an enclosed useful duration.

There is no inherited return, repeated compatible stage sequence, singular
solution, or full Navier–Stokes blowup proof.

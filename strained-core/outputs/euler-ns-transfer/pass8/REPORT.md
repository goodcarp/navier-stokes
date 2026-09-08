# Eighth pass: diagnose depletion and construct an initially energy-gaining replacement

8 September 2026. **The new leading-envelope candidate simultaneously grows
the receiving mean maximum and its perturbation energy, while retaining
favorable full central pressure feedback, including viscosity.**
The full Navier–Stokes goal remains open: no useful quantitative stage
duration, inherited return or singular solution has been proved.

## The selected candidate

The [same-datum pressure theorem](leading-full-pressure-feedback.md) and
[leading-envelope construction](leading-envelope-alternative.md) use a
compact radial envelope centered at 7/50 with width 1/10, angular mode four,
radial carrier -20, seed amplitude between 1/8 and 1/4, and exact outer-swirl
retuning. Viscosity remains 1/1000.

| Actual initial quantity for this new family | Proved lower or upper bound |
|---|---|
| Receiving mean angular-momentum derivative | G_t >8871/800 (>11.088) |
| Fractional nonaxisymmetric kinetic-energy derivative | K'/K >12853/2205 (>5.829) |
| Complete central pressure-feedback test | p_zz'(0)+32 <-12493177/1120000 (<-11.15) |
| Central ratio acceleration | (b/Omega)''(0)>12493177/2240000 (>5.57) |

The declining spatial envelope permits positive receiving torque even though
the covariance is positive. Its reversed chirality extracts energy from the
decreasing mean swirl. Total kinetic energy still dissipates.
The full nonlocal pressure and axial viscous costs are retained.

The [independent full-family audit](leading-family-independent-audit.md) and
[independent envelope/viscosity audit](leading-envelope-independent-audit.md)
accepted these inequalities. Compactness and smooth local dependence give
some common positive interval of simultaneous gain. Its useful numerical
duration has not been enclosed.

## What changed the direction

For the previous trailing seed, an [exact actual energy identity](actual-fluctuation-energy.md)
shows K'(0)/K(0)<-43.66. Its initial transfer drains the supplying fluctuation.
The pressure-corrected [Kelvin model](kelvin-covariance.md) also shows that
favorable winding depletes covariance instead of amplifying it; that
all-time model statement is kept separate from actual NS dynamics.

The [full-pressure initial jet certificate](initial-receiver-deceleration.md)
proves that the old packet-center receiver strongly decelerates. Even after
following its radial maximum at fixed axial height, its value has second
derivative below -240000 and the radius initially moves outward. The
unrestricted global maximum can select other axial locations, and no
turning time follows from an initial second derivative alone.

Those old-data statements do not apply automatically to the new leading
envelope. The replacement is a new initial datum, not a regenerated endpoint
of the old solution. The steep envelope is essential; the old constant-
amplitude Kelvin calculation is not its exact evolution.

## The next mathematical task

Prove a useful interval for the actual coupled new solution, carrying the
spatial covariance gradient, pressure-generated vertical velocity, changing
mean geometry, viscosity, core curvature and the old field. Then prove a
controlled inherited endpoint and a compatible continuation. Rebuilding a
fresh seed at the endpoint would not satisfy this requirement.

The [actual covariance identities](covariance-evolution.md), [whole-space
pressure-gradient bound](covariance-evolution-pressure-slab.md), and
[endcap pressure formulas](exterior-dynamics-pressure.md) provide tools
for this next estimate. Exploratory component quadratures have not been
used as certified PDE time evolution.

Run `python3 checks/run_checks.py` to execute all eleven programs, including
recomputation of three new interval certificates. [Verification](VERIFICATION.md)
states the exact scope and trust base.

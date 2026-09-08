#!/usr/bin/env python3
from pathlib import Path
import shutil,json,hashlib,sys
import mpmath,sympy
root=Path(__file__).resolve().parents[2];source=root/'work/pass8'
top=root/'outputs/euler-ns-transfer';dest=top/'pass8'
(dest/'checks').mkdir(parents=True,exist_ok=True)
(dest/'experiments').mkdir(exist_ok=True)
for p in source.glob('*.md'):shutil.copy2(p,dest/p.name)
checks=['verify_covariance_evolution.py','verify_covariance_pressure_slab.py','certify_fixed_circle_acceleration.py','verify_fixed_circle_acceleration_audit.py','verify_actual_fluctuation_energy.py','kelvin-covariance-check.py','verify_endcap_gradient_bound.py','exterior-dynamics-verify.py','certify_leading_envelope.py','certify_leading_pressure.py','verify_leading_family_independent_audit.py','run_checks.py']
certificates=['fixed-circle-acceleration-certificate.json','leading-envelope-certificate.json','leading-pressure-certificate.json']
for name in checks+certificates:shutil.copy2(source/name,dest/'checks'/name)
for name in ['exterior-dynamics-pressure.py','exterior-dynamics-pressure-48.json','exterior-dynamics-pressure-96.json']:
    shutil.copy2(source/name,dest/'experiments'/name)
(dest/'REPORT.md').write_text("""# Eighth pass: diagnose depletion and construct an initially energy-gaining replacement

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

Run python3 checks/run_checks.py to recompute the three new interval
certificates and the eleven finite programs. [Verification](VERIFICATION.md)
states the exact scope and trust base.
""")
(dest/'VERIFICATION.md').write_text("""# Verification scope

Eleven finite programs cover the actual covariance/torque algebra, the
harmonic-slab full-pressure bound, a newly recomputed old-family initial-jet
certificate, its independent audit, the actual energy identity, the separate
Kelvin and endcap models, the exact exterior-pressure identities, and the
two newly recomputed leading-envelope certificates with an independent
same-family audit.

Run python3 checks/run_checks.py. It regenerates:

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
""")
sources=dict(date='2026-09-08',inherited_manifests=[dict(path=f'../pass{i}/SOURCES.json',sha256=hashlib.sha256((top/f'pass{i}/SOURCES.json').read_bytes()).hexdigest()) for i in [5,6,7]],
    primary_sources=[dict(title='NIST DLMF 10.9.4, Poisson integral for Bessel J',url='https://dlmf.nist.gov/10.9.E4',checked='Used for the self-contained positive-weight derivative bound; independently verified by the endcap-bound agent.'),dict(title='NIST DLMF 10.22(v), Hankel transform and inversion',url='https://dlmf.nist.gov/10.22#v',checked='Normalization compared with the explicitly stated weighted radial convention.')],
    new_derivations='Actual covariance/receiver jets, whole-space harmonic-slab pressure-gradient enclosure, old-family energy drain and deceleration, and a leading-envelope family with three compatible initial gains.',
    verification=dict(python=sys.version,mpmath_version=mpmath.__version__,sympy_version=sympy.__version__,finite_check_programs_passed=None,new_interval_certificates_recomputed=False,formal_kernel_replay=False,NS_time_evolution=False,new_family_mean_initial_gain=True,new_family_fluctuation_energy_initial_gain=True,new_family_full_pressure_initial_gain=True,useful_numerical_stage_duration_certified=False,same_solution_return_proved=False,NS_blowup_proved=False))
(dest/'SOURCES.json').write_text(json.dumps(sources,indent=2)+'\n')
print(dest)

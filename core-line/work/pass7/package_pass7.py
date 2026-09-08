#!/usr/bin/env python3
"""Package pass7; leave historical pass2--pass6 unchanged."""
from pathlib import Path
import shutil,json,hashlib,sys
import mpmath,sympy
root=Path(__file__).resolve().parents[2]
source=root/'work/pass7';top=root/'outputs/euler-ns-transfer';dest=top/'pass7'
(dest/'checks').mkdir(parents=True,exist_ok=True)
notes=['common-family-feedback.md','D-pressure-kernel.md','D-pressure-independent-audit.md','E-interaction-reduction.md','E-certificate-independent-audit.md','outer-E-interval-independent-audit.md']
for name in notes:shutil.copy2(source/name,dest/name)
for name in ['next-stage-shear-diagnostic.md']:
    if (source/name).exists():shutil.copy2(source/name,dest/name)
checks=['verify_D_pressure_kernels.py','certify_D_viscosity.py','D-viscosity-certificate.json','verify_E_interaction_reduction.py','certify_outer_E.py','outer-E-certificate.json','verify_E_certificate_audit.py','verify_common_family.py','verify_next_stage_shear.py','run_checks.py']
for name in checks:shutil.copy2(source/name,dest/'checks'/name)
(dest/'REPORT.md').write_text("""# Seventh pass: the same three-dimensional family passes both initial tests

8 September 2026. **The finite-amplitude compatibility gap from pass6 is closed.**
The same compact, smooth three-dimensional Navier–Stokes datum now has certified
favorable central pressure feedback and positive initial growth of the global
azimuthal-mean angular-momentum maximum, including viscosity.
**A useful quantitative stage, inherited return and blowup remain unproved.**

The [common-family theorem](common-family-feedback.md) applies to every seed
amplitude 3/4 <= lambda <= 1, with the earlier exact pressure retuning:

| Quantity for the actual datum, at viscosity 1/1000 | Certified result |
|---|---|
| Central pressure | p_zz(0) = -8 exactly |
| Complete pressure-feedback test | p_zz'(0)+32 < -435297/70000 (below -6.2185) |
| Initial ratio acceleration | (b/Omega)''(0) > 435297/140000 (above 3.109) |
| Mean angular momentum at an initial maximizing circle | time derivative >626/35 (above 17.885) |

These requirements now hold on **one and the same family**. The earlier tiny
inner perturbation is not being combined with a different outer datum.
The fixed-point mean estimate and smooth local existence give simultaneous
strict gains for a common short interval, but no useful numerical interval
length or prescribed finite gain has been enclosed.

## What resolved the full pressure calculation

The [D kernel](D-pressure-kernel.md) includes all core and meridional pressure
responses and proves D_nu < 4-(96/5)Cw. Its viscosity estimate requires one
new one-dimensional outward range integral.

The [E reduction](E-interaction-reduction.md) retains the nonlocal cross
pressure of the strong outer swirl and the seed. Its complete bound is
E<1/60. The proof uses an exact radial Green trial and an energy estimate
for its whole-space residual, with both radial tails included.
The 2048-panel range calculation bounds E above by less than 0.016616141.
That is an upper bound, not an approximation to the coefficient's value.

Keeping the favorable Cw term while substituting these two estimates closes
the inequality for the full amplitude interval. Independent audits check
[D](D-pressure-independent-audit.md), [E and the common interval](E-certificate-independent-audit.md),
and [the actual interval implementation](outer-E-interval-independent-audit.md).

## What remains

The axisymmetric return obstruction from pass6 remains valid.
The new candidate permits initial three-dimensional transfer, but its evolving
correlation, viscosity, pressure, core deformation and outer geometry must
still be controlled together. Positive initial derivatives do not justify
integrating a favorable Riccati inequality to blowup.

The next task is a useful quantitative interval that retains gain and leaves
a controlled inherited endpoint, followed by a genuine return or another
compatible continuation mechanism. There is no packet reset, artificial
prescribed stress, singular solution, or formal proof-kernel replay here.

Run python3 checks/run_checks.py to recompute both decisive interval
calculations and the seven finite checks. [Verification](VERIFICATION.md)
distinguishes the analytic proof, interval evidence and open obligations.
""")
if (dest/'next-stage-shear-diagnostic.md').exists():
    with (dest/'REPORT.md').open('a') as f:f.write('\nThe [next-stage shear diagnostic](next-stage-shear-diagnostic.md) separates a frozen-background winding estimate from the coupled evolution that remains to be proved.\n')
(dest/'VERIFICATION.md').write_text("""# Verification scope

The new result is the common-family initial feedback theorem and its qualitative
simultaneous local-time consequence. It does not prove an inherited return, a
fixed useful gain, an infinite sequence or blowup.

Run python3 checks/run_checks.py. It runs seven finite programs:

1. Exact full D tensors, cone signs, and energy constants.
2. Recomputed 1024-panel D viscosity integral. PASS requires R1<120000 and
   nu*d_w<4, using exact endpoint comparisons.
3. Exact E local/source/Green/residual identities, with both radial tails.
4. Recomputed 2048-panel E range integrals. PASS requires E<1/60 and the
   sufficient compatibility bound; amplitude scaling handles either sign
   of a future E upper expression.
5. Independent exact-dyadic replay and common amplitude arithmetic.
6. Both actual certificates joined with the central jets and local-time constants.
7. Auxiliary exact initial shear and frozen-model comparisons, not actual persistence.

The interval engine is inherited unchanged from pass5/certificate and uses
mpmath.iv with rational cells and outward elementary-function arithmetic.
The result also depends on the documented analytic integration identities,
local-existence facts and earlier certified profile and retuning bounds.
Displayed intervals for upper-bound expressions are not two-sided intervals
for D, E, or the exact pressure.

Three independent audits inspect the D and E derivations and the range code;
an independent E run reproduced its upper endpoint. Symbolic checks verify
finite identities. No Lean or other formal proof replay was performed.

Exploratory finite-box pressure calculations in work/pass7 helped select the
certificate, but their values and convergence are not premises of the proved
inequalities. No numerical Navier–Stokes evolution was run.

The common existence interval follows qualitatively from strict uniform
initial margins and smooth dependence on a compact parameter family.
No useful numerical duration or return-map bound has been certified.
""")
sources=dict(date='2026-09-08',inherited_manifests=[dict(path=f'../{p}/SOURCES.json',sha256=hashlib.sha256((top/p/'SOURCES.json').read_bytes()).hexdigest()) for p in ['pass5','pass6']],
    new_derivations='Full horizontal D kernel, m=4 mixed-pressure Green/residual bound, outward interval estimates and common finite-amplitude local feedback.',
    source_scope='No new external theorem or contemporary announcement is used in this pressure certificate. Earlier source attributions remain in the unchanged manifests.',
    verification=dict(python=sys.version,mpmath_version=mpmath.__version__,sympy_version=sympy.__version__,finite_check_programs_passed=None,decisive_interval_calculations_recomputed=False,formal_kernel_replay=False,NS_time_evolution=False,finite_outer_family_full_initial_gate_certified=True,finite_outer_family_mean_initial_growth_certified=True,qualitative_common_local_time_consequence=True,useful_numerical_stage_duration_certified=False,same_solution_return_proved=False,NS_blowup_proved=False))
(dest/'SOURCES.json').write_text(json.dumps(sources,indent=2)+'\n')
print(dest)

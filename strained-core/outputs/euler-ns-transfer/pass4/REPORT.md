# Strained core, meridional feedback, and a concrete numerical candidate

September 8, 2026. **The full initial-pressure test now has a favorable numerical example at positive viscosity. It is not yet a certified inequality, a reusable stage, or a Navier–Stokes blowup proof.**

The previous pass proved one finite receiving-velocity gain but exposed a leading strain that prevented a rigid-core reset. This pass enlarges the initial class to include strain, computes its actual pressure response, and tests an added outer meridional flow. No evolved solution is replaced by a fresh profile in these calculations.

## Exact results

For the specified compact radial extension of the affine core

    L=diag(-b,-b,2b)+Omega J,

the axial core pressure is exactly

    p_zz,core=-(18/7)b^2+(2/5)Omega^2.

Disjoint outer fields can tune b'=2b^2 and Omega'=2bOmega initially. Retaining the subsequent pressure response, the exact ratio test at this neutral insertion is

    (b/Omega)''(0)=-[p_zz'(0)+32b^3]/(2Omega).

The pressure profile cannot be assumed to retain its initial radial form. An explicit defect term shows why an outer pressure contribution alone does not determine this derivative. These identities and the contact term were checked independently.

A passive outer swirl loses pressure contribution relative to the growing core rotation. An explicit compact meridional flow reverses the initial swirl-channel loss by moving the outer swirl inward. Its own pressure response must be included. The full derivative reduces to eight cubic coefficients and two potentially nonzero viscous coefficients for the chosen profiles: the two core viscous coefficients vanish exactly. The core cubic coefficients also admit exact one-dimensional integral formulas.

See [the full pressure test](full-pressure-feedback-gate.md), [the reservoir and pump construction](outer-pressure-reservoir.md), and [the independent coefficient audit](coefficient-expansion-audit.md).

## The numerical lead

For one explicit C-infinity compact datum with b=Omega=1, pump amplitude c=2 and viscosity nu=0.01, the refined calculation gives

    p_zz'(0)+32 approximately -241.35.

Negative is favorable. A second equivalent integral representation gives approximately -241.42. The refined split representation changes by about 0.00105 between the two largest resolutions. The positive-viscosity contribution and meridional backreaction are included.

This is exploratory evidence. Angular truncation and radial quadrature still need rigorous error bounds; the neutral amplitude is tuned with numerically estimated profile constants. The [numerical evidence note](numerical-pressure-evidence.md) records both resolutions, controls, residuals, the misleading coarse calculation, and raw reproducible data.

## What has not been achieved

Selecting ever-smaller scalar receivers cannot manufacture singularity on a regular solution: the selected scale relation would itself force an unbounded gradient. An actual inherited-state return is still required. The enlarged strained-core class has not been shown to persist or return, and its exterior pressure supply has not been regenerated over repeated stages.

The next immediate task is to certify the full initial pressure sign for the explicit pump datum. A successful certificate would then support a quantitative finite persistence argument. Neither step would by itself prove a cascade. The broader stage requirements remain in [NEXT_TARGET.md](../NEXT_TARGET.md).

An [exact angular-error estimate](angular-residual-certificate.md) now supports that certification step: pressure modes above degree four have zero effect on the selected core functional, and their outer contribution has an explicit energy bound. The required residual norms and low-mode errors have not yet been rigorously enclosed.

The analytical notes, independent audits, symbolic checks, numerical controls and raw experiments are packaged together. [VERIFICATION.md](VERIFICATION.md) separates what each check establishes. No formal kernel replay or NS time simulation was performed, and no novelty claim is made for the local identities.

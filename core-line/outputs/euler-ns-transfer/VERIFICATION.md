# Verification record

This file records the first pass. The latest [ninth-pass verification](pass9/VERIFICATION.md) records eight successful finite programs, including two recomputed interval certificates, explicit-amplitude exact algebra and scoped nonlinear slice checks. Pass2–pass8 remain unchanged historical packages. Actual initial identities, exact reduced-model statements, numerical slice diagnostics and conditional full-field validation are distinguished. No full compact 3D time integration, useful actual duration, inherited return or blowup is established.

## Reproduce

Use Python 3.9 or newer with SymPy 1.14.0, then run from this folder:

```sh
python3 checks/run_checks.py
```

The check programs contain our finite calculations. They do not download or execute the released Lean project. They do not solve the Navier–Stokes equations numerically.

The successful packaged run is saved in [CHECK_RESULTS.txt](CHECK_RESULTS.txt). Exact package hashes are in [SHA256SUMS.txt](SHA256SUMS.txt); hashes identify bytes and do not verify mathematical claims.

## What was checked

| Work | Verification | Scope |
|---|---|---|
| Phase heat factorization | Analytic proof, independent review, exact residual of a three-mode example with noncommuting time-dependent component matrices | Auxiliary equation with phase-independent coefficients |
| Reset-versus-loss control | Exact ODE solution and endpoint identities; explicit strict analytic amplitude inequality | One finite model pulse |
| Viscous circulation/reduced-vorticity operators | Two separately written symbolic derivations and analytical review | Coordinate identities |
| Principal damping and scale conversion | Characteristic polynomial, ratio identities, two-generation exponent, independent source/limit review | Retained Euler coefficient history and published scales |
| Thin-ring condition | Independent dimensional and algebraic review | Controlled-profile necessary condition |
| Nested-core proposal | Rational exponent checks and independent energy/dissipation/residual-scaling review | A scale budget, not a PDE stage |
| Source theorem scope | Selected files pinned to a public commit; exported propositions compared textually | No type checking or proof replay |

## Corrections made during review

- Required continuous diffusion coefficients for the classical phase equation.
- Restricted full Sobolev exponential decay to mean-zero phase data; differentiated seminorms need no such restriction.
- Described the exact control as a contracting reset pulse, not a demonstrated repeated cycle.
- Explicitly retained heat-induced distortion of the source's linear phase cell and the separate need for pointwise lower bounds.
- Required assembled-force residuals to include joining and cross-stage interactions.
- Retained potentially large carrier self-interactions; the strain/diffusion scale comparison alone does not justify linearization.
- Kept smooth forced Clay C/D constructions within scope.

No blocking error was found in these scoped analytical reviews. They are internal reviews of the stated calculations, not external mathematical certification.

## Source handling

The Euler PDF was inspected through extracted text and rendered pages where superscripts or geometry mattered. In particular, equation (12.1)'s parameter thresholds are powers \(2^{18}\) and \(2^{20}\), not the integers 218 and 220 produced by flattened text extraction. The asymptotic proof does not depend on sampling those thresholds.

Nine selected public proof files were downloaded at commit `d0124689230b58b4f86e7b90ac59de06404b3b6b`. The main Challenge and Solution propositions match after whitespace normalization. This does not establish equality of all imported definitions or a closed proof dependency graph. The public repository's build/axiom claims remain the maintainers' claims for purposes of this package.

The paper's stronger geometric details and the exported Lean existence statement are kept distinct. No ordinary Navier–Stokes theorem is inferred from an Euler proof.

See [SOURCES.json](SOURCES.json) for source URLs, inspected-byte hashes and local research references.

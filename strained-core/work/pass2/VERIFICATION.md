# Second-pass verification scope

Four standalone Python programs accompany this pass. Three use SymPy; the correction-window check uses the Python standard library. Run `python3 checks/run_checks.py` from any directory. The supplied CHECK_RESULTS.txt records the packaged run.

- `verify_active_core_pressure.py`: core divergence and nonlinear pressure source; distributional angular coefficient; fourth-derivative pressure tensor; sign on the cone; initial vorticity identities.
- `verify_local_pressure_audit.py`: an independent derivation of the same pressure signs and constants.
- `verify_affine_stretch_gate.py`: exact Gaussian and Kelvin model residuals, gain formulas and the Reynolds/localization identity.
- `verify_correction_window.py`: exact dimensional exponents, the displayed finite scalar recursion, Catalan majorant, conditional asymptotic schedules and a model mean-residual falsifier.

The active-core local-time conclusion also uses the analytical argument in active-core-pressure.md and standard local smooth NS existence for compact smooth data. The software checks its algebra; they are not a formal proof of local existence or of an evolution estimate. The two pressure calculations were independently reviewed in local-pressure-independent-audit.md. The correction argument was independently reviewed in correction-window-independent-audit.md. The short exact stress/cross-frequency formulas in mean-feedback.md received an independent algebraic review; they are diagnostics rather than solution estimates.

The Jeong–Yoneda note extracts a published PDE result and derives an explicit scaling conversion. This pass does not audit its entire proof. The current-publication search is bounded, with access limitations recorded in source-check.md.

No time-dependent numerical simulation, complete Lean build, new formal kernel proof, quantified NS return map, full-force viscous estimate, or infinite cascade was established in this pass. Conditional hypotheses and exact model identities are labeled at their point of use. SHA256SUMS.txt records the delivered file contents; SOURCES.json records the downloaded source and inherited pinned provenance.

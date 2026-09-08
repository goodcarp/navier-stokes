# Verification and limitations of pass four

Run `python3 checks/run_checks.py` from this directory or by absolute path. SymPy, NumPy and SciPy are required. CHECK_RESULTS.txt records the packaged run.

The symbolic checks cover the compact vector-potential extension; divergence and symmetry; exact pressure coefficients and distributional contact; normalized core/receiver identities; passive and pumped outer pressure derivatives; the full neutral feedback formula; coefficient polarization and support/parity cancellations; exact core viscosity zeros; radial core pressure and cubic reductions; and the independent second-jet viscosity identity. These are finite algebraic checks combined with the analytical arguments in the notes.

The numerical core-control check compares the general whole-space Legendre evaluator with a separately derived one-dimensional pressure calculation. It checks refinement and finite-harmonic controls. It does not place rigorous bounds on the outer-pump quadrature.

The tenth check verifies the degree-four cutoff in the core's pressure-error coupling and the prefactors in the high/low angular residual estimates. The accompanying analytic energy argument does not evaluate or enclose the numerical residual norms. Its initial checker needed canonicalization of an exactly zero symbolic rational sum; the corrected checker uses exact monomial integration and rational simplification. No mathematical coefficient changed.

The full pump experiments are saved separately under experiments/. They evaluate the initial derivative only. Two mathematically equivalent source/stress representations share an approximate pressure solve. Refinement, agreement, and residual checks are evidence of numerical reliability, not certified sign bounds. There was no NS time evolution and no independent kernel proof replay.

The following remain unverified:

- An interval or analytical bound certifying the explicit pump datum's favorable total pressure sign at positive viscosity.
- Quantitative persistence of that sign and a useful finite transfer interval for this enlarged class.
- Full terminal profile and exterior-environment control on inherited data.
- A same-solution return, infinite compatible stage sequence, admissible assembled force if used, and a singular solution.

File hashes identify artifacts; successful checker exits do not certify these missing claims. Earlier pass2/pass3 packages are retained unchanged.

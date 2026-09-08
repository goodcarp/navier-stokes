# Reproduce the claim you intend to check

Verification is layered. Start with the archive's integrity, then select
the finite calculation or experiment relevant to the claim. Passing one
layer does not silently discharge the others.

For the [pass12 addendum](../core-line/outputs/euler-ns-transfer/pass12/VERIFICATION.md),
`python3 tools/check_pass12_addendum.py` verifies the package manifests
and runs its four bounded programs in a temporary sibling-package copy.
It leaves the archived sources/results unchanged. The programs check
exact finite algebra and manufactured floating reconstruction cases;
they do not enclose an actual NS residual. Download the data assets first
to include their hashes in the package-integrity check.

## 1. Check the collected record

From the repository root:

```sh
python3 tools/check_archive.py
```

The checker verifies the preserved original files, every copied source
hash and local links in the new reader guides. It also runs a narrow
credential-pattern scan. Its JSON report records historical portability
warnings without modifying those source documents. It does not verify
remote links, Markdown anchors or mathematical claims.

[The source manifest](../manifests/collected-sources.json) distinguishes
files stored in Git from large local/release data assets. Download the
[data assets](data-assets.md) to the listed paths before a complete
integrity check if they are not already present in a fresh clone.

## 2. Run a packaged finite check

The current [core-line verification guide](../core-line/outputs/euler-ns-transfer/VERIFICATION.md)
links each pass's runner and exact scope. Prefer these packaged scripts
over working drafts. For example:

```sh
python3 core-line/outputs/euler-ns-transfer/pass9/checks/run_checks.py
python3 core-line/outputs/euler-ns-transfer/pass10/checks/run_checks.py
python3 core-line/outputs/euler-ns-transfer/pass11/checks/run_checks.py
```

Read each pass's `VERIFICATION.md` before running it. Some checks replay
stored rational interval endpoints; others recompute integrals or test
manufactured numerical fields. Expensive recomputations and production
evolutions are separately identified. Follow the packaged runner's use
of temporary copies where it regenerates output, so the hashed originals
remain unchanged.

The recorded Python work uses Python 3 with NumPy, SciPy, SymPy, mpmath
and, for figures, Matplotlib. Exact versions and the numerical environment
are in the relevant check or reproduction records; there is no claim that
every historical draft has identical dependency requirements. A recorded
mpmath interval certificate and a symbolic algebra check are not the same
trust object as a Lean theorem.

The [historical 94-run report](../strained-core/REPLICATION.md) records a
separate reproduction of the old pass2–9 snapshot. Its original included
and excluded runs, timing fields and discrepancies remain visible. Newer
files do not inherit that reproduction claim.

## 3. Inspect the Lean cores at their precise scope

The [Lean README](../lean/README.md) supplies its pinned Lean/Mathlib
versions, build commands and a list of what is outside the formalization.
The stored `AXIOMS.txt` corresponds to 91 elementary statements. This
archive publication preserves that evidence; it does not claim a new
independent build or a formalized Navier–Stokes PDE theorem.

The [external Euler replication note](../replication/README.md) is a
different project. Its recorded in-progress status must not be replaced
by an assumed successful build or comparator replay.

## 4. Reproduce a numerical trajectory without changing its question

Use the exact parameter block, solver source hash and checkpoint metadata
recorded with the pass11 experiment. The selected full datum has amplitude
1020 and leading-seed amplitude 1/4; it is not the earlier pressure-neutral
family or a planar slice. The solver's finite radial wall and axial period
are part of the numerical problem and must remain in the report.

Compare full complex endpoint fields, generated-mode energies, receiving
observables and core extraction as described in the
[trajectory criteria](../core-line/work/pass11/trajectory-diagnostic-criteria.md).
Time agreement alone leaves spatial, domain and whole-space residual
errors unresolved. A rational error-propagator test with manufactured
inputs is not a propagation certificate for these numerical fields.

## 5. Treat a failed check as information

A missing dependency, changed hash, failed assertion or unresolved
convergence comparison should remain in the record with its actual scope.
Do not repair it by weakening the scientific target, substituting a
different datum, resetting inherited errors, or reclassifying a model
calculation as the full PDE. New mathematical repairs belong in new
dated revisions with their own evidence.

# Pass12 verification and reproduction

The mathematical results are written conditional proofs with independent
review and finite algebra checks. No Lean kernel replay, actual NS
residual enclosure, useful whole-space stage or return is claimed.

The source is in [experiments](experiments/). Run scripts from that
directory with Python, NumPy, SciPy and SymPy. The adjacent pass11 source
package supplies the unchanged initial-field evaluator and MAC operator.

```
python3 check_c4_adapted_reference.py
python3 verify_comoving_gradient_decay.py
python3 check_compact_potential_reconstruction.py
python3 check_continuous_potential_fit.py
```

The first two programs check selected exact identities and constants, not
the complete PDE arguments. The potential tests compare manufactured
fields, dense QR and compact banded systems, verify axis/zero-extension
identities and the axial-window commutator, and construct a counterexample
to inference from sampled errors. They do not validate the stored actual
trajectory or certify high-derivative quadrature. The continuous test
includes a different radial mapping to catch hidden mapping assumptions.

To reproduce a corrected endpoint reconstruction using the full pass11
half-time-step checkpoint:

```
python3 continuous_potential_fit.py --checkpoint ../../pass11/data/resolved-192-2049-dt50.npz --intervals 144 --output continuous-fit-144.json
python3 observe_compact_reconstruction.py continuous-fit-144.npz --output continuous-observables-144.json
python3 observe_compact_reconstruction.py continuous-fit-144.npz --nz 8199 --gauss 40 --output continuous-observables-144-refined.json
```

Use `--intervals 96` for the radial trial-space comparison. The original
`compact_potential_reconstruction.py` command uses the rejected sampled
least-squares selector and is retained for reproducing that failure. Do not
promote it because its manufactured checks pass.

Four compressed potential datasets are stored in [data](data/), with a
separate SHA-256 manifest. Both sampled fits and both corrected continuous
fits are preserved. They are reconstruction coefficients, not additional
time-integrated NS fields. Their input checkpoint SHA-256 is
`603663cf57011ece8f13028fc6175f74a4dfd9cb867fb90e3b791c72aacd1940`.

The initial sampled-fit runs used a status string saying “smooth”. That
metadata was corrected to **C6/H7**, with no numerical changes. Exact
executed source versions are preserved in [provenance](provenance/), and
the reports record their hashes. Later changes added mapping propagation
and portable sibling-package imports. Numerical reports identify the
source revision actually used; the archive manifest covers the final
portable source. No complete run is claimed to have used a later revision.

All reported norms remain floating quadrature/linear-algebra quantities.
Polynomial exactness describes the mathematical quadrature rule, not an
interval enclosure of its floating evaluation. The refined quadrature
comparison concerns one fixed reconstruction only. Old pass2–pass11
archive contents and manifests are preserved.

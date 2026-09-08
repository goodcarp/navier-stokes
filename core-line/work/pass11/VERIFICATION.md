# Verification and reproduction

Run `python3 checks/run_checks.py`. It executes nine bounded programs on
temporary copies, preserving the published archive:

1. Exact poloidal–toroidal potential/operator identities.
2. Manufactured shared-frequency Hankel projection and reconstruction.
3. Compatible MAC divergence, pressure, viscosity, curl and nonlinearity.
4. Fast/reference equivalence and independent signed-frequency Stokes solves.
5. Exact-rational five-component error propagation (44 tests).
6. Independent rational propagator audit (161 assertions).
7. Exterior harmonic energy, stiffness and pressure metric.
8. Implementable exterior MAC blocks, trace matching and axis reduction.
9. Phase-preserving endpoint comparator and incompatible-input rejection.

The exterior and potential designs are not secretly enabled in the
production MAC model. That model still has the stated no-slip radial wall
and periodic axial direction. The checkers establish exact algebra or
bounded floating-point consistency as specified by each audit. They do
not enclose a whole-space space-time residual or establish a PDE endpoint.

## Reproduce the principal numerical runs

From the directory containing `experiments/`:

```text
python3 experiments/run_resolved_mac.py --nr 192 --nz 2049 --J 4 --T .001 --dt .0001 --output baseline.json --checkpoint baseline.npz
python3 experiments/run_resolved_mac.py --nr 192 --nz 2049 --J 4 --T .001 --dt .00005 --output fine-time.json --checkpoint fine-time.npz
python3 experiments/compare_mac_endpoints.py --a-json baseline.json --a-npz baseline.npz --b-json fine-time.json --b-npz fine-time.npz --output comparison.json
python3 experiments/probe_initial_resolution.py --nr 384 --nz 4097 --output initial-refined.json
```

These are substantial computations, deliberately excluded from the short
checker runner. `run_resolved_mac.py` records source hashes, all parameters,
each time-step observation, initial representation changes, generated-mode
energies and whether the requested interval completed. Its NPZ checkpoint
contains the final complex coefficient field and grid; it is not a smooth
continuous reconstruction. The comparator rejects mismatched metadata,
grids, source hashes and initial observations, and retains phase differences.

NumPy, SciPy and SymPy are required; Matplotlib is used only for the summary
figure. FFT worker counts affect performance and are recorded explicitly.
Wall-time measurements are machine/load dependent.

## Evidence boundaries

- The old short filtered and even-grid histories are retained as numerical
  failure/calibration records, not accepted core trajectories.
- The complete MAC initial state is the sampled field after a discrete
  projection; its difference from the selected exact datum is not zero.
- Time-step agreement does not estimate radial, axial, angular, domain or
  whole-space reconstruction errors.
- The raw fast nonlinear array omits excluded Fourier outputs and cannot
  serve as a complete unresolved-product residual.
- Rational comparison examples use synthetic coefficients. No actual
  fluid-trajectory error radius is certified by them.
- Exterior Bessel integrals and matrix tests are numerical/symbolic checks,
  not interval-certified special-function or evolved-field bounds.
- Classical PDE identities, Sobolev estimates, earlier initial certificates
  and the stated analytic arguments remain in the mathematical trust base.
  No new Lean proof-kernel replay is performed in this pass.

`SHA256SUMS.txt` covers the text/code/result archive. Large endpoint arrays
have a separate manifest under `data/`; a GitHub distribution may provide
them as release assets. Their role and limitations are the same on disk
and on GitHub. Historical pass2–pass10 archives are preserved unchanged.

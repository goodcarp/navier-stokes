# Actual transport and the logarithmic record clock

The Clay goal remains unresolved. This package contains a local
estimate for actual smooth finite-energy, unforced Navier–Stokes
solutions on physical R³. It does not prove an a priori bound through
a singular time or a nonsummable record charge.

Read [LOGARITHMIC_RECORD_CLOCK.md](LOGARITHMIC_RECORD_CLOCK.md)
and the separately derived
[ENTROPY_KERNEL_AUDIT.md](ENTROPY_KERNEL_AUDIT.md).
[PRIOR_ART_SCOPE.md](PRIOR_ART_SCOPE.md) distinguishes classical
inputs, related results, and the complete argument supplied here.
No novelty or cross-vendor certification is claimed.

Reproduce the exact controls with Python and SymPy:

```bash
python3 entropy_clock_gate.py
python3 transport_time_control.py
```

The root gate checks Gaussian relative entropy, failure when its
displacement term is omitted, Nash scaling algebra, logarithmic time
integration and original Navier–Stokes scaling. The independently
written second script verifies the time-dependent translated heat
kernel, and rejects wrong drift signs and time arguments at an exact
interior witness. Its constant-in-space drift has infinite kinetic
energy and is only a scalar kernel control.

These scripts do not prove BMO, John–Nirenberg, Nash smoothing,
entropy inequalities, or continuation. The main proof and analytical
audit establish their use under the stated hypotheses. Fresh
source-bound receipts record versions, command, exit status, hashes,
stdout and stderr. SHA256SUMS binds every package file except itself;
the top-level manifest binds the return, archive and manifest.

Earlier frozen returns remain unchanged. Bridge delivery remains
OWED after automatic approval review rejected recipient/data
authorization. Local manifests are not external receipts.

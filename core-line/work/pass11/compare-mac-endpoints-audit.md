# Saved MAC endpoint comparator

`compare_mac_endpoints.py` compares two completed `run_resolved_mac.py` JSON/NPZ pairs without evolving a field or constructing a new MAC solver. Example:

```bash
python3 work/pass11/compare_mac_endpoints.py \
  --a-json RUN_A.json --a-npz RUN_A.npz \
  --b-json RUN_B.json --b-npz RUN_B.npz \
  --output COMPARISON.json
```

B supplies the denominator of relative errors; signed observable differences are B minus A. The full and nonzero-angular-mode L2 differences use `B_j-A_j` as a complex coefficient before taking its modulus. All positive angular modes have kinetic weight two, and m=0 has weight one. The reported seed norm includes generated nonzero harmonics; separate per-mode results distinguish m=4 from the others. Both absolute state norms and absolute/relative differences are returned.

The receiver is reconstructed at physical z=4 using `exp(i*k*(4-z[0]))`, including on odd axial grids. Core derivatives and the radial cubic branch observable reproduce the driver. The latter is only a branch chosen in .18<r<.3 at z=4, not a global maximum. Recomputed endpoint observables must agree with their JSON history; the comparator rejects a wrong checkpoint/JSON pairing.

Comparisons reject different grid/physical/datum parameters, allowing only dt, steps and worker-count metadata to differ. The source hashes and initial observables must match; changing either would prevent interpreting the result as an isolated time-step comparison. Each checkpoint's parameter JSON must exactly match its run report, and its actual coordinates must match the specified grid, including z origin. Completed status, final time, step-count consistency, finite arrays and real mode zero are checked.

For each endpoint the report includes K and mean-circulation gain differences against its own saved discrete initial values. It also includes comparisons to K0=6.4437717665751, meanM0=41.87025495048, b0=Omega0=1, with explicit numerical-reference labels. These are not interval bounds. Core b, Omega, beta, and b−Omega are reported as well. None of these differences bounds the whole-space residual, domain/spatial errors, or an inherited successor class.

`check_compare_mac_endpoints.py` creates temporary manufactured JSON/NPZ files on Nz=32 and Nz=33, never reads the research checkpoint, and deletes its fixtures afterward. It checks the exact phase identity `||exp(i*theta)w-w||^2=2(1-cos(theta))||w||^2`, mode weights, generated-mode energy, a known core rotation increment and the analytic receiver increment at physical z=4. The phase test keeps m=4 energy unchanged while producing a nonzero L2 state difference. There are also 28 mismatch-rejection checks across the two grids. The test passed; exact source hashes and numerical discrepancies are saved in `compare-mac-endpoints-checks.json`.

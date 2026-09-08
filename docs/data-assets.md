# Full numerical data

Large numerical arrays stay in the organized local copy and are distributed
as separate GitHub release assets. The original snapshot and later
addenda keep their own immutable asset URLs:

- [research-archive-2026-09-08](https://github.com/goodcarp/forced-route-2026-09/releases/tag/research-archive-2026-09-08)
- [research-pass12-2026-09-08](https://github.com/goodcarp/forced-route-2026-09/releases/tag/research-pass12-2026-09-08)

They are excluded from Git blobs; their exact bytes remain part of the
[source inventory](../manifests/collected-sources.json). The release must
exist and its assets must be downloaded before the large-data part of an
integrity check can pass in a fresh clone.

| Asset | Role | Local placement relative to repository root | Bytes |
|---|---|---|---:|
| [pass11-resolved-192-2049-dt100.npz](https://github.com/goodcarp/forced-route-2026-09/releases/download/research-archive-2026-09-08/pass11-resolved-192-2049-dt100.npz) | Finite-cylinder Fourier/MAC endpoint | `core-line/outputs/euler-ns-transfer/pass11/data/resolved-192-2049-dt100.npz`<br>`core-line/work/pass11/resolved-192-2049-dt100.npz` | 76614004 |
| [pass11-resolved-192-2049-dt50.npz](https://github.com/goodcarp/forced-route-2026-09/releases/download/research-archive-2026-09-08/pass11-resolved-192-2049-dt50.npz) | Finite-cylinder Fourier/MAC endpoint | `core-line/outputs/euler-ns-transfer/pass11/data/resolved-192-2049-dt50.npz`<br>`core-line/work/pass11/resolved-192-2049-dt50.npz` | 76688623 |
| [pass12-compact-fit-144.npz](https://github.com/goodcarp/forced-route-2026-09/releases/download/research-pass12-2026-09-08/pass12-compact-fit-144.npz) | Rejected sampled-fit potential coefficients | `core-line/outputs/euler-ns-transfer/pass12/data/compact-fit-144.npz`<br>`core-line/work/pass12/compact-fit-144.npz` | 44919162 |
| [pass12-compact-fit-96.npz](https://github.com/goodcarp/forced-route-2026-09/releases/download/research-pass12-2026-09-08/pass12-compact-fit-96.npz) | Rejected sampled-fit potential coefficients | `core-line/outputs/euler-ns-transfer/pass12/data/compact-fit-96.npz`<br>`core-line/work/pass12/compact-fit-96.npz` | 29871042 |
| [pass12-continuous-fit-144.npz](https://github.com/goodcarp/forced-route-2026-09/releases/download/research-pass12-2026-09-08/pass12-continuous-fit-144.npz) | Continuous-L2-fit potential coefficients | `core-line/outputs/euler-ns-transfer/pass12/data/continuous-fit-144.npz`<br>`core-line/work/pass12/continuous-fit-144.npz` | 44903899 |
| [pass12-continuous-fit-96.npz](https://github.com/goodcarp/forced-route-2026-09/releases/download/research-pass12-2026-09-08/pass12-continuous-fit-96.npz) | Continuous-L2-fit potential coefficients | `core-line/outputs/euler-ns-transfer/pass12/data/continuous-fit-96.npz`<br>`core-line/work/pass12/continuous-fit-96.npz` | 29872421 |

Use [the data manifest](../manifests/data-assets.json) for SHA-256 values.
Pass11 files are full Fourier/MAC endpoint fields of the finite-cylinder
experiment. Pass12 files are compact potential reconstructions of an
increment from the same endpoint; both rejected sampled fits and corrected
continuous-L2 fits are preserved. Their induced velocities are C6/H7.
None is a validated whole-space NS state. Keep the accompanying
JSON parameters, solver source hashes and comparison report with them.
No interpolated or independently phase-aligned replacement is the original
checkpoint.

For a private repository, download through an authenticated GitHub session
or `gh release download`, then place each file at its listed local path.
Run `python3 tools/check_archive.py` afterward. The archive checker validates
the stored hashes; it does not execute arrays or certify the PDE.

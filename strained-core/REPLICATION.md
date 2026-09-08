# REPLICATION: second-host re-run of every strained-core check script

Seat: REPLICATION, Codex half of the 2026-09-08 sitting. Date: 2026-09-08.

This file records a re-run, on a second host, of every check program in `work/pass2` through `work/pass9` and the
top-level `work/*.py`. It replicates the checks as written. It does not grade the mathematics, and a green run here
certifies nothing beyond what each pass's own VERIFICATION.md already scopes: these are finite algebra programs,
interval enclosures and exploratory quadratures, not a Navier-Stokes existence or blowup proof.

## How this was run

- Host: darwin 24.6.0, arm64 Mac, 12 cores. Python 3.9.6 (Apple, /usr/bin/python3).
- Packages present: numpy 1.26.4, scipy 1.13.1, sympy 1.14.0, mpmath 1.3.0, matplotlib 3.9.4. Nothing had to be installed.
- A Lean build was running throughout; load average reached 270 and the machine was swapping. All runs used nice -n 10, one process at a time.
- Everything ran in a scratch copy of this tree. The repository copy and the originals under `~/<core-line-working-tree>` were
  never touched; this file and `REPLICATION_RESULTS.json` are the only two files written into the repository.
- Runner: one script at a time, nice -n 10, cwd set to the pass folder, stdout/stderr/exit code/wall time captured.
- work/pass2 through work/pass6 contain no run_checks.py; the byte-identical copies from outputs/euler-ns-transfer/pass{2,3,4,6}/checks and pass5/certificate were placed in the scratch pass folders and run there. work/pass7, pass8 and pass9 have their own run_checks.py.
- Every packaged check script under outputs/euler-ns-transfer/ was diffed against its work/ counterpart before running: all identical, run_checks.py excepted. The one exception outside the pass folders is the top-level checker, where the packaged `viscosity_check.py` and the working `viscosity-check.py` differ only in the closing line's cross reference (`../REPORT.md` against `viscosity-audit.md`).
- Raw stdout, stderr, exit codes and wall times for all 94 runs were captured per script; the machine-readable summary is in `REPLICATION_RESULTS.json` alongside this file.

## Result in one line

**94 runs, 94 PASS, 0 FAIL, 16 scripts not re-run.** Every check program in the tree that is a check program exited zero. Every stored JSON certificate that a script regenerates was reproduced field for field, exact dyadic endpoints included, with the differences confined to elapsed-time fields, two prose strings and one dependency path.

## Top level, work/*.py

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `phase_heat_check.py` | Exact auxiliary phase-heat and return-model algebra: three-mode PDE residual with noncommuting B(t), damped-pulse ODE, return time, retained amplitude. | 0 | PASS | 4.1 s |  |
| `viscosity-check.py` | Viscosity-audit identities: circulation and reduced-vorticity diffusion in volume coordinates, damping in the projective ratio, two-generation frequency conversion, ordinary viscosity exponent. | 0 | PASS | 1.5 s |  |
| `verify_direction_scout.py` | Direction-scout algebra: both exact reduced diffusion operators, physical-coordinate growth exponent, frozen viscous characteristic polynomial, summability exponents. | 0 | PASS | 1.7 s |  |

## Pass 2

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `run_checks.py (from the packaged pass2/checks)` | Aggregate runner for the four pass2 programs; PASS is a zero exit and the closing line. | 0 | PASS | 19.0 s | stdout byte-identical to the packaged pass2/CHECK_RESULTS.txt |
| `verify_active_core_pressure.py` | Core divergence and nonlinear pressure source, distributional angular coefficient, fourth-derivative pressure tensor, sign on the cone, initial vorticity identities. | 0 | PASS | 3.7 s |  |
| `verify_local_pressure_audit.py` | Independent derivation of the same pressure signs and constants (core pressure Hessian coefficient, horizontal fourth derivatives, cone eigenvalues). | 0 | PASS | 2.2 s |  |
| `verify_affine_stretch_gate.py` | Exact Gaussian and Kelvin model residuals, gain formulas, Reynolds/localization identity. | 0 | PASS | 3.2 s |  |
| `verify_correction_window.py` | Exact dimensional exponents, the displayed finite scalar recursion, Catalan majorant, conditional asymptotic schedules, model mean-residual falsifier. | 0 | PASS | 10.6 s |  |

## Pass 3

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `run_checks.py (from the packaged pass3/checks)` | Aggregate runner for the four pass3 programs. | 0 | PASS | 18.3 s | stdout byte-identical to the packaged pass3/CHECK_RESULTS.txt |
| `verify_quantitative_core_stage.py` | Explicit Sobolev counting constants, derivative-norm Fourier comparison, time-derivative bound scaling, finite remainder/gain inequalities including g = k0*tau^2/4. | 0 | PASS | 1.0 s |  |
| `verify_rotational_receiver.py` | Local projected acceleration, radial harmonic moment examples, exact receiver scaling and uniform gain/radius arithmetic (15 exact checks). | 0 | PASS | 8.7 s |  |
| `verify_material_core_transfer.py` | Material deformation/area identities, rigid-rotation/Laplacian commutator, harmonic initial-jet relations used in the circulation argument. | 0 | PASS | 8.2 s |  |
| `verify_receiver_flux.py` | Receiver divergence, nonlinear angular-momentum flux, viscous weight Laplacian, harmonic quadratic coefficient, matched energy exponent. | 0 | PASS | 2.5 s |  |

## Pass 4

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `run_checks.py (from the packaged pass4/checks)` | Aggregate runner for the ten pass4 programs; runs verify_pressure_viscosity.py with --symbolic-only. | 0 | PASS | 77.4 s | stdout byte-identical to the packaged pass4/CHECK_RESULTS.txt |
| `derive_affine_core_pressure.py` | Compact vector-potential extension, divergence and symmetry, exact core pressure p_zz = -(18/7)b^2 + (2/5)Omega^2, tuned initial b' and Omega'. | 0 | PASS | 4.9 s |  |
| `verify_affine_pressure_independent.py` | Independent cylindrical rederivation of the strain pressure including the contact term, and the exterior-cone profile defect. | 0 | PASS | 2.2 s |  |
| `verify_outer_reservoir.py` | Exact swirl transport factor, cone lower bound, divergence-free pump, pressure-limit coefficient, tuned initial core ratios, neutral amplitude A^2 = (38/7 b^2 + 2/5 Omega^2 - c^2 C_P)/C_v. | 0 | PASS | 2.8 s |  |
| `verify_iteration_audit.py` | Two-step Reynolds telescoping, matched velocity power, strain-to-rotation evolution, normalization exponents. | 0 | PASS | 2.7 s |  |
| `verify_full_feedback_gate.py` | The full neutral feedback formula beta'' = -(p_zz' + 32 b^3)/(2 Omega), including the contact coefficient, and the profile-defect decomposition. | 0 | PASS | 2.2 s |  |
| `verify_coefficient_expansion.py` | Eight cubic and four viscous coefficients, support/parity cancellations, exact core viscosity zeros dS = dW = 0. | 0 | PASS | 2.8 s |  |
| `derive_core_pressure_derivative.py` | Radial core pressure and cubic reductions; the H2/H4 ODE identities. | 0 | PASS | 33.0 s |  |
| `verify_pressure_quadrature.py` | Analytic core pressure coefficients, divergence, full-harmonic trace, and convergence of the general Legendre evaluator to the independent one-dimensional cubic reduction. | 0 | PASS | 9.2 s | reproduces the recorded core control values t300 = -6.734827516735 at nr=800 and -6.734553125174 at nr=1600 |
| `verify_angular_residual_certificate.py` | Three core weights, boundary cancellation, degree-four cutoff, angular Poincare eigenvalues, outer high/low mode factors. | 0 | PASS | 10.3 s |  |
| `evaluate_core_pressure_radial.py` | Independent one-dimensional core control: t300, tWb, dS, dW at n = 2001..16001, and source-vs-reduced agreement. | 0 | PASS | 1.1 s | t300 = -6.734463311168844 and tWb = 2.563122131694385 at n=2001, converging to -6.7344633111671 by n=16001; the note records -6.73446331117 and 2.563122131694. source_vs_reduced agreement about 4e-11, the note says about 5e-11 |
| `verify_pressure_viscosity.py (no --symbolic-only)` | Independent second-jet viscosity identity plus the exploratory same-grid numerical cross-check at nr=400, nmu=512, lmax=2. | 0 | PASS | 11.5 s | run without --symbolic-only, so it also performs the exploratory numerical cross-check the packaged runner skips; the mixed-amplitude additivity assertion holds to rtol 1e-9 |
| `evaluate_pressure_gate.py (defaults nr=600 nmu=256 lmax=128)` | The full whole-space initial pressure-gate evaluator at its script defaults; this is the coarse, angularly underresolved setting the note warns about. | 0 | PASS | 14.0 s | reproduces the coarse, angularly underresolved setting; poisson_trace_relative 0.5669, representation_difference up to -18981 at c=4, exactly the misleading regime the note documents |
| `evaluate_pressure_gate.py --nr 1200 --nmu 1024 --lmax 768 --c 0 1 2` | Coarser of the two refined resolutions used for the numerical lead. | 0 | PASS | 41.7 s | every field matches the stored pressure-split-1200.jsonl to all printed digits |
| `evaluate_pressure_gate.py --nr 1600 --nmu 1536 --lmax 1280 --c 0 1 2` | The documented reproduce command for the pass4 numerical lead p_zz'(0)+32 at c=2, nu=0.01. | 0 | PASS | 1155.9 s | every field matches the stored pressure-split-1600.jsonl to all printed digits; 19.3 min wall time because the host was swap-bound |

## Pass 5

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `run_checks.py (from the packaged pass5/certificate)` | Aggregate runner for the twelve pass5 verify programs. | 0 | PASS | 55.6 s | every per-check line identical to the packaged pass5/CHECK_RESULTS.txt; only the runner's own closing summary line differs (the packaged file was written by the packaging script) |
| `verify_actual_return_next_lemma.py` | Exact quartic-pressure lower bound 3693184/105861 > 34; RMS normalization and cubic core jet. | 0 | PASS | 3.3 s |  |
| `verify_angular_tail_budget.py` | Cartesian fourth-derivative kernel to spherical block, radial/tangent degree structure, angular Sobolev weights. | 0 | PASS | 10.7 s |  |
| `verify_core_exterior_pressure_response.py` | Angular contraction and contact terms by degree, radial moments, and the derivative form (18/7) b p_zz - (b moment/15) p_zzzz. | 0 | PASS | 5.6 s |  |
| `verify_finite_feedback_corollary.py` | beta and excess-strain initial jets, quotient third-derivative majorant, uniform receiving-mode restrictions. | 0 | PASS | 1.7 s |  |
| `verify_full_gate_bound.py` | Aggregation of the saved exact dyadic endpoints into the full gate strict upper bound -290649/8750 and beta'' lower bound 290649/17500. | 0 | PASS | 0.1 s |  |
| `verify_inherited_geometry.py` | Azimuthal equation and adjoint, moment rates, moving-axis pump rates, initial packet margins, outer-stress rate > 226249/51000. | 0 | PASS | 3.2 s |  |
| `verify_interval_cutoff.py` | Interval consistency: plateau and midpoint cutoff jets, cutoff mass 5/8 enclosed, derivative integral -1 enclosed. | 0 | PASS | 3.6 s |  |
| `verify_local_swirl_feedback_kernel.py` | Signs and coefficients of Qtheta, K_loc, S_theta; K/Q and core/Cv ranges; exact annular cancellation tWc = 0. | 0 | PASS | 2.9 s |  |
| `verify_plateau_pump_bound.py` | Local centripetal and swirl-advection kernel, exact multiplier upper bound -15.53532849417749 < -31/2. | 0 | PASS | 4.4 s |  |
| `verify_radial_annular_pressure_cutoff.py` | Fourth-derivative Newton kernel, exact q0/q2/q4, degree-four polynomial identity, radial functional coefficient. | 0 | PASS | 4.9 s |  |
| `verify_robust_initial_data.py` | Non-neutral central second derivatives and the rational A = 1100 margins read from local-interval-results.json. | 0 | PASS | 1.8 s |  |
| `verify_swirl_adjoint_kernel.py` | Three exact adjoint radial Poisson equations, radial Green symmetry, plateau adjoints and w^2 kernel coefficients. | 0 | PASS | 4.3 s |  |
| `certify_core_twb.py (default 1024 panels, 35 dps)` | Interval enclosure of the tWb coefficient; PASS requires tWb strictly below four. | 0 | PASS | 5.3 s | 26 of 26 fields identical to the stored core-twb-interval.json |
| `certify_radial_strain.py (default c=7/5, 1024 panels)` | Interval enclosure of the combined radial meridional cubic coefficient in three stages (core, inner, outer). | 0 | PASS | 9.7 s | 25 of 25 fields identical to the stored radial-strain-interval-1024.json apart from the timing field |
| `certify_full_initial_gate.py` | Aggregation only of the stored component bounds into the full initial gate. | 0 | PASS | 0.2 s |  |
| `certify_full_initial_gate.py --recompute` | Independent regeneration of all four integral enclosures (radial meridional cubic, core strain/swirl, positive outer swirl stress, its viscous coefficient) and reassembly of the gate. | 0 | PASS | 64.2 s | 103 of 103 certificate fields identical to the stored FULL_INITIAL_GATE_CERTIFICATE.json, including every exact dyadic mantissa and exponent; only the three per-integral timing fields differ |
| `evaluate_radial_pump.py (defaults nr=1200 nmu=1024 lmax=4)` | Exploratory radial annular pump evaluation; not an input to the certificate. | 0 | PASS | 53.4 s | run at the script defaults nr=1200; the stored radial-pump-800.jsonl was produced at nr=800, so the values differ by resolution as expected. Core pressures reproduce -18/7 and 2/5 and C_v agrees with the 1200-node pressure-gate value to 10 digits |

## Pass 6

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `run_checks.py (from the packaged pass6/checks)` | Aggregate runner for the seven pass6 verify programs. | 0 | PASS | 16.5 s | every per-check line identical to the packaged pass6/CHECK_RESULTS.txt; only the closing summary line differs |
| `verify_circulation_return_obstruction.py` | Circulation equation and scaling, receiver lower bounds, forcing scaling, conditional meridional recurrence and normalization. | 0 | PASS | 2.5 s |  |
| `verify_escape_route_audit.py` | The escape-route arithmetic for the normalized axisymmetric classes. | 0 | PASS | 2.6 s |  |
| `verify_nonaxisymmetric_torque.py` | Fully three-dimensional angular equation and exact Fourier-mode covariance. | 0 | PASS | 4.2 s |  |
| `verify_outer_mean_maximum.py` | Recomputes the 64-panel cutoff-jet bounds with the pass5 interval engine and emits the outer mean-maximum certificate. | 0 | PASS | 0.7 s | 25 of 25 fields identical to the stored outer-mean-maximum-certificate.json |
| `verify_outer_pressure_compatibility.py` | The amplitude-polynomial compatibility reduction for the outer finite-amplitude family. | 0 | PASS | 3.6 s |  |
| `verify_rotational_receiver_torque.py` | Selected-sign receiver torque and its frequency-independent energy bound. | 0 | PASS | 5.9 s |  |
| `verify_three_dimensional_initial_gate.py` | H4 pressure-pairing counting constant, pressure-functional Lipschitz polynomial, exact neutral retuning, discrete-symmetry central identities. | 0 | PASS | 2.1 s |  |

## Pass 7

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `run_checks.py` | The seven pass7 programs in order, recomputing the 1024-panel D and 2048-panel E interval certificates and asserting their PASS status. | 0 | PASS | 72.2 s | every per-check line identical to the packaged pass7/CHECK_RESULTS.txt; differs only in the elapsed-seconds field and the closing summary wording (the packaged file predates the seventh program being added to the list) |
| `verify_D_pressure_kernels.py` | Exact full D tensors, cone signs and energy constants; conservative horizontal bounds core < -12/5 and pump < -12. | 0 | PASS | 10.3 s |  |
| `certify_D_viscosity.py --panels 1024` | Recomputed 1024-panel D viscosity integral. PASS requires R1 < 120000 and nu*d_w < 4 by exact endpoint comparison. | 0 | PASS | 5.2 s | 27 of 27 fields identical to the stored D-viscosity-certificate.json |
| `verify_E_interaction_reduction.py` | Exact E local/source/Green/residual identities with both radial tails. | 0 | PASS | 4.6 s |  |
| `certify_outer_E.py --panels 2048` | Recomputed 2048-panel E range integrals. PASS requires E < 1/60 and the sufficient compatibility bound. | 0 | PASS | 18.2 s | 147 of 147 fields identical to the stored outer-E-certificate.json apart from the timing field |
| `verify_E_certificate_audit.py` | Independent exact-dyadic replay of the E certificate and the common amplitude arithmetic. | 0 | PASS | 1.3 s |  |
| `verify_common_family.py` | Joins both actual certificates with the central jets and local-time constants: T_lambda < -435297/70000 for every lambda in [3/4,1]. | 0 | PASS | 2.6 s |  |
| `verify_next_stage_shear.py` | Auxiliary initial-shear and frozen-model comparisons, kept separate from actual-solution persistence. | 0 | PASS | 3.0 s |  |
| `verify_alternative_direction.py` | Reflection signs, Q kernels, local E, projection split, carrier optimum and mode selection for the alternative direction. | 0 | PASS | 5.0 s |  |
| `certify_outer_adjoint_norm.py --panels 18000` | Interval enclosure of nd_squared = \|\|r div(Qw)\|\|^2. Not part of run_checks; its recorded status is VALID_ENCLOSURE_TARGET_NOT_MET both stored and here. | 0 | PASS | 111.9 s | 29 of 30 fields identical to the stored outer-adjoint-norm-certificate.json; only elapsed_seconds differs. Run at 18000 panels to match the stored certificate rather than the script default of 24000. Its status string is VALID_ENCLOSURE_TARGET_NOT_MET in both, an expected outcome, not a failure of the run |
| `evaluate_alternative_budget.py` | Exploratory midpoint quadrature of the alternative pressure budget; explicitly not an enclosure. | 0 | PASS | 3.6 s |  |
| `evaluate_horizontal_trial.py` | Exploratory separable whole-space pressure trial; explicitly not a certified bound. | 0 | PASS | 7.8 s |  |
| `evaluate_outer_E.py` | Exploratory finite-box E quadrature; not a premise of the proved inequalities. | 0 | PASS | 3.3 s |  |

## Pass 8

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `run_checks.py` | The eleven pass8 programs, regenerating the 2048-panel fixed-circle certificate and the two 1024-panel leading-envelope certificates and asserting PASS on each. | 0 | PASS | 94.3 s | every per-check line identical to the packaged pass8/CHECK_RESULTS.txt except the elapsed-seconds field and one certificate scope string (see the certificate section) |
| `verify_covariance_evolution.py` | Covariance production and viscosity, pressure/axial flux collapse, fixed-circle Gtt, exact global fluctuation-energy production. | 0 | PASS | 7.8 s |  |
| `verify_covariance_pressure_slab.py` | Harmonic-slab full-pressure bound: saved residual ne < 85/4, whole-space \|pressure_r - P'\| < 5, exact derivative-reduced Green gradient. | 0 | PASS | 3.6 s |  |
| `certify_fixed_circle_acceleration.py --panels 2048` | The 2048-panel radial-root / Green-moment / fourth-jet certificate. PASS requires the moving-radial-maximum Gtt upper below -240000, Gt in (85,186), radial speed above 1/5. | 0 | PASS | 43.5 s | 138 of 138 fields identical to the packaged pass8/checks copy apart from timing; against the work/pass8 copy one prose scope string differs because that file predates the final wording |
| `verify_fixed_circle_acceleration_audit.py` | Independent audit of the logistic/g jets through order four and the saved exact initial-jet inequalities. | 0 | PASS | 16.0 s |  |
| `verify_actual_fluctuation_energy.py` | Full initial strain contraction and the rational fractional fluctuation-energy drain. | 0 | PASS | 4.9 s |  |
| `kelvin-covariance-check.py` | The separate Kelvin model: exact rotating Kelvin equations, polarization, viscosity, stress depletion, integrated budget. | 0 | PASS | 11.4 s |  |
| `verify_endcap_gradient_bound.py` | Exact endcap estimate: positive Bessel weight, Hankel coefficient, mode Green jump, tails, slab prefactor. | 0 | PASS | 5.9 s |  |
| `exterior-dynamics-verify.py` | Exact exterior-pressure identities: Newton derivative kernels, harmonicity, radial density, reflected-packet bounds. | 0 | PASS | 4.0 s |  |
| `certify_leading_envelope.py --panels 1024` | The 1024-panel leading-envelope torque / extraction / viscosity certificate. | 0 | PASS | 17.3 s | 65 of 66 fields identical to the packaged pass8/checks copy; the one difference is the relative path recorded in cutoff_derivative_certificate_dependency, which reflects the work tree layout rather than the package layout |
| `certify_leading_pressure.py --panels 1024` | The 1024-panel leading-envelope full-pressure norm certificate; PASS requires E_absolute_upper strictly below 1/3 and nu*d below 7. | 0 | PASS | 16.3 s | 76 of 76 fields identical to both the work and packaged copies |
| `verify_leading_family_independent_audit.py` | Independent same-family audit: non-flat leading covariance/torque, exact retuning and energy margins, full-pressure E < 1/3, common actual initial T < -11. | 0 | PASS | 5.3 s |  |
| `exterior-dynamics-pressure.py (default n=64)` | Exploratory whole-space endcap quadrature feeding exterior-dynamics-verify; explicitly not an error enclosure. | 0 | PASS | 4.5 s | run at the script default n=64; the stored artifacts are n=48 and n=96, so this is an additional grid rather than a comparison |
| `evaluate_mean_maximum_acceleration.py` | Exploratory finite-box pressure evaluation of the mean-maximum acceleration; no actual finite duration follows. | 0 | PASS | 11.5 s |  |

## Pass 9

| Script | What it establishes, per VERIFICATION.md | Exit | Verdict | Wall | Note |
|---|---|---|---|---|---|
| `run_checks.py` | The eight pass9 programs, recomputing the general-envelope jet interval certificate and the affine-envelope budget, and asserting PASS on the three stored certificates. | 0 | PASS | 86.9 s | every per-check line identical to the packaged pass9/CHECK_RESULTS.txt except the elapsed-seconds field |
| `verify_general_envelope_jet.py` | General-envelope covariance, pressure-flux and mean-jet algebra, and the exact flat-envelope limit. | 0 | PASS | 12.3 s |  |
| `certify_general_envelope_initial_jet.py --panels 2048` | Newly recomputed full-pressure initial-jet interval enclosure. | 0 | PASS | 9.7 s | 227 of 227 fields identical to the stored general-envelope-initial-jet-certificate.json apart from the timing field |
| `verify_general_envelope_initial_jet_certificate.py` | Exact rational replay of the jet certificate bounds and dependency paths: Gt > 11, actual Gtt < -13000, tracked Gtt < -12000, pressure error < 1.3. | 0 | PASS | 0.1 s |  |
| `verify_evolution_validation_bridge.py` | H4/H7 residual majorants, core formulas, endpoint-error algebra, explicit Riccati supersolution. | 0 | PASS | 2.2 s |  |
| `verify_affine_envelope_budget.py --panels 2048` | Exact affine-envelope finite-time gain algebra and its new interval moment; emits the affine-envelope budget certificate. | 0 | PASS | 4.1 s | 27 of 27 fields identical to the stored affine-envelope-budget-certificate.json |
| `verify_fixed_amplitude_initialization.py` | Fixed-amplitude off-neutral core and initial-gain rational inequalities at A = 1020, lam = 1/4. | 0 | PASS | 2.0 s | 17 of 17 fields identical to the stored fixed-amplitude-initialization-certificate.json |
| `check_slice_solver.py` | Analytic-solution, incompressibility and circulation checks for the slice solver at nr = 192, 384, 768, requiring second-order convergence ratios in (3.5, 4.5). | 0 | PASS | 20.3 s | 17 of 17 fields identical to the stored slice-solver-checks.json |
| `compare_slice_runs.py` | Sensitivity comparison of the stored nonlinear slice runs (radial/time refinement, angular and outer-boundary comparisons, selected amplitude). | 0 | PASS | 0.5 s | 142 of 142 fields identical to the stored slice-comparison.json |
| `affine-envelope-spectrum.py (defaults n=8192, length=2)` | Exploratory Fourier quadrature of the exact linear rotating-affine envelope model. | 0 | PASS | 11.9 s | 70 of 70 fields identical to the stored affine-envelope-spectrum-8192.json |
| `evolve_circular_slice.py (defaults nr=768, nphi=48, dt=1e-6, A=950, T=0.002)` | Nonlinear planar slice diagnostic in the prescribed affine strain; explicitly not compact 3D NS. | 0 | PASS | 61.7 s | run at the script defaults (T=0.002). Every observable at every sampled time it shares with the stored slice-768.json (21 times, t=0 through t=0.002) is bit-identical, including mean_maximum, fluctuation energy, total energy and torque |
| `affine-envelope-spectrum.py --n 16384 (length left at default 2)` | Resolution variant; superseded by the --length 4 run that matches the stored artifact. | 0 | PASS | 7.5 s | run at the script default length=2, so it is not comparable with the stored length=4 artifact; kept only as a resolution variant |
| `affine-envelope-spectrum.py --n 16384 --length 4` | The parameters of the stored affine-envelope-spectrum-16384.json artifact. | 0 | PASS | 8.8 s | 63 of 70 fields identical to the stored affine-envelope-spectrum-16384.json. The seven differences are two grid-index selections and the values they carry: the stored file records the sample at t=0.00024 where this run records t=0.00025, and the first joint-sign break at t=0.0008 where this run gives t=0.00079. Every parameter, energy, enstrophy, cap, peak ratio and peak time matches bit for bit |

## Scripts not re-run

| Script | Pass | Reason |
|---|---|---|
| `work/fetch_proof_scope.py` | top level | Not a check script. It downloads pinned source text over the network and writes into work/sources. Skipped as out of scope and network dependent. |
| `work/package_transfer.py` | top level | Packaging script; writes into outputs/ and an absolute archive path outside the tree. Skipped per instruction. |
| `work/pass2/package_pass2.py` | pass2 | Packaging script. Skipped per instruction. |
| `work/pass3/package_pass3.py` | pass3 | Packaging script. Skipped per instruction. |
| `work/pass4/reduce_core_cubic_scratch.py` | pass4 | Not one of the named check patterns; a scratch derivation that loads derive_core_pressure_derivative.py through the repo-root relative path work/pass4/..., so it only runs from the tree root. Skipped. |
| `work/pass5/interval_local_integrals.py` | pass5 | Library module holding the interval engine. It is imported and exercised by certify_full_initial_gate.py, certify_radial_strain.py, certify_core_twb.py and the pass6 to pass9 certifiers, all of which were run. |
| `work/pass5/package_pass5.py` | pass5 | Packaging script. Skipped per instruction. |
| `work/pass6/package_pass6.py` | pass6 | Packaging script. Skipped per instruction. |
| `work/pass7/finalize_pass7.py` | pass7 | Finalize script. Skipped per instruction. |
| `work/pass7/package_pass7.py` | pass7 | Packaging script. Skipped per instruction. |
| `work/pass7/update_authoritative_state.py` | pass7 | Skipped per instruction. |
| `work/pass8/finalize_pass8.py` | pass8 | Finalize script. Skipped per instruction. |
| `work/pass8/package_pass8.py` | pass8 | Packaging script. Skipped per instruction. |
| `work/pass9/finalize_pass9.py` | pass9 | Finalize script. Skipped per instruction. |
| `work/pass9/package_pass9.py` | pass9 | Packaging script. Skipped per instruction. |
| `work/pass9/plot_slice_pulse.py` | pass9 | Figure generation, not a check; writes into pass9/figures. |

The brief flagged `evolve_circular_slice.py` and `affine-envelope-spectrum.py` as candidates for exclusion on runtime.
A dry read put both well inside the ten-minute budget, so both were run: the slice evolution took 61.7 s at its
defaults and the spectrum quadrature 11.9 s. Neither needed to be skipped.

## Failures

None. Every executed program exited zero.

Two results are worth naming so they are not mistaken for failures:

- `pass7/certify_outer_adjoint_norm.py` reports `status: VALID_ENCLOSURE_TARGET_NOT_MET`. That is the stored
  status too. The enclosure is valid; the norm target `9/2500` is not reached, and the script records the weaker
  rational bound `773/200000` it does prove. It is not part of `run_checks.py` and its failure to reach the
  target is a recorded fact of the pass, not a replication defect.
- `pass4/evaluate_pressure_gate.py` at its script defaults (nr=600, nmu=256, lmax=128) produces the
  angularly underresolved numbers the pass explicitly disowns: a relative Poisson trace residual of 0.5669 and a
  representation difference of -18981 at c=4. This reproduces the misleading coarse regime documented in
  `numerical-pressure-evidence.md`; it is the expected behaviour of that setting, not a discrepancy.

## Certificate comparisons

Where a script recomputes a stored JSON certificate, the recomputed output was compared field by field against the
pristine repository copy. Comparison is on every leaf, which for these certificates includes the exact dyadic
mantissa and exponent of every interval endpoint, so a match is bit-level, not a match to some number of digits.

| Stored certificate | Regenerated by | Identical / total | Differences | Verdict |
|---|---|---|---|---|
| `pass5/FULL_INITIAL_GATE_CERTIFICATE.json` | `certify_full_initial_gate.py --recompute` | 103 / 103 | three per-integral seconds fields only | MATCH |
| `pass5/core-twb-interval.json` | `certify_core_twb.py` | 26 / 26 | none | MATCH |
| `pass5/radial-strain-interval-1024.json` | `certify_radial_strain.py` | 25 / 25 | seconds only | MATCH |
| `pass6/outer-mean-maximum-certificate.json` | `verify_outer_mean_maximum.py` | 25 / 25 | none | MATCH |
| `pass7/D-viscosity-certificate.json` | `certify_D_viscosity.py --panels 1024` | 27 / 27 | none | MATCH |
| `pass7/outer-E-certificate.json` | `certify_outer_E.py --panels 2048` | 147 / 147 | seconds only | MATCH |
| `pass7/outer-adjoint-norm-certificate.json` | `certify_outer_adjoint_norm.py --panels 18000` | 29 / 30 | elapsed_seconds only | MATCH |
| `pass8/fixed-circle-acceleration-certificate.json` | `certify_fixed_circle_acceleration.py --panels 2048` | 138 / 138 | seconds only, against the packaged copy; the work/pass8 copy carries an older prose scope string | MATCH |
| `pass8/leading-envelope-certificate.json` | `certify_leading_envelope.py --panels 1024` | 65 / 66 | cutoff_derivative_certificate_dependency records the work-tree relative path instead of the package path | MATCH on all numbers |
| `pass8/leading-pressure-certificate.json` | `certify_leading_pressure.py --panels 1024` | 76 / 76 | none | MATCH |
| `pass9/general-envelope-initial-jet-certificate.json` | `certify_general_envelope_initial_jet.py --panels 2048` | 227 / 227 | seconds only | MATCH |
| `pass9/affine-envelope-budget-certificate.json` | `verify_affine_envelope_budget.py --panels 2048` | 27 / 27 | none | MATCH |
| `pass9/fixed-amplitude-initialization-certificate.json` | `verify_fixed_amplitude_initialization.py` | 17 / 17 | none | MATCH |
| `pass9/slice-comparison.json` | `compare_slice_runs.py` | 142 / 142 | none | MATCH |
| `pass9/slice-solver-checks.json` | `check_slice_solver.py` | 17 / 17 | none | MATCH |
| `pass9/affine-envelope-spectrum-8192.json` | `affine-envelope-spectrum.py (defaults)` | 70 / 70 | none | MATCH |
| `pass9/affine-envelope-spectrum-16384.json` | `affine-envelope-spectrum.py --n 16384 --length 4` | 63 / 70 | two grid-index selections near a sign crossing and the five values they carry | PARTIAL |
| `pass9/slice-768.json` | `evolve_circular_slice.py (defaults, T=0.002)` | all shared fields | none over the 21 shared sample times | MATCH on every shared observable |
| `pass4/pressure-split-1200.jsonl` | `evaluate_pressure_gate.py --nr 1200 --nmu 1024 --lmax 768 --c 0 1 2` | all shared fields | none to all printed digits | MATCH |
| `pass4/pressure-split-1600.jsonl` | `evaluate_pressure_gate.py --nr 1600 --nmu 1536 --lmax 1280 --c 0 1 2` | all shared fields | none to all printed digits | MATCH |

### The three differences worth stating plainly

1. `pass8/fixed-circle-acceleration-certificate.json` and `pass8/leading-envelope-certificate.json` as stored in
   `work/pass8` carry an earlier `scope` sentence than the scripts now emit. The packaged copies under
   `outputs/euler-ns-transfer/pass8/checks/` carry the current wording and match the recomputation exactly. So the
   `work/` JSONs are one wording revision behind their own scripts. No number is affected.
2. `leading-envelope-certificate.json` records `cutoff_derivative_certificate_dependency` as
   `../../pass6/checks/outer-mean-maximum-certificate.json` in the package and `../pass6/outer-mean-maximum-certificate.json`
   when run from `work/`. That is the script resolving its dependency against whichever layout it finds, by design.
3. `pass9/affine-envelope-spectrum-16384.json` reproduces 63 of 70 fields bit for bit under
   `--n 16384 --length 4`. The seven that differ are two grid-index selections and the values they carry: the stored
   file reports the sample at `t = 0.00024` where the current script selects `t = 0.00025`, and reports the first
   joint sign break at `t = 0.0008` where this run gives `t = 0.00079`. Both are index picks near a sign crossing;
   the stored file appears to come from a marginally earlier version of the sample list. Every physical quantity in
   that file, including initial energy, initial enstrophy, the enstrophy-energy cap, the maximum sampled energy
   ratio and the sampled peak time, matches exactly. The default `n = 8192` artifact matches on all 70 fields.

## Headline numbers independently confirmed

**pass2: p_zz,core = -(18/7) b^2 + (2/5) Omega^2** (partly here, fully in pass4)

pass2 verify_active_core_pressure.py asserts ang = -4/15, rad = -1/4 and core = 2/3 - 4*ang*rad = 2/5, the Omega^2 coefficient; verify_local_pressure_audit.py rederives the same core pressure Hessian coefficient independently. The complete formula including -18/7 b^2 is asserted in pass4 derive_affine_core_pressure.py, whose output line reads 'Core p_zz = 2*Om**2/5 - 18*b**2/7', and is rederived in pass4 verify_affine_pressure_independent.py.

**pass3: finite gain g = k0 tau^2 / 4** (confirmed)

verify_quantitative_core_stage.py computes gain = k0*tau*tau/4 over k0 in {1/10, 1, 7} with tau = min(TE, k0/M2, omega*k0/M3) and asserts the modal lower bound 1 + k0 t^2/4; the program passed.

**pass3: receiver exponents at a = 1/4** (confirmed)

verify_rotational_receiver.py passed all 15 checks including 'RMS velocity scale match implies Reynolds scale match', 'weight mass and RMS scaling give exact energy ratio' and the scale exponents; REPORT.md records that at a = 1/4 these are exactly the candidate velocity, Reynolds and energy exponents.

**pass4: p_zz'(0) + 32 approximately -241.35 (split) and -241.42 (source) at c=2, nu=0.01** (confirmed)

recomputed at nr=1600, nmu=1536, lmax=1280: T_split = -271.61999921482595, T_source = -271.6869334279583, viscous coefficient 3027.118162213184, so T_0.01 split = -241.3488175926941 and T_0.01 source = -241.41575180582646. The note records -241.348818 and -241.415752.

**pass4: resolution change about 0.00105 in the split representation** (confirmed)

split T_0 at 1200 nodes = -271.6210449882172, at 1600 nodes = -271.61999921482595; the difference is 0.0010457733912403455.

**pass4: the two representations differ by about 0.06693 at the finer resolution** (confirmed)

representation_difference = 0.06693421313235604 in the recomputed 1600 run.

**pass4: C_P about -0.42697138636, C_v about 0.00000627156153, A about 1096.21507043, dP about -170.89870411, dv about 0.00308791848** (confirmed)

recomputed 1600 run: Cp = -0.42697138636192755, Cv = 6.2715615295074624e-06, A at c=2 = 1096.215070427095, dP = -170.89870411328232, dv = 0.0030879184800330544.

**pass4: core controls t300 about -6.73446331117 and tWb about 2.563122131694** (confirmed)

evaluate_core_pressure_radial.py: t300 = -6.734463311167136 and tWb = 2.563122131694385 at n=16001, with source_vs_reduced about -3.8e-11.

**pass5: the full initial gate certificate** (confirmed)

certify_full_initial_gate.py --recompute regenerated all four integral enclosures and reproduced every field of FULL_INITIAL_GATE_CERTIFICATE.json, including every exact dyadic mantissa and exponent: status INITIAL_FEEDBACK_CERTIFIED_USING_ANALYTICAL_REDUCTIONS_AND_MPMATH_IV, full_gate_upper_exact = -290649/8750 = -33.21702857142857, beta_second_derivative_lower_exact = 290649/17500 = 16.608514285714286, reservoir_log_derivative_lower_exact = 226249/51000, coarse bounds C300_Phi < 66, tWb < 4, dv/Cv < 901, core kernel < -12/5, pump kernel < -31/2. verify_full_gate_bound.py independently reprinted the same two exact rationals.

**pass7: the D certificate** (confirmed)

certify_D_viscosity.py --panels 1024 reproduced D-viscosity-certificate.json field for field: status PASS, radial second-derivative energy interval [100952.1793303089676983218159257833218569, 113430.7471049831419909679864742410035166] with radial_upper_below_120000 true, viscous upper bound interval [3.263029119503281114289078677167430340491, 3.658516517079472127305937188046058041206] with viscous_upper_below_four true, sufficient threshold 5093339/210000.

**pass7: the E certificate** (confirmed)

certify_outer_E.py --panels 2048 reproduced outer-E-certificate.json field for field: status PASS, E_upper_strictly_below_one_sixtieth true, E_upper_expression [0.01308016428837384303488356186423805542445, 0.01661614072972984853428293641397043671335], compatibility_upper_expression [17.34176757414131989558123310152281651597, 20.94846354432444550496859514224984546097] against threshold 5093339/210000. verify_E_certificate_audit.py and verify_common_family.py then replayed the exact dyadic endpoints to T_lambda < -435297/70000 and beta_second > 435297/140000 for every lambda in [3/4, 1].

**pass8: the fixed-circle acceleration certificate** (confirmed)

certify_fixed_circle_acceleration.py --panels 2048 reproduced fixed-circle-acceleration-certificate.json: status PASS, unique_global_radial_maximum_interval ['1077377876667/5000000000000', '269344469167/1250000000000'], Gtt_upper_expression [-304786.594933592446230202556004150034858572597, -293333.188604612339424175039523497704572402962] with Gtt_upper_strictly_below_minus_250000 true, moving_radial_maximum_Gtt_upper_below_minus_240000 true, Gt_interval [85.6465125470615141262294034129150102935027932, 171.279549323398698324019039360862838176865749] with initial_Gt_between_85_and_186 true, radial_speed_above_one_fifth true.

**pass8: the leading-envelope and leading-pressure certificates** (confirmed)

certify_leading_envelope.py --panels 1024 reproduced R1 interval [154278.6422641070164361793701028804422644, 175257.3100341986565244852148553171403894], radial_extraction [0.02196178344685765685795724681084006528738, 0.02438365289314620584089211100546093441205], mean_derivative_lower 8871/800 and fluctuation_energy_fractional_derivative_lower 12853/2205. certify_leading_pressure.py --panels 1024 reproduced E_absolute_upper_expression [0.2896117095151673399599848780222081378255, 0.3036125064520005326017319565070058626543] with E_absolute_upper_strictly_below_one_third true, nu_d_rational_upper 9078400/1419857 and common_T_negative_margin 12493177/1120000.

**pass9: the affine-envelope budget** (confirmed)

verify_affine_envelope_budget.py --panels 2048 reproduced affine-envelope-budget-certificate.json exactly: status PASS, envelope_second_derivative_moment interval [480868.53632052795729, 641080.35446054721251] over 314 panels, A2_coarse_upper 800000, inviscid_additive_gain_lower 16, viscosity_loss_upper 2062985492/1406514375, E0_upper 366140/3969, net_additive_gain_lower 27927869/1984500, relative_gain_lower 3/20, time_upper 7/5000.

**pass9: the general-envelope initial-jet certificate and the fixed-amplitude initialization** (confirmed)

certify_general_envelope_initial_jet.py --panels 2048 reproduced all 227 fields; verify_general_envelope_initial_jet_certificate.py replayed Gt > 11, actual Gtt < -13000, tracked Gtt < -12000 and pressure error < 1.3. verify_fixed_amplitude_initialization.py reproduced X_lower 293913/50000, beta_prime_lower 17391/700000, complete_pressure_T_upper -634112687/50000000 and beta_second_lower 4264878809/700000000 at A = 1020, lam = 1/4.

**pass9: slice solver validation** (confirmed)

check_slice_solver.py reproduced slice-solver-checks.json exactly: discrete divergence below 1.8e-15, nonlinear circulation error below 1.8e-14, and Poisson relative max errors 1.316768e-04, 3.290126e-05, 8.225802e-06 and Gaussian relative max errors 1.712974e-06, 4.286274e-07, 1.071808e-07 across nr = 192, 384, 768, that is refinement ratios of 4.00 in both, inside the 3.5 to 4.5 window the script asserts. evolve_circular_slice.py at defaults reproduced every observable of slice-768.json bit for bit at all 21 shared sample times.

## What this replication does not say

Nothing here is a mathematical judgement. Every VERIFICATION.md in the tree already states that its programs check
finite algebra, interval enclosures and exploratory quadratures, and that no Navier-Stokes time integration, formal
proof-kernel replay, useful stage duration, inherited return or blowup is certified. Re-running those programs on a
second host confirms that they run and that they produce the recorded numbers. It does not enlarge a single claim.


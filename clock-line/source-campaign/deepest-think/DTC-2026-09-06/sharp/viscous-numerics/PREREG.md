# PREREG — sharp seat `viscous-numerics` — DTC-2026-09-06
### written BEFORE any production run (L-09: deviations logged, thresholds never re-cut)

Tier: receipts, single seat (L-08). Nothing here is promoted. Numerics falsify, never prove.

## 0. Question
Does the ACTUAL viscous axisymmetric no-swirl NS evolution of a dyadic N-ring datum realise a
doubling time `t_d <= C/(M log Re_E)` — i.e. is Astra's continuation clock
`H0 = c/(M(1+log_+ Re))` SHARP in the log — or does something (viscous decay of the inner
rings, deformation of the stack, the O(M) corrections) prevent it?

## 1. System solved (no model, no closure)
Lift variables: `eta = omega^theta / r`, `L5 = d_rr + (3/r) d_r + d_zz`, `L5 psi1 = -eta`,
`Psi = r^2 psi1`, `u^r = -(1/r) d_z Psi = r a` with `a = -d_z psi1`, `u^z = (1/r) d_r Psi`,
and the EXACT no-swirl axisymmetric NS transport `D_t eta = nu L5 eta`.
Inviscidly `eta` is a material invariant, so `omega^theta = r eta` grows EXACTLY like the
material radius: `d ln(omega^theta)/dt = a` along parcels. Viscosity is carried, not modelled.

## 2. Data (fixed now)
Units `rho0 = 1`, `M = 1`; every time below is therefore `t*M`.

**Datum A (primary, "smoothed bang-bang shell over N octaves")**
`eta_0 = -A * (2z/s^2) * Theta(s)`, `s = |x| = sqrt(r^2+z^2)`,
`Theta(s) = (1/2)[tanh((s-rho0)/w0) - tanh((s-R)/w1)]`, `w0 = 0.25 rho0`, `w1 = 0.10 R`,
`R = 2^N rho0`. Then `omega^theta = -A sin(2 phi)`, so `|omega^theta| <= A`, vanishing on the
axis (admissible smooth axisymmetric data), maximal on the 45-degree cones at EVERY octave.
`A` is rescaled so the DISCRETE `sup|omega^theta|` at `t = 0` equals exactly `M = 1`.
This is the brief's "bang-bang field, smoothed", made axis-regular; `eta_0` is odd in `z`,
which the dynamics preserves exactly, so the run is done on `z >= 0` with `psi1(z=0) = 0`.

**Datum B (variant, "N discrete rings")** `N` Gaussian ring-pairs centred at
`(r_k, z_k) = s_k(sin 45, cos 45)`, `s_k = rho0 2^k`, `k = 0..N-1`, widths `sigma_k = 0.20 s_k`,
amplitudes `C_k = sqrt(2) M / s_k` so each ring alone has `sup|omega^theta| = M`; same rescale.

`nu = M rho0^2 / Re0` with **Re0 = 100** (primary). Sensitivity run at `Re0 = 400`.
`Re_E = M l^2 / nu`, `l = E^{1/5} M^{-2/5}`, `E = ||u||_2^2` (no 1/2) — E computed
numerically from the t=0 velocity field over the FULL space, never assumed.

`N in {1,2,3,4,5,6}` primary, `7` if it fits the budget.
Grid: uniform, `h = rho0/8` primary; box `rmax = zmax = max(1.5 R, 4 rho0)`.
Convergence: `h = rho0/6` and `rho0/12` at `N <= 4`, `rho0/6` at the largest N run.
Domain check: `lambda = 2.25` at `N = 3`.
Horizon `T_max = max(0.6, 4/N) / M`; the run stops early once `2M` is reached.

## 3. Measured quantities
`T32(N)` = first `t` with `sup|omega^theta| >= 1.5 M`;  `T2(N)` = first `t` with `>= 2 M`
(linear interpolation between recorded steps). Also: `a` at the argmax, the spherical radius
`s*` of the argmax (reported as octave index `log2(s*/rho0)`), `sup|eta|` (must be
NON-INCREASING — exact maximum principle for `D_t eta = nu L5 eta`), and `Re_E`.

## 4. PRE-REGISTERED GATES
Define `P(N) = T32(N) * M * log Re_E(N)` and `Q(N) = T32(N) * M * log(R/rho0)`.
Fit `log(T32(N)*M) = alpha - beta * log N` by least squares over every `N >= 2` that reached
`1.5 M` within `T_max`.

- **SURVIVE** iff `beta >= 0.7` AND `P(N) <= 8` for every N run.
- **KILL** iff `beta <= 0.3` (i.e. `T32*M` does not fall as octaves are added), or if the
  largest N run fails to reach `1.5 M` within `T_max` while `N = 2` reaches it.
- **INCONCLUSIVE** iff `0.3 < beta < 0.7`.
The first-order prediction is `beta = 1`, `Q -> ln(1.5)/kappa` with `kappa` the measured
strain-per-e-fold constant of the datum.

## 5. Controls named in advance (L-88) — each must fire
- **C1 sign control.** The sign-reversed datum (`A -> -A`) has `a < 0` at the argmax; its
  `sup|omega^theta|` must DECREASE monotonically and never reach `1.5 M`. If both signs grow,
  the code is wrong and every number here is void.
- **C2 operator control.** Replacing `3/r` by `1/r` in `L5` (the plain 3D Laplacian) must
  break the exact-Hill validation by >100%.
- **C3 exact-Hill validation.** Full-domain solve on Hill's spherical vortex must return
  `psi1` to a few 1e-3, `a_nose -> 1/5`, axis peak `/U -> 5/2`, Eulerian toll `-> 5/4`.
- **C4 5D diffusion law.** Pure diffusion (u = 0) must give `d<|X|^2>/dt = 10 nu`.
- **C5 Lagrangian identity.** Tracers must satisfy `omega^theta(t)/omega^theta(0) = r(t)/r(0)`
  to the inviscid accuracy of the run.
- **C6 single-octave control.** `N = 1` must show NO log gain: its measured `a` at the inner
  cone must equal `kappa M ln 2 + O(M)`, and the fitted `a(N)` must be linear in `N`.

## 6. Discipline
Every number below comes from a script in this folder that I wrote and ran; nothing is typed
from memory. Exponents/constants claimed exact are derived by hand and only CHECKED
numerically. Grid convergence is reported for every headline number. If a gate cannot resolve
what it was set to resolve, that is logged as a deviation and the gate still stands as written.

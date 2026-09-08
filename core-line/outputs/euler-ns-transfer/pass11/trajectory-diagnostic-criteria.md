# A bounded endpoint diagnostic for the selected full field

Target one uninterrupted evolution of the fixed pass9 datum
`u0=M+1020 v+(1/4) w_L`, with `nu=.001`, through `T=.001`. Keep the
amplitude, profiles, phases and all generated modes fixed by that initial
value problem. This note specifies useful **finite-cylinder numerical
evidence** and the later whole-space validation gates. It reports no new
trajectory or error enclosure.

## 1. Measure retained endpoint gains directly

Use the same physical coordinates and initial datum in every comparison.
Record these quantities at common times, including exactly zero and T:

| Quantity | Useful endpoint comparison | Qualification |
|---|---|---|
| Mean angular momentum `G=r mean(u_theta)` | A selected receiver value above the exact initial global mean maximum `M0` | A single receiver suffices to lower-bound the later global maximum. |
| Nonaxisymmetric energy `K=||u−mean(u)||2²/2` | `sqrt(2K(T))−sqrt(2K0)>0` | This amplitude margin converts directly to an L2 error budget. |
| Core strain and rotation | `b(T)−Omega(T)>0`, plus `Omega(T)>0` | Prefer the linear cone margin to a quotient alone. Also report `Omega(T)−1` if rotational growth is desired. |
| Specified full-velocity weighted RMS `V_L` | `L V_L(T)/nu` exceeds the chosen receiving Reynolds target | Mean circulation or K does not substitute for full receiving velocity. |

Report each endpoint's absolute value and its change from that run's
represented initial field. Also compare with the independently evaluated
analytic initial references. Subtracting a numerically biased initial
value is a useful change diagnostic; it does not remove initialization
error from a claim about the selected exact datum.

The current `observe` routine supplies G, K, b and Omega, but **does not
supply full receiving RMS, mode-resolved data, or an inherited-state
test**. Its `mean_G_branch` selects the largest spline critical value in
`.18<r<.3` at **z=4** independently at each observation. This is a
radial-window stationary-point selector, not necessarily a continuously
tracked branch or the global maximum in `(r,z)`. It does not compare the
window endpoints and falls back to the old radius if there is no selected
critical point. Store the receiver radius and value; a
search over a declared `(r,z)` region can locate an additional receiver.
An initial axial plateau makes a unique global maximizing point
unavailable without further analysis. A changed maximizer is permissible
for an endpoint lower bound, provided the location and interpolation error
are retained.

No endpoint gate requires `G_t`, `K_t`, or `beta''` to remain positive at
every intermediate time. A pulse can still leave a useful positive gain.
Conversely, initial signs or a quadratic jet cannot establish gains at T.
The selected outer rotation has angular phase scale `m A T=4.08`, so
short Taylor polynomials are particularly unsuitable here.

## 2. Reject unresolved initialization before interpreting core dynamics

The odd-grid initial reference at `nr=192` gives:

| nz | K(0) | Discrete Omega'(0) |
|---|---:|---:|
| 1025 | 6.4563975452 | −4.7541692 |
| 2049 | 6.4563975452 | 1.9968441 |

The actual insertion identity is `Omega'(0)=2`. Thus converged energy,
small discrete divergence and a stable initial Omega do not establish
resolved core evolution. The 2049 result is a promising calibration point,
not proof that later axial derivatives are resolved.

Before a long refinement, compare initial `b,Omega`, the full generator's
`b',Omega'`, initial G and K, initial nonaxisymmetric energy work, and the
sampling/projection correction. Use `b=Omega=1`, `Omega'=2` and the
independent pressure references. For b', compare to the chosen finite
domain's pressure response separately from the whole-space value. Check
core extraction with independent parity-aware radial fits and local axial
derivative stencils; a single four-cell extrapolation is not an error
estimate. The current nr192 initialization correction is about 1.80 in
discrete L2, which cannot be treated as zero merely because projection
removes discrete divergence.

The older even-grid/filter pilot histories have corrupted initial core
rotation and should remain implementation diagnostics. They do not count
as a convergence sequence toward an accepted core endpoint.

## 3. Minimum comparisons for a credible finite-cylinder endpoint

Start from one baseline with resolved initial jets. Refine **one source
of error at a time**, re-evolving the same analytic datum each time:

1. **Time:** compare `dt` with `dt/2` at the same grid and T. The symmetric
   CN/RK4/CN composition is generally second order in time, despite its
   fourth-order nonlinear substep. If the observed difference consumes
   the gain margin, add the next halving or reject that endpoint. A
   two-grid difference is not a rigorous remainder bound.
2. **Radial:** compare 192 with 384 cells, or the next comparable
   refinement, at fixed axial/angular resolution and time step. Check
   the core separately from the outer packet joins and RMS receiver.
3. **Axial:** compare 2049 with 4097 points at fixed period, initially
   using the cheap initial-jet gate. An endpoint axial refinement is
   still required before accepting axial convergence. Odd grids remove
   the even Nyquist ambiguity; they do not remove Fourier truncation.
4. **Angular:** compare `J=4` with at least `J=6` (maximum physical modes
   16 and 24), preserving dealiasing. Check low-mode changes as well as
   the newly resolved modes. If terminal-mode content is significant or
   the endpoint changes appreciably, increase J further.

The same fine baseline can serve all four comparisons. If discrepancies
interact or are close to a target margin, require a jointly refined run;
do not infer a rate from one diagonal refinement. A practical screening
rule is that the sum of the observed endpoint discrepancies and
initial-reference offsets be below one quarter of every claimed gain
margin. This is a conservative **numerical prioritization rule**, not a
certified error bar; cancellation or unresolved tails can defeat it.

Use energy balance, dissipation, discrete divergence, reality and C4/odd
symmetry as integrity checks throughout. In particular compare the
energy-balance defect to both the viscous loss and the much smaller K
change, rather than only to the large total energy. Excellent total-energy
balance can coexist with an incorrect core derivative, as the old pilots
demonstrate.

## 4. Energies alone miss the relevant phase and transport

Save complex Fourier coefficients of the complete endpoint, its initial
representation and enough intermediate snapshots to inspect a pulse.
For genuine angular coefficients `a_j`, report

`K_j = 2 pi integral |a_j|² r dr dz`, `j>0`,

with the corresponding MAC mass quadrature. Then `K=sum K_j`. Track
energy and differentiated energy near the highest retained angular and
axial modes, plus radial join-resolution diagnostics. A small final mode
is helpful but cannot certify a small omitted tail.

Compare complex coefficient fields after evaluation on a common physical
grid, respecting face versus cell locations, Fourier normalization and
the axial origin phase. Include absolute L2 differences for the full
field and separately for the nonaxisymmetric field; division only by the
large pump's norm would hide a poor seed approximation. Add local C1/core
and receiver discrepancies. For a future H4 validation target, compare
ordered Cartesian derivatives through order four from a smooth
axis-regular reconstruction; raw cylindrical component differences
without the geometric terms are not that norm.

Do not independently align each mode's phase or refit the packet before
comparison. Those operations can conceal differential phase and alter
quadratic stress. A separately reported global-rotation fit is a useful
symmetry diagnostic, but retain the unaligned same-coordinate difference.
As a mechanism diagnostic also save the actual covariances

`R=2 Re sum_{j>0} a_{j,r} conjugate(a_{j,theta})`,
`Z=2 Re sum_{j>0} a_{j,z} conjugate(a_{j,theta})`,

and the mean torque `−2R−r R_r−r Z_z` near the receiver. These distinguish
stored seed energy from useful transfer; they retain generated vertical
motion and do not impose frozen polarization.

## 5. Domain error and eventual certification remain separate

At period 16 the pass10 initial axial-image contribution to p_zz is about
`+0.48013`, versus about `+0.008383` at period 32. Its instantaneous b'
effect is about `−0.24007`, comparable over a short interval with the
small desired core-cone gain. It is not legitimate to multiply that
initial correction by T and call it the actual later domain error.
After the fixed-domain refinements, enlarge the axial period and radial
wall radius while preserving local physical resolution. Periodic copies
and the no-slip wall are distinct effects. No finite number of such
comparisons provides a whole-space error enclosure.

For a later smooth divergence-free whole-space approximation v, the
pass11 hierarchy must receive certified full-slab coefficients, residuals
and inherited initial errors. If its endpoint radius is `rho_ell` and
its zeroth component is `z0`, the pass10 conversions give
`epsilon_infty=ell^(−3/2)rho_ell`,
`epsilon_grad=ell^(−5/2)rho_ell`, and `||u−v||2<=z0`.
The usable budgets are then explicit:

| Desired gain | Sufficient positive margin after errors |
|---|---|
| Mean gain eta_G | `G(v,T,r,z)−(1+eta_G)M0_up−r epsilon_infty>0` |
| Seed-energy gain eta_K | `sqrt(2K(v,T))−sqrt(2(1+eta_K)K0_up)−z0>0` |
| Core ratio at least 1+eta_beta | `b_v−(1+eta_beta)Omega_v−(2+eta_beta)epsilon_grad>0` and `Omega_v>epsilon_grad` |
| Rotation above Omega_min | `Omega_v−Omega_min−epsilon_grad>0` |
| Full weighted RMS | `V_L(v)−sqrt(||chi_L||infty/M_L) z0` exceeds the specified target |

Coordinate uncertainty, initial-reference enclosure width, reconstruction
jumps and exterior continuation also consume these budgets. The core
margin is likely the tightest: a change of order `10^−4` requires core
gradient accuracy appreciably better than `10^−4`, irrespective of total
energy accuracy. This is a scale warning, not a predicted endpoint value.

A validated positive endpoint would establish a retained finite gain.
A successor result additionally needs a named full-field class with a
proved stable stage lemma and the actual normalized inherited endpoint
inside it. The present observers do not test this condition. The evolved
outer field, core deformation, phases and tails must be retained; the
initial pressure-jet identities cannot be reapplied to a freshly reset
template.

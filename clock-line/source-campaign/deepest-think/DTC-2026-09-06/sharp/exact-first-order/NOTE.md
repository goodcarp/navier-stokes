# exact-first-order — N-ring datum, exact strain, energy, c1, second order
Deepest-Think-Claude sub-seat, DTC-2026-09-06. All numbers below come from scripts in this folder that
I wrote and ran; every script is listed in `SHA256SUMS`. Numerics falsify, they never prove; where a
statement is an identity it is marked EXACT and derived, not fitted.

---

## 0. The instrument (`ns5d.py`, calibrated in `c1_calibrate.py`)

Re-derived from scratch, not quoted. For axisymmetric no-swirl flow, with `eta = omega^theta / r` and
`Psi = r^2 psi`:

    Delta_5 psi = -eta ,  Delta_5 = d_rr + (3/r) d_r + d_zz   (Laplacian on R^5 = R^4_r x R_z)
    a := u^r/r = -d_z psi ,   u^z = 2 psi + r d_r psi ,   G_5 = 1/(8 pi^2 |x|^3)

so with `d^4x' = r'^3 dr' sin^2(th) dth dOmega_2`:

    a(r,z) = (3/(2 pi)) int int eta(r',z') r'^3 (z-z') J5 dr'dz' ,  J5 = int_0^pi sin^2 th (A-B cos th)^{-5/2} dth
    psi(r,z) = (1/(2 pi)) int int eta(r',z') r'^3      J3 dr'dz' ,  J3 = same with exponent -3/2
    A = r^2+r'^2+(z-z')^2 ,  B = 2 r r' .

**EXACT closed forms for the theta integrals** (derived in the header of `ns5d.py` via u = 1-2w, the Beta
integral, and a Pfaff transformation), in terms of the classical elliptic modulus
`k2 = 4 r r' / ((r+r')^2 + (z-z')^2)`:

    J5 = (pi/2) 2F1(1/2, 3/2; 3; k2) / [ ((r-r')^2+(z-z')^2) * ((r+r')^2+(z-z')^2)^{3/2} ]
    J3 = (pi/2) 2F1(3/2, 3/2; 3; k2) /                          ((r+r')^2+(z-z')^2)^{3/2}

This form is the reason the whole computation is clean: `2F1(1/2,3/2;3;.)` is **finite** at k2 = 1
(c-a-b = 1 > 0; value 16/(3 pi)), so the entire singularity of the strain kernel is the explicit factor
`1/((r-r')^2+(z-z')^2)`, and in polar coordinates about the evaluation point the integrand is bounded.
`2F1(3/2,3/2;3;.)` has c-a-b = 0, i.e. exactly the log divergence of the ring self-energy; scipy returns
`inf` there, so `ns5d._F3` splices in the exact leading behaviour
`(8/pi)[-log(1-x) + 4 log 2 - 4]` below 1-x = 1e-9 (verified against mpmath to 0 ulp at 1-x = 1e-16).

### Calibration (`c1_log.txt`)
| check | result |
|---|---|
| J5, J3 closed forms vs mpmath, 40 random points | worst rel 8.9e-15 / 1.3e-14 |
| J5 near coincidence, t = 1e-1 … 1e-8 | `J5 * 3 r^3 t^2 -> 1.000000` (exact flat-space codim-2 limit) |
| **Hill's vortex `a = A z/5`**, 3 off-axis + 2 on-axis interior points | rel 7.5e-15 … 1.8e-11 |
| Hill potential `psi = A(R^2/6 - \|x\|^2/10)`, 3 interior points | rel 1.2e-9 … 1.6e-9 |
| `E_Hill = 16 pi A^2 R^7/315`, i.e. `KE = (10/7) pi U^2 R^3` | exact (symbolic) |
| bang-bang shell axis identity `a(0,0) = (M/2)log(R/rho0)`, R/rho0 = 4, 4096, 2^20 | rel < 2.2e-15 |
| sign control (omega -> -omega) | exact negation |

`c1_calibrate.py` exits 1 on exactly one row: the t = 1e-8 near-coincidence comparison lands at
1.22e-8 relative against my hard-coded 1e-10 threshold. There the **mpmath reference** is the weaker
side (mp.quad on a peak of width 1e-8), and the same row's exact-limit ratio is 1.000000. Every
physical calibration passes. I left the threshold as written rather than loosening it after the fact.

---

## 1. The datum (explicit)

Fix `M > 0`, `rho0 > 0`, aspect `s in (0, 1/6]`, polar angle `phi*`. For `k = 0..N-1`:

    rho_k   = rho0 2^k
    centre  (r_k, z_k) = rho_k (sin phi*, cos phi*)
    sigma_k = s rho_k
    eta_k(r,z) = -(M/r_k) exp( -((r-r_k)^2 + (z-z_k)^2) / (2 sigma_k^2) )
    omega^theta = r * sum_k eta_k        (so omega^theta -> 0 at the axis; eta is the transported scalar)

All cores in `z > 0` with `omega^theta < 0` — the bang-bang sign that the estate's extremal identity
`sup a(axis) = (M/2) log(R/rho0)` selects (there `omega^theta = -M sgn z`).

**phi\* is not a choice.** The exact reduction of the axis strain is
`a(0,0) = (3/4) int int omega^theta (-sin^2 phi cos phi) d log rho dphi`, so the per-e-fold weight is
`sin^2 phi |cos phi|`, maximised at `phi* = arccos(1/sqrt 3)` with value `2/(3 sqrt 3) = 0.3849001795`
(numeric argmax 0.95531313 vs exact 0.95531662).

**Normalisation.** `||omega_0||_inf` is computed, never assumed, and every reported quantity is divided
by it. `Z := ||omega_0||_inf / M` = 1.0011978, 1.0026892, 1.0047660, 1.0074200, 1.0117387 for
s = 0.04, 0.06, 0.08, 0.10, 0.125 (the excess is the `r`-factor in `omega = r eta`, `~ 1 + s^2/(2 sin^2 phi*)`,
plus neighbour tails).

**Smoothness caveat, stated exactly:** each ring's Gaussian has amplitude `exp(-r_k^2/(2 sigma_k^2)) =
exp(-0.3333/s^2)` on the axis (2.4e-4 at s = 0.125, 1e-91 at s = 0.04), so the datum is within that of a
genuinely smooth axisymmetric field; mollifying changes nothing reported here at the quoted digits for
s <= 0.125. I kept s <= 0.125 throughout so that `r_k - 6 sigma_k > 0`.

---

## 2. Exact strain at the innermost core

**Scale invariance.** `a` is invariant under `x -> lambda x` at fixed M, so the contribution of the ring
at `rho_0 2^m` to `a` at the ring-0 core depends only on `(m, s, phi*)`. Call it `M A_m(s)`. Then exactly

    a_j(N)/M = sum_{m=-j}^{N-1-j} A_m  ,      a_inner(N)/M = a_0(N)/M = sum_{m=0}^{N-1} A_m .

Quadrature convergence (`c2_log.txt`): triples (nt,nal) = (20,256),(32,512),(44,768) agree to
relative 1.5e-16 … 5.5e-16 for m = 1, 3, 10.

### THEOREM A (EXACT, not numerics) — the self-strain of a ring vanishes identically
`J5` depends on `(z-z')` only through `(z-z')^2`, and the strain integrand carries the odd factor
`(z_0 - z')`. If a core is symmetric in z about its own centre, the integrand is **odd** under
`z' -> 2 z_k - z'` and the integral is **exactly zero**. Hence `A_0 = 0` for every s and every
`phi*`. Measured: `|A_0| <= 8e-19` (machine zero) at all five s. So the N = 1 control gives not merely
`O(M)` and no log — it gives exactly 0.

### THEOREM B (EXACT to O(s^6), `c7_kappa_analytic.py`) — the per-octave increment
For `m -> infinity` the evaluation point is at the origin relative to the source ring, where
`J5 = (pi/2)/rho'^5` exactly, so with `f := r'^3 z'/rho'^5` and the Gaussian moments
`int G = 2 pi sigma^2`, `int G d_i d_j = 2 pi sigma^4 delta_ij`:

    kappa(s)/M = (3 pi/2)(sigma^2/r_k)[ f + (sigma^2/2) Laplacian f ]_{(r_k,z_k)} + O(s^6) .

Sympy gives, **at phi\* and only there**, `Laplacian f / f = -6` exactly, hence

    kappa(s) = (pi / sqrt 3) s^2 (1 - 3 s^2) + O(s^6) .

Independently confirmed by the single-ring far-field check in `c1` (ratios 0.998800, 0.992512, 0.970187
at s = 0.02, 0.05, 0.10 — i.e. 1-3s^2 = 0.9988, 0.992500, 0.9700) and by the c2 quadrature (rel 4.8e-6
at s = 0.04 rising to 4.8e-4 at s = 0.125, consistent with the O(s^6) truncation).

### The answer to "is `a_inner = (M/2)(N ln 2) + O(M)`?" — **NO.**
`A_m/kappa` reaches 1 geometrically (s = 0.10: m = 1,2,3,4,6,8 give 1.3583, 1.1470, 1.0736, 1.0375,
1.0096, 1.0024), and `A_m < 0` for m < 0 (inner rings **compress** outer ones; `A_-1/kappa = -0.0957`,
decaying by ~1/16 per step). So `a_inner(N)/M` is **exactly affine in N**:

| s | Z | kappa (normalised) | b | kappa / ((1/2)ln2) | 5.2340 s^2(1-3s^2) | c1_asym = 2(ln2)^2/kappa |
|---|---|---|---|---|---|---|
| 0.040 | 1.0011978 | 2.884710e-03 | -8.2156e-04 | 0.008324 | 0.008333 | 333.10 |
| 0.060 | 1.0026892 | 6.441998e-03 | -1.9273e-03 | 0.018588 | 0.018637 | 149.16 |
| 0.080 | 1.0047660 | 1.133233e-02 | -3.6191e-03 | 0.032698 | 0.032851 |  84.79 |
| 0.100 | 1.0074200 | 1.746766e-02 | -6.0340e-03 | 0.050401 | 0.050765 |  55.01 |
| 0.125 | 1.0117387 | 2.671157e-02 | -1.0327e-02 | 0.077073 | 0.077941 |  35.97 |

`b/kappa = -0.345` at s = 0.10, i.e. `a_inner = M kappa (N - 0.345)`; the -1 from `A_0 = 0` is partly
repaid by the near-neighbour excess (`sum_m (A_m/kappa - 1) = -1 + 0.652 = -0.348`, matching the fit).

**The rate is `(M/2)ln2` per octave times the VOLUME FRACTION, not per octave.**

    kappa(s) / ((M/2) ln 2) = (2 pi / (sqrt3 ln 2)) s^2 (1 - 3 s^2) = 5.2340 s^2 (1 - 3 s^2),

which equals 1 only at s ~ 0.437 — where `sigma_k = 0.44 rho_k` and the "rings" have merged into the
bang-bang shell (at s = 0.437 each ring's Gaussian already has amplitude `exp(-1/(2s^2)) = 0.073 M` at
its neighbours' centres). A genuine, well-separated ring stack loses a factor `5.234 s^2`: at s = 0.10
it achieves **5.0%** of the extremal per-octave rate.

Inviscid first-order growth of the sup: `d/dt sup|omega| |_{t=0} = a_inner * M` (all cores carry the same
peak at t = 0, and `a_j` is decreasing in j, so the innermost core is the maximiser).

### Controls (`c2_log.txt`, section 4)
* **(a) N = 1:** `a/M = 1.95e-19` — exactly zero by Theorem A, no N and no log. PASS (stronger than asked).
* **(b) sign-reversed stack (N = 8, s = 0.10):** bang-bang `+1.3461499788e-01`, reversed
  `-1.3461499788e-01`, sum `0.00e+00` (exact negation). `a_inner > 0` for bang-bang, `< 0` reversed. PASS.
* **(c)** `A_m < 0` for m < 0 with ratio 0.0512, 0.0580, 0.0602 -> 1/16, as the far-field
  `rho_j^4/rho_k^4` dipole scaling requires. PASS.

---

## 3. Energy, Reynolds number, and c1

**EXACT energy identity** (derived, Hill-verified in `c1`):

    E := ||u||_2^2 (no 1/2) = 2 pi int int omega^theta Psi dr dz = 2 pi int int eta psi r^3 dr dz
       = (1/pi) int_{R^5} eta psi = (1/pi) || grad psi ||^2_{L^2(R^5)} >= 0 .

For the ring stack (`c3_energy.py`), scale invariance of the bilinear form gives
`E_{k,l} = M^2 rho_min(k,l)^5 h_{|k-l|}(s)`, hence

    E(N) = M^2 rho0^5 [ sum_{k<N} 32^k h_0 + 2 sum_k sum_{j=1}^{N-1-k} 32^k h_j ] .

`h_j` grows like ~2^j (the interaction energy `E_{0,j} ~ M^2 s^4 rho_0^4 rho_j`), but the pair (k, k+j)
contributes `32^k 2^j`, so the (k = N-1-j, j) family falls as `32^{N-1} 2^{-4j}`: the j-sum converges
geometrically and truncating at j = 6 costs < 2^-28. **E is dominated by the outermost ring's own
self-energy** — as it must be, since `E ~ M^2 R^5`.

**Two viscous-floor conventions, both reported.** The floor must sit on the *smallest* length in the
datum, which for a ring stack is the core, not the innermost radius:
* **(I) core-limited (physical):** `sigma_0 = sqrt(nu/M)`, i.e. `nu = M s^2 rho0^2`.
* **(II) the brief's:** `rho_0 = sqrt(nu/M)`, i.e. `nu = M rho0^2`.

`Re_E = E^{2/5} M^{1/5}/nu`, `t_d = ln2 / a_inner`, `c1 := t_d * M * log Re_E`.

**Asymptotics, exact.** `E(N) = C(s) M^2 rho0^5 32^{N-1}(1+o(1))` so `log Re_E = 2(N-1)ln2 + O(1)`,
while `a_inner = M kappa(N + b/kappa)`. Hence

    c1(N) -> 2 (ln 2)^2 / kappa(s)     as N -> infinity ,

with the O(1/N) correction fixed by the additive constants. Values in the table of §2.

Quadrature convergence of the energy at s = 0.10: `h_0` = 5.9571887570e-03, 5.9571890770e-03,
5.9571890772e-03 for outer Gauss-Hermite n = 12, 18, 24, and 5.9571890770e-03 with a finer inner polar
grid (nt 26->38, nal 256->384) — 10 digits stable. `h_0/h_thin-ring-uniform-core` = 0.809, 0.788,
0.772, 0.761, 0.752 for s = 0.04..0.125 (a Gaussian core is lighter than a uniform one, as expected).

### `c1(N, s)` for the ring stack (`c3_log.txt`), viscous floor (I) `sigma_0 = sqrt(nu/M)`

| N | s=0.04 | s=0.06 | s=0.08 | s=0.10 | s=0.125 |
|---|---|---|---|---|---|
| 4  | 480.76 | 209.49 | 116.95 | 74.95 | 48.57 |
| 8  | 399.42 | 176.10 |  99.08 | 63.83 | 41.52 |
| 12 | 376.61 | 166.82 |  94.15 | 60.78 | 39.59 |
| 20 | 358.95 | 159.65 |  90.34 | 58.43 | 38.12 |
| 26 | 352.92 | 157.20 |  89.05 | 57.63 | 37.61 |
| **N -> inf** | **333.10** | **149.16** | **84.79** | **55.01** | **35.97** |

`c1` descends monotonically to `2(ln2)^2/kappa(s)` from above, as the additive constants predict.
Convention (II) (`nu = M rho0^2`) shifts `log Re_E` by the constant `-2 log(1/s)` and so gives the same
limit with a larger O(1/N) correction (e.g. s = 0.10: 50.51 at N = 26); at small N it makes `Re_E < 1`,
which is the sign that (II) is the wrong floor for a ring stack — the core, not `rho0`, is the smallest
length.

At N = 26, s = 0.125 the datum has `E = 6.23e35 M^2 rho0^5`, `Re_E = e^{37.1}`, `a_inner = 0.684 M`,
`t_d = 1.013/M`.


### The extremal endpoint, computed exactly (`c6_shell_exact.py`)

The `s -> space-filling` end of the family is the bang-bang shell itself, and there everything is exact.
Independent route (no singular quadrature anywhere): 5D zonal-harmonic separation in Gegenbauer
`C_l^{3/2}`, `H_l` exact rationals, C^1 matching solved in adaptive precision
(`dps = 40 + (2l+3) K log10 2 + 20`), radial integrals in closed form, `l` up to 201.

    eta = -M sgn(t)/(rho sqrt(1-t^2)) = (1/rho) sum_{l odd} H_l C_l(t),  H_1 = -5/6, H_3 = 3/40, H_5 = -247/1680, ...
    E = 2 pi sum_l N_l H_l int_{rho0}^R rho^3 A_l drho ,   N_l = 2(l+1)(l+2)/(2l+3)

Checks: `H_l` exact-vs-quadrature agree to 14 digits at l = 1, 5, 21, 61; **every mode's contribution is
non-negative** (it is that mode's Dirichlet energy) — 0 negative modes at all four R; partial sums reach
1 - 1.3e-6 by l = 201.

| R/rho0 | E / (M^2 R^5) | a_inner = (M/2)log(R/rho0) | log Re_E | 2 log(R/rho0) | t_d M | **c1** |
|---|---|---|---|---|---|---|
| 2^8  | 1.72403978e-01 | 2.772589 | 10.387189 | 11.090355 | 0.25000000 | 2.596797 |
| 2^12 | 1.72403978e-01 | 4.158883 | 15.932366 | 16.635532 | 0.16666667 | 2.655394 |
| 2^16 | 1.72403978e-01 | 5.545177 | 21.477544 | 22.180710 | 0.12500000 | 2.684693 |
| 2^20 | 1.72403978e-01 | 6.931472 | 27.022721 | 27.725887 | 0.10000000 | 2.702272 |

`E/(M^2 R^5)` is **constant to 9 digits** across four decades of R/rho0 — the inner cutoff contributes
nothing at leading order, and there is no log enhancement (the l = 1 mode's `rho log rho` particular
solution is exactly cancelled by its homogeneous matching). So

    E_shell = 0.172403978 M^2 R^5 ,   log Re_E = 2 log(R/rho0) - 0.7031660 ,

    ****  c1 = 4 ln 2  -  (2 ln 2)(0.7031660) / log(R/rho0)   ->   4 ln 2 = 2.7725887222   ****

approached **from below**. Reproduced by all four rows to 6 digits.

---

## 4. Second order at t = 0, and the window on which first order is valid

### EXACT ALGEBRA (derived in the header of `c4_second_order.py`, then confirmed numerically)
Following a fluid particle, `dr/dt = u^r`, so `D_t a = D_t(u^r)/r - a^2`. For axisymmetric **no-swirl**
flow the radial momentum equation carries no centrifugal term, `D_t u^r = -d_r p + nu(Delta - 1/r^2)u^r`.
Hence at a core, where `omega^theta = r eta` with `eta` materially conserved,

    d^2/dt^2 omega^theta = omega^theta ( D_t a + a^2 ) = omega^theta * D_t(u^r)/r
                         = omega^theta * [ -(1/r) d_r p  +  (nu/r)(Delta - 1/r^2) u^r ] .

**The `a^2` compounding term is exactly cancelled by the `-(u^r)^2/r^2` inside `D_t a`.** The genuine
second-order driver is `beta := D_t a + a^2 = -(1/r) d_r p` (inviscid) — the radial pressure gradient,
not `a^2`. This is an identity, and it is the cleanest thing this sub-seat found.

### Numerics (self-consistent evolution of eta under its own flow)
`eta` is transported, so each core keeps its peak, translates with `u(x_k)`, and its covariance evolves
as `C -> (I + tau G) C (I + tau G)^T`, `G = grad u` at the core — exact to first order in tau and
leading order in sigma/rho. `D_t a = [a(tau) - a(-tau)]/(2 tau)` at the **moving** inner core.

Internal checks at N = 4, s = 0.10: incompressibility residual `d_r u^r + u^r/r + d_z u^z = 1.7e-6`
against a scale of 8.0e-2; and `r_0 (D_t a + a^2)` matches `D_t(u^r)` computed the same way to
rel 4.9e-7 (deforming cores) / 7.9e-10 (translation only) — i.e. the exact algebra above is reproduced
numerically. tau-independence over tau = 4e-3, 2e-3, 1e-3 at 8-9 digits.

**Result at N = 4, s = 0.10** (un-normalised `a = 6.29796899e-02 M`):

| N | variant | `a/M` | `D_t a / M^2` | `beta = D_t a + a^2` | `beta/a^2` | `t*M := 2a/\|beta\|` | `t_d M = ln2/a` | `t*/t_d` |
|---|---|---|---|---|---|---|---|---|
| 4 | translation only | 6.29797e-02 | -1.3503815e-03 | +2.6160598e-03 | 0.6596 | 48.15 | 11.006 | 4.37 |
| 4 | + core deformation | 6.29797e-02 | -1.9667318e-03 | +1.9997095e-03 | 0.5042 | 62.99 | 11.006 | **5.72** |
| 8 | translation only | 1.34615e-01 | -4.0853103e-03 | +1.4035887e-02 | 0.7746 | 19.18 | 5.149 | 3.73 |
| 8 | + core deformation | 1.34615e-01 | -1.0373129e-02 | +7.7480683e-03 | 0.4276 | 34.75 | 5.149 | **6.75** |
| 12 | translation only | 2.05083e-01 | -2.3863912e-03 | +3.9672778e-02 | 0.9432 |  10.34 | 3.380 | 3.06 |
| 12 | + core deformation | 2.05083e-01 | -1.9160456e-02 | +2.2898713e-02 | 0.5444 |  17.91 | 3.380 | **5.30** |

Two things to read off:

1. **`D_t a < 0` but `beta > 0`.** The strain itself *decays* at t = 0 (the stack starts to disassemble),
   but the second-order coefficient of `omega` is positive: `omega(t) = omega_0 (1 + a t + (beta/2) t^2 + ...)`
   with `beta > 0`, so the true growth is **faster** than first order, not slower. The first-order
   doubling time is an **upper** bound on the true doubling time over this window — which is the
   direction a lower-bound construction needs.
2. **`t*/t_d > 5` at every N tested.** `beta/a^2` = 0.504, 0.428, 0.544 and `t*/t_d` = 5.72, 6.75,
   5.30 at N = 4, 8, 12 — no drift with N, consistent with `beta/a^2 ~ 1/2` and `t*/t_d = 2/( (beta/a^2) ln2)
   ~ 5.8`. So first order is valid for **~5 doubling times** and the first-order doubling time is a
   legitimate estimate of the actual one. Core deformation matters at the 25-100% level in `beta` and
   *widens* the window (translation-only would give 4.37, 3.73, 3.06); including it is the more accurate
   and here the more favourable choice.
3. **The stack does move a lot, and that is already accounted for.** At N = 12 the inner core is swept
   at `u^z = -86 M rho0` — over one doubling time that is 291 rho0, i.e. 14% of the stack's outer
   radius `2^11 rho0`. Most of it is a near-uniform sweep that does not change the relative geometry:
   the quantity that matters, `D_t a / a = -0.093 M`, says the strain itself e-folds *down* only on the
   timescale 10.7/M, about 3 doubling times. That is the honest statement of how long the configuration
   holds together, and it is consistent with `t*/t_d = 5.3`.

**Caveat, stated plainly.** This `D_t a` is a *self-consistent estimate*, not an exact evaluation: it
uses the exact velocity and velocity gradient at each core and the exact linearised flow map, but it
truncates the core's evolution at its first two moments (translation + affine deformation), dropping
`O(sigma^2/rho^2)` shape corrections and the ring's own log-divergent self-propulsion beyond the
translation captured by `u(x_k)`. The internal check `r_0(D_t a + a^2) = D_t u^r` is passed by
construction of the same evolution, so it validates the *algebra*, not the truncation. The independent
statement that does not depend on the truncation is the identity `beta = -(1/r) d_r p`.

Self-propulsion timescale, for the record: `Gamma_k = 2 pi sigma_k^2 M`, self-induced speed
`~ (Gamma/(4 pi rho))(log(8 rho/sigma))`, so the time for a ring to translate its own radius is
`~ 2/(M s^2 log(8/s))`; the ratio to `t_d = ln2/(kappa M N)` is `~ 0.19 log(8/s)/N`, i.e. **< 1 for
N >~ 5** — the doubling happens before the stack rearranges. Consistent with `t*/t_d = 5.7`.

---

## 5. Does viscosity prevent it? — the decisive comparison (`c5_viscous.py`, `c8_viscous_shell.py`)

`D_t eta = nu Delta_5 eta`, so the relative first-order rate of the peak of `|omega^theta| = |r eta|` at a
material point is exactly

    (d/dt log|omega^theta|)|_{t=0} = a + nu (Delta_5 eta)/eta .

Both terms below are **exact symbolic** results (`c8`), not scalings:

* **Gaussian ring core**, `eta = -(M/r_k)exp(-|d|^2/(2 sigma^2))`: at the core centre
  `nu Delta_5 eta/eta = -2 nu/sigma^2` exactly. With the floor `sigma_0 = sqrt(nu/M)` the innermost
  core pays a penalty of exactly `2M` and ring j pays `2M 4^{-j}`.
* **Bang-bang shell**, `eta = -M sgn(z)/r`: `Delta_5(1/r) = 2/r^3 - 3/r^3 = -1/r^3`, so
  `nu Delta_5 eta/eta = -nu/r^2` exactly at every interior point (cross-checked in omega form:
  `(Delta - 1/r^2)omega^theta/omega^theta = -1/r^2` for locally constant `omega^theta`).
  With `rho0 = sqrt(nu/M)` the penalty at the **innermost** radius is `M/sin^2 phi = 1.5 M` at phi*,
  and it falls as `rho^{-2}` outward.

**A plateau has no curvature at its own peak; a ring core has `1/sigma^2`.** That is the whole
difference, and it is worth an enormous factor:

| datum | condition for net first-order growth at the inner scale | threshold |
|---|---|---|
| bang-bang shell | `(M/2)log(R/rho0) > 1.5 M` | `log Re_E > 5.297`, i.e. **Re_E > 200** |
| ring stack s=0.125 | `kappa N > 2` | `log Re_E > 103.8` (Re_E > 10^45) |
| ring stack s=0.10 | `kappa N > 2` | `log Re_E > 158.7` (Re_E > 10^69) |
| ring stack s=0.04 | `kappa N > 2` | `log Re_E > 961.1` (Re_E > 10^417) |

For the ring stack the *innermost* core therefore never wins at any physical Reynolds number. But the
question asked is about `sup|omega|`, and the sup simply moves outward. With
`g_j(N) = a_j(N) - 2 M 4^{-j}` (`c5_viscous.py`) the maximiser is interior:

    dg/dj = 0  =>  4^{-j*} = kappa/(2 ln 4) ,  j* = log_4(2 ln 4 / kappa) ,
    max_j g_j = M kappa ( N - j* - 1/ln 4 + b/kappa ) .

At s = 0.10, `j* = 3.655`, so the loss is `j* + 1/ln4 = 4.38` octaves. Checked against the table:
N = 26, s = 0.10 gives `max_j g_j/M = 0.36868` at j = 4, versus `kappa(N - 0.345 - 4.376) = 0.37170` —
agreeing to the integer-vs-continuous-j discrepancy. Measured argmax j: 5 (s = 0.04), 4 (s = 0.06,
0.08, 0.10), 3 (s = 0.125); net growth of the sup starts at N ~ 5-8 for all five s.

**So: viscosity does not destroy the mechanism — it costs `j* + 1/ln 4 = O(log(1/kappa))` octaves,
which is ADDITIVE in N and therefore `O(1/N)` in `c1`.** The law `t_d = c1/(M log Re_E)` survives with
the same `c1`. For the shell the cost is a bounded O(1) penalty from the first octave only.

---

## 6. Answer to the brief

**Q: is `a_inner = (M/2)(N ln 2) + O(M)`?**
No, not for a genuine (well-separated) ring stack. Exactly:
`a_inner(N) = M[kappa(s) N + b(s)]` with `kappa(s) = (pi/sqrt3) s^2 (1 - 3 s^2) + O(s^6)` and
`b/kappa ~ -0.35`. The `(M/2) ln 2` per octave is a **space-filling** rate; a ring stack recovers it
only as `s -> 0.437`, where the rings have merged into the bang-bang shell. The self-strain of any
z-symmetric core is **exactly zero** (Theorem A), so a single ring gives 0, not O(M).

**Q: c1 in `t_d = c1/(M log Re_E)`?**

    ring stack:   c1(N) -> 2 (ln 2)^2 / kappa(s)   = 333.1, 149.2, 84.8, 55.0, 36.0 for s = 0.04..0.125
    bang-bang shell (the extremal endpoint of the same family, computed exactly):
                  c1 = 4 ln 2 - (2 ln 2)(0.7031660)/log(R/rho0)  ->  ****  c1 = 4 ln 2 = 2.7725887  ****

Both use `E` computed from the exact identity `E = (1/pi)||grad psi||^2_{L^2(R^5)}`; for the shell
`E = 0.172403978 M^2 R^5` to 9 digits, independent of `R/rho0`.

**Q: does the viscous NS evolution realise a doubling time `<= C/(M log Re_E)`, or does something
prevent it?**

At **first order in t**, and modulo the second-order window being wide (`t*/t_d = 5.7`, and the
second-order correction has the *favourable* sign, `beta > 0`), the answer is **yes, with the upper
bound's form matched exactly**, and nothing in this computation prevents it:

* the strain is `a = (M/4)(log Re_E + 0.703)` for the extremal datum — the same `M log Re` as the
  upper bound `H_0 = c/(M_0(1 + log_+ Re_0))`, with no gap in the power of the log;
* viscosity costs a bounded additive number of octaves (O(1) for the shell, `log_4(2 ln4/kappa) + 1/ln4`
  for the ring stack) and therefore does not change `c1`;
* the deformation of the stack shows up in `D_t a < 0` — the strain *does* start to decay — but the
  second-order coefficient of `omega` itself is `beta = -(1/r)d_r p > 0`, so `omega` accelerates.

**The quotable consequence.** If the first-order rate is actually realised, the universal constant in
the upper bound cannot exceed

    a = (M/4)(log Re_E + 0.703)  =>  time to reach (3/2)M_0  is  4 log(3/2)/(M log Re_E) + O(1/(M log^2))
    ****  c <= 4 log(3/2) = 1.6218604  ****

in `H_0 = c/(M_0(1 + log_+ Re_0))`. Astra's upper bound and this lower-bound construction differ only
in the constant, and the constant is pinned between 0 and 1.6219.

**What is NOT established here.** Everything above is a statement about `t = 0` (first and second
derivative) for a *specific* datum. It is not a proof that the solution's `||omega||_inf` actually
doubles in time `c1/(M log Re_E)`: that needs control of the evolution over the whole interval, and the
honest status is that the first-order rate holds, the second-order correction has the right sign and a
window of ~5.7 doubling times, and the O(M) corrections and core deformation are bounded but not
controlled to all orders. The ring stack is also the *wrong* extremiser — the shell beats it on both
counts (rate by `1/(5.234 s^2)`, viscous penalty by `1/sin^2 phi` versus `2/s^2` in units of M).

**Recommendation for the campaign.** Drop the discrete-ring realisation. The object that carries the
mechanism is the **mollified bang-bang plateau**: it is extremal for the strain (`(M/2) log(R/rho0)`
exactly), it is viscously neutral at O(1) Reynolds number because `Delta omega = 0` on a plateau, and it
gives `c1 = 4 ln 2` exactly. Rings are a smooth *approximation* to it that pays the volume fraction
`5.234 s^2` twice over.

---

## 7. Files

| file | what it establishes |
|---|---|
| `ns5d.py` | the 5D-lift Biot-Savart library; exact closed forms for J5, J3 in terms of 2F1 and the elliptic modulus |
| `c1_calibrate.py` / `c1_log.txt` / `c1_results.json` | kernel vs mpmath; Hill `a = Az/5` and `psi`; Hill energy; shell axis identity; sign control |
| `c2_strain.py` / `c2_log.txt` / `c2_results.json` | `A_m(s)`, `kappa(s)`, `a_inner(N)`, `\|\|omega\|\|_inf`, all three controls |
| `c3_energy.py` / `c3_log.txt` / `c3_results.json` | `h_j(s)`, `E(N)`, `Re_E(N)`, `c1(N,s)` for both viscous conventions |
| `c4_second_order.py` / `c4_log.txt` / `c4_results.json` | `D_t a`, `beta = D_t a + a^2 = -(1/r)d_r p`, validity window `t*/t_d` |
| `c5_viscous.py` / `c5_results.json` | `g_j = a_j - 2M 4^{-j}`, the interior maximiser `j*`, thresholds |
| `c6_shell_exact.py` / `c6_log.txt` / `c6_results.json` | the extremal endpoint exactly: `E = 0.172403978 M^2 R^5`, `c1 -> 4 ln 2` |
| `c7_kappa_analytic.py` / `c7_results.json` | `kappa(s) = (pi/sqrt3)s^2(1-3s^2)`, from `Laplacian f/f = -6` at phi* |
| `c8_viscous_shell.py` | exact symbolic viscous rates: `-nu/r^2` (plateau) vs `-2nu/sigma^2` (ring core) |

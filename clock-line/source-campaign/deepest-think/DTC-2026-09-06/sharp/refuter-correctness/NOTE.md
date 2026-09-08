# REFUTER (mathematical-correctness lens) — the two-sided doubling clock

Seat label: `refuter-correctness`.  Everything below was produced by scripts in
`scripts/`; reruns of the two construct seats are in `rerun-exact-first-order/` and
`rerun-viscous-numerics/`.  Nothing outside this folder was written.

## 0. Verdict

**REFUTED as written — MAJOR.**  Three defects, in decreasing order of severity:

1. the statement's quantifiers make its left half **false** (§1);
2. its right half is **not a theorem** — it is a t=0 derivative plus a six-point
   simulation (§2);
3. a **refereed published theorem contradicts the right half** (§3), and until that is
   settled in print the upper bound cannot be asserted.

The content survives in corrected form (§7).  Two of my own attempted refutations of the
construct seats **failed** — the mechanism's constant κ = 1/2 per e-fold is right where the
vorticity actually is (§5), and the admissible class reaches it (§6).

## 1. The quantifier defect (script `r1_quantifier_scaling.py`)

Candidate: `c1/(M(1+log_+Re_E)) <= inf_data T_d <= c2/(M log Re_E)`.  The middle term is a
single number; the outer terms depend on a free datum through `M` and `Re_E`.

Under the NS scaling `u_lam(x,t) = lam u(lam x, lam^2 t)` at fixed nu, measured on an
explicit smooth divergence-free field (div control 9.7e-12):

| lam | E/(E1/lam) | M/(lam^2 M1) | Re_E |
|---|---|---|---|
| 0.5 | 0.9939 | 1.000000000 | 0.9703929 |
| 1.0 | 1.0000 | 1.000000000 | 0.9727595 |
| 2.0 | 1.0000 | 1.000000000 | 0.9727595 |
| 3.0 | 1.0000 | 1.000000000 | 0.9727595 |

`Re_E` is scale invariant (10 digits for lam >= 1; the lam = 0.5 row is box truncation) and
`M T_d` is scale invariant because `T_d ~ lam^-2`.  Hence at **any** fixed `Re_E`,
`inf_data T_d = 0` (send lam -> infinity).  So, read literally, the left inequality is
false for every `c1 > 0` and the right inequality is true but empty.  This is a
book-keeping fault, not a fluid-mechanical one; the fix is §7.

Also checked: `Re_E` uses the **initial** energy, and that choice is harmless — on `[0,T_d]`
energy only decreases and `M` grows by at most 3/2, so `Re_E(t) <= (3/2)^{1/5} Re_E(0)`,
i.e. `log Re_E` moves by at most **+0.0811** (script `r6`).

## 2. The right half is not established

* `exact-first-order` states it itself (its caveat 1): everything is the **first and second
  derivative at t = 0** for one datum.
* `viscous-numerics` is six points, `log Re_E` from 6.46 to 13.40 — a factor 2.07.  Numerics
  falsify, never prove.

Worse, the pre-registered gate is **not diagnostic for the hypothesis it was written for**
(script `r6_gate_diagnosticity.py`).  Fitting the same six points with a *saturating* model
`T32*M = A + c/(logRe - d)` — the exact negation of "the log is necessary":

| forced floor A | c | d | rms | rms / rms(A=0) | beta |
|---|---|---|---|---|---|
| 0.000 | 1.8659 | 5.0730 | 2.66e-03 | 1.00 | 1.4185 |
| **0.011 (best fit)** | 1.7990 | 5.1127 | **6.51e-04** | **0.24** | – |
| 0.020 | 1.7438 | 5.1456 | 2.30e-03 | 0.86 | 1.3931 |
| 0.050 | 1.5690 | 5.2519 | 9.71e-03 | 3.65 | 1.3539 |
| 0.100 | 1.2993 | 5.4216 | 2.30e-02 | 8.64 | 1.2855 |
| 0.150 | 1.0557 | 5.5823 | 3.74e-02 | 14.06 | 1.2137 |

The unconstrained fit **prefers a positive floor** (A = +0.0109, rms four times better),
and every saturating model clears the gate's `beta >= 0.7` survive line.  The data
therefore *bound c2 from above*; they do not show `T32*M -> 0`.  A floor of 0.05, i.e. a
log-free clock with `c(3/2) = 20`, is not excluded by this evidence.  (Reaching
`T32*M = 1/20` along the fitted log model needs `log Re_E = 42.4`, ~28 octaves.)

`beta = 1.4059` reproduced exactly from the raw histories; it exceeds the first-order
prediction 1 only because of the additive offset `d`, and is not a shape discriminator.

## 3. The published contradiction (scripts `r3_bfg_gap.py`; source read at
`../literature-short-time/txt/bfg.txt`)

Bradshaw–Farhat–Grujic arXiv:1704.05546v4 = ARMA 2019, Thm 10: for `omega_0 in L^2 ∩ L^inf`
and `M > 1` there is `c(M)` with a mild solution on `T >= 1/(c(M)||omega_0||_inf)` and
`||omega(t)||_{L^inf(Omega_t)} <= M ||omega_0||_inf`.

**I verified at source that `Omega_t` is a complex strip** `{x+iy in C^3 : |y| < sqrt(t/c(M))}`,
not a shrinking spatial subdomain.  It contains `R^3`, so the bound is a genuine global sup
bound on the real vorticity, and with `M = 3/2` it says `T_d >= 1/(c(3/2)M_0)` with **no
Reynolds number**.  That contradicts the candidate's right half for large `Re_E`.  The
literature seat's reading is correct and the conflict is real.

Their sketch (p.10) turns on `int |f(x)|/(|x|+1)^4 dx <= c ||f||_BMO`, which they
parenthetically concede needs the local average subtracted, and they advertise "no
quantitative dependence on `||omega_0||_2`" (p.8, verbatim, line 411).

**That step fails on this campaign's own extremiser.**  For the bang-bang shell (`||omega||_inf
= M = 1` for every `R`), `grad u` is smooth in the empty hole and has eigenvalues
`(a, a, -2a)` with `a = u^r/r`, so `|grad u| = sqrt(6) a`:

| R/rho0 | a(0,0) | (M/2)log R | min a on \|x\|<=1/2 | lower bound for int \|grad u\|/(1+\|x\|)^4 |
|---|---|---|---|---|
| 4 | 0.693150 | 0.693147 | 0.681594 | 0.259016 |
| 16 | 1.386293 | 1.386294 | 1.372979 | 0.521752 |
| 64 | 2.079888 | 2.079442 | 2.066462 | 0.785286 |
| 256 | 2.777843 | 2.772589 | 2.765041 | 1.050757 |
| 1024 | 3.446639 | 3.465736 | 3.433371 | 1.304732 |

The lower bound grows **linearly in log(R/rho0)** (slope 0.18902 per e-fold, intercept
-0.00182) while `||grad u||_BMO <= C||omega||_inf = CM` is bounded uniformly in `R`.  So the
inequality is false with an absolute constant, and the deficit is exactly one power of the
clock's logarithm.

Restoring the omitted term with `|f_{B(x,L)}| <= CM(1 + log_+(ell/L))` (the `L >= ell` case
is Astra pass 8 eq. (6); the `L < ell` case is its own telescoping) turns BFG's own string
into `tau[3/2 + (1/2)log(Re/tau)] <= 1/(2C)`, whose solutions have `tau log Re` tending to a
constant: 0.465, 0.618, 0.698, 0.783, 0.869, 0.924, 0.957 at `Re = 1e2 … 1e80` — a
logarithmic clock.  Dropping it gives `tau <= 1/(2C)`, log-free.  **The whole discrepancy
between BFG Thm 8/10 and Astra pass 8 is that one term.**

This is my reading of a sketch, not an erratum, and I found no published challenge.  It is
enough to say the campaign must not present the upper bound as settled while BFG stands.

## 4. Reproduction of the two construct seats

* `exact-first-order`: all 10 `.py` hash-verified against its `SHA256SUMS`.  Re-run from a
  copy — `c2_strain.py` reproduces its `c2_results.json` **bit for bit** (no field differs by
  more than 1e-12 relative; scan clean); `c6_shell_exact.py`, `c7`, `c8` reproduce every
  printed digit; `c1_calibrate.py` reproduces including its single documented failing row
  (`FAILURES: ['near_coincidence']`, exit 1) — the caveat is accurate.
* `viscous-numerics`: 16 files (`.py` + `results_*.json`) hash-verified.  `analyse.py`
  re-run: output **byte-identical** to the stored `analysis_log.txt`.  `validate.py` re-run:
  `validate.json` identical.
* Independent re-extraction of `T32` from the raw `hist` arrays by linear interpolation of
  the first upcrossing of `1.5*om0` (script `r2`): agrees with every stored `T32` to `<= 7e-16`
  relative, and `Q' = 0.9880 +- 0.0184` reproduces the seat's headline exactly.

## 5. An attempted refutation that FAILED — the strain at the empty origin (scripts
`r4_shell_strain_independent.py`, `r4b_conv.py`, `r5b_axis_diag.py`)

`exact-first-order`'s `c1 = 4 ln 2` rests on `a_inner = (M/2)log(R/rho0)`, which is the strain
**at the origin** — the centre of the empty hole, where there is no vorticity to stretch.
What sets the growth rate of `||omega||_inf` is `a` where `|omega|` is maximal.

I built an **independent instrument** for this: the classical elliptic-integral vortex-ring
velocity (Lamb), integrated in polar coordinates about the evaluation point so the
principal-value singularity becomes bounded.  Controls, all of which fire:

* vs a direct 3D Biot–Savart line integral of the same filament: worst rel **3.1e-13**;
* Hill's spherical vortex `omega^theta = A r` must give `a = A z/5`: worst rel **4.7e-09** at 5 points;
* the shell axis identity `a(0+,0) = (M/2)log(R/rho0)`: rel **2.1e-06** (R=64), **1.7e-06** (R=4096).

Result, after fixing an under-resolution in my own first pass (geometric radial subdivision;
`r4b`), at the inner edge and phi = 10 deg:

| R/rho0 | a(material) | (M/2)log R | difference |
|---|---|---|---|
| 64 | 2.296241 | 2.079442 | 0.216800 |
| 256 | 2.989364 | 2.772589 | 0.216775 |
| 1024 | 3.682508 | 3.465736 | 0.216772 |
| 4096 | 4.375657 | 4.158883 | 0.216774 |

`a_material = (M/2)log(R/rho0) + 0.216773 +- 2e-6`, a **constant** offset over six octaves.
So kappa = 1/2 per e-fold is correct on the material, and the O(M) correction has the
favourable sign.  **The construct's headline constant survives this test.**

Recorded deviation: my first pass (`r4`, nsplit = 1) reported `kappa_mat = 0.473` from an
under-resolved radial quadrature; `r4b` supersedes it.  Likewise the `kappa_axis` column of
`r5` is void — `r5b` shows the `r0 = 1e-5` evaluation point carries 31% error at R = 4096
while `r0 in [1e-3, 1e-1]` reproduces the exact identity to 2e-6.  Only the material column
of `r5` is used.

## 6. A second attempted refutation that FAILED — admissibility (script `r5_admissible_taper.py`)

The bang-bang shell is **not** admissible data: `omega^theta` does not vanish on the axis, so
`omega = omega^theta e_theta` is discontinuous there and `eta = omega^theta/r` is unbounded.
`exact-first-order`'s "extremal endpoint" is therefore outside the candidate's class, and
`viscous-numerics` reports its own repaired datum's `kappa = 0.41` as the price of
admissibility (its caveat 6).

Tested a family that vanishes only near the axis, `omega^theta = -M sgn(z) min(1, phi_ax/delta)`
(`|omega| <= M`, `eta <= M/(rho0 sin delta)`, admissible after mollification).  kappa measured as an
increment over six octaves (R = 64 -> 4096), so additive offsets cancel exactly:

| datum | kappa on the material | c2(3/2) = 2 ln(3/2)/kappa |
|---|---|---|
| bang-bang (inadmissible) | **0.50000** | 1.6219 |
| taper delta = 30 deg | 0.48363 | 1.6768 |
| taper delta = 15 deg | 0.49781 | 1.6290 |
| taper delta = 7.5 deg | 0.49972 | 1.6228 |
| `-M sin 2phi` (viscous-numerics' datum) | **0.40000** | 2.0273 |

So the admissible class reaches `kappa = 1/2` to any accuracy.  The seat's 0.41 is the cost
of **also** killing the equator, not the cost of admissibility — its caveat 6 is
pessimistic, and the first-order `c2(3/2) = 4 ln(3/2) = 1.6219` stands for admissible data.
(The taper does raise the datum's smallest feature to `rho0 sin delta`, which the viscous floor
must respect; that is an O(log(1/delta)) bookkeeping cost, not evaluated here.)

## 7. Corrected statement

Fix `nu > 0`.  For smooth finite-energy `u_0` let `M_0 = ||omega_0||_inf`,
`E_0 = ||u_0||_2^2`, `Re_E = E_0^{2/5} M_0^{1/5}/nu`, and let `T_d` be the first time
`||omega||_inf = (3/2)M_0` (`+infinity` if never).  Both `M_0 T_d` and `Re_E` are invariant
under the NS scaling, so

    T(Re) := inf { M_0 T_d(u_0) : u_0 smooth, finite energy, Re_E(u_0) <= Re }

is well defined, dimensionless and non-increasing in `Re`.  Then

**(i) THEOREM** (Astra pass 8, referee-confirmed end to end).  `T(Re) >= c1/(1 + log_+ Re)`
for a universal `c1 > 0`.  Equivalently every smooth finite-energy solution keeps
`||omega||_inf <= (3/2)M_0` for a time `c1/(M_0(1+log_+ Re_E))`.

**(ii) CONJECTURE**, not proved by anything in this campaign.  `T(Re) <= c2/log Re` for all
`Re >= Re_*`, realised by mollified axisymmetric no-swirl plateau data: `omega^theta = -M sgn z`
on `rho0 < |x| < R`, tapered to zero within angle `delta` of the axis, `rho0 >~ sqrt(nu/M)`.
Evidence, at the stated strength and no more:
  * *exact*, and independently confirmed here on the material rather than at the empty
    origin: the initial strain is `a = (M/2)log(R/rho0) + O(M)`, whence a first-order
    `c2(3/2) = 4 ln(3/2) = 1.62186`;
  * *numerical*, six points, `Re_E <= 6.6e5`: `M T_d log Re_E` falls monotonically 8.69 -> 3.04
    and is consistent with an asymptote near 2, with the invariant `M T_d log(R/s*) = 0.988 +- 0.018`;
  * the named spoilers (viscous death of the inner rings, stack deformation, O(M)
    corrections) do not stop the mechanism over the range tested.
The evidence bounds `c2`; it does **not** establish that `T(Re) -> 0` (§2).

**(iii) OPEN CONFLICT.**  Bradshaw–Farhat–Grujic (ARMA 2019) Theorem 10 asserts
`T(Re) >= 1/c(3/2)` uniformly in `Re`, contradicting (ii).  Their sketch's weighted-BMO step
fails on the family of (ii) by a factor `log(R/rho0)` (§3), and restoring the omitted term
returns the shape of (i).  (ii) may not be asserted, and "sharp up to constants" may not be
claimed, until that is resolved in print.

Also: `exact-first-order`'s "`c1 = 4 ln 2`" is the constant of the **upper** bound in the
candidate's labelling (a `c2`), and it is the factor-**2** doubling constant; the candidate's
threshold is 3/2, for which the first-order value is `4 ln(3/2) = 1.62186`.

## 8. Smaller findings

* `exact-first-order` control (c) declares "A_m for m<0 must decay ~4^m", prints the expected
  ratio 0.25, measures 0.051/0.058/0.060, and the summary reinterprets the target as 1/16
  after the fact.  The reinterpretation is very likely right; the control as pre-declared
  failed and was re-cut post hoc.
* `literature-short-time` says the only "Navier–Stokes" strings in Kim–Jeong are three
  bibliography lines; there is also a body mention (intro, line 63 of the extraction).  The
  Euler-only verdict stands: **zero** occurrences of "viscous"/"viscosity" in that paper, and
  Bang–Cheskidov has 2 "Navier" (both bibliography) and 1 "viscous".
* Arithmetic tie-out: `0.172403978^{2/5} = 0.495015631`, `log = -0.7031659` vs the quoted
  `-0.7031660`.

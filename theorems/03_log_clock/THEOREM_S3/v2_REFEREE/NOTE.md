# REFEREE REPORT on THEOREM (S3) version 2, and on `fix4`

Seat `s3close/round2/refute-THEOREM_S3_v2`, single referee, sitting of 2026-09-08.
Laws: `TORMENT NEXUS/LAWS.md` first 120 lines; seat header `TEMPLATES/SEAT_HEADER_2026-09-08.md`.
Posture: FL-000 will fall; find the move, and find the hole. This report finds the hole.

**VERDICT: NOT PROVED as stated.** One FATAL, five MAJOR, five MINOR. The fatal item is a
hypothesis of an input theorem that the mollification of `fix4` destroys and that no sheet
recomputed. It is repairable, and the repair costs nothing in the budget, but it is not made
and the theorem as written is false in its claim that "every analytic item that
`THEOREM_S3.md` sec.3 carried as unproved is discharged".

Every number below came out of a script in `scripts/`, run on COPIES in a scratch directory,
with results in `results/`. Nothing under `round2/THEOREM_S3/` was read-write; the target
files are untouched. `SHA256SUMS` computed with `shasum -a 256`.

**FL-000 stands. Nothing here touches the headline problem.**

---

## 0. WHAT WAS RE-RUN, AND WHAT REPRODUCED

`fix4/check_fix4.py` on a byte copy: **251 CHECKS, 251 PASS, 0 FAIL**. The sheet is
internally consistent and its transcription is clean. The findings below are not transcription
errors.

An instrument written here from the mathematics of `THEOREM_S3_v2.md` sec.1.1 alone
(`s1_profile_indep.py`: adaptive-quadrature Gaussian convolution with the kink breakpoints
passed explicitly, plus rigorous derivative bounds from
`||W_sigma^{(k+1)}||_inf <= ||Wtil'||_inf ||chi_sigma^{(k)}||_1`; no shared code with
`fix4/p1_profile.py`, which is an FFT on a 600001-point grid) reproduces:

| constant | `fix4` | here | relative |
|---|---|---|---|
| `N_sigma` | `1.0124508488` | enclosure `[1.0124508487, 1.0124508495]` | contains |
| `kappa_delta` | `0.4917868801` | `[0.4917868796, 0.4917868800]` | contains |
| `c_*` | `0.8244732108` | `0.8244732112` | `5.2e-10` |
| `c_2` | `1.6489464217` | `1.6489464227` | `5.6e-10` |
| clock floor | `1.6700567265e-02` | `1.6700567884e-02` | `3.7e-08` |
| `lam_max(c_*)` | `1.8558724485` | `1.8558724492` | `3.7e-10` |
| `E_0` | `7.5523076942` | `[7.5523076888, 7.5523076951]` | contains |
| `kappa_delta` kinked control | `0.4978224383` | `0.4978224382` | `1.2e-10` |
| `r_h` grid minimum at `c_*` | `0.9035285738` | `0.9035285735` | `3.3e-10` |
| `P_h'` at `lam_max(c_*)` | `0.3549367307` | `0.3549090746` | `7.8e-05` |
| `lam_mono` | `2.0693384635` | `2.0693179283` | `9.9e-06` |
| monotone cap `c/c_*` | `1.1760705119` | `1.1760544630` | `1.4e-05` |
| `int_{-1}^1 w^2 dt` | `1.8239293735` | `1.8239293727` | `4.4e-10` |
| `||omega_0||_2^2/(M^2R^3)` | `1.4381082627` | `1.4381082621` | `4.2e-10` |
| `log Re_E` shift | `2.9135781820` | `2.9135781816` | `1.3e-10` |
| `L_*` proved, `eps<=cap` | `1424610.4953` | `1424610.4953` | exact re-run |
| `L_*` proved, `eps<=0.1` | `1555469.4004` | `1555469.4004` | exact re-run |
| `log Lambda_*` | `2849223.9041` | `2849223.9041` | exact re-run |

The arithmetic is right. The chain is not.

---

## 1. FINDING F-1, **FATAL**: (H3'_vartheta) is FALSE for the datum of sec.1.1

Script `scripts/s3_vartheta.py`, results `results/s3_results.json`.

### What is claimed

Dependency row **B4** of `THEOREM_S3_v2.md` sec.2:

> `B4 | Theorem V.4' with (H3'_vartheta) for k = 0..4; vartheta(f = 4) = 0.0556 |
> round2/pmax-h3v Part B (ADDENDUM_1) | PROVED at f >= 4`

and `fix4/FIX4.md` sec.7 closes with "**M1 (P-max) and M2 (H3-V)** were discharged in
`ADDENDUM_1` and are untouched here."

`(H3'_vartheta)` (`pmax-h3v/PROOF.md` lines 674-680) asks for a number `vartheta` in `[0,1)`
with

```
    sup_{|v|=1} | d_v^k (eta_0 - eta_P)(x) |  <=  vartheta k! M / r(x)^{k+1} ,
    k = 0,1,2,3,4 ,  for every x in B(x_*, d) ,      eta_P := -M/r ,
    d = min( f , rho_* sin(phi_0 - delta) , rho_* sin(pi/2 - delta_m - phi_0) ) rho_0 .
```

Theorem V.4' needs `vartheta < 1` twice: through `Q = (1+vartheta)/(1-vartheta)` and through
`vt = vartheta/(1-vartheta)`. At `vartheta >= 1` the theorem has no content.

### Why `vartheta = 0.0556` was small, and why that reason is gone

`pmax-h3v` sec.B1 and sec.B3 item 9 say it in as many words:

> "the two angular mollifiers are `min(1, .)`, so on a ball that misses both kinks the angular
> factor is **exactly** `1` and contributes nothing; but the radial ramp is `tanh` ..."
> "On the ball the angular profile is identically `1` (that is what `d` is chosen for), so the
> entire deviation is the `tanh` radial ramp".

`vartheta(f=4) = 0.055612` is therefore a statement about `psi = 1 - Theta` alone. It is a
number for the KINKED angular profile and for no other.

`fix4` replaces that angular profile by its Gaussian mollification. Three things follow at
once:

1. The mollified profile is **nowhere** identically `1`. On the bulk it exceeds `1` by up to
   `1.2450848776 %` before normalisation and falls `1.0908879024 %` short after it, which
   `FIX4.md` sec.1.1 states itself: `A_sigma(phi_0)/N_sigma = 0.9890911210`.
2. At `f = 4` the binding branch of `d` is `rho_* sin(phi_0 - delta) = 5 sin(22.5 deg)
   = 1.9134171618`, so **the ball is exactly tangent to the cone `phi = delta`**. That is not
   an accident of `f`: the same branch binds at `f = 8, 16, 32`, the whole grid `FS_F4` the
   budget minimises over. The tangency is structural in `f`.
3. Mollifying at scale `sigma` puts the smoothing layer `|phi - delta| <~ 3 sigma` **inside**
   the ball, and in that layer the angular derivatives of the profile are `O(sigma^{-(k-1)})`.
   Exactly: with `Delta := W'(delta^+) - W'(delta^-) = -58.527913`,
   `W_sigma^{(k+1)} = Delta chi_sigma^{(k-1)}(phi - delta) + smooth`, so
   `W_sigma'''' (delta) = Delta chi_sigma''(0) = 2.9187e+06`.

### The measurement

`s3_vartheta.py` computes `r^{k+1}|d_v^k(eta_0 - eta_P)|/k!` exactly along straight lines
`x + t v` by truncated Taylor arithmetic to order 4 (`sqrt`, `log`, `atan`, reciprocal and
composition series), scanning `B(x_*,d)` and meridian directions. It returns a **lower** bound
on the supremum, which is all a refutation needs.

**Control (L-14).** The same instrument on the kinked datum at `f = 4`:

```
   vartheta_k , k = 0..4 :  0.0008923  0.00432952  0.0126102  0.0285659  0.0555009
   vartheta = 0.0555009        against  pmax-h3v's  0.055612       (2.0e-03 relative)
   vartheta_0 = 0.0008923      against  pmax-h3v's  8.9614e-04     (exact clause psi(rho_*-d))
```

The instrument is the right instrument.

**The test.** The `sigma = 0.02` mollified datum of `THEOREM_S3_v2.md` sec.1.1:

| `f` | `d` | `vartheta_0` | `vartheta_1` | `vartheta_2` | `vartheta_3` | `vartheta_4` | `vartheta` |
|---|---|---|---|---|---|---|---|
| **4** | `1.913417` | `0.0631978` | `0.597707` | `1.92543` | `4.14053` | **`8.76183`** | **`8.76183`** |
| 8 | `3.444151` | `0.0631647` | `0.597691` | `1.92545` | `4.14071` | `8.76295` | `8.76295` |
| 16 | `6.505618` | `0.0631644` | `0.597691` | `1.92545` | `4.14071` | `8.76296` | `8.76296` |
| 32 | `12.628553` | `0.0631644` | `0.597691` | `1.92545` | `4.14071` | `8.76296` | `8.76296` |

The supremum is attained at `phi = 9.6458 deg`, i.e. `1.87 sigma` above the taper cone, inside
the smoothing layer, at the binding order `k = 4` (the same order `pmax-h3v` reports for the
kinked datum at `f >= 0.25`).

**`vartheta >= 8.76` at every `f` the theorem admits.** The hypothesis fails at `k = 2` already
(`1.925 > 1`), and it fails at `k = 0` by a factor `70` against the kinked datum. There is no
admissible `vartheta`, so Theorem V.4' does not apply, `M2 (H3-V)` is not discharged for this
datum, and row `F4` (`eps_v = eps_bulk + E_hess + E_4 + E_tail`, "PROVED given B4") has no
warrant. `ADDENDUM_1`'s "CONSEQUENCE: the theorem must be stated at `f >= 4`" is a repair for
the kinked datum only; the failure here is angular, not radial, and larger `f` does not touch
it.

Note the exact shape of the trade `fix4` made and did not price. `fix3` chose to mollify the
angle in order to make `Psi(0)` finite (M3). `pmax-h3v` chose `min(1,.)` angular mollifiers
precisely **so that** the plateau hypothesis would hold on the ball (its sec.B1 says so). The
two choices are incompatible, and `fix4` took one of them without checking the other.

### The repair, and its price

Shrinking the Theorem V.4' ball so it misses the smoothing layer restores the hypothesis
(`results/s3_results.json`, `d_shrink_scan_f4`):

| `d` as a fraction of ASSEMBLY's `d` | `d/rho_0` | `vartheta` (lower bound) |
|---|---|---|
| `1.00` | `1.913417` | `8.712` |
| `0.90` | `1.722075` | `8.325` |
| **`0.80`** | **`1.530734`** | **`0.100621`** |
| `0.70` | `1.339392` | `0.0861577` |
| `0.50` | `0.956709` | `0.0422732` |

`0.80` is `d = rho_* sin(phi_0 - delta - 4 sigma)` to two digits: the ball has to stand about
four mollification widths off the cone. At `SAFETY = 2`, which is `pmax-h3v` sec.B5's own
convention for a grid lower bound on a supremum, `2 vartheta = 0.2012 < 1`, so the repaired
hypothesis is admissible.

**The budget does not notice** (`s4_budget_rerun.py`, `ball_shrink_cost`; the viscous budget is
`t3_budget.viscous_budget` transcribed with `d` exposed, checked identical to the imported
function at `dscale = 1` to `0.000e+00`):

| variant | `L_*` (proved, `eps <= cap`) | `eps_v` at `L = 2e+06` |
|---|---|---|
| stock | `1424610.4953` | `3.1295e-09` |
| `delta_eff = delta + 4 sigma` | `1424610.4953` | `3.3852e-09` |
| `delta_eff = delta + 6 sigma` | `1424610.4953` | `3.6931e-09` |
| `dscale = 0.8` | `1424610.4953` | `3.3951e-09` |
| `dscale = 0.5` | `1424610.4953` | `4.7928e-09` |

So the hole is cheap to fill in `L_*` and expensive in status: the theorem as written is not
proved, and the sheet that would prove it has not been run.

---

## 2. FINDING F-2, **MAJOR**: the budget runs at `vartheta = 0`, i.e. on V.4 and not V.4'

`pmax-h3v` sec.B5 states the corrected budget line explicitly:

```
   eps_v(vartheta) = eps_bulk [ 1 + vt (9 + 2 sigma_z/sigma_y) ] + Q ( E_hess + E_4 + E_tail )
```

and supplies `q4_budget_theta.py` to carry it. `t3_budget.viscous_budget`, which `fix4/p4`
imports unchanged and which produces every `eps_v` in `FIX4.md` sec.6.2 and in
`THEOREM_S3_v2.md`, returns

```
   eps_v = eps_bulk + E_hess + E_4 + E_tail
```

with no `Q` and no `E_vartheta`. Row `F4` of the dependency graph therefore quotes V.4, not
V.4', while row `B4` claims V.4'. Even granting the kinked datum's `vartheta = 0.0556`, the
"all proved" `eps_v` is the `vartheta = 0` number. Numerically the correction is invisible
(`eps_v ~ 3e-09` against `eps ~ 6.7e-02`); structurally it means the discharge of M2 was
recorded in `ADDENDUM_1` and then not wired into the instrument.

---

## 3. FINDING F-3, **MAJOR**: Proposition K's hypothesis is not satisfied by the datum

`THEOREM_S3_v2.md` sec.1.3 opens: "let `u` be its maximal strong solution, which exists and is
unique by Proposition K of `forced-route-2026-09/theorems/02_ADDENDUM_uniqueness_2026-09-08.md`
sec.3". Proposition K's hypothesis is "`u_0` smooth, divergence free and axisymmetric
satisfying (4)", and (4) is Fefferman's decay condition, quantified over every multi-index and
every `K`. Smooth plus (4) is the Schwartz condition, and the note says so.

The datum is not Schwartz and is not smooth on `R^3`: `omega_0` decays like `rho^{-8}` and is
`C^{7,1}` at the origin, so `omega_0 in H^s` only for `s < 9.5`. `FIX4.md` sec.7 M6 concedes
this ("*The correction.* The brief says the family's data are Schwartz-class. They are not")
and asserts the repair, "so K(b)-(f) apply at, say, `s = 9`", closing with "Recorded rather
than repaired". But K(a) through K(e) as written are quantified over **every** `s >= 0`:
K(a) concludes `u_0 in H^s` for every `s`, K(c)'s blow-up characterisation is justified by
"(A1) holds for every `s` by (a)", and K(e) states (H\*) as "`u in C([0,T');H^s)` for every
`s >= 0`". None of those is true for this datum above `s = 9.5`.

The substance is very likely fine (every consumer in the chain needs a finite `s`, and `s = 9`
covers them), but the theorem's existence and uniqueness clause currently cites a proposition
whose hypothesis its own datum fails. That is a hypothesis used in the chain and not in the
statement. What is owed is a restated Proposition K with hypothesis "`u_0 in H^{s_0}` for a
fixed `s_0 > 5/2`, divergence free, axisymmetric" and conclusions quantified over `s <= s_0`.

I verified the two tail claims the correction rests on (`s5_results.json`, `tail_check`):
`Theta(u = -3, L = 40) = 5.109079826e-12` against the closed form
`rho^8/(e^2 + rho^8) = 5.109089028e-12`, and the outer tail `e^{-8(u - L + eps_r)}` to the same
place. The `rho^{-8}` / `rho^{8}` statement is right, and `H^s` for `s < 9.5` follows.

---

## 4. FINDING F-4, **MAJOR**: (Gamma-off) is PROVED only under an a priori posture, and the
bootstrap's closing step is not written

Row `D7` reads "`(Gamma-off)` `Gamma(s) <= 2 a(0,s) + C'' M`; `C''` the fixed point; existence
by sign change | `u2` sec.8 + `fix2/f5` + `fix4/p4` | **PROVED** for `L >= L_Gamma^exist`".

`u2`'s theorem carries hypothesis **(H4)**, "the bootstrap posture of `B(t)`:
`Gamma(sigma) := ||grad u(.,sigma)||_{L^inf(R^5)} <= Gamma-bar` for `sigma <= s`", as an
assumption. Section 8 then shows a self-consistent `Gamma-bar` **exists** for `L >= L_*` and
bisects for the existence boundary. Existence of a fixed point is not the same statement as
the a priori assumption propagating: what closes a bootstrap is a continuity / first-crossing
argument (`Gamma` continuous in `s`, `Gamma(0) < Gamma-bar` strictly, so the first crossing
cannot occur). `u2` does not write it; its own sec.12 item 6 only asks a referee to check that
the bisection finds the smallest branch. By the brief's rule, a step marked PROVED whose
closing argument is a sketch is MAJOR. The step is routine. It is not written, and `D7` should
read PROVED-modulo until it is.

---

## 5. FINDING F-5, **MAJOR**: Lemma T' is invoked outside its stated range, and one of its
hypotheses is asserted rather than proved

(a) **Range.** `write/lemma-T-shell-dependent/PROOF.md` line 129 states Lemma T' for
`lambda in C^1([rho_0,R];[1,3/2])`, and its Corollary 2 line 305 defines
`r_h := inf_{lambda in [1,3/2]} P_{h_delta}(lambda)/lambda`. The chain invokes Corollary 2
with `r_h` on `[1, lam_max]`, and the self-consistent window puts `lam_max` well above `3/2`:
my re-run of the proved column at `L = 2e+06` uses `lam_max = 1.9342710419` and
`r_h = 0.8779222518` (`s6_results.json`), and at the certified cap `lam_max = 2.0570755611`.
The source lemma's own parenthetical (lines 138-139) disclaims range-dependence, and its
pivotal identity `int_0^pi sin^2 phi / g^4 dphi = pi/(2 lambda)` is exact for every
`lambda > 0`, so the extension is very probably sound. But no seat states the extended-range
Corollary 2 as a lemma, and `THEOREM_S3_v2.md` row `B3` marks the unextended source PROVED.
The margin is not large: `P_h` stops increasing at `lam_mono = 2.0693`, which is what caps the
window, so the invocation runs to within `0.6 %` of a real obstruction.

(b) **`lambda(.,s) in C^1`.** `u1/PROOF.md` line 363 asserts it bare. `(P2)` supplies only the
derivative bound `|rho lambda'/lambda| <= kappa_s/L`; that `partial frak_a / partial log rho
= -F(2 rho, s)` is continuous in `rho` at each `s`, which is what `C^1` means and what Lemma
T' hypothesis (H1)/(S) consume, is nowhere argued. `F(2 rho, s)` is an `l = 1` spherical
coefficient of the transported vorticity.

---

## 6. FINDING F-6, **MAJOR**: the ENCLOSED list contains grid values, not certified interval
computations

`THEOREM_S3_v2.md` sec.3 lists eight items. Taking them in order, with the brief's test
(certified outward-rounded interval computation, or a grid value called an enclosure):

1. **`r_h` and `P_h` increasing.** The certificate is a genuine Lipschitz-plus-grid argument
   with closed-form primitives, and I confirm it is conservative: the true infimum at `c_*` is
   `0.9035285735` against the certified `0.9033435309`, and the true `min P_h'` is
   `0.3549090746` against the certified `0.3335513465`. It is **not** interval arithmetic, and
   one of its inputs is mislabelled: the "proved Lipschitz bound `|h_sigma''| <=
   225.9466903722`" is computed by `p1_profile.lipschitz_hpp` as `np.max` over the same
   600001-point grid, i.e. it is a scan value, not a proved bound. My scan of the same quantity
   returns `225.94670713`, so the number is right; the **rigorous** a priori bound, from
   `|W_sigma''| <= ||Wtil'||_inf ||chi_sigma'||_1 = 58.19333257 x 39.89422804`, is
   `2415.5503891825`, a factor `10.7` larger. Substituting it inflates the grid pad on
   `sup|h'| = 11.3248852` from `1.1831e-03` to `1.2648e-02`, which moves the Lipschitz
   correction inside `r_h` from `1.850e-04` to `1.852e-04`. **The direction and the size are
   safe, and I supply the replacement.** MAJOR by the rule, repaired here.
   Separately, `FIX4.md` sec.2.1's stated reason for testing the majorant only at `lam_max`
   ("`lam^9` up, every `I_i` down") is a non-sequitur: an increasing factor times decreasing
   factors need not increase. The conclusion is nevertheless true, checked here at 441 values
   of `lambda` on `[1, 2.10]` with zero decreasing steps (`s2_results.json`,
   `majorant_monotone_claim`). MINOR, listed as F-9 below.
2. **`E_0`, `Gfrak_0`, `Psi_sigma(0)`, `N_sigma`.** These are described as "Suprema of
   real-analytic functions over the same grid" with **no** pad at all. They are grid values
   labelled ENCLOSED. I supply the enclosures: with the rigorous bounds
   `|A_sigma''| <= 2445.626` and `|W_sigma''| <= 2321.578` and the grid step `5.236e-06`, the
   interior-maximum pad is `8.381e-09` and `7.956e-09` respectively, giving
   `N_sigma in [1.0124508487, 1.0124508495]` and `E_0 in [7.5523076888, 7.5523076951]`.
   Both contain `fix4`'s values. Note that for `N_sigma` the direction of a grid undershoot is
   the unsafe one, since `N_sigma` is a denominator and the theorem's normalisation is
   `||omega_0||_inf = M`: the true supremum can exceed `M` by at most `7.5e-10` relative.
   MAJOR by the rule, repaired here, no consequence.
3. **`C_E`.** Described as "converged to `7.3e-06` in `lmax` and `4.7e-07` in `du`". A
   convergence statement is not an enclosure. It enters only as `(2/5) log C_E` in the
   dictionary, so `7.3e-06` relative moves `log Lambda_*` by `3e-06` out of `2.8e+06`.
   MINOR in effect. The arithmetic checks: `2 log s + (2/5) log(0.0551906693) = 2.9135781816`
   against the quoted `2.9135781820`, and `2 x 1424610.4953 + 2.9135781820 = 2849223.9042`
   against the quoted `2849223.9041`.
4. **`C''`, `Ghat`, `C_R`, `chat_a`, `L_Gamma^exist`.** The `(Gamma-off)` root's existence is a
   genuine intermediate-value argument. `C_R`'s supremum over the angle is a sampled search
   (`C_R_sup_fast`: a 160-point coarse scan and a 41-point refinement). I checked it against
   the imported `T2.C_R_sup` at `n = 200, 400, 800` and against the fast search at
   `n = 80, 160, 320`: all six agree at `642.71110` to `3e-08` relative, with the maximum at
   `phi = 32.93 deg`, an interior smooth maximum. **No finding.**
5. **`K_2`.** See F-7.
6. **The far/near constants.** Truncated-series quadratures in the source seat, not intervals.
   Their three hypotheses are the load-bearing part and they are genuinely profile-free; I
   re-verified all three independently (sec.8 below). MINOR.
7. **Interior parabolic Schauder**, cited not proved, declared. Correctly listed.
8. **The CKN attribution not used.** Correct as far as it goes, but see F-3: the alternative
   route, Proposition K, has its own unmet hypothesis.

---

## 7. THE PROFILE-CONSISTENCY LEDGER (attack 1)

Every constant in the chain, and how `fix4` treated it. (a) recomputed for the mollified
profile, (b) proved profile-independent, (c) imported without re-derivation.

**(a) Recomputed.** `kappa_delta`, `c_*`, `c_2`, the clock floor, `lam_max`, `r_h`,
`min P_h'`, `lam_mono`, the window cap, `N_sigma`, `E_0`, `Gfrak_0`, `||grad eta_0||_inf`,
`Psi_sigma(0)`, `C_kink`, `Psi_2(0)`, `eps_a`, `C_E` and the `log Re_E` shift,
`||omega_0||_2^2`, the shell density `J`, `K2hat` at four `lambda`, `C''`, `Ghat`, `chat_a`,
`c_G`, `p`, `C_R`, `L_Gamma^exist`, `ELL_RAMP`. All checked above or re-run; all correct.

**(b) Proved profile-independent, and the proof checked.**
`C_K = 3/(8 pi^2)`, `|S^4|`, `C_K|S^4| = 1`, the Riesz constant `pi^4/2`, `3 pi^2/16`, `R_A`:
properties of the kernel and of `R^5`. `fix2` sec.4(c)'s branch-bound structure and its
`YMAX = 6.4072265625` (a statement about `Theta`, not about the angle). Lemma T' itself, whose
only datum hypothesis (D) is `|eta_0| <= M/r` with support in the shell and which needs no
regularity of the profile at all. Theorem V.1, whose `B1` depends on `c, q, n` only.
The far/near constants `0.291999`, `0.014754`, `3.999218/sin phi`, `pi/8`: their derivation
eliminates the Gegenbauer coefficients by Cauchy-Schwarz against Parseval, so the profile
enters only through `int_{-1}^1 w^2 dt <= 2 M^2`, and their source script takes no profile
argument. I re-verified all three hypotheses for the mollified normalised profile with my own
instrument: `sup|A_sigma|/N_sigma = 0.9999999943` (i.e. `1`), z-oddness residual `1.97e-15`,
`int_{-1}^1 w^2 dt = 1.8239293727 <= 2`. Genuine category (b).

**(c) Imported without re-derivation.**
1. **`vartheta(f=4) = 0.0556`**. FATAL, sec.1 above. Never recomputed anywhere: `vartheta`
   appears in the whole `s3close` tree only in `pmax-h3v/` and in the two `THEOREM_S3*.md`
   files.
2. **`K_2 <= 161.7735 M/rho_0`**. Proved in `hk2` for its datum (D-B), which is `tanh` in the
   **angle** as well as in the radius, at `f = 0`, over `lambda in [1,3/2]`
   (`t3_budget.py` line 47 says so in its own comment). `fix4/p3` re-checks the mollified
   datum's `K2hat` at **four** `lambda` values only (`1.0`, `1.5`, `lam_max(c_*)`,
   `lam_max(cap)`) at `L = 10`, and reports the majorant clearing by `3.9889277977`. See F-7.
3. **The exclusion geometry** (the four cones `delta`, `pi/2 +- delta_m`, `pi - delta`) is kept
   as `fix2` had it. For `K_2` that is conservative, as `FIX4.md` sec.3.3 argues. For the
   Theorem V.4' ball it is **not**: it is the mechanism of F-1.
4. **`sigma_z/sigma_y = 3.5130747` and the constant `16.026` in `E_vartheta`** were computed at
   the kinked clock's `c_*`; `c_*` has moved. Not load-bearing only because `E_vartheta` is not
   carried at all (F-2). MINOR.
5. `s = 1/sin delta = 7.66129757554039` inside `viscous_budget` (`nutau = c/(s^2 L)`) is
   **definitional**, since sec.1.1 sets `nu := M rho_0^2/s^2`, and is correctly not moved. But
   note that for the kinked datum `s` and `E_0` coincided and now do not (`7.6613` against
   `7.5523`). `viscous_budget`'s `N = r_0/sin(delta)` is the `||eta_0||_inf r_*/M` of Theorem
   V.4' and, at `7.6613` against the datum's true `7.5523`, is conservative. Safe, and worth a
   line in the sheet.
6. The MEASURED entries of `t3_budget.RECORD` (`Cprime_measured = 11.74`,
   `Ghat_measured = 0.98142`, `material_offset_M = 0.069`, `K2_measured_window = 5.3854`) are
   **not** touched by the "all proved" column: I read `COLUMNS["proved"] = ("u2_fixedpoint",
   "far_near", "proved", "K2_proved_window", 1.0, True)` and confirmed each branch. The only
   record entry the proved column reads is `K2_proved_window`. **The "all proved" column is
   clean of measured constants**, with the single qualification of F-7.

**Stale-constant leak test.** `p4_budget.bind` rebinds `t3_budget`'s `KAPPA`, `CSTAR`, `C2`,
`E0`, `G0`, `SHIFT` but not its `LAM_MAX` and `R_H`, which still hold the kinked datum's
`1.8420112114` and `0.9186365333`. I poisoned both with `NaN` and re-ran the proved column at
`L = 2e+06`: every output is bit-identical (`eps = 0.06691234177450323` both times). No stale
constant reaches the budget; `lam_max = e^{3c/4}` and `r_h` are taken at the self-consistent
`c` everywhere. **No finding.**

---

## 8. THE WINDOW AND THE CAP (attack 2)

`s2_window.py` rebuilds `Q(lambda) = P_h(lambda)/lambda` from the independent profile with
300-point Gauss panels. The probability normalisation `int_0^1 3 lambda^9 v^2 D^{-5/2} dv = 1`
holds to `2.3e-15` at every `lambda` tested, so the measure `mu_lam` is right.

| `lambda_max` | at | `r_h` (true inf) | `min P_h'` (true) | `fix4` certified |
|---|---|---|---|---|
| `1.8558724491` | `c = c_*` | `0.9035285735` | `0.3549090746` | `0.9033435309` / `0.3335513465` |
| `2.0570755611` | certified cap `1.1664585076` | `0.8332958213` | `0.0192685010` | `> 0` |
| `2.0693384642` | monotone cap `1.1760705119` | `0.8284145899` | `-0.0000321223` | boundary |

Both `fix4` certificates are genuine lower bounds and both are conservative. `lam_mono` here is
`2.0693179283`, giving a monotone cap `c/c_* = 1.1760544630` against `fix4`'s `1.1760705119`;
`fix4`'s is `1.4e-05` high, which is inside its own certified cap `1.1664585076` by a wide
margin, so nothing turns on it. The self-consistent window is imposed at one `c` everywhere,
and `lam_max = e^{3c/4}` is evaluated at that `c` (leak test, sec.7). **The window and the cap
check out.**

One presentational defect, MINOR (F-8): `THEOREM_S3_v2.md` sec.1.2 lists
`r_h = inf_{[1,lam_max]} P_h/lambda >= 0.9033435309` and row A3 repeats it, but that is the
value at `lam_max(c_*)`. The proof does not confine itself to `c = c_*`: the theorem quantifies
over `c in [c_*, 1.1664585076 c_*]`, and at `L = 2e+06` the number actually used is
`0.8779222518`, at the cap `0.8332958213`. The constant table in the statement is not the
constant the proof uses.

---

## 9. HYPOTHESES USED BUT NOT IN THE STATEMENT (attack 3)

The statement says: datum of sec.1.1, `nu > 0`, `f >= 4`, `L >= L_Gamma^exist`, and "**No
further hypothesis.**" Walking the graph:

| hypothesis | source | in the statement? | status |
|---|---|---|---|
| `(H3'_vartheta)`, `vartheta in [0,1)` | Theorem V.4', row B4 | no | **FATAL**, and FALSE for this datum (F-1) |
| Prop K's (4) (Schwartz) | row M6 / sec.1.3 | no, and the datum fails it | **MAJOR** (F-3) |
| `(H4)`, u2's a priori `Gamma <= Gamma-bar` | row D7 | no | **MAJOR** (F-4) |
| Lemma T' `lambda in [1,3/2]`; `lambda(.,s) in C^1` | row B3, F3 | no | **MAJOR** (F-5) |
| the contradiction hypothesis `||omega^theta|| <= (3/2)M` on the window | u1 (P2),(P3); u2 (H3); hk2 (D1') | not stated | **MINOR**: legitimate as the contradiction posture, but "No further hypothesis" hides that the theorem is a proof by contradiction and every constant is conditioned on it |
| pmax-h3v (B1) `Gamma-bar, K_2 < infinity`; (B2) `b(0,s)=0`; (B3) `x.b <= Gamma_rad|x|^2` | rows D1, D2, D3 | no | **NONE**: (B2),(B3) are proved in u2; (B1)'s two finiteness clauses follow qualitatively from Prop K on `[0,T']`, and `pmax-h3v` uses `K_2` only qualitatively |
| pmax-h3v (D1)-(D5) on the datum | row D3 | implied by sec.1.1 | **NONE**: `W_sigma/N_sigma` is real-analytic and even at both poles, so (D5) holds; the `rho^{-9}` tail (D4) holds with `Afrak = O(e^{-L})` |
| hk2 (D1') and the transport of the shell density (2.1) over the window | row B6 | no | **MINOR**: `hk2`'s own status table records "transport of the shell structure (2.1) over the window: NOT PROVED, inherited from prove-lagrangian's GAP T". Defused because (B1) needs only `K_2 < infinity`, which Prop K supplies |
| interior parabolic Schauder | row D3 | declared in sec.3 item 7 | **NONE** |
| `eps(L)` monotone decreasing to the floor | sec.1.3 | asserted | **MINOR**: tabulated at five values of `L`, not proved |

**On the `H^s` question the brief raises.** With `s < 9.5` the datum is short of Schwartz but
long past everything the analytic machinery consumes. The P1/P2 maximum principles need only
`eta_0 in W^{1,inf} cap C` and `z`-oddness (P1) plus (D1)-(D5) (P2); `pmax-h3v` line 178 states
"no decay hypothesis and no `L^1` hypothesis is assumed", the decay conditions having become
Lemma F and Corollary F'. Nothing in P1/P2 breaks. The `eps`-regularity escape is not used
(sec.3 item 8). The only casualty is Proposition K's literal hypothesis, F-3.

---

## 10. THE BUDGET (attack 4)

Re-run from a byte copy (`s4_budget_rerun.py`):

```
   L_*  proved, eps <= 0.1664585076 :  1424610.4953     (FIX4: 1424610.4953)
   L_*  proved, eps <= 0.1          :  1555469.4004     (FIX4: 1555469.4004)
   log Lambda_* = 2 L_* + 2.9135781820  =  2849223.9041 (FIX4: 2849223.9041)
   shift  = 2 log(7.66129757554039) + (2/5) log(0.0551906693) = 2.9135781816
```

Exact reproduction. The "all proved" column reads no MEASURED record entry (sec.7 item 6).
`C_E` recomputed from `p2_results.json` is `0.05519066934750944`, and the shift arithmetic
closes.

**F-7, MINOR: the `K_2` import.** `t3_budget.RECORD["K2_proved_window"] = 161.7735` is
commented "PROVED over lambda in [1,3/2]" and was proved for `hk2`'s (D-B) datum at `f = 0`;
`hk2`'s own sec.11 grades it "PROVED (modulo the two rows above)", those rows being a collar
quadrature that is "not interval-certified" and near-ball suprema "computed by dense
sampling". `fix4/p3` re-checks the mollified datum at four `lambda`. So the majorant is
imported from a different datum, is itself partly a sampled bound, and is re-verified only
pointwise. **The direction is safe and the item is not load-bearing**, which I established
rather than took on trust by sweeping `C_K` far past the majorant:

```
   C_K =        1 :  L_*(eps<=0.1) = 1555474.3954
   C_K = 161.7735 :  L_*(eps<=0.1) = 1555474.3954     (identical)
   C_K =     1000 :  L_*(eps<=0.1) = 1555474.3954     (identical)
   C_K =    10000 :  L_*(eps<=0.1) = 1555484.3854     (+6.4e-06 relative)
   C_K =  1000000 :  L_*(eps<=0.1) = 1556593.6799     (+7.2e-04 relative)
```

`K_2` would have to be wrong by four orders of magnitude to move the fourth digit of `L_*`.
`FIX4.md`'s claim that `(H-K2)` has left the budget is correct.

**F-10, MINOR.** The theorem's displayed conclusion carries an all-measured column
(`L_* = 1757.8779`, `log Lambda_* = 3518.6693`). It is labelled, but a measured number inside a
theorem statement invites exactly the misreading the campaign's own ledger keeps recording.

---

## 11. THE BFG CONFRONTATION (attack 5)

BFG, arXiv:1704.05546v4, Theorem 10, at line 582 of
`sharp/literature-short-time/txt/bfg.txt`. Its hypotheses, as stated there: the initial datum
`omega_0` in `L^2 cap L^inf`, and a constant `M` larger than 1. Its conclusion: a constant
`c(M) > 1` and a unique mild solution `omega` in `C_w([0,T], L^inf)` with
`T >= (1/c(M))(1/||omega_0||_inf)`, holomorphic on a complex strip of width proportional to
`sqrt(t/c(M))`, and, in its own words, "`kw(t)kL infinity (Ω t ) ≤ M kω0 k∞`".

**The family satisfies each hypothesis, re-checked here for the mollified profile**
(`s5_bfg_datum.py`, independent of `fix4/p2` and `fix4/p3`):

| clause | requirement | measured here |
|---|---|---|
| `omega_0 in L^inf` | finite | `sup|A_sigma|/N_sigma = 0.9999999943`, i.e. `1` exactly by construction; `||omega_0||_inf = M sup Theta` |
| `omega_0 in L^2` | finite | `||omega_0||_2^2/(M^2 R^3) = 1.4381082621`, finite at every `L` (`fix4`: `1.4381082627`) |
| decay | "suitable" | `O(rho^{-8})` at infinity, `O(rho^8)` at the origin, both verified against the closed forms to `1.8e-06` relative |
| the class | `C_w([0,T], L^inf)` | the maximal strong solution of Prop K lies in it |
| `M > 1` | `M = 3/2` | yes |

**So the confrontation is well posed and, granting the chain, the clause is violated for
`L >= L_*`**: the theorem gives `M T_d <= c_2 (1 + eps)/(2L) -> 0` at fixed `M`, while BFG's
`c(3/2)` is absolute (no dependence on `||omega_0||_2`, on `nu` by the scaling
`u -> lam u(lam x, lam^2 t)`, or on geometry).

**But the chain is not granted.** What stands between the family and the clause is no longer
"the enclosure list and two citations", as sec.4 says. It is:

1. `(H3'_vartheta)`, which is **false** for the datum as stated (F-1). Until the ball is
   shrunk and `vartheta` recomputed on it, Theorem V.4' does not apply and `eps_v` is unproved.
2. Proposition K's hypothesis, which the datum fails as literally written (F-3), so the
   object the two statements are supposed to be about is not yet pinned down. This is the
   clause `THEOREM_S3_v2.md` sec.4 leans on hardest: "the two statements are about the same
   object".
3. The unwritten bootstrap closure in `(Gamma-off)` (F-4) and the unstated extended-range
   Corollary 2 (F-5).
4. Then the enclosures, of which two are grid values (F-6), and the two citations.

The estate's standing rule that the refutation is not asserted until a seat that did not build
the chain has refereed it has now been exercised, and the answer is: not yet.

---

## 12. THE CORRECTED STATEMENT

The theorem of sec.1.3 should read, until the repairs are made:

> ... Assume, on `[0, c/(ML)]`: `L >= L_Gamma^exist`; **and (H3'_vartheta) for the datum of
> sec.1.1 on the ball `B(x_*, d)` with `d = rho_* sin(phi_0 - delta - 4 sigma)`, which holds
> with `vartheta <= 0.1006` measured** (grid lower bound; at `pmax-h3v` sec.B5's
> `SAFETY = 2`, `2 vartheta = 0.2012 < 1`). Then, for every `L >= L_*` ...

with the budget re-run through `pmax-h3v/q4_budget_theta.py` so that
`eps_v = eps_bulk[1 + vt(9 + 2 sigma_z/sigma_y)] + Q(E_hess + E_4 + E_tail)` is what the
columns report, and with row `M6` restated on a finite-`s` Proposition K. On the evidence of
sec.1 the reach does not move: `L_* = 1424610.4953` and `log Lambda_* = 2849223.9041` are
unchanged to the last displayed digit at every ball radius down to half of ASSEMBLY's.

**Overall grade: NOT PROVED as stated. PROVED-MODULO-(H3'_vartheta on a shrunken ball; a
finite-`s` Proposition K; the bootstrap continuity step in (Gamma-off); an extended-range
Corollary 2 and `lambda(.,s) in C^1`; two grid values re-certified; two citations)** once
those five are written, and the numerical content of the sheet survives all five.

## 13. THE REMAINING GAP, IN ONE SENTENCE

`fix4` mollified the angular profile to make `Psi(0)` finite and did not check the one place
where `pmax-h3v` had chosen the kinked profile on purpose, so the theorem now needs its
Theorem V.4' ball moved four mollification widths off the taper cone and `vartheta` recomputed
there, which costs nothing in `L_*` and has not been done.

**FL-000 stands. Nothing here touches the headline problem.**

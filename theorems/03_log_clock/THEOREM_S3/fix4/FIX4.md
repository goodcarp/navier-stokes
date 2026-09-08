# FIX4: the kinked angular profile replaced by its `sigma`-mollification, and every constant that depends on it recomputed

Seat `s3close/round2/THEOREM_S3/fix4`, sitting of 2026-09-08.
Laws: `TORMENT NEXUS/LAWS.md`, first 120 lines, read before any work.
Posture: FL-000 will fall; find the move. Nothing below moves it.

Every number in this sheet came out of a script in **this** folder that I wrote and ran
(`p1_profile.py`, `p2_energy.py`, `p3_kernel_k2.py`, `p4_budget.py`, `p5_epsa.py`,
`p6_terms.py`, `p7_derived.py`), re-asserted by `check_fix4.py`; `SHA256SUMS` is computed with
`shasum -a 256`, never typed. Nothing outside `fix4/` was written, except the one new file
`round2/THEOREM_S3/THEOREM_S3_v2.md` the brief asks for. `THEOREM_S3.md`, its three addenda,
`fix2/`, `fix3/`, `hk2/`, `pmax-h3v/`, `u1/`, `u2/` and the `forced-route-2026-09/` tree are
read-only here.

**Read** (L-14 declaration): `THEOREM_S3.md` sec.1 (the datum), sec.3 (the modulo list M1-M8),
sec.4, sec.5; `ADDENDUM_1_2026-09-08.md`; `ADDENDUM_2_2026-09-08.md`;
`ADDENDUM_3_2026-09-08.md`; `fix2/FIX2.md` in full and `f1_window.py`, `f2_ramp_k2.py`,
`f2_datum_s3.py`, `f4_scans.py`, `f5_gamma_exist.py`; `fix3/FIX3.md` sec.0-sec.2 and
`g1_vstrain.py`; `round2/u2/PROOF.md` sec.8 and sec.9.1; `round2/u1/PROOF.md` sec.3;
`round2/pmax-h3v/PROOF.md` sec.A6 (Proposition 3) and sec.A7;
`rebuild/far-near-kernel-lemma/NOTE.md` sec.1-sec.3;
`forced-route-2026-09/theorems/02_ADDENDUM_uniqueness_2026-09-08.md` sec.3 (Proposition K) and
`02_ADDENDUM_uniqueness_XVENDOR_2026-09-08.md`.

**Code imported** (byte copies in `imported/`, hashes in `SHA256SUMS`; all nineteen verified
equal to the entry in their source seat's own `SHA256SUMS` by `check_fix4.py`): from `fix2/`,
`f1_window.py`, `f2_datum_s3.py`, `f2_ramp_k2.py`, `f4_scans.py`, `f5_gamma_exist.py`,
`f6_controls.py`, `f1_results.json`, `f4_results.json`; from `fix2/imported/`,
`t1_results.json`, `t2_gamma_CR.py`, `t3_budget.py`, `hk2lib.py`, `k1_algebra.py`,
`k2_datum.py`, `k3_bound.py`, `k4_direct.py`; from `fix3/`, `g1_vstrain.py`, `g1_results.json`;
from `THEOREM_S3/`, `t1_datum.py`.
L-14's inverse convention, declared. The objects under **test** are imported unchanged, because
the point is to run the same object on a new datum: `t3_budget`'s viscous half
(`viscous_budget`, `gauss_tail_aniso`, `J_pow`, `eps_Tprime`, the `RECORD` block) and its
`bootstrap`; `t2_gamma_CR`'s `chat_a` and `C_R_of_phi`; `f5_gamma_exist`'s direct
`(Gamma-off)` root; `hk2`'s `k3_bound.lambda_bulk / assemble / sups_on_ball_B`; `fix2`'s
`DatumS3` and `StrainedS3` field objects. Everything I am testing a **claim** about (the
angular profile and every constant it carries, `kappa_delta`, `P_h`, `r_h`, `E_0`, `Gfrak_0`,
`Psi_sigma(0)`, `C_E`, the shell density, the far/near hypotheses, the self-consistent window)
is re-implemented here from the mathematics.

`t3_budget`'s module-level datum constants (`KAPPA`, `CSTAR`, `C2`, `E0`, `G0`, `SHIFT`) are
rebound to the mollified datum's before any call. That rebinding **is** the substitution under
test; it is stated here so it is not mistaken for an edit of the imported file, which is
untouched.

---

## 0. STATUS TABLE

| unit | item | status | what it says |
|---|---|---|---|
| 1 | the mollified datum, stated exactly | **PROVED** | `W_sigma := (W even-reflected) * chi_sigma`, then divided by `N_sigma = sup|W_sigma sin phi|`. The normalisation is **not** optional and was missing from `fix3`'s table. |
| 2 | `kappa_delta`, `P_h`, `r_h`, `P_h` increasing, the window cap | **ENCLOSED** (certificates) | `kappa_delta = 0.4917868801`, `c_* = 0.8244732108`, `c_2 = 1.6489464217`, floor `1.6700567e-02`; `r_h >= 0.9033435309` on `[1, 1.8558724485]`; `P_h' >= 0.3335513465 > 0` there; the certified window cap is `eps <= 0.1664585076`, down from `fix2`'s `0.1943662079`. |
| 3 | `E_0`, `Gfrak_0`, `\|grad eta_0\|_inf` | **ENCLOSED** (certified branch bound) | `E_0 = 7.5523076942`, `Gfrak_0 = 36.0795702396` (against `7.6612975755` and `58.6954805410` kinked), `\|grad eta_0\|_inf = 12.6002066238`. |
| 4 | `Psi_sigma(0)`, `C_kink`, `eps_a` (**M3**) | **DISCHARGED** | `Psi_sigma(0) = 489.0279393839` in units `M/rho_0^2`, `C_kink = sigma Psi_sigma(0) = 9.7805587877`; `eps_a = 9.3574e-07` at `L = 1e+07` in the proved column, `3.9169e-05` of `eps`. For the kinked datum `Psi(0) = +infinity`. |
| 5 | what does **not** move | **PROVED** | the exact kernel constants; the far/near constants (their three hypotheses are re-verified for the mollified profile); `K_2` (the imported majorant `161.7735` still dominates, by `3.9889277977`); the `ELL_RAMP` count. |
| 6 | `C_E` and the `log Re_E` shift | **PROVED** | `C_E = 0.0551906693`, shift `2.9135781820` against the kinked `2.9234858948`, a move of `-9.9077128e-03`. |
| 7 | `(Gamma-off)`, `C_R`, `C''(L)` | **PROVED** given M1 | `C''(1e+09) = 1002.5338`, `Ghat = 792.0570`, `C_R = 546.8521`, `L_Gamma^exist = 9298.1782` at `c = c_*`. |
| 8 | the budget | **COMPUTED**, conditional on the modulo list | `L_* = 1424610.4953` (all proved, `eps <= 0.1664585076`), `1555469.4004` (all proved, `eps <= 0.1`); `log Lambda_* = 2849223.9041` and `3110941.7143`. |

**The headline.** Against `ADDENDUM_3`'s kinked-datum table, `L_*` in the proved column falls
by a factor `1.3251248106` at each sheet's own certified cap (`1887786.7127` to
`1424610.4953`) and by `1.2711978138` at `eps <= 0.1` (`1977309.3012` to `1555469.4004`), and
`log Lambda_*` falls by the same `1.3251244816`. The all-measured column moves the other way,
by `1.1426886044` and `1.2029353604`, because its `L_*` sits near the clock floor and the floor
rose.

**FL-000 stands. Nothing here touches the headline problem.**

---

## 1. THE DATUM

### 1.1 The definition, exactly

Write `phi_ax := min(phi, pi - phi)`, `delta = 7.5 deg`, `delta_m = 5.0 deg`,

```
    A(phi) = sgn(cos phi) min(1, phi_ax/delta) min(1, |phi - pi/2|/delta_m) ,
    W(phi) = A(phi)/sin phi .
```

`A` is odd across each pole and `sin phi` is odd there too, so `W` is **even** about `phi = 0`
and about `phi = pi`; that is the geometrically correct extension, and it is the one
`pmax-h3v` Proposition 3 and `fix3` sec.1.6 use. Let `Wtil` be `W` extended by even reflection
at both poles (`2 pi` periodic) and let

```
    chi_sigma(y) = exp(-y^2/(2 sigma^2))/(sigma sqrt(2 pi)) ,
    W_sigma := (Wtil * chi_sigma)|_{[0,pi]} ,     A_sigma := W_sigma(phi) sin phi ,
    N_sigma := sup_{[0,pi]} |A_sigma| .
```

> **THE MOLLIFIED DATUM (fix4).** With `Theta`, `rho_0`, `nu`, `R`, `L`, `eps_r = 0.25` exactly
> as in `THEOREM_S3` sec.1.1,
> ```
>     omega_0^theta(rho, phi)  :=  - M Theta(rho) A_sigma(phi) / N_sigma ,
>     eta_0                     =  omega_0^theta / r
>                               =  - (M/rho) Theta(u) W_sigma(phi)/N_sigma ,      u = log(rho/rho_0),
> ```
> axisymmetric, no swirl, odd in `z`, at `sigma = 0.02` rad (`1.146 deg`).

`W_sigma` is **real-analytic** on `[0, pi]` (a Gaussian convolution of a bounded function is
real-analytic), and even at both poles, so `W_sigma'(0) = W_sigma'(pi) = 0` and `eta_0` is
`C^infinity` in the angle including the axis. The four kink cones `phi = delta`,
`pi/2 - delta_m`, `pi/2 + delta_m`, `pi - delta` are gone.

The tracked point `x_* = ((1+f) rho_0, phi_0)`, `phi_0 = 30 deg`, sits `19.6349540849 sigma`
from the nearest of those four angles, so the Gaussian weight reaching it from any of them is
`e^{-193}`: the profile at the tracked point is unchanged by the smoothing. What does change
there is the normalisation: `A_sigma(phi_0)/N_sigma = 0.9890911210`, i.e. the tracked point
carries `1.0908879024 %` less than the nominal amplitude, against `0 %` for the kinked datum.

### 1.2 Why the normalisation is not optional (a correction to `fix3`'s table)

On the bulk `delta < phi < pi/2 - delta_m` the kinked `W` is exactly `1/sin phi`, which is
**convex**. Convolution with a probability kernel therefore raises it: `W_sigma >= W` there by
Jensen, so `A_sigma = W_sigma sin phi >= 1` on the bulk. **Mollifying `W` overshoots the
amplitude.** Measured, at `sigma = 0.02`:

```
     N_sigma = 1.0124508488 ,   attained at  phi = 9.9885000509 deg ,
     i.e. the unnormalised profile exceeds M by 1.2450848776 % .
```

The theorem's entire normalisation is `||omega_0^theta||_inf = M`: `T_d` is the first time
`||omega||_inf` reaches `(3/2) M`, and the contradiction hypothesis is `lam_om = M_s/M <= 3/2`.
So the overshoot has to be divided out, and `kappa_delta = (1/2) P_h(1)` is computed for
`h = |A_sigma|/N_sigma`, not for `|A_sigma|`.

`fix3` sec.1.6's table reports the **unnormalised** profile. Its
`kappa_delta(0.05) = 0.4983883 > 0.4978224` and the "free lunch, recorded not claimed" that
`ADDENDUM_3` sec.1 M5 draws from it are the overshoot, not a gain. Normalised, `kappa_delta`
goes the other way at every `sigma`:

| `sigma` | `N_sigma` | `kappa_delta` unnormalised (`fix3`) | `kappa_delta` normalised (here) | floor `(1-2k)/(2k)` |
|---|---|---|---|---|
| kinked | `1` | `0.4978224383` | `0.4978224383` | `4.3741734e-03` |
| `0.002` | `1.0002090379` | `0.4978232887` | `0.4977192705` | `4.5823613e-03` |
| `0.005` | `1.0011591563` | `0.4978278606` | `0.4972514885` | `5.5274073e-03` |
| `0.01` | `1.0039591496` | `0.4978442141` | `0.4958809666` | `8.3064962e-03` |
| **`0.02`** | **`1.0124508488`** | **`0.4979100241`** | **`0.4917868801`** | **`1.6700567e-02`** |
| `0.05` | `1.0484593277` | `0.4983882748` | `0.4753530044` | `5.1849879e-02` |

The unnormalised column reproduces `fix3` `g1`'s own numbers to `1.7e-11` or better at every
`sigma`, so the two instruments agree and the disagreement is entirely about which profile the
theorem's `M` normalises.

**This is why `sigma = 0.05` is not the right choice**, and the brief's default is refused
(L-17). `sigma = 0.05` costs `5.18e-02` of clock floor against a target of `0.1`, leaves a
certified monotone window of only `eps <= 0.0853498988`, and lands at
`L_* = 3238884.2136936598` in the coarse scan of sec.6.3, `1.6362795943` times **worse** than
the kinked datum. The scan's minimum over `{0.05, 0.02, 0.01, 0.005, 0.002}` is at
`sigma = 0.02`, and that is the value this sheet carries.

---

## 2. THE DATUM CONSTANTS

### 2.1 `kappa_delta`, `P_h`, `r_h`, and the window cap

`fix2` sec.4's exact structure is unchanged and is used unchanged: with
`D(v, lam) = 1 + (lam^6 - 1)v^2`, `F_lam(v) = lam^9 v^3 D^{-3/2}` is a probability distribution
function on `[0,1]`, `Q(lam) := P_h(lam)/lam = E_{mu_lam}[H]` with `H(v) = h(arcsin v)`, and

```
     Q'(lam) = -9 lam^8 int_0^1 h'(arcsin v) v^3 sqrt(1-v^2) D^{-5/2} dv .
```

What changes is `h`: `h_sigma = |A_sigma|/N_sigma` is smooth, so `fix2`'s closed-form
`Abar`/`Bbar` split of `H'` into a taper piece and an equatorial piece no longer applies. The
certificate here is the cell-wise version of the same idea, with the same closed-form primitive

```
     int_0^V v^3 (1 + K v^2)^{-5/2} dv = K^{-2}[ 2/3 + (1/3) S^{-3/2} - S^{-1/2} ] ,
     S = 1 + K V^2 ,  K = lam^6 - 1  (exact; the K -> 0 limit is V^4/4) ,
```

and the weight `v^3 sqrt(1-v^2) D^{-5/2} <= v^3 D^{-5/2}` majorised cell by cell. The two
certified statements are

```
     r_h  = inf_{[1, lam_max]} Q  >=  0.9033435309       (grid minimum 0.9035285738,
                                                          max Lipschitz bound 25.4923791809)
     min_{[1, lam_max]} P_h'      >=  0.3335513465  > 0  (the value at lam_max is 0.3549367307)
```

at `lam_max = e^{3 c_*/4} = 1.8558724485`. Because the majorant `9 lam^9 sum_i max(M_i,0) I_i`
is increasing in `lam` (`lam^9` up, every `I_i` down), testing it at `lam_max` certifies the
whole interval.

`P_h'` vanishes at `lam_mono = 2.0693384635`, so the window may not exceed
`c/c_* = 1.1760705119`; the largest window for which the **certificate** above still closes is

```
     c/c_*  <=  1.1664585076 ,      i.e.   eps  <=  0.1664585076 ,
```

against `fix2`'s `0.1943662079` for the kinked profile. **The cap moves, and it moves down**,
by a factor `0.8637588763`: the mollified `h` is slightly flatter near the axis, so `Q` falls
away a little sooner. `lam_max` at that cap is `2.0570755604`.

Status: **ENCLOSED**, with a certificate whose only non-exact input is the evaluation of the
real-analytic `h_sigma` on a `600001`-point grid; the grid pad added to `sup |h_sigma'|` is
`1.1831e-03`, from the proved Lipschitz bound `|h_sigma''| <= 225.9466903722`.

Derived: `c_* = log(3/2)/kappa_delta = 0.8244732108`, `c_2 = 2 c_* = 1.6489464217`,
`1 - 2 kappa_delta = 1.6426239743e-02`, floor `= 1.6700567265e-02`.

### 2.2 `E_0`, `Gfrak_0` and `||grad eta_0||_inf`

`fix2` sec.4(c)'s branch bound is profile-independent and is used verbatim: with `P = W^2`,
`S = W'^2` and `Theta = sig1 - sig2`,

```
     (rho^2 |grad eta_0|/M)^2  <=  max( P + S            if 81 P <= 8 S ,
                                        6.4072265625 P + S  otherwise )        pointwise in phi,
```
maximised over `phi`. That is a rigorous upper bound over **all** `u` at once. Its value and
the direct `(phi, u)` scan agree here to the last digit, so the enclosure is tight:

```
     E_0     = sup rho|eta_0|/M      =  7.5523076942     (kinked: 7.6612975755 = 1/sin delta)
     Gfrak_0 = sup rho^2|grad eta_0|/M = 36.0795702396   (kinked: 58.6954805410 = 1/sin^2 delta)
                                                         attained at phi = 9.0474000515 deg
     ||grad eta_0||_inf rho_0^2/M    = 12.6002066238     (kinked: 20.3344667336)
```

`Gfrak_0` falls by a factor `1.6268342260`. That is the whole reason the sheet is worth
running: `Gfrak_0` is the binding datum constant of `Ghat`, `C''`, `C_R` and hence `L_*`
(`THEOREM_S3` finding 1).

The kinked column reproduces `fix2` sec.4(c)'s closed forms `1/sin delta` and `1/sin^2 delta`
to `7.0e-09` and `1.4e-08`, and `||grad eta_0||_inf = 20.3344667336` is a **third** independent
instrument agreeing with `fix3`'s `20.3339` and `20.3261` and disagreeing with
`pmax-h3v` sec.A7's `31.093556949918312`. That discrepancy stands as `fix3` recorded it;
nothing in the chain reads it.

### 2.3 `Psi_sigma(0)`, `C_kink`, and `Psi_2(0)`

With `eta_0 = G(rho) W(phi)`, `rho^3 Lap_5 eta_0/M = -W(Theta_uu + Theta_u - 2 Theta)
- Theta(W'' + 3 cot(phi) W')` (`fix3` sec.1.6), and `Psi(0) = sup rho|Lap_5 eta_0|`:

```
     Psi_sigma(0) = 489.0279393839   (units M/rho_0^2),  attained at u = 0.39 ,
     C_kink := sigma Psi_sigma(0) = 9.7805587877 ,
     Psi_2(0) := sup rho |grad^2 eta_0|_F = 306.4888158429 ,  attained at u = 0.40 .
```

For the **kinked** datum `Psi(0) = +infinity`: `A` is piecewise linear, so `A''` carries a Dirac
mass on each of the four kink cones. Its absolutely continuous part is `171.4566296577`,
reproducing `fix3`'s `171.4257` to `1.80e-04`. Unnormalised, this sheet's `Psi_sigma(0)`
reproduces `fix3`'s whole `sigma` table to `1.6e-04` or better.

---

## 3. WHAT DOES NOT MOVE, CHECKED RATHER THAN ASSERTED

### 3.1 The exact kernel constants

`K(w) = -d_z G_5(w)`, `C_K = 3/(8 pi^2) = 0.0379954439`, `|S^4| = 8 pi^2/3 = 26.3189450696`,
`C_K |S^4| = 1` exactly; the Riesz composition `int |x-x'|^{-4}|x'|^{-2} dx' = (pi^4/2)/|x|` in
`R^5`, `pi^4/2 = 48.7045455170`, confirmed by an independent polar quadrature at
`48.7126360810` (relative `1.6612e-04`); `3 pi^2/16 = C_K pi^4/2 = 1.8505508252`;
`2 R_A = 3.9992184446`. These are properties of the kernel and of the dimension. No angular
profile enters.

### 3.2 The far/near constants: the hypotheses re-verified

`far-near-kernel-lemma` sec.2-3 gives `C1 = 0.29199853` (far), `C2in = 0.01475427` (inner) and
`|a_collar| <= (3.999218/sin phi + pi/8) M`. Their hypotheses are exactly three, and none of
them mentions the taper: `|omega^theta| <= M` pointwise, `omega^theta` odd in `z` (only odd `l`
survive), and the Parseval bound `sum_l g_l^2 N_l = int_{-1}^1 w(t)^2 dt <= 2 M^2`. For the
mollified, normalised profile:

```
     sup |A_sigma|/N_sigma          =  1.0                    (by construction, exactly)
     max |A_sigma(pi-phi) + A_sigma(phi)|/N_sigma = 5.995e-15  (z-oddness, to machine precision)
     int_{-1}^{1} w(t)^2 dt         =  1.8239293735   <=  2    (kinked: 1.8751740891)
```

So the four constants are unchanged, and they are unchanged for a reason that has been checked
rather than assumed.

### 3.3 The shell density and `K_2`

`hk2` Lemma 5.1's `J`, with the ramp kept and the supremum over `rho'` taken as (2.1) requires:

| datum | `J` supremum over `rho'` | `J` plateau |
|---|---|---|
| `THEOREM_S3`'s kinked datum (control) | `96.2564876464` | `75.2494579308` |
| **the mollified datum** | **`93.6629034481`** | **`73.3641888444`** |

The control reproduces `fix2`'s `96.249844` and `75.251104` to `6.9e-05` and `2.2e-05`. The
mollified datum's shell density is **lower**, by `2.7 %`.

`K_2` is then recomputed with `hk2`'s own instrument (`k3_bound.lambda_bulk`, `assemble`,
`sups_on_ball_B`, all unchanged) on `fix2`'s `DatumS3` field object with only its angular factor
`W(phi, k)` replaced by the mollified profile, at `f = 4` (ADDENDUM_1), `phi_0 = 30 deg`,
`L = 10`, over the full corrected `lambda` range:

```
  K2hat = K_2 rho_0/M ,  f = 4
  ------------------------------------------------------------------------------------
   lambda            1.0       1.5     1.8558725 (lam_max)   2.0570756 (lam at the cap)
   mollified      18.9649   25.0428     33.5336                40.5556
   kinked (ctl)   19.1685   25.4176     33.7096 (at 1.8420112)  42.0083 (at 2.0741)
  ------------------------------------------------------------------------------------
```

The kinked control reproduces `fix2`'s `19.169 / 25.418 / 33.710 / 42.008` to `2.4e-05` or
better at every `lambda`. **The budget's imported majorant `161.7735` still dominates**, over
the whole corrected range, by `3.9889277977`. The exclusion geometry is kept as `fix2` had it
(the same four cones, now layer centres rather than kinks), which is conservative: a smooth
datum admits a larger ball, and a larger ball can only lower `K_2` through `sups_on_ball_B`.

### 3.4 The `ELL_RAMP` count

`int_{-inf}^{inf} Theta du = L - 2 eps_r = L - 0.5` exactly, confirmed at `L = 12`
(`11.500000000000005`) and `L = 40` (`39.5`); each tail is `log(1 + e^{-2})/8 = 0.0158660014`.
So the `tanh` ramp costs `0.5` e-folds and `t3_budget.py:65`'s `ELL_RAMP = 0.25` undercharges by
`0.25`. This sheet charges `0.5`. `ADDENDUM_3` owed that correction to the budget and it is now
paid.

Its size, measured directly rather than through a bisection: at `L = 2e+06` and the
self-consistent window `c = 1.0669124225 c_*`, moving `ELL_RAMP` from `0.25` to `0.5` moves
`eps` from `0.0669122018` to `0.0669123418`, a difference of `1.3995047e-07`. That is below the
`L_*` bisection's own resolution, which is why the sensitivity block of `p4` returns the same
`L_*` for both. Recorded as a null result with its size, not as an absence.

`C_K` is the same kind of null: at the same point, `C_K = 1` gives `eps = 0.0669123357` and
`C_K = 161.7735` gives `0.0669123418`, a difference of `6.0847127e-09`. `(H-K2)` has left the
budget entirely.

---

## 4. THE ENERGY CONSTANT AND THE `log Re_E` SHIFT

`THEOREM_S3` sec.1.4's dictionary is unchanged:
`log Re_E = 2L + 2 log s + (2/5) log C_E`, `C_E = ||u_0||_2^2/(M^2 R^5)`, with
`C_E = 2 pi sum_{l odd} 2 N_l Hhat_l^2 D_l/(2l+3)`. `Theta` does not change, so every `D_l` is
unchanged and the whole move is in the Gegenbauer coefficients of the angular profile.

| datum | `C_E` | shift `2 log s + (2/5) log C_E` |
|---|---|---|
| bang-bang cap (control) | `0.1724039760` (record `0.172403978`) | |
| sharp radial edges (control) | `0.1703256330` (`ASSEMBLY`: `0.17032563295410327`) | `3.3643454569` |
| `tanh` ramp, kinked angle (control) | `0.0565747735` (`THEOREM_S3`: `0.05657477`) | `2.9234858948` |
| **the mollified datum** | **`0.0551906693`** | **`2.9135781820`** |

The shift moves by `-9.9077128e-03`. Two independent controls on the coefficients (L-98):
Gegenbauer orthogonality `int_{-1}^1 (1-t^2)[C_l^{3/2}]^2 dt = N_l` reproduced to `3.1e-11`;
Parseval `sum_l Hhat_l^2 N_l = int_{-1}^1 g_odd^2 dt` reproduced to `4.8409e-09`; and the
synthesis `sum_l Hhat_l C_l^{3/2}(t) = -W_sigma(arccos t)/N_sigma` to `4.690e-04` on
`|t| < 0.98` at `LMAX = 121`. Convergence: `4.7e-07` in `du`, `7.3e-06` in `lmax`.

The `L^2` size of the datum, which BFG's hypothesis needs finite:
`||omega_0||_2^2/(M^2 R^3) = 1.4381082627` for the mollified datum, `1.4785130340` for the
kinked one (reproducing `THEOREM_S3`'s `1.4785137` to `4.5e-07`).

---

## 5. `(Gamma-off)`, `C_R` AND `C''`

`u2` sec.8's closure and `t2_gamma_CR`'s implementation are unchanged; only `E_0` and `Gfrak_0`
are the mollified datum's, and the root is `fix2` `f5`'s **direct** root rather than the Picard
iteration. Proved column, frozen window `c = c_*`:

| `L` | `C''` | `c_G` | `Ghat` | `C_R` (sup over all `phi`) | `chat_a` sup | `chat_a(phi_0)` |
|---|---|---|---|---|---|---|
| `1e+04` | `1829.2356` | `1.38753` | `1584.8911` | `1072.2304` | `121.4223` | `13.046833` |
| `3e+04` | `1125.0971` | `1.26763` | `908.0867` | `623.7327` | `107.7552` | `13.046833` |
| `1e+05` | `1034.4930` | `1.24524` | `822.2342` | `566.8468` | `105.3794` | `13.046833` |
| `3e+05` | `1012.8147` | `1.23949` | `801.7581` | `553.2798` | `104.7783` | `13.046833` |
| `1e+06` | `1005.5792` | `1.23754` | `794.9300` | `548.7556` | `104.5746` | `13.046833` |
| `1e+09` | `1002.5338` | `1.23671` | `792.0570` | `546.8521` | `104.4884` | `13.046833` |

Against `THEOREM_S3` sec.4.3's kinked-datum row at `1e+09` (`C'' = 1460.8194`,
`Ghat = 1250.4810`, `C_R = 850.620`), every drive constant falls by about the factor
`Gfrak_0` falls by. `chat_a(phi_0) = 13.046833` is unchanged, as it must be: `chat_a` reads
`E_0` only through the collar branch `2 R_A e^{c_G} E_0`, and at `phi_0` the other branch binds.

`(Gamma-off)`'s existence boundary, by `fix2` `f5`'s sign-change certificate:

```
     L_Gamma^exist  =  9298.1782   at c = c_* ,      16231.8518  at the certified cap ,
     L_Gamma^picard =  9298.4988   at c = c_* .
```

Both are far below `L_*`, so M8's threshold clause is not binding. The kinked datum's boundary
was `14259.948`; the smaller `Gfrak_0` flattens `F` and moves it down.

---

## 6. THE BUDGET

### 6.1 `eps(L)`, self-consistent window, `f >= 4`

The window is `fix2`'s: `c >= c_*(1 + eps(c,L))`, with `lam_max = e^{3c/4}`, `r_h`, `ell_loss`,
the majorant ODE, `c_G` and the whole viscous budget re-read at that `c`, and `c/c_*` capped at
the certified monotone window `1.1664585076`.

| column | `eps(1e+03)` | `eps(1e+04)` | `eps(1e+05)` | `eps(1e+06)` | `eps(1e+07)` |
|---|---|---|---|---|---|
| all proved | `inf` | `inf` | `inf` | `inf` | `0.0238897079` |
| proved, `C_a` measured | `inf` | `inf` | `inf` | `inf` | `0.0235580850` |
| all measured | `inf` | `0.0267974804` | `0.0176414180` | `0.0167941841` | `0.0167100481` |

`inf` means no `c` inside the cap satisfies the self-consistency condition at that `L`, so the
theorem says nothing there. The proved column first closes at `L_* = 1424610.4953`.

### 6.2 The term table, where the proved column is alive

| column | `L` | `c/c_*` | `mu` | `mu_J` | `eps_T'` | `eps_v` | `eps_a` | `eps` | `C_R` | `C''` | `r_h` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| proved | `2e+06` | `1.0669124` | `4.264720e-03` | `1.00734e-04` | `0.04704700` | `3.12955e-09` | `0.04706270` | `0.0669123` | `642.7111` | `1165.2701` | `0.877922` |
| proved | `1e+07` | `1.0238897` | `6.592570e-04` | `1.83300e-05` | `7.01814e-03` | `2.31754e-10` | `7.02126e-03` | `0.0238896` | `579.0127` | `1057.3298` | `0.894421` |
| proved, `C_a` meas. | `1e+07` | `1.0235581` | `6.337820e-04` | `1.71914e-07` | `6.69911e-03` | `2.31599e-10` | `6.69959e-03` | `0.0235580` | `557.2293` | `1056.5523` | `0.894541` |
| measured | `1e+04` | `1.0267975` | `5.190590e-03` | `1.73212e-04` | `9.34088e-03` | `2.72186e-06` | `9.82669e-03` | `0.0267974` | `2.6192` | `11.7400` | `0.893364` |
| measured | `1e+06` | `1.0167941` | `4.984540e-05` | `1.69697e-06` | `8.70858e-05` | `1.80006e-09` | `9.19150e-05` | `0.0167940` | `2.6192` | `11.7400` | `0.896965` |
| measured | `1e+07` | `1.0167100` | `4.982840e-06` | `1.69667e-07` | `8.70338e-06` | `1.58064e-10` | `9.18628e-06` | `0.0167099` | `2.6192` | `11.7400` | `0.896995` |

The `eps` column is the budget's own value at that window, which sits a unit in the last
displayed digit below the window's `c/c_* - 1` because the self-consistency condition is met
with slack at the bisection's resolution. The binding term is `eps_T'` through `mu`, in every
column, exactly as in `fix2`; `eps_v` and `eps_ell` are invisible; `eps` is `eps_T'` plus the
floor.

### 6.3 The choice of `sigma`

A coarse scan (`n = 10`, `nbis = 10`, `tol = 1e-3`, `nphi_CR = 120`; used only to choose
`sigma`, never quoted as a constant), proved column:

| `sigma` | `L_*` (`eps <= 0.1`) | against the kinked datum | `L_*` (`eps <= own cap`) | against the kinked datum |
|---|---|---|---|---|
| kinked | `1979419.79` | `1` | `1890360.14` | `1` |
| `0.05` | `3238884.21` | `1.6362795943` | `3238884.21` | `1.7133688675` |
| **`0.02`** | **`1555714.17`** | **`0.7859445383`** | **`1424715.72`** | **`0.7536742276`** |
| `0.01` | `1599810.83` | `0.8082221047` | `1509106.58` | `0.7983169689` |
| `0.005` | `1712780.24` | `0.8652940857` | `1630347.53` | `0.8624533996` |
| `0.002` | `1838254.88` | `0.9286836919` | `1755546.63` | `0.9286836919` |

The trade is two-sided: larger `sigma` lowers `Gfrak_0` (good, linearly in the drive) and raises
the clock floor and lowers the certified window cap (bad, and the floor enters the denominator
`target - floor`). The minimum on this grid is at `sigma = 0.02`, and the full-precision run
below uses it.

### 6.4 `L_*` and `log Lambda_*`, three columns

`log Lambda_* = 2 L_* + 2.9135781820`.

| column | `L_*` (`eps <= 0.1664585076`) | `log Lambda_*` | `L_*` (`eps <= 0.1`) | `log Lambda_*` |
|---|---|---|---|---|
| **all proved** | **`1424610.4953`** | **`2849223.9041`** | **`1555469.4004`** | **`3110941.7143`** |
| proved, `C_a` measured | `1380034.8075` | `2760072.5285` | `1500454.5140` | `3000911.9416` |
| all measured | `1757.8779` | `3518.6693` | `2046.8908` | `4096.6951` |

Against `ADDENDUM_3`'s table (kinked datum, `eps <= 0.1943662` and `eps <= 0.1`):
proved `1887786.7127 / 1977309.3012`, `C_a` measured `1848687.4 / 1932180.6`, measured
`1538.37 / 1701.58`. The proved column improves by `1.3251248106` and `1.2711978138`; the
`C_a`-measured column by `1.3395947624` and `1.2877302057`; the all-measured column gets worse
by `1.1426886044` and `1.2029353604`, because it sits close to the floor and the floor rose by
`3.8179938670`.

### 6.5 Control against `fix2` (L-11, L-14)

The same instrument, driven with the **kinked** datum, `fix2`'s cap `0.1943662079` and
`fix2`'s `ELL_RAMP = 0.25`:

```
     L_*(eps <= 0.1943662)  =  1889578.7372     against fix2's 1887786.7127   (+0.0949273 %)
     L_*(eps <= 0.1)        =  1978766.7792     against fix2's 1977309.3012   (+0.0737102 %)
```

Both are high by less than a tenth of a percent, in the conservative direction, and the cause is
identified: this sheet's `r_h` certificate is coarser than `fix2`'s closed-form one, returning
`0.9178761522` at `c = c_*` against `fix2`'s `0.9186364112`, a relative deficit of `8.276e-04`,
which enters `eps_T'` as `1/r_h`. `kappa_delta` agrees with `fix2`'s to `4.0e-08`.

---

## 7. THE MODULO LIST AFTER THIS SHEET

**M3 (V-strain): DISCHARGED.**
*What it needed.* `fix3` sec.1 proved the reduction in full: Lemma 1 (the strain kernel is
`Kcal = d_z G_5 = -3z/(8 pi^2 rho^5)` and its weight is exactly `3/8` per e-fold), Lemma 2 (the
difference `w = eta - eta^tr` obeys `D_t w = nu Lap_5 eta`, `w(0) = 0`, so
`sup rho|w(s)| <= nu e^{c_R} int_0^s Psi`), Lemma 3 (the harmonic identity), the propagation
lemma, and the theorem `eps_a = (3/8)(1 + 1/(4L)) Psibar c e^{c_R}/(s^2 L)`. The single open
input was `(H-Delta-eta)`, whose only content is `Psibar = sup rho|Lap_5 eta| < infinity`. For
`THEOREM_S3`'s kinked datum that is **false**: `Psi(0) = +infinity`.
*What supplies it now.* `Psi_sigma(0) = 489.0279393839` at `sigma = 0.02` (sec.2.3). With
`e^{c_R} <= e^{c_G}` (conservative, as `fix3` `g1` sec.D does) and the self-consistent window:

| `L` | column | `c_G` | `eps_a` | `eps` | `eps_a/eps` |
|---|---|---|---|---|---|
| `1e+04` | measured | `1.27084` | `9.4266e-04` | `0.0267975` | `3.5177e-02` |
| `1e+05` | measured | `1.25863` | `9.2288e-05` | `0.0176414` | `5.2314e-03` |
| `1e+06` | measured | `1.25749` | `9.2107e-06` | `0.0167941` | `5.4845e-04` |
| `1e+07` | measured | `1.25738` | `9.2089e-07` | `0.0167100` | `5.5110e-05` |
| `1e+07` | **proved** | `1.26634` | **`9.3574e-07`** | `0.0238897` | **`3.9169e-05`** |

The propagation half is not binding either: `tau K_2 sqrt(G_inf N_2) = 0.2843828760` at
`L = 1e+04` and `2.8438288e-03` at `L = 1e+06`, against `Psi_2(0) = 306.4888158429`, i.e.
`9.2787e-04` and `9.2787e-06` of the datum term.

**M5 (H5-Lip): DISCHARGED.**
*What it needed.* `u2`'s (H5) asks for `Gfrak_0 < infinity` and (H1) for a bounded classical
solution; `THEOREM_S3`'s datum has kinks at four cones, so `Gfrak_0 = 58.6954805410` was an
**essential** supremum and `eta_0` was not `C^1` in the angle. `ADDENDUM_1` discharged the
P-max half of M5 by `pmax-h3v` Proposition 3, which uses the same `sigma`-family as an
**approximating sequence** for the linear problem; the rest of the chain still faced a datum
that is only Lipschitz.
*What supplies it now.* The datum **is** one member of that family. `W_sigma` is real-analytic
in the angle and even at both poles, so `eta_0` is `C^infinity` in the angle, `Gfrak_0` is an
attained supremum, and there are no kink cones for any later step to avoid.
*What it costs, stated.* `kappa_delta` falls by `1.2123917551 %`, `c_2` rises by a factor
`1.0122727109`, and the clock floor rises from `4.3741734e-03` to `1.6700567e-02`. What is
bought is `Gfrak_0` down by `1.6268342260` and `Psi_sigma(0)` finite, and the net on `L_*` is
the `1.2711978138` of sec.6.4. **M3 and M5 close together, exactly as `fix3` said they would.**

**M6 (H\*): DISCHARGED for the constructed family, with one correction to the brief's premise.**
Proposition K of `forced-route-2026-09/theorems/02_ADDENDUM_uniqueness_2026-09-08.md` sec.3
gives, for a smooth axisymmetric divergence-free datum and a force satisfying (5) (here `f = 0`,
which satisfies it trivially): a unique maximal strong solution `u in C([0,T'];H^s) cap
L^2(0,T';H^{s+1})` for every `T' < T_*` (K(b)); `T_*` characterised by
`sup_{t<T_*}||u(t)||_{L^inf} = infinity` (K(c)); axisymmetry and swirl-freeness preserved
(K(d)); and (H\*) **as literally written**, with no escape of the `H^s` norm to spatial infinity
on `[0,T']`, because membership in `C([0,T'];H^s)` is how `T_*` is defined (K(e)). Its
uniqueness corollary C1 is PROVED there, and the whole note carries a second-vendor blind
confirmation (`02_ADDENDUM_uniqueness_XVENDOR_2026-09-08.md`: items (1), (2) and (3) all
returned CONFIRMED, item (1) with an independent nested-ball proof; the residue is the
attribution of the forced epsilon-regularity criterion to CKN 1982, which neither line verified
at source).

*The correction.* The brief says "the family's data are Schwartz-class". They are not, and the
statement does not need them to be. `omega_0` decays like `rho^{-8}`, not faster than every
polynomial, so it is not Schwartz; and at the origin `Theta(rho) = rho^8/(rho^8 + e^2)` exactly
(because `eps_r = 1/4` makes `2/eps_r = 8`), so `eta_0 = -M Theta W_sigma/(N_sigma rho)`
behaves like `|x|^7` and `omega_0` like `|x|^8`: the datum is `C^infinity` away from the origin
and `C^{7,1}` at it, and `omega_0 in H^s` for `s < 9.5`, hence `u_0 in H^s` for `s < 10.5`.
That is what Proposition K's local theory consumes (it needs `s > 5/2` and propagation upward),
so K(b)-(f) apply at, say, `s = 9`, and (H\*) holds. Recorded rather than repaired; the repair
is cheap and unexercised: replacing the inner `tanh` edge by any profile vanishing to infinite
order at the origin makes `omega_0` genuinely `C^infinity` on `R^3` and changes `Theta` only
inside the inner tail, whose entire contribution to every constant is bounded by
`ADDENDUM_3` sec.3's tail table (`<= 1.7e-04` relative, or exactly zero).

**M4 (lam-cap): SLACK, unchanged.** (P3) proves `lam <= e^{3c/4}` from the contradiction
hypothesis alone, and the proved column uses it. A self-consistent argument should give
`(3/2)(1 + O(1/L))`, which would improve `lam_max` from `1.8558724485` and therefore `r_h`,
`mu` and `L_*`. It is a slack in the constants, not a gap in the validity, and it is not
written.

**M8: the threshold clause of the statement, unchanged in kind, smaller in size.**
`L >= L_Gamma^exist = 9298.1782` at `c = c_*` and `16231.8518` at the certified cap. Not
binding: `L_* = 1424610.4953`.

**M1 (P-max) and M2 (H3-V)** were discharged in `ADDENDUM_1` and are untouched here.
**M7 (Theta-tail)** was discharged in `ADDENDUM_3` and is untouched: the tails are a property of
`Theta`, which this sheet does not change.

**Carried, unchanged:** `refute-u2` MAJOR-1 (the P1/P2 "sharpness" tautology, about `14x` of
`C''` available but unproved) and MAJOR-2 (the datum choice). ADDENDUM_1's `f >= 4` constraint
stands and every number above is computed under it.

---

## 8. THE DECORRELATION TABLE

Every control this sheet runs, with the number it had to hit and the number it returned.

| control | reference | here | relative |
|---|---|---|---|
| `kappa_delta`, kinked | `fix2` `f4` `0.4978224182` | `0.4978224383` | `4.04e-08` |
| `E_0`, kinked | closed form `1/sin delta` | `7.5523076942` at `sigma`, `7.6612975222` kinked | `6.96e-09` |
| `Gfrak_0`, kinked | closed form `1/sin^2 delta` | `58.6954797236` | `1.39e-08` |
| `lam_mono`, kinked | `fix2` `2.0769162` | `2.0770339680` | `5.67e-05` |
| `Psi(0)` a.c. part, kinked | `fix3` `171.4257` | `171.4566296577` | `1.80e-04` |
| `Gfrak_0(sigma)` unnormalised | `fix3`'s five rows | five rows | `<= 6.24e-05` |
| `Psi_sigma(0)` unnormalised | `fix3`'s five rows | five rows | `<= 1.59e-04` |
| `kappa(sigma)` unnormalised | `fix3`'s five rows | five rows | `<= 1.23e-10` |
| `C_E` sharp edges | `ASSEMBLY` `0.17032563295410327` | `0.1703256330` | `1.63e-16` |
| `C_E` `tanh` ramp, kinked | `THEOREM_S3` `0.05657477` | `0.0565747735` | `6.13e-08` |
| `log Re_E` shift, kinked | `THEOREM_S3` `2.9234859` | `2.9234858948` | `1.79e-09` |
| `\|\|omega_0\|\|_2^2/(M^2R^3)`, kinked | `THEOREM_S3` `1.4785137` | `1.4785130340` | `4.50e-07` |
| `J` supremum, kinked | `fix2` `96.249844` | `96.2564876464` | `6.90e-05` |
| `J` plateau, kinked | `fix2` `75.251104` | `75.2494579308` | `2.19e-05` |
| `K2hat`, kinked, `f = 4` | `fix2`'s four `lambda` | four `lambda` | `<= 2.36e-05` |
| Riesz `pi^4/2` | `48.7045455170` | quadrature `48.7126360810` | `1.66e-04` |
| Gegenbauer orthogonality | `N_l` | three `l` | `<= 3.1e-11` |
| Parseval, mollified | `int g_odd^2 dt` | `sum Hhat_l^2 N_l` | `4.84e-09` |
| `L_*` proved, kinked, `fix2` settings | `fix2` `1887786.7127` | `1889578.7372` | `9.49e-04` |
| `||grad eta_0||_inf`, kinked | `fix3` `20.3339`, `20.3261` | `20.3344667336` | third instrument, agrees |

---

## 9. WHAT THIS SHEET CANNOT CATCH

It changes one input, the angular profile, and recomputes everything downstream of it. It does
not check that `C_R` is the right constant, that `u2`'s P1 and P2 are the right instruments,
that the modulo list is complete, or that the clock argument is right. The `r_h` and `P_h`
certificates are enclosures whose only non-exact input is the evaluation of a real-analytic
function on a fixed grid; a finer grid can only move them in the fourth decimal. `K_2` and the
measured column are measurements, and numerics falsify but never prove.

One thing it does catch that no earlier sheet did: **`fix3`'s `sigma` table is unnormalised**,
and `ADDENDUM_3` sec.1 M5's "both moves are favourable" reads it as a free lunch. Normalised, the
`kappa_delta` move is unfavourable at every `sigma`, and at `sigma = 0.05` it is unfavourable
enough to lose more than the `Gfrak_0` gain wins. The net at `sigma = 0.02` is still a gain, and
that is the sheet's result.

---

## GATE, AND FILES

```
251 CHECKS, 251 PASS, 0 FAIL
```

`check_fix4.py` re-reads every displayed constant from the stored JSONs, checks that the string
actually appears in `FIX4.md` or `THEOREM_S3_v2.md`, re-derives the ratios and factors quoted
between them rather than copying, and verifies the nineteen byte copies in `imported/` against
their source seats' own `SHA256SUMS`.

| file | what it establishes |
|---|---|
| `p1_profile.py` / `p1_results.json` | the mollified profile and `N_sigma`; `kappa_delta`, `P_h`, the certified `r_h`, `P_h` increasing, `lam_mono`, the window cap; `E_0`, `Gfrak_0`, `\|grad eta_0\|_inf`, `Psi_sigma(0)`, `C_kink`; the kinked controls and the unnormalised comparison with `fix3` |
| `p2_energy.py` / `p2_results.json` | `C_E` and the `log Re_E` shift for the mollified datum; the bang-bang, sharp-edge and kinked controls; Gegenbauer orthogonality, Parseval and synthesis controls; the BFG `L^2` check |
| `p3_kernel_k2.py` / `p3_results.json` | the exact kernel constants and the Riesz control; the three far/near hypotheses re-verified; `J`; `K_2` on the mollified datum over the corrected `lambda` range with the kinked control; the `ELL_RAMP` count |
| `p4_budget.py` / `p4_results.json`, `p4_scan.json` | the self-consistent window with the mollified datum: `eps(L)`, `L_*`, `log Lambda_*`, three columns; the constants against `L`; `L_Gamma^exist`; the `sigma` scan; the kinked control against `fix2` |
| `p5_epsa.py` / `p5_results.json` | `eps_a` (M3) at the self-consistent window, `Psi_2(0)`, the propagation term |
| `p6_terms.py` / `p6_results.json` | the term table where the proved column is alive, and the direct `ELL_RAMP` and `C_K` sensitivities the bisection cannot resolve |
| `p7_derived.py` / `p7_results.json` | every ratio, percentage and difference quoted in the prose |
| `imported/` | byte copies, with hashes, of every file imported from another seat |
| `check_fix4.py` | the gate |
| `SHA256SUMS` | computed with `shasum -a 256`, never typed |

**FL-000 stands. Nothing here touches the headline problem.**

# FIX5: the six findings of the round-2 referee report, worked

Seat `s3close/round2/THEOREM_S3/fix5`, sitting of 2026-09-08.
Laws: `TORMENT NEXUS/LAWS.md`, first 120 lines, read before any work; seat header
`TORMENT NEXUS/TEMPLATES/SEAT_HEADER_2026-09-08.md`.
Posture: FL-000 will fall; find the move. Nothing below moves it.

Every number in this sheet came out of a script in **this** folder that I wrote and ran
(`fx_profile.py`, `x1_vartheta.py`, `x2_budget.py`, `x3_propK.py`, `x4_bootstrap.py`,
`x5_lemmaT.py`, `x6_certified.py`), re-asserted by `check_fix5.py`; `SHA256SUMS` is computed
with `shasum -a 256`, never typed. Nothing outside `fix5/` was written except the one new file
`round2/THEOREM_S3/THEOREM_S3_v3.md` the brief asks for. `THEOREM_S3.md`, its three addenda,
`THEOREM_S3_v2.md`, `fix2/`, `fix3/`, `fix4/`, `hk2/`, `pmax-h3v/`, `u1/`, `u2/`,
`refute-THEOREM_S3_v2/`, `write/lemma-T-shell-dependent/` and the `forced-route-2026-09/` tree
are read-only here and are untouched.

**Read** (L-14 declaration): `round2/refute-THEOREM_S3_v2/NOTE.md` in full, with its six scripts
and six result files; `round2/THEOREM_S3/THEOREM_S3_v2.md`; `fix4/FIX4.md` with `p1_profile.py`
and `p4_budget.py`; `round2/pmax-h3v/PROOF.md` Part B (B1 to B5) and `q3_theta.py`,
`q4_budget_theta.py`; `round2/u2/PROOF.md` sec.1, sec.8, sec.10 to sec.12; `round2/u1/PROOF.md`
sec.1 to sec.5; `write/lemma-T-shell-dependent/PROOF.md` sec.1 to sec.2 with Corollaries 1 to 3;
`forced-route-2026-09/theorems/02_ADDENDUM_uniqueness_2026-09-08.md` sec.3 and sec.4;
`forced-route-2026-09/theorems/01_two_fences.md` sec.1.1 to sec.1.4; `hk2/PROOF.md` sec.0.

**Code copied** (byte copies in `copies/`, hashes in `SHA256SUMS`): the referee's six scripts and
six result files; `fix4/p1_profile.py`, `p2_energy.py`, `p3_kernel_k2.py`, `p4_budget.py`,
`p5_epsa.py`, `p6_terms.py`, `p7_derived.py` and their JSONs; `fix4/imported/` entire (nineteen
files, including `t3_budget.py`, `t2_gamma_CR.py`, `f5_gamma_exist.py`, `t1_results.json`);
`pmax-h3v/q3_theta.py`, `q4_budget_theta.py`, `q3_results.json`.
L-14's inverse convention, declared. The objects under **test** are imported unchanged, because
the point is to run the same object on a corrected input: `fix4`'s window search
(`p4_budget.Datum`, `bind`, `column_constants`, `assemble_c`, `self_consistent`, `Lstar`),
`t3_budget`'s `bootstrap`, `J_pow`, `eps_Tprime`, `gauss_tail_aniso` and `RECORD`,
`t2_gamma_CR`'s `chat_a` and `C_R_of_phi`, `f5_gamma_exist`'s direct `(Gamma-off)` root, and
`p1_profile`'s `certified_r_h` and `Ph_prime_lower`. Everything I am testing a **claim** about
(the profile and its derivatives to order four, `vartheta`, the viscous budget's `eps_v`, the
certified suprema, `Proposition K`, the bootstrap closure, Lemma T's range) is re-implemented
here from the mathematics.

**FL-000 stands. Nothing here touches the headline problem.**

---

## 0. STATUS TABLE

| unit | finding | grade in the report | status here |
|---|---|---|---|
| 1 | **F-1** `(H3'_vartheta)` is false for the mollified datum on ASSEMBLY's ball | FATAL | **CLOSED**. The ball is moved four mollification widths off the taper cone, `d = rho_* sin(phi_0 - delta - 4 sigma)`, and `vartheta = 0.1305648864` is measured there by the five-parameter search; at `SAFETY = 2`, `2 vartheta = 0.2611297728 < 1`. |
| 2 | **F-2** the budget runs at `vartheta = 0`, i.e. on V.4 and not V.4' | MAJOR | **CLOSED**. `eps_v = eps_bulk[1 + vt(9 + 2 sigma_z/sigma_y)] + Q(E_hess + E_4 + E_tail)` is carried, with the reduced `d` propagated into `E_hess`, `E_4` and `E_tail`. `L_*` in the proved column does not move. |
| 3 | **F-3** Proposition K's hypothesis is not satisfied by the datum | MAJOR | **CLOSED**. Proposition K' is stated at finite Sobolev index and proved from the same two sources; every consumer in the chain is listed with the highest derivative it uses, the largest is five, and `s < 10.5` supplies eight. |
| 4 | **F-4** `(Gamma-off)` is proved only under an a priori posture; the bootstrap's closing step is not written | MAJOR | **CLOSED**. The closed set, the strict margin at `t = 0`, the propagation with a strict margin, the three continuity inputs and the connectedness step are written, with every constant computed. |
| 5 | **F-5** Lemma T' is invoked outside its stated range; `lambda(.,s) in C^1` is asserted | MAJOR | **CLOSED**. Lemma T'' is stated on `[1, lam_max]`, every step of the proof is checked to use no upper bound on `lambda`, and the `C^1` claim is proved from the chain's own P1 and P2 constants. |
| 6 | **F-6** the ENCLOSED list contains grid values | MAJOR | **CLOSED**, and it costs something: the certified window cap moves down from `0.1664585076` to `0.1647266909`, and `fix4`'s `Gfrak_0` is found to carry a discretisation error of its own. `L_*` moves by four parts in ten thousand, in the favourable direction. |

**The headline.** All six findings close. `L_*` in the proved column is unchanged to the last
displayed digit at `1424610.4953` (`eps <= 0.1664585076`) and `1555469.4004` (`eps <= 0.1`) when
the repair is carried on `fix4`'s own constants, and moves to `1424207.9721` and `1554565.5693`
when the certified constants of unit 6 are substituted. The one number that moves materially is
the certified monotone window cap, and it moves down.

---

## 1. UNIT F-1: `(H3'_vartheta)` ON A BALL THAT MISSES THE SMOOTHING LAYER

### 1.1 What went wrong, restated

Dependency row `B4` of `THEOREM_S3_v2.md` imports

```
    vartheta(f = 4) = 0.055612 ,    d = min(f, rho_* sin(phi_0 - delta),
                                            rho_* sin(pi/2 - delta_m - phi_0)) rho_0
```

from `pmax-h3v` Part B. That number is a statement about the **kinked** angular profile.
`pmax-h3v` sec.B1 chose `min(1, .)` angular mollifiers precisely so that on the ball the angular
factor is exactly `1` and the whole deviation `eta_0 - eta_P` is the radial `tanh` ramp. At
`f >= 4` the binding branch of `d` is `rho_* sin(phi_0 - delta)`, so the ball is exactly tangent
to the cone `phi = delta`. `fix4` replaced the kinked angle by its Gaussian mollification at
`sigma = 0.02`, which puts a smoothing layer of width a few `sigma` **inside** that ball, and in
that layer the angular derivatives are `O(sigma^{-(k-1)})`.

### 1.2 The instrument, and its two controls

`x1_vartheta.py` implements `pmax-h3v` sec.B3's five-parameter search for the mollified datum:
three parameters for the point (a radius in `(0,d]` and a direction on `S^2`, writing a ball
point as `x_* + (w1, w2, 0, 0, w5)` by the `SO(4)` symmetry) and two for the unit direction
`v = a Yhat + b e_perp + c e_z`. Along `x + s v` everything is composed by truncated Taylor
arithmetic to order four:

```
    e := eta_0 - eta_P  =  1/r  -  Theta(u) W_sigma(phi)/(rho N_sigma)      (M = rho_0 = 1, z > 0)
```

with `r(s)`, `z(s)`, `rho(s)`, `u(s) = log rho`, `phi(s) = atan(r/z)` all series in `s`. The
angular derivatives `W_sigma^{(k)}(phi)`, `k = 0..4`, come from `fx_profile.py`, a third
instrument for the profile written here: it evaluates

```
    W_sigma^{(k)}(phi) = int Wtil(t) chi_sigma^{(k)}(phi - t) dt
```

in the integration variable `t`, so the breakpoints of `Wtil` sit at fixed abscissae independent
of `phi` and one panel structure serves every `phi`, with 24-point Gauss-Legendre on panels of
width at most `sigma/2`. Its controls are in `fx_controls.json`: against the referee's
adaptive-quadrature instrument it reproduces `W_sigma`, `W_sigma'`, `W_sigma''` and
`W_sigma'''` at ten angles to `2.0272e-10` relative in the worst case; its `W_sigma''''` agrees with a
fourth-order finite difference of its own `W_sigma'''` to `2.8988e-09`, and with the referee's spline
table at `phi = delta` to `3.8052e-09`.

**Control 1 (L-14).** The same search on the **kinked** datum at ASSEMBLY's own `d`, `f = 4`:

```
    vartheta = 0.0552690     against pmax-h3v's 0.055612      relative 6.168e-03, binding k = 4
```

Below `pmax-h3v`'s value, which is what a coarser grid lower bound on a supremum must be. The
instrument is the right instrument.

**Control 2.** The same search on the **mollified** datum at ASSEMBLY's own `d`, `f = 4`:

```
    vartheta = 8.78252       against the referee's 8.76183 ,  binding k = 4 at phi = 9.6147 deg
```

The refutation reproduces, and slightly strengthens.

### 1.3 The corrected ball, and `vartheta` on it

```
    d  :=  rho_* min( sin(phi_0 - delta - 4 sigma) ,
                      sin(pi/2 - delta_m - phi_0 - 4 sigma) ,  f/rho_* ) ,     sigma = 0.02 .
```

The taper branch binds at every `f >= 4`, so `d = rho_* sin(phi_0 - delta - 4 sigma)` and the
smallest polar angle attained on the ball is exactly `delta + 4 sigma`, four mollification widths
clear of the taper cone.

| `f` | `d` (corrected) | `d` (ASSEMBLY) | `R_- = r_* - d` | `rho_min` on the ball |
|---|---|---|---|---|
| **4** | **`1.5381397413`** | `1.9134171618` | `0.9618602587` | `3.4618602587` |
| 8 | `2.7686515343` | `3.4441508913` | `1.7313484657` | `6.2313484657` |
| 16 | `5.2296751203` | `6.5056183502` | `3.2703248797` | `11.7703248797` |
| 32 | `10.1517222923` | `12.6285532680` | `6.3482777077` | `22.8482777077` |

Per-order values at `f = 4`, from the fine search with six rounds of local polish:

```
    vartheta_0 = 0.0118175206      vartheta_1 = 0.0153696032      vartheta_2 = 0.0422876051
    vartheta_3 = 0.0706575208      vartheta_4 = 0.1305648864  <-- binding
    vartheta   = 0.1305648864      at phi = 12.0892249189 deg = delta + 4.005 sigma
```

and at the other admissible `f`:

| `f` | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| `vartheta` | `0.1305648864` | `0.1311569011` | `0.1311567764` | `0.1311567756` |
| `2 vartheta` | `0.2611297728` | `0.2623138022` | `0.2623135528` | `0.2623135512` |
| admissible | yes | yes | yes | yes |

**The corrected `(H3'_vartheta)` instance.** With `x_* = ((1+f) rho_0, phi_0)`, `phi_0 = 30 deg`,
`f >= 4`, `sigma = 0.02` and `d = rho_* sin(phi_0 - delta - 4 sigma)`,

```
   sup_{|v| = 1} | d_v^k ( eta_0 - eta_P )(x) |  <=  vartheta k! M / r(x)^{k+1} ,
   k = 0, 1, 2, 3, 4 ,  for every x in B(x_*, d) ,   eta_P := -M/r ,
   vartheta = 0.2623138022  (the largest of the four values above, at SAFETY = 2) < 1 ,
   Q = (1 + vartheta)/(1 - vartheta) = 1.7111799108 ,  vt = vartheta/(1 - vartheta) = 0.3555899554 ,
```

and at the budget's own minimiser `f = 4`, `vartheta = 0.2611297728`, `Q = 1.7068352824`,
`vt = 0.3534176412`. Theorem V.4' applies, so `M2 (H3-V)` is discharged for this datum and row
`F4` has its warrant back.

### 1.4 How far off the cone the ball has to stand, and how thin that is

`x1`'s offset scan at `f = 4`, `SAFETY = 2`:

| offset | `d` | `vartheta` | `2 vartheta` | admissible |
|---|---|---|---|---|
| `0` | `1.9134171618` | `8.78252` | `17.56504` | no |
| `1 sigma` | `1.8206526970` | `8.78248` | `17.56495` | no |
| `2 sigma` | `1.7271599953` | `8.57131` | `17.14263` | no |
| `3 sigma` | `1.6329764527` | `2.80969` | `5.619373` | no |
| `3.5 sigma` | `1.5856373782` | `0.883225` | `1.76645` | no |
| **`4 sigma`** | **`1.5381397413`** | **`0.131171`** | **`0.2623412`** | **yes** |
| `5 sigma` | `1.4426877944` | `0.102464` | `0.2049283` | yes |
| `6 sigma` | `1.3466587917` | `0.0884151` | `0.1768303` | yes |
| `8 sigma` | `1.1530234754` | `0.0602701` | `0.1205402` | yes |

The transition between `3.5 sigma` and `4 sigma` is a factor `6.8` in half a mollification
width. That is the shape of the cliff the mollification created, and it is worth saying plainly:
at `4 sigma` the hypothesis holds with a factor `3.8` of room against `vartheta < 1`, but the
room in the **geometry** is about a quarter of a `sigma`. The `6 sigma` ball is a strictly safer
statement at no cost in `L_*` (sec.2.4), and a referee who does not like standing this close to a
cliff should read the theorem at `6 sigma`.

**Grid stability.** Coarse settings give `0.1311706039` against the fine `0.1305648864`, a
relative change of `4.639e-03`; the polish improves the coarse grid maximum by at most a factor
`1.0004927432` at any order. `vartheta` is a grid lower bound on a supremum, inflated by
`SAFETY = 2`; that is `pmax-h3v` sec.B5's own convention, and it is not an interval enclosure.
It is listed as such in sec.7.

---

## 2. UNIT F-2: THE BUDGET WITH `E_vartheta` AND `Q` CARRIED

### 2.1 What was wrong

`pmax-h3v` Corollary V.5' is

```
   eps_v(vartheta) = eps_bulk [ 1 + vt (9 + 2 sigma_z/sigma_y) ] + Q ( E_hess + E_4 + E_tail ) .
```

`t3_budget.viscous_budget`, which `fix4/p4` imports unchanged, returns
`eps_v = eps_bulk + E_hess + E_4 + E_tail`: the `vartheta = 0` line, which is Theorem V.4 and not
V.4'. Row `B4` claimed V.4' while row `F4` quoted V.4.

### 2.2 The instrument

`x2_budget.py`'s `viscous_budget_theta` is `t3_budget.viscous_budget` transcribed with exactly
two changes: the ball radius `d` is exposed, so unit 1's repair can be propagated into
`E_hess`'s `r_*/R_-^2` and `4 N/d` terms, into `E_4`'s `r_*/R_-^5`, and into `E_tail`'s three
Gaussian exit probabilities; and `eps_v` is assembled by Corollary V.5'. At `d` equal to
ASSEMBLY's and `vartheta = 0` the transcription reproduces the imported function to
**`0.000e+00`** at four `(L, c, f)` points. Everything else in the budget is imported and
untouched, and the control reproduces `fix4`'s own reach:

```
   L_* proved, eps <= cap :  1424610.4953   (FIX4: 1424610.4953)
   L_* proved, eps <= 0.1 :  1555469.4004   (FIX4: 1555469.4004)
```

### 2.3 `eps_v` and its share

Proved column, `f = 4`, `vartheta = 0.2611297728` at `SAFETY = 2`, `d = 1.5381397413`:

| `L` | window | `eps_bulk` | `E_vartheta` | `E_hess` | `E_4` | `E_tail` | `eps_v` | `eps_v` in `fix4` | ratio | `eps` | `eps_v/eps` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `1e+05` | `c = c_*` | `1.5397e-08` | `8.8766e-08` | `7.5932e-07` | `1.0500e-11` | `0` | `1.400217e-06` | `6.995716e-07` | `2.0015` | `5.8353147148` | `2.3996e-07` |
| `1e+06` | `c = c_*` | `1.5397e-09` | `8.8766e-09` | `7.4555e-09` | `1.0500e-13` | `0` | `2.314178e-08` | `8.257398e-09` | `2.8026` | `0.0840729869` | `2.7526e-07` |
| `2e+06` | `c = c_*` | `7.6985e-10` | `4.4383e-09` | `1.8620e-09` | `2.6251e-14` | `0` | `8.386373e-09` | `2.447619e-09` | `3.4263` | `0.0487599461` | `1.7199e-07` |
| `2e+06` | self-consistent, `c/c_* = 1.0669124225` | `8.0238e-10` | `4.6252e-09` | `2.4028e-09` | `2.8452e-14` | `0` | `1.003684e-08` | `3.129548e-09` | `3.2071` | `0.0669123599` | `1.5000e-07` |

At `L = 1e+05` and `1e+06` the proved column has no self-consistent window (`eps` there is the
budget's value at the frozen `c = c_*`, and the self-consistency condition `c >= c_*(1 + eps)` is
not met), so those two rows are reported at the frozen window and labelled. At `L = 2e+06` both
are reported.

**The reading.** Carrying `vartheta` and shrinking the ball costs a factor between `2.0` and
`3.4` on `eps_v`, and `eps_v` is between `1.5e-07` and `2.8e-07` of `eps`. The two effects pull
in opposite directions and neither is visible: `Q = 1.7068` and the new `E_vartheta` term push
`eps_v` up, while the smaller `d` raises `R_-` from `0.5865828382` to `0.9618602587` and pushes
`E_hess`'s and `E_4`'s `r_*/R_-^2` and `r_*/R_-^5` down hard. `E_tail` is `0` at every `L` the
theorem uses, at both radii, because `d/2` exceeds the anisotropic Gaussian's standard deviation
by many orders of magnitude at every `L >= 1e+05`.

### 2.4 `L_*`, three columns, with the repair carried

`log Lambda_* = 2 L_* + 2.9135781820`. Constants are `fix4`'s, so that the only change against
`FIX4.md` sec.6.4 is unit 1 and unit 2.

| column | `L_*` (`eps <= 0.1664585076`) | `log Lambda_*` | `L_*` (`eps <= 0.1`) | `log Lambda_*` |
|---|---|---|---|---|
| **all proved** | **`1424610.4953`** | **`2849223.9041`** | **`1555469.4004`** | **`3110941.7143`** |
| proved, `C_a` measured | `1380035.9154` | `2760074.7444` | `1500454.5140` | `3000911.9416` |
| all measured | `1761.3926` | `3525.6987` | `2051.2352` | `4105.3840` |

Against `FIX4.md`: the proved column does not move at all; the `C_a`-measured column moves by
`1.1079` in `L_*`, eight parts in ten million; the measured column moves by `3.5147` and
`4.3444`, two parts in a thousand, because it sits close to the clock floor where `eps_v`'s
share is largest.

**The `6 sigma` ball.** Repeating the whole run with the ball at `d = rho_* sin(phi_0 - delta -
6 sigma) = 1.3466591247` and `vartheta = 0.0884151` gives `L_* = 1424610.4953` and
`1555469.4004`, identical to the last displayed digit. The safer geometry is free.

---

## 3. UNIT F-3: PROPOSITION K AT FINITE SOBOLEV INDEX

### 3.1 What was wrong

`THEOREM_S3_v2.md` sec.1.3 cites Proposition K of
`forced-route-2026-09/theorems/02_ADDENDUM_uniqueness_2026-09-08.md` sec.3. Its hypothesis is
"`u_0` smooth, divergence free and axisymmetric satisfying (4)", and (4) is Fefferman's decay
condition quantified over every multi-index and every `K`; smooth plus (4) is the Schwartz
condition, and K(a) says so. The datum is not Schwartz: `omega_0` decays like `rho^{-8}`, and at
the origin `Theta(rho) = rho^8/(rho^8 + e^2)` exactly, so `omega_0` is `C^{7,1}` there. K(a),
K(c) and K(e) are quantified over every `s >= 0` and are false for this datum above `s = 9.5`.
`fix4` sec.7 M6 records the defect and marks it "recorded rather than repaired".

### 3.2 The datum's exact regularity, checked

`x3_propK.py`:

1. **The inner edge is exactly `rho^8/(rho^8 + e^2)`.** With `eps_r = 1/4`, `2/eps_r = 8`, and
   `(1/2)(tanh(4u - 1) + 1) = e^{8u}/(e^{8u} + e^2)`; the residual of that identity in exact
   symbolic arithmetic is `0`, and the two sides agree numerically to `1.12e-12` at eight radii
   between `0.3` and `20`. (Below `rho = 0.3` the double-precision `tanh` saturates at `-1` and
   `Theta` underflows to zero while the closed form is still representable; the identity is
   exact and is checked symbolically there rather than numerically.)
2. **The structure at the origin.** In Cartesian components,
   `omega_0 = -(M/N_sigma) [Theta(rho)/rho] W_sigma(phi) (-x_2, x_1, 0)`, and
   `Theta(rho)/rho = rho^7/e^2 (1 + O(rho^8))`, so the leading term of every component is
   `c |x|^7 x_j W_sigma(phi)`, homogeneous of degree `8` and smooth off the origin, and the next
   term is homogeneous of degree `16`.
3. **`W_sigma` is smooth on the sphere.** It is even at both poles, so every odd derivative
   vanishes there: measured, `W_sigma'(0) = 2.547e-09`, `W_sigma'''(0) = 3.380e-09`,
   `W_sigma'(pi) = 2.547e-09`, against `W_sigma(0) = 7.64` and `|W_sigma'''| ~ 1e+04` a few
   degrees away. The angular factor contributes no singularity at the origin or on the axis.
4. **`|x|^7 x_1` is `C^{7,1}` and not `C^8`.** In exact arithmetic,
   `d^8/dx_1^8 (|x|^7 x_1)` is homogeneous of degree `0` and takes the three distinct values
   `40320` along `e_1`, `0` along `e_3`, and `903105 sqrt(2)/32` along `(1,1,0)/sqrt 2`; so it is
   bounded and discontinuous at the origin, and `d^9` is homogeneous of degree `-1`, with
   `t d^9` having a finite nonzero limit along each ray. A compactly supported function whose
   ninth derivatives are `O(|x|^{-1})` lies in `H^s_loc(R^3)` exactly for `s < 9 + 1/2`.
5. **The Fourier exponent.** The Fourier transform of a degree-`8` homogeneous function on `R^3`
   that is smooth off the origin is homogeneous of degree `-11`, and for the `l = 1` angular mode
   the coefficient is `2^{a+n} pi^{n/2} Gamma((a+n+l)/2)/Gamma((l-a)/2)` at `(n,l,a) = (3,1,8)`,
   whose denominator `Gamma(-7/2)` has no pole, so the constant is `5066760.6` and not zero.
   Hence `int (1+|xi|^2)^s |omega_0hat|^2 dxi ~ int k^{2s-20} dk`, which converges exactly for
   `s < 9.5`. Confirmed on the exact model `f(x) = |x|^7 x_1 e^{-|x|^2}`, whose `l = 1` radial
   transform has the closed form
   `G(k) = sqrt(pi/(2k)) [Gamma(6)/(2 Gamma(5/2))] (k/2)^{3/2} 1F1(6; 5/2; -k^2/4)`; evaluated in
   60-digit arithmetic the successive decay exponents are

   ```
      -11.04729   -11.01172   -11.00292   -11.00073   -11.00018      at k = 50 .. 1600 .
   ```

   (A compactly supported `C^infinity` cutoff was tried first and is useless here: its own
   transform decays like `exp(-c sqrt k)`, which dominates `k^{-11}` over every `k` a quadrature
   reaches. Recorded rather than hidden.)

So `omega_0 in H^s(R^3)` for every `s < 9.5` and for no larger `s`, and `u_0 in H^s` for every
`s < 10.5`.

### 3.3 Proposition K'

> **PROPOSITION K'.** Let `u_0` be divergence free and axisymmetric on `R^3` with
> `u_0 in H^{s}(R^3)` for every `s < S_1`, `S_1 > 5/2` (here `S_1 = 10.5`), and let `f` be smooth
> and axisymmetric satisfying Fefferman's (5) with `div f = 0` (here `f = 0`, which satisfies (5)
> trivially). Then:
>
> **(a')** `u_0 in H^s` for every `s in (5/2, S_1)`, and `f in L^1(0,inf;H^s) cap L^2(0,inf;H^s)`
> for every `s`, so hypothesis (A1) of `[FENCES]` Lemma A holds at every index.
> **(b')** There is `T_* in (0,inf]` and a unique `u` with
> `u in C([0,T'];H^s) cap L^2(0,T';H^{s+1})` for every `T' < T_*` and every `s in (5/2, S_1)`,
> maximal among such; and `T_*` does not depend on `s` in that range.
> **(c')** If `T_* < inf` then `sup_{t < T_*} ||u(t)||_{L^inf} = inf`.
> **(d')** `u` is axisymmetric, and swirl free if `u_0` and `f` are.
> **(e')** On every closed `[0,T']` with `T' < T_*`,
> `u in L^inf([0,T']; L^2 cap L^6 cap L^inf)` and `grad u in L^inf(R^3 x [0,T'])`; `u` satisfies
> (H\*) in the finite-index form `u in C([0,T');H^s)` for every `s in (5/2, S_1)`.
> **(f')** `u in L^inf(0,T_*;L^2) cap L^2(0,T_*;H^1)` up to `T_*`.

*Proof.* Local existence and uniqueness are `[FENCES]` sec.1.3 Step 3, which is written out there
rather than cited (the citation it originally carried was read at source and found not to say
what was attributed to it): for `nu > 0` and a force obeying (A1), pairing with `Lambda^s u`
makes the dissipative term nonnegative and gives
`d/dt ||u||_{H^s} <= C ||u||_{H^s}^2 + ||f||_{H^s}` for `s > 5/2` (using `H^s subset W^{1,inf}`),
hence existence up to `t_0 + c/(||u(t_0)||_{H^s} + ||f||_{L^1(t_0,t_0+1;H^s)})` by a
regularisation-plus-compactness scheme, and uniqueness in `L^inf_t H^s` from the `L^2` estimate on
the difference. Higher regularity propagates: `u in L^2(t_0,T;H^{s+1})` gives an admissible `t_1`
with `u(t_1) in H^{s+1}`, and the same estimate at `s+1` from `t_1` on, iterated, yields every
`s' in (s, S_1)`. The iteration stops at `S_1` because the **datum** stops there, not because the
argument does.

`T_*` is independent of `s` in the range: `T_*(s)` is non-increasing in `s`; and if
`T_*(s) < T_*(s')` for some `s > s'` then on `[0, T_*(s))` the `L^inf` norm is bounded, since
`H^{s'} subset L^inf`, so `[FENCES]` Lemma A at `alpha = 2`, whose hypothesis (A1) holds at index
`s` by (a'), gives `sup_{t < T_*(s)} ||u(t)||_{H^s} < inf` and the local theory continues past
`T_*(s)`, a contradiction. (c') is that Lemma read contrapositively at any single admissible `s`.
(d') is K(d) unchanged: rotations about `e_z` and the reflection `(x_1,x_2,x_3) -> (x_1,-x_2,x_3)`
preserve the equation, the class of (b'), the datum and the force, and uniqueness does the rest.
(e') is Sobolev embedding on a compact time interval, needing only `s > 5/2`. (f') is the energy
identity of `[FENCES]` Lemma A Step 1, which uses only `u in L^2` and the force in `L^1_t L^2`. ∎

Corollary C1 of the addendum is unchanged: its Lemma U argument needs the reference in
`L^inf_t(L^2 cap L^6 cap L^inf)` with `grad u in L^inf`, which is K'(e'), and no higher index.

### 3.4 Every consumer, and the highest derivative it uses

| consumer | object | highest order in `u` |
|---|---|---|
| `u2` sec.2, `\|\|grad u\|\|_op <= 2\|a\| + r\|grad a\| + \|omega^theta\|` | `grad u`, `omega` | 1 |
| `u2` sec.4, `Gamma_rad` and `x . b <= Gamma_rad \|x\|^2` | `grad b` | 1 |
| `u2` sec.5 **P1**, `sup rho\|eta\| <= e^{c_R} E_0 M` | `eta` a classical solution of `D_t eta = nu Lap_5 eta` | 3 |
| `u2` sec.5 **P2**, `sup rho^2\|grad eta\| <= e^{p c_G} Gfrak_0 M` | `grad eta` a classical solution of its own equation | **4** |
| `u2` sec.6, `r\|grad a\|` by the Riesz composition | `grad_5 eta` | 2 |
| `pmax-h3v` Part A (the P1/P2 proofs, barriers, comparison) | as `u2` sec.5; interior parabolic Schauder supplies the classical regularity | **4** |
| `pmax-h3v` (D1) to (D5) on the datum | `eta_0 in W^{1,inf} cap C`, `z`-odd, real-analytic at the poles, `rho^{-9}` tail | 1 |
| Theorem V.1, the Feynman-Kac representation | `eta` a bounded classical solution | 3 |
| Theorem V.4' hypothesis `(H3'_vartheta)` | `d_v^k(eta_0 - eta_P)`, `k <= 4`, on `B(x_*,d)` | **5**, on the datum |
| `hk2` `(H-K2)`, `K_2 = \|\|grad^2_5 b\|\|_{L^inf(N_tau)}` | `grad^2 b`, and `sup_B \|\|grad^2 eta\|\|_F` | 3 |
| Lemma T'/T'' hypotheses (D), (H1), (H2), (H3) | `eta_0` measurable with `\|eta_0\| <= M/r`; `Phi` a `C^1` diffeomorphism; `lambda(.,s) in C^1` | 2 |
| `u1` (P1), (P2), (P3) and the majorant | the `l = 1` coefficient of `omega^theta` and its `rho`-derivative | 2 |
| `(Gamma-off)`, `u2` sec.8 | `Gamma = \|\|grad u\|\|_inf` | 1 |
| Corollary C1 of the addendum | the reference in `L^inf_t(L^2 cap L^6 cap L^inf)`, `grad u in L^inf` | 1 |

The largest is **five**, and it is a condition on the datum at `t = 0`, on a ball whose closest
point to the origin is `rho_* - d = 3.4618602587 rho_0`, where `eta_0` is real-analytic; the
origin's `C^{7,1}` corner is nowhere near it. `H^s subset C^k` for `s > k + 3/2`, so `s < 10.5`
gives `u_0` and `u(t)` in `C^8` in space at every `t < T_*`, and for `t > 0` with `nu > 0`
interior parabolic regularity gives `C^infinity` anyway. **There is no gap.**

---

## 4. UNIT F-4: THE CONTINUITY AND FIRST-CROSSING STEP

### 4.1 What was wrong

`u2`'s theorem carries **(H4)**, "the bootstrap posture of `B(t)`:
`Gamma(sigma) := ||grad u(.,sigma)||_{L^inf} <= Gammabar` for `sigma <= s`", as an assumption.
Its sec.8 shows that a self-consistent `Gammabar` **exists** for `L >= L_Gamma^exist` and bisects
for the existence boundary. Existence of a fixed point is not the same statement as the a priori
assumption propagating. Row `D7` recorded PROVED. The step is routine. It was not written.

### 4.2 The argument

Throughout, the contradiction hypothesis of (S3) is in force: `T_d(u_0) > tau`, `tau = c/(ML)`,
and `T_d <= T_*` by definition, so by Proposition K' the solution is strong on `[0,tau]`.

**The three functions.**

```
    Gamma(s) := ||grad u(.,s)||_{L^inf(R^5)} ,
    m(s)     := sup_x |Phi_s(x) - Lambda_s(x)|/|x| ,        mu(s)   := lam_max^2 m(s) ,
    w(s)     := sup_x |log( J_Phi(x,s)/lambda(|x|,s)^2 )| ,  mu_J(s) := e^{w(s)} - 1 .
```

Each is continuous on `[0,tau]`. `Gamma`: `H^s subset W^{1,inf}` for `s > 5/2`, and
`u in C([0,tau];H^s)` by K'(b'), so `s -> grad u(.,s)` is continuous into `L^inf`. `m` and `w`:
`Phi_s` is the flow of `b`, Lipschitz in `x` uniformly on `[0,tau]` with constant `Gammabar` and
continuous in `s`; `Lambda_s(x) = T_{lambda(|x|,s)} x` with
`lambda(rho,s) = exp int_0^s frak_a(rho,sigma) dsigma` continuous in `s` and, by unit 5, `C^1` in
`rho`.

**The closed set.**

```
    Bfrak := { s in [0,tau] : for all sigma in [0,s] ,
                 (i)   Gamma(sigma) <= Gammabar := (3/2) M L + C'' M ,
                 (ii)  mu(sigma)    <= mubar    := 1/2 ,
                 (iii) mu_J(sigma)  <= mu_Jbar  := 1/2 } .
```

`Bfrak` is closed by continuity. The strain lower bound
`a_ref(sigma) >= (M/2) r_h A(sigma) > 0` is a **consequence** of (P1) `lambda >= 1` and `h >= 0`,
not a bootstrap hypothesis; it is recorded here so the reader can see the argument is not
circular.

**`s = 0`, with a strict margin of order `M L`.** By Consequence A and the evenness of
`A_sigma(phi) cos phi sin^2 phi` about `phi = pi/2`, `F(rho',0) = M Theta(rho') kappa_delta`, and
`int Theta dlog rho = L - 2 eps_r = L - 1/2` exactly, so `a(0,0) = M kappa_delta (L - 1/2)`
exactly. At `s = 0`, `c_G = int_0^0 Gamma = 0`, so (H4) is vacuous and `u2`'s chain applies
unconditionally, giving `Gamma(0) <= 2 a(0,0) + C''_0 M` with `C''_0` the `(Gamma-off)` constant
at `c_G = 0`. Hence

```
    Gammabar - Gamma(0)  >=  (3/2 - 2 kappa_delta) M L  +  kappa_delta M  +  (C'' - C''_0) M ,
    3/2 - 2 kappa_delta  =  0.5164262397 > 0 .
```

`mu(0) = mu_J(0) = 0 < 1/2` because `Phi_0 = Lambda_0 = id` and `J_Phi(.,0) = lambda(.,0)^2 = 1`.

**Propagation, with the same strict margin.** Let `s_*` in `Bfrak`. On `[0,s_*]` hypothesis (H4)
holds with `c_G = int_0^{s_*} Gamma <= Gammabar tau = (3/2 + C''/L) c`, which is the `c_G` at
which `C''` was computed, so `u2`'s theorem applies at every `sigma <= s_*` and gives
`Gamma(sigma) <= 2 a(0,sigma) + C'' M <= lamhat M L + C'' M`, where
`lamhat := sup_{[0,tau]} ||omega(sigma)||_inf/M`. Now `omega = curl u` is continuous into `L^inf`
(`H^{s-1} subset L^inf` for `s > 5/2`), `[0,tau]` is compact, and the contradiction hypothesis
gives `||omega(sigma)||_inf < (3/2) M` at every `sigma in [0,tau]`; so the supremum is attained
and `lamhat < 3/2`. With `theta_0 := (3/2 - lamhat)/2 > 0`,

```
    Gamma(sigma)  <=  Gammabar - 2 theta_0 M L      for every sigma <= s_* .
```

**This is where strictness comes from, and it comes from the contradiction hypothesis, not from
the constants:** `theta_0` depends on the solution, `C''` does not. Simultaneously, on `[0,s_*]`
every hypothesis of `u1` sec.2.3's majorant holds, so
`dm/dtheta = G m + lam_max lam_om (C_R + 2 log lam_max)/L`, `m(0) = 0`, `G = Gammabar/(ML)`, whose
solution is `m(theta) = (drive/G)(e^{G theta} - 1)`, and `mu = lam_max^2 m`; and (2.4) gives
`mu_J`.

**Conclusion.** `Bfrak` is nonempty, closed, and relatively open in `[0,tau]`: at any `s_*` in
`Bfrak` all three bounds hold on `[0,s_*]` with a strict margin, and all three functions are
continuous, so they persist on `[s_*, s_* + delta]` for some `delta > 0`. `[0,tau]` is connected,
so `Bfrak = [0,tau]`, and (H4) holds on the whole window without being assumed. ∎

### 4.3 Every constant, and every margin

Proved column, `sigma_* = 0`, `f` immaterial:

| `L` | window | `C''` | `C''` at `c_G = 0` | `Gammabar/(ML)` | `Gamma(0)/(ML)` bound | margin | `mu(c)` | `mu_J(c)` |
|---|---|---|---|---|---|---|---|---|
| `1e+05` | `c = c_*` | `1034.4930` | `129.5940` | `1.5103449` | `0.9848648` | `0.5254801` | `5.9540e-02` | `1.7548e-03` |
| `1e+05` | cap | `1539.5613` | `129.5940` | `1.5153956` | `0.9848648` | `0.5305308` | `1.6481e-01` | `2.5288e-03` |
| `1e+06` | `c = c_*` | `1005.5792` | `129.5940` | `1.5010056` | `0.9837029` | `0.5173027` | `5.7380e-03` | `1.7399e-04` |
| `2e+06` | `c = c_*` | `1004.0511` | `129.5940` | `1.5005020` | `0.9836383` | `0.5168637` | `2.8633e-03` | `8.6953e-05` |
| `2e+06` | `c/c_* = 1.0669124` | `1165.2701` | `129.5940` | `1.5005826` | `0.9836383` | `0.5169443` | `4.2647e-03` | `1.0073e-04` |
| `2e+06` | cap | `1457.5018` | `129.5940` | `1.5007288` | `0.9836383` | `0.5170904` | `7.6814e-03` | `1.2449e-04` |
| `1e+07` | `c = c_*` | `1002.8342` | `129.5940` | `1.5001003` | `0.9835867` | `0.5165136` | `5.7176e-04` | `1.7384e-05` |

Every margin is strict in every row. `mu(c) <= 1.65e-01 < 1/2` and `mu_J(c) <= 2.53e-03 < 1/2`
across the whole table, so `(1 - mu)^5 > 0` and `lam_max^{-2} - m > 0` throughout, which is what
Lemma T'' (H2) and `u1` (2.4) need. `mu` at `L = 2e+06` and the self-consistent window is
`4.2647e-03`, which is `FIX4.md` sec.6.2's own `4.264720e-03`.

**Where strictness is obtained, one line each.** (i) from `lamhat < 3/2`, i.e. from the
contradiction hypothesis on a compact interval; (ii) from `mu(c) <= 1.65e-01`; (iii) from
`mu_J(c) <= 2.53e-03`; at `s = 0`, from `3/2 - 2 kappa_delta = 0.5164262397 > 0`.

Row `D7` is now **PROVED**, not PROVED-modulo.

---

## 5. UNIT F-5: LEMMA T'' ON `[1, lam_max]`, AND `lambda(.,s) in C^1`

### 5.1 Where `lambda <= 3/2` is used, step by step

`write/lemma-T-shell-dependent/PROOF.md` line 129 states Lemma T' for
`lambda in C^1([rho_0,R];[1,3/2])`, and its own parenthetical at lines 138 to 139 says the range
enters only in Step 0's condition and in Corollary 2's `r_h`. `x5_lemmaT.py` checks that
quantitatively at `lambda` well above `3/2`:

| step | what it uses | range-dependent? | check |
|---|---|---|---|
| 0(a) | `J_Lambda = lambda^2[1 + (rho lambda'/lambda)(1 - 3 cos^2 phi)]`, `1 - 3 cos^2 phi in [-2,1]` | no | identity reproduced to `7.63e-15` over `125` cases with `lambda` to `2.5` and `rho lambda'/lambda` free; the bracket's range measured as `[-2, 1]` |
| 0(b) | injectivity by contraction, `\|d log ghat/d log lambda\| <= 2` | no | supremum over `psi in [0,pi]` and `lambda in [1, 2.5]` is `2.000000000000`, attained on the axis; it is a convex combination of `-1` and `+2` for every `lambda > 0` |
| 0 | the standing condition `2 kappa_s/L < 1` | no | `kappa_s <= (3/4) c` by `u1` (P2), from the contradiction hypothesis; at the certified cap `2 kappa_s = 1.5 c = 1.4404`, against `L >= 1.4e+06` |
| 2 | `sup\|w\|^4\|K\| = 3/(8 pi^2)`, `sup\|w\|^5\|grad K\| = 3/(2 pi^2)` | no | properties of the kernel in `R^5`; no `lambda` appears |
| 3 | `int_0^pi sin^2 phi/g(phi;lambda)^4 dphi = pi/(2 lambda)` | no | 30-digit quadrature at seven `lambda` up to `3`, maximum relative error `1.88e-31`; the source proves it for every `lambda > 0` |
| 4 | `(2.7)`'s `(1-mu)^5`, from `\|w\| >= (1-mu)\|Lambda x\|` on the segment | no | a statement about `mu`; the `\|Lambda x\|^{-4}` it multiplies is integrated by (2.6), which uses Step 3 at the same `lambda` |
| 5 | `\|J_Phi - lambda^2\| <= mu_J lambda^2` | no | hypothesis (H3) |
| Cor. 2 | `r_h := inf_{lambda in [1,3/2]} P_h(lambda)/lambda` | **yes** | the one place; recomputed on `[1, lam_max]` |

### 5.2 Lemma T''

> **LEMMA T''.** Let `lam_max >= 1` and `lambda in C^1([rho_0,R];[1,lam_max])` satisfy (S) with
> `2 kappa_s/L < 1`, let `eta_0` satisfy (D), and let `Phi` satisfy (H1), (H2) with `mu` and (H3)
> with `mu_J`. Then `a[eta_0 o Phi^{-1}](0)` and `a_ref[eta_0;lambda]` both converge absolutely,
> `|a_ref| <= (3 pi/8) M A`, `|a[eta_0 o Phi^{-1}](0)| <= (3 pi/8)(1+mu_J)(1-mu)^{-4} M A`, and
> ```
>    | a[eta_0 o Phi^{-1}](0) - a_ref[eta_0;lambda] |
>              <=  pi M A [ 3 mu (1 + mu_J)/(2(1-mu)^5) + 3 mu_J/8 ] .
> ```
>
> **COROLLARY 2''.** For `omega_0^theta = -M sgn(z) h(phi)`, with
> `r_h := inf_{lambda in [1, lam_max]} P_h(lambda)/lambda`,
> ```
>    |Delta|/a_ref  <=  (2 pi/r_h) [ 3 mu (1 + mu_J)/(2(1-mu)^5) + 3 mu_J/8 ] .
> ```

*Proof.* Lemma T's proof verbatim. Steps 0(a), 0(b), 2, 3, 4, 5 and 6 use no upper bound on
`lambda`, as the table above checks; Step 0's standing condition is `2 kappa_s/L < 1`, and
`kappa_s` comes from `u1` (P2) and the contradiction hypothesis, not from the range of `lambda`;
Corollary 2's `r_h` is the only range-dependent object, and it is an infimum over a larger
interval. Every constant is identical to Lemma T's and Corollary 2's. ∎

The `lam_max` the chain uses is `e^{3c/4}`: `1.8558724485` at `c = c_*` and `2.0548738638` at the
certified cap of unit 6. `P_h` stops increasing at `lam_mono = 2.0693384635`, so the invocation
runs to within `0.70 %` of a real obstruction; that margin is the window cap, and unit 6 tightens
it. `r_h` on `[1, lam_max]` is unit 6's certified value, `0.9033433307` at `c = c_*` and
`0.8784367504` at the window `L = 2e+06` actually uses. The source lemma's own value on
`[1, 3/2]` is `0.981822`.

**Instrument control.** `P_h(lambda)/lambda` recomputed here from the Gauss-Legendre profile,
against `fix4`'s FFT-grid instrument, at four `lambda`: relative agreement `2.37e-10` or better.

### 5.3 `lambda(.,s) in C^1`

`u1/PROOF.md` line 363 asserts this bare. Here it is, from the chain's own constants.
`u1` sec.1.1 (P2) proves the identity `d frak_a/d log rho = -F(2 rho, s)` exactly, and by
Consequence A, `F(rho',s) = -(3/4) int_0^pi omega^theta(rho',phi,s) cos phi sin^2 phi dphi`.
Differentiating under the integral, with `omega^theta = r eta`, `r = rho sin phi`, so
`rho d_rho(r eta)|_phi = r eta + r rho d_rho eta`, and `int_0^pi |cos phi| sin^2 phi dphi = 2/3`,

```
    | dF/d log rho' |  <=  (3/4)(2/3) sup |rho d_rho omega^theta|
                       <=  (1/2) ( sup rho|eta| + sup rho^2|grad eta| )
                       <=  (1/2) ( e^{c_R} E_0 + e^{p c_G} Gfrak_0 ) M ,
```

by `u2`'s P1 and P2, using `e^{c_R} <= e^{c_G}` (conservative, as `fix3` `g1` sec.D does). Hence
`F(2 . , s)` is Lipschitz in `log rho`, uniformly on `[0,tau]`; `frak_a(.,s)` is `C^{1,1}`; and
`lambda(rho,s) = exp int_0^s frak_a(rho,sigma) dsigma` is `C^1` in `rho` with
`rho d_rho lambda/lambda = - int_0^s F(2 rho, sigma) dsigma`, continuous in `rho`. That is exactly
Lemma T's (H1) and (S), with the `kappa_s` that (P2) already supplies. ∎

Measured, proved column, `p = 2.1107` (`u2` sec.8.2's largest tabulated value):

| `L` | window | `c_G` | `sup rho\|eta\|/M` | `sup rho^2\|grad eta\|/M` | `\|dF/d log rho\|/M` |
|---|---|---|---|---|---|
| `1e+05` | `c = c_*` | `1.24524` | `26.2349` | `499.7223` | `262.9786` |
| `1e+05` | cap | `1.45515` | `32.3626` | `778.3004` | `405.3315` |
| `1e+06` | `c = c_*` | `1.23754` | `26.0337` | `491.6662` | `258.8500` |
| `2e+06` | `c = c_*` | `1.23712` | `26.0229` | `491.2356` | `258.6292` |
| `2e+06` | cap | `1.44113` | `31.9119` | `755.5995` | `393.7557` |
| `1e+07` | `c = c_*` | `1.23679` | `26.0143` | `490.8923` | `258.4533` |

Finite at every `L` and every window in the theorem's range, so `lambda(.,s) in C^1` and row `F3`
is discharged as stated.

---

## 6. UNIT F-6: THE ENCLOSED LIST, CERTIFIED

### 6.1 The proved inputs

`fx_profile.W_bounds` supplies, in closed form,

```
    ||W_sigma||_inf   <= ||Wtil||_inf = 1/sin delta = 7.6612975755 ,
    ||W_sigma^{(k)}||_inf <= ||Wtil'||_inf ||chi_sigma^{(k-1)}||_1 ,  k >= 1 ,
    ||Wtil'||_inf = cos delta/sin^2 delta = 58.1933325682 ,
    ||chi||_1 = 1.0000000000 ,  ||chi'||_1 = 2/(sigma sqrt(2 pi)) = 39.8942280401 ,
    ||chi''||_1 = 4 e^{-1/2}/(sigma^2 sqrt(2 pi)) = 2419.7072451914 ,
    ||chi'''||_1 = (2 + 8 e^{-3/2})/(sigma^3 sqrt(2 pi)) = 188751.6250163096 ,
```

each `L^1` norm being the total variation of the previous derivative, evaluated at its critical
points; hence

```
    B_W = [7.6612975755, 58.1933325682, 2321.5780798926, 140810.8284371617, 10984086.0873664636]
    B_A = [7.6612975755, 65.8546301438, 2445.6260426045, 147957.8039721196, 11561499.3042223137]
```

with `B_A(k) = sum_j C(k,j) B_W(j)` from `|sin^{(m)}| <= 1`. For a smooth interior maximum on a
uniform grid of step `h` the pad is `(1/8) h^2 sup|f''|`, because the nearest grid point is within
`h/2` of the argmax and `f'` vanishes there; where the maximised function is only piecewise smooth
or the maximum can sit at a region boundary, the pad is `(h/2) sup|f'|`, which needs no interior
hypothesis. Upper bounds below are always the **global** grid maximum plus the pad; local
refinements are used only to raise lower bounds, since a refinement of one window says nothing
about the rest of `[0,pi]`.

### 6.2 The certified values, at `h = 3.1416e-05`

| quantity | certified enclosure | `fix4` | inside? |
|---|---|---|---|
| `N_sigma` | `[1.0124508487, 1.0124511502]` | `1.0124508488` | yes |
| `sup \|W_sigma\|` | `[7.6463403352, 7.6463406216]` | (implicit) | -- |
| `E_0 = sup rho\|eta_0\|/M` | `[7.5523054461, 7.5523079773]` | `7.5523076942` | yes |
| `Gfrak_0` (branch bound) | `[36.0783630557, 36.0784195640]` | `36.0795702396` | **no**, see sec.6.3 |
| `sup \|h_sigma'\|` | `<= 11.3628286061` (grid `11.3248852293`, pad `3.7943e-02`) | grid `11.3248852` | -- |
| `\|h_sigma''\|` | `<= 2415.5503900804` (proved) | `225.9466903722`, a grid scan | -- |
| `Psi_sigma(0)` | `<= 495.1064102419` | `489.0279393839` | yes, below |

`Gfrak_0`'s branch bound needs a word. The maximised function
`max(P + S if 81 P <= 8 S else YMAX P + S)`, `P = (W_sigma/N_sigma)^2`, `S = (W_sigma'/N_sigma)^2`,
`YMAX = 6.4072265625`, is discontinuous where the branch switches, so its cells are split: cells
where `81 P - 8 S` cannot change sign inside the cell (`|81P - 8S|` above
`margin + (h/2) Lip_cond`, `Lip_cond = 2179227.06`, `margin = 34.23`) carry the interior pad
`(1/8) h^2 (2.6769e+07) = 3.3024e-03`, applied to the grid maximum `1301.6491`; the `2734` cells that could
contain a switch carry the Lipschitz pad `(h/2)(269169.38) = 4.2281`. The maximum sits at
`phi = 170.9532 deg` with `8 S - 81 P = 6844.05`, far outside any switch cell, and the switch
cells' own bound is `860.3189`, well below the interior cells' `1301.6524`.
So the certificate closes on the interior
branch alone.

`Psi_sigma(0)` is certified by the separated bound
`Psi <= sup_u e^{-2u}( |rad(u)| sup|w_0| + |Theta(u)| sup|angL| )`, with
`sup|angL| <= 1371.7827305119` on `phi in [0.02, pi - 0.02]` (grid plus the proved Lipschitz pad)
and `<= 4 sup_{[0,0.02]}|w_2| = 10.3435246793` on each polar cap (there `W_sigma` is even, so
`|w_1(phi)| <= phi sup|w_2|` and `3 cot(phi)|w_1| <= 3 sup|w_2|` since `phi cot phi <= 1`), and
with `u` split at `-1` and `8`: for `u <= -1`, `Theta(u) <= e^{8(u - 1/4)}` and `|rad| <= 74 Theta`
give `Psi <= e^{6u-2}(74 E_0 + sup|angL|) <= 0.6476621037` at `u = -1`; for `u >= 8`,
`e^{-2u} <= e^{-16}` gives `1.6994e-04`; on `[-1,8]` a uniform cell bound with `hu = 2.0000e-06` and
the proved bounds on `|rad'|` and `|Theta'|` gives `495.1064102419`, attained at `u = 0.3915`.
That is `1.24 %` above `fix4`'s tight scan value, and `Psi_sigma(0)` enters only `eps_a`, which is
`3.9e-05` of `eps`, so the looseness costs nothing.

### 6.3 A discretisation error inside `fix4`, found on the way

`fix4`'s `Gfrak_0 = 36.0795702396` lies **above** the certified interval, by `1.15e-03`
(`3.2e-05` relative). The cause is located: `fix4/p1_profile.Profile` computes
`W_sigma' = fftconvolve(Wp, chi)*h`, a uniform Riemann sum for a convolution whose first factor
`Wtil'` **jumps** at the four kink cones. A jump inside a cell makes the composite rule first
order, with local error about `(h/2) x jump x chi(distance)`; at the `Gfrak_0` argmax
`phi = 0.15790692`, which is `1.35 sigma` from the cone `phi = delta`, that predicts
`(5.236e-06/2)(58.528)(8.03) = 1.23e-03`, and the measured discrepancy is `0.0012276591`. Three
checks confirm the direction (`fx_controls.json`):

```
   phi_G = 0.1579069196  (the Gfrak_0 argmax)
   referee's adaptive quadrature, W_sigma'(phi_G)  =  -35.9605572617
   this seat's Gauss-Legendre,    W_sigma'(phi_G)  =  -35.9605572617
   fix4's FFT at nphi =  600001                    =  -35.9617849209   (h = 5.236e-06)
   fix4's FFT at nphi = 2400001                    =  -35.9608638645   (h = 1.309e-06)
   error ratio on the fourfold refinement          =  4.0041
```

The error falls by the factor `4` that a first-order rule and a fourfold refinement predict, and
the two breakpoint-aware instruments agree exactly. `W_sigma` itself (order `0`) is unaffected,
since `Wtil` is continuous: `fix4`'s `W_sigma(phi_G) = 6.4111336495` against this seat's `6.4111336484`,
relative `1.6700e-10`. So `kappa_delta`, `P_h`, `r_h` and the window cap, which read only
`W_sigma`, are untouched, and the only casualty is `Gfrak_0`, which `fix4` reports `3.2e-05` too
**large**. That is the conservative direction: `Gfrak_0` is the binding drive constant of `Ghat`,
`C''`, `C_R` and hence `L_*`. Graded MINOR, direction safe, reported because a sheet that hides
its own instrument's error has not been checked.

### 6.4 `r_h` and `min P_h'` with the proved Lipschitz constant

`fix4`'s certificate is re-run unchanged with `lip_pad = |h_sigma''| x h_grid` and the scan value
`225.9466903722` replaced by the proved `2415.5503900804`, i.e. the pad on `sup|h_sigma'|` moved
from `1.1830541035e-03` to `1.2647792258e-02` (`fix4` quotes `1.1831e-03`, the referee
`1.2648e-02`; both reproduce). At `nsub = 20000`, `fix4`'s own value:

| window | `lam_max` | `r_h`, scan pad | `r_h`, **proved** pad | `min P_h'`, scan pad | `min P_h'`, **proved** pad |
|---|---|---|---|---|---|
| `c = c_*` | `1.8558724485` | `0.9033435309` | **`0.9033433307`** | `0.3335513465` | **`0.3294825889`** |
| `c/c_* = 1.0669124225` | `1.9342710419` | `0.8784369398` | **`0.8784367504`** | `0.2016642746` | **`0.1978911613`** |
| `c/c_* = 1.0238897079` | `1.8834914830` | `0.8949627711` | **`0.8949625747`** | `0.2872563497` | **`0.2833025588`** |
| `c/c_* = 1.1664585076` (fix4's cap) | `2.0570755604` | `0.8331460368` | **`0.8331458640`** | `0.0004499697` | **`-0.0029985709`** |

The scan-pad column reproduces `THEOREM_S3_v2.md`'s `r_h >= 0.9033435309` and
`P_h' >= 0.3335513465` exactly, so the substitution is the only change. `r_h` moves by `2.0e-07`.
`min P_h'` moves by `4.1e-03`, and at `fix4`'s certified cap it goes **negative**: the certificate
no longer closes there. Bisecting for the largest window at which it does:

```
     certified window cap, proved pad :   c/c_* <= 1.1647266909 ,  i.e.  eps <= 0.1647266909
     (fix4's, scan pad :                  c/c_* <= 1.1664585076 ,        eps <= 0.1664585076 )
     lam_max at the certified cap     :   2.0548738638      (fix4: 2.0570755604)
```

This is the one place where unit 6 costs something. It costs a factor `0.9985153` on the window.

### 6.5 The budget on the certified constants

Re-running the whole budget with `E_0`, `Gfrak_0` at the certified **upper** ends, `r_h` with the
proved pad, and the certified cap, alongside unit 1's ball and unit 2's `eps_v`:

| column | `L_*` (`eps <= 0.1647266909`) | `log Lambda_*` | `L_*` (`eps <= 0.1`) | `log Lambda_*` |
|---|---|---|---|---|
| **all proved** | **`1424207.9721`** | **`2848418.8578`** | **`1554565.5693`** | **`3109134.0522`** |
| proved, `C_a` measured | `1379464.3531` | `2758931.6198` | `1499581.4464` | `2999165.8064` |
| all measured | `1761.8423` | `3526.5981` | `2050.2770` | `4103.4676` |

with `E_0 = 7.5523079773`, `Gfrak_0 = 36.0784195640`, `r_h` from the proved pad, and the certified
cap `eps <= 0.1647266909`. Against the same run on `fix4`'s constants (sec.2.4):

| column | `L_*` move at the cap | relative | `L_*` move at `eps <= 0.1` | relative |
|---|---|---|---|---|
| all proved | `402.5232` | `2.8255e-04` | `903.8310` | `5.8107e-04` |
| proved, `C_a` measured | `571.5623` | `4.1416e-04` | `873.0676` | `5.8187e-04` |
| all measured | `-0.4497` | `-2.5533e-04` | `0.9582` | `4.6713e-04` |

**Does `L_*` move?** In the proved column, by `402.5232` at the cap and `903.8310` at `eps <= 0.1`, in
each case **downwards**, i.e. `2.8255e-04` and `5.8107e-04` relative, in the favourable direction: the
corrected `Gfrak_0` more than pays for the tightened window cap. The answer to the brief's
question is yes, by about three parts in ten thousand at the cap and six at `eps <= 0.1`, and the
sign is the good one.

---

## 7. WHAT IS STILL NOT A CERTIFIED ENCLOSURE

Listed so that `THEOREM_S3_v3.md` sec.3 can quote it rather than paraphrase it.

1. **`vartheta`** is a grid-plus-polish lower bound on a supremum over a five-dimensional set,
   inflated by `SAFETY = 2`. That is `pmax-h3v` sec.B5's convention and this sheet keeps it. Its
   margin: `2 vartheta = 0.2623138022` against the requirement `< 1`, a factor `3.81`; and the
   coarse-to-fine change is `4.6e-03`. It is not an interval computation.
2. **`C_E`** is a convergence statement, not an enclosure: `7.3e-06` in `lmax`, `4.7e-07` in
   `du`. It enters only as `(2/5) log C_E` in the dictionary, so `7.3e-06` relative moves
   `log Lambda_*` by `3e-06` out of `2.8e+06`. Unchanged from `fix4`; not attempted here.
3. **`C''`, `Ghat`, `C_R`, `chat_a`, `L_Gamma^exist`** are quadratures with a sign-change
   existence proof. `C_R`'s supremum over the angle is a sampled search; the referee checked it
   at six resolutions and found agreement to `3e-08` at an interior smooth maximum. Unchanged
   here.
4. **`K_2 <= 161.7735 M/rho_0`** is imported from `hk2` for its own datum and re-checked
   pointwise at four `lambda` by `fix4/p3`. The referee swept `C_K` past it by four orders of
   magnitude and `L_*` moved in the fourth digit only at `C_K = 1e+04`. Not load-bearing.
5. **The far and near constants** are Parseval-based numerical evaluations in the source seat,
   not intervals. Their three hypotheses are profile-free and were re-verified for this datum by
   `fix4/p3` and independently by the referee.
6. **Interior parabolic Schauder** is cited, not proved, in `pmax-h3v` Part A. That is the one
   textbook citation the P1/P2 chain rests on, and it is the only genuinely unproved analytic
   item left in the graph.
7. **`eps(L)` monotone decreasing to the floor** is tabulated at five values of `L`, not proved.
8. **The geometry margin in unit 1.** `4 sigma` is a quarter of a `sigma` past the cliff of
   sec.1.4. The `6 sigma` statement costs nothing and is safer; the brief specified `4 sigma` and
   that is what the headline carries.

---

## 8. THE DECORRELATION TABLE

Every control this sheet runs, with the number it had to hit and the number it returned.

| control | reference | here | relative |
|---|---|---|---|
| `W_sigma`, `W_sigma'`, `W_sigma''`, `W_sigma'''` at ten angles | referee's adaptive quadrature | Gauss-Legendre panels | `2.0272e-10` |
| `W_sigma''''` | fourth-order finite difference of `W_sigma'''` | direct | `2.8988e-09` |
| `W_sigma''''(delta)` | referee's spline table | direct | `3.8052e-09` |
| `vartheta`, kinked, `f = 4`, ASSEMBLY's `d` | `pmax-h3v` `0.055612` | `0.0552690` | `6.2e-03`, below |
| `vartheta`, mollified, `f = 4`, ASSEMBLY's `d` | referee `8.76183` | `8.78252` | `2.4e-03`, above |
| `viscous_budget` transcription at `vartheta = 0` | imported function | transcribed | `0.000e+00` at four points |
| `L_*` proved, `eps <= cap` | `fix4` `1424610.4953` | `1424610.4953` | exact re-run |
| `L_*` proved, `eps <= 0.1` | `fix4` `1555469.4004` | `1555469.4004` | exact re-run |
| `r_h` certified, scan pad, `c = c_*` | `THEOREM_S3_v2` `0.9033435309` | `0.9033435309` | exact |
| `min P_h'` certified, scan pad, `c = c_*` | `THEOREM_S3_v2` `0.3335513465` | `0.3335513465` | exact |
| `N_sigma` | `fix4` `1.0124508488` | `[1.0124508487, 1.0124511502]` | contains |
| `E_0` | `fix4` `7.5523076942` | `[7.5523054461, 7.5523079773]` | contains |
| `Psi_sigma(0)` | `fix4` `489.0279393839` | `<= 495.1064102419` | contains |
| `Q(lambda) = P_h/lambda` at four `lambda` | `fix4`'s FFT grid instrument | Gauss-Legendre profile | `<= 2.4e-10` |
| `mu` at `L = 2e+06`, self-consistent window | `FIX4` `4.264720e-03` | `4.2647e-03` | exact |
| `J(lambda) = pi/(2 lambda)` | closed form | 30-digit quadrature, seven `lambda` | `<= 1.9e-31` |
| `sup \|d log ghat/d log lambda\|` | `2` | `2.000000000000` | attained |
| Fourier decay exponent of `\|x\|^7 x_1 e^{-\|x\|^2}` | `-11` | `-11.00018` at `k = 1600` | `1.6e-05` |
| inner edge closed form | `rho^8/(rho^8 + e^2)` | `Theta` at eight radii | `1.12e-12`, and exact in symbolic arithmetic |
| `W_sigma'` at the `Gfrak_0` argmax | referee `-35.9605572617` | `-35.9605572617` | exact agreement |
| `W_sigma` at the `Gfrak_0` argmax | `fix4`'s FFT `6.4111336495` | `6.4111336484` | `1.6700e-10` |

---

## 9. WHAT THIS SHEET CANNOT CATCH

It repairs six named findings and recomputes what they touch. It does not check that `C_R` is the
right constant, that `u2`'s P1 and P2 are the right instruments, that the modulo list is
complete, that the clock argument is right, or that the two remaining citations say what they are
said to say. `vartheta` and `C_R` remain sampled suprema with safety factors rather than
interval computations. `K_2` and the measured column are measurements, and numerics falsify but
never prove.

One thing it does catch that no earlier sheet did: **`fix4`'s own `W_sigma'` carries a
first-order discretisation error at the kink layer**, because an FFT convolution of a function
with a jump is a Riemann sum through that jump. It is `3.2e-05` of `Gfrak_0`, in the
conservative direction, and it would have been invisible to any instrument built the same way.

---

## GATE, AND FILES

```
424 CHECKS, 424 PASS, 0 FAIL
```

`check_fix5.py` re-reads every displayed constant from the stored JSONs, checks that the string
actually appears in `FIX5.md` or `THEOREM_S3_v3.md`, re-derives the ratios and factors quoted
between them rather than copying, and verifies the byte copies in `copies/` against their source
files.

| file | what it establishes |
|---|---|
| `fx_profile.py` / `fx_controls.json` | the third profile instrument: `W_sigma^{(k)}`, `k = 0..4`, by fixed-node Gauss-Legendre convolution in the integration variable, with the proved `L^1` norms of the Gaussian derivatives and the derived `B_W`, `B_A` |
| `x1_vartheta.py` / `x1_results.json` | unit F-1: the five-parameter `vartheta` search on the corrected ball, the two controls, the offset scan, the grid stability |
| `x2_budget.py` / `x2_results.json` | unit F-2: `eps_v` with `E_vartheta` and `Q`, the reduced `d` propagated, `eps_v`'s share, `L_*` three columns at both targets, the `6 sigma` variant, and unit F-6's budget re-run on the certified constants |
| `x3_propK.py` / `x3_results.json` | unit F-3: Proposition K', the consumer table, the inner-edge closed form, the origin singularity in exact arithmetic, the Fourier exponent |
| `x4_bootstrap.py` / `x4_results.json` | unit F-4: the closed set, every constant, every strict margin, the continuity inputs |
| `x5_lemmaT.py` / `x5_results.json` | unit F-5: Lemma T'' and the step-by-step range check, Corollary 2'' with its `r_h`, and the `C^1` proof with its constants |
| `x6_certified.py` / `x6_results.json` | unit F-6: the proved derivative bounds, the certified `N_sigma`, `E_0`, `Gfrak_0`, `Psi_sigma(0)`, `sup\|h'\|`, `\|h''\|`, `r_h`, `min P_h'`, the certified window cap, and the `fix4` `W_sigma'` diagnosis |
| `copies/` | byte copies of the referee's scripts and results, of `fix4`'s scripts, results and `imported/` tree, and of `pmax-h3v`'s `q3`/`q4` |
| `fx_stdout.txt`, `x1_stdout.txt` .. `x6_stdout.txt` | the console output of each run, kept and hashed |
| `check_fix5.py` | the gate |
| `SHA256SUMS` | computed with `shasum -a 256`, never typed |

**FL-000 stands. Nothing here touches the headline problem.**

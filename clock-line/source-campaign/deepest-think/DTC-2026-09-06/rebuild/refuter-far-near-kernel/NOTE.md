# REFUTER REPORT — `far-near-kernel-lemma` (GATE H-c)
Seat `refuter-far-near-kernel`, DTC-2026-09-06.  PREREG.md written before any script here ran (mtime 13:12:08;
first run 13:1x).  Target folder re-run from a COPY in `rerun/`.  Every number below from a script in this folder.

## Verdict
**REFUTED — MAJOR.**  The mathematical core survives an independent re-derivation in full.  Four defects:
one reported *constant* is not a constant, the gate verdict contradicts its own registered decision rule,
one additive find is false, and the assembled inequality is false as written.  None of the corrections
reverses a conclusion; the constant correction **strengthens** §H.3(a).

## R1 PROVENANCE — did not fire
`PREREG.md` mtime 2026-09-06 12:34:13; every script 12:37–13:03, every output 12:40–13:05, NOTE 13:05.
PREREG precedes every run.  All kills reported in the build result (K1,K2,K3,K3',K4,K5,K5',K6) are in
PREREG with the thresholds actually used.  No registered kill was silently dropped.  `shasum -a 256 -c
SHA256SUMS` in the target folder: 17/17 OK.

## R2 REPRODUCTION — did not fire
Re-ran s1–s6 from the copy.  Every headline number reproduces to the printed precision, including
C1=0.756945/0.291999, C2in=0.026093/0.014754, 3.999218/s+pi/8, K1 diffs −1.13e−5/−1.93e−6/−3.69e−6,
axis control 4.15888308 (rel 2.1e−16), Parseval 1.997888, exponents 5.0000/5.0000 and 4.0009/4.0009,
the far/collar/inner table, near-flatness slope −0.00099, c_edge(10°)=+0.216773 and the finite-R table,
and the whole s6 profile.  (s5's c_edge(90°) prints +0.115901 vs s3's +0.115737 — Cesàro tail, immaterial.)

## R3 INDEPENDENT DERIVATION — did not fire (`r1_independent_algebra.py`)
By my own route (mpmath dps=40, numeric z-differentiation of the solid harmonics — not sympy `simplify`):
* `d_z[rho^l C_l^{3/2}(z/rho)] = (l+2) rho^{l-1} C_{l-1}^{3/2}`, l=1..20, 3 points each: ALL AGREE.
  At l=2 the measured constant is **4.0**, so PREREG's `(2l+1)=5` is false and the repair `(l+2)` is right.
  Closed form independently: `c_l = l C_l^{3/2}(1)/C_{l-1}^{3/2}(1) = l·[(l+1)(l+2)/2]/[l(l+1)/2] = l+2`.
* `d_z[rho^{-l-3}C_l] = -(l+1) rho^{-l-4} C_{l+1}`: ALL AGREE (and by hand `(l+3)C_l(1) = (l+1)C_{l+1}(1)`).
* `N_l = (l+1)(l+2)/(l+3/2)` by quadrature, l=0..6: agrees to 1e−25.
* `alpha_l = g_l/(2l+3)`: from the single-layer jump `d_rho psi|_- − d_rho psi|_+ = g_l` ⟹ `(l + (l+3))alpha = g_l`.
* `Phi(0) = −(3/5)g_1 = −(3/4)INT_0^pi w cos sin^2 dphi`; bang-bang `g_1 = −5/6`, `Phi(0)=1/2`;
  `sup|Phi(0)|/M = (3/4)INT|cos|sin^2 = 1/2` ⟹ **kappa_0 = 1/2 is sharp**.  Confirmed.
* `C1 = 0.7569449261 / 0.2919985325`, `C2in = 0.02609298143 / 0.01475427158` — agree to 7 digits.
  Leading coefficients checked in closed form: `Q(u)/u → (4/7)·3/sqrt(24/7) = 0.925820`;
  `S(v)/v^4 → sqrt(3)/2 = 0.866025`; `S(v)/v^5 → sqrt(2.4) = 1.549193`.
* Collar: `R_A = (2^5−2^{-5})^{1/5} = 1.999609`; coefficient of `1/s` is `2 R_A = 3.999218`;
  `INT dz/(z^2+1/4)^2 = 4pi` and the second piece is exactly `pi/8 = 0.3926991`.  Confirmed.
* Instrument prefactor `3/(2pi)` re-derived from `dOmega_3 = 4pi sin^2 th dth`, `K = 3w_z/(8pi^2|w|^5)`.

## R4 THIRD INSTRUMENT — did not fire (`r2_psi_fd_instrument.py`)
Both of the build's numeric instruments (`toys/g2` and s4's `a_direct`) use the *analytic* `-d_z G_5`;
a wrong derivative constant would corrupt both together.  I therefore built a third instrument that
quadratures the POTENTIAL `psi1 = G_5 * eta` and takes `a = -d_z psi1` by central finite difference.
* axis identity `a(0,0) = 4.158883` vs `(1/2)log 4096`, rel **1.0e−08**;
* interior points vs the build's series: `+9.2e−4`, `−5.0e−4`, `+3.5e−4` (my exclusion ball / FD step);
* single-octave decay exponents **4.9955** (z-odd) and **4.0121** (constant) — the powers 5 and 4 confirmed
  without using any analytic z-derivative.  The brief's guessed `2^{-3j/2}` is indeed wrong.

## R5 FL-043 TWO-SIDEDNESS — FIRES (partial, undisclosed)
* K1 (series vs `g2`), K2, K3 (sympy identities), K6 (vs the sharp-clock refuter's elliptic-integral
  instrument): genuinely two-sided.  K2 in fact fired.  K6 is declared non-failing IN PREREG — disclosed.
* **K4 and the inner half of K5 are self-confirming.**  The "direct numerical evaluation" in `s5` is the
  *same* modal series with the *same* `g_l`; `C1` and `C2in` are the Cauchy–Schwarz/triangle majorants of
  exactly that series.  A majorant can fail to dominate its own series only through a coding typo.  These
  kills cannot falsify the shell representation, the kernel, `alpha_l`, or the derivative constants.  The
  measured slack (21×, 2–5×) makes them weak even as arithmetic checks.  NOTE §4's "every one was checked
  against an exact evaluation" overstates the independence; the word "exact" reads as an outside check.
  What K4/K5 *do* test: the Parseval normalisation and `|C_m| <= C_m(1)`.
* The **collar half of K5 is genuine** (real-space bathtub estimate vs the modal series — different routes).
* K3' is partially correlated (series and `a_direct` share the analytic `-d_z G_5`); my R4 closes that gap.
Resolution: rerun the split-check against an instrument that does not share the expansion — R4 above is
one; it agrees, so the substance is unharmed.  The disclosure, not the result, is the defect.

## R6 DECISION-RULE COMPLIANCE — FIRES (MAJOR)
PREREG: *"GATE H-c PASSES iff K1,K2,K3,K3',K4,K5,K5' all fail to fire.  K6 firing does not fail the gate."*
**K2 FIRED.**  The reported verdict is "GATE H-c PASSES on its substance".  By the registered rule the
gate FAILED.  What actually stands is a post-hoc repaired, *unregistered* re-run — legitimate science,
transparently ledgered, but it is not the pre-registered gate.  Verdict of record must read:
**GATE H-c FAILED as registered (K2); the repaired claim was re-verified outside pre-registration and
independently confirmed by this seat.**  (Estate laws: pin the rule, then obey it.)

## R7 THE CONSTANT `D` IS NOT A CONSTANT — FIRES (MAJOR) (`r3_D_constant.py`, `r3_out.txt`)
`s4` reports `D(all l) = 124.3991` for `|inner-octave strain| <= D M_in (rho_in/|x|)^5`, "valid for
`rho >= 4 rho_in`".  It is the l1 sum
`SUM_{l odd} |(l+1)g_l/(2l+3)| ((2^{l+4}−1)/(l+4)) ((l+2)(l+3)/2) * 0.5^{l-1}` truncated at LMAX=120.
Two errors compound:
1. the halving factor `0.5^{l-1}` corresponds to `rho >= 2 rho_in`, not the stated `rho >= 4 rho_in`;
2. with `0.5^{l-1}` the factor `2^{l+4}·0.5^{l-1} = 32` is **l-independent**, so the summand is `~8 l |g_l|`.
   Measured `|g_l| ~ l^{-1.4461}` (fit over l=400..2400) ⟹ summand `~ l^{-0.45}` ⟹ **the sum diverges**.
Measured `Dall(LMAX)`: **65.4 (40), 123.2 (120), 180.6 (240), 268.1 (500), 387.1 (1000), 560.3 (2000),
834.7 (4000)** — growing like `LMAX^0.55`.  `Dall(2000)/Dall(120) = 4.55`.  **124.40 is a truncation artifact.**
Corrected constants (both convergent, computed here):
* at the stated range `rho >= 4 rho_in` the correct l1 sum (factor `0.25^{l-1}`) is **`D = 13.4477`**,
  LMAX-independent from LMAX=40 to 4000;
* a **datum-independent** Cauchy–Schwarz constant, valid for ANY z-odd `|omega^theta| <= M_in` (the form the
  lemma actually needs — `D=124.4` was computed from the bang-bang `g_l` alone and is not universal):
  `D_CS = sqrt(2)·2^{5j}·INT_{2^{-j}}^{2^{1-j}} S(v) dv/v` = **14.67 (j=2), 13.65 (3), 13.59 (4), 13.58 (j>=5)**.
* exact measured worst value for the bang-bang octave: **12.2737** (at `phi -> 0`), consistent with
  `D(l=1 alone) = 12.4001`, which the build reports correctly.
Direction: 124.4 is *over*-conservative, so the lemma's inequality remains TRUE as stated — it is the
number and its derivation that are wrong.  Effect on §7's dynamical reading: `2^{5j} <= D` becomes
`2^{5j} <= 13.6`, i.e. `j <= 0.75` instead of `j <= 1`.  **The conclusion is strengthened, not weakened.**
(Had s4 run at LMAX=4000 it would have printed 835 and got `j <= 1.94` — still `j <= 1`.  So §H.3(a) is
robust to the error in both directions; only the printed constant is wrong.)

## R8 STATEMENT AUDIT — FIRES (two MINORs and one false additive find)
**(a) The ASSEMBLED SPLIT is false as written.**
`| a(x) − (M/2) log(R/(2 rho)) | <= C_tot(phi) M  for |omega^theta| <= M` — take `omega^theta ≡ 0`, which
satisfies the hypothesis for every `M>0`: LHS `= (M/2)log(R/2rho)`, unbounded, RHS `= C_tot M`.  NOTE §4
repairs this in the very next sentence ("with `(M/2)L_out` replaced by the exact
`INT_{2rho}^R (−(3/5)g_1) drho'/rho'`"); the BUILD RESULT's ASSEMBLED SPLIT paragraph drops the repair.
Correct form: `| a(x) − INT_{2rho}^{R} (−(3/5) g_1(rho')) drho'/rho' | <= C_tot(phi) M`, and separately
`| that integral | <= (M/2) log(R/(2rho))` with equality only for the bang-bang cap.  MINOR (the operative
"CONCRETE HANDOFF" form is correct).
**(b) exponent typo.**  NOTE §7 displays `D M_in 2^{-(j-1)(p)}` and then, two lines down, `124.40 M_in 2^{-5j}`.
`2^{-(j-1)p}` is 32× larger, so the display is the weaker (still true) statement, but it contradicts both
the following line and the BUILD RESULT's `2^{-pj}`.  MINOR.
**(c) additive find (3) is FALSE** (`r4_collar_phi.py`, `r4_out.txt`).  The build asserts: *"the 1/sin(phi)
in the proved collar bound is REAL, not an artifact of crude estimation — the exact c_edge genuinely blows
up at the axis."*  Two things are conflated.  `c_edge` is the INNER-EDGE plateau offset; the collar is a
different piece.  Measured, at `rho = 64` (the deep interior the lemma's `C_tot(phi)` bounds), the exact
two-sided collar **saturates** as `phi -> 0`:
| phi (deg) | 45 | 10 | 3 | 1 | 0.3 | 0.1 | 0.03 |
|---|---|---|---|---|---|---|---|
| exact collar | 0.1665 | 0.0582 | 0.0528 | 0.0523 | 0.0522 | 0.0522 | 0.0520 |
| bound `4/s+pi/8` | 6.05 | 23.4 | 76.8 | 229.6 | 764 | 2292 | 7640 |
| slack | 36× | 403× | 1455× | 4389× | 14630× | **146940×** |
The exact far remainder is likewise flat (`0.0134` at every `phi`), inside `C1=0.292` by 22×.
Even where growth IS real — the one-sided collar at `rho = rho0+`, which is what carries `c_edge` — the
exact rate is **logarithmic**: `0.577 (10°), 0.873 (3°), 1.147 (1°), 1.448 (0.3°), 1.726 (0.1°)`, fitting
`0.21702 log(1/phi)` (≈ `1/4 = kappa_0/2`, matching the build's own s6 fit `A = 0.25171`).  A `log(1/phi)`
is not a `1/sin(phi)`.  So: the `1/s` **is** an artifact of the bathtub estimate; the collar bound is
sound but over by up to 1.5e5; and there is no obstruction to a `phi`-uniform collar bound.
This *removes* what would otherwise have been the lemma's worst usability problem (a near-axis error
budget of `23.7 M`, exceeding the whole `(M/2)L` term until `L > 47` e-folds) — the object is small there;
only the estimate was bad.

## What survives, in full
Shell representation `(*)`; the interior/exterior modal expansions with `(l+2)` and `-(l+1)`;
`alpha_l = g_l/(2l+3)`; `N_l`; `kappa_0 = 1/2` per e-fold as the exact, sharp `l=1` interior mode with
equality only at `w = -M sgn z`; `C1 = 0.756945 / 0.291999`; `C2in = 0.026093 / 0.014754`; the collar bound
`(3.999218/s + pi/8) M` (true, but very slack, and the `1/s` removable); the multipole powers **5** (z-odd)
and **4** (general), now confirmed by a third instrument; `D(l=1) = 12.4001`; the closed form for
`c_edge(phi)` and its value `+0.2167733` at 10°, reproducing the sharp-clock refuter's `+0.216773` to
2.7e−7 and its whole finite-R table; the `c_edge` profile, its `(1/4)log(1/phi)` divergence for the
inadmissible bang-bang datum, and its saturation for admissible tapers; `kappa` = 0.497808 / 0.499721.
The build's own scoping ("NOT proved here: any dynamical statement") is correct and should be kept.

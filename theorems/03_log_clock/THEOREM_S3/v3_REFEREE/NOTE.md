# REFEREE REPORT on THEOREM (S3) version 3, and on `fix5`

Seat `s3close/round2/refute-THEOREM_S3_v3`, single referee, sitting of 2026-09-09.
Laws: `TORMENT NEXUS/LAWS.md` first 120 lines; seat header `TEMPLATES/SEAT_HEADER_2026-09-08.md`,
including the P1 checklist appended to it.
Posture: FL-000 will fall; find the move, and find the hole.

**VERDICT: NOT PROVED AS STATED.** Two MAJOR defects in the proof as written, one MAJOR defect in
the statement's quantifier, six MINOR. None of them is a mathematical obstruction, and only one
(R-5) moves a displayed number, in its seventh digit. With the three named repairs made, the grade
would be
**PROVED-MODULO-(interior parabolic Schauder; and the quadrature and sampled constants `C''`,
`Ghat`, `C_R`, `chat_a`, `L_Gamma^exist`, `C_E`, `K_2`)**, which is what the reach numbers
`L_* = 1424207.9721` and `log Lambda_* = 2848418.8578` actually rest on.

**The one thing that got better.** The report the fast vendor and this seat both attacked first,
`vartheta`, is now **settled in the theorem's favour**, by a certified upper bound rather than by
a convention. On the `4 sigma` ball,

```
     vartheta  <=  0.2030001      (f = 4;  0.2005878 / 0.2005659 / 0.2005657 at f = 8, 16, 32)
     vartheta  <=  0.2030156      uniformly over every f in [4, 32] and every L >= 40
```

against the theorem's stated `0.2623138022`. So `(H3'_vartheta)` at the stated `vartheta` is
**true**, and row `B4` may be graded PROVED for the right reason. The sheet's own reason -- a grid
lower bound multiplied by `SAFETY = 2` -- remains no reason at all, and it is what this report
grades.

Every number below came out of a script in `scripts/`, written from scratch here and run on
COPIES in a scratch directory; results in `results/`. Nothing under `round2/THEOREM_S3/` was
opened read-write; the target files are untouched. `SHA256SUMS` computed with `shasum -a 256`.

**FL-000 stands. Nothing here touches the headline problem.**

---

## 0. WHAT WAS RE-RUN, AND WHAT REPRODUCED

`fix5/x6_certified.py` and `fix5/x2_budget.py` were re-run from byte copies in the scratch tree.
Both reproduce to the last digit the sheet prints.

| quantity | `FIX5.md` / `THEOREM_S3_v3.md` | re-run here | |
|---|---|---|---|
| `N_sigma` | `[1.0124508487, 1.0124511502]` | `[1.012450848737, 1.012451150151]` | exact |
| `E_0` | `[7.5523054461, 7.5523079773]` | `[7.5523054461, 7.5523079773]` | exact |
| `Gfrak_0` | `[36.0783630557, 36.0784195640]` | `[36.0783630557, 36.0784195640]` | exact |
| `\|h_sigma''\|` proved | `2415.5503900804` | `2415.5503900804` | exact |
| `Psi_sigma(0)` | `<= 495.1064102419` | `495.106410` | exact |
| `r_h`, proved pad, `c = c_*` | `0.9033433307` | `0.9033433307` | exact |
| `min P_h'`, proved pad, at `fix4`'s cap | `-0.0029985709` | `-0.0029985709` | exact |
| certified window cap | `1.1647266909` | `1.1647266909` | exact |
| `L_*` proved, `eps <= cap` | `1424207.9721` | `1424207.9721111932` | exact |
| `L_*` proved, `eps <= 0.1` | `1554565.5693` | `1554565.5693347803` | exact |
| `L_*` measured, `eps <= cap` | `1761.8423` | `1761.8422827244044` | exact |
| `L_*` proved `C_a` measured, cap | `1379464.3531` | `1379464.3531055540` | exact |
| `kappa_delta`, `c_*`, `c_2`, `lam_max`, clock floor | sec.1.3 | recomputed in `r5` | agree to `<= 1e-9` |

The transcription is clean and the arithmetic is right. What follows is not a transcription
finding.

---

## 1. FINDING R-1, **MAJOR** (method; the conclusion survives): `vartheta` is a lower bound times two, and this report replaces it with a certified upper bound

### 1.1 The defect

`THEOREM_S3_v3.md` sec.1.4 assumes

> `(H3'_vartheta)` ... which holds with `vartheta = 0.2623138022 < 1` measured (a five-parameter
> grid lower bound on the supremum, inflated by `pmax-h3v` sec.B5's `SAFETY = 2`)

and row `B4` grades the row **PROVED** with `vartheta` **ENCLOSED-sampled**. `fix5/x1_vartheta.py`
computes `vartheta_k = sup r^{k+1}|d_v^k(eta_0 - eta_P)|/k!` by a grid search with local polish
over a five-parameter set. A grid search returns a value **attained**, hence a lower bound on a
supremum. Two times a lower bound is a larger lower bound. It is not an upper bound, and
`(H3'_vartheta)` asks for an upper bound. The sheet says so itself in `FIX5.md` sec.7 item 1 and
in `THEOREM_S3_v3.md` sec.3 item 3; saying so does not repair it. **This is the hypothesis of an
imported theorem (V.4'), and Theorem V.4' has no content without it**, so the defect is not
cosmetic: without an upper bound, row `F4`'s `eps_v` and therefore the whole budget carry no
warrant.

Symmetrically, the theorem's cliff narrative rests on the same convention in the other direction.
`THEOREM_S3_v3.md` sec.3 item 3 states "at an offset of `3.5 sigma` ... `2 vartheta = 1.76645` and
the hypothesis fails". What `x1` measured at `3.5 sigma` is `0.883225`, a **lower** bound, and
`0.883225 < 1`. A lower bound below one refutes nothing. The convention cannot both certify a
bound and refute one.

### 1.2 The certified upper bound, and how it is obtained

`scripts/r0_iv.py`, `scripts/r1_profile_cert.py`, `scripts/r2_theta_cert.py`, results
`results/r1_profile_results.json`, `results/r2_results.json`, `results/r4_results.json`.
No code is shared with `fix5` (L-14); the objects under test are re-implemented from the
mathematics.

Three reductions make the problem small enough to enclose.

1. `e := eta_0 - eta_P` depends on the point only through `rho = |x|` and the polar angle `phi`,
   and along the line `x + s v` it depends on the direction only through `a = <v, Yhat>` and
   `c = <v, e_z>`, because `|v| = 1` gives `rho(s)^2 = rho^2 + 2s(r a + z c) + s^2` and
   `z(s) = z + s c`. So the five-parameter search of `x1` is a four-parameter search.
2. **Rescale `s = rho tau`.** Then `rho(s) = rho Pp(tau)`, `r(s) = rho Rr(tau)`,
   `z(s) = rho Z(tau)` with
   ```
       Pp^2 = 1 + 2 p tau + tau^2 ,  p = a sin phi + c cos phi ,   Z = cos phi + c tau ,
       Rr^2 = sin^2 phi + 2 a sin phi tau + (1 - c^2) tau^2 ,
   ```
   all independent of `rho`, and
   ```
       e(s) = (1/rho) E(tau) ,  E = 1/Rr - Theta(log rho + log Pp) W_sigma(Phi)/(Pp N_sigma) ,
       vartheta_k = sup  sin^{k+1}(phi) | [tau^k] E | .
   ```
   `rho` survives only inside `Theta`, whose deviation from `1` on the ball is below `5e-05`, so
   the whole radial extent of the ball is carried in one interval and the search is
   **one-dimensional, in `phi`**. That is the step that makes an enclosure affordable: the
   dependency blow-up of interval arithmetic in the radial variable disappears because the
   radial variable is not really there.
3. The `tau`-coefficient of order `k` is a polynomial of degree `<= k` in `(a, c)`, carried
   exactly as a polynomial with interval coefficients, so the direction is not gridded in the
   main pass; the maximum over the closed unit disc is bounded by
   `sum |coefficient| x max |a|^i |c|^j` and then sharpened by interval evaluation over a
   `64 x 64` grid of `(a, c)` cells.

The one non-elementary input is the mollified profile. `r1_profile_cert.py` encloses it. In the
window `t in [-0.125, 1.25]` the reflected profile `Wtil` has exactly **one** corner, at
`t = delta`; on `(-delta, delta)` it is `t/(delta sin t)` and on `(delta, pi/2 - delta_m)` it is
`1/sin t`. Distributionally,

```
    W_sigma^{(m)}(phi) = sum_{i=1}^{m-1} J_i chi_sigma^{(m-1-i)}(phi - delta)
                         + int {Wtil^{(m)}}(t) chi_sigma(phi - t) dt ,
    J_i = P^{(i)}(delta) - Q^{(i)}(delta) ,
    J_1 = -58.528002 , J_2 = 889.128627 , J_3 = -20436.391287 , J_4 = 624475.727173 .
```

The jump sum is closed form. The integral is banded on a fixed `t`-grid of width `sigma/32`; on
the `1/sin` side the band's extremes are its endpoint values, because the derivatives of `csc`
alternate in sign with monotone magnitude on `(0, pi/2)` -- **proved symbolically here**, by
writing `(1/sin)^{(m)} = (-1)^m N_m(cos t)/sin^{m+1} t` and checking that every coefficient of
`N_m` is non-negative, for `m = 1..7` (`results/r1_profile_results.json`,
`sign_alternation_certificate`, all `True`). On the `t/(delta sin t)` side a Cauchy estimate on
the disc `|t| <= 1/2` bounds every derivative. Gaussian band masses come from the error function;
the tail outside the window is bounded by `sup|{Wtil^{(m)}}| x Phi(-16.7) < 1e-50`. A `phi`
INTERVAL enclosure is (value at the midpoint) `+/-` (halfwidth) `x` (the same bound one order up).

**Control (L-98).** Every enclosure contains `fix5/fx_profile.py`'s Gauss-Legendre value, an
instrument built on a different principle, at eleven test angles and every order `k = 0..4`:
`all_contained = true`. Relative enclosure widths at `phi = delta + 4 sigma` are
`3.0e-03` at `k = 0` rising to `2.1e-02` at `k = 4`.

### 1.3 The result

Outward-rounded interval arithmetic (one ulp on `+ - * /`, four on `sqrt`, `log`, `exp`,
`arctan`, `tanh`), a two-sided geometrically graded `phi` partition refined to `14800` boxes, and
`64 x 64` direction cells:

| offset | certified `vartheta` upper bound | `x1`'s grid lower bound | binding `k` | argmax `phi` |
|---|---|---|---|---|
| `3.5 sigma` | **`1.0558700`** | `0.883225` | 4 | `11.5107 deg` |
| **`4 sigma`** | **`0.2030001`** (f = 4) | `0.1305649` | 4 | `12.0837 deg` |
| `4 sigma` | `0.2005878 / 0.2005659 / 0.2005657` at `f = 8, 16, 32` | `0.1311569` | 4 | `12.0837 deg` |
| `4 sigma` | **`0.2030156` uniformly over `f in [4, 32]`** | -- | 4 | `12.0837 deg` |
| `5 sigma` | `0.1782419` | `0.102464` | 4 | `13.2322 deg` |
| `6 sigma` | `0.1558262` | `0.0884151` | 4 | `14.3755 deg` |

Per order at `4 sigma`, `f = 4`: `0.0129, 0.0194, 0.0541, 0.1075, 0.2030`. The lower and upper
bounds bracket consistently (`0.1306 <= 0.2030`), the ratio `1.55` being the cost of the
direction cells and the profile enclosure.

**Answers to the three questions the brief asks about `vartheta`.**

* **Does a certified `vartheta` exist, and does `2 vartheta < 1` survive?** A certified upper
  bound exists at `4 sigma`, `5 sigma` and `6 sigma`. There is no `2` to survive: with a genuine
  upper bound the `SAFETY` factor is not needed and must not be applied. The correct reading is
  `vartheta <= 0.2030156 < 1` at `4 sigma` uniformly in `f in [4, 32]`, hence
  `Q <= 1.5094606667` and `vt <= 0.2547303334` -- both **smaller** than the sheet's `1.7111799108`
  and `0.3555899554`, so the sheet's budget is conservative and nothing has to be recomputed.
* **At which offset?** `4 sigma` closes with a factor `4.93` of room to `1`. `5 sigma` and
  `6 sigma` close with more. **`3.5 sigma` does not close**: the certified bound there is
  `1.0559 > 1`, so at `3.5 sigma` the hypothesis is neither established nor refuted. The theorem's
  claim that it *fails* at `3.5 sigma` is therefore still not established -- what is established
  is that the certified bound crosses `1` between `3.5 sigma` and `4 sigma`. The cliff is real;
  the sentence describing it is not proved.
* **The third branch of `d`.** See R-4. It never binds, and no number changed.

### 1.4 What the certificate itself rests on, said plainly

So that the next referee has the same list this report demanded of `fix5`:

* the interval arithmetic is outward-rounded by `numpy.nextafter`, one ulp on `+ - * /` and four
  on `sqrt`, `log`, `exp`, `arctan`, `tanh`. It is not a formally verified library;
* the band extremes of `{Wtil^{(m)}}` on the `1/sin` side are endpoint values, justified by the
  symbolic sign-alternation certificate, and on the `t/(delta sin t)` side by a Cauchy estimate on
  `|t| <= 1/2`; both are evaluated in double precision from `sympy`-generated closed forms, with a
  relative pad of `1e-12` and an absolute pad of `1e-12` on every band sum;
* Gaussian band masses come from `scipy.special.ndtr`, whose error is far inside those pads;
* the truncation of the `t` window to `[-0.125, 1.25]` is covered by a tail bound below `1e-50`.

So the certificate is `PROVED-modulo-(double-precision evaluation of four explicit elementary
functions, padded at 1e-12)`. That is three orders of magnitude tighter than the margin it is used
at (`0.2030` against `1`, and against the stated `0.2623`), and it is the same standard the rest of
the sheet works to. An `mpmath.iv` re-run at 30 digits would remove even that; it was not needed to
decide the question.

### 1.5 Grade

The number in the theorem is right. The argument for it is not an argument. Because a certified
bound now exists and is below the stated value, the finding is **MAJOR (method)** and not FATAL:
`(H3'_vartheta)` at `vartheta = 0.2623138022` is true, and the correct sentence to write is
"proved by the interval computation of `refute-THEOREM_S3_v3/r2_theta_cert.py`", not "measured,
inflated by `SAFETY = 2`". `fix5` sec.7 item 1 and `THEOREM_S3_v3.md` sec.3 item 3 should be
struck and replaced.

---

## 2. FINDING R-2, **MAJOR**: the continuation step is asserted by a false sentence

### 2.1 The defect

`FIX5.md` sec.4.2 opens:

> Throughout, the contradiction hypothesis of (S3) is in force: `T_d(u_0) > tau`, `tau = c/(ML)`,
> and **`T_d <= T_*` by definition**, so by Proposition K' the solution is strong on `[0,tau]`.

and `THEOREM_S3_v3.md` sec.5 item 5 repeats it:

> a referee should satisfy themselves that `T_d > tau` does entail `tau < T_*`, which it does only
> because `T_d <= T_*` by definition.

`T_d(u_0)` is defined in sec.1.1 as *the first time `||omega(t)||_inf >= (3/2) M`*. If the
vorticity maximum never reaches `(3/2) M` before the maximal time `T_*`, there is no such time and
`T_d = +infinity`. `T_d <= T_*` is then false. And that is exactly the case the argument is in:
the contradiction hypothesis **is** the no-hit case. The step is vacuous precisely where it is
used, so the existence of the strong solution on the window `[0, tau]` -- which every estimate in
`u1`, `u2`, `hk2`, Lemma T'' and Theorem V.1 presumes -- is asserted and not proved. Under the
brief's rule that a PROVED which is a sketch is MAJOR, row `D7` and `FIX5.md` sec.4's closing
"Row `D7` is now **PROVED**, not PROVED-modulo" are overstated.

### 2.2 The dichotomy that has to replace it, written out

Let `T_* in (0, infinity]` be the maximal time of Proposition K'(b'), `tau := c/(ML)`, and
`tau' := min(tau, T_*)`.

1. **Contradiction hypothesis.** `||omega(t)||_inf < (3/2) M` for every `t in [0, tau')`.
2. **Strictness on compacts.** For each `s < tau'`, `[0,s]` is compact and
   `t -> ||omega(t)||_{L^inf}` is continuous, so `lamhat_s := sup_{[0,s]} ||omega||_inf/M < 3/2`
   and `theta_0(s) := (3/2 - lamhat_s)/2 > 0`. The supremum must be taken over a **closed
   subinterval**, not over `[0, tau]`, since `tau` is not yet known to be reached.
3. **Bootstrap.** `fix5` sec.4.2's closed-open-connected argument, run on `[0,s]` with
   `lamhat_s`, gives `Gamma(sigma) <= Gammabar - 2 theta_0(s) M L` for `sigma <= s`. `s < tau'`
   being arbitrary, `Gamma <= Gammabar` on `[0, tau')`.
4. **Continuation.** With `f = 0`, `[FENCES] sec.1.3` Step 3's estimate
   `d/dt ||u||_{H^s} <= C ||grad u||_{L^inf} ||u||_{H^s}` and Gronwall give
   `sup_{t < tau'} ||u(t)||_{H^s} <= ||u_0||_{H^s} e^{C Gammabar tau} < infinity` for every
   admissible `s in (5/2, 10.5)`; in particular `sup_{t<tau'} ||u(t)||_{L^inf} < infinity`
   (equivalently, from the energy identity K'(f') and
   `||u||_inf <= C ||u||_2^{2/5} ||grad u||_inf^{3/5}`). If `tau' = T_* < infinity` this
   contradicts K'(c'). Hence `tau' = tau` and `T_* > tau`: the strong solution exists on the
   **closed** window `[0, tau]`, and the three bounds extend to `t = tau` by continuity.
5. **The clock fires.** On the closed `[0, tau]`, `C6`'s quasimonotone comparison and `C7`'s
   conservative clock apply and force `lambda(t_*) >= 3/2` at some `t_* <= tau`, i.e.
   `||omega(t_*)||_inf >= (3/2) M`, contradicting 1.
6. **Conclusion.** The first hitting time exists and `T_d <= tau = t_*`. Note that the conclusion
   is at the **closed** right endpoint, so step 4's extension to `[0, tau]` is load-bearing and
   not a courtesy.

### 2.3 Where each continuity comes from, and what is still asserted

The brief asks that every strict inequality be established for the scalar suprema themselves, with
continuity in `t` of those scalars rather than pointwise continuity under a sup. Three of the four
are fine and one is asserted.

| scalar | continuity in `t` | status |
|---|---|---|
| `Gamma(s) = \|\|grad u(.,s)\|\|_{L^inf}` | `u in C([0,s];H^sigma)` by K'(b') and `H^sigma subset W^{1,inf}` for `sigma > 5/2`: the norm of a continuous curve | **proved**, and `fix5` says so |
| `\|\|omega(s)\|\|_{L^inf}` | same, via `H^{sigma-1} subset L^inf` for `sigma - 1 > 3/2` | proved; `fix5` uses it without naming the index condition |
| `m(s) = sup_x \|Phi_s(x) - Lambda_s(x)\|/\|x\|` | `fix5` gives a one-line reason (Lipschitz in `x`, continuous in `s`) | **asserted**; see below |
| `w(s) = sup_x \|log(J_Phi/lambda^2)\|` | same one-line reason | **asserted**; see below |

What is missing for the last two is uniformity in `x`, and it is available. `b(0) = 0` and
`||grad b||_inf <= Gammabar` give `|Phi_s(x)| <= |x| e^{Gammabar s}` and
`|d_s(Phi_s(x) - Lambda_s(x))| <= Gammabar|Phi_s(x)| + |x| |d_s lambda|/lambda <= C |x|`, so the
quotient is **Lipschitz in `s` with a constant independent of `x`**, and the supremum of a family
that is uniformly Lipschitz in `s` is Lipschitz in `s`. Likewise
`d_s log J_Phi = (div b)(Phi_s(x),s)`, bounded by `5 Gammabar`, and
`d_s log lambda^2 = 2 frak_a`, bounded by `u1` (P2); so `w` is Lipschitz in `s` uniformly in `x`.
Two lines, and they are the difference between "continuous because each term is" and "continuous".
Graded MINOR inside this MAJOR: R-9.

The strictness itself is fine. `fix5` sec.4.3 establishes `Gammabar - Gamma(0) >= 0.5164262397 M L`
at `s = 0`, `mu(c) <= 1.65e-01 < 1/2` and `mu_J(c) <= 2.53e-03 < 1/2` across the whole table, and
the propagation margin `2 theta_0 M L > 0` from `lamhat < 3/2`. The vendor's worry that the closed
inequalities are never made strict is **not** realised; the worry that the interval is not known to
exist **is**.

---

## 3. FINDING R-3, **MAJOR**: the statement quantifies `(H3'_vartheta)` over every `f >= 4`, and it is false there

`scripts/r3_checklist.py`, results `results/r3_results.json` section C.

The statement of sec.1.4 reads "Let `u_0` be the datum of sec.1.1, with `nu > 0` and `f >= 4`",
and sec.1.3 lists `vartheta = 0.2623138022` as "the largest over `f >= 4`", on the evidence of
four values `f = 4, 8, 16, 32`. The budget minimises over exactly those four
(`p4_budget.FS_F4 = [f for f in t3_budget.FS if f >= 4]`), so nothing in the numerics needs more.
But the hypothesis as quantified is false. At `x = x_*` the `k = 0` clause reads

```
      | 1 - Theta(log rho_*) A_sigma(phi_0)/N_sigma |  <=  vartheta ,
      A_sigma(phi_0)/N_sigma = 0.9890911210 ,   rho_* = (1+f) rho_0 .
```

`Theta` is the radial ramp, and it turns off at the outer edge `u = L`. At the theorem's own
`L = L_* = 1424207.9721`:

| `u = log rho_*` | `Theta` | `vartheta_0` at `x_*` | exceeds the stated `vartheta` |
|---|---|---|---|
| `L - 1` | `0.997527` | `0.013355` | no |
| `L - 0.5` | `0.880797` | `0.128811` | no |
| `L` | `0.119203` | `0.882097` | **yes** |
| `L + 1` | `4.53979e-05` | `0.999955` | **yes** |

The threshold is `log(1+f) > L - 0.384557`. Beyond it `(H3'_vartheta)` fails at the stated
`vartheta`, and as `f` grows further `vartheta_0 -> 1`, at which point Theorem V.4' has no content
at all. The tracked point has simply left the shell.

**Correction.** The statement must bound `f` above. The honest form is
`f in {4, 8, 16, 32}` (what the budget uses), or `4 <= f <= 32` together with the uniform
certificate of R-1 (`vartheta <= 0.2030156` for every `f in [4, 32]`, proved here), or a named
`F_max(L)` with `log(1 + F_max) <= L - 1`. As written, the theorem asserts a hypothesis over a
range on which it is false.

---

## 4. FINDING R-4, MINOR: the third branch of `d` is dimensionally wrong as printed

`THEOREM_S3_v3.md` sec.1.2 and `FIX5.md` sec.1.3 print

```
    d  :=  rho_* min( sin(phi_0 - delta - 4 sigma), sin(pi/2 - delta_m - phi_0 - 4 sigma), f/rho_* ) .
```

The first two entries are dimensionless; the third is not, and `rho_* x (f/rho_*) = f` is a pure
number where a length is wanted. `pmax-h3v` sec.B1 and ASSEMBLY have
`d = min(f, rho_* sin(phi_0-delta), rho_* sin(pi/2-delta_m-phi_0)) rho_0`, so the inner branch is
`f rho_0` and the nondimensional entry is **`f rho_0/rho_*`**, as the fast vendor says.

**No number moves**, and `results/r3_results.json` section A says why: `x1.d_rule` and
`x2.d_of` both compute `min(rho_* sin(.), rho_* sin(.), f)` in the units `rho_0 = 1` that every
script uses, which is the correct `f rho_0`. The taper branch binds whenever
`f >= sin(phi_0-delta-4 sigma)/(1 - sin(phi_0-delta-4 sigma)) = 0.4443101761`, hence at every
`f >= 4`, so the inner branch never binds and the misprint is inert. Graded MINOR, correction to
the printed formula only.

---

## 5. FINDING R-5, MINOR: the energy-Reynolds display is stronger than the bound that is proved

`scripts/r5_statement_checks.py`, results `results/r5_results.json` section a.

The theorem proves `M T_d <= c_*(1 + eps)/L` and then displays

```
     T(Lambda) <= c_2 (1 + eps(Lambda))/log Lambda ,   c_2 = 2 c_* ,  log Lambda = 2 L + 2.9135781820 .
```

Since `c_*(1+eps)/L = c_2(1+eps)/(2L)` and `2L = log Lambda - 2.9135781820 < log Lambda`, the
displayed inequality is **strictly stronger** than the proved one, by the factor
`log Lambda/(log Lambda - 2.9135781820)`. At `log Lambda_* = 2848418.8578` that factor is
`1.0000010229`, so the displayed `eps <= 0.1647266909` has to be read as
`eps <= 0.1647278823`. At the measured column's `log Lambda_* = 3526.5981` the factor is
`1.0008268` and `eps <= 0.1647266909` becomes `eps <= 0.1656898`. The correction is below the
precision anyone will use, but the inequality as displayed does not follow from the inequality that
is proved, and the fix is one symbol: divide by `log Lambda - 2 log s - (2/5) log C_E`, or raise
`eps` by `(1+eps) x 2.9135781820/(log Lambda - 2.9135781820)`.

---

## 6. FINDING R-6, MINOR: the Sobolev bookkeeping, re-derived

The fast vendor's arithmetic is right and the chain does not depend on it. Recorded so nobody has
to redo it.

**The embedding.** On `R^5`, `H^s subset C^{k,alpha}` needs `s > k + alpha + 5/2`, so
`H^s(R^5)` with `s < 9.5` gives `C^{6,alpha}` for every `alpha < 1` and **not** `C^7`. The vendor
is right about that, and right that `eta_0`, the five-dimensional field, is homogeneous of degree
`7` at the origin (`omega_0^theta ~ rho^8 sin phi`, divided by `r = rho sin phi`), whose Fourier
transform on `R^5` is homogeneous of degree `-12`, giving `H^s(R^5)` exactly for `s < 9.5`.

**But that is not the statement `fix5` makes.** `FIX5.md` sec.3.2 and sec.3.4 work on `R^3`:
`omega_0 in H^s(R^3)` for `s < 9.5`, `u_0 in H^s(R^3)` for `s < 10.5`, and
`H^s(R^3) subset C^k` for `s > k + 3/2`, so `u_0, u(t) in C^{8,alpha}`, `alpha < 1`. That is
correct as written; there is no `C^7` claim from `R^5` in the sheet to correct.

**The highest derivative any consumer uses, re-derived from the differentiated hypothesis.**
Read at source rather than from `fix5`'s table:

| consumer | object | order in `u` |
|---|---|---|
| `u2` sec.2, `\|\|grad u\|\|_op` bound | `grad u`, `omega` | 1 |
| `u2` sec.4, `Gamma_rad` | `grad b` | 1 |
| `u2` sec.5 **P1** | `eta` a classical solution of `D_t eta = nu Lap_5 eta` | 3 |
| `u2` sec.5 **P2** | `grad eta` a classical solution of its own equation | **4** |
| `u2` sec.6, Riesz composition | `grad_5 eta` | 2 |
| `pmax-h3v` Part A barriers and comparison | as P1/P2 | **4** |
| `hk2` `(H-K2)`, `K_2 = \|\|grad^2_5 b\|\|_inf`, `sup_B \|\|grad^2 eta\|\|_F` | `grad^2 b`, `grad^2 eta` | 3 |
| Lemma T''/(H1),(S), and `lambda(.,s) in C^1` | `sup rho\|eta\|`, `sup rho^2\|grad eta\|` | 2 |
| Theorem V.1, Feynman-Kac | `eta` a bounded classical solution | 3 |
| the continuity/first-crossing step | `Gamma`, `\|\|omega\|\|_inf` | 1 |
| **Theorem V.4' hypothesis `(H3'_vartheta)`** | `d_v^k(eta_0 - eta_P)`, `k <= 4`, **on the datum** | **5** |

The vendor's caution -- that a fourth-order remainder applied to `D^k grad eta` would reach `k+5`
-- does **not** apply: `pmax-h3v` sec.B2's proof of V.4' applies the fourth-order Taylor remainder
to `eta_0` itself (step (iii): `k = 4` gives `sup|d_v^4 eta_0| <= (1+vartheta) 4! M/R_-^5`), never
to `grad eta_0`. Read at source. So the ceiling is **four on the reference solution** and **five
on the datum**, and both are conditions at a distance `rho_* - d = 3.4618602587 rho_0` from the
origin where `eta_0` is real-analytic. `C^{6,alpha}` on the window would suffice with two orders
to spare; `C^{8,alpha}` from `s < 10.5` on `R^3` suffices with four.

**Is the reference solution's `H^s` norm bounded on `[0, tau]`?** Yes, by K'(b') plus compactness
of `[0, tau]` -- *provided* `tau < T_*`, which is R-2. The vendor is right that analyticity of the
datum on the ball supplies nothing about `t > 0`; what supplies it is K'(b'), and K'(b') is only
usable on a window that has been shown to exist.

**What is genuinely unstated.** Every consumer above lives on `R^5` (`b`, `eta = omega^theta/r`,
`Lap_5`), while the regularity is established on `R^3` for `u_0` and `u(t)`. The transfer -- that
the five-dimensional lift of an axisymmetric no-swirl field, and `omega^theta/r`, inherit the
required regularity **across the axis** `r = 0` -- is standard (it is the evenness of the
axisymmetric profiles in `r`, already used in sec.1.1 to say `eta_0` is `C^infinity` in the angle
including the axis) but it is nowhere written as a step, and the consumer table's column heading
"highest order in `u`" quietly performs it. One paragraph, in `FIX5.md` sec.3.4. Graded MINOR.

---

## 7. FINDING R-7, MINOR: `K_2` is re-checked at the wrong `L` and at the superseded ball

`fix4/p3_results.json` `C_K2_mollified_f4_L10`: the mollified-datum re-check of the imported
majorant `K_2 <= 161.7735 M/rho_0` was run at `f = 4`, **`L = 10`**, four values of `lambda`, and
with `d = 0.5842755976`, `d_layers = 1.9134171618` -- ASSEMBLY's **tangent** ball, the one unit
F-1 replaced. `max over lambda = 40.5556`, slack `3.9889277977`. The theorem runs at
`L >= 1424207.9721` and on the ball `d = 1.5381397413`.

The direction is safe on the geometry (`K_2` is a supremum over a neighbourhood built from `d`, and
the corrected `d` is smaller, so the corrected `K_2` cannot exceed the one that was checked), and
the referee's own `C_K` sweep past four orders of magnitude moves `L_*` only in the fourth digit.
But "re-checked in `fix4/p3`" in row `B6` is doing more work than the check supports: it is a
four-point sample at an `L` five orders below the theorem's, on the superseded ball. Category, for
the P1 checklist: **imported from another datum, sampled re-check at a non-representative `L`,
sensitivity-swept.**

Related mislabel: row `B6` ends "not load-bearing". The *value* is not load-bearing. The
*finiteness* is: V-b's `(H1)` asks for `b in C^{1,1}`, which is `K_2 < infinity`, and without it
Theorem V.4' does not apply at all. Say which.

---

## 8. FINDING R-8, MINOR: "`eps(L)` decreasing to the clock floor" is in the statement and is a tabulation

`results/r5_results.json` section b. The displayed conclusion of sec.1.4 asserts that `eps(L)`
decreases to `(1 - 2 kappa_delta)/(2 kappa_delta) = 1.6700567265e-02`; sec.3 item 7 records that
this is tabulated at five values of `L`. Tabulated here at fifteen more, with `fix5/x2`'s own
instrument on the certified constants:

```
     L = 2e+06  eps = 0.0668630   L = 1e+07  eps = 0.0238848   L = 1e+09  eps = 0.0167685
     L = 3e+06  eps = 0.0449111   L = 3e+07  eps = 0.0190057   L = 1e+10  eps = 0.0167074
     L = 5e+06  eps = 0.0320098   L = 1e+08  eps = 0.0173834   L >= 1e+12 eps = 0.0167006777
```

Monotone decreasing throughout, and the plateau `1.6700677660e-02` sits `6.6e-06` in relative terms
above the stated floor `1.6700567265e-02` -- inside the window search's own bisection resolution
(`nbis = 16` over a window of width `0.1647 c_*` resolves `2.5e-06` in `c/c_*`). So the claim is
supported and is still not proved, and the floor is printed to eleven significant figures by an
instrument that resolves six. Either prove the monotonicity or weaken the sentence to
"`eps(L) -> (1-2 kappa_delta)/(2 kappa_delta)`, tabulated".

---

## 9. FINDING R-9, MINOR: the two flow suprema are asserted continuous

Folded into R-2 sec.2.3 above, with the two lines that repair it.

---

## 10. THE P1 CHECKLIST, ITEM BY ITEM

### 10.1 Every constant recomputed for this datum, or proved profile-independent

| constant | category |
|---|---|
| `kappa_delta`, `c_*`, `c_2`, `lam_max`, clock floor | **recomputed** for the mollified profile (`fix4/p1`); reproduced here to `1e-9` |
| `N_sigma`, `sup\|W_sigma\|`, `E_0`, `Gfrak_0`, `sup\|h'\|`, `\|h''\|`, `Psi_sigma(0)`, `r_h`, `min P_h'`, window cap | **recomputed and certified** (grid maximum plus a proved analytic pad), `fix5/x6`; re-run here exactly |
| `vartheta` | **recomputed** for the mollified datum and the corrected ball, but as a **grid lower bound**; certified here (R-1) |
| `C_E` | recomputed (`fix4/p2`); a convergence statement, `7.3e-06` in `lmax`, not an enclosure |
| `N = \|\|eta_0\|\|_inf r_*/M` inside the viscous budget | **not recomputed**: `t3_budget` uses the kinked `1/sin delta = 7.6612975755` against the certified `E_0 <= 7.5523079773`. **Proved conservative**: `E_hess` and `E_tail` are affine increasing in `N` |
| `C_R`, `chat_a`, `Ghat`, `C''`, `L_Gamma^exist` | quadratures and sampled suprema, recomputed on the mollified profile but **not enclosed** |
| `K_2` | **imported for another datum**; sampled re-check at `L = 10` on the superseded ball (R-7) |
| far and near kernel constants | Parseval evaluations in the source seat; the three hypotheses are **profile-free** and re-verified |
| `vartheta`'s independence of `L` | **proved profile-independent**: `Theta_L - Theta_inf` and all its `u`-derivatives are bounded by `4 e^{2(u - L + eps_r)/eps_r}`, which at `f = 32`, `L = 40` is `3.77e-125` and at `L = L_*` underflows to zero in 40-digit arithmetic. The same estimate is what fails at `f ~ e^L` (R-3) |
| interior parabolic Schauder | **cited, not proved**; the one textbook citation |

### 10.2 Every extremum analytic or an outward-rounded enclosure

`fix5/x6`'s pads are **not** bare `np.max`: each is a grid maximum plus an explicitly proved pad,
`(1/8) h^2 sup|f''|` at an interior critical point and `(h/2) sup|f'|` where the maximised function
is only piecewise smooth, with the `sup` taken from the proved `||W_sigma^{(k)}||_inf <=
||Wtil'||_inf ||chi_sigma^{(k-1)}||_1` chain and the closed-form `L^1` norms. That is the right
shape of argument and it checks out. Three observations:

* the `d2` clause needs the maximiser to be an interior critical point. It is: `N_sigma`'s argmax
  is at `9.9885 deg`, `sup|W_sigma|`'s at `175.5720 deg`, `Gfrak_0`'s branch maximum at
  `170.9532 deg` with `8S - 81P = 6844.05`, about a hundred grid steps from the nearest cell where
  the branch could switch. The prose says so; the script asserts `argmax_in_switch_cell = False`
  but does not assert the distance. Add the assertion.
* the arithmetic is plain double precision with no outward rounding. The pads are many orders
  above one ulp, so this is not a defect, but "certified" and "outward-rounded" are not the same
  word and the sheet uses the first.
* the profile values that the pads are applied to come from a 24-point Gauss-Legendre convolution
  with no quadrature-error term. The independent enclosures of `r1_profile_cert.py` contain them
  at every test angle, which is the missing check, supplied here.

`vartheta` was the one extremum that was a grid value where a certified bound was claimed. R-1.

### 10.3 Imported hypotheses matched to the statement

Theorem V.4' (`pmax-h3v` sec.B2) assumes V-b's `(H1)`, `(H2)`, `(H4)` and `(H3'_vartheta)`.

| hypothesis | what it says | discharged by | in the statement? |
|---|---|---|---|
| V-b `(H1)` | `b` is `C^{1,1}` in `x`, continuous in `t`, `\|b\| <= B_0(1+\|x\|)` | `(H-K2)` for `grad^2 b` bounded, `u2` sec.4 `(D2)` and `Gamma <= Gammabar` for the growth | no, and it should be named: it is the qualitative half of `K_2` |
| V-b `(H2)` | `eta` a bounded classical solution on `R^5 x [0,tau]` | Theorem V.1 plus Proposition K'(e') **and R-2's continuation** | no |
| V-b `(H4)` | `nu > 0`, `tau > 0` | trivially | yes, `nu > 0` |
| `(H3'_vartheta)` | the five clauses `k = 0..4` | R-1, now certified | yes, as an assumption |
| `u2`'s `(H4)` (the bootstrap posture, a different `(H4)`) | `Gamma <= Gammabar` on `[0,s]` | `FIX5.md` sec.4, modulo R-2 | declared not assumed |
| Prop K' | `u_0 in H^s`, `s < S_1 = 10.5`, `div f = 0`, `f` obeying (5) | sec.3.2, `f = 0` | yes, by citation |
| `[FENCES]` Lemma A / sec.1.3 Step 3 | `f in L^1 cap L^2(H^s)` -- (A1) | `f = 0` | not named, needed by R-2 step 4 |
| interior parabolic Schauder | classical regularity for P1/P2 | cited | named in sec.3 item 1 |

Two `(H4)`s with different meanings, four sheets apart, is a trap for the next reader. Rename one.

### 10.4 The window of the estimates versus the interval of the conclusion

The estimates hold on `[0, tau]`, `tau = c/(ML)`; the conclusion asserts `T_d <= t_* = c/(ML)`.
They are the same interval, and the conclusion is at its **closed right endpoint**. That is
admissible, and it is exactly why R-2 step 4 has to deliver existence on the closed interval and
not merely on `[0, tau)`. `c` ranges over `[c_*, 1.1647266909 c_*]`, inside the certified monotone
cap, consistently in the statement, in `FIX5.md` sec.6.4 and in the re-run.

### 10.5 Geometric containments, checked by script at the actual parameters

`results/r3_results.json` section B. At `f = 4, 8, 16, 32`, and scale-free for every `f`:

```
    d/rho_*        = sin(phi_0 - delta - 4 sigma) = 0.3076279483
    R_-/rho_*      = sin(phi_0) - sin(phi_0 - delta - 4 sigma) = 0.1923720517   > 0
    d < r_*        <=>  0.3076279483 < 0.5                                       true, every f
    z_min > 0      <=>  cos(phi_0) = 0.8660254038 > 0.3076279483                 true, every f
    phi_min on the ball = delta + 4 sigma exactly                                every f
    phi_max on the ball = 2 phi_0 - delta - 4 sigma = 47.9163 deg                every f
    clearance to the equatorial kink pi/2 - delta_m = 85 deg : 32.4 sigma        every f
    equatorial branch never binds: sin(pi/2-delta_m-phi_0-4 sigma) = 0.7706949701 > 0.3076279483
    inner branch never binds for f >= 0.4443101761
```

At `f = 4`: `d = 1.5381397413`, `R_- = 0.9618602587`, `rho_min = 3.4618602587`,
`rho_max = 6.5381397413`, `z_min = 2.7919872777`. All the containments the chain needs hold, and
they hold for the algebraic reason, not just at the four sampled `f`. The one containment that
fails is the one nobody checked: the ball must also sit inside the shell's radial plateau, and
that is R-3.

---

## 11. THE RESIDUAL LIST'S HONESTY, AND THE GRADE OF THE STATEMENT

The brief asks whether anything labelled `enclosed` or `proved` is a quadrature, and whether the
statement depends on it at PROVED level.

**Yes, on both counts, and the theorem's own tables mostly say so.** `THEOREM_S3_v3.md` sec.2 marks
`A5` and `E4` and `D7`'s `C''` **ENCLOSED-sampled**, `G1` **COMPUTED**, and `G3` **CONDITIONAL**;
sec.3 items 2, 3, 4, 5, 6 list the offenders. That is honest bookkeeping and it is better than
`v2`'s. What is not honest is the boxed statement of sec.1.4, which displays

```
     L_* = 1424207.9721 (all proved) ,  log Lambda_* = 2848418.8578
```

with no tag, while `L_*` is a monotone function of `C_R`, `C''`, `Ghat`, `chat_a`,
`L_Gamma^exist`, `C_E` and `K_2`, of which

* `C_R` is a sampled supremum over the angle (six resolutions, `3e-08`), not an enclosure;
* `C''`, `Ghat`, `chat_a`, `L_Gamma^exist` are quadratures with a sign-change existence proof;
* `C_E` is a convergence statement;
* `K_2` is imported for another datum (R-7);

and the whole P1/P2 route rests on a cited interior parabolic Schauder estimate. The column label
"all proved" refers to which *branch* of the budget is taken (proved rather than measured `C_a`),
not to the grade of the constants, and a reader will not know that.

**So the statement's grade is `PROVED-MODULO-(interior parabolic Schauder; C'', Ghat, C_R,
chat_a, L_Gamma^exist, C_E, K_2)`,** and the reach numbers should carry that tag wherever they are
displayed, including in the abstract and in sec.4's BFG paragraph. Sec.4 already says "conditional
on the residual list of sec.3 and on nothing else"; the statement should match it.

**Does the certified `vartheta` change any of this?** No. `scripts/r4_uniform_and_budget.py` sweeps
`vartheta` through the budget on the certified constants with `fix5/x2`'s own instrument, imported
unchanged:

```
     vartheta = 0        L_* proved, certified cap 1424207.9721111932  eps<=0.1 1554565.5693347803  measured 1758.7135233
     vartheta = 0.203    L_* proved, certified cap 1424207.9721111932  eps<=0.1 1554565.5693347803  measured 1760.9683838
     vartheta = 0.262314 L_* proved, certified cap 1424207.9721111932  eps<=0.1 1554565.5693347803  measured 1761.8606704
```

(the last row uses `0.262314` at every `f`; `fix5/x2`'s own run, which uses `0.2611298` at `f = 4`
and `0.2623138` above, gives `1761.8422827` for the same column, and the difference is the
per-`f` value, not the instrument. The sweep above `vartheta = 0.2623` did not finish inside this
sitting and is not claimed; `stdout/r4_stdout.txt` is the receipt for the three rows that did.)

The proved column does not move at all between `vartheta = 0` and the stated value; only the
measured column moves, by `1.8e-03` relative, because it sits near the clock floor where `eps_v`'s
share is largest. `eps_v` is between `1.5e-07` and `2.8e-07` of `eps`. Carrying `Q` and `vt` was
the right thing to do and it changes nothing in the reach, and neither does replacing the stated
`vartheta` by the certified one; the certified value being the smaller, the substitution is in the
safe direction in any case.

---

## 12. CORRECTED STATEMENT

Only the changed clauses are given; everything else is as `THEOREM_S3_v3.md` sec.1.4.

> **THEOREM (S3), version 3, corrected.** Let `u_0` be the datum of sec.1.1 with `nu > 0` and
> **`f in {4, 8, 16, 32}`** (or any `f` with `4 <= f <= 32`), and let `u` be its maximal strong
> solution, which exists and is unique by Proposition K'. Then, on the window
> `[0, tau]`, `tau = c/(ML)`:
>
> * **`(H3'_vartheta)` is not assumed. It is proved**, with
>   `vartheta = 0.2030156 < 1` on the ball `B(x_*, d)`,
>   `d = rho_* min( sin(phi_0-delta-4 sigma), sin(pi/2-delta_m-phi_0-4 sigma), f rho_0/rho_* )`,
>   by the interval computation of `refute-THEOREM_S3_v3/scripts/r2_theta_cert.py`, uniformly
>   over `f in [4,32]` and `L >= 40`. Hence `Q <= 1.5094606667` and `vt <= 0.2547303334`; the
>   budget's `Q = 1.7111799108`, `vt = 0.3555899554` are conservative and may be kept.
> * `L >= L_Gamma^exist` as before.
> * The solution exists on the closed window: this is the dichotomy of sec.2.2 of this report
>   (bootstrap on `[0, min(tau, T_*))`, Gronwall, Proposition K'(c')), **not** the sentence
>   "`T_d <= T_*` by definition", which is false for the case at hand.
>
> Then for every `L >= L_*` there is `c in [c_*, 1.1647266909 c_*]` with `c >= c_*(1 + eps(c,L))`
> and
> ```
>     T_d(u_0)  <=  t_*  =  c/(M L) ,
>     L_*  =  1424207.9721  (proved branch) ,  1379464.3531 (C_a measured) ,  1761.8423 (measured)
>                                                      at  eps <= 0.1647266909 ,
> ```
> **all three grades being PROVED-MODULO-(interior parabolic Schauder; `C''`, `Ghat`, `C_R`,
> `chat_a`, `L_Gamma^exist`, `C_E`, `K_2`)**. In energy Reynolds number the same bound reads
> ```
>     M T_d  <=  c_2 (1 + eps)/(log Lambda - 2.9135781820) ,
> ```
> equivalently `c_2 (1 + eps')/log Lambda` with `eps' = (1+eps) log Lambda/(log Lambda -
> 2.9135781820) - 1 = 0.1647278823` at `log Lambda_* = 2848418.8578`.
> `eps(L)` is tabulated, not proved, to decrease to `(1 - 2 kappa_delta)/(2 kappa_delta)`.

---

## 13. THE SINGLE NEXT UNIT

**Write the continuation dichotomy of sec.2.2 as a numbered lemma in the `fix5` line, with the
two uniform-Lipschitz continuity arguments of sec.2.3, and re-issue sec.4 of `FIX5.md` with the
sentence "`T_d <= T_*` by definition" struck.** It is the only item on this list that is a hole in
the proof rather than a defect in the writing or in a quantifier: `(H3'_vartheta)` is now certified
(R-1), the `f` range is a one-line edit (R-3), and everything else is bookkeeping. Nothing in the
budget moves, so the unit is pure prose plus one script that exhibits `Gammabar`,
`sup_{t<tau'}||u||_{H^s}` and the Gronwall constant at the theorem's own `L`. Estimated cost: one
sitting.

Second, if there is appetite: state the theorem at `6 sigma`. The certified `vartheta` there is
`0.1558262`, `L_*` is identical to the last displayed digit, and the geometric margin stops being
a quarter of a `sigma`. `FIX5.md` sec.7 item 8 already says the brief is the only reason it is not.

---

## 14. WHAT THIS REPORT CANNOT CATCH

It certifies one constant, re-runs two scripts, checks the geometry and the quantifiers, and reads
four hypotheses at source. It does not check that `C_R` is the right constant, that `u2`'s P1 and
P2 are the right instruments, that the far and near kernel constants are right, that the clock
argument is right, or that BFG's Theorem 10 says what sec.4 says it says. `C_R`, `C''`, `Ghat`,
`chat_a` and `L_Gamma^exist` remain sampled suprema and quadratures; this report did not enclose
them, and the same objection that R-1 makes about `vartheta` applies to `C_R` word for word --
`C_R`'s supremum over the angle is a sampled search, and the theorem uses it as a bound. That is
the obvious next `vartheta`-shaped target, and the machinery in `scripts/` is most of what it
needs.

**FL-000 stands. Nothing here touches the headline problem.**

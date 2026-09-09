# THEOREM (S3), version 3: the sharp logarithmic doubling clock, on a smooth datum, with the Theorem V.4' ball off the taper cone

Seat `s3close/round2/THEOREM_S3/fix5`, sitting of 2026-09-08. This file restates the theorem after
the six findings of `round2/refute-THEOREM_S3_v2/NOTE.md`. It supersedes nothing:
`THEOREM_S3.md`, its three addenda and `THEOREM_S3_v2.md` are unchanged and remain the record of
the earlier versions. Every constant below comes from a script in `fix5/`, listed in
`fix5/FIX5.md` and re-asserted by `fix5/check_fix5.py`.

**What changed against `THEOREM_S3_v2.md`.** Five things, one of them structural.

1. The Theorem V.4' ball is moved four mollification widths off the taper cone:
   `d = rho_* sin(phi_0 - delta - 4 sigma)` in place of `rho_* sin(phi_0 - delta)`. On the tangent
   ball that `v2` used, `(H3'_vartheta)` is **false** for the mollified datum, with
   `vartheta >= 8.78` at every admissible `f`; on the corrected ball it holds with
   `vartheta = 0.1305648864` at `f = 4`, and `2 vartheta = 0.2611297728 < 1` at `SAFETY = 2`.
2. The budget carries `E_vartheta` and `Q = (1 + vartheta)/(1 - vartheta)`, so the `eps_v` the
   columns report is Theorem V.4's successor V.4' and not V.4, with the reduced `d` propagated
   everywhere `d` enters.
3. Proposition K is replaced by Proposition K', stated at finite Sobolev index, which the datum
   satisfies; the Schwartz hypothesis it carried is not available here and is not needed.
4. The bootstrap posture `(H4)` is closed by a continuity and first-crossing argument rather than
   assumed, and Lemma T' is restated as Lemma T'' on `[1, lam_max]` with the `C^1` regularity of
   the reference profile proved rather than asserted.
5. Every entry of `v2` sec.3's ENCLOSED list is replaced by a certified bound. That costs one
   thing: the certified monotone window cap falls from `eps <= 0.1664585076` to
   `eps <= 0.1647266909`. It also uncovers a first-order discretisation error in `fix4`'s own
   `W_sigma'`, which made `Gfrak_0` `3.2e-05` too large; the correction is favourable.

**FL-000 stands.** This is not a proof of the headline problem and does not bear on it.

---

## 1. THE STATEMENT

### 1.1 The datum

Unchanged from `THEOREM_S3_v2.md` sec.1.1. Fix `M > 0`, `nu > 0`, `delta = 7.5 deg`,
`delta_m = 5.0 deg`, `phi_0 = 30.0 deg`, `s := 1/sin delta = 7.66129757554039`, `eps_r = 0.25`,
`sigma = 0.02` rad, `f >= 4`, `L > 0`. Write `u := log(rho/rho_0)`,
`phi_ax := min(phi, pi - phi)`. Put

```
    nu     := M rho_0^2/s^2 ,     R := rho_0 e^L ,
    Theta(rho) := (1/2)[ tanh((u - eps_r)/eps_r) - tanh((u - L + eps_r)/eps_r) ] ,
    A(phi)  := sgn(cos phi) min(1, phi_ax/delta) min(1, |phi - pi/2|/delta_m) ,
    W(phi)  := A(phi)/sin phi ,   Wtil := W even-reflected at 0 and pi, 2 pi periodic ,
    chi_sigma(y) := exp(-y^2/(2 sigma^2))/(sigma sqrt(2 pi)) ,
    W_sigma := (Wtil * chi_sigma)|_{[0,pi]} ,   A_sigma := W_sigma(phi) sin phi ,
    N_sigma := sup_{[0,pi]} |A_sigma| in [1.0124508487, 1.0124511502] ,
    omega_0^theta(rho, phi) := - M Theta(rho) A_sigma(phi)/N_sigma ,
    eta_0 = omega_0^theta/r = -(M/rho) Theta(u) W_sigma(phi)/N_sigma ,
```

axisymmetric, no swirl, odd in `z`. The tracked material point is `x_* = (rho_*, phi_0)` with
`rho_* := (1+f) rho_0`; `T_d(u_0) :=` the first time `||omega(t)||_inf >= (3/2) M`.

`W_sigma` is real-analytic on `[0, pi]` and even at both poles, so `eta_0` is `C^infinity` in the
angle including the axis. At the origin `Theta(rho) = rho^8/(rho^8 + e^2)` **exactly**, so every
Cartesian component of `omega_0` is `c |x|^7 x_j W_sigma(phi)` plus a term homogeneous of degree
`16`: `omega_0` is `C^infinity` on `R^3 \ {0}`, `C^{7,1}` at the origin, lies in
`L^2 cap L^inf` with `||omega_0||_inf = M sup Theta`, decays like `rho^{-8}` at infinity, and
belongs to `H^s(R^3)` **for `s < 9.5` and for no larger `s`**; hence `u_0 in H^s` for `s < 10.5`.

### 1.2 The Theorem V.4' ball

```
    d  :=  rho_* min( sin(phi_0 - delta - 4 sigma) ,
                      sin(pi/2 - delta_m - phi_0 - 4 sigma) ,  f/rho_* ) .
```

The first branch binds at every `f >= 4`, so `d = rho_* sin(phi_0 - delta - 4 sigma)` and the
smallest polar angle attained on `B(x_*, d)` is exactly `delta + 4 sigma`, four
mollification widths clear of the taper cone. At `f = 4`, `d = 1.5381397413 rho_0`,
`R_- = r_* - d = 0.9618602587 rho_0`, and the ball's inner radius is `3.4618602587 rho_0`.

### 1.3 The constants of the datum

```
    kappa_delta  = 0.4917868801        (kinked: 0.4978224383)
    c_*          = log(3/2)/kappa_delta = 0.8244732108
    c_2          = 2 log(3/2)/kappa_delta = 1.6489464217
    lam_max      = e^{3 c/4}           = 1.8558724485 at c = c_*, 2.0548738638 at the cap
    r_h          = inf_{[1, lam_max]} P_h(lam)/lam  >=  0.9033433307 at c = c_*
    E_0          = sup rho|eta_0|/M    in [7.5523054461, 7.5523079773]
    Gfrak_0      = sup rho^2|grad eta_0|/M in [36.0783630557, 36.0784195640]
    Psi_sigma(0) = sup rho|Lap_5 eta_0|/M  <=  495.1064102419
    C_E          = ||u_0||_2^2/(M^2 R^5) = 0.0551906693
    vartheta     = 0.2623138022  (SAFETY = 2, the largest over f >= 4; 0.2611297728 at f = 4)
```

### 1.4 The theorem

> **THEOREM (S3), version 3.** Let `u_0` be the datum of sec.1.1, with `nu > 0` and `f >= 4`, and
> let `u` be its maximal strong solution, which exists and is unique by **Proposition K'** of
> `fix5/FIX5.md` sec.3.3 (finite Sobolev index; the datum satisfies its hypothesis and does not
> satisfy Proposition K's). Assume, on `[0, c/(ML)]` with `c` the self-consistent window fixed
> below:
>
> * **(H3'_vartheta)** for the datum of sec.1.1 on the ball `B(x_*, d)` of sec.1.2, which holds
>   with `vartheta = 0.2623138022 < 1` measured (a five-parameter grid lower bound on the
>   supremum, inflated by `pmax-h3v` sec.B5's `SAFETY = 2`); hence
>   `Q = 1.7111799108` and `vt = 0.3555899554`.
> * `L >= L_Gamma^exist`, so that `(Gamma-off)`'s own fixed point exists:
>   `L_Gamma^exist = 9298.1782` at `c = c_*` and `16231.8518` at the certified window cap.
> * **(lam-cap)** nothing. `lam <= e^{3c/4}` is proved from the contradiction hypothesis alone.
>
> **No further hypothesis.** In particular `(H4)`, the bootstrap posture of `u2`, is not assumed:
> it is closed by the continuity and first-crossing argument of `FIX5.md` sec.4.
>
> Then, for every `L >= L_*`, there is `c` in `[c_*, 1.1647266909 c_*]` with
> `c >= c_*(1 + eps(c, L))`, and with `eps(L) := c/c_* - 1` for the smallest such `c`,
> ```
>     T_d(u_0)  <=  t_*  :=  c/(M L)  =  c_* (1 + eps(L))/(M L)
>                        =   2 log(3/2) (1 + eps(L))/(2 kappa_delta M L) ,
> ```
> with `c_* = 0.8244732108` and `c_2 = 1.6489464217`, and `eps(L)` decreasing to the clock floor
> `(1 - 2 kappa_delta)/(2 kappa_delta) = 1.6700567265e-02`. The reach is
> ```
>     L_*  =  1424207.9721  (all proved) ,  1379464.3531  (proved, C_a measured) ,
>             1761.8423     (all measured)                    at  eps <= 0.1647266909 ,
>     L_*  =  1554565.5693  (all proved) ,  1499581.4464  (proved, C_a measured) ,
>             2050.2770     (all measured)                    at  eps <= 0.1 .
> ```
>
> **In energy Reynolds number.** With `Re_E = E_0^{2/5} M_0^{1/5}/nu`, `E_0 = ||u_0||_2^2`, and
> `log Re_E = 2L + 2 log s + (2/5) log C_E = 2L + 2.9135781820`,
> ```
>     T(Lambda)  <=  c_2 (1 + eps(Lambda))/log Lambda ,   c_2 = 1.6489464217 ,
>     for  log Lambda >= log Lambda_* = 2 L_* + 2.9135781820 ,
>      log Lambda_*  =  2848418.8578 / 2758931.6198 / 3526.5981   at eps <= 0.1647266909 ,
>                       3109134.0522 / 2999165.8064 / 4103.4676   at eps <= 0.1 .
> ```

**The window cap.** A2 and C6 need `P_h` increasing on `[1, lam_max]`. With the **proved**
Lipschitz constant `|h_sigma''| <= 2415.5503900804` in place of `fix4`'s grid scan
`225.9466903722`, the largest window for which the cell-wise certificate closes is
`c/c_* <= 1.1647266909`, i.e. `eps <= 0.1647266909`, against `fix4`'s `0.1664585076`;
`min P_h'` at `fix4`'s cap is `-0.0029985709` with the proved pad. `P_h'` still vanishes at
`lam_mono = 2.0693384635`, so the certified cap sits `0.70 %` below the real obstruction.

---

## 2. THE DEPENDENCY GRAPH

`PROVED` means proved in the named file with its hypotheses met here. `ENCLOSED-certified` means
the statement is proved and the constant carries an outward bound with an explicitly proved pad.
`ENCLOSED-sampled` means the constant is a sampled supremum with a safety factor, not an interval
computation. `SLACK` means a proved but improvable constant, not a gap.

| # | statement | file | status |
|---|---|---|---|
| A1 | the datum; `N_sigma`; `omega_0` in `L^2 cap L^inf`; the decay; `omega_0 in H^s` exactly for `s < 9.5` | `fix4/p1`, `fix4/p2`, `fix5/x3`, `fix5/x6` | **ENCLOSED-certified** |
| A2 | `kappa_delta`, `c_*`, `c_2`, `lam_max`; `P_1(lam) = lam`; `mu_lam` a probability measure | `fix2/f4` (structure) + `fix4/p1` | **PROVED** structure, **ENCLOSED-certified** value |
| A3 | `r_h >= 0.9033433307` and `P_h' >= 0.3294825889 > 0` on `[1, lam_max(c_*)]`; the cap `0.1647266909` | `fix4/p1` certificate, proved Lipschitz input from `fix5/x6` | **ENCLOSED-certified** |
| A4 | `E_0`, `Gfrak_0` by `fix2` sec.4(c)'s branch bound | `fix2/f4` (structure) + `fix5/x6` | **PROVED** structure, **ENCLOSED-certified** value |
| A5 | the energy identity, `G_l`, `D_l`, `C_E`, the `L` to `log Re_E` dictionary | `THEOREM_S3/t1` + `fix4/p2` | **ENCLOSED-sampled** (a convergence statement, see sec.3 item 2) |
| B1 | Lemma 2: `eta <= 0` on `{z > 0}` for all `t`, all `nu >= 0` | `lower/prove-lagrangian` sec.3 | **PROVED** |
| B2 | far and near kernel lemma; Prop. 1, Prop. 2, Consequence A | `rebuild/far-near-kernel-lemma` | **PROVED**; hypotheses re-verified for this datum in `fix4/p3` |
| B3 | **Lemma T''** and Corollary 2'' on `[1, lam_max]`; the sensitivity `eps_T'` | `write/lemma-T-shell-dependent` + `fix5/x5` | **PROVED** |
| B4 | Theorem V.4' with `(H3'_vartheta)` for `k = 0..4` on the ball of sec.1.2; `vartheta = 0.2623138022` | `round2/pmax-h3v` Part B + `fix5/x1` | **PROVED** at `f >= 4`, `vartheta` **ENCLOSED-sampled** |
| B5 | Theorem V.1, `B1` | `gaps/gap-V-aronson` | **PROVED** |
| B6 | `(H-K2)`: `K_2 <= 161.7735 M/rho_0`, no `log` | `s3close/hk2` sec.0 with `(D1')` | **PROVED**; re-checked in `fix4/p3`, slack `3.9889277977`; not load-bearing |
| C1 | the slaved reference `frak_a`; (P1), (P2), (P3) | `round2/u1` sec.1 | **PROVED** |
| C2 | `Lambda` is the flow of `frak_a . D X`; `J_Lambda` | `u1` sec.1.2 | **PROVED** |
| C3 | the `l = 1` rate error is identically zero | `u1` sec.2.1 | **PROVED** |
| C4 | `frak_a` and Lemma T's `a_ref` are the same functional | `u1` sec.4.1 | **PROVED** |
| C5 | the restriction step (positivity plus `z`-oddness) | `u1` sec.4.2 step 2 | **PROVED** |
| C6 | the integro-differential inequality and the quasimonotone comparison | `u1` sec.5.1 + A3 | **PROVED** given A3 |
| C7 | the conservative clock `c_2 = 2 log(3/2)/kappa_delta` | `u1` sec.5.3 | **PROVED** |
| D1 | `\|\|grad u\|\|_op <= 2\|a\| + r\|grad a\| + \|omega^theta\|` | `u2` sec.2 | **PROVED** |
| D2 | `b(0) = 0`; `x . b <= Gamma_rad \|x\|^2` | `u2` sec.4 | **PROVED** |
| D3 | **P1**, **P2** | `round2/pmax-h3v` Part A | **PROVED** (interior parabolic Schauder cited, not proved) |
| D4 | `r\|grad a\| <= (3 pi^2/16) e^{p c_G} Gfrak_0 M sin phi`; Riesz `pi^4/2` | `u2` sec.6 | **PROVED** |
| D5 | `C_1 = 0` | `u2` sec.7(a) | **PROVED** |
| D6 | the axis-safe collar | `u2` sec.7(b) | **PROVED** |
| D7 | `(Gamma-off)`, and **(H4) closed by continuity and first crossing** | `u2` sec.8 + `fix2/f5` + `fix4/p4` + `fix5/x4` | **PROVED** for `L >= L_Gamma^exist`; `C''` **ENCLOSED-sampled** |
| E1 | `\|R_y\| <= chat_a(phi) M rho sin phi` | `refute-u1` F3(a) + B2 + D6 | **PROVED** |
| E2 | the `z`-integral identity | `refute-u1` F3(b) | **PROVED** |
| E3 | `\|R_z\|`, with the `r\|grad a\|` integral | `THEOREM_S3/t2` + D4 | **PROVED** |
| E4 | `C_R = sup_phi`, finite over every angle | `THEOREM_S3/t2` + `fix4/p4` | **PROVED**, value **ENCLOSED-sampled** |
| F1 | the majorant ODE, `mu = lam_max^2 m` | `u1` sec.2.3 | **PROVED** given D7, E4 |
| F2 | the Jacobian equation, `mu_J`, Lemma T'' (H3) | `u1` sec.2.3 | **PROVED** |
| F3 | Lemma T'' applied once, with `r_h` on `[1, lam_max]`; `lambda(.,s) in C^1` | `u1` sec.4.2 step 4 + `fix5/x5` | **PROVED** given A3 |
| F4 | `eps_v = eps_bulk[1 + vt(9 + 2 sigma_z/sigma_y)] + Q(E_hess + E_4 + E_tail)` | `u1` sec.4.2 step 5 + `fix5/x2` | **PROVED** given B4 |
| F5 | **(V-strain)**: the viscous defect of the strain functional, `eps_a` | `fix3` sec.1 + `fix4/p5` + `fix5/x6` | **PROVED**; `Psi_sigma(0)` **ENCLOSED-certified** |
| G1 | the self-consistent window; `eps(L)`, `L_*`, `log Lambda_*` | `fix2/f1` (instrument) + `fix5/x2` | **COMPUTED**, conditional on the ENCLOSED rows |
| G2 | the family satisfies BFG's hypotheses | sec.4 below | **PROVED** |
| G3 | BFG Theorem 10 is false as stated | sec.4 below | **CONDITIONAL** on G1 and on the ENCLOSED rows |
| M4 | `(lam-cap)`: the a priori cap `e^{3c/4}` | `u1` (P3) | **SLACK** |
| M6 | existence, uniqueness and (H\*) for this datum | Proposition K', `fix5/x3` | **PROVED** |

---

## 3. THE RESIDUAL LIST

Nothing in the chain is unproved except item 1. What else remains is a list of constants that are
sampled suprema with safety factors rather than interval computations.

1. **Interior parabolic Schauder** is cited, not proved, in `round2/pmax-h3v` Part A. That is the
   one textbook citation the P1/P2 chain rests on, and the only genuinely unproved analytic item
   in the graph.
2. **`C_E`** is a convergence statement (`7.3e-06` in `lmax`, `4.7e-07` in `du`), not an
   enclosure. It enters only as `(2/5) log C_E`, so `7.3e-06` relative moves `log Lambda_*` by
   `3e-06` out of `2.8e+06`.
3. **`vartheta`** is a grid-plus-polish lower bound on a supremum over a five-dimensional set,
   inflated by `SAFETY = 2` per `pmax-h3v` sec.B5. Margin: `0.2623138022` against `1`, a factor
   `3.81`; coarse-to-fine change `4.6e-03`. It is not an interval computation. And the geometric
   margin is thin: at an offset of `3.5 sigma` instead of `4 sigma`, `2 vartheta = 1.76645` and the
   hypothesis fails. The `6 sigma` ball gives `2 vartheta = 0.1768303` at no cost in `L_*` and is
   the safer statement.
4. **`C''`, `Ghat`, `C_R`, `chat_a`, `L_Gamma^exist`** are quadratures with a sign-change
   existence proof; `C_R`'s supremum over the angle is a sampled search, checked at six
   resolutions to `3e-08`.
5. **`K_2`** is imported from `hk2` for its own datum and re-checked pointwise at four `lambda`.
   Sweeping `C_K` past it by four orders of magnitude moves `L_*` in the fourth digit only.
6. **The far and near constants** are Parseval-based numerical evaluations in the source seat,
   not intervals; their three hypotheses are profile-free and were re-verified for this datum.
7. **`eps(L)` monotone decreasing to the floor** is tabulated at five values of `L`, not proved.
8. **The CKN attribution is NOT used by this theorem.** Proposition E and Corollary C2 of
   `02_ADDENDUM_uniqueness_2026-09-08.md` are needed only to rule out a general Clay-class
   competitor. This theorem is a statement about the maximal strong solution supplied by
   Proposition K', which is proved by in-estate citations; if a referee insists on the Clay-class
   formulation, then Proposition E enters and that attribution becomes an ENCLOSED step.
9. **A discretisation error in `fix4`**, found here and corrected: `fix4/p1_profile` computes
   `W_sigma'` by an FFT convolution of the discontinuous `Wtil'`, which is first order through
   the jump, so `Gfrak_0 = 36.0795702396` is `3.2e-05` too large. The direction is conservative
   and the corrected value is used above.

---

## 4. THE BFG CONSEQUENCE, FOR THE RECORD

Bradshaw, Farhat and Grujic, arXiv:1704.05546v4, Theorem 10, p.11. The clause, stated rather than
quoted: for a datum `omega_0` in `L^2 cap L^inf` and any constant `M > 1`, there is `c(M) > 1`
such that a unique mild solution exists in `C_w([0,T], L^inf)` on a time interval of length at
least `(1/c(M))(1/||omega_0||_inf)`, the solution is on `(0,T]` the restriction to `R^3` of a
holomorphic function, and its `L^inf` norm on the corresponding domain stays below
`M ||omega_0||_inf`. Taking `M = 3/2` gives

```
     M_0 T_d(u_0)  >=  1/c(3/2)  >  0     for EVERY datum with omega_0 in L^2 cap L^inf,
```

with `c(3/2)` absolute: no dependence on `||omega_0||_2` (BFG's own sentence on p.8), none on
`nu` (both sides are invariant under `u -> lam u(lam x, lam^2 t)` at fixed `nu`), none on
geometry.

**The family satisfies each of BFG's hypotheses, re-checked for the mollified profile** by
`fix4/p2`, `fix4/p3` and independently by the round-2 referee's `s5_bfg_datum.py`:

| clause | requirement | this datum |
|---|---|---|
| `omega_0 in L^inf` | finite | `sup\|A_sigma\|/N_sigma = 1` exactly by construction; `\|\|omega_0\|\|_inf = M sup Theta` |
| `omega_0 in L^2` | finite | `\|\|omega_0\|\|_2^2/(M^2 R^3) = 1.4381082627`, finite at every `L` |
| decay | "suitable" | `O(rho^{-8})` at infinity and `O(rho^8)` at the origin, both against the closed forms |
| the class | `C_w([0,T], L^inf)` | the maximal strong solution of **Proposition K'** lies in it, by K'(e') |
| `M > 1` | `M = 3/2` | yes |

The fourth row is the one that changed. In `v2` it appealed to Proposition K, whose hypothesis
this datum fails; K'(e') gives the same conclusion from a hypothesis the datum satisfies, so the
two statements are now about the same object for a reason that has been checked.

The theorem of sec.1.4 gives `M T_d <= c_2 (1 + eps)/(2L) -> 0` as `L -> infinity` at fixed `M`.
The two cannot both hold. **Conditional on the residual list of sec.3 and on nothing else, BFG
Theorem 10 is false as stated**, and the located defect in its displayed proof (the p.10
weighted-BMO step without the local average, which two campaign referees called false and which
nobody has re-verified) is where to look.

*What is different from `THEOREM_S3_v2.md` sec.4.* There, the distance was "the list of
enclosures in sec.3, the two citations in it, and a referee's reading of the chain", and the
referee's reading found a false hypothesis. The distance is now one textbook citation, a list of
sampled constants with safety factors, and a second referee's reading. It is not a public claim:
the estate's standing rule is that the refutation is not asserted until the chain has been
refereed by a seat that did not build it, and this version has been refereed once.

---

## 5. WHAT A REFEREE ATTACKS FIRST

1. **The `4 sigma` offset, and whether it is enough.** `vartheta` is a grid lower bound on a
   supremum over a five-dimensional set, inflated by a factor two, and the geometry it sits in has
   a cliff: at `3.5 sigma` the hypothesis fails at `SAFETY = 2` and at `4 sigma` it holds with a
   factor `3.81` of room. A referee should ask (a) whether the five-parameter search is dense
   enough near `phi = delta + 4 sigma`, where the binding order `k = 4` lives and where the
   fourth angular derivative is still `5.0e+04`; (b) whether `SAFETY = 2` is the right convention
   for a supremum over a set of that dimension; and (c) why the statement is not simply made at
   `6 sigma`, where `2 vartheta = 0.1768303` and `L_*` is identical. The honest answer to (c) is
   that the brief asked for `4 sigma`. This is where I would spend the first day.

2. **The normalisation, and whether `kappa_delta` is still the right strain constant.** Unchanged
   from `v2` sec.5 item 1, and unaffected by anything here. The whole cost of the smooth datum is
   that `kappa_delta` falls from `0.4978224383` to `0.4917868801`, which raises the clock floor by
   a factor `3.8179938670`, and the whole gain is that `Gfrak_0` falls by `1.6268342260`. A
   referee should check that mollifying `W = A/sin phi` rather than `A` is the right operation and
   that dividing by `N_sigma` is the right response to the convexity overshoot.

3. **`C_R`, and whether it is really the supremum**, given that `mu` and `mu_J` are suprema over
   the shell while the pieces of `R` are estimated at a field point. Unchanged from
   `THEOREM_S3.md` sec.7 item 2. The answer changes `C_R` linearly, hence `L_*` linearly.

4. **The window cap, which is now tighter and is the only place unit 6 cost anything.** Beyond
   `eps = 0.1647266909` the certificate that `P_h` is increasing does not close, and C6's
   quasimonotone comparison then has no hypothesis at all. The cap moved because the Lipschitz
   constant it rests on was a grid scan and is now proved, a factor `10.69` larger. A referee
   should check that the cell-wise certificate is genuinely a bound, and that the `9 lam^9 sum_i
   max(M_i,0) I_i` majorant really is increasing in `lam` (`FIX4.md` sec.2.1's stated reason for
   testing only at `lam_max` is a non-sequitur, as the referee's F-9 records; the conclusion was
   checked directly at 441 values of `lambda`).

5. **`theta_0` in the bootstrap.** The strictness that closes `(H4)` comes from
   `lamhat := sup_{[0,tau]} ||omega||_inf/M < 3/2`, which is the contradiction hypothesis on a
   compact interval, not from any constant in the sheet. That is legitimate and it is also the
   whole of the argument's contact with the contradiction posture; a referee should satisfy
   themselves that `T_d > tau` does entail `tau < T_*`, which it does only because `T_d <= T_*`
   by definition.

---

## 6. ABSTRACT

For axisymmetric no-swirl Navier-Stokes we bound the vorticity doubling time of a family of shell
data. The family is a smooth, equatorially tapered and radially mollified vortex ring of amplitude
`M`, shell aspect ratio `L = log(R/rho_0)` and inner scale `rho_0 = s sqrt(nu/M)`; its angular
profile is the Gaussian mollification, at scale `sigma = 0.02` rad, of a piecewise linear taper,
renormalised so that the vorticity maximum is exactly `M`. For `L` above a threshold, the first
time the vorticity maximum reaches `(3/2)M` satisfies `M T_d <= c_2 (1 + eps(L))/(2L)` with
`c_2 = 2 log(3/2)/kappa_delta = 1.6489464217`, and in energy Reynolds number
`T(Lambda) <= c_2 (1 + eps)/log Lambda` for `log Lambda` above `2848418.8578`. The error
`eps(L)` decreases to `1.6700567265e-02`, the strain deficit of the mollified profile. Every
analytic hypothesis the earlier assembly carried is discharged, the Theorem V.4' ball now stands
four mollification widths off the taper cone so that its own hypothesis holds, the existence and
uniqueness clause rests on a proposition whose hypothesis this datum satisfies, the bootstrap
posture is closed rather than assumed, and the numerical enclosures carry proved pads. What
remains is one textbook citation and a list of sampled constants with safety factors. Granting
them, the bound contradicts Bradshaw, Farhat and Grujic's Theorem 10, whose lower bound on the
same quantity is absolute.

**FL-000 stands. Nothing here touches the headline problem.**

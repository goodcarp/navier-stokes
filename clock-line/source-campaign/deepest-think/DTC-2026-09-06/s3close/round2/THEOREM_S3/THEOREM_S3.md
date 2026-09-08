# THEOREM (S3), assembled: the sharp logarithmic doubling clock

Seat `s3close/round2/THEOREM_S3`, sitting of 2026-09-08.  Re-assembler.
Laws: `TORMENT NEXUS/LAWS.md` (first 120 lines read).

Every number below came out of a script in **this** folder that I wrote and ran
(`t1_datum.py`, `t2_gamma_CR.py`, `t3_budget.py`, `t4_derived.py`), re-asserted by
`check_constants.py`; `SHA256SUMS` is computed with `shasum`, never typed.  Nothing outside
this folder was written.

**Read, in the order the brief set** (L-14): `s3close/assembly/ASSEMBLY.md` in full;
`s3close/hk2/PROOF.md` §0 and §1; `s3close/round2/u1/PROOF.md` in full;
`s3close/round2/refute-u1/NOTE.md` in full; `s3close/round2/u2/PROOF.md` in full.
**Code imported:** `s3close/round2/u1/u5_budget.py` (itself a copy of
`s3close/assembly/a2_budget.py`) is the base of `t3_budget.py`; its viscous half
(`viscous_budget`, `gauss_tail_aniso`, `J_pow`, `eps_Tprime`, the `RECORD` block) is carried
over **unchanged**, so the viscous column is literally the assembly's own instrument and the
comparison is apples to apples.  Byte copies of every imported file, with their hashes, are in
`copies/`.  This is L-14's inverse convention, declared: the object under test is imported;
everything I am testing a *claim* about (the datum constants, the energy, U2's fixed point,
`C_R`) was re-implemented here from the mathematics.

---

## 0. HEADLINE

**(S3) is assembled, modulo a named list of eight items, and the reach is `L_* = 311511.9`.**

```
   L_*(eps <= 1/2)        log Lambda_*        route
   4.5253e14              9.0506e14           ASSEMBLY  (a priori reference, e^{34} feedback)
   16161.986              32327.336           u1        (slaved reference, at phi_0, G_collar = 0)
   52437.352              104878.07           refute-u1 (the same, honestly on the cone)
   311511.9               623026.8            THIS NOTE (all proved, (R-z) and (A-cone) DISCHARGED)
      422.4                  847.8            THIS NOTE (all measured)
```

The number goes **up** by a factor `5.94` against `refute-u1`'s corrected figure, and it goes up
for a reason worth stating plainly: `u1` reached `16161.986` while carrying two unproved
hypotheses, **(R-z)** (the collar gradient constant `G_collar`, set to `0` in every headline
number) and **(A-cone)** (the majorant established on a cone that excludes the axis).  This note
**discharges both**, using `u2`'s proved global gradient bound
`r|grad a| <= (3 pi^2/16) e^{p c_G} Gfrak_0 M sin phi` and `u2`'s axis-safe collar
`|a_collar| <= 2 R_A e^{c_G} E_0 M`.  Setting an unproved constant to zero is cheaper than
proving it.  The exchange is `L_*` from `52437.352` to `311511.9` and a theorem with two fewer
hypotheses.

**Four findings.**

1. **The binding constant has moved.**  `u1`'s conclusion was that what binds is the far/near
   constant `C_a`, whose measured value is `0.069`: replacing it alone took `L_*` from `16161.986`
   to `2389.4697` there.  Here, replacing `C_a` by `0.069` and changing nothing else moves `L_*`
   by **`2.8 %`** (`311511.9` to `302704.3`).  What binds now is `Ghat`, the proved pointwise
   bound on `r|grad a|`: it supplies **`75.6 %`** of `C_R`.  Through `Ghat` the binding *datum*
   quantity is `Gfrak_0 = sup |x|^2 |grad eta_0| / M`.
2. **`(Gamma-off)` has an existence threshold of its own, and nobody has stated it.**  `C''` is
   not a constant; it is the solution of a fixed point that runs away below
   `L_Gamma_* = 14260.5` for this datum (`u2`'s own `L_* = 5313.8` is for `u2`'s datum).  Below
   that `L` the theorem is empty for a reason that has nothing to do with the clock.
3. **The whole viscous column is now numerically irrelevant and still hypothetically load-bearing.**
   At `L = 1e+06`, `eps_v = 1.7343e-08` and `E_tail` is exactly `0.0`.  `(H-K2)` is worth
   `1.4e-06` of `L_*`.  But `(H3-V)` is still a *hypothesis* the theorem's own datum does not
   satisfy, and a small number is not a discharged hypothesis.
4. **The mandatory radial mollification costs `0.44086` in `log Lambda_*`.**  `hk2` forces a
   `tanh` radial ramp (a sharp radial edge makes `K_2` infinite).  `ASSEMBLY`'s energy constant
   `3.364345456912714` is for sharp edges.  For the theorem's actual datum
   `C_E = 0.05657477`, not `0.17032563295410327`, and the shift is **`2.9234859`**.

**FL-000 stands.**  (S3) is not proved: the modulo list of §3 is not empty.  Nothing here
touches the headline problem.

### 0.1 Grade

**FILED.**  No Solid, no Major.  The content is (a) the theorem stated with every hypothesis
named, (b) the discharge of `(R-z)` and `(A-cone)`, (c) the budget recomputed with the surviving
constants and a datum-dependence that was not previously tracked, (d) four findings.

---

## 1. THE THEOREM

### 1.1 The datum, stated fully

Fix `M > 0`, `delta = 7.5 deg`, `delta_m = 5.0 deg`, `phi_0 = 30.0 deg`,
`s := 1/sin delta = 7.66129757554039`, `eps_r = 0.25`, `f > 0`, `L > 0`.  Write
`u := log(rho/rho_0)`, `phi_ax := min(phi, pi - phi)`.  Put

```
    nu   := M rho_0^2 / s^2        (equivalently rho_0 = s sqrt(nu/M)) ,      R := rho_0 e^L ,

    Theta(rho) := (1/2)[ tanh( (u - eps_r)/eps_r ) - tanh( (u - L + eps_r)/eps_r ) ] ,

    omega_0^theta(rho,phi) := - M * Theta(rho) * sgn(cos phi)
                                  * min(1, phi_ax/delta) * min(1, |phi - pi/2|/delta_m) ,

    omega_0 = omega_0^theta e_theta ,   u_0 = Biot-Savart of omega_0 ,
```
axisymmetric with no swirl, odd in `z`.  `eta_0 := omega_0^theta / r`.

This is `ASSEMBLY` §1.1's **angular** profile with `hk2`'s **radial** mollification.  The radial
ramp is not a choice: `hk2` §0 proves `K_2 = +infinity` for a sharp radial edge, and its
sharpening of the V-b refuter's Finding 4a says the coupling lemma needs the *global* `K_2`, so
an inset does not rescue it.  Written in the log variable, `rho Theta'(rho_0) = 1/(2 eps_r) = 2`,
which is exactly `hk2`'s `rho_0/(2 w_0)` at `w_0 = 0.25 rho_0`: the inner slope agrees with
`hk2` (D-B) to the digit.

The **tracked material point** is `x_* = (rho_*, phi_0)` with `rho_* := (1+f) rho_0`;
`L := log(R/rho_0)`; `T_d(u_0) :=` the first time `||omega(t)||_inf >= (3/2) M`.

Two properties of this datum, computed (`t1`):

* `||omega_0||_inf = M sup_rho Theta = M (1 - O(e^{-2L/eps_r}))`; the nominal scale is `M`.
* `E_0 := sup |x| |eta_0| / M = 7.6612976`, attained at the taper edge `phi = delta`, where the
  profile's `min(1, phi_ax/delta)/sin phi` peaks at `1/sin delta = 7.6612976`.
  `Gfrak_0 := sup |x|^2 |grad eta_0| / M = 58.695476`, attained at `phi = pi - delta`, the
  taper kink, where `d/dphi [ 1/sin phi ]` jumps to `-58.193328`.
  (For the fully-`tanh` campaign datum (D-C) the same instrument returns `7.6606020` and
  `20.919714`, reproducing `u2`'s `7.66060196` and `20.9197070`.  The theorem's datum has
  `Gfrak_0` larger by `2.8057`, and that factor goes straight into `Ghat`, `C_R` and `L_*`.)

### 1.2 The strain constants

`kappa_delta := (1/2) P_h(1)`, `P_h(lam) = 3 int_0^1 h(arcsin v) v^2 (A v^2 + B)^{-5/2} dv`,
`A = lam^2 - lam^{-4}`, `B = lam^{-4}`, with `h` the **full** angular profile
(taper *and* equatorial mollifier), which is what the datum has and what `refute-u1` F4 and
`u2` §9.1 both say `ASSEMBLY` §1.2 got wrong by using the taper alone:

| profile | `kappa` | `1 - 2 kappa` | clock floor `(1-2kappa)/(2kappa)` |
|---|---|---|---|
| bang-bang | `0.5` | `0` | `0` |
| axis taper `7.5 deg` only (`ASSEMBLY` §1.2) | `0.49972123` | `5.5754e-04` | `5.5785e-04` |
| **the theorem's datum (`7.5 deg` and `5.0 deg`)** | **`0.49782244`** | **`4.3551236e-03`** | **`4.3741737e-03`** |
| campaign (D-C), `tanh` taper and `tanh` equator | `0.47473436` | `5.0531e-02` | `5.3220594e-02` |

`0.49782244` reproduces `u2` §9.1's `0.4978224382` to sixteen digits and `refute-u1` F4's
`0.4978226672` to `4.6e-07`; the taper-only entry reproduces `ASSEMBLY`'s
`0.4997212305210886`.  Consequently

```
   c_*  := log(3/2)/kappa_delta            = 0.8144774        (u1's taper-only value: 0.8113826)
   c_2  := 2 log(3/2)/kappa_delta = 2 c_*  = 1.6289547        (refute-u1 F4: 1.6289540)
   lam_max := e^{3 c_*/4}                  = 1.8420112
   r_h  := inf_{[1,lam_max]} P_h(lam)/lam  = 0.91863653       (attained at lam_max)
```
`P_h` is increasing on `[1, lam_max]` (minimum increment `0.0040871` on a `0.01` grid; the first
decrease is at `lam = 2.09 > lam_max`), which is what makes the shell system quasimonotone.
`P_1(lam) = lam` exactly, at `lam = 1.5` and at `lam_max`, as a control.

### 1.3 The statement

> **THEOREM (S3), assembled form.**  Let `u_0` be the datum of §1.1 and let `u` be its `H^s`
> strong solution.  Assume, on `[0, tau]` with `tau := c_*/(M L)`:
>
> * **(P-max)** the approximate-maximum step of `u2`'s two weighted maximum principles P1 and
>   P2 (§3, item M1);
> * **(H3-V)** `ASSEMBLY` BLOCK 5's hypothesis for Theorem V.4 (§3, item M2);
> * **(V-strain)** the viscous defect of the strain *functional* (§3, item M3);
> * **(lam-cap)** nothing (this is the a priori cap `lam <= e^{3 c_*/4}`, which is proved; it is
>   listed in §3 only because a self-consistent argument would improve it);
> * **(H5-Lip)** `u2`'s (H5) for a datum that is Lipschitz, not `C^1`, in the angle (§3, M5);
> * **(H\*)** the solution is the `H^s` strong solution, unique in `C_w([0,T], L^inf)` (§3, M6);
> * **(Theta-tail)** the shell support statement of (D1) holds up to the `tanh` tails (§3, M7);
> * `L >= L_Gamma_* = 14260.5`, so that `(Gamma-off)`'s own fixed point exists (§3, M8).
>
> Then there are `L_* < infinity` and an explicit `eps(L)` decreasing to
> `(1 - 2 kappa_delta)/(2 kappa_delta) = 4.3741737e-03` such that for every `L >= L_*`
> ```
>     T_d(u_0)  <=  t_*  :=  c_* (1 + eps(L)) / (M L)
>                        =   2 log(3/2) (1 + eps(L)) / (2 kappa_delta M L) ,
> ```
> with `c_* = 0.8144774` and `2 log(3/2)/kappa_delta = c_2 = 1.6289547`, and
> ```
>     L_* = 311511.9   (all proved) ,   302704.3   (proved, C_a measured) ,   422.4  (all measured)
> ```
> for `eps <= 1/2`, and `1056242.9 / 1025043.0 / 1119.3` for `eps <= 0.1`.

The argument is `u1`'s, by contradiction: assume `T_d > tau`.  Then `||omega(t)||_inf < (3/2) M`
on `[0, tau)`, so the range hypothesis `M_s <= (3/2) M` that every record theorem carries is not
an assumption but the contradiction hypothesis, and the a priori cap
`lam <= e^{3 c_*/4} = 1.8420112` follows from it alone.  The clock closes when
`kappa_delta c (1 - eps_a) + log(1 - eps_v) >= log(3/2)`, i.e.

```
    1 + eps  =  [ 1 + log(1/(1-eps_v))/log(3/2) ] / [ (1 - eps_a)(1 - eps_delta) ] ,
    eps_delta := 1 - 2 kappa_delta = 4.3551236e-03 ,
    eps_a     := (ell_loss + ell_ramp)/L + eps_T' + chat_a(phi_0)/(kappa_delta L) .
```

### 1.4 The conversion to `T(Lambda)`

`T(Lambda) := inf{ M_0 T_d(u_0) : Re_E(u_0) <= Lambda }`, `Re_E = E_0^{2/5} M_0^{1/5}/nu`,
`E_0 = ||u_0||_2^2` (`ASSEMBLY` §1.4's convention; the brief's `E_0 = (1/2)\int|u_0|^2` shifts
`log Re_E` by the additive constant `(2/5) log(1/2)` and nothing else).  Exactly,

```
    log Re_E = 2L + 2 log s + (2/5) log C_E ,      C_E := ||u_0||_2^2/(M^2 R^5) ,
```
and `kappa_delta` does **not** enter this shift.  `C_E` is derived here, not quoted: with
`||u_0||_2^2 = (1/pi) \int psi_1 eta dx_5`, the radial Green function
`G_l(rho,rho') = rho_<^l rho_>^{-(l+3)}/(2l+3)` against `rho'^3 drho'`, and
`W(t) = sum_l H_l C_l^{3/2}(t)` with `N_l = (l+1)(l+2)/(l+3/2)`,

```
    C_E = 2 pi sum_{l odd} 2 N_l Hhat_l^2 D_l/(2l+3) ,
    D_l = \int\int_{x<y} Theta(x) Theta(y) e^{(l+4)x - (l-1)y} dx dy ,   x = log(rho/R) .
```

| datum | `C_E` | shift `2 log s + (2/5) log C_E` |
|---|---|---|
| bang-bang cap (control) | `0.17240397` (record `0.172403978`, rel `2.0e-08`) | |
| sharp radial edges, `delta = 7.5`, `delta_m = 5.0` | `0.17032563` (`ASSEMBLY`: `0.17032563295410327`) | `3.364345456912714` |
| **the theorem's datum (`tanh` ramp, `eps_r = 0.25`)** | **`0.05657477`** | **`2.9234859`** |

The `D_l` recursion agrees with the sharp closed form `D_l = 1/(5(l+4))` (in the `L -> infinity`
limit) to `5.0e-05`, the `du`-refinement moves `C_E` by `4.8e-07` and `l_max` by `9.5e-06`.

> **`T(Lambda) <= c_2 (1 + eps(Lambda))/log Lambda`, `c_2 = 1.6289547`, for `Lambda >= Lambda_*`,
> `log Lambda_* = 2 L_* + 2.9234859`.**

So `log Lambda_* = 623026.8` (all proved, `eps <= 1/2`), `605411.5` (proved, `C_a` measured),
`847.8` (all measured).  With `ASSEMBLY`'s sharp-edge shift the same `L_*` gives `623027.3`,
which is the size of the correction: `0.44086`, invisible against `L_*` and not against
anything else.

---

## 2. THE DEPENDENCY GRAPH, DATUM TO `T_d`

Each row: the statement, the file it comes from, and its status **in this assembly**.
`PROVED` means proved there and used here with its hypotheses met.

| # | statement | file | status |
|---|---|---|---|
| A1 | the datum, `E_0`, `Gfrak_0`, the trajectory margins | this folder, `t1` | **PROVED** |
| A2 | `kappa_delta`, `r_h[1,lam_max]`, `P_h` increasing on `[1,lam_max]`, `P_1(lam)=lam` | this folder, `t1` (`lower/prove-lagrangian` §2 for the identity) | **PROVED** |
| A3 | the energy identity, `G_l`, `D_l`, `C_E`, the `L <-> log Re_E` dictionary | this folder, `t1` (`ASSEMBLY` §1.4 derived the sharp case) | **PROVED** |
| B1 | **Lemma 2**: `eta <= 0` on `{z>0}` for all `t`, all `nu >= 0` | `lower/prove-lagrangian` §3 | **PROVED** |
| B2 | far/near kernel lemma, Prop. 1, Prop. 2, Consequence A (`z`-odd constants `0.291999`, `0.014754`, `3.999218`, `pi/8`) | `rebuild/far-near-kernel-lemma` | **PROVED** |
| B3 | **Lemma T'** and its Corollary 2; the sensitivity `eps_T'` | `write/lemma-T-shell-dependent` (refuter: STANDS) | **PROVED** |
| B4 | Theorem V.4 / Cor. V.5 with the refuter's sizing (`87.166`) | `write/V-b-...` + `write/refute-V-b-...` | **PROVED-modulo-(H3-V)** |
| B5 | Theorem V.1, `B1` | `gaps/gap-V-aronson` | **PROVED** |
| B6 | **(H-K2)**: `K_2 <= 161.7735 M/rho_0` over `lam` in `[1,3/2]`, no `log`, for the `tanh`-radially-mollified datum | `s3close/hk2` §0 | **PROVED** |
| C1 | the slaved reference `frak_a(rho,s) := \int_{2rho}^inf -(3/5) g_1 dlog rho'`; (P1) `frak_a >= 0`; (P2) `\|d frak_a/dlog rho\| <= M_s/2`; (P3) `lam <= e^{3c/4}`, decreasing | `s3close/round2/u1` §1 (refuter: re-derived, holds) | **PROVED** |
| C2 | `Lambda` is the flow of `frak_a . D X`; `J_Lambda = lam^2[1 + (rho lam'/lam)(1 - 3 cos^2 phi)]` | `u1` §1.2 | **PROVED** |
| C3 | the driver: the `l = 1` rate error is **identically zero** | `u1` §2.1 (refuter: re-derived) | **PROVED** |
| C4 | `frak_a` and Lemma T's `a_ref` are the same functional up to the inner-shell integral | `u1` §4.1 (refuter: exact, nothing hidden) | **PROVED** |
| C5 | the restriction step (positivity plus `z`-oddness) | `u1` §4.2 step 2 | **PROVED** |
| C6 | the integro-differential inequality and the quasimonotone comparison | `u1` §5.1 + A2's monotonicity | **PROVED** |
| C7 | the conservative clock `c_2 = 2 log(3/2)/kappa_delta` | `u1` §5.3 | **PROVED** |
| D1 | `||grad u||_op = ||grad_5 b||_op <= 2|a| + r|grad a| + |omega^theta|` | `u2` §2 | **PROVED** |
| D2 | `b(0)=0`; `x . b <= Gamma_rad |x|^2`; the `Gamma_rad` envelope | `u2` §4 | **PROVED** |
| D3 | **P1** `sup rho|eta| <= e^{c_R} E_0 M`, **P2** `sup rho^2|grad eta| <= e^{p c_G} Gfrak_0 M`, `p <= 3` | `u2` §5 | **PROVED-modulo-(P-max)** |
| D4 | `r|grad a| <= (3 pi^2/16) e^{p c_G} Gfrak_0 M sin phi`, Riesz composition `pi^4/2` in `R^5` | `u2` §6 | **PROVED-modulo-(P-max)** |
| D5 | `C_1 = 0` (`a_far(0,s) <= a(0,s)`, by B1) | `u2` §7(a) | **PROVED** |
| D6 | the axis-safe collar `\|a_collar\| <= 2 R_A e^{c_G} E_0 M`, no angular dependence | `u2` §7(b) | **PROVED-modulo-(P-max)** |
| D7 | **(Gamma-off)** `Gamma(s) <= 2 a(0,s) + C'' M` on all of `R^5`, `C''` the fixed point | `u2` §8 | **PROVED-modulo-(P-max)**, and only for `L >= L_Gamma_*` |
| E1 | `\|R_y\| <= chat_a(phi) M rho sin phi` (the collar's `1/sin phi` cancels against `\|Y\|`) | `refute-u1` F3(a) + B2 + D6 | **PROVED-modulo-(P-max)** |
| E2 | `\int_0^{\|z\|} log(rho/rho_zeta) dzeta = \|z\| - r arctan(\|z\|/r)`, so there is no `log(1/sin phi)` | `refute-u1` F3(b) | **PROVED** |
| E3 | `\|R_z\|`, with `\int_0^{\|z\|} r\|grad a\| dzeta <= Ghat M r log((\|z\|+rho)/r)` | this folder, `t2` (from D4) | **PROVED-modulo-(P-max)**.  **(R-z) is DISCHARGED**: `G_collar` is not used |
| E4 | `C_R = sup_phi` of the above, finite over **every** angle | this folder, `t2` | **PROVED-modulo-(P-max)**.  **(A-cone) is DISCHARGED** |
| F1 | the majorant `dm/dtheta = G m + lam_max(C_R + 2 lam_om log lam_max)/L`, `mu = lam_max^2 m` | `u1` §2.3, re-parametrised against `M` here | **PROVED** given D7, E4 |
| F2 | the Jacobian equation, `mu_J`, Lemma T' (H3) obtained without Corollary 3 | `u1` §2.3 | **PROVED** |
| F3 | Lemma T' applied **once**, to `eta_0 1_{S^out}`, with `r_h` on `[1, lam_max]` | `u1` §4.2 step 4 | **PROVED** (B3; the shell supremum is now finite, E4) |
| F4 | `eps_v = eps_bulk + E_hess + E_4 + E_tail` | `u1` §4.2 step 5, instrument from `a2_budget.py` | **PROVED-modulo-(H3-V)**, carried gap **(V-strain)** |
| G1 | the budget, `eps(L)`, `L_*`, `log Lambda_*` | this folder, `t3` | **COMPUTED**, conditional on the modulo list |
| G2 | BFG hypothesis check for this family | §5 below | **PROVED** |
| G3 | BFG Theorem 10 is false as stated | §5 below | **CONDITIONAL** on the modulo list |

**Off the critical path, stated so it is not looked for.**

* **BLOCK 2, the restart lemma, is NO LONGER NEEDED.**  `u1`'s slaved reference removes the
  `l = 1` rate error from the driver identically, so the feedback exponent is `G c = O(1)`
  rather than `34.0573`, and the one-window argument closes at a finite `L`.  `ASSEMBLY` §3.4
  priced a restart at twelve orders of magnitude; that price is no longer payable because the
  thing it was buying is already had.  BLOCK 2 remains open and is no longer binding.
* **Theorem A, Corollary A1, the L3v sum, `(B.2)`, `(B.3)`, `(B.4)`, `(B.5)`, `(B.6)`,
  homogeneity, the Gegenbauer machinery and the slab are all off the critical path** (`u2` §10).
  `u2` proves hypothesis (i) of `B(t)` with no modal expansion: the chain is
  `||grad u|| <= 2|a| + r|grad a| + |omega|`, then P1/P2, then Riesz, then far/near.  They remain
  the sharp route for the exact plateau (`C' <= 151.15` there against `C'' ~ 1461` here), and
  the L3v refuter's MAJOR-1 therefore no longer propagates into the bootstrap.
* **BLOCK 1** (no algebraic contraction) is not fatal and never was: the majorant is a
  differential comparison.  **BLOCK 4** (origin versus material shell) was solved by
  `ASSEMBLY` §2.4.  **BLOCK 6** is `(H-K2)`, closed by `hk2`.

---

## 3. THE MODULO LIST, COMPLETE

Eight items.  Nothing else in the chain above is unproved.

**M1.  (P-max): the approximate-maximum step in `u2`'s P1 and P2.**  *What it needs.*  P1 and
P2 run the parabolic maximum principle on `V := |x| |eta|` and `W := |x|^2 |grad eta|`.  Neither
is `C^2`: `|.|` is Lipschitz and not differentiable at its zero set, `|x|` is not smooth at the
origin, and the argument uses the upper Dini derivative of `sup V` (Hamilton's trick) rather
than a classical interior maximum.  The standard fix is to mollify `|.|` to `sqrt(. ^2 + eps^2)`,
work on `R^5 \ B_eps`, use `V(0) = W(0) = 0` (which holds because `eta` is bounded) and the
decay at infinity to exclude escape of the supremum, and let `eps -> 0`.  It is routine and it
is **not written out anywhere in the chain**.  *What it carries.*  Everything: P1 gives the
axis-safe collar and hence `Chat_a` and the discharge of (A-cone); P2 gives `Ghat` and hence
`C''`, `(Gamma-off)`, and the discharge of (R-z).  If M1 fails, D3 through D7 and E1, E3, E4 all
fail and nothing in the budget survives.  This is the single most load-bearing item on the list
and the cheapest to close.

**M2.  (H3-V): `ASSEMBLY` BLOCK 5.**  Theorem V.4 and Cor. V.5 require
`eta_0 = -M sgn(z)/r` **exactly** on `B(x_*, d)`.  The theorem's datum never equals that
anywhere: the `tanh` radial ramp is never identically `1`, and the two angular mollifiers
change `eta_0` off the plateau.  `ASSEMBLY` §2.6 BLOCK 5 gives the repair (restate V.4 with
`|Delta_5 eta_0 + eta_0/r^2| <= vartheta M/r^3` on the ball, a change to Step 2's Taylor
remainder, not written) and notes that `min(1,.)` was chosen over `tanh` precisely to make (H3)
literally true away from the kinks; `hk2` then forced `tanh` back into the radial direction.
The two requirements are in direct conflict and the conflict is unresolved.  *Cost if closed at
the obvious price:* nothing measurable, since `eps_v <= 1.2e-07` at every `L >= 1e+05` in the
proved column.  A small number is not a discharged hypothesis.

**M3.  (V-strain): the viscous defect of the strain functional.**  Theorem V.4 is pointwise, on
a tube around the tracked point; `a[.](0)` is a weighted integral of `eta` over the whole shell,
and the defect at the axis taper, the equatorial layer and the two radial edges is not covered.
`refute-u1` F10 is right that `u1` §6.7 prices this by assuming the answer.  Not done in the
record, not done here.

**M4.  (lam-cap): the a priori cap.**  (P3) proves `lam <= e^{3 c_*/4} = 1.8420112` from the
contradiction hypothesis alone.  A self-consistent argument should give `(3/2)(1 + O(1/L))`
(material on a shell where the taper is flat has `|omega^theta| = M lam_material <= (3/2)M`, and
`lam_material >= lam e^{-C_a c/L}`).  Not written.  *Priced:* `L_*` falls from `311511.9` to
`171527.1`, a factor `1.8161`, and `r_h` rises from `0.91863653` to `0.98148716`.

**M5.  (H5-Lip): the datum is Lipschitz in the angle, not `C^1`.**  `u2`'s (H5) asks for
`Gfrak_0 = sup |x|^2 |grad eta_0| < infinity` and (H1) for a bounded classical solution;
`u2` proves them for a `C^infinity` datum (D-C).  The theorem's datum has kinks at
`phi = delta`, `pi/2 - delta_m`, `pi/2` and their mirrors, so `Gfrak_0 = 58.695476` is an
**essential** supremum and `eta_0` is not `C^1`.  The repair is parabolic smoothing for `s > 0`
plus an approximation argument, or a change of datum.  *The change of datum is a real trade and
it is not free:* going to the fully-`tanh` (D-C) profile takes `Gfrak_0` from `58.695476` to
`20.919714` (worth about `2.8057` in `L_*`) and takes `kappa` from `0.49782244` to `0.47473436`,
i.e. `c_2` from `1.6289547` to `1.7081768` and the clock floor from `4.3741737e-03` to
`5.3220594e-02`.  The theorem's constant gets worse to make its hypothesis cleaner.

**M6.  (H\*): the solution class.**  The theorem is a statement about the `H^s` strong solution
on `[0,tau]`.  For this datum that is automatic (`omega_0` in `L^2` and `L^inf` with
`O(rho^{-8})` decay, so the solution is classical for `t > 0` and unique in
`C_w([0,T], L^inf)`), and it is the same class BFG's Theorem 10 asserts uniqueness in, so the
confrontation of §5 is between two statements about the same object.  Listed because it is a
hypothesis and not a theorem of this note.

**M7.  (Theta-tail): the shell support.**  The `tanh` ramps have tails, so `supp omega_0 = R^3`
with `O(rho^{-8})` decay at infinity and `O(rho^{8})` at the origin, not compact support.  This
is fine for BFG (`L^2 \cap L^inf` plus "a suitable decay") and it is fine for the Riesz route
(`u2` §6 needs only `grad eta` in `L^1 \cap L^inf`).  It is **not** literally the hypothesis
(D1) of `hk2` §2, which asks for support in `{rho_-(t) <= rho <= rho_+(t)}`, nor the shell
splitting of the far/near lemma at `2 rho`.  The tails are exponentially small in `L` and
nobody has written the estimate.

**M8.  `(Gamma-off)`'s own existence threshold.**  `C''` is a fixed point, not a constant:
`C'' = 2 Chat_a + (3 pi^2/16) e^{p c_G} Gfrak_0 + lambda` with `c_G = (lambda + C''/L) c_*`.  It
runs away below `L_Gamma_* = 14260.5` for this datum (`13888.0` on `{sin phi >= 1/2}`), against
`5273.2` for `u2`'s own datum, which reproduces `u2`'s reported `5313.8` to `0.8 %`.  Below
`L_Gamma_*` there is no `C''` and the theorem is empty.  This is a computed restriction rather
than an unproved statement, but `ASSEMBLY` and `u1` both carry `C'` as an `L`-independent
constant, so it has to be said.

**Explicitly discharged here, and no longer on the list:**
**(R-z)** (`u1`'s `G_collar`), by `u2` §6's global `r|grad a|` bound;
**(A-cone)** (`u1`'s cone restriction, `refute-u1` F1 and F2), by `u2` §7(b)'s axis-safe collar
together with the geometry `refute-u1` F3 keeps.  The supremum of `C_R` over **all** angles is
attained at `phi = 33.084 deg`, an interior point, so the cone restriction now costs exactly
nothing: `855.1552` over all `phi` against `855.1552` over `{phi >= delta}` and `851.1165` at
`phi_0 = 30.0 deg` itself.  `refute-u1`'s factor `3.2445` has evaporated.

---

## 4. THE BUDGET, RECOMPUTED

`t3_budget.py` is `u1`'s `u5_budget.py` (sha256 `4c2cdbc32cb0866a...`, itself a copy of
`ASSEMBLY`'s `a2_budget.py`) with `bootstrap()`, `COLUMNS` and `assemble()` replaced and the
datum constants re-read from `t1_results.json`.  `viscous_budget()`, `gauss_tail_aniso()`,
`J_pow()`, `eps_Tprime()` and the `RECORD` block are carried over unchanged.

### 4.1 What changed, term by term

| quantity | `u1` | here | why |
|---|---|---|---|
| `kappa_delta` | `0.49972123` (taper only) | `0.49782244` (full profile) | `refute-u1` F4, `u2` §9.1 |
| `c_*` | `0.8113826` | `0.8144774` | follows |
| `c_2` | `1.6227652` | `1.6289547` | follows |
| `lam_max` | `1.8377407` | `1.8420112` | follows |
| `r_h[1,lam_max]` | `0.9199394` | `0.91863653` | full profile, at the new `lam_max` |
| `C'` / `C''` | `151.15`, `L`-independent, `C_1 = 1` | fixed point, `1537.3` at `1e+05` and `1460.8` at `L -> infinity`, `C_1 = 0` | `u2` §7(a), §8 |
| `G_collar` in `C_R` | unproved, set to `0` | not used; `Ghat = 1250.5` at `L -> infinity` | `u2` §6 |
| angular sup | at `phi_0` (`u1`) or on the cone (`refute-u1`) | over **all** `phi` | `u2` §7(b) + `refute-u1` F3 |
| `C_R` | `29.024165` (`phi_0`), `98.289489` (cone), `100.6122` (F3, sharpened) | `855.1552` | above |
| `K_2` | `161.7735` | `161.7735` | `hk2`, unchanged |
| `log Lambda_*` shift | `3.364345456912714` | `2.9234859` | the `tanh` radial ramp's `C_E` |

`C_R` is dominated by the `r|grad a|` term: setting `Ghat = 0` takes it from `865.4507` to
`210.9797` at `L = L_*`, so **`75.6 %` of the drive is `Ghat`**, and `Ghat` is
`(3 pi^2/16) e^{p c_G} Gfrak_0` with `3 pi^2/16 = 1.8505508`, `p -> 2.0`, `c_G -> 1.2217172`.

### 4.2 The three columns

| column | `(Gamma-off)` | `C_a` | `Ghat` | `K_2/(M/rho_0)` | `P(|Z|>d/2)` | Lemma T' sensitivity |
|---|---|---|---|---|---|---|
| `proved` | `u2` fixed point, `sigma_* = 0` | far/near, `chat_a` | proved, `1257.3` at `1e+06` | `161.7735` | Theorem V.1 `B1` | as proved |
| `proved_Ca_measured` | same | `0.069` (R9, measured) | proved | `161.7735` | Theorem V.1 `B1` | as proved |
| `measured` | `C'' = 11.74` | `0.069` | `0.98142` (measured) | `5.3854` | exact affine Gaussian | `/ 6.0668328` |

All three at `delta = 7.5`, `delta_m = 5.0`, `phi_0 = 30.0`, `c = c_*`,
`lam_max = 1.8420112`, `r_h = 0.91863653`, `lam_om = M_s/M <= 1.5`, `C_1 = 0`, `f` minimised
over `{0.05, ..., 32}`.

### 4.3 `C''`, `Ghat` and `C_R` against `L` (proved column)

| `L` | `C''` | `c_G` | `p` | `Ghat` | `C_R` (sup over all `phi`) | `chat_a` sup | `chat_a(phi_0)` |
|---|---|---|---|---|---|---|---|
| `1e+03` | no fixed point | | | | | | |
| `1e+04` | no fixed point | | | | | | |
| `3e+04` | `1785.7897` | `1.2702` | `2.1004` | `1565.1227` | `1059.130` | `109.5835` | `13.046833` |
| `1e+05` | `1537.3117` | `1.2342` | `2.0262` | `1324.3538` | `899.574` | `105.7290` | `13.046833` |
| `3e+05` | `1484.9214` | `1.2257` | `2.0085` | `1273.7433` | `866.036` | `104.8391` | `13.046833` |
| `1e+06` | `1467.9116` | `1.2229` | `2.0025` | `1257.3246` | `855.155` | `104.5435` | `13.046833` |
| `1e+09` | `1460.8194` | `1.2217` | `2.0000` | `1250.4810` | `850.620` | `104.4192` | `13.046833` |

*Control (L-14, decorrelated).*  Driven with `u2`'s own datum norms `E_0 = 7.66060196`,
`Gfrak_0 = 20.9197070` and `u2`'s `c = 2 log(3/2)`, this instrument returns
`C''(1e+04) = 582.2400` on `{sin phi >= 1/2}` and `801.0451` on `R^5`, against `u2` §8.2's
`582.24` and `801.05` (relative `1.3e-08` and `6.1e-06`); `c_G = 1.2636109` and `1.2813545`
against `1.2636` and `1.2814`; `p = 2.1067843` and `2.1097300` against `2.1068` and `2.1097`.
`Chat_a` at `lambda = 1`, `sigma_* = 1/2` returns `8.6978888`, reproducing `ASSEMBLY` §2.4's
`8.697888`.  The Riesz composition constant `pi^4/2` is confirmed by an independent quadrature
to `1.7e-04`, and `C_K |S^4| = 1` exactly.

### 4.4 The term table at `c = c_*`, `f = 1`

| column | `L` | `eps_ell` | `eps_Ca` | `mu` | `mu_J` | `eps_T'` | `eps_v` | `eps_a` | `eps` |
|---|---|---|---|---|---|---|---|---|---|
| proved | `1e+03` | (no `(Gamma-off)` fixed point) | | | | | | | `inf` |
| proved | `1e+04` | (no `(Gamma-off)` fixed point) | | | | | | | `inf` |
| proved | `1e+05` | `3.6505e-05` | `2.6208e-04` | `0.090555` | `1.7392e-03` | `1.5004` | `9.0e-07` | `1.5007` | `inf` |
| proved | `3e+05` | `1.1756e-05` | `8.7359e-05` | `0.028914` | `5.7445e-04` | `0.34519` | `1.19e-07` | `0.34529` | `0.53407` |
| proved | `1e+06` | `3.4860e-06` | `2.6208e-05` | `8.551e-03` | `1.7181e-04` | `0.092034` | `1.73e-08` | `0.092064` | `0.10622` |
| proved_Ca_measured | `1e+05` | `3.6462e-05` | `1.386e-06` | `0.088417` | `1.6494e-05` | `1.4411` | `9.0e-07` | `1.4412` | `inf` |
| proved_Ca_measured | `1e+06` | `3.4855e-06` | `1.386e-07` | `8.3386e-03` | `1.6091e-06` | `0.089212` | `1.73e-08` | `0.089216` | `0.10276` |
| measured | `1e+03` | `3.5582e-03` | `1.386e-04` | `0.044643` | `1.6284e-03` | `0.095704` | `3.2147e-04` | `0.099401` | `0.11611` |
| measured | `1e+04` | `3.4778e-04` | `1.386e-05` | `4.4413e-03` | `1.6073e-04` | `7.7488e-03` | `4.0168e-06` | `8.1105e-03` | `0.012597` |
| measured | `1e+05` | `3.4698e-05` | `1.386e-06` | `4.439e-04` | `1.6052e-05` | `7.5914e-04` | `1.2567e-07` | `7.9523e-04` | `5.1738e-03` |
| measured | `1e+06` | `3.4690e-06` | `1.386e-07` | `4.4388e-05` | `1.605e-06` | `7.5759e-05` | `9.8125e-09` | `7.9366e-05` | `4.4539e-03` |

`E_tail` is exactly `0.0` in every row of this table.  The `120`-orders-of-magnitude swing
`ASSEMBLY` §5.3 item 4 flagged is a small-`L` artefact and it has left the budget entirely at
these `L`.  The binding term is `eps_T'`, through `mu`, in every column: `mu(c) L = 8499.62` at
`L -> infinity` in the proved column, from the closed form
`mu(c) L = lam_max^3 (C_R + 2 lam_om log lam_max)(lam_max^2 - 1)/(3/2)`, which the ODE solver
reproduces to `1e-06` at `L = 1e+06`.

### 4.5 `eps(L)` and `L_*`, the three columns

`eps` at `c = c_*`, minimised over `f`:

| column | `eps(1e+03)` | `eps(1e+04)` | `eps(1e+05)` | `eps(1e+06)` | best `f` at the first finite entry |
|---|---|---|---|---|---|
| all proved | `inf` | `inf` | `inf` | `0.10622` | `0.1` |
| proved, `C_a` measured | `inf` | `inf` | `inf` | `0.10276` | `0.1` |
| all measured | `0.11572` | `0.012553` | `5.1682e-03` | `4.4533e-03` | `0.5` |

`L_*` and `log Lambda_* = 2 L_* + 2.9234859`:

| column | `L_*` (`eps <= 1/2`) | `log Lambda_*` | `L_*` (`eps <= 0.1`) | `log Lambda_*` |
|---|---|---|---|---|
| **all proved** | **`311511.9`** | **`623026.8`** | `1056242.9` | `2112488.8` |
| proved, `C_a` measured | `302704.3` | `605411.5` | `1025043.0` | `2050088.9` |
| all measured | `422.4` | `847.8` | `1119.3` | `2241.5` |

`311511.9` is inside the `1e+04` to `1e+06` band the brief asks for.  The all-measured column
sits at `422.4`, below the `679 ... 1435` band `u2` quotes for the exact-plateau route's own
crossover, which means that on measured constants the clock is no longer the binding piece of
the campaign.

### 4.6 What dominates, and what a factor buys

| perturbation | `L_*` (proved, `eps <= 1/2`) | factor against the headline |
|---|---|---|
| headline | `311511.9` | `1` |
| `C_a -> 0.069` (measured), nothing else | `302704.3` | `0.9717261` |
| `lam_max = 3/2` (the self-consistent cap, M4) | `171527.1` | `0.5506277` |
| `C_K = 1` instead of `hk2`'s proved `161.7735` | `311511.5` | `0.9999986` |

For `C_R` the sensitivity is run against a base with `C_R` frozen at its `L = 1e+06` value
(`307839.3`), so that the multiplier is the only thing changing:

| `C_R` | `L_*` | factor against that base |
|---|---|---|
| times `0.5` | `155319.5` | `0.5045472` |
| times `1` | `307839.3` | `1` |
| times `2` | `612899.7` | `1.9909730` |
| times `5` | `1528097.1` | `4.9639442` |

Two things this says.  **`C_R` is in the drive, not in the exponent**: a factor `5` in it is a
factor `4.9639442` in `L_*`, not `e^5`.  That is `u1`'s structural gain and it survives intact.
And **`(H-K2)` is now worth `1.4e-06`**: `ASSEMBLY` §3.6 said closing it was necessary and not
sufficient because `e^{34}` swamped it; `u1` §6.6 said the reason had changed to "it is simply
small"; here it has become invisible.

---

## 5. THE BFG CONSEQUENCE

BFG Theorem 10 (arXiv:1704.05546v4, p.11, quoted in `ASSEMBLY` §4 from the text layer of the
PDF): *"Let the initial datum omega0 be in L2 cap Linf, and M a constant larger than 1.  Then
there is a constant c(M) > 1 such that there exists a unique mild solution omega in Cw[0,T],
Linf where `T >= (1/c(M))(1/||omega_0||_inf)`, and for any t in (0,T] the solution omega is the
R3-restriction of a holomorphic function ... moreover `||omega(t)||_{L^inf(Omega_t)} <= M
||omega_0||_inf`."*  Take `M = 3/2`: since `R^3` is contained in `Omega_t`, the clause gives
`||omega(t)||_{L^inf(R^3)} <= (3/2)||omega_0||_inf` for all `t <= T`, hence
`M_0 T_d(u_0) >= 1/c(3/2) > 0` for **every** datum with `omega_0` in `L^2` and `L^inf`, with
`c(3/2)` absolute: no dependence on `||omega_0||_2` (BFG's own sentence, p.8, says the `L^2`
assumption is soft and carries no quantitative dependence), none on `nu` (both sides are
invariant under `u -> lam u(lam x, lam^2 t)` at fixed `nu`), none on geometry.  The family of
§1.1 satisfies every hypothesis BFG state (`ASSEMBLY` §4.2, re-checked here: `||omega_0||_inf`
is finite and equal to `M` up to `O(e^{-2L/eps_r})`; `||omega_0||_2^2/(M^2 R^3) = 1.4785137`, so
`omega_0` is in `L^2` at every `L`; the decay at infinity is `O(rho^{-8})`; the solution class is
the one BFG assert uniqueness in, M6).  The assembled theorem of §1.3 gives
`M T_d <= c_2 (1+eps)/(2L) -> 0` as `L -> infinity` at fixed `M`.  The two cannot both hold.
**Therefore, conditional on exactly the modulo list of §3 and on nothing else, BFG Theorem 10 is
false as stated**, and the located defect in its displayed proof (the p.10 weighted-BMO step
without the local average, which two campaign referees called false and which nobody has
re-verified) is where to look.  Not triggered until M1 is written out; M1 is a routine argument;
that is the whole distance between this note and a refutation.

---

## 6. ABSTRACT

For axisymmetric no-swirl Navier-Stokes we assemble an upper bound on the vorticity doubling
time of a family of shell data. The family is a tapered, equatorially and radially mollified
vortex ring of amplitude M with shell aspect ratio L = log(R/rho_0) and inner scale
rho_0 = s sqrt(nu/M). For L above a threshold, the first time the vorticity maximum reaches
(3/2)M satisfies M T_d <= c_2 (1 + eps(L))/(2L), with c_2 = 2 log(3/2)/kappa_delta = 1.6289547,
and in energy Reynolds number T(Lambda) <= c_2 (1 + eps)/log Lambda for log Lambda above
623026.8. The error eps(L) decreases to 4.3741737e-03, the strain deficit of the mollified
profile. The assembly is complete modulo eight named items, of which one, the
approximate-maximum step in two weighted maximum principles, is routine and carries all the
others. Granting them, the bound contradicts Bradshaw, Farhat and Grujic's Theorem 10, whose
lower bound on the same quantity is absolute.

---

## 7. WHAT A REFEREE SHOULD ATTACK FIRST

1. **M1, the approximate-maximum step, and specifically whether the supremum can escape.**  P1
   and P2 are stated as suprema over all of `R^5` of `rho|eta|` and `rho^2|grad eta|`.  The
   Grönwall step needs the supremum to be attained, or the Dini derivative argument to be run on
   a family that exhausts it; `eta` is only bounded, and the weights `rho` and `rho^2` grow.
   `u2`'s proof asserts the decay `Theta(rho) = O(rho^{-8})` is preserved "up to Gaussian tails"
   by the heat-transport semigroup, which is exactly the point at issue: the transport term
   `b . grad_5` moves mass outward at rate `Gamma_rad ~ (3/4) M L`, and `rho^2` is being asked to
   stay integrable against it.  If the supremum escapes to infinity, P2 fails, `Ghat` fails,
   `(Gamma-off)` fails, and nothing here stands.  This is where I would spend the first day.

2. **Whether `C_R` is really the supremum, given that `mu` and `mu_J` are suprema over the
   *shell* while the pieces of `R` are estimated at a *field point*.**  §4.3's `chat_a` supremum
   is `104.5435`, taken at the axis where the collar branch caps; `chat_a(phi_0) = 13.046833`.
   The Jacobian equation uses the shell supremum (`refute-u1` F14 is right that `u1` did not),
   and the drive uses `C_R`'s supremum over `phi`.  But the `z`-integrals in `R_z` run along a
   segment from the equator to the field point, sweeping the whole angular range, and I bound
   them by integrating `chat_a(phi_zeta)` along that segment.  A referee should check that this
   is the right object and not an over-count of the collar, which is where all of `chat_a`'s
   size lives.  The answer changes `C_R` linearly, hence `L_*` linearly, so it is worth a factor
   and not an exponential.

3. **M5 and M2 together, because they pull in opposite directions and the theorem's datum sits
   between them.**  `hk2` forces a `tanh` radial ramp or `K_2` is infinite.  `ASSEMBLY` BLOCK 5
   chose `min(1,.)` angular profiles so that (H3) holds exactly away from the kinks.  `u2`'s P1
   and P2 want a `C^infinity` datum.  The `kappa_delta` of the fully-`tanh` datum is
   `0.47473436` against `0.49782244`, which moves the theorem's own constant `c_2` from
   `1.6289547` to `1.7081768` and the floor from `4.3741737e-03` to `5.3220594e-02`.  There may
   be no datum that satisfies all three at once, and if there is not, the theorem has to name
   which one it gives up.  Nobody in the chain has checked this, and I have not either.

---

## GATE, AND FILES

```
89 CHECKS, 89 PASS, 0 FAIL
```

`check_constants.py` rebuilds every displayed quantity from the stored JSONs and from
mathematics (never by comparing a literal with itself) and then audits this document token by
token.  Constants quoted from a named seat are listed in its `FOREIGN` dictionary with their
origin and appear only in comparison columns.

What the gate **cannot** catch, stated so it is not mistaken for coverage: it checks that the
numbers displayed are the numbers the scripts produced and that the algebra between them is
consistent.  It does not check that (P-max) is true, that `C_R` is the right constant, that the
modulo list is complete, or that `u2`'s P1 and P2 are the right instruments.  That is what §7 is
for.

| file | what it establishes |
|---|---|
| `t1_datum.py` / `t1_results.json` | the datum; `kappa` for four profiles; `c_*`, `c_2`, `lam_max`, `r_h`, `P_h` monotone; `E_0` and `Gfrak_0` for the theorem's datum and for (D-C); the energy identity, `H_l`, `D_l`, `C_E` for sharp and `tanh` radial ramps and the `log Re_E` shift; the BFG `L^2` check |
| `t2_gamma_CR.py` / `t2_results.json` | the exact kernel constants (`C_K|S^4| = 1`, Riesz `pi^4/2` two ways, `3 pi^2/16`, `2 R_A`); `u2`'s `(Gamma-off)` fixed point re-implemented and controlled against `u2`'s own numbers; `L_Gamma_*`; `C_R(phi)` with the geometry kept and the axis-safe collar, and its supremum over every angle |
| `t3_budget.py` / `t3_results.json` | the budget: `u5_budget.py` with `bootstrap`, `COLUMNS`, `assemble` replaced; the term table, `eps(L)`, `L_*`, `log Lambda_*`, the dominance decomposition, the sensitivities |
| `t4_derived.py` / `t4_results.json` | the ratios, percentages and differences quoted inline in this document's prose, so that each is traceable |
| `check_constants.py` | the gate and the document audit |
| `copies/` | byte copies, with hashes, of every file imported from another seat |
| `SHA256SUMS` | computed with `shasum -a 256`, never typed |

**FL-000 stands.  Nothing here touches the headline problem.**

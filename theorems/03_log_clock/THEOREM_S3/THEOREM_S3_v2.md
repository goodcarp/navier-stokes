# THEOREM (S3), version 2: the sharp logarithmic doubling clock, on a smooth datum

Seat `s3close/round2/THEOREM_S3/fix4`, sitting of 2026-09-08. This file restates the theorem
for the `sigma`-mollified datum. It supersedes nothing: `THEOREM_S3.md` and its three addenda
are unchanged and remain the record of the kinked-datum version. Every constant below comes
from a script in `fix4/`, listed in `fix4/FIX4.md` and re-asserted by `fix4/check_fix4.py`.

**What changed against `THEOREM_S3.md` plus `ADDENDUM_1/2/3`.** The angular profile
`min(1, phi_ax/delta) min(1, |phi - pi/2|/delta_m)` is replaced by its Gaussian mollification at
`sigma = 0.02` rad, normalised so that `||omega_0^theta||_inf` is still `M`. That closes M3 and
M5 together and leaves the modulo list with no unproved analytic item. It costs `kappa_delta`
and buys `Gfrak_0`; the net on the reach is a factor `1.2711978138` in the right direction.

**FL-000 stands.** This is not a proof of the headline problem and does not bear on it.

---

## 1. THE STATEMENT

### 1.1 The datum

Fix `M > 0`, `nu > 0`, `delta = 7.5 deg`, `delta_m = 5.0 deg`, `phi_0 = 30.0 deg`,
`s := 1/sin delta = 7.66129757554039`, `eps_r = 0.25`, `sigma = 0.02` rad, `f >= 4`, `L > 0`.
Write `u := log(rho/rho_0)`, `phi_ax := min(phi, pi - phi)`. Put

```
    nu     := M rho_0^2/s^2       (equivalently rho_0 = s sqrt(nu/M)) ,     R := rho_0 e^L ,

    Theta(rho) := (1/2)[ tanh((u - eps_r)/eps_r) - tanh((u - L + eps_r)/eps_r) ] ,

    A(phi)  := sgn(cos phi) min(1, phi_ax/delta) min(1, |phi - pi/2|/delta_m) ,
    W(phi)  := A(phi)/sin phi ,
    Wtil    := W extended to R by even reflection at phi = 0 and phi = pi (2 pi periodic) ,
    chi_sigma(y) := exp(-y^2/(2 sigma^2))/(sigma sqrt(2 pi)) ,
    W_sigma := (Wtil * chi_sigma)|_{[0,pi]} ,
    A_sigma := W_sigma(phi) sin phi ,        N_sigma := sup_{[0,pi]} |A_sigma| = 1.0124508488 ,

    omega_0^theta(rho, phi) := - M Theta(rho) A_sigma(phi)/N_sigma ,
    omega_0 = omega_0^theta e_theta ,        u_0 = Biot-Savart of omega_0 ,
    eta_0   = omega_0^theta/r = -(M/rho) Theta(u) W_sigma(phi)/N_sigma ,
```

axisymmetric, no swirl, odd in `z`. The tracked material point is `x_* = (rho_*, phi_0)` with
`rho_* := (1+f) rho_0`; `L := log(R/rho_0)`; `T_d(u_0) :=` the first time
`||omega(t)||_inf >= (3/2) M`.

`W_sigma` is real-analytic on `[0, pi]` and even at both poles, so `eta_0` is `C^infinity` in
the angle including the axis; `omega_0` is `C^infinity` on `R^3 \ {0}` and `C^{7,1}` at the
origin, lies in `L^2 cap L^inf` with `||omega_0||_inf = M sup Theta = M(1 - O(e^{-2L/eps_r}))`,
decays like `rho^{-8}` at infinity and like `rho^8` at the origin, and belongs to `H^s(R^3)` for
every `s < 9.5`; hence `u_0 in H^s` for every `s < 10.5`.

The normalisation by `N_sigma` is not cosmetic: on the bulk the kinked `W` equals `1/sin phi`,
which is convex, so `W_sigma >= W` there and the unnormalised `A_sigma` exceeds `1` by
`1.2450848776 %`. Dividing by `N_sigma` is what makes `M` the amplitude the clock is about.

### 1.2 The constants of the datum

```
    kappa_delta  = 0.4917868801        (kinked: 0.4978224383)
    c_*          = log(3/2)/kappa_delta = 0.8244732108
    c_2          = 2 log(3/2)/kappa_delta = 1.6489464217
    lam_max      = e^{3 c_*/4}         = 1.8558724485
    r_h          = inf_{[1, lam_max]} P_h(lam)/lam  >=  0.9033435309
    E_0          = sup rho|eta_0|/M    = 7.5523076942        (kinked: 7.6612975755)
    Gfrak_0      = sup rho^2|grad eta_0|/M = 36.0795702396   (kinked: 58.6954805410)
    Psi_sigma(0) = sup rho|Lap_5 eta_0|/M  = 489.0279393839  (kinked: +infinity)
    C_kink       = sigma Psi_sigma(0)  = 9.7805587877
    C_E          = ||u_0||_2^2/(M^2 R^5) = 0.0551906693
    ||omega_0||_2^2/(M^2 R^3)          = 1.4381082627
```

### 1.3 The theorem

> **THEOREM (S3), version 2.** Let `u_0` be the datum of sec.1.1, with `nu > 0` and `f >= 4`,
> and let `u` be its maximal strong solution, which exists and is unique by Proposition K of
> `forced-route-2026-09/theorems/02_ADDENDUM_uniqueness_2026-09-08.md` sec.3 (item M6 below).
> Assume, on `[0, c/(ML)]` with `c` the self-consistent window fixed below:
>
> * **(lam-cap)** nothing. `lam <= e^{3c/4}` is proved from the contradiction hypothesis alone;
>   it is listed in sec.2 only because a self-consistent argument would improve it.
> * `L >= L_Gamma^exist`, so that `(Gamma-off)`'s own fixed point exists:
>   `L_Gamma^exist = 9298.1782` at `c = c_*` and `16231.8518` at the certified window cap.
>
> **No further hypothesis.** Every analytic item that `THEOREM_S3.md` sec.3 carried as unproved
> is discharged: (P-max) and (H3-V) in `ADDENDUM_1`, (Theta-tail) in `ADDENDUM_3`, (V-strain)
> and (H5-Lip) here.
>
> Then, for every `L >= L_*`, there is `c` in `[c_*, 1.1664585076 c_*]` with
> `c >= c_*(1 + eps(c, L))`, and with `eps(L) := c/c_* - 1` for the smallest such `c`,
> ```
>     T_d(u_0)  <=  t_*  :=  c/(M L)  =  c_* (1 + eps(L))/(M L)
>                        =   2 log(3/2) (1 + eps(L))/(2 kappa_delta M L) ,
> ```
> with `c_* = 0.8244732108` and `c_2 = 1.6489464217`, and `eps(L)` decreasing to the clock floor
> `(1 - 2 kappa_delta)/(2 kappa_delta) = 1.6700567265e-02`. The reach is
> ```
>     L_*  =  1424610.4953  (all proved) ,  1380034.8075  (proved, C_a measured) ,
>             1757.8779     (all measured)                    at  eps <= 0.1664585076 ,
>     L_*  =  1555469.4004  (all proved) ,  1500454.5140  (proved, C_a measured) ,
>             2046.8908     (all measured)                    at  eps <= 0.1 .
> ```
>
> **In energy Reynolds number.** With `Re_E = E_0^{2/5} M_0^{1/5}/nu`, `E_0 = ||u_0||_2^2`, and
> `log Re_E = 2L + 2 log s + (2/5) log C_E = 2L + 2.9135781820`,
> ```
>     T(Lambda)  <=  c_2 (1 + eps(Lambda))/log Lambda ,   c_2 = 1.6489464217 ,
>     for  log Lambda >= log Lambda_* = 2 L_* + 2.9135781820 ,
>     log Lambda_*  =  2849223.9041 / 2760072.5285 / 3518.6693   at eps <= 0.1664585076 ,
>                      3110941.7143 / 3000911.9416 / 4096.6951   at eps <= 0.1 .
> ```

**The window cap.** A2/C6 need `P_h` increasing on `[1, lam_max]`. `P_h'` vanishes at
`lam_mono = 2.0693384635`, so the window may not exceed `c/c_* = 1.1760705119`; the largest
window for which the certificate of `fix4` sec.2.1 closes is `c/c_* = 1.1664585076`, i.e.
`eps <= 0.1664585076`. `fix2`'s cap for the kinked profile was `0.1943662079`; the cap moves
down with the datum, and the `eps <= 1/2` column withdrawn by `ADDENDUM_2` stays withdrawn.

---

## 2. THE DEPENDENCY GRAPH

`PROVED` means proved in the named file with its hypotheses met here. `ENCLOSED` means the
statement is proved and the constant is a certified numerical enclosure, not a closed form.
`SLACK` means the item is a proved but improvable constant, not a gap.

| # | statement | file | status |
|---|---|---|---|
| A1 | the datum; `N_sigma`; `omega_0` in `L^2 cap L^inf`; the decay and regularity | `fix4/p1`, `fix4/p2` | **ENCLOSED** |
| A2 | `kappa_delta`, `c_*`, `c_2`, `lam_max`; `P_1(lam) = lam`; `mu_lam` a probability measure | `fix2/f4` (structure) + `fix4/p1` | **PROVED** structure, **ENCLOSED** value |
| A3 | `r_h >= 0.9033435309` and `P_h' >= 0.3335513465 > 0` on `[1, lam_max]`; the cap `0.1664585076` | `fix4/p1` | **ENCLOSED** (cell-wise Lipschitz certificate with closed-form primitives) |
| A4 | `E_0`, `Gfrak_0` by `fix2` sec.4(c)'s branch bound | `fix2/f4` (structure) + `fix4/p1` | **PROVED** structure, **ENCLOSED** value |
| A5 | the energy identity, `G_l`, `D_l`, `C_E`, the `L` to `log Re_E` dictionary | `THEOREM_S3/t1` (structure) + `fix4/p2` | **ENCLOSED** (three controls: orthogonality, Parseval, synthesis) |
| B1 | Lemma 2: `eta <= 0` on `{z > 0}` for all `t`, all `nu >= 0` | `lower/prove-lagrangian` sec.3 | **PROVED** |
| B2 | far/near kernel lemma; Prop. 1, Prop. 2, Consequence A; `0.291999`, `0.014754`, `3.999218`, `pi/8` | `rebuild/far-near-kernel-lemma` | **PROVED**; hypotheses re-verified for this datum in `fix4/p3` |
| B3 | Lemma T' and Corollary 2; the sensitivity `eps_T'` | `write/lemma-T-shell-dependent` | **PROVED** |
| B4 | Theorem V.4' with `(H3'_vartheta)` for `k = 0..4`; `vartheta(f = 4) = 0.0556` | `round2/pmax-h3v` Part B (`ADDENDUM_1`) | **PROVED** at `f >= 4` |
| B5 | Theorem V.1, `B1` | `gaps/gap-V-aronson` | **PROVED** |
| B6 | `(H-K2)`: `K_2 <= 161.7735 M/rho_0`, no `log`, for the `tanh`-radially-mollified datum | `s3close/hk2` sec.0, with `(D1')` of `ADDENDUM_3` | **PROVED**; the majorant re-checked for this datum in `fix4/p3`, slack `3.9889277977` |
| C1 | the slaved reference `frak_a`; (P1), (P2), (P3) | `round2/u1` sec.1 | **PROVED** |
| C2 | `Lambda` is the flow of `frak_a . D X`; `J_Lambda` | `u1` sec.1.2 | **PROVED** |
| C3 | the `l = 1` rate error is identically zero | `u1` sec.2.1 | **PROVED** |
| C4 | `frak_a` and Lemma T's `a_ref` are the same functional | `u1` sec.4.1 | **PROVED** |
| C5 | the restriction step (positivity plus `z`-oddness) | `u1` sec.4.2 step 2 | **PROVED** |
| C6 | the integro-differential inequality and the quasimonotone comparison | `u1` sec.5.1 + A3 | **PROVED** given A3 |
| C7 | the conservative clock `c_2 = 2 log(3/2)/kappa_delta` | `u1` sec.5.3 | **PROVED** |
| D1 | `\|\|grad u\|\|_op <= 2\|a\| + r\|grad a\| + \|omega^theta\|` | `u2` sec.2 | **PROVED** |
| D2 | `b(0) = 0`; `x . b <= Gamma_rad \|x\|^2` | `u2` sec.4 | **PROVED** |
| D3 | **P1** `sup rho\|eta\| <= e^{c_R} E_0 M`, **P2** `sup rho^2\|grad eta\| <= e^{p c_G} Gfrak_0 M` | `round2/pmax-h3v` Part A (`ADDENDUM_1`) | **PROVED** (interior parabolic Schauder cited, not proved) |
| D4 | `r\|grad a\| <= (3 pi^2/16) e^{p c_G} Gfrak_0 M sin phi`; Riesz `pi^4/2` | `u2` sec.6 | **PROVED**; `pi^4/2` cross-checked to `1.66e-04` in `fix4/p3` |
| D5 | `C_1 = 0` | `u2` sec.7(a) | **PROVED** |
| D6 | the axis-safe collar `\|a_collar\| <= 2 R_A e^{c_G} E_0 M` | `u2` sec.7(b) | **PROVED** |
| D7 | `(Gamma-off)` `Gamma(s) <= 2 a(0,s) + C'' M`; `C''` the fixed point; existence by sign change | `u2` sec.8 + `fix2/f5` + `fix4/p4` | **PROVED** for `L >= L_Gamma^exist`; `C''` **ENCLOSED** |
| E1 | `\|R_y\| <= chat_a(phi) M rho sin phi` | `refute-u1` F3(a) + B2 + D6 | **PROVED** |
| E2 | `int_0^{\|z\|} log(rho/rho_zeta) dzeta = \|z\| - r arctan(\|z\|/r)` | `refute-u1` F3(b) | **PROVED** |
| E3 | `\|R_z\|`, with the `r\|grad a\|` integral; **(R-z) discharged** | `THEOREM_S3/t2` + D4 | **PROVED** |
| E4 | `C_R = sup_phi`, finite over every angle; **(A-cone) discharged** | `THEOREM_S3/t2` + `fix4/p4` | **PROVED**, value **ENCLOSED** |
| F1 | the majorant ODE, `mu = lam_max^2 m` | `u1` sec.2.3 | **PROVED** given D7, E4 |
| F2 | the Jacobian equation, `mu_J`, Lemma T' (H3) | `u1` sec.2.3 | **PROVED** |
| F3 | Lemma T' applied once, with `r_h` on `[1, lam_max]` | `u1` sec.4.2 step 4 | **PROVED** given A3 |
| F4 | `eps_v = eps_bulk + E_hess + E_4 + E_tail` | `u1` sec.4.2 step 5, instrument from `a2_budget.py` | **PROVED** given B4 |
| F5 | **(V-strain)**: the viscous defect of the strain functional, `eps_a` | `fix3` sec.1 + `fix4/p1`, `fix4/p5` | **PROVED**; `Psi_sigma(0)` **ENCLOSED** |
| G1 | the self-consistent window; `eps(L)`, `L_*`, `log Lambda_*` | `fix2/f1` (instrument) + `fix4/p4` | **COMPUTED**, conditional on the ENCLOSED rows |
| G2 | the family satisfies BFG's hypotheses | sec.4 below, `fix4/p2` | **PROVED** |
| G3 | BFG Theorem 10 is false as stated | sec.4 below | **CONDITIONAL** on G1 and on the ENCLOSED rows |
| M4 | `(lam-cap)`: the a priori cap `e^{3c/4}` | `u1` (P3) | **SLACK** (proved; a self-consistent argument would improve it) |

---

## 3. THE REMAINING ENCLOSED STEPS, LISTED

Nothing in the chain is unproved. What remains is a list of constants that are certified
numerical enclosures rather than closed forms, plus two citations.

1. **`r_h` and `P_h` increasing** (`fix4/p1`). Proved statements; the certificate is cell-wise,
   with closed-form primitives, and its only non-exact input is the evaluation of the
   real-analytic `h_sigma` on a `600001`-point grid, padded by the proved Lipschitz bound
   `|h_sigma''| <= 225.9466903722`. `r_h >= 0.9033435309` against a grid minimum
   `0.9035285738`.
2. **`E_0`, `Gfrak_0`, `Psi_sigma(0)`, `N_sigma`** (`fix4/p1`). Suprema of real-analytic
   functions over the same grid; `Gfrak_0` additionally carries `fix2` sec.4(c)'s branch bound,
   which is a rigorous upper bound over all `u` at once and which agrees with the direct scan
   here to the last digit.
3. **`C_E`** (`fix4/p2`). A Gegenbauer sum plus a radial recursion, converged to `7.3e-06` in
   `lmax` and `4.7e-07` in `du`, with three independent controls.
4. **`C''`, `Ghat`, `C_R`, `chat_a`, `L_Gamma^exist`** (`fix4/p4`). The existence of the
   `(Gamma-off)` root is proved by a sign change (intermediate value theorem); its value and
   `C_R`'s supremum over the angle are quadratures.
5. **`K_2`** (`fix4/p3`). A quadrature, used only through the imported majorant `161.7735`,
   which it clears by `3.9889277977`; and `(H-K2)` is worth `6.0847127e-09` of `eps`, so the
   enclosure is not load-bearing.
6. **The far/near constants** `C1 = 0.291999`, `C2in = 0.014754`, `C_collar = 3.999218/sin phi
   + pi/8` are Parseval-based numerical evaluations in `rebuild/far-near-kernel-lemma`; their
   three hypotheses are re-verified for this datum in `fix4/p3`.
7. **Interior parabolic Schauder** is cited, not proved, in `round2/pmax-h3v` Part A. That is
   the one textbook citation the P1/P2 chain rests on.
8. **The CKN attribution is NOT used by this theorem.** Proposition E and Corollary C2 of
   `02_ADDENDUM_uniqueness_2026-09-08.md`, which carry the unverified attribution of the forced
   epsilon-regularity criterion to CKN 1982, are needed only to rule out a *general Clay-class*
   competitor. This theorem is a statement about the maximal strong solution supplied by
   Proposition K, which is PROVED there by in-estate citations, and the BFG confrontation of
   sec.4 is between two statements about that same object. If a referee insists on the
   Clay-class formulation, then Proposition E enters and that attribution becomes an ENCLOSED
   step; it is recorded here so the choice is visible.

---

## 4. THE BFG CONSEQUENCE, FOR THE RECORD

Bradshaw, Farhat and Grujic, arXiv:1704.05546v4, Theorem 10, p.11. The clause, stated rather
than quoted (the verbatim text is at `THEOREM_S3.md` sec.5, transcribed there from the text
layer of the PDF): for a datum `omega_0` in `L^2 cap L^inf` and any constant `M > 1`, there is
`c(M) > 1` such that a unique mild solution exists in `C_w([0,T], L^inf)` on a time interval of
length at least `(1/c(M)) (1/||omega_0||_inf)`, the solution is on `(0,T]` the restriction to
`R^3` of a holomorphic function, and its `L^inf` norm on the corresponding domain stays below
`M ||omega_0||_inf`. Taking `M = 3/2` and using that `R^3` sits inside that domain gives
`||omega(t)||_{L^inf(R^3)} <= (3/2)||omega_0||_inf` for all `t <= T`, hence

```
     M_0 T_d(u_0)  >=  1/c(3/2)  >  0     for EVERY datum with omega_0 in L^2 cap L^inf,
```

with `c(3/2)` absolute: no dependence on `||omega_0||_2` (BFG's own sentence on p.8 says the
`L^2` assumption is soft and carries no quantitative dependence), none on `nu` (both sides are
invariant under `u -> lam u(lam x, lam^2 t)` at fixed `nu`), none on geometry.

**The family of sec.1.1 satisfies each of BFG's hypotheses, re-checked for the mollified
profile** (`fix4/p2`, `fix4/p3`):

* `||omega_0||_inf = M sup Theta = M(1 - O(e^{-2L/eps_r}))`, finite and equal to `M` up to that
  correction, **by construction**: the normalisation by `N_sigma = 1.0124508488` is exactly what
  makes this true, and `sup |A_sigma|/N_sigma = 1.0` exactly.
* `omega_0 in L^2`: `||omega_0||_2^2/(M^2 R^3) = 1.4381082627`, finite at every `L`.
* decay: `O(rho^{-8})` at infinity and `O(rho^8)` at the origin, so the "suitable decay" BFG ask
  for is present with room.
* the solution class: BFG assert uniqueness in `C_w([0,T], L^inf)`; the solution here is the
  maximal strong solution of Proposition K, which lies in that class, so the two statements are
  about the same object (item M6 of `fix4/FIX4.md` sec.7).

The theorem of sec.1.3 gives `M T_d <= c_2 (1 + eps)/(2L) -> 0` as `L -> infinity` at fixed `M`.
The two cannot both hold. **Conditional on the ENCLOSED list of sec.3 and on nothing else, BFG
Theorem 10 is false as stated**, and the located defect in its displayed proof (the p.10
weighted-BMO step without the local average, which two campaign referees called false and which
nobody has re-verified) is where to look.

*What is different from `THEOREM_S3.md` sec.5.* There, the trigger was "M1 is a routine argument
and that is the whole distance". M1 has since been written (`ADDENDUM_1`), M2 and M7 with it,
and M3 and M5 are closed here. The distance is now not an unproved analytic step but the list of
enclosures in sec.3, the two citations in it, and a referee's reading of the chain. That is a
different kind of gap and a smaller one, and it is the reason this file exists. It is not a
public claim: the estate's standing rule is that the refutation is not asserted until the chain
has been refereed by a seat that did not build it.

---

## 5. WHAT A REFEREE ATTACKS FIRST

1. **The normalisation, and whether `kappa_delta` is still the right strain constant.** The
   whole cost of this version is that `kappa_delta` falls from `0.4978224383` to `0.4917868801`,
   which raises the clock floor by a factor `3.8179938670`, and the whole gain is that
   `Gfrak_0` falls by `1.6268342260`. Both come from the same operation. A referee should check
   (a) that mollifying `W = A/sin phi` rather than `A` is the right operation, given that it is
   `W` that must be even at the poles for `eta_0` to be smooth there, and (b) that dividing by
   `N_sigma` is the right response to the convexity overshoot rather than, say, mollifying `A`
   with an odd reflection, which would keep `|A| <= 1` automatically but has to be checked for
   axis smoothness. The two choices give different `kappa_delta` and the difference is
   comparable to `1 - 2 kappa_delta` itself, which is the most delicate number in the theorem.
   This is where I would spend the first day.

2. **Whether `C_R` is really the supremum, given that `mu` and `mu_J` are suprema over the
   shell while the pieces of `R` are estimated at a field point.** Unchanged from
   `THEOREM_S3.md` sec.7 item 2 and unaffected by the mollification: `chat_a`'s size lives at
   the axis (`104.4884` against `13.046833` at `phi_0`), the `z`-integrals in `R_z` sweep the
   whole angular range, and the bound integrates `chat_a(phi_zeta)` along that segment. The
   answer changes `C_R` linearly, hence `L_*` linearly.

3. **The two enclosures that are not merely tight but structural: `r_h` and the window cap.**
   `r_h` enters `eps_T'` as `1/r_h` and the certified value here, `0.9033435309`, is `1.7 %`
   below the kinked datum's `0.9186364112` because `lam_max` grew with `c_*`. The cap
   `eps <= 0.1664585076` is not a convenience: beyond it `P_h` is not increasing and C6's
   quasimonotone comparison has no hypothesis at all. A referee should check that the cell-wise
   certificate is genuinely a bound and not a sampled statement, and that the window
   self-consistency of `fix2` item 1 is imposed at the same `c` everywhere it is used.

---

## 6. ABSTRACT

For axisymmetric no-swirl Navier-Stokes we bound the vorticity doubling time of a family of
shell data. The family is a smooth, equatorially tapered and radially mollified vortex ring of
amplitude `M`, shell aspect ratio `L = log(R/rho_0)` and inner scale `rho_0 = s sqrt(nu/M)`; its
angular profile is the Gaussian mollification, at scale `sigma = 0.02` rad, of a piecewise
linear taper, renormalised so that the vorticity maximum is exactly `M`. For `L` above a
threshold, the first time the vorticity maximum reaches `(3/2)M` satisfies
`M T_d <= c_2 (1 + eps(L))/(2L)` with `c_2 = 2 log(3/2)/kappa_delta = 1.6489464217`, and in
energy Reynolds number `T(Lambda) <= c_2 (1 + eps)/log Lambda` for `log Lambda` above
`2849223.9041`. The error `eps(L)` decreases to `1.6700567265e-02`, the strain deficit of the
mollified profile. Every analytic hypothesis the earlier assembly carried is now discharged;
what remains is a list of certified numerical enclosures and two citations. Granting them, the
bound contradicts Bradshaw, Farhat and Grujic's Theorem 10, whose lower bound on the same
quantity is absolute.

**FL-000 stands. Nothing here touches the headline problem.**

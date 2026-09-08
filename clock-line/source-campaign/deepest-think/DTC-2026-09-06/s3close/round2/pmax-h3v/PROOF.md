# (P-max) and (H3-V): items M1 and M2 of THEOREM_S3's modulo list

Seat `s3close/round2/pmax-h3v`, DTC-2026-09-06, sitting of 2026-09-08.
Laws: `TORMENT NEXUS/LAWS.md`, first 120 lines read.
FL-000 will fall; nothing in this note touches it.

Every number below came out of a script in **this** folder that I wrote and ran
(`q1_identities.py`, `q2_datum.py`, `q3_theta.py`, `q4_budget_theta.py`), re-asserted by
`check_constants.py`. `SHA256SUMS` is computed with `shasum -a 256`, never typed. Nothing
outside this folder was written.

**Read, in the order the brief set** (L-14): `round2/THEOREM_S3/THEOREM_S3.md` in full;
`round2/u2/PROOF.md` sections 0 to 6; `round2/refute-u2/NOTE.md` in full;
`write/V-b-bulk-viscous-loss/PROOF.md` sections 1 to 5;
`write/refute-V-b-bulk-viscous-loss/PROOF.md` section 4 and section 7;
`s3close/assembly/ASSEMBLY.md` sections 2.5 and 2.6.

**Code imported, L-14 inverse convention declared.** For TASK B the object under test is
THEOREM_S3's own budget instrument, so it is imported byte for byte and not re-implemented:
`copies/t3_budget.py`, `copies/t2_gamma_CR.py`, `copies/t1_results.json`, with hashes in
`copies/IMPORTED_SHA256SUMS.txt`. `copies/u5_budget.py` is there as a receipt: its sha256 is
`4c2cdbc32cb0866a894f6e1ba8853a2e5e81ba47677d3c6bea084cd40854c7b0`, which is the value
THEOREM_S3 section 4 quotes for its own provenance chain. Nothing in any imported file is
edited; the `vartheta` terms are added on top, in `q4_budget_theta.py`.
Everything I am testing a *claim* about (the datum constants `E_0` and `Gfrak_0`, the weight
identities, the barrier, `vartheta`) is re-derived here from the mathematics.

---

## 0. HEADLINE

**M1 (P-max) closes.** P1 and P2 are theorems with no gap, and the referee's first attack
point closes with them.

1. The **escape of the supremum to infinity** (refute-u2 section 7 item 1: *"this is where I
   would spend the first day"*) is not a gap and does not need the datum's tails. The
   comparison runs against an explicit barrier `h = e^{Ks}(rho^2+1)^{beta/2}`, and the only
   growth information it needs is `rho|eta| <= ||eta_0||_inf rho` for P1 and
   `rho^2|grad eta| <= C rho^2` for P2, both of which are free. Lemma E case (ii).
2. The **two unstated hypotheses** of refute-u2 MINOR-4 stop being hypotheses. The decay
   `rho|eta| -> 0`, `rho^2|grad eta| -> 0` is Lemma F, proved from the datum by a second
   barrier with explicit rate `Lambda_9 = 9 Gamma-bar + 99 nu/R^2`. `grad eta in L^1` is
   Corollary F', with an explicit bound; it is used in exactly one place in the whole chain,
   `u2` section 6's `grad a = K * grad_5 eta`, and nowhere in P1 or P2.
3. The **mollification is exact**. Replacing `|.|` by `sqrt(.^2 + eps^2) - eps` and `rho` by
   `R_eps = sqrt(rho^2 + eps^2)` produces the two differential inequalities with **no error
   term at all** in P1 and a single `Gamma eps R_eps^2` source in P2, and the zeroth-order
   weight coefficients are `-(2 rho^2 + 5 eps^2)/R_eps^4` and `-(2 rho^2 + 10 eps^2)/R_eps^4`,
   both `<= 0` for every `eps`. At `eps = 0` these are `k(k-3)/rho^2 = -2/rho^2` for `k = 1`
   and `k = 2`, so the brief's sign check holds, and it holds with room: `k(k-n+2) <= 0` for
   every `0 <= k <= n-2`, which at `n = 5` is `k <= 3`. Verified symbolically for
   `n = 2 ... 7` and `k = 0 ... 4` (`q1`, 49 checks, 0 failures).
4. **The datum's two constants have closed forms.**
   ```
      E_0     = 1/sin(delta)    = 7.661297575540389        (= s, the campaign's own s)
      Gfrak_0 = 1/sin(delta)^2  = 58.69548054098104  = E_0^2
   ```
   Both exact, both attained at the taper kinks `phi = delta` and `phi = pi - delta`. On the
   plateau arc `F^2 + F'^2 = (sin^2 + cos^2)/sin^4 = 1/sin^4 phi`, and the arc's closest
   approach to the axis is `phi = delta`; that identity is the whole derivation.
   THEOREM_S3's `58.695476` is this number to `7.2e-08` and refute-u2's `58.6948961` to
   `1.0e-05`; both are grid values of the same closed form.
5. **Brief-layer correction (L-17).** The brief expects `E_0 ~ 4.09 M`. That number is
   `||eta_0||_inf rho_0/M = 4.093544431740673`, a **different functional** (it is
   refute-u2 MINOR-2's quantity, computed here for the theorem's datum rather than for
   (D-C), where refute-u2 got `4.093173`). `E_0 = sup rho|eta_0|/M = 7.6613`. Both are
   reported.
6. **A plausible sharpening of P2's rate is worth exactly nothing.** The rate `Gamma` in P2
   can be replaced by the sharp `-Lambda_min(sym grad_5 b)`, by the same computation as
   Lemma 4.2. On the reference strain `b = a diag(1,1,1,1,-2)` that quantity is `2a`, and
   `Gamma = ||grad_5 b||_op = 2a` as well. No gain. Recorded so the next seat does not spend
   a day on it. (refute-u2 MAJOR-1's `4.18x` sits in *where the supremum is attained*, not in
   the rate, and is untouched by this note.)

**M2 (H3-V) is restated and priced, and the price is a constraint on `f`, not on `L_*`.**

7. **The repair ASSEMBLY BLOCK 5 promises is not a one-line change.** BLOCK 5 and the brief
   propose replacing (H3) by two clauses, `|Delta_5 eta_0 + eta_0/r^2| <= vartheta M/r^3` and
   `|eta_0 + M sgn(z)/r| <= vartheta M/r`. Step 2 of Theorem V.4 contracts the Hessian
   against `C_tau = diag(2 sigma_y I_4, 2 sigma_z)`, which is **anisotropic**:
   `sigma_z/sigma_y = I_4/I_2 = 3.5130747`. Writing
   `(1/2)C:grad^2 e = sigma_y Delta_5 e + (sigma_z - sigma_y) d_z^2 e` shows the Laplacian
   clause leaves `d_z^2 e` uncontrolled with a coefficient `(sigma_z - sigma_y)/sigma_y =
   2.513`. Step 2's remainder is moreover a **fourth**-derivative statement. The correct
   hypothesis is the scale-free family `k = 0,1,2,3,4` of section B2, and it does contain the
   brief's two clauses (`|Delta_5 eta_0 + eta_0/r^2| <= 11 vartheta M/r^3`).
8. **Theorem V.4' is stated with every hypothesis** (section B2), with the three additive
   costs, and reduces to V.4 verbatim at `vartheta = 0`:
   ```
     |eta(X(tau),tau)/eta_0(x_*) - (1 - s_C)|  <=  E_vartheta + Q ( E_hess + E_4 + E_tail ) ,
     Q := (1+vartheta)/(1-vartheta) ,  vt := vartheta/(1-vartheta) ,
     E_vartheta := vt (9 + 2 sigma_z/sigma_y) s_C  =  16.026 vt s_C .
   ```
9. **`vartheta` for the theorem's datum is `O(1)` at the `f` the budget picks, and the
   hypothesis is therefore vacuous there.** On the ball the angular profile is identically
   `1` (that is what `d` is chosen for), so the entire deviation is the `tanh` radial ramp:
   `eta_0 - eta_P = M psi(rho)/r` with `psi = 1 - Theta` and, at `eps_r = 1/4`,
   `psi_in(rho) = e^2/(e^2 + (rho/rho_0)^8)` exactly. Then `vartheta_0 = psi(rho_* - d)`
   exactly, and
   ```
      f          0.05     0.1      0.25     0.5      1        2        4         8
      rho_min    1.000    1.000    1.000    1.000    1.235    1.852    3.087     5.556
      vartheta   1.3007   1.4593   5.0734   9.6389   6.8754   2.4235   0.055612  5.0751e-04
   ```
   The budget's own minimiser is `f = 0.1`, where `vartheta = 1.4593 > 1` and (H3), (H3')
   are both empty. **`vartheta` for `f = 1` is `6.87194`.** The hypothesis first becomes
   usable at `f = 4`.
10. **A second defect at the same place, independent of (H3-V).** At `f = 0.1` the tracked
    material point sits at `u_* = log 1.1 = 0.0953`, below the ramp's half-height
    `u = eps_r = 0.25`: `Theta(rho_*) = 0.224868`. The tracked point carries `22.5 %` of the
    nominal amplitude `M`. Nothing in `ell_loss(f) = log 2 + log(1+f) + 3 log lam_max +
    log((1+mu)/(1-mu)) + eps_r` forces `Theta(rho_*) ~ 1`, so the `f`-minimisation is free to
    put the tracked point inside the inner ramp, and it does.
11. **The cost of fixing both is small.** `f` is nearly free in the proved column
    (`eps_ell = 1.18e-05` against `eps_a = 0.345` at `L = 3e+05`), so forcing `f >= 4`
    changes `L_*` by `{{LSTAR_CHANGE}}`. Numbers in section B5.

### 0.1 Grade

**FILED.** No Solid, no Major. The content is (a) P1 and P2 written out with no gap and with
two of their hypotheses converted into proved lemmas, (b) closed forms for `E_0` and
`Gfrak_0`, (c) Theorem V.4' with `vartheta` carried through, (d) the measurement that the
budget's own `f` makes (H3-V) vacuous, and the price of not doing that.

**FL-000 stands.**

---

# PART A. (P-max)

## A1. The setting, and every hypothesis

`R^5 = R^4_y x R_z`, `r = |y|`, `rho = |x|`, `phi` the polar angle from the `z` axis,
`eta = omega^theta/r`, `b` the 5-D lift, `D_t = d_s + b.grad_5`. The equation is the one
`u2` section 5 works with:

```
      D_t eta  =  nu Delta_5 eta        on  R^5 x (0, tau] .                       (A.1)
```

**Hypotheses on the drift.** These are `u2`'s (H1) and (H4) plus the regularity that
`(H-K2)` already supplies.

* **(B1)** `b` is continuous on `R^5 x [0,tau]`, `b(.,s) in C^{1,1}(R^5)` for each `s`, with
  `Gamma-bar := sup_{[0,tau]} ||grad_5 b(.,s)||_{L^inf, op} < infinity` and
  `K_2 := sup_{[0,tau]} ||grad_5^2 b(.,s)||_{L^inf} < infinity`.
  *Where it comes from:* `Gamma-bar` is the bootstrap posture (H4), and `K_2 < infinity` is
  `(H-K2)`, **PROVED** in the chain by `s3close/hk2` (`K_2 <= 161.7735 M/rho_0`, no `log`,
  for the `tanh`-radially-mollified datum, which is the theorem's datum). So (B1) costs the
  chain nothing new. This is the only place in the P-max argument where `K_2` is used, and it
  is used only qualitatively.
* **(B2)** `b(0,s) = 0`. **PROVED**, `u2` Lemma 4.1, from `z`-oddness.
* **(B3)** `x.b(x,s) <= Gamma_rad(s)|x|^2` and `|b(x,s)| <= Gamma(s)|x|`, with
  `Gamma_rad(s) = sup_x Lambda_max(sym grad_5 b(x,s))`. **PROVED**, `u2` Lemma 4.2 and
  Corollary 4.3.
* **(B4)** `div_5 b = 2a`. Recorded because it is the thing that does **not** appear; see
  Lemma D.

**Hypotheses on the datum.** `M > 0`, `rho_0 = s sqrt(nu/M)`, `R = rho_0 e^L`,
`u = log(rho/rho_0)`.

* **(D1)** `eta_0 in W^{1,inf}(R^5) cap C(R^5)`, `z`-odd.
* **(D2)** `E_0 M := sup_x |x| |eta_0(x)| < infinity`.
* **(D3)** `Gfrak_0 M := ess sup_x |x|^2 |grad eta_0(x)| < infinity`.
* **(D4)** the tail: `|eta_0(x)| <= Afrak (M/rho_0) Xi_9(rho)`, `Xi_m(rho) := (1+(rho/R)^2)^{-m/2}`.
* **(D5)** the structure `eta_0 = -(M/rho) F(phi) Theta(u)` with `F in W^{1,inf}(0,pi)`, even
  about `phi = 0` and `phi = pi`. Used **only** by the approximation step, Proposition 3.

For the theorem's datum all five hold, with the constants of section A6.

**Hypothesis on the solution.**

* **(S)** `eta` is a bounded solution of (A.1) on `R^5 x (0,tau]`, continuous on
  `R^5 x [0,tau]` with `eta(.,0) = eta_0`. It is unique in that class: Lemma E with
  `gamma = 0` applied to `+-(eta - eta')` gives uniqueness, so (S) is an existence statement
  only. That existence is THEOREM_S3's item M6, `(H*)`, and is not claimed here.

Nothing else is assumed. In particular **no decay hypothesis and no `L^1` hypothesis is
assumed**; refute-u2 MINOR-4's two items are proved below from (D4) and (B1).

## A2. Step 0: the solution is smooth for `s > 0`

> **LEMMA A (interior regularity). PROVED.** Under (B1) and (S), for every
> `alpha in (0,1)`, `eta in C^{3+alpha}_{loc}(R^5 x (0,tau])`. In particular `eta(.,s)` is
> `C^3` in `x` and `C^1` in `s` on compact subsets of `R^5 x (0,tau]`, and so are `grad_5 eta`
> up to order 2.

*Proof.* Write (A.1) as `d_s eta - nu Delta_5 eta + b.grad_5 eta = 0`, a uniformly parabolic
non-divergence equation whose second-order coefficient is the constant matrix `nu I` and whose
drift `b` is, by (B1), Lipschitz in `x` and continuous in `t`, hence `C^{alpha}` on every
compact cylinder. `eta` is bounded. The interior Schauder estimate for bounded solutions of a
uniformly parabolic equation with `C^{alpha}` coefficients (Lieberman, *Second Order Parabolic
Differential Equations*, Theorem 4.9; equivalently Friedman chapter 3, or LSU chapter IV)
gives `eta in C^{2+alpha}_{loc}` on `(0,tau]` with `|eta|_{2+alpha;Q'} <= C |eta|_{0;Q}` for
`Q' \subset\subset Q`.

Bootstrap once. By Lemma D below, `v := d_j eta` solves
`d_s v - nu Delta_5 v + b.grad_5 v = -(d_j b).grad_5 eta` on the same cylinder. The right side
lies in `C^{alpha}_{loc}`, because `d_j b in C^{0,1}` by (B1) and `grad_5 eta in
C^{1+alpha}_{loc}` by the previous paragraph, and `v` is locally bounded for the same reason.
The same interior estimate gives `v in C^{2+alpha}_{loc}`, i.e. `eta in C^{3+alpha}_{loc}`. ∎

*Remark.* The linear growth of `b` at infinity is irrelevant: interior Schauder is a local
statement and `b` is bounded on every compact cylinder. This is why `(B1)`'s global `K_2` is
not needed quantitatively here.

*Consequence.* For every `eps > 0` the functions `u_eps`, `G_eps`, `V_eps`, `W_eps` built in
A3 are genuinely `C^{2,1}` on `R^5 x (0,tau]`: `sqrt(t + eps^2)` is smooth on `t >= 0` and
`R_eps = sqrt(rho^2 + eps^2)` is smooth everywhere including the origin. **No viscosity
solutions, no semi-jets, no distributional maximum principle is needed anywhere in this note.**

## A3. The exact regularised identities

> **LEMMA B (weights). EXACT, symbolic.** In `R^n`, for a `C^2` function `g` and every `k`,
> ```
>     rho^k Delta g  =  Delta(rho^k g) - (2k/rho^2) x.grad(rho^k g) + k(k-n+2)(rho^k g)/rho^2 .
> ```
> The zeroth-order coefficient `k(k-n+2)` is `<= 0` exactly for `0 <= k <= n-2`; at `n = 5`
> that is `k <= 3`, and `k(k-3) = -2` for both `k = 1` and `k = 2`.
> With `R_eps := sqrt(rho^2 + eps^2)`, `V := R_eps g`, `W := R_eps^2 g`, in `n = 5`,
> ```
>     R_eps   Delta g = Delta V - (2/R_eps^2) x.grad V - (2 rho^2 +  5 eps^2) V / R_eps^4 ,  (A.2)
>     R_eps^2 Delta g = Delta W - (4/R_eps^2) x.grad W - (2 rho^2 + 10 eps^2) W / R_eps^4 .  (A.3)
> ```
> Both zeroth-order coefficients are `<= 0` for every `rho` and every `eps >= 0`, and both
> reduce to the `eps = 0` rows at `eps = 0`.

*Proof.* Cartesian expansion. `q1` verifies the general-`n` identity as a sympy residual `0`
for a generic `sympy.Function` at `n = 2,...,7` and `k = 0,...,4` (30 checks), and (A.2),
(A.3) and their `eps -> 0` limits as four further residuals `0`, together with
`|grad R_eps|^2 = rho^2/R_eps^2`, `Delta R_eps = (5 eps^2 + 4 rho^2)/R_eps^3` and
`Delta R_eps^2 = 10`. ∎

This is the point of the `eps`-regularisation of the **weight**, and it is worth saying
plainly: mollifying `rho` costs nothing. It does not introduce an error term, it does not
spoil the sign, and it makes the first-order coefficient `2k x/R_eps^2` a **bounded** vector
field (`|2k x/R_eps^2| <= k/eps`), which is what removes the singularity at the origin that
`u2` section 5 leaves implicit.

> **LEMMA C (regularised Kato). EXACT, symbolic.** Let `eta` be `C^2` and `g := grad_5 eta`.
> For `eps > 0` put `u_eps := sqrt(eta^2 + eps^2) - eps` and `G_eps := sqrt(|g|^2+eps^2) - eps`.
> Then `0 <= u_eps <= |eta|`, `0 <= G_eps <= |g|`, both increase to their limits as `eps` falls,
> and
> ```
>     D_t u_eps  <=  nu Delta_5 u_eps ,                                              (A.4)
>     D_t G_eps  <=  Gamma(s) (G_eps + eps) + nu Delta_5 G_eps .                      (A.5)
> ```

*Proof.* Write `u~ := sqrt(eta^2+eps^2)`. Then
`Delta u~ - eta Delta eta/u~ = eps^2 |grad eta|^2 / u~^3 >= 0` (residual `0`, `q1` C), and
`D_t u~ = eta D_t eta/u~ = nu eta Delta eta / u~ <= nu Delta u~`. Subtracting the constant
`eps` changes neither `D_t` nor `Delta`, which gives (A.4) with **no** error term.
For the vector, `G~ := sqrt(|g|^2+eps^2)` has
`Delta G~ - g.Delta g/G~ = |grad g|^2/G~ - |g.grad g|^2/G~^3` (residual `0`, `q1` C), and this
is `>= 0` because `sum_j |g . d_j g|^2 <= |g|^2 |grad g|^2 <= G~^2 |grad g|^2` by
Cauchy-Schwarz. By Lemma D, `g . D_t g = -g^T (grad b) g + nu g.Delta g <= Gamma |g|^2 +
nu g.Delta g`, so `D_t G~ <= Gamma |g|^2/G~ + nu Delta G~ <= Gamma G~ + nu Delta G~`, which is
(A.5) after subtracting `eps`. ∎

The single `Gamma eps` on the right of (A.5) is the whole price of the shift, and it is the
only inexactness anywhere in A3. It is removed in the limit in A5.

> **LEMMA D (the gradient equation, and the absence of a zeroth-order term). EXACT, symbolic.**
> From (A.1),
> ```
>     d_s(d_j eta) + b.grad_5(d_j eta) + (d_j b).grad_5 eta  =  nu Delta_5 (d_j eta) .   (A.6)
> ```
> `div_5 b = 2a` appears nowhere in (A.6). It appears only when (A.1) is written in
> conservation form, `d_s eta + div_5(b eta) = nu Delta_5 eta + 2 a eta`, which the argument
> never uses.

*Proof.* Differentiate (A.1) in `x_j`. `q1` D verifies (A.6) for `j = 1,...,5` as five
residuals `0` for generic `sympy.Function`s `eta` and `b_1,...,b_5`, and verifies the
advective/conservative identity as a sixth. ∎

This is the same answer `u2` Lemma 5.0 and `gap-V-aronson` section 2.3 give, obtained here
independently.

## A4. The two differential inequalities

> **PROPOSITION 1 (P1's inequality). PROVED.** Under (B1) to (B3), (S), for every `eps > 0`,
> with `V_eps := R_eps u_eps`,
> ```
>   D_t V_eps  <=  Gamma_rad(s) V_eps
>                  + nu [ Delta V_eps - (2/R_eps^2) x.grad V_eps
>                                     - (2 rho^2 + 5 eps^2) V_eps / R_eps^4 ]        (A.7)
> ```
> pointwise on `R^5 x (0,tau]`.

*Proof.* `D_t R_eps = (x.b)/R_eps <= Gamma_rad rho^2/R_eps <= Gamma_rad R_eps` by (B3) and
`rho <= R_eps`. Hence `D_t V_eps = u_eps D_t R_eps + R_eps D_t u_eps <= Gamma_rad V_eps +
nu R_eps Delta u_eps` by (A.4), and (A.2) rewrites `R_eps Delta u_eps`. The zeroth-order term
is `<= 0` because `V_eps >= 0`. ∎

> **PROPOSITION 2 (P2's inequality). PROVED.** Under the same hypotheses, with
> `W_eps := R_eps^2 G_eps`,
> ```
>   D_t W_eps  <=  (Gamma(s) + 2 Gamma_rad(s)) W_eps  +  Gamma(s) eps R_eps^2
>                  + nu [ Delta W_eps - (4/R_eps^2) x.grad W_eps
>                                     - (2 rho^2 + 10 eps^2) W_eps / R_eps^4 ] .      (A.8)
> ```

*Proof.* `D_t(R_eps^2) = 2 x.b <= 2 Gamma_rad rho^2 <= 2 Gamma_rad R_eps^2`. Then
`D_t W_eps = G_eps D_t(R_eps^2) + R_eps^2 D_t G_eps <= 2 Gamma_rad W_eps + R_eps^2[Gamma(G_eps
+ eps) + nu Delta G_eps]` by (A.5), and (A.3) rewrites `R_eps^2 Delta G_eps`. ∎

Both are exact identities plus one signed inequality each, and both have the *good* sign on
the zeroth-order weight term, which is the `n = 5` fact the brief asks to be checked.

## A5. The comparison lemma: two routes, neither of which has a gap

> **LEMMA E (comparison on `R^5`). PROVED.** Let `gamma, sigma in L^inf(0,tau)`, `sigma >= 0`,
> `k in {0,1,2}`, `eps in (0,1]`, and let `Z in C(R^5 x [0,tau]) cap C^{2,1}(R^5 x (0,tau])`
> satisfy, **at every point where `Z > 0`**,
> ```
>   D_t Z  <=  gamma(s) Z + sigma(s)
>              + nu [ Delta Z - (2k/R_eps^2) x.grad Z - q(x) Z ] ,      q >= 0 .      (A.9)
> ```
> Suppose **either**
> * **(i)** `Z(x,s) -> 0` as `rho -> infinity`, uniformly for `s` in compact subsets of
>   `(0,tau]`;  **or**
> * **(ii)** `Z(x,s) <= C (1 + rho^{beta - theta})` on `R^5 x [0,tau]` for some
>   `2 <= beta <= 3` and `theta > 0`.
>
> Then for `0 <= s_0 <= s <= tau`,
> ```
>   sup_x Z(x,s)^+  <=  e^{int_{s_0}^s gamma} [ sup_x Z(x,s_0)^+ + int_{s_0}^s sigma ] .
> ```

*Proof.* Write `m(s) := sup_x Z(x,s)^+`.

*Case (i).* `Z^+` is continuous, vanishes at infinity uniformly on `[s_0, tau]`, so the
supremum is attained and all maximisers lie in a fixed compact set. Where `m(s) > 0`, at any
maximiser `x_s` we have `grad Z = 0` and `Delta Z <= 0`, hence
`d_s Z(x_s,s) = D_t Z(x_s,s) - b.grad Z(x_s,s) = D_t Z(x_s,s) <= gamma m(s) + sigma(s)`, using
`q Z >= 0`. The upper right Dini derivative of `m` satisfies `D^+ m(s) <= d_s Z(x_s,s)`
(Hamilton's trick, valid for the supremum of a family of `C^1` functions attained on a compact
set), and where `m(s) = 0` the Dini derivative is `<= sigma(s)`. Gronwall for the Dini
derivative gives the claim.

*Case (ii).* Put `h(x,s) := e^{Ks} (rho^2+1)^{beta/2}` with
```
        K  :=  beta Gamma-bar  +  nu beta(beta+3)  +  ||gamma||_inf  +  1 .
```
Then `h > 0` and, pointwise,
`D_t h - nu[Delta h - (2k/R_eps^2) x.grad h - q h] - gamma h >= h > 0`. Indeed
`D_t h = K h + e^{Ks} phi'(rho)(x.b)/rho >= K h - Gamma-bar rho phi' e^{Ks} >= (K -
beta Gamma-bar) h` using `rho phi' <= beta phi` (`q1` E, symbolic); `- nu Delta h =
-nu e^{Ks}(phi'' + 4 phi'/rho) >= -nu beta(beta+3) h` using
`phi'' + 4 phi'/rho = 5 beta (rho^2+1)^{beta/2-1} + beta(beta-2) rho^2 (rho^2+1)^{beta/2-2}
<= beta(beta+3)(rho^2+1)^{beta/2-1}` for `beta >= 2` (`q1` E, symbolic identity plus a grid
check whose minimum gap is `0`, attained at `beta = 2`); and the two remaining terms,
`+nu(2k/R_eps^2) x.grad h >= 0` and `+ nu q h >= 0`, are dropped.
For `delta > 0` set `Z_delta := Z - delta h`. By (ii) and `beta - theta < beta`,
`Z_delta(x,s) -> -infinity` as `rho -> infinity`, uniformly in `s in [0,tau]`, so
`m_delta(s) := sup_x Z_delta(x,s)^+` is attained on a compact set whenever it is positive, and
where `Z_delta > 0` we have `Z > 0` so (A.9) applies; subtracting the barrier inequality,
`Z_delta` satisfies (A.9) with the same `gamma, sigma` and a strictly negative slack. The
Case (i) argument now runs verbatim on `Z_delta`, giving
`m_delta(s) <= e^{int gamma}[m_delta(s_0) + int sigma] <= e^{int gamma}[m(s_0) + int sigma]`.
Since `Z(x,s) <= Z_delta(x,s) + delta h(x,s)` and `h` is finite, letting `delta -> 0`
pointwise gives the claim. ∎

**Case (ii) is the answer to the referee's first attack point.** It needs no decay of the
solution at all: only that `Z` grows more slowly than `rho^beta` for some `beta <= 3`, and the
`beta <= 3` ceiling is exactly Lemma B's `k(k-n+2) <= 0` ceiling at `n = 5`. For P1,
`V_eps <= (rho + eps)|eta| <= (rho+1)||eta_0||_inf`, so `beta = 2, theta = 1`. For P2,
`W_eps <= (rho^2+1)(|grad eta| + 1)`, and `sup_x |grad eta(.,s)| < infinity` for each
`s > 0` by Corollary F' below, so `beta = 3, theta = 1`. The supremum cannot escape, because
the barrier catches it before it gets anywhere.

> **LEMMA F (decay, with explicit constants). PROVED.** Assume (B1) to (B3), (D4), (S). Then
> for every `s in [0,tau]` and every `x`,
> ```
>   |eta(x,s)|  <=  Afrak (M/rho_0) e^{Lambda_9 s} Xi_9(rho) ,
>   Lambda_9 := 9 Gamma-bar + 99 nu/R^2 ,      Xi_9(rho) = (1 + (rho/R)^2)^{-9/2} .
> ```
> In particular `e^{Lambda_9 tau} = e^{9 c_G + 99 nu tau/R^2}` and, since
> `nu tau/R^2 = c e^{-2L}/(s^2 L)`, the second term is below `1e-300` at every `L` the theorem
> uses, so the constant is `e^{9 c_G}` for all practical purposes.

*Proof.* Let `psi(x,s) := Afrak (M/rho_0) e^{Lambda_9 s} R^9 (rho^2+R^2)^{-9/2}`. With
`phi := (rho^2+R^2)^{-9/2}`: `rho|phi'| = 9 rho^2 (rho^2+R^2)^{-11/2} <= 9 phi`, so
`b.grad psi >= -Gamma-bar rho |phi'| (...) >= -9 Gamma-bar psi` by (B3); and
`Delta phi = -45 (rho^2+R^2)^{-11/2} + 99 rho^2 (rho^2+R^2)^{-13/2} <= 99 phi/(rho^2+R^2)
<= 99 phi/R^2`, so `-nu Delta psi >= -99 nu psi/R^2`. Hence
`d_s psi + b.grad psi - nu Delta psi >= (Lambda_9 - 9 Gamma-bar - 99 nu/R^2) psi = 0`, i.e.
`psi` is a supersolution. Now put `Z_eps := u_eps - psi` with `u_eps` of Lemma C. It is continuous on
`R^5 x [0,tau]` and `C^{2,1}` on `R^5 x (0,tau]` (Lemma A), satisfies
`D_t Z_eps <= nu Delta Z_eps` by (A.4) and the supersolution property, has
`Z_eps(.,0) <= |eta_0| - psi(.,0) <= 0` by (D4), and is bounded above by `||eta_0||_inf`.
Lemma E case (ii) with `k = 0` (no weight term), `beta = 2`, `theta = 2`,
`gamma = sigma = 0` gives `sup Z_eps(.,s)^+ <= sup Z_eps(.,0)^+ = 0`, i.e. `u_eps <= psi`;
let `eps -> 0`. ∎

> **COROLLARY F' (gradient decay, global gradient bound, and `L^1`). PROVED.** Under the same
> hypotheses there is `C_* = C_*(nu, Gamma-bar, K_2, rho_0)` such that for
> `s in [s_0, tau]`, `s_0 > 0`, and `rho >= 2 rho_0`,
> ```
>   |grad eta(x,s)|  <=  2^9 C_*(s_0) Afrak (M/rho_0^2) e^{Lambda_9 tau} Xi_9(rho) .
> ```
> Consequently, for each `s in (0,tau]`,
> (a) `sup_x |grad eta(.,s)| < infinity`;
> (b) `rho^2|grad eta(x,s)| -> 0` and `rho|eta(x,s)| -> 0` as `rho -> infinity`, uniformly on
>     `[s_0,tau]`;
> (c) `grad eta(.,s) in L^1(R^5)`, with
>     `||grad eta(.,s)||_{L^1} <= |S^4|[ Gfrak_0 M (2 rho_0)^3/3 + C (M/rho_0^2) int_{2rho_0}^infinity Xi_9 rho^4 drho ] < infinity`.

*Proof.* Apply the interior gradient estimate of Lemma A on the cylinder
`Q_x := B(x,rho_0) x (s - min(s_0, rho_0^2/nu)/2, s]`, on which `b` is bounded (by
`Gamma-bar(rho+rho_0)`) and Lipschitz: `|grad eta(x,s)| <= C_*(s_0) rho_0^{-1}
||eta||_{L^inf(Q_x)}`. For `rho >= 2 rho_0` every point of `Q_x` has radius `rho' >= rho/2`,
and `Xi_9(rho') <= Xi_9(rho/2) <= 2^9 Xi_9(rho)`; Lemma F bounds `||eta||_{L^inf(Q_x)}`.
(a) follows since `grad eta` is continuous and bounded on `{rho <= 2 rho_0} x [s_0,tau]` by
Lemma A. (b) follows from `rho^2 Xi_9(rho) -> 0` and `rho Xi_9(rho) -> 0`.
(c) split at `rho = 2 rho_0`, use `|grad eta| <= Gfrak_0 M/rho^2` (Theorem P2' below) inside
and the display outside; `int Xi_9 rho^4 drho < infinity` because `9 > 5`. ∎

*Where the `L^1` is used (the brief's item (d)).* Exactly once in the whole chain: `u2`
section 6, first display, *"Since `grad eta in L^1 cap L^inf`, `grad a = K * grad_5 eta`
(differentiation under the convolution)"*. It is used nowhere in P1, nowhere in P2, nowhere in
section 7's far/near split, and nowhere in section 8's fixed point. refute-u2 MINOR-4(b) is
right that P2 alone gives only `|grad eta| <~ rho^{-2}`, which is not `L^1` in `R^5`; the
tails are what supply it, and Corollary F'(c) supplies them with a constant. The alternative
route, `u2` section 5.4(c)'s `TV(eta(s)) <= e^{3 c_G} TV(eta_0)`, also works and the refuter
re-derived it; it needs an integration by parts on `R^5` that itself wants the decay, so
F'(c) is the shorter road. The datum's own value is
`TV(eta_0) = 31.138 M R^3` (`q2`, stable in `L` to `4e-06` between `L = 20` and `L = 40`).
The `L^1` norm is therefore enormous, of order `M R^3`; only its finiteness is used, and no
constant in the budget depends on it.

## A6. The two theorems

> **THEOREM P1'. PROVED.** Assume (B1) to (B3), (D1), (D2), (S). Then for every
> `s in [0,tau]`
> ```
>     sup_{x in R^5} |x| |eta(x,s)|  <=  e^{c_R(s)} E_0 M ,      c_R(s) := int_0^s Gamma_rad .
> ```
> No decay hypothesis and no `L^1` hypothesis is used.

*Proof.* Fix `eps in (0,1]`. By Lemma A and the consequence noted there,
`V_eps in C(R^5 x [0,tau]) cap C^{2,1}(R^5 x (0,tau])` and `V_eps >= 0`. Proposition 1 is
(A.9) with `k = 1`, `gamma = Gamma_rad`, `sigma = 0`,
`q = (2 rho^2 + 5 eps^2)/R_eps^4 >= 0`. Growth: `V_eps = R_eps u_eps <= (rho + eps)|eta| <=
(rho+1) ||eta_0||_inf`, so hypothesis (ii) of Lemma E holds with `beta = 2`, `theta = 1`.
Hence
```
   sup_x V_eps(.,s)  <=  e^{c_R(s)} sup_x V_eps(.,0)  <=  e^{c_R(s)} [ E_0 M + eps ||eta_0||_inf ] ,
```
using `sqrt(rho^2+eps^2) <= rho + eps` and `u_eps <= |eta_0|` at `s = 0`. Let `eps -> 0`:
`R_eps -> rho` and `u_eps -> |eta|` pointwise, so `rho|eta(x,s)| <= e^{c_R(s)} E_0 M` for
every `x`. ∎

> **THEOREM P2'. PROVED.** Assume (B1) to (B4) [(B4) is used only through Lemma D], (D1) to
> (D5), (S). Then for every `s in (0,tau]`
> ```
>     sup_{x in R^5} |x|^2 |grad eta(x,s)|  <=  e^{c_G(s) + 2 c_R(s)} Gfrak_0 M
>                                            =  e^{p c_G} Gfrak_0 M ,   p := 1 + 2 c_R/c_G <= 3 .
> ```

*Proof, for a datum with `grad eta_0 in C^{alpha}`.* Fix `eps in (0,1]` and `s_0 in (0,s)`.
`W_eps in C^{2,1}(R^5 x (0,tau])`, `W_eps >= 0`. By Corollary F'(b), `W_eps <= (rho^2+1)
|grad eta| -> 0` as `rho -> infinity` uniformly on `[s_0, tau]`, with a bound **independent of
`eps`**; hence there is `D < infinity`, independent of `eps in (0,1]` and of `s in [s_0,tau]`,
such that every maximiser of `W_eps(.,s)` lies in `{rho <= D}`. Proposition 2 is then (A.9)
with `k = 2`, `gamma = Gamma + 2 Gamma_rad`, `q = (2 rho^2+10 eps^2)/R_eps^4 >= 0`, and
`sigma(s) = Gamma(s) eps (D^2 + 1)`, since only points with `rho <= D` matter for the Dini
step. Lemma E case (i) gives
```
   sup_x W_eps(.,s)  <=  e^{int_{s_0}^s (Gamma + 2 Gamma_rad)} [ sup_x W_eps(.,s_0)
                                                                 + Gamma-bar eps (D^2+1) tau ] .
```
Let `eps -> 0`, then `s_0 -> 0`. The second term vanishes with `eps`. For the first,
`sup_x rho^2 |grad eta(.,s_0)| -> Gfrak_0 M` as `s_0 -> 0`: locally uniform convergence
`grad eta(.,s_0) -> grad eta_0` follows from the interior Schauder estimate applied on
cylinders reaching `t = 0`, which is legitimate because `grad eta_0 in C^{alpha}`, and
uniform smallness outside a compact set follows from Lemma F together with the same interior
estimate applied at `s_0`. Hence
`sup rho^2|grad eta(.,s)| <= e^{int_0^s(Gamma + 2 Gamma_rad)} Gfrak_0 M`, which is the claim.
Case (ii) of Lemma E gives the same conclusion without Corollary F'(b), at the cost of
carrying the `eps`-source through the barrier; both routes were written and they agree. ∎

> **PROPOSITION 3 (the datum is only Lipschitz in the angle; the approximation). PROVED.**
> Let `eta_0` satisfy (D5) and let `chi_sigma` be a smooth even mollifier at scale `sigma`.
> Put `F_sigma := F * chi_sigma` (after even reflection at `phi = 0` and `phi = pi`) and
> `eta_0^{(sigma)} := -(M/rho) F_sigma(phi) Theta(u)`. Then
> * (a) `eta_0^{(sigma)}` is `z`-odd and `C^infinity` away from the origin, and `C^{1,alpha}`
>   on `R^5` (near the origin `Theta = O(rho^8)` so `eta_0^{(sigma)} = O(rho^7)`; near the
>   axis `F_sigma` is even about `phi = 0`, hence a smooth function of `r^2`);
> * (b) `sup rho|eta_0^{(sigma)}| <= E_0 M` and `ess sup rho^2|grad eta_0^{(sigma)}| <=
>   Gfrak_0 M`, **for every `sigma > 0`**;
> * (c) `eta_0^{(sigma)} -> eta_0` uniformly.
> Let `eta^{(sigma)}` solve the **same linear equation, with the same `b`**, from
> `eta_0^{(sigma)}`. Then `||eta^{(sigma)}(.,s) - eta(.,s)||_inf <= ||eta_0^{(sigma)} -
> eta_0||_inf -> 0` (Lemma E, `gamma = sigma = 0`, applied to `+-` the difference), and by
> Lemma A `grad eta^{(sigma)}(.,s) -> grad eta(.,s)` locally uniformly for each `s > 0`.
> Therefore P2' holds for `eta` itself, for every `s in (0,tau]`.

*Proof of (b).* `||F_sigma||_inf <= ||F||_inf` and `F_sigma' = F' * chi_sigma` so
`||F_sigma'||_inf <= ||F'||_inf`. With `A(u) := Theta - dTheta/du` and `B(u) := Theta`,
`rho|eta_0| = M|F| B` and `rho^2|grad eta_0| = M sqrt(F^2 A^2 + F'^2 B^2)`, so both weighted
suprema are non-decreasing functions of the pair `(||F||_inf, ||F'||_inf)` at fixed `u`, and
for this datum both suprema are attained at the **same** angle `phi = pi - delta`, where
`|F| = 1/sin delta` and `|F'| = cos delta/sin^2 delta`. `q2` confirms it numerically at four
mollification scales:
```
   sigma      E_0(sigma)/E_0    Gfrak_0(sigma)/Gfrak_0     ||F_sigma'||_inf
   0.02       0.99931989        0.79723209                 46.283387
   0.01       0.99964616        0.88422792                 51.394955
   0.005      0.99981964        0.93725313                 54.509148
   0.002      0.99992702        0.97329425                 56.625546
```
every ratio below `1`, both rising to `1` as `sigma -> 0`. ∎

*Remark on what the approximation is and is not.* `b` is **not** regenerated from
`eta^{(sigma)}`; `eta^{(sigma)}` is not a Navier-Stokes solution. P1 and P2 are statements
about the **linear** equation (A.1) with a given drift, which is exactly the posture of `u2`'s
(H1) and (H4), so the approximation is legitimate and changes nothing about `b`, `Gamma`,
`Gamma_rad` or `c_G`. This is also why THEOREM_S3's item M5 (H5-Lip) is, for P1 and P2,
discharged: the datum's kinks are removed at no cost in either constant, and the classical
regularity that M5 asks for is Lemma A.

## A7. The two constants, for the theorem's datum

The datum is THEOREM_S3 section 1.1: ASSEMBLY section 1.1's angular profile with `hk2`'s
`tanh` radial ramp. In the form (D5),

```
   F(phi) = sgn(cos phi) g(phi)/sin phi ,
   g(phi) = min(1, phi_ax/delta) min(1, |phi - pi/2|/delta_m) ,  delta = 7.5 deg, delta_m = 5 deg,
   Theta(u) = (1/2)[ tanh((u-eps_r)/eps_r) - tanh((u-L+eps_r)/eps_r) ] ,  eps_r = 0.25 .
```

> **CLOSED FORMS.** In the limit `L -> infinity` (at finite `L` both carry a factor
> `sup Theta = 1 - O(e^{-2(L-2eps_r)/eps_r})`),
> ```
>    E_0     = 1/sin delta    = 7.661297575540389      attained at phi = delta and pi - delta
>    Gfrak_0 = 1/sin^2 delta  = 58.69548054098104      attained at phi = pi - delta
>            = E_0^2 .
> ```

*Derivation.* On the plateau arc `delta < phi < pi/2 - delta_m` and its mirror, `g == 1`, so
`F = sgn(cos phi)/sin phi`, `F' = -sgn(cos phi) cos phi/sin^2 phi`, and
```
     F^2 + F'^2  =  (sin^2 phi + cos^2 phi)/sin^4 phi  =  1/sin^4 phi ,
```
maximised on that arc at its closest approach to the axis, `phi = delta`. In the taper arc
`phi < delta`, `|F| = phi/(delta sin phi) <= 1/sin delta` and `|F'| = |sin phi - phi cos phi|/
(delta sin^2 phi) <= 0.3345`; in the equatorial arc `|F'| <= 1/delta_m = 11.459`. Deep in the
plateau `Theta = 1`, `dTheta/du = 0`, so `A = B = 1`, and the joint supremum over `u` is
attained there (`q2` computes the joint supremum exactly, by reducing the `u`-scan to the
`21874` vertices of the upper convex hull of `{(A(u)^2, B(u)^2)}`, which is legitimate because
the objective is linear in that pair at fixed `phi`). ∎

**Confirmation and comparison.**

| quantity | closed form | this seat, scan | rel. gap | record |
|---|---|---|---|---|
| `E_0` | `7.661297575540389` | `7.661297575537042` | `4.4e-13` | `7.6612975721936944` (`THEOREM_S3` `t1`), rel `4.4e-10` |
| `Gfrak_0` | `58.69548054098104` | `58.695480532064344` | `1.5e-10` | `58.695476310955364` (`t1`), rel `7.2e-08`; `58.6948961` (refute-u2), rel `1.0e-05` |
| `ess sup|F'|` | `cos delta/sin^2 delta = 58.19333256822213` | `58.19333255930511` | `1.5e-10` | `58.193328` (`t1`) |
| `sup_u e^{-u}Theta` | | `0.534314767358705` at `rho = 1.63763 rho_0` | | `0.53431477` (refute-u2), rel `4.9e-09` |
| `||eta_0||_inf rho_0/M` | `E_0 sup_u e^{-u}Theta` | **`4.093544431740673`** | | `4.093173` for (D-C) (refute-u2 MINOR-2), rel `9.1e-05` |
| `||grad eta_0||_inf rho_0^2/M` | | `31.093556949918312` | | not in the record |
| `Afrak` (tail, `m = 9`) | | `4.0935444` (`L = 10, 20, 40`, stable to `4e-08`) | | not in the record |
| `TV(eta_0)/(M R^3)` | | `31.138` (`L = 20, 40`, stable to `4e-06`) | | not in the record |

The two end caps `phi < 1e-03` and `phi > pi - 1e-03` are excluded from the scan and handled
analytically instead, where `|F| <= 7.639438541650668` and `|F'| <= 2.549e-03`, both far below
the suprema; a floating-point evaluation there suffers catastrophic cancellation
(`sin t - t cos t ~ t^3/3` against `sin^2 t ~ t^2`) and that, not the mathematics, is the only
reason the record's grid values sit `1e-05` low. The angular derivative is taken from the
exact per-arc closed form, not by finite differences; the two agree to `1.4e-08` away from the
kinks and the end caps.

**`Gfrak_0 = E_0^2 = 1/sin^2 delta`, and `E_0 = s`,** the same `s` that fixes
`rho_0 = s sqrt(nu/M)`. THEOREM_S3 finding 1 says the binding datum quantity is `Gfrak_0`;
this says `Gfrak_0` is `1/sin^2 delta` and nothing else, so the entire datum dependence of
`Ghat`, `C''` and `L_*` in the proved column is a function of the taper angle alone, through
`sin^{-2} delta`. Widening `delta` from `7.5 deg` to `15 deg` would divide `Gfrak_0` by
`3.9556` at once. What that costs in `kappa_delta` is not computed here.

## A8. What is still a hypothesis after this note

For P1 and P2:

| | status |
|---|---|
| `b in C^{1,1}` with `Gamma-bar`, `K_2` finite | (H4) plus `(H-K2)`, both already in the chain; `K_2` used only qualitatively, in Lemma A |
| `b(0,s) = 0`, `x.b <= Gamma_rad rho^2` | PROVED, `u2` Lemmas 4.1, 4.2 |
| `eta` is **the** bounded solution with datum `eta_0` | existence is THEOREM_S3's M6 `(H*)`; uniqueness is proved here (Lemma E, `gamma = 0`) |
| the contradiction hypothesis `||omega|| <= lambda M` | (S3)'s own, not an extra assumption |
| `rho|eta| -> 0`, `rho^2|grad eta| -> 0` at infinity | **no longer a hypothesis**: not needed at all (Lemma E case (ii)), and proved anyway (Lemma F, Corollary F') |
| `grad eta in L^1` | **no longer a hypothesis**: proved (Corollary F'(c)); used only in `u2` section 6 |
| `eta_0 in C^1` in the angle | **no longer a hypothesis**: Proposition 3, at no cost in `E_0` or `Gfrak_0` |

So the answer the brief asks for: **P1 and P2 are unconditional theorems about the linear
problem (A.1), given (B1) to (B3) and (D1) to (D5). The remaining hypotheses are the ones
(S3) already carries: the bootstrap posture that supplies `Gamma-bar`, `(H-K2)` for `K_2`,
existence in the bounded class (M6), and the contradiction hypothesis. Nothing new is
assumed.**

Two honest caveats, stated so they are not mistaken for coverage.

* Lemma A cites interior parabolic Schauder rather than proving it. That is a textbook
  theorem with a name and a reference, not a gap in the chain, but it is a citation.
* P2' at `s = 0` needs `grad eta_0 in C^{alpha}`, which the theorem's datum does not have;
  Proposition 3 supplies it by mollifying the angle, and the two weighted norms provably do
  not increase. The conclusion of P2' is therefore stated for `s in (0,tau]`, not `[0,tau]`.
  Since `c_G(0) = c_R(0) = 0`, nothing is lost.

## A9. One negative result, recorded so it is not re-attempted

`u2`'s P2 uses `Gamma = ||grad_5 b||_op` as the exponential rate for `|grad eta|`. The sharp
rate is not `Gamma` but `-Lambda_min(sym grad_5 b) = max(-a, (a + sqrt((3a+2Q)^2 +
(2P-omega)^2))/2)`, by exactly the eigenvalue computation of `u2` Lemma 4.2 (the
`(e_r, e_z)` block has trace `-a` and discriminant `(3a+2Q)^2 + (2P-omega)^2`; the transverse
`R^3` carries `a` with multiplicity three). This looks like free headroom and it is not. On
the reference strain `b = a diag(1,1,1,1,-2)`, which is the `O(M log)` part of `u`, the
spectrum of `sym grad_5 b` is `{a (x4), -2a}`, so `-Lambda_min = 2a` while
`||grad_5 b||_op = 2a` as well: the two coincide exactly, and `Gamma_rad = a`, reproducing
`u5`'s `Gamma_rad/Gamma = 1/2`. The sharpening buys **nothing** in the regime that dominates.
refute-u2 MAJOR-1's `4.18x` of headroom lives in *where the weighted supremum is attained*
(on the axis for `E_0`, in the taper transition for `Gfrak_0`), not in the rate, and this note
does not touch it.

---

# PART B. (H3-V)

## B1. What is wrong with (H3), precisely

Theorem V.4 (V-b section 5) assumes

> **(H3)** there are `x_0` with `z_0 > 0`, `r_0 = r(x_0)`, `0 < d < r_0` such that on
> `B(x_0,d)` the datum coincides with the bare plateau `eta_P(x) := -M sgn(z)/r`.

The V-b refuter's Finding 4b and ASSEMBLY section 2.6 BLOCK 5 both record that no member of
the campaign's family satisfies this at any `d > 0`. For the theorem's datum the reason is
sharp and computable: the two angular mollifiers are `min(1, .)`, so on a ball that misses
both kinks the angular factor is **exactly** `1` and contributes nothing; but the radial ramp
is `tanh` and is never identically `1`. So on `B(x_*,d)`,

```
     eta_0  =  eta_P Theta(u) ,      e := eta_0 - eta_P = -eta_P psi(rho) = M psi(rho)/r ,
     psi := 1 - Theta = psi_in + psi_out ,
     psi_in(rho)  = 1/(1 + e^{2(u-eps_r)/eps_r}) = e^2/(e^2 + (rho/rho_0)^8)   at eps_r = 1/4,
     psi_out(rho) = 1/(1 + e^{-2(u-L+eps_r)/eps_r}) <= e^{-8(L - eps_r - u)} ,
```

`psi_out` underflows at every `L` the theorem uses. **The whole of the (H3) violation is the
inner `tanh` ramp, and it is an explicit rational function.** This is exactly the conflict
ASSEMBLY BLOCK 5 and BLOCK 6 leave unresolved: BLOCK 5 chose `min(1,.)` angular profiles
*so that* (H3) would hold, and BLOCK 6 then forced `tanh` back into the radial direction
because a sharp radial edge makes `K_2` infinite.

## B2. THEOREM V.4', with every hypothesis

> **HYPOTHESIS (H3'_vartheta).** There are `x_*` with `z_* > 0`, `r_* := r(x_*)`, and
> `0 < d < r_*`, and a number `vartheta in [0,1)`, such that on `B(x_*,d)`
> ```
>      sup_{|v| = 1} | d_v^k ( eta_0 - eta_P )(x) |  <=  vartheta * k! * M / r(x)^{k+1}
>      for k = 0, 1, 2, 3, 4 and every x in B(x_*, d) .
> ```
> Put `R_- := r_* - d`. At `vartheta = 0` this is (H3).

The normalisation is the exact value of `sup_{|v|=1}|d_v^k eta_P|` (V-b identity (2.4),
`k!/r^{k+1}` from the Legendre generating function), so `vartheta` is a pure relative
deviation and is scale free.

**Why five clauses and not the two BLOCK 5 proposes.** BLOCK 5 and the brief propose
`|Delta_5 eta_0 + eta_0/r^2| <= vartheta M/r^3` together with
`|eta_0 + M sgn(z)/r| <= vartheta M/r`. Two things go wrong.

* Step 2 does not contract the Laplacian. It contracts the Hessian against
  `C_tau = diag(2 sigma_y I_4, 2 sigma_z)`, and
  `(1/2) C_tau : grad^2 e = sigma_y Delta_5 e + (sigma_z - sigma_y) d_z^2 e`. The Laplacian
  clause controls the first term and leaves the second free. The coefficient is not small:
  by Corollary V.5, `sigma_z/sigma_y = I_4/I_2 = 1.7905793948266/0.5096901045590 =
  3.5130747`, so `(sigma_z - sigma_y)/sigma_y = 2.5130747`. The exact cancellation
  `d_z^2 eta_P = 0` that makes `sigma_z` drop out of V.4 identically (Lemma V.3') is a
  property of `eta_P`, and it is precisely the property the perturbation `e` need not have.
* Step 2's remainder is a fourth-derivative statement,
  `|eta_P(x_*+zeta) - P_3(zeta)| <= (1/4!)sup|d_v^4 eta_P| |zeta|^4`, and Step 1 uses
  `|grad eta_0| <= M/R_-^2`. Neither is implied by a value bound and a Laplacian bound.

(H3'_vartheta) does contain the brief's two clauses: `k = 0` is the second verbatim, and from
`k = 0` and `k = 2`, using `Delta_5 eta_P = -eta_P/r^2`,
```
     |Delta_5 eta_0 + eta_0/r^2| = |Delta_5 e + e/r^2| <= 5 (2 vartheta M/r^3) + vartheta M/r^3
                                 = 11 vartheta M/r^3 .
```

> **THEOREM V.4'. PROVED (given V.4's own Steps 0 and 1 machinery).** Assume V-b's (H1),
> (H2), (H4) and (H3'_vartheta) with `vartheta < 1`. Write, as in V.4,
> `s_C := nu int_0^tau dt/r(t)^2 = sigma_y/r_*^2` (Corollary V.5, unchanged and exact),
> `Pfrak := P(|Z_tau| > d/2) + P(|Z^a_tau| > d/2) + P(|Z^a_tau| >= d)`, `N := ||eta_0||_inf r_*/M`,
> `W` as in the coupling lemma (5.2), and put
> ```
>     Q  := (1 + vartheta)/(1 - vartheta) ,      vt := vartheta/(1 - vartheta) .
> ```
> Then
> ```
>   | eta(X(tau),tau)/eta_0(x_*) - (1 - s_C) |  <=  E_vartheta + Q ( E_hess + E_4 + E_tail ) ,
>
>   E_vartheta := vt ( 9 + 2 sigma_z/sigma_y ) s_C ,
>   E_hess     := (r_*/R_-^2) (1/2) K_2 e^{c} tau V_1  +  4 N (1/2) K_2 e^{c} tau V_1 / d ,
>   E_4        := (r_*/R_-^5) [ (tr C_tau)^2 + 2 tr(C_tau^2) ] ,
>   E_tail     := 2 N Pfrak + Pfrak + sum_{k=1}^{3} r_*^{-k} (E|Z^a_tau|^{2k})^{1/2} Pfrak^{1/2} ,
> ```
> with `E_hess`, `E_4`, `E_tail` exactly V.4's in the V-b refuter's corrected form. At
> `vartheta = 0`, `Q = 1`, `vt = 0`, `E_vartheta = 0`, and the statement is V.4 verbatim.

*Proof.* Only Step 2's algebra and the normalisation change; Step 0 (the Feynman-Kac
representation `eta(X(tau),tau) = E[eta_0(x_* + Z_tau)]`, from Theorem V.1) and the coupling
lemma (5.2) are untouched, because neither mentions `eta_0`.

*(i) The normalisation.* By `k = 0`, `|eta_0(x_*) - eta_P(x_*)| <= vartheta M/r_*`, so
`|eta_0(x_*)| >= (1-vartheta) M/r_*` and `|eta_P(x_*)/eta_0(x_*) - 1| <= vt`.

*(ii) Step 2's Hessian term.* `Z^a` is exactly Gaussian with mean zero, so for the third-order
Taylor polynomial `P_3^{eta_0}` of `eta_0` at `x_*`,
`E[P_3^{eta_0}(Z^a)] = eta_0(x_*) + (1/2) C_tau : grad^2 eta_0(x_*)` exactly. Split
`grad^2 eta_0 = grad^2 eta_P + grad^2 e`. Lemma V.3' gives
`(1/2)C_tau : grad^2 eta_P(x_*) = sigma_y Delta_y eta_P(x_*) = -eta_P(x_*) sigma_y/r_*^2 =
-eta_P(x_*) s_C`, with `sigma_z` dropping out identically. For the perturbation, `k = 2` gives
`||grad^2 e||_op = sup_{|v|=1}|d_v^2 e| <= 2 vartheta M/r_*^3`, so, `C_tau` being positive
semidefinite,
```
   | (1/2) C_tau : grad^2 e(x_*) |  <=  (1/2) tr(C_tau) ||grad^2 e||_op
                                    =  (1/2)(8 sigma_y + 2 sigma_z)(2 vartheta M/r_*^3) .
```
Dividing by `|eta_0(x_*)| >= (1-vartheta)M/r_*` and adding the normalisation error of (i),
```
   | E[P_3^{eta_0}(Z^a)]/eta_0(x_*) - (1 - s_C) |  <=  vt s_C  +  vt (8 + 2 sigma_z/sigma_y) s_C
                                                    =  E_vartheta .
```

*(iii) Steps 1 and 2's remainders.* On `B(x_*,d)`, `k = 1` gives `|grad eta_0| <=
(1+vartheta) M/R_-^2`, `k = 4` gives `sup|d_v^4 eta_0| <= (1+vartheta) 4! M/R_-^5`, and
`k <= 3` gives `sup|d_v^k eta_0| <= (1+vartheta) k! M/r_*^{k+1}`. Each of V.4's estimates is
divided by `|eta_0(x_*)| >= (1-vartheta)M/r_*` rather than by `M/r_*`, and each has a
numerator that has grown by at most `(1+vartheta)`. Every one of `E_hess`, `E_4`, `E_tail`
therefore acquires at most the single factor `(1+vartheta)/(1-vartheta) = Q`: Step 1's first
piece through `|grad eta_0| <= (1+vartheta)M/R_-^2` and its second through
`N <= ||eta_0||_inf r_*/((1-vartheta)M)`; `E_4` through `sup|d_v^4 eta_0| <=
(1+vartheta) 4! M/R_-^5`; `E_tail` through both. ∎

**Corollary V.5' (the budget's `eps_v`).** With the notation of ASSEMBLY section 2.5,

```
   eps_v(vartheta)  =  eps_bulk [ 1 + vt (9 + 2 sigma_z/sigma_y) ]
                       +  Q ( E_hess + E_4 + E_tail ) .
```

At the theorem's `c = c_*` the ratio `sigma_z/sigma_y` is computed by the imported instrument
(it is `J_pow(c,4)/J_pow(c,-2)`), so `E_vartheta = 16.026 vt eps_bulk` at Corollary V.5's own
`sigma_z/sigma_y = 3.5131`.

## B3. `vartheta` for the theorem's datum

`d` is ASSEMBLY section 2.5's rule as coded in `THEOREM_S3/t3_budget.py`:
`d = min(f, rho_* sin(phi_0 - delta), rho_* sin(pi/2 - delta_m - phi_0))`, with
`x_* = ((1+f)rho_0, phi_0)`, `phi_0 = 30 deg`. The first entry is the **sharp**-edge inner
constraint `rho_* - d >= rho_0`; it has no meaning for the `tanh` ramp, and the `tanh` ramp's
cost is exactly what `vartheta` measures. Both variants are computed.

`e = M psi(rho)/r` restricted to a line `x + s v` is
`M e^2 (R(s)^2)^{-1/2} (e^2 + (P(s)^2)^4)^{-1}` with `R(s)^2 = r^2 + 2 r a s + (a^2+b^2)s^2`
and `P(s)^2 = rho^2 + 2 (x.v) s + s^2` both quadratic, so every `d_v^k e` is obtained
**exactly** by formal power series (Miller recurrence for `(R^2)^{-1/2}`, polynomial powers,
reciprocal recurrence, one convolution). Against `sympy` at a test point the agreement is
`0`, `1.4e-16`, `0`, `1.0e-15`, `2.6e-15` for `k = 0,...,4`. The supremum over
`B(x_*,d) x S^4` is then a five-parameter search (three for the point, two for the direction,
after using the `SO(4)` symmetry), done on a grid with six rounds of vectorised local
refinement; the coarse-to-fine change is `5.1e-04` at `f = 1` and `6.7e-03` at `f = 8`. As an
exact control, `vartheta_0 = psi(rho_* - d)` identically, and the computed value reproduces
it to `0`.

**`vartheta` at ASSEMBLY's `d`:**

| `f` | `d` | `rho_min = rho_* - d` | `Theta(rho_*)` | `psi(rho_min) = vartheta_0` | `vartheta` | binding `k` |
|---|---|---|---|---|---|---|
| `0.05` | `0.0500` | `1.0000` | `0.166633` | `8.8080e-01` | `1.3007` | 1 |
| **`0.1`** | `0.1000` | `1.0000` | `0.224868` | `8.8080e-01` | **`1.4593`** | 1 |
| `0.25` | `0.2500` | `1.0000` | `0.446493` | `8.8080e-01` | `5.0734` | 4 |
| `0.5` | `0.5000` | `1.0000` | `0.776211` | `8.8080e-01` | `9.6389` | 4 |
| **`1`** | `0.7654` | `1.2346` | `0.971946` | `5.7781e-01` | **`6.8754`** | 4 |
| `2` | `1.1481` | `1.8519` | `0.998875` | `5.0695e-02` | `2.4235` | 4 |
| `4` | `1.9134` | `3.0866` | `0.999981` | `8.9614e-04` | `0.055612` | 4 |
| `8` | `3.4442` | `5.5558` | `1.000000` | `8.1392e-06` | `5.0751e-04` | 4 |
| `16` | `6.5056` | `10.4944` | `1.000000` | `5.0227e-08` | `3.1320e-06` | 4 |
| `32` | `12.6286` | `20.3714` | `1.000000` | `2.4912e-10` | `1.5535e-08` | 4 |

The brief's three values, with the per-order breakdown at the finer grid:

| `f` | `vartheta_0` | `vartheta_1` | `vartheta_2` | `vartheta_3` | `vartheta_4` | `vartheta` |
|---|---|---|---|---|---|---|
| `0.5` | `0.880797` | `1.922101` | `3.275772` | `4.434305` | `9.641744` | `9.64174` |
| `1` | `0.577815` | `1.663234` | `3.164467` | `5.001180` | `6.871935` | **`6.87194`** |
| `2` | `0.050695` | `0.233885` | `0.652847` | `1.372662` | `2.420944` | `2.42094` |

**Three readings.**

* **(H3'_vartheta) is vacuous at the `f` the budget picks.** At `f = 0.1`, `vartheta = 1.4593`,
  and the hypothesis requires `vartheta < 1` for `Q` and `vt` to be finite. At `f = 1`,
  `vartheta = 6.87`. The hypothesis first has content at `f = 4`.
* **The mechanism is one line.** `psi ~ e^2 (rho_min/rho_0)^{-8}` and, when `d` is
  taper-limited, `rho_min = (1 - sin(phi_0 - delta)) rho_* = 0.61732 rho_*`, so `vartheta`
  falls like `f^{-8}`. To get `psi` below `1e-02`, `1e-03`, `1e-04`, `1e-06` needs
  `rho_min/rho_0` above `2.2834`, `3.0449`, `4.0604`, `7.2206` respectively.
* **The binding order is `k = 4`**, except at very small `f` where `k = 1` binds because
  `d <= f` shrinks the ball faster than the ramp decays. So the fourth-derivative clause, the
  one the two-clause repair does not contain, is the one that actually decides `vartheta`.

## B4. The second defect at the same place

Independently of (H3-V): at the budget's minimiser `f = 0.1` the tracked material point sits
at `u_* = log 1.1 = 0.0953`, **below** the `tanh` ramp's half-height `u = eps_r = 0.25`, and
`Theta(rho_*) = 0.224868`. The tracked point carries `22.5 %` of the nominal amplitude `M`.
`ell_loss(f) = log 2 + log(1+f) + 3 log lam_max + log((1+mu)/(1-mu)) + eps_r` charges
`log(1+f)` for the inset and a flat `eps_r` for the outer ramp, and nothing in it forces
`Theta(rho_*)` anywhere near `1`; the `f`-minimisation is therefore free to put the tracked
point inside the inner ramp, and it does. `Theta(rho_*)` first exceeds `0.999` at `f = 2`.

This is not an item on THEOREM_S3's modulo list and I am not adding one. It is a bookkeeping
question about `ell_loss` for the `tanh` datum, it points the same way as (H3-V), and it is
fixed by the same move: require `f` large enough that the tracked ball is in the plateau.

## B5. The budget with `vartheta` carried

{{B5}}

---

## C. FINDINGS

{{FINDINGS}}

---

## GATE, AND FILES

{{GATE}}

| file | what it establishes |
|---|---|
| `q1_identities.py` / `q1_results.json` | every algebraic identity and sign: the weight identity in general `n` and `k`, the `eps`-regularised identities in `n = 5` and the sign of their zeroth-order terms, the regularised Kato remainders, the gradient equation and the absence of a zeroth-order term, the barrier coefficient bounds, the `sym grad_5 b` spectrum |
| `q2_datum.py` / `q2_results.json` | `E_0` and `Gfrak_0` with their closed forms and a scan that confirms them; `||eta_0||_inf`, `||grad eta_0||_inf`; the decay-lemma tail constants; `TV(eta_0)`; the mollification control that licenses Proposition 3 |
| `q3_theta.py` / `q3_results.json` | `vartheta` for the theorem's datum, per order `k` and per `f`, with the series-vs-sympy control and the grid-refinement control |
| `q4_budget_theta.py` / `q4_results.json` | THEOREM_S3's budget re-run with `vartheta` carried; the control that reproduces THEOREM_S3's own `L_*` from the imported instrument |
| `copies/` | byte copies, with hashes, of every file imported from another seat |
| `check_constants.py` | the gate and the document audit |
| `SHA256SUMS` | computed with `shasum -a 256`, never typed |

**FL-000 stands. Nothing here touches the headline problem.**

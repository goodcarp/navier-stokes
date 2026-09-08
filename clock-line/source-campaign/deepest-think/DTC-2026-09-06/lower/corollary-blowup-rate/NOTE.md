# NOTE — sub-seat `corollary-blowup-rate` — DTC-2026-09-06 / lower

**Task.** State and prove the blow-up-rate corollary of theorem (i) (Astra pass 8, the
logarithmic record clock), then sweep the literature *at source* for any published lower
bound on `||omega(t)||_inf` near a putative Navier–Stokes singularity, and compare.

**Register.** The corollary is a five-line consequence of a theorem I did not prove and did
not re-verify; it is exactly as strong as its input. What I did verify myself: the algebra of
the three forms, the scaling invariance of all three, the bootstrap that makes form (c)
explicit, the classical comparator's exponents, and the source text of every paper cited.
Everything numeric came from `scripts/s1_corollary_algebra.py` and
`scripts/s2_leray_nu_and_invariance.py`, both listed in `SHA256SUMS`; logs are `scripts/s1_log.txt`,
`scripts/s2_log.txt`. Numerics falsify, never prove — blocks (B) and (E) are counterexample
searches, not proofs; blocks (A), (C), (D) and all of `s2` are symbolic identities.

---

## 0. Bottom line

1. **The corollary is proved** (§2–§3), in three forms: implicit (a), energy-monotone (b),
   and fully explicit in `(T-t)` (c). All three are Navier–Stokes-scaling invariant
   (`s2` block 2, three exact zeros).
2. **It beats the only classical route to a vorticity lower bound** — Leray/Giga in `L^inf`
   pushed through the interpolation `||u||_inf <= C E^(1/5) M^(3/5)` — which gives exponent
   `5/6`, not `1`. The ratio grows like `(T-t)^(-1/6)/log` (§5, fitted slope 0.131–0.154
   towards 1/6).
3. **Novelty sweep verdict: NOVEL-CONDITIONAL.** I located **no** published lower bound on
   `||omega(t)||_inf` for Navier–Stokes, conditional or otherwise, **except** one that is
   strictly stronger and implies ours: **Bradshaw–Farhat–Grujić, ARMA 2019, Theorem 8/10**,
   whose proof is a self-described *sketch* with a defective step already located by two
   referees of this campaign. Until that conflict is settled in print, the honest statement is:
   *the corollary is new relative to everything else found, and is implied by a published
   theorem whose proof does not currently support it.*
4. The closest true analogue in the literature is **Ingimarson–Kukavica, arXiv:2603.17431v2
   (21 Mar 2026), Theorem 2.1** — same `1/((T-t) log)` shape, but for **Euler**, only in
   `lim sup`, with a constant depending on `||u_0||_{H^3}`. Our log carries coefficient
   **1/5** where theirs carries **1**, holds **for every** `t < T`, and its constant depends
   only on `E_0` and `nu` through a scaling-invariant Reynolds number.
5. **Correction to a machine-summarised source.** A tool summary of arXiv:1503.03063 reported
   "Theorem 1.1 (Robinson, Sadowski, Silva): `||omega(t)||_inf >= C(T-t)^(-1)` for 3D
   Navier–Stokes". That is a fabrication. The paper is **Cortissoz–Montero**, is about
   `Ḣ^s` norms of the **velocity**, and the strings `vorticity`, `omega` and `curl` do not
   occur in it at all (`grep`, §6.1). Do not propagate it.

---

## 1. Setting, and the exact class the input theorem lives in

`nu > 0` fixed, unforced 3D Navier–Stokes on `R^3`. Initial datum `u_0` divergence-free,
finite energy, `u_0 in H^m(R^3)` with integer `m >= 4` (this is verbatim the class of the
input theorem: *"the maximal strong solution with smooth finite-energy initial data, for
example divergence-free `H^m` data with integer `m >= 4`"*, LOGARITHMIC_RECORD_CLOCK.md §1).
Let

    u in C([0, T_max); H^m_sigma(R^3)) ,   T_max in (0, infinity]

be **the** maximal strong solution, `T_max` its maximal existence time in that class, and

    E(t) = ||u(t)||_2^2 ,  M(t) = ||omega(t)||_inf ,  Re_E(t) = E(t)^(2/5) M(t)^(1/5) / nu ,
    E_0 = E(0) .

Two facts about the class that the argument uses and that `m >= 4` supplies:

* `H^{m-1} ⊂ L^inf` for `m-1 > 3/2`, so `M(t) < infinity` and `t -> M(t)` is continuous on
  `[0, T_max)` (needed by the input theorem's own first-exit step, not by mine);
* the energy equality `d/dt E = -2 nu ||grad u||_2^2 <= 0` holds, so `E(t) <= E_0`.

**Blow-up hypothesis.** *The solution blows up at `T`* means `T := T_max < infinity`.

**The input (theorem (i), Astra pass 8, referee-confirmed end to end).** There is a
universal constant `c_1 in (0,1)` such that: for every `t_0 < T_max` with `E(t_0), M(t_0) > 0`,
writing `A(t_0) = 1 + log_+ Re_E(t_0)` and

        H(t_0) = c_1 / ( M(t_0) A(t_0) ) ,                                  (H)

the solution continues on `[t_0, t_0 + H(t_0)]` inside the class above, and
`sup_{[t_0, t_0+H(t_0)]} ||omega||_inf <= (3/2) M(t_0)`.
(That `c_1 < 1` is explicit in the source: *"Let `A = 1+log_+ R_0` and `H = c/(M_0 A)` with
`0 < c < 1` … Choose a universal `c` small enough"*, §5. I use `c_1 < 1` only to drop a
`min(c_1,1)` in form (c); block (B) shows the guard is otherwise needed.)

---

## 2. The corollary

> **COROLLARY (blow-up rate of the vorticity maximum).**
> Let `nu > 0` and let `u` be the maximal strong solution with divergence-free finite-energy
> `u_0 in H^m(R^3)`, `m >= 4`, and suppose it blows up at `T < infinity`. Let `c_1 in (0,1)`
> be the constant of (H). Then for every `t in [0,T)`:
>
> **(a) implicit form.**
>
>       ||omega(t)||_inf  >=  c_1 / ( (T-t) * (1 + log_+ Re_E(t)) ) ,
>       Re_E(t) = E(t)^(2/5) ||omega(t)||_inf^(1/5) / nu .
>
> **(b) energy-monotone form.** Since `E(t) <= E_0`,
>
>       ||omega(t)||_inf  >=  c_1 / ( (T-t) * (1 + log_+( E_0^(2/5) ||omega(t)||_inf^(1/5) / nu )) ) .
>
> **(c) explicit form.** With `K := E_0^(2/5) / nu`,
>
>       ||omega(t)||_inf  >=  c_1 / ( (T-t) * (1 + log_+( K (T-t)^(-1/5) )) )
>                          =  c_1 / ( (T-t) * (1 + [ (2/5)log E_0 - log nu + (1/5)log(1/(T-t)) ]_+ ) ) .
>
> All three are invariant under `u -> lambda u(lambda x, lambda^2 t)` at fixed `nu`.

Form (c) is the usable one: its right-hand side involves only `E_0`, `nu` and `T-t`.
Asymptotically `(c)` reads

       ||omega(t)||_inf  >=  ( 5 c_1 + o(1) ) / ( (T-t) log(1/(T-t)) )    as t -> T^- .

---

## 3. Proof

### 3.1 The five-line continuation argument (form (a))

Fix `t in [0,T)`.

*Line 0 (non-degeneracy).* If `M(t) = 0` then `omega(t) = 0`; a finite-energy divergence-free
field on `R^3` with vanishing curl is `0` (Biot–Savart, or: `u(t)` is harmonic in each component
and in `L^2`), so `u(t) = 0`, and by uniqueness `u ≡ 0` on `[t, infinity)`, contradicting
`T_max = T < infinity`. Likewise `E(t) > 0`. So (H) applies at `t_0 = t`.

*Line 1.* Suppose, for contradiction, `T - t < H(t)`.

*Line 2.* By theorem (i) applied at `t_0 = t`, the strong solution with datum `u(t)` at time
`t` exists in the class `C([t, t+H(t)]; H^m_sigma)` — an interval that strictly contains
`[t, T]` by Line 1.

*Line 3.* By uniqueness of the maximal strong solution in that class, that solution coincides
with `u` on `[t, T)`; hence `u` extends to a strong solution on `[0, t + H(t)]`, an interval
strictly containing `[0,T)`.

*Line 4.* This contradicts `T = T_max` being the maximal existence time.

*Line 5.* Therefore `T - t >= H(t) = c_1/(M(t)(1+log_+ Re_E(t)))`, i.e.
`M(t) (T-t) (1 + log_+ Re_E(t)) >= c_1`, which is (a). ∎

Two things this argument does **not** do, worth stating because they are the usual failure
modes: it never uses the cap `(3/2)M_0` (only the continuation half of theorem (i)), and
`H(t)` never refers to `T`, so there is no circularity.

### 3.2 Form (b)

`E(t) <= E_0` by the energy equality, and `r -> 1 + log_+ r` is non-decreasing, so
`1 + log_+ Re_E(t) <= 1 + log_+(E_0^(2/5) M(t)^(1/5)/nu)`. Dividing the (positive) constant
`c_1/(T-t)` by the larger denominator weakens the bound, so (a) implies (b). ∎

### 3.3 Form (c) — removing `M(t)` from the logarithm

Write `s = T-t`, `Re_* = K s^(-1/5) = E_0^(2/5)/(nu s^(1/5))`, and
`m := c_1 / ( s (1 + log_+ Re_*) )`. Suppose `M(t) < m`. Since `1 + log_+ Re_* >= 1` and
`c_1 < 1`, we get `M(t) < c_1/s < 1/s`, hence

        E_0^(2/5) M(t)^(1/5) / nu  <  E_0^(2/5) (c_1/s)^(1/5) / nu  =  c_1^(1/5) Re_*  <=  Re_* ,

so `log_+( E_0^(2/5) M(t)^(1/5)/nu ) <= log_+ Re_*`. Feeding that into (b),

        M(t)  >=  c_1 / ( s (1 + log_+(E_0^(2/5)M(t)^(1/5)/nu)) )  >=  c_1 / ( s (1 + log_+ Re_*) )  =  m ,

contradicting `M(t) < m`. Hence `M(t) >= m`, which is (c). ∎
(For a general `c_1 > 0` the same proof gives (c) with `c_1` replaced by `min(c_1,1)`; the
guard is not cosmetic — block (B) below exhibits violations without it.)

### 3.4 Scaling

Under `u -> lambda u(lambda x, lambda^2 t)` at fixed `nu`: `E -> E/lambda`,
`M -> lambda^2 M`, `T-t -> (T-t)/lambda^2`. Then `M(T-t)` is invariant, `Re_E` is invariant
(`s1` block A), and `E_0^(2/5)/(nu (T-t)^(1/5))` is invariant. `s2` block 2 checks all three
combinations symbolically: **scaled minus original = 0, 0, 0.**

### 3.5 Verification receipts

`s1` block (A):

    ell_lam/ell = 1/lambda        Re_E(lam)/Re_E = 1
    Re_E = E**(2/5)*M**(1/5)/nu   (M*T)_lam/(M*T) = 1

`s1` block (B) — 4,000,000 random `(E_0, nu, s, c_1, M)` over 12–22 decades each,
`c_1` sampled up to `10^0.6 > 1` on purpose:

    (b) true                       = 2,053,207
    (b) true and (c) false         = 0            <-- the §3.3 step survives
    (b) true and (c-without-min)   = 11           <-- the min(c_1,1) guard is load-bearing

---

## 4. Consequence: a quantitative BKM accumulation bound

Integrating (c) in `t` on the branch where the `log_+` is active (`s1` block (D) verifies the
antiderivative `F(r) = -5 c log(1 + log K - (1/5) log r)` by differentiation: `F' - integrand = 0`):

        int_0^t ||omega(s)||_inf ds  >=  5 c_1 * log log ( 1/(T-t) )  -  C(E_0, nu, T) .

So the corollary reproves the divergence `int_0^T ||omega||_inf = infinity` (Beale–Kato–Majda /
Serrin-type) with an explicit `log log` rate whose constant `5 c_1` is universal and whose
error term depends only on `E_0, nu, T`. This is the exact viscous counterpart of
Ingimarson–Kukavica's Euler bound `int_0^t ||omega||_inf >= (1/C) log log(1/(T_*-t))`
(their (2.1), C depending on `||u_0||_{H^3}`) — see §7.

---

## 5. Comparison with the only classical route (Leray/Giga through the interpolation)

Leray's `L^q` lower bounds (`q in (3, inf]`, proved by Giga) read, at `nu = 1`,
`||u(t)||_q >= c_q (T-t)^{-(1-3/q)/2}`. `s2` block 1 restores the viscosity by the exact
normalisation `w(x, sigma) = u(x, sigma/nu)/nu`, which solves the `nu=1` equations and blows
up at `nu T`:

        ||u(t)||_q  >=  c_q * nu^((1+3/q)/2) * (T-t)^(-(1-3/q)/2) ,
        q = inf :  ||u(t)||_inf >= c nu^(1/2) (T-t)^(-1/2) .

Feed that through the input theorem's own velocity interpolation
`||u||_inf <= C E^(1/5) M^(3/5)` (LOGARITHMIC_RECORD_CLOCK.md eq. (4)). `s1` block (C):

        M(t)  >=  c' * nu^(5/6) * E_0^(-1/3) * (T-t)^(-5/6) ,
        check:  E^(1/5) * (nu^(5/6) E^(-1/3) s^(-5/6))^(3/5)  =  sqrt(nu)/sqrt(s)   ✓

**So the best pre-existing route to a vorticity lower bound gives exponent `5/6`.**
The corollary gives exponent `1` minus a logarithm. With both unknown constants set to 1
(`s1` block (E)) the ratio corollary/classical is:

| `E_0` | `nu` | `T-t = 1e-6` | `1e-9` | `1e-12` | `1e-15` |
|---|---|---|---|---|---|
| 1 | 1 | 2.66 | 6.15 | 15.3 | 40.0 |
| 1 | 1e-2 | 55.5 | 150 | 417 | 1.17e3 |
| 1e3 | 1e-3 | 2.35e3 | 6.75e3 | 1.95e4 | 5.69e4 |
| 1e-3 | 1e-1 | 2.06 | 4.60 | 11.2 | 28.9 |

Fitted slope of `log(ratio)` against `log(1/(T-t))` on `[1e-15, 1e-6]`: **0.1314** (`E_0=nu=1`)
and **0.1537** (`E_0=1e3, nu=1e-3`), rising towards the asymptotic **1/6 = 0.1667** as the log
factor's relative weight decays. The crossover *location* depends on the two unknown constants
and is not a claim; the exponent gap `1/6` is.

---

## 6. Novelty sweep at source

### 6.1 Method and correction

Queries run (WebSearch, Sept 2026): *lower bound vorticity blow-up rate Navier-Stokes L
infinity norm*; *"lower bounds" blow-up Navier-Stokes "vorticity" "(T-t)" rate*; *Navier-Stokes
"||omega(t)||" lower bound "1/(T-t)" vorticity maximum norm singularity local existence bounded
vorticity*; *Ingimarson Kukavica … viscous analogue quantitative BKM*; *Giga 1986 … Leray*;
*Bradshaw Farhat Grujić … citing papers*; *arXiv 2026 Navier-Stokes lower bound vorticity
L^infty blow-up rate logarithm Reynolds energy dependent lifespan*. Every hit that could
plausibly contain a vorticity lower bound was **downloaded as a PDF and read** (`pdf/`, text in
`txt/`, page numbers located with `scripts/pagefind.sh`).

**A machine summary of a source was false and is corrected here.** A `WebFetch` of
`arxiv.org/pdf/1503.03063` returned "*Theorem 1.1 (Page 1, Robinson, Sadowski, Silva):
`||omega(t)||_inf >= C(T-t)^{-1}` … Theorem 2.1 (Page 3) … Euler … exponent −1*". None of that
is in the paper. The paper is **Cortissoz–Montero**, *Lower bounds for possible singular
solutions for the Navier–Stokes and Euler equations revisited*, arXiv:1503.03063v2 (2 Sep 2016)
= JMFM 19 (2017); it is about `Ḣ^s` norms of the **velocity** on `T^3` at `nu = 1`; and

    grep -c -i -E "vortic|omega|curl" txt/rss-1503.03063.txt  ->  0

Its actual Theorem 2.1 (PDF p.3) is `C_s t^{-(1/2)(s-1/2)} <= ||u(T-s)||_{Ḣ^s(T^3)}`,
`1/2 < s < 5/2`; its Corollary 2.1 (PDF p.5) is the matching lifespan bound
`K_s ||u_0||_s^{-4/(2s-1)} <= T`. Not a vorticity statement.

### 6.2 What exists, at source

* **Bradshaw–Farhat–Grujić**, *An algebraic reduction of the 'scaling gap' in the Navier–Stokes
  regularity problem*, arXiv:1704.05546v4, ARMA 231 (2019) 1983–2005, DOI
  10.1007/s00205-018-1314-5. **Theorem 8** (arXiv v4 **p.8**, verified by me in `txt/`):
  *"Let the initial datum `omega_0` be in `L^2 ∩ L^inf`. Then there exists a unique mild
  solution `omega` in `C_w([0,T], L^inf)` where `T >= (1/c)(1/||omega_0||_inf)` for an absolute
  constant `c > 0`."* **Theorem 10** (**p.11**) adds `||omega(t)||_{L^inf(Omega_t)} <= M ||omega_0||_inf`
  on that interval. The paragraph immediately before Theorem 8 (**p.8**) says the `L^2`
  hypothesis is *"a 'soft assumption', i.e., there will be no quantitative dependence on
  `||omega_0||_2` in the proof."* Applied at `t` with `omega_0 = omega(t)`, Theorem 8 gives
  `||omega(t)||_inf >= 1/(c(T-t))` — **log-free, energy-free, `nu`-uniform**, i.e. strictly
  stronger than the corollary. The proof is introduced as *"we present a sketch here"* (p.8);
  two referees of this campaign located the defective step (the BMO estimate
  `int |f|/(|x|+1)^4 dx <= c ||f||_BMO` used without the local average, pp.9–10, whose deficit
  on the plateau family is exactly `log(R/rho_0)`). No published erratum; no citing paper
  contesting it (searched).
* **Ingimarson–Kukavica**, *Lower bounds on the blowup rate of vorticity in the Euler
  equations*, arXiv:2603.17431v2, 21 Mar 2026. **Euler only.** Abstract (**p.1**): *"not much
  is known regarding the rate at which either of the two quantities … blow up"* (of BKM's
  `int ||omega||_inf` and `limsup ||omega||_inf`). **Theorem 2.1** (**p.3**):
  `int_0^t ||omega(s)||_inf ds >= (1/C) log log(1/(T_*-t))` for `T_*-t < 1/C`, and
  `limsup_{t->T_*^-} (T_*-t) log(1/(T_*-t)) ||omega(t)||_inf >= 1/C`, with **`C` depending on
  `||u_0||_{H^3}`**. Their **Theorem 2.5** (p.3) is pointwise-in-time but for `||D omega||_inf`
  (`>= 1/(C(T_*-t)^{7/5})`), not for `omega`. Their route is BKM + Kozono–Taniuchi
  `||Du||_inf <= C(1 + ||omega||_inf(1 + log_+ ||u||_{H^3}))` (their (3.3), p.4) + Gronwall;
  the `H^3` norm inside the logarithm is what forces the `log log` and the non-invariant
  constant.
* **Luo**, *Sharp norm inflation for 3D Navier-Stokes equations in supercritical spaces*,
  arXiv:2504.08288, **p.3** (verified verbatim in `txt/`): *"no rigorous lower bounds for (1.5)
  `[|omega(t)|_inf/|omega_0|_inf]` have been established for solutions of the Navier-Stokes
  equations to date. In other words, there remains no proof of growth—let alone sustained
  growth or blowup—for solutions of (1.1)."* His Theorem 1.3 (p.3) is an unconditional growth
  *ratio* at a long clock (`M_0 t_* -> infinity`), not a rate at a singularity.
* **Cortissoz–Montero**, arXiv:1503.03063v2 — `Ḣ^s` velocity, §6.1 above. Their §4 records
  that the `Ḣ^{1/2}` blow-up rate is still unknown, and credits the `L^p` bounds to Leray
  (stated without proof) and Giga (proved).
* **Leray**, Acta Math. 63 (1934) 193–248 (`L^p` lower bounds stated) / **Giga**, JDE 62 (1986)
  186–212 (proved): `||u(t)||_q >= c_q nu^{(1+3/q)/2}(T-t)^{-(1-3/q)/2}`, `q in (3, inf]`
  (viscosity power restored in `s2`).
* **Tao**, *Quantitative bounds for critically bounded solutions to the Navier–Stokes
  equations* (2019): `limsup_{t->T_*} ||u(t)||_{L^3} / (log log log(1/(T_*-t)))^c = +infinity`.
  Critical **velocity** norm; his Theorem 1 gives an *upper* bound on `|omega|` in terms of `A`.
* **Barker**, *Quantitative classification of potential Navier-Stokes singularities beyond the
  blow-up time*, arXiv:2510.20757 (v3, 11 Aug 2026): quantitative lower bounds near a blow-up
  time via Tao/Carleman, for approximately axisymmetric data; localized `L^2` vorticity
  concentration, not `||omega(t)||_inf`.
* **Chen–Strain–Tsai–Yau**, arXiv:math/0701796 and arXiv:0709.4230: axisymmetric-with-swirl
  Type-I **velocity** criteria (`|v| <= C_*(r^2-t)^{-1/2}` implies regularity). Different
  quantity, different class.
* Not lower bounds, checked and excluded: Kozono–Ogawa–Taniuchi (Kyushu J. Math. 57 (2003)
  303–324) — a `Ḃ^0_{∞,∞}` continuation criterion, no rate; Grujić arXiv:2607.08866 — a
  conditional regularity result; Kim–Jeong JFA 283 (2022), Bang–Cheskidov arXiv:2605.16502,
  Lim–Jeong ARMA 249 (2025), Egamberganov–Yao arXiv:2512.13456 — all **Euler** upper bounds or
  ill-posedness, per the campaign's `sharp/literature-short-time` seat, which I did not re-run.

---

## 7. Comparison table

| # | Statement | Equation | Quantity | Rate | Constant depends on | For which `t` | Source (verified at source) |
|---|---|---|---|---|---|---|---|
| **0** | **This corollary (c)** | **NS, `nu>0`** | **`\|\|omega(t)\|\|_inf`** | **`c_1 /((T-t)(1+log_+(E_0^{2/5} nu^{-1}(T-t)^{-1/5})))`, i.e. `~ 5c_1/((T-t)log(1/(T-t)))`** | **universal `c_1`; `E_0, nu` only inside a scaling-invariant `log`** | **every `t<T`** | **§2–3 here; input = LOGARITHMIC_RECORD_CLOCK.md** |
| 1 | `T >= (1/c)(1/\|\|omega_0\|\|_inf)`, `c` absolute ⟹ `\|\|omega(t)\|\|_inf >= 1/(c(T-t))` | NS | `\|\|omega\|\|_inf` | `(T-t)^{-1}`, **no log** | absolute; **no `E`, no `nu`** | every `t<T` | Bradshaw–Farhat–Grujić, arXiv:1704.05546v4 **Thm 8, p.8**; Thm 10 p.11; ARMA 231 (2019) 1983–2005. ⚠ proof is a *sketch*; defective BMO step pp.9–10 |
| 2 | `limsup (T_*-t) log(1/(T_*-t)) \|\|omega(t)\|\|_inf >= 1/C` | **Euler** | `\|\|omega\|\|_inf` | `(T-t)^{-1}/log`, **log coefficient 1** | `\|\|u_0\|\|_{H^3}` (not scaling-invariant) | **`limsup` only** | Ingimarson–Kukavica, arXiv:2603.17431v2 **Thm 2.1, p.3** |
| 3 | `int_0^t \|\|omega\|\|_inf ds >= (1/C) loglog(1/(T_*-t))` | **Euler** | accumulation | `loglog` | `\|\|u_0\|\|_{H^3}` | `T_*-t < 1/C` | same, **(2.1), p.3** |
| 3' | *(implied by row 0)* `int_0^t \|\|omega\|\|_inf ds >= 5c_1 loglog(1/(T-t)) - C(E_0,nu,T)` | **NS** | accumulation | `loglog`, coefficient **`5c_1` universal** | `E_0, nu, T` in the additive error only | `t` near `T` | §4 here |
| 4 | `\|\|D omega(t)\|\|_inf >= 1/(C(T_*-t)^{7/5})` | **Euler** | `\|\|D omega\|\|_inf` | `7/5` | `\|\|u_0\|\|_{L^2}` | every `t<T_*` | Ingimarson–Kukavica **Thm 2.5, p.3** |
| 5 | `\|\|u(t)\|\|_q >= c_q nu^{(1+3/q)/2}(T-t)^{-(1-3/q)/2}`, `q in (3,inf]` | NS | **velocity** `L^q` | `1/2` at `q=inf` | `c_q`; `nu` power restored in `s2` | every `t<T` | Leray, Acta Math. 63 (1934); proved by Giga, JDE 62 (1986) 186–212 |
| 5' | *(row 5 + `\|\|u\|\|_inf <= C E^{1/5}M^{3/5}`)* `\|\|omega(t)\|\|_inf >= c nu^{5/6}E_0^{-1/3}(T-t)^{-5/6}` | NS | `\|\|omega\|\|_inf` | **`5/6`** | `c`, `E_0`, `nu` **algebraically** | every `t<T` | §5 here (`s1` block C) — the classical baseline row 0 beats |
| 6 | `C_s t^{-(1/2)(s-1/2)} <= \|\|u(T-t)\|\|_{Ḣ^s(T^3)}`, `1/2<s<5/2` | NS (`nu=1`, `T^3`) | **velocity** `Ḣ^s` | `(2s-1)/4` | `C_s` | every `t<T` | Cortissoz–Montero, arXiv:1503.03063v2 **Thm 2.1, p.3** (**not** a vorticity bound) |
| 7 | `limsup \|\|u(t)\|\|_{L^3}/(logloglog(1/(T_*-t)))^c = +infinity` | NS | **velocity** `L^3` (critical) | `logloglog` | absolute | `limsup` | Tao (2019), *Quantitative bounds for critically bounded solutions* |
| 8 | *"no rigorous lower bounds for `\|omega(t)\|_inf/\|omega_0\|_inf` have been established for solutions of the Navier-Stokes equations to date"* | NS | — | — | — | — | Luo, arXiv:2504.08288, **p.3** |

### How row 0 compares, precisely

* **vs row 1 (the only stronger NS statement).** Row 1 implies row 0 and drops the logarithm
  entirely. Row 0 is therefore *not* new if BFG Theorem 8 stands. The campaign's two referees
  located a defective step in its sketch; that is not a published erratum, so the correct
  register is **conflict unresolved**, not "refuted". Note the direction of the stake: the
  campaign's conjecture (ii) — a matching *upper* bound `c_2/log Re` — would refute row 1, and
  row 0 is the theorem row 1 would subsume.
* **vs row 2 (the closest analogue).** Same `1/((T-t) log)` shape. Four differences, all in
  row 0's favour except the equation: (i) row 2 is **Euler**, row 0 is **Navier–Stokes with
  `nu > 0`**, so neither implies the other; (ii) row 2 is a **`limsup`**, row 0 holds at
  **every** `t < T`; (iii) row 2's log is `log(1/(T_*-t))`, row 0's is
  `(1/5)log(1/(T-t)) + O(1)` — a factor **5** in the constant of the asymptotic rate;
  (iv) row 2's constant depends on `||u_0||_{H^3}`, a subcritical non-invariant norm, whereas
  row 0's constant is universal and all data-dependence sits inside a **scaling-invariant**
  Reynolds number `E_0^{2/5}M^{1/5}/nu`. Row 0 degenerates as `nu -> 0` (the log diverges), so
  it does **not** recover row 2 in the inviscid limit.
* **vs row 5' (the classical baseline).** Row 0 is strictly stronger for small `T-t`: exponent
  `1` against `5/6`, ratio growing like `(T-t)^{-1/6}/log` (fitted slopes 0.131 and 0.154
  towards 1/6). Row 5' is algebraically better in `E_0` and `nu` for large `T-t`; the crossover
  depends on constants neither statement pins down.
* **vs rows 6, 7, 8.** Different quantities (velocity `Ḣ^s`, velocity `L^3`) or no statement.
  Row 8 is the field's own assessment, published Apr 2025, that NS vorticity lower bounds do
  not exist — for the *unconditional* growth question, which is strictly harder than row 0's
  conditional one.

---

## 8. What is established, and what is not

**Established here.** (a), (b), (c) and the `log log` accumulation bound §4, given theorem (i);
their mutual equivalences and scaling invariance, verified symbolically; the classical
comparator's exponent `5/6`; the source text of every citation.

**Not established here.** Theorem (i) itself (Astra pass 8, separately referee-confirmed;
I did not re-verify it). Sharpness: nothing here says the logarithm is necessary — that is
conjecture (ii), and the campaign's own numerics do not exclude a log-free clock. The status of
BFG Theorem 8/10. And the corollary is vacuous unless a blow-up occurs.

**One-line honest claim.** *Conditional on Astra pass 8, a Navier–Stokes solution that blows up
at `T` satisfies `||omega(t)||_inf >= c_1/((T-t)(1+log_+(E_0^{2/5} nu^{-1}(T-t)^{-1/5})))` for
every `t < T`, with `c_1` universal — a pointwise-in-time, scaling-invariant, energy-and-viscosity
bound for which I found no published counterpart, the nearest being an Euler `limsup` result
(Ingimarson–Kukavica 2026) and a stronger but sketch-proved NS claim (Bradshaw–Farhat–Grujić 2019).*

---

## 9. Files

    NOTE.md                                  this note
    SHA256SUMS                               recomputed, never typed
    scripts/s1_corollary_algebra.py          blocks (A)-(E); log s1_log.txt
    scripts/s2_leray_nu_and_invariance.py    viscosity power + invariance; log s2_log.txt
    scripts/pagefind.sh                      locates the PDF page containing a string
    pdf/  txt/                               the four primary sources, as downloaded and extracted:
      ik-2603.17431.pdf     Ingimarson-Kukavica v2
      bfg-1704.05546v4.pdf  Bradshaw-Farhat-Grujic v4
      luo-2504.08288.pdf    Luo
      rss-1503.03063.pdf    Cortissoz-Montero (the misattributed one)

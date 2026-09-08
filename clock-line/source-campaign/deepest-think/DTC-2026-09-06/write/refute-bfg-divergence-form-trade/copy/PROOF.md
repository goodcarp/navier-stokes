# The divergence-form trade in BFG's vorticity mild formulation

**Seat** `write/bfg-divergence-form-trade`, DTC-2026-09-06.
**Task** BFG (Bradshaw–Farhat–Grujić, arXiv:1704.05546v4, ARMA 231 (2019)) prove their
Theorem 8 / Theorem 10 clock `T >= 1/(c ||omega_0||_inf)` from the mild vorticity
equation (8) with the **undifferentiated** heat kernel `G`, through a weighted-BMO
inequality whose local-average subtraction they drop. Because `div u = div omega = 0`,
both nonlinear terms are in divergence form; moving the derivative onto the kernel
replaces `G` (mean one) by `grad G` (mean zero), for which the subtraction is free.
Question: what clock does that repair actually produce, and how does it compare with the
logarithmic repair and with the Giga–Inui–Matsui / Kukavica velocity clock?

**Answer, in one line.** The divergence-form repair produces
`M_0 T_d >= c_div / Re_E` with `c_div = 3·3^{1/5} π^{11/5}/12800 = 3.622974321832e-03`
and the Reynolds exponent **exactly −1** (not −1/2, as the brief's guess had it); this is
*the same clock as* Giga–Inui–Matsui/Kukavica's `T >= c nu/||u_0||_inf^2` once that is fed
the same velocity bound, with `C_u` cancelling identically; and it is weaker than the
logarithmic clock at every `Re_E >= 1` for every `c_1 >= c_div`. Neither repair is
energy-free, and by a scaling argument (§8) **none can be**.

Every number below came out of a script in this folder that I wrote and ran
(`s1`…`s7`, re-asserted by `check_constants.py`: 50 checks pass, 50/50 single-number
mutations caught; `s7_audit.py` re-derives all 9 closed forms quoted below, 0 mismatches). `SHA256SUMS` is computed, never typed. Nothing outside this folder was
written or modified.

---

## 0. Verdict table

| item of the brief | verdict |
|---|---|
| (1) the a-priori bound is `||omega_0||_inf + C ∫ (t−s)^{−1/2}||omega||_inf||u||_inf ds` | **ESTABLISHED**, with `C = 2/sqrt(pi nu)` and the constant `2` shown sharp for the pointwise pairing (§3–§4) |
| (2) the resulting clock, exact exponent | **ESTABLISHED**: `M_0 T_d >= c_div/Re_E`, `c_div = 3·3^{1/5}π^{11/5}/12800 = 3.622974321832e-03`. Exponent of `Re_E` is **−1**. The brief's `Re_E^{-1/2}` is **wrong** (§6) |
| (3) two honest repairs, neither energy-free | **ESTABLISHED**, and strengthened: no repair through `||u||_inf` can be energy-free, because `||u||_inf` is not a function of `||omega||_inf` (§8, claim N1) |
| (4) comparison with the GIM/Kukavica velocity clock | **ESTABLISHED**: identical functional form and identical exponent; the ratio is `pi·2^{−46/5}/c_G = 5.341631355917e-03/c_G` and `C_u` cancels (§7) |
| does the divergence-form route recover BFG Thm 8/10 as stated? | **NO** — it is weaker by a factor `∝ Re_E` (§7). Both available repairs are strictly weaker than the published statement |
| does either repair contradict (S3) `T(Λ) <= c_2/log Λ`? | **NO**. Only BFG's `Re`-free statement does (§7.3) |

**Novelty, stated honestly.** The *exponent* `−1` is already on the estate record: it is
`verify-astra/refuter-A1-entropy-clock/r5_final.py` item (d) ("`T <~ nu/(M^2 ell^2)
= 1/(M Re)`") and `r3_closure.py` line 47. What is new here is (i) the mean-zero
(divergence-form) route specifically — BFG's own kernel, not a heat-route paraphrase —
with the BMO step removed entirely rather than repaired; (ii) every constant explicit and
derived (`C(nu)`, `Q`, `C_u`, `c_div`), where the literature's are unquantified; (iii) the
exact identification with the GIM/Kukavica clock, including the cancellation of `C_u`; and
(iv) the impossibility claim N1. Grade: **receipts**, not Solid; nothing here touches
FL-000.

---

## 1. Setting, notation, standing hypotheses

3D incompressible Navier–Stokes on `R^3`, viscosity `nu > 0`:

```
    d_t u + (u.grad)u = nu Δu − grad p ,   div u = 0 ,
    omega := curl u ,
    d_t omega + (u.grad)omega = nu Δ omega + (omega.grad)u .            (V)
```

Notation, fixed once:

* `M_0 := ||omega_0||_inf`, `M(t) := ||omega(t)||_inf`;
* `E := ||u||_{L^2}^2` — **the estate convention, no 1/2**
  (`sharp/exact-first-order/NOTE.md` §3, line "`E := ||u||_2^2 (no 1/2)`");
* `Re_E := E^{2/5} M^{1/5} / nu`, likewise that seat's definition (dimensionless: check in
  §6.3);
* `T_d` := the first time `||omega(t)||_inf = (3/2)M_0` (the estate's doubling threshold,
  `RETURN_ADDENDUM_8_CLOCK_LINE.md`, "Object");
* `G_nu(x,t) := (4 pi nu t)^{-3/2} exp(−|x|^2/(4 nu t))`.

**(H1) Regularity.** `u` is a smooth solution of NS on `R^3 × [0,T_max)` with
`u(0) = u_0`, and on each `[0,T] ⊂ [0,T_max)`
`sup_t (||u||_inf + ||omega||_inf + ||grad u||_inf) < ∞`, with enough spatial decay of `u`,
`omega`, `grad u` that (a) the Duhamel representation (§3) holds pointwise and (b) the
integration by parts moving `d_i` from the nonlinearity onto `G_nu` produces no boundary
term. (Schwartz-class data, or `u_0 ∈ H^s`, `s > 5/2`, with compactly supported `omega_0`,
suffice; the campaign's mollified plateau is in this class.)

**(H2) Finite energy.** `u_0 ∈ L^2(R^3)`, `E_0 := ||u_0||_2^2 < ∞`, and the solution obeys
the energy inequality `||u(t)||_2^2 <= E_0` for `t >= 0`.

**(H2) is strictly stronger than BFG's hypothesis.** BFG assume `omega_0 ∈ L^2 ∩ L^inf`
(their Thm 8, p.8), which gives `u ∈ Ḣ^1 ⊂ L^6` but **not** `u ∈ L^2`. So the
divergence-form repair is not a drop-in replacement for their argument: it needs finite
energy, which their statement deliberately does not assume ("no quantitative dependence on
`||omega_0||_2`", p.8). §8 shows this is unavoidable, not an artefact.

All estimates below are *a priori* estimates for a solution satisfying (H1)–(H2). To turn
them into a local existence theorem one runs the identical chain on the Picard iterates
`omega^{(n)}`, exactly as Kukavica (Prop. 3.2) and BFG (their (9)) do; nothing in the
constants changes.

---

## 2. LEMMA D (the divergence form) — exact

**Statement.** If `div u = 0` and `div omega = 0`, then with
`A_{ij} := u_i omega_j − omega_i u_j` (an **antisymmetric** matrix field),

```
    u_i d_i omega_j − omega_i d_i u_j  =  d_i A_{ij} ,   j = 1,2,3.        (D)
```

*Proof.* `d_i(u_i omega_j) = (d_i u_i)omega_j + u_i d_i omega_j = u_i d_i omega_j` and
`d_i(omega_i u_j) = (d_i omega_i)u_j + omega_i d_i u_j = omega_i d_i u_j`. ∎

Verified in `s2_algebra.py`: antisymmetry residual `0` (symbolic), and (D) checked as an
identity of explicit polynomial divergence-free fields, residual `[0,0,0]` in all three
components.

**LEMMA D′ (the pointwise constant).** For every `v ∈ R^3`,
`A^T v = −A v = −[(v·omega)u − (v·u)omega]`, and

```
    Q := sup { |(n·u)omega − (n·omega)u| / (|u||omega|) : |n| = 1, u,omega ≠ 0 }  =  1 ,
```

**exactly**, and the supremum is attained (`u = e_1`, `omega = e_2`, `n = (e_1+e_2)/√2`).

*Proof.* Put `a = u/|u|`, `b = omega/|omega|`, `c = a·b`, `alpha = n·a`, `beta = n·b`. Then

```
    |(n·u)omega − (n·omega)u|^2 / (|u|^2|omega|^2) = alpha^2 + beta^2 − 2c·alpha·beta .
```

For a **unit** `n` the triple `(a,b,n)` has positive-semidefinite Gram matrix
`[[1,c,alpha],[c,1,beta],[alpha,beta,1]]`, whose determinant is
`(1−c^2) − (alpha^2 − 2c·alpha·beta + beta^2)`; and in dimension `>= 3` every
`(alpha,beta)` with non-negative determinant is realised by some unit `n`. So the
admissible set is exactly `{alpha^2 − 2c·alpha·beta + beta^2 <= 1 − c^2}`, and the
objective **is** that same quadratic form. Hence the supremum over `n` is `1 − c^2 <= 1`,
with equality at `c = 0` on the boundary of the ellipse. (Degenerate case `c = ±1`, i.e.
`omega` parallel to `u`: then `A = 0` and both sides are `0`, consistent with `1 − c^2 = 0`.) ∎

Verified in `s2_algebra.py`: the Gram determinant identity has symbolic residual `0`;
4 000 000 random `(u,omega,n)` give max `0.999999980775060 <= 1`; the stated configuration
attains `1.000000000000000`.

> The triangle inequality would give `2` here. The antisymmetry is worth a factor `2` in
> `||u||_inf||omega||_inf`, hence a factor `4` in the clock: `c_div` would be
> `9.057435804579e-04` instead of `3.622974321832e-03` (`s4_log.txt`).

---

## 3. LEMMA K (kernel facts) — exact, `s1_kernel.py`

All four are derived symbolically and confirmed by independent quadrature.

```
 K1   ∫_{R^3} G_nu(x,t) dx        = 1
 K2   ∫_{R^3} d_1 G_nu(x,t) dx    = 0                    (grad G is MEAN ZERO)
 K3   ∫_{R^3} |grad G_nu(x,t)| dx = 2 / sqrt(pi nu t)     EXACTLY
 K4   || G(·,tau) ||_{L^2(R^3)}   = (8 pi tau)^{-3/4}
 K5   ∫_0^t (t−s)^{-1/2} ds       = 2 sqrt(t)
```

K3 route: `|grad G_nu| = (|x|/(2 nu t)) G_nu`, so the integral is `E|X|/(2 nu t)` for
`X ~ N(0, 2 nu t · I_3)`; `E|X| = 2σ√(2/π) = 4 sqrt(nu t/pi)`. Independent radial
quadrature at `(nu,t) = (1,1), (1,0.01), (7,3), (0.001,5)` returns
`1.12837916709551256`, `11.2837916709551251`, `0.246232521229829082`,
`15.9576912160573077` against the closed form, relative error `0` to `1e-31`
(`s1_log.txt`). This reproduces `[K1]` of `gaps/kukavica-at-source` (`2/sqrt(pi tau)` at
`nu = 1`) by a different route.

Two further constants, recomputed here as cross-checks of that seat's `[K2]`,`[K3]`
(not used in the proof):

```
 K6   sup_x |grad G(x,1)|(|x|+1)^4 = 0.9231098781   at |x| = 2.810821027   (their 0.9231098779)
 K7   sup_x G(x,1)(|x|+1)^4        = 0.7109819927   at |x| = 2.372281323   (their 0.7109819926)
 K8   ∫_0^t∫ (|z|+√(t−s))^{-4} dz ds       = (8 pi/3)√t = 8.377580409573·√t
      ∫_0^t∫ √(t−s)(|z|+√(t−s))^{-4} dz ds = (4 pi/3) t = 4.188790204786·t
```

K8 is the exponent fingerprint: the **mean-zero** kernel produces `t^{1/2}`, the
**mean-one** kernel `t^1`. That single difference is what turns an inverse-linear clock in
`M_0` into an inverse-quadratic one in the datum, and it is the reason the repair below
lands on a power of `Re` (§6.4).

---

## 4. PROPOSITION 1 — the a-priori bound (item (1) of the brief)

**Statement.** Under (H1), for all `t ∈ [0,T_max)`,

```
   ||omega(t)||_inf  <=  ||omega_0||_inf
                       + (2/sqrt(pi nu)) ∫_0^t (t−s)^{-1/2} ||omega(s)||_inf ||u(s)||_inf ds .
                                                                                   (P1)
```

*Proof.* Write (V) as `d_t omega_j − nu Δ omega_j = −d_i A_{ij}` using Lemma D, and apply
the heat semigroup:

```
   omega_j(x,t) = (G_nu(·,t) * omega_{0j})(x) − ∫_0^t ∫ G_nu(x−y,t−s) d_{y_i}A_{ij}(y,s) dy ds .
```

Integrating by parts in `y` (no boundary term, by (H1)) and using
`d_{y_i}[G_nu(x−y,τ)] = −(d_iG_nu)(x−y,τ)`,

```
   omega_j(x,t) = (G_nu(·,t) * omega_{0j})(x) − ∫_0^t ∫ (d_iG_nu)(x−y,t−s) A_{ij}(y,s) dy ds .
                                                                                   (MILD)
```

This is BFG's (8) with `G` replaced by `grad G` and the integrand `u ⊗ omega` replaced by
the antisymmetric `A = u⊗omega − omega⊗u`; **both** nonlinear terms of (8) are absorbed
into the single term, and no BMO inequality — subtracted or not — occurs anywhere.

First term: `||G_nu(·,t)*omega_0||_inf <= ||G_nu(·,t)||_{L^1}||omega_0||_inf =
||omega_0||_inf` by K1 and `G_nu >= 0`; constant exactly `1`.

Second term: the `R^3`-vector with components `∑_i (d_iG_nu)A_{ij}` is `A^T v` with
`v = (grad G_nu)(x−y,t−s)`, so by Lemma D′ its magnitude is
`<= Q |grad G_nu(x−y,t−s)| |u(y,s)||omega(y,s)| = |grad G_nu| |u||omega|`. Hence

```
   |Duhamel| <= ∫_0^t ||u(s)||_inf ||omega(s)||_inf ∫|grad G_nu(z,t−s)|dz ds
             =  (2/sqrt(pi nu)) ∫_0^t (t−s)^{-1/2}||omega(s)||_inf||u(s)||_inf ds
```

by K3. ∎

So the brief's `C` is `C = 2 Q/sqrt(pi nu) = 2/sqrt(pi nu)`, and the numeral `2` is
**sharp for the pointwise pairing** (Lemma D′ attains `Q = 1`; K3 is an identity, not a
majorant). No claim is made that `2/sqrt(pi nu)` is the norm of the Duhamel operator
itself.

**Where the local average went.** BFG need the weighted-BMO inequality because their
kernel has mean one and their second factor is `|d_i u_j|`, which is not in `L^inf`. In
(MILD) the second factor is `|u||omega|`, and both are in `L^inf` by hypothesis; the
kernel is mean zero and integrable. There is no step at which a local average is dropped —
the objection of `sharp/literature-short-time/NOTE.md` §"the step where the energy/log was
dropped" simply does not arise. The cost is that `||u||_inf` has appeared. That is the
trade.

---

## 5. LEMMA U — the velocity bound, with its exact constant (`s3_velocity.py`)

**Statement.** Let `u` be smooth, `div u = 0`, `u ∈ L^2(R^3)`, `omega = curl u ∈ L^inf`.
With `E = ||u||_2^2` and `M = ||omega||_inf`,

```
   ||u||_inf  <=  C_u · E^{1/5} M^{3/5} ,
   C_u = ((3/2)^{2/5}+(2/3)^{3/5}) · (4/√pi)^{3/5} · (8 pi)^{-3/10}
       = 5·2^{9/10}·3^{2/5} / (6 pi^{3/5})
       = 1.214239420800536 .
```

*Proof (the heat split at `ℓ = √tau`).* For any `tau > 0` write
`u = e^{tau Δ}u + (u − e^{tau Δ}u)` (an auxiliary, non-dynamical mollification; the
`Δ` here carries no `nu`).

*Far/low-frequency half.* By Cauchy–Schwarz and K4,
`||e^{tau Δ}u||_inf <= ||G(·,tau)||_{L^2}||u||_{L^2} = (8 pi tau)^{-3/4} E^{1/2}`.

*Near/high-frequency half.* Since `div u = 0`, `Δu = grad(div u) − curl curl u = −curl omega`,
so `u − e^{tau Δ}u = −∫_0^tau Δ e^{σΔ}u dσ = ∫_0^tau curl(G_σ * omega) dσ`, and by K3 (at
`nu = 1`)

```
   ||u − e^{tau Δ}u||_inf <= M ∫_0^tau (2/√(pi σ)) dσ = (4/√pi) M √tau .
```

(The limit `e^{σΔ}u → u` in `L^inf` as `σ → 0` uses uniform continuity of `u`, which (H1)
gives; it is used qualitatively only.)

*Optimisation.* With `A = (4/√pi)M`, `B = (8 pi)^{-3/4}E^{1/2}`, minimise
`f(tau) = A tau^{1/2} + B tau^{-3/4}` over `tau > 0`:

```
   tau_* = (3B/(2A))^{4/5} ,     f(tau_*) = gamma · A^{3/5} B^{2/5} ,
   gamma := (3/2)^{2/5} + (2/3)^{3/5} = 5·2^{-2/5}3^{-3/5} = 1.960131704207789 .
```

Substituting `A`, `B` gives `C_u` as stated. ∎

Both the symbolic optimisation and an independent golden-section minimisation at
`M = E = 1` return `C_u = 1.214239420800536` (agreeing to `1e-25`), at
`tau_* = 0.1042175951` (`s3_log.txt`). The scaling residual is `0`: `E^{1/5}M^{3/5}` has
exactly the dimension of a velocity (`s3_velocity.py`).

**The enstrophy variant, and why it is useless here.** The same split against
`||omega||_{L^2}` (using `||e^{tau Δ}·BS||_{L^2} = 2^{3/4}/(4 pi^{3/4} tau^{1/4})`) gives

```
   ||u||_inf <= C_w M^{1/3} ||omega||_2^{2/3} ,   C_w = 3·2^{1/6}/(2 pi^{2/3}) = 0.784927737928 ,
```

which matches BFG's hypothesis class (`omega_0 ∈ L^2 ∩ L^inf`) better than Lemma U does.
(`C_w` here is computed by identifying the Biot–Savart symbol with the scalar `1/|xi|`; the
**exponents** `1/3, 2/3` are exact and scaling-forced, the numeral is indicative — I did not
verify the matrix operator norm, and nothing below depends on it.)
It cannot be used in a clock: enstrophy is **not** non-increasing along NS, so
`||omega(s)||_2` does not propagate over the window, whereas `E(s) <= E_0` does by (H2).
This is the precise sense in which BFG's remark that there is "no quantitative dependence
on `||omega_0||_2`" is correct and unhelpful: `||omega_0||_2` is the wrong `L^2` quantity.

---

## 6. THEOREM C — the divergence-form clock (item (2))

**Statement.** Assume (H1)–(H2), `M_0 = ||omega_0||_inf > 0`, `E_0 = ||u_0||_2^2 < ∞`,
`Re_E = E_0^{2/5}M_0^{1/5}/nu`. Then

```
   ||omega(t)||_inf <= (3/2) M_0    for all   0 <= t <= t_* ,
   t_* = 3·3^{1/5} pi^{11/5} · nu / (12800 · E_0^{2/5} M_0^{6/5}) ,
```

equivalently, in scaling-invariant form,

```
   M_0 · T_d  >=  M_0 t_*  =  c_div / Re_E ,
   c_div = 3·3^{1/5} pi^{11/5}/12800 = pi/(2^{46/5} C_u^2) = 3.622974321832e-03 .
                                                                             (C)
```

### 6.1 Proof

Let `T* := sup{ t ∈ [0,T_max) : ||omega(s)||_inf <= 2 M_0 for all s <= t }`; `T* > 0` by
continuity. On `[0,T*]`, Lemma U with (H2) (`E(s) <= E_0`) gives

```
   ||u(s)||_inf <= C_u E_0^{1/5}(2M_0)^{3/5} .
```

Feed this and `||omega(s)||_inf <= 2M_0` into (P1) and use K5:

```
   ||omega(t)||_inf <= M_0 + (2/√(pi nu))·2√t·(2M_0)·C_u E_0^{1/5}(2M_0)^{3/5}
                     = M_0 + 2^{18/5} C_u √t M_0^{8/5} E_0^{1/5} / √(pi nu)
                     = M_0 + (40√2·3^{2/5}/(3 pi^{11/10})) · E_0^{1/5}M_0^{8/5}√t/√nu .
```

Setting the Duhamel term `= M_0/2` and solving for `t`:

```
   √t_* = √(pi nu) / (2^{23/5} C_u M_0^{3/5} E_0^{1/5}) ,
    t_* = pi nu / (2^{46/5} C_u^2 M_0^{6/5} E_0^{2/5}) = 3·3^{1/5}pi^{11/5}nu/(12800 E_0^{2/5}M_0^{6/5}) .
```

For `t <= t_*` we then have `||omega(t)||_inf <= (3/2)M_0 < 2M_0`, so `T*` cannot be
attained before `t_*` (continuity of `t ↦ ||omega(t)||_inf` under (H1)), i.e. the bootstrap
closes on `[0, t_* ∧ T_max)`. In particular the first time `||omega||_inf` reaches
`(3/2)M_0` satisfies `T_d >= t_*`. ∎

Numerically checked end to end in `check_constants.py`: at `M_0 = 3.7`, `E_0 = 0.41`,
`nu = 0.013` (`Re_E = 69.9536750335`) the formula gives `t_* = 1.3997581237509200e-05`,
at which the Duhamel majorant equals `1.8500000000000001 = M_0/2` to relative error `0`,
and `M_0 Re_E t_* = 3.6229743218317728e-03 = c_div`.

### 6.2 The exponent, extracted rather than guessed

`t_*` has exponents `M_0^{-6/5} E_0^{-2/5} nu^{+1}` (read off by logarithmic
differentiation in `s4_clock.py`). Ask for `t_* = c' M_0^{alpha} Re_E^{beta}` with
`Re_E = E^{2/5}M^{1/5}nu^{-1}`. Matching the three exponents gives a system of **three**
equations in **two** unknowns:

```
   nu :   +1   = −beta
   E  :  −2/5  = (2/5) beta
   M  :  −6/5  = alpha + (1/5) beta
```

which is overdetermined and **consistent**, with the unique solution
`alpha = −1, beta = −1` (sympy `solve` returns `[{alpha: -1, beta: -1}]`). The
consistency of the `nu` and `E` rows is the real check: a clock of the shape
`c/(M Re_E^{1/2})` would need `beta = −1` from `E` and `beta = −2` from `nu` simultaneously.

> **The brief's guess `tau >= c/(M·Re_E^{1/2})` is therefore wrong.** The exponent is
> exactly `1`. The brief's own sketch line confirms it once carried through:
> `t^{1/2}·M·(E^{1/5}M^{3/5})/√nu <= 1/2` gives `t <= nu/(4E^{2/5}M^{6/5})`, i.e.
> `M t <= 1/(4 Re_E)`. Stronger: `M t_* = nu E^{-2/5}M^{-1/5} = Re_E^{-1}` **identically**,
> so a law of the form `M t_* = c·Re_E^{-1/2}` is not merely false but dimensionally
> impossible — `M t_*/Re_E^{-1/2} = √nu/(E^{1/5}M^{1/10})` is not a pure number
> (`s7_audit.py`).

### 6.3 `Re_E` is dimensionless (consistency of the estate's definition)

`[E] = L^5 T^{-2}`, `[M] = T^{-1}`, `[nu] = L^2 T^{-1}`, so
`[E^{2/5}M^{1/5}/nu] = (L^2 T^{-4/5})(T^{-1/5})(L^{-2}T) = 1`. ✔

### 6.4 Where the exponent comes from

It is the kernel's, not the argument's. K8: the mean-zero kernel gives `t^{1/2}` in the
Duhamel term, hence a bound quadratic in the data (`√t · X^2 <= X` ⟹ `t <~ X^{-2}`), and
the datum `X` must carry the velocity. The mean-one kernel gives `t^1`, hence a bound
linear in the data (`t·M^2 <= M` ⟹ `t <~ M^{-1}`) — BFG's Theorem 8. Kukavica's original
carries `∂_jG` throughout and lands on `T >= 1/(C||u_0||_inf^2)`; BFG's modification
carries `G` and lands on `T >= 1/(c||omega_0||_inf)`. Moving the derivative onto the
kernel is not presentational: it is exactly the trade between those two clocks.
(The `t^{1/2}`/`t` dichotomy is `[K3]` of `gaps/kukavica-at-source`; reproduced here as K8.)

---

## 7. The three clocks compared (items (3) and (4))

All windows written in the scaling-invariant dimensionless form `M_0 · T`.

| | clock | `M_0 T >=` | `Re_E` dependence | status |
|---|---|---|---|---|
| **(a)** | log clock — retain the local average (Astra, (S1)) | `c_1/(1+log_+ Re_E)` | logarithmic | PROVED (Astra pass 8, campaign-internal referee; `RETURN_ADDENDUM_8` (S1)). `c_1` not numerically pinned |
| **(b)** | **divergence form — this note** | `3.622974321832e-03 / Re_E` | **power, exponent −1** | proved here (§6), constants explicit |
| **(c)** | GIM / Kukavica velocity clock, fed Lemma U | `(c_G/C_u^2)/Re_E = 0.678252406508·c_G/Re_E` | power, exponent −1 | Kukavica Thm 2.2 / Prop 3.2 `T >= 1/(C||u_0||_inf^2)`, `C = 1/c_G` unquantified |
| **(d)** | BFG Thm 8/10 **as stated** | `1/c`, `c` absolute | **none** | published; its displayed proof is a sketch with the local-average step dropped |

### 7.1 (b) and (c) are the same clock

Kukavica's `T >= 1/(C ||u_0||_inf^2)` is stated at `nu = 1`; restoring `nu` by the NS
scaling `u ↦ λ u(λx, λ^2 t)` (which fixes `nu`) the invariant form is
`T >= c_G nu/||u_0||_inf^2`. Feeding the *same* Lemma U bound,

```
   T_GIM >= c_G nu / (C_u E_0^{1/5}M_0^{3/5})^2 = (c_G/C_u^2)·nu·M_0^{-6/5}E_0^{-2/5}
          = [3·6^{1/5}pi^{6/5}c_G/25] / (M_0 Re_E) .
```

Same functional form, same exponent. The ratio of the two constants is

```
   c_div / c_GIM = pi·2^{-46/5} / c_G = 5.341631355917e-03 / c_G ,
```

and **`C_u` cancels identically** (`s4_clock.py`, symbolic residual `0`). That is the
sharpest way to say it: the divergence-form repair of BFG's vorticity argument *is* the
GIM/Kukavica velocity clock transcribed into vorticity variables. It contributes an
explicit numeral, `pi·2^{-46/5} = 5.341631355917e-03`, where the source leaves `C`
unquantified — and nothing else.

### 7.2 (a) beats (b) everywhere that matters

`s4_clock.py`:

| `Re_E` | (a)`/c_1` | (b) | (a)/(b) |
|---|---|---|---|
| `1e0`  | `1.000000e+00` | `3.622974e-03` | `2.760163e+02 · c_1` |
| `1e2`  | `1.784067e-01` | `3.622974e-05` | `4.924316e+03 · c_1` |
| `1e5`  | `7.991736e-02` | `3.622974e-08` | `2.205850e+06 · c_1` |
| `1e10` | `4.162183e-02` | `3.622974e-13` | `1.148831e+11 · c_1` |
| `1e20` | `2.125322e-02` | `3.622974e-23` | `5.866234e+20 · c_1` |
| `1e40` | `1.074075e-02` | `3.622974e-43` | `2.964621e+40 · c_1` |

Crossover: `(a) = (b)` iff `Re_E/(1+log Re_E) = c_div/c_1`. Since
`Re_E/(1+log Re_E) >= 1` for `Re_E >= 1`, the log clock is the larger of the two at
**every** `Re_E >= 1` as soon as `c_1 >= c_div = 3.62e-03`. Only for a very small `c_1`
does the power clock win in a bounded window: `c_1 = 1e-3` ⟹ crossover at `Re_E = 12.883`;
`c_1 = 1e-6` ⟹ `Re_E = 42209.1` (`s4_log.txt`).

### 7.3 What each does to the campaign's own statements

* **Neither (a) nor (b) contradicts (S3)** `T(Λ) <= c_2/log Λ`: a lower bound `c/Re` and an
  upper bound `c_2/log Re` never cross for large `Re`.
* **(d) does.** BFG's `M_0 T >= 1/c` with `c` absolute is incompatible with
  `T(Λ) → 0`. The stake recorded in `RETURN_ADDENDUM_8` is unchanged by this note.
* **The divergence-form route does not recover (d).** (b) is weaker than (d) by a factor
  `∝ Re_E`. So the honest summary is: *both* available repairs of BFG's displayed argument
  — retaining the local average, or moving to divergence form — produce statements
  strictly weaker than Theorem 8/10 as published, one by a logarithm and one by a power.
  The published statement is supported by neither.

---

## 8. Neither repair is energy-free — and none can be (`s6_no_energy_free.py`)

**CLAIM N1.** There is no function `F : [0,∞) → [0,∞)` with
`||u||_inf <= F(||omega||_inf)` for all smooth divergence-free finite-energy `u` on `R^3`.

*Proof.* Fix such a `u` with `omega = curl u ≢ 0` and set `u_λ(x) := λ^{-1}u(λx)` for
`λ > 0`. Then `div u_λ = 0` and `curl u_λ(x) = omega(λx)`, so `||curl u_λ||_inf` is
**constant** in `λ` while `||u_λ||_inf = λ^{-1}||u||_inf → ∞` as `λ → 0`. ∎

Both identities are verified symbolically on an explicit Gaussian-localised field
(`div u_λ = 0`; `curl u_λ − omega(λ·) = [0,0,0]`), and the ratio
`||u_λ||_inf/||omega_λ||_inf` is exhibited at `λ = 1, 0.5, 0.1, 0.01` as
`0.517409, 1.03482, 5.17409, 51.7409` (`s6_log.txt`).

Under that same family (`||omega||_inf` fixed) the candidate control quantities scale as

```
   ||u||_inf ~ λ^{-1} ,  E = ||u||_2^2 ~ λ^{-5} ,  enstrophy ~ λ^{-3} ,  ||omega||_inf ~ λ^0 .
```

**Consequence.** Proposition 1 pairs `|omega|` against `|u|`; by N1 the `||u||_inf` factor
cannot be eliminated in favour of `||omega||_inf`. Some second quantity must be paid.
Of the candidates only the energy is non-increasing along NS, so only the energy
propagates over the window — Lemma U is not one choice among many. Therefore:

> **The divergence-form repair is energy-dependent by necessity.** It removes the
> local-average defect completely (there is no BMO step left) and pays for it with an
> `Re_E`-dependence that no rearrangement of the same estimate can remove.

The logarithmic repair pays the same debt in a different currency: there the dropped term
is `|avg_{B(x,√(t−s))} grad u|`, and the standard Biot–Savart splitting bounds it by
`||omega||_inf log(L/ρ) + (large-scale term controlled by the energy)` with `L` the energy
scale — i.e. `log Re` (`sharp/literature-short-time/NOTE.md` §"the step where the
energy/log was dropped"; that diagnosis is that seat's, quoted, not re-derived here).
Neither currency is free. **This is the trade of the title, in its exact form:
`log Re_E` (mean-one kernel, local average retained) versus `Re_E^{1}` (mean-zero kernel,
no BMO at all).**

---

## 9. What the two clocks say for the campaign's own datum (`s5_datum.py`)

Quoted from `sharp/exact-first-order/NOTE.md` §3 (`c6_shell_exact.py`), **not re-derived
here**: for the bang-bang shell, `E_shell = 0.172403978 M^2 R^5` and
`log Re_E = 2L − 0.7031660` with `L = log(R/ρ_0)`; and the frozen-frame doubling window is
`M t_d = 2 ln2/L` (since `a_inner = (M/2)L`, `t_d = ln2/a_inner`). Using the shell rather
than the `δ`-tapered value is the choice most favourable to the power clock (the taper adds
`2log(1/δ)` to `log Re_E`, per `gaps/refute-gap-T-lipschitz`), so the table below
*understates* how bad the power clock is.

| `L` | `log Re_E` | `Re_E` | (b) power window `M T` | (a) log window `M T /c_1` | model `M t_d` |
|---|---|---|---|---|---|
| 8.3178 | 15.932434 | `8.305546e+06` | `4.362115e-10` | `5.905825e-02` | 0.166666 |
| 20 | 39.296834 | `1.165194e+17` | `3.109332e-20` | `2.481585e-02` | 0.069315 |
| 50 | 99.296834 | `1.330660e+43` | `2.722690e-46` | `9.970404e-03` | 0.027726 |
| 100 | 199.296834 | `3.576970e+86` | `1.012861e-89` | `4.992590e-03` | 0.013863 |

Ratio of the model window to each proved lower bound:

| `L` | model / (b) | model / (a) |
|---|---|---|
| 8.3178 | `3.820761e+08` | `2.822061 / c_1` |
| 20 | `2.229248e+18` | `2.793164 / c_1` |
| 50 | `1.018327e+44` | `2.780819 / c_1` |
| 100 | `1.368691e+87` | `2.776704 / c_1` |

Two readings, both worth recording.

1. **The power clock is vacuous for this datum.** At `L = 20` it certifies a window `18`
   orders of magnitude below the one the frozen-frame model predicts, and the gap grows
   like `e^{2L}`. It is a correct theorem and a useless instrument here.
2. **The log clock is essentially sharp for this datum**, and the ratio is not merely
   bounded — it converges. `model/(a) → 4 ln 2 = 2.7725887222`, which is exactly the
   constant `c_1 = 4 ln 2` computed independently in `sharp/exact-first-order` §3 for the
   shell (`c1 = 4 ln2 − (2 ln2)(0.7031660)/L`, approached from below). The two numbers
   agree because both are `(M t_d)·log Re_E` in the limit; that they agree is a
   consistency check on this note against that seat, and it is the quantitative version of
   "(a) is the right shape and (b) is not".

---

## 10. Honest limits of this note

1. **Nothing here is a new theorem about NS.** (b) is a re-derivation, with constants, of
   the standard mild-theory baseline; the estate already had the exponent
   (`verify-astra/refuter-A1-entropy-clock`, `r5_final.py` item (d), `r3_closure.py` l.47).
   The contribution is exactness and the identification with GIM/Kukavica, plus N1.
2. **`c_1` is not pinned.** (S1) is proved with an unquantified `c_1`; every comparison
   with (a) above is therefore reported as a multiple of `c_1`, never as a number.
   Pinning `c_1` is the obvious next unit of work, and it is what would let §7.2 be stated
   without a free parameter.
3. **`c_G` is not pinned either.** Kukavica's `C` (Thm 2.2, Prop 3.2) is not made explicit
   in the source, so §7.1's ratio carries `c_G`. I did **not** obtain Kukavica or GIM at
   source myself: the statements and page references are taken from
   `gaps/kukavica-at-source/NOTE.md` §3.2, which read them by eye from the author's
   accepted manuscript, with that seat's own caveat that it is not the Elsevier version of
   record.
4. **Hypotheses are strictly stronger than BFG's** ((H2), finite energy). A repair
   operating inside BFG's stated hypothesis class `omega_0 ∈ L^2 ∩ L^inf` is not
   constructed here, and by N1 + §5 (enstrophy does not propagate) I do not believe the
   divergence-form route supplies one.
5. **`c_div` is not claimed optimal.** The bootstrap cap `2M_0` and the target `M_0/2` are
   convenient, not optimised; `Q = 1` and K3 are sharp, but the cap costs a factor
   `2^{18/5}` that a variable-cap Gronwall argument would improve. The *exponent* is
   unaffected, and it is the exponent the brief asks for.
6. **BFG Theorem 8/10 is not adjudicated here.** This note establishes what the
   divergence-form route yields; it does not establish that no other route yields (d). The
   campaign's verdict stays where `RETURN_ADDENDUM_8` left it: (b), unsupported by the
   displayed proof, uncontradicted by any published NS result.

---

## 11. Files and reproduction

| file | what it is |
|---|---|
| `s1_kernel.py` / `s1_log.txt` / `s1_results.json` | Lemma K: K1–K8, symbolic + quadrature |
| `s2_algebra.py` / `s2_log.txt` / `s2_results.json` | Lemma D and D′: divergence form, Gram argument, `Q = 1`, 4e6-sample search |
| `s3_velocity.py` / `s3_log.txt` / `s3_results.json` | Lemma U: `gamma`, `C_u` (two routes), `tau_*`, the enstrophy variant `C_w` |
| `s4_clock.py` / `s4_log.txt` / `s4_results.json` | Theorem C, the overdetermined exponent solve, `c_div`, the GIM ratio, the clock table, crossovers |
| `s5_datum.py` / `s5_log.txt` / `s5_results.json` | the campaign datum: both windows at `L = 8.3178, 20, 50, 100` |
| `s6_no_energy_free.py` / `s6_log.txt` / `s6_results.json` | claim N1 and the scaling exponents |
| `s7_audit.py` / `s7_log.txt` / `s7_results.json` | adversarial re-check of every closed form and numeral quoted in this file (9 identities, 0 mismatches) |
| `check_constants.py` / `check_log.txt` | the gate: 50 independent recomputations + a 50/50 mutation test + an end-to-end numeric closure of the bootstrap |
| `SHA256SUMS` | computed over all of the above |

Reproduce with `python3 s1_kernel.py; python3 s2_algebra.py; python3 s3_velocity.py;
python3 s4_clock.py; python3 s5_datum.py; python3 s6_no_energy_free.py; python3 s7_audit.py;
python3 check_constants.py` (sympy 1.14.0, mpmath 1.3.0, numpy 1.26.4). `check_constants.py`
exits `0` only if all 50 checks pass **and** all 50 single-number mutations are caught; the
last line of `check_log.txt` prints the real counters, not a literal.

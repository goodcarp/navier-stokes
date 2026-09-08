# Scout: the dissipation ceiling of the layered cascade

Scope: how far the Alpoge-Buckmaster (A-B) layered mechanism can be pushed if a term
`-nu |grad|^alpha` is added to the momentum equation. Sources read: the extracted text of
the two papers released 2026-09-08 (Boussinesq, 76 pp; Euler, 112 pp) and the arXiv HTML
rendering of Cordoba-Martinez-Zoroa-Zheng (CMZZ), arXiv:2407.06776. Nothing below is a
theorem. Every step is labelled *read* (taken from the papers), *derived* (my algebra),
*checked* (numerically verified in the toy), or *conjectured*.

## Assumptions

1. Damping enters the layer amplitude system as a diagonal decay at the layer's own
   frequency: rate `D_q = nu lambda_q^alpha |zeta_q|^alpha`, with `|zeta_q|` of order one.
   *Read*: Theorem 3.2(ii) of the Boussinesq paper gives `1 - C eps <= r_j <= zeta_+`.
   This is the single-harmonic (`n = 1`) model. The paper's profile (4.1),
   `F(s) = sin s + chi(s)(s - sin s)`, has higher harmonics, each damped at
   `(n lambda |zeta|)^alpha`, so one amplitude pair no longer closes the residual identity of
   Lemma 3.1 ("no sinusoidal identity was used" is what makes the undamped case exact); and
   `(-Delta)^{alpha/2}` is nonlocal, so the affine region where the next layer lives receives a
   smooth residual of size about `lambda^alpha Theta` from the oscillating part. `F` is smooth, so
   the harmonics decay fast and the damping rate is `(lambda |zeta|)^alpha` up to constants: the
   exponent conditions below are unaffected, but this is the first place "nothing else changes"
   fails, before the steering root named in section 8.
2. Nothing else changes: corrections, localization errors, the steering feedback and the
   root selection are assumed unaffected. This is the main idealization, and section 8
   argues it is exactly where the released proof would break.
3. `nu = 1` (it scales away).
4. A stage "grows" iff its realized logarithmic gain over its own time interval is at least
   the prescribed gain `L_q`.

## 1. The growth rate does not see the frequency

*Derived.* Lemma 3.1 gives `zeta' = -D^T zeta`, `Theta' = -(J zeta . G) Omega / (lambda |zeta|^2)`,
`Omega' = lambda zeta_1 Theta`. The two off-diagonal coefficients carry `1/lambda` and
`lambda`, so their product, and hence the growing eigenvalue, is

    mu = sqrt( - zeta_1 (J zeta . G) / |zeta|^2 ) ,

free of `lambda`. In the frozen model of Subsection 1.2 (`G = -A e_2`, `zeta = e(phi)`) this is
`sqrt(A) sin phi = sigma sin phi` with `sigma = |G|^{1/2}`. So the guess in the brief is right in
form but needs one correction: the rate is `sqrt(|G|)` **times the insertion angle**, and the
angle is not free. Euler equation (3.4) has the identical structure, with
`d_m = (J zeta_m) . grad Gamma_<m` and `c_m = 2 (W Gamma_<m)`; `N_m` again cancels between the
two entries. The papers' own rate symbol is `Gamma_q = sigma_{q-1} sin s_q`.

## 2. The damped stage condition

*Derived.* With `y = (Theta, sigma Omega / lambda)` and damping `D` on the vorticity slot and
`delta D` on the scalar slot (`delta = 1` if the scalar is also damped):

    y' = [[-delta D, Gamma], [Gamma, -D]] y .

* `delta = 0`: `mu_+ = (-D + sqrt(D^2 + 4 Gamma^2))/2`. This is positive for every `D`.
  Momentum dissipation alone never quenches the instability; it degrades the rate from
  `Gamma` to `Gamma^2 / D` when `D >> Gamma`.
* `delta = 1`: `mu_+ = Gamma - D`, a hard threshold at `D = Gamma`.

Stage `q` grows iff `mu_+ tau_q >= L_q`. When the scalar is damped there is a second
condition: a finished layer's gradient `A_j` decays at `nu lambda_j^alpha`, so all later
stages must fit inside time `(nu lambda_j^alpha)^{-1}`. Call this retention.

## 3. What the released schedules give as written

*Read.* Both papers use `beta = 1/8`. The amplitude is amplified from the seed
`lambda_q^{-k_q-6}` to `lambda_q^{-(1-beta)}`, so `A_q = lambda_q |zeta_q| Theta_q ~ lambda_q^beta`,
`sigma_q^2 = |G_{<q+1}| ~ lambda_q^beta` (Theorem 3.2(ii): `c_sigma lambda^{1/8} < sigma^2 < C_sigma lambda^{1/8}`),
and the growth-rate exponent is `beta/2`, the exponent of the growth scale `sigma_{q-1}`. The
vorticity amplitude is a different quantity: on the growing eigenline `Omega = (lambda r / sqrt(A_old)) Theta`
(Subsection 1.2, and the normalization `y_q = (Theta_q, sigma_{q-1} Omega_q / lambda_q)` of (3.17)), so the
peak vorticity of a growing layer is `A_q / sigma_{q-1} ~ lambda_q^{beta (1 - 1/(2Q))}`, essentially
`lambda_q^beta`; the value `~ sigma_q` is the retained post-steering vorticity of Theorem 3.2(iii)
(`[c_W, C_W] sigma_q` for `t >= t^b_{q+1}`), reached only after the pulse has returned `Omega_q` to zero.
The ceiling is governed by the rate exponent `beta/2`, not by either vorticity exponent.
Frequencies: `lambda_q = lambda_{q-1}^{Q_q}`, `Q_q = Q_* + q` with `Q_* >= 200` (Boussinesq (3.7));
the Euler paper takes `Q_i = q_0 + i` with `q_0 >= 2^20` (12.1)-(12.3). Insertion angles
`s_q = L_q sigma_{q-2}/sigma_{q-1}` (an angle, not its sine, (3.10); the Euler paper writes
`s_i = Lambda_i sigma_{i-2}/sigma_{i-1}` with `Lambda_i` its growth length, since its `L_i` is the
derivative allowance of (11.1)), so `Gamma_q = L_q sigma_{q-2} (1 + O(eps))` (Lemma 3.13 proof) and stage length
`tau_q ~ 1/sigma_{q-2}`, consistent with `|[t_q^in, t_q^hat]| < C_S / sigma_{q-2}`.

*Derived.* The light-damping condition `D_q <~ Gamma_q` then reads
`lambda_q^alpha <~ L_q lambda_{q-2}^{beta/2}`, that is

    alpha < beta / (2 Q_{q-1} Q_q)  ~  beta / (2 Q^2) .

Since both schedules hard-wire `Q_q -> infinity` (`Q_q = Q_* + q` in Boussinesq (3.7), `Q_i = q_0 + i`
in Euler Section 12), the stage-`q` ceiling `beta/(2 Q_{q-1} Q_q)` tends to zero and no fixed
`alpha > 0` survives all stages, before any appeal to `k_q`. The numbers `1.5e-6` (`Q_1 = 201`) and
`5.7e-14` (`q_0 = 2^20`) are `beta/(2Q^2)` at the earliest damped stage only. *Checked*: the toy with
`Q_q = 200 + q` gives `alpha* = 1.77e-6, 1.51e-6, 1.48e-6, 1.45e-6` at `K = 3, 4, 6, 8`, each equal to
`beta/(2 Q_{K-1} Q_K)`.
The `Q^2` (rather than `Q`) is entirely due to the choice of a very small insertion angle,
which ties stage `q`'s rate to `sigma_{q-2}` rather than `sigma_{q-1}`.

## 4. The ceiling after re-optimization

*Derived, checked numerically, in the amplitude ODE.* Three levers can be moved without
leaving the amplitude system; whether the released steering admits them is a separate question:
the insertion angle (*conjectured, requires re-choosing `Lambda_q`*: the cone tail
`sum_{j>k} s_j <= eps s_k` is not the only constraint; (3.52b) also imposes `L_j s_j <= eps_sch s_{j-1}`
and (3.52c) `2(C_phi + C_theta) L_j s_j / s_{j-1} <= eps_sch`, and the shooting family of
Subsection 3.7.6 needs `Lambda_q sin s_q <~ 1` through `|Z| <= (c_p H_2 + C eps) L sin s` and
`|d_2| <= C (L s)^2`, with `c_p Lambda_1 H_2 sin s_* <= 1/4` in (3.9). With the released
`Lambda_q = L_q`, `s_j ~ 1/j^2` sends `L_j s_j` to infinity, and the widest angle the released
control admits is `s_j ~ eps^j / prod_{i<=j} L_i`, which still gives `Gamma_q = sigma_{q-1}` up to
logarithmic factors, so the exponent below stands; but a schedule with `s_j` of order one needs a
pulse compression `Lambda_q = O(1)`, a different steering design from the released one), the stage length (only `sum tau_q < infinity`
is needed, so `tau_q ~ 1/q^2` is admissible and buys the heavy-damping regime), and `beta`
(bounded scalar forces `Theta_q = lambda_q^{-(1-beta)} -> 0`, hence `beta < 1`). Then:

* **Both fields damped** (the case for hypodissipative Navier-Stokes, where the swirl `Gamma`
  and the vorticity `omega` are both parts of `u`): `alpha < beta / (2Q)`. Supremum `1/2`,
  approached only as `beta -> 1` and `Q -> 1`; this is an amplitude-ODE supremum. Inside the
  released correction scheme the expansion parameter is `delta = Pi_q^5 lambda_q^{-7/8}`
  (Proposition 7.2), which with `Pi_q = lambda_{q-1}^4 (log lambda_q)^{10}` from (3.7) and `7/8`
  read as `1 - beta` is `lambda_q^{20/Q - (1 - beta)}` up to logs, so `delta <= 1/2` needs
  `beta < 1 - 20/Q` (`beta < 0.9` at `Q = 201`), and (7.15) needs `Q >= 9k + 42`: the joint limit
  lies outside the released budget on two counts besides the `k_q -> infinity` count of section 5.
* **Momentum only, scalar undamped, long stages** (scalar-gradient blowup only, overdamped
  regime): `alpha < beta / Q`, supremum `1`. This figure lives entirely in `D_q >> Gamma_q`
  (rate `Gamma^2/D`), where the wave vorticity is quasi-static, `Omega ~ b Theta / D`, and the
  velocity is slaved to buoyancy, an IPM-like fractional-Darcy regime. With `tau_q = T_0/q^2` the
  older layers' vorticities decay during each stage (`lambda_{q-1}^alpha tau_q >> 1`), so the
  rotation feedback `F_q` of (3.4), the vorticity lower bound (3.19) and the `limsup ||omega|| = infinity`
  conclusion of (1.5) are lost. If the
  scalar carries any dissipation, retention destroys the long-stage option and the ceiling
  falls back to `beta/(2Q)`.

*Interpretation, derived.* The factor two between `1/2` and the bounded-velocity wall at `1`
is the square root in the mechanism: the rate is `sqrt(A) <= lambda^{1/2}`, because a
bounded scalar caps `A = lambda Theta` at `lambda` and the two-field coupling then takes a
square root. A rate proportional to the full `lambda U` with `U` of order one would reach `1`.

## 5. What actually forces `Q` large, and the consequence

*Read.* Boussinesq (7.4)-(7.15): with `d = k_q` controlled derivatives and `J_q = 2 k_q + 8`
correction levels, the required order is `r_* = d + 2(2 J_q + 4) + 1 = 9d + 41`, and the
construction needs `r_* <= Q - 1`; (7.15) alone needs `Q >= 9k + 42`. The stated hypothesis
`120 d <= Q` with `Q >= 201` is not only slack for (7.15): it is also load-bearing in the gain-ratio
bound (3.52a), `100 <= L_j/L_{j-1}`, whose proof bounds the ratio below by
`(41/8) Q_j / (Q_{j-1}/120 + 41/8) >= 100`; with `k <= Q/c` in place of `Q/120` that bound needs
`c >= ~20` for large `Q` and `c >= ~39` at `Q = 201`, and `a_L = Q_1/120 + 41/8` (Section 8) uses it
too. So within the released scheme the derivative budget costs at least `Q >~ 20 k_q`, unless the
constant 100 in (3.52a) is itself renegotiated. Euler: `c_Q >= 2^18` and `Q_i = q_0 + i`,
`q_0 >= 2^20`; its (7.15) analogue is `B + 3 = 9k + 51 <= Q_m` ((10.3) with `r_* = k + 4J_m + 9`),
encoded as `9/c_Q + 51/(q_0 + 1) < 1` in (12.2), but the binding use of `c_Q` is the majorant
`mu_j = 2048 k_j + 9216 <= Q_j/4` (Section 10), i.e. `Q_j >= 8192 k_j + 36864`.
Lemma 7.4 and Theorem 8.2 deliver `k_q -> infinity` (from `Q_q = Q_* + q` and increasing
`log lambda_q`); Lemma 9.3 and Proposition 9.4(b) via (9.6) require it; and it is structurally
necessary, because uniform `C^k` control of layer `q` exists only for `k <= k_q` (Definition 7.1:
"Beyond these ranges the last hypothesis gives finiteness for each fixed layer, not the same uniform
constants"; Euler Section 11 says the same), and a finite correction count `J_q = 2k_q + 8` leaves a
residual whose `k`-th derivative gains `lambda_q^k`, so bounded `k_q` yields a `C^{k-bar}` force,
not `C^infinity`. That is the `C^infinity` force.

*Derived, and the main finding.* For the released schedules this is redundant: `Q_q -> infinity`
is hard-wired and section 3 already kills every fixed `alpha`. The chain is the operative
argument for a hypothetical bounded-`Q` re-optimization: there `k_q <= Q_q/c` (`c = 120` as fixed,
`c >= ~20` the least value the released Boussinesq proof visibly tolerates, `c ~ 2^13` for Euler)
caps `k_q`, so a `C^infinity` force forces `Q_q -> infinity` and the ceiling `beta/(2 Q_q) -> 0`.
Within this architecture a fixed `alpha > 0` and a `C^infinity` force are incompatible. A `C^k`
force costs `Q >~ c k` and buys at most `alpha <~ beta/(2ck)`, i.e. `beta/(240k)` at the released
constant and `beta/(40k)` at best for Boussinesq. This is *conjectured* as a structural statement
(it is a statement about A-B's correction scheme, not about all possible schemes), but the
arithmetic behind it is *read* and *derived*.

## 6. CMZZ's heuristic, and where `alpha_0` actually comes from

*Read* (arXiv HTML rendering, since checked against the PDF via pdftotext; Section 1.2.4 matches).
The four constraints below are CMZZ's own back-of-the-envelope Section 1.2.4, which ends with
"all we need is s(alpha) > 0, which holds if alpha is a small positive number" and the disclaimer that
"in reality there are other sources of error getting in the way". That section does not derive
`alpha_0`; `alpha_0` is introduced in Theorem 1 and fixed in Section 4. Their heuristic constraints are
`M_n^alpha <= A_{n-1}`, `A_n M_n^s / K_n <= 1`, `K_n / L_n <= 1/M_{n-1}`,
`A_n^2 L_n M_n^{s-1} / A_{n-1} <= 1`, with `A_n` the vorticity amplitude, `M_n` the frequency,
`K_n` the stretching factor, `L_n` the concentration scale, `s` the force regularity index.
The first is exactly "damping rate at most growth rate", the growth rate being the previous
layer's amplitude. *Derived*: eliminating `K_n, L_n` gives
`A_n^3 M_n^{2s-1} M_{n-1} <= A_{n-1}`; setting `A_n = M_{n+1}^alpha` (the binding value from the
first constraint) gives the three-term log-recurrence

    3 alpha log M_{n+1} + (2s - 1 - alpha) log M_n + log M_{n-1} <= 0 ,

and a geometric solution `log M_n = R^n` (so `M_{n+1} = M_n^R`, i.e. `R` is the scale-separation
exponent) exists iff `3 alpha R^2 + (2s-1-alpha) R + 1 <= 0` for some `R > 1`, whose
discriminant condition is `s <= s(alpha) := (1+alpha)/2 - sqrt(3 alpha)`.

*Checked.* The heuristic optimum at `alpha_0 = (22-8 sqrt 7)/9 = 0.092665` would be
`R = 1/sqrt(3 alpha_0) = 1.8966` and amplitude exponent `alpha_0 R = 0.17575`, but those are not
CMZZ's parameters. Section 4 of the paper takes `R = sqrt(2/(7 alpha)) = 1.7559` and
`A = N^{sqrt(2 alpha/7)}`, amplitude exponent `a = 0.16271` at `alpha_0`. In both the heuristic and
the proof `alpha = a/R` holds by construction (equality in the first constraint,
`M_n^alpha <= A_{n-1}`): it is the damping-at-most-growth constraint itself, the same constraint as
assumption 1 above, so **admissible `alpha` = (growth-rate exponent) / (scale-separation exponent)**
is a restatement of that assumption, not an independent check. The content of `alpha_0` lies in
what caps `a`. CMZZ reach 0.0927 with `a = 0.163` and separation `R = 1.756`; A-B sit at a
growth-rate exponent `beta/2 = 0.0625` and separation at least 201.

*Resolved (2026-09-08 referee pass).* The bare heuristic condition `s(alpha) > 0` gives
`5 - 2 sqrt 6 = 0.10102`, above `alpha_0 = 0.09267`. The tightening constraint is the quantitative
outer-velocity linearization remainder of CMZZ Section 4.3.4 (Lemma 11): the term
`(Omega_n * K - U_n) . grad omega_n` carries three powers of `1/L_n` and the outer layer's `C^2`
norm `M_{n-1}^2`, and after the lifetime factor the base exponent is
`N^{8 sqrt(2 alpha/7) + r + 1} L^{-3} = N^{2 sqrt(14 alpha) - 2 - 3 alpha + r + 3s}`, closed by
`r + 3s < 4s = 3 alpha + 2 - 2 sqrt(14 alpha)`. That is the defining equation of their `s`. In
heuristic exponents (`a = log_N A`, `l = log_N L`, per `R^n`) it reads `a + 1 + s + 2/R <= 3l`,
whereas the heuristic third constraint `K_n/L_n <= 1/M_{n-1}` reads `a + s + 1/R <= l`. Corrected
derivation, *checked numerically*: keep constraint 1 tight (`R = a/alpha`) and constraint 4 tight
(`l = 1 + alpha - s - 2a`, the self-interaction bound of Section 4.3.2, which is why
`L = N^{1 + alpha - s - 2 sqrt(2 alpha/7)}`), replace constraint 3 by the remainder constraint. Then
`4s <= 2 + 3 alpha - 7a - 2 alpha/a`, maximized at `a = sqrt(2 alpha/7)`, giving `R = sqrt(2/(7 alpha))`
and `4s = 2 + 3 alpha - 2 sqrt(14 alpha)`: exactly the Section 4 choices, with `s > 0` iff
`alpha < (22 - 8 sqrt 7)/9` (numeric root 0.0926655). The heuristic third constraint with the same
optimization gives `a = sqrt(alpha/3)`, `2s = 1 + alpha - 2 sqrt(3 alpha)`, root `5 - 2 sqrt 6`, the
number flagged above; it is not binding in the proof (exponent margin 0.035 at `alpha_0`). The
side conditions `alpha < 2/7` (Section 4.3.5) and `L < N^{5/6}` (Lemma 5 hypothesis) do not bind.

## 7. The ladder, and what the gaps mean

| level | value | status |
|---|---|---|
| CMZZ, proved, force in `L^1_t C^{1,eps} ∩ L^inf_t L^2` | `alpha < 0.0927` | theorem |
| A-B schedules as written, dissipation added naively | no fixed `alpha > 0` survives all stages (`Q_q = Q_* + q`, `Q_i = q_0 + i`); earliest-stage ceilings `1.5e-6` (`Q_1 = 201`), `5.7e-14` (`q_0 = 2^20`); stage-`q` ceiling `beta/(2 Q_{q-1} Q_q) -> 0` | schedule inequality, derived + checked |
| A-B mechanism re-optimized, both fields damped | `alpha < 1/2` | derived + checked |
| A-B mechanism re-optimized, momentum only, scalar clean | `alpha < 1` | derived + checked |
| bounded velocity wall | `alpha = 1` | derived (rate `U lambda >> nu lambda^alpha` with `U` bounded) |
| Navier-Stokes | `alpha = 2` | - |

The interesting gap is 0.0927 to 1/2, a factor 5.4, and it is a gap in optimization, not in
mechanism: both endpoints obey the same inequality `alpha < (growth-rate exponent)/(separation
exponent)`, which in both is the damping-at-most-growth constraint (assumption 1). The gap 1/2 to 1 is a gap in mechanism: it needs a rate proportional to `lambda`
rather than `lambda^{1/2}`, i.e. an instability that is first order in the background rather
than a two-field Rayleigh-Taylor pair. The gap 1 to 2 needs unbounded velocity, which no
construction in this family has.

## 8. Pricing a "sharpen the hypodissipative threshold" campaign

*Conjectured, with reasons.*

* **Not a re-parameterization of the released construction.** The schedules are tuned so that
  each stage's prescribed gain equals its time budget almost exactly (`Gamma_q tau_q = 1 + L_q`
  against a required `L_q`); any damping consumes more than the slack of 1. Adding slack is
  cheap, but the amplitude matrix then stops being antidiagonal, and two load-bearing steps
  depend on that: the steering step that returns `Omega_q` to zero by selecting a pulse
  amplitude (Theorem 3.2(i), a sign-change and continuity argument in `h`), and the propagator
  bound (3.21) `||U(t,s)|| <= C_R |y(t)|/|y(s)|`, which controls every correction relative to
  the growing amplitude. Estimate: a new Section 3 and a reworked Section 5, not a re-run of
  Section 8.
* **Cheapest real gain is on the CMZZ side.** The four-inequality optimization with an explicit
  discriminant is their heuristic, which the paper says is not the real constraint set. In the proof
  `alpha_0` is set jointly by two estimates: the self-interaction bound (Section 4.3.2, Lemma 7),
  which fixes `l`, and the outer-velocity linearization remainder (Section 4.3.4, Lemma 11), which
  fixes `s`; the threshold function at the optimum is `4s = 2 + 3 alpha - 7a - 2 alpha/a`. Improving
  either estimate moves `a` and `R` and changes the threshold through a square root, and the
  remainder estimate is at least as load-bearing as the cubic term. Still a parameter
  re-optimization rather than a new mechanism.
* **Passing 1/2 needs a new mechanism**, by section 4.
* **The `C^infinity` force is the expensive requirement.** If section 5 is right, asking for a
  smooth force and a fixed `alpha > 0` at once asks for `Q` bounded and `k_q` unbounded. The
  first question for anyone holding the unreleased hypodissipative manuscript is which of the
  two was given up, or what replaces the `r_* <= Q - 1` correction budget.
* **Next toy after this one**: put the steering pulse back in and ask whether a root `h` with
  `Omega_q(t_q^hat; h) = 0` still exists when `Omega` is damped. That is the first place the
  released control argument breaks (the profile-harmonic and nonlocal-residual issues under
  assumption 1 break earlier, at the exact-wave step), and it is one-dimensional root finding on
  the same ODE.

## Table: empirical `alpha` threshold from the toy

Toy: per stage `k`, integrate `y' = [[-delta D, Gamma],[Gamma, -D]] y` with
`Gamma_k = sigma_{k-1} sin s_k`, `sigma_k = exp(beta x_k / 2)`, `x_k = log lambda_k = Q x_{k-1}`,
`D_k = nu lambda_k^alpha`, required gain `L_k = 6 x_k`; a schedule passes iff every stage
`k = 3..K` reaches its gain within its budget, and `alpha*` is found by bisection.
`x_1 = 2000`, `K = 6`, `nu = 1`, time-budget slack factor 2. "paper" uses A-B's own
`s_k = L_k sigma_{k-2}/sigma_{k-1}` (the toy writes it as `sin s_k`; (3.52b) forces `s_j <= eps_sch`,
so the difference is harmless) and `tau_k = c(1+L_k)/Gamma_k`; "wide-tight" uses
`sin s_k = 1/k^2` with the same minimal `tau`; "wide-slow" uses `sin s_k = 1/k^2` and a
summable budget `tau_k = 1/k^2` independent of the gain. `delta = 1` also enforces retention.
Verified: the closed-form gain agrees with `scipy` Radau integration of the log-radial system
to `<= 6e-15` relative error over seven regimes, including `D/Gamma = 10^3`.

| schedule | delta | Q=1 (geom.) | Q=1.5 | Q=2 | Q=3 | Q=5 | Q=10 | Q=20 | Q=201 | fitted law |
|---|---|---|---|---|---|---|---|---|---|---|
| paper | 0 | 0.2500 | 0.1119 | 0.0627 | 0.02781 | 0.01000 | 0.002500 | 6.25e-4 | 6.19e-6 | `beta/(2Q^2)` |
| paper | 1 | 0.2499 | 0.1118 | 0.0627 | 0.02781 | 0.01000 | 0.002500 | 6.25e-4 | 6.19e-6 | `beta/(2Q^2)` |
| wide-tight | 0 | 0.2492 | 0.1663 | 0.1248 | 0.08323 | 0.04996 | 0.02499 | 0.01250 | 1.24e-3 | `beta/(2Q)` |
| wide-tight | 1 | 0.2469 | 0.1660 | 0.1246 | 0.08317 | 0.04994 | 0.02499 | 0.01250 | 1.24e-3 | `beta/(2Q)` |
| wide-slow | 0 | 0.4954 | 0.3298 | 0.2479 | 0.16570 | 0.09963 | 0.04990 | 0.02497 | 2.49e-3 | `beta/Q` |
| wide-slow | 1 | 3.7e-4 | 3.2e-4 | 1.0e-4 | 2.0e-5 | 2.6e-6 | 1.6e-7 | 1.0e-8 | 9.9e-13 | retention-limited, `-> 0` |

All entries for `beta = 1/2`. The `Q=1` column is `x_k = x_1 + 2(k-1)` with `x_1 = 5000`,
`K = 8`. Scaling in `beta` is exactly linear: repeating at `beta = 1/8` and `beta = 0.9`
reproduces `alpha* Q^2 = beta/2` (paper), `alpha* Q = beta/2` (wide-tight),
`alpha* Q = beta` (wide-slow, `delta = 0`) to four digits. `alpha*` is stable in `K`
(varies by 1% from `K=4` to `K=8`; the retention-limited `delta = 1` wide-slow row is the exception,
and its `Q = 201` cell was originally run at `K = 5`, giving `1.8e-10`; at the stated `K = 6` it is
`9.9e-13`. That row's retention test checks `D_{k-1} tau_k <= log 2` for the next stage only rather
than the full remaining time `sum_{j >= k} tau_j`; the difference is a polynomial factor with no
effect on exponents) and independent of `nu` and of the gain constant to within
0.1% over `nu` in `[1e-3, 1e3]` and gain constants `6` to `600`: it is a pure exponent
condition, as the derivation predicts.

Re-run at the released values `beta = 1/8`, `Q = 201`, the same three schedules give
`alpha* = 1.55e-6` (paper), `3.11e-4` (wide-tight) and `6.22e-4` (wide-slow). The fixed
`Q = 201` here is a stand-in for the earliest stage: the released Boussinesq schedule has
`Q_q = 200 + q`, and its stage-`q` ceiling `beta/(2 Q_{q-1} Q_q)` decreases without bound (section 3),
so it tolerates no fixed `alpha > 0`. The most permissive re-optimization inside the same amplitude
ODE at `Q = 201` reaches `6.2e-4`. Neither is near 0.0927, and the reason is `Q`, not `beta`.

## Reproduction

Toy and sweeps: `/private/tmp/claude-501/-Users-spaceman-Desktop/3b835b56-c69b-49b8-a172-7b94c2a0802b/scratchpad/scout-ceiling/toy.py`
(numpy + scipy; `toy.alpha_threshold(...)` is the entry point, `toy.log_gain_numeric` is the
integration cross-check). Patched 2026-09-08: `run_schedule` now passes `log(Gamma tau)` directly to
`log_gain_exact` instead of forming `logGamma + logtau`; the sum cancelled two logs of order `1e18`
and returned `0.0` for `K >= 10` at `Q ~ 200`, so `alpha*` read `0` there. Table values are unchanged
(error at `K = 6` was about `1e-2` in a gain of about 36); `K = 10, 12` now reproduce `6.188e-6`.

## Referee corrections (2026-09-08)

Sixteen findings from a two-lens referee pass (growth-rate law and ceiling formula; reading of
the released schedules; CMZZ `alpha_0` derivation). Each was re-checked against the note, the
extracted paper texts (`boussinesq.txt`, `euler.txt`, `cmzz.txt` via pdftotext) and the toy before
being applied. All sixteen were confirmed; none rejected. Edits are in place above with minimal
diffs; the pre-correction note is preserved at
`scratchpad/SCOUT_before_referee.md`.

1. **Vorticity amplitude confused with growth-rate exponent** (MAJOR, section 3, applied). The
   eigenline relation `Omega = (lambda r / sqrt(A_old)) Theta` (Subsection 1.2) and the
   normalization (3.17) give peak growing vorticity `A_q/sigma_{q-1} ~ lambda_q^{beta(1-1/(2Q))}`;
   `~ sigma_q` is the post-steering value of Theorem 3.2(iii). The governing exponent is the rate
   exponent `beta/2`. Section 4's "rate is `Omega = sqrt(A)`" renamed. Ceiling arithmetic unchanged.
2. **Insertion-angle lever not admitted by the released steering** (MAJOR, section 4, applied).
   Subsection 3.7.6 needs `Lambda_q sin s_q <~ 1` (`|Z| <= (c_p H_2 + C eps) L sin s`,
   `|d_2| <= C(Ls)^2`) and (3.9) needs `c_p Lambda_1 H_2 sin s_* <= 1/4`; with `Lambda_q = L_q`,
   `s_j ~ 1/j^2` sends `L_j s_j` to infinity. Lever relabelled conjectured, requiring `Lambda_q = O(1)`.
   The `beta/(2Q)` and `beta/Q` ceilings are ceilings of the amplitude ODE.
3. **`beta -> 1, Q -> 1` outside the correction budget** (MINOR, section 4, applied).
   `delta = Pi_q^5 lambda_q^{-7/8}` (Proposition 7.2) with `Pi_q = lambda_{q-1}^4 (log)^10` gives
   `beta < 1 - 20/Q` for `delta <= 1/2`; (7.15) gives `Q >= 9k + 42`. `1/2` marked as an
   amplitude-ODE supremum.
4. **Assumption 1 is the `n = 1` harmonic model** (MINOR, assumptions, applied). The profile (4.1)
   has higher harmonics and the fractional Laplacian is nonlocal; exponent conditions unaffected,
   but this is the first structural casualty, earlier than the steering root.
5. **`beta/Q` is scalar-gradient blowup only** (MINOR, section 4, applied). Overdamped regime,
   vorticity quasi-static; (3.4), (3.19) and the `limsup ||omega||` conclusion of (1.5) lost.
6. **Table cell wide-slow / `delta = 1` / `Q = 201`** (MINOR, table, applied). Toy at the stated
   `K = 6` gives `9.855e-13`; the printed `1.8e-10` matches the `K = 5` run (`1.756e-10`). The rest
   of that row reproduces at `K = 6`. Cell replaced; retention-test caveat noted.
7. **Angle rule is an angle, not a sine; Euler symbol is `Lambda_i`** (MINOR, section 3 and
   table caption, applied). (3.10) `s_q = L_q sigma_{q-2}/sigma_{q-1}`; Euler `s_i = Lambda_i
   sigma_{i-2}/sigma_{i-1}`; `Gamma_q in [(1 - eps^2/6) L_q sigma_{q-2}, L_q sigma_{q-2}]`.
8. **`s_j ~ 1/j^2` violates (3.52b)-(3.52c)** (MINOR, section 4, applied, merged with item 2).
   `L_j s_j <= eps_sch s_{j-1}` forces `s_j ~ eps^j / prod L_i`; only logarithmic factors lost,
   `beta/(2Q)` stands.
9. **The 120 in `120 k_q <= Q_q` is load-bearing** (MAJOR, section 5, applied). (3.52a)'s proof
   bounds `L_j/L_{j-1} >= (41/8) Q_j/(Q_{j-1}/120 + 41/8) >= 100`; with `Q/c` this needs `c >= ~20`
   (large `Q`) or `c >= ~39` (`Q = 201`); `a_L = Q_1/120 + 41/8` also uses it. Euler: (7.15)
   analogue `9k + 51 <= Q_m` and (12.2) confirmed, but the binding use is `mu_j = 2048 k_j + 9216
   <= Q_j/4`, i.e. `Q_j >= 8192 k_j + 36864`. `C^k` cost restated as `Q >~ c k`, ceiling
   `beta/(2ck)`.
10. **Direction of the `k_q -> infinity` statement** (MINOR, section 5, applied). Lemma 7.4 and
    Theorem 8.2 deliver it; Lemma 9.3 and Proposition 9.4(b) via (9.6) require it; structurally
    necessary because uniform `C^k` control exists only for `k <= k_q` (Definition 7.1).
11. **`1.5e-6` / `5.7e-14` are first-stage numbers, not schedule ceilings** (MAJOR, sections 3, 7,
    8, applied). Both schedules hard-wire `Q_q -> infinity`, so the stage-`q` ceiling
    `beta/(2 Q_{q-1} Q_q) -> 0` and no fixed `alpha > 0` survives all stages. Toy check with
    `Q_q = 200 + q`: `alpha* = 1.773e-6, 1.511e-6, 1.480e-6, 1.452e-6` at `K = 3, 4, 6, 8`, each
    equal to the closed form. Section 5 re-presented as the argument for a bounded-`Q` variant.
12. **Toy cancellation at `K >= 10`** (MINOR, toy, applied). `logGamma + logtau` cancelled two
    logs of order `1e18` and returned `0.0`; `run_schedule` now passes `log(Gamma tau)` directly.
    Table unchanged; `K = 10, 12` reproduce `6.188e-6`. Original saved as `toy.py.bak`.
13. **CMZZ Section 1.2.4 does not derive `alpha_0`** (MAJOR, section 6, applied). The section
    ends with "s(alpha) > 0 ... holds if alpha is a small positive number" and the disclaimer
    about other sources of error; `alpha_0` is fixed in Section 4 by `R = sqrt(2/(7 alpha))`,
    `A = N^{sqrt(2 alpha/7)}`, `L = N^{1+alpha-s-2 sqrt(2 alpha/7)}`, `s = (3 alpha + 2 - 2 sqrt(14
    alpha))/4`, and is the root of `s = 0`.
14. **The `5 - 2 sqrt 6` vs `0.0927` discrepancy** (MAJOR, section 6, applied, resolved). The
    tightening constraint is the Section 4.3.4 outer-velocity linearization remainder
    (`a + 1 + s + 2/R <= 3l`), not the heuristic `K_n/L_n <= 1/M_{n-1}`. Re-derived and checked
    numerically: root `0.0926655`, heuristic root `0.101021`, heuristic third-constraint margin
    `0.035` at `alpha_0`; `alpha < 2/7` and `L < N^{5/6}` do not bind.
15. **`0.176` and `1.90` are heuristic values, not CMZZ's** (MINOR, sections 6 and 7, applied).
    CMZZ use `a = 0.163`, `R = 1.756`; `alpha = a/R` is the damping-at-most-growth constraint by
    construction, so the "same shape as `beta/(2Q)`" remark restates assumption 1.
16. **"Raises `alpha_0` linearly" via the cubic constraint** (MINOR, section 8, applied). The
    threshold is set jointly by the self-interaction bound (Lemma 7) and the remainder (Lemma 11)
    and moves through a square root; "linearly" dropped.

**What survives.** The mechanism ceiling formula stands: under assumption 1 (single-harmonic
diagonal damping at the layer frequency) and assumption 2 (nothing else changes), a stage grows
iff `mu_+ tau_q >= L_q` with `mu_+` the top eigenvalue of `[[-delta D, Gamma],[Gamma, -D]]`, the
rate is `Gamma_q = sigma_{q-1} sin s_q` with rate exponent `beta/2`, and the resulting exponent
conditions are `alpha < beta/(2 Q_{q-1} Q_q)` for the released angle rule, `alpha < beta/(2Q)` for
the widest admissible angle with both fields damped (amplitude-ODE supremum `1/2`), and
`alpha < beta/Q` for momentum-only damping with long stages (scalar-gradient blowup only,
overdamped regime, supremum `1`); all three are reproduced by the toy to four digits and are
exponent conditions independent of `nu` and of the gain constant. The `C^infinity`-force-vs-fixed-
`alpha` incompatibility claim survives with this exact scope: it is a statement about the A-B
correction scheme under assumptions 1 and 2, not about all schemes; for the released schedules it
is redundant, because `Q_q = Q_* + q` and `Q_i = q_0 + i` already send the stage ceiling to zero
(`1.5e-6` and `5.7e-14` are earliest-stage values); for a bounded-`Q` re-optimization it is the
operative argument, since the papers' own budgets `120 k_q <= Q_q` (Boussinesq, visibly
relaxable to about `20 k_q <= Q_q`, not `9 k_q`) and `k_j <= Q_j/c_Q` with `c_Q >= 2^18` (Euler,
binding at `Q_j >= 8192 k_j + 36864`) cap `k_q`, while `C^infinity` needs `k_q -> infinity`
(Lemma 9.3, Definition 7.1), so a `C^k` force buys at most `alpha <~ beta/(2ck)`. The CMZZ
`alpha_0` derivation status: Section 1.2.4 is a heuristic yielding only "alpha small";
`alpha_0 = (22 - 8 sqrt 7)/9` is fixed in Section 4 and is the root of `s = 0` for the exponent set
jointly by the self-interaction bound (Section 4.3.2) and the outer-velocity linearization
remainder (Section 4.3.4); the note's re-derivation of that pair reproduces `alpha_0` exactly,
and the earlier flagged discrepancy is closed. The ladder's endpoints, CMZZ at `0.0927` and the
amplitude-ODE ceiling at `1/2`, both obey `alpha < (growth-rate exponent)/(separation exponent)`,
which in both is the damping-at-most-growth constraint itself.

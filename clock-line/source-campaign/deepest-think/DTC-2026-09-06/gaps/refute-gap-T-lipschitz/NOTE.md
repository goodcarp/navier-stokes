# REFUTATION of `gaps/gap-T-lipschitz` — the lemma stands, the bridge to GAP T does not

Seat `gaps/refute-gap-T-lipschitz`, DTC-2026-09-06.  Adversarial re-check of
`gaps/gap-T-lipschitz`.  Every number below came out of a script in **this** folder that I wrote
and ran (`r1`…`r4`, re-asserted by `check_constants.py`); `SHA256SUMS` is computed, never typed.
Nothing outside this folder was written.  The target seat's scripts were re-run **from a copy** in
a scratch directory, never in place.

---

## 0. Verdict

| claim under test | verdict |
|---|---|
| LEMMA T (the estimate `(T)` with hypotheses (D),(H1),(H2),(H3)) | **STANDS.**  Every constant re-derived independently and confirmed. |
| the sharp kernel constants `3/(8pi^2)`, `3/(2pi^2)`, the shell constant `pi^3` | **STANDS** (re-derived from scratch; a new exact identity supplied, below) |
| necessity of (H3): the control C1 is a genuine `C^1` `SO(4)`-equivariant diffeomorphism | **STANDS** (I attacked the diffeomorphism property directly; it survives) |
| Corollary 1 `|Delta|/a_T <= 11.7809725 mu` | **STANDS** (arithmetic correct) |
| **"GAP T is downgraded from GAP to PROVED-MODULO L3v"** | **REFUTED** |
| **"With prove-lagrangian's own §4(2) figures both hypotheses hold with `mu, mu_J = O(c/L)`"** | **REFUTED — measured `mu = 0.3853`, exactly `L`-independent** |
| "the `O(c/L)` shear between shells … is a radial rearrangement — the least dangerous of the three error terms" | **REFUTED — it is an angular shear; up to `17.16°` of polar-angle shift, `tan`-ratio factor `3.375`** |
| "`Delta = 0` … for every amplitude `mu` (however large)" (§3(N)) | **OVERCLAIM** — true only inside the diffeomorphism range; the seat's own `s5` `mu=3` row contradicts it and `check_constants.py` filters that row out |

Severity: **MAJOR**, not fatal.  Lemma T is a correct lemma about a different object from the one
`prove-lagrangian` §4(2) names, and the repair is available inside the seat's own machinery — I ran
it (§4) and it works, yielding relative error `O(1/L)` after all.  But the seat's headline status
downgrade is not earned by what it proved.

---

## 1. Reproduction (`re-run from a copy`)

`shasum -a 256 -c SHA256SUMS` on the target folder: **14/14 OK**.  `s1`, `s2`, `s3`, `s4`, `s5`
re-run in a scratch copy reproduce every printed digit; `check_constants.py` prints
`ALL ASSERTIONS PASSED`.  No arithmetic error was found anywhere in the seat's own numbers.

## 2. What I independently confirmed (`r1_algebra.py`)

* `grad K(w) = -(3/(8 pi^2))|w|^{-5}(e_z - 5 c what)`, `|e_z - 5 c what|^2 = 1 + 15 c^2 <= 16`
  (symbolic residual `0`), so `sup |w|^5|grad K| = 3/(2 pi^2) = 0.15198177546350666`, and
  `sup |w|^4 |K| = 3/(8 pi^2) = 0.037995443865876666`.  **The seat's Step 2 is right.**
* **A new exact identity, which the seat needed but never wrote down.**  With
  `g(phi;lambda) = |T_lambda x|/|x| = sqrt(lambda^2 sin^2 phi + lambda^{-4} cos^2 phi)`,

  ```
  J(lambda) := int_0^pi sin^2(phi) / g(phi;lambda)^4 dphi = pi/(2 lambda)   EXACTLY.
  ```
  Route (each step symbolic, residual `0`): `d/dt[atan(t sqrt(a/b))/sqrt(ab)] = 1/(b+a t^2)`;
  `t = tan phi` gives `int_0^pi dphi/(b+(a-b)sin^2 phi) = pi/sqrt(ab)`; differentiate in `a`;
  set `a = lambda^2`, `b = lambda^{-4}`.  30-digit `mpmath` confirmation at six `lambda`,
  max error `3.0e-17`.  This is what makes the seat's Step-3 constant `pi^3 M L/lambda` true in
  the `x`-variables (`s1` proved only the `u`-substitution route; `s3` only checked it numerically).
* `int_0^{pi/2} cos phi sin^2 phi / g^5 dphi = lambda/3` exactly ⟹ **the strain density is
  `(M/2) lambda(rho)` per unit `log rho`, for the bang-bang cap, for ANY profile `lambda(.)`.**
* For the **shell-dependent** strain `Lambda(x) = T_{lambda(|x|)} x`,
  `det D Lambda = lambda^2 [ 1 + (dlog lambda/dlog rho)(sin^2 phi - 2 cos^2 phi) ]`; checked against
  a finite-difference `5x5` determinant of the actual map at 400 random points, max relative error
  `3.2e-10`.
* Instrument cross-checks against the neighbouring seats (`r2`, my own quadrature, nothing imported):
  `a[T_1](0)/(ML) = 0.500000000000011`, `a[T_{3/2}](0)/(ML) = 0.750000000000017` (a fourth
  independent confirmation of `P_1(lambda) = lambda`), and `kappa_delta(7.5°) = 0.4997212305210992`,
  agreeing with the target seat's `0.4997212305210886` to `1.1e-14` and with
  `prove-lagrangian` §2 to `3e-8`.

## 3. THE REFUTATION — Lemma T compares against the wrong map (`r2_shell_dependent.py`)

**Verbatim at source.**  `lower/prove-lagrangian` §5, status row **T**:

> `flow map on [0,τ] = T_{λ(s)} up to relative O(c/L); and (1.1) is stable under that`

and §4(2), third bullet:

> `the shear between shells: ρ ∂_ρλ/λ = −τM/2 = −c/(2L), so the 5D Jacobian and the log-radial measure are preserved to 1 + O(c/L)`

`T_{λ(s)}` is indexed by the **material shell** `s`.  Its closed form is that seat's own (4.1),
`λ(σ,θ) = (1 − (1−σ)κθ/2)^{-2}`, `σ = log(ρ/ρ0)/L`, and step **P** consumes it shell-by-shell
(`H(σ,θ) = ∫_σ^1 e^{F(σ',θ)} dσ'`).  At the terminal time `κθ = 2(1−√(2/3)) = 0.3670068381`,

```
lambda(sigma=0) = 3/2   ...   lambda(sigma=1) = 1 ,   log-mean = sqrt(3/2) = 1.224744871392 .
```

The bullet `ρ∂_ρλ/λ = −c/(2L)` is a **rate per log-radian**.  Integrated across the `L`
log-radians of the shell it gives a total spread of `log(3/2)` — **`Theta(1)`, uniformly in `L`**.
The seat's §6 reads the rate as the total and concludes `mu = O(c/L)` against a **single** `T_lambda`.

Measured (own quadrature; grid `4001 x 2001`; the comparison `lambdabar` chosen to *minimise* the
sup relative displacement, so these are the best possible constants for the seat's hypotheses):

| `L` | best `lambdabar` | `mu` of (H2) | `mu_J` of (H3) | Lemma-T bound `/ a(0)` |
|---|---|---|---|---|
| 8.318 | 1.1770 | **0.385329** | 0.799702 | **70.44** |
| 20 | 1.1770 | **0.385329** | 0.697167 | **67.90** |
| 50 | 1.1770 | **0.385329** | 0.653365 | **66.78** |
| 100 | 1.1770 | **0.385329** | 0.638764 | **66.40** |
| 400 | 1.1770 | **0.385329** | 0.627813 | **66.11** |

`mu` is **identical to twelve digits at every `L`** (spread `< 1e-12` across `L = 8.3 … 400`): the
displacement depends only on the ratio `lambda/lambdabar`, and that profile is `L`-independent.
So (H2) against a single `T_lambda` holds with `mu = 0.3853`, never with `O(c/L)`; and Lemma T,
fed its own honest constants, certifies a relative error of **6600 % – 7000 %** — vacuous at every
`L`.  Taking the single `lambda` at either end of the range `prove-lagrangian` names (`lambda = 1`
or `3/2`) misstates `a(0)` by `18.4 % – 22.4 %`, again with no decay in `L`.

**The seat's §6 sentence — "With prove-lagrangian's own §4(2) figures … both hypotheses hold with
`mu, mu_J = O(c/L)`" — is false**, and it is the only sentence that connects Lemma T to GAP T.

### 3b. The seat's "radial rearrangement" consequence is backwards (`r3_radial_claim.py`)

§3(N) says the shear between shells "is a radial rearrangement — the least dangerous of the three
error terms".  A radial rearrangement is `Psi(u) = s(|u|) u`: polar angle fixed.  `T_lambda` moves
the polar angle by `tan phi' = lambda^3 tan phi`, and `lambda^3` runs from `3.375` to `1` across
the shell.  Measured against `lambdabar = sqrt(3/2)`: **maximum polar-angle shift `17.16°`**,
sup relative angular displacement `0.3175`, and at some points the displacement is `100 %` angular
(`max angular fraction 0.99999999999957`).  The shear between shells lives in exactly the
**angular** direction the seat's own control C1 identifies as the dangerous one — not in the null
direction.

## 4. THE REPAIR — and it works (`r2_shell_dependent.py`, part E)

Nothing is lost mathematically; only the statement has to change.  Two routes, both measured:

**(R1) State Lemma T against the shell-dependent reference.**  Steps 1, 2, 4, 5 of the seat's proof
never use that the comparison map is a single `T_lambda`; Step 3 generalises by my §2 identity,

```
int_S |eta_0| |Lambda x|^{-4} dx_5  <=  pi^3 M int_0^L dlog(rho)/lambda(rho) ,
```
so Lemma T holds verbatim with `T_lambda -> Lambda` and `lambda M L -> M int lambda dlog rho`.
Against that reference (H2) is *vacuous* (`mu = 0`) and only (H3) enters, with

```
mu_J(local) = 0.8990 / L    (measured: 0.10808, 0.044949, 0.017980, 0.008990, 0.002247 at
                             L = 8.318, 20, 50, 100, 400 -- exactly 0.8990/L at all five)
```
— i.e. `1 + O(c/L)` with constant `1.109 c`, `c = 2 log(3/2)`, which is the §4(2) bullet as
written.  Lemma T then certifies

```
|Delta| / a(0)  <=  2.03/L, 2.08/L, 2.10/L, 2.11/L, 2.12/L   -> O(1/L).
```

**(R2) Or keep the single-`lambda` lemma and decompose.**  Split `S` into `L` shells of unit
log-width, apply Lemma T on each against its own `T_{lambda(mid)}`, and sum (legitimate: `a(0)` is
linear in `eta`).  Measured per-shell constants and summed bound:

| `L` | shells | `mu` max | `mu_J` max | `sum(bound)/a(0)` |
|---|---|---|---|---|
| 8.318 | 8 | 0.056570 | 0.171662 | 1.06807 |
| 20 | 20 | 0.022537 | 0.068632 | 0.36963 |
| 50 | 50 | 0.009000 | 0.027162 | 0.13967 |
| 100 | 100 | 0.004497 | 0.013533 | 0.06856 |
| 400 | 400 | 0.001124 | 0.003374 | 0.01691 |

`sum(bound)/a(0) ~ 6.8/L` — `O(1/L)`, as needed.  **So the seat's conclusion survives; its proof of
that conclusion does not.**  The missing paragraph is the decomposition (or the restatement), not
new analysis.

**What the shear actually costs.**  The true error of the shell-resolved reference, i.e.
`| a[eta_0 o Lambda^{-1}](0) − (M/2)∫lambda dlog rho |/a(0)`, is entirely the Jacobian-shear term
and measures `0.3351/L, 0.3431/L, 0.3467/L, 0.3479/L, 0.3488/L` — converging to `≈ 0.349/L`.
(Quadrature: 200 vs 400 nodes agree to `1.6e-14`.)

**Honest reach of the constant.**  The proved constants only bite for very large `L`: route (R1)
needs `L >~ 21` and route (R2) `L >~ 68` for a 10 % relative error, i.e. `log Re_E >~ 45` and
`>~ 139` respectively via `log Re_E = 2L + 2log(1/delta) − 0.7031659`.  The seat's remark that
"sharpening is not worth doing" should be read against that: the `o(1)` is genuine but the
crossover is astronomically far out.

## 5. What survived my attack

* **The control C1 is a genuine diffeomorphism** (`r4_control_diffeo.py`).  I expected the seat's
  parenthetical `1 + beta' >= 1 − mu_beta > 0` to be false, because `beta'` carries terms of size
  `mu_beta/phi_c` (up to `1.35e+12` at `phi_c = 1e-12`).  It is nonetheless **exactly right**:
  on a grid refined inside both ramps (`200001` points per ramp), `min(1 + beta') = 1 − mu_beta`
  to `1e-9` at every one of the 15 `(mu_beta, phi_c)` pairs, `phi + beta` is monotone, and it maps
  `[0,pi]` onto `[0,pi]`.  So (H1) holds and the necessity of (H3) is established as claimed.
* Every constant of Lemma T, re-derived independently (§2).

## 6. Smaller findings

1. **`Delta = 0` "for every amplitude" is an overclaim.**  §3(N)/§3(A)/§5 assert the log-periodic
   null direction holds "for every … amplitude `mu` (however large)" and "unchanged past the
   diffeomorphism limit `mu = 0.9342`".  Outside `mu sqrt(1+(k pi/L)^2) < 1` the map `s(|u|)u` is
   not injective (and for `mu > 1`, `s` changes sign), so `a[eta_0 o Phi^{-1}](0)` is not defined by
   the change-of-variables formula at all.  The seat's own `s5` run reports `rel = 1.076` for
   `k = 2, 4` at `mu = 3.0`, and stores it as `null_max_rel_k_even = 1.076e+00` — a number that
   appears in no table in `NOTE.md`; `check_constants.py` asserts the null property only over
   `mu <= 0.5`, which is the sub-range where it is true.  A self-check filtered to the passing
   rows is not a control.  (My own re-derivation of that integrand differs from the seat's by a
   factor `sign(s)` once `s < 0`, which is exactly the ill-definedness.)
2. **`mu_J > mu` in all twelve verification rows** (by `2.6x` to `7.5x`), so the clean corollary
   constant `C(0+) = 15 pi/8`, which presumes `mu_J <= mu`, is in force in none of them.  Not an
   error; but the table does not test the regime the corollary is quoted in.
3. **Corollary 2 may be repairing the wrong statement.**  `prove-lagrangian` §4(2) writes the
   metric as `|δX|/|X|` with `X` the (current, deformed) position, and its own bullet measures the
   remainder displacement `C_3 M rho tau` against the current `rho`.  The `max(lambda^{-1},lambda^2)`
   conversion is therefore probably not needed against that source; I cannot check the orchestrator
   brief the seat quotes, so I record this as unresolved rather than wrong.
4. **The status label is too generous even on the seat's own account.**  §6 calls GAP T
   "PROVED-MODULO L3v", but what it leaves open is bullet (a) — the trajectory closure — which needs
   L3v *and* the integration of the velocity remainder, and (per §3 above) also needs the
   shell-dependent reference.  "PROVED-MODULO L3v" understates it.

## 7. Corrected statement

> **LEMMA T (as proved).**  Unchanged — it is correct as stated, for a **single** `T_lambda`
> comparison map, under (D),(H1),(H2),(H3), with `C(0+) = 15 pi/8` and relative form `15 pi/4`.
>
> **What it does NOT give.**  It does not apply to the map `prove-lagrangian` §4(2) names.  The
> flow map is compared there to `T_{lambda(s)}` with `lambda` running over `[1, 3/2]` across the
> shell; against any single `T_lambdabar` the best achievable hypothesis constant is
> `mu = 0.385329`, independent of `L`, and (T) then certifies only a `66x` relative error.
>
> **Repaired bridge (verified, not proved here).**  Either (R1) restate Lemma T with the
> shell-dependent reference `Lambda(x) = T_{lambda(|x|)}x` — Steps 1,2,4,5 are unchanged and Step 3
> generalises by `int_0^pi sin^2/g(.;lambda)^4 = pi/(2 lambda)` to
> `int_S |eta_0||Lambda x|^{-4} <= pi^3 M int dlog rho/lambda(rho)` — whereupon `mu = 0`,
> `mu_J = 0.899/L`, and the relative error is `<= 2.12/L`; or (R2) decompose `S` into unit-log-width
> shells and sum, giving `<= 6.8/L`.  Both are `O(1/L)`, the order `prove-lagrangian` already
> carries.  Writing out (R1) or (R2) is a short exercise; it has not been written.

## 8. Remaining gap

GAP T is **still a gap**, in two pieces:

1. **the restatement/decomposition above** — cheap, verified numerically here, unwritten;
2. **the PDE half, untouched by both seats**: that the true Navier–Stokes flow map on `[0,tau]`
   satisfies, against the shell-dependent `T_{lambda(s)}`, a displacement bound `mu = O(c/L)`
   *and* a Jacobian bound `mu_J = O(c/L)`, uniformly over the whole shell `rho_0 < |x| < R` and
   including both boundary layers (where `prove-lagrangian` §4(2) itself records
   `sup|a − A(rho)| = 0.7047 M` at the outer edge, against `0.3187 M` one octave in).  Nothing in
   either seat bears on this.  L3v remains SKETCH.

## 9. Files

| file | what it establishes |
|---|---|
| `r1_algebra.py` / `r1_results.json` | independent re-derivation of every Lemma-T constant; the exact identity `int_0^pi sin^2/g^4 = pi/(2 lambda)`; `det D Lambda`; the per-shell core `lambda/3` |
| `r2_shell_dependent.py` / `r2_results.json` | the refutation (`mu = 0.385329`, `L`-independent; bound `66x`) and both repairs; instrument cross-checks |
| `r3_radial_claim.py` / `r3_results.json` | the shear between shells is angular (`17.16°`, `lambda^3` ratio `3.375`), not radial; the null direction past the diffeomorphism limit |
| `r4_control_diffeo.py` / `r4_results.json` | the seat's control C1 IS a diffeomorphism — `min(1+beta') = 1-mu_beta` exactly at all 15 `(mu_beta, phi_c)` |
| `check_constants.py` | re-asserts every number above against the JSONs |
| `SHA256SUMS` | computed, never typed |

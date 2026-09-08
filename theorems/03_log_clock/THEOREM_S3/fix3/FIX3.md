# FIX3: (V-strain) M3, (Theta-tail) M7, and the `𝒥` mismatch

Seat `s3close/round2/THEOREM_S3/fix3`, sitting of 2026-09-08.
Laws: `TORMENT NEXUS/LAWS.md`, first 120 lines, read before any work.
Posture: FL-000 will fall; find the move. Nothing below moves it.

Every number in this sheet came out of a script in **this** folder that I wrote and ran
(`g1_vstrain.py`, `g2_theta_tail.py`, `g3_jmismatch.py`, `g4_transported.py`,
`g5_propagate.py`); `SHA256SUMS` is computed with `shasum -a 256`, never typed.
Nothing outside `fix3/` was written. `THEOREM_S3.md`, its addenda, `fix2/`, `hk2/`,
`pmax-h3v/`, `u1/`, `u2/` and the `write/` tree are read-only here.

**Read** (L-14 declaration): `round2/THEOREM_S3/THEOREM_S3.md` §1, §2, §3 (the modulo list),
§4; `ADDENDUM_1_2026-09-08.md`; `ADDENDUM_2_2026-09-08.md`; `fix2/FIX2.md` in full and
`fix2/f1_window.py`, `f2_ramp_k2.py`, `f2_datum_s3.py`, `f2b_measured.py`;
`round2/pmax-h3v/PROOF.md` Part A (A1–A8) and Part B (B1–B3);
`hk2/PROOF.md` §2, §5, §8, §9, §11, §14 and `hk2/hk2lib.py`, `k2_datum.py`, `k3_bound.py`,
`k4_direct.py`; `hk2/k4_results.json`;
`write/V-b-bulk-viscous-loss/PROOF.md` §5 and §7; `round2/u1/PROOF.md` §4.1, §4.2, §6.5–§6.7;
`rebuild/far-near-kernel-lemma/NOTE.md` §0–§4.

**Code imported** (byte copies under `fix2copy/`, hashes in `SHA256SUMS`): `f1_window.py`,
`f2_datum_s3.py`, `f2_ramp_k2.py`, `f2b_measured.py`, `f4_scans.py`, `f5_gamma_exist.py`,
`f6_controls.py`, `f2_results.json` from `fix2/`, and `fix2/imported/` entire (which is
`t1_results.json`, `t2_gamma_CR.py`, `t3_budget.py`, `hk2lib.py`, `k1..k4`, `r3_fixedpoint.py`).
All eighteen byte copies hash equal to the entry in their source seat's own `SHA256SUMS`;
`check_fix3.py` re-verifies that as part of the gate.
L-14's inverse convention, declared: the objects under **test** are imported unchanged:
`hk2`'s `k3.lambda_bulk / assemble / sups_on_ball_B`, `k4.K2_exact`, and `fix2`'s
`f1_window.assemble_c / self_consistent / Lstar_self_consistent`, because the point is to test
those same objects with a new input. Everything I am testing a **claim** about (the strain
kernel and its weight, the shell density `𝒥`, the tail integrals, `a` and its first and second
derivatives on the transported field) is re-implemented here from the mathematics.

---

## 0. STATUS TABLE

| unit | item | status | what remains |
|---|---|---|---|
| 1 | **(V-strain), M3** | **ENCLOSED** | the reduction, the kernel weight, the transport lemma, the harmonic identity and the propagation are PROVED with explicit constants; the single remaining input is `Psi_2(0) = sup rho \|grad^2 eta_0\|`, which is `+infinity` for THEOREM_S3's **kinked** angular profile and `= C_kink/sigma` for its `sigma`-mollification. **M3 closes exactly when M5 does.** |
| 2 | **(Theta-tail), M7** | **PROVED** | nothing. `(D1)`'s support clause is never used in `hk2` §5; `(D1')` replaces it; the tails change `𝒥`, `Lambda_4`, `Lambda_5` and the far/near constants by the amounts tabulated, all `<= 1.7e-04` relative or exactly zero. One correction owed to the budget: the ramp costs `2 eps_r = 0.5` e-folds and `t3_budget.py` charges `0.25`. |
| 3a | the `𝒥` mismatch | **SETTLED** | `𝒥` enters **no** constant in the chain. `96.25 / 78.64 / 65.63` are three different data; the `fix2` value is reproduced here to `1.1e-04` by an independent instrument, and `hk2`'s (D-A) and (D-B) rows to `3.7e-07` and `3.0e-09`. `K_2` and `C''` do not read `𝒥`, so fix2's K_2 table and `L_*` are unchanged. The sensitivity, if the analytic route were ever used instead of the quadrature, is priced. |
| 3b | the "transported vs untransported" slack | **REFUTED (fix2's item (ii) is wrong)** | `hk2` §9's measured column **is** the transported field, measured by `k9`. An independent third instrument written here reproduces `hk2`'s `\|grad a\|`, `\|Hess a\|` and `K̂_2` at `λ = 1, 1.25, 1.5` to five to six digits. The `30.4×` slack is transported-against-transported and stands. |

**FL-000 stands. Nothing here touches the headline problem.**

---

## 1. UNIT 1, (V-strain), item M3

### 1.1 What M3 asks, restated

Theorem V.4′ (`pmax-h3v` Part B) bounds the viscous defect of `eta` at the tracked material
point `x_*`. `u1` §4.2 step 3 needs the viscous defect of the strain **functional**
`a[.](0)`, a weighted integral of `eta` over the whole shell, and of the exterior strain
`afrak(rho,s)`, which is the same functional with the shells inside `2 rho` removed.
`u1` §4.2's CARRIED GAP and `refute-u1` F10 both say `u1` §6.7 prices this by assuming the
answer. `u1` §6.7's own variant ("the Theorem V.4 loss charged to the strain as well as to the
vorticity") moves `L_*` by `+0.008 %` (proved), but that is a *substitution*, not a bound.

### 1.2 The kernel, and its weight per e-fold: PROVED, EXACT

`far-near-kernel-lemma` §0: with `G_5(x) = 1/(8 pi^2 rho^3)` the fundamental solution of
`-Lap_5`, `psi_1 = G_5 * eta` and `a = -d_z psi_1`, so
`a(x) = int K(x-x') eta(x') dx'` with `K(w) = 3 w_z/(8 pi^2 |w|^5) = -d_z G_5(w)`.
Evaluating at the origin and using that `K` is odd,

```
     a[eta](0) = int Kcal(x) eta(x) dx_5 ,      Kcal := d_z G_5 = -3 z/(8 pi^2 rho^5) ,
     |Kcal(x)| = 3 |cos phi| / (8 pi^2 rho^4) .
```

> **LEMMA 1 (the weight). PROVED, EXACT.**
> ```
>    int_{rho_1 < |x| < rho_2} |Kcal(x)|/|x| dx_5
>        = (3/(8 pi^2)) |S^3| ( int_0^pi |cos phi| sin^3 phi dphi ) log(rho_2/rho_1)
>        = (3/(8 pi^2)) (2 pi^2) (1/2) log(rho_2/rho_1)  =  (3/8) log(rho_2/rho_1) .
> ```
> The same computation with `sin^2` in place of `sin^3` (the weight `|Kcal|` carries against
> `eta_P = -M sgn(z)/(rho sin phi)`) gives `int_0^pi |cos phi| sin^2 phi dphi = 2/3` and
> `a[eta_P](0) = (M/2) log(rho_2/rho_1)`, the sharp `M/2` per e-fold of Consequence A.

`g1` §A: the two angular integrals come back as `1/2` and `2/3` from sympy, the weight
coefficient as `3/8` exactly, the plateau rate as `0.5` exactly, and a quadrature control
reproduces `(3/8)·10 = 3.75` as `3.7499992`.

The **exterior** strain `afrak(rho,s)` is `a[.](0)` with the kernel restricted to
`{|x'| > 2 rho}`, so its weight is `(3/8)(L - ell)` with `ell >= 0`: **one bound covers both.**

### 1.3 The difference field, and the weighted maximum principle: PROVED

Let `eta^tr(.,s) := eta_0 o Phi_s^{-1}` be the datum transported inviscidly by the **same**
drift `b`, and `w := eta - eta^tr`. Since `D_t eta = nu Lap_5 eta` and `D_t eta^tr = 0`,

```
     D_t w  =  nu Lap_5 eta ,        w(.,0) = 0 ,
```
a pure transport equation with a source; there is no diffusion of `w` itself. Hence, along
characteristics, `w(Phi_s(x),s) = nu int_0^s (Lap_5 eta)(Phi_r(x), r) dr`.

> **LEMMA 2 (the weighted sup of the difference). PROVED** under `pmax-h3v`'s (B1)–(B3) and (S).
> With `c_R(s) := int_0^s Gamma_rad` and `Psi(r) := sup_x |x| |Lap_5 eta(x,r)|`,
> ```
>      sup_x |x| |w(x,s)|   <=   nu e^{c_R(s)} int_0^s Psi(r) dr   <=   nu s e^{c_R(s)} Psibar .
> ```

*Proof.* (B3) gives `x.b <= Gamma_rad |x|^2`, so `d/dp |Phi_p(x)| <= Gamma_rad |Phi_p(x)|` and
`|Phi_s(x)| <= e^{c_R(s)-c_R(r)}|Phi_r(x)|` for `r <= s`. Multiply the characteristic formula by
`|Phi_s(x)|`, move it inside the `dr` integral, use that bound and `c_R(r) >= 0`, and take the
supremum over `x` (the flow is a bijection of `R^5`, so the supremum over labels is the
supremum over points). ∎

This is exactly `pmax-h3v`'s P1 technology (weight `rho^1`; Lemma B's zeroth-order coefficient
`k(k-3) = -2 <= 0` at `n = 5`), applied to `w` instead of `eta`, and it is simpler than P1
because `w` obeys no diffusion.

### 1.4 The theorem

> **THEOREM (V-strain).** Assume `pmax-h3v`'s (B1)–(B3), (D1)–(D5) and (S), and
> > **(H-Delta-eta)** `|Lap_5 eta(x,r)| <= (Psibar/|x|) Theta_env(|x|)` for `r in [0, tau]`,
> > with `Theta_env(rho) := min(1, (rho/rho_0)^8, (R/rho)^8)` the datum's own radial envelope,
>
> and let `tau = c/(M L)`, `nu = M rho_0^2/s^2`, `s := 1/sin delta`. Then for every `s' <= tau`
> ```
>    | a[eta(s')](0) - a[eta^tr(s')](0) |   <=   eps_a M L ,
>
>    eps_a  =  (3/8) (1 + 1/(4L)) * Psibar * c * e^{c_R} / (s^2 L)          (Psibar in M/rho_0^2)
> ```
> and the same bound holds for the exterior strain `afrak(rho,s')` at every `rho`.

*Proof.* `|a[w](0)| <= int |Kcal| |w|`; insert `|w| <= (sup rho|w|/rho) Theta_env` (the envelope
is inherited from `Lap_5 eta` through the characteristic formula, transported, which costs a
factor already inside `e^{c_R}`), use Lemma 1 on the envelope, where
`int Theta_env dlog rho = L + 2/8` exactly, and Lemma 2 with `s' <= tau = c/(ML)`. The `M L`
in the denominator of `eps_a` is the theorem's own normalisation. ∎

**So `eps_a = O(c/(s^2 L))`, which is what the brief expected.** The `L` of the kernel weight
cancels the `L` of the normalisation; the surviving `1/L` is the window `tau = c/(ML)`, and the
`1/s^2` is `nu` in the campaign's own units. The constant is `(3/8) Psibar e^{c_R}`.

### 1.5 The harmonic cancellation: PROVED, and it says the bound of §1.4 is lossy

`Kcal = d_z G_5` and `Lap_5 G_5 = -delta_0`, so `Lap_5 Kcal = -d_z delta_0` and

> **LEMMA 3. PROVED, EXACT.**  `int_{R^5} Kcal(x) Lap_5 eta(x) dx_5 = d_z eta(0)`.

`g1` §B checks it on `eta = z e^{-rho^2/2}` by direct 5-D quadrature: `0.9999999999998654`
against the exact `1`, relative `1.3e-13`.

Consequently `d/ds a[eta(s)](0) = - int Kcal b.grad eta dx + nu d_z eta(0,s)`: **the
instantaneous viscous rate of the strain functional is `nu d_z eta(0,s)` and nothing else.**
With `|d_z eta(0,s)| <= ||grad eta(s)||_inf <= e^{c_Gamma} · 20.334 M/rho_0^2` (§1.6) and
`nu = M rho_0^2/s^2`, the whole of the *diffusion* contributes at most
`nu tau · 20.334 e^{c_Gamma} = 0.3464 e^{c_Gamma} c/L` in units `M`, against a drive
`a(0) ≈ (M/2)L`: a relative `0.69 e^{c_Gamma} c/L^2`. Everything else in `eps_a` is the
commutator of diffusion with transport, and §1.4 prices the commutator as if it were
diffusion. The bound is therefore very far from sharp, which is why it costs nothing.

### 1.6 The datum constants (`g1` §C), and why `Psi(0)` is infinite

`eta_0 = G(rho) W(phi)`, `G = -M Theta/rho`, `W = A/sin phi`, and in `R^5`

```
    rho^3 Lap_5 eta_0 / M  =  - W (Theta_uu + Theta_u - 2 Theta)  -  Theta ( W'' + 3 cot(phi) W' )
```
(checked against the bare plateau, where it returns `1/sin^3 phi`, i.e. `Lap_5(-M/r) = M/r^3`).
`A` is piecewise **linear**, so `A'' = 0` on each piece and a **Dirac mass at each kink**:
`phi = delta, pi/2 - delta_m, pi/2 + delta_m, pi - delta`. Hence

> **`Psi(0) = sup_x |x| |Lap_5 eta_0| = +infinity` for THEOREM_S3's datum.**
> Its absolutely continuous part is `171.4257` (units `M/rho_0^2`), attained at `u = 0.392`.

Mollifying the angular profile as `pmax-h3v` Proposition 3 does, namely `F_sigma = F * chi_sigma`
with `chi_sigma` a Gaussian of standard deviation `sigma`, `F` even-reflected at `phi = 0` and
`phi = pi`, makes it finite, with the closed form

```
    Psi_sigma(0)  ~  [W'] * max_u( Theta e^{-2u} ) / (sigma sqrt(2 pi))  =  C_kink/sigma ,
    [W'] = (1/delta)/sin delta = 58.52800        (the taper kink; the equator kink gives 11.503)
    max_u Theta e^{-2u} = 0.34565  at u = 0.38733
    C_kink = 8.07063  predicted ,   8.11753  measured at sigma = 5e-04 .
```

| `sigma` (rad) | `Psi_sigma(0)` | `\|grad eta_0\|_inf` | `Gfrak_0(sigma)` | `kappa_delta(sigma)` |
|---|---|---|---|---|
| kinked (a.c. part only) | `171.43` | `20.3339` | `58.6946` | `0.4978224182` |
| `0.05` | `244.07` | `9.1779` | `25.2888` | `0.4983882748` |
| `0.02` | `495.16` | `12.7523` | `36.5308` | `0.4979100241` |
| `0.01` | `901.42` | `15.2711` | `43.9502` | `0.4978442142` |
| `0.005` | `1709.06` | `17.1846` | `49.5547` | `0.4978278606` |
| `0.002` | `4130.11` | `18.7895` | `54.2132` | `0.4978232887` |
| `0.001` | `8165.29` | `19.4701` | `56.1885` | `0.4978226358` |
| `0.0005` | `16235.07` | `19.8625` | `57.3270` | `0.4978224726` |

Two things worth carrying out of that table.

* **A free lunch, recorded not claimed.** Mollifying the angle at `sigma = 0.05` rad (`2.9°`)
  divides `Gfrak_0` by `2.32` (`58.69 -> 25.29`) **and** raises `kappa_delta` from `0.4978224`
  to `0.4983883`. Both moves are favourable: `Gfrak_0` is the binding datum constant of
  `Ghat`, `C''`, `C_R` and hence `L_*` (THEOREM_S3 finding 1), and a larger `kappa` lowers both
  `c_* = log(3/2)/kappa` and the clock floor `(1-2 kappa)/(2 kappa)`. This is the *opposite* of
  the trade THEOREM_S3 M5 describes for the fully-`tanh` datum (D-C), which buys
  `Gfrak_0 = 20.92` at the price of `kappa = 0.47473`; the difference is that (D-C) also widens
  the equatorial layer from `delta_m = 5°` to `w = 0.20` in `cos phi`, and that is what costs
  `kappa`. **A rounding of the two corners at fixed `delta` and `delta_m` costs neither.**
  Whether to change the datum is Alexander's call; this sheet only prices it.
* **A number that does not reproduce.** `pmax-h3v` §A7's table reports
  `||grad eta_0||_inf rho_0^2/M = 31.093556949918312` ("not in the record"). Two independent
  instruments here, my own `(u, phi)` scan and `fix2`'s `DatumS3.grad_norm` on an `(r,z)`
  grid, give `20.3339` and `20.3261`. The closed form supports the smaller value:
  `rho^2|grad eta_0|/M = sqrt(W^2(Theta-Theta_u)^2 + W'^2 Theta^2)`, whose supremum is
  `Gfrak_0 = 1/sin^2 delta = 58.6955` (agreed by all three), and `||grad eta_0||_inf` is the
  supremum of the same expression times `e^{-2u}`, which the taper kink attains at
  `u = 0.368`. `31.09` would need `e^{-2u} = 0.53` **and** the full `58.70` at the same `u`,
  and those are attained `0.36` apart in `u`. Recorded as a discrepancy against `pmax-h3v`;
  nothing in the chain depends on it (it appears in no bound).

### 1.7 The propagation (`g5`), and the size of `eps_a` (`g1` §D)

`Psi <= sqrt(5) Psi_2`, `Psi_2(s) := sup rho |grad^2 eta|_F`. Differentiating the equation
twice (`pmax-h3v` Lemma D) and running the same weighted maximum principle with `k = 1`,

```
   Psi_2(s)  <=  e^{c_R + 2 c_Gamma} [ Psi_2(0)  +  s K_2 sqrt( G_inf N_2 ) ] ,
   G_inf(s) = ||grad eta||_inf <= e^{c_Gamma} · 20.334 M/rho_0^2 ,   N_2(s) <= e^{p c_G} Gfrak_0 M ,
```
the source being closed by `rho K_2 |grad eta| <= K_2 min(rho G_inf, N_2/rho)
<= K_2 sqrt(G_inf N_2)`, and neither weighted bound alone works (one blows up at the origin, the
other at infinity). `g5`: `tau K_2 sqrt(G_inf N_2)` is `5.25` at `L = 1e4` (frozen window),
`0.0525` at `L = 1e6`, against `Psi_2(0) = 1847` at `sigma = 0.002`: the propagation is
`0.28 %` at `L = 1e4` and `2.8e-05` at `L = 1e6`. **The propagation is not the binding half of
(H-Delta-eta); the datum is.**

`eps_a(L)`, proved column, `f >= 4`, `c_G` read from `fix2`'s `f1_window.assemble_c` at that
window (`L = 1e4` has no `(Gamma-off)` root in the proved column, `L_Gamma^exist = 14259.9`,
so `c_G` there is the measured column's):

| `L` | window | `c_G` | `sigma = 0.05` | `sigma = 0.01` | `sigma = 0.002` | `sigma = 0.001` | kinked, `sqrt` route | `eps(L)` of the budget |
|---|---|---|---|---|---|---|---|---|
| `1e4` | frozen `c_*` | `1.22267` | `4.31e-04` | `1.59e-03` | `7.30e-03` | `1.44e-02` | `3.44e-02` | `0.012686` (meas.) |
| `1e4` | cap `1.19437 c_*` | `1.46032` | `6.53e-04` | `2.41e-03` | `1.11e-02` | `2.19e-02` | `4.77e-02` | `0.023813` (meas.) |
| `1e5` | frozen `c_*` | `1.23424` | `4.36e-05` | `1.61e-04` | `7.38e-04` | `1.46e-03` | `1.10e-02` | `inf` (proved) |
| `1e5` | cap | `1.48377` | `6.69e-05` | `2.47e-04` | `1.13e-03` | `2.24e-03` | `1.54e-02` | `inf` (proved) |
| **`1e6`** | frozen `c_*` | `1.22291` | **`4.31e-06`** | **`1.59e-05`** | **`7.30e-05`** | **`1.44e-04`** | **`3.44e-03`** | `0.106218` (proved) |
| `1e6` | cap | `1.46141` | `6.54e-06` | `2.42e-05` | `1.11e-04` | `2.19e-04` | `4.77e-03` | `0.541568` (proved) |

**Relative to `fix2`'s `eps(L)`**: at `L = 1e6`, frozen window, proved column, `eps_a` is
`0.0041 %` of `eps` at `sigma = 0.05`, `0.069 %` at `sigma = 0.002`, and `3.24 %` on the
kinked `sqrt` route. **`eps_a` is not a load-bearing term in any of these regimes.**

The "kinked, `sqrt` route" column is the alternative if the datum's kinks are *kept*: then
`Lap_5 eta_0` is a surface measure, parabolic smoothing gives `Psi(r) ~ C [A'] rho_0/sqrt(nu r)`
and `nu int_0^s Psi = 2 C [A'] sqrt(nu s)`, so

```
    eps_a  =  (3/4)(1 + 1/(4L)) C_kink' [A'] e^{c_G} sqrt(c/L) / s ,
```
i.e. `O(sqrt(c/L))` rather than `O(c/(s^2 L))`. It is reported with `C_kink' = 1` and
`[A'] = 1/delta_m = 11.459`; the constant is **not** proved: it needs an explicit
parabolic-smoothing estimate with drift, which is not written here or anywhere in the chain.
That column is the reason this unit is ENCLOSED and not PROVED.

### 1.8 Status

> **UNIT 1: ENCLOSED.**
> PROVED: Lemma 1 (the weight `3/8` per e-fold, exact), Lemma 2 (the transport-with-source
> weighted sup), Lemma 3 (the harmonic identity, and with it the statement that pure diffusion
> costs `O(1/L^2)` of the drive), the propagation lemma with its constants, and the theorem of
> §1.4 conditionally on (H-Delta-eta).
> OPEN: (H-Delta-eta) itself. For THEOREM_S3's datum as written it is **false**,
> `Psi(0) = +infinity`, and it becomes true, with `Psi_sigma(0) = 8.07/sigma`, exactly when
> M5's angular kinks are mollified. **M3 and M5 are one item, not two.**
> The `sqrt` route that avoids mollification needs a parabolic-smoothing lemma nobody has
> written, and gives `eps_a = 3.44e-03` at `L = 1e6` rather than `7.30e-05`.

---

## 2. UNIT 2, (Theta-tail), item M7

### 2.1 What (D1) is used for, and what it is not used for

`hk2` §2 (D1): *`|omega^theta(.,t)| <= M_t <= (3/2)M` and `eta(.,t)` is supported in
`{rho_-(t) <= rho <= rho_+(t)}`.* Reading `hk2` §5's proof: Lemma 5.2's four integrals are

```
  int_{rho'>2rho} |x-x'|^{-4} d|D eta| <= 16 M 𝒥 int_{2rho}^{infinity} rho'^{-2}drho' = 8 M 𝒥/rho ,
  int_{rho'<rho/2} |x-x'|^{-4} d|D eta| <= 16 rho^{-4} M 𝒥 (rho/2)^3/3 = (2/3) M 𝒥/rho ,
```
and their `Lambda_5` mirrors. **They already run over `(2 rho, infinity)` and `(0, rho/2)`.**
The support hypothesis is never invoked in their proof; what is invoked is that the shell
density (2.1) holds with one finite `𝒥` at every `rho'`. And `hk2`'s numerical half never used
the support either: `k3.lambda_bulk` integrates `|grad eta|` over the actual field out to
`s = 3R`, tails included.

> **(D1′), the corrected hypothesis.**  `|omega^theta(.,t)| <= M_t <= (3/2)M` on `R^5`, and
> ```
>    |D eta|( {rho' < |x'| < rho' + drho'} )  <=  M rho'^2 𝒥 drho'     for every rho' in (0, infinity)
> ```
> with one finite `𝒥`.  No support statement is used anywhere in `hk2` §§4–5.
> `(D3)`(b) (the radial edge mollified over a width `w_0 >= kappa_0 rho_0`) already covers the
> `tanh` ramp, and the far/near lemma's Propositions 1 and 2 need only `|omega^theta| <= M` on
> `{rho' >= 2 rho}` and `{rho' <= rho/2}`, which holds on all of `R^5`.

### 2.2 The tails, priced with the actual `psi_in` (`g2`)

`Theta(u) = (1/2)[tanh((u-eps_r)/eps_r) - tanh((u-L+eps_r)/eps_r)]`, `u = log(rho/rho_0)`,
`eps_r = 1/4`; `psi_in = e^2/(e^2 + (rho/rho_0)^8)` and `Theta = 1 - psi_in` near the inner end.

> **(a) The e-fold count, EXACTLY.**  `int_{-infinity}^{infinity} Theta(u) du = L - 2 eps_r`,
> and each tail is
> ```
>     int_{-infinity}^{0} Theta du  =  (eps_r/2) log(1 + e^{-2 eps_r/eps_r})
>                                   =  log(1 + e^{-2})/8  =  0.0158660 ,
> ```
> the outer tail the same.  Numeric control at `L = 12`: whole line `11.500000000000005`
> against `L - 0.5 = 11.5`; inner tail `0.015865584`; shell only `11.468267` against
> `L - 0.5 - 2(0.0158660) = 11.4682680`.

This is the `C e^{-2(...)}` the brief expected, in closed form: **the `e^{-2}` is exactly the
`e^2` of `psi_in`**, and the `1/8` is `eps_r/2` with `eps_r = 1/4`.

> **(b) The local shell density in the tails.**  With `Theta_u = 8 Theta + O(Theta^2)` inside
> and `-8 Theta + O(Theta^2)` outside,
> ```
>    𝒥(rho')  <=  Theta(rho') * 𝒥coef_in  = 236.512 Theta   (u < 0) ,
>    𝒥(rho')  <=  Theta(rho') * 𝒥coef_out = 292.662 Theta   (u > L) ,
>    Theta <= e^{-2}(rho/rho_0)^8  (u<0) ,   Theta <= e^{-2}(R/rho)^8  (u>L)   EXACTLY
> ```
> (`Theta = sigma(8u-2)` at the inner end and `sigma(x) <= e^x`; numeric ratio `1.0000005`).
> The envelope is sharp: at `u = -1` and `u = L+1` the measured `𝒥` is `0.99991` of it.

Measured (`g2` §B, `L = 10` and `L = 40`, both):

| | `𝒥` |
|---|---|
| supremum over `rho'`, what (2.1) asks | `96.2602`, attained at `u = L - 0.336`, in the **outer ramp** |
| plateau (`Theta ≡ 1`) | `75.2511` |
| at `u = 0` (the inner shell edge) | `25.0233` |
| at `u = L` (the outer shell edge) | `31.6880` |
| supremum over the whole inner tail `u < 0` | `25.0233` |
| supremum over the whole outer tail `u > L` | `31.6880` |
| at `u = -1` / `u = L+1` | `0.010737` / `0.013286` |

**The tails do not raise `𝒥`.** The supremum lives in the outer ramp, a factor `3.85` above the
inner-tail supremum and `3.04` above the outer-tail one, and `𝒥` decays like `e^{-8|u|}` beyond.
`(D1')`'s `𝒥` is `96.2602`, the same number `(D1)` would have had, because the sup was never
in the tails.

> **(c) `Lambda_4` and `Lambda_5`.**  Analytically, the tails' share of `hk2`'s own uniform-`𝒥`
> majorants over the same regions is
> ```
>    inner:  int_0^{rho_0} 𝒥(rho') rho'^2 drho'  <=  𝒥coef_in e^{-2} rho_0^3/11  =  2.9099 ,
>            against  𝒥 (rho_*/2)^3/3 = 501.36  at f = 4   ->   0.58 %   (9.07 % at f = 1) ,
>    outer:  int_R^{infinity} 𝒥(rho') rho'^{-2} drho' <= 𝒥coef_out e^{-2}/(9R) = 2.0e-04 ,
>            against  𝒥/(2 rho_*) = 9.626   ->   2.1e-05 = (2/9) e^{-2} (𝒥coef_out/𝒥)(rho_*/R) ,
> ```
> i.e. `O(e^{-L})` at the outer end.  These are **shares of an existing majorant, not additions
> to it**: `hk2` integrates `rho'` from `0` and to `infinity` already.
>
> Measured, with `hk2`'s own `k3.lambda_bulk` run twice, once on the actual field and once with
> the field masked outside the transported shell, at `f = 4`, `L = 10`, the tracked point:
> ```
>    lam = 1.0000 :  Lambda_4  26.574429 full  vs  26.570094 masked   ->  tails 1.63e-04
>                    Lambda_5   7.192621 full  vs   7.191770 masked   ->  tails 1.18e-04
>    lam = 1.8420 :  Lambda_4  47.470194 full  vs  47.470145 masked   ->  tails 1.04e-06
>                    Lambda_5  11.875768 full  vs  11.875764 masked   ->  tails 3.46e-07
> ```

### 2.3 The corrected constants

Nothing moves. `𝒥 = 96.2602` (the sup, tails included); `Lambda_4`, `Lambda_5` and hence the
`K̂_2` table change by `<= 1.7e-04` relative if the tails are removed, and `hk2`'s quadrature
already includes them; the far/near constants `C1 = 0.291999`, `C2in = 0.014754`,
`C_collar = 3.999218/sin phi + pi/8` are unchanged, because their hypotheses are pointwise
bounds on `|omega^theta|` that hold on all of `R^5`.

**One correction is owed, and it is not a tail.** The exact e-fold content of the datum is
`L - 2 eps_r = L - 0.5`; `t3_budget.py:65` sets `ELL_RAMP = 0.25` ("e-folds lost at the outer
tanh edge"), i.e. it charges one edge. The missing `0.25` e-folds enter `eps_a` as `0.25/L`:
`2.5e-05` at `L = 1e4`, `2.5e-06` at `1e5`, `2.5e-07` at `1e6`, below `eps_Tprime` by four to
five orders at every `L` the theorem uses, and conservative to fix. Recorded, not repaired
(repairing it means editing `t3_budget.py`, which is out of scope here).

### 2.4 Status

> **UNIT 2: PROVED.**  `(D1')` is the correct hypothesis and `hk2` §5 needs nothing more;
> the tail contributions are `<= 1.7e-04` relative on `Lambda_4`, `Lambda_5` and zero on `𝒥`
> (the supremum is not in the tails); the far/near constants are unchanged; the e-fold
> correction is `+log(1+e^{-2})/8 = +0.0158660` per end, exactly, with the drive's sign.
> M7 can be struck from the modulo list, replaced by the `(D1')` clause.
> Owed to the budget, separately: `ELL_RAMP` should be `2 eps_r = 0.5`, not `0.25`.

---

## 3. UNIT 3, the `𝒥` mismatch, and the transported measurement

### 3.1 `𝒥` for the three data, one instrument (`g3` §A)

`hk2` Lemma 5.1: `|D eta|({rho'<|x'|<rho'+drho'}) = M rho'^2 𝒥(rho') drho'` with
`𝒥(rho') = |S^3| int_0^pi rho'^2 |grad eta|(rho',phi) sin^3 phi dphi / M` plus the jump terms.
The instrument written here takes `rho'^2|grad eta|` from each datum object's **own** derivative
fields, so one instrument serves all three profiles.

| datum | `𝒥` plateau (`Theta ≡ 1`) | `𝒥` supremum over `rho'`, what (2.1) asks | reported |
|---|---|---|---|
| **(D-A)** `delta = 7.5°`, sharp `sgn(z)` | `39.162759` a.c. `+ 39.478418` jump `= 78.641176` | `+infinity` at the two sharp radial edges (carried by `hk2` (D3)(a) as surface terms, not by `𝒥`) | `hk2` §8(a) `78.641147`, rel `3.7e-07` |
| **(D-B)** `delta = 7.5°`, `w = 0.20`, `w_0 = 0.25 rho_0`, `w_1 = 0.10 R` | `65.625910` | **`148.828`**, at `u = 9.9954` (`rho ≈ R`), in **its own outer ramp** | `hk2` §8(a) `65.625910`, rel `3.0e-09` |
| **(S3)** THEOREM_S3's datum | `75.26229` | `96.25539`, at `u = L - 0.34`, in its outer ramp | `fix2` §2.2 `75.251104` / `96.249844`, rel `1.5e-04` / `5.8e-05` |

*(`phi` grid `40001`; the `𝒥` values are quoted to five digits. `g2` §B, on a `400001` grid,
gives the (S3) supremum as `96.2602`, so the fourth digit is grid, not instrument.)*

> **THE RECONCILIATION.** `96.25` and `65.63` are not the same quantity.
> `96.25` is a **supremum over `rho'`**, which is what hypothesis (2.1) asks for; `65.63` and
> `78.64` are **plateau values**, quoted by `hk2` §8(a) as reference rows. On like for like:
> * supremum against supremum, `(S3) 96.26` against `(D-B) 148.83`, so **the theorem's datum is
>   `1.55×` better**, because its outer ramp is `eps_r = 0.25` wide in `log rho` while `hk2`'s
>   is `w_1 = 0.10 R` wide in `rho`, i.e. `0.10` wide in `log rho`, which is sharper;
> * plateau against plateau, `75.26` against `65.63`, so the theorem's datum is `1.15×` worse,
>   not `1.47×`.
>
> `fix2`'s `1.47×` and `1.22×` compared a supremum with two plateau values. The direction of
> the finding survives (the theorem's datum has a larger plateau density than `hk2`'s, because
> its `min(1,·)` equatorial mollifier at `delta_m = 5°` is steeper than `hk2`'s `w = 0.20`);
> the size does not.
>
> The same reading applies to `hk2`'s own sheet: **`hk2` §8(a)'s `𝒥` rows are plateau values,
> and (2.1) as stated needs the supremum, `148.83` for (D-B).** Nothing downstream moves
> (see §3.2), but the row is mislabelled in exactly the way `fix2` caught for the theorem's datum.

### 3.2 Where `𝒥` is used: nowhere that produces a constant (`g3` §B)

A file-by-file count of any shell-density symbol in the scripts that produce constants:

```
   hk2/k3_bound.py  0     hk2/k4_direct.py 0     hk2/k5_consequence.py 0    hk2/k6_confront.py 0
   hk2/hk2lib.py    0     THEOREM_S3/t2_gamma_CR.py 0     t3_budget.py 0    fix2/f1_window.py 0
```

`hk2`'s **PROVED** `K̂_2` table (§8(c)) comes from `k3.lambda_bulk`, a direct quadrature of
`int_{|x-x'|>d}|x-x'|^{-p}|grad eta|dx'` on the actual field; `C''` comes from
`t2_gamma_CR.gamma_off`, whose inputs are `E_0`, `Gfrak_0`, `C_far`, `C_inner`, `C_collar`.
**Neither reads `𝒥`.** `𝒥` appears only in Lemma 5.2's analytic statement, the *qualitative*
"no `log(R/rho_0)`" result, and in the reference table of §8(a).

So: **`fix2`'s K_2 table is unchanged, and `L_*` is unchanged.** Re-running a copy of
`fix2/f2_ramp_k2.py`'s instrument (§3.3) and of `fix2/f1_window.py` confirms it. The `L_*`
re-run reproduces `fix2`'s headline to the digit:

| `C_K` fed to `f1_window.assemble_c` | `L_*` (`eps <= 0.1943662`) | `L_*` (`eps <= 0.1`) |
|---|---|---|
| `161.7735` (the budget's import; `fix2` ADDENDUM_2's row) | `1887786.7127` | `1977309.3012` |
| `33.710` (`fix2` §2.3's quadrature value, S3 datum, `lam_max`, `f = 4`) | `1887786.6179` | `1977309.2020` |
| `102.6625` (the analytic route at `𝒥 = 96.26`, §3.3) | `1887786.7127` | `1977309.3012` |
| `176.4010` (the analytic route at `𝒥 = 192.52`, i.e. `𝒥` **doubled**) | `1887786.7127` | `1977309.4004` |

`fix2` ADDENDUM_2 reports `1887786.7` and `1977309.3`: **reproduced.** A `5.2×` range in `K̂_2`
moves `L_*` by `9.5e-08` relative, and `eps(1e7, proved) = 0.014141132567409676` is identical
to sixteen digits at every one of the four `C_K`: `fix2`'s `eps(1e7)` row, reproduced.

### 3.3 What `𝒥` would cost if the analytic route were used (`g3` §C)

The one place `𝒥` is load-bearing is the fully analytic `Lambda`: Lemma 5.2's far and near
bounds in `𝒥`, plus the collar `rho/2 < rho' < 2 rho` by quadrature. At `f = 4`, `L = 10`,
`phi_0 = 30°`, THEOREM_S3's datum, `d` optimised as `fix2` does:

```
    J = 96.2602 : lam=1.0000  L4=186.727 (collar 19.876)  L5=29.256   K2hat = 66.056
                  lam=1.8420  L4=205.488 (collar 30.906)  L5=32.146   K2hat = 102.663
    J = 78.6411 : lam=1.8420                                          K2hat =  89.166
    J = 65.6259 : lam=1.8420                                          K2hat =  79.196
    J = 192.520 : lam=1.8420  L4=380.069                L5=54.628     K2hat = 176.401
    (fix2's DIRECT-QUADRATURE value at the same point: K2hat = 33.710)
```

so `dlog K_2/dlog 𝒥 = 0.781` at `𝒥 = 96.26`, and the analytic route is `3.05×` weaker than the
quadrature one. Even so, doubling `𝒥` moves `L_*` by `0.1` in `1.98e+06`.

> **`𝒥` is worth nothing in `L_*`, by a margin of seven orders of magnitude.** `fix2`'s
> "the number to watch if the `K_2` margin is ever spent" is correct as a caution and is not a
> live cost: `(H-K2)` is already worth `1.4e-06` of `L_*` (THEOREM_S3 §4.6), and `𝒥` sits one
> level below `K_2`.

### 3.4 The transported measurement, and the refutation of `fix2`'s item (ii) (`g4`)

`fix2` ADDENDUM_2, "Also open": *"`hk2`'s quoted measured slack over the window compares a
transported bound with an untransported measurement (`hk2` §14 item 7 declares the untransported
evaluation wrong for the bound but the measured column still uses it)."*

**That is wrong.** `hk2` §9's confrontation table says, in its own parenthesis, *"measurement by
`k9` on the transported field; `k4`'s zonal series is the cross-check at `lambda = 1`"*. What
`fix2` read is `hk2/k4_results.json`'s `measured` block, which is the **withdrawn**
untransported `k4` run, the one §14 item 7 retires, not the table `hk2` §9 publishes.

Settled here with a **third** instrument, written from the mathematics and sharing no code with
either `k4`'s zonal series or `k9`: a direct 5-D kernel quadrature that puts every derivative on
`eta` and none on the kernel,

```
   a      = int K(w) eta(x-w) dw ,                    K(w) = 3 w_z/(8 pi^2|w|^5)
   a_r    = int K eta_r mu ,   a_z = int K eta_z ,    mu = yhat.yhat' = (r - s sinT cos X)/r'
   a_rr   = int K [ eta_rr mu^2 + (eta_r/r')(1-mu^2) ] ,   a_rz = int K eta_rz mu ,  a_zz = int K eta_zz
   a_perp = int K [ eta_rr q + (eta_r/r')(1-q) ] ,    q = s^2 sin^2 T sin^2 X/(3 r'^2)
```
absolutely convergent for the real-analytic (D-B), and **needing no separability**, so it works
on `eta_lam = eta_0 o T_lam^{-1}` at every `lambda`. Two controls that do not mention the
answer: `a_perp` must equal `a_r/r` (axisymmetry of `a`), and `Lap_5 a = a_rr + 3a_r/r + a_zz`
must equal `d_z eta(x)` (the PDE). At the finest grid `(ns,nT,nc) = (1400,260,260)` they hold to
`3.3e-16 … 7.7e-15` and `5.4e-08 … 1.4e-05` respectively.

| `lambda` | `a` | `\|grad a\|` | `\|\|Hess a\|\|` | `K̂_2` (transported, here) | `hk2` §9's measured | `hk2` bound | slack |
|---|---|---|---|---|---|---|---|
| `1.00` | `4.651525` | `0.516404` | `1.554688` | **`2.190342`** | `2.1903` | `66.6622` | `30.43` |
| `1.25` | `5.993196` | `0.637229` | `3.662912` | **`3.382506`** | `3.3827` | `107.6738` | `31.83` |
| `1.50` | `6.922827` | `0.728751` | `6.247641` | **`5.385261`** | `5.3854` | `161.7735` | `30.04` |

**Every entry of `hk2` §9's measured column is reproduced to five or six digits by an
independent instrument, on the transported field.** The withdrawn `k4` numbers, for contrast,
are `2.190342 / 1.543781 / 1.477726`: they agree at `lambda = 1` (where transport is the
identity) and are `2.19×` and `3.64×` too small at `lambda = 1.25` and `3/2`, which is exactly
the error §14 item 7 describes. The slack quoted by `hk2`, `30.4×`, `31.8×`, `30.0×`, is
transported against transported and **stands**; `t3_budget.py`'s
`K2_measured_window = 5.3854` is the correct transported value at `lambda = 3/2`.

**A side find in `hk2`.** `hk2` §9's control table says *"`a(x)` vs a direct 5-D quadrature
… `4.651534585` vs `4.651534585`, relative `~1e-11`"*, but `hk2/k4_results.json`'s
`control_kernel` row records `a_series = 4.651534585`, `a_kernel = 0.730625060`,
`rel = 0.843`. The instrument here gives `4.651525` at three grid refinements, agreeing with
the series and with `hk2` §9's prose: **the PROOF's claim is right and the stored JSON row is
the artefact.** It is a control, used in no bound; recorded because a displayed "passed at
`1e-11`" and a stored `0.84` cannot both be the same run.

### 3.5 Status

> **UNIT 3(a): SETTLED.** The three `𝒥` values are three different quantities of three
> different data, one supremum and two plateau values, and the like-for-like comparisons are
> `96.26` against `148.83` (supremum, the theorem's datum better by `1.55×`) and `75.26`
> against `65.63` (plateau, worse by `1.15×`). `𝒥` enters no constant in the chain: `K_2` is a
> direct quadrature and `C''` is a fixed point in `E_0`, `Gfrak_0` and the far/near constants.
> `fix2`'s `K_2` table and `L_*` are unchanged, and the re-run reproduces `1887786.7127` /
> `1977309.3012`. If the analytic route were used instead, `dlog K_2/dlog 𝒥 = 0.781` and even a
> doubled `𝒥` moves `L_*` by `5e-08` relative.
>
> **UNIT 3(b): REFUTED.** `fix2`'s open item (ii) rests on a misreading of
> `hk2/k4_results.json`; `hk2` §9's measured column is the transported field, and a third
> instrument reproduces it at `lambda = 1, 1.25, 3/2` to five to six digits. The item is
> closed, not repaired.

---

## 4. WHAT THIS SHEET CANNOT CATCH

It closes M7, encloses M3 (down to M5), and settles two open items of `fix2`. It does **not**
touch M1 (discharged in ADDENDUM_1), M2 (discharged), M4, M5, M6 or M8, it does not check that
`C_R` is the right constant or that `u2`'s P1 and P2 are the right instruments, and it does not
check that the modulo list is complete. Three of its own numbers are measurements
(`Psi_sigma(0)`, the `𝒥` values, the transported `K̂_2`): numerics falsify, they never prove.
The `sqrt` route of §1.7 is a *scaling*, not a bound: its constant is set to `1`.

**FL-000 stands. Nothing here touches the headline problem.**


---

## GATE, AND FILES

```
208 CHECKS, 208 PASS, 0 FAIL          (check_fix3.py)
```

`check_fix3.py` re-reads every displayed constant from the stored JSONs, checks the string
actually appears in `FIX3.md` or `ADDENDUM_3_2026-09-08.md`, and **re-derives** every ratio and
factor quoted between two numbers (`3.85×`, `3.04×`, `1.55×`, `1.15×`, `3.05×`, `0.069 %`,
`3.24 %`, the three slacks, `C_kink`'s closed form, `log Λ_*`) instead of copying it. It also
checks the scalings that carry the argument: `ε_a` falling like `1/L` on the smooth route and
like `1/√L` on the kinked one, `Ψ_σ(0)σ` flat, `𝒥`'s supremum independent of `L`, and
`ε(10⁷)` identical at every `C_K`.

| file | what it establishes |
|---|---|
| `g1_vstrain.py` / `g1_results.json` / `g1_stdout.txt` | UNIT 1: the kernel `𝒦` and its exact `3/8` weight per e-fold; the plateau control `a[η_P](0) = (M/2)L`; the harmonic identity `∫𝒦Δ₅η = ∂_zη(0)`; the datum's `Ψ_σ(0)`, `𝔊₀(σ)`, `κ_δ(σ)` and `‖∇η₀‖_∞`; `ε_a(L)` at three `L` and two windows |
| `g2_theta_tail.py` / `g2_results.json` / `g2_stdout.txt` | UNIT 2: `∫Θ du = L − 2ε_r` exactly, the tail e-folds `log(1+e^{−2})/8`, the tail envelopes, `𝒥(ρ′)` across the whole line, and `Λ₄`, `Λ₅` with and without the tails through `hk2`'s own `k3.lambda_bulk` |
| `g3_jmismatch.py` / `g3_results.json` / `g3_stdout.txt` | UNIT 3(a): `𝒥` for (D-A), (D-B), (S3) on one instrument, with `hk2`'s two rows reproduced; the file-by-file count showing `𝒥` is read by no script; the analytic-`Λ` route priced against `𝒥`; `L_*` re-run at four `C_K` |
| `g4_transported.py` / `g4_results.json` / `g4_stdout.txt` | UNIT 3(b): the third instrument, a direct 5-D kernel quadrature with every derivative on `η`, measuring `a` and its first and second derivatives on the **transported** field at `λ = 1, 1.25, 3/2`, with the axisymmetry and PDE controls |
| `g5_propagate.py` / `g5_results.json` / `g5_stdout.txt` | UNIT 1: the propagation of `Ψ_2` over the window, and its size against `Ψ_2(0)` |
| `fix2copy/` | byte copies, with hashes, of `fix2`'s scripts and of `fix2/imported/` |
| `check_fix3.py` | the gate: 208 checks over the two documents, plus the eighteen copy hashes |
| `SHA256SUMS` | computed with `shasum -a 256`, never typed |

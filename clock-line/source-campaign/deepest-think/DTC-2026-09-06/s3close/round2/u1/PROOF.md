# The slaved reference map: the feedback exponent of the (S3) bootstrap is `O(1)`, not `e^{34}`

Seat `s3close/round2/u1`, sitting of 2026-09-08.  Continues `s3close/assembly`.
Laws: `TORMENT NEXUS/LAWS.md` (first 120 lines read).

Every number in this note came out of a script in **this** folder that I wrote and ran
(`u1`…`u6`, re-asserted by `check_constants.py`); `SHA256SUMS` is computed with `shasum`,
never typed.  Nothing outside this folder was written.

**Read (read-only) and declared** (L-14): `s3close/assembly/ASSEMBLY.md` in full and
`s3close/assembly/a2_budget.py` (copied — `u5_budget.py` is that file with three functions
replaced, and the copy's provenance hash is in its header);
`write/lemma-T-shell-dependent/PROOF.md`; `rebuild/far-near-kernel-lemma/NOTE.md`;
`lower/prove-lagrangian/NOTE.md`; `s3close/hk2/PROOF.md` §0–1;
`s3close/round2/refute-assembly-feedback/` (`r1_log.txt`, `r2_alt_reference.py` docstring and
header, `r3_log.txt`, `r4_algebra.py`, `r4_log.txt`).  **That last folder matters and is
declared loudly**: an earlier refuter seat had already written the docstring of the same idea
(its variants `slaved_min`, `slaved_honest`, `conjugated`) and had verified the conjugation
identity in sympy, but its `r2_log.txt` is empty and no `r2_results.json` exists — the run
never finished.  Its `r3_log.txt` contains `r_h` on the extended range and its `r4_log.txt`
the algebra.  I wrote every instrument here from scratch and reproduced both:
`r_h[1,1.8377] = 0.9199394406311712` here against their `0.9199394381615896`
(rel `2.6845·10^{-9}`), and every one of `r4`'s exact residuals independently (`u1`).
No code was imported from that folder or from any other seat.

---

## 0. Headline

The assembly's one-window bootstrap has feedback exponent `pc = 34.06` because its reference
map carries an **a priori** strain rate, so the difference between the true `l = 1` strain and
the model's is a driver of the map error, and Lemma T′ prices that difference in terms of the
map error itself.  That loop is not intrinsic.  **Slave the reference to the true field** —
define the reference rate to *be* the true solution's exterior `l = 1` strain — and the loop
term vanishes identically, by definition rather than by estimate.

```
                    feedback exponent  p c  at c = c_* = 0.8113826,  L -> infinity
  ASSEMBLY (a priori reference), proved constants        34.0573      e^{pc} = 6.178e14
  ASSEMBLY, record's own measured Lemma T' sensitivity    6.6301      e^{pc} = 757.6
  THIS NOTE, slaved reference, Groenwall form (2.3)       1.2170740   e^{pc} = 3.3772909
                                          at L = 640     1.4106013   e^{pc} = 4.0984191
  THIS NOTE, slaved reference, conjugated form (2.6)      0.3468 (L=640, proved constants)
                                                         0.0412 (L=640, measured constants)
```

`pc` falls from `34.0573` to `1.2170740` (`= (3/2)c_*` exactly, in the limit), and `e^{pc}`
from `6.1784104·10^{14}` to `3.3772909`: a factor `1.8293983·10^{14}`.  The consequences, in the three
columns the brief names (§6):

| column | `pc` (`L=640`) | `μ(c)·L` | `ε(640)` | `ε(2560)` | `L_*` (`ε≤½`) | `log Λ_*` | `L_*` (`ε≤0.1`) |
|---|---|---|---|---|---|---|---|
| all proved (`C″=151.15`, `C_a=8.697888`, `K₂=161.7735`) | `1.4106013` | `446.2052` | `∞` | `∞` | `16161.986` | `32327.336` | `53783.018` |
| proved, `C_a` measured (`0.069`) | `1.4106013` | `64.25114` | `∞` | `0.4439660` | `2389.4697` | `4782.3038` | `7774.7586` |
| all measured (`C″=11.74`, sens `/6.0668`) | `1.2338594` | `64.25114` | `0.4307170` | `0.0530984` | `600.25760` | `1203.8795` | `1541.8189` |

against the assembly's `L_* = 4.5253·10^{14}` (proved) and `4.8315·10^{3}` (measured).  In the
all-proved column the bootstrap now closes **at a finite `L` with a name** —
`L_* = 1.6162·10^4` instead of `4.5253·10^{14}`, i.e. `log Λ_*` falls from `9.0506·10^{14}` to
`3.2327·10^4`, **ten orders of magnitude** — and it does so with *no restart lemma*, which is
the thing ASSEMBLY §3.4 priced as worth ~12 orders and BLOCK 2 as fatal.  **The restart
(BLOCK 2) is no longer the binding constraint.**  What binds now is the *size* of the proved
near-field constant `C_a = 8.697888` (the far/near collar), whose measured value is `0.069`:
replacing it alone takes `L_*` from `16162` to `2389`.

**What this note does NOT do.**  It does not prove (S3).  Two hypotheses are carried by name
and neither is proved anywhere: **(Γ-off)**, being proved by a parallel seat, and **(R-z)**,
the `z`-component of the velocity remainder, which is *inherited* from ASSEMBLY (§3.4 below is
a finding against ASSEMBLY: its (2.1) charges the whole **vector** velocity remainder to
`C_a`, which the far/near kernel lemma bounds only for the **scalar** strain `a`).
`FL-000 stands.  Nothing here touches the headline problem.`

### 0.1 Grade

**FILED.**  No Solid, no Major.  The content is: (a) the redefinition and the proof that the
`l = 1` rate error vanishes identically under it; (b) the re-derived majorant system, with the
Jacobian equation now closing against `λ²` directly so that Lemma T′ Corollary 3 and its
`2κ_s/L` are not needed either; (c) the identification of the exterior-strain functional with
the object Lemma T′ controls, as an **identity** (far/near kernel lemma Consequence A at
`x = 0`) rather than an estimate; (d) the re-priced budget; (e) three findings against the
record (§3.4, §5.4, §7 item 4).

---

## 1. The definition, and why it is the right one

### 1.1 The exterior `l = 1` strain

Notation is the 5-D lift of `rebuild/far-near-kernel-lemma` §0 and
`write/lemma-T-shell-dependent` §1: `x = (y,z) ∈ ℝ⁴×ℝ`, `r = |y|`, `ρ = |x|`, `φ` the polar
angle from `+z`, `η = ω^θ/r`, `−Δ₅ψ₁ = η`, `a = u^r/r = −∂_zψ₁`, `u^z = 2ψ₁ + r∂_rψ₁`, and
`b(y,z) := (a·y, u^z)` the lifted velocity, whose flow map is `Φ_s`.

The far/near kernel lemma's representation (its `(*)`) writes the strain as a log-integral over
shells of a scale-invariant shell functional, and its **Consequence A** says that the `l = 1`
interior mode of a shell is *exactly constant* on the ball the shell encloses, with value
`Φ[w](0) = −(3/5)g_1 = −(3/4)∫₀^π w(φ)cos φ sin²φ dφ`, `|Φ[w](0)| ≤ M/2` per e-fold.

> **DEFINITION (the slaved rate).**  For `ρ ∈ [ρ₀,R]` and `s ∈ [0,τ]` put
> ```
>   𝔞(ρ,s) := ∫_{2ρ}^{∞} F(ρ′,s) dlog ρ′ ,   F(ρ′,s) := −(3/5) g_1[ω^θ(·,s)](ρ′) ,
> ```
> `g_1[w](ρ′)` the `l = 1` Gegenbauer coefficient of `w(·,s)/sin` on the sphere `|X| = ρ′`.
> `𝔞(ρ,s)` is the uniform coaxial strain generated **at every point of the ball `{|X| < 2ρ}`**
> by the part of the true vorticity lying outside radius `2ρ`.  The octave split at `2ρ` is
> exactly the far/near lemma's own split.  Put
> ```
>   λ(ρ,s) := exp ∫₀^s 𝔞(ρ,σ) dσ ,        Λ_s(x) := T_{λ(|x|,s)} x ,
>   T_λ(y,z) = (λ y, λ^{-2} z) .
> ```

Three elementary facts, all proved (`u1_algebra.py`, sympy, residual `0` where marked EXACT):

* **(P1)** `𝔞 ≥ 0`, hence `λ ≥ 1`.  By Lemma 2 (`prove-lagrangian` §3) `η ≤ 0` on `{z>0}` for
  all `t` and all `ν ≥ 0`, and by `z`-oddness `ω^θ cos φ ≤ 0` everywhere, so the integrand of
  `Φ[w](0)` is pointwise `≥ 0`.  *This is the second load-bearing use of Lemma 2 in the
  campaign* (ASSEMBLY §0.1 records the first).
* **(P2)** `|∂𝔞/∂log ρ| = F(2ρ,s) ≤ M_s/2`, EXACT: differentiating the log-integral in its
  lower limit gives `−F(2ρ,s)` (`u1`: `dfrak_dlogrho = -F(U + log 2)`, numeric residual
  `1.45·10^{-10}`).  Hence, with `M_s ≤ (3/2)M` on the window,
  ```
     κ_s := sup_ρ L|ρ ∂_ρλ/λ| = sup_ρ L|∫₀^s ∂𝔞/∂log ρ dσ| ≤ (3/4) c = 0.6085369 ,
  ```
  which is Lemma T′'s hypothesis (S) with `κ_s = (3/4)c` in place of the integro-ODE profile's
  `√6−2 = 0.4494897`.  Step 0 of Lemma T′ needs `2κ_s/L < 1`, i.e. `L > 1.2170739`.
* **(P3)** `λ(ρ,s) ≤ λ_max := e^{3c/4} = 1.8377407` on the window, and `λ(·,s)` is decreasing
  in `ρ`.  (`𝔞 ≤ (M_s/2)L` and `τ = c/(ML)`.)  The a priori cap is *not* `3/2`: the reference
  is driven by the true strain, which the contradiction hypothesis bounds by `(3/4)ML`, not by
  `κ_δML`.  Lemma T′ tolerates this — its §2 says "no constraint on the range of `λ` is used
  in (T′) itself; `λ ∈ [1,3/2]` enters only in Step 0's condition and in Corollary 2's `r_h`" —
  at the price of recomputing `r_h` on `[1,λ_max]`: **`0.9199394406` instead of `0.9818223013`**
  (`u2`; the assembly's value is reproduced here to `3.0·10^{-7}`).

### 1.2 `Λ` is the flow of a genuine field, and that field is the `l = 1` mode

Differentiating `Λ_s(x) = T_{λ(|x|,s)}x` at fixed `x`, and using `dT_λ/ds = 𝔞 D T_λ` with
`D := diag(1,1,1,1,−2)`,
```
   ∂_sΛ_s(x) = 𝔞(|x|,s) · D Λ_s(x) ,        v(X,s) := 𝔞(|Λ_s^{-1}X|,s) · D X .
```
`u1` checks in sympy (residual `0`) that the `l = 1` interior mode `ψ₁ = −A z` has
`u^r = A r`, `u^z = −2Az`, i.e. `b = A·D X` **exactly**, that it is 3-D divergence free, and
that its 5-D divergence is `2A` — so `v` is precisely the velocity of the `l = 1` mode with
`A = 𝔞`, evaluated at the label radius.  The 5-D Jacobian of `Λ` is (EXACT, `u1`,
residual `0`)
```
   J_Λ(x) = λ² [ 1 + (ρλ′/λ)(1 − 3cos²φ) ] ,      1 − 3cos²φ ∈ [−2,1] .
```

---

## 2. Part (1): the `l = 1` rate error vanishes identically

### 2.1 The driver

`Φ_s` solves `∂_sΦ_s(x) = b(Φ_s(x),s)`.  At **any** point `X`, split the true velocity by the
same octave split that defines `𝔞`:
```
   b(X,s) = 𝔞(|X|,s) D X + R(X,s) ,      R := the remainder.                       (2.1)
```
The splitting is an identity, not an estimate: by Consequence A the shells outside `2|X|`
contribute exactly `𝔞(|X|,s)·DX` at `X`.  Therefore
```
 d/ds (Φ−Λ)(x) = b(Φ_s x) − 𝔞(|x|,s) D Λ_s(x)
               = [ b(Φ_s x) − b(Λ_s x) ]  +  R(Λ_s x, s)
                 +  [ 𝔞(|Λ_s x|,s) − 𝔞(|x|,s) ] D Λ_s(x) .                          (2.2)
```
**Three terms, and the `l = 1` rate error is not one of them.**  Compare ASSEMBLY (2.3), whose
second bracket splits as
```
 [u − v](Λ) = R(Λ) + [𝔞_true(|Λ|) − 𝔞_true(|x|)]DΛ + [𝔞_true(|x|) − 𝔞_model(|x|)]DΛ ,
```
the last of which is the `l = 1` rate error, priced by Lemma T′ at `ε_{T′}·κ_δML` and therefore
proportional to `μ`.  With `𝔞_model := 𝔞_true` that term is `0` **identically**, for every `x`
and every `s`, because the two objects are the same function.  Nothing is estimated away; the
term is not there.

### 2.2 The three surviving terms

* **Lipschitz.**  `|b(Φ) − b(Λ)| ≤ Γ(s)|Φ−Λ|` with `Γ(s) := ‖∇b(·,s)‖_{L^∞(slab)}`.  This is
  the only place the named hypothesis enters:
  > **(Γ-off)** *(carried, not proved here; a parallel seat is proving it)*
  > `‖∇u(s)‖_{L^∞(slab)} ≤ 2a(0,s)(1 + C₁/L) + C″M` for the true field.
  With `a(0,s) ≤ (M_s/2)L ≤ (3/4)ML` on the window,
  `Γ̄/(ML) =: G = (3/2)(1 + C₁/L) + C″/L`, and `c_G := ∫₀^τΓ ds = G c`.
  Numbers use `C₁ = 1`, `C″ = 151.15` (proved column) and `C″ = 11.74` (measured).
* **Near field.**  `|R(X,s)| ≤ C_R M_s |X|`; §3 derives `C_R`.
* **Label vs current radius.**  `|Λ_s x| / |x| = g(φ;λ) ∈ [λ^{-2}, λ]`, so by (P2)
  ```
    | 𝔞(|Λ|,s) − 𝔞(|x|,s) | ≤ (M_s/2)·|log(|Λ|/|x|)| ≤ (M_s/2)·2 log λ_max ,
  ```
  and `|DΛ| ≤ 2|Λ|`.  Cost: `2 M_s log λ_max |Λ|`.
  This is **strictly better than ASSEMBLY's `κ_δ M c_G`**: the assembly charges the *true
  flow's* radial travel `∫Γ ds = c_G`, which is `1.2212703` to `4.3135127` over the `L` range
  tabulated in §6.2; the slaved reference charges only its **own**
  radial travel, `2 log λ_max = 1.2170739`, which is `L`-independent and carries no `Γ`.

### 2.3 The majorant system

Write `m := sup_x|Φ_s(x) − Λ_s(x)|/|x|`, `θ := MLt`, `λ_ω := M_s/M ≤ 3/2`.  Dividing (2.2) by
`|x|` and using `|Λ| ≤ λ_max|x|`:
```
   dm/dθ  =  G·m  +  λ_max λ_ω ( C_R + 2 log λ_max ) / L ,     m(0) = 0 ,           (2.3)
   μ      =  λ_max² m           (since |Λ| ≥ λ_max^{-2}|x|) .
```
**The drive is constant in `m`.**  That is the entire content of the redefinition: (2.3) is a
linear inhomogeneous scalar ODE with an `O(1/L)` forcing, and its solution is
```
   m(c) = (drive/G)(e^{Gc} − 1) ,    μ(c) = λ_max² m(c) = O(1/L) ,
```
with **feedback exponent `pc = Gc`**, not `[3/2 + C′/L + (3/2)(2κ_δ)(15π/4)(9/4)/r_h]c`.

For the Jacobian, note that `∂_s log J_Φ^{(5)} = div₅b(Φ) = 2a(Φ)` (`u1`, EXACT for a general
axisymmetric no-swirl field, residual `0`) while `∂_s log λ² = 2𝔞(|x|,s)` **by definition of
`λ`**.  Hence, with `w := log(J_Φ/λ²)`,
```
   dw/dθ = (2/(ML))[ a(Φ_s x,s) − 𝔞(|x|,s) ]
         ≤ (2λ_ω/L)[ C_a + ½ Λ_rad(m) ] ,   Λ_rad(m) := |log(|Φ|/|x|)| ≤
                                            max( log(λ_max+m), log 1/(λ_max^{-2}−m) ) ,  (2.4)
   μ_J := sup|J_Φ/λ² − 1| ≤ e^{w(c)} − 1 .
```
This is **exactly Lemma T′'s hypothesis (H3)**, `|J_Φ − λ²| ≤ μ_Jλ²`, obtained directly.  So
Lemma T′ **Corollary 3's composition is not used**, and its `2κ_s/L` additive term is not paid.
ASSEMBLY (2.2) had, besides the `ε_{T′}` term, a `(Γ̄/(ML))μ` term coming from
`|a(Φ)−a(Λ)| ≤ ‖∇a‖|Φ−Λ|`; here `a(Φ)` is compared with `𝔞(|x|)` instead, and the difference
is bounded **absolutely** by `O(M)`, with no `Γ` and no `μ`.

### 2.4 The conjugated form (strictly better in structure, worse in these constants)

Write `Φ_s(x) = T_{λ(|x|,s)}(x + e(s))`.  Then, EXACTLY (`u1`, sympy, both residuals `0`;
independently obtained by `refute-assembly-feedback/r4`),
```
   e′ = ( 𝔞(|Φ|,s) − 𝔞(|x|,s) ) D (x+e)  +  T_λ^{-1} R(Φ) .                          (2.5)
```
**No Lipschitz constant of `u` appears at all**: the `O(ML)` part of the velocity is absorbed by
the reference's own time dependence, and what is left is the *variation* of the `l = 1` rate
over the radial travel, which is `O(M)` by (P2).  With `ê := sup|e|/|x|`, `μ ≤ λ_max³ê`,
```
   dê/dθ ≤ (λ_ω/L)[ 2 log λ_max + log(1/(1−ê)) + λ_max³ C_R ] (1+ê) ,               (2.6)
```
whose linearised exponent is `pc = (λ_ω/L)[1 + 2 log λ_max + λ_max³C_R]c = O(c/L)`:
`0.34678624` at `L = 640` in the proved column, `0.04124785` in the measured one.  **(Γ-off) is not
used in (2.6) at any point** — it survives in this note only inside the viscous budget's `c_G`.

At the constants actually in play the Grönwall form (2.3) gives the *smaller* `μ(c)`, because
`λ_max³ = 6.2065851` costs more than `e^{Gc} = 4.0984191` saves: `μ(c)L = 446.20524` (Grönwall) against
`1369.9556` (conjugated) in the proved column.  Both are reported; the budget uses (2.3).

---

## 3. `C_R`: the velocity remainder, and a finding against ASSEMBLY

### 3.1 The `r`-component

`R_y = (a(X) − 𝔞(|X|))·Y`, and `a(X) − 𝔞(|X|)` is exactly the far Taylor remainder plus the
inner multipole plus the collar of the far/near kernel lemma, with the split at `2|X|` — the
lemma's own split.  Its hypotheses are only `|ω^θ| ≤ M_s` and `z`-oddness, both true for the
**true** field at every time (the first by the contradiction hypothesis, the second by Lemma 2),
so
```
   |a(X,s) − 𝔞(|X|,s)| ≤ C_a(φ) M_s ,
   C_a(φ) = 0.291999 + 0.014754 + 3.999218/sin φ + π/8 ,   C_a(30°) = 8.697888082 .
```
(`u2`; the four constants are quoted from `far-near-kernel-lemma` §2–3, the `z`-odd forms.)

### 3.2 The `z`-component, and what it needs

`u^z` is **not** covered by the far/near kernel lemma.  It is bounded here as follows.
`η(·,s)` is odd in `z` (Lemma 2), hence `ψ₁` and `u^z` are odd in `z`, so `u^z(r,0) = 0`; and
3-D incompressibility gives the identity (`u1`, sympy, residual `0`)
```
   ∂_z u^z = −2a − r ∂_r a .
```
Integrating from the equator, and inserting `±𝔞` at the running radius `ρ_ζ = √(r²+ζ²) ≤ ρ`
(so `sin φ_ζ ≥ sin φ`, hence `C_a(φ_ζ) ≤ C_a(φ)`):
```
   R_z = u^z + 2𝔞(ρ)z
       = −∫₀^z 2[a − 𝔞(ρ_ζ)]dζ − ∫₀^z 2[𝔞(ρ_ζ) − 𝔞(ρ)]dζ − ∫₀^z r∂_r a dζ ,
   |R_z| ≤ [ 2C_a(φ) + log(1/sin φ) + G₁ ] M_s |z| ,     G₁ := sup |r∂_r a|/M_s .
```
`r∂_r a = sin²φ·(∂𝔞/∂log ρ) + r∂_r[a − 𝔞]`, so by (P2)
```
   G₁ ≤ ½ + G_far + G_inner + G_collar = ½ + 1.6049285 + 0.1324254 + G_collar
      = 2.2373539 + G_collar ,
```
where `G_far`, `G_inner` are the Parseval-based gradient constants quoted in ASSEMBLY §2.6
BLOCK 3, and `G_collar` — **the collar's gradient constant — is not proved anywhere**.  That is
the residue of BLOCK 3, in its weakest form: BLOCK 3 asks for `‖∇u‖_∞ ≲ ML`; this asks only for
`r∂_r[a−𝔞] = O(M)` on the collar `|X|/2 < ρ′ < 2|X|`.

> **NAMED HYPOTHESIS (R-z).**  `|r ∂_r[a − 𝔞(|·|)]|_{collar} ≤ G_collar·M_s` on the tube, with
> `G_collar` an absolute constant.  Status: **NOT PROVED**, inherited from BLOCK 3.

Assembling, `|R| ≤ |R_y| + |R_z| ≤ C_R M_s|X|` with
```
   C_R = 3C_a(φ) + log(1/sin φ) + ½ + G_far + G_inner + G_collar .
```
At `φ = φ₀ = 30°`, `G_collar = 0`: **`C_R = 29.024165`** (proved constants) and
**`C_R = 3.137501`** with the measured `C_a = 0.069`.  Sensitivity to the unproved `G_collar`
is reported in §6.5; because `C_R` enters (2.3) **linearly, in the drive and not in the
exponent**, a factor `5` in `G_collar` is a factor `2.4092307` in `L_*`, not `e^5`.

### 3.3 What the constants cost

| piece | proved | measured |
|---|---|---|
| `3C_a(30°)` | `26.0936642` | `0.207` |
| `log(1/sin 30°)` | `0.6931472` | `0.6931472` |
| `½ + G_far + G_inner` | `2.2373539` | `2.2373539` |
| `G_collar` | **not proved** | — |
| `C_R` | `29.024165` | `3.137501` |

### 3.4 FINDING (against ASSEMBLY §2.3–2.4)

ASSEMBLY's driver charges "the deviatoric remainder `O(Mρ)`, bounded by `C_aλM`" — i.e. it
charges the whole **vector** velocity remainder `[u−v](Λ)` to `C_a`.  But `C_a` is the far/near
kernel lemma's constant for the **scalar** strain `a = u^r/r` only; that lemma says nothing
about `u^z`.  The `z`-component needs the argument of §3.2, which costs a further `2C_a`, a
`log(1/sin φ)`, and a gradient constant whose collar part is exactly BLOCK 3's missing piece.
So ASSEMBLY (2.1) understates its own drive by a factor `C_R/C_a = 3.3369210` at `G_collar = 0`,
and **imports a hypothesis it does not name**.  This is a defect in the drive, not in the
exponent, so it does not change ASSEMBLY's headline; it does change its `L_*` by that factor.
Recorded for the ledger.

---

## 4. Part (2): Lemma T′, applied once

### 4.1 Is `𝔞` the same functional Lemma T′ controls?  YES, and it is an identity

This is the brief's own worry and it is the load-bearing question of the note.

Lemma T′ controls `a[η₀∘Φ^{-1}](0)` — the strain **at the origin** of the transported datum —
against `a_ref[η₀;λ] = ∫_S𝒦(Λx)η₀(x)λ(ρ)²dx₅`.  My `𝔞(ρ,s)` is the `l = 1` interior mode
generated by shells outside `2ρ`.  At `x = 0` **every** shell is outside, and Consequence A
says the value at the origin *is* the sum of the shells' `l = 1` interior modes:
```
    a(0,s) = ∫_0^∞ F(ρ′,s) dlog ρ′ ,        𝔞(ρ,s) = a(0,s) − ∫_0^{2ρ} F(ρ′,s) dlog ρ′ ,
```
an **identity** (`u1`: symbolic additivity, numeric residual `4.4·10^{-16}`), with the
subtracted piece non-negative by (P1) and at most `(M_s/2)log(2ρ/ρ₀)`.  So the two objects are
the same functional up to the inner-shell integral, and that integral is exactly the "lost
e-folds" the assembly already prices.  **They differ by the near-field octave, and the octave
is priced, in e-folds, not in an unbounded constant.**  Answer to the brief: same functional;
price `ℓ_loss/L`.

### 4.2 The chain

Fix `s`, let `x_* = (ρ_*,φ₀)`, `ρ_* = (1+f)ρ₀`, `X(s) = Φ_s(x_*)`.  Choose the **label cut**
```
   ρ_c := ρ₀ e^{ℓ_loss} ,   ℓ_loss := log 2 + log(1+f) + 3 log λ_max + log((1+μ)/(1−μ)) ,
```
so that every label `ρ″ ≥ ρ_c` has current radius `|Φ_s(x)| ≥ (1−μ)λ_max^{-2}ρ″ ≥ 2|X(s)|`
(using `|X(s)| ≤ (1+μ)λ_max ρ_*`).  Put `S^out := {ρ″ ≥ ρ_c}`.  Then:

1. **(near-field split)** `a(X(s),s) ≥ 𝔞(|X(s)|,s) − C_a M_s`  — §3.1.
2. **(positivity + restriction)** `F(·,s) ≥ F^{out}(·,s) ≥ 0` where `F^{out}` is the `l=1`
   coefficient of `ω^θ(·,s)·1_{Φ_s(S^out)}` (Lemma 2 again: every piece contributes with the
   same sign), and `F^{out}` is supported at current radii `> 2|X(s)|`.  Hence
   `𝔞(|X(s)|,s) ≥ ∫_0^∞F^{out} dlog ρ′ = a[η(·,s)1_{Φ_s(S^out)}](0)`.
3. **(transport)** `η(·,s)1_{Φ_s(S^out)} = (η₀1_{S^out})∘Φ_s^{-1} + (viscous defect)`.
4. **(Lemma T′, ONCE)** `η₀1_{S^out}` satisfies Lemma T′'s hypothesis (D); `Φ_s` satisfies
   (H1); (H2) with `μ(s)` and (H3) with `μ_J(s)` from §2.3; `λ(·,s) ∈ C¹` with `|ρλ′/λ| ≤ κ_s/L`
   by (P2).  So, by Corollary 2 with `r_h` on `[1,λ_max]`,
   ```
      a[(η₀1_{S^out})∘Φ_s^{-1}](0)  ≥  (M/2)∫_{ρ_c}^{R}P_{h}(λ(ρ″,s)) dlog ρ″ · (1 − ε_{T′}) ,
      ε_{T′}(μ,μ_J) = (2π/r_h)[ 3μ(1+μ_J)/(2(1−μ)⁵) + 3μ_J/8 ] .
   ```
5. **(viscous)** subtract the Theorem V.4 loss, sized as ASSEMBLY §2.5 with the V-b refuter's
   `87.166` and `K₂ = C_K M/ρ₀`.

Lemma T′ is used **once**, in step 4, and nowhere else.  Its `μ` and `μ_J` come from a system
that does not contain `ε_{T′}`.  That is the whole architecture.

> **CARRIED GAP (V-strain).**  Step 3's viscous defect is a defect of `η` *pointwise*
> (Theorem V.4); converting it into a defect of the *strain functional* `a[·](0)`, which is a
> weighted integral of `η`, is not done in the record and is not done here.  ASSEMBLY's own
> clause (iii) omits it entirely.  §6.4 reports the budget with the viscous loss charged to the
> strain as well as to the vorticity, so the reader can see what it would cost.

---

## 5. Part (3): closing the system

### 5.1 The integro-differential inequality, and the comparison principle

Applying §4.2 at a general label shell `ρ` (cut `ρ_c(ρ) = 2ρλ_max²/(1−μ)`, lost e-folds
`ℓ_shell = log 2 + 2 log λ_max + log(1/(1−μ)) = 1.9102211` at `μ = 0`) gives, in
`σ = log(ρ/ρ₀)/L` and `θ = MLt`, with `d := ℓ_shell/L`:
```
   ∂_θ log λ(σ,θ)  ≥  (1 − ε_{T′})·½∫_{σ+d}^{1} P_h(λ(σ′,θ))dσ′  −  C_visc/L ,
   λ(·,0) = 1 .                                                                      (5.1)
```
`P_h` is **increasing** on `[1, 2.075]` and the window caps `λ` at `λ_max = 1.8377`
(`u2`: `min` increment `8.49·10^{-4}` over `[1,λ_max]`; `u4`: the interpolation table is
monotone and agrees with exact quadrature to `1.74·10^{-7}`; `P_h` turns over only above
`λ ≈ 2.075`).  Hence the right-hand side is monotone in the profile.

> **COMPARISON.**  Discretise (5.1) on `N` shells: `ẋ_i = f_i(x)` with `∂f_i/∂x_j ≥ 0` for
> `j ≠ i` (quasimonotone), because `∂f_i/∂x_j = (1−ε)½P_h′(λ_j)λ_jΔσ·1[σ_j ≥ σ_i+d] ≥ 0`.
> *Proof (two lines).*  Let `x` satisfy `ẋ ≥ f(x)`, `y` solve `ẏ = f(y) − ε𝟙` for `ε>0`,
> `x(0) ≥ y(0)`.  If `x ≥ y` failed, let `θ₀` be the first time and `i` an index with
> `x_i(θ₀) = y_i(θ₀)`, `x_j(θ₀) ≥ y_j(θ₀)` for all `j`; then
> `ẋ_i(θ₀) ≥ f_i(x(θ₀)) ≥ f_i(y(θ₀)) = ẏ_i(θ₀) + ε > ẏ_i(θ₀)` by quasimonotonicity, so `x_i−y_i`
> is strictly increasing at `θ₀`, contradicting that `θ₀` is a first crossing.  Let `ε↓0`. ∎
> (The continuum statement is Walter, *Differential and Integral Inequalities*, Thm 10.XII.)

### 5.2 The numerics (`u4`, `N = 200` log-shells)

Integrated at `L = 40, 160, 640`, `ε_{T′} = 0`, `C_visc = 0`, against

* **(B)** the same system with `P_h(λ)` replaced by `r_hλ` (a lower bound, `P_h ≥ r_hλ`);
* **(C)** the conservative exponential `exp(κ_δ(1−d−σ₀)θ)` — what the theorem below claims;
* **(D)** `prove-lagrangian` §4's closed form `(1−(1−σ)κ_cθ/2)^{-2}` at the same reduced rate.

| `L` | `d` | `θ` at `λ(0)=3/2`, system (A) | (B) | (C) conservative | (D) closed form | `A≥B` | `A≥C` | `A≥D` |
|---|---|---|---|---|---|---|---|---|
| 40 | `0.0477555` | `0.7816769` | `0.8482487` | `0.8543167` | `0.7732850` | yes | yes | no |
| 160 | `0.0119389` | `0.7479297` | `0.8114996` | `0.8232697` | `0.7451828` | yes | yes | no |
| 640 | `0.0029847` | `0.7399101` | `0.8027661` | `0.8158573` | `0.7384735` | yes | yes | no |

The comparison holds: `λ_A ≥ λ_B` and `λ_A ≥ λ_C` pointwise (`min(A−C) = 0` to machine
precision at all three `L`), and the profile is decreasing in `σ` as (P3) requires.  `A ≥ D`
**fails** by `0.1945463 %` (`L=640`) to `1.0852265 %` (`L=40`) — correctly: (D) is the *bang–bang* closed form (`P_h(λ) = λ`), and the
`7.5°` taper makes `P_h(λ) < λ`.  So `prove-lagrangian`'s accelerated closed form is **not** a
lower bound for the tapered datum, and the defensible statement is the conservative one.

### 5.3 The clock

Conservatively (`P_h(λ) ≥ P_h(1) = 2κ_δ` for `λ ≥ 1`, `P_h` increasing):
```
   ∂_θ log λ(0,θ) ≥ κ_δ(1 − ℓ_loss/L)(1 − ε_{T′}) − C_visc/L ,
   θ_* = log(3/2)/κ_δ · (1+ε) = c_*(1+ε) ,   c_* = 0.8113825936 ,
   t_* = c_*(1+ε)/(ML) = 2 log(3/2)(1+ε)/(ML) · (1/(2κ_δ)) ,
   c₂ = 2 log(3/2)/κ_δ = 1.6227651872   (ASSEMBLY's 4 log(3/2) = 1.6218604 up to 2κ_δ) .
```

### 5.4 FINDING (a refusal of the brief, L-17)

The brief asks me to conclude `t_* = (2/(ML))∫₁^{3/2}dλ/(λP_h(λ))·(1+ε)`.  That formula is not
what the system yields, and it is **optimistic**.  Computed (`u4`):
```
   2∫₁^{3/2}dλ/(λ P_h(λ))  =  0.6698262   (the brief's θ_*)
   4(1−√(2/3))             =  0.7340137   (prove-lagrangian's closed form, bang-bang)
   log(3/2)/κ_δ            =  0.8113826   (the conservative value, and what is proved)
   system (A) at L = 640   =  0.7399101   (what the tapered integro-ODE actually gives)
```
The brief's expression is `8.7447189 %` below the closed form and `9.4719469 %` below the
true discretised answer, because it replaces the integro-ODE's shell integral `H(0,θ) = ∫₀¹λ dσ′ = √λ` by
`P_h(λ(0))` — an over-estimate of the rate, hence an under-estimate of `t_*`.  It is not used.
The theorem below uses `c_*`; the true acceleration of the tapered integro-ODE over the
conservative bound is `10.264382 %` at `L = 640` (`u4`), and that is what a sharper clock could
claim — not the brief's `17.4463134 %`.

---

## 6. Part (5): the budget

`u5_budget.py` is `s3close/assembly/a2_budget.py` (sha256 `96c1ee0b10be43c3…`) with
`bootstrap()`, `COLUMNS` and `assemble()` replaced.  `viscous_budget()`, `gauss_tail_aniso()`,
`J_pow()`, `C_a_proved()`, `eps_Tprime()` and the `RECORD` block are carried over unchanged, so
the viscous half of the budget is literally the assembly's instrument and the comparison is
apples to apples.  Datum constants are re-read from `u2_results.json` (computed here) instead
of `a1_results.json`.

### 6.1 The columns

| column | `C″` | `C_a` | `C_R` | `K₂/(M/ρ₀)` | `P(|Z_τ|>d/2)` | Lemma T′ sensitivity |
|---|---|---|---|---|---|---|
| `proved` | `151.15` | `8.697888` | `29.024165` | `161.7735` | Theorem V.1 `B1` | as proved |
| `proved_Ca_measured` | `151.15` | `0.069` | `3.137501` | `161.7735` | Theorem V.1 `B1` | as proved |
| `measured` | `11.74` | `0.069` | `3.137501` | `5.3854` | exact affine Gaussian | `× 1/6.0668328` |

All three use `C₁ = 1`, `δ = 7.5°`, `δ_m = 5°`, `φ₀ = 30°`, `G_collar = 0`,
`λ_max = e^{3c/4} = 1.8377407`, `r_h = 0.9199394` on `[1,λ_max]`, `κ_δ = 0.4997212305210886`,
`c = c_* = 0.8113825936`.

### 6.2 The term table at `c = c_*`, `f = 1`

| column | `L` | `ε_ℓ` | `ε_{C_a}` | `μ` | `μ_J` | `ε_{T′}` | `c_G` | `ε_v` | `ε` |
|---|---|---|---|---|---|---|---|---|---|
| proved | 40 | `0.270308` | `0.652706` | `∞` | `∞` | `∞` | `4.31351` | `1.5961e4` | `∞` |
| proved | 160 | `0.0675769` | `0.163176` | `∞` | `∞` | `∞` | `1.99118` | `18.3529` | `∞` |
| proved | 640 | `0.0083177` | `0.0407941` | `0.784016` | `0.0368686` | `1.77198e4` | `1.41060` | `0.0297174` | `∞` |
| proved | 2560 | `0.00139637` | `0.0101985` | `0.179440` | `0.00892439` | `5.00872` | `1.26546` | `0.00131599` | `∞` |
| proved_Ca_measured | 160 | `0.029830` | `0.00129447` | `0.652964` | `0.0125752` | `1345.77` | `1.99118` | `18.3529` | `∞` |
| proved_Ca_measured | 640 | `0.00537291` | `0.000323618` | `0.112894` | `0.00266631` | `2.11770` | `1.41060` | `0.0297174` | `∞` |
| proved_Ca_measured | 2560 | `0.00127484` | `8.09045e-5` | `0.0258384` | `0.000649358` | `0.303593` | `1.26546` | `0.00131599` | `0.444220` |
| measured | 160 | `0.0256421` | `0.00129447` | `0.418122` | `0.0118359` | `10.7153` | `1.28422` | `0.0135302` | `∞` |
| measured | 640 | `0.00533659` | `0.000323618` | `0.101408` | `0.00266000` | `0.294187` | `1.23386` | `0.000774592` | `0.431788` |
| measured | 2560 | `0.00127431` | `8.09045e-5` | `0.0251613` | `0.000649269` | `0.0485689` | `1.22127` | `5.00003e-5` | `0.0532646` |

The binding term is `ε_{T′}`, through `μ = 446.21/L` (proved) or `64.251/L` (measured `C_a`):
`ε_{T′} ≈ (15π/4)μ/r_h = 12.8062478 μ` (proved sensitivity), so the window closes when
`μ ≲ 0.04`.  Everything else is `O(1/L)` with an `O(1)` constant.  Compare ASSEMBLY's own term
table, where `ε_{T′} = ∞` in *every* row at `c = c_*` and every `L ≤ 2560` in every column.

### 6.3 `ε(L)` at `c = c_*`, minimised over `f ∈ {0.05,…,32}`

| column | `ε(40)` | `ε(160)` | `ε(640)` | `ε(2560)` | best `f` at the first finite entry |
|---|---|---|---|---|---|
| proved | `∞` | `∞` | `∞` | `∞` | — |
| proved_Ca_measured | `∞` | `∞` | `∞` | `0.4439660` | `0.5` |
| measured | `∞` | `∞` | `0.4307170` | `0.0530984` | `0.5` / `0.25` |

### 6.4 `L_*` and `Λ_*`

`log Λ_* = 2L_* + 3.3643455` (ASSEMBLY (1.2); the shift is reproduced here to `1.1·10^{-7}`
in `C_E` (rel `1.2506208·10^{-7}`) from `s = 1/sin δ = 7.66129757554`).

| column | `L_*` (`ε≤½`) | `log Λ_*` | `L_*` (`ε≤0.1`) | `log Λ_*` | ASSEMBLY's `L_*` (`ε≤½`) |
|---|---|---|---|---|---|
| proved | **`16161.986`** | **`32327.336`** | `53783.018` | `107569.40` | `4.5253·10^{14}` |
| proved_Ca_measured | **`2389.4697`** | **`4782.3038`** | `7774.7586` | `15552.882` | — |
| measured | **`600.25760`** | **`1203.8795`** | `1541.8189` | `3087.0021` | `4.8315·10^{3}` |

The improvement is `2.7999653·10^{10}` in `L_*` in the proved column and `8.0490443` in the
measured one,
at `N = 1` sub-windows — i.e. **without** the restart lemma ASSEMBLY §3.4 prices at `N = 16`.
For reference, ASSEMBLY's `N = 16` restart gives `L_* = 366.96` (proved, `C_K = 66.6622`); the
slaved reference alone gives `16162` at `N = 1`, and the two are independent, so a restart on
top of the slaved reference would multiply the gains.

### 6.5 Sensitivity to the two things that are not proved

**The collar gradient `G_collar` (hypothesis (R-z)).**  It enters `C_R` additively and `C_R`
enters the drive linearly:

| `G_collar` | `L_*` proved (`ε≤½`) | `L_*` measured (`ε≤½`) |
|---|---|---|
| `0` | `16161.986` | `600.25760` |
| `C_a(30°) = 8.697888` | `20717.086` | `1765.1714` |
| `5C_a(30°) = 43.489440` | `38937.953` | `6428.8103` |

A factor `5` in the unproved constant costs a factor `2.4092307` in `L_*` (proved) — **not an
exponential**.  That is the structural gain: `C_R` is in the drive, not in the exponent.

**The a priori cap `λ_max`.**  §1.1 (P3) proves `λ ≤ e^{3c/4} = 1.8377407` with no extra
hypothesis.  A self-consistent argument would give `λ_max = (3/2)(1+O(1/L))` (material on a
shell where `h_δ = 1` has `|ω^θ| = Mλ_material ≤ (3/2)M`, and `λ_material ≥ λ e^{-C_a c/L}`);
that is not written here.  At `λ_max = 3/2`:

| | `L_*` proved (`ε≤½`) | `L_*` measured (`ε≤½`) | `ε(640)` measured |
|---|---|---|---|
| `λ_max = e^{3c/4} = 1.8377407` | `16161.986` | `600.25760` | `0.4307170` |
| `λ_max = 3/2` | `8338.8577` | `293.75883` | `0.1239899` |

So proving the self-consistent cap is worth a factor `1.9381535` (proved) to `2.0433687`
(measured) in `L_*`.

### 6.6 `(H-K2)` costs almost nothing, again

| `C_K` | `L_*` proved | `L_*` proved_Ca_measured | `L_*` measured |
|---|---|---|---|
| `1` (ASSEMBLY's `a2` baseline) | `16159.492` | `2374.4690` | `598.82938` |
| `3.0202` (hk2 measured, global probe) | `16159.528` | `2374.6584` | `599.53487` |
| `5.3854` (hk2 measured, over the window) | `16159.569` | `2374.8801` | `600.25760` |
| `66.6622` (hk2 proved, `λ=1`) | `16160.661` | `2380.6144` | `614.46961` |
| `161.7735` (hk2 proved, over the window) | `16161.986` | `2389.4697` | `632.15997` |

The whole of `(H-K2)` moves `L_*` by `+0.015 %` (proved) and `+5.6 %` (measured).  ASSEMBLY's
§3.6 conclusion — closing `(H-K2)` was necessary and is not sufficient — is unchanged, and the
reason is now different: it is no longer swamped by `e^{34}`, it is simply small.

### 6.7 Two variants that must be reported (`u6`)

**(V-strain).**  The Theorem V.4 viscous loss charged to the strain lower bound as well as to
the vorticity (§4.2's carried gap).  It costs almost nothing:

| column | `L_*` (`ε≤½`), headline | `L_*` with (V-strain) | `ε(2560)` headline | `ε(2560)` with (V-strain) |
|---|---|---|---|---|
| proved | `16161.986` | `16163.288` | `∞` | `∞` |
| proved_Ca_measured | `2389.4697` | `2397.0839` | `0.4439660` | `0.4462687` |
| measured | `600.25760` | `601.18392` | `0.0530984` | `0.0531698` |

So the carried gap, if it were closed at the obvious cost, would move `L_*` by `+0.0080543 %`
(proved) to `+0.3186568 %` (proved, `C_a` measured).  It is not a load-bearing gap; it is bookkeeping.

**(conj).**  The budget with the conjugated bootstrap (2.6) in place of (2.3).  Its feedback
exponent is `0.34679` (proved) and `0.04125` (measured) at `L = 640` — an order of magnitude
below (2.3)'s — but its geometric factor `λ_max³ = 6.2066` in the drive costs more than the
exponent saves at these constants:

| column | `L_*` (`ε≤½`), Grönwall (2.3) | `L_*`, conjugated (2.6) | `pc` (`L=640`) conj |
|---|---|---|---|
| proved | `16161.986` | `48986.447` | `0.34678624` |
| proved_Ca_measured | `2389.4697` | `5587.9475` | `0.04124785` |
| measured | `600.25760` | `1442.1759` | `0.04124785` |

The conjugated form is nonetheless the structurally better estimate and is the one to sharpen:
it uses **no** gradient bound on `u` at all, so it does not need **(Γ-off)** in the map error,
and its exponent is `O(c/L)` rather than `O(c)`.  What it needs is a component-wise treatment
of `T_λ^{-1}R` in place of the crude `|T_λ^{-1}R| ≤ λ²|R|`.

### 6.8 FINDING: the angular supremum (`u7`)

`μ` is a supremum over **all** `x ∈ S` — Lemma T′'s (H2) and (H3) are stated that way — but the
far/near collar constant `C_a(φ) = … + 3.999218/sin φ + π/8` **diverges on the axis**.
ASSEMBLY §2.4 evaluates it "at the worst angle over the window, `φ = φ₀ = 30°`", which is the
worst angle along the *trajectory*, not over the *shell*.  Neither note therefore establishes
the majorant on the whole shell.

| `φ` | `C_a(φ)` | `C_R(φ)` (`G_collar=0`) |
|---|---|---|
| `90°` | `4.6986701` | `16.333364` |
| `30°` (both notes' choice) | `8.6978881` | `29.024165` |
| `10°` | `23.730031` | `75.178170` |
| `7.5°` (`= δ`, the taper edge) | `31.338651` | `98.289489` |

Restricting honestly to the cone `{δ ≤ φ ≤ π−δ}` and taking `C_a(δ)`:

| column | `L_*` (`ε≤½`) at `φ₀` | `L_*` on the cone at `δ` |
|---|---|---|
| proved | `16161.986` | `52437.352` |
| measured | `600.25760` | `9884.1636` |

A factor `3.2444869` (proved) and `16.466536` (measured), **not an exponential** — again because `C_a`
sits in the drive.  Inside the taper cone `{φ < δ}` the bound is `+∞` as it stands, and
closing that needs `far-near-kernel-lemma` §6's observation that admissibility buys the
divergence back (its `c_edge` saturates for a `δ`-tapered datum) turned into a proved constant.
That is a real hole, it is **inherited**, and it is where I would attack this note first.

---

## 7. The theorem, its hypotheses, and the status of every step

> **THEOREM (S3, slaved-reference form).**  Fix `δ = 7.5°`, `δ_m = 5°`, `φ₀ = 30°`,
> `s = 1/sin δ`, `f > 0`, and the datum of `ASSEMBLY` §1.1 **with the hk2 correction**: the
> radial cutoff `tanh`-mollified at `w₀ = 0.25ρ₀` (a sharp radial edge makes `K₂ = +∞`).  Put
> `ν = Mρ₀²/s²`, `R = ρ₀e^L`, `τ = c_*/(ML)`, `c_* = log(3/2)/κ_δ = 0.8113825936`.
> Assume, on `[0,τ]`:
>
> * **(Γ-off)** `‖∇u(σ)‖_{L^∞(slab)} ≤ 2a(0,σ)(1 + C₁/L) + C″M`,  `C₁ = 1`, `C″ = 151.15`.
>   *Not proved here; a parallel seat is proving it.*
> * **(R-z)** `|r∂_r[a − 𝔞(|·|)]| ≤ G_collar·M_σ` on the collar `|X|/2 < ρ′ < 2|X|`,
>   `G_collar` absolute.  *Not proved anywhere; it is the gradient half of ASSEMBLY BLOCK 3
>   restricted to the collar, and ASSEMBLY needs it too without naming it (§3.4).*
> * **(A-cone)** the majorant is established on `{δ ≤ φ ≤ π−δ}` only (§6.8).
> * **(H3-V)** ASSEMBLY BLOCK 5's hypothesis for Theorem V.4, which §1.1's datum satisfies away
>   from its kinks.
> * *carried gap* **(V-strain)**: the viscous defect of the strain functional (§4.2), worth
>   `≤ 0.3186568 %` in `L_*` (§6.7).
>
> Then, for every `L ≥ L_*`,
> ```
>    T_d(u₀) ≤ t_* = c_*(1+ε(L))/(ML) = 2 log(3/2)(1+ε)/(2κ_δ ML) ,
>    𝒯(Λ) ≤ c₂(1+ε)/log Λ ,  c₂ = 2 log(3/2)/κ_δ = 1.6227652 ,  log Λ ≥ 2L_* + 3.3643455 ,
> ```
> with `L_* = 16161.986` on all-proved constants (`G_collar = 0`, `φ = φ₀`), `2389.4697` with
> the measured `C_a`, and `600.25760` on all-measured constants.

### 7.1 Status per step

| # | step | status |
|---|---|---|
| T1 | the datum, the trajectory margins at the reference's own cap (`φ(λ_max) = 74.407°`, equator margin `15.593° > δ_m = 5°`) | **PROVED** (`u2`) |
| T2 | `κ_δ = 0.4997212305210886`, `r_h[1,λ_max] = 0.9199394406`, `P_h` increasing on `[1,2.075]`, `P₁(λ) = λ` | **PROVED** (`u2`, `u4`; reproduces the record to `0.0`, `3.0e-7`, `1.6e-8`) |
| T3 | the definition of `𝔞`; (P1) `𝔞 ≥ 0`; (P2) `\|∂𝔞/∂logρ\| ≤ M_s/2`; (P3) `λ ≤ e^{3c/4}`, decreasing | **PROVED** (`u1`, `u2`; (P1) is the second load-bearing use of Lemma 2) |
| T4 | `Λ` is the flow of the `l=1` field `𝔞·DX`; `J_Λ = λ²[1+(ρλ′/λ)(1−3cos²φ)]` | **PROVED** (sympy, residual `0`) |
| T5 | the driver (2.2); **the `l = 1` rate error is identically zero** | **PROVED** (identity, not estimate) |
| T6 | `\|R_y\| ≤ C_a M_s\|X\|` for the true field, no homogeneity used | **PROVED** (far/near lemma + Lemma 2) |
| T7 | `\|R_z\| ≤ [2C_a + log(1/sinφ) + G₁]M_s\|X\|` | **PROVED-modulo-(R-z)** |
| T8 | the Lipschitz term `Γ\|Φ−Λ\|` | **PROVED-modulo-(Γ-off)** |
| T9 | the majorant (2.3) and its solution `μ(c) = O(1/L)`, `pc = Gc` | **PROVED** given T6–T8 (linear ODE; closed form re-derived independently in the gate and agreeing with the solver to `1e-13`) |
| T10 | `μ_J` from `d/dθ log(J_Φ/λ²)` directly; Lemma T′ (H3) obtained without Corollary 3 | **PROVED** |
| T11 | `𝔞` and Lemma T′'s `a_ref` are the same functional, up to the inner-shell integral `≤ (M_s/2)ℓ_loss` | **PROVED** (Consequence A at `x=0`; an identity) |
| T12 | Lemma T′ applied once, to `η₀1_{S^out}`, with `r_h` on `[1,λ_max]` | **PROVED** (R1) |
| T13 | the restriction step (positivity + `z`-oddness of the restricted piece) | **PROVED** (Lemma 2 + reflection-invariance of `S^out`) |
| T14 | the integro-differential inequality (5.1) and the quasimonotone comparison | **PROVED** (two-line proof, §5.1; numerics `u4` confirm `A ≥ B`, `A ≥ C` at `L=40,160,640`) |
| T15 | the conservative clock, `c₂ = 2log(3/2)/κ_δ` | **PROVED** given T14 |
| T16 | `ε_v` (Theorem V.4 + V.1, refuter's sizing, `K₂` from hk2) | **PROVED-modulo-(H3-V)**, plus the carried gap (V-strain) |
| T17 | the angular supremum | **NOT ESTABLISHED on the taper cone — (A-cone)**, §6.8 |
| T18 | the budget, `L_*`, `log Λ_*` | **COMPUTED** (conditional on T7, T8, T17) |
| T19 | BFG confrontation | **NOT TRIGGERED** — (S3) is still not proved, so ASSEMBLY §4.3 stands verbatim |
| — | ASSEMBLY (2.1) charges the vector remainder to the scalar constant `C_a` | **FINDING** (§3.4) |
| — | the brief's `t_*` formula is `9.5 %` optimistic | **FINDING** (§5.4) |
| — | the angular sup is taken along the trajectory, not over the shell | **FINDING** (§6.8) |

### 7.2 Where a referee should attack

1. **(the brief's own worry, and mine) Is the exterior-strain functional that defines `𝔞` the
   same object Lemma T′'s `a_ref` controls?**  §4.1 answers yes and does so by an *identity*:
   at `x = 0` every shell is exterior, and Consequence A says the origin value **is** the sum
   of the shells' `l=1` interior modes, so `𝔞(ρ,s) = a(0,s) − ∫_{ρ₀}^{2ρ}F dlog ρ′` with the
   subtracted piece in `[0,(M_s/2)log(2ρ/ρ₀)]`.  **The place to attack is not that identity but
   step 2 of §4.2**: I lower-bound `𝔞(|X(s)|,s)` by the contribution of the *image of a label
   sub-shell*, which needs (i) `F ≥ 0` for each piece separately — true because `Φ_s` commutes
   with `z ↦ −z` (the datum and NS preserve the symmetry, Lemma 2) and `S^out` is
   reflection-invariant — and (ii) the label cut `ℓ_loss` to be large enough that the image is
   outside `2|X(s)|`, which uses `|Φ| ≥ (1−μ)λ_max^{-2}|x|` and therefore feeds `μ` back into
   `ℓ_loss`.  That feedback is *logarithmic* (`log((1+μ)/(1−μ))`), not exponential, so it does
   not restore the loop; but it is a feedback and a referee should check it.
2. **Is the near-field velocity at the point really `≤ C_a M` for the true field at all times
   with no homogeneity?**  For the `r`-component (the strain `a`): **yes**, and that is the
   whole point of `far-near-kernel-lemma` — its hypotheses are `|ω^θ| ≤ M_s` and `z`-oddness,
   nothing else.  For the `z`-component: **no**, §3.2, and this is the sharpest hole.  It is
   inherited from ASSEMBLY, which charges the whole vector remainder to `C_a` (§3.4).
3. **The angular supremum, §6.8.**  `C_a(φ) → ∞` on the axis and both notes evaluate at `φ₀`.
   The cone version costs `×3.2444869`; the taper cone itself is uncovered.
4. **The a priori cap `λ_max = e^{3c/4} = 1.8377407` instead of `3/2`.**  It is what the
   contradiction hypothesis alone gives, and it costs `×1.9381535`–`×2.0433687` in `L_*` and takes `r_h`
   from `0.9818223` to `0.9199394`.  A self-consistent argument should recover `3/2(1+O(1/L))`.
5. **`ε_{T′}` is still evaluated at the END of the window** (as in ASSEMBLY), i.e. the largest
   `μ`, and it multiplies the *whole* strain integral.  Evaluating it under the time integral
   would gain a further factor `≈ (e^{Gc}−1)/(Gc)`.  I did not do it.
6. **The conjugated form (2.6) is the better estimate and is not the one used.**  Its
   `λ_max³ = 6.2066` comes from the crude `|T_λ^{-1}R| ≤ λ²|R|` and `|Φ| ≤ λ(1+ê)|x|`; a
   component-wise treatment would remove most of it and make the exponent `O(c/L)` decisive.
7. **The window is still one window.**  Nothing here is a restart lemma; BLOCK 2 remains open.
   It is simply no longer the binding constraint.

---

## 8. Gate and files

`check_constants.py` rebuilds every load-bearing quantity from mathematics and asserts it
against the stored JSONs, then audits this document token by token.

| file | what it establishes |
|---|---|
| `u1_algebra.py` / `u1_results.json` | the `l=1` mode is exactly `𝔞·DX` (3-D divergence free, 5-D divergence `2𝔞`); `∂_zu^z = −2a − r∂_ra`; `J_Λ`; the conjugation identity (2.5); `T_λ` commutes with `D` and `λ^{-2}\|v\| ≤ \|T_λv\| ≤ λ\|v\|`; `∂𝔞/∂logρ = −F(2ρ)`; the origin identity |
| `u2_constants.py` / `u2_results.json` | `P_h`, `κ_δ`, `r_h` on `[1,3/2]` and `[1,λ_max]`, `Φ_h(3/2)`, `P₁(λ)=λ`, monotonicity of `P_h`; `λ_max`, `κ_s`; `C_a(φ)`, `C_R`; `ℓ_loss`, `ℓ_shell`; trajectory margins; the record cross-checks |
| `u3_bootstrap.py` / `u3_results.json` | the three feedback exponents; `μ(c)`, `μ_J(c)`, `ε_{T′}` on an `L`-grid, in both the Grönwall and the conjugated form |
| `u4_shells.py` / `u4_results.json` | the 200-shell integro-differential system, the comparison against (B) and (C), the acceleration ratio, and the arithmetic refuting the brief's `t_*` |
| `u5_budget.py` / `u5_results.json` | the budget: a copy of `a2_budget.py` with `bootstrap`, `COLUMNS`, `assemble` replaced; the term table, `ε(L)`, `L_*`, `log Λ_*`, the `G_collar` and `λ_max` sensitivities, the `C_K` sweep |
| `u6_variants.py` / `u6_results.json` | (V-strain) and the conjugated budget |
| `u7_cone.py` / `u7_results.json` | the angular supremum: `C_a(φ)`, `C_R(φ)`, and `L_*` on the cone `{δ ≤ φ ≤ π−δ}` |
| `check_constants.py` / `gate_stats.json` | the gate |
| `SHA256SUMS` | computed with `shasum -a 256`, never typed |

```
95 CHECKS, 95 PASS, 0 FAIL
DOC AUDIT: 468 numeric tokens in PROOF.md (the gate's own report box excluded),
           415 traced to this folder's JSONs, 53 quoted from a named source, 0 untraced
```

The pass count is `len(CHECKS)`, printed, never a literal.  No check compares a literal against
itself: the five that did on the first run (the `C_R/C_a` factor, the `ε_{T′}` slope, the ratio
of the two exponentials, and the two improvement factors) were caught by the gate as FAILs
because the values I had typed into the document were rounded, and they were replaced by the
computed ones.  The `FOREIGN` dictionary in `check_constants.py` lists every quoted constant
with the seat each comes from.

What the gate **cannot** catch, stated so it is not mistaken for coverage: it checks that the
numbers displayed are the numbers the scripts produced and that the algebra between them is
consistent.  It does **not** check that (2.3)–(2.4) is the right majorant, that `C_R` is the
right constant, that (R-z) or (Γ-off) are true, or that §4.2's restriction step is airtight.
That is what §7.2 is for.

**`FL-000` stands.  Nothing here touches the headline problem.**

# refuter-lagrangian-b — adversarial review of `lower/prove-lagrangian`

Sub-seat `lower/refuter-lagrangian-b`, DTC-2026-09-06. Every number below came from a script in
this folder that I wrote and ran (`x1`–`x5`), or from re-running the attempt's own scripts from a
copy (`copy/`). Nothing outside this folder was written. Numerics falsify, never prove.

## 0. Verdict

**REFUTED — MAJOR, on the `key` claim, not on the theorem's arithmetic.**

The attempt's mathematics is right where it says it is exact. I confirmed LEMMA 1 by two routes it
did not use, and I ran the one test it skipped (the strain at the **material** point of the
deformed configuration, not the empty origin) — that test passes too, which *strengthens* the
attempt. What fails is the `key`:

> "The step (2) of the brief (perturbative persistence of the strain) is not needed … the
> perturbative O(c) step the brief asked for is replaced by an identity."

It is not replaced. It is relocated into GAP T and made **strictly stronger**, and the attempt's
own §5 status table says so. Separately, the identity is not the cancellation it is advertised to
be, the accelerated constant `c₂ = 1.4680` has no empirical support and the only measurement of it
sits on the value it claims to improve, and the gate contains checks that cannot fail.

## 1. What I re-ran and what reproduced

`copy/` holds the attempt's `s1`–`s7` + `check_constants.py`, re-run by me. All exit 0.
`check_constants.py` run in place in the attempt's own folder: **rc = 0, 88 `PASS` lines, 0 `FAIL`**
(`copy/gate_out.txt`). Its stored `s3['mol_vs_closed']` max is `2.537e-9` against a tolerance of
`3e-9` — a 16 % margin on the tightest check in the gate.

## 2. FIND 1 (MAJOR) — the perturbative step is relocated, not removed

LEMMA 1 is an identity about the *hypothetical* field `η₀∘T_λ^{-1}`. Using it needs, in order:

1. **L3v** — the deformed field's velocity is still a per-shell pure strain + `O(Mρ)`
   (attempt's status: **SKETCH**, "computed, not proved");
2. **GAP T** — the true flow map equals `T_{λ(σ)}` up to relative `O(c/L)`, *and* the functional
   `a[η](0)` is Lipschitz-stable under an `O(c/L)` relative perturbation of the flow map in the
   scale-invariant metric, with constant `C·M·L` (attempt's status: **GENUINE GAP — SKETCH**;
   "Everything downstream is conditional on it").

The brief's step (2) asked for a **one-sided lower** bound, `a ≥ κML(1−O(c))` on `[0,τ]`. GAP T is
a **two-sided Lipschitz estimate for a Calderón–Zygmund-type kernel written in log-radial
coordinates** — strictly stronger than what it replaces, and unproved. So the `key` sentence
inverts the direction of the difficulty. The theorem statement in the attempt is honest
(`PROVED-MODULO`); the `key` is not.

## 3. FIND 2 (MAJOR) — LEMMA 1 contains no cancellation (deflation)

The attempt's §6.3 explains LEMMA 1 as an exact balance: "rotation toward the equator, where the
weight `|cos φ| sin²φ` vanishes, versus inward collapse near the axis … exactly neutral. This is
the content of `∫₀¹ v²(Av²+B)^{-5/2}dv = λ/3`." That is a misdescription of its own identity.

In Eulerian variables (my `x1_deflate_lemma1.py`):

* `η₀ = −M sgn(z) h(φ)/r` is a function of `r` **alone** (times `sgn z` and the taper). `T_λ`
  scales `r` by the single constant `λ`, so as a function of the Eulerian point
  `η_λ(w) = η₀(T_λ^{-1}w) = −Mλ sgn(z_w) h(φ'(w))/r_w` — literally `λ` times the same `1/r`
  profile.
* `T_λ(shell) = {ρ₀ < |w| G(t) < R}`, `G(t) = √((1−t²)/λ² + λ⁴t²)`. **Both** shell radii scale by
  `1/G(t)`, so the log-thickness of the support is `log(R/ρ₀) = L` in **every** direction `t`.
* `a(0)` is linear in `η` with a scale-free kernel ⟹ `a_λ(0) = λ a₀(0)` for `h ≡ 1`. Two lines,
  no geometry.

Numerical demonstration that the geometry never enters: deleting `G(t)` from the Eulerian
integrand leaves the answer unchanged to 9 digits.

| λ | Eulerian `a_λ(0)/(ML)` | same with `G(t)` deleted | `λ/2` |
|---|---|---|---|
| 1.0 | 0.49999999955921703 | 0.49999999955921703 | 0.5 |
| 1.25 | 0.6249999994490216 | 0.6249999994490216 | 0.625 |
| 1.5 | 0.7499999993388258 | 0.7499999993388258 | 0.75 |
| 5.0 | 2.4999999977960865 | 2.4999999977960865 | 2.5 |

(Also: sympy antiderivative residual `0`, `P₁(λ)−λ ≡ 0`, mpmath 30-dps rel `0.0` — the attempt's
Lagrangian route is *correct*, just an unnecessary change of variables.)

Consequences. (a) The identity is an algebraic property of **this** datum (`η` a function of `r`
alone) under **exactly** the uniform coaxial strain; it generalises to nothing, which is precisely
why GAP T is load-bearing. (b) "The strain scales precisely like the stretched vorticity" is
circular inside the ansatz: `λ` *is* the assumed stretching factor. (c) §7.1's novelty claim is
over-sold.

## 4. FIND 3 (MAJOR, FL-043) — the gate cannot fail on part of what it reports, and is blind to
77 % of the numbers it reads

`x4_gate_mutation.py`. Baseline: the gate exits 0, prints 87 check lines (80 `ck()` results + 7
plain-`print` lines from the `Phi>=1` asserts), and its last line is the hardcoded
`print(("ALL %d CHECKS PASS" % 0) ...)` → **`ALL 0 CHECKS PASS`**. The NOTE's "88 assertions …
exits 0 (`ALL CHECKS PASS`)" is one off on the count, and the count the gate itself reports is a
literal `0`.

**Checks that cannot fail for any outcome (6, each verified by reading the line):**

* *masked to zero* —
  `ck('sum |H| ||C||/l^2 (lam=1)', s5['sum_absH_over_l2_tail_from_3']*0+0.3958, 0.3958, 1e-9)`.
  The stored value is multiplied by `0`. This is the **only** gate line for the NOTE's L3
  remainder constant, and the `s5` key it names holds **0.011311110612413619** — a different
  quantity from the `0.3958` the NOTE attributes to it (`0.3958` actually comes from `s6`'s
  `worstcase_Lemma3_bound`). A 35× mismatch between the cited key and the displayed value passes.
* *no stored input at all (5)* — each recomputes the constant it is compared against:
  `ck('model Dt a / (M^2 L^2) = 1/8', 1/8, 0.125, 0)`;
  `ck('beta/(M^2L^2) = 3/8', 1/8+1/4, 0.375, 0)`;
  `ck('log-rho excursion 2 log(3/2)', 2*np.log(1.5), 0.8109302162163288, 1e-14)`;
  `ck('phi0=30deg -> after lam=1.5', np.rad2deg(np.arctan(np.tan(np.deg2rad(30))*1.5**3)), 62.833, 1e-4)`;
  `ck('measured T32 M a* / log(3/2)', 0.479/np.log(1.5), 1.1814, 1e-4)` — the last from a hardcoded
  `0.479`, not from the `viscous-numerics` data it purports to check.

**Dynamic mutation.** Every stored number in `s1`–`s7_results.json` perturbed one at a time
(`v → 1.5v + 0.37`), gate re-run for each:

| stored numbers | gate flips to non-zero exit | gate does not notice |
|---|---|---|
| **416** | **95** | **321 (77.2 %)** |

So the gate returns the same verdict for 321 of the 416 outcomes it reads. (Some blindness is
benign — stored diagnostics the NOTE never quotes — but it includes displayed numbers, e.g.
`s1/Pprime1/flat/0 = 1.0000000000192764` (the NOTE's `Φ'(1) = +1`, the sign of the acceleration),
`s2/taper_table/30.0deg/kappa = 0.48362523497974136`, and every entry of
`s1/route3_eulerian_taper7.5`.) `x4_out.txt` lists the first 40; the full list is in
`x4_gate_mutation_results.json`. *(My static regex in T1 over-fires — it flags any `ck` body
without a literal `s1[`…`s7[` token, so lines reading through the local aliases `t`, `v`, `r`,
`mx`, `val` are wrongly listed. The six above are hand-verified; the dynamic count is the sound
measurement.)*

## 5. FIND 4 (MINOR) — `0.3958` is a truncated partial sum

The NOTE displays `Σ_{l≥3}|H_l|·‖C_l‖_∞/(l(l+3)−4) = 0.3958` as the L3 remainder constant. `s6`
truncates at `LMAX = 601`. My independent computation (`x1`, Gegenbauer streamed to `l = 2001`,
60001-point θ-trapezoid; `H₁ = −0.83333333276`, `H₃ = 0.074999999`, `H₅ = −0.147023808` against the
exact `−5/6, 3/40, −247/1680`):

| truncation | 101 | 401 | **601** | 801 | 1201 | 1601 | 2001 |
|---|---|---|---|---|---|---|---|
| partial sum | 0.34990 | 0.38858 | **0.39582** | 0.40015 | 0.40529 | 0.40837 | 0.41046 |

Terms behave like `0.8254 · l^{-3/2}` (`term·l^{3/2}` = 0.8776, 0.8430, 0.8307, 0.8212 at
`l` = 101, 401, 801, 1601), so the series **does** converge — L3's boundedness claim survives — but
the true sum is ≈ **0.4295**, and the displayed `0.3958` understates it by ~8.5 % with no statement
of its truncation dependence. This is exactly the number the gate is blind to (FIND 3).

## 6. FIND 5 (MAJOR) — the accelerated constant has no empirical support, and the only data sit on
the value it claims to beat

`x5_empirical_status.py`, using `viscous-numerics`' measured `Q' = T₃₂·M·log(R/s*) = 0.988 ± 0.018`
(N = 2…7) and its `t = 0` strain law `κ = 0.437 / 0.415 / 0.411` at N = 3 / 5 / 6:

| κ | frozen model `log(3/2)/κ` | z | accelerated `2(1−√(2/3))/κ` | z |
|---|---|---|---|---|
| 0.437 (N=3) | 0.92784 | −3.34 | 0.83983 | −8.23 |
| 0.415 (N=5) | 0.97702 | −0.61 | 0.88435 | −5.76 |
| **0.411 (N=6)** | **0.98653** | **−0.08** | **0.89296** | **−5.28** |

At the finest resolution the **frozen** clock lands on the measurement (z = −0.08); the
**accelerated** clock is 5.3 sd below it. The advertised improvement is
`1 − 1.4680274/1.6218604 = 9.485 %` of `c₂`; the accelerated model's shortfall against the data is
`1 − 0.892961/0.988 = 9.619 %`. **The entire gain the `key` sells is the size of the model's
disagreement with the only measurement of it.** The attempt states this in §6.6 and then still
carries `1.4680274` in the headline table and in clause (b) of the theorem.

## 7. FIND 6 (MINOR) — LEMMA 2's novelty claim is wrong at source, and it is inert

§7.5: "the first non-perturbative statement in this line about `t > 0`". Sign preservation of
`η = ω^θ/r` under axisymmetric no-swirl Navier–Stokes is the classical maximum principle for the
transported-diffused `η` — it is what makes `‖ω^θ/r‖_∞` non-increasing and it underlies the
standard axisymmetric global-regularity theory (Ukhovskii–Yudovich line). Only the packaging (odd
in `z` ⟹ signed integrand ⟹ `a(0,t) ≥ 0`) is new. And it is **never used**: no step of §4 invokes
LEMMA 2, and `a(0,t) ≥ 0` is a statement at the centre of the empty hole, where `ω ≡ 0` by
construction, so it cannot bound the growth of `‖ω‖_∞`.

## 8. What survives — including one step the attempt did not test, which I ran and it passes

I built an independent instrument (`gegen5d.py`): zonal Gegenbauer `C_l^{3/2}` on `S⁴` with the
radial Green's function
`ψ_l(ρ) = ρ/(2l+3)[∫_0^∞ e^{(1-l)w}q_l(ρe^w)dw + ∫_{-∞}^0 e^{(l+4)w}q_l dw]`,
`a = −Σ_l[ψ_l' t C_l + (ψ_l/ρ)(1−t²)C_l']`, applied to the strained field
`η_λ(ρ,t) = Ξ(t)/ρ` on `{ρ₀ < ρG(t) < R}`. Controls: Hill's ball (`η = A` on `ρ<b`) gives
`a = Az/5` analytically from the same formula; `a(0,0) − (M/2)L` = 4.5e-4 … −9.6e-4 over
`R = 64…4096`; material offset at `(ρ₀, 10°)` = **+0.21423** (R = 4096) against the campaign's
`+0.216773 ± 2e-6` — 2.5e-3, my grid.

*(Recorded deviation: my first pass mis-shifted the Gegenbauer recurrence and gave −0.113 for that
offset. Fixed to `n C_n = 2t(n+α−1)C_{n−1} − (n+2α−2)C_{n−2}`; verified against
`C₃ = 17.5t³−7.5t`, `C₅ = 86.625t⁵−78.75t³+13.125t`, `C_l(1) = (l+1)(l+2)/2`. The attempt's
`geg_all` was correct all along.)*

**A fifth route to LEMMA 1, and the material-point test the attempt skipped** (`x2`; `M = 1`,
tracked point starts at `ρ₀`, `φ₀ = 30°`, moved by `T_λ`):

| R | λ | `a(0) − ½λL` | `a(material) − ½λL` | ratio |
|---|---|---|---|---|
| 65536 | 1.0 | −6.3e-4 | −0.01808 | 0.99674 |
| 65536 | 1.1 | −6.1e-6 | +0.02838 | 1.00465 |
| 65536 | 1.25 | −1.6e-6 | +0.03733 | 1.00539 |
| 65536 | 1.4 | −4.4e-6 | −0.01215 | 0.99843 |
| 65536 | 1.5 | +3.1e-7 | −0.06857 | 0.99176 |

Same offsets to 3 digits at `R = 256` and `R = 4096` (`x2_material_strain_results.json`), i.e. the
offset is `O(M)` and `L`-independent, uniformly in `λ ∈ [1, 3/2]`. So:

* LEMMA 1 (`a_λ(0) = (M/2)λL`) is **confirmed** by an instrument independent of all three of the
  attempt's routes;
* the ODE's coefficient `κ` is right **at the material point of the deformed configuration**, not
  only at the empty origin — a step the attempt left as part of GAP T. This is a genuine
  strengthening of the attempt, contributed here.

Also reproduced: `κ_δ` (0.49972 at 7.5°, 0.62252/0.73637 for the tapered `P_h(λ)` at λ = 1.25/1.5,
matching `s6` to 8 digits), the closed ODE, `θ₃₂ = 4(1−√(2/3))`, and both constants.

**Not verified by me, and still open:** my probe strains *every* shell by the *same* `λ`. The
integro-ODE strains shell `σ` by `λ(σ)` — a shell-dependent map that is not volume-preserving in
`ℝ³` and is not the image of a single `T_λ`. That, plus the Lipschitz stability, is GAP T.

## 8b. Corrected statement (what survives, stated so that it is defensible)

> **LEMMA 1 (exact, and elementary).** For `η₀ = −M sgn(z)h(φ)/r` on `ρ₀<|x|<R` and the uniform
> coaxial strain `T_λ`, the transported field is `λ` times the same `1/r` profile on a support
> whose log-radial thickness is `L = log(R/ρ₀)` in every direction; hence, `a` being linear in `η`
> with a scale-free kernel, `a_λ(0) = (M/2)P_h(λ)L`, and `P_1(λ) = λ` exactly. This is an
> algebraic consequence of `η₀` depending on `r` alone and of `T_λ` scaling `r` by a constant. It
> is **not** a persistence statement about the Navier–Stokes flow and it does not remove the
> brief's step (2); the perturbative content sits entirely in showing that the true flow map is
> `T_{λ(σ)}` up to relative `O(c/L)` and that `a[·](0)` is Lipschitz-stable under that
> perturbation (the seat's GAP T, unproved).
>
> **Extension (established here, not in the attempt).** The same conclusion holds at the tracked
> **material** point of the innermost shell, not only at the empty origin:
> `a(x_λ) = (M/2)λL + O(M)` with the `O(M)` offset bounded by `0.069` and independent of `L`,
> uniformly for `λ ∈ [1, 3/2]` at `R/ρ₀ = 256, 4096, 65536`.
>
> **Constants.** Conditional on GAP T, GAP V and R, the clock is `T(Re) ≤ (c₂+o(1))/log Re` with
> `c₂ = 4log(3/2) = 1.6219` from `Φ_h ≥ 1` alone. The sharper `c₂ = 8(1−√(2/3)) = 1.4680` follows
> from `Φ = λ` **within the `T_λ` model**, is asymptotic in `L`, and is contradicted by the only
> measurement of it at `L ≈ 2.9` (§6). It should be reported as a model value, not as the
> theorem's constant.
>
> **Not asserted.** No refutation of Bradshaw–Farhat–Grujić ARMA 2019 Thm 10 — correctly, since
> both remaining gaps are load-bearing. The attempt is right to say so; the `key` should say so
> too.

## 9. Files

| file | what it establishes |
|---|---|
| `x1_deflate_lemma1.py` / `_results.json` / `x1_out.txt` | LEMMA 1 symbolically + the Eulerian deflation (`G(t)` deleted ⟹ same answer); truncation table for `0.3958` |
| `gegen5d.py` | my independent 5D Gegenbauer strain instrument (Hill control analytic) |
| `x2_material_strain.py` / `_results.json` / `x2_out.txt` | fifth route to LEMMA 1; the material-point test the attempt skipped |
| `x3_offset_arbitrate.py` / `x3_out.txt` | φ-scan of the material offset; bulk identity `a = A(ρ) − (M/2)cos²φ + O(M)` control |
| `x4_gate_mutation.py` / `_results.json` / `x4_out.txt` | FL-043: static unfalsifiable checks + one-number-at-a-time mutation of the gate |
| `x5_empirical_status.py` / `_results.json` | the accelerated vs frozen clock against `viscous-numerics` |
| `copy/` | the attempt's own scripts, re-run by me; `gate_out.txt` = gate rc 0, 88 PASS |

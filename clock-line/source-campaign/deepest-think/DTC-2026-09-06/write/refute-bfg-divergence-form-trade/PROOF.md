# REFUTER report — `write/bfg-divergence-form-trade`

**Seat** `write/refute-bfg-divergence-form-trade`, DTC-2026-09-06.
**Target** `write/bfg-divergence-form-trade/PROOF.md` + its 8 scripts, re-run from a copy at
`./copy/` (never from the original tree; nothing outside this folder was written).
**Method** (i) hash-verify and re-execute every script of the target from a copy and diff the
output against its committed logs; (ii) re-derive every load-bearing constant by a route of
my own choosing that reads nothing from the target; (iii) run a falsification probe against
the one inequality in the note that is not an identity (Lemma U and its two halves);
(iv) verify every quotation against the cited source file; (v) test each stated *control* for
whether it can fail.

---

## 0. Verdict

| item | verdict |
|---|---|
| reproducibility (24/24 hashes, 8/8 scripts, logs byte-identical) | **CLEAN** |
| Lemma D (divergence form) | **SURVIVES** — re-proved here by hand and symbolically |
| Lemma D′ (`Q = 1`, Gram/PSD) | **SURVIVES** — re-proved; 12 000 000-sample independent search max `0.99999995107541062` |
| Lemma K (`K1`–`K8`) | **SURVIVES** — `∫|∇G_ν| = 2/√(πνt)` reproduced by two routes of mine (radial mpmath quad, rel err `< 1e-30`; a raw Cartesian Riemann sum, rel err `1.3e-06`) |
| Proposition 1 (`C = 2/√(πν)`) | **SURVIVES** — sign, index and integration-by-parts chain checked line by line |
| Lemma U (`C_u = 1.214239420800536`) | **SURVIVES** — re-derived symbolically and by my own bisection; **and survived a 37-shape falsification probe** (worst `Φ = 0.372777 ≪ C_u`) |
| Theorem C (`c_div = 3.622974321832e-03`, exponent `−1`) | **SURVIVES** — re-derived from scratch, symbolic residual `0`; end-to-end closure reproduced |
| claim N1 (no energy-free repair) | **SURVIVES** — the dilation family is correct |
| §7 comparison table, §7.2 crossovers, §9 datum table | **SURVIVES** — every number reproduced by my own arithmetic |
| **§6.2 / headline: "the brief's `Re_E^{-1/2}` is *dimensionally impossible*"** | **REFUTED — F1, MAJOR** |
| **§6.2 / `s7_audit.py` "the exponent claim re-derived a THIRD way"** | **FL-043 — F2**, a control that cannot fail, and it drops `c_div` by a factor `276` |
| **§9 "an unplanned cross-seat consistency check"; §7.1 "`C_u` cancels identically"** | **FL-043 — F3**, two controls that cannot fail |
| §10.5's enumeration of what makes `c_div` non-optimal | **INCOMPLETE — F4, MINOR** |
| §0 novelty item (i) | **OVERSTATED — F5, MINOR** |

**Overall: `refuted = true`, severity MAJOR.** The mathematics of the note — every lemma,
every constant, the theorem, and the answer to all four items of the brief — stands
unaltered under adversarial recomputation. What is refuted is a *false supporting claim
carried in the note's own headline*, together with the three "verifications" that were
built so that they could not have caught it.

---

## 1. What I re-ran, and that it reproduces

```
shasum -a 256 -c SHA256SUMS            -> 24/24 OK   (PROOF.md included)
python3 s1..s7, check_constants.py     -> exit 0 for all 8
diff <fresh stdout> <committed *.txt>  -> IDENTICAL for all 8
check_constants.py                     -> 50 PASS / 0 FAIL, 50/50 mutations caught
s7_audit.py                            -> 9 closed forms, 0 mismatches
```

The gate `check_constants.py` is honestly built for the part it covers: every stored number
is the *left* side, the recomputation route differs from the producing script's, and the
mutation test `v -> 1.5v+0.37` is dynamic. Two of its fifty "recomputations" are typed
literals (`s2 Q exact -> 1`, `s2 Q triangle -> 2`); those two are re-proved independently
below, so nothing turns on it. Its coverage is the JSON, not the prose: a numeral mistyped
into `PROOF.md` would not be caught by any script in the folder.

## 2. Independent re-derivation of every load-bearing constant (`r1_independent.py`)

My routes, none of which reads the target's files:

* `∫_{R³}|∇G_ν(·,t)|dx`: radial `mpmath.quad` on the exact density and, separately, a raw
  241³ Cartesian Riemann sum. Closed form `2/√(πνt)` confirmed at
  `(ν,t) = (1,1),(1,0.01),(7,3),(0.001,5)`:
  `1.12837916709551256`, `11.2837916709551251`, `0.246232521229829082`, `15.9576912160573077`
  (radial rel err `≤ 2.4e-51`; Cartesian rel err `1.3e-06`, its own discretisation floor).
  `‖G(·,1)‖₂ = 0.0890881837943737764 = (8π)^{-3/4}` confirmed numerically.
* `Q`: the Gram identity re-derived symbolically here (`objective − (1−c²)` on the
  constraint boundary `= 0`), plus **12 000 000** fresh uniform samples on `S²×S²×S²`,
  max `0.99999995107541062 ≤ 1`, and attainment `1.000000000000000` at `(e₁,e₂,(e₁+e₂)/√2)`.
* `γ`: my symbolic optimisation returns `5·2^{3/5}3^{2/5}/6 = 1.96013170420778926`, which is
  the note's `5·2^{-2/5}3^{-3/5}` (same number, different reduced form).
* `C_u`: symbolic `5·2^{9/10}3^{2/5}/(6π^{3/5}) = 1.21423942080053582`; independent
  bisection on `f′` gives the same to `< 1e-14`, at `τ_* = 0.104217595142273536`.
  Difference from the note's closed form: symbolic `0`.
* `c_div`: I re-ran the bootstrap symbolically from the definitions (`Q = 1`, `K3`, `K5`,
  cap `2M₀`, target `M₀/2`) and obtained
  `t_* = 3·3^{1/5}π^{11/5}ν/(12800 E₀^{2/5}M₀^{6/5})` and
  `c_div = M₀ Re_E t_* = 3·3^{1/5}π^{11/5}/12800 = 0.00362297432183177276`.
  Difference from the note: symbolic `0`. Also `π/(2^{46/5}C_u²)` agrees,
  `1/C_u² = 0.678252406508449757`, `π·2^{-46/5} = 0.00534163135591710948`,
  `c_div(Q=2) = 0.000905743580457943191`, ratio exactly `4`.
* Datum table (`L = 8.3178, 20, 50, 100`) and both crossovers (`Re_E = 12.883`, `42209.1`)
  reproduced digit for digit.
* End-to-end closure at the note's own test point: `Re_E = 69.9536750335`,
  `t_* = 1.3997581237509200e-05`, `Duhamel(t_*) = 1.8500000000000001 = M₀/2`,
  `M₀ Re_E t_* = 3.6229743218317728e-03`.

**Every constant tracks.** I found no arithmetic error anywhere in the folder.

## 3. Falsification probe of Lemma U (`r2_lemmaU_probe.py`)

Lemma U is the only genuinely *inequality*-shaped step (everything else is an identity or a
pointwise algebraic bound), so it is the one worth attacking. The functional
`Φ[u] := ‖u‖_∞/(E^{1/5}M^{3/5})` is invariant under `u ↦ au` and `u ↦ u(λ·)` (verified:
`Φ(u) = Φ(3.7u) = 0.2838505214`), so one probe per *shape* suffices. I built 37 smooth
compactly-concentrated divergence-free fields as `u = curl A` on a `128³` box (relative
`|div u| ≤ 1.3e-13`): 7 Gaussian shapes (isotropic, oblate `3,3,1`, prolate `1,1,4`, `4,1,1`,
`3,1,2`, `1,1,1`, `4,4,4`) and 30 random multi-blob fields.

```
worst  Phi = 0.372777   (random blobs #17)      vs   C_u = 1.214239   ->  NO VIOLATION
near-half  ||u - e^{tau Lap}u||_inf / ((4/sqrt pi) M sqrt tau)   max ratio 0.201327   (15 tests)
far-half   ||e^{tau Lap}u||_inf / ((8 pi tau)^{-3/4} E^{1/2})    max ratio 0.428882   ( 8 tests)
best tightness  ||u||_inf / (Lemma-U bound)                      0.2912   (12 shapes)
```

Lemma U and both of its halves survive. The probe also measures something the note does not
report — see F4.

## 4. Sources checked at their cited location

* `sharp/exact-first-order/NOTE.md` §3 — `E := ‖u‖₂² (no 1/2)` (l.157),
  `Re_E = E^{2/5}M^{1/5}/ν`, `c1 := t_d·M·log Re_E` (l.175),
  `E_shell = 0.172403978 M²R⁵`, `log Re_E = 2 log(R/ρ₀) − 0.7031660`,
  `c1 = 4 ln 2 − (2 ln 2)(0.7031660)/log(R/ρ₀) → 4 ln 2 = 2.7725887222` (ll.235,237).
  **All four quotations are faithful.**
* `RETURN_ADDENDUM_8_CLOCK_LINE.md` — `(S1) PROVED: 𝒯(Λ) ≥ c₁/(1+log₊Λ)` (l.5);
  `(S3) CONJECTURE 𝒯(Λ) ≤ c₂/log Λ: SKETCH` (l.7); the BFG verdict `(b)` and the symmetric
  stake (l.10). **Faithful.**
* `gaps/kukavica-at-source/NOTE.md` — `T ≥ 1/(C‖u₀‖²_{L^∞})` (Thm 2.2 / Prop 3.2) at l.37;
  `sup|∇G|(|x|+√t)⁴ = 0.9231098779`, `sup G(|x|+√t)⁴/√t = 0.7109819926` at l.274;
  `(4π/3)t` mean-one vs mean-zero at l.45; the accepted-manuscript caveat at ll.21,288.
  **Faithful, and the note carries the caveat forward as required.**
* `verify-astra/refuter-A1-entropy-clock/r5_final.py` l.38 —
  `T <~ nu/(M^2 ell^2) = 1/(M*Re)`; `r3_closure.py` ll.47–50, including the same
  `Re/(1+log Re)` ratio. **Faithful; the prior-art disclosure in §0 is correct.**
* The *brief* for this key is not on disk anywhere under `DTC-2026-09-06`, so the note's
  quotation of the brief's guess `τ ≥ c/(M·Re_E^{1/2})` could not be checked at source.
  This does not affect the mathematics; it is recorded because F1 is *about* that guess.

---

## 5. F1 — MAJOR. The headline "dimensionally impossible" claim is false

The note asserts, in its one-line answer, in the §0 verdict table, and as a blockquote in
§6.2:

> "The brief's `Re_E^{-1/2}` is dimensionally impossible: `M t_* = ν E^{-2/5}M^{-1/5} = Re_E^{-1}`
> identically, and `M t_*/Re_E^{-1/2} = √ν/(E^{1/5}M^{1/10})` is **not a pure number**."

Both halves of that sentence are wrong.

**(a) `√ν/(E^{1/5}M^{1/10})` *is* a pure number — it is `Re_E^{-1/2}` itself.**
Symbolically (`r3_findings.py`):

```
    sqrt(nu)/(E^{1/5} M^{1/10})  -  Re_E^{-1/2}   =   0        (sympy, exact)
```

and with `[E] = L⁵T^{-2}`, `[M] = T^{-1}`, `[ν] = L²T^{-1}` its dimension evaluates to `1`.
The note's own §6.3 proves `Re_E` is dimensionless; `Re_E^{-1/2}` is therefore dimensionless
too, and the note contradicts itself two subsections apart.

**(b) A law `M T ≥ c·Re_E^{-1/2}` is dimensionally admissible.** `[M·T] = 1` and
`[Re_E^{-1/2}] = 1` (both computed, not asserted). Nothing forbids such a law on dimensional
grounds; it is simply a *stronger* statement than the one this route proves, for every
`Re_E > (c/c_div)²`:

```
    Re_E = 1e2   c_div/Re_E = 3.622974e-05   c_div/sqrt(Re_E) = 3.622974e-04
    Re_E = 1e6   c_div/Re_E = 3.622974e-09   c_div/sqrt(Re_E) = 3.622974e-06
    Re_E = 1e17  c_div/Re_E = 3.622974e-20   c_div/sqrt(Re_E) = 1.145685e-11
```

**(c) The displayed identity is off by `c_div`.** The note writes
`M t_* = ν E^{-2/5}M^{-1/5} = Re_E^{-1}` *identically*, three lines after Theorem C states
`M₀ t_* = c_div/Re_E` with `c_div = 3.622974321832e-03`. The two displays disagree by a
factor `1/c_div = 276.0`. A constant that does not track, inside one subsection.

**What survives.** The *conclusion* is correct and is independently re-derived here: within
this scheme the exponent is exactly `−1`, because `√t·M^{8/5}E^{1/5}/√ν ≲ M` forces
`t ≲ ν M^{-6/5}E^{-2/5}`. The brief's `−1/2` does not follow from this route. That is all
that may be said.

**Corrected statement for §6.2 and the headline:**

> The Reynolds exponent produced by this route is exactly `−1`: `M₀ t_* = c_div/Re_E`.
> The brief's `Re_E^{-1/2}` is *not obtainable from this estimate* — it is strictly stronger
> than what the divergence-form bootstrap yields, by a factor `Re_E^{1/2}`. It is **not**
> ruled out on dimensional grounds: `Re_E^{-1/2}` is dimensionless, and
> `√ν/(E^{1/5}M^{1/10})` *is* `Re_E^{-1/2}`. Ruling out a `−1/2` clock for NS in general is
> not attempted here and does not follow from anything in this note.

## 6. F2 — FL-043. `s7_audit.py`'s "third way" cannot fail

`s7_audit.py` §"the exponent claim, re-derived a THIRD way (pure dimensional analysis)"
begins

```python
    t = nu*M**Rational(-6,5)*E**Rational(-2,5)      # <- the answer, typed in
    assert sp.simplify(M*t - Re**-1) == 0
    print(" M t_* / Re^{-1/2} = ...  (not a constant -> the -1/2 law is dimensionally impossible)")
```

Three defects, all confirmed in `r3_findings.py`:

1. It does not read `t_*` from Theorem C; it re-types the exponents whose correctness is the
   thing under test. If Theorem C's exponents were wrong the assertion would still pass.
   **A control that cannot fail.**
2. Its input differs from Theorem C's `t_*` by exactly `c_div = 3.622974321832e-03`
   (`ratio t_*/t_s7 = 3·3^{1/5}π^{11/5}/12800`, verified symbolically), which is how the false
   identity of F1(c) got into the prose.
3. The printed clause is a non-sequitur: the computation establishes "not a **constant**",
   the printed conclusion says "dimensionally impossible". No `assert` guards that clause —
   it is a string, and it is the one sentence of the section that reached the headline.

Note also that the *sound* exponent extraction, `s4_clock.py`'s overdetermined `3×2` solve,
is not the independent check the note calls "the real check": with three dimensional inputs
`(E,M,ν)` and two dimensions, Buckingham gives exactly one dimensionless group, so *any*
dimensionally correct `t_*` is automatically of the form `M^{-1}·c·Re_E^{β}` and the system
is automatically consistent. It does detect the brief's `−1/2`, so it is not empty; but it
is dimensional analysis, not a second derivation.

## 7. F3 — FL-043. Two more controls that cannot fail

**(a) §9's "unplanned cross-seat consistency check."** The note reports
`model/(a) → 4 ln 2 = 2.7725887222` and calls the agreement with
`sharp/exact-first-order`'s `c1` "a consistency check on this note against that seat". It is
not. The ratio is

```
    model/(a) · c_1  =  (M t_d)(1 + log Re_E)  =  (2 ln2 / L)(1 + 2L - 0.7031660)
```

— built from **two numbers quoted from that one seat** and nothing else. `c_div` does not
appear; perturbing `c_div` by `10^6` leaves every entry of the column unchanged. And
`sharp/exact-first-order/NOTE.md` l.175 *defines* `c1 := t_d·M·log Re_E`, so the limit is
that definition evaluated on its own inputs. Reproduced (`r3_findings.py`):
`2.8220608543 / 2.7931636873 / 2.7808187082 / 2.7767037152 / 2.7725891337` at
`L = 8.3178, 20, 50, 100, 10⁶`. Nothing about this note could have made it disagree.
There is also a name collision: the `c_1` of clock (a) (Astra, unpinned) and the `c1` of
`exact-first-order` (`= 4 ln 2`) are different constants, and §9's prose
"`model/(a) → 4 ln 2`" silently drops the `/c_1` its own table carries.

**(b) §7.1's "`C_u` cancels identically (symbolic residual 0)."** Both sides were
constructed by dividing by `C_u²` — `c_div = π/(2^{46/5}C_u²)` and `c_GIM := c_G/C_u²` — so
the cancellation is structural. With an arbitrary symbol `X` in place of `C_u`,
`c_div(X)/c_GIM(X) = 2^{4/5}π/(1024 c_G)` and `∂_X` of it is `0` for every `X`
(`r3_findings.py`). The cancellation is a consequence of feeding both clocks the same
velocity bound, not evidence that the two arguments are the same clock. The functional-form
and exponent match is real; "IS the GIM/Kukavica clock" remains an interpretation, and with
`c_G` unpinned it cannot be checked at the level of constants.

## 8. F4 — MINOR. The list of what makes `c_div` non-optimal is incomplete

§10.5: *"`Q = 1` and `K3` are sharp, but the cap costs a factor `2^{18/5}`"*. My probe (§3)
measures `‖u‖_∞ / (Lemma-U bound) ≤ 0.2912` over 12 shapes — Lemma U is loose by a factor
`≥ 3.43` on everything I could build, and `c_div` carries `C_u^{-2}`, so **at least a factor
`≈ 11.8`** of the shortfall sits in Lemma U, which §10.5 does not list. The heat-split
constant `C_u` is *exact for that split*; it is not the sharp constant of the inequality, and
the note's phrase "Lemma U — the velocity bound, with its exact constant" invites the
opposite reading.

## 9. F5 — MINOR. Novelty item (i) is overstated

§0 claims as new "(i) the mean-zero (divergence-form) route specifically". But
`RETURN_ADDENDUM_8_CLOCK_LINE.md` l.10 already records it, together with the velocity bound
it needs: *"Its divergence-form repair (div u = div ω = 0) trades the logarithm for an energy
dependence (‖u‖_∞ ≤ CE^{1/5}M^{3/5})"*. The route and the shape of Lemma U were both on the
estate record; what is new is the execution with constants, which is items (ii)–(iv) and is
correctly claimed. §10.1 is honest ("nothing here is a new theorem"); §0 item (i) should be
narrowed to match it. Similarly, §7.2's `Re/(1+log Re)` ratio is already in
`verify-astra/refuter-A1-entropy-clock/r3_closure.py` ll.47–50.

---

## 10. What I did **not** find

* No error in Lemma D, Lemma D′, Lemma K, Proposition 1, Lemma U, Theorem C, or N1.
* No constant that fails to track. Every numeral in the note's `constants` block that is
  claimed as derived-in-folder reproduced under my own route, most of them exactly.
* No sign or index error in the mild representation. I checked
  `∂_t ω_j − νΔω_j = −∂_i A_{ij}`, `∂_{y_i}[G(x−y)] = −(∂_iG)(x−y)`, and
  `(A^T v)_j = (v·u)ω_j − (v·ω)u_j`, all by hand and symbolically.
* No misquotation of any source file (four sources checked at their cited lines).
* No hidden hypothesis: (H1), (H2), the `e^{σΔ}u → u` qualitative step, the stronger-than-BFG
  energy assumption, the un-obtained Kukavica source, and the Picard-iterate step are all
  stated in §1 or §10.
* The reproduction discipline is exemplary: 24/24 hashes, 8/8 scripts, byte-identical logs.

## 11. Remaining gap after this refutation

The note's four brief items stand as written *except* for the justification of item (2)'s
second half. What remains open is unchanged from the note's own §10 (`c_1` unpinned, `c_G`
unpinned, hypotheses stronger than BFG's, `c_div` not optimal, BFG Thm 8/10 not
adjudicated), plus one item this refutation adds: **no lower bound on the achievable
Reynolds exponent has been established.** The note proves `−1` for *this* estimate and
believed it had also excluded `−1/2` on dimensional grounds; it has not. Whether some other
route inside finite-energy hypotheses reaches `M₀T ≳ Re_E^{-1/2}` — or `Re_E^{-θ}` for any
`θ < 1` — is open, and is the natural successor unit of work alongside pinning `c_1`.

## 12. Files

| file | what it is |
|---|---|
| `PROOF.md` | this report |
| `r1_independent.py` / `r1_log.txt` / `r1_results.json` | my own route to `∫|∇G|`, `Q`, `γ`, `C_u`, `τ_*`, `c_div`, `1/C_u²`, `π2^{-46/5}`, the datum table, the crossovers, and the dimension of `√ν/(E^{1/5}M^{1/10})` |
| `r2_lemmaU_probe.py` / `r2_log.txt` / `r2_results.json` | 37-shape falsification probe of Lemma U and of both halves of the heat split |
| `r3_findings.py` / `r3_log.txt` / `r3_results.json` | F1, F2, F3a, F3b pinned symbolically; end-to-end closure of Theorem C |
| `copy/` | byte copy of the target folder, from which all re-runs were made |
| `SHA256SUMS` | computed over the above |

Reproduce with `python3 r1_independent.py; python3 r2_lemmaU_probe.py; python3 r3_findings.py`
(sympy 1.14.0, mpmath 1.3.0, numpy 1.26.4); target re-run with
`cd copy && shasum -a 256 -c SHA256SUMS && for f in s1_kernel s2_algebra s3_velocity s4_clock s5_datum s6_no_energy_free s7_audit check_constants; do python3 $f.py; done`.

---

## 13. Second refuter pass (`r4_recheck.py`, `r4_log.txt`, `r4_results.json`, `copy2/`)

An independent re-run of this refutation, from a *fresh* copy of the target (`copy2/`),
with a script that reads nothing from either the target folder or from `r1`–`r3`.

**Reproduction, again.** `shasum -a 256 -c SHA256SUMS` in `copy2/` → **24/24 OK**; all eight
target scripts exit `0`; fresh stdout is **byte-identical** to all eight committed logs;
`check_constants.py` → 50 PASS / 0 FAIL, 50/50 mutations caught; `s7_audit.py` → 9 / 0
mismatches.

**Constants, re-derived a third time by my own route** (`r4_log.txt` §§A–D):
`∫|∇G_ν| = 2/√(πνt)` by radial integration *and* by the probabilistic route `E|X|/(2νt)`
with `X ∼ N(0,2νt·I₃)` — symbolic residual `0`; `‖G(·,τ)‖₂ = 2^{3/4}/(8π^{3/4}τ^{3/4})
= (8πτ)^{-3/4}`; the Gram identity `det + objective − (1−c²) = 0` exactly, with a fresh
3 000 000-sample search maxing at `0.999999896043568 ≤ 1`; `γ = 5·2^{3/5}3^{2/5}/6
= 1.9601317042077893`; `C_u = 5·2^{9/10}3^{2/5}/(6π^{3/5}) = 1.2142394208005358`, matched by
my own **bisection on `f′`** (not golden section) at `τ_* = 0.104217595142273536`, agreeing
to `< 1e-25`; `t_* = 3·3^{1/5}π^{11/5}ν/(12800 E^{2/5}M^{6/5})`;
`c_div = 3·3^{1/5}π^{11/5}/12800 = 3.622974321831772763e-03`, symbolic residual `0` against
both `π/(2^{46/5}C_u²)` and the reduced form, and `∂_M c_div = ∂_E c_div = 0`.
End-to-end closure at **a test point of my own choosing** (`M=0.83, E=11.9, ν=0.207`,
`Re_E = 12.533260244131`): `t_* = 3.482756456330e-04`, `Duhamel(t_*) = 0.415000000000000
= M/2` (residual `< 1e-40`), `M·Re_E·t_* = 3.622974321831772763e-03 = c_div`. **Every
constant in the note tracks.**

**F1 re-confirmed, harder.** Symbolically `√ν/(E^{1/5}M^{1/10}) − Re_E^{-1/2} = 0`, and
under `[E]=L⁵T^{-2}`, `[M]=T^{-1}`, `[ν]=L²T^{-1}` the substitution returns, with no `L`
or `T` surviving: `[Re_E] = 1`, `[Re_E^{-1/2}] = 1`, `[√ν/(E^{1/5}M^{1/10})] = 1`, and
`[M·t_*] = 3·3^{1/5}π^{11/5}/12800` — i.e. `M t_*` reduces to the *pure number* `c_div`,
not to `Re_E^{-1}`. That single substitution refutes both halves of the note's blockquote at
once: the ratio it calls "not a pure number" is dimensionless (it is `Re_E^{-1/2}` itself),
and the display it calls an identity is off by `1/c_div = 276.0163090238`.

**F6 — FL-043, new in this pass.** `s7_audit.py`'s ninth "closed form" is

```python
eq("C_w = 3*2^(1/6)/(2 pi^(2/3))",
   3*2**sp.Rational(1,6)/(2*sp.pi**sp.Rational(2,3)),
   3*2**sp.Rational(1,6)/(2*sp.pi**sp.Rational(2,3)))
```

— the **same sympy expression on both sides** (`srepr` identical, verified). It reports `OK`
for any value whatever; it is a control that cannot fail. Classifying all nine:
**five substantive** (`γ`'s reduced form, `C_u` via the `γ` route, the Duhamel prefactor,
`c_div`'s two forms, `1/C_u²`), one notational (`π/2^{46/5}` vs `π·2^{-46/5}`), two trivial
exponent arithmetic (`2·2^{18/5} = 2^{23/5}`, `(2^{23/5})² = 2^{46/5}`), one tautology.
So the status line *"`s7_audit.py` re-derives all 9 closed forms quoted in `PROOF.md`,
0 mismatches"* overcounts by nearly two to one — and the one constant it audits
tautologically is `C_w`, the very numeral §5 and gap 5 flag as unverified. (`C_w`'s numeral
*is* covered elsewhere: `check_constants.py`'s `rc_Cw` re-minimises it independently. Nothing
mathematical turns on this; the overstated audit count does.)

**F2, F3a, F3b re-confirmed.** `t_*(Theorem C) / t(s7's "third way") = c_div` exactly, so
that section is structurally blind to Theorem C's constant. With an arbitrary symbol `X` in
place of `C_u`, `c_div(X)/c_GIM(X) = 2^{4/5}π/(1024 c_G)` with `∂_X ≡ 0` — the cancellation
is structural. The cross-seat "check" `(2 ln2/L)(1+2L−0.7031660)` gives
`2.8220608543 / 2.7931636873 / 2.7808187082 / 2.7767037152` at `L = 8.3178,20,50,100` and
`2.7725891337` at `L = 10⁶`, from the two quoted formulas alone; and the quoted offset is
itself `(2/5)·ln(0.172403978) = −0.7031659387`, i.e. the second quoted number is a function
of the first. `sharp/exact-first-order/NOTE.md` l.175 *defines* `c1 := t_d·M·log Re_E`, so
the limit is that definition evaluated on its own inputs.

**F5 re-confirmed at source.** `RETURN_ADDENDUM_8_CLOCK_LINE.md` l.10, verbatim:
*"Its divergence-form repair (div u = div ω = 0) trades the logarithm for an energy
dependence (‖u‖_∞ ≤ CE^{1/5}M^{3/5})"*. Both the route and the shape of Lemma U were on the
estate record before this note; the execution with constants was not.

**Sources re-checked in this pass at their cited lines:** `RETURN_ADDENDUM_8` ll.4–10
((S1), (S3), the BFG verdict (b), the divergence-form sentence); `sharp/exact-first-order`
ll.175, 235, 237 (`c1` definition, `E_shell`, `log Re_E`, `4 ln 2`);
`gaps/kukavica-at-source` ll.21, 37, 274, 288 (`T ≥ 1/(C‖u₀‖²_∞)`, the two majorant
numerals, the accepted-manuscript caveat). **All faithful.**

**Verdict unchanged: `refuted = true`, severity MAJOR** — one false claim carried in the
note's headline (F1), four controls that cannot fail (F2, F3a, F3b, F6), two minor
overstatements (F4, F5), and no error anywhere in the mathematics.

# REFUTER of `s3close/round2/u2/PROOF.md` — BLOCK 3, the hypothesis (Γ-off)

Single refuting seat, sitting of 2026-09-08. Laws: `TORMENT NEXUS/LAWS.md` (first 120 lines).
Everything below was re-derived or re-run in a **copy** of `u2/` at
`/private/tmp/claude-501/-Users-spaceman-Desktop/e34d31b2-1bd6-4688-bd1c-674bd591026e/scratchpad/refute-u2/`.
`u2/` was not edited. Where I checked a constant I wrote my own instrument from the
definition in the source lemma (`rebuild/far-near-kernel-lemma/NOTE.md`,
`gaps/gap-V-aronson/NOTE.md`, `lower/prove-lagrangian/NOTE.md`,
`s3close/assembly/ASSEMBLY.md`), not from `u2`'s script.

---

## 0. Verdict

**NOT REFUTED. (Γ-off) with `C₁ = 0` and a finite, `L`-free `C″` SURVIVES.**

The chain `(2.1) → P1/P2 → Riesz → far/near` is sound. Every exact identity, every kernel
constant and every displayed table entry that I could recompute reproduces. The two genuinely
new moves — `C₁ = 0` from `κ ≥ 0`, and the axis-safe collar `3.9992184 e^{c_R}E₀M` from P1 —
are correct, and they do what the note says they do: they remove the `1/sin φ` divergence and
make the statement global on `ℝ⁵`.

What does not survive is a set of **claims about the proof**, not the proof:

* **2 MAJOR** — the sharpness claim for P1/P2 (which is a tautology in the script and is false
  for the note's own datum), and a datum-scope gap (the quoted `C″`/`L_*` do not apply to the
  datum `(S3)` is actually stated for).
* **8 MINOR** — an `L_*` that is the Picard-convergence boundary rather than the existence
  boundary, a false `‖η₀‖_∞` row, a support argument that covers only the un-mollified datum,
  two unstated hypotheses, an `a(0,s) ≤ (λ/2)ML` that is not exact for `s > 0`, an inconsistent
  crossing angle, and numerics that never probe the region the note's novelty is about.
* **0 FATAL.**

Grade unchanged in substance: the note is **PROVED modulo the approximate-maximum argument**,
with the two extra regularity/decay hypotheses of MINOR-4 added to (H1)–(H5).

---

## 1. What I re-derived and confirmed (NONE)

| item | `u2` | this seat | agree |
|---|---|---|---|
| `∫_{ℝ⁵}|x−x′|^{-4}|x′|^{-2}dx′ = (π⁴/2)/|x|` | `48.70454551700121` | `48.704545517001218618` (40 dps, closed-form angular integral + `mpmath`) | `9.8e-36` |
| `|S⁴| = 8π²/3`, `C_K|S⁴| = 1`, `3π²/16 = C_K·π⁴/2` | `1.8505508252042546` | same | exact |
| `div₅b = 2a`; `∇₅b = a·diag(1,1,1,1,−2)+E`; `(2.1)` | residual `0`, worst ratio `0.998932` | re-derived by hand from the 3-D `∇u` in `(ê_r,ê_θ,ê_z)`; `u1` re-run: residuals `0` | ✓ |
| `Δ₅(ρ^k g)` rearrangement, sign of the `W/ρ²` term | `−2/ρ²` both `k = 1, 2` | general `n`: coefficient is `k(k−n+2)`; at `n = 5` that is `k(k−3) = −2` for `k = 1` **and** `k = 2` | ✓ and it is an `n=5` fact |
| `Λ_max(sym∇₅b) = max(a, ½[−a+√((3a+2Q)²+(2P−ω)²)])` | rel `1.0e-15` | derived: block trace `−a`, `(tr)²−4det = (3a+2Q)²+(2P−ω)²` | ✓ |
| envelope (4.1), third branch `2Ĉ_aM+2Ĝ+λM/2` | `0` violations | needs `a ≥ −Ĉ_aM`, which **is** available: `a(x) ≥ a_far(0,s) − Ĉ_aM ≥ −Ĉ_aM` since `κ ≥ 0`. Without it the branch is `2A+2Ĝ+…` and `p → 3`, so this step is load-bearing for the `p ≈ 2` column and it holds | ✓ |
| far Taylor remainder, `z`-odd | `0.29199853254882` | `0.291998532549` | ✓ |
| inner multipole, `z`-odd | `0.014754271582666` | `0.0147542715827` | ✓ |
| `2R_A = 2(2⁵−2^{-5})^{1/5}` | `3.9992184446453` | `3.99921844465` | ✓ |
| `Ĉ_a(λ=1,σ_*=1/2)` | `8.697888775` | `8.69788877` = `ASSEMBLY §2.4`'s `C_a` | ✓ |
| `E₀ = sup ρ|η₀|/M` | `7.66060196237754` | `7.66060196237754` (`= tanh(1/w)/sin δ`, attained on the axis) | ✓ |
| `𝔊₀ = sup ρ²|∇η₀|/M` | `20.919707021704` | **`20.9197121979`** — the true sup is the bulk value `max_φ √(h²+h′²)` at `φ = 6.6641°`; `u2`'s grid value is `6.5e-8` low and its reported argmax radius (`160.75` / `1.66e15`) is an artefact of a maximum that is flat in `ρ` | ✓ (immaterial) |
| collar branch crossing | `sin φ = 0.05834932`, `3.345°` | `0.05834932`, `3.3451°` (at `e^{c_G}`) | ✓ arithmetically; see MINOR-6 |
| `§8.1` frozen table (`Ĉ_a`, `Ĝ`, `C″`) | all rows | reproduced from an independently written fixed-point solver: `214.38 / 316.97 / 468.56`; axis `335.46 / 464.84 / 650.18`; `Ĝ(p=3) = 440.97/810.10/1488.26`; `Ĝ(map μ=.05) = 96.03/183.77/324.11` | ✓ exact |
| `§8.2` `C″(10⁴)` every row | `232.82 361.48 582.24 / 362.85 526.76 801.05 / 3604.69 / 329.53 342.05 361.34 / 658.12 / 431.16 / 557.90` | identical to 2 dp | ✓ |
| `§9.1` `κ` feature table | 7 rows | all 7 reproduce to 12 digits; `κ_δ(7.5°) = 0.4997212305210886` via `½P_{h_δ}(1)` | ✓ |
| `𝒥_ac` | `65.62590996` | (`u2`'s own value; `hk2` reports `65.625910`) | ✓ |
| `u4` re-run from a copy | — | **bit-identical** `u4_results.json` and log | ✓ |
| `check_constants.py` | 287 checks | re-run in the copy: `ALL CHECKS PASS` | ✓ |

**`§9`'s two structural claims hold in the re-run.** `max Γ − 2a(0) = −0.52347` at `L = 10`
and at `L = 40` (`λ=1`), `−0.81629`/`−0.81629` (`λ=3/2`) — `L`-independent to five digits;
`max r|∇a|` moves `3.5e-5` (`λ=1`) and `1.3e-8` (`λ=3/2`) between `L = 10` and `L = 40`.
The `457.5×` slack is `440.97/0.96387` and is arithmetically correct; its decomposition
`e^{2c_G}=11.39 × |S⁴|𝔊₀/𝒥=8.39 × 4.79` is consistent, but the third factor is a **residual
defined a posteriori**, not a measured quantity — and see MAJOR-1, which shows the first
factor is not the "sharp" one the note says it is.

**The by-product finding of `§9.1` STANDS.** `ASSEMBLY §1.3`'s `ε_δ = 5.5754·10^{-4}` is
`1 − 2κ` for the **taper alone**; `§1.1`'s datum also carries `min(1,|φ−π/2|/δ_m)` with
`δ_m = 5°`. Recomputed here from the definition `κ = (3/2)∫₀^{π/2}h(φ)sin²φ cos φ dφ`:

```
  taper 7.5° alone            kappa = 0.499721230521   1-2k = 5.5753896e-4   floor 5.5784998e-4
  linear equatorial 5° alone  kappa = 0.498101207666   1-2k = 3.7975847e-3
  ASSEMBLY 1.1 datum (both)   kappa = 0.497822438187   1-2k = 4.3551236e-3   floor 4.3741737e-3
  campaign (D-B)/(D-C)        kappa = 0.474734355388   1-2k = 5.0531289e-2   floor 5.3220594e-2
```
Factor `7.841` above the displayed floor for the sealed datum, `95.4` for the run datum.
Confirmed; `δ_m` is the term that sets the floor, not `δ`.

---

## 2. MAJOR-1 — "P1, P2 exactly sharp on the reference strain" is a tautology, and is false
for the note's own datum

`§0` item 4, `§5` Sharpness, and the `§11` row *"P1, P2 exactly sharp on the reference strain —
**PROVED** (rel. err `0.0`, both)"*.

**What the script does.** `u5` block (c) sets `P1_maxprinciple = exp(∫Γ_rad) = λ` and
`P1_exact = lam` — the latter **hard-coded** — then reports `|λ−λ|/λ = 0.0`. Same for
`P2_exact = lam**4`. No field is constructed and nothing is transported. The `0.0` is the
identity `e^{log λ} = λ`, not a measurement. (Block (d) does compute `‖T_λ‖ = λ` and
`‖T_λ^{-1}‖ = λ²` correctly, but those are the two factors in the *upper* bound.)

**Why the claim is wrong as applied.** Equality in P1 needs `sup|α||η₀(α)|` to be attained
where `T_λ` expands, i.e. in the `y`-plane. For datum (D-C) it is attained **on the axis**,
which `T_λ` contracts by `λ^{-2}`. Equality in P2 needs `sup|α|²|∇η₀|` at the equator (where
`∇η₀ ∥ ê_z` and `α` lies in the `y`-plane); `u2`'s own `§12` item 3 records that `𝔊₀` is
attained at `φ = 6.660°`, in the taper transition. Measured here on the exact transported
field `η_λ = η₀∘T_λ^{-1}` for datum (D-C):

| | `λ = 1.25` | `λ = 3/2` |
|---|---|---|
| exact `sup ρ|η_λ|` | `4.902785 M` | `3.404712 M` |
| P1 bound `e^{c_R}E₀M = λE₀M` | `9.575752 M` | `11.490903 M` |
| **P1 loss** | **`1.9531 = λ³`** | **`3.3750 = λ³`** |
| exact `sup ρ²|∇η_λ|` (`= 5λ⁴M`, attained at the equator) | `12.207026 M` | `25.312489 M` |
| P2 bound `e^{c_G+2c_R}𝔊₀M = λ⁴𝔊₀M` | `51.073516 M` | `105.906043 M` |
| **P2 loss** | **`4.1839`** | **`4.1839`** |

The `λ³` is exact and structural: the maximiser stays on the axis, so the true factor is
`λ^{-2}` where the principle allows `λ`.

**Consequences inside the note.** `§5`'s *"all of the loss in the theorem sits in
`Γ_rad ≤ Γ`, none in the weights"* is false — at `λ = 3/2` the weights lose `3.375×` and
`4.184×` on the reference strain itself. `§10` item 1's slack attribution and `§9`'s
decomposition of the `457×` are correspondingly misallocated by about `14×`.

**Corrected statement.** *P1 and P2 have the correct exponential rate on the pure
axisymmetric strain: `e^{c_R}` and `e^{c_G+2c_R}` are attained by a datum whose weighted sup
sits in the expanding directions (for P1) and at the equator (for P2). For datum (D-C) neither
is attained; the constants are lossy by `λ³` and `4.184` respectively.* The `§11` row should be
**PROVED-for-a-witness-datum / NOT-SHARP-for-(D-C)**, and the `u5` (c) check should be struck
as circular.

Not fatal: P1 and P2 remain valid upper bounds, and the direction of the loss is safe.
It is, however, an identified `≈14×` of headroom in `C″`.

---

## 3. MAJOR-2 — the quoted `C″` and `L_*` are for a datum `(S3)` is not stated for

`u2` works with **(D-C)**: `tanh` taper, `tanh` equatorial mollifier, `tanh` radial ramp in
`log ρ`. `(S3)` as stated in `ASSEMBLY §1.1` uses
`min(1, φ_ax/δ)·min(1, |φ−π/2|/δ_m)` with `δ = 7.5°`, `δ_m = 5°`. `BLOCK 6` forces the
**radial** cutoff to be mollified; it says nothing about the two **angular** profiles.

`C″` is dominated by `Ĝ = (3π²/16)e^{pc_G}𝔊₀` (`440.97` of `468.56` at `λ = 3/2`, `p = 2`,
`σ_* = 1/2`), and `𝔊₀` is a pure functional of the angular profile:

```
  (D-C)  tanh taper + tanh equator :  G0 = 20.9197122   (attained at phi = 6.664 deg)
  ASSEMBLY 1.1 angular family      :  G0 = 58.6948961   (attained at the kink phi = delta = 7.5 deg,
                                                         where h = 1/sin d = 7.6613 and
                                                         h' = -cos d/sin^2 d = -58.13)
```
a factor **`2.806`**. So for the datum the assembly's theorem is about, `Ĝ` and hence `C″` are
`≈2.7×` larger (`C″ ≈ 1.3·10³` at `λ = 3/2`, `σ_* = 1/2`, frozen feedback), and `L_*` moves up
with it. `E₀` barely moves (`7.6613` vs `7.6606`), so `Ĉ_a` is unaffected.

The note nowhere states that its headline numbers do not transfer to `§1.1`'s family. `§12`
item 3 comes close ("a datum with a wider taper would lower `𝔊₀`") but does not price the
family actually in the seal. **Corrected statement:** *the theorem is datum-uniform in
structure; its constant is not. `C″ = 582` is `(D-C)`'s number. For `ASSEMBLY §1.1`'s
piecewise-linear angular profiles `𝔊₀ = 58.695` and `C″ ≈ 1.3·10³`; a `C¹` angular profile is
required if the `582` is to be quoted, and `§1.1` must be amended to a `tanh` (or otherwise
`C¹`) taper as well as a `tanh` radial ramp.*

---

## 4. MINOR findings

**MINOR-1 — `L_*` is the Picard-convergence boundary, not the existence boundary.**
`§8` says `u3` "bisects for the smallest `L` at which it exists". `u3.solve` bisects on whether
its simultaneous iteration in `(Ĝ, c_G, p)` **converges**. Solving the actual closure condition
directly — find `Γ̄` with `F(Γ̄) := λL + C″(c_G(Γ̄)) ≤ Γ̄`, `c_G = Γ̄c/L` — gives

| column (`λ = 3/2`) | `u2`'s `L_*` | existence boundary | over-statement |
|---|---|---|---|
| proved, `σ_* = 1/2` | `4997.6` | **`4904.7`** | `+1.89 %` |
| proved, `σ_* = 0` | `5313.8` | **`5273.0`** | `+0.77 %` |
| proved, `σ_* = sin 7.5°` | `5088.2` | **`5012.7`** | `+1.51 %` |

All `λ = 1`, `λ = 1.25`, `p3` and `map` rows agree exactly (`2189.7`, `3277.2`, `9909.2`,
`670.0/696.2/736.5`). At `L = 4997.6` the smallest fixed point is `C″ = 1035.64`,
`c_G = 1.3844`, `p = 2.3544` — matching `§9`'s boundary quote, so when the iteration converges
it converges to the smallest root, as `§12` item 6 hopes. Direction: **conservative**. The
wording of `§8` and the row `L_* = 4997.6` should say *smallest `L` at which the iteration
converges*, and the existence boundary is `4904.7`.

**MINOR-2 — the `‖η₀‖_∞` row is false for (D-C).**
`§3` and the `§11` row *"brief-layer correction: `‖η₀‖_∞ = M/(ρ₀ sin δ)`, not `M/ρ₀` —
**PROVED** (`u2`)"*. What `u2` computes is `E₀ = sup ρ|η₀|/M = 7.66060196`, a different
functional. The actual sup norm for (D-C) is
```
  ||eta_0||_inf = M/rho0 * sup_rho (rho0 Theta(rho)/rho) * sup_phi h(phi)
                = 0.53431477 * 7.66060196 = 4.093173 M/rho0     (attained at rho = 1.63761 rho0, on the axis)
```
against `M/(ρ₀ sin δ) = 7.661298 M/ρ₀`. The correction's *direction* is right (the taper, not
the radial edge, sets the constant; the brief's `M/ρ₀` is too small by `4.09`), and
`M/(ρ₀ sin δ)` is a valid majorant — for the **sharp-edge** datum it is the value. For (D-C)
it is not. Row should read `‖η₀‖_∞ ≤ M/(ρ₀ sin δ)`, with the (D-C) value `4.0932 M/ρ₀`.

**MINOR-3 — `§5.4(b)`'s support argument covers only the un-mollified datum.**
The argument is: `Φ_{s#}|Dη₀|` lives on `Φ_s(supp|Dη₀|)`, the equatorial jump lives on the
null set `{z=0}`, the flow preserves it, and for `ν > 0` the solution's `|Dη(s)|` is a.c. off
it. That needs `|Dη₀|` to carry a singular part on a null set. It does **not** for (D-C)
(real-analytic, and its `tanh` radial ramp never vanishes, so `supp|Dη₀| = ℝ⁵` and
`Φ_s(supp) = ℝ⁵`), and it does not for `ASSEMBLY §1.1`'s datum either: `min(1,|φ−π/2|/δ_m)`
cancels the `sgn(cos φ)` flip and the profile crosses the equator **linearly and
continuously**, so `|Dη₀|` is a.c. there too. The stated proof therefore applies only at
`δ_m = 0`. The conclusion is likely still true (diffusion spreads `|Dη|` off any transported
support), but the `§11` row *"the brief's pushforward inequality is FALSE for `ν > 0` —
**PROVED** (support argument)"* is a sketch outside `δ_m = 0`.
**Does it matter?** No. `§5.4(c)` supplies a correct global TV statement
(`TV(η(s)) ≤ e^{3c_G}TV(η₀)`, which I re-derived: `∫b·∇|g| = −∫2a|g|`, `|a| ≤ Γ`,
`∫Δ₅|g| = 0`), and the note's own route never uses TV.

**MINOR-4 — two hypotheses used but not stated.**
(a) P1/P2 need `V = ρ|η| → 0` and `W = ρ²|∇η| → 0` at infinity so that the sup is an interior
maximum. `§5` argues this from (D-C)'s `Θ = O(ρ^{-8})` "up to Gaussian tails"; nothing in
(H1)–(H5) gives it, and (H1) alone (bounded classical solution) certainly does not.
(b) `§6` needs `∇η ∈ L¹∩L^∞` to write `∇a = K*∇₅η`. P2 gives only `|∇η| ≲ ρ^{-2}`, which is
**not** `L¹` at infinity in `ℝ⁵` (`ρ^{-2}·ρ⁴dρ`). The composition integral itself converges
(`|x−x′|^{-4}|x′|^{-2} ∼ |x′|^{-6}`), so the *bound* is fine; the differentiation-under-the-
convolution step needs the decay hypothesis. Both hold for (D-C). Add them to (H5).
(The choice of `K*∇η` over `∇K*η` is right and should stay: `∇K` is degree `−5` in `ℝ⁵`, not
locally integrable, and would force a principal value and a Calderón–Zygmund constant. No
integration by parts is performed, so there is nothing to justify there.)

**MINOR-5 — `a(0,s) ≤ (λ/2)ML` is not exact for `s > 0`.**
`§8` reads it off `κ ≤ M_s/2`, but `a(0,s) = ∫₀^∞ κ(ρ′,s)dρ′/ρ′` and for `s > 0` the
log-support of `ω^θ` is not `[log ρ₀, log R]`: transport moves `log ρ` by at most
`c_G` at each end, and viscosity adds tails. The honest majorant is
`a(0,s) ≤ (λ/2)M(L + 2c_G) + tail`. At `L ≈ 5·10³`, `2c_G/L ≈ 5·10^{-4}`, so `L_*` moves
`≲0.05 %`. Used only in the fixed point; the theorem `Γ ≤ 2a(0,s) + C″M` is untouched.
(`u4` measures `a(0)/(ML/2) = 0.902`–`1.392`, always inside `λ`, consistent.)

**MINOR-6 — the crossing angle is computed with the wrong exponent.**
The theorem writes the collar branch as `3.9992184 e^{c_R}E₀`; `u3` computes it with
`e^{c_G} ≥ e^{c_R}` (`§12` item 2 flags the substitution as conservative, which it is for
`Ĉ_a`). But the **crossing angle** inherits it in the unsafe direction:
`σ_cross = 2R_Aλ/(2R_A e^{X}E₀ − πλ/8)` is decreasing in `X`, so
```
  X = c_G = 1.2164  ->  sin phi = 0.0583493   (phi = 3.345 deg)   <- what sec.1 displays
  X = c_R ~ c_G/2   ->  sin phi = 0.1077      (phi = 6.18 deg)
```
i.e. the region in which the axis branch is the binding one is about twice as wide in angle as
`§1` states. Harmless — `σ_* = 0` is admissible and the second branch is angle-free — but the
displayed `3.345°` is not the crossing of the theorem's own two branches.

**MINOR-7 — `{sin φ ≥ 1/2}` does not contain `N_τ` at `λ = 1`.**
`§1`: "on `{sin φ ≥ 1/2}` — which contains the tube `N_τ` at every `λ ∈ [1,3/2]`". At `λ = 1`
the tracked point sits at `φ = 30°` **exactly** (`hk2 §8(b)`), i.e. on the boundary
`sin φ = 1/2`, so a ball of radius `d` around it is not contained. Harmless (the `σ_* = 0`
column is the one that is quoted as global), but the containment statement is false as written.

**MINOR-8 — the datum note about the inner edge.**
`§1` (D-C): "the inner edge agrees with `hk2`'s `tanh` ramp exactly
(`ρΘ′(ρ₀) = 1/(2ε_r) = 2 = ρ₀/(2w₀)`)". For (D-C),
`ρΘ′(ρ₀) = ½sech²(−1)/ε_r = 0.83995`, not `2`; the value `2` is the **maximum** log-slope,
attained at `ρ = ρ₀e^{ε_r} = 1.284ρ₀`, where also `Θ = ½`. `Θ(ρ₀) = 0.1192`. So (D-C)'s edges
sit at `ρ₀e^{ε_r}` and `Re^{-ε_r}` and its effective log-width is `L − 2ε_r`, not `L`. The two
ramps share a maximum log-slope, not an inner-edge value. Cosmetic at `L ~ 5·10³`; it is a
`5 %` effect at the `L = 10` used in `§9`.

**MINOR-9 — the numerical confrontation never probes the new content.**
`§9` samples 4 radii in `[2ρ₀, R/2]` × 3 angles, all with `sin φ ≥ 1/2`, and the `max Γ` is
always at the innermost radius `2ρ₀` (which is why `max Γ − 2a(0)` is `L`-independent — the
binding geometry is pinned to `ρ₀` and does not move with `L`; the `L`-independence is
therefore a weaker confirmation than `§9` reads it as). Nothing is measured near the axis, in
the collar of an axis point, or inside either radial ramp. The note's headline novelty —
the axis-safe collar (branch B), and with it the claim that the bound is global on `ℝ⁵` —
carries **no** measurement. Not an error; a hole in the controls.

---

## 5. Answers to the specific attack points

1. **STEP 3b.** `D_t(ρ|η|)` and `D_t(ρ²|∇η|)`: re-derived independently. `Kato` is used in the
   right direction (`νΔ|η| ≥ sgn(η)νΔη = D_t|η|`). `div₅b = 2a` produces **no** zeroth-order
   term — correct, the term `−(div b)η` exists only in conservation form, and `(5.1)` follows
   from `∂_j` of the advective equation. The weight rearrangement is right and its good sign is
   an `n = 5` fact: `ρ^kΔg = Δ(ρ^kg) − (2k/ρ²)x·∇(ρ^kg) + k(k−n+2)(ρ^kg)/ρ²`, and at `n = 5`,
   `k(k−3) = −2` for both `k = 1` and `k = 2`. Hamilton's trick: at an interior maximum
   `∇V = 0` kills `b·∇V`, so `∂_sV = D_tV ≤ Γ_rad V`; the only defect is the Lipschitz
   regularity of `|·|` and of `ρ` at `0`, which is the note's own declared "modulo", plus the
   decay hypothesis of MINOR-4(a). `Γ_rad` (Lemma 4.2) is correct. **Sharpness: no — see
   MAJOR-1.**
2. **STEP 4.** `∫_{ℝ⁵}|x−x′|^{-4}|x′|^{-2}dx′ = (π⁴/2)/|x|` **confirmed to 36 digits** by an
   instrument written here (closed-form angular integral
   `J(t) = [(1+t²)log((1+t)/|1−t|) − 2t]/(2t³)`, then `|S³|∫t³J dlog t`; series continuation at
   `t → 0, ∞` to kill the cancellation). The kernel used is `K*∇₅η` with `|K| ≤ C_K|w|^{-4}` —
   the right choice, see MINOR-4(b). The bound `r|∇a| ≤ (3π²/16)e^{pc_G}𝔊₀M sin φ` holds at
   every `x`: near the axis `sin φ → 0` makes it smaller, at small `ρ` the `1/ρ` in `|∇a|`
   cancels against `r = ρ sin φ`, and there is no layer-dependence because the estimate is a
   global sup bound on `ρ²|∇η|`. No defect found.
3. **STEP 5.** The far/near split at an arbitrary field point is legitimate for the **true**
   field: Propositions 1 and 2 need only `|ω^θ| ≤ λM` on the respective region (Parseval per
   shell) and `z`-oddness — no support hypothesis — so they apply at `ρ < ρ₀` and `ρ > R/2`
   with the empty region contributing `0`, as `§12` item 4 says. `a_far(0,s)` must be taken to
   `∞` (not to `R`) for `(7.1)` to be an inequality between the same objects; `§7(a)` does
   write `∫_{2ρ}^{∞}`, correctly. `C₁ = 0` is right: `κ ≥ 0` at every radius (Lemma 2 gives
   `ω^θ cos φ ≤ 0`), hence `0 ≤ a_far(0,s) ≤ a(0,s)`, and the same fact supplies `a ≥ −Ĉ_aM`
   for Corollary 4.3. The axis-safe collar `(7.2)` is airtight: bathtub/rearrangement gives
   `∫_A|x−x′|^{-4}dx′ ≤ |S⁴|R_A`, `R_A = 1.99960922ρ`, and `C_K|S⁴| = 1`, so
   `|a_collar| ≤ 2R_A·(2e^{c_R}E₀M/ρ)/1 = 3.99921844 e^{c_R}E₀M`. **Branch crossing: `3.345°`
   is computed at `e^{c_G}`; at the theorem's own `e^{c_R}` it is `≈6.18°` — MINOR-6.**
4. **STEP 6.** `C″ = 2Ĉ_a + (3π²/16)e^{pc_G}𝔊₀ + λ` follows from `(2.1)`, `(6.1)`, `(7.3)`,
   (H3) — correct. `p = 1 + 2c_R/c_G ≤ 3` is correct (`Γ_rad ≤ ‖sym∇₅b‖ ≤ ‖∇₅b‖ = Γ`). The
   fixed point `c_G = (λ + C″/L)c` is consistent and I reproduced every `C″(L)`. `λ_max`
   (`= Λ_max(sym∇₅b)`) enters **only** through `Γ_rad`, i.e. through `c_R`, i.e. through `p`
   and through the collar's `e^{c_R}`; it never enters `Γ̄` or `C″` directly. The `λ ≤ 3/2`
   ceiling enters through (H3) in three places: the `λ` term of `C″`, the `λ`-scaling of
   `Ĉ_a`, and `a(0,s) ≤ (λ/2)ML`. See MINOR-1 for the `L_*` definition and MINOR-5 for the
   `a(0,s)` bound.
5. **STEP 7.** `u4` re-run from a copy: **bit-identical**. `457.5×` confirmed; its
   decomposition is arithmetic on an a-posteriori residual. `Γ − 2a(0)` is `L`-independent to
   five digits, but see MINOR-9 for why that is weaker evidence than `§9` reads it as.
6. **By-products.** `ε_δ` floor `4.3741737e-3`: **confirmed** (7 profiles, 12 digits).
   `‖η₀‖_∞`: **wrong for (D-C)** — MINOR-2.
7. **Transported-TV at `ν > 0`.** The claim is probably right but is proved only for
   `δ_m = 0` — MINOR-3. It does not matter: the note's route does not use TV.
8. **`§11` PROVED rows that are sketches.** Three: the sharpness row (MAJOR-1), the `‖η₀‖_∞`
   row (MINOR-2), the `§5.4(b)` row (MINOR-3). The `§8` row inherits MINOR-1 and MINOR-5. The
   `§5` P1/P2 rows are honestly labelled `PROVED-modulo`. Everything else in `§11` I could
   check, I confirmed.

---

## 6. Corrected statement of the theorem

Unchanged in structure. Two hypotheses added, one constant qualified:

> **THEOREM (Γ-off), corrected.** Let `η` be a bounded classical solution of `D_tη = νΔ₅η` on
> `ℝ⁵×[0,τ]` for the 5-D lift of axisymmetric no-swirl NS, with `z`-odd datum `η₀` satisfying
> `E₀ := sup|x||η₀|/M < ∞`, `𝔊₀ := sup|x|²|∇η₀|/M < ∞`, **and, for every `s ∈ [0,τ]`,
> `ρ|η(·,s)| → 0` and `ρ²|∇η(·,s)| → 0` as `ρ → ∞`, with `∇η(·,s) ∈ L¹(ℝ⁵)`**. Suppose on
> `[0,s]` `‖ω^θ‖_∞ ≤ λM` (`λ ≤ 3/2`) and `Γ(σ) := ‖∇u(·,σ)‖_{L^∞(ℝ⁵)} ≤ Γ̄`. Put
> `c_G := ∫₀^sΓ`, `c_R := ∫₀^sΓ_rad`, `p := 1 + 2c_R/c_G ≤ 3`. Then for **every** `x ∈ ℝ⁵`
> ```
>   ‖∇u(x,s)‖_op  ≤  2 a(0,s)  +  C″ M ,          (C₁ = 0)
>   C″ = 2 Ĉ_a  +  (3π²/16) e^{p c_G} 𝔊₀  +  λ ,
>   Ĉ_a = 0.30675280 λ + min( 3.99921844 λ/σ_* + πλ/8 ,  3.99921844 e^{c_R} E₀ ) ,
> ```
> `σ_* = 0` admissible. `3π²/16 = 1.8505508252042546 = C_K·π⁴/2`.
>
> For datum (D-C): `E₀ = 7.66060196`, `𝔊₀ = 20.9197122`. For `ASSEMBLY §1.1`'s angular family:
> `E₀ = 7.66129753`, **`𝔊₀ = 58.6948961`**, and every `C″`/`L_*` below must be recomputed.
>
> The window fixed point `c_G = (λ + C″/L)c`, `c = 2log(3/2)`, closes for
> `L ≥ L_* = 4904.7` (`σ_* = 1/2`, `λ = 3/2`) and `L ≥ 5273.0` (`σ_* = 0`), with
> `C″(10⁴) = 582.24` and `801.05`. The `a(0,s) ≤ (λ/2)ML` used there should read
> `a(0,s) ≤ (λ/2)M(L + 2c_G) + O(tail)`.
>
> **Status:** PROVED modulo (i) the standard approximate-maximum argument for `V = ρ|η|` and
> `W = ρ²|∇η|`, and (ii) the contradiction hypothesis `‖ω‖ ≤ λM`, which is (S3)'s own.
> P1 and P2 are **not** sharp on the reference strain for (D-C); they are lossy by `λ³`
> and `4.18`.

---

## 7. Remaining gaps, in order of what they are worth

1. **The `≈14×` in the weights (MAJOR-1).** P1/P2 are lossy by `λ³ · 4.18` on the reference
   strain for (D-C) because `E₀` sits on the axis and `𝔊₀` in the taper transition. A
   direction-aware version — carrying `sup |T_λ x̂|·ρ|η|` rather than `sup ρ|η|`, or splitting
   the weight by polar angle — would recover most of it. This is the cheapest available
   reduction of `C″`.
2. **The `8.39×` pointwise-vs-TV factor** (`§5.4`, `|S⁴|𝔊₀/𝒥 = 8.39`) is the note's own
   identified headroom and stands.
3. **The datum (MAJOR-2).** Either amend `ASSEMBLY §1.1` to `C¹` angular profiles, or
   recompute the whole `§8` table at `𝔊₀ = 58.695`.
4. **The approximate-maximum argument** (the note's declared modulo) plus the two decay
   hypotheses of MINOR-4. All routine; none written.
5. **`a(0,s)` at `s > 0`** (MINOR-5): a clean statement of the log-support spreading would
   close the last unstated step in the fixed point.
6. **No measurement near the axis** (MINOR-9): the cheapest missing control is `u4` re-run at
   `φ = 1°, 3.3°, 6.2°` and at `ρ ∈ (ρ₀, 2ρ₀)`, which is where branch B is binding.

**FL-000 stands.** Nothing here bears on it; BLOCK 3's hypothesis (i) is closed as `u2` claims,
with the corrections above.

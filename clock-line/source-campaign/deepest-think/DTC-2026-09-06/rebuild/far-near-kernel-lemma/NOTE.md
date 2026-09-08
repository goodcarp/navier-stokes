# GATE H-c — the far/near split of the 5D-lift strain kernel, with explicit constants
Seat `far-near-kernel-lemma`, DTC-2026-09-06.  Pre-registration in `PREREG.md` (written before any run).
Every number below comes from a script in this folder; SHA256 in `SHA256SUMS`.
Numerics falsify, they never prove: the algebra below is exact (sympy), the numerics are refuters.

---

## 0. Setting and the one representation everything follows from

5D lift of axisymmetric no-swirl NS.  `x = (y,z)` with `y in R^4`, `r = |y|`, `rho = |x|`, `t = z/rho = cos phi`
(`phi` measured from the `+z` axis).  `eta = omega^theta/r`; `-(d_rr + (3/r) d_r + d_zz) psi1 = eta` is
exactly `-Laplacian_5 psi1 = eta`; `a = u^r/r = -d_z psi1`.  With `G_5(x) = 1/(8 pi^2 |x|^3)`,

    a(x) = INT_{R^5} K(x - x') eta(x') dx',      K(w) = 3 w_z / (8 pi^2 |w|^5),

homogeneous of degree `-4`.  (This reproduces the campaign's own instrument exactly: writing
`dx' = r'^3 dr' dOmega_3 dz'`, `|S^3| = 2 pi^2`, `dOmega_3 = 4 pi sin^2(th) dth`, gives the
`(3/(2 pi)) INT INT eta r'^3 (z-z') J dr' dz'` of `toys/g2_selfstretch_biot_savart.py`.)

**Shell representation.**  Put `w := omega^theta`, `s(sigma) := sin` of the polar angle of `sigma in S^4`, and

    Phi[w](xi) := INT_{S^4} K(xi - sigma) w(sigma)/s(sigma) dOmega_4(sigma).

Because `K` is homogeneous of degree `-4` and `eta dx' = w rho'^3 sin^2(phi') drho' dphi' dOmega_3`,

    **(*)   a(x) = INT_{rho0}^{R} Phi[ w(rho' .) ](x/rho')  drho'/rho'.**

The strain is a **log-integral over shells** of a scale-invariant shell functional evaluated at the
rescaled evaluation point.  "e-fold count" is literally the measure `drho'/rho'` in (*).

## 1. The exact modal identities (script `s1_identities.py`, sympy, exact rationals)

Expand the shell density in 5D zonal harmonics, `w(sigma)/s(sigma) = SUM_l g_l C_l^{3/2}(t)`, weight `(1-t^2)`,
`N_l := INT_{-1}^1 (C_l^{3/2})^2 (1-t^2) dt = (l+1)(l+2)/(l+3/2)` (verified exactly, `l <= 8`).
Matching across the shell gives the interior/exterior coefficient `alpha_l = g_l/(2l+3)` (exact).
The two derivative identities, verified **exactly for every `l <= 14`**:

    I1   d_z [ rho^l      C_l^{3/2}(z/rho) ] =  (l+2) rho^{l-1}  C_{l-1}^{3/2}(z/rho)      (l >= 1)
    I2   d_z [ rho^{-l-3} C_l^{3/2}(z/rho) ] = -(l+1) rho^{-l-4} C_{l+1}^{3/2}(z/rho)      (l >= 0)

> **KILL K2 FIRED, and was repaired.**  PREREG registered I1 with the constant `(2l+1)`.  That is FALSE
> from `l = 2` on (`d_z[rho^2 C_2^{3/2}] = 12 z`, not `15 z`).  The error was matching the coefficient of
> `t^l` in `C_l^{3/2}` instead of the coefficient of `z^l` in `rho^l C_l^{3/2}(z/rho)`; the correct constant
> is `l C_l^{3/2}(1)/C_{l-1}^{3/2}(1) = l+2`.  Both forms are tested side by side in `s1`.  Every constant
> downstream is computed from the repaired `(l+2)`.  I2 (registered `-(l+1)`) was correct as registered.

Hence, with `xi = x/rho'`:

    interior (|xi| < 1):   Phi[w](xi) = - SUM_{l>=1} ((l+2) g_l/(2l+3)) |xi|^{l-1} C_{l-1}^{3/2}(t)
    exterior (|xi| > 1):   Phi[w](xi) = + SUM_{l>=0} ((l+1) g_l/(2l+3)) |xi|^{-(l+4)} C_{l+1}^{3/2}(t)

**Consequence A (the uniform-strain fact, now a one-line exact statement).**  The `l = 1` interior term
carries `|xi|^0 C_0^{3/2} = 1`: it is **exactly constant inside the shell**.  That constant is

    Phi[w](0) = -(3/5) g_1 = -(3/4) INT_0^pi w(phi) cos(phi) sin^2(phi) dphi,

and `|Phi[w](0)| <= (3/4) M INT_0^pi |cos phi| sin^2 phi dphi = M/2`, **with equality iff `w = -M sgn(cos phi)`**
(the near-cap bang-bang shell).  So `kappa_0 = 1/2` per e-fold is not just the axis value: it is the exact,
sharp, `x`-independent `l=1` interior mode, and the "uniform strain `(r,z) -> (r/2, -z)`" of the brief is
precisely this `l=1` solid harmonic.

**Consequence B (the multipole power).**  The leading exterior term is `l = 0`: `|xi|^{-4}`.  Its coefficient
`g_0` is proportional to the shell mean of `omega^theta`, so **`g_0 = 0` whenever `omega^theta` is odd in `z`** —
which is exactly the near-cap stack.  For a z-odd shell the leading exterior term is `l = 1`: `|xi|^{-5}`.

## 2. GATE H-c, PROPOSITION 1 (FAR).  `s2_constants.py`, checked by `s5_split_check.py`

Let `|omega^theta| <= M` on `{|x'| >= 2|x|}` and put `rho = |x|`,
`a_far(x) := INT_{2 rho}^{R} Phi[w(rho' .)](x/rho') drho'/rho'`.  Then

    a_far(0) = INT_{2 rho}^{R} ( -(3/5) g_1(rho') ) drho'/rho'      (the "M-weighted e-fold count", rate <= M/2)
    |a_far(x) - a_far(0)| <= C1 M,      C1 = sqrt(2) INT_0^{1/2} Q(u) du/u,
    Q(u)^2 = SUM_{l>=2} ((l+2)/(2l+3))^2 (l(l+1)/2)^2 u^{2(l-1)} / N_l .

Proof: interior expansion above; `|C_{l-1}^{3/2}(t)| <= C_{l-1}^{3/2}(1) = l(l+1)/2`; Cauchy-Schwarz in `l`
against Parseval `SUM_l g_l^2 N_l = INT_{-1}^1 w(t)^2 dt <= 2 M^2`; then `u = rho/rho'`, `drho'/rho' = du/u`,
and `Q(u) = O(u)` makes `INT_0^{1/2} Q(u) du/u` finite.  Numerically

    **C1 = 0.756945   (general bounded omega^theta)**
    **C1 = 0.291999   (omega^theta odd in z -- the near-cap stack; only odd l, remainder starts at u^2)**

Measured (exact, `s5`, bang-bang cap, `R/rho0 = 4096`, nine `(rho,phi)` points): `|a_far(x) - a_far(0)|`
lies in **`0.0098 .. 0.0141 M`** -- inside `C1 = 0.292` by a factor `~21..30`.  **KILL K4 did not fire.**

## 3. GATE H-c, PROPOSITION 2 (NEAR).  The near part is O(M) with NO log.

Split `{rho' < 2 rho}` into the **inner** part `{rho' <= rho/2}` and the **collar** `{rho/2 < rho' < 2 rho}`.

**(a) Inner shells -- exterior multipole, geometric, no log.**  From the exterior expansion, with `v = rho'/rho`,

    |a_inner(x)| <= sqrt(2) M INT_0^{1/2} S(v) dv/v,
    S(v)^2 = SUM_{l>=0} ((l+1)/(2l+3))^2 ((l+2)(l+3)/2)^2 v^{2(l+4)} / N_l ,   S(v) = O(v^4)  (O(v^5) if z-odd)

    **C2in = 0.026093 (general)      C2in = 0.014754 (z-odd)**

The point of principle: because the 5D kernel `K` is homogeneous of degree `-4`, an inner shell's field at
`rho >> rho'` falls like `(rho'/rho)^4` -- so the inner shells' total is a CONVERGENT GEOMETRIC SUM, not a
log.  All the log lives in the `l=1` interior mode of the OUTER shells.  Measured `|a_inner|` (exact, `s5`):
`3.2e-3 .. 1.2e-2 M`, inside `C2in`.  **KILL K5 did not fire.**

**(b) The two-octave collar -- a bounded, `phi`-dependent constant.**  Bounding absolutely in `R^5`,
`|a_collar| <= (3/(8 pi^2)) INT_{A} |omega^theta| /(r' |x-x'|^4) dx'`, `A = {rho/2 < |x'| < 2 rho}`.  Split `A` at
`r' = r/2`.  On `{r' >= r/2}`: `1/r' <= 2/r` and the bathtub/rearrangement bound gives
`INT_A |x-x'|^{-4} dx' <= |S^4| R_A` with `R_A = (2^5 - 2^{-5})^{1/5} rho = 1.999609 rho`.  On `{r' < r/2}`:
`|x-x'|^2 >= (z-z')^2 + r^2/4`, and `(1/r') dx' = 2 pi^2 r'^2 dr' dz'`, so
`INT_0^{r/2} r'^2 dr' INT_R dz'/((z-z')^2 + r^2/4)^2 = (r^3/24)(4 pi/r^3) = pi/6`, a pure number.  Hence

    **|a_collar(x)| <= ( 3.999218 / s  +  pi/8 ) M,     s = sin(phi) at the evaluation point.**

Measured (exact, `s5`): `0.058 M` (`phi=10 deg`), `0.167 M` (`45 deg`), `0.609 M` (`90 deg`) -- the bound holds
everywhere, crudely (`x 7 .. x 400`).  The `1/s` is NOT removable: see the additive find in section 6.

**(c) The near part carries no log (KILL K5').**  Exact evaluation of `collar + inner` at `phi = 45 deg` over
six octaves `rho/rho0 = 2 .. 64`:  `0.16650, 0.16196, 0.16181, 0.161810, 0.161809, 0.161809`.
Fitted `d(near)/d log(rho/rho0) = -0.00099` (threshold `0.02`).  **KILL K5' did not fire.**

## 4. GATE H-c, THE SPLIT (assembled)

For `|omega^theta| <= M` supported in `rho0 < |x'| < R`, at any `x` with `rho0 <= |x| = rho <= R/2`:

    | a(x)  -  (M/2) * L_out(x) |  <=  C_tot(phi) M,     L_out(x) := log( R / (2 rho) ),
    C_tot(phi) = C1 + 4/sin(phi) + pi/8 + C2in
               = 5.175 (general) / 4.699 (z-odd)  at phi = 90 deg
               = 6.831 (general) / 6.355 (z-odd)  at phi = 45 deg
               = 24.21 (general) / 23.73 (z-odd)  at phi = 10 deg

with `(M/2) L_out` replaced by the exact `INT_{2rho}^{R} (-(3/5) g_1(rho')) drho'/rho'` when the outer shells
are not at the cap (`|that| <= (M/2) L_out`, equality only for the bang-bang cap).
**Every constant is explicit and every one was checked against an exact evaluation; none was violated.**

## 5. Verification against the validated instrument (KILL K1)

`s3_series_vs_instrument.py` builds `a(rho,phi)` for the bang-bang cap by term-by-term e-fold integration
of (*) -- a route that shares no code with the campaign's quadrature -- and compares with
`toys/g2_selfstretch_biot_savart.py` (`rho0 = 1`, `R = 4096`):

| `rho` | `phi/pi` | series | instrument | diff |
|---|---|---|---|---|
| 8 | 1/6 | 2.85169537 | 2.85170670 | -1.13e-05 |
| 128 | 0.45 | 1.86901346 | 1.86901538 | -1.93e-06 |
| 512 | 1/3 | 0.94693501 | 0.94693870 | -3.69e-06 |

Instrument axis control: `a(0,0) = 4.15888308` vs exact `(M/2) log 4096 = 4.15888308`, rel `2.1e-16`.
`g_l` from quadrature vs the exact sympy rationals of `s1`: worst `3.9e-14`.  **KILL K1 did not fire.**

## 6. The plateau constant +0.216773 -- REPRODUCED, and then corrected in scope

The shell representation gives a **closed form** for the inner-edge offset of the bang-bang cap:

    c_edge(phi) := lim_{R/rho0 -> inf} [ a(rho0+, phi) - (M/2) log(R/rho0) ]
                 = - SUM_{l odd >= 3} ( (l+2) g_l / ((2l+3)(l-1)) ) C_{l-1}^{3/2}(cos phi).

    c_edge(10 deg) = **+0.2167733**    vs the sharp-clock refuter's measured **+0.216773**   (diff 2.7e-7)

and at the refuter's own finite geometry (`rho = rho0(1+1e-7)`, `phi = 10 deg`) the series reproduces its
whole convergence table -- `R = 64: 2.296244` (refuter `2.296241/2.296244`), `256: 2.989363` (`2.989364`),
`1024: 3.682509` (`3.682508`), `4096: 4.375656` (`4.375657`).  Two independent instruments (their elliptic-
integral ring kernel, my 5D zonal series) now agree to `~1e-6`.  **KILL K6 did not fire: the constant stands.**

**But it is not a constant of the geometry -- it is the value at one angle.**  The exact profile:

| `phi` (deg) | 1 | 2 | 3 | 5 | 7.5 | 10 | 15 | 20 | 30 | 45 | 60 | 75 | 90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `c_edge/M` | +0.787 | +0.614 | +0.513 | +0.386 | +0.286 | **+0.2168** | +0.122 | +0.059 | -0.017 | -0.062 | -0.050 | +0.008 | +0.116 |

It changes sign, and it **diverges logarithmically as `phi -> 0`**.  `s6_plateau_profile.py` fits
`c_edge ~ A log(1/phi) + B` over `phi <= 5 deg` and settles what the divergence is:

| datum | log rate `kappa = -(3/5) g_1` | `A` | `B` | `c_edge` at `phi = 0.25 deg` |
|---|---|---|---|---|
| bang-bang `w = -M sgn z` (INADMISSIBLE) | 0.500000 | **0.25171** (`~ 1/4 = kappa_0/2`) | -0.2299 | +1.1430 (still growing) |
| taper `delta = 15 deg` (admissible) | 0.497808 | 0.00056 | +0.0990 | +0.1017 (SATURATED) |
| taper `delta = 7.5 deg` (admissible) | 0.499721 | 0.00459 | +0.2720 | +0.2939 (SATURATED) |

(The `kappa` column reproduces the sharp-clock refuter's admissibility table independently: it reports
`0.49781` at `delta = 15 deg` and `0.49972` at `delta = 7.5 deg`; I get `0.497808` and `0.499721`.)

So for ADMISSIBLE data the inner-edge offset is bounded, but its ceiling is set by the taper width:
`sup_phi c_edge ~ (1/4) log(1/delta) + O(1)` (`0.1017` at `delta = 15 deg`, `0.2939` at `7.5 deg`,
difference `0.192` vs `(1/4) log 2 = 0.173`).  Three consequences the estate should carry:
* G17's `c_self = 0.2168` is the self-strain coefficient **at `phi = 10 deg` for the bang-bang datum only**;
  on the same material it runs over `[-0.062, +inf)`, and for an admissible `delta`-tapered cap its ceiling is
  `~(1/4) log(1/delta)`.  The octave model's `c_self` is a free parameter of the datum, not a constant of the
  geometry.  (G17's verdict is a comparison of spreads and is unlikely to move; but any argument that prices
  `c_self` numerically has to price `delta` too.)
* the `1/sin(phi)` in the proved collar bound is REAL, not an artifact of a crude estimate: the exact
  `c_edge` genuinely blows up at the axis for the extremal (inadmissible) datum.
* the divergence rate `A = 1/4 = kappa_0/2` is itself a clean number and is what admissibility has to buy back.

## 7. THE LEMMA (the rebuild-time kernel input)

**Lemma (stack strain / inner-octave decay).**  Let `omega^theta` be axisymmetric, no swirl, supported in
`rho0 <= |x| <= R`, with `|omega^theta| <= M` on `{|x| >= 2 rho_in}` and `|omega^theta| <= M_in` on the inner
octave `{rho_in <= |x| <= 2 rho_in}` (`M_in >> M`).  Let `x` be a point of the `j`-th octave,
`|x| = 2^j rho_in`, `j >= 2`, at polar angle `phi`, `s = sin(phi)`.  Then

    a(x)  =  (1/2) * ( M-weighted e-fold count outside 2|x| )  +  E,
    |E|  <=  ( C1 + 4/s + pi/8 + C2in ) M   +   D * M_in * 2^{-(j-1)(p)} ,

with `p = 5` when the inner octave is odd in `z` (the near-cap stack: its shell mean of `omega^theta` vanishes,
so its 5D monopole `l=0` is absent) and `p = 4` otherwise.  For the z-odd near-cap octave the explicit
constants are `D(l=1 alone) = 12.400`, `D(all l, valid for |x| >= 4 rho_in) = 124.40`, i.e.

    | inner-octave contribution |  <=  124.40 * M_in * (rho_in/|x|)^5   =   124.40 * M_in * 2^{-5j}.

*The exact power is 5 (z-odd) / 4 (general).*  The brief's guess `2^{-3j/2}` is **corrected**: the decay is
far faster than a half-power -- a quadrupole, not a "geometric factor".

Verification (`s4_multipole_decay.py`, series vs an independent direct 5D quadrature, `phi = 40 deg`,
`rho/rho_in = 4 .. 128`): fitted decay exponents **5.0000 / 5.0000** (z-odd) and **4.0009 / 4.0009**
(`omega^theta = -M` constant, `g_0 = -3 pi/8 = -1.178097` exactly).  **KILL K3' did not fire.**

**Dynamical reading (this is the H-b input, and it strengthens H-b).**  For octave `j` to be rebuilt to a new
cap level `M_new` at the new cap's own rate, it must feel strain `~ M_new`.  The lemma says the strain it
feels is at most `(M/2) L_out + C(phi) M + 124.4 M_new 2^{-5j}`.  So an amplified inner cap can drive octave
`j` at rate `M_new` only while `2^{5j} <= 124.4`, i.e. `j <= 1`: **beyond the first octave the inner cap is
dynamically invisible** and the octave is rebuilt only at the ambient rate `(M/2) L + O(M)`.  The rebuild
time of octave `j` is therefore `>= min( 2^{5j}/(124.4 M_new), 2/(M L) )` -- not `2^{j/2}/M_new`.
The "outer log cannot be rebuilt in O(1) turnovers of the new cap" conclusion of §H.3(a) survives with a
much larger margin than the brief assumed.

## 8. Status, and what this does and does not establish

* PROVED (exact, sympy-checked, explicit constants): the shell representation (*), the interior/exterior
  modal expansions, `kappa_0 = 1/2` per e-fold as the sharp `l=1` interior mode, the far Taylor remainder
  `C1`, the inner multipole bound `C2in`, the collar bound `4/s + pi/8`, and the multipole powers `5`/`4`.
* NOT proved here: any dynamical statement.  Section 7's "dynamical reading" is arithmetic on the kernel
  bound plus the frozen-position hypothesis of G17; it is not a PDE theorem.  The strain bound is
  instantaneous and says nothing about transport of the octaves.
* The constants `4/s`, `C1`, `D` are crude by factors `7-400` against exact evaluation.  Sharpening them is
  worth doing only if a downstream argument needs a specific numerical value; the STRUCTURE (log only from
  outside; geometric from inside; `1/s` at the collar) is what the rebuild lemma consumes.
* One pre-registered claim was FALSE and was repaired: the interior derivative constant is `(l+2)`, not
  `(2l+1)`.  Everything downstream uses the repaired constant.

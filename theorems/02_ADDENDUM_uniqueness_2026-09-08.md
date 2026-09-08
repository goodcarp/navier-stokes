# Addendum to Theorem 02: discharging (H*) in the Clay class

Date: 2026-09-08. Seat: sub-seat of the HOLD UP SHIPS campaign, forced-route line.
Scope: the standing hypothesis (H*) of `02_forced_axisymmetric_on_axis_typeII.md` (below [02]),
its referee finding R1 in `02_REFEREE_REPORT.md`, and the uniqueness caveat of
`01_two_fences.md` section 3 (below [FENCES]).
Script beside this note: `02_ADDENDUM_flux_exponents.py`, results in
`02_ADDENDUM_flux_exponents.json`, 50/50 checks pass. Digests in `SHA256SUMS`.

**Result.** (H*) is discharged. For Fefferman data (4) and force (5) with `div f = 0`, every
solution in the Clay class (6), (7) coincides with the maximal strong `H^s` solution, and the
two maximal times are equal. Theorem 02 therefore holds with (H*) deleted from its hypotheses.
The chain has four links, three of them written out here and one cited; the cited one is the
Caffarelli-Kohn-Nirenberg epsilon-regularity proposition at unit scale, verified against the
same secondary source [FENCES] already relies on for Lemma B, with CKN 1982 itself not accessed.
This is a single-seat result and is a candidate, not a promotion (L-16).

Status lines, one per unit:

| Unit | Statement | Status |
|---|---|---|
| Lemma U | uniqueness of a bounded-energy smooth competitor against a strong reference | PROVED |
| Prop K | the reference solution exists, is the `H^s` solution, and satisfies (H*) | PROVED (in-estate citations) |
| Cor C1 | every Clay-class solution equals the strong solution up to `T_*` | PROVED |
| Prop D | every Clay-class solution has finite dissipation on every finite interval | PROVED |
| Prop E | every Clay-class solution is bounded on `R^3 x [0,T]`, `T < T_v` | ENCLOSED (rests on the cited CKN proposition) |
| Cor C2 | `T_v = T_*`; (H*) discharged; Theorem 02' | ENCLOSED (inherits Prop E) |

Notation follows [02] throughout: (4), (5), (6), (7) are Fefferman's conditions as quoted in
[FENCES] section 3; `(H)`, `(H*)`, `(Dec)`, `Sigma` are as in [02] sections 1.1 and 2.

---

## 1. The comparison lemma at source

Source: the unpublished forced Navier-Stokes finite-time blowup preprint, read on this machine
as the local text extraction `navier-stokes.txt` held in this session's scratchpad (below
[PAPER]); its sha256 is in section 10, which is how to identify the exact bytes the line numbers
below refer to. Lemma 10.5 at lines 6277 to 6395,
printed pages 121 to 123, in subsection 10.3 "Energy and comparison"; the surrounding
Lemma 10.4 at lines 6255 to 6275; the closing remark used in section 8 below at lines 6404 to
6427, printed page 124. Unpublished; nothing below depends on the correctness of the rest of
that preprint, only on this one lemma and its proof, which is self-contained.

### 1.1 The statement as printed

> **[PAPER] Lemma 10.5.** Fix `T < 1`. If `v, P` is a smooth solution of (1.1) at viscosity one
> on `R^3 x [0,T]`, with the force `f` of Lemma 10.3, zero initial velocity, and
> `v in L^inf([0,T]; L^2(R^3))`, then `v = u` on that interval, where `u` is the localized
> velocity of Proposition 10.1.

The reference solution `u` is not a general object. By [PAPER] Proposition 10.1 (line 6069) it
is smooth on `R^3 x [0,1)` with `supp u(.,t) ∪ supp p(.,t)` contained in **one fixed compact
set `K`** for every `t`, and it vanishes for all small `t`. The force `f` of Lemma 10.3 (line
6196) is in `C_c^inf(R^3 x (0,inf); R^3)`.

**What is assumed of the reference `u`.** Smoothness on the closed slab; fixed compact spatial
support. The proof uses that support in exactly one place (section 2.2 below). Everything else
it uses is a consequence: `||u(t)||_2`, `||u(t)||_6`, `||grad u(t)||_inf` bounded on `[0,T]`.

**What is assumed of the competitor `v`.** Smoothness on the closed slab, and
`v in L^inf([0,T];L^2)`. That is all. No decay of `v` at spatial infinity, no bound on
`grad v`, no enstrophy, no energy inequality, no integrability in time beyond `L^inf_t L^2_x`,
no relation between `v` and any weak-solution class.

**What is assumed of the pressures.** Nothing beyond smoothness and that each pair solves the
equations. The lemma's own preamble (line 6272) says so: "These hypotheses leave the growth of
spatial derivatives at infinity unrestricted. We recover the pressure gradient from the
equation before passing to the limit in the localized energy identity." `P` is not assumed to
be the Leray pressure and is not assumed to decay or to be integrable.

**What is assumed of the force.** Only that both solutions carry the same `f`. It cancels in
the difference equation and never appears again.

### 1.2 The proof, in four moves

Write `w = v - u`, `pi = P - p`, `g_ij = v_i v_j - u_i u_j = w_i w_j + w_i u_j + u_i w_j`. The
difference solves

    d_t w + (v.grad)w + (w.grad)u = Laplace w - grad pi,    div w = 0.       (10.15)

On `[0,T]`, `||w(t)||_2 <= C_T` (from `v in L^inf_t L^2` and `u in L^inf_t L^2`) and
`sum_ij ||g_ij(t)||_1 <= C_T` (each summand is a product of two `L^2` functions).

**(a) The pressure gradient is determined.** Set `pi_* = sum_ij R_i R_j g_ij`, `R_i` the Riesz
transform. Since `g in L^inf_t L^1_x`, its Fourier transform is bounded, so `pi_*` lies
uniformly in `H^{-s}` for each `s > 3/2`, and `Laplace pi_* = - sum_ij d_i d_j g_ij`. For
`a in C_c^inf(0,T)`, the conservative form of (10.15) writes `int a grad pi dt` as the sum of a
term in `H^{-2}`, a term in `L^2` and a term in `H^{-3}`, so
`H_a := int a (grad pi - grad pi_*) dt` lies in `H^{-3}`. Taking the divergence of (10.15) gives
`Laplace pi = - sum_ij d_i d_j g_ij = Laplace pi_*`, so `Laplace H_a = 0`; the Fourier transform
of `H_a` is a weighted `L^2` function supported at the origin, hence zero. So
`grad pi = grad pi_*` in the interior, **with the spatial growth of `P` unrestricted**. This is
the move that replaces a decay assumption on the pressure.

**(b) The pressure flux is bounded by the dissipation inside the ball.** Fix `phi` smooth,
`0 <= phi <= 1`, compactly supported, equal to one on the unit ball, `phi_R(x) = phi(x/R)`, and
`chi_R = phi_R^8`. Set

    E_R = int chi_R |w|^2,   A_R = ( int chi_R |grad w|^2 )^{1/2},   B_R = || phi_R^4 w ||_6.

Sobolev and the product rule give `B_R <= C(A_R + R^{-1}||w||_2)` (their (10.17)). Multiplying
`pi_*` by `phi_R^4` and commuting,

    phi_R^4 pi_* = sum R_i R_j (phi_R^4 g_ij) + sum [phi_R^4, R_i R_j] g_ij.       (10.18)

The first sum is `O(B_R + 1)` in `L^{3/2}`, using `||phi^4 w_i w_j||_{3/2} <= B_R ||w||_2` and
`||phi^4 w_i u_j||_{3/2} <= ||w||_2 ||u||_6`. **This second inequality is the only place the
reference solution's own integrability enters the pressure estimate, and it asks for `L^6`, not
for compact support.** The commutator kernel is bounded by
`K_R(z) = C|z|^{-3} min{|z|/R, 1}`, whose `L^{4/3}` norm is `C R^{-3/4}` (the radial integration
is reproduced and checked two ways in the script, symbolically and by quadrature: the exact
value of `R . ||K_R||_{4/3}^{4/3}` is `4C^{4/3}`, independent of `R`). Young's convolution
inequality then bounds the second sum in `L^{4/3}` by `C_T R^{-3/4}`. Since
`|grad chi_R| <= C R^{-1} phi_R^7 = C R^{-1} phi_R^4 . phi_R^3` and
`||phi^3 w||_3 <= B_R^{1/2}||w||_2^{1/2}`, `||phi^3 w||_4 <= B_R^{3/4}||w||_2^{1/4}`,

    | int pi w . grad chi_R | <= C_T R^{-1} (B_R + 1) [ B_R^{1/2} + R^{-3/4} B_R^{3/4} ].   (10.19)

The step from `int pi w.grad chi_R` to `int pi_* w.grad chi_R` is move (a) tested against the
compactly supported field `a(t) chi_R w`, whose divergence is `w . grad chi_R`.

**(c) The difference energy.** Pairing (10.15) with `chi_R w`, all integrations having compact
support,

    (1/2) E_R' + A_R^2 = - int chi_R (w.grad)u . w + (1/2) int |w|^2 Laplace chi_R
                         + (1/2) int |w|^2 v . grad chi_R + int pi w . grad chi_R.

**Here the compact support of the reference is used, once:** "We take `R` large enough that
`chi_R = 1` on a neighborhood of `supp u`. Then `v = w` on `supp grad chi_R`", so the transport
flux is `C R^{-1} int phi_R^6 |w|^3 <= C_T R^{-1} B_R^{3/2}`.

**(d) Grönwall.** Every power of `A_R` appearing in the flux bounds is at most `3/2 < 2`, so
Young's inequality absorbs each into `A_R^2` and leaves an inverse power of `R`. What remains is

    (1/2) E_R' + (1/2) A_R^2 <= ||grad u||_inf E_R + C_T / R,

with `E_R(0) = 0` because the two solutions share their initial datum. Grönwall gives
`E_R(t) <= C_T' / R`, and `chi_R = 1` on each fixed ball for large `R`, so `w = 0`.

The coefficient `||grad u||_inf` is the second and last demand on the reference solution.

---

## 2. Does it apply in our setting

### 2.1 The audit

The competitor side transfers with nothing to check. A Clay-class solution is `C^inf` on
`R^3 x [0,T_v)` and satisfies (7), so on every closed `[0,T]` with `T < T_v` it is smooth with
`sup_t ||v(t)||_2^2 < C`. Its pressure is only smooth, which is exactly what move (a) is built
for. The two solutions share `u_0` and `f`, so `w(0) = 0` and the force cancels.

The reference side needs, and only needs:

| Demand of the proof | Where it is used | Supplied by (H*)? |
|---|---|---|
| `u` smooth on the closed slab | throughout | yes |
| `u in L^inf_t L^2` | `||w||_2 <= C_T`, `||g||_1 <= C_T` | yes, `H^0` |
| `u in L^inf_t L^6` | `||phi^4 w_i u_j||_{3/2} <= ||w||_2 ||u||_6` in (10.18) | yes, `H^1 -> L^6` |
| `grad u in L^inf_{t,x}` | the Grönwall coefficient | yes, `H^s -> W^{1,inf}` for `s > 5/2` |
| `u in L^inf_{t,x}` | **not in the source proof**; needed once compact support is dropped | yes |
| fixed compact spatial support | move (c), transport flux, once | **no** |

So the lemma does not apply as printed. One step has to be replaced.

### 2.2 The replacement, and the extended lemma

**Lemma U (PROVED-HERE).** Let `nu > 0` and `T > 0`. Let `(u,p)` and `(v,P)` both solve (NS) on
`R^3 x [0,T]` with the same force `f` and the same initial datum `u_0`, both smooth on the
closed slab, and suppose

    sup_{[0,T]} ( ||u(t)||_2 + ||u(t)||_6 + ||u(t)||_inf + ||grad u(t)||_inf ) < inf,
    sup_{[0,T]} ||v(t)||_2 < inf.

Then `v = u` on `R^3 x [0,T]`. No decay, integrability or growth condition is imposed on `P`,
on `p`, or on `grad v`.

*Proof.* Moves (a), (b), (d) of section 1.2 are unchanged: (a) and (b) never used the support of
`u`, only `||u||_6` and the `L^1` bound on `g`, and (d) is arithmetic. Only the transport flux in
move (c) changes. Since `v = w + u` everywhere rather than only outside `supp u`,

    (1/2) int |w|^2 v . grad chi_R = (1/2) int |w|^2 w . grad chi_R + (1/2) int |w|^2 u . grad chi_R.

The first term is the source's own, bounded by `C R^{-1} B_R^{3/2} ||w||_2^{3/2}`. The second is
bounded crudely, without any cancellation:

    (1/2) | int |w|^2 u . grad chi_R | <= C R^{-1} ||u(t)||_inf int phi_R^7 |w|^2
                                       <= C R^{-1} ||u(t)||_inf ||w(t)||_2^2.

This carries no power of `A_R` at all, so nothing needs absorbing, and it contributes a term of
order `R^{-1}` to the Grönwall inhomogeneity. Collecting,

    (1/2) E_R' + (nu/2) A_R^2 <= ||grad u(t)||_inf E_R + C_T / R,

with `C_T` depending on `T`, `nu`, and the four suprema above, and not on `R`. The viscosity is
carried through unchanged: `nu` multiplies `A_R^2` on the left, and Young's inequality absorbs
`R^{-k} A_R^m` with `m < 2` into `(nu/4) A_R^2` at the cost of `C(nu) R^{-2k/(2-m)}`. Since
`E_R(0) = 0`, Grönwall gives `E_R(t) <= C_T'/R` and `R -> inf` gives `w = 0`. ∎

The script tabulates the residual exponents. With compact support the slowest-decaying residual
is `R^{-4/3}`; the replacement term above decays only like `R^{-1}`, so the extension lands
exactly on the order `C_T/R` the source itself states, and does not degrade it.

---

## 3. The reference solution is the strong `H^s` solution

**Proposition K.** Let `u_0` be smooth, divergence free and axisymmetric on `R^3` satisfying (4),
and let `f` be smooth and axisymmetric satisfying (5), with `div f = 0` (see the note below).
Then:

**(a)** `u_0` is Schwartz, so `u_0 in H^s(R^3)` for every `s >= 0`; and `f in L^1(0,inf;H^s) ∩
L^2(0,inf;H^s)` for every `s`, since taking `K` large in (5) gives
`||f(t)||_{H^s} <= C_s (1+t)^{-2}`. Also `f in L^q(R^3 x [0,T])` for every `q in [1,inf]`
([FENCES] section 2.1, taking `K = 4`).

**(b)** There is `T_* in (0,inf]` and a unique `u` with `u in C([0,T'];H^s) ∩ L^2(0,T';H^{s+1})`
for every `s >= 0` and every `T' < T_*`, solving (NS) with datum `u_0` and force `f`, maximal
among such. `T_*` does not depend on `s`.
*Source:* [FENCES] section 1.3 Step 3, which gives local existence in `H^s`, `s > 5/2`, for the
forced equation with existence time bounded below in terms of `||u(t_0)||_{H^s}` and
`||f||_{L^1(t_0,t_0+1;H^s)}`, uniqueness in `L^inf_t H^s`, and propagation of higher regularity to
every `s' > s`. That step is written out there rather than cited, after the citation it
originally carried was read at source and found not to say what was attributed to it. The
standard reference for the `alpha = 2` case is Kato's `H^s` theory; I did not read it at source
in this pass and do not lean on it.

**(c)** `T_*` is characterised by the velocity: if `T_* < inf` then
`sup_{t < T_*} ||u(t)||_{L^inf} = inf`.
*Source:* [FENCES] Lemma A at `alpha = 2`, whose hypothesis (A1) holds for every `s` by (a) and
whose solution class is exactly the class in (b). Lemma A says a finite `L^inf` bound forces a
finite `H^s` bound and hence continuation. This is also why the maximal time in (b) is the same
for every `s`.

**(d)** `u` is axisymmetric; and if `u_0` and `f` are swirl free then so is `u`. Rotations about
`e_z` preserve the equation, the class in (b), the datum and the force, so uniqueness gives
axisymmetry. The reflection `x -> (x_1, -x_2, x_3)` is orthogonal, preserves the equation, and
sends `(u_r, u_theta, u_z)` to `(u_r, -u_theta, u_z)`; if the datum and force are swirl free it
fixes both, so uniqueness gives `u_theta = -u_theta = 0`.

**(e)** On every closed `[0,T']` with `T' < T_*`, `u in L^inf([0,T']; L^2 ∩ L^6 ∩ L^inf)` and
`grad u in L^inf(R^3 x [0,T'])`, by Sobolev embedding and continuity in time on a compact
interval. In particular `u` satisfies every demand of Lemma U as a reference solution, and
satisfies (H*) as literally written in [02] section 1.1, namely `u in C([0,T');H^s)` for every
`s >= 0`. **No escape of the `H^s` norm to spatial infinity is possible on `[0,T']`, `T' < T_*`,
because membership in `C([0,T'];H^s)` is how `T_*` is defined.** The escape question, which is a
real one, is not about the reference solution at all; it is about whether a competitor can
outlive `T_*`, and that is section 6.

**(f)** `u in L^inf(0,T_*;L^2) ∩ L^2(0,T_*;H^1)` **up to `T_*`**, not merely on closed
subintervals: the energy identity of [FENCES] Lemma A Step 1 gives
`||u(t)||_2 <= ||u_0||_2 + int_0^inf ||f||_2 =: M_2` for all `t < T_*`, and then
`nu int_0^{T'} ||grad u||_2^2 dt <= (1/2)||u_0||_2^2 + M_2 int_0^inf ||f||_2 dt` uniformly in
`T' < T_*`, so monotone convergence gives the integral up to `T_*`. This is the exact energy
consequence part (i) of Theorem 02 needs.

*Note on `div f`.* [02] section 3.1 reduces a general force to a divergence-free one by Leray
projection, absorbing the gradient part into the pressure; the projection preserves (5) up to
constants and preserves (A1). Everything below is stated for `div f = 0`; a general `f` is
handled by that reduction, which changes the pressure and not the velocity.

---

## 4. Uniqueness in the Clay class

**Corollary C1 (PROVED).** Let `u_0`, `f` be as in Proposition K, `u` the maximal strong solution,
`T_*` its maximal time. Let `(v,P)` be any pair with `v, P in C^inf(R^3 x [0,T_v))` solving (NS)
with the same datum and force and satisfying `sup_{t<T_v} int |v(x,t)|^2 dx < C`, and let `T_v` be
maximal in the sense of [02] section 1 (no extension to a `C^inf` solution with (7) on a longer
closed slab). Then:

**(a)** `v = u` on `R^3 x [0, min(T_*, T_v))`.
Apply Lemma U on `[0,T]` for each `T < min(T_*,T_v)`, with reference `u` by Proposition K(e) and
competitor `v`; then exhaust. No uniformity as `T` rises is needed, since the conclusion is
pointwise equality on an increasing union.

**(b)** `T_v >= T_*`. If `T_v < T_*` then `v = u` on `[0,T_v)` by (a), and `u` restricted to
`[0,T']` for `T_v < T' < T_*` is `C^inf` on the closed slab with (7), so `v` extends and `T_v` was
not maximal.

**(c)** If `T_* = inf` then `v = u` on `R^3 x [0,inf)`: the Clay-class solution is unique and
global.

**(d)** If `T_* < inf` and `T_v > T_*` then `v` is `C^inf` on the closed slab `R^3 x [0,T_*]`,
satisfies (7) there, and, since `v = u` on `[0,T_*)` and Proposition K(c) gives
`sup_{t<T_*}||u(t)||_inf = inf`, satisfies

    sup_{ R^3 x [0,T_*] } |v| = inf   with `v` smooth on that closed slab.        (Escape)

So the only way a Clay-class solution can outlive the strong solution is (Escape): a smooth,
bounded-energy velocity field whose supremum over a compact time slab is infinite, necessarily
through spatial infinity. That is precisely the failure of the zeroth-order half of hypothesis
(H) for the competitor, and it is the whole of what remains of (H*). Sections 5 and 6 close it.

Two things are worth recording about what C1 already does. First, it is stated for arbitrary
data satisfying (4), (5); **axisymmetry is nowhere used**. Second, it removes the reason
[FENCES] section 3 gave for the caveat: the energy proof of uniqueness there needed the
difference in `L^2_t H^1`, needed `u in L^inf` to kill `int_{|x|=R} u_n |w|^2`, and needed a
pressure normalisation to kill `int_{|x|=R}(p - p~) w_n`. Lemma U needs none of the three. The
first is replaced by absorbing the flux into the dissipation inside the ball rather than
controlling it by a global gradient bound; the second by the crude `||u||_inf` bound of section
2.2 (the reference is in `L^inf`, the competitor need not be); the third by move (a), which
determines `grad pi` from the equation instead of normalising `P`.

---

## 5. Every Clay-class solution has finite dissipation

**Proposition D (PROVED-HERE).** Let `nu > 0`, let `f` be smooth with `div f = 0` and
`f in L^1(0,T;L^2)`, and let `(v,P)` be `C^inf` on `R^3 x [0,T]` solving (NS) with
`M_2 := sup_{[0,T]} ||v(t)||_2 < inf`. Then

    nu int_0^T ||grad v(t)||_2^2 dt <= (1/2)||v(0)||_2^2 + M_2 int_0^T ||f||_2 dt + C(nu, M_2),

and the global energy identity holds on `[0,T]`. Consequently `v in L^inf(0,T;L^2) ∩ L^2(0,T;H^1)`,
`v in L^{10/3}(R^3 x (0,T))`, and the pressure normalised as `P_* = sum R_i R_j (v_i v_j)`, which
differs from `P` by a function of `t` alone, lies in `L^{5/3}(R^3 x (0,T))`.

*Proof.* This is the machinery of section 1.2 run on a single solution, `u = 0`, `g_ij = v_i v_j`.

Taking the divergence of the equation with `div f = 0` gives `Laplace P = - d_i d_j (v_i v_j)`.
Since `v(t) in L^2`, `v ⊗ v in L^inf_t L^1_x`, so `P_* := sum R_i R_j (v_i v_j)` is defined,
uniformly in `H^{-s}` for `s > 3/2`, and satisfies the same Poisson equation. Move (a) applies
verbatim, with `- d_t v - div(v ⊗ v) + nu Laplace v + f` in place of the difference equation:
`int a grad P dt` lies in `H^{-3}`, `Laplace(int a (grad P - grad P_*) dt) = 0`, and the
distribution vanishes. So `grad P = grad P_*` and `P = P_* + c(t)`.

Pair the equation with `chi_R v`, `chi_R = phi_R^8` as before, and write
`E_R = int chi_R |v|^2`, `A_R = (int chi_R |grad v|^2)^{1/2}`, `B_R = ||phi_R^4 v||_6`:

    (1/2) E_R' + nu A_R^2 = (nu/2) int |v|^2 Laplace chi_R + (1/2) int |v|^2 v . grad chi_R
                            + int P_* v . grad chi_R + int chi_R f . v.

The additive `c(t)` drops out of the pressure flux because `int v . grad chi_R = - int chi_R div v = 0`.
Move (b) applies with `w` replaced by `v` and with the `u_j` cross terms absent, so the
right-hand side is bounded by

    C nu R^{-2} M_2^2 + C R^{-1} B_R^{3/2} M_2^{3/2}
      + C R^{-1} [ M_2^{3/2} B_R^{3/2} + R^{-3/4} M_2^{9/4} B_R^{3/4} ] + ||f(t)||_2 M_2.

Every power of `B_R`, hence of `A_R`, is at most `3/2 < 2`, so Young absorbs each into
`(nu/2) A_R^2` with residuals `R^{-4}`, `R^{-14/5}`, `R^{-2}` (the script tabulates them; all are
strictly negative powers). Integrating on `[0,T]` and using `E_R(0) <= ||v(0)||_2^2` gives a bound
on `int_0^T A_R^2 dt` uniform in `R >= 1`; monotone convergence as `chi_R -> 1` gives the stated
inequality. The same estimate shows every flux tends to zero in `L^1(0,T)`, so the identity, not
merely the inequality, survives the limit.

The interpolation `||v||_{10/3} <= C ||v||_2^{2/5} ||grad v||_2^{3/5}` then gives
`int_0^T ||v||_{10/3}^{10/3} dt <= C M_2^{4/3} int_0^T ||grad v||_2^2 dt < inf`, so
`v ⊗ v in L^{5/3}(R^3 x (0,T))` and, Riesz transforms being bounded on `L^{5/3}`,
`P_* in L^{5/3}(R^3 x (0,T))`. ∎

Two remarks. First, Proposition D says the Clay class (6), (7) sits inside the Leray energy class
automatically, which is the thing [FENCES] section 3 said "(6), (7) do not give for free". It is
free after all, for smooth solutions, once the pressure is recovered from the equation rather
than assumed to decay. Second, it gives a second and largely decorrelated route to Corollary C1:
a Clay-class solution is now a Leray-Hopf solution satisfying the energy identity, and the strong
solution lies in `L^inf_t L^inf_x`, so Prodi-Serrin weak-strong uniqueness applies. I record that
route but do not lean on it; Prodi 1959 and Serrin 1962 are marked "not verified directly" in the
[FENCES] verification record, and Lemma U is the load-bearing argument.

I did not search the literature for a prior statement of Proposition D. If it is known, the
credit is not mine; the argument is written out here so it can be checked without a citation.

---

## 6. No escape: the Clay-class lifetime equals the strong lifetime

**Cited input.** Caffarelli-Kohn-Nirenberg epsilon-regularity at one scale, in the form of
Proposition 1.5 of the Ulm seminar notes on CKN theory (`uni-ulm.de`, `skipper/SeminarCKN207.pdf`,
fetched and text-extracted 2026-09-08), attributed there to Proposition 1 of [CKN82]:

> **Proposition 1.5.** There are absolute constants `eps_1, C_1 > 0` and a constant `eps_2(q) > 0`
> with the following properties. If `(u,p)` is a suitable weak solution of the NSE on `Q_1(0,0)`
> with force `f in L^q` for some `q > 5/2` and
>
>     ∫∫_{Q_1} (|u|^3 + |u||p|) dx dt + ∫_{-1}^0 ( ∫_{B_1} |p| dx )^{5/4} dt <= eps_1
>     and  ∫∫_{Q_1} |f|^q dx dt <= eps_2,
>
> then `u in L^inf(Q_{1/2}(0,0))` with `||u||_{L^inf(Q_{1/2}(0,0))} <= C_1`.

The same notes state (their sentence following the proposition) that it may be shifted and
rescaled to any `Q_r(x,t)`, and their Definition 1.1 fixes "suitable weak solution" as:
`f in L^q`, `q > 5/2`, `div f = 0`; `p in L^{5/4}` (they note the modern `L^{5/3}`);
`u in L^inf(0,T;L^2) ∩ L^2(0,T;H^1)`; `- Laplace p = d_i d_j (u_i u_j)`; and the local energy
inequality with the `2(u.f)phi` term. This is the same secondary source, and the same
definition, that [FENCES] section 2.1 already relies on for Lemma B. **CKN 1982 itself was not
accessed.** That is the single external dependency of this section, and it is the reason
Proposition E and Corollary C2 are graded ENCLOSED rather than PROVED.

**Proposition E (ENCLOSED).** Let `u_0`, `f` be as in Proposition K and let `(v,P)` be Clay-class
on `[0,T_v)` as in Corollary C1. Then for every `T < T_v`, `sup_{R^3 x [0,T]} |v| < inf`.

*Proof.* Normalise the viscosity: `v~(x,t) = nu^{-1} v(x, t/nu)`, `P~ = nu^{-2}P(x,t/nu)`,
`f~ = nu^{-2}f(x,t/nu)` solves at viscosity one, preserves smoothness, preserves (7) and preserves
the hypotheses on `f`. So take `nu = 1`.

Fix `T < T_v` and `t_1 in (0, min(T, T_*))`. On `R^3 x [0,t_1]`, `v = u` by Corollary C1(a) and
`u in C([0,t_1];H^s)`, so `v` is bounded there. It remains to bound `v` on `R^3 x [t_1,T]`.

By Proposition D on `[0,T]`, `v in L^inf_t L^2 ∩ L^2_t H^1`, `v in L^{10/3}(R^3 x (0,T))`, and
`P_* in L^{5/3}(R^3 x (0,T))` with `P = P_* + c(t)`. The pair `(v, P_*)` is a suitable weak
solution on `R^3 x (0,T)` in the sense of the definition above: the integrability items are
these three plus `f in L^q` for every `q` by Proposition K(a); `- Laplace P_* = d_i d_j (v_iv_j)`
by construction; and the local energy inequality holds with equality, since `v` is smooth and
the additive `c(t)` contributes `2c(t) ∫ v . grad phi = - 2c(t) ∫ phi div v = 0`. Restriction to
any subcylinder preserves all of it.

Put `rho = min(1, sqrt(t_1)) > 0`, fixed. For every `x_0 in R^3` and every `t_0 in [t_1, T]` the cylinder
`Q_rho(x_0,t_0)` lies in `R^3 x (0,T]`. Under the parabolic rescaling that carries
`Q_rho(x_0,t_0)` to `Q_1(0,0)`, the three quantities in Proposition 1.5 pick up the fixed factors
`rho^{-2}`, `rho^{-2}`, `rho^{-13/4}`, and the force quantity picks up `rho^{3q-5}` with
`3q - 5 > 5/2 > 0` (all four exponents are computed in the script). Since `rho` is fixed, it
suffices that the unrescaled integrals over `Q_rho(x_0,t_0)` tend to zero as `|x_0| -> inf`,
uniformly in `t_0 in [t_1,T]`. They do, by absolute continuity of a finite integral: with
`S(x_0) := B_rho(x_0) x (0,T)`,

    ∫∫_{S} |v|^{10/3},   ∫∫_{S} |P_*|^{5/3},   ∫∫_{S} |f|^q     all tend to 0 as |x_0| -> inf,

each integrand being globally integrable on `R^3 x (0,T)`. Hölder on the finite-measure set
`Q_rho` converts these into the quantities the proposition asks for:
`∫∫|v|^3 <= |Q|^{1/10} (∫∫|v|^{10/3})^{9/10}`; `∫∫|v||P_*| <= ||v||_{L^{10/3}} ||P_*||_{L^{10/7}}`
with `L^{5/3} ⊂ L^{10/7}` on finite measure; and
`∫_B |P_*| <= |B|^{2/5} ||P_*||_{L^{5/3}(B)}` before raising to the power `5/4` and integrating in
time. Every exponent is checked in the script.

So there is `rho_0` with: for all `|x_0| >= rho_0` and all `t_0 in [t_1,T]`, the hypotheses of
Proposition 1.5 hold on `Q_rho(x_0,t_0)`, giving `|v| <= C_1/rho` on `Q_{rho/2}(x_0,t_0)`. That
bound is an essential supremum, and `v` is continuous, so it holds at every point of the closure,
including the top centre `(x_0,t_0)`. Hence `|v| <= C_1/rho` on `{|x| >= rho_0} x [t_1,T]`. On the
complementary set `{|x| <= rho_0} x [0,T]`, which is compact, `v` is continuous, hence bounded.
Together with the bound on `R^3 x [0,t_1]`, `v` is bounded on `R^3 x [0,T]`. ∎

**Corollary C2 (ENCLOSED).** With `u_0`, `f` as in Proposition K and `(v,P)` Clay-class as above:

**(a)** `T_v = T_*`, and `v = u` on `R^3 x [0,T_*)`.
If `T_v > T_*` then pick `T` with `T_* < T < T_v`; Proposition E bounds `v` on `R^3 x [0,T]`, so
`sup_{t<T_*}||u(t)||_inf < inf` by Corollary C1(a), contradicting Proposition K(c). With C1(b),
`T_v = T_*`.

**(b) (H*) is discharged.** Every Clay-class solution with data (4), (5) satisfies
`v in C([0,T');H^s)` for every `s >= 0` and every `T' < T_v`, and satisfies
`v in L^inf(0,T_v;L^2) ∩ L^2(0,T_v;H^1)` up to `T_v` by Proposition K(f). Hypothesis (H) holds
for it as well.

**(c) The Clay reduction.** `(u_0,f)` witnesses Fefferman's (C) if and only if `T_* < inf`.
If `T_* = inf` then `u` itself is a global Clay-class solution: `u in C([0,inf);H^s)` for every
`s` gives (6), its Leray pressure is smooth, and `||u(t)||_2 <= M_2` gives (7). If `T_* < inf`
then any global Clay-class solution would have `T_v = inf > T_*`, contradicting (a).

---

## 7. Theorem 02 restated

**Theorem 02' (Clay-class form).** Let `nu > 0`. Let `u_0` be smooth, divergence free and
axisymmetric on `R^3` satisfying (4), and let `f` be smooth and axisymmetric on `R^3 x [0,inf)`
satisfying (5). Let `(u,p)` be a solution of (NS) with this datum and force which is `C^inf` on
`R^3 x [0,T)` and satisfies (7) there, with `T < inf` its first singular time. Define

    Sigma := { x_0 in R^3 : sup_{ B_rho(x_0) x (T - rho, T) } |u| = inf  for every rho > 0 }.

Then:

**(i)** `Sigma ⊂ {r = 0}`. If moreover `sup_{t<T} |u(x,t)| -> 0` as `|x| -> inf`, then
`Sigma != ∅`.

**(ii-a)** There is no `C_* < inf` with `|u(x,t)| <= C_*(r^2 + (T-t))^{-1/2}` on `R^3 x (0,T)`.

**(ii-b)** If in addition (Dec) holds, there is no `C_* < inf` with
`|u(x,t)| <= C_*(T-t)^{-1/2}` on `R^3 x (0,T)`; equivalently the singularity is Type II,
`limsup_{t->T-} sqrt(T-t) sup_x |u(x,t)| = inf`.

**(iii)** If `u_0` and `f` are swirl free then no such `T` exists: the solution is global and
smooth, and it is the only Clay-class solution with that datum and force.

The hypotheses are (4), (5), axisymmetry, membership in the Clay class, and, for (ii-b) only,
(Dec). **(H*) and (H) are gone.** The `T` of the statement is the strong solution's maximal time
`T_*` by Corollary C2(a), so the two readings of "first singular time" agree.

*Proof.* By Proposition K there is a maximal strong solution `u~` with maximal time `T_*`, and it
is axisymmetric. By Corollary C2(a), `u = u~` on `[0,T)` and `T = T_*`. So `u` satisfies (H*) and
(H), and Theorem 02 of [02] applies verbatim, giving (i), (ii-a), (ii-b). For (iii), Proposition
K(d) makes `u~` swirl free, so [02] section 7.2 gives `T_* = inf`, and Corollary C1(c) gives
uniqueness in the Clay class. ∎

### 7.1 Which steps still need the strong class of the reference

All of them. Nothing in [02] is rewritten; what changes is that the class is now supplied rather
than assumed. The dependency, step by step, using [02]'s own labels:

| Step | What it needs | Supplied by |
|---|---|---|
| S3 ([FENCES] Lemma A, forced) | `u in C([0,T);H^s) ∩ L^2_loc H^{s+1}`, (A1) | Prop K(b), K(a) |
| part (i) (Lemma B step (ii)) | `u in L^inf_t L^2 ∩ L^2_t H^1` and `int_0^T ||grad u||_2^2 < inf` up to `T` | Prop K(f) |
| part (i) (CKN suitability) | Leray pressure in `L^{5/3}_loc`, `f in L^q`, `q > 5/2` | Prop K(a), K(b) |
| S19.1, S19.3, S20 | `||omega(t)||_2` and `||u(t)||_{H^s}` finite for the Grönwall arguments | Prop K(b) |
| S16, part (iii) closure | `H^s` continuation from `||u||_{L^6}` and `||f||_{H^s}` | Prop K(b) |
| identification of the Clay solution with the reference | Lemma U | sections 2.2, 4 |

The separating field of [02] section 1.1, uniformly smooth with finite `L^2` norm and infinite
`||grad u||_{L^2}`, is not affected by anything here: it shows (H) does not imply (H*) and that
remains true. What is shown here is that a Clay-class **solution** of this equation with this
data cannot be that field, which is a different statement.

---

## 8. The remark on the source's page 124, and whether we inherit it

At lines 6404 to 6427 the preprint proves its Theorem 1.1 and, in passing, writes: "The embedding
`H^3(R^3) -> L^inf(R^3)` excludes a classical `H^3` continuation through time one. In the
classical `H^3` class, the same difference energy calculation is justified by Sobolev
approximation and `H^3 -> W^{1,inf}`; Gronwall gives uniqueness on every shorter interval."

That half-sentence is doing real work, and it marks a limit. Uniqueness for merely classical
solutions with no decay class is false, not merely unproved. The standard witness: for unforced
(NS) at any viscosity, and any smooth `b` with `b(0) = 0`,

    u(x,t) = b(t),    p(x,t) = - b'(t) . x

is a smooth solution with zero initial datum, since `d_t u = b'`, `(u.grad)u = 0`,
`Laplace u = 0`, `grad p = -b'`. There is a continuum of them. Two features kill them, and they
are exactly the two hypotheses the comparison lemma keeps: the velocity is not in `L^2`, and the
pressure grows linearly in `x`. So the source's lemma is not being careless about the pressure;
the pressure is where the counterexamples live, and move (a) of section 1.2 is the answer to
them. The bounded-energy condition (7) is not a technical convenience in Fefferman's problem
statement, it is what makes the question well posed at all.

Do we inherit the remark's weakness? No, for one reason and with one caveat.

The reason: the remark is about uniqueness **within** the classical `H^3` class, which is what
that proof needed for its own maximal-lifespan claim. Lemma U is the stronger statement, and it
is the one we use: the competitor is required to be smooth with bounded energy and nothing else,
and it is the **reference** that must sit in a class where the difference energy identity is
justified. Our reference is the strong `H^s` solution for every `s`, which is strictly inside
`H^3`, so the half-sentence's own justification applies to it a fortiori.

The caveat: the asymmetry is real and permanent. Nothing here gives uniqueness between two
merely classical bounded-energy solutions when neither is known to be strong. The whole
construction is anchored at `t = 0`, where Proposition K produces a strong solution from Schwartz
data, and it propagates forward only as far as that strong solution lives. Corollary C2 closes
the forward end by an entirely different mechanism (epsilon-regularity, not energy), and that is
why section 6 is a separate argument rather than a corollary of section 4.

---

## 9. Self-refutation

One pass, distinct attack surface: not the analysis of sections 1 to 6 but the logic of the
discharge claim and the standing of its weakest citation.

**The strongest objection.** *You have not discharged (H*); you have relabelled it. Theorem 02
was a statement about a solution in the Clay class, hypothesised to be strong. Theorem 02' is a
statement about the strong solution, asserted to be the only Clay-class one. If the identification
fails anywhere, the restated theorem is about an object that is not the obstruction to Fefferman's
(C), and the referee's finding R1 stands untouched. The identification rests, at its one load-
bearing external point, on a proposition from a 1982 paper you did not read, taken from a student
seminar's lecture notes, whose statement in those notes carries a typographical inconsistency
(the constant is introduced as `eps` and then used as `eps_1`). Grade the chain by its weakest
link and this is ENCLOSED at best, and the honest headline is "(H*) reduced, not discharged".*

**Response, in three parts.**

*Part one, conceded.* The grading is right and is already applied: Proposition E and Corollary C2
are marked ENCLOSED, not PROVED, and the summary at the head of this note names the CKN citation
as the single external dependency. The `eps`/`eps_1` slip in the notes is real; it is a naming
slip, the quantifier and the two-sided structure of the statement are unambiguous, and the same
notes are the source [FENCES] section 2.1 already uses for Lemma B, so this note does not lower
the estate's standard, it inherits it. What would raise the grade to PROVED is one thing: CKN
1982 Proposition 1 read at source, or an equivalent statement in Lin 1998 or
Ladyzhenskaya-Seregin 1999 read at source, with its force hypothesis checked. That is a bounded,
cheap piece of work and it is the first thing a second seat should do.

*Part two, the reduction is not a relabelling, and is safe even if the objection lands.* Strike
section 6 entirely and what survives is still strictly stronger than the referee's position.
Corollary C1 is PROVED and independent of CKN: it gives `v = u` on `[0, min(T_*,T_v))` and
`T_v >= T_*`, so a Clay-class solution can never be **shorter** than the strong one, and the
entire residue is the single scenario (Escape) of C1(d). That residue is the zeroth-order half of
(H) for the competitor: one `L^inf` statement, about the competitor only, on a closed slab. (H*)
was an infinite family of `L^2`-based statements plus a uniqueness theorem plus a
no-escape-to-infinity statement. Trading that for one `L^inf` bound is not a relabelling, and it
is what the referee asked for minus one step. Moreover part (iii) of Theorem 02' has no residue
even without section 6: `T_* = inf` in the swirl-free case, so (Escape) cannot arise and
Corollary C1(c) closes it outright.

*Part three, where the objection would actually bite.* The place to attack is not the citation, it
is Proposition D, which is mine and is load-bearing for section 6 in a way Lemma U is not. If
Proposition D is wrong then `(v,P_*)` is not known to be a suitable weak solution, the CKN input
cannot be applied, and section 6 collapses. Three probes, run:

1. *Does the pressure recovery survive when only one solution is in play?* The difference
   equation's `g = v ⊗ v - u ⊗ u` becomes `v ⊗ v`, still `L^inf_t L^1_x` by (7) alone, and the
   `H^{-3}` bookkeeping in move (a) uses only `d_t v`, `Laplace v`, `div(v ⊗ v)` and `f`, each of
   which is available. It survives. But note what it needs: `div f = 0`, else
   `Laplace P = - d_i d_j (v_iv_j) + div f` and `P_*` is not the Riesz form. The reduction is in
   Proposition K's closing note and is [02]'s own section 3.1; it is stated, not skipped.
2. *Do the flux exponents still absorb with `u = 0`?* The `(B_R + 1)` in (10.19) came from the
   `w_i u_j` cross terms; without them the bound is `C R^{-1} B_R^{3/2}` and the `1` branch, which
   was the slowest residual at `R^{-4/3}`, disappears. Absorption gets easier, not harder. The
   script checks the single-solution table separately (`D_transport`, `D_pressure`,
   `D_commutator`, `D_laplacian`); worst residual `R^{-2}`.
3. *Is the `R -> inf` limit legitimate?* `int_0^T A_R^2 dt` is bounded uniformly in `R >= 1` and
   `chi_R` increases to `1` pointwise, so monotone convergence applies to
   `int_0^T int chi_R |grad v|^2`. No uniform-in-time control is needed and none is claimed.

A fourth probe I could not close, and record as the residual gap of this note: Proposition D is
plausible enough that it may well be known, in which case the right move is a citation rather than
a proof, and if a published version carries a hypothesis I have not noticed, my derivation is where
that hypothesis is hiding. I did not search. That is a literature debt, not a mathematical one,
and it does not change the grade.

**Two further probes on sections 1 to 4, both dissolved.**

*The Grönwall constant `C_T` and the coefficient `||grad u||_inf` blow up as `T` rises to `T_*`,
so the identification `v = u` on `[0,T_*)` is an unjustified limit.* No: Lemma U gives equality on
each closed `[0,T]`, `T < T_*`, and `[0,T_*)` is the increasing union of those. Pointwise equality
on an increasing union needs no uniformity.

*Lemma U requires the reference in `L^inf`, so it cannot be applied on `[0,T]` as `T` approaches a
time where `||u||_inf` diverges.* Correct, and not needed: `||u(t)||_inf` is finite on every closed
`[0,T]` with `T < T_*` by Proposition K(e), and only diverges in the limit `t -> T_*`, which is the
limit taken outside the lemma.

---

## 10. Verification record

| Claim | Source | Status |
|---|---|---|
| [PAPER] text extraction `navier-stokes.txt`, sha256 `423b0a6eb4d189496d2663057ae79a9f2a18267b1a489d8a7919255b3168218e` | local copy, this session's scratchpad | Digest recorded so the line numbers below are checkable against exact bytes |
| Comparison lemma statement, hypotheses on `u`, `v`, pressures, force | [PAPER] lines 6277-6279 | Read verbatim from the extracted text |
| Proof moves (a)-(d), including the `H^{-3}` argument, the commutator kernel and the Grönwall | [PAPER] lines 6280-6395 | Read in full |
| Reference solution has fixed compact spatial support | [PAPER] Proposition 10.1, line 6069 | Read verbatim |
| Force of the lemma is `C_c^inf(R^3 x (0,inf))` | [PAPER] Lemma 10.3, line 6196 | Read verbatim |
| The `H^3 -> W^{1,inf}` remark and the classical-class uniqueness sentence | [PAPER] lines 6404-6427, printed p. 124 | Read verbatim |
| `||K_R||_{4/3}^{4/3} = 4 C^{4/3} R^{-1}` | recomputed | Symbolic (sympy) and numeric (quadrature) agreement, script checks `kernel_*` |
| All flux powers of `A_R` are `< 2`; residual exponents; extension lands at `R^{-1}` | recomputed | Script, 50/50 |
| CKN epsilon-regularity Proposition 1 at one scale, with the force hypothesis and the `C_1` conclusion | Ulm seminar notes `skipper/SeminarCKN207.pdf`, Proposition 1.5, fetched and text-extracted 2026-09-08 | Verified in secondary source only; **CKN 1982 not accessed** |
| Definition of suitable weak solution, `f in L^q`, `q > 5/2`, `div f = 0`, local energy inequality | same notes, Definitions 1.1 and 2.1 | Verified in secondary source only |
| Shift and rescale of the epsilon-regularity proposition to `Q_r(x,t)` | same notes, sentence following Proposition 1.5 | Verified in secondary source only; the four scaling exponents recomputed in the script |
| Local existence in `H^s`, `s > 5/2`, forced, existence time bounded below | [FENCES] section 1.3 Step 3 | In-estate, written out there; Kato's `H^s` theory **not read at source in this pass** |
| Lemma A at `alpha = 2` and the `L^inf` characterisation of `T_*` | [FENCES] section 1.1, 1.3 | In-estate |
| Fefferman (4), (5), (6), (7), statement (C) | [FENCES] section 3, from `claymath.org navierstokes.pdf` | Quoted there verbatim; not re-fetched here |
| Prodi-Serrin weak-strong uniqueness (second route in section 5) | Prodi 1959, Serrin 1962 | **Not verified directly**; recorded, not relied on |
| Proposition D has a prior statement in the literature | none | **Not searched.** Recorded as a literature debt in section 9 |

**What a second seat should do first,** in order of payoff: read CKN 1982 Proposition 1 (or Lin
1998, or Ladyzhenskaya-Seregin 1999) at source and check the force hypothesis, which promotes
Proposition E and Corollary C2 from ENCLOSED to PROVED; then re-derive Proposition D independently,
without reading section 5, since it is the only new mathematics here that nothing else checks; then
search for a prior statement of Proposition D. Per L-10 the second seat should be a different
vendor, and per L-14 should write its own instrument rather than run this note's script.

## Second-vendor verdict (2026-09-08, added after the blind check)

A second, independent line re-derived Proposition D from scratch by a nested-ball local energy argument (and obtained the energy equality, so the additive constant in the statement is unnecessary), confirmed the modified transport-flux step of Lemma U at order R⁻¹ with the Grönwall closing on the global L² bound, and confirmed that a Clay-class solution is a suitable weak solution once finite dissipation holds (v ∈ L^{10/3}, p₀ ∈ L^{5/3}, local energy equality), so the ε-regularity escape argument stands. The one residue is the attribution of the forced ε-regularity criterion to the 1982 paper, which neither line read at source. Capture: `02_ADDENDUM_uniqueness_XVENDOR_2026-09-08.md`. Status of Theorem 02′: two-engine confirmed modulo that citation check.

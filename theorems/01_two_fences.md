# Two fences that any forced full Navier-Stokes blowup must pass

Date: 2026-09-08. Scope: 3D incompressible Navier-Stokes with a smooth external force, on
`R^3` and `T^3`. Purpose: give two checkable necessary conditions, with proofs or precise
citations, so that a claimed proof of Fefferman's (C) or (D) can be tested against them
quickly.

Throughout, `Lambda = |grad| = (-Delta)^{1/2}`, and we study

    d_t u + (u.grad)u + grad p = -nu Lambda^alpha u + f,    div u = 0,    nu > 0,   (NS_alpha)

with `alpha = 2` the full Navier-Stokes case. Write `q_c(alpha) = 3/(alpha-1)`.

**What is proved here versus cited.** Lemma A is *proved in full* below (scaling remark plus a
complete energy argument valid for every `1 < alpha <= 2`); it uses one external ingredient,
the Kato-Ponce product estimate, quoted with verified hypotheses (Grafakos-Oh on `R^3`,
Bényi-Oh-Zhao on `T^3`). Lemma B is *assembled from cited theorems* (the
Caffarelli-Kohn-Nirenberg epsilon-regularity criterion applied to the classical solution on
backward cylinders), with the rotation-invariance step proved here; its failure for
`alpha < 2` is cited (Tang-Yu). All verification notes are collected in section 5.

---

## 0. Context and one honesty note

Alpöge and Buckmaster are reported (2026-09-08) to have proved finite-time blowup for 3D
incompressible **Euler** on `R^3` with a smooth compactly supported force, axisymmetric with
swirl, data and force supported in a solid torus away from the axis, with velocity and
circulation bounded and vorticity and circulation gradient blowing up. The mechanism is the
Córdoba-Martínez-Zoroa layered cascade.

**I could not verify this preprint.** The arXiv API index is current through 2026-09-02 for
the author `Alpöge` (21 entries, most recent 2609.02606, 2026-09-02) and returns no such
paper; searches returned nothing matching. Everything I say about the Alpöge-Buckmaster
result below is therefore *as reported to me*, not verified, and is used only to identify a
class of mechanisms, not as a cited theorem.

The dissipative analogue is verified:

> **[CMZ-Z]** D. Córdoba, L. Martínez-Zoroa, F. Zheng, *Finite Time Blow-Up for the
> Hypodissipative Navier Stokes Equations with a Force in `L^1_t C_x^{1,eps} ∩ L^inf_t L^2_x`*,
> Arch. Ration. Mech. Anal. **250** (2026), art. 38, doi:10.1007/s00205-026-02198-0;
> arXiv:2407.06776 (submitted 2024-07-09).

Verified from the abstract and the open-access full text: dissipation is `|grad|^alpha` with
symbol `|xi|^alpha`; the range is `alpha in [0, alpha_0)` with `alpha_0 = (22 - 8 sqrt 7)/9`,
which is **approximately 0.0927** (not 0.118 or 0.547, both of which appear in automated
summaries of this paper; `8 sqrt 7 = 21.1660`, so `alpha_0 = 0.8340/9`). The force is rough in
time. Crucially for what follows, their solution satisfies `u in L^inf_{t,x} ∩ L^inf_t L^2_x`,
and what diverges is `int_0^t ||grad u(s)||_{L^inf} ds` as `t -> T`. The paper does not
identify `alpha = 1` as a threshold; `alpha_0` comes from parameter balancing in the
construction.

So both the Euler result and the hypodissipative result belong to the same class: **the
velocity stays bounded and only derivatives blow up.** Lemma A says that class stops at
`alpha = 1`.

---

## 1. Lemma A: bounded velocity is subcritical above order one

### 1.1 Statement

**Lemma A.** Let `1 < alpha <= 2`, `nu > 0`, and let `u` be the strong solution of (NS_alpha)
on `[0,T) x R^3` (or `[0,T) x T^3`), that is

    u in C([0,T); H^s) ∩ L^2_loc([0,T); H^{s + alpha/2}),   p the Leray pressure,

with smooth divergence-free finite-energy data `u_0 in H^s`, `s > 5/2`, and force `f`
satisfying, for some `delta > 0`,

    f in L^1(0,T+delta; L^2) ∩ L^2(0,T+delta; H^s).                              (A1)

(The solution class is what Steps 1 and 2 use, and it removes the ambiguity
`p -> p + c(t).x`, `u -> u + shift` of a merely classical solution on `R^3`. If the lemma is
to be applied to a Clay-class solution, `C^inf` with bounded energy, one first identifies it
with this strong solution on `[0,T)`; for `alpha = 2` this is the energy argument on the
difference using (A2), which needs the difference to lie in the energy class, see the
uniqueness caveat at the end of section 3. The interval `[0,T+delta)` in (A1) is what Step 3
needs to continue past `T`; under Fefferman's (5) or (9) the force is global and this is
automatic.)

If

    M := sup_{0 <= t < T} ||u(t)||_{L^inf} < inf,                                 (A2)

then `sup_{t<T} ||u(t)||_{H^s} < inf`, and consequently `u` extends as a classical solution
past `T`. In particular no singularity forms at `T`.

**Scope condition on the force.** (A1) is what the proof uses. "Smooth and bounded with all
derivatives on `[0,T] x R^3`" does **not** imply (A1) on `R^3` (a constant force is bounded
with all derivatives but is not in `L^2`). On `T^3` it does, since the torus has finite
measure. On `R^3`, Fefferman's decay condition (5) below implies (A1) for every `s`, and so
does any smooth compactly supported force. This is the honest hypothesis: (A1), not "smooth
and bounded".

### 1.2 The scaling argument: `L^inf` is critical exactly at `alpha = 1`

If `(u, p, f)` solves (NS_alpha), then for `lambda > 0` so does

    u_lambda(x,t) = lambda^{alpha-1} u(lambda x, lambda^alpha t),
    p_lambda(x,t) = lambda^{2(alpha-1)} p(lambda x, lambda^alpha t),
    f_lambda(x,t) = lambda^{2 alpha - 1} f(lambda x, lambda^alpha t).

(Direct check: every term of (NS_alpha) picks up exactly `lambda^{2 alpha - 1}`.) Then

    ||u_lambda(.,t)||_{L^q(R^3)} = lambda^{alpha - 1 - 3/q} ||u(., lambda^alpha t)||_{L^q},

so `L^q` is scaling invariant precisely when `q = q_c(alpha) = 3/(alpha - 1)`. Thus
`q_c(2) = 3` (the `L^3` endpoint of Escauriaza-Seregin-Šverák), and `q_c(alpha) -> inf` as
`alpha -> 1+`, so **`L^inf` is the critical Lebesgue space exactly at `alpha = 1`**. For
`alpha > 1` the exponent `alpha - 1 - 3/inf = alpha - 1` is strictly positive, so zooming in
(`lambda -> 0`) sends `||u_lambda||_{L^inf} = lambda^{alpha-1} ||u||_{L^inf} -> 0`: a bounded
velocity carries no scale-invariant content at a concentrating singularity. This is a
heuristic; the proof follows.

### 1.3 Proof of Lemma A

Fix `s > 5/2` (for instance `s = 3`).

**Ingredient (cited, hypotheses verified).** Kato-Ponce / fractional Leibniz rule in the
`L^inf` endpoint form. From L. Grafakos and S. Oh, *The Kato-Ponce inequality*, Comm. Partial
Differential Equations **39** (2014), 1128-1157 (arXiv:1303.5144), Theorem 1: for
`1/2 < r < inf`, `1 < p_1, p_2, q_1, q_2 <= inf` with `1/r = 1/p_1 + 1/q_1 = 1/p_2 + 1/q_2`,
and `sigma > max(0, n/r - n)` or `sigma in 2N`,

    ||D^sigma (F G)||_{L^r} <= C ( ||D^sigma F||_{L^{p_1}} ||G||_{L^{q_1}}
                                 + ||F||_{L^{p_2}} ||D^sigma G||_{L^{q_2}} ).

The endpoint value `inf` is admitted for `p_i, q_i`. Take `n = 3`, `r = 2` (so the condition
is `sigma > 0`), `(p_1,q_1) = (2, inf)`, `(p_2,q_2) = (inf, 2)`, `F = G = u` componentwise:

    ||Lambda^sigma (u ⊗ u)||_{L^2} <= C ||u||_{L^inf} ||Lambda^sigma u||_{L^2},    sigma > 0.   (KP)

Grafakos-Oh Theorem 1 is stated on `R^n` for Schwartz functions. On `T^3` the same estimate
with identical index conditions is Proposition 1 of Á. Bényi, T. Oh, T. Zhao, *Fractional
Leibniz rule on the torus*, arXiv:2311.07998 (their condition (i) is the Grafakos-Oh
condition; `D^s = (-Delta)^{s/2}` on `T^d`, `f, g in C^inf(T^d)`); the authors state that
there is otherwise no standard reference for the periodic fractional Leibniz rule, so the
citation is not optional. On `T^3`, `Lambda^s` annihilates the mean mode; the mean of `u`
evolves by `d/dt mean(u) = mean(f)` (the nonlinear, pressure and dissipative terms have zero
mean) and is controlled by Step 1, so `Ḣ^s + L^2 = H^s` closes as on `R^3`.

**Step 1 (energy).** Pairing (NS_alpha) with `u`:

    (1/2) d/dt ||u||_{L^2}^2 + nu ||Lambda^{alpha/2} u||_{L^2}^2 = <f, u> <= ||f||_{L^2} ||u||_{L^2},

so `d/dt ||u||_{L^2} <= ||f||_{L^2}` and `||u(t)||_{L^2} <= ||u_0||_{L^2} + int_0^t ||f||_{L^2}`,
bounded on `[0,T)` by (A1).

**Step 2 (homogeneous `H^s`).** Apply `Lambda^s`, pair with `Lambda^s u`. The pressure term
drops because `Lambda^s u` is divergence free. Using `(u.grad)u = div(u ⊗ u)`:

    (1/2) d/dt ||Lambda^s u||_2^2 + nu ||Lambda^{s + alpha/2} u||_2^2 = -N + <Lambda^s f, Lambda^s u>,

where `N = <Lambda^s div(u ⊗ u), Lambda^s u>`. Integrating by parts once and splitting the
`2s+1` derivatives by Plancherel (`R_j = d_j Lambda^{-1}` is the Riesz transform, bounded on
`L^2`; the symbol identity is `|xi|^{s+1-alpha/2} . |xi|^{s+alpha/2} . (i xi_j/|xi|) = |xi|^{2s} i xi_j`):

    |N| <= ||Lambda^{s+1-alpha/2}(u ⊗ u)||_{L^2} . ||Lambda^{s+alpha/2} u||_{L^2}.

Set `X := ||Lambda^{s+alpha/2} u||_{L^2}` and `Y := ||Lambda^s u||_{L^2}`. Since
`sigma := s + 1 - alpha/2 > 0` (as `s > 1` and `alpha <= 2`), (KP) gives

    |N| <= C ||u||_{L^inf} ||Lambda^{s+1-alpha/2} u||_{L^2} . X.

**The exponent check.** Write `s + 1 - alpha/2 = s + theta (alpha/2)` with

    theta = (2 - alpha)/alpha.

Then `theta in [0,1]` if and only if `1 <= alpha <= 2`, which is exactly our range. By
Plancherel and Hölder with exponents `1/(1-theta), 1/theta`,

    ||Lambda^{s + theta alpha/2} u||_{L^2} <= Y^{1-theta} X^{theta}.

Hence `|N| <= C ||u||_{L^inf} Y^{1-theta} X^{1+theta}`. Now Young's inequality with conjugate
exponents `2/(1+theta)` and `2/(1-theta)` (**legitimate precisely because `theta < 1`, i.e.
`alpha > 1`**) gives, for any `nu > 0`,

    |N| <= (nu/2) X^2 + C(nu, alpha) ||u||_{L^inf}^{2/(1-theta)} Y^2,

and since `1 - theta = 2(alpha-1)/alpha`,

    2/(1 - theta) = alpha/(alpha - 1).

(Sanity checks: `alpha = 2` gives exponent `2`, the classical `||u||_inf^2`; `alpha = 3/2`
gives `3`; `alpha -> 1+` sends the exponent to `+inf`, and at `alpha = 1` exactly,
`theta = 1` and Young's inequality is unavailable. The threshold in the proof sits at the
same place as the threshold in the scaling.)

The force term is `<Lambda^s f, Lambda^s u> <= (1/2)||Lambda^s f||_2^2 + (1/2) Y^2`.
Collecting, with `y(t) = Y(t)^2`,

    y' + nu X^2 <= ( 1 + C(nu,alpha) M^{alpha/(alpha-1)} ) y + ||Lambda^s f||_{L^2}^2.

Grönwall on `[0,T)` with (A1) and (A2) gives `sup_{t<T} y(t) < inf`. With Step 1,
`sup_{t<T} ||u(t)||_{H^s} < inf`.

**Step 3 (continuation).** What is needed is local existence in `H^s`, `s > 5/2`, for
(NS_alpha) with `nu > 0` and force, on an interval whose length is bounded below in terms of
`||u(t_0)||_{H^s}` and the force; applied at `t_0` close to `T`, with (A1) holding on
`[0,T+delta)`, it extends the solution past `T`. For `alpha = 2` this is Kato's `H^s` theory.
For `1 < alpha < 2` it is the same Kato-type argument, and no citation is relied on: after
pairing with `Lambda^s u` the dissipative term is nonnegative, so the Euler estimate

    d/dt ||u||_{H^s} <= C ||grad u||_{L^inf} ||u||_{H^s} + ||f||_{H^s} <= C ||u||_{H^s}^2 + ||f||_{H^s}

(using `H^s ⊂ W^{1,inf}` for `s > 5/2`) holds with a constant independent of `nu`, and a
continuity argument gives existence up to time
`t_0 + c / ( ||u(t_0)||_{H^s} + ||f||_{L^1(t_0, t_0+1; H^s)} )`. Existence follows from this
a priori bound by a standard regularization or Galerkin scheme (for instance the
`eps(-Delta)^a` regularization with `eps`-uniform bounds and Aubin-Lions compactness used in
the proof of Wu's Theorem 3.1 cited below); uniqueness in `L^inf_t H^s` follows from the
`L^2` estimate on the difference `w` of two solutions,
`d/dt ||w||_2^2 <= C ||grad v||_{L^inf} ||w||_2^2`, the dissipative term again having the
good sign. Higher regularity propagates: `u in L^2(t_0, T; H^{s+alpha/2})` gives a.e. `t_1`
with `u(t_1) in H^{s+alpha/2}`, and the estimate of Step 2 at `s' = s + alpha/2` from `t_1`
on, iterated, yields every `s' > s`.

*Correction of the citation (2026-09-08, Wu 2003 read directly from the author PDF).* The
earlier version of this step cited J. Wu, *Generalized MHD equations*, J. Differential
Equations **195** (2003), 284-312, via Ma-Wu, Filomat 39:1 (2025), 239-247, p. 239, as
containing local strong solutions for `alpha > 0` from `H^1` data. It does not. Wu's only
local classical existence theorem is Theorem 3.1 (with Lemma 3.2): for the system
(3.1)-(3.2) with **no velocity dissipation** (`nu = 0`), `eta >= 0`, no external force, data
`(u_0, b_0) in H^m` with `m > max{2, beta} + d/2` (so `m > 7/2` in 3D, above the `s > 5/2`
used here), on `R^d` or `T^d`, with existence time depending only on the data. Its proof
regularizes with `eps(-Delta)^a u` (eq. (3.3)) with bounds uniform in `eps >= 0`, so the same
estimate covers a fixed `nu > 0`, but that is a re-derivation, not Wu's statement. For
`nu > 0`, `alpha > 0` Wu gives global weak solutions from `L^2` data (Theorem 2.2) and global
classical solutions for `alpha >= 1/2 + d/4` from `H^s` data, `s >= max{2 alpha, 2 beta}`
(Theorem 2.3). The Ma-Wu paraphrase is inaccurate on this point. Wu is therefore not relied
on here. ∎

### 1.4 Corollary and what it forbids

**Corollary A.** No blowup mechanism for (NS_alpha) with `1 < alpha <= 2` can have bounded
velocity. In particular the Córdoba-Martínez-Zoroa layered cascade, in the form used by
[CMZ-Z] (where `u in L^inf_{t,x}` is part of the conclusion) and, as reported, by
Alpöge-Buckmaster for Euler (`nu = 0`, velocity and circulation bounded), **cannot be pushed
to any `alpha > 1`, and a fortiori not to `alpha = 2`**, no matter how the force is chosen,
as long as the force satisfies (A1).

This closes the door on the naive extrapolation "Euler with force, then hypodissipative NS
with force, then full NS with force". The gap `[alpha_0, 1]` with `alpha_0 ≈ 0.0927` may be
technical; the gap `(1, 2]` is not, for this class of mechanism. Any claimed full-NS blowup
must have `sup_{t<T} ||u(t)||_{L^inf} = inf` (the limit statement `||u(t)||_{L^inf} -> inf`
needs in addition a lower bound on the local existence time in terms of `||u(t_0)||_{L^inf}`
for the forced equation, which is not supplied here).

### 1.5 What the `L^3` endpoint gives at `alpha = 2`, and an honest caveat

> **[ESS]** L. Escauriaza, G. A. Seregin, V. Šverák, *`L_{3,inf}`-solutions of the
> Navier-Stokes equations and backward uniqueness*, Russian Math. Surveys **58**:2 (2003),
> 211-250.

> **[Ser12]** G. Seregin, *A certain necessary condition of potential blow up for
> Navier-Stokes equations*, Comm. Math. Phys. **312** (2012), 833-845; arXiv:1104.3615.
> Result: if `T` is a blow-up time then `lim_{t↑T} ||v(.,t)||_{L^3} = inf`.

**Caveat, stated plainly.** The system as written in [Ser12], equation (1.1), is
`d_t v + v.grad v - Delta v = -grad q, div v = 0`. **There is no force term.** [ESS] is
likewise stated for the unforced Cauchy problem. So the `L^3` fence is verified only for
`f = 0`. It is plausible that it survives a force satisfying Fefferman's (5), because the
force is subcritical under the blow-up rescaling (`f_lambda = lambda^3 f(lambda x, lambda^2 t) -> 0`
as `lambda -> 0`), so the ancient solution produced by the blow-up procedure still solves the
unforced equations. This is a sketch only: the [ESS] argument needs the rescaled solutions to
be suitable weak solutions with a uniform local energy inequality including the force term,
and the limit to be an ancient unforced solution; none of that is carried out here. **I did
not find a published forced version and I am not asserting one.**

---

## 2. Lemma B: axisymmetric full-NS singularities lie on the axis

### 2.1 The partial regularity input

> **[CKN]** L. Caffarelli, R. Kohn, L. Nirenberg, *Partial regularity of suitable weak
> solutions of the Navier-Stokes equations*, Comm. Pure Appl. Math. **35** (1982), 771-831.

Hypothesis on the force, as used in the definition of a suitable weak solution and in the
epsilon-regularity proposition:

    f in L^q(D) for some q > 5/2,   div f = 0,

together with `u in L^inf(0,T;L^2) ∩ L^2(0,T;H^1)`, `p in L^{5/4}(D)`, and the local energy
inequality

    2 ∫∫ |grad u|^2 phi <= ∫∫ [ |u|^2 (phi_t + Delta phi) + (|u|^2 + 2p) u.grad phi + 2 (u.f) phi ],
    phi >= 0, phi in C_c^inf.

Conclusion (CKN Main Theorem B): the singular set `S` (the complement of the points having a
space-time neighbourhood on which `u` is essentially bounded) satisfies `P^1(S) = 0`, where
`P^1` is one-dimensional **parabolic** Hausdorff measure. CKN also prove existence of a
suitable weak solution for every such `f`.

*Verification note: I could not access CKN 1982 directly. I verified the force hypothesis
`f in L^q, q > 5/2, div f = 0`, the local energy inequality with its `2(u.f)phi` term, the
existence statement, and `P^1(S) = 0` against the Ulm seminar notes on CKN theory
(uni-ulm.de, `SeminarCKN207.pdf`, Definitions 1.1, 2.1, 4.1, Theorem 1.2, Proposition 1.5),
a secondary source that attributes each of these to [CKN82]. Alternative formulations exist
(F.-H. Lin, Comm. Pure Appl. Math. **51** (1998), 241-257; O. Ladyzhenskaya and G. Seregin,
J. Math. Fluid Mech. **1** (1999), 356-387); I did not check their force hypotheses.*

A smooth force that is compactly supported in space-time, or that satisfies Fefferman's (5)
below, lies in `L^q(R^3 x [0,T])` for every `q in [1, inf]`: taking `K = 4` in (5) gives
`|f| <= C(1+|x|)^{-4}`. So `q > 5/2` is satisfied with room to spare. (If the intended force
is not divergence free, replace `f` by its Leray projection and absorb the gradient part into
the pressure; this preserves smoothness and decay.)

### 2.2 Statement and proof

**Lemma B.** Let `alpha = 2`. Let `u` be an axisymmetric classical solution of (NS_2) on
`[0,T) x R^3` (or `T^3`) with smooth axisymmetric divergence-free data of finite energy and
smooth axisymmetric force `f` satisfying (A1) and `f in L^q`, `q > 5/2`, `div f = 0`, and
suppose `T` is a first singular time. Define the blow-up set

    Sigma := { x_0 : sup_{ B_rho(x_0) x (T - rho, T) } |u| = inf   for every rho > 0 }.

Then `Sigma` is contained in the axis of symmetry `{r = 0}`.

*Proof.* Three steps. The localisation of axisymmetric Navier-Stokes singularities to the
axis via CKN is a standard remark, not new here: the abstract of [CSTY] part I
(arXiv:math/0701796) attributes to [CKN] that axisymmetric strong solutions "could only blow
up on the axis of symmetry", and the proof of [CSTY] II Theorem 1.1 (section 6) reduces to
showing that every point on the axis is regular. Step (iii) records the measure-theoretic
detail.

(i) *`Sigma` is nonempty and rotation invariant.* Rotation invariance is immediate: `u` is
axisymmetric, so `|u|` is invariant under rotations about the axis, hence so is the defining
condition. Nonemptiness: if `Sigma = ∅` then every point has a neighbourhood on which `u` is
bounded near `T`; on `T^3` compactness upgrades this to `sup_{T-delta < t < T} ||u(t)||_{L^inf} < inf`
and Lemma A (with `alpha = 2`) contradicts `T` being singular. **On `R^3` this step needs an
extra input at spatial infinity** (a uniform decay statement for `u`, e.g. from decaying data
and force), because a countable cover of `R^3` gives no uniform bound. This is a genuine
scope condition and I flag it rather than paper over it.

(ii) *The CKN criterion applies to `u` itself on backward cylinders.* For `x_0` and small
`r > 0` the backward cylinder `Q_r(x_0,T) = B_r(x_0) x (T - r^2, T)` lies in `R^3 x [0,T)`
(or `T^3 x [0,T)`), where `u` is classical. On it `u` is a suitable weak solution:
`u in L^inf_t L^2 ∩ L^2_t H^1` by Step 1 of Lemma A (the energy inequality bounds
`nu ∫_0^T ||grad u||_2^2 dt`), `p in L^{5/3}_{t,x}` locally (Leray pressure; `u in L^{10/3}_{t,x}`
by interpolation, then Calderón-Zygmund), `f in L^q` with `q > 5/2`, and the local energy
inequality holds with equality for a classical solution. CKN Proposition 2 (Ulm notes
Proposition 5.1: if `limsup_{r->0} r^{-1} ∫∫_{Q_r(x_0,T)} |grad u|^2 <= eps_3` then `u` is
bounded on a smaller backward cylinder at `(x_0,T)`) therefore gives, for every
`x_0 in Sigma`, `limsup_{r->0} r^{-1} ∫∫_{Q_r(x_0,T)} |grad u|^2 > eps_3`. The Vitali covering
argument of CKN Theorem B (Ulm notes Theorem 11.1 with Lemma 11.3), run on these top-anchored
cylinders with `∫_0^T ||grad u||_2^2 dt < inf`, gives `P^1(Sigma x {T}) = 0`. No solution past
`T` and no uniqueness theorem is needed. (*The earlier version of this step extended `u` by a
CKN suitable weak solution `v` with the same data and force and identified `v = u` on `[0,T)`
by Prodi-Serrin weak-strong uniqueness, G. Prodi, Ann. Mat. Pura Appl. **48** (1959), 173-182,
J. Serrin, Arch. Ration. Mech. Anal. **9** (1962), 187-195. That route needs `v` to satisfy
the global energy inequality, which the CKN construction provides but the definition of a
suitable weak solution does not, and I had not verified a forced statement page-by-page; it
is replaced by the argument above.*)

(iii) *Parabolic measure on a single time slice, and the circle.* Covering `Sigma x {T}` by
parabolic cylinders `Q_{r_i} = B_{r_i} x (t_i - r_i^2, t_i)` and intersecting with the slice
`{t = T}` produces a cover of `Sigma` by Euclidean balls of the same radii `r_i`. Hence
`P^1(Sigma x {T}) = 0` if and only if `H^1(Sigma) = 0`, `H^1` being one-dimensional Hausdorff
measure in `R^3`. Now suppose `x_0 in Sigma` with `r_0 = dist(x_0, axis) > 0`. By (i) the
entire circle `C(x_0)` of radius `r_0` about the axis through `x_0` lies in `Sigma`, and
`H^1(C(x_0)) = 2 pi r_0 > 0`, contradicting `H^1(Sigma) = 0`. So every point of `Sigma` has
`r = 0`. ∎

### 2.3 Why this fails below `alpha = 2`

> **[TY]** L. Tang, Y. Yu, *Partial regularity of suitable weak solutions to the fractional
> Navier-Stokes equations*, Comm. Math. Phys. **334** (2015), 1455-1482 (erratum: Comm. Math.
> Phys. **335** (2015), 1057-1063). Their convention is `(-Delta)^a` with `3/4 < a < 1`, and
> the conclusion is that the suitable weak solution is regular away from a relatively closed
> singular set whose `(5 - 4a)`-dimensional Hausdorff measure is zero.

**Conversion to the `|grad|^alpha` convention.** `(-Delta)^a = |grad|^{2a}`, so `alpha = 2a`,
`a = alpha/2`, and the Tang-Yu bound reads

    dim(S) <= 5 - 4 (alpha/2) = 5 - 2 alpha,    valid for   3/2 <= alpha < 2

(endpoint `alpha = 3/2`, i.e. `a = 3/4`: W. Ren, Y. Wang, G. Wu, arXiv:1505.05585, Remark 1.2,
for the parabolic Hausdorff measure of the singular set of smooth solutions at the first
blow-up time). Tang-Yu assume a force `f in L^{q'}_{t,x}` with `q' > (9 + 6a)/(4a + 1)`
(Ren-Wang-Wu Remark 1.3); the measure is the parabolic Hausdorff measure (Ren-Wang-Wu
Remark 1.2).

At `alpha = 2` this returns `1`, recovering CKN. For every `alpha < 2` it exceeds `1`, so a
one-dimensional set such as a circle is **not** excluded. Strictly below `alpha = 3/2` no
theorem of this kind is available at all, and in particular none covers the [CMZ-Z] range
`alpha < 0.0927`. (For comparison on the other side: N. Katz and N. Pavlović, Geom. Funct.
Anal. **12** (2002), 355-379, arXiv:math/0104199, give `dim <= 5 - 4a` for the Hausdorff
dimension of the *spatial* singular set at the first blow-up time, for hyperdissipative
`(-Delta)^a`, `1 < a < 5/4`, i.e. `dim < 1` for `alpha > 2`; this is not the same object as
the space-time `P^k` statements of CKN and Tang-Yu.)

**Consequence.** An off-axis torus mechanism, in which the singular set is by construction a
circle around the axis, is compatible with the hypodissipative regime and with Euler, and is
*incompatible with `alpha = 2`* even setting aside Lemma A. The Alpöge-Buckmaster mechanism,
as reported, is excluded at `alpha = 2` by Lemma A because the velocity is bounded; Lemma B
excludes in addition any version of it whose singular set is an off-axis circle. Off-axis
support of the data and force does not by itself fix the singular set for `nu > 0`: the
support of `u` spreads instantly, and the singular point may sit on the axis (for Euler,
transport keeps the singular set inside the support geometry; for Navier-Stokes it need not).
So the second exclusion applies to the singular set, not to the support of the data.

### 2.4 Type I is excluded in the axisymmetric case (cited; unforced)

> **[KNSS]** G. Koch, N. Nadirashvili, G. Seregin, V. Šverák, *Liouville theorems for the
> Navier-Stokes equations and applications*, Acta Math. **203** (2009), 83-105;
> arXiv:0709.3599. Their definition (verbatim from the preprint): a singularity at time `T`
> is **type I** if `sup_x |u(x,t)| <= C/sqrt(T-t)`; type II is any singularity that is not
> type I. Theorem 6.2: let `u` be axisymmetric on `R^3 x (0,T)`, in `L^inf_{x,t}(R^3 x (0,T'))`
> for each `T' < T`, a weak solution satisfying `|u| <= C/sqrt(T-t)` on `R^3 x (0,T)` and, for
> some `R_0 > 0`, `|u(x,t)| <= C/sqrt(x_1^2 + x_2^2)` whenever `sqrt(x_1^2+x_2^2) >= R_0`
> (as holds for mild solutions with sufficiently fast decaying data). Then `|u| <= M(C)` on
> `R^3 x (0,T)`, so `u` is regular. **The system in [KNSS] is `u_t + u.grad u + grad p - Delta u = 0`,
> equation (1.1): unforced.**

> **[CSTY]** C.-C. Chen, R. M. Strain, T.-P. Tsai, H.-T. Yau, *Lower bound on the blow-up
> rate of the axisymmetric Navier-Stokes equations*, Int. Math. Res. Not. (2008), and *II*,
> Comm. Partial Differential Equations **34**:3 (2009), 203-232; arXiv:0709.4230. Part II,
> Theorem 1.1 (paraphrased from the PDF): let `(v,p)` be an axisymmetric strong solution of
> (N-S) in `D = R^3 x (-T_0, 0)` with `v|_{t=-T_0} = v^0 in H^{1/2}`, `r v_theta^0 in L^inf`,
> `p in L^{5/3}(D)`, and suppose either `|v(x,t)| <= C_* |t|^{-1/2}` on `D` or, for some
> `eps in [0,1]`, `|v(x,t)| <= C_* r^{-1+eps}|t|^{-eps/2}` on `D`; then
> `v in L^inf(B_R x [-T_0, 0])` for every `R > 0` (the case `eps = 0`, `|v| <= C_* r^{-1}`, is
> handled in their appendix). Part I (arXiv:math/0701796, abstract) proves the same
> conclusion under the bound `|v(x,t)| <= C_*(r^2 - t)^{-1/2}` for `-T_0 <= t < 0`. The
> system (N-S) is `d_t v + (v.grad)v + grad p = Delta v`, `div v = 0`; **no external force
> appears**.

So: **both Type I exclusion results are for `f = 0`**, and as stated they are not fences for a
*forced* blowup. The subcriticality remark of section 1.5 applies again, which makes a forced
version likely, but I did not find one published and do not assert it.

**Combined statement (with the above caveats).** Any axisymmetric blowup for full
Navier-Stokes must occur **on the axis** (Lemma B, proved above modulo the `R^3` decay input)
and, if Type I exclusion extends to forced solutions, must be **Type II**, which by [KNSS]'s
definition means that the bound `sup_x |u(x,t)| <= C/sqrt(T-t)` fails for *every* `C`, i.e.

    limsup_{t -> T} sqrt(T - t) . sup_x |u(x,t)| = inf.

So the velocity must blow up strictly faster than the self-similar rate. Combined with Leray's
lower bound `sup_x |u(x,t)| >= eps_1/sqrt(T-t)` (unforced; not re-proved here for `f != 0`),
the admissible window for an axisymmetric forced blowup is narrow on both sides.

---

## 3. Fefferman's statements, quoted

From C. L. Fefferman, *Existence and smoothness of the Navier-Stokes equation*, official Clay
Mathematics Institute problem description (claymath.org, `navierstokes.pdf`), fetched and text
extracted 2026-09-08. Verbatim:

Decay conditions:

    (4)  |d_x^a u_0(x)| <= C_{aK} (1 + |x|)^{-K}   on R^n, for any a and K
    (5)  |d_x^a d_t^m f(x,t)| <= C_{amK} (1 + |x| + t)^{-K}   on R^n x [0,inf), for any a, m, K

Acceptable-solution conditions on `R^3`:

    (6)  p, u in C^inf(R^n x [0,inf))
    (7)  ∫_{R^n} |u(x,t)|^2 dx < C   for all t >= 0   (bounded energy)

Periodic case: `u_0(x + e_j) = u_0(x)`, `f(x + e_j, t) = f(x,t)` (8); the Errata block at the
end of the same pdf states that the further condition `p(x + e_j, t) = p(x,t)` should be made
explicit in (8). The pdf also states that in the periodic case, in place of (4) and (5), `u_0`
is assumed smooth (periodic data are required only to be smooth and to satisfy (8)) and `f`
satisfies

    (9)  |d_x^a d_t^m f(x,t)| <= C_{amK} (1 + |t|)^{-K}   on R^3 x [0,inf), for any a, m, K
    (10) u(x,t) = u(x + e_j, t)   on R^3 x [0,inf) for 1 <= j <= n
    (11) p, u in C^inf(R^n x [0,inf))

The two breakdown statements:

> **(C) Breakdown of Navier-Stokes solutions on `R^3`.** Take `nu > 0` and `n = 3`. Then there
> exist a smooth, divergence-free vector field `u_0(x)` on `R^3` and a smooth `f(x,t)` on
> `R^3 x [0,inf)`, satisfying (4), (5), for which there exist no solutions `(p,u)` of (1),
> (2), (3), (6), (7) on `R^3 x [0,inf)`.

> **(D) Breakdown of Navier-Stokes Solutions on `R^3/Z^3`.** Take `nu > 0` and `n = 3`. Then
> there exist a smooth, divergence-free vector field `u_0(x)` on `R^3` and a smooth `f(x,t)`
> on `R^3 x [0,inf)`, satisfying (8), (9), for which there exist no solutions `(p,u)` of (1),
> (2), (3), (10), (11) on `R^3 x [0,inf)`.

Two things to note. First, (5) and (9) are strong: **every** space and time derivative of the
force must decay faster than **every** polynomial rate, uniformly. A force built by summing
layers of a cascade with growing frequencies and shrinking supports has to be checked against
this, not merely against smoothness and compact support at each fixed time. Second, (C) and
(D) ask for *nonexistence in the class*, not merely for one solution that breaks down. A
constructed blowup delivers this only together with a uniqueness statement in the Clay class.
A uniqueness theorem covering the class `C^inf ∩ L^inf_t L^2_x` of (6), (7) is needed;
weak-strong uniqueness (Prodi, Serrin) suffices only if the hypothetical competitor is shown
to be Leray-Hopf, with `grad u in L^2_{t,x}` and the energy inequality, which (6), (7) do not
give for free. The Clay pdf itself notes that uniqueness of weak solutions is not known.

---

## 4. Checklist for a claimed forced full Navier-Stokes blowup

Apply in order. A failure of item 1, of item 3 (axisymmetric case, with the `R^3` decay input
of 2.2(i)), of item 5 or of item 6 is decisive. Items 2 and 4 are expected but unproved for
`f != 0`: a failure there is a red flag, not a refutation.

1. **Is `sup_{t<T} ||u(t)||_{L^inf}` finite?** If yes, the claim is false, by Lemma A
   (`alpha = 2`, proved above), provided the force satisfies (A1) (implied by Fefferman's
   (5) or (9)). This is the fence that a Córdoba-Martínez-Zoroa layered cascade, and the
   reported Alpöge-Buckmaster Euler mechanism, fail. Ask for the claimed growth rate of
   `||u(t)||_{L^inf}` explicitly. (Lemma A gives `sup_{t<T} ||u(t)||_{L^inf} = inf`, not the
   limit. For `alpha = 2` the conclusion is classical, Leray 1934, and the Clay pdf itself
   states on page 3 that at a finite blowup time the velocity becomes unbounded near the
   blowup time; Lemma A's contribution is the fractional range `1 < alpha < 2` and the explicit
   force class (A1).)

2. **Does `||u(t)||_{L^3} -> inf` as `t -> T`?** Necessary for `f = 0` on `R^3` by [ESS] and
   [Ser12]; unproved here for `f != 0` and for `T^3` (both results are stated on `R^3`; no
   periodic version is cited). For a forced problem this is expected but, as far as I could
   verify, not published; treat a proof that produces bounded `L^3` as a red flag and ask
   the authors to address it. For `f = 0` on `R^3` the same holds for the critical Besov
   norms `Ḃ^{-1+3/p}_{p,q}`, `3 < p, q < inf` (I. Gallagher, G. Koch, F. Planchon, Comm. Math.
   Phys. **343** (2016), 39-82, arXiv:1407.4156; D. Albritton, Analysis & PDE **11** (2018),
   1415-1456, arXiv:1612.04439), with the same `f = 0` caveat.

3. **Where is the singular point, and is the flow axisymmetric?** If axisymmetric, every
   blow-up point lies on the axis `{r = 0}` (Lemma B; on `R^3` subject to the decay input of
   2.2(i)). A claim whose singular set contains an off-axis point, for instance a circle
   around the axis, is excluded at `alpha = 2` (`P^1` of a circle is positive). Off-axis
   support of the data or the force is not by itself excluded at `alpha = 2`: for `nu > 0` the
   support of `u` spreads instantly and the singular point may sit on the axis, so ask where
   the singular point is. Check that the force too is compatible: `f in L^q`,
   `q > 5/2`, `div f = 0` after Leray projection. If the flow is not axisymmetric, this fence
   does not apply and item 4 does not either.

4. **What is the rate?** In the axisymmetric case, Type I (`sup_x|u| <= C/sqrt(T-t)` for some
   `C`) is excluded by [KNSS] Theorem 6.2 and [CSTY], **for the unforced equations**; a claimed
   axisymmetric forced blowup at a Type I rate must either extend those theorems or explain why
   the force breaks them. If Type I exclusion extends to the forced axisymmetric problem on
   the relevant domain (both results are on `R^3`, not `T^3`, and [KNSS] Theorem 6.2 also
   assumes `|u(x,t)| <= C/sqrt(x_1^2 + x_2^2)` for `sqrt(x_1^2 + x_2^2) >= R_0`, which must be
   checked for the claimed solution), then the rate satisfies
   `limsup_{t->T} sqrt(T-t) sup_x|u| = inf` from above. Leray's
   `sup_x|u(x,t)| >= eps_1/sqrt(T-t)` from below is the unforced bound, stated here without a
   forced proof or citation.

5. **Energy bookkeeping.** With `f in L^1_t L^2_x`, `d/dt (1/2)||u||_2^2 <= ||f||_2 ||u||_2`,
   so `||u(t)||_{L^2} <= ||u_0||_2 + int_0^t ||f||_2`, which stays bounded on `[0,T)`, and for
   any Leray-Hopf continuation by the energy inequality (and uniformly on `[0,inf)` under (5),
   which makes `||f(t)||_{L^2}` integrable in `t`: (5) with exponent `K` gives
   `||f(t)||_2 <= C(1+t)^{3/2-K}`, integrable for `K > 5/2`). Fefferman's condition (7) is
   therefore never the obstruction. A scaling heuristic (at the Leray scale `l ~ sqrt(T-t)`
   with `|u| ~ 1/sqrt(T-t)` the energy in the concentrating ball is `~ (T-t)^{1/2}`) suggests
   that the concentrating scales carry vanishing kinetic energy, but a bounded `L^2` norm is
   compatible with a fixed amount of energy collapsing to a point, and no theorem forces this
   for a general blowup. The checkable request: give the energy budget of the cascade layer
   by layer and show that the sum of the layer energies stays bounded.

6. **Is it actually the Clay problem?** (a) Does the force satisfy (5) on `R^3` (or (9) on
   `T^3`) with the quantifier "for any `a`, `m`, `K`" honoured jointly, not just for each
   layer separately, with the pressure periodic on `T^3` (errata to (8))? (b) Is the data in
   the Clay class: (4) on `R^3`, or smooth and periodic per (8) on `T^3`? (c) Is the asserted
   conclusion nonexistence in the class `C^inf(R^3 x [0,inf))` with bounded energy (on `T^3`:
   periodic `u` and `p`), and is the uniqueness step that converts "this solution breaks down"
   into "no solution exists" supplied? (d) Is `nu > 0` fixed throughout (no `nu -> 0` limit),
   and is the final force a single smooth function of `(x,t)` satisfying (5) (or (9)) with the
   constants `C_{amK}` finite? Dependence of `f` and `u_0` on `nu` is allowed: (C) and (D)
   choose the data and force after `nu`, and `nu` can be normalised to `1` by
   `u(x,t) = nu w(x, nu t)`, which preserves (4), (5), (8), (9). Every constructed blowup
   defines `f` as the residual of an ansatz; what the Clay class requires is only that the
   resulting `f` be a fixed function of `(x,t)` in the class.

---

## 5. Verification record

| Claim | Source | Status |
|---|---|---|
| Clay statements (C), (D), conditions (4)-(11) | claymath.org `navierstokes.pdf` | Fetched, text extracted, quoted verbatim; errata on (8) (pressure periodicity) and the periodic-data sentence added 2026-09-08 |
| [CMZ-Z] `alpha_0 = (22-8 sqrt 7)/9`, dissipation `\|grad\|^alpha`, force in `L^1_t C^{1,eps}_x ∩ L^inf_t L^2_x` | arXiv:2407.06776 abstract via arXiv API | Verified verbatim. `alpha_0 ≈ 0.0927` computed directly |
| [CMZ-Z] journal reference ARMA **250** (2026) art. 38 | Springer link, doi:10.1007/s00205-026-02198-0 | Verified |
| [CMZ-Z] `u in L^inf_{t,x}`, blowup of `int_0^t \|\|grad u\|\|_inf` | Open-access full text (PMC13161257) | Verified |
| Kato-Ponce endpoint, Thm 1 index conditions | Grafakos-Oh CPDE **39** (2014) 1128-1157, author PDF | Fetched, theorem read verbatim |
| CKN force hypothesis `f in L^q, q>5/2, div f = 0`; local energy inequality; `P^1(S)=0`; existence | Ulm seminar notes on CKN (secondary) | Verified in secondary source only; CKN 1982 not accessed |
| [TY] `(-Delta)^a`, `3/4<a<1`, `(5-4a)`-dim parabolic measure zero, force `f in L^{q'}`, `q' > (9+6a)/(4a+1)`; CMP **334** (2015) 1455-1482 | Ma-Wu Filomat 39:1 (2025) p. 240; Ren-Wang-Wu arXiv:1505.05585 Remarks 1.2-1.3; Kwon-Ożański arXiv:2010.12105 intro | Verified through three independent secondary sources; Tang-Yu itself not accessed (paywall) |
| [KNSS] Type I definition, Thm 6.2, unforced system (1.1) | arXiv:0709.3599 PDF, text extracted | Verified verbatim |
| [CSTY] II Theorem 1.1 hypotheses (`eps in [0,1]`, `v^0 in H^{1/2}`, `r v_theta^0 in L^inf`, `p in L^{5/3}(D)`); part I bound `(r^2 - t)^{-1/2}`; CPDE **34**:3 (2009) 203-232 | arXiv:0709.4230 PDF; arXiv:math/0701796 abstract | Theorem 1.1 read from the PDF and paraphrased; no force in the system |
| [Ser12] `lim \|\|v(t)\|\|_{L^3} = inf`; system (1.1) unforced | arXiv:1104.3615 PDF, text extracted | Verified verbatim |
| Fractional Serrin criterion `2a/p + 3/q <= 2a-1`, `3/(2a-1) < q <= inf`, `1 <= a <= 3/2` for `(-Delta)^a` | Y. Zhou, Ann. IHP (C) **24**:3 (2007), 491-505, quoted in Ma-Wu Filomat 39:1 (2025) p. 240 | Form verified. **Note: Zhou's range `a in [1,3/2]` is `\|grad\|^beta` with `beta in [2,3]`, hyperdissipative only. It does NOT cover `1 < beta < 2`.** This is why Lemma A is proved here rather than cited for that range |
| Local well-posedness in `H^s` for (NS_alpha), `alpha > 0` | J. Wu, JDE **195** (2003) 284-312, author PDF read directly 2026-09-08 | Wu Thm 3.1 has `nu = 0`, no force, `H^m` with `m > max{2,beta} + d/2`; it does not state the result previously attributed to it. No longer relied on: Step 3 now gives the Kato-type argument |
| Weak-strong uniqueness with forcing | Prodi 1959, Serrin 1962 | **Not verified directly**; no longer used in Lemma B (step (ii) replaced by the backward-cylinder argument 2026-09-08); still relevant to the section 3 uniqueness caveat |
| Alpöge-Buckmaster 2026-09-08 forced Euler blowup | arXiv API (author `Alpoge`, 21 entries through 2026-09-02), web search | **Not found, not verified.** Treated throughout as reported, never as cited |
| OpenAI claim of full-NS blowup | none | Unpublished and unseen; nothing here depends on it |

**Bottom line.** Fence 1 (Lemma A) is proved here and is unconditional for `1 < alpha <= 2`
given (A1): a forced full-NS blowup must have unbounded velocity, so the layered-cascade
family, which keeps the velocity bounded, cannot reach `alpha = 2`. Fence 2 (Lemma B) is
proved here from CKN and confines an axisymmetric forced full-NS singular set to the axis,
which rules out an off-axis singular set such as a singular circle (not, by itself, off-axis
support of data and force), and it degrades continuously as `alpha`
drops below `2` (Tang-Yu bound `5 - 2 alpha`), which explains why the same geometry is
available in the hypodissipative and Euler regimes. The Type I exclusion and the `L^3`
endpoint remain verified only for `f = 0`, and I have not asserted forced versions.

---

## 6. Referee corrections (2026-09-08)

Each finding was re-checked against the note and, where a citation is at stake, against the
source (Wu 2003 author PDF; claymath.org `navierstokes.pdf`; arXiv PDFs 0709.4230, 1505.05585,
2010.12105, 2311.07998; arXiv abstracts math/0104199, math/0701796, 1407.4156, 1612.04439;
Ulm `SeminarCKN207.pdf`). Applied unless stated otherwise.

1. *1.3 Step 2, pressure term.* Confirmed: `Lambda^s` commutes with `div`, so the pressure
   pairing vanishes, and `Lambda^s p in L^2` on `R^3` for the Leray pressure by (KP). The only
   loose end was the solution class. Applied: Lemma A now states the strong solution class
   `u in C([0,T);H^s) ∩ L^2_loc([0,T);H^{s+alpha/2})` with `p` the Leray pressure, plus a
   parenthetical on identifying a Clay-class solution with it.
2. *Kato-Ponce on `T^3`.* Confirmed: Grafakos-Oh Theorem 1 is on `R^n`; Bényi-Oh-Zhao
   arXiv:2311.07998 Proposition 1 gives the periodic statement with the same index conditions
   and says there is otherwise no standard reference. Applied, with the mean-mode remark.
3. *1.3 Step 3, Wu 2003.* Confirmed MAJOR: Wu's Theorem 3.1 has `nu = 0`, no force, data in
   `H^m` with `m > max{2, beta} + d/2`; Theorems 2.2 and 2.3 are global results. The earlier
   attribution (via Ma-Wu) was wrong. Applied as option (b): `s > 5/2` kept, the Kato-type
   local existence argument supplied in Step 3, Wu's actual content stated, Wu no longer
   relied on. Table row corrected.
4. *(A1) on `(0,T)` only; smoothing for `s' > s`.* Confirmed. Applied: (A1) now on
   `[0,T+delta)`; one sentence on parabolic smoothing added to Step 3.
5. *2.2 step (ii), weak-strong uniqueness.* Confirmed: the CKN definition of a suitable weak
   solution carries only the local energy inequality, and the step is not needed. Applied:
   step (ii) replaced by the backward-cylinder argument (CKN Proposition 2 applied to the
   classical solution on `Q_r(x_0,T) ⊂ R^3 x [0,T)`, then the Vitali cover); the old route is
   recorded in a parenthetical. Table row and header updated.
6. *2.3 Tang-Yu.* Confirmed: Ren-Wang-Wu Remark 1.3 gives the force hypothesis
   `q' > (9+6a)/(4a+1)` and Remark 1.2 gives the parabolic measure and the endpoint `a = 3/4`
   for smooth solutions at the first blow-up time; Kwon-Ożański arXiv:2010.12105 states the
   `(5-4s)` result. Applied: range `3/2 <= alpha < 2` with the endpoint attribution, force
   hypothesis added, parabolic-vs-Euclidean caveat replaced, third source added to the table.
7. *Katz-Pavlović.* Confirmed from the arXiv abstract ("singular set at time of first blow
   up"). Applied: "spatial singular set at the first blow-up time", with the remark that this
   is not the CKN/Tang-Yu space-time object.
8. *2.4 [CSTY] "abstract (verbatim)".* Confirmed: the quoted text was a paraphrase of
   Theorem 1.1, with `eps in [0,1]` and the hypotheses `v^0 in H^{1/2}`, `r v_theta^0 in L^inf`,
   `p in L^{5/3}(D)` dropped; part I's bound is `(r^2 - t)^{-1/2}`. Applied: blockquote
   rewritten as a labelled paraphrase of Theorem 1.1 with the part I bound stated separately;
   table row corrected.
9. *Wu 2003 misattributed (second lens).* Same defect as finding 3; applied there.
10. *Section 3, Clay errata and periodic-data sentence.* Confirmed from the pdf text
    ("The further condition p(x + e_j, t) = p(x, t) should be made explicit in Eqn (8)";
    "In place of (4) and (5), we assume that u_0 is smooth"). Applied to section 3, items
    6(a), 6(b), 6(c) and the table row.
11. *Energy bookkeeping wording.* Confirmed (`(5)` gives `||f(t)||_2 <= C(1+t)^{3/2-K}`).
    Applied: "up to and beyond `T`" replaced; the standing solution class added to Lemma A
    (see 1).
12. *Section 4 header.* Confirmed: items 2 and 4 are verified only for `f = 0`, so "the first
    failure is decisive" contradicted the note's own caveats. Applied: decisive restricted to
    items 1, 3, 5, 6; items 2 and 4 labelled red flags.
13. *Item 1, sup versus limit; Leray 1934.* Confirmed: Lemma A yields `sup = inf`, and the
    Clay pdf page 3 states that the velocity becomes unbounded near a finite blowup time.
    Applied in 1.4 and item 1; the forced Leray lower bound is stated as not supplied rather
    than added.
14. *Item 2, `R^3` only, Besov sentence.* Confirmed. Applied: item 2 marked necessary for
    `f = 0` on `R^3` and unproved for `f != 0` and `T^3`; Besov sentence now cites
    Gallagher-Koch-Planchon (arXiv:1407.4156, CMP 343 (2016) 39-82, `3 < p, q < inf`) and
    Albritton (arXiv:1612.04439, APDE 11 (2018) 1415-1456), both `f = 0`, journal data from
    the arXiv records; the 1.5 rescaling sketch is labelled as a sketch.
15. *Item 3 and Lemma B: `R^3` caveat, off-axis support, "new".* Confirmed on all three
    points; the [CSTY] part I abstract attributes the on-axis localisation to [CKN], and
    [CSTY] II section 6 reduces to axis points. Applied: item 3 rewritten as proposed, the
    2.3 "Consequence" paragraph and the bottom line adjusted from "geometry excluded" to
    "off-axis singular set excluded", and "only the third is new" replaced by the citation
    check. Leonardi-Málek-Nečas-Pokorný 1999 was not accessed and is not cited.
16. *Item 4, unconditional "must".* Confirmed. Applied: item 4 made conditional, the [KNSS]
    decay hypothesis listed as a check, Leray's lower bound marked unforced (also in 2.4).
17. *Item 5, "must carry vanishing kinetic energy".* Confirmed heuristic. Applied: labelled as
    a scaling heuristic and replaced by the checkable request.
18. *Item 6(d), force independent of `nu`.* Confirmed: (C), (D) choose data and force after
    `nu`, and `u(x,t) = nu w(x, nu t)` normalises `nu` while preserving (4), (5), (8), (9).
    Applied as proposed.
19. *Section 3, "weak-strong uniqueness suffices".* Confirmed over-claim. Applied: replaced by
    the hedged statement. No Galdi citation added, since I did not check one.
20. *Table row, Clay quotations.* Confirmed. Applied.

Not touched: everything outside the passages named above.

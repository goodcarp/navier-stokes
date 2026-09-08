# Forced axisymmetric Navier-Stokes in the Clay (C) class: singularities are on the axis and are Type II

Date: 2026-09-08. Seat: sub-seat of the HOLD UP SHIPS campaign.
Companion: `FENCES_forced_NS_blowup.md` in this directory (referred to below as [FENCES]).
Scope: 3D incompressible Navier-Stokes with `nu > 0` and a smooth external force, on `R^3`,
in Fefferman's class (C); one section on `T^3` (class (D)).

This note closes the gap [FENCES] section 2.4 left open. That section recorded that both
Type I exclusion theorems for axisymmetric Navier-Stokes, Koch-Nadirashvili-Seregin-Šverák
2009 and Chen-Strain-Tsai-Yau 2008/2009, are stated for the **unforced** equations, and
declined to assert a forced version. Here the forced version is proved, by going through the
KNSS route step by step. Every step carries a label:

* **CITED-VERBATIM** - the statement is used exactly as it appears in the source, with the
  source read at the URL recorded in section 10.
* **PROVED-HERE** - the argument is written out in this note.
* **ROUTINE-EXTENSION** - the source's argument goes through with a force, and the reason
  is stated. Where the reason is not obvious the extension is written out and relabelled
  PROVED-HERE.

Section 8 lists the steps that do **not** go through, and what is assumed instead.

---

## 1. Setting

Write `r = sqrt(x_1^2 + x_2^2)`, `z = x_3`, and use cylindrical coordinates `(r, theta, z)`
with the frame `e_r, e_theta, e_z`. A vector field `u` on `R^3` is *axisymmetric* if
`u(Qx) = Q u(x)` for every rotation `Q` about the `x_3` axis; equivalently
`u = u_r(r,z,t) e_r + u_theta(r,z,t) e_theta + u_z(r,z,t) e_z`. It is *swirl free* if
`u_theta = 0`. The *poloidal part* is `b = u_r e_r + u_z e_z`.

The system is

    d_t u + (u . grad) u + grad p = nu Laplace u + f,   div u = 0,   u(.,0) = u_0.   (NS)

Fefferman's conditions, quoted in [FENCES] section 3 from the Clay problem description:

    (4)  |d_x^a u_0(x)| <= C_{aK} (1 + |x|)^{-K}                on R^3, for any a, K
    (5)  |d_x^a d_t^m f(x,t)| <= C_{amK} (1 + |x| + t)^{-K}     on R^3 x [0,inf), for any a, m, K
    (6)  p, u in C^inf(R^3 x [0,inf))
    (7)  int_{R^3} |u(x,t)|^2 dx < C   for all t >= 0

**(C)** asks for `u_0` and `f` satisfying (4), (5) for which no `(p,u)` satisfying (6), (7)
solves (NS) on `R^3 x [0,inf)`.

Throughout, `T < inf` is the *first singular time* of a solution `(u,p)` of (NS) which is
`C^inf` on `R^3 x [0,T)` and satisfies (7) there: `(u,p)` does not extend to a solution of
(NS) that is `C^inf` on `R^3 x [0,T']` and satisfies (7), for any `T' > T`.

### 1.1 The standing hypothesis that (4)-(7) do not supply

    (H)  for every T' < T and every integer k >= 0,  sup_{R^3 x [0,T']} |grad^k u| < inf.

(H) says the solution is *uniformly* smooth on every slab strictly before `T`, including at
spatial infinity. Smoothness plus (7) does not give this: a countable cover of `R^3` yields
no uniform bound, and (7) controls only an integral. This is the same spatial-infinity input
that [FENCES] section 2.2 step (i) flagged and refused to paper over, and it is what makes
`R^3` harder here than `T^3` or a bounded domain, where compactness supplies it for free.

**(H) is not enough, and is not what is carried below** (referee finding R1, 2026-09-08). Three
steps need `u(t)` to lie in an `L^2`-based space, which no `L^inf` bound supplies: S3, which is
[FENCES] Lemma A and is stated there for the strong class; part (i), whose Lemma B proof uses
`u in L^inf_t L^2 ∩ L^2_t H^1` and `int_0^T ||grad u||_2^2 dt < inf` up to `T`; and S19.3, S20,
whose Grönwall arguments presuppose `||omega(t)||_{L^2}` and `||u(t)||_{H^s}` finite. The two
classes really differ: with `psi in C_c^inf(B_1)`, points `x_j -> inf` pairwise `4`-separated,
`w_j = 1/log j`, `a_j = ((log j)/j)^{1/2}`, the divergence-free field
`u = curl( sum_j a_j w_j psi((x-x_j)/w_j) e )` has every derivative uniformly bounded and
`||u||_{L^2}^2 ~ sum 1/(j (log j)^2) < inf`, but `||grad u||_{L^2}^2 ~ sum 1/j = inf`. So the
hypothesis actually carried below is

    (H*)  u in C([0,T'); H^s(R^3))  for every s >= 0 and every T' <= T.

(H*) implies (H) by Sobolev embedding (`H^{k+2} -> W^{k,inf}` on `R^3`), and implies
`int_0^T ||grad u||_{L^2}^2 dt < inf` through the energy identity ([FENCES] Lemma A Step 1). It
is exactly the solution class of [FENCES] Lemmas A and B. It is **not** implied by (4)-(7):
identifying a Clay-class solution with the strong solution needs a uniqueness theorem in the
class `C^inf ∩ L^inf_t L^2_x`, which [FENCES] section 3 records as unavailable (the energy proof
of uniqueness needs the difference in `L^2_t H^1`, needs `u in L^inf` to kill
`int_{|x|=R} u_n |w|^2`, and needs a pressure normalisation to kill `int_{|x|=R}(p - p~)w_n`,
since (6) imposes no decay on `p`), and it also needs the `H^s` norm not to escape to spatial
infinity before `T`. Neither is closed here. (H*) is carried explicitly in every statement
below, **including part (i)**, which needs its energy consequence though not (H) itself.

A second, separate decay hypothesis is needed for one half of part (ii) only, and is quoted
verbatim from KNSS, where it is their (6.6):

    (Dec)  there exist R_0 > 0 and C_0 < inf with |u(x,t)| <= C_0 / r
           whenever r >= R_0 and 0 < t < T.

KNSS remark after Theorem 6.2 that (Dec) holds "when `u_0` decays sufficiently fast at
infinity", citing Brandolese; that remark is for the unforced problem and is not re-derived
here for `f != 0`. (Dec) is *implied* by the CSTY-form Type I bound, which is why part (ii)
splits into an unconditional half and a conditional half.

---

## 2. The theorem

> **Theorem.** Let `nu > 0`. Let `u_0` be smooth, divergence free and axisymmetric on `R^3`
> and satisfy (4). Let `f` be smooth and axisymmetric on `R^3 x [0,inf)` and satisfy (5).
> Let `(u,p)` solve (NS), be `C^inf` on `R^3 x [0,T)` with bounded energy (7) there, let
> `T < inf` be its first singular time, and assume (H*) of section 1.1. Define the blow-up set
>
>     Sigma := { x_0 in R^3 : sup_{ B_rho(x_0) x (T - rho, T) } |u| = inf  for every rho > 0 }.
>
> Then:
>
> **(i) On the axis.** `Sigma` is contained in the axis of symmetry `{r = 0}`. This uses (H*)
> only through its energy consequence `u in L^inf(0,T;L^2) ∩ L^2(0,T;H^1)`, and does not use
> (H) itself or (Dec). Under (H*) together with `sup_{t<T} |u(x,t)| -> 0` as `|x| -> inf`,
> `Sigma` is also nonempty.
>
> **(ii) Type II.** Assume (H*). Then:
>
> **(ii-a)** There is no constant `C_* < inf` with
>
>     |u(x,t)| <= C_* ( r^2 + (T - t) )^{-1/2}    on R^3 x (0,T).      (T1-CSTY)
>
> This is the bound of Chen-Strain-Tsai-Yau part I, Theorem 1.1. No decay hypothesis beyond
> (H*) is needed, because (T1-CSTY) implies `|u| <= C_*/r` on all of `R^3 x (0,T)`.
>
> **(ii-b)** Assume in addition (Dec). Then there is no constant `C_* < inf` with
>
>     |u(x,t)| <= C_* (T - t)^{-1/2}              on R^3 x (0,T).      (T1-KNSS)
>
> Equivalently, in the KNSS terminology the singularity at `T` is of **type II**:
>
>     limsup_{t -> T-}  sqrt(T - t) . sup_{x in R^3} |u(x,t)|  =  inf.
>
> Since (T1-CSTY) implies (T1-KNSS), (ii-b) is the stronger exclusion where (Dec) is
> available, and (ii-a) is the exclusion available without it.
>
> **(iii) Swirl.** Assume (H*). If in addition `u_0` and `f` are swirl free, that is
> `u_{0,theta} = 0` and `f_theta = 0` identically, then no such `T` exists: the solution is
> global and smooth. Consequently, inside the axisymmetric sub-class of (C), any construction
> that breaks down in finite time must have `u_{0,theta} != 0` or `f_theta != 0`.

**Reading.** Restricted to the axisymmetric sub-class of the Clay (C) problem, a winning
construction must place its singularity on the axis, must blow up strictly faster than the
self-similar rate `(T-t)^{-1/2}` (unconditionally in the `(r^2 + (T-t))^{-1/2}` form), and
must carry swirl in the data or in the force, all of it under (H*). Combined with [FENCES]
Lemma A, which says the velocity must be unbounded, this is a narrow corridor: `sup_x |u(x,t)|` must diverge, and
faster than `(T-t)^{-1/2}`, at a point on the axis, in a flow with nonzero swirl.

---

## 3. The force under the blow-up rescaling

### 3.1 Reduction to a divergence-free force

**Step S0 (PROVED-HERE).** Write `f = g + grad phi` with `phi = Laplace^{-1} div f` and
`g = P f` the Leray projection. Then `(u, p)` solves (NS) with force `f` if and only if
`(u, p - phi)` solves (NS) with force `g`, and `div g = 0`. Under (5), `div f(.,t)` is
smooth with all derivatives decaying faster than every polynomial rate, so `phi` is smooth
and `p - phi` is smooth: the Clay condition (6) survives the replacement. So we may and do
assume `div f = 0`, writing `g` for the force from here on. Note that `g` itself does **not**
satisfy (5): by S1 its decay stops at `(1+|x|)^{-3}`, and generically no faster. Nothing below
needs it to; (A1) survives because the Riesz transforms are bounded on `L^2` and on `H^s`, and
S6 needs only the `(1+|x|)^{-3}` bound. (This is the same replacement
[FENCES] section 2.1 makes for the CKN hypothesis.)

### 3.2 What (5) buys, and why boundedness would not have been enough

**Step S1 (PROVED-HERE).** Let `f` satisfy (5) and `g = P f`. Then for all `k >= 0` and all
`K`,

    |grad^k g(x,t)| <= C_{kK} (1 + |x|)^{-3} (1 + t)^{-K},           (S1a)
    ||grad^k g(.,t)||_{L^q(R^3)} <= C_{kKq} (1 + t)^{-K},  1 < q <= inf.   (S1b)

*Proof.* `g_i = f_i + R_i R_j f_j` with `R_j` the Riesz transforms. Fix `t`. The kernel of
`R_i R_j` is `c delta + k(x)` with `k` smooth away from `0`, homogeneous of degree `-3` and
of mean zero on spheres. For `|x| >= 1`, split `int k(x-y) f_j(y) dy` at `|y| = |x|/2`. On
`|y| <= |x|/2` we have `|k(x-y)| <= C|x|^{-3}` and the integral is bounded by
`C|x|^{-3} ||f_j(.,t)||_{L^1}`. On `|y| > |x|/2` the kernel is still singular at `y = x`, so that region must be split once
more (referee finding R5, 2026-09-08): cut out `|y - x| <= |x|/4`, on which `f 1_{|y-x|<=|x|/4}`
is Schwartz with every seminorm `O_N(|x|^{-N})` by (5), and use
`||R_iR_j h||_{L^inf} <= C ||hhat||_{L^1} <= C ||h||_{W^{4,1}}`; on the rest, `k` is
non-singular and the integral converges absolutely against the decay of `f`. (The earlier
version of this line said "`k` is integrable against the remaining decay", which is false:
`k` is homogeneous of degree `-3` on `R^3` and is not locally integrable.) Hence
`|R_iR_jf_j(x,t)| <= C(1+|x|)^{-3}` with the `t` decay inherited from (5). Derivatives commute
with `R_iR_j`. (S1b) follows by integration; note `q > 1` is sharp, `(1+|x|)^{-3}` is not in
`L^1(R^3)`. `∎`

**Why this matters.** KNSS section 3 record explicitly that solutions of the Stokes system
with a right-hand side in `L^inf_{x,t}` that is *not* in divergence form "are not well
defined", because `P` maps `L^inf(R^n)` into `BMO(R^n)` only modulo constants. Their entire
mild-solution framework is built on right-hand sides in divergence form (`f_k = -u_k u`)
precisely to avoid this. A merely smooth and bounded force therefore does **not** fit the
framework. (5) is exactly what repairs it: by (S1b) `g(.,t)` lies in every `L^q`, `1 < q`,
and the Leray projection is an honest function with no additive ambiguity. This is a genuine
scope condition of the theorem, and it is the force hypothesis that is load-bearing. It is
the exact analogue of the scope condition (A1) that [FENCES] Lemma A carries.

### 3.3 The scaling

**Step S2 (PROVED-HERE; machine-checked).** If `(u,p,g)` solves (NS) then so does

    u_L(x,t) = L u(Lx, L^2 t + T),  p_L(x,t) = L^2 p(Lx, L^2 t + T),
    g_L(x,t) = L^3 g(Lx, L^2 t + T),

for every `L > 0`, and `div u_L = 0`. Hence for every `j >= 0`

    || grad_x^j g_L ||_{L^inf} = L^{3+j} || grad^j g ||_{L^inf}.      (S2a)

The exponent `3` on the force is what makes `(S2a)` a *gain* as `L -> 0`: the force is
strongly subcritical under the blow-up scaling. The covariance and the exponent are verified
symbolically in `check_forced_axisym.py`, check C1, together with a control showing that the
exponent `2` does not work.

### 3.4 The velocity is unbounded at `T`

**Step S3 (CITED-VERBATIM, from the companion note).** Under (4), (5), (7) and `nu > 0`, at
the first singular time

    sup_{0 <= t < T} || u(t) ||_{L^inf(R^3)} = inf.

This is [FENCES] Lemma A at `alpha = 2`, whose hypothesis (A1) is implied by (5), proved
there **with the force present**. The label is CITED-VERBATIM only once (H*) is in force:
Lemma A is stated for the strong class `u in C([0,T);H^s) ∩ L^2_loc([0,T);H^{s+1})`, and its own
parenthetical requires a Clay-class solution to be identified with the strong solution first.
(H*) is that identification, assumed rather than proved; see section 1.1 and F-2. It replaces Leray's lower bound `sup_x|u| >= eps_1/sqrt(T-t)`,
which KNSS quote at the head of their blow-up construction and which is an unforced result
(recorded as such in [FENCES] item 4). The construction below needs only `M_k -> inf`, and S3
supplies it. This is one place where the two fences fit together: Lemma A supplies exactly
the input the KNSS rescaling needs, and Lemma A is already forced.

---

## 4. Forced mild solutions on `R^3`

All of this section is the KNSS sections 3 and 4 machinery with the extra Duhamel term
carried along. Write `S(t)` for the heat semigroup, `K_{ijk}` for the kernel in KNSS (3.5),
and

    G(t) := int_0^t S(t - s) g(s) ds        (g divergence free, so P g = g).

**Step S4 (ROUTINE-EXTENSION).** *Definition and estimates.* A bounded measurable
`u : R^3 x (0,T) -> R^3` is a **forced mild solution** with datum `u_0 in L^inf` and force
`g` if

    u_i(x,t) = int Gamma(x-y,t) u_{0i}(y) dy
             + int_0^t int K_{ijk}(x-y, t-s) ( - u_k u_j )(y,s) dy ds
             + G_i(x,t).                                                   (4.1)

*Reason the extension is routine.* The added term is an explicit linear functional of `g`
that is not singular: since `g` is smooth with bounded derivatives (S1),

    || grad^k G ||_{L^inf(R^3 x (0,T))} <= T || grad^k g ||_{L^inf},        (4.2)
    || d_t G ||_{L^inf(R^3 x (0,T))} <= || g ||_{L^inf} + T || Laplace g ||_{L^inf},  (4.3)

both by differentiating under the integral sign. In particular no smoothing estimate is
consumed by the force term, and the bilinear estimate KNSS (4.5),
`||B(u,v)||_inf <= C sqrt(T) ||u||_inf ||v||_inf`, is untouched. Consequently:

* *(S4a) Local existence and uniqueness.* `u = U + B(u,u) + G` is solved by the contraction
  argument of KNSS section 4 on `[0,T_1]` with `T_1 = T_1( ||u_0||_inf, ||g||_inf ) > 0`
  bounded below in terms of those two quantities only. Uniqueness is immediate from the
  formula, as in KNSS.
* *(S4b) Blow-up criterion.* If `T_max < inf` then `||u(t)||_{L^inf} -> inf` as `t -> T_max`,
  since otherwise (S4a) would continue the solution. (This is an independent route to S3;
  we use S3, which is a statement about the Clay-class solution.)
* *(S4c) Smoothing.* KNSS Proposition 4.1 holds with `u_0` replaced by `(u_0, g)` and the
  right-hand side of (4.6) replaced by `C(k,l)(||u_0||_inf + T ||grad^{k} g||_inf + ...)`,
  by (4.2)-(4.3). Only boundedness of finitely many derivatives of `g` is used.

**Step S5 (PROVED-HERE).** *Forced decomposition (KNSS Lemma 3.1 with a force).* Let
`u in L^inf(R^3 x (0,T))` be a weak solution of the forced Stokes system

    d_t u + grad p - Laplace u = d_{x_k} F_k + g,   div u = 0,

with `F in L^inf_{x,t}` and `g` as in S1, in the sense that
`int int u (phi_t + Laplace phi) = int int F_k d_{x_k} phi - int int g . phi` for all
smooth compactly supported divergence-free `phi`. Then

    u = v + G + w + b(t),

where `v` is the mild solution with datum `0` and divergence-form data `F`, `w` solves the
heat equation, and `b` depends on `t` only, with `||w||_inf + ||b||_inf <= C(T) ||u||_inf`.

*Proof.* `G` is divergence free, bounded by (4.2), and satisfies
`int int G (phi_t + Laplace phi) = - int int g . phi` for divergence-free `phi`, since
`d_t G - Laplace G = g` classically and `G(.,0) = 0`. Hence `u - G` is a bounded weak
solution of the *unforced* Stokes system with the same divergence-form data, and KNSS
Lemma 3.1 applies to it verbatim. `∎`

**Step S6 (PROVED-HERE).** *Decay of `G`.* For `0 <= t <= T`,

    |G(x,t)| <= C(T) (1 + |x|)^{-3},   and hence   |G(x,t)| <= C'(T) / r  on R^3 x (0,T).

*Proof.* By (S1a), `|g(y,s)| <= A (1+|y|)^{-3}`. For `0 < tau <= T`, split
`int Gamma(x-y,tau) (1+|y|)^{-3} dy` at `|y| = |x|/2`. The outer piece is at most
`(1+|x|/2)^{-3}`. On the inner piece `|x-y| >= |x|/2`, so the heat kernel is at most
`sup_{0<tau<=T} (4 pi tau)^{-3/2} e^{-|x|^2/(16 tau)}`; writing `sigma = |x|^2/tau >= |x|^2/T`
this equals `c |x|^{-3} sup_{sigma >= |x|^2/T} sigma^{3/2} e^{-sigma/16}`, which for `|x|`
large is `<= C_N |x|^{-N}` for every `N`, and the inner integral of `(1+|y|)^{-3}` grows only
logarithmically. So `|S(tau) g(.,s)(x)| <= C (1+|x|)^{-3}` uniformly for `tau <= T`, and
`|G(x,t)| <= T sup <= C(T)(1+|x|)^{-3}`. The second bound follows since `(1+|x|)^{-3} <= r^{-3} <= r^{-1}`
for `r >= 1` and `|x| >= r`, while for `r <= 1` the bound `C'/r` exceeds `||G||_inf`. `∎`

**Step S7 (PROVED-HERE).** *Compactness with vanishing force (KNSS Lemmas 4.1 and 6.1 with a
force).* Let `u^{(k)} in L^inf_loc(R^3 x (A_k, 0))` be forced mild solutions with forces
`g^{(k)}`, with `A_k -> -inf`, and

    || grad^j g^{(k)} ||_{L^inf}  ->  0   for j = 0, 1, 2.                    (4.4)

**(a)** If `|u^{(k)}| <= C` uniformly on `R^3 x (A_k, 0)`, then a subsequence converges locally
uniformly on `R^3 x (-inf, 0]`, endpoint included, to an **ancient mild solution of the
unforced Navier-Stokes equations** `u` with `|u| <= C`.
**(b)** If instead only `sup_{R^3 x (A_k,-delta)} |u^{(k)}| <= C(delta)` for every `delta > 0`,
the same holds with convergence locally uniform on the open set `R^3 x (-inf, 0)`.

(The two forms are used respectively in S8/S10 and in S13, where `w^{(k)}` obeys only
`|w^{(k)}| <= C/sqrt(-tau)`; referee finding R3, 2026-09-08.)

*Proof.* Fix `S > 0` and take `k` with `A_k < -S`. Restarting the Duhamel formula at time
`-S`,

    u^{(k)}(t) = S(t+S) u^{(k)}(-S)
               + int_{-S}^t S(t-s) P div( - u^{(k)} (x) u^{(k)} )(s) ds
               + int_{-S}^t S(t-s) g^{(k)}(s) ds,   -S < t < 0,

and the last term is bounded in `L^inf` by `S ||g^{(k)}||_inf -> 0`, with its first two
spatial derivatives bounded by `S ||grad^j g^{(k)}||_inf -> 0` by (4.2). The first two terms
are handled exactly as in KNSS Lemma 4.1: KNSS (3.10) gives uniform parabolic Hölder bounds
on the bilinear term from `||u^{(k)}||_inf <= C`, and heat smoothing gives equicontinuity of
the first term on `t >= -S + delta`. **This reaches `t = 0`**: KNSS (3.10) bounds
`||u||_{C^alpha_par(Q(z_0,R))}` on every backward parabolic cylinder
`Q(z_0,R) = B(x_0,R) x (t_0 - R^2, t_0)` contained in the time interval, and such a cylinder may
end at the final time, so in case (a) the equicontinuity holds on `R^3 x [-S+delta, 0]` and the
limit is attained at `s = 0` as well. Extract a locally uniformly convergent subsequence,
pass to the limit in the identity above; the force term vanishes in the limit, so the limit
satisfies the unforced Duhamel identity on `(-S, 0)`. Diagonalise over `S = 1, 2, 3, ...`.
The limit is therefore an ancient mild solution of the unforced equations. `∎`

**This is the structural point of the whole note.** The Liouville theorems are applied to the
limit, and the limit is unforced. No Liouville theorem has to be extended.

**Step S8 (PROVED-HERE).** *Forced Proposition 6.1.* Let `(u,p)` be as in the Theorem and
assume (H*). Then `u` restricted to `R^3 x (0,T)` is a bounded weak solution of forced (NS) on
every `R^3 x (0,T')`, `T' < T`, and, *provided the hypothesis at hand forces `b = 0` in S5*, a
finite-time singularity at `T` generates a **nonzero bounded ancient mild solution of the
unforced Navier-Stokes equations**. (The proviso is not decoration: (6), (7) and (H*) do not by
themselves kill the `b(t)` mode, and KNSS Remark 3.1 records that `w` and `b` are determined
only up to a constant. The two applications below supply it, from (6.4) and from (Dec).)

*Proof.* Set `h(t) = sup_x |u(x,t)|`, finite for `t < T` by (H), and `H(t) = sup_{s<=t} h(s)`.
By S3, `H(t) -> inf` as `t -> T`. Choose `t_k` increasing to `T` with `h(t_k) = H(t_k)`,
`gamma_k` decreasing to `1`, `N_k = H(t_k)`, and `x_k` with `M_k = |u(x_k,t_k)| >= N_k/gamma_k`;
then `M_k -> inf`. Set

    v^{(k)}(y,s) = (1/M_k) u( x_k + y/M_k , t_k + s/M_k^2 ),

defined for `s in (A_k, B_k)`, `A_k = -M_k^2 t_k -> -inf`. By S2 with `L = 1/M_k`, `v^{(k)}`
is a forced mild solution with force

    g^{(k)}(y,s) = M_k^{-3} g( x_k + y/M_k , t_k + s/M_k^2 ),
    || grad_y^j g^{(k)} ||_inf = M_k^{-3-j} || grad^j g ||_inf  ->  0,

so (4.4) holds. Also `|v^{(k)}| <= gamma_k` on `R^3 x (A_k, 0)` and `|v^{(k)}(0,0)| = 1`. That
`v^{(k)}` is a forced *mild* solution, rather than merely weak, on `(A_k, 0)` is the rescaling
of the corresponding statement for `u`, which is S5 plus S6: `u = v + G + w + b` and, in the
two applications below, the hypothesis forces `b = 0`. Apply S7 in form (a). The limit `v`
satisfies `|v| <= 1` on `R^3 x (-inf,0]` and `|v(0,0)| = 1`.

The value at the endpoint `s = 0` is the whole content of the contradiction below, so it is
worth saying where it comes from (referee finding R2, 2026-09-08). KNSS get it from Leray's
lower bound (6.1), which makes `B_k = M_k^2(T - t_k)` bounded below by a positive constant and
so puts `s = 0` in the interior of the interval on which the `v^{(k)}` live. That route is
unforced. Here it is supplied instead by the closed endpoint of S7(a), which rests only on
KNSS (3.10) on backward cylinders ending at the final time. (A second forced route: S4a gives a
local existence time from `t_k` bounded below by `c(gamma_1, ||g||_inf) M_k^{-2}`, hence
`B_k >= c > 0`, which reproduces the KNSS geometry with the force present.) `∎`

---

## 5. Proof of part (ii-a): the CSTY form is excluded

This is KNSS Theorem 6.1 with a force.

**Step S9 (CITED-VERBATIM).** *KNSS Theorem 5.3.* Let `v` be a bounded weak solution of the
**unforced** Navier-Stokes equations in `R^3 x (-inf, 0)`, axisymmetric, with
`|v(x,t)| <= C/sqrt(x_1^2+x_2^2)` there. Then `v = 0`. (Read from arXiv:0709.3599, Theorem
5.3, proof pages 16-18; the equation there is (1.1), `u_t + u grad u + grad p - Laplace u = 0`,
unforced.)

**Step S10 (PROVED-HERE).** *Forced KNSS Theorem 6.1.* Let `u` be as in the Theorem, assume
(H*), and suppose

    |u(x,t)| <= C / r     on R^3 x (0,T).                                (6.4)

Then (6.4) is impossible. (KNSS's Theorem 6.1 concludes the quantitative `|u| <= M(C)`; the
argument below does not, and does not need to. The earlier version of this line asserted
`sup |u| < inf`, which the proof never establishes; referee finding R8, 2026-09-08.)

*Proof.* Suppose it does. First, `u` is a forced mild solution on `R^3 x (0,T)` for a suitable
datum. Indeed, `u` is a bounded weak solution on each `R^3 x (0,T')` by (H), so S5 applies
with `F_k = -u_k u`, giving `u = v + G + w + b`. The bound (6.4) is inherited by `v` and by
`w` from the decay of the kernels `K_{ijk}` (KNSS (3.7)) and of the heat kernel, exactly as in
KNSS's proof of Theorem 6.1, and by `G` from S6. Hence `b`, which is constant in `x`, obeys
`|b(t)| <= C/r` for every `r`, so `b = 0`.

Now run S8. Since `|x'_k| <= C/M_k` by (6.4) (writing `x' = (x_1,x_2)`), the fields `v^{(k)}`
are axisymmetric about axes parallel to the `y_3` axis at distance at most `C` from it; pass
to a subsequence so that these axes converge, and the limit `v` is axisymmetric about a fixed
axis. The bound (6.4) is scale invariant, so `|v| <= C/r` in the limiting coordinates. By S7
the limit is an ancient mild solution of the **unforced** equations, hence in particular a
bounded ancient weak solution of the unforced equations, so S9 applies and `v = 0`. This
contradicts `|v(0,0)| = 1`. `∎`

**Step S11 (PROVED-HERE).** *(ii-a).* Suppose (T1-CSTY) held:
`|u(x,t)| <= C_*(r^2 + (T-t))^{-1/2}`. Since `r^2 + (T-t) >= r^2`, this gives
`|u(x,t)| <= C_*/r` on `R^3 x (0,T)`, which is (6.4). By S10 this is impossible. `∎`

*Remark.* CSTY part I, Theorem 1.1 (arXiv:math/0701796) states exactly the bound
`|v(x,t)| <= C_*(r^2 - t)^{-1/2}` on `D = R^3 x (-T_0, 0)`, singular time `0`; with the
singular time called `T` this is (T1-CSTY). Their hypotheses also include `p in L^{5/3}(D)`
and `v` smooth in `x`, Hölder in `t`. The route above does not use their proof, which is a
De Giorgi-Nash-Moser argument, and does not need their pressure hypothesis; it uses only
that their bound implies KNSS's (6.4). KNSS themselves note that their Theorems 6.1 and 6.2
"can be thought of as generalizations" of CSTY's results.

---

## 6. Proof of part (ii-b): the KNSS form is excluded under (Dec)

This is KNSS Theorem 6.2 with a force. The rescaling is done in two stages, and the force
behaves differently in the two stages; that is the only new bookkeeping.

**Step S12 (CITED-VERBATIM).** *KNSS Theorem 5.1 and Remark 6.1.* A bounded weak solution of
the unforced Navier-Stokes equations in `R^2 x (-inf,0)` has the form `u(x,t) = b(t)`; a
bounded ancient **mild** solution of that form is constant. (arXiv:0709.3599, Theorem 5.1
page 12 and Remark 6.1 page 18.)

**Step S13 (PROVED-HERE).** *Forced KNSS Theorem 6.2.* Let `u` be as in the Theorem, assume
(H*) and (Dec), and suppose

    |u| <= C / sqrt(T - t)     on R^3 x (0,T).                            (6.5)

Then (6.5) is impossible. (Same correction as S10; referee finding R8.)

*Proof.* By (Dec) and the argument of S10, `u` is a forced mild solution for a suitable
datum, hence smooth on `R^3 x (0,T)` by (S4c). Put `F(x,t) = r |u(x,t)|`. By S10 it suffices
to show `F` is bounded, since `F <= C` is (6.4).

Suppose not. Choose `t_k -> T` and `x_k` with
`M_k = F(x_k,t_k) = sup_{R^3} F(.,t_k) = sup_{s<=t_k} sup_{R^3} F(.,s) -> inf`. Put
`lambda_k = |x'_k|`. By (Dec), `lambda_k < R_0` for `k` large, since `r >= R_0` forces
`F <= C_0`. Define the first rescaling

    v^{(k)}(y,s) = lambda_k u( lambda_k y', lambda_k y_3 + x_{3k}, T + lambda_k^2 s ),

with `s_k = -(T-t_k) lambda_k^{-2}`. By S2 with `L = lambda_k`, `v^{(k)}` is a forced mild
solution with force `lambda_k^3 g(lambda_k y', lambda_k y_3 + x_{3k}, T + lambda_k^2 s)`,
whose sup norm is at most `R_0^3 ||g||_inf`. **This force is bounded but not small**, so no
limit is taken at this stage; the sequence `v^{(k)}` is only an intermediate object, exactly
as in KNSS. The scale-invariant bounds KNSS (6.9), (6.10), (6.11) are unaffected by the
force, being consequences of (6.5) and of the definition of `M_k`.

Define the second rescaling, as in KNSS (6.12),

    w^{(k)}(x,tau) = (1/M_k) v^{(k)}( e_1 + x/M_k , s_k + tau/M_k^2 ),  tau in (A_k, 0],

so that the total scaling factor applied to `u` is `lambda_k / M_k`. Then `w^{(k)}` is a
forced mild solution with force

    g_w^{(k)}, ||grad^j g_w^{(k)}||_inf = (lambda_k/M_k)^{3+j} ||grad^j g||_inf
                                        <= (R_0/M_k)^{3+j} ||grad^j g||_inf  ->  0,

because `lambda_k <= R_0` and `M_k -> inf`. So (4.4) holds and S7 applies, in form (b): the
`w^{(k)}` are **not** uniformly bounded up to `tau = 0`, since (6.16) gives only
`|w^{(k)}| <= C/sqrt(-tau)` globally and (6.14) gives `|w^{(k)}| <= 2` only off the cylinder
`C_k`, which escapes to infinity. A subsequence converges locally uniformly on the open set
`R^3 x (-inf,0)` to a bounded ancient mild solution `w` of the **unforced** equations, with
`|w| <= 2` by KNSS (6.14) and (6.16). Since the `v^{(k)}` are axisymmetric and `M_k -> inf`, the
zoom is centred on a point of the unit circle and the limit `w` is independent of `x_2`, so
`(w_1,w_3)` is a two-dimensional divergence-free field; by S12 it is a constant. **The constant
is zero** because (6.16) passes to the limit: `|w(x,tau)| <= C/sqrt(-tau) -> 0` as
`tau -> -inf`. (S12 by itself gives only "constant", not "zero"; referee finding R4,
2026-09-08.) With zero drift the remaining component solves the heat equation and inherits the
same decay, so `w = 0`, as in KNSS.

The final step is KNSS's: `|w^{(k)}(0,0)| = 1`, and one must show `w^{(k)}(0,0) -> w(0,0)`.
KNSS do this by applying the representation formula on `R^3 x (-1,0)` with `w^{(k)}(.,-1)` as
datum and showing that the contribution of the cylinder `C_k` where `|w^{(k)}|` may be large
is negligible, via their integral `I(M) -> 0`. With a force the representation formula
carries the extra term `int_{-1}^0 S(-s) g_w^{(k)}(s) ds`, whose sup norm is at most
`||g_w^{(k)}||_inf -> 0` by (4.2). It is therefore negligible, and the rest of KNSS's estimate
is unchanged. So `|w(0,0)| = 1`, contradicting `w = 0`. `∎`

**Step S14 (PROVED-HERE).** *(ii-b).* (T1-KNSS) is (6.5). By S13 it cannot hold, for any
`C_*`. By KNSS's definition, a singularity that is not type I is type II, and "not type I"
means precisely that `sup_x |u(x,t)| <= C/sqrt(T-t)` fails for every `C`, that is
`limsup_{t->T} sqrt(T-t) sup_x |u(x,t)| = inf`. `∎`

*Remark on CSTY II.* CSTY part II, Theorem 1.1 (arXiv:0709.4230) excludes both
`|v| <= C_*|t|^{-1/2}` (their (1.2), the KNSS form) and
`|v| <= C_* r^{-1+eps}|t|^{-eps/2}` for some `eps in [0,1]` (their (1.3)), under the
hypotheses `v^0 in H^{1/2}`, `r v_theta^0 in L^inf`, `p in L^{5/3}(D)`. Their (1.3) with
`eps = 0` is `|v| <= C_*/r`, which is S10 above; with `eps = 1` it is (1.2). The intermediate
range `0 < eps < 1` is **not** covered by the route of this note: it is not implied by (6.4)
and is not KNSS's hypothesis. Their proof of the intermediate range runs through De
Giorgi-Moser and Nash local regularity applied to `Gamma = r v_theta` together with local
energy estimates, and I did not extend it to a forced setting. See section 8, item F-5.

---

## 7. Part (i), part (iii), and `T^3`

### 7.1 Part (i): on the axis

**Step S15 (CITED-VERBATIM, from the companion note).** [FENCES] Lemma B, proved there with
the force present, states: for `alpha = 2`, an axisymmetric classical solution on
`[0,T) x R^3` with smooth axisymmetric divergence-free finite-energy data and smooth
axisymmetric force `f` satisfying (A1) and `f in L^q` for some `q > 5/2` with `div f = 0`
has `Sigma` contained in `{r = 0}`. The proof there is: rotation invariance of `Sigma`; the
Caffarelli-Kohn-Nirenberg epsilon-regularity criterion applied to `u` itself on backward
cylinders `Q_r(x_0,T) subset R^3 x [0,T)`; the Vitali covering argument giving
`P^1(Sigma x {T}) = 0`; and the observation that an off-axis point of `Sigma` would drag a
whole circle of positive `H^1` measure into `Sigma`.

**Step S16 (PROVED-HERE).** The force hypothesis of S15 is satisfied here: by (5) with
`K = 4`, `|f| <= C(1+|x|)^{-4}`, so `f in L^q(R^3 x [0,T])` for every `q in [1, inf]`; after the
S0 replacement the object actually used is `g = Pf`, which lies in `L^q` for `1 < q <= inf` by
(S1b) but **not** in `L^1`, and `q > 5/2` still holds with room to spare. (A1) is implied by (5),
as [FENCES] section 1.1 records. So part (i) holds. It does not use (H) itself and does not use
(Dec): the CKN criterion is applied on backward cylinders inside the region where `u` is already
classical, and the covering argument is local. It does, however, use (H*), through the one
global input of [FENCES] Lemma B step (ii): `u in L^inf_t L^2 ∩ L^2_t H^1` with
`int_0^T ||grad u||_{L^2}^2 dt < inf`, needed up to `T` because the cylinders `Q_r(x_0,T)` reach
`T`. The earlier version of this step said part (i) was free of every standing hypothesis;
that was wrong (referee finding R1, 2026-09-08). Nonemptiness of `Sigma` is a separate matter and
is the one place where a decay statement at spatial infinity enters, exactly as [FENCES] section
2.2 step (i) flagged. `∎`

### 7.2 Part (iii): swirl-free with a swirl-free force

**Step S17 (PROVED-HERE; machine-checked).** *The forced swirl equation.* In cylindrical
coordinates the `e_theta` component of forced (NS) reads, for `Gamma = r u_theta`,

    d_t Gamma + b . grad Gamma = nu ( Laplace Gamma - (2/r) d_r Gamma ) + r f_theta.

If `f_theta = 0` and `Gamma(.,0) = 0` then `Gamma = 0` for as long as the solution is
classical. The justification must not be phrased as "uniqueness for a linear drift-diffusion
equation with bounded coefficients": the coefficient `2 nu / r` is unbounded at the axis, and
`Gamma = r u_theta` grows linearly in `r` under (H), so a growth class would be needed as well
(referee finding R7, 2026-09-08). Argue on `u_theta` instead. Its equation is

    d_t u_theta + b . grad u_theta + (u_r/r) u_theta = nu ( Laplace - 1/r^2 ) u_theta + f_theta,

in which the potential `- nu/r^2` has the favourable sign for a maximum principle, `u_r/r` is
bounded under (H) because `u_r` vanishes on the axis and `u` is smooth with bounded second
derivatives, and `u_theta` is bounded on each slab by (H). With `f_theta = 0` and
`u_theta(.,0) = 0` the maximum principle gives
`||u_theta(t)||_inf <= e^{Ct} ||u_theta(0)||_inf = 0`. So the flow stays swirl free, and this
step uses (H), not only `f_theta = 0`. (Verified against Nowakowski-Zajaczkowski arXiv:2302.00730
equation (1.3), which gives the same equation with source `r f_phi`; conversely, a force with
`f_theta != 0` generates swirl out of swirl-free data, so this hypothesis is not removable.)

**Step S18 (PROVED-HERE; machine-checked).** *The forced vorticity equations, no swirl.* With
`omega = curl u = omega_theta e_theta` and `F_theta = (curl f) . e_theta`,

    d_t omega_theta + b . grad omega_theta
        = nu ( Laplace - 1/r^2 ) omega_theta + (u_r / r) omega_theta + F_theta,   (7.1)

and, setting `Omega = omega_theta / r`,

    d_t Omega + b . grad Omega = nu ( d_r^2 + (3/r) d_r + d_z^2 ) Omega + F_theta / r.  (7.2)

The operator in (7.2) is the five-dimensional Laplacian acting on `SO(4)`-invariant
functions, as in KNSS (5.15)-(5.16). Both identities are verified symbolically in
`check_forced_axisym.py`, checks C2 and C3, using a Stokes stream function so that
`div u = 0` holds identically and the pressure drops out of the poloidal curl. They agree
with the unforced forms recorded in Q. S. Zhang's survey (arXiv:2101.04905) equations (2.4),
(2.5), and with the forced `omega_phi` equation of Nowakowski-Zajaczkowski (1.9), whose
source is `F_phi = rot f . e_phi`.

Note that (7.2) has **no vortex stretching term**: the term `-2 v_r J` of Zhang (2.4)
carries `J = omega_r / r` and `omega_r = -d_z u_theta = 0` without swirl. The force enters
only as the source `F_theta/r`, and for smooth axisymmetric swirl-free `f` this source is a
smooth function across the axis (both `d_z f_r` and `d_r f_z` vanish to first order at
`r = 0`, exactly as `omega_theta` does), with

    || F_theta / r ||_{L^p(R^3)} <= C_p (1 + t)^{-K}   for every p in [1, inf], every K,  (7.3)

by (5): for `r >= 1` use `|F_theta|/r <= |grad f|`, and for `r <= 1` use
`|F_theta|/r <= sup_{r<=1} |d_r F_theta| <= C (1+|z|)^{-K}`.

**Step S19 (PROVED-HERE).** *A priori bounds.* Assume (H*) and `f_theta = 0`, `u_{0,theta} = 0`.
Items 1 and 3, and S20, are Grönwall arguments on `||u(t)||_{L^2}`, `||omega(t)||_{L^2}` and
`||Lambda^s u(t)||_{L^2}`; each presupposes that the quantity differentiated is finite and
continuous on `[0,T)`, which is (H*) and is **not** supplied by (H) (see section 1.1; (H) gives
`L^inf` bounds on derivatives, not `L^2` bounds). Then on `[0,T)`:

1. `||u(t)||_{L^2} <= ||u_0||_{L^2} + int_0^t ||f||_{L^2}`, bounded, by the forced energy
   estimate ([FENCES] Lemma A Step 1).
2. `||Omega(t)||_{L^inf} <= ||Omega(0)||_{L^inf} + int_0^t ||F_theta/r||_{L^inf} ds`,
   bounded by (7.3). This is the parabolic maximum principle applied to (7.2), read as an
   equation on `R^5 x (0,T')` with bounded drift, exactly as KNSS apply their Lemma 2.1 to
   their (5.16). It requires `Omega` bounded on each slab `R^3 x [0,T']`, which (H) supplies:
   `|omega_theta / r| <= C sup |grad^2 u|` near the axis, since `omega_theta` vanishes on the
   axis, and `|omega_theta / r| <= |grad u| / r` away from it.
3. `d/dt ||omega(t)||_{L^2} <= ||Omega||_{L^inf} ||u||_{L^2} + ||curl f||_{L^2}`, hence
   `||omega(t)||_{L^2}` is bounded on `[0,T)`.

   *Proof of 3.* Pair (7.1) with `omega_theta` over `R^3`. Since
   `Laplace(omega_theta e_theta) = ((Laplace - 1/r^2) omega_theta) e_theta`, the dissipative
   term is `-nu ||grad omega||_{L^2}^2 <= 0`. The drift term vanishes because `div b = 0`.
   The stretching term is
   `int (u_r/r) omega_theta^2 dx = int u_r Omega omega_theta dx <= ||Omega||_inf ||u||_{L^2} ||omega||_{L^2}`.
   The force term is `int F_theta omega_theta dx <= ||curl f||_{L^2} ||omega||_{L^2}`. `∎`

**Step S20 (PROVED-HERE).** *Closure.* From 1 and 3, `||u(t)||_{H^1}` is bounded on `[0,T)`,
hence so is `||u(t)||_{L^6}` by Sobolev. Now repeat [FENCES] Lemma A Step 2 at `alpha = 2`
with `L^6` in place of `L^inf`. Fix `s = 3`. The Kato-Ponce estimate of Grafakos-Oh quoted in
[FENCES], with `n = 3`, `r = 2`, `(p_1,q_1) = (6,3)`, `(p_2,q_2) = (3,6)` (both admissible:
`1/2 = 1/6 + 1/3`, and `sigma = s > max(0, n/r - n) = 0`), gives

    || Lambda^s (u (x) u) ||_{L^2} <= C || u ||_{L^6} || Lambda^s u ||_{L^3}.

With `X = ||Lambda^{s+1} u||_{L^2}`, `Y = ||Lambda^s u||_{L^2}`, Sobolev
`||Lambda^s u||_{L^3} <= C ||Lambda^{s+1/2} u||_{L^2}` (since `dot H^{1/2}(R^3) -> L^3`), and
Plancherel interpolation `||Lambda^{s+1/2}u||_{L^2} <= Y^{1/2} X^{1/2}`,

    |N| <= || Lambda^s (u (x) u) ||_{L^2} X <= C || u ||_{L^6} Y^{1/2} X^{3/2}
        <= (nu/2) X^2 + C(nu) || u ||_{L^6}^4 Y^2,

by Young with exponents `4/3` and `4`. Adding the force term
`<Lambda^s f, Lambda^s u> <= (1/2)||Lambda^s f||^2 + (1/2) Y^2` and applying Grönwall on
`[0,T)` with `||u||_{L^6}` bounded and `||f||_{H^s}` integrable by (5) gives
`sup_{t<T} ||u(t)||_{H^s} < inf`. By [FENCES] Lemma A Step 3 (Kato-type local existence in
`H^s`, `s > 5/2`, for the forced equation, with existence time bounded below in terms of
`||u(t_0)||_{H^s}` and the force), the solution continues past `T`. So `T` was not a singular
time. `∎`

This proves part (iii). Note where each hypothesis was used: `f_theta = 0` in S17, (H) in S17
and S19.2, the `H^s` half of (H*) in S19.1, S19.3 and S20 (the Grönwall arguments need the
quantities to be finite before they can be estimated), and (5) in S19 and S20.

**What the sources say.** The unforced statement is classical: O. A. Ladyzhenskaya, Zap.
Nauchn. Sem. LOMI **7** (1968), 155-177, and M. R. Ukhovskii, V. I. Yudovich, Prikl. Mat.
Mekh. **32** (1968), 59-69 (English: J. Appl. Math. Mech. **32** (1968), 52-61), with a
later and simpler proof by S. Leonardi, J. Málek, J. Nečas, M. Pokorný, Z. Anal. Anwendungen
**18** (1999), 639-649. I could **not** obtain the text of any of these three. What I could
verify at source is that the standard modern presentations state them **without a force**:
Q. S. Zhang's survey (arXiv:2101.04905) writes the axisymmetric system as its (1.5), which
has no force term, and says of that system that "if the swirl `v_theta = 0`, then it is known
for long time (Ladyzhenskaya, Ukhovskii and Yudovich) that finite energy solutions to (1.5)
are smooth for all time", citing Leonardi-Málek-Nečas-Pokorný as well; and Seregin-Šverák
(arXiv:0804.1803) Theorem 1.3, which they describe as "a local version of the corresponding
global regularity result in [11] and [30]" (Ladyzhenskaya; Ukhovskii-Yudovich), is likewise
stated for their (1.1), which is unforced. So: **the classical no-swirl theorem is stated
without a force in every source I could read, and the forced version is supplied here** as
S17-S20. The answer to the question posed in the brief is therefore yes, it survives, but
with two conditions that must be stated: the force itself must be swirl free, and on `R^3`
the maximum-principle step needs (H).

### 7.3 `T^3` and Fefferman (D)

**Step S21 (PROVED-HERE; machine-checked).** *Periodic axisymmetry is trivial.* Let
`u : R^3 -> R^3` be continuous, `Z^3`-periodic (Fefferman's (8): `u(x + e_j) = u(x)`), and
axisymmetric about the `x_3` axis. Then `u` does not depend on `x_1, x_2`. If moreover
`div u = 0`, then `u` is constant in space and of the form `(0, 0, c)`.

*Proof.* Let `Q` be the rotation by `alpha` about the `x_3` axis and, for `a in Z^3`, let
`Q_a(x) = Q(x - a) + a` be the rotation by `alpha` about the vertical line through `a`.
Periodicity gives `u(Q_a x) = u(Q(x-a)) = Q u(x - a) = Q u(x)`, so `u` is axisymmetric about
that line too. Then `Q_0^{-1} Q_a (x) = x + (Q^{-1} - I) a`, a pure translation (check C4 in
the script), and `u(Q_0^{-1} Q_a x) = Q^{-1} u(Q_a x) = Q^{-1} Q u(x) = u(x)`. So `u` is
invariant under translation by `(Q^{-1} - I)a` for every `alpha` and every `a in Z^3`. Taking
`a = e_1` and `a = e_2` and letting `alpha -> 0`, these vectors are
`(cos alpha - 1, -sin alpha) ~ (0, -alpha)` and `(sin alpha, cos alpha - 1) ~ (alpha, 0)`:
arbitrarily short, in two independent directions (check C4'). The group they generate is
dense in the horizontal plane, so by continuity `u` is independent of `x_1, x_2`. Axisymmetry
of a field depending only on `x_3` forces `u(x) = Q u(x)` for all `Q`, hence `u_r = u_theta = 0`,
so `u = (0,0,u_3(x_3))`; then `div u = d_3 u_3 = 0` gives `u_3` constant. `∎`

**Consequence for (D).** In Fefferman's periodic problem, the axisymmetric sub-class consists
of the fields `u(x,t) = (0,0,c(t))` with `c' = mean(f_3)`, which are global and smooth. (The
same argument applied to `f`, which uses only continuity, periodicity and axisymmetry, gives
`f = (0,0,f_3(x_3,t))`; periodicity of `p`, the Clay errata to (8), then forces
`c' = mean(f_3)`.) So the theorem above is **vacuous on `T^3`**: there is nothing in the
axisymmetric sub-class of (D) to constrain. The correct statement about (D) is not "the
argument transfers" but "the class is trivial".

**Scope of S21, stated precisely** (referee finding R10, 2026-09-08). The proof uses the **full
continuous** rotation group about a fixed vertical line, through the limit `alpha -> 0`. That is
the notion of axisymmetry the KNSS and CSTY theorems are about, and on `T^3` it is trivial. The
only rotations compatible with the lattice are the multiples of `pi/2`, and invariance under
that finite subgroup forces nothing: four-fold-symmetric periodic fields are a large class. So
S21 refutes the transfer of *this* notion of axisymmetry to (D), and says nothing about
lattice-symmetric constructions, which no theorem cited here covers either.

**What is and is not covered on `T^3`.** The brief's suggestion that "the rescaling is local
so the Liouville step runs" is half right and the half that fails is worth stating exactly.
Under the blow-up rescaling with `L -> 0` a periodic solution on `R^3/Z^3` becomes a solution
on `R^3/(L^{-1} Z^3)`, whose fundamental domain exhausts `R^3`, so the *limit object* does
live on `R^3` and the Liouville theorems S9, S12 would be available for it. What does not
transfer is the machinery of section 4: KNSS's sections 3 and 4 are built on the explicit
`R^3` Stokes kernels `K_{ij}, K_{ijk}` and the Duhamel formula (4.1) on `R^3`, which are
global objects, and the periodic Stokes kernel is a different function. (The earlier version of
this sentence also blamed the Liouville theorem for bounded curl-free divergence-free fields
used in KNSS Lemma 3.1. That is wrong: on `T^n` such a field is harmonic, hence constant by
compactness, so that input is easier on the torus, not harder. Referee finding R9, 2026-09-08.)
A `T^3` version would have to be built on a local theory instead. The natural tool is Seregin-Šverák, *On Type I singularities of the local
axi-symmetric solutions of the Navier-Stokes equations*, Comm. PDE **34** (2009), 171-201,
arXiv:0804.1803, whose Theorem 1.1 assumes only `v in L^3(Q(z_0,R))`, an associated pressure
`q in L^{3/2}(Q(z_0,R))`, axial symmetry, and the *local* Type I condition
`ess sup_{Q(z_0,R)} sqrt(t_0 - t) |v-bar(x,t)| < inf` on the poloidal part `v-bar`, and
concludes that `z_0` is a regular point. Their Theorem 1.2 does the same for the local
version of `r|v-bar| < inf`, and their Theorem 1.3 is the local no-swirl statement. Those
theorems would remove hypothesis (Dec) on `R^3` as well. **I did not extend Seregin-Šverák to
a forced setting.** Their route goes through suitable weak solutions, the local energy
inequality, and CKN-type epsilon-regularity, so the forced version would have to carry the
`2 (u . f) phi` term in the local energy inequality (as [FENCES] section 2.1 already does for
CKN), check that the scaled quantities (1.2)-(1.8) of their paper are still controlled with a
force, and re-run their Theorems 2.4 and 2.8. That is a separate piece of work and is listed
in section 8 as F-4.

---

## 8. What fails, and what is assumed

**F-1. A merely smooth and bounded force is not enough (genuine failure, and the fix).**
Decay does two separate jobs, and only the second is the operative one (referee finding R6,
2026-09-08).

*(a)* KNSS section 3 state that solutions of the Stokes system with a right-hand side in
`L^inf_{x,t}` not in divergence form are not well defined, because `P` is defined on
`L^inf(R^n)` only modulo constants. Decay makes `Pf` an honest function. This job exists only
when `div f != 0`.

*(b)* The step that fails for **every** bounded force, divergence free or not, is S6 and its use
in S10. Without `|G| <= C/r` the `b(t)` mode of the decomposition S5 is not killed, so `u` is
not shown to be a mild solution, so S8 does not run. This is exactly the failure KNSS flag after
their Theorem 6.2: "the statement fails, for trivial reasons, if we drop assumption (6.6).
(Consider `u(x,t) = b(t)`.)" A spatially constant divergence-free force generates precisely that
mode, and for it (a) is vacuous.

What (5) supplies is therefore: `g = Pf` in every `L^q`, `1 < q <= inf`, with pointwise decay
`(1+|x|)^{-3}` (S1), hence `|G| <= C/r` (S6). **Corrected statement:** the theorem holds for
forces in the Clay class (5); it does not hold as stated for forces that are merely smooth and
bounded on `R^3`, and this is true even for smooth bounded forces that are already divergence
free. This mirrors [FENCES] Lemma A's scope condition (A1), which likewise fails for a constant
force.

**F-2. Hypothesis (H*) is not supplied by the Clay class, and (H) alone is not enough.** (6)
and (7) give smoothness and bounded energy, not uniform bounds on slabs, so (H) is assumed. But
(H) is also not what the proofs need: S3, part (i), and the Grönwall arguments of S19 and S20 all
need `L^2`-based control, and an `L^inf` bound on every derivative does not give it (section 1.1
carries an explicit divergence-free field in `L^2` with all derivatives bounded and
`grad u` not in `L^2`). The hypothesis carried is (H*), `u in C([0,T');H^s)` for every `s`, the
solution class of [FENCES] Lemmas A and B. On `T^3` or a bounded domain compactness supplies the
`L^inf` half for free; on `R^3` neither half is free, and closing (H*) needs both a uniqueness
theorem in the class `C^inf ∩ L^inf_t L^2_x` (the gap [FENCES] section 3 records) and a
no-escape-to-spatial-infinity statement for the `H^s` norm. I did not close either. **Part (i)
is not free of (H*)**, contrary to the earlier version of this note; it is free of (H) and of
(Dec) but uses the energy consequence of (H*) (referee finding R1, 2026-09-08).

**F-3. Hypothesis (Dec) is KNSS's own and is not supplied here.** KNSS Theorem 6.2 carries
their (6.6) and remark that it holds for mild solutions with sufficiently fast decaying data,
citing Brandolese (2001, 2004) for the unforced problem. I did not check whether the forced
Duhamel term preserves the relevant decay class, and (Dec) is therefore a hypothesis of
(ii-b). Part (ii-a) does not need it.

**F-4. Seregin-Šverák's local theorems are not extended.** See section 7.3. Extending them
would (a) remove (Dec) from (ii-b), (b) give a `T^3` statement for *locally* axisymmetric
solutions, which is the only meaningful periodic statement given S21, and (c) require
carrying the force through a local energy inequality and epsilon-regularity, which is real
work, not bookkeeping.

**F-5. The intermediate CSTY II range is not covered.** CSTY II Theorem 1.1 also excludes
`|v| <= C_* r^{-1+eps} |t|^{-eps/2}` for `0 < eps < 1`. Neither KNSS Theorem 6.1 nor 6.2
covers that family, so the route of this note does not either. Extending it would require
force-extending CSTY's own De Giorgi-Moser and Nash arguments on `Gamma = r v_theta`, which
use local energy estimates; not attempted.

**F-6. Steps that use exact scaling invariance: none fail.** This is worth stating because it
was the anticipated failure mode. The blow-up rescaling is exactly covariant for the forced
system (S2, machine-checked), with the force picking up `L^3` rather than `L^1`. Because
`L -> 0`, the force is *strongly subcritical*: it does not merely survive the rescaling, it
vanishes in the limit together with all its derivatives. The Liouville theorems S9 and S12
are therefore applied to a genuinely unforced ancient solution and are used verbatim, with no
hypothesis of theirs modified. Two places where the force does not become small: the
intermediate rescaling `v^{(k)}` in S13, where no limit is taken; and the `b = 0` step of S10 and
S13, where the force enters through `G` at a fixed scale and its decay, not its smallness, is
what is used (F-1(b); referee finding R6, 2026-09-08).

**F-7. Steps that use a global energy identity: none, on the KNSS route.** The KNSS route
uses no energy inequality at all, local or global; it is built on the Duhamel formula and
maximum principles for scalar quantities. The energy inequality enters only through part (i)
(the CKN local energy inequality, which [FENCES] already carries with its `2(u.f)phi` term)
and through part (iii) step S19.

**F-8. Leray's lower bound is unforced and is not used, but the reason matters.** KNSS quote
`sup_x|u(x,t)| >= eps_1/sqrt(T-t)` before their construction. It is an unforced result and I did
not extend it. The earlier version of this item said the construction "needs only `M_k -> inf`".
That is not an accurate account of KNSS: `M_k -> inf` is not what they use (6.1) for. They use it
to make `B_k = M_k^2(T - t_k)` bounded below, which puts `s = 0` in the *interior* of the time
interval on which the `v^{(k)}` are defined, and `|v(0,0)| = 1` is the entire content of the
contradiction in S10 and S13. Here the endpoint is closed instead by S7(a), whose equicontinuity
reaches `t = 0` because KNSS (3.10) holds on backward cylinders ending at the final time; and a
second, purely forced route is S4a, whose local existence time gives
`B_k >= c(gamma_1, ||g||_inf) > 0` directly. Either way Leray's bound is not needed, but it had
to be replaced, not merely dropped (referee finding R2, 2026-09-08).

---

## 9. What this says about the Clay problem

Inside the axisymmetric sub-class of Fefferman (C), with force in the class (5) and under the
two standing hypotheses (H*) and (Dec), a breakdown at a finite time `T` must satisfy all of:

1. `sup_{t<T} ||u(t)||_{L^inf} = inf` ([FENCES] Lemma A, given (5) and (H*); *not*
   unconditional, see F-2);
2. every blow-up point lies on the axis `{r = 0}` (part (i), given (5) and the energy
   consequence of (H*); *not* unconditional, see S16 and F-2);
3. `limsup_{t->T} sqrt(T-t) sup_x |u(x,t)| = inf`, that is, blow-up strictly faster than the
   self-similar rate (part (ii-b) under (Dec); part (ii-a) unconditionally excludes the
   stronger `(r^2 + (T-t))^{-1/2}` form);
4. nonzero swirl in the data or in the force (part (iii)).

Item 4 is the one that most directly touches the reported Alpöge-Buckmaster mechanism, which
is (as reported, and unverified, see [FENCES] section 0) axisymmetric **with** swirl, so it
passes item 4. That result is for Euler, `nu = 0`, which is outside the range of every fence
here; what [FENCES] Corollary A says is that the mechanism, as reported, keeps the velocity
bounded and therefore cannot be transplanted to any `alpha > 1`, in particular not to
`alpha = 2`, where item 1 would be violated (referee finding R11, 2026-09-08). Items 2 and 3 are
new fences for any successor construction that repairs the velocity bound.

The remaining corridor for an axisymmetric forced blowup is narrow on both sides: from below
by Leray's `sup_x|u| >= eps_1/sqrt(T-t)` (unforced, not re-proved), from above by item 3, and
in geometry by item 2. Nothing here says the corridor is empty. FL-000 still stands.

---

## 10. Verification record

Everything in this table was fetched and read during the session that produced this note.
PDFs were downloaded from the URLs shown and converted with `pdftotext -layout`; line numbers
refer to that conversion.

| Claim | Source, read at | Status |
|---|---|---|
| KNSS system (1.1) is unforced; scaling `u -> L u(Lx, L^2 t)` | arXiv:0709.3599 PDF, section 1 | Verified verbatim |
| KNSS section 3: `P` on `L^inf` only modulo constants, so non-divergence-form `L^inf` right-hand sides give ill-defined solutions | arXiv:0709.3599, section 3, paragraph after (3.7) | Verified verbatim; this is the source of F-1 |
| KNSS mild solution definition (3.5); kernel bounds (3.6), (3.7); estimates (3.10)-(3.14) | arXiv:0709.3599, section 3 | Verified verbatim |
| KNSS Lemma 3.1 (decomposition `u = v + w + b`) and Remark 3.1 | arXiv:0709.3599, section 3, end | Verified verbatim |
| KNSS bilinear estimate (4.5); Proposition 4.1; Lemma 4.1; vorticity bootstrap (4.7)-(4.11) | arXiv:0709.3599, section 4 | Verified verbatim |
| KNSS Theorem 5.1 (2D), Theorem 5.2 (no swirl), Theorem 5.3 (`|u| <= C/r`) | arXiv:0709.3599, section 5 | Verified verbatim, with proofs read |
| KNSS equations (5.12), (5.14), (5.15), (5.16); Remark 5.1 (`omega_theta/r` smooth across the axis) | arXiv:0709.3599, section 5 | Verified verbatim |
| KNSS Lemma 6.1, Remark 6.1, Proposition 6.1 and the rescaling (6.2), (6.3) | arXiv:0709.3599, section 6 | Verified verbatim |
| KNSS Theorem 6.1 (`|u| <= C/r`) and Theorem 6.2 (Type I plus (6.6)), with both proofs including the two-step rescaling (6.8), (6.12) and the final `I(M) -> 0` estimate | arXiv:0709.3599, section 6 | Verified verbatim, proofs read line by line |
| KNSS Type I definition; "type II is any singularity that is not type I" | arXiv:0709.3599, section 1 | Verified verbatim |
| KNSS citation of Brandolese for (6.6) under fast-decaying data | arXiv:0709.3599, section 6, remark after Theorem 6.2, and refs [1], [2] | Verified verbatim |
| CSTY I Theorem 1.1: `\|v(x,t)\| <= C_*(r^2 - t)^{-1/2}` on `D = R^3 x (-T_0,0)`, `p in L^{5/3}(D)`, conclusion `v in L^inf(B_R x [-T_0,0])`; system (N-S) unforced | arXiv:math/0701796 PDF, section 1 | Verified verbatim |
| CSTY II Theorem 1.1: hypotheses `v^0 in H^{1/2}`, `r v_theta^0 in L^inf`, `p in L^{5/3}(D)`, bounds (1.2) `\|v\| <= C_*\|t\|^{-1/2}` and (1.3) `\|v\| <= C_* r^{-1+eps}\|t\|^{-eps/2}`, `eps in [0,1]`; proof route is De Giorgi-Moser and Nash, not Liouville; system unforced | arXiv:0709.4230 PDF, section 1 and section 2 opening | Verified verbatim |
| CSTY II acknowledge KNSS proved "results similar to Theorem 1.1 using a different approach based on Liouville theorems" | arXiv:0709.4230, section 1, last paragraph | Verified verbatim |
| Seregin-Šverák local Theorems 1.1, 1.2, 1.3; local Type I conditions (1.9), (1.10) on the poloidal part; system (1.1) unforced; Theorem 1.3 described as the local version of Ladyzhenskaya and Ukhovskii-Yudovich | arXiv:0804.1803 PDF, section 1 and refs [11], [16], [30] | Verified verbatim |
| Axisymmetric system stated **without a force** in the standard survey; no-swirl regularity attributed to Ladyzhenskaya, Ukhovskii-Yudovich, Leonardi-Málek-Nečas-Pokorný; unforced `Omega = omega_theta/r` equation (2.4) with stretching term `-2 v_r J`, `J = omega_r/r`; unforced `omega_theta` equation (2.5); summary of the KNSS proof route | Q. S. Zhang, arXiv:2101.04905 PDF, sections 1 and 2 | Verified verbatim |
| Forced axisymmetric equations: swirl equation with source `r f_phi` (1.3); `omega_phi` equation with source `F_phi = rot f . e_phi` (1.9) | Nowakowski-Zajaczkowski, arXiv:2302.00730 PDF, section 1 | Verified verbatim. Note: their domain is a bounded cylinder, so this is used only to confirm the form of the source terms, not as a theorem |
| Leonardi-Málek-Nečas-Pokorný, Z. Anal. Anwendungen **18** (1999), 639-649, title, authors, journal data, abstract | ems.press/journals/zaa/articles/11075 | Bibliographic data and abstract verified; **full text not obtained**, so whether their system carries a force is unverified |
| Ladyzhenskaya, Zap. Nauchn. Sem. LOMI **7** (1968), 155-177; Ukhovskii-Yudovich, Prikl. Mat. Mekh. **32** (1968), 59-69 / J. Appl. Math. Mech. **32** (1968), 52-61 | Bibliographies of arXiv:0709.3599 [18], [29], arXiv:0804.1803 [11], [30], arXiv:2101.04905 [58], [112] | Bibliographic data cross-checked across three independent bibliographies; **neither paper obtained** |
| Clay conditions (4)-(11) and statements (C), (D) | [FENCES] section 3, which records them as fetched from claymath.org `navierstokes.pdf` and quoted verbatim | Relied on, not re-fetched |
| [FENCES] Lemma A (`alpha = 2`, forced) and Lemma B (forced) | `FENCES_forced_NS_blowup.md` in this directory, sections 1 and 2, refereed with 20 corrections applied | Relied on as a companion result |
| Alpöge-Buckmaster preprint | not found | Same status as in [FENCES]: reported, never cited |

**Machine checks.** `check_forced_axisym.py` in this directory, `python3 check_forced_axisym.py`,
all seven checks pass. It verifies: C1, exact scaling covariance of the forced system under
`(u,p,f) -> (L u, L^2 p, L^3 f)` with a control showing the exponent `2` fails; C2a, the
stream-function ansatz is divergence free; C2, the poloidal curl of the forced momentum
equations equals the forced `omega_theta` equation (7.1); C3, the change of variable
`omega_theta = r Omega` produces (7.2) with the five-dimensional Laplacian and source
`F_theta/r`; C4 and C4', the translation identity and density used in S21. Hashes in
`SHA256SUMS`.

---

## 11. New versus assembled

**Assembled.** The Liouville theorems (KNSS 5.1, 5.2, 5.3), the blow-up rescaling and the two
axisymmetric applications (KNSS 6.1, 6.2), the mild-solution framework (KNSS 3, 4), the
Type I/II vocabulary, the CKN partial regularity input, the classical swirl-free regularity
theorem, and Fefferman's class conditions are all taken from the sources named. [FENCES]
Lemma A and Lemma B, both already proved with a force, are used as given.

**New here.** Four things. First, the forced statement itself: no forced version of the
axisymmetric Type I exclusion appears in any source I could read, and [FENCES] section 2.4
explicitly declined to assert one; sections 4 to 6 supply it, with the extra Duhamel term
carried through every step of KNSS sections 3, 4 and 6. Second, the observation that the
force enters the *Liouville* step only through terms that vanish in the rescaling limit, so
that the Liouville theorems are applied to an unforced ancient solution and need no extension
at all; this is why the exclusion is *unconditional in the force* once (5) is assumed, rather
than requiring smallness or structure. (It is not true that the force enters *only* through
vanishing terms: it also enters at fixed scale through `G` in the `b = 0` step, where decay
rather than smallness is what is used. See F-1(b) and F-6.) Third, the identification of which force hypothesis is
load-bearing: not smoothness and not boundedness but decay. The reason is not primarily the
`L^inf`-to-`BMO` failure of `P` (which is vacuous for a force that is already divergence free)
but the `b(t)` mode: without `|G| <= C/r` the decomposition S5 is not shown to have `b = 0`, so
`u` is not shown to be a mild solution. A smooth bounded force on `R^3` is not covered even when
it is divergence free, and Fefferman's (5) is exactly the repair (F-1, S1, S6). Fourth, the observation
that the axisymmetric sub-class of Fefferman (D) contains only spatially constant vertical
fields (S21), so that no axisymmetric statement about the periodic Clay problem has content;
this replaces the expected "the argument transfers to `T^3`" with a sharper negative.

The forced swirl-free result (S17-S20) sits between the two categories: the unforced theorem
is classical, and the extension is routine in structure, but every source I could read states
it without a force, so the proof is written out here, and it turns up two conditions worth
stating, namely that the force must itself be swirl free and that the maximum-principle step
needs (H) on `R^3`.

---

## 12. Step-label index

| Step | Content | Label |
|---|---|---|
| S0 | Leray projection of the force; pressure absorbs the gradient part | PROVED-HERE |
| S1 | `g = Pf` decays like `(1+\|x\|)^{-3}`, lies in every `L^q`, `q > 1` | PROVED-HERE |
| S2 | Scaling covariance of forced (NS); `\|\|grad^j f_L\|\|_inf = L^{3+j} \|\|grad^j f\|\|_inf` | PROVED-HERE (machine-checked) |
| S3 | `sup_{t<T} \|\|u(t)\|\|_inf = inf` at a first singular time | [FENCES] Lemma A, forced, **plus (H\*)**: Lemma A is stated for the strong class, so this is not CITED-VERBATIM on its own |
| S4 | Forced mild solutions: definition, estimates (4.2)-(4.3), local existence, blow-up criterion, smoothing | ROUTINE-EXTENSION (force term is an additive non-singular Duhamel term; bilinear estimate untouched) |
| S5 | Forced decomposition `u = v + G + w + b` | PROVED-HERE (reduces to KNSS Lemma 3.1 for `u - G`) |
| S6 | `\|G\| <= C(1+\|x\|)^{-3}`, hence `<= C/r` | PROVED-HERE |
| S7 | Compactness with vanishing force; the ancient limit is unforced; two forms, closed and open endpoint | PROVED-HERE |
| S8 | Forced Proposition 6.1: a singularity generates a nonzero bounded unforced ancient mild solution | PROVED-HERE |
| S9 | KNSS Theorem 5.3, applied to the unforced limit | CITED-VERBATIM |
| S10 | Forced KNSS Theorem 6.1: `\|u\| <= C/r` excluded | PROVED-HERE |
| S11 | Part (ii-a): the CSTY form is excluded | PROVED-HERE |
| S12 | KNSS Theorem 5.1 and Remark 6.1, applied to the unforced limit | CITED-VERBATIM |
| S13 | Forced KNSS Theorem 6.2: `\|u\| <= C/sqrt(T-t)` excluded under (Dec) | PROVED-HERE |
| S14 | Part (ii-b): Type II in the KNSS sense | PROVED-HERE |
| S15 | On-axis localisation with force | CITED-VERBATIM ([FENCES] Lemma B, forced; its step (ii) needs the energy consequence of (H\*)) |
| S16 | (5) supplies `f in L^q`, `q > 5/2`, `div f = 0` | PROVED-HERE |
| S17 | Forced swirl equation; `f_theta = 0` preserves swirl-freeness | PROVED-HERE (machine-checked; maximum principle on `u_theta`, uses (H)) |
| S18 | Forced `omega_theta` and `Omega = omega_theta/r` equations, no swirl; no stretching term | PROVED-HERE (machine-checked) |
| S19 | A priori bounds: energy, `\|\|Omega\|\|_inf` by maximum principle, `\|\|omega\|\|_{L^2}` | PROVED-HERE (uses (H) and the `H^s` half of (H\*)) |
| S20 | `H^1` bound closes to `H^s` via Kato-Ponce at `(6,3)`; part (iii) | PROVED-HERE |
| S21 | Periodic axisymmetric divergence-free fields are constant vertical fields | PROVED-HERE (machine-checked) |

Steps that failed or were not attempted: F-1 (a merely bounded smooth force is not covered,
even when divergence free; corrected statement given), F-4 (Seregin-Šverák not force-extended,
so (Dec) is not removed and `T^3` gets no local statement), F-5 (the intermediate CSTY II range
`0 < eps < 1` is not covered), plus the two carried hypotheses (H*) and (Dec) recorded as F-2
and F-3.

---

## 13. Referee corrections (2026-09-08)

Single refuter pass. Full report: `REFUTE_forced_axisymmetric_on_axis_typeII.md` in this
directory. The script was re-run from a copy (seven PASS, exit 0). KNSS arXiv:0709.3599,
CSTY I arXiv:math/0701796, CSTY II arXiv:0709.4230 and Seregin-Šverák arXiv:0804.1803 were
re-fetched and read; every CITED-VERBATIM statement drawn from them was checked at source and is
accurate. Verdict on the note as it stood: **refuted as stated**, on R1; the conclusions survive
under the corrected hypothesis.

**Applied.**

1. **R1 (MAJOR).** *The standing hypothesis was understated.* (H) is not the only thing (4)-(7)
   fail to supply, and part (i) is not free of it. S3 rests on [FENCES] Lemma A, stated for the
   strong class; part (i) rests on Lemma B step (ii), which uses
   `int_0^T ||grad u||_{L^2}^2 dt < inf` up to `T`; S19.3 and S20 are Grönwall arguments needing
   `||omega(t)||_{L^2}` and `||u(t)||_{H^s}` finite a priori. (H) gives none of these: an explicit
   divergence-free field with all derivatives bounded, in `L^2`, with `grad u` not in `L^2`,
   separates the classes. Applied: section 1.1 rewritten around
   `(H*) u in C([0,T');H^s)` for every `s`, with the separating field written out and with the
   two reasons (H*) cannot itself be discharged; theorem statement now assumes (H*) throughout,
   with part (i) attributed to its energy consequence; S3, S16, S19 annotated; F-2 rewritten;
   section 9 items 1 and 2 no longer say "unconditional"; step-index rows for S3, S15, S19
   corrected.
2. **R2 (MINOR, load-bearing).** *`|v(0,0)| = 1` was not delivered by S7 as stated.* S7 concluded
   convergence on the open set `R^3 x (-inf,0)`; S8 used the value at `s = 0`. KNSS obtain the
   endpoint from Leray's (6.1), which bounds `B_k = M_k^2(T-t_k)` below, and F-8 declared that
   bound unnecessary without replacing its actual function. Applied: S7 restated with a closed
   endpoint in form (a), its proof annotated with the reason (KNSS (3.10) holds on backward
   cylinders ending at the final time), S8 given the justifying paragraph plus the alternative
   forced route through S4a, and F-8 rewritten.
3. **R3 (MINOR).** *S7's uniform bound is not satisfied in its second application.* `w^{(k)}` in
   S13 obeys only `|w^{(k)}| <= C/sqrt(-tau)` globally and `<= 2` off the escaping cylinder `C_k`.
   Applied: S7 split into forms (a) and (b), and S13 now invokes form (b) explicitly.
4. **R4 (MINOR).** *`(w_1,w_3) = 0` did not follow from S12.* S12 gives "constant", not "zero".
   Applied: S13 now kills the constant with KNSS (6.16), `|w| <= C/sqrt(-tau) -> 0` as
   `tau -> -inf`, and spells out why `w = 0` follows.
5. **R5 (MINOR).** *S1's proof mishandled the Riesz kernel.* The outer region of the split
   contains the singularity `y = x`, and `k`, homogeneous of degree `-3` on `R^3`, is not locally
   integrable. Applied: the cutoff-at-`|y-x| <= |x|/4` argument written in, with the
   `||R_iR_j h||_inf <= C||h||_{W^{4,1}}` bound; the `q > 1` sharpness noted; S16 corrected to
   attribute `L^1` to `f` and `L^q`, `q > 1`, to `g`; S0 now says `g` does not itself satisfy (5).
6. **R6 (MINOR).** *F-1 and section 11 misidentified which property of (5) is load-bearing.* A
   smooth bounded force that is already divergence free needs no Leray projection, so the
   `L^inf`-to-`BMO` objection is vacuous for it, yet the theorem still fails, at S6 and the
   `b = 0` step. Applied: F-1 split into (a) and (b) with the KNSS `u = b(t)` remark quoted; F-6
   now lists the `b = 0` step as a second place where the force does not vanish; section 11's
   second and third "new" items corrected.
7. **R7 (MINOR).** *S17's justification was wrong.* "Bounded coefficients" is false (`2 nu / r`)
   and `Gamma = r u_theta` is unbounded on `R^3`. Applied: the argument replaced by a maximum
   principle on `u_theta`, and the record that S17 uses (H); section 7.2's closing attribution
   corrected.
8. **R8 (MINOR).** *S10 and S13 stated a conclusion their proofs do not establish.* Neither proof
   bounds `sup |u|`; each shows directly that the hypothesised bound is impossible. Applied: both
   statements restated, with a note that KNSS's quantitative `|u| <= M(C)` is not reproduced.
9. **R9 (MINOR).** *Section 7.3's Liouville clause was wrong.* Bounded curl-free divergence-free
   fields on `T^n` are constant by compactness, so that input transfers trivially; the genuine
   obstruction is the explicit `R^3` Stokes kernels. Applied: the clause deleted with the reason
   recorded.
10. **R10 (MINOR).** *S21 is correct but its scope was not stated.* The proof needs the full
    continuous rotation group; the lattice-compatible subgroup (multiples of `pi/2`) forces
    nothing. Applied: a "Scope of S21" paragraph added, the `f` and pressure bookkeeping for
    `c' = mean(f_3)` written out, and "the class is empty" corrected to "the class is trivial".
11. **R11 (MINOR).** *Section 9's Alpöge-Buckmaster sentence was flatter than the evidence.* The
    reported result is for Euler, outside the range of every fence here. Applied: rewritten as a
    statement about transplanting the mechanism to `alpha > 1`.

**Not applied, and why.**

* The referee's observation that S13 uses (Dec) only through `lambda_k = |x'_k|` bounded, so that
  the weaker `limsup_{r -> inf} sup_{t<T} r|u| < inf` would suffice, is correct but is a
  strengthening rather than a correction; F-3 already carries (Dec) honestly and is left as is.
* The referee's note that the `C3` docstring in `check_forced_axisym.py` reads as though the
  source term `F_theta/r` is verified when it is only implied (the check verifies the operator
  identity `L_omega[r Om] = r L_Om[Om]`, from which the source divides through) is accurate but
  is not an error in the note or in the check. Left as is.
* No change to the section 10 verification table: every row re-checkable at source was re-checked
  and is accurate. One addition the referee flagged and I record here rather than in the table:
  Seregin-Šverák Theorems 1.2 and 1.3, unlike their Theorem 1.1, also assume `v` essentially
  bounded on `B(x_0,R) x (t_0-R^2, t')` for each `t' < t_0`, a local form of (H). The Theorem 1.1
  route named in section 7.3 as the way to remove (Dec) does not carry that hypothesis, so the
  claim there stands.

**Not touched:** everything outside the passages named above.

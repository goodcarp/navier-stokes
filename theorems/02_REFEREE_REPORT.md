# Referee report on `THEOREM_forced_axisymmetric_on_axis_typeII.md`

Date: 2026-09-08. Single refuter seat, HOLD UP SHIPS campaign.
Target: `THEOREM_forced_axisymmetric_on_axis_typeII.md` (809 lines) and `check_forced_axisym.py`,
in `Solve Navier Stokes/campaign/external/alpoge-buckmaster-2026-09-08/`.
Companion read: `FENCES_forced_NS_blowup.md` (Lemmas A and B, Clay conditions quoted in its section 3).

**Verdict: REFUTED AS STATED, one MAJOR and ten MINOR findings, no FATAL.**

The mathematical content survives. What does not survive is the hypothesis bookkeeping. The note's
headline claim, that `(H)` is "the one standing hypothesis that (4)-(7) do not supply" and that part
(i) and items 1 and 2 of section 9 are "unconditional given (5)", is false. Every part of the theorem,
including part (i), needs the solution to lie in the strong/energy class that `[FENCES]` Lemma A
already assumes, and that class is strictly stronger than (4)-(7) plus `(H)`. With `(H)` replaced by
`(H*)` (below) the theorem is proved, modulo the ten MINOR repairs, all of which are one-sentence
fixes and none of which changes a conclusion.

The script was re-run from a copy in a scratch directory: all seven checks PASS in 7.3 s. The
symbolic content of C1, C2a, C2, C3, C4, C4' was read line by line and is correct; see finding N-0.

Sources fetched at source and read: KNSS arXiv:0709.3599 (PDF, `pdftotext -layout`, 26 pp),
CSTY I arXiv:math/0701796, CSTY II arXiv:0709.4230, Seregin-Sverak arXiv:0804.1803. Every
CITED-VERBATIM label in the note that touches those four papers was checked against the source and
is accurate; see section "Citations that check out".

---

## Findings

### R1 (MAJOR). The standing hypothesis is understated: every part needs the strong/energy class, not `(H)`

**What the note says.** Section 1.1: "`(H)` ... is the one standing hypothesis that (4)-(7) do not
supply". Theorem part (i): "This uses neither `(H)` nor `(Dec)`". Section 8 F-2: "Part (i) is free of
it." Section 9: item 1 "unconditional given (5)", item 2 "unconditional given (5)".

**Why it is wrong.** Three independent places need `u(t)` to lie in an `L^2`-based space, and `(H)`
supplies only `L^inf` bounds.

1. **S3.** S3 is labelled CITED-VERBATIM from `[FENCES]` Lemma A. But Lemma A is stated for the
   *strong solution class* `u in C([0,T); H^s) ∩ L^2_loc([0,T); H^{s+alpha/2})`, and Lemma A's own
   parenthetical says so explicitly: "If the lemma is to be applied to a Clay-class solution, `C^inf`
   with bounded energy, one first identifies it with this strong solution on `[0,T)` ... see the
   uniqueness caveat at the end of section 3." `[FENCES]` section 3 then records that the needed
   uniqueness theorem in the class `C^inf ∩ L^inf_t L^2_x` is not available, and that Prodi-Serrin
   needs `grad u in L^2_{t,x}` and the energy inequality, "which (6), (7) do not give for free". The
   theorem note drops that caveat entirely. So S3 is not CITED-VERBATIM; it is Lemma A plus an
   unstated hypothesis.

2. **Part (i).** `[FENCES]` Lemma B step (ii) uses `u in L^inf_t L^2 ∩ L^2_t H^1` "by Step 1 of
   Lemma A (the energy inequality bounds `nu int_0^T ||grad u||_2^2 dt`)", and the Vitali covering is
   run "with `int_0^T ||grad u||_2^2 dt < inf`". The backward cylinders `Q_r(x_0,T)` reach up to `T`,
   so the dissipation bound is needed *up to* `T`, and `(H)`, which is a statement about slabs
   `[0,T']` with `T' < T`, does not give it. Part (i) is therefore not free of the class hypothesis.

3. **Part (iii), S19.3 and S20.** S19.3 writes `d/dt ||omega(t)||_{L^2} <= ...` and S20 runs Gronwall
   on `Y = ||Lambda^s u||_{L^2}`. Both presuppose that the quantity being differentiated is finite,
   that is `u(t) in H^{s+1}` on `[0,T)`. `(H)` does not give even `u(t) in H^1`.

**Explicit separation of the two classes.** `(H)` plus (7) does not imply `grad u(t) in L^2`. Let
`psi in C_c^inf(B_1)`, pick points `x_j` with `|x_j| -> inf` and pairwise distance `> 4`, and set
`w_j = 1/log j`, `a_j = ((log j)/j)^{1/2}`, `A = sum_j a_j w_j psi((x - x_j)/w_j) e`, `u = curl A`.
Then `u` is smooth and divergence free, `|grad^k u| ~ a_j w_j^{-k} = ((log j)^{2k+1}/j)^{1/2} -> 0`
so every derivative is uniformly bounded; `||u||_{L^2}^2 ~ sum a_j^2 w_j^3 = sum 1/(j (log j)^2) < inf`;
and `||grad u||_{L^2}^2 ~ sum a_j^2 w_j = sum 1/j = inf`. So a field satisfying (7) and the `t`-slice
of `(H)` for every `k` need not lie in `H^1`. Whether the *Navier-Stokes evolution from Schwartz data*
can reach such a state is open, which is exactly why it has to be assumed rather than asserted.

**Secondary point.** The decomposition S5 leaves a mode `b(t)` that is constant in `x`. The note kills
it only under (6.4) or `(Dec)` (S10, S13). Neither (6), (7) nor `(H)` kills it, and KNSS Remark 3.1
records that `w` and `b` are determined only up to a constant. So "`u` is a mild solution" is not
available from the Clay conditions either.

**Corrected hypothesis.** Replace `(H)` by

    (H*)  u in C([0,T'); H^s(R^3))  for every s >= 0 and every T' <= T.

`(H*)` implies `(H)` by Sobolev embedding (`H^{k+2} -> W^{k,inf}` on `R^3`), implies
`int_0^T ||grad u||_{L^2}^2 dt < inf` through the energy identity (Lemma A Step 1), and makes the
Gronwall arguments of S19.3 and S20 legitimate. It is exactly the class `[FENCES]` Lemma A and Lemma B
already work in. It is **not** implied by (4)-(7).

**Answer to the attack "is `(H)` automatic?" (the requested two-line argument, and why it fails).**
The attempted argument is: `u_0 in ∩_s H^s` by (4) and `f in L^1_loc([0,inf); H^s)` by (5), so Kato's
local theory gives a unique maximal strong solution `u~ in C([0,T_*); H^s)` for every `s`, and Sobolev
gives `(H)` on every `[0,T']`, `T' < T_*`. It fails at two points, neither of them cosmetic.

* It produces `(H)` for `u~`, not for the given Clay solution `(u,p)`. The identification needs
  uniqueness in the class `C^inf ∩ L^inf_t L^2_x` of (6), (7). The energy proof of uniqueness needs
  the difference to lie in `L^2_t H^1` (see the separation above), needs `u in L^inf` to kill the
  boundary term `int_{|x|=R} u_n |w|^2`, and needs a pressure normalisation to kill
  `int_{|x|=R} (p - p~) w_n`, since Fefferman's (6) imposes no decay on `p`. `[FENCES]` section 3
  records precisely this gap and refuses to close it.
* Even granting the identification, it gives `(H)` only on `[0, min(T,T_*))`. Ruling out `T_* < T`
  means excluding `||u(t)||_{H^3} -> inf` as `t` rises to `T_*` while `u` stays `C^inf` on
  `R^3 x [0,T)`, that is, `H^3` mass escaping to spatial infinity. Smoothness plus (7) does not
  forbid that; the field above shows the two classes differ.

So `(H)` cannot be discharged, and neither can `(H*)`. Both must be carried. The note is right to
carry a hypothesis; it is wrong about which one, and wrong that part (i) escapes it.

**Severity: MAJOR.** The theorem as stated is not established. With `(H*)` it is.

---

### R2 (MINOR, load-bearing). S8's `|v(0,0)| = 1` is not delivered by S7 as stated, and F-8 is right for the wrong reason

S7 concludes convergence "locally uniformly on `R^3 x (-inf, 0)`", an open set. S8 then asserts
"The limit `v` satisfies `|v| <= 1` on `R^3 x (-inf,0)` and `|v(0,0)| = 1`", using the value at the
endpoint `s = 0`. Locally uniform convergence on an open set says nothing at its boundary, and
`|v(0,0)| = 1` is the entire content of the contradiction in S10 and S13.

KNSS do not have this problem, because they use Leray's lower bound (6.1) exactly here: it gives
`B_k = M_k^2 (T - t_k)` bounded below by a positive constant, so `s = 0` is an interior point of the
time interval on which the `v^{(k)}` are defined. Section 8 F-8 of the note says Leray's bound "is not
used ... It is not needed: the construction needs only `M_k -> inf`". The first half is true of the
note's own route and the second half is false as an account of KNSS: `M_k -> inf` is what Leray's
bound is *not* used for, and the forward interval is what it *is* used for.

The note's route can close the endpoint without Leray, and its own S7 proof already does the work:
KNSS (3.10) bounds `||u||_{C^alpha_par(Q(z_0,R))}` on every backward parabolic cylinder
`Q(z_0,R) = B(x_0,R) x (t_0 - R^2, t_0)` contained in `R^n x (0,T)`, and such a cylinder may end at the
final time. So the bilinear term is uniformly parabolic-Holder up to `s = 0`, the heat term is smooth
there, and (4.2)-(4.3) handle the force term. Equicontinuity therefore holds on `R^3 x [-S+delta, 0]`.

**Fix.** S7 should conclude "converges locally uniformly on `R^3 x (-inf, 0]`", and F-8 should say
that KNSS use (6.1) for the forward interval and that the note replaces it by closed-endpoint
compactness (or, alternatively, by S4a: the forced local existence time from `t_k` is bounded below by
`c(gamma_1, ||g||_inf) M_k^{-2}`, which reproduces `B_k >= c > 0` with the force present).

---

### R3 (MINOR). S7's uniform-boundedness hypothesis is not satisfied in its second application

S7 assumes `|u^{(k)}| <= C` uniformly on `R^3 x (A_k, 0)`. In S13 it is applied to `w^{(k)}`, which
satisfies only KNSS (6.16), `|w^{(k)}| <= C/sqrt(-tau)` on `R^3 x (A_k,0)`, unbounded as `tau -> 0`,
together with (6.14), `|w^{(k)}| <= 2` outside the cylinder `C_k`. The cylinder escapes to infinity, so
the bound is uniform on compact sets and on `R^3 x (A_k, -delta)` for each `delta > 0`, but not on
`R^3 x (A_k, 0)`.

This is why KNSS need their separate `I(M) -> 0` estimate at the end of Theorem 6.2, which the note
does quote and does carry the force through correctly. But the invocation "So (4.4) holds and S7
applies" is not licensed by S7 as stated.

**Fix.** State S7 in two forms, or in one form with the hypothesis
"for every `delta > 0`, `sup_{R^3 x (A_k, -delta)} |u^{(k)}| <= C(delta)`", concluding convergence
locally uniformly on `R^3 x (-inf,0)`; and note that in the S8/S10 application the bound is uniform up
to `s = 0` (so R2's closed-endpoint version applies) whereas in the S13 application the endpoint is
handled separately by the `I(M)` estimate.

---

### R4 (MINOR). S13: `(w_1,w_3) = 0` does not follow from S12 as stated

S12 as quoted gives: a bounded weak solution on `R^2 x (-inf,0)` is `b(t)`, and a bounded ancient
*mild* solution of that form is *constant*. A constant is not zero, and `|w(0,0)| = 1` would be
perfectly consistent with a unit constant vector. S13 writes "applying S12 to the two-dimensional
field `(w_1,w_3)` gives `(w_1,w_3) = 0`", which skips the step that kills the constant.

The missing ingredient is in the note already, one line earlier: KNSS (6.16), `|w^{(k)}| <= C/sqrt(-tau)`,
passes to the limit, so `|w(x,tau)| <= C/sqrt(-tau) -> 0` as `tau -> -inf`. A constant obeying that is
zero. (This is also what makes `w = 0` follow from `(w_1,w_3) = 0`: with zero drift the remaining
component solves the heat equation and inherits the same decay.)

**Fix.** One sentence in S13.

---

### R5 (MINOR). S1's proof mishandles the singularity of the Riesz kernel; and `g` is not in `L^1`

(S1a) and (S1b) are true. The proof is not. It splits `int k(x-y) f_j(y) dy` at `|y| = |x|/2` and says
of the outer piece "`k` is integrable against the remaining decay". It is not: `k` is homogeneous of
degree `-3` on `R^3`, hence not locally integrable, and the outer region contains the singularity
`y = x`. The composition `R_i R_j` is a Calderon-Zygmund operator and needs its principal value or a
cutoff argument there.

**Fix.** Cut off instead at `|y - x| <= |x|/4`: on that piece `f 1_{|y-x|<=|x|/4}` is Schwartz with
every seminorm `O_N(|x|^{-N})` by (5), and `||R_iR_j h||_{L^inf} <= C ||h^||_{L^1} <= C ||h||_{W^{4,1}}`;
on the complement the kernel is non-singular and the two-region split as written is absolutely
convergent.

**Related.** (S1b) is correct to say `1 < q`: `(1+|x|)^{-3}` is not in `L^1(R^3)`, and generically
`R_iR_j f` decays exactly like `|x|^{-3}` (no faster, unless `int f = 0`), so `g = Pf` is genuinely not
in `L^1`. S16 then says "by (5) with `K = 4`, `|f| <= C(1+|x|)^{-4}`, so `f in L^q` for every
`q in [1,inf]` ... and by S0 we may take `div f = 0`". The `L^q` claim is for `f`; after the S0
replacement the correct claim is `g in L^q` for `1 < q <= inf`. Harmless, since Lemma B needs only
`q > 5/2`, but the sentence as written attributes an `L^1` bound to the object actually used.

Note also that `g = Pf` does **not** satisfy (5) (its decay stops at `|x|^{-3}`). Nothing downstream
needs it to: (A1) survives (Riesz transforms are bounded on `L^2` and on `H^s`), and S6 needs only the
`(1+|x|)^{-3}` bound. Worth one sentence in S0 so the reader does not assume the Clay class is
preserved for the force.

---

### R6 (MINOR). F-1 and section 11 misidentify which property of (5) is load-bearing

The note says the load-bearing point is that `P` maps `L^inf` to `BMO` only modulo constants, so a
bounded smooth force does not fit KNSS's framework, and that (5) is "exactly the repair".

Half of that is right and the operative half is missing. Take `f` smooth, bounded, with all
derivatives bounded, and **already divergence free**. Then S0 is vacuous, `g = f`, the Leray
projection is never applied, and KNSS's `L^inf`/`BMO` objection does not arise. The theorem still
fails for such an `f`, and it fails at S6: without decay one has no `|G| <= C/r`, so the `b(t)` mode in
S5 is not killed, so `u` is not shown to be a mild solution, so S8 does not run. That is exactly the
failure KNSS flag after their Theorem 6.2: "the statement fails, for trivial reasons, if we drop
assumption (6.6). (Consider `u(x,t) = b(t)`.)" A spatially constant divergence-free force generates
precisely that mode.

Relatedly, section 11's second "new" item, "the force enters the KNSS route only through terms that
vanish in the rescaling limit", is false as stated. The force also enters through `G` in the `b = 0`
step of S10 and S13, and that entry does not vanish under rescaling; it is a fixed-scale decay
statement. F-6 correctly notes the intermediate rescaling in S13 as one exception; the `b = 0` step is
a second one.

**Fix.** State both roles of (5): (a) `Pf` is an honest function, needed only when `div f != 0`;
(b) `|G| <= C(1+|x|)^{-3} <= C'/r`, needed always, and the one a bounded divergence-free force fails.

**Bookkeeping question answered.** The rescaling itself needs *only* `||grad^j g||_{L^inf} < inf` for
`j = 0,1,2`. Under `g_L(x,t) = L^3 g(Lx, L^2 t + T)` one has
`||grad^j g_L||_inf = L^{3+j} ||grad^j g||_inf`, so the gain is `L^{3+j}` with no reference to decay,
and the mild-solution norms used in (4.2), (4.3) and (4.4) are `L^inf` norms of `g` and its first two
derivatives. Decay does no work there. It does work only at S0 and S6.

---

### R7 (MINOR). S17's uniqueness justification is wrong, and S17 uses `(H)`

S17 concludes `Gamma = r u_theta = 0` "by uniqueness for this linear drift-diffusion equation with
bounded coefficients". The coefficients are not bounded: the equation is
`d_t Gamma + b.grad Gamma = nu(Laplace Gamma - (2/r) d_r Gamma) + r f_theta`, and `2 nu / r` blows up
at the axis. Nor is `Gamma` bounded: under `(H)`, `|Gamma| <= r sup|u|` grows linearly in `r`, so
uniqueness needs a growth class as well.

**Fix.** Run the maximum principle on `u_theta` instead. Its equation is
`d_t u_theta + b.grad u_theta + (u_r/r) u_theta = nu(Laplace - 1/r^2) u_theta + f_theta`;
the potential `-nu/r^2` has the favourable sign, `u_r/r` is bounded under `(H)` because `u_r` vanishes
on the axis and `u` is smooth with bounded second derivatives, and `u_theta` is bounded under `(H)`.
With `f_theta = 0` and `u_theta(.,0) = 0` the maximum principle gives
`||u_theta(t)||_inf <= e^{Ct} ||u_theta(0)||_inf = 0`.

**Consequence.** The closing sentence of section 7.2, "Note where each hypothesis was used:
`f_theta = 0` in S17, `(H)` in S19.2, and (5) in S19 and S20", is wrong: `(H)` is used in S17 too.

---

### R8 (MINOR). S10 and S13 state a conclusion their proofs do not establish

S10: "Then `sup_{R^3 x (0,T)} |u| < inf`, contradicting S3. Hence (6.4) cannot hold." The proof never
establishes a bound on `sup |u|`; it assumes (6.4), invokes S3 to get `M_k -> inf`, builds the ancient
limit, and contradicts `|v(0,0)| = 1`. KNSS's Theorem 6.1 does conclude `|u| <= M(C)` quantitatively;
the note's argument does not, and does not need to. S13 has the same mismatch.

**Fix.** State both as "then (6.4) [resp. (6.5)] is impossible", which is what is proved and what is
used.

---

### R9 (MINOR). Section 7.3's stated reason why S5 has no periodic analogue is wrong

The note says the decomposition S5 "relies on the Liouville theorem for bounded curl-free
divergence-free fields on `R^n`" and "has no periodic analogue in the form used". The Liouville input
does have a periodic analogue, and a trivial one: a curl-free divergence-free field on `T^n` is
harmonic, hence constant by compactness. KNSS's Lemma 3.1 uses exactly that statement, and it is
easier on the torus than on `R^n`.

The genuine obstruction is the one the note names first, the explicit `R^3` Stokes kernels
`K_{ij}, K_{ijk}` and the Duhamel formula (4.1). That part stands.

**Fix.** Delete the Liouville clause.

---

### R10 (MINOR). S21 is correct, but its scope is not stated, and "empty" should be "trivial"

The proof of S21 is correct. Every step was rechecked:
`u(Q_a x) = u(Q(x-a)) = Q u(x-a) = Q u(x)`; `Q_0^{-1} Q_a(x) = x + (Q^{-1} - I)a`;
`(Q^{-1} - I)e_1 = (cos alpha - 1, -sin alpha, 0)`, `(Q^{-1} - I)e_2 = (sin alpha, cos alpha - 1, 0)`;
these are arbitrarily short and span the horizontal plane as `alpha -> 0`; continuity then gives
independence of `x_1, x_2`; axisymmetry of a field of `x_3` alone forces `u = (0,0,u_3(x_3))`;
`div u = 0` gives `u_3` constant. Script checks C4 and C4' verify the translation identity and the
spanning. The "Consequence for (D)" is also correct: `f` continuous, periodic and axisymmetric is
`(0,0,f_3(x_3,t))` by the same argument, and periodicity of `p` (Fefferman's errata to (8)) forces
`c' = mean(f_3)`.

Two things are missing.

* The argument uses the **full continuous** rotation group about a fixed vertical line, through the
  limit `alpha -> 0`. That is essential and should be said. The only lattice-compatible rotations are
  the multiples of `pi/2`, and invariance under that finite subgroup does not force triviality at all;
  the four-fold-symmetric periodic fields are a large class. So the sharp statement is: *the notion of
  axisymmetry that makes KNSS-type theorems apply is trivial on `T^3`; the notion that is compatible
  with the lattice is not, and is not what any of the cited theorems are about.*
* "the class is empty" (end of the Consequence paragraph) contradicts the same paragraph's
  "consists of the fields `u(x,t) = (0,0,c(t))`". The class is trivial, not empty.

---

### R11 (MINOR). Section 9's Alpoge-Buckmaster sentence is flatter than the evidence

"It fails item 1 (`[FENCES]` Corollary A: velocity bounded)." The reported Alpoge-Buckmaster result is
for **Euler**, `nu = 0`, and `[FENCES]` Corollary A covers `1 < alpha <= 2` only, so it says nothing
about the Euler result itself; what it says is that the *mechanism* cannot be pushed to any
`alpha > 1`. And the preprint was not found (`[FENCES]` section 0 and the note's own section 10 row:
"not found ... reported, never cited"). The hedge "as reported, and unverified" appears earlier in the
sentence and then the conclusion is stated flatly.

**Fix.** "The mechanism, as reported, keeps the velocity bounded, so by `[FENCES]` Corollary A it
cannot be transplanted to any `alpha > 1`, in particular not to `alpha = 2`."

---

### N-0 (NONE). Machine checks, and the two decay hypotheses

**Script.** Re-run from a copy: seven PASS, exit 0, 7.3 s. Reading the source:

* C1 composes the scaling map explicitly and uses `sp.Subs(...).doit()`, which substitutes correctly
  inside `Derivative` objects; it verifies `F_L = L^3 (transported F)` for all three momentum
  components and `div u_L = L^2 (transported div u)`. Correct, and `nu` is carried symbolically, so
  the covariance is `nu`-independent as the note needs. C1' is a weak control (it only checks
  `L^3 != L^2`), but it is honestly labelled a control.
* C2a, C2 are correct: `(Delta u)_r = Delta u_r - u_r/r^2` in the swirl-free case is used in `Rr`,
  `(curl A)_theta = d_z A_r - d_r A_z` in `curlR`, and `omega_theta = d_z u_r - d_r u_z`. The target
  matches (7.1) sign for sign, including the stretching term `+ (u_r/r) omega_theta` on the right.
* C3 verifies the operator identity `L_omega[r Om] = r L_Om[Om]` and does not carry the source term
  explicitly. That is sufficient: the `omega` equation is `L_omega[omega] = F_theta`, so
  `r L_Om[Om] = F_theta`, that is `L_Om[Om] = F_theta/r`. The docstring's phrasing suggests the source
  is checked; it is implied, not checked. Not an error.
* C4, C4' as described in R10.

**(Dec).** F-3 is right to carry it. (4) and (5) do not give it. KNSS's own remark attributes (6.6) to
Brandolese for mild solutions with fast-decaying data in the **unforced** problem, and even there it is
a theorem about a decay class, not a consequence of the Clay conditions; the statement needed is
uniform in `t` up to the singular time, which is the hard direction. One sharpening the note could
take for free: S13 uses `(Dec)` only to get `lambda_k = |x'_k|` bounded, so the weaker hypothesis
`limsup_{r -> inf} sup_{t<T} r |u| < inf` suffices.

**(7.3).** Checked and correct, including the vanishing of `F_theta` on the axis (for smooth
axisymmetric swirl-free `f`, `f_r` is odd in `r` and `f_z` even, so `d_z f_r` and `d_r f_z` are both
`O(r)`) and the two-region bound.

**S6, S4, S5, S20.** Checked in detail and correct. In particular the Kato-Ponce step of S20 is a
legitimate instance of Grafakos-Oh: with `r = 2`, `n = 3`, `sigma = s = 3 > max(0, n/r - n) = 0`, and
the pairs `(p_1,q_1) = (3,6)`, `(p_2,q_2) = (6,3)`, both terms of the fractional Leibniz rule collapse
to `||u||_{L^6} ||Lambda^s u||_{L^3}` because `u (x) u` is symmetric. The note prints the pairs in the
opposite order from its own conclusion, which reads as an error but is only a labelling slip; the
estimate is right. The subsequent chain `||Lambda^s u||_{L^3} <= C||Lambda^{s+1/2}u||_{L^2} <= Y^{1/2}X^{1/2}`
and the Young step with exponents `4/3` and `4` are correct.

---

## Citations that check out

Read at source during this review, not taken from the note.

| Claim in the note | Source, read | Status |
|---|---|---|
| S9 = KNSS Theorem 5.3: bounded weak solution on `R^3 x (-inf,0)`, axisymmetric, `\|u\| <= C/sqrt(x_1^2+x_2^2)`, then `u = 0` | arXiv:0709.3599 p. 15, (5.17) | Verbatim, accurate |
| S12 = KNSS Theorem 5.1 (`R^2`, bounded weak, `u = b(t)`) and Remark 6.1 (bounded ancient mild of that form is constant) | arXiv:0709.3599 pp. 12, 18 | Verbatim, accurate |
| KNSS weak-solution test identity and the sign convention used in S5 | arXiv:0709.3599 section 3 | Signs recomputed and match |
| KNSS (3.10) parabolic Holder bound on any `Q(z_0,R) ⊂ R^n x (0,T)` | arXiv:0709.3599 (3.10) | Verbatim; this is what repairs R2 |
| KNSS Lemma 3.1 decomposition `u = v + w + b`, Remark 3.1 | arXiv:0709.3599 section 3 end | Verbatim, accurate |
| KNSS Theorem 6.1 hypotheses (`u in L^inf(R^3 x (0,T'))` for each `T' < T`, weak, (6.4)) and the `\|x'_k\| <= C/M_k` step | arXiv:0709.3599 p. 20 | Verbatim; note that KNSS's own `L^inf` hypothesis is `(H)` at `k = 0` |
| KNSS Theorem 6.2 hypotheses (6.5), (6.6), the two-stage rescaling (6.8), (6.12), bounds (6.9)-(6.11), (6.14)-(6.16), `I(M) -> 0` | arXiv:0709.3599 pp. 21-24 | Verbatim; the note's transcription is faithful except R3 and R4 |
| KNSS remark that Theorem 6.2 "fails, for trivial reasons, if we drop (6.6). (Consider `u(x,t) = b(t)`)" | arXiv:0709.3599 p. 22 | Verbatim; this is the evidence for R6 |
| KNSS use of Leray (6.1) to give `B_k = M_k^2(T-t_k)` bounded below | arXiv:0709.3599 pp. 19-20 | Verbatim; this is the evidence for R2 |
| CSTY I Theorem 1.1: axisymmetric `(v,p)` on `D = R^3 x (-T_0,0)`, `v` smooth in `x` and Holder in `t`, `p in L^{5/3}(D)`, `\|v\| <= C_*(r^2-t)^{-1/2}`, conclusion `v in L^inf(B_R x [-T_0,0])` | arXiv:math/0701796 p. 2 | Verbatim; the note's section 5 Remark is accurate |
| CSTY II Theorem 1.1: `v^0 in H^{1/2}`, `r v_theta^0 in L^inf`, `p in L^{5/3}(D)`, bounds (1.2) and (1.3) with `eps in [0,1]`, `eps = 1` is (1.2), `eps = 0` in the appendix | arXiv:0709.4230 p. 2 | Verbatim; the note's F-5 is accurate |
| Seregin-Sverak Theorems 1.1-1.3, conditions (1.9), (1.10) on the poloidal part `v-bar`, `v in L^3(Q)`, `q in L^{3/2}(Q)`, Theorem 1.3 the local no-swirl statement | arXiv:0804.1803 pp. 4-5 | Verbatim. One omission: Theorems 1.2 and 1.3 also assume `v` essentially bounded on `B x (t_0-R^2,t')` for each `t' < t_0`, a local `(H)`. Theorem 1.1, the one the note leans on to remove `(Dec)`, does not |

The note's own section 10 table is accurate on every row I could re-check.

---

## Corrected theorem statement

> **Theorem (corrected).** Let `nu > 0`. Let `u_0` be smooth, divergence free and axisymmetric on
> `R^3` and satisfy (4). Let `f` be smooth and axisymmetric on `R^3 x [0,inf)` and satisfy (5). Let
> `(u,p)` solve (NS), be `C^inf` on `R^3 x [0,T)` with bounded energy (7) there, let `T < inf` be its
> first singular time, and assume
>
>     (H*)  u in C([0,T'); H^s(R^3))  for every s >= 0 and every T' <= T.
>
> `(H*)` implies `(H)` by Sobolev embedding and implies `int_0^T ||grad u||_{L^2}^2 dt < inf` through
> the energy identity. It is the solution class of `[FENCES]` Lemmas A and B. It is not implied by
> (4)-(7); see R1.
>
> Define `Sigma := { x_0 : sup_{B_rho(x_0) x (T-rho,T)} |u| = inf for every rho > 0 }`. Then:
>
> **(i)** `Sigma ⊂ {r = 0}`. This uses `(H*)` only through the energy consequence
> `u in L^inf(0,T;L^2) ∩ L^2(0,T;H^1)`, and does not use `(Dec)`. If in addition
> `sup_{t<T} |u(x,t)| -> 0` as `|x| -> inf`, then `Sigma` is nonempty.
>
> **(ii-a)** There is no `C_* < inf` with `|u(x,t)| <= C_*(r^2 + (T-t))^{-1/2}` on `R^3 x (0,T)`.
>
> **(ii-b)** If in addition `(Dec)` holds (KNSS (6.6); implied by (ii-a)'s bound, and implied by
> `limsup_{r -> inf} sup_{t<T} r|u| < inf`), there is no `C_* < inf` with
> `|u(x,t)| <= C_*(T-t)^{-1/2}` on `R^3 x (0,T)`; equivalently
> `limsup_{t -> T-} sqrt(T-t) sup_x |u(x,t)| = inf`, the singularity is type II in KNSS's sense.
>
> **(iii)** If `u_{0,theta} = 0` and `f_theta = 0` identically, then no such `T` exists: the solution
> is global and smooth.

Everything else in the note's section 9 reading stands with "unconditional given (5)" struck from
items 1 and 2 and replaced by "given (5) and `(H*)`".

---

## What remains

1. `(H*)` is not discharged, and this review does not discharge it. Discharging it means a uniqueness
   theorem in the Clay class `C^inf ∩ L^inf_t L^2_x` together with a no-escape-to-spatial-infinity
   statement for the `H^s` norm. `[FENCES]` section 3 already flags the first half. This is the single
   largest hole in the whole two-note stack, and it is the same hole in both notes.
2. `(Dec)` remains a hypothesis of (ii-b). F-3 is honest about it. The route to removing it is
   Seregin-Sverak Theorem 1.1, force-extended; F-4 is honest that this was not done.
3. The intermediate CSTY II range `0 < eps < 1` remains uncovered (F-5). Confirmed against
   arXiv:0709.4230.
4. On `T^3` there is still no statement, and after R10 the reason is sharper: the axisymmetric class in
   the sense the theorems need is trivial, and the lattice-compatible symmetry class is untouched by
   any of this.
5. Nothing here narrows the corridor further, and nothing here opens it. FL-000 stands.

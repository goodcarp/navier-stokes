# Actual NS time polynomials, exact residuals, and a first-order validation obstruction

The first-order time polynomial is an exactly divergence-free whole-space
approximation, but it is a poor candidate for the current H4 core-cone
validator. A rigorous residual lower bound prevents that validator from
certifying even `b(T)>=Omega(T)` with this approximation at `T=0.001`.
This is a limitation of the approximation **and the specified norm-radius
certificate**, not a lower bound on the actual solution error or a failure
theorem for the NS evolution. The quadratic polynomial removes the leading
residual; its usable interval still requires actual coefficient norms.

## 1. One fixed datum and full NS coefficients

Fix the explicit pass9 datum

\[
 a_0=u_0=M+1020v+\tfrac14w_L,\qquad \nu=1/1000.
\]

Every profile is the unchanged smooth compact profile; `u0` is supported
in the ball of radius six. No pressure retuning is performed at a later
time. Use the full whole-space Leray bilinear operator

\[
 B(a,b)=\mathbb P(a\cdot\nabla b),\quad
 N(u)=\nu\Delta u-B(u,u),
\]

and define the actual initial time derivatives

\[
 a_1=N(u_0)=N_0,
\]
\[
 a_2=DN(u_0)a_1
 =\nu\Delta a_1-B(a_0,a_1)-B(a_1,a_0)=N_1,
\]
\[
 a_3=\nu\Delta a_2-B(a_0,a_2)-B(a_2,a_0)-2B(a_1,a_1).
 \tag{1}
\]

These are obtained from the actual full field, including the core, axial
endcaps, seed envelope, and all pressure interactions. They are not the
coefficients of the planar or affine model. Each is smooth, solenoidal,
odd, and C4-equivariant. Thus both

\[
 v_1=a_0+ta_1,\qquad v_2=a_0+ta_1+\tfrac12t^2a_2
 \tag{2}
\]

are smooth, exactly divergence-free full-space fields with exactly the
chosen initial datum. They are not asserted to solve unforced NS.

## 2. Exact vector residuals, without an omitted pressure term

Define `R(v)=v_t−nu Delta v+B(v,v)`. Bilinearity gives identically

\[
 \boxed{\mathcal R(v_1)=-ta_2+t^2B(a_1,a_1),}
 \tag{3}
\]

\[
 \boxed{\mathcal R(v_2)=
 -\tfrac12t^2a_3
 +\tfrac12t^3[B(a_1,a_2)+B(a_2,a_1)]
 +\tfrac14t^4B(a_2,a_2).}
 \tag{4}
\]

No Taylor remainder is hidden in these residual formulas: they are exact
finite polynomials for the chosen approximate fields. In particular the
`−2B(a1,a1)` in `a3` is necessary. A numerical pressure construction must
represent the Leray operator in (1), (3), and (4) on the whole space, or
charge its error separately.

## 3. A quantitative obstruction for the linear candidate

Let `ell(a)=r_* overline(a_theta)(r_*,4)` be the fixed-circle mean
angular-momentum functional at the previously isolated exact initial
maximum. The audited full-NS initial jet includes the explicit amplitude
`A=1020` and gives

\[
 \ell(a_2)=G_{tt}(0,r_*,4)<-13000,\qquad 0<r_*<1/4.
 \tag{5}
\]

The archive's scalar embedding implies
`|ell(f)|<=r_* ||f||H4`. Write `a=−ell(a2)>13000` and
`b=ell(B(a1,a1))`, without making any assumption on the size or sign
of `b`. Equation (3) gives `ell(R(v1))=at+bt²`.

There is a cancellation-resistant lower bound valid for **every** `b`:

\[
 \int_0^T|at+bt^2|\,dt\ge aT^2/8.
 \tag{6}
\]

For an exact rational proof, take the test function `phi(s)=61/64` on
`[0,4/5]` and `phi(s)=−1` on `(4/5,1]`. It has `|phi|<=1`,
`integral phi(s)s² ds=0`, and `integral phi(s)s ds=1/8`.
Testing `at+bt²` against `phi(t/T)` proves (6). The optimal constant
is `2^(−2/3)−1/2>1/8`, but the rational bound is enough.

Consequently any valid H4 residual majorant for this exact polynomial
must satisfy

\[
 \boxed{\int_0^T\delta_4(t)\,dt
 \ge\int_0^T\|\mathcal R(v_1)(t)\|_4\,dt>6500T^2.}
 \tag{7}
\]

This cannot be evaded by hoping the unknown quadratic residual coefficient
cancels the linear coefficient over the whole interval.

The H4 validator from pass9 uses the nonnegative-growth majorant
`rho′>=560||v1||5 rho+270rho²+delta4`, or its displayed explicit
supersolution. Hence `rho(T)>=integral delta4`. On the other hand,
write the exact initial off-neutral value

\[
 d=\tfrac12(1020^2C_v+C_L/16-204/35).
\]

The existing bounds `C_v<279/40000000` and `C_L<3/2` give
`0<d<4/5`. The linear polynomial's central gradient is exactly

\[
 b_{v_1}(t)=1+(2+d)t,\qquad
 \Omega_{v_1}(t)=1+2t,
\]

so its apparent core-cone margin is only `dT`. The pass9 H4 endpoint
test for `b>=Omega` requires

\[
 b_{v_1}(T)-\Omega_{v_1}(T)-2\rho(T)\ge0,
 \quad\hbox{hence}\quad \rho(T)<(2/5)T.
 \tag{8}
\]

Equations (7) and (8) are incompatible whenever

\[
 \boxed{T\ge1/16250\simeq6.154\times10^{-5}.}
 \tag{9}
\]

At the proposed diagnostic time `T=1/1000`, the residual forces the
validator's radius above `13/2000=0.0065`, while (8) permits less than
`1/2500=0.0004`. A strictly growing ratio cone is more demanding still.

This rules out this **linear approximation + nonnegative H4 residual
majorant + isotropic H4 cone-error test** for that endpoint. It does not
rule out a sharper signed/linearized propagator estimate, an observable-
specific error certificate, a different successor class, the quadratic
approximation, or the actual NS stage. In particular (7) is not a lower
bound on `||u(T)−v1(T)||4`; parabolic evolution and cancellation can make
the actual error smaller than an absolute residual accumulation.

## 4. Practical recursive upper bounds for the quadratic candidate

For the derivative-sum Sobolev norms, put `M_j,s>=||a_j||Hs`. For
`s>=3`, elementary product counting gives

\[
 \|B(a,b)\|_s\le K_s\|a\|_{s+1}\|b\|_{s+1},\qquad
 K_s=3\,2^s\sqrt{\binom{s+3}{3}},
 \quad \|\Delta a\|_s\le3\|a\|_{s+2}.
 \tag{10}
\]

In particular one may take `K4=288` and `K5=720`. The lower derivative
factor in each differentiated product fits the scalar H2 embedding;
the other has order at most `s+1`.

More generally, the exact NS time-jet recurrence is

\[
 a_{n+1}=\nu\Delta a_n-
       \sum_{j=0}^n\binom nj B(a_j,a_{n-j}),
\]

and a directly computable choice of successive upper bounds is

\[
 M_{n+1,s}:=3\nu M_{n,s+2}
 +K_s\sum_{j=0}^n\binom nj M_{j,s+1}M_{n-j,s+1}.
 \tag{11}
\]

Starting with valid initial bounds, this recurrence inductively gives
valid bounds for each actual coefficient. Direct certified norms of the actual coefficients
can substantially improve this elementary recursion.

Equations (3)–(4) yield the usable residual majorants

\[
 \delta_{4,1}(t)=tM_{2,4}+288t^2M_{1,5}^2,
\]
\[
 \boxed{\delta_{4,2}(t)=
 \tfrac12t^2M_{3,4}
 +288t^3M_{1,5}M_{2,5}
 +72t^4M_{2,5}^2.}
 \tag{12}
\]

The coefficient required by the error-growth factor is bounded by
`||v2(t)||5<=M0,5+tM1,5+(t²/2)M2,5`. Initial derivative bounds through
H8 suffice for the first polynomial's crude residual bound; H10 suffices
for the quadratic bound through (11). The latter includes the worst
`nu³ Delta³u0` term in `a3`. These norms are finite for the actual cutoff
but may be very large. Finiteness alone does not make (12) useful at
`T=0.001`.

A concrete next computation is to enclose `M1,5`, `M2,5`, and `M3,4`
from the actual fields, or to bound the complete residual polynomial
directly with interval correlations retained. The positive endpoint
margins must then survive the resulting radius. No current file certifies
these three norms or that outcome.

## 5. Harmonic pressure tails are part of the polynomial

Let `N=(-Delta)^(-1)` denote the Newtonian potential. The initial pressure
and its actual first time derivative are

\[
 p_0=\partial_{ij}N*(a_{0i}a_{0j}),\qquad
 p_1=\partial_{ij}N*(a_{0i}a_{1j}+a_{1i}a_{0j}).
 \tag{13}
\]

Both tensors in (13) are compactly supported in `|x|<=6`, even though
`a1` itself is not compact. Outside this ball, exactly

\[
 a_1=-\nabla p_0,\qquad a_2=-\nabla p_1.
 \tag{14}
\]

The Laplacian of `a1` vanishes there and every local mixed term with
`a0` vanishes. Thus (14) gives genuine harmonic velocity tails, generally
of order `|x|^(−4)`, not zero boundary data. Oddness and C4 symmetry do
not imply that their anisotropic leading stress moment vanishes.

For explicit tail estimates let `T_l` be the two compact tensors in
(13) and `Q_l=sum_ij ||T_l,ij||L1`. One has
`Q0<=3||a0||L2²` and `Q1<=6||a0||L2||a1||L2`. Define

\[
 S_n=\sum_{p=0}^{\lfloor n/2\rfloor}
 \frac{n!}{2^pp!(n-2p)!}(2n-2p-1)!!,
 \qquad(-1)!!=1.
\]

The Cartesian Newtonian derivative formula gives
`|D^alpha (4pi|x|)^(-1)|<=S_n |x|^(−n−1)/(4pi)` for `|alpha|=n`.
For `l=0,1`, `a_(l+1)` in (14), and `R>=12`, this implies the
explicit whole exterior bound

\[
 \sum_{|\alpha|\le4}\int_{|x|\ge R}|D^\alpha a_{l+1}|^2dx
 \le\frac{3Q_l^2}{4\pi}
 \sum_{j=0}^4\binom{j+2}{2}
 \frac{S_{j+3}^2\,2^{2j+8}}{2j+5}R^{-(2j+5)}.
 \tag{15}
\]

The estimate uses `|x−y|>=|x|/2`; its constants are conservative.
Tail bounds for `v1,v2` follow by the triangle inequality with their
time coefficients. These estimates supplement the interior calculation;
they do not remove the pressure contact term inside the support.

The higher nonlinear coefficients also admit compact-source evaluations.
For divergence-free fields with compactly supported vorticity,

\[
 B(a,a)=-\mathbb P(a\times\operatorname{curl}a),
\]
\[
 B(a,b)+B(b,a)=
 -\mathbb P[a\times\operatorname{curl}b+
             b\times\operatorname{curl}a].
 \tag{16}
\]

The gradients of kinetic-energy products have been projected away.
All finite NS time coefficients here have vorticity supported inside
the original ball, as follows inductively from the differentiated
vorticity equation. Consequently the Lamb-vector sources needed in
(4) and (16) are compact. Their whole-space Leray tails must still be
retained. This provides a practical alternative to multiplying truncated
nonlocal velocity fields on a periodic box.

## 6. A heat-resummed alternative and its exact remaining residual

If the high cutoff derivatives make the polynomial diffusion terms
too costly, retain the exact heat semigroup `H(t)=exp(nu t Delta)`.
With `S01=B(a0,a1)+B(a1,a0)`, define

\[
 v_{H,1}(t)=H(t)a_0-\int_0^tH(t-s)B(a_0,a_0)\,ds,
\]
\[
 v_{H,2}(t)=H(t)a_0-\int_0^tH(t-s)
                   [B(a_0,a_0)+sS_{01}]\,ds.
 \tag{17}
\]

They are exactly divergence free, have the same initial datum, and
match `a1` and, for the second field, `a2` at zero. Their exact residuals
are respectively

\[
 B(v_{H,1},v_{H,1})-B(a_0,a_0),
\]
\[
 B(v_{H,2},v_{H,2})-B(a_0,a_0)-tS_{01}.
 \tag{18}
\]

These fields pay viscosity from the outset. They are still approximations
until (18), reconstruction error, and all tails are bounded. In particular
the obstruction (7) is for the literal polynomial `v1`, not automatically
for (17).

No heat Taylor remainder may be discarded when replacing these fields
by polynomials. The exact identity is

\[
 H(t)f-\sum_{j=0}^p\frac{(\nu t)^j}{j!}\Delta^jf
 =\frac{\nu^{p+1}}{p!}\int_0^t(t-s)^pH(s)\Delta^{p+1}f\,ds,
\]

with Hs norm at most `(nu t)^(p+1)||Delta^(p+1)f||Hs/(p+1)!`.
For a compact `f` supported in radius `R0`, a weighted heat-energy
estimate also gives

\[
 \|1_{|x|\ge R}D^\alpha H(t)f\|_2
 \le e^{-(R-R_0)^2/(4\nu t)}\|D^\alpha f\|_2.
 \tag{19}
\]

For example, apply the heat energy estimate with a Lipschitz distance
weight and optimize `exp(−a(R−R0)+a²nu t)` in `a`. Formula (19) applies
to compact sources. It must not be applied as though the pressure tails
in (14) were compact. A smooth local/tail split must charge its cutoff
derivatives; heat contraction then controls the already bounded tail.

## Scope

The practical next ansatz should include at least the full quadratic
coefficient or a justified heat-resummed evolution. The first-order
polynomial is explicitly excluded only for the stated H4 core-cone
certificate at the proposed time; the actual full NS stage remains open.
No model has been substituted for the full initial coefficients, no
endpoint has been reset, and no useful actual duration is claimed.
Using a residual as a chosen physical force would be a separate permitted
forced construction requiring its own pressure choice, spatial behavior,
and assembly estimates; it is not an unforced error validation.

# A center-independent gradient-decay repair

This is a finite conditional PDE lemma for the scalar drift–diffusion
equation. It repairs the location-dependent constant in the live Claude
`pmax-h3v/PROOF.md`, Corollary F′. It does not establish the log-clock
assembly, a Navier–Stokes lifespan, or the drift bootstrap assumed there.
The source review and exact live-file pin are recorded in the companion
checker output. No Claude source or published archive was changed.

## 1. Precise statement

Let \(\nu>0\), \(n\ge1\), and let

\[
 \partial_t u+b(x,t)\cdot\nabla u=\nu\Delta u
 \quad\hbox{on }\mathbb R^n\times(0,T].                 \tag{1}
\]

Assume that \(b\) is continuous in time, continuously differentiable in
space, and globally Lipschitz in space, uniformly on the time interval:
\(\|D_xb\|_{\mathrm{op}}\le L<\infty\). It suffices that these properties
hold locally with the displayed global bound and that \(b(0,t)\) is
bounded on finite intervals. This gives a complete forward/backward ODE
flow. For the direct proof, assume that \(u\), \(\nabla u\), their needed
time derivatives, and spatial derivatives through order three are
continuous on positive-time compact cylinders, and that (1) holds there.
No a priori *global* bound on \(\nabla u\) is assumed.
Here \(L\) denotes the drift Lipschitz bound, not the log-shell parameter
called \(L\) in the Claude budget.

Fix \(x\), \(t>0\), \(0<h<t\), and \(r>0\). Let the backward-centered
trajectory solve

\[
 X'(s)=b(X(s),s),\qquad X(t)=x,\qquad t-h\le s\le t,
\]

and let
\(M_Q=\sup\{|u(X(s)+y,s)|:|y|\le r,\ t-h\le s\le t\}\).
Set \(C_n=\max(4n,24)\). Then

\[
 \boxed{\quad
 |\nabla u(x,t)|\le M_Q
 \sqrt{\frac{1+2Lh+C_n\nu h/r^2}{2\nu h}}.
 \quad}                                                     \tag{2}
\]

There is no dependence on \(|b(x,t)|\), the trajectory's position,
\(D^2b\), a time-Hölder modulus, or the length \(T\) beyond the chosen
interval's displayed bounds. The norm in \(L\) is the operator norm;
\(|D^2u|\) below is the Frobenius norm. This is an a priori estimate, not
an existence assertion under the listed drift assumptions.

## 2. Proof with a translating ball and a Bernstein cutoff

Write \(\sigma=s-(t-h)\),
\(w(y,\sigma)=u(X(s)+y,s)\), and
\(a(y,\sigma)=b(X(s)+y,s)-b(X(s),s)\). Then

\[
 \mathcal D w=0,\qquad
 \mathcal D=\partial_\sigma+a\cdot\nabla-\nu\Delta,
 \qquad |a(y,\sigma)|\le L|y|,
 \qquad \|D_ya\|_{\mathrm{op}}\le L.                 \tag{3}
\]

This coordinate change is a translation. It introduces neither a
deformed diffusion tensor nor a Jacobian factor. In particular it does
not assume the full drift is affine.

Let \(q=|\nabla w|^2\), \(H=D^2w\),
\(d=|y|^2/r^2\), and \(\zeta=1-d\) on the closed ball. Differentiating
(3) gives

\[
 \mathcal Dq=-2\nu|H|^2
   -2\nabla w\cdot(Da)^T\nabla w
 \le-2\nu|H|^2+2Lq,\qquad
 \mathcal D(w^2)=-2\nu q.                             \tag{4}
\]

The cutoff has
\(|\nabla\zeta|^2=4d/r^2\),
\(\Delta\zeta=-2n/r^2\). Its drift term satisfies
\(a\cdot\nabla(\zeta^2)\le4L\zeta d\). The potentially dangerous
product term is bounded by

\[
 \begin{split}
 -2\nu\sigma\nabla(\zeta^2)\cdot\nabla q
 &\le8\nu\sigma\zeta|\nabla\zeta|\,|H|\sqrt q\\
 &\le2\nu\sigma\zeta^2|H|^2
       +8\nu\sigma|\nabla\zeta|^2q.                  \tag{5}
 \end{split}
\]

The first term cancels the negative Hessian term in (4). Therefore for
\(F=\sigma\zeta^2q+A w^2\),

\[
 \begin{split}
 \mathcal DF\le\bigg[\zeta^2
  +2L\sigma(1-d^2)
  +\frac{\nu\sigma}{r^2}\{4n(1-d)+24d\}
  -2\nu A\bigg]q.                                   \tag{6}
 \end{split}
\]

Indeed the diffusion cutoff terms are
\(-\nu\Delta\zeta^2+8\nu|\nabla\zeta|^2
=\nu\{4n(1-d)+24d\}/r^2\), and the drift terms combine as
\(2L\zeta^2+4L\zeta d=2L(1-d^2)\).
For \(0\le d\le1\), \(0\le\sigma\le h\), choose

\[
 A=\frac{1+2Lh+C_n\nu h/r^2}{2\nu}.
\]

Then \(\mathcal DF\le0\). On the initial face \(\sigma=0\) and lateral
boundary \(|y|=r\), \(F=A w^2\le A M_Q^2\). The elementary parabolic
maximum principle on this bounded cylinder gives \(F\le A M_Q^2\).
At \((y,\sigma)=(0,h)\), \(h q\le F\), proving (2). Each cylinder is
compact and stays at positive times, so this argument assumes no
unproved uniform bound on gradients at infinity.

As a useful independent consequence, if \(|u|\le B\) throughout the
time slab, let \(r\to\infty\) in (2). This yields

\[
 \|\nabla u(\cdot,t)\|_\infty
 \le B\sqrt{\frac{1+2Lh}{2\nu h}}.                    \tag{7}
\]

Thus even the *global* gradient bound required before the weighted
maximum-principle argument follows without presupposing gradient decay.

## 3. Polynomial tails remain polynomial gradient tails

Assume now \(b(0,s)=0\) and, throughout the slab,

\[
 |u(z,s)|\le D\,\Xi_p(z),\qquad
 \Xi_p(z)=(1+|z|^2/R^2)^{-p/2},\quad p>0.             \tag{8}
\]

For \(z=X(s)+y\), \(|y|\le r\), the flow growth bound implies
\(|x|\le e^{Lh}(|z|+r)\). For nonnegative \(a,c\),

\[
 \sqrt{1+(a+c)^2}\le\sqrt{1+a^2}+c
                  \le(1+c)\sqrt{1+a^2}.
\]

Consequently, **at every spatial point**, not just in an exterior region,

\[
 \Xi_p(z)\le e^{pLh}(1+r/R)^p\Xi_p(x).
\]

Inserting this into (2) proves

\[
 \boxed{\quad
 |\nabla u(x,t)|\le
 D e^{pLh}(1+r/R)^p
 \sqrt{\frac{1+2Lh+C_n\nu h/r^2}{2\nu h}}\,\Xi_p(x).
 \quad}                                                     \tag{9}
\]

If \(b(0,s)\ne0\) but \(|b(0,s)|\le B_0\), the same proof replaces
\(1+r/R\) by \(1+(r+B_0h)/R\). A spatially constant drift can translate a
tail, so its magnitude must enter a tail bound about a *fixed* origin;
it still does not enter the local gradient constant (2).

For the Claude application take \(n=5\), \(p=9\),
\(D=\mathfrak A(M/\rho_0)e^{\Lambda_9T}\) from its scalar-tail lemma,
and \(R\ge\rho_0>0\). For any \(t\in[t_0,T]\), choose

\[
 h=\min\left(\frac{t_0}{2},\frac1L,\frac{\rho_0^2}{\nu}\right),
 \qquad r=\sqrt{\nu h};                               \tag{10}
\]

omit \(1/L\) when \(L=0\). These choices give \(C_5=24\),
\(Lh\le1\), \(r\le\rho_0\le R\), and the explicit coarse version

\[
 |\nabla u(x,t)|\le
 D(2e)^9\sqrt{\frac{27}{2\nu h}}\,\Xi_9(x).            \tag{11}
\]

The sharper constant is the one in (9), not the universal coarse factor
in (11). In the notation of the source's Corollary F′, (11) has the form
\(2^9 C_*(t_0)\mathfrak A(M/\rho_0^2)e^{\Lambda_9T}\Xi_9\) with
\(C_*(t_0)=e^9\rho_0\sqrt{27/(2\nu h)}\).
It is independent of the center and of \(K_2\); its growth as
\(t_0\downarrow0\) or \(\nu\downarrow0\) is explicit.
In particular, this bound alone does not justify a weighted-gradient
initial trace or passing a weighted maximum estimate from \(t_0>0\) to
\(t_0=0\); that passage needs its own datum/solution argument.

In particular \(|x|^2|\nabla u|\to0\) uniformly on \([t_0,T]\).
For \(p>n\), integrating (9) gives

\[
 \|\nabla u(\cdot,t)\|_1\le
 D e^{pLh}(1+r/R)^p
 \sqrt{\frac{1+2Lh+C_n\nu h/r^2}{2\nu h}}
 \frac{\pi^{n/2}R^n\Gamma((p-n)/2)}{\Gamma(p/2)}.       \tag{12}
\]

For \((n,p)=(5,9)\), the last factor equals
\(16\pi^2R^5/105\). This supplies finiteness without an integration by
parts that already assumes the desired gradient tail. It does not claim
the favorable initial-data scaling \(MR^3\).

The source's displayed inner-region term in F′(c) reuses
\(\mathfrak G_0\), an **initial** weighted-gradient bound, at positive
times without a propagation estimate. That particular constant is not
justified by the stated proof. Replace the whole display by (12), or use
the interior volume times the valid global bound (7). Only finiteness was
said to be used downstream, so this repair introduces no certified
quantitative budget improvement.

## 4. What this does and does not repair about regularity

The proof above applies directly to actual classical smooth solutions on
compact positive-time intervals. Its constants use only the displayed
Lipschitz bound. It does **not** infer smoothness of such a solution from
the draft's undefined phrase “bounded solution.”

Spatial Lipschitz regularity and time continuity do not imply parabolic
Hölder regularity. For example choose a continuous bounded scalar
\(\alpha(t)\) with a logarithmic modulus at an interior time and put
\(b(x,t)=\alpha(t)x\). This drift satisfies the spatial bounds with
\(D^2b=0\), while it need not have any positive Hölder exponent in time.
The bounded one-dimensional solution

\[
 u(x,t)=e^{-\nu\Sigma(t)}\cos(e^{-A(t)}x),\quad
 A(t)=\int_0^t\alpha(s)\,ds,\quad
 \Sigma(t)=\int_0^t e^{-2A(s)}\,ds
\]

shows that its time derivative can inherit that non-Hölder modulus at
positive times. Thus the source's blanket Schauder step is not repaired
by changing cylinders. The Bernstein estimate needs less than that
blanket regularity statement.

A precise stable extension is available: if classical approximants
\(u_j,b_j\) converge locally uniformly to \(u,b\), their drift Lipschitz
constants are uniformly bounded by \(L\), \(b_j(0,t)=0\), and the scalar tails (8) have
uniform \(D,R,p\), then (9) passes to an a.e. weak-gradient bound for
\(u\). To see this, the uniform estimate gives weak-* compactness of
\(\nabla u_j\) on each compact cylinder; integration against a compactly
supported test function identifies every subsequential limit with the
distributional gradient of \(u\). The closed pointwise bound is preserved
under weak-* limits. If \(u\) has a continuous spatial gradient, it holds
everywhere. This extension requires **the approximation/convergence
premise**; the present note does not silently obtain it from (S).

For a general weak/mild interpretation of (S), one must therefore either
provide that standard construction/uniqueness bridge with its solution
class, or cite and check an appropriate regularization theorem. For an
actual smooth NS solution, positive-time differentiability is already
available on its classical existence interval. One must still verify
that its lifted drift has the asserted *global* Lipschitz bound and that
the scalar equation (1) and scalar tail premise hold. Smoothness on
spatial compact sets alone does not prove a global derivative bound, and
no estimate here establishes a uniform bound up to an unknown singular
endpoint.

For orientation, Menozzi–Pesce–Zhang study gradient estimates with
linearly growing drifts using an auxiliary flow in
[their primary paper](https://arxiv.org/abs/2006.07158). No theorem from
that paper is imported into (2)–(12); the proof and constants above are
explicit. Its broader regularity theory should not be cited as an
automatic verification of the draft's exact assumptions.

## 5. Bounded conclusion

The spatial-center issue is repairable: translate with the actual drift,
apply the elementary cutoff estimate, and transfer the scalar tail along
the trajectory. This gives the global gradient bound, vanishing weighted
gradient tail, and \(L^1\) finiteness needed in that part of the argument,
conditional on the precise classical or approximation hypotheses above.
It leaves the separate regularity interpretation, quantitative drift
bootstrap, affine-covariance premise, receiver placement, and S3 assembly
obligations untouched.

## 6. An initial-trace-uniform extension under a representation premise

The positive-time factor \(h^{-1/2}\) is unnecessary when a bounded
initial gradient is available and the solution is the drift–diffusion
evolution of that datum. Here is a separate lemma that makes this premise
explicit and supplies the weighted initial trace.

Let \(u_0\in C_b^1(\mathbb R^n)\), meaning that \(u_0\) and its continuous
first derivatives are bounded. Keep \(\|Db\|_{\mathrm{op}}\le L\) and
\(|b(0,s)|\le B_0\) on \([0,T]\), with \(b\) and \(Db\) continuous.
For each fixed final \(t\), let \(W\) be standard \(n\)-dimensional
Brownian motion, \(\sigma=\sqrt{2\nu}\), and solve

\[
 dX_s=-b(X_s,t-s)\,ds+\sigma\,dW_s,\quad X_0=x,
 \quad 0\le s\le t.                                  \tag{13}
\]

Assume that the chosen solution is represented by

\[
 u(x,t)=\mathbb E[u_0(X_t)].                           \tag{14}
\]

This is a hypothesis identifying the solution; it is not supplied merely
by labeling an unspecified bounded object a solution. For a **bounded
classical** solution, continuous to its initial datum, (14) follows by
applying Itô's formula to \(u(X_s,t-s)\), stopped in a ball and before
the initial-time endpoint. Its drift vanishes by (1). Global Lipschitz
growth gives no explosion; boundedness permits removal of the spatial
stopping and then the initial-time cutoff by dominated convergence.
This argument also gives uniqueness within that bounded classical
class. For a weak/mild solution, (14) or a corresponding proved
uniqueness/approximation bridge must be stated separately. This lemma
does not provide existence for the coupled NS problem.

### Derivative formula without a time-Hölder assumption

Subtracting the additive Brownian path turns (13) into an ordinary
integral equation with a globally Lipschitz, spatially \(C^1\) vector
field. Its derivative with respect to \(x\) obeys the pathwise ODE

\[
 J_s'=-Db(X_s,t-s)J_s,\qquad J_0=I,
 \qquad\|J_s\|\le e^{Ls}.                             \tag{15}
\]

Bounded \(\nabla u_0\) and (15) justify differentiating (14) under the
expectation. Therefore

\[
\nabla u(x,t)=\mathbb E[J_t^T\nabla u_0(X_t)].         \tag{16}
\]

In particular \(\|\nabla u(\cdot,t)\|_\infty
\le e^{Lt}\|\nabla u_0\|_\infty\), uniformly down to time zero.

No Hessian of the drift, derivative of the drift in time, Bismut formula,
or Gaussian covariance ansatz is used. In an application where \(b\)
comes from the actual NS solution, it is that entire actual coefficient
field in (13), rather than a frozen affine replacement.

Let \(\phi_s\) solve the noiseless ODE with the same initial point and
time-reversed drift. With \(S_t=\sup_{s\le t}|W_s|\), pathwise Gronwall
gives \(|X_t-\phi_t|\le\sigma e^{Lt}S_t\). The backward deterministic
flow bound gives \(|x|\le e^{Lt}(|\phi_t|+B_0t)\). Thus

\[
 |x|\le e^{Lt}(|X_t|+Z_t),\qquad
 Z_t=B_0t+\sigma e^{Lt}S_t.                           \tag{17}
\]

### The source's homogeneous quadratic gradient weight

Suppose
\(G_0=\sup_x|x|^2|\nabla u_0(x)|<\infty\), and put
\(B=\|\nabla u_0\|_\infty\). Then
\(|y||\nabla u_0(y)|\le\sqrt{G_0B}\).
Using (15)–(17) and \(\|S_t\|_2\le2\sqrt{nt}\) gives

\[
 \boxed{\quad
 \sup_x|x|^2|\nabla u(x,t)|
 \le e^{3Lt}\left[
   \sqrt{G_0}+
   \left(B_0t+2e^{Lt}\sqrt{2n\nu t}\right)\sqrt B
 \right]^2.
 \quad}                                               \tag{18}
\]

Indeed, before taking the Brownian moment bound, the right side is at
most
\(e^{3Lt}\{G_0+2\sqrt{G_0B}\mathbb EZ_t+B\mathbb EZ_t^2\}
\le e^{3Lt}(\sqrt{G_0}+\sqrt B\|Z_t\|_2)^2\).
For \(b(0,s)=0\), set \(B_0=0\); in dimension five the noise term is
\(2e^{Lt}\sqrt{10\nu t}\).

The bound is uniform for \(0\le t\le T\) and tends to \(G_0\) as
\(t\downarrow0\). For each fixed \(x\), (13) implies \(X_t\to x\) in
probability, while \(\|J_t-I\|\le e^{Lt}-1\). Hence (16) and continuity
of \(\nabla u_0\) give pointwise \(\nabla u(x,t)\to\nabla u_0(x)\).
Choosing a point arbitrarily close to the initial supremum then yields

\[
 \lim_{t\downarrow0}\sup_x|x|^2|\nabla u(x,t)|=G_0.     \tag{19}
\]

This is convergence of the **weighted norm values**. It is not a claim
of convergence of the weighted difference without an additional
uniform-continuity/tail argument. Norm convergence (19) is sufficient
for the initial-supremum passage in the proposed weighted estimate when
its other premises are in place. Smooth compactly supported data satisfy
all the initial assumptions in this lemma; the source's angularly
Lipschitz datum still needs its stated approximation passage.

### Uniform integrable tails from an initial weighted gradient

More generally let \(p\ge1\), \(R>0\), and

\[
 H_{p,0}=\sup_x(1+|x|^2/R^2)^{p/2}|\nabla u_0(x)|<\infty.
\]

For an even integer \(q\ge\max(2,p)\), define

\[
 c_{n,q}=\frac{q}{q-1}
 \left[\prod_{j=0}^{q/2-1}(n+2j)\right]^{1/q}.
\]

Doob's inequality and the Gaussian radial moment give
\(\|S_t\|_q\le c_{n,q}\sqrt t\). Applying the same weight comparison as
in Section 3 to (17), followed by Minkowski in \(L^p\), proves

\[
 \boxed{\quad
 |\nabla u(x,t)|\le H_{p,0}\,e^{(p+1)Lt}
 \left[1+\frac{B_0t+\sqrt{2\nu}\,e^{Lt}c_{n,q}\sqrt t}{R}\right]^p
 \Xi_p(x).
 \quad}                                               \tag{20}
\]

The coefficient tends to \(H_{p,0}\) as \(t\downarrow0\). For \(p>n\)
this gives the gradient's \(L^1\) bound uniformly down to the initial
time by the same integral in (12). For \((n,p,q)=(5,9,10)\),
\(c_{5,10}=(10/9)45045^{1/10}<10/3\), so a convenient explicit version is

\[
 |\nabla u(x,t)|\le H_{9,0}\,e^{10Lt}
 \left[1+\frac{B_0t+(10/3)e^{Lt}\sqrt{2\nu t}}R\right]^9\Xi_9(x).
                                                               \tag{21}
\]

For completeness, the Brownian supremum estimate used here follows
from the nonnegative submartingale \(|W_s|\): stopping at its first
crossing of \(\lambda\) gives
\(\lambda\Pr(S_t\ge\lambda)
\le\mathbb E[|W_t|\mathbf1_{S_t\ge\lambda}]\).
Multiply by \(q\lambda^{q-2}\), integrate, and apply Hölder to obtain
\(\|S_t\|_q\le(q/(q-1))\|W_t\|_q\); truncation removes any preliminary
integrability issue. The radial Gaussian integral gives
\(\mathbb E|W_t|^{2k}=t^k\prod_{j=0}^{k-1}(n+2j)\).

Thus the extra representation premise supplies the formerly missing
initial-trace control for the stated smooth datum class. It does not
repair an unverified drift bound or extend the actual solution beyond
its known classical interval, and does not complete S3.

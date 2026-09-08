# Separate prospective audit of the entropy-kernel route

The proposed kernel estimate is sound under the smooth, bounded,
divergence-free drift hypotheses specified below. An explicit first-exit
argument then gives the proposed energy-dependent logarithmic doubling
delay for a smooth, finite-energy, unforced Navier–Stokes parent.
The kernel must be oriented correctly: in elapsed backward time its
particle drift is **minus** the actual velocity at reversed original
time. This note is a separate Pass 8 analytical result, not part of the
completed Pass 7 Picard argument. It is not a Clay solution, a
terminal-uniform adjoint estimate, or a novelty claim.

## 1. Kernel orientation and hypotheses

Work on physical $\mathbb R^3$. Let $u(r,y)$ be a smooth bounded
divergence-free drift on a fixed slab, and write $\mathcal U(t,s)$
for the scalar evolution of

\[
\partial_t v+u(t)\cdot\nabla v-\nu\Delta v=0.
\]

Its positive kernel is defined by

\[
[\mathcal U(t,s)g](x)=\int K(t,x;s,y)g(y)\,dy.
\tag{1}
\]

The original-time generator on scalar functions is
$L_s=\nu\Delta-u(s)\cdot\nabla$. For fixed terminal $(t,x)$,
the source-variable equation is

\[
-\partial_sK=L_s^*K
=\nu\Delta_yK+\nabla_y\cdot(u(s,y)K).
\]

Consequently, for $\tau=t-s$, the density
$p(\tau,y)=K(t,x;t-\tau,y)$ solves

\[
\partial_\tau p=\nu\Delta p+
\nabla\cdot(u(t-\tau,y)p),\qquad p(0)=\delta_x.
\tag{2}
\]

It is the Fokker–Planck density of

\[
dY_\tau=-u(t-\tau,Y_\tau)\,d\tau+\sqrt{2\nu}\,dW_\tau,
\qquad Y_0=x.
\tag{3}
\]

Thus the sign is $-u$, and the coefficient time is reversed. Treating
this as the forward physical tracer with drift $+u(s)$ would give the
wrong kernel for (1). Boundedness supplies nonexplosion, and
$\int K(t,x;s,y)dy=1$. For a smooth actual NS parent these statements
are used only on slabs preceding its maximal smooth time.

## 2. The kernel bounds are independent of strain

Let $\|u(r)\|_\infty\le U$ throughout $[s,t]$. Then

\[
0\le K(t,x;s,y)\le C(\nu\tau)^{-3/2},
\qquad \tau=t-s>0,
\tag{4}
\]

with $C$ dimensional and independent of $U$ or derivatives of $u$.
Here incompressibility is essential. For the forward equation and its
time-reversed adjoint, drift work vanishes from the $L^2$ identity,
and $L^1$ is contractive. The Nash inequality

\[
\|v\|_2^{10/3}\le C\|\nabla v\|_2^2\|v\|_1^{4/3}
\]

then yields the $L^1\to L^2$ bound $C(\nu\tau)^{-3/4}$.
Duality for the time-dependent adjoint gives the same $L^2\to
L^\infty$ bound. Composition at the midpoint gives (4). This argument
does not invoke an Aronson constant depending exponentially on drift.
The classical Nash inequality originates in
[Nash's 1958 paper](https://www.karlin.mff.cuni.cz/~kaplicky/pages/pages/2011z/Nash1958.pdf);
the cancellation and nonautonomous application above are the argument
needed here.

From (3) and $|\int_0^\tau u(t-a,Y_a)da|\le U\tau$,

\[
\int|y-x|^2K(t,x;s,y)dy
\le2U^2\tau^2+12\nu\tau.
\tag{5}
\]

The diffusion contribution is $2\nu\mathbb E|W_\tau|^2=6\nu\tau$
before the elementary factor two used in (5). Both (4) and (5) apply
to the correctly oriented kernel; no pointwise Gaussian tail bound
with a drift-independent center was assumed.

## 3. Gaussian exponential averaging of BMO oscillation

Assume that at the source time, $f=\nabla u(s)$ and that the physical
velocity and vorticity have upper bounds

\[
\|u(s)\|_2^2\le E,\qquad
\|\operatorname{curl}u(s)\|_\infty\le M,
\qquad E,M>0.
\]

Use the Frobenius magnitude for matrices. The Calderón–Zygmund
oscillation bound is $[f]_{\mathrm{BMO}}\le CM$. For a ball
$B=B(x,L)$ let $m_B=f_B$, and let

\[
G_L(y-x)=(4\pi L^2)^{-3/2}e^{-|y-x|^2/(4L^2)}.
\]

There are dimensional $a,C>0$ such that

\[
\int G_L(y-x)
\exp\!\left(\frac{a|f(y)-m_B|}{M}\right)dy\le C.
\tag{6}
\]

This is a consequence of the ordinary ball John–Nirenberg inequality,
not a claim that arbitrary normalized weights preserve BMO estimates.
On $B(x,2^{j+1}L)$, the mean shifts by at most $CM(j+1)$ from
$m_B$. John–Nirenberg controls the exponential of the remaining
oscillation if $a$ is sufficiently small. Multiplication by the
Gaussian bound on the annulus leaves summable factors of the form
$C2^{3j}e^{Caj}e^{-c4^j}$. The inner ball is handled by the same
inequality. For a matrix, apply the scalar inequality to its finitely
many components and decrease $a$ as needed. The classical input is
[John–Nirenberg's exponential oscillation inequality](https://onlinelibrary.wiley.com/doi/10.1002/cpa.3160140317).

The earlier mean-retaining heat argument also gives, uniformly in
$x$,

\[
|m_B|\le CM[1+\log_+(\ell_E/L)],\qquad
\ell_E=(E^{1/2}/M)^{2/5}.
\tag{7}
\]

Indeed $|m_B-e^{L^2\Delta}f(x)|\le CM$, and the signed heat
average is bounded using energy at length $\ell_E$ and integrating
$\nabla\operatorname{curl}e^{a\Delta}\omega$ between the two
scales. This retains the large-scale mean that failed in the scalar
seminorm-only argument.

## 4. Relative entropy gives the proposed average

For probability densities $p$ and $G_L$, define
$D(p\Vert G_L)=\int p\log(p/G_L)$. Equations (4)–(5) imply

\[
D(p\Vert G_L)
\le C+3\log\frac{L}{\sqrt{\nu\tau}}
 +C\frac{\nu\tau+U^2\tau^2}{L^2}.
\tag{8}
\]

For example, use $\log p\le\log\|p\|_\infty$ and the exact
quadratic logarithm of $G_L$. The negative part in the entropy
definition is integrable, and the displayed upper bound is finite;
one need not assume a separately finite differential entropy.

For any nonnegative measurable $F$, the elementary entropy inequality
is

\[
\int pF\le D(p\Vert G_L)+\log\int G_Le^F.
\tag{9}
\]

It follows from nonnegativity of relative entropy against the tilted
probability density $G_Le^F/\int G_Le^F$; truncation justifies it
for an unbounded $F$. Apply this to $F=a|f-m_B|/M$ and use (6).

Now suppose the energy and vorticity upper caps $E,M$ hold throughout
the intervening drift slab. The independently checked interpolation
estimate gives

\[
U\le CE^{1/5}M^{3/5}=CM\ell_E.
\tag{10}
\]

Choose $L=\max(\ell_E,\sqrt{\nu\tau})$. Then (7) costs at
most $CM$, the diffusion moment term in (8) is bounded, and
$U^2\tau^2/L^2\le C(M\tau)^2$. Thus

\[
\boxed{\displaystyle
\sup_x\int K(t,x;s,y)|\nabla u(s,y)|\,dy
\le CM\left[1+\log_+\frac{\ell_E}{\sqrt{\nu(t-s)}}
 +(M(t-s))^2\right].}
\tag{11}
\]

The drift cap in (10) is needed throughout $[s,t]$, whereas the
gradient BMO and mean estimates are at source time $s$. Conflating
those two quantifiers would be a gap. Fixed caps on the whole slab
supply both. Smooth finite-energy NS solutions have the required
regularity; the assertions here do not address an arbitrary rough
drift coupled to a merely distributional vector solution.

## 5. A noncircular first-doubling bootstrap

Let a smooth finite-energy unforced NS parent start at $t_0$ with
$E_0=\|u(t_0)\|_2^2>0$ and
$M_0=\|\omega(t_0)\|_\infty>0$. Define

\[
\ell_0=(E_0^{1/2}/M_0)^{2/5},\qquad
\mathrm{Re}_E=\frac{M_0\ell_0^2}{\nu},\qquad
A_0=1+\log_+\mathrm{Re}_E.
\tag{12}
\]

Take the first-exit interval on which the solution is smooth and
$\|\omega\|_\infty<2M_0$. The energy inequality supplies
$E\le E_0$, and (10) supplies a uniform velocity cap there. In
(11), take the cap $M=2M_0$; its associated energy length is
$2^{-2/5}\ell_0$, which may be enlarged to $\ell_0$ in the
upper bound, with universal constants unchanged in nature.

The actual vector equation has the scalar-kernel Duhamel formula

\[
\omega(t,x)=\int K(t,x;t_0,y)\omega(t_0,y)dy
 +\int_{t_0}^t\int K(t,x;s,y)
 [(\omega(s,y)\cdot\nabla)u(s,y)]\,dy\,ds.
\tag{13}
\]

The kernel in (13) incorporates the actual advection. Its positivity
and unit mass control the initial vector magnitude by $M_0$.
Using the bootstrap cap on the factor $|\omega|$ and (11), for
$T=t-t_0$ within the first-exit interval,

\[
\|\omega(t)\|_\infty
\le M_0+CM_0^2T
\left[1+\log_+\frac{\ell_0}{\sqrt{\nu T}}
 +(M_0T)^2\right].
\tag{14}
\]

This uses the same integrable logarithm calculation as the pure-heat
case and integrates the squared-time term directly.

Choose

\[
H=\frac{c}{M_0A_0},\qquad 0<c<1.
\tag{15}
\]

To check the constants, write $a=M_0H=c/A_0$ and
$r=\log_+\mathrm{Re}_E$. Then

\[
\log_+\frac{\ell_0}{\sqrt{\nu H}}
=\tfrac12\log_+\frac{\mathrm{Re}_E}{a}
\le\tfrac12\left[r+\log(A_0/c)\right].
\]

Since $\log A_0\le A_0$ and $A_0\ge1$, the increment in (14),
divided by $M_0$, is at most
$C[c(1+\log(1/c))+c^3]$. A sufficiently small universal $c$
makes this at most $1/2$. The integral bound is increasing in its
terminal time, so the same improvement holds for all $T\le H$.
Consequently the bootstrap improves $2M_0$ to $3M_0/2$.

A first vorticity exit before $H$ is excluded by continuity. If the
maximal smooth interval ended earlier without such an exit, the
already-established energy/vorticity interpolation would give a
uniform velocity supremum. Restarting the standard velocity mild
theory at times approaching that endpoint gives a uniform positive
extension time, and uniqueness continues the actual solution. This
contradicts maximality. The kernel was only used on already-smooth
subslabs; no future cap was assumed to construct it.

Thus this route proves the scoped delay

\[
\boxed{\displaystyle
\sup_{0\le T\le c/[M_0(1+\log_+\mathrm{Re}_E)]}
\|\omega(t_0+T)\|_\infty\le\frac32M_0.}
\tag{16}
\]

In particular, the first subsequent doubling cannot occur sooner.
The scale check is consistent: under the original NS scaling,
$\ell_0\mapsto\lambda^{-1}\ell_0$,
$\mathrm{Re}_E$ is invariant, and this delay scales as
$\lambda^{-2}$.

## 6. Scope of the positive audit

The argument depends on physical incompressibility, finite energy,
bounded vorticity on the bootstrap slab, and the actual scalar
advection-diffusion kernel. It applies neither to the five-dimensional
lift with nonzero divergence and reaction nor to its weighted adjoint
without further analysis. It does not justify an unqualified $c/M_0$
clock independent of the energy Reynolds number.

The delay remains a local statement whose size changes with the
current norms. It supplies no bound on those norms up to a possible
singular time, no prescribed detector ancestry, no repeated-epoch
estimate, and no resolution of the Clay problem. The proof is written
separately so that its classical analytical inputs and kernel
orientation can receive further scrutiny before any later promotion.

## 7. Final read of the independently organized root proof

I read `work/pass8/LOGARITHMIC_RECORD_CLOCK.md` after the preceding
derivation and the independent nonautonomous kernel control. I found
no substantive estimate or scope discrepancy. This final read checks
the mathematical argument; it does not independently re-audit the
additional source agent's literature comparisons or execute root code.

The direct Gaussian bound on $\exp(a|\nabla u|/M)$ at $L\ge\ell_E$
is valid: the previously retained ball mean is $O(M)$, so restoring it
only multiplies the exponential oscillation integral by a dimensional
constant. It is equivalent to the mean-subtracted route used above.

The entropy wording is justified. The negative part of
$p\log(p/G_L)$ integrates to at most $1/e$ after expressing it as
$G_L r\log r$. Bounded $p$ controls its positive logarithmic part,
and the finite second moment controls the quadratic term in
$-\log G_L$. Thus the upper bound and the truncated entropy
variational argument involve finite quantities; no unproved lower
Gaussian bound for $p$ is needed.

The upper bound used throughout $0<T\le H$ is increasing. On the
branch $T<b^2$,

\[
\frac d{dT}\{T[1+\log(b/\sqrt T)]\}
=\frac12+\log(b/\sqrt T)>0,
\]

and on $T>b^2$ its derivative is 1. The cubic term is also increasing.
Alternatively, monotonicity follows before estimating the positive
time integral. Hence evaluating the bootstrap increment at $H$ does
control every preceding time.

The $H^m$ continuation is a valid use of bounded-velocity mild theory,
with its standard persistence step understood as follows. The restart
keeps finite $L^2$ energy and is smooth at positive elapsed times.
Choose a positive elapsed time still before the hypothesized old
endpoint, where uniqueness identifies it with the original $H^m$
solution. On the remaining interval away from the restart, the mild
solution has bounded spatial derivatives, in particular a finite
integral of $\|\nabla u\|_\infty$. The usual commutator energy
estimate

\[
\frac d{dt}\|u\|_{H^m}^2
\le C_m\|\nabla u\|_\infty\|u\|_{H^m}^2
\]

with the dissipative term retained on the left if desired, propagates
the same $H^m$ solution beyond the old endpoint. Thus the restart
does not merely substitute a lower-regularity bounded trajectory for
a maximal strong solution. This is a persistence/uniqueness use of
the cited local theory, not an extra norm cap derived from entropy.

Finally, both declared comparisons are correct. Replacing the current
energy by a fixed initial upper cap only shortens the guaranteed
clock. At attained levels $M_j=2^jM_0$, the resulting lower delays
are of order $2^{-j}/(1+j)$ at large $j$ and have finite sum. They
therefore cannot exclude finite-time record accumulation. In the
specified normalization $R(M)=M^{-11/20}$ and at fixed energy cap
and viscosity, the ratio of this guaranteed clock to $R(M)^2/\nu$
is of order $\nu M^{1/10}/(1+\log_+(CM^{1/5}))$, which tends to
infinity. This is only a duration comparison and entails none of
the occupation, ancestry, overlap, or adjoint conclusions expressly
left open in the root proof.

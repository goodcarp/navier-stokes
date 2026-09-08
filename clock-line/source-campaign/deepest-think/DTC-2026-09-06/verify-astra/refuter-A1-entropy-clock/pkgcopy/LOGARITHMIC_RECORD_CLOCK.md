# An energy-dependent logarithmic record clock

The Clay problem remains unresolved. The result below is a local
bound for an actual smooth solution of the unforced three-dimensional
Navier–Stokes equation on physical $\mathbb R^3$. It does not assert
a uniform bound through a possible singular time. The proof uses
classical kernel, BMO, and entropy tools; no novelty is claimed.

## 1. Statement and conventions

Let $u$ be the maximal strong solution with smooth finite-energy
initial data, for example divergence-free $H^m$ data with integer
$m\ge4$, and viscosity $\nu>0$. At a time $t_0$ before its maximal
smooth endpoint, put

\[
E_0=\|u(t_0)\|_2^2,\quad M_0=\|\omega(t_0)\|_\infty,
\quad \ell_0=E_0^{1/5}M_0^{-2/5},\quad
\mathcal R_0=\frac{M_0\ell_0^2}{\nu}
=\frac{E_0^{2/5}M_0^{1/5}}{\nu}.
\tag{1}
\]

There is a dimensional constant $c>0$ such that the solution
continues through the interval of length

\[
\boxed{H_0=\frac{c}{M_0(1+\log_+\mathcal R_0)}}
\tag{2}
\]

and satisfies

\[
\sup_{0\le T\le H_0}\|\omega(t_0+T)\|_\infty
\le\tfrac32M_0.
\tag{3}
\]

Here $\log_+r=\max(\log r,0)$. If $E_0=0$ or $M_0=0$, finite
energy and incompressibility imply $u=0$, so that case is trivial
and (1) need not be evaluated. Energy has no factor of one half.
Vorticity uses Euclidean magnitude and gradients use Frobenius
magnitude. No numerical best constant is claimed.

## 2. Three velocity facts with quantitative means retained

Suppose a smooth divergence-free snapshot has
$\|u\|_2^2\le E$ and $\|\operatorname{curl}u\|_\infty\le M$,
with positive caps $E,M$. Set $\ell=E^{1/5}M^{-2/5}$.

First, the heat split

\[
u=e^{L^2\Delta}u+
 \int_0^{L^2}\operatorname{curl}e^{a\Delta}\omega\,da
\]

gives $\|u\|_\infty\le C(E^{1/2}L^{-3/2}+ML)$.
At $L=\ell$ this is

\[
\|u\|_\infty\le CE^{1/5}M^{3/5}=CM\ell.
\tag{4}
\]

Second, the Calderón–Zygmund representation gives
$[\nabla u]_{\mathrm{BMO}}\le CM$. This is a homogeneous
oscillation estimate, not a bound on a local average.

Third, for $f=\nabla u$, comparison of its Gaussian average at
radius $L$ with its ball mean $f_{B(x,L)}$ yields

\[
|f_{B(x,L)}-e^{L^2\Delta}f(x)|\le CM.
\tag{5}
\]

To see this, telescope ball averages on annuli of radii $2^jL$.
The oscillation grows at most as $C(j+1)M$ and the Gaussian
weights sum as $\sum_j2^{3j}(j+1)e^{-c4^j}<\infty$.
The signed heat average obeys
$\|\nabla e^{L^2\Delta}u\|_\infty\le CE^{1/2}L^{-5/2}$.
Thus, whenever $L\ge\ell$,

\[
|f_{B(x,L)}|\le CM.
\tag{6}
\]

These are the same mean-retaining tools used in the preceding
pure-heat bound. Only (6), at large enough radius, is needed below.
Qualitative decay or finite enstrophy has not replaced its
quantitative energy input.

The ball John–Nirenberg inequality and (6) imply the following
Gaussian exponential estimate, for dimensional $a,C>0$:

\[
\boxed{\int G_L(y-x)
 \exp\left(\frac{a|\nabla u(y)|}{M}\right)dy\le C
 \quad\text{for every }L\ge\ell,}
\quad
G_L(z)=(4\pi L^2)^{-3/2}e^{-|z|^2/(4L^2)}.
\tag{7}
\]

For clarity, apply John–Nirenberg to each of the finitely many
matrix components and then Hölder's inequality, decreasing $a$.
On a ball of radius $2^{j+1}L$, the mean shift from $f_{B(x,L)}$
is at most $CM(j+1)$. Multiplying the resulting exponential
integral by the Gaussian bound on the annulus leaves summable
factors $C2^{3j}e^{Caj}e^{-c4^j}$. The inner ball follows directly
from John–Nirenberg. Finally (6) restores the mean at a bounded
multiplicative cost. This proof is specific to the Gaussian
reference measure; no arbitrary-weight BMO assertion is used.

## 3. The actual scalar transport kernel

On a fixed smooth slab let $\mathcal U(t,s)$ solve
$\partial_t v+u(t,x)\cdot\nabla v=\nu\Delta v$, and define

\[
[\mathcal U(t,s)g](x)=\int K(t,x;s,y)g(y)\,dy.
\tag{8}
\]

For fixed $(t,x)$, set $\tau=t-s$ and
$p(\tau,y)=K(t,x;t-\tau,y)$. The source-variable equation is

\[
\partial_\tau p=\nu\Delta_y p+
 \nabla_y\cdot(u(t-\tau,y)p),\qquad p(0)=\delta_x.
\tag{9}
\]

Equivalently, this is the density of
$dY_\tau=-u(t-\tau,Y_\tau)d\tau+\sqrt{2\nu}\,dW_\tau$,
with $Y_0=x$. Both the minus sign and the reversed coefficient
time matter. Positivity and bounded drift imply nonexplosion and
$\int p\,dy=1$. Physical incompressibility also gives the other
kernel mass identity, by integrating the original scalar equation.

There is a dimensional, drift-independent bound

\[
\|p(\tau)\|_\infty\le C(\nu\tau)^{-3/2}.
\tag{10}
\]

Here is the required argument, so that no hidden drift-dependent
Gaussian constant enters. For both scalar evolution and its
reversed adjoint, incompressibility cancels the drift in the
$L^2$ energy identity. Both are $L^1$ contractions. The Nash
inequality
$\|v\|_2^{10/3}\le C\|\nabla v\|_2^2\|v\|_1^{4/3}$
then yields the $L^1\to L^2$ norm $C(\nu\tau)^{-3/4}$.
Duality gives the same $L^2\to L^\infty$ bound for the other
half of the interval. Composition at the midpoint proves (10).
Smooth approximations of a point mass justify the kernel version.
The Nash inequality itself follows by splitting its Fourier
integral at a radius and optimizing the low-frequency $L^1$ and
high-frequency gradient bounds.

As a primary-source check, a stronger concentration comparison for
bounded measurable, time-dependent divergence-free drifts implies
(10), even with the heat-kernel supremum constant; see
[Hess-Childs–Raquépas–Rowan, Theorem 1.7 and Corollary 1.10](https://arxiv.org/pdf/2503.16723),
pp. 2–3. We only need the elementary dimensional bound above.
[Qian–Xi](https://arxiv.org/pdf/1612.07727), pp. 3–4 and 15,
records the adjoint kernel framework and both mass identities.

If $\sup_{s\le r\le t}\|u(r)\|_\infty\le U$, the stochastic
representation in (9) gives the separate moment bound

\[
\int|y-x|^2p(\tau,y)\,dy\le2U^2\tau^2+12\nu\tau.
\tag{11}
\]

This uses $|\int_0^\tau u(t-r,Y_r)dr|\le U\tau$,
$\mathbb E|W_\tau|^2=3\tau$, and
$|v+w|^2\le2|v|^2+2|w|^2$. It is not a drift-independent
pointwise Gaussian tail estimate centered at $x$.

## 4. Relative entropy pays for transport

Assume the snapshot caps $E,M$ hold throughout $[s,t]$, and put
$\ell=E^{1/5}M^{-2/5}$. They bound the drift by (4) on that whole
slab and bound the source gradient $\nabla u(s)$ by (7). These are
different time quantifiers, both supplied by the same bootstrap.

Choose $L=\max(\ell,\sqrt{\nu\tau})$ and use $G_L(y-x)$ as
a reference probability density. Equations (10)–(11) give

\[
\begin{aligned}
D(p\Vert G_L)
&=\int p\log(p/G_L)\\
&\le C+3\log\frac{L}{\sqrt{\nu\tau}}
 +C\frac{\nu\tau+U^2\tau^2}{L^2}\\
&\le C\left[1+\log_+\frac{\ell}{\sqrt{\nu\tau}}
 +(M\tau)^2\right].
\end{aligned}
\tag{12}
\]

Indeed $\log p\le\log\|p\|_\infty$ and $-\log G_L$ is
an explicit quadratic. Bounded density and the finite second
moment justify the entropy integrals; the negative part of
$p\log(p/G_L)$ is integrable since $r\log r\ge-1/e$ relative
to $G_L$. In the last line, $U/\ell\le CM$ and $L\ge\ell$.

For any nonnegative $F$, the entropy inequality is

\[
\int pF\le D(p\Vert G_L)+\log\int G_Le^F.
\tag{13}
\]

It follows from nonnegativity of entropy relative to the tilted
probability density $G_Le^F/\int G_Le^F$; first truncate $F$ and
then pass monotonically. Apply (13) with
$F=a|\nabla u(s)|/M$ and use (7). We obtain the actual-kernel
estimate

\[
\boxed{
\sup_x\int K(t,x;s,y)|\nabla u(s,y)|\,dy
\le CM\left[1+\log_+\frac{\ell}{\sqrt{\nu(t-s)}}
 +(M(t-s))^2\right].}
\tag{14}
\]

The kernel is generated by the actual advecting velocity. Formula
(14) has been derived on an already smooth capped slab; it has
not assumed that such a slab continues indefinitely.

## 5. First exit and continuation

Start at $t_0$ and consider times before both the maximal smooth
endpoint and the first exit from $\|\omega\|_\infty<2M_0$.
Energy supplies the cap $E_0$. On each such slab, apply (14) with
$M=2M_0$. Its length $\ell=2^{-2/5}\ell_0$ may be enlarged to
$\ell_0$ in the upper bound, changing only dimensional constants.

The original vorticity equation has the exact scalar-kernel
Duhamel representation

\[
\omega(t)=\mathcal U(t,t_0)\omega(t_0)
 +\int_{t_0}^t\mathcal U(t,s)[(\omega(s)\cdot\nabla)u(s)]\,ds.
\tag{15}
\]

The advection is already included in $\mathcal U$. Positivity,
unit mass, the cap $|\omega|<2M_0$, and (14) yield, for
$T=t-t_0$ in the first-exit interval,

\[
\|\omega(t)\|_\infty
\le M_0+CM_0^2T\left[1+\log_+\frac{\ell_0}{\sqrt{\nu T}}
 +(M_0T)^2\right].
\tag{16}
\]

The time integration uses
$\int_0^T\log_+(\ell_0/\sqrt{\nu a})da
\le T\log_+(\ell_0/\sqrt{\nu T})+T/2$ and
$\int_0^T a^2da=T^3/3$.

Let $A=1+\log_+\mathcal R_0$ and $H=c/(M_0A)$ with $0<c<1$.
At that time,

\[
\log_+\frac{\ell_0}{\sqrt{\nu H}}
=\tfrac12\log_+\frac{\mathcal R_0 A}{c}
\le\tfrac12\left[\log_+\mathcal R_0+\log A+\log(1/c)\right].
\]

Since $\log A\le A$, the increment in (16), divided by $M_0$,
is at most $C[c(1+\log(1/c))+c^3]$. Choose a universal $c$
small enough to make it at most $1/2$. The integrated bound is
increasing with $T$ (also directly check the derivative of
$T[1+\log_+(b/\sqrt T)]$). Hence this improves the cap to
$3M_0/2$ for every $T\le H$ still in the first-exit interval.
Continuity excludes a first exit at or before $H$.

If a maximal smooth endpoint came first, (4) would give a uniform
velocity supremum on the interval approaching it. The standard
bounded-velocity mild theory supplies a positive restart duration
$c_1\nu/\|u(s)\|_\infty^2$, uniformly bounded below under
that cap. Restarting sufficiently close to the endpoint and using
uniqueness extends the same solution beyond it. This contradicts
maximality. More explicitly, the mild restart retains finite $L^2$
norm and has derivative bounds at positive times. Away from the
restart instant these bounds and the usual $H^m$ energy estimate
preserve the original strong regularity across the former endpoint.
The input is the checked local theory of
[Giga–Inui–Matsui](https://eprints.lib.hokudai.ac.jp/repo/huscap/all/69160/re410.pdf).
Thus the actual solution reaches $H$ and (2)–(3) follow. No future
cap was used to construct a kernel before this continuation step.

## 6. What this improves, and what remains missing

The preceding conservative clock was proportional to
$1/(M_0\mathcal R_0)$. Formula (2) is longer at large
$\mathcal R_0$. Under original Navier–Stokes scaling,
$E_0\mapsto\lambda^{-1}E_0$, $M_0\mapsto\lambda^2M_0$,
$\ell_0\mapsto\lambda^{-1}\ell_0$, and $\mathcal R_0$ is
invariant. The duration scales as $\lambda^{-2}$, as required.

For actual attained records $M_j=2^jM_0$, restarting the argument
gives a lower bound on the next doubling delay. With the fixed
initial energy as an upper cap, these lower bounds have the form
$c/[2^jM_0(1+\log_+(C2^{j/5}))]$. Their sum is finite. Nothing
here prevents infinitely many records from accumulating at a
finite time, nor does it establish a dissipation charge for them.

In the previously considered normalization $R(M)=M^{-11/20}$,
fixed $E_0,\nu$ give
$H(M)/(R(M)^2/\nu)\asymp\nu M^{1/10}/(1+\log_+(C M^{1/5}))$
as $M\to\infty$, which tends to infinity. This addresses a
duration mismatch only. It provides no lower swirl occupation,
record ancestry, cutoff overlap, or terminal-uniform adjoint trace.

The physical divergence-free kernel facts do not apply unchanged
to the five-dimensional lift with nonzero divergence and a
reaction term. Its earlier gauge obstruction remains. The result
also retains a logarithmic energy dependence and makes no claim
of a uniform $c/M_0$ clock.

The [independent analytical audit](ENTROPY_KERNEL_AUDIT.md)
checks the same argument. The nonautonomous scalar time control
tests both reversed-time choices separately. Algebraic checks do
not prove John–Nirenberg, Nash, entropy inequalities, or PDE
continuation; the proof above and the stated classical inputs are
essential. No current twin material was used.

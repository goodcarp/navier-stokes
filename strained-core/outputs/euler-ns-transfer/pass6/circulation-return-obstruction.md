# A circulation obstruction to repeated normalized rotating-core returns

**Result.** The proposed supercritical return cannot repeat indefinitely while retaining a uniformly nontrivial normalized rotating core. The scalar circulation maximum principle gives a decreasing quantity under every proposed return:
\[
G(V):=\|\rho V^\theta\|_\infty,\qquad
\boxed{G(V_{n+1})\le q_n^\alpha G(V_n),\quad\alpha>0.}
\tag{1}
\]
If the admitted states have any fixed positive lower bound \(G(V)\ge g_*>0\), this rules out an infinite sequence whose cumulative radius tends to zero. A fixed-radius rotational receiver lower bound is enough; an exactly affine template is not required.

This conclusion is about a particular swirl-bearing return class. It is not a global-regularity proof for axisymmetric Navier–Stokes, nor an exclusion of arbitrary meridional or RMS concentration scenarios. The previously certified finite episode remains valid.

## 1. Published maximum principle and its applicable interval

For an axisymmetric velocity \(u=u^\rho e_\rho+u^\theta e_\theta+u^ze_z\), put \(\Gamma=\rho u^\theta\). Lei and Zhang state its scalar equation and maximum principle as equations (1-3)–(1-4), page 171, in [Criticality of the axially symmetric Navier–Stokes equations, Pacific Journal of Mathematics 289 (2017)](https://msp.org/pjm/2017/289-1/pjm-v289-n1-p06-s.pdf). The [authors' arXiv version](https://arxiv.org/abs/1505.02628) states the same bound. With an explicit viscosity coefficient, the identities are
\[
\partial_t\Gamma+u^\rho\partial_\rho\Gamma+u^z\partial_z\Gamma
=\nu\left(\partial_{\rho\rho}\Gamma-\frac1\rho\partial_\rho\Gamma
                         +\partial_{zz}\Gamma\right),
\qquad
\|\Gamma(t)\|_\infty\le\|\Gamma(0)\|_\infty.
\tag{2}
\]
We use this only on intervals where the actual solution is smooth. It does not assume that the solution stays smooth for all physical time.

For completeness, multiplying the azimuthal NS equation
\[
\partial_tu^\theta+u^\rho\partial_\rho u^\theta
 +u^z\partial_z u^\theta+\frac{u^\rho u^\theta}{\rho}
=\nu\left(\partial_{\rho\rho}u^\theta+\frac1\rho\partial_\rho u^\theta
 +\partial_{zz}u^\theta-\frac{u^\theta}{\rho^2}\right)
\]
by \(\rho\) gives (2). Smoothness on the symmetry axis implies \(\Gamma=O(\rho^2)\), hence zero axis boundary value. At an interior maximum its gradient vanishes and its diffusion is nonpositive. This supplies the scalar comparison argument.

There is no unnoticed assumption of compact support at positive time. On each closed smooth interval of the present finite-energy solution, \(u\) is bounded, and \(|\Gamma|\) grows at most linearly at spatial infinity. A barrier
\(\epsilon e^{ct}(1+\rho^2+z^2)\) dominates it there. The drift is bounded, and the diffusion operator in (2) maps \(1+\rho^2+z^2\) to \(2\), so \(c\) can be chosen to make this a strict supersolution. Apply comparison to both signs, then let \(\epsilon\downarrow0\). The compact initial datum has finite \(G(u_0)\).

## 2. Exact normalized no-return theorem

Let \(V_n\) be complete axisymmetric incoming velocities with finite \(G(V_n)\). Evolve each by its actual unforced NS equation, at any positive viscosity \(\nu_n\), for a smooth Sobolev interval ending at \(\tau_n\), with bounded velocity on that closed interval. The actual \(H^{10}\) solutions in the prior passes meet this hypothesis. Suppose the exact successor is
\[
V_{n+1}(y)=q_n^{1+\alpha}u_n(\tau_n,q_ny),
\qquad 0<q_n<1,\quad\alpha>0.
\tag{3}
\]
The whole field is retained. No swirl is added or replaced between stages. An axial translation is harmless, but the common symmetry axis must be preserved.

At corresponding points \(\rho_x=q_n\rho_y\),
\[
\rho_y V_{n+1}^\theta(y)
=q_n^\alpha[\rho_xu_n^\theta(\tau_n,x)].
\]
Taking the global supremum and applying (2) proves (1). This holds for every positive \(\nu_n\), including the varying normalized viscosities
\(\nu_{n+1}=q_n^\alpha\nu_n\).

Let \(Q_n=\prod_{j=0}^{n-1}q_j\). Induction gives
\[
\boxed{G(V_n)\le Q_n^\alpha G(V_0).}
\tag{4}
\]
If all admitted states satisfy \(G(V)\ge g_*>0\), then every admitted \(n\) obeys
\[
\boxed{
Q_n^\alpha\ge\frac{g_*}{G(V_0)},\qquad
\alpha\sum_{j=0}^{n-1}\log(1/q_j)
\le\log\frac{G(V_0)}{g_*}.}
\tag{5}
\]
Thus no infinite admitted sequence can have \(Q_n\to0\). In particular, if \(q_n\le q_+<1\), its number \(N\) of successful returns satisfies
\[
\boxed{N\le
\left\lfloor\frac{\log[G(V_0)/g_*]}
                   {\alpha\log(1/q_+)}\right\rfloor.}
\tag{6}
\]
If \(G(V_0)<g_*\), the initial state is not admitted; if \(G(V_0)=g_*\), even one strict shrink cannot remain admitted. A class need not have a uniform upper bound on \(G\) for this theorem: the finite value of the actual initial state is enough.

For variable positive exponents \(\alpha_j\), the exact bound instead concerns
\(\sum_j\alpha_j\log(1/q_j)\). A conclusion from \(Q_n\to0\) alone requires a fixed exponent or a positive lower bound on the exponents. Infinite episodes with \(q_j\to1\) and \(Q_n\) bounded away from zero are not excluded by (5), but they do not realize indefinite spatial shrinking.

## 3. Minimal quantitative rotation assumptions

Fix a nonnegative nonzero radial receiver weight
\(\chi_R(x)=\chi(|x|^2/R^2)\), of finite positive mass and second moment. Define
\[
\Omega_R(V)=
\frac{\int\chi_R V\cdot Jx\,dx}{\int\chi_R|Jx|^2\,dx},
\qquad Jx=(-x_2,x_1,0).
\]
Because \(V\cdot Jx=\rho V^\theta=\Gamma_V\), this is exactly
\[
\Omega_R(V)=\frac{\int\chi_R\Gamma_V}{D_R},
\qquad D_R=\int\chi_R\rho^2.
\]
Let
\[
d_\chi=\frac{\int\chi_1\rho^2}{\int\chi_1}>0.
\]
Then
\[
\boxed{R^2|\Omega_R(V)|\le \frac{G(V)}{d_\chi},
\qquad G(V)\ge d_\chi R^2|\Omega_R(V)|.}
\tag{7}
\]
No sign assumption on the actual swirl is needed. Positivity of the receiver weight and the triangle inequality suffice.

Therefore a uniform receiver condition
\[
R_n\ge R_*>0,\qquad|\Omega_{R_n}(V_n)|\ge\omega_*>0
\]
supplies \(g_*=d_\chi R_*^2\omega_*\). This is weaker than requiring return to an entire core profile.

Alternatively, central rotation and uniform local shape control suffice. Suppose every state has \(\Omega(V,0)\ge\omega_*>0\), and
\(|\partial_{x_1x_1}V_2|\le M\) on a fixed ball \(B_{R_*}\). At \(x=(a,0,0)\),
\[
V^\theta(a,0)=V_2(a,0)\ge\omega_*a-\frac12Ma^2.
\]
Choose
\[
a=\min\{R_*/2,\omega_*/M\}>0,
\]
omitting the second restriction if \(M=0\). Then
\[
\boxed{G(V)\ge\frac{\omega_*a^2}{2}>0.}
\tag{8}
\]
A fixed local \(C^2\) bound follows from the proposed bounded template parameter class and a fixed weighted \(H^{10}\) residual bound, since the weight is bounded below on a fixed ball. The pass5 proposal also explicitly demands local derivative bounds. Consequently adding a bounded cubic-core coordinate or bounded residual does not evade the obstruction while the normalized central rotation remains uniformly positive.

Central rotation **alone** is insufficient. The compact solenoidal radial swirls
\[
V_\epsilon(x)=\chi(|x|^2/\epsilon^2)Jx,\qquad \chi=1\text{ near }0,
\]
all satisfy \(\Omega(V_\epsilon,0)=1\), but \(G(V_\epsilon)=O(\epsilon^2)\). Their core radius shrinks and their higher local derivative norms are not uniform. The theorem must retain a quantitative profile radius, a receiver lower bound, or equivalent derivative control.

## 4. A global ceiling on the rotational receiving Reynolds observable

For the physical solution at fixed viscosity \(\nu>0\), (2) and (7) imply, at every smooth time and every receiver radius,
\[
\boxed{
\frac{R^2|\Omega_R(u(t))|}{\nu}
\le\frac{G(u_0)}{d_\chi\nu}.}
\tag{9}
\]
The normalized circulation and normalized viscosity both acquire the same \(q^\alpha\) factor, so \(G(V_n)/\nu_n\le G(V_0)/\nu_0\). There is no contradiction with a finite rotational receiving gain: the initial outer swirl may provide a large finite circulation budget, and the initial receiver need not saturate (9).

Equation (9) controls a particular linear rotational receiver, not arbitrary meridional velocity or the full RMS norm. Even a swirl RMS estimate cannot simply be deduced by squaring \(|u^\theta|\le G/\rho\) and integrating across the axis: the resulting upper integral is singular. Thus this argument does not exclude every conceivable RMS cascade. Its decisive implication is that a return class with uniformly nonzero normalized rotation on a fixed normalized scale cannot support the proposed unbounded rotational growth.

## 5. Smooth forcing with a finite weighted torque budget

For axisymmetric forcing \(f\), the exact scalar equation acquires the source \(\rho f^\theta\). Comparison gives
\[
\boxed{G(u(t))\le G(u_0)+
 \int_0^t\|\rho f^\theta(s)\|_\infty\,ds.}
\tag{10}
\]
This extension assumes the displayed weighted source is integrable in time. Smoothness alone on all of space, or smoothness only before a terminal time, does not automatically supply that global finite bound. Compactly supported or appropriately decaying smooth forcing on a closed finite time interval does.

Suppose a complete physical cascade occurs before \(T<\infty\), with
\[
B_T=G(u_0)+\int_0^T\|\rho f^\theta(s)\|_\infty\,ds<\infty.
\]
Its cumulative normalized states satisfy
\[
G(V_n)\le Q_n^\alpha B_T.
\]
Equations (5)–(6) therefore remain valid with \(B_T\) replacing \(G(V_0)\). A finite smooth torque budget cannot sustain a uniformly nontrivial normalized rotating core under indefinite supercritical shrinkage.

The scaling is consistent with the actual transformed force. At cumulative radius \(Q_n\), its normalized force is
\(F_n=Q_n^{3+2\alpha}f\) at the corresponding physical points and times. Multiplying by normalized radius and integrating over normalized time makes its torque budget \(Q_n^\alpha\) times the physical torque budget. Assigning a fresh forcing independently at each normalized stage would be a different problem and need not obey a fixed finite physical budget.

## 6. What the obstruction changes

The requested pass5 return class combines \(\alpha>0\), shrink ratios bounded above by a number below one, uniformly positive normalized rotation, and uniform local/profile bounds. Those requirements imply (8), and so they cannot all persist through infinitely many exact unforced returns.

An alternative construction would have to alter at least one relevant requirement: permit normalized rotation/profile mass to vanish or its core width to collapse, cease indefinite supercritical shrinking, or leave the common-axis setting. This statement does not establish that any alternative succeeds. In particular, \(G(V_n)\to0\) by normalization does not itself prove regularity: its ratio to the simultaneously shrinking \(\nu_n\) need not tend to zero.

The initial pressure certificate and finite-feedback corollary are unchanged. This is a rigorous rejection of their proposed indefinitely recurring normalized rotating-core class, not a rejection of the finite episode or a resolution of the general Navier–Stokes problem.

The companion verify_circulation_return_obstruction.py checks the scalar PDE algebra, normalization powers, receiver-radius factors, a finite-stage example, the local Taylor factor, and the forcing-budget scaling. The comparison and no-return proofs are supplied above.

# The finite receiver gain does not return to a rigid rotating profile

This is a local analytical obstruction to treating the new finite gain as an already repeatable copy of the initial rigid core. It does not rule out a larger terminal class, a genuine return phase, or other Navier–Stokes mechanisms.

## A first-order profile change accompanies a second-order gain

Use the compact datum, notation and explicit interval from quantitative-core-stage.md. In the initial rigid neighborhood,

    u_0(x)=Omega Jx,
    u_t(0,x)=grad H(x),
    Delta H=0,
    H_zz(0)=K>0.

Let S(t)=sym grad u(t,0). Then

    S(0)=0,
    S(t)=t Hess H(0)+O(t^2).

The quantified theorem gives the stronger one-sided bound

    ||S(t)||op >= |S_zz(t)| >= k0 t/2,    0<=t<=tau.

Every rigid rotation, about any axis and with any translation, has skew velocity gradient and hence zero symmetric gradient. Thus its C1 distance from the evolved core is at least k0 t/2 when measured by the gradient operator norm at the origin. An orthogonal coordinate change preserves the norm of this symmetric part and cannot remove it. The estimate only asserts a local gradient mismatch; it is not a lower bound in a weaker velocity norm that does not control derivatives.

Meanwhile the rotational receiver amplitude is

    a_r(t)=1+(K/2)t^2+O(t^3).

The pressure-generated strain is therefore first order in time, while the first increase of the rotational coefficient is second order. Making the interval shorter does not make this profile discrepancy small relative to the gain. In dimensionless terms its size relative to the starting rotation is O(t ||Hess H||/Omega), whereas the relative rotational gain starts at K t^2/2.

## The radius selected from full RMS velocity has the same issue

Fix an allowed receiver radius r and exponent 0<a<1/2. Define the exact RMS selection function from rotational-receiver.md:

    F(q,t)=q^(a+1) U_(qr)(t)/U_r(0).

At t=0, F(q,0)=q^(a+2). Hence at (q,t)=(1,0),

    F=1, partial_q F=a+2>0, partial_t F=0.

The positive initial receiver energy and smooth solution permit differentiating the RMS near this point. The implicit function theorem gives a unique local branch q(t) near 1 with F(q(t),t)=1. This is an additional local choice rule; the more general finite-time intermediate-value argument did not require uniqueness.

Write

    D_r=integral chi_r|Jx|^2 dx,
    Q_r=integral chi_r|grad H|^2 dx,
    beta_r=(1/2)[K+Q_r/(Omega^2 D_r)]>0.

The exact energy derivatives in receiver-flux-budget.md imply

    U_r(t)/U_r(0)=1+beta_r t^2+O(t^3),
    q(t)=1-[beta_r/(a+2)]t^2+O(t^3).

In particular the isotropic change of scale and its matched amplitude change begin only at second order. Normalize the selected velocity back to the original coordinates/amplitude by

    v(t,y)=q(t)^(1+a) u(t,q(t)y).

This is a comparison of states, not a time-dependent change of variables claimed to solve the same equation. At the origin,

    sym grad_y v(t,0)=q(t)^(a+2) S(t)
                     =t Hess H(0)+O(t^2).

Thus the selected scaling does not cancel the first-order strain mismatch. Along the local branch, for sufficiently small positive t it remains at least k0 t/4 in operator norm. The full RMS gain and exact energy/Re scale relations are valid, but they do not return the normalized velocity to the initial rigid-rotation class at accuracy of the second-order gain.

## Consequence for the next mathematical step

A repeatable argument must either evolve a class that includes the growing strain and prove its invariant bounds, or add a true dynamical return which reduces that strain while retaining the useful receiving gain. It must also replenish the outer packet geometry and track the supplying material reservoir. Assigning the terminal state to a fresh rigid core would discard a leading term of the actual velocity. The local finite transfer should therefore be used as a measured gain component inside a future return construction, not as its completed return map.

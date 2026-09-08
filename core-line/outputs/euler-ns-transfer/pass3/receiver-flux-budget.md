# Exact flux and energy budgets for a fixed rotating receiver

This calculation identifies how the finite receiving-mode gain in rotational-receiver.md is supplied. It does not prescribe the pressure or an external strain, and it supplies no infinite return map.

Let Jx=(-y,x,0), let h(s) be a nonnegative smooth radial weight supported strictly inside the initially rigid core, and set phi(x)=h(|x|^2)Jx. Define

    I(t)=integral u(t,x) dot phi(x) dx,
    D=integral h(|x|^2)|Jx|^2 dx,
    Omega_h(t)=I(t)/D,
    E_h(t)=(1/2)integral h(|x|^2)|u(t,x)|^2 dx.

For the smooth unforced NS solution, compact support of the test field permits every integration by parts below. The solution itself need not remain compactly supported.

## 1. Angular momentum enters through advective and viscous flux

The field phi is divergence free. Consequently pressure makes no contribution to I', and the exact identity is

    I'(t) = 2 integral h'(|x|^2)(u dot x)(u dot Jx) dx
            + nu integral [4|x|^2 h''(|x|^2)+10h'(|x|^2)](u dot Jx) dx.

The two elementary identities used here are

    u_i u_j partial_j phi_i = 2h' (u dot x)(u dot Jx),
    Delta phi = (4|x|^2 h''+10h')Jx.

The hJ part of grad phi is skew and hence has zero contraction with u tensor u. A radius-scaled weight is included by taking h(s)=chi(s/r^2); its chain-rule factors must be retained.

At time zero the velocity in the support of h is Omega Jx. Its radial velocity is zero, and Delta u_0=0. Therefore I'(0)=0, including the viscous term. This is stronger than ignoring viscosity in a model: the term vanishes exactly at this endpoint.

Put H(x)=Omega^2(x^2+y^2)/2-p_0(x). The actual initial acceleration in the core is u_t(0,x)=grad H(x). The function H is harmonic there because Delta p_0=2Omega^2. Since Delta u_t(0)=0, differentiation gives

    I''(0) = 2Omega integral h'(|x|^2)|Jx|^2 x dot grad H(x) dx.

Radial harmonic moments yield I''(0)=Omega D H_zz(0)=Omega D K. This agrees with the independent receiver-mode calculation. Thus pressure first sets up an irrotational radial/strain velocity, and its coupling to the existing rotation then supplies angular-momentum flux into the fixed receiver. The pressure Hessian and the receiving-mode gain are consistent with the fact that pressure annihilates a divergence-free receiver test.

This receiver is a fixed spatial region. It can acquire material from a surrounding annulus and is not the material disk considered in material-core-transfer.md. A positive receiving-mode gain therefore does not contradict the cancellation of vorticity gain against material-area contraction.

## 2. Local kinetic energy has an additional potential-flow contribution

Initially I(0)=Omega D and E_h(0)=Omega^2 D/2. Since h u_0=Omega phi is divergence free,

    E_h'(0)=0,
    E_h''(0)=integral h|grad H|^2 dx + Omega^2 D K.

In particular, for K>0 the actual weighted kinetic energy grows at its first nonzero time order. Its increase is at least the one obtained from rotational projection, but it can be larger: the first-order harmonic-gradient velocity also contributes energy at second order. It must not be omitted when selecting a scale from the full velocity field.

For completeness the exact localized energy identity for arbitrary smooth times is

    E_h'(t) = integral (|u|^2/2+p) u dot grad h dx
              + (nu/2) integral |u|^2 Delta h dx
              - nu integral h|grad u|^2 dx.

This includes pressure work through the weighted boundary and viscous transport as well as dissipation. Local energy growth is compatible with nonincrease of total energy. It is supplied by flux from the rest of the same solution.

## 3. Consequence for the continuing stage target

The fixed receiving region can have both angular-momentum and kinetic-energy gain during the controlled short interval. The scale-matching construction in rotational-receiver.md uses the full RMS velocity, so its selected energy relation includes the potential-flow contribution above.

This proves only a finite transfer into a chosen receiver. To repeat it, the state outside and inside that receiver must generate a suitable successor without resetting the field. In particular, the annular and outer-packet reservoir that paid the flux cannot be treated as unchanged or as a free new copy. The material-circulation calculation and this fixed-region flux calculation are complementary tests for any proposed restart.

No new literature priority is claimed for these integration-by-parts identities. verify_receiver_flux.py checks the exact contractions and harmonic-moment coefficient used in this note; the local smooth solution is supplied by the analytical stage argument.

# The outer pressure reservoir and an explicit initial replenishment mechanism

This note analyzes the actual initial derivatives of compact smooth axisymmetric NS data on R3. It first rejects a fixed reservoir-ratio boundary for the passive outer-swirl design, then constructs an added meridional velocity with a favorable **swirl-channel** derivative. It does not establish positivity of the derivative of the total outer pressure, persistence, or a return map.

## 1. Definitions and a passive-reservoir loss

Use cylindrical coordinates (r,theta,z), R=sqrt(r^2+z^2), N(x)=1/(4pi R), and the pressure kernel

    Q_ij(x)=-partial_zzij N(x).

The rotating and straining core has initial velocity Lx near the origin, with

    L=diag(-b,-b,2b)+Omega J,  b>0, Omega>0.

Choose the circular, axisymmetric version of the disjoint outer swirl packets, so their velocity is A v=A v_theta(r,z)e_theta and their support lies strictly inside r<=|z|/4, away from the core. All fields are smooth on the axis when expressed in Cartesian coordinates. The radial compact extension of the core is specified in affine-core-pressure.md.

Let zeta be a smooth fixed cutoff equal to one on a neighborhood of the initial outer support and zero near the core. Define a contribution of the actual velocity's stress to central axial strain by

    C_F(t)=integral zeta(x) Q_ij(x) u_i(t,x)u_j(t,x) dx.

At insertion C_F(0)=A^2 C_v>0. This is a fixed-region quadratic stress contribution; no support separation is asserted for positive time.

The azimuthal direction is an eigenvector of Q, with no azimuthal-poloidal cross component. At time zero, the nondiffusive acceleration on the pure-swirl outer support is poloidal: self-advection is radial and the full axisymmetric pressure gradient has no azimuthal component. The core velocity is zero on that support. Therefore

    C_F'(0)=2nu A^2 integral Q_ij v_i Delta v_j dx
           =-2nu A^2 integral Q_ij partial_k v_i partial_k v_j dx.

The integration by parts uses Delta Q_ij=0 away from the origin and symmetry of Q. Derivatives of zeta vanish wherever the initial velocity or its derivatives are nonzero. Every partial_k v is horizontal. The horizontal block of Q is positive definite on the support cone, with lower bound 6/(4pi R^5). Hence C_F'(0)<0 at nu>0; at the Euler endpoint it is zero.

The actual central angular speed satisfies Omega'(0)=2bOmega. Thus the dimensionless reservoir ratio R_F=C_F/Omega^2 obeys

    R_F'(0)=C_F'(0)/Omega^2-4bR_F(0)<0

for every nu>=0 when b>0. A state rescaling u(x)->q^(1+a)u(qx) multiplies C_F and Omega^2 by the same factor q^(4+2a), when the cutoff/outer geometry is rescaled with the state. It does not restore this ratio.

The compact-core pressure calculation gives

    b'=C_F/2-(5/7)b^2-(1/5)Omega^2,
    (b/Omega)'=Omega[R_F/2-(19/7)(b/Omega)^2-1/5]

at insertion. For the natural margin

    F=R_F-(38/7)(b/Omega)^2-2/5,

the boundary F=0 has (b/Omega)'=0 and F'=R_F'<0. Thus that fixed boundary is crossed in the unfavorable direction by the unmodified compact-core/pure-swirl datum. This is a tangent failure of this specific proposed cone; interior slack and more general evolving classes are not ruled out.

## 2. An outer meridional velocity can reverse the swirl-channel derivative

Add a compact smooth axisymmetric meridional velocity cP_1, supported away from the core and equal to

    c(r e_r-2z e_z)

on a neighborhood of the entire swirl support. The full initial datum is now

    u_0=u_core+cP_1+A v.

The core and the two outer fields have disjoint supports from each other only in the core-versus-outer sense. The meridional and swirl outer fields deliberately overlap. Axisymmetry is retained.

Let w=u_theta. Its exact NS equation is

    w_t+u_r w_r+u_z w_z+(u_r/r)w=nu(Delta-1/r^2)w.

The swirl contribution to the same far pressure functional is

    C_theta(t)=integral zeta Q_theta w(t)^2 dx,
    Q_theta=3(4z^2-r^2)/(4pi(r^2+z^2)^(7/2)).

At time zero integration by parts, using div u=0, gives the exact derivative

Here the unweighted derivative integrals use only the initial outer swirl `w_out,0=A v_theta`, extended by zero away from the packets. They do not use the full initial azimuthal velocity, which also contains the core rotation. The cutoff defining C_theta is zero near that core and constant near every initial support, so no cutoff derivative contributes at this endpoint.

    C_theta'(0)
      =c integral Q_theta w_out,0^2 A(t) dx -2nu D_theta,
    D_theta=integral Q_ij partial_k(w_out,0 e_theta)_i
                           partial_k(w_out,0 e_theta)_j dx>0,
    t=r^2/z^2,
    A(t)=8-6t/(4-t)-21t/(1+t).

The variable t in A(t) is a dimensionless spatial ratio, not evolution time. One obtains it from

    [(r partial_r-2z partial_z)Q_theta-2Q_theta]/Q_theta.

The diffusion term retains the full Cartesian derivatives, including the variation of e_theta. Replacing it by the scalar gradient of w alone would omit a term. On 0<=t<=1/16,

    A(t)>=gamma:=2381/357>6.

Consequently

    (C_theta/Omega^2)'(0)
      >=(C_theta(0)/Omega^2)
         [gamma c-4b-2nu D_theta/C_theta(0)].

The quotient D_theta/C_theta(0) is independent of the swirl amplitude A. For every fixed nonzero smooth swirl profile and fixed finite nu,b, it is finite. Choosing

    c>[4b+2nu D_theta/C_theta(0)]/gamma

therefore reverses the previously unfavorable **swirl-pressure ratio** derivative. This is an actual initial NS calculation, including meridional advection and swirl diffusion, rather than a prescribed pressure control.

## 3. A compact meridional field with the required local velocity

Choose a smooth alpha(s), equal to one near s=0 and zero outside a compact interval, and a smooth even beta(z), supported in two intervals around +/-d away from zero, and equal to one on neighborhoods containing the axial part of the swirl support. For a small sigma>0, define the axisymmetric streamfunction

    Psi(r,z)=-c z r^2 alpha(r^2/sigma^2) beta(z).

Using u_r=-(1/r)Psi_z and u_z=(1/r)Psi_r gives

    (cP_1)_r=c r alpha(r^2/sigma^2)[beta+z beta'],
    (cP_1)_z=-2c z[alpha(r^2/sigma^2)
                  +(r^2/sigma^2)alpha'(r^2/sigma^2)]beta.

This is a smooth compact divergence-free vector field, including at r=0. Where both cutoffs are identically one it equals c(r e_r-2z e_z). Choose the swirl support strictly inside those plateaus. With sigma small enough the full outer support also lies in the required cone and remains disjoint from the compact core.

The pressure coefficient of this unit meridional field has a useful sign for sufficiently thin support:

    C_P=-partial_zz p[P_1](0),
    C_P/sigma^2 -> -48 I_alpha I_beta <0,
    I_alpha=integral_0^infinity s[alpha(s^2)+s^2alpha'(s^2)]^2 ds,
    I_beta=integral_R |z|^(-3) beta(z)^2 dz.

The limit follows by r=sigma s. The vertical velocity is order one, the radial velocity is order sigma, and Q_zz tends to -24/(4pi|z|^5); the axisymmetric volume measure is 2pi sigma^2 s ds dz. Dominated convergence applies because beta is supported away from z=0 and all profiles are compact. Both integrals are positive for the stated nontrivial profiles. With fixed profiles the error in C_P is O(sigma^4).

If sigma is changed, the swirl profile must continue to fit inside the plateau. Its diffusion quotient can grow as its transverse scale shrinks. No uniform gain or finite-c bound as sigma->0 is asserted.

## 4. A nonempty set of favorable initial conditions, with a crucial limit

Axisymmetric meridional and swirl velocities have identically zero cross contribution to the pressure Poisson source. This applies to the overlapping P_1 and v. Core/outer stresses are initially disjoint. Hence the total outer axial pressure coefficient at insertion is exactly

    C_total=A^2 C_v+c^2 C_P.

For a sufficiently thin fixed meridional support, C_P<0. Choose b,Omega>0, then choose c to make the swirl-ratio derivative positive as above. Finally set

    A^2=[(38/7)b^2+(2/5)Omega^2-c^2 C_P]/C_v>0.

The actual initial core then satisfies

    Omega'=2bOmega>0,
    b'=2b^2,
    (b/Omega)'=0,
    (C_theta/Omega^2)'>0.

The strained-receiver identities imply positive initial rotational and RMS velocity growth. These conditions can be achieved by one smooth compact datum at a fixed positive viscosity. They remove the specific passive-swirl tangent failure; they do not close a return map.

**The missing term must remain explicit:** C_total(t) includes the evolving meridional stress and any newly populated region between supports. Positivity of C_theta'(0) does not imply positivity of C_total'(0). The meridional field responds to the same swirl and global pressure that drive the core. Its derivative may be of the same order as the favorable swirl term. No rigorous sign for that full response follows from this channel calculation. The full cubic pressure derivative is treated separately in full-pressure-feedback-gate.md, and a favorable but uncertified numerical example is recorded in numerical-pressure-evidence.md. Initial freedom to choose c is not a license to keep prescribing it during evolution.

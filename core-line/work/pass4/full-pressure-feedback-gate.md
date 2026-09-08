# The full pressure derivative needed to preserve the strained-core ratio

The meridional pump in outer-pressure-reservoir.md improves one actual outer-swirl contribution. The following exact test keeps the missing pressure response explicit. It is the next scalar sign problem for this candidate, not a solved invariant-cone result.

## 1. Actual central equations and a neutral initial boundary

For an odd axisymmetric smooth NS solution, let

    b(t)=partial_z u_z(t,0)/2,
    Omega(t)=omega_z(t,0)/2,
    beta(t)=b(t)/Omega(t),
    B(t)=partial_z Delta u_z(t,0)/2,
    W(t)=Delta omega_z(t,0)/2.

While Omega is nonzero, the exact central equations are

    b'=-2b^2-p_zz/2+nu B,
    Omega'=2bOmega+nu W,
    beta'=-(p_zz+8b^2)/(2Omega)
             +nu[B/Omega-bW/Omega^2].

Take the compact initial datum with an affine core Lx and tune its total initial pressure to p_zz(0)=-8b^2. This is precisely b'(0)=2b^2 and beta'(0)=0. The full central matrix then has L'(0)=2bL.

On the initial affine neighborhood Delta u_0=0 and Delta p_0=-tr(L^2) is constant. The exact initial acceleration satisfies

    u_t(0)=-L^2 x-grad p_0,
    Delta u_t(0)=0.

Consequently B(0)=W(0)=B'(0)=W'(0)=0. These are local endpoint identities, not neglect of viscosity during the later evolution.

Differentiating the exact ratio equation therefore gives

    beta''(0)=-[p_zz'(0)+32b^3]/(2Omega).

Thus a strictly favorable bend of beta at this neutral boundary requires

    p_zz'(0)<-32b^3.

Positive swirl-pressure growth alone is not this inequality.

## 2. An exact evaluation target with every pressure contribution retained

Let G_0=grad u_0 and compute the initial acceleration from the full datum on R3:

    u_1=nu Delta u_0-P div(u_0 tensor u_0),
    g_1=2 tr(G_0 grad u_1),
    p_1=N*g_1,   N(x)=1/(4pi|x|).

Here P is the whole-space Leray projection. Then p_zz'(0)=partial_zz p_1(0) exactly. With the distributional Hessian of N this is

    p_zz'(0)=-g_1(0)/3
       +PV integral [(3z^2-|x|^2)/(4pi|x|^5)] g_1(x) dx.

At the neutral insertion, L'=2bL and tr(L^2)=6b^2-2Omega^2, so

    g_1(0)=24b^3-8bOmega^2,
    -g_1(0)/3=-8b^3+(8/3)bOmega^2.

The contact term must be included. The remaining source is supported where grad u_0 is nonzero, although u_1 itself has its full nonlocal pressure tail. Near the origin g_1 is smooth; the principal-value prescription is part of the formula.

For fixed profiles this target is a cubic functional of the velocity amplitudes plus viscosity times a quadratic functional. It can be evaluated without a full time simulation, but requires the actual nonlocal u_1. Differentiating the initial compact-core pressure formula while pretending the evolved field stays in the same radial-profile family would calculate a different quantity.

## 3. Why the previous outer-reservoir test is insufficient by itself

For the initially disjoint core and outer field, define the fixed-region outer stress contribution C_F as in outer-pressure-reservoir.md. At insertion

    p_zz=-(18/7)b^2+(2/5)Omega^2-C_F.

To track the possible failure of that profile identity later, define the exact defect

    D(t)=p_zz(t,0)+C_F(t)+(18/7)b(t)^2-(2/5)Omega(t)^2.

D(0)=0, but D'(0) is not fixed by the initial identity. Writing

    F=C_F/Omega^2-(38/7)beta^2-2/5,

the actual ratio equation is

    beta'=Omega F/2-D/(2Omega)
             +nu[B/Omega-bW/Omega^2].

At the neutral insertion,

    beta''(0)=Omega F'(0)/2-D'(0)/(2Omega).

For passive pure-swirl outer packets, the proved F'(0)<0 is an actual exit from that explicitly defined observable cone. It does not by itself determine beta'' without the evolving core/profile defect. For the new pump, growth of C_theta/Omega^2 is likewise only one part of the full test. The formula above explains exactly which additional response is missing.

## 4. For fixed profiles, the gate has a finite coefficient expansion

Write the proposed datum as

    u_0=bS+Omega W+cP+Av,

where S,W are the unit core strain/rotation profiles, and P,v are the unit outer meridional/swirl profiles. Define the exact bilinear operators

    B(a,b)=(1/2)P_Leray[(a dot grad)b+(b dot grad)a],
    Pi(a,b)=partial_zz p[a,b](0),
    p[a,b]=N*tr(grad a grad b).

Pi and B are symmetric bilinear. The initial pressure derivative is exactly

    p_zz'=2Pi(u_0,nu Delta u_0-B(u_0,u_0)).

Disjoint core/outer supports make B(S,P), B(S,v), B(W,P), and B(W,v) zero. Axisymmetric poloidal/azimuthal fields have zero bilinear pressure source. These facts reduce the derivative to

    p_zz'=t300 b^3+t210 b^2c+t120 bc^2+t030 c^3
          +Omega^2(tWb b+tWc c)+A^2(tvb b+tvc c)
          +nu(dS b^2+dW Omega^2+dP c^2+dv A^2).

All coefficients depend only on the four chosen profiles. Their exact definitions are:

| Coefficient | Full nonlocal expression |
|---|---|
| t300 | -2Pi(S,B(S,S)) |
| t210 | -2Pi(P,B(S,S)) |
| t120 | -2Pi(S,B(P,P)) |
| t030 | -2Pi(P,B(P,P)) |
| tWb | -2Pi(S,B(W,W))-4Pi(W,B(S,W)) |
| tWc | -2Pi(P,B(W,W)) |
| tvb | -2Pi(S,B(v,v)) |
| tvc | -2Pi(P,B(v,v))-4Pi(v,B(P,v)) |
| dj, j=S,W,P,v | 2Pi(j,Delta j) |

In particular, the coefficient tvc includes both the swirl advection and the meridional response to swirl. Keeping only its second term would omit exactly the backreaction that the new pump must pass. These expressions use the whole-space Leray projection, not a truncation to the four profiles; the generated velocity fields B(a,b) are retained in full.

For completeness, B(S,W) and B(P,v) are azimuthal and already divergence free, so Leray leaves them unchanged and they retain their respective compact supports. This proves the further cancellations Pi(v,B(S,W))=Pi(W,B(P,v))=0. The generated poloidal fields generally have nonlocal tails; the table keeps their cross-region effects.

The [independent coefficient audit](coefficient-expansion-audit.md) sharpens the twelve-slot formula for these particular profiles. Exactly dS=dW=0, because Laplacian variations stay within the radial profile classes and do not change their central linear amplitudes. Also dv>0 for the selected horizontal cone-supported swirl. Thus ten coefficients can be nonzero. There is no corresponding established sign for dP.

Substitute the neutral-boundary amplitude

    A^2=[(38/7)b^2+(2/5)Omega^2-c^2 C_P]/C_v.

After choosing theta=Omega/b, mu=c/b and epsilon=nu/b in fixed profile units, write H=[38/7+(2/5)theta^2-C_P mu^2]/C_v. The full gate divided by b^3 is

    G=32+t300+t210 mu+t120 mu^2+t030 mu^3
      +theta^2(tWb+tWc mu)+H(tvb+tvc mu)
      +epsilon[dP mu^2+dv H].

The required sign is G<0. This is an exact finite reduction of the initial-time test, retaining the full generated velocity. Apart from the exact core zeros and the positive sign of dv, rigorous values or bounds for these coefficients have not yet been established here. Exploratory quadrature is reported separately and does not certify this inequality.

## Next calculation

For the compact meridional-pump datum, rigorously bound

    T(u_0,nu)=partial_zz[N*(2tr(grad u_0 grad u_1))](0)+32b^3.

The [exploratory evaluations](numerical-pressure-evidence.md) now give a concrete candidate with a favorable sign after the positive-viscosity term is included. They are not interval-certified. Bound the source projection, radial integration, reconstructed pressure, amplitude tuning and arithmetic errors to prove a nonempty parameter range with T<0, or reject the candidate if validated bounds fail. Even a certified favorable sign would be a further initial-time compatibility test. It would still need a quantitative persistence estimate, receiver gain, full profile control and an actual return on inherited data before it could be iterated.

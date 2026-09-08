# A compact active-vortex datum with positive initial strain feedback

This is a local-in-time lemma for the ordinary, unforced three-dimensional Navier–Stokes equation on R3, at every fixed viscosity nu>0. It is not a blowup theorem, a quantitative large-amplification theorem, or a regeneration lemma. The calculation strengthens the earlier irrotational-core diagnostic by putting nonzero vorticity in the core. No literature novelty is claimed.

## Statement and explicit data

Choose Omega>0 and a smooth radial cutoff psi(s) that is one near s=0 and zero for s>=R^2. Put

    c(x,y,z) = Omega psi(x^2+y^2+z^2) (-y,x,0).

This field is smooth, compactly supported, divergence free, and a rigid rotation near the origin. Its central vorticity is 2 Omega e_z.

Choose d>0 and 0<epsilon<d/5 with R<d-epsilon. For a nonzero smooth bump eta supported in [0,1), set

    chi(x,y,z) = eta((x^2+2y^2+(z-d)^2)/epsilon^2)
               +eta((x^2+2y^2+(z+d)^2)/epsilon^2),
    v = curl(chi e_z) = (chi_y,-chi_x,0).

The bump can be chosen smooth as a function of its nonnegative argument, with all derivatives zero at 1. Then v is a nonzero, smooth, compact, horizontal, divergence-free field. Its support is disjoint from c and lies in the double cone

    x^2+y^2 <= z^2/16.

The unequal x/y coefficients make this example nonaxisymmetric. Both c and v are odd under full spatial inversion. Let p[v] be the decaying pressure defined by

    -Delta p[v] = partial_i partial_j(v_i v_j),
    C_v = -partial_zz p[v](0).

We prove below that

    C_v >= (6/(4 pi)) integral |v(y)|^2 / |y|^5 dy > 0.

For initial data u_0=c+A v, choose

    K := A^2 C_v - (2/5) Omega^2 > 0.

The unique local smooth NS solution has zero velocity at the origin (by inversion symmetry), and

    partial_z u_z(0,0)=0,
    partial_t partial_z u_z(0,0)=K,
    omega_z(0,0)=2 Omega,
    partial_t omega(0,0)=0,
    partial_tt omega_z(0,0)=2 Omega K.

Consequently, as t decreases to zero through positive times,

    partial_z u_z(t,0) = K t + o(t),
    omega_z(t,0) = 2 Omega + Omega K t^2 + o(t^2).

The axial strain and central vorticity therefore both increase for a sufficiently short positive interval. The data have finite kinetic energy and force is exactly zero. The derivatives displayed here are independent of nu because the initial field is affine in an open neighborhood of the core; the interval of validity need not be uniform in nu or the parameters.

## 1. Exact pressure of the compact rotating core

Use N(x)=1/(4 pi |x|), so -Delta N=delta. For c above, the pressure source is

    g=-Delta p[c] = -Omega^2 (2 psi^2 +4 (x^2+y^2) psi psi'),

where prime means derivative with respect to s=|x|^2. The distributional identity

    partial_zz N = PV [(3z^2-|x|^2)/(4 pi |x|^5)] - delta/3

must retain its delta term. The isotropic term 2 psi^2 has zero principal-value angular contribution. With mu=z/|x|,

    average_sphere [(3 mu^2-1)(1-mu^2)] = -4/15,
    integral_0^infinity r psi(r^2) psi'(r^2) dr = -1/4.

It follows that

    partial_zz p[c](0) = (2/3) Omega^2 -(4/15) Omega^2
                       = (2/5) Omega^2.

Thus the isolated compact rotating core initially compresses axially. The outer packets have to overcome this term; it cannot be omitted.

## 2. A sign-definite pressure kernel for horizontal outer packets

Because the origin is outside the support of v, no principal value or distributional contact term is needed for its contribution:

    partial_zz p[v](0) = integral partial_zzij N(y) v_i(y) v_j(y) dy.

Restrict this fourth-derivative tensor to horizontal vectors. Writing r=|y|, y_perp=(x,y), it is the matrix

    (1/(4 pi r^9)) [
        (105 z^2-15 r^2) y_perp y_perp^T
        + (3 r^4-15 r^2 z^2) I_2 ].

For t=|y_perp|^2/z^2 in [0,1/16], its radial and tangential eigenvalues, after multiplication by 4 pi r^5, are respectively

    (-12+81t-12t^2)/(1+t)^2,
    (-12-9t+3t^2)/(1+t)^2.

Both are at most -6 throughout that interval. For example, adding 6 to the first numerator gives -6+93t-6t^2, which is increasing on the interval and still negative at 1/16; the second is even smaller. This proves the stated lower bound for C_v. This estimate needs no oscillatory approximation, phase average, or large carrier frequency.

## 3. Exact disjoint-support and vorticity calculation

Pointwise c_i v_j=v_i c_j=0, so the full pressure at time zero satisfies

    p[u_0] = p[c] + A^2 p[v]

up to an irrelevant constant. This equality is only asserted at time zero: the subsequent velocity fields need not remain disjoint.

Near the origin u_0=Omega(-y,x,0). Its Laplacian is zero, its vorticity is the constant 2 Omega e_z, and its vorticity time derivative is identically zero there by the NS vorticity equation. The initial nonlinear acceleration has zero z component. Therefore

    partial_t partial_z u_z(0,0) = -partial_zz p[u_0](0) = K.

Differentiate

    omega_t = -u dot grad omega + omega dot grad u + nu Delta omega.

At the origin all differentiated terms except omega_0 dot grad u_t vanish, since u_0(0)=0, grad omega_0=0, and omega_t is identically zero in an initial core neighborhood. In particular Delta omega_t(0,0)=0. This yields

    partial_tt omega_z(0,0) = 2 Omega partial_z u_t,z(0,0) = 2 Omega K.

Standard local smooth NS well-posedness for this compact smooth datum justifies these right time derivatives and the Taylor conclusions. The symmetry u_0(-x)=-u_0(x) is preserved by uniqueness, so the origin is also a material trajectory.

## What this does and does not supply

This provides an active central vortex and an actual pressure-mediated onset of increased stretching, with explicit sign and amplitude threshold. Unlike a homogeneous affine strain, the inducing velocity has finite energy. It does not quantify a large retained gain, keep the outer packets favorable for a full strain time, restore a next-stage envelope, or show increasing Reynolds number. Those require control of the evolving, coupled field beyond its initial time derivatives. The absence of those estimates is the remaining mathematical gap, not a force-admissibility problem in this local lemma.

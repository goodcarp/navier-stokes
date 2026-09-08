# Axial images explain the initial core-pressure offset

For the full compact datum at A=1020, lambda=1/4, c=7/5, the axial image
correction accounts for nearly all of the difference between the periodic
pilot and the whole-space initial pressure. Using 32 positive/negative image
pairs, the independent retained-image quadratures give:

| Axial period | Image correction | Whole-space reference plus images |
|---:|---:|---:|
| 16 | +0.480133465314626 | -8.22890173732355 |
| 32 | +0.008383184525667 | -8.70065201811251 |
| 64 | +0.000209482248515 | -8.70882572038966 |

The whole-space reference is p_zz=-8.709035202638178. The fd768 pilot's
reported values, -8.228850200576058 at period 16 and -8.700602790969402 at
period 32, exceed these retained-image references by approximately
5.15367e-5 and 4.92271e-5 respectively. Those similar residuals are consistent
with a remaining discretization/boundary-representation error. This
comparison does not identify its exact cause or provide a certified solver
error bound.

The retained-image calculations are numerical quadratures. The omitted-image
bound below is analytic and separately verified, but does not enclose the
quadrature error of the retained images.

## 1. Exact initial periodic-image formula

The initial datum is supported in the ball of radius six. For L>12 its
periodic copies have disjoint initial supports, so the initial periodic
pressure source is exactly the periodization of the whole-space source.
The pressure Hessian image series is absolutely convergent:

    delta_L p_zz(0)
      = sum_(n nonzero) integral D_zz N(r,z+nL) g(x) dx,
    N(x)=1/(4pi|x|),  g=tr((grad u0)^2).

Only the angular-zero source is needed for this scalar at the axis. No
finite radial pressure boundary is used in this image computation. There
are no image contact terms: every image center lies outside the compact
source support.

For a compact divergence-free velocity,

    g=partial_i partial_j(u_i u_j),   integral g dx=0.

Consequently a numerically better source kernel for the nth positive/
negative pair is

    K_n(r,z)=D_zz N(r,z+nL)+D_zz N(r,z-nL)-1/[pi(nL)^3].

The subtracted constant is the exact two-image kernel value at the source
origin. This uses the exact zero-monopole identity; it is not an arbitrary
numerical adjustment of the pressure. The pair kernel has no linear term
at the origin and now visibly begins at order (nL)^-5. The raw order-72
quadrature of the full source monopole was 3.5e-11; leaving that residual in
an uncorrected source sum would introduce a spurious order n^-3 term.

An independent nonsingular stress formula is

    delta_L p_zz(0)
      = sum_(n nonzero) integral D_zzij N(r,z+nL) u_i u_j dx.

The script uses the full angularly averaged stress, including the seed's
mode-zero quadratic covariance. The radial M contribution is integrated
over the entire sphere R<=6. The outer/seed additions use cylindrical
quadrature and both axial packets. Source and stress sums agree to 2.6e-15
at period 16 and more closely at periods 32 and 64.

## 2. Image-count convergence

The following are positive/negative image pairs: N=8 means all images
n=-8,...,-1,1,...,8.

| L | N=8 correction | N=16 correction | N=32 correction |
|---:|---:|---:|---:|
| 16 | 0.480124410315586 | 0.480132868110961 | 0.480133465314626 |
| 32 | 0.008382904442204 | 0.008383165917005 | 0.008383184525667 |
| 64 | 0.000209473518451 | 0.000209481667416 | 0.000209482248515 |

For very distant image pairs the leading term is

    [6/(pi(nL)^5)] integral (2u_z^2-u_r^2-u_theta^2) dx.

The numerical stress moment here is 104061.661771480, making its
coefficient 198743.134287456/(nL)^5. The first image at L=16 is not in a
regime where retaining only this leading moment is accurate: higher spatial
moments of the large compact annular field matter. Thus simply dividing
the period-16 correction by 32 does not accurately predict period 32.

## 3. A rigorous operator-norm bound for every omitted image

Let H=4pi R^5 (D_zzij N), and t=z^2/R^2. In the cylindrical basis,

    H_rr=105t(1-t)-12,  H_tt=3-15t,
    H_zz=105t^2-90t+9,
    H_rz^2=225t(1-t)(7t-3)^2.

For its r,z block,

    det(24I-H)=180(1-t)(t+3),
    det(24I+H)=36[16-5(t-1)^2].

Both determinants are nonnegative on 0<=t<=1. Their first diagonal pivots
are positive:

    24-H_rr=105(t-1/2)^2+39/4,
    24+H_rr=12+105t(1-t).

Together with -12<=H_tt<=3, this proves

    ||D_zzij N(x)||op <= 24/(4pi |x|^5).

Writing E=||u0||_2^2/2 and using |x+nLe_z|>=nL-6 on the full support,
the exact tail after N pairs satisfies

    |tail_(L,N)| <= (24E/pi) sum_(n>N) (nL-6)^-5
                 <= 6E/[pi L(LN-6)^4].                         (1)

The integral test applies to the decreasing positive summand. Unlike a raw
source-L1 estimate, (1) retains the correct fifth-power image decay without
requiring numerical cancellation of source moments.

## 4. An exact coarse energy bound removes quadrature from the tail

The prior cutoff certificate gives 0<=chi<=1, -4<chi'<=0. For the compact
radial meridional extension, direct spherical averaging and integration by
parts give the exact identity

    E(S_Phi)=(8pi/15) integral_0^infinity s^(7/2) Phi'(s)^2 ds.

On the disjoint core, inner and outer transitions, use the maximum
derivative times the total variation. With c=7/5 this yields

    integral s^(7/2) Phi'^2 ds
      <= 4+(5/2)^7(16c^2/25)+6^7(3c^2/11).

For the core swirl, its coefficient has absolute value at most 11/3 in
the unit ball, so E(W_psi)<=484pi/135. For the outer swirl and seed,

    E(v)<= (9pi/20)(63/200)^4,
    K <= (pi/32)(392/5)(12/5).

The last inequality uses the existing Pass8 R0 and Z0 bounds and
lambda=1/4. Angular averaging, pointwise orthogonality of meridional and
azimuthal components, and disjoint core/outer swirl supports remove the
cross-energy terms. Using pi<22/7 therefore gives

    E < 2015718896182553/7560000000 < 270000.

Together with pi>3, (1) has the wholly analytic bound

    |tail_(L,N)| < 540000/[L(LN-6)^4].                          (2)

| L | Tail bound after N=8 | after N=16 | after N=32 |
|---:|---:|---:|---:|
| 16 | 1.524e-4 | 8.640e-6 | 5.149e-7 |
| 32 | 4.320e-6 | 2.575e-7 | 1.572e-8 |
| 64 | 1.288e-7 | 7.857e-9 | 4.853e-10 |

The table rounds upward. These bounds concern the exact omitted image
integrals only. The script also reports a tighter value obtained by
inserting the numerical energy into (1); that latter value is explicitly
not the certified tail bound.

## 5. Files and scope

initial_periodic_images.py contains both integral formulations and the
source-monopole subtraction. initial-periodic-images.json stores the
order-72 reference and all 8/16/32-pair comparisons; the separate order-48
file supplies a quadrature comparison. verify_periodic_image_tail.py
checks the matrix identities, radial energy identity and rational energy/
tail bounds; its result is periodic-image-tail-certificate.json.

Changing the per-piece quadrature order from 48 to 72 changed the 32-pair
stress correction by less than 4.8e-14 at L=16, 4.2e-16 at L=32, and
4.8e-18 at L=64. The independently differentiated source formulation changed
by at most 4.4e-12. These observed changes remain distinct from the exact
omitted-image bounds.

This isolates the initial periodic-versus-whole-space pressure difference.
It does not correct a later periodic trajectory into a whole-space NS
trajectory, certify the retained-image quadrature, or remove radial,
temporal, projection or unresolved-mode errors from the pilot.

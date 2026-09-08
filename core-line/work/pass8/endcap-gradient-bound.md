# A direct whole-space radial-gradient bound for the m=4 endcap residual

This is an exact pressure estimate. It uses all radial tails and all axial endcaps, and does not impose an artificial radial or axial boundary. The original proposed Hankel constants are correct; a sharper Bessel derivative bound removes its second and third terms.

## 1. Residual sign and mode convention

Use complex mode amplitudes, with the physical pressure equal to `Re[e^(im theta) p_m(r,z)]`. Let

`A_m=-partial_rr-r^(-1)partial_r+m^2/r^2`, `A_m P=S`,

where the actual mixed source is `Re[e^(im theta) S(r)F(z)]`. The separated horizontal trial is `P(r)F(z)`. If p_m is the true decaying whole-space pressure amplitude, then

`-Delta_m[P F]=S F-P F''`,

so the error `e_m=p_m-P F` has the **positive** source

`-Delta_m e_m=P F''`, `e=N*[P(r)e^(im theta)F''(z)]`.       (1)

Any factor two in the definition of the mixed pressure source must already be present in S and hence P. There is no additional factor two in the estimate below. The bound on the complex amplitude controls the radial derivative of its real part at every angle.

## 2. Hankel normalization and the factor one-half

Take the order-m Hankel convention

`P_hat(k)=integral_0^infinity P(rho) J_m(k rho) rho d rho`,

with inverse measure `k dk`. Hankel transformation of the three-dimensional modal Poisson equation gives the one-dimensional operator `-partial_zz+k^2`, whose whole-line Green function is `exp(-k|z-zeta|)/(2k)`. Therefore, for `d=|z0-zeta|>0`,

`partial_r e_m(r,z0)=integral F''(zeta) integral_0^infinity P(rho) K_(m,r)(r,rho,d) rho d rho d zeta`,

`K_(m,r)(r,rho,d)=(1/2) integral_0^infinity exp(-kd) k J_m'(kr)J_m(k rho) dk`.       (2)

Hankel Plancherel gives exactly

`integral_0^infinity rho |K_(m,r)|^2 d rho`

` =(1/4) integral_0^infinity exp(-2kd) k |J_m'(kr)|^2 dk`.       (3)

No pi or two-pi factor is missing: (2)--(3) are radial mode-amplitude formulas. Such factors appear only when converting amplitudes to complete angular L2 norms. The standard Hankel transform and its inversion are recorded in [DLMF §10.22(v)](https://dlmf.nist.gov/10.22#v); the chosen weighted convention is equivalent to its square-root-weight convention.

## 3. A sharper global derivative bound from a positive weight

For integer m>=1 and real x>=0, Poisson's integral gives

`J_m(x)=c_m x^m integral_(-1)^1 exp(ixt)(1-t^2)^(m-1/2) dt`,

`c_m=1/[2^m sqrt(pi) Gamma(m+1/2)]`.

This representation is [DLMF 10.9.4](https://dlmf.nist.gov/10.9.E4). Differentiate and integrate the term `t partial_t exp(ixt)` by parts. The endpoint term vanishes, including for m=1. With `w=(1-t^2)^(m-1/2)`, the resulting weight is

`(m-1)w-tw'=[(m-1)+m t^2](1-t^2)^(m-3/2)>=0`.

Its integral is `m integral w`. The beta integral therefore proves the global bound

`|J_m'(x)|<=x^(m-1)/[2^m(m-1)!]`.       (4)

At x=0 this holds by continuity. For m=4 it is simply `|J_4'(x)|<=x^3/96`. The weaker recurrence bound `x^3/96+x^5/7680` is also valid, but the second term is unnecessary.

Substituting (4) into (3), and using `integral_0^infinity exp(-2dk) k^7 dk=7!/(2d)^8`, gives

`||K_(4,r)(r,.,d)||_(L2(rho d rho))^2 <= 35 r^6/(65536 d^8)`,

or

`||K_(4,r)|| <= sqrt(35) r^3/(256 d^4)`.       (5)

For comparison, expanding the weaker bound reproduces the originally proposed expression exactly:

`(1/4)[a^2 r^6 7!/(2d)^8+2ab r^8 9!/(2d)^10+b^2 r^10 11!/(2d)^12]`,

with `a=1/96`, `b=1/7680`. Thus the sharpened result is a reduction of a valid bound, not a normalization correction.

## 4. A computable one-dimensional error budget

Define the full radial norm and the weighted axial endcap norm

`P2=integral_0^infinity rho |P(rho)|^2 d rho`,

`W4(z0)=integral_R |F''(zeta)|/|z0-zeta|^4 d zeta`.

Cauchy--Schwarz in rho and the triangle inequality in zeta give the direct whole-space bound

`|partial_r e_m(r,z0)| <= [sqrt(35)/256] r^3 sqrt(P2) W4(z0)`       (6)

for m=4. In squared form, suitable for rational interval certificates,

`|partial_r e_m|^2 <= (35/65536) r^6 P2 W4^2`.

At the maximizing circle, `r_*^2<=3969/64000`, so a uniform rational prefactor is

`|partial_r e_m(r_*,z0)|^2 <= [437664515463/3435973836800000000] P2 W4(z0)^2`.       (7)

For the actual product axial cutoff at `z0=4`, F is constant on `|z-4|<=9/40`; its F'' source is separated by at least `9/40=0.225`. Thus W4 is finite. It must include the reflected packet near z=-4 as well as both endcaps of the packet near z=4. The simpler bound `W4<=(40/9)^4 ||F''||_1` is valid but may lose useful distance information.

P2 must include the horizontal trial's radial tails. If S is supported in `[rho_-,rho_+]`, regularity at zero and decay at infinity give

`P(rho)=P(rho_-)(rho/rho_-)^4` below rho_-,

`P(rho)=P(rho_+)(rho_+/rho)^4` above rho_+.

The omitted tail integrals are therefore exactly

`|P(rho_-)|^2 rho_-^2/10` and `|P(rho_+)|^2 rho_+^2/6`.

These formulas permit a finite radial quadrature plus exact tails. An estimate using only the source annulus would not certify P2.

Equations (6)--(7) reduce this method to rigorous one-dimensional enclosures for P2 and W4. No numerical values for those two integrals are certified in this note. The separate harmonic-slab energy estimate is already sufficient for the current parent certificate, so no costly additional computation is required here.

## 5. Independent check of the alternative slab normalization

For a real mode `Re[e^(im theta)e_m]` harmonic in `|z-z0|<h`, write its Hankel amplitude as `A(k)cosh(k(z-z0))+B(k)sinh(k(z-z0))`. Exact integration over the symmetric slab gives

`||grad e_real||_2^2 >= pi integral_0^infinity k^2 sinh(2kh)|A(k)|^2 dk`.

Applying Cauchy--Schwarz to `e_m,r(r,z0)=integral A(k)k^2 J_m'(kr)dk` gives

`|e_m,r|^2 <= [||grad e_real||_2^2/pi] integral_0^infinity k^2 J_m'(kr)^2/sinh(2kh) dk`.

The factors in the independently developed slab certificate therefore agree with the same Hankel convention. This note does not substitute a fixed-circle acceleration bound for the acceleration of a moving global maximum, particularly when the initial axial maximum has a flat plateau.

`verify_endcap_gradient_bound.py` checks the positive-weight algebra, coefficient and factorial identities, mode Green-function jump, radial tails, and slab energy prefactor.

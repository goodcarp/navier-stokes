# The exact Kelvin model preserves the stress sign but depletes its magnitude

2026-09-08. For the existing outer seed, incorporating incompressible velocity polarization reverses the naive inference from passive phase winding: the radial carrier increases, but the useful negative Reynolds covariance decreases strictly. The model has an exact finite integrated-stress budget. This is not a persistence theorem or a no-go theorem for the actual compact NS datum.

## 1. An exact affine model with the correct initial local gradient

At an initial mean-angular-momentum maximizing circle, set

`h=m/r_*`, `k0=20`, `m=4`, `c=7/5`, `Omega_s=A B_*>0`.

The proved circle bounds give

`1000/63<h<1600/63`, `390<Omega_s<1020`.

The initial mean gradient is, in its cylindrical basis,

`A0=[[c,-Omega_s,0],[-Omega_s,c,0],[0,0,-2c]]`.

Let `R(t)` rotate about z at the constant rate `Omega_s`, and prescribe the Cartesian affine gradient `A(t)=R(t) A0 R(t)^T`. This is an exact divergence-free affine NS background: A is symmetric, so `A'+A^2` is symmetric and its acceleration is balanced by a quadratic pressure. It is unbounded in space and has infinite energy. Its prescribed later pressure and gradient are **not** the later pressure and gradient of the compact PASS7 datum.

For this model, one horizontal plane wave is an exact Kelvin perturbation. Its self-advection vanishes because its velocity polarization is perpendicular to its wave vector. In the rotating frame put `J=[[0,-1],[1,0]]`, let `kappa` be the horizontal wave vector and `a` its real velocity amplitude. The exact equations are

`kappa'=-(A0+Omega_s J)kappa`,

`a'=-(A0+Omega_s J)a`

`    +2 kappa[kappa dot A0 a]/|kappa|^2-nu|kappa|^2 a`,

`kappa dot a=0`.       (1)

The frame-rotation terms and the pressure projection are both included. Initially the flat-envelope part of the actual seed gives `kappa(0)=(k0,h)` and `a(0)=lambda(-h,k0)`. This matches the local phase gradient and velocity polarization, not the entire spatial jet of the curved, localized seed.

## 2. Closed solution, including viscosity and pump dilation

Write

`S=2Omega_s h`, `K(t)=k0+St`, `Q(t)=K(t)^2+h^2`, `Q0=k0^2+h^2`.

Then

`kappa(t)=exp(-ct)(K(t),h)`.

The scalar vertical perturbation-vorticity amplitude `zeta=kappa_r a_theta-kappa_theta a_r` obeys

`zeta'=-(2c+nu|kappa|^2)zeta`.

Consequently, with

`I(t)=integral_0^t exp(-2cs)Q(s) ds`, `D(t)=exp[-nu I(t)]`,

the exact polarization is

`a(t)=lambda Q0 exp(-ct)D(t) (-h,K(t))/Q(t)`.       (2)

For c>0 the integral is explicit:

`I=Q0 J0+2k0 S J1+S^2 J2`,

`J0=(1-exp(-2ct))/(2c)`,

`J1=[1-exp(-2ct)(1+2ct)]/(4c^2)`,

`J2=[1-exp(-2ct)(1+2ct+2c^2t^2)]/(4c^3)`.

At c=0 this reduces continuously to `I=Q0 t+k0 S t^2+S^2 t^3/3`. With c>0 its infinite-time limit is finite, because horizontal pump expansion eventually reduces physical wave numbers. The stress nevertheless vanishes through pump dilution and polarization. The short-time cubic damping diagnostic should not be extrapolated unchanged to arbitrarily large times when c>0.

## 3. Exact covariance and its monotone depletion for this seed

Average over the plane-wave phase. Define the useful positive shear-stress magnitude

`sigma(t)=-<w_r w_theta>=-a_r a_theta/2`.

Equation (2) gives

`sigma(t)=(lambda^2/2) Q0^2 exp(-2ct-2nu I) h K/Q^2>0`,

`sigma(t)/sigma(0)=exp(-2ct-2nu I) [K/k0][Q0/Q]^2`.       (3)

Its exact logarithmic derivative is

`sigma'/sigma=-2c-2nu|kappa|^2`

`             +S [h^2-3K^2]/[K(K^2+h^2)]`.       (4)

The actual range satisfies `h^2<(1600/63)^2<3k0^2`. Since K increases, every term on the right side of (4) is negative for all t>=0. **The stress never changes sign, but its magnitude decreases strictly.** At c=nu=0 it decays asymptotically like `t^(-3)`.

When the reduced radial carrier has doubled, `K=2k0`,

`sigma/sigma(0) <= 2[(k0^2+h^2)/(4k0^2+h^2)]^2 < 7/16`.

Thus more than 9/16 of this model's initial covariance magnitude is lost during the first winding, even before attributing any loss to viscosity. The wave energy is then less than half its initial value. This is the missing polarization effect in the passive-phase calculation.

The conclusion is specific to the initial carrier ratio. An inviscid pump-free wave with `0<k0<h/sqrt(3)` initially increases its covariance until `K=h/sqrt(3)`, then loses it. With pump and viscosity, initial increase requires the positive shear term in (4) to exceed `2c+2nu Q0`. This observation permits a different transient seed design, but changes the already certified datum and does not remove the integrated budget below.

## 4. An exact integrated-stress ceiling

Let the phase-averaged perturbation kinetic energy density be

`E(t)=|a(t)|^2/4=lambda^2 Q0^2 exp(-2ct-2nu I)/(4Q)`.

Directly from the pressure-projected amplitude equation,

`E'=-2c E-2Omega_s sigma-2nu|kappa|^2 E`.       (5)

Hence, for every T>0,

`2Omega_s integral_0^T sigma dt`

` +2c integral_0^T E dt+2nu integral_0^T |kappa|^2 E dt=E(0)-E(T)`.

In particular,

`integral_0^infinity sigma dt <= E(0)/(2Omega_s)=lambda^2 Q0/(8Omega_s)`.       (6)

In the pump-free inviscid case the finite-time identity is explicitly

`integral_0^T sigma dt=lambda^2 Q0/(8Omega_s)[1-Q0/Q(T)]`,

and the bound becomes equality as T tends to infinity. With the actual c>0 and nu>0 the bound is strict. The stress in this model consumes a finite initial perturbation-energy reservoir; it is not renewed by phase winding.

For the actual allowed h and Omega_s, the total budget measured in units of the initial stress obeys the conservative rational bound

`[integral_0^infinity sigma dt]/sigma(0)`

` <= (k0^2+h^2)/(4Omega_s k0 h)`

` < 10369/7862400 < 33/25000 = 0.00132`.       (7)

This is an infinite-time ceiling **inside the affine Kelvin model**, not a bound on the actual compact solution.

## 5. Connection to the actual initial NS energy identity, and what is still missing

There is an independently checkable actual-solution counterpart at t=0. Let `P0` be azimuthal averaging and `K_fluct=||u-P0u||_2^2/2`. For the actual retuned datum with outer seed `lambda w`, orthogonality, incompressibility, and integration by parts give

`K_fluct'(0)=-2c K_fluct(0)-nu lambda^2||grad w||_2^2`

` +pi lambda^2 m k0 A [integral r C'(r)e(r)^2 dr]`

`                         *[integral q(z)^2 q_v(z) dz]`,       (8)

where `v_theta=r C(r)q_v(z)`. This last term is negative. The exact envelope plateaus give

`-integral r C'e^2 dr>63/800`, `integral q^2q_v dz>=9/10`,

and the known energy bound gives

`K_fluct(0)<=pi lambda^2(6764/75)(12/5)/2`.

Therefore the actual initial fractional energy drain satisfies

`K_fluct'(0)/K_fluct(0)<-2953517/67640<-43.66`,

even after omitting the additional viscous drain. This validates the **initial direction** of the Kelvin depletion using the complete finite-energy NS datum. It does not propagate the rate or prove any finite-time drain fraction.

Finally, covariance and torque are different observables. The actual azimuthal-mean angular-momentum source contains the divergence of stress,

`-(1/r)partial_r[r^2 overline(v_r v_theta)]-partial_z[r overline(v_z v_theta)]`.

A single homogeneous Kelvin wave does not supply these spatial-envelope derivatives. Decreasing stress at one frozen ray does not prove decreasing actual mean torque, and the actual perturbation can acquire vertical velocity and other modes. A subsequent actual-stage theorem must control this spatial stress divergence and show how the perturbation-energy reservoir is replenished, or else budget its finite use together with central pressure feedback. Equations (3)--(7) identify the polarization/depletion term that such a theorem must retain.

`kelvin-covariance-check.py` verifies the exact ODE solution, energy identity, inviscid stress integral, and rational inequalities. It does not verify localization, a return, or a singularity.

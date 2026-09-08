# The outer carrier initially winds in the favorable direction, on a fast time scale

This is a **passive phase / local linearized diagnostic**, not a persistence proof for the finite-amplitude NS seed. The initial stationary-circle identities and the rational bounds below are exact. Extending their coefficients to positive time is an explicitly frozen model.

The PASS7 datum has `c=7/5`, `nu=1/1000`, `m=4`, initial radial carrier `k0=20`, and `780<A<1020`. At either axial packet center `z=±4`, write the unit swirl as `v_theta=r B(r,z)`. At every positive outer mean-angular-momentum maximizing circle,

`Gamma_mean=A r^2 B`, `2B+r B_r=0`,

`63/400<r_*<63/250`, `1/2<B_*<=1`, `B_z=0`.       (1)

The bound on B follows from the already proved maximum location `1/4<r_*^2/a^2<5/8` and the monotone cutoff, whose value at `5/8` is `1/2`.

## 1. Initial phase winding

For a phase advected by the **axisymmetric mean background** `U=c(r e_r-2z e_z)+A rB e_theta`, let `k=partial_r phi`, `ell=partial_z phi`, and `partial_theta phi=m`. Along a mean-flow ray the exact eikonal equations are

`D k=-c k-m A B_r`, `D ell=2c ell-m A B_z`, `D m=0`.

Thus the contribution specifically from swirl shear is

`S=-m partial_r(U_theta/r)=2m A B_*/r_*>0`,

and the complete initial radial derivative includes the pump dilation:

`k'(0)=S-c k0=S-28`, `ell'(0)=0`.       (2)

The angular wave-number component `m/r` has initial relative rate `-c`. From (1),

`260000/21<S<1088000/21`,

`259412/21<k'(0)<1087412/21`.       (3)

The initial winding therefore **reinforces** the carrier sign `mk>0` responsible for the proved positive mean-maximum derivative. Freezing the full initial derivative `k'(0)` gives a linearized radial-carrier doubling time between

`105/271853` and `105/64853`,

or conservatively between `1/2600` and `1/600`. This is roughly `0.000386–0.001620` in the datum's time units. It is not a proved time to doubling in the actual solution.

For reference, if the prescribed mean field were held steady, the exact passive formula would be

`k(t)=exp(-ct)[k0-m integral_0^t exp(cs) A B_r(r_*exp(cs),z_*exp(-2cs)) ds]`.

Even this formula evolves the mean coefficients along the ray; constant S is a further simplification. The actual mean NS field also evolves, and the finite seed contributes to the full particle velocity.

## 2. Viscosity initially lags phase winding, then becomes relevant

In the frozen shear diagnostic with fixed radius and `k(t)=k0+St`, scalar Fourier damping has exponent

`D_nu(t)=nu[(k0^2+(m/r_*)^2)t+k0 S t^2+S^2 t^3/3]`.       (4)

Here `exp(-D_nu)` is the passive viscous amplitude factor. The pump and later background evolution are omitted from (4). Using the largest allowed S and smallest radius proves the conservative exact model bound

`D_nu(1/600)<22357/2551500<1/100`.

Thus throughout the initial winding window up to `1/600`, this frozen damping model has less than a one-percent viscous amplitude loss, even with this deliberately unfavorable combination of bounds. This is a favorable early-time diagnostic, not a lower bound for the actual seed amplitude. The linearized time estimate above freezes `S-ck0`; formula (4) deliberately uses the larger pump-free S in its damping comparison.

The same frozen model eventually accumulates the cubic phase-mixing cost. Exact rational substitutions give

`D_nu(9/1000)<20550697/27562500<1`,

`D_nu(3/100)>536479/330750>1`.

Consequently its first one-e-fold damping time lies between `0.009` and `0.03`. The positivity of winding cannot be treated as free persistent gain: increasing radial frequency incurs a growing viscous cost. These two endpoint statements belong only to the frozen model; they do not bound the actual NS damping time.

## 3. Polarization and pressure act on the same fast scale

At a maximizing circle, (1) gives the exact initial **mean** gradient in the moving orthonormal cylindrical frame:

`grad U = [[c,-Omega_s,0],[-Omega_s,c,0],[0,0,-2c]]`,

where `Omega_s=A B_*` and `390<Omega_s<1020`. Its transverse part is a strong symmetric strain, not a rigid rotation. The apparent discrepancy with the slower central core strain is real: the outer angular-momentum maximum has `partial_r(rU_theta)=0`.

For an affine-background Kelvin wave, the wave vector and velocity polarization obey

`kappa'=-[grad U]^T kappa`,

`a'=-(grad U)a+2 kappa[kappa dot (grad U)a]/|kappa|^2-nu|kappa|^2 a`,

`kappa dot a=0`.

The second term is the pressure projection required to preserve incompressibility. These are exact for the affine Kelvin model; here they describe the local leading-order linearized diagnostic only. In the rotating cylindrical frame the rotation of the basis supplies the additional terms that reduce the radial equation to (2).

Since the strain is of size `Omega_s`, polarization can change by order one on the same time scale as the first winding. In fact, ignoring the small pump correction, `Omega_s*(k0/S)=k0 r_*/(2m)` lies between `63/160` and `63/100`. There is no small parameter justifying a frozen polarization during the predicted doubling. The actual perturbation can also develop vertical velocity through its nonlocal pressure; the compact form `curl[f e_z]` is an initial representation, not a preserved ansatz. Therefore the initial plateau identity `mean torque=mk e^2q^2/(2r)` cannot simply be propagated by replacing k with `k(t)`; outside that flat plateau the radial derivative of the covariance is also required.

## 4. The next actual-stage requirement

The combined PASS7 initial coefficient bound gives `beta''(0)>435297/140000`, while `b'(0)=|Omega_core|'(0)=2`; smooth existence yields some positive interval of favorable central feedback and mean-maximum gain. No enclosed duration has yet been shown comparable to the winding or damping times above. The initial derivative signs do not establish this comparison.

A concrete next lemma must control the **actual seed covariance and pressure** through a time window of order `10^(-3)`: prove that the transported negative covariance `overline(v_r v_theta)` still supplies positive mean torque at a tracked receiver, while the full central pressure feedback remains favorable and the pressure-generated vertical component is controlled. It must retain the strong background-shear polarization and the accumulated viscous term. A scalar phase calculation alone cannot supply either pressure or covariance persistence.

The diagnostic offers a favorable sign but a demanding time-scale requirement. It does not provide an obstruction theorem, an actual finite-duration estimate, or a repeated stage.

`verify_next_stage_shear.py` checks the exact maximum-circle shear bounds, gradient matrix, and rational frozen-model comparisons.

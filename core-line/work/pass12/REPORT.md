# Pass12: adapted geometry, a repaired gradient passage, and a whole-space reconstruction

8 September 2026. **The Navier–Stokes objective remains open.** This pass
supplies conditional mathematical estimates and an implemented H7
whole-space endpoint approximation. It does not validate an NS stage,
retained actual gains, an inherited return, or a singular construction.

The fixed datum remains `u0=M+1020v+(1/4)w_L`, viscosity `1/1000`.
The input numerical endpoint is the completed pass11 half-time-step run
at `T=.001`. No new datum, packet, phase alignment, or amplitude retuning
is used. This pass contains no new NS time integration.

## 1. What transfers from Claude's reference idea

The [C4 reference note](c4-adapted-reference.md) matches the **actual**
central gradient `A(t)=grad u(0,t)` and evolves `F'=AF`. For the actual
material flow, `Y=F^-1 Phi` satisfies exactly

`Y'=F^-1 [u(FY,t)-A(t)FY]`.

This removes a rate-mismatch term. Odd C4 symmetry makes the central
remainder cubic and the central matrices commute, while preserving the
whole nonaxisymmetric flow and its pressure/mean feedback. The outer
receiver has strong symmetric shear and still requires full matrix
evolution with time ordering.

A useful new connection avoids asking the H4 error bound to control a
third derivative in L-infinity. The explicit endpoint Sobolev estimate is

`[D²e]_(C^(0,1/2),op) <= (4 sqrt(3)/pi) ||e||_D4`.

For odd error this bounds its non-affine part by
`(16 sqrt(3)/(15 pi)) ||e||_D4 |x|^(5/2)`. Combining it with the
reconstruction's cubic remainder gives an explicit mixed-power material
radius and C1 bound. Central gradient errors also enclose F and its phase.
This is compatible with the existing H4 validation bridge; no actual H4
radius or useful finite-slab inputs have yet been certified.

The transformed field equation retains pressure, anisotropic diffusion,
the `2Bv` coupling, and the evolved viscous linear jet. A convenient
material coordinate system does not remove those terms from NS.

## 2. A specific repair for the parallel Claude argument

The [gradient note](comoving-gradient-decay.md) and its
[independent audit](comoving-gradient-independent-audit.md) repair the
location-dependent gradient estimate in the reviewed live draft.
For `u_t+b.grad u=nu Delta u`, translate a ball along the actual drift.
If `||Db||op<=L`, a local Bernstein argument gives

`|grad u(x,t)| <= M_Q sqrt((1+2Lh+max(4n,24)nu h/r²)/(2nu h))`.

The constant does not depend on the ball's spatial center, the magnitude
of the drift there, a drift Hessian, or a time-Hölder modulus. It transfers
a scalar polynomial tail to a gradient tail and supplies L1 finiteness.

The note also treats the initial-time passage for smooth bounded-gradient
data. The backward drift–diffusion representation and its spatial flow
derivative give a weighted gradient bound uniform down to zero. For
`G0=sup |x|²|grad u0|` and `B=||grad u0||infinity`,

`sup |x|²|grad u(t)| <= exp(3Lt) [sqrt(G0)+(B0 t+2 exp(Lt)sqrt(2n nu t))sqrt(B)]²`.

Here `B0` bounds the drift at the origin. The solution class,
representation/uniqueness premise and global Lipschitz bound are explicit.
The proof gives convergence of weighted norm values to G0, not an
unqualified weighted difference estimate. Nonsmooth angular data still
need their approximation passage. The drift bootstrap, general versus
affine covariance issue, receiver placement and complete S3 assembly
remain separate. These lemmas do not extend a classical NS solution
beyond an unknown singular endpoint.

## 3. An actual whole-space reconstruction, with the initial field retained

The [reconstruction audit](whole-space-reconstruction-audit.md) describes
the implementation and its scope. The velocity increment is represented
by `curl(T e_z)+curl curl(P e_z)`. Each positive mode has potentials
`r^m B(r²)`, using degree-nine B-splines and compact outer zero extension.
An axial C-infinity window is applied **inside** the potentials, including
its exact velocity contribution `q' grad_h P`.

These potentials are C8/H9; the induced velocity is C6/H7. It is exactly
divergence-free and regular on the Cartesian axis, but is not C-infinity
and cannot serve as an all-order force assembly without further work.
The support cylinder `r<=8, |z|<=8` fits in a ball of radius `8 sqrt(2)`.
Its exterior pressure residual generally does not vanish.

The fit uses `U(T)-U_discrete(0)` from pass11. The approximate endpoint is
the **exact analytic u0 plus this reconstructed increment**. Thus no fit
of u0 is substituted. This choice does not cancel the initial MAC
projection's consequences for the later trajectory; they must appear in
the complete evolution/residual error analysis.

All axial Fourier coefficients and angular modes 0,4,8,12,16 are retained.
The angular mode set enforces exact C4 equivariance. Stored coefficients
preserve inversion oddness to floating roundoff. An exact odd projection
is available when applying the symmetric comparison lemma; its residual
must be evaluated for that projected reconstruction.

## 4. The failed fit is evidence, not a successful endpoint

A fit minimizing weighted errors only at MAC nodes passed its linear
algebra checks but concealed between-node oscillations. Increasing the
radial trial space from 96 to 144 intervals reduced its sampled error
while making its continuous energy and central jets much worse.

| Endpoint reconstruction | Sampled increment fit error L2 | Fluctuation energy K | Core rotation Omega | Core strain b |
|---|---:|---:|---:|---:|
| Sampled least squares, 96 intervals | .745874 | 6.913278 | 1.001992 | .938180 |
| Sampled least squares, 144 intervals | .471716 | 13.210425 | .214239 | 45.383929 |
| Continuous L2 fit, 96 intervals | .745918 | 6.911374 | 1.001989 | 1.002167 |
| Continuous L2 fit, 144 intervals | .473013 | 6.901887 | 1.001995 | 1.002126 |

These are approximation diagnostics, not values of an actual whole-space
NS solution. The failed fits and their full coefficient files are retained.
A separate exact compact-potential counterexample has zero velocity at
every sampled radius but nonzero core rotation and energy. It explains
why a small sampled error cannot control the quantities being claimed.

The corrected selector uses the **continuous** velocity L2 Gram matrix
against a specified parity-preserving cubic radial interpolant of the MAC
increment. Poloidal and toroidal fields decouple by horizontal integration
by parts. Gauss quadrature on the union of potential and target knots
integrates these polynomial forms exactly in exact arithmetic; the solves
still use floating point. The interpolant itself is only an L2 target,
not an axis-regular Cartesian reference for higher modes.

## 5. What survives the corrected reconstruction, and what does not follow

For the continuous 144-interval reconstruction, numerical quadrature gives
`K=6.9018872459`, compared with `K(u0)=6.4437717666` evaluated by quadrature
from the exact initial field, and
`b-Omega=.0001300998`. The original receiving circle loses mean angular
momentum (`G=41.70834`); a selected circle near radius `.2198551` gives
`G=41.88696`. The selector is a bounded numerical optimization, not a
certificate of the global maximum.

Increasing axial quadrature from 4099 to 8199 nodes and radial Gauss order
from 36 to 40 changes K by about `7.1e-15`. This checks integration of one
fixed reconstruction. It does **not** measure its error against NS.

Changing the corrected radial trial space from 96 to 144 changes K by
about `.00949`, selected G by `.01103`, and the core margin by
`4.85e-5`. That last difference is comparable to the proposed small
margin. The axial localization alone changes the sampled increment by
about `.26044` in L2, predominantly in the mean. No H4 representation
error or full residual has been enclosed. Neither corrected fit supplies
a retained-gain certificate.

The reconstruction is an endpoint object. A time-slab reconstruction,
its time derivatives, full pressure, noncompact projected residual and
tail estimates are still required. A small normal-equation residual is
not a small NS residual. A finite positive pulse would also need a stable
successor class and a supplied reservoir before it could be iterated.

## 6. Next action

Use the continuous potential representation to construct a time-slab
approximation with controlled higher derivatives, preserving the exact
initial datum and the inherited endpoint. Price all five residual
components for the existing H4 hierarchy, including axial/radial joins,
omitted angular products and noncompact pressure. Evaluate the new
mixed-power core deformation test only with those certified inputs.

For the Claude path, apply the repaired gradient lemma under its exact
datum/solution hypotheses, then resolve the covariance and quantitative
strain-functional defect without restoring a sampled constant to a
proved column. Neither path currently supplies a completed S3 or a
repeatable ordinary-NS blowup mechanism.

See [verification and reproduction](VERIFICATION.md) for check scope.

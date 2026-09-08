# Independent audit: compact core plus remote pressure-induced stretching

2026-09-08. This audits the root agent's proposed local pressure calculation. It is a statement about the initial second time derivative of an actual smooth NS solution, not about a finite amplification factor or a repeated cascade.

**Verdict:** the signs and constants are correct under the stated definitions. The remote-packet lower bound has a small amount of spare margin. The second time derivative is independent of viscosity for this initial local geometry.

## Core pressure

Let `R^2=x^2+y^2+z^2`, `r_perp^2=x^2+y^2`, and

`u_c=Omega psi(R^2)(-y,x,0)`,

where `psi` is smooth, compactly supported as a radial profile, and equals one near zero. A prime on `psi` below means derivative in its argument `R^2`. The field is divergence-free and

`g_c=-Delta p_c=tr[(grad u_c)^2]=-Omega^2[2 psi^2+4 r_perp^2 psi psi']`.

For the decaying Newtonian pressure with `N=1/(4 pi R)`, the distributional kernel is

`partial_zz N = PV[(3 z^2-R^2)/(4 pi R^5)] - delta_0/3`.

The delta contribution is `-g_c(0)/3=2 Omega^2/3`. The radial `psi^2` term has zero angular principal value. With `mu=z/R`,

`average[(3 mu^2-1)(1-mu^2)]=-4/15`,

and `integral_0^infinity R psi(R^2) psi'(R^2) dR=-1/4`. Hence the remaining term is `-4 Omega^2/15` and

`p_c,zz(0)=2 Omega^2/5`.

No monotonicity of the cutoff is needed: its boundary values alone give the radial integral. Omitting the delta term would give the wrong answer and sign.

## Remote horizontal packet

Take a nonzero smooth compactly supported divergence-free `v=(v_x,v_y,0)` supported away from the core and inside `r_perp <= |z|/4`. Such packets exist, for example by taking a compact scalar bump in the interior of that region and setting `v=(partial_y chi,-partial_x chi,0)`.

Because the support avoids the evaluation point, integration by parts gives

`p_v,zz(0)=integral partial_zzij N(-x) v_i(x)v_j(x) dx`, with `i,j` horizontal.

The horizontal matrix of this kernel, multiplied by `4 pi R^5`, is

`(3-15 mu^2) I_2 + (105 mu^2-15) n_perp tensor n_perp`,

where `n_perp=(x,y)/R`. Its eigenvalues are

`lambda_tangent=3-15 mu^2`,

`lambda_radial=-12+105 mu^2(1-mu^2)`.

The support cone implies `mu^2>=16/17`. On that interval their respective maxima are `-189/17` and `-1788/289`; both are less than `-6`. It follows that

`c_v := -p_v,zz(0) >= [6/(4 pi)] integral |v|^2/R^5 dx > 0`.

The radial/tangential terminology refers to horizontal directions; at `r_perp=0` the rank-one term vanishes and the formula extends continuously.

## Actual NS second time derivative

For `u_0=u_c+A v`, disjoint supports imply `u_0 tensor u_0=u_c tensor u_c+A^2 v tensor v` exactly. Thus the decaying pressure satisfies

`p_0,zz(0)=2 Omega^2/5-A^2 c_v`.

In a whole neighborhood of the origin the initial velocity is rigid rotation, the initial vorticity is `(0,0,2 Omega)`, and `omega_t=0`. Differentiating the actual unforced vorticity equation at `(0,0)` therefore gives

`omega_tt(0,0)=2 Omega partial_z u_t(0,0)`.

The local nonlinear acceleration has zero z-component and the local Laplacian of the initial velocity vanishes. Consequently

`omega_z,tt(0,0)=-2 Omega p_0,zz(0)=2 Omega[A^2 c_v-2 Omega^2/5]`.

For `Omega>0` and `A^2 c_v>2 Omega^2/5`, the actual Eulerian axial vorticity initially increases at second order. A remote pressure gradient can move the origin's fluid particle later; it does not change this Eulerian second-jet calculation because the initial vorticity is spatially constant and its first time derivative vanishes locally.

Smooth compact initial data lie in every Sobolev space. Standard local smooth NS existence gives a genuine solution for some positive time for each finite `A`; enough space/time derivatives exist at the initial endpoint to justify the calculation and Taylor expansion. The guaranteed interval depends on the full data and can shrink as `A` grows. Therefore making this initial acceleration arbitrarily large does **not** by itself prove a uniform finite-factor vorticity gain, much less a same-solution regenerative sequence.

## Verification

`python3 work/pass2/verify_local_pressure_audit.py` verifies the pressure source, distributional-kernel regular part, angular constant, radial boundary identity, fourth-derivative horizontal kernel, and exact cone constants. The local time-jet derivation uses the actual PDE identities displayed above. No numerical pressure approximation or singular solution is used.

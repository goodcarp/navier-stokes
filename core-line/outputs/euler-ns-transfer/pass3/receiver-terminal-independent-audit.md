# Independent final audit: receiver flux and terminal profile

2026-09-08. Reviewed `receiver-flux-budget.md`, `terminal-profile-gap.md`, the definitions in `rotational-receiver.md`, and consistency with `material-core-transfer.md`. No original file was changed.

**Verdict:** no blocking correction is needed. The weighted flux/energy signs, harmonic coefficient, RMS expansion, implicit radius branch and first-order strain mismatch are correct in the stated local smooth setting.

## Weighted angular momentum and energy

For `phi=h(|x|^2)Jx`, direct differentiation verifies `div phi=0`,

`u_i u_j partial_j phi_i=2h'(u·x)(u·Jx)`,

`Delta phi=(4|x|^2 h''+10h')Jx`.

Integrating `u_t=-u·grad u-grad p+nu Delta u` against `phi` gives the positive advective sign and positive viscous transport sign printed in the note. The pressure term vanishes because the test is divergence free. At the initial endpoint, both the radial advective factor and the viscous pairing vanish exactly. Differentiating once leaves

`I''(0)=2Omega integral h'|Jx|^2 x·grad H`.

The harmonic radial-moment argument gives `I''(0)=Omega D H_zz(0)=Omega D K`. The verifier checks its quadratic coefficient; validity for arbitrary harmonic `H` follows from the full radial harmonic first-moment identity proved in `rotational-receiver.md` (or orthogonality of spherical harmonic degrees). Higher harmonic terms are not being approximated away.

The localized energy identity has the correct flux and dissipation signs:

`E_h'=integral (|u|^2/2+p)u·grad h +(nu/2)integral |u|^2 Delta h -nu integral h|grad u|^2`.

At zero, differentiating the actual weighted energy directly gives

`E_h'(0)=0`,

`E_h''(0)=Q_h+Omega I''(0)=Q_h+Omega^2 D K`,

with `Q_h=integral h|grad H|^2`. Thus the extra potential-flow energy term has the correct coefficient. It should indeed be included when the selected velocity scale is full RMS velocity.

## RMS and implicit radius branch

Since `E_r(0)=Omega^2D_r/2`, Taylor expansion gives

`E_r(t)/E_r(0)=1+[K+Q_r/(Omega^2D_r)]t^2+O(t^3)`.

Taking the positive square root yields exactly

`beta_r=(1/2)[K+Q_r/(Omega^2D_r)]`.

For `F(q,t)=q^(a+1)U_(qr)(t)/U_r(0)`, the initial rigid rotation gives `F(q,0)=q^(a+2)`. Therefore `F_q(1,0)=a+2`, `F_t(1,0)=0`, and `F_tt(1,0)=2beta_r`. The implicit branch satisfies

`q'(0)=0`, `q''(0)=-2beta_r/(a+2)`.

Hence the displayed `q(t)=1-beta_r t^2/(a+2)+O(t^3)` is correct. The argument is a one-sided local implicit-function argument at the initial time; it does not require solving NS backward. The receiver support is strictly inside the initial rigid neighborhood, which also permits a small open interval of radii around the chosen `r` for applying that argument.

## Terminal strain mismatch

The initial local acceleration is a harmonic gradient, so its gradient is symmetric. Consequently

`sym grad u(t,0)=t Hess H(0)+O(t^2)`.

Every rigid rotation has zero symmetric velocity gradient. Taking the symmetric part contracts the operator norm, so the comparison with **any** rigid rotation has gradient error at least `||sym grad u(t,0)||_op >= k0 t/2` on the quantified interval. Orthogonal coordinate changes preserve this norm; translations do not affect velocity gradients.

For the matched state `v(t,y)=q(t)^(1+a)u(t,q(t)y)`,

`sym grad_y v(t,0)=q(t)^(a+2) sym grad u(t,0)`.

Since `q(t)=1+O(t^2)`, the leading `t Hess H(0)` remains. The local lower bound `k0 t/4` for sufficiently small positive time is justified by `q(t)^(a+2)>=1/2`. The state comparison is expressly not claimed to be a time-dependent NS coordinate transformation. Exact RMS/energy/Re matching therefore does not establish return of the velocity profile.

## Consistency with material circulation

The receiver is fixed in space, whereas `material-core-transfer.md` transports its initial disk with the fluid. Weighted angular momentum/RMS gain can arise from changing material content and irrotational energy without increasing a transported disk's circulation. The material result uses the full three-dimensional Laplacian and proves only the first four circulation derivatives vanish. These claims are compatible; neither should be restated as a general no-go for receiver gain or as exact positive-time Kelvin conservation.

## Checks

`python3 work/pass3/verify_receiver_flux.py` passed. An independent SymPy calculation also checked the RMS square-root coefficient, the implicit radius coefficient, and the unchanged first-order normalized strain. The existing material-jet verifier had already passed. These are algebra/identity checks supporting the analytical arguments, not an infinite return-map verification.

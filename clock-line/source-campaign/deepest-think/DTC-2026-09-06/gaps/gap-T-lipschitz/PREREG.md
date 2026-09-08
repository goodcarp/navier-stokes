# PREREG — gap-T-lipschitz (written before any script was run)

Seat: `gaps/gap-T-lipschitz`, DTC-2026-09-06.  Target: GAP T of `lower/prove-lagrangian` §4(2).

## Claim to be proved
With `S = {rho0<|x|<R} ⊂ R^5`, `K(w) = -(3/8pi^2) w_z/|w|^5`, `T_lambda(y,z)=(lambda y, lambda^-2 z)`,
`a[f](0) = ∫ K(w) f(w) dw`, and `eta_0` supported in `S` with `|eta_0| <= M/r`:

  |a[eta_0∘Phi^-1](0) − a[eta_0∘T_lambda^-1](0)| <= pi*lambda*M*L*[ 3mu(1+muJ)/(2(1-mu)^5) + 3muJ/8 ]

under (H1) Phi a C^1 diffeo of S onto its image with J_Phi>0, (H2) |Phi−T_lambda| <= mu|T_lambda|, mu<1,
(H3) |J_Phi − lambda^2| <= muJ lambda^2.

## Registered kill criteria (a KILL FIRES if the stated thing happens)
* **K1** the shell integral `I := ∫_S |eta_0| |T_lambda x|^-4 dx_5` does NOT converge, or is not `<= pi^3 M L/lambda`.
  (Independent recomputation by 2-D quadrature must match the claimed closed form to <1e-8 relative.)
* **K2** the gradient constant: `sup_w |w|^5 |grad K(w)| != 3/(2 pi^2)`.
* **K3** on any of the 2 map families x 2 data x 3 (lambda,mu) points, the measured
  `|Delta|` EXCEEDS the bound computed with the MEASURED `(mu,muJ)`.
* **K4** the quadrature is not converged: doubling the node count moves any reported `Delta` by more
  than 1e-6 relative.
* **K5** the identity check `a[eta_0∘T_lambda^-1](0) = (M/2) lambda L` for the bang-bang cap fails at 1e-10.
* **K6 (control, must fire)** the linear small-mu corollary `|Delta| <= C0 mu lambda M L`, `C0 = 15pi/8`,
  must be VIOLATED by the registered control family (angular shear with an axis ramp of width phi_c,
  phi_c -> 0, at fixed mu of order 1).  If it is never violated the control has failed and I say so.

## Registered controls
* **C1** `Psi_c`: half-plane angular shear by `beta(phi) = mu_beta * m(phi_ax/phi_c)`, `m` the C^1 ramp
  `t^2(3-2t)` on [0,1] then 1; `phi_ax = min(phi,pi-phi)`.  `mu = 2 sin(mu_beta/2)`.
  Prediction: `|Delta| ~ (3/4) M L mu_beta^3 log(mu_beta/phi_c)` -> unbounded, while `mu` is fixed.
* **C2** `mu >= 1` (hypothesis (H2) void): the segment used by the mean-value step contains `w=0`.

## Fixed numerical settings (registered before running)
rho0 = 1, R = 4096 (L = log 4096 = 8.317766...), M = 1.
(lambda,mu) grid: (1.00, 0.05), (1.25, 0.10), (1.50, 0.20).
Data: D1 = bang-bang cap `omega^theta = -M sgn(z)`; D2 = 7.5-degree axis taper `-M sgn(z) min(1, phi_ax/delta)`.
Maps: A = radial ripple `u -> u(1+mu sin(2 pi log(|u|/rho0)/L))`; B = angular shear `beta = mu_beta sin^2 phi`,
mu_beta chosen so that `2 sin(mu_beta/2) = mu`.
Quadrature: Gauss-Legendre, panels split at phi = pi/2 (and at the taper corners for D2), n and 2n.

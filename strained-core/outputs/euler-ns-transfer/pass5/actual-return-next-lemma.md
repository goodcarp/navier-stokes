# The actual inherited-state return test after the certified episode

The certified radial datum now has a genuine finite interval of favorable feedback. It still does not return to its initial template under receiver shrinking and isotropic normalization. The new neutral tuning removes the old **first-order central-matrix** mismatch, but a definite **first-order cubic core deformation** remains. The exterior packets and pump also acquire shape changes which the same normalization cannot erase.

## 1. The successor must be the full evolved field

Fix the receiving weight and initial receiver radius `r0` inside the affine core. Let `U(r,t)` be its actual weighted RMS velocity. For a chosen exponent `0<alpha<1/2`, the local RMS matching branch satisfies

\[
 q(t)^{1+\alpha}\frac{U(q(t)r_0,t)}{U(r_0,0)}=1,
 \qquad q(0)=1.
\]

For the certified datum `U_t(r,0)/U(r,0)=2` at every allowed radius. The implicit function theorem therefore gives

\[
 \boxed{q(t)=1-\kappa t+O(t^2),\qquad
 \kappa=\frac2{2+\alpha}\in(4/5,1).}
 \tag{1}
\]

At a selected endpoint `tau`, the successor is exactly

\[
 \boxed{V_+(y)=q^{1+\alpha}u(\tau,qy).}
 \tag{2}
\]

This formula includes the entire inherited pump, swirl packets, pressure-generated velocity, diffusive tails and earlier field. No component is discarded or planted again.

For continuation, freeze this chosen `q`. The full rescaled solution

\[
 V_+(s,y)=q^{1+\alpha}u(\tau+q^{2+\alpha}s,qy)
\]

solves NS with

\[
 \boxed{\nu_+=q^\alpha\nu.}
 \tag{3}
\]

The physical viscosity of `u` is unchanged. Equation (3) is its dimensionless viscosity after a constant stage change of variables. A continuously time-dependent substitution of `q(t)` would introduce additional terms and is not the equation used here.

In particular, endpoint quantities obey

\[
 b_+=q^{2+\alpha}b(\tau),\quad
 \Omega_+=q^{2+\alpha}\Omega(\tau),\quad
 \beta_+=\beta(\tau).
\]

The initial certificate is at `nu=.001` and its viscosity bounds extend downwards to `0<nu<=.001`, but this does not make a general inherited field another certified radial template.

## 2. What the new neutral tuning really fixes

Let `M=D+J`. The exact central insertion jet is

\[
 \nabla u_0(0)=M,\qquad\partial_t\nabla u(0,0)=2M.
\]

Together with (1), this gives

\[
 \boxed{\nabla V_+(0)=M+O(\tau^2).}
 \tag{4}
\]

Thus isotropic RMS normalization now cancels the complete first-order central matrix change. This is an improvement over the earlier rigid-rotation candidate, whose first-order strain could not be removed that way.

Nevertheless, the certified strain-ratio acceleration is scale invariant:

\[
 \boxed{\beta_+=1+\frac12\beta''(0)\tau^2+O(\tau^3),
 \qquad\beta''(0)\ge\frac{290649}{17500}.}
 \tag{5}
\]

An exact return to `beta=1` already fails at second order. A useful enlarged class can admit a range of strain ratios, but it must track that range rather than reset the ratio to one.

## 3. A decisive first-order higher-jet obstruction

Write the actual radial strain profile as

\[
 \Phi(s)=\psi(s)+\frac75\eta(s),\qquad s=|x|^2.
\]

The core and annular strain together are `S_Phi`. The exact radial pressure formula from [the pass4 radial pressure identity](../pass4/core-pressure-derivative.md) gives

\[
 p[S_\Phi]_{zzzz}(0)=\frac{256}{7}\int_0^\infty(\Phi'(s))^2\,ds.
 \tag{6}
\]

The core radial swirl contributes no quartic pressure harmonic in the affine neighborhood. The outer swirl contribution satisfies, using its exact kernel from `local-swirl-feedback-kernel.md`,

\[
 p[Av]_{zzzz}(0)\ge-\frac{30 A^2C_v}{z_{\min}^2},
 \qquad z_{\min}=\frac{71}{20}.
\]

On the core transition `[1/4,1]`, `Phi=psi`, with endpoint change from one to zero. Cauchy–Schwarz gives the exact lower bound

\[
 \int(\Phi')^2\ge\int_{1/4}^1(\psi')^2
 \ge\frac1{1-1/4}=\frac43.
\]

All poloidal/azimuthal pressure cross sources vanish; core and outer swirl supports are disjoint. With exact tuning `A²Cv=204/35`, the **full initial pressure** therefore satisfies

\[
\boxed{
 p_{0,zzzz}(0)\ge K_4:=\frac{256}{7}\frac43
        -\frac{204}{35}\frac{30}{(71/20)^2}
      =\frac{3693184}{105861}>34.}
 \tag{7}
\]

This is an exact rational bound, independent of numerical pressure reconstruction.

Inside the initially affine core, both `nu Delta u0` and the third spatial derivative of `(u0·grad)u0=M²x` vanish. Hence the actual NS equation gives

\[
 \partial_z^3u_z(0,0)=0,\qquad
 \partial_t\partial_z^3u_z(0,0)=-p_{0,zzzz}(0).
\]

For the successor (2),

\[
\boxed{
 \partial_{y_z}^3(V_+)_z(0)
 =q^{4+\alpha}\partial_z^3u_z(\tau,0)
 =-p_{0,zzzz}(0)\tau+O(\tau^2).}
 \tag{8}
\]

Thus for sufficiently small positive `tau`, its absolute value is at least `K4 tau/2`.

Every fresh template in the original family, even after changing `b,Omega,c,A`, pump radii, or packet widths and distance, is exactly affine on a neighborhood of the origin. Its third core jet is zero. Therefore

\[
 \boxed{\inf_{U_{\rm template}}
 \|V_+-U_{\rm template}\|_{C^3(B_{r_*})}
       \ge\frac{K_4}{2}\tau}
 \tag{9}
\]

for a common smaller affine neighborhood and sufficiently small `tau`. In particular there is no return to that family with error `o(tau)`. Isotropic rescaling or a different affine strain amplitude cannot remove this cubic jet.

This does not rule out a class that includes the inherited deformation. The leading deformation is a harmonic cubic velocity field, so it is consistent with `Delta u_t(0)=0`: the initial curvature cancellations do not imply profile preservation. One concrete additional core coordinate is the coefficient of `−grad H4`, where

\[
 H_4=(35z^4-30|x|^2z^2+3|x|^4)/8.
\]

A compact divergence-free extension can be defined by taking the curl of `−chi(|x|²)x×(−grad H4)/5`, with `chi=1` near the core. This only supplies a candidate profile coordinate. Its changed pressure and interactions would need new estimates; higher harmonic terms still belong in the full residual.

## 4. The inherited exterior is a separate first-order condition

The companion `inherited-geometry-first-variation.md` tracks the actual positive upper swirl tag through the linear swirl equation driven by the **actual** meridional velocity. Its conserved angular-momentum measure gives an axial center `d`, transverse RMS radius `a`, and axial RMS width `h` with exact initial laws

\[
 d'=-2cd,\qquad (a^2)'=2ca^2+8\nu,
 \qquad(h^2)'=-4ch^2+2\nu.
\]

Consequently

\[
 \boxed{(\log(a/d))'=3c+4\nu/a^2>0,
 \qquad(\log(h/d))'=\nu/h^2>0.}
 \tag{10}
\]

These tagged shape ratios are unchanged by every isotropic normalization (2). They avoid treating viscous support as a transported compact set. Exact return to the same packet shape is therefore also impossible on the first small-time branch.

The same companion calculation finds meridional vorticity generated on the packet's axial transition:
`partial_t omega_theta=(1/rho)partial_z(w0²)`, which is nonzero. A single new pump amplitude does not reproduce that deformation. A return criterion must admit and control this inherited pump response as well as the changing packet moments.

The available cone and plateau margins can survive for a small interval; this is consistent with the finite-feedback corollary. A subsequent true return would require the integrated actual geometry rates to restore their allowed range, not just a new name for the outer scale.

## 5. A concrete complete-state return criterion

A proposed return class should consist of **whole fields with tags**, rather than a list of fresh radial profiles. One usable starting specification is:

- A finite parameter box for central rotation and strain, the cubic core coordinate, pump geometry/amplitude, and tagged packet centers and second moments. It must contain the motions in (5), (8), (10), with room on the appropriate sides.
- A profile residual `e=V−U_theta` measured on all of space, for example

\[
 \|e\|_X^2=\sum_{|\gamma|\le10}
       \int_{\mathbb R^3}(1+|y|^2)^{-2}|\partial^\gamma e|^2dy
       \le\epsilon^2.
 \tag{11}
\]

  The template may include the extra cubic coordinate, and the residual includes every remaining inherited field. A local norm alone would hide the exterior that determines pressure.
- Explicit local derivative bounds and tail-stress bounds sufficient to enclose the actual central pressure and its time derivative. These should be verified from the full field; defining the weighted norm (11) is not itself a proof of every needed pressure or lifespan estimate.
- Bounds on the actual feedback quantities and receiver coefficients, including the actual viscosity parameter `0<nu<=.001`. No later zero-curvature condition should be imposed merely because it held at insertion.

For clarity, the exact feedback tests can be evaluated on an incoming full field without pretending it is radial. Define

\[
 a_\nu(V)=\nu\Delta V-\mathbb P\operatorname{div}(V\otimes V),
\]

\[
 a^{[2]}_\nu(V)=\nu\Delta a_\nu
 -\mathbb P[(a_\nu\cdot\nabla)V+(V\cdot\nabla)a_\nu].
\]

Extract `(b,Omega)`, `(b1,Omega1)` and `(b2,Omega2)` from the central gradients of `V,a_nu,a_nu^[2]` respectively. The full tests are

\[
 \mathcal F_\nu(V)=\frac{b_1\Omega-b\Omega_1}{\Omega^2},
\]

\[
 \mathcal G_\nu(V)=\frac{b_2\Omega-b\Omega_2}{\Omega^2}
                   -2\frac{\Omega_1}{\Omega}\mathcal F_\nu(V).
 \tag{12}
\]

These are the actual `beta'` and `beta''`. At a neutral radial insertion they reduce to the certified formulas. At an inherited endpoint, using `−(pzz'+32b³)/(2Omega)` in place of (12) would omit terms that no longer vanish.

For example, an enlarged class could require central `Omega in [1/2,3/2]`, `beta in [1,3/2]`, `F_nu>=0`, and `G_nu>=delta/4`, together with positive receiver and swirl-ratio margins and explicit geometric/residual bounds. These are proposed conditions to test, not an assertion that this chosen box or its boundary is invariant.

The actual endpoint condition is then

\[
 \boxed{(q^{1+\alpha}u_V(\tau,q\,\cdot),\ q^\alpha\nu,
             \text{the rescaled actual tags})\in\mathcal K,}
 \tag{13}
\]

including (11), all pressure/curvature tests, and the inherited exterior bounds. The old field is not replaced by the best-fitting template; it remains `U_theta+e` in the successor state.

Global energy needs its own ledger. Under (2),

\[
 E[V_+]=q^{2\alpha-1}E[u(\tau)]
       \le q^{2\alpha-1}E[V].
\]

Since `alpha<1/2`, a uniform global normalized energy bound does not follow from physical energy dissipation. This is one reason to track remote history by a weighted norm and pressure tails, while preserving the exact physical finite-energy budget.

## 6. The next verifiable lemma

The next substantive target is a **robust return on inherited data**: choose an explicit class `K`, constants `tau_->0`, `tau_+<infinity`, and `0<q_-<=q_+<1`, and prove that **every** admitted complete incoming state evolves under its actual NS equation to (13) for some `tau in [tau_-,tau_+]`, `q in [q_-,q_+]`, with the same residual and geometry bounds.

Before attempting that theorem, the present calculations give two inexpensive rejection tests: a candidate class that insists on an affine core with residual `o(tau)` fails (9), and one that fixes the tagged packet aspect ratios fails (10). Enlarging the class must include those explicit deformations and then prove control at its boundaries; first-time membership in a small tube by continuity does not establish repeated return.

The present certified interval supplies initial positive feedback and small receiving gain. It does not yet supply a uniform contraction, a restoring geometry phase, or an invariant defect budget. The higher-jet lower bound (7) identifies one concrete new core variable that any next return calculation must carry.

`verify_actual_return_next_lemma.py` checks the exact higher-jet lower bound, normalization powers, first-order matrix cancellation, strain-ratio invariance, and the constant NS scaling law. The geometric identities are checked in the companion inherited-geometry verifier.

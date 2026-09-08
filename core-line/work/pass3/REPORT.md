# Third pass: a quantitative finite transfer into a smaller receiver

8 September 2026. **Ordinary Navier–Stokes blowup remains unproved.** This pass upgrades the preceding pressure-feedback onset to a finite velocity-gain estimate for one actual solution, and identifies why that estimate does not yet repeat.

## The new finite result

Use the compact rotating core and disjoint outer vortex packets from the [active-core pressure lemma](../pass2/active-core-pressure.md), with positive pressure margin

\[
K=A^2C_v-\frac25\Omega^2>0.
\]

The equation is ordinary three-dimensional unforced NS on R3, with fixed positive viscosity and smooth finite-energy initial data. No modified nonlinearity, prescribed affine strain, or evolving support separation is used.

The [quantitative core theorem](quantitative-core-stage.md) gives explicit profile-norm constants M2, M3 and an energy-method lifespan TE. For any known lower margin `0<k0<=K`, set

\[
\tau=\min\{T_E,k_0/M_2,\Omega k_0/M_3\},\qquad
g=\frac{k_0\tau^2}{4}>0.
\]

At time tau, central axial vorticity has grown by at least `1+g`. More substantially, every radial rotating test mode supported inside the initially rigid core has normalized velocity projection at least `1+g`. The remainder bound is uniform in receiver radius. Weighted Cauchy–Schwarz therefore gives a corresponding increase in actual local RMS velocity and weighted kinetic energy.

The constants explicitly include viscosity, all evolved core/packet interactions and the high derivatives of the localized profiles. They can certify a very small gain. This is a rigorous finite local estimate assembled from the analytical argument and standard local existence; the accompanying programs check its algebra and constants, not a formal proof assistant build.

## One smaller receiver can match the proposed scale relations

Let `U_r(t)` be the velocity RMS with a fixed nonnegative radial weight scaled to radius r, `Re_r=r U_r/nu`, and `E_r` the corresponding weighted kinetic energy. For any exponent `0<a<1/2`, the [receiver lemma](rotational-receiver.md) selects a radius ratio `0<q<1` at the receiving time such that the actual solution obeys

\[
\frac{\mathrm{Re}_{qr}(\tau)}{\mathrm{Re}_r(0)}=q^{-a},\qquad
\frac{U_{qr}(\tau)}{U_r(0)}=q^{-(1+a)},\qquad
\frac{E_{qr}(\tau)}{E_r(0)}=q^{1-2a}.
\]

At `a=1/4`, these are exactly the finite velocity/Re/energy exponents of our candidate budget. A fixed normalized profile/parameter class also gives explicit bounds `0<q_-<=q<=q_+<1` and a common positive lower gain.

The radius is selected from one evolved velocity field by a continuous scalar equation. This corollary describes one finite transfer; it is not a constructed sequence of shrinking solutions or an output-profile theorem. Its energy is the weighted energy of the receiving region, not a claim that the entire solution has been confined there.

## Where the receiving gain comes from

The [exact receiver budget](receiver-flux-budget.md) shows that pressure first generates a harmonic-gradient strain. Its interaction with the existing rotation then transports angular momentum into the fixed receiving region. Local kinetic energy also receives a contribution from the new strain velocity. All of this is compatible with the ordinary global energy inequality.

The receiving region is not material. For the original transported disk, the independent [material-core calculation](material-core-transfer.md) instead gives

\[
\Gamma_{\mathrm{material}}(t)=\Gamma_{\mathrm{material}}(0)+O(t^5).
\]

The first four circulation derivatives vanish for this particular datum. Its second-order vorticity increase is exactly cancelled by material-area contraction. The fixed receiver can gain by sampling additional neighboring material; that reservoir must be tracked in any restart.

## The remaining return problem is now explicit

The [terminal-profile calculation](terminal-profile-gap.md) shows that the normalized output has a symmetric strain of first order in time, while its rotational gain and selected radius change begin at second order. Isotropic rescaling and a coordinate rotation cannot remove that strain. Consequently the finite gain does not return the actual field to the original rigid rotating core at an error comparable to the gain.

A next proof must retain and control a larger class of core/outer states, or supply a genuine dynamical return that reduces this leading strain mismatch while retaining gain. It must regenerate the surrounding packets, control inherited fields and account for the material/energy reservoir. Smooth force may be used if the total all-order admissibility estimates close; merely forcing a fresh copy into place does not establish those estimates.

| Continuing requirement | Evidence after this pass |
|---|---|
| Actual coupled ordinary NS field | Established for the explicit finite family |
| Positive finite receiver velocity and Reynolds gain | Established, with explicit interval and a selected smaller receiver |
| Velocity/Re/energy scale arithmetic | Realized by actual weighted observables for one receiving time |
| Full usable terminal geometry and repeated transfer | Unproved; a leading rigid-core profile mismatch is identified |
| Smooth force for the finite episode | Exactly zero |
| One compatible infinite cascade, finite terminal time and singular continuation quantity | Unproved |

Four verification programs and independent reviews support the finite derivations; see [VERIFICATION.md](VERIFICATION.md). The full objective remains active. A finite receiver gain does not establish a Navier–Stokes singularity.

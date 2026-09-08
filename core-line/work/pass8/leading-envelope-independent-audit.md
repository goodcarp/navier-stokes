# Independent audit of the leading-envelope initial datum

**Accepted for the two initial growth claims and exact central pressure
neutrality.** I reviewed `leading-envelope-alternative.md`,
`certify_leading_envelope.py`, and `leading-envelope-certificate.json`,
independently checked the covariance, energy, and full viscosity algebra,
and reran the 1024-panel certificate successfully. No mathematical
correction was identified. The full central pressure-derivative gate is
outside this audit and is being checked separately.

## Torque and energy signs

For `psi=e(r)q(z)cos(m theta+kr)` and `w=curl(psi e_z)`, the exact components are

\[
 w_r=-\frac{meq}{r}\sin\phi,\qquad
 w_\theta=-e'q\cos\phi+keq\sin\phi,\qquad w_z=0.
\]

Thus `overline(w_r w_theta)=−mk e²q²/(2r)` even for the nonconstant
envelope: the derivative cross term averages to zero. The torque is
`mk q²(e²+2ree′)/(2r)`. With `m=4,k=−20`, covariance is positive and
the decreasing envelope makes torque positive at the previously isolated
mean-maximizing radius. The same envelope and chirality give positive
integrated work from the mean swirl, since `C′≤0` and

\[
 -\int w_iw_j\partial_j(Av)_i
 =\pi mkA\left(\int rC'e^2\,dr\right)
                \left(\int q^2q_v\,dz\right)>0.
\]

The mean meridional field is exactly `c(x,y,−2z)` throughout the seed
support. Its energy work is `−c||w||²`, or `−2cK` after inserting the
amplitude. The compact, solenoidal fluctuation has zero integrated
pressure work and self-transport work. These cancellations apply to
integrated energy, not to the local pressure or covariance equations.

## Complete viscosity, including cylindrical basis derivatives

The claimed norm identities are exact:

\[
 \|w\|_2^2=\pi R_0Z_0,\qquad
 \|\nabla w\|_2^2=\pi(R_1Z_0+R_0Z_1).
\]

An independent check used the full horizontal vector-gradient density

\[
 H=2r\overline{(\partial_rw_r)^2+(\partial_rw_\theta)^2
 +r^{-2}[(\partial_\theta w_r-w_\theta)^2
          +(\partial_\theta w_\theta+w_r)^2]}
\]

with `q=1`. If `J` denotes the proposed `R1` integrand, then symbolic
differentiation gives exactly

\[
 H-J=-\frac{d}{dr}\left[(e')^2+k^2e^2
                  -\frac{2m^2ee'}r+\frac{m^2e^2}{r^2}\right].
\]

The boundary term vanishes because the smooth envelope is compactly
supported away from the axis. The axial-gradient contribution is
`pi R0 Z1`. Therefore the certificate does not omit angular-basis
derivatives or axial dissipation. The `k²=400` and `m²=16` coefficients
in its `R1` are correct.

## Enclosures, pressure retuning, and symmetry

The new integration panels cover the full radial support `[1/25,6/25]`.
All envelope derivatives, including its two transition regions, use the
unchanged cutoff engine with endpoint-tail enclosures. Its adaptive
stopping criterion affects sharpness only; each contribution remains
interval width times a range enclosure. The rerun reproduced
`I_r>1/50`, `R1<200000`, and unit torque `>1950`.

The upper bound `R0≤392/5` correctly uses `||e′||∞≤80`, total variation
two, and `log6<9/5`. Two disjoint axial bumps give `Z0≤12/5` and
`Z1≤160/3`; their common plateau with the swirl gives `I_z≥9/10`.
Consequently all three certified margins are valid:

\[
 G_t(0,r_*,\pm4)>\frac{8871}{800}>11,
\]
\[
 \frac{K'(0)}{\pi\lambda^2}>\frac{205648}{375}>548,
 \qquad \frac{K'(0)}{K(0)}>\frac{12853}{2205}>5.
\]

The support lies strictly in the affine pump plateau and away from both
the axis and the central core. The horizontal pressure tensor is positive
on its cone. Its operator upper bound, applied to the **full** seed
velocity norm, gives

\[
 0<C_{\rm new}\le\frac{1764000}{1419857}<\frac32.
\]

With `lambda²≤1/16`, the certified `Cv` bounds imply `900<A_lambda<1020`
for the named exact amplitude. All base/seed bilinear central pressure
terms vanish by angular Fourier selection. Hence
`A_lambda² Cv+lambda² C_new=204/35` preserves `p_zz(0)=−8` exactly.

The scalar potential is even under full spatial inversion and invariant
under a quarter turn. Its curl is odd and C4-equivariant. Since it
vanishes near the core, the original affine central velocity is
unchanged. The symmetry makes the central pressure Hessian isotropic
in the horizontal plane; the neutral axial pressure therefore preserves
the stated central first-jet identities. No assertion about the new
second central jet follows from neutrality alone.

## Scope

For each datum with `1/8≤lambda≤1/4`, smooth local existence and the
strict initial margins imply simultaneous mean-maximum and fluctuation-
energy growth on some sufficiently short positive interval. Evaluating
the mean at an initially maximizing circle is sufficient; no moving
maximizer differentiability is needed. No numerical duration is proved.

This is a constructive **replacement initial datum**. It is not shown
to arise from the earlier trailing seed's actual evolution. Its steep
envelope cannot be propagated using the flat-envelope plane-wave model
without additional analysis. Persistence, full-state return, and
regeneration remain separate questions even if the independently checked
central pressure-derivative gate also passes.

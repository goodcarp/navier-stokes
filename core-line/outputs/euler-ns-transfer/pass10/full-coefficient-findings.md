# Full initial 3D coefficients: useful calibration, no evolved endpoint

The selected datum is still `M+1020v+(1/4)w_L`, with viscosity `1/1000`.
The new evaluator includes the complete compact radial pump, both axial
endcaps, the core, and all initial angular interactions. The root program
computes numerical approximations to `a1=N(u0)` and `a2=DN(u0)a1`.
It does **not** integrate the full NS solution in time. Its finite-cylinder
Taylor values are not accepted endpoint gains.

## Independent whole-space reference

Piecewise Gaussian quadrature at orders 32, 48 and 72 gives the following
reference values. Source/contact and stress calculations of pressure are
independent; tensor and separated-profile energy identities are independent.
These comparisons are numerical evidence, not interval error certificates.

| Actual initial quantity | Whole-space quadrature reference |
|---|---:|
| Central `p_zz` | −8.70903520263818 |
| `beta'(0)=−4−p_zz/2` | 0.35451760131909 |
| Total kinetic energy | 81037.4265803467 |
| Total energy derivative | −4611.6270902433 |
| Fluctuation energy `K` | 6.44377176657510 |
| Fluctuation energy derivative `K'` | 447.4261166273 |
| `K'/K` | 69.43543825499 |

The compact radial pump alone carries about 77623 units of kinetic energy.
Leaving it out would miss most of the full field's energy and a substantial
part of its nonlocal pressure. The pass9 rational initial-sign inequalities
remain the certificates; the decimals above do not replace them.

## Pressure solver and axial-domain comparison

The scalar pressure solve uses Fourier-periodic axial coordinates and the
decaying Bessel exterior condition in radius. For the initial pressure and
its first time derivative, the exact source is compact inside radius six;
placing the Robin boundary beyond that support is valid for the continuous
periodic problem. Finite-cylinder **velocity norms** still omit harmonic
tails. This boundary argument does not apply automatically to later evolved
pressure sources.

The legacy finite-volume solver was followed by separately discretized
gradients. The new sixth-order radial collocation solver improves pressure
accuracy and removes the deliberate source compatibility repair. Neither
is a commuting discrete Leray projection or a smooth full-space
reconstruction. Manufactured solutions expose that difference even when
the scalar solve residual is tiny. Axis parity also does not enforce every
required higher-order vanishing condition.

At radial resolution 768 and axial spacing `1/128`:

| Axial period | Numerical central `p_zz` | Independent image contribution |
|---|---:|---:|
| 16 | −8.22885020057606 | +0.48013346531463 |
| 32 | −8.70060279096940 | +0.00838318452567 |

The image calculation predicts the periodic value by adding its contribution
to the whole-space reference. Residual discrepancies are about `5.15e-5`
and `4.92e-5`, respectively. The image note separately bounds the omitted
image **tail**; it does not certify the Gaussian quadrature of the included
images or these discretization errors. Refining a fixed axial period cannot
remove the image contribution.

The refined initial mixed-mode radial pressure gradient is approximately
`−46.86124+55.05736i` after division by the genuine Fourier factor
`lambda*A/2`. The new full-gradient certificate from pass9 bounds the
whole-space leading-pressure approximation error by 1.3. This numerical
agreement is a consistency test, not a new enclosure.

## Derivative failure and its diagnosis

On the 384×1024 grid, direct extraction from the quadratic field gives
`beta''≈−553.75`, whereas the pressure identity gives `+57.06` for the
same periodic computation. A read-only independent audit finds the cause:
the global axial derivative of the large outer convection product leaks
into the central derivative. The convection contribution to `b''` is about
`−619.24`; a local derivative estimates the underlying derivative at 16.89,
close to the central identity's 16.91. The two viscous contributions are
only about −0.073 and −0.045.

Refinement improves this mismatch but does not eliminate it. The
768×4096, period-32 run has pressure-based `beta''≈38.96` and directly
extracted `beta''≈50.57`. Its numerical `div(a1)` L2 norm is still about
920. No sign or error enclosure for the evolved core is inferred from the
raw reconstructed polynomial.

The finite-volume pilot has an additional structural inconsistency: it
repairs the source but then uses the unrepaired source in the later
Laplacian identity. Also `2 tr(grad u0 grad a1)` equals the divergence of
the linearized convection only when `div(a1)=0`. That pilot is retained
as a rejected diagnostic, not an accepted NS coefficient computation.

## Rotation and the quadratic field

At time 0.001, the literal quadratic polynomial has fluctuation energy
about 261, starting from about 6.44. This is not evidence of fortyfold
growth. The phase `m*A*t` can be about four, outside a small-angle
quadratic approximation to rotation. Even the scalar quadratic polynomial
of a unit complex exponential has spurious amplitude growth.

Applying the exact global rotating-frame jet identities reduces the
period-32 polynomial's fluctuation energy to about 69.87 at rotation rate
1020, or 24.04 at the finite-cylinder L2-optimal rate about 755.16. The
remaining spread is substantial. One global rotation does not remove
differential rotation or meridional dynamics; these values are still
uncertified polynomials. Mean observables are unchanged by this rotation.

The computed initial radial mean acceleration at the receiver is about
6635 and the initial axial mean strain derivative about −31493. These
large full-field terms are absent from the prescribed-strain slice model.
That model cannot be promoted to the present compact 3D evolution.

## Required next step

Use a pressure and velocity representation with enforced axis regularity
and a controlled divergence-free reconstruction. Integrate the full viscous
evolution, or construct a sufficiently accurate rotation-aware approximation,
including generated modes and evolving axial and meridional fields. Bound
the full space-time residual, reconstruction, domain and tail errors in a
validator whose growth constants are usable. Then test a retained endpoint
and its successor-class membership. None of those endpoint obligations has
been completed here.

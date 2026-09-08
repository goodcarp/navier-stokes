# Independent audit of the full initial Fourier field

Accepted as a floating-point evaluator of the stated analytic initial datum.
This audit does not turn the evaluator into an interval certificate or a
time integrator.

The independent comparison with `cylindrical-NS-representation.md` checks:

- The radial strain extension and compact core swirl retain all cutoff
  derivative terms. The annular profile is the original compact profile;
  its difference-of-cutoffs expression is valid because the transitions are
  disjoint. Both axial packets and every radial/axial transition are present.
- The leading seed uses the actual phase `4 theta - 20 r` and a genuine
  Fourier coefficient `lambda/2`. Its radial and axial jets include the
  envelope derivatives and the inverse-radius derivatives.
- The gradient is the full cylindrical covariant gradient, with output
  component rows and differentiation-direction columns. Multiplying it by
  the velocity gives full convection, including basis-curvature terms.
- The pressure source is `trace(G_a G_b)`. The mode-zero mixed contribution
  is `2 trace(G_4 conjugate(G_4))`, the positive mode-four term is
  `2 trace(G_0 G_4)`, and mode eight is `trace(G_4 G_4)`. These are the
  correct Fourier factors; the trace is not a Frobenius norm.
- The vector Laplacian has the correct radial/azimuthal coupling signs.
  The axis limits of the mean and initially supported-away-from-axis mode
  four are consistent. More generally, smooth C4 fields have helical
  orders three and five in mode four, so their gradient and vector
  Laplacian also vanish exactly at the axis. Evolution does not preserve
  the initial vanishing of mode four throughout an axis neighborhood.
- Componentwise conjugation reconstructs the negative cylindrical modes;
  summing twice the real part of each positive coefficient is correct.

The symbolic representation checker independently verifies the general
divergence, Laplacian and helical identities used in these conclusions.
The module's separate Cartesian finite-difference checks remain numerical
diagnostics. Floating evaluation, including its explicitly documented
underflow treatment at flat cutoff endpoints, is not outward-rounded
evaluation. Subsequent pressure and evolution calculations must generate and
retain additional angular modes as required; the initial set `{0,4,8}` for
quadratic quantities is not a closed dynamical truncation.

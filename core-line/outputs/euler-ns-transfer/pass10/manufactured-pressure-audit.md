# Independent manufactured pressure check

This is a numerical consistency test of `Cylinder` and `HighOrderCylinder`, not a rigorous NS evolution certificate. The independent program is `check_manufactured_pressure.py`; its four-resolution output is `manufactured-pressure-checks.json`. It completed in about seven seconds. No solver file was changed by this check.

For angular modes 0, 4, and 8, the exact pressure is

\[
p_m(r,z)=(r/1.4)^m f((r/1.4)^2)
 [1+0.3\cos(2\pi z/6)+0.2\sin(4\pi z/6)],
\]

with an additional complex multiplier \(1+0.5i\) for nonzero modes, and
\(f(s)=\exp[-s/(1-s)]\) for \(s<1\), zero otherwise. Its analytic source is \(-\Delta_m p_m\). Thus the test does not manufacture the source by applying either tested discrete operator. The source and exact gradient use independent elementary derivatives. The field is smooth at the axis, periodic in \(z\), and zero well before the outer radial boundary. Both the Bessel boundary condition and the zero-mode gauge therefore agree with the exact solution.

The runs use \(R=4,L=6\), the solver's mapped grid with mapping parameter 4, 64 axial points, and 96, 192, 384, and 768 radial points. Axial derivatives of the deliberately low-frequency manufactured field are spectrally resolved. These runs isolate the radial pressure solve and derivative handling; they do not test axial truncation of the research datum or replacement of whole space by a periodic cylinder.

At 768 radial points:

| Solver | Mode | Relative pressure L2 error | Relative gradient L2 error | Relative L2 of div(grad p)+source |
|---|---:|---:|---:|---:|
| FV | 0 | 2.58e-5 | 1.88e-4 | 1.23e-3 |
| FV | 4 | 1.16e-4 | 3.89e-4 | 1.51e-3 |
| FV | 8 | 2.17e-4 | 5.75e-4 | 1.75e-3 |
| High order | 0 | 2.04e-9 | 3.06e-7 | 2.24e-5 |
| High order | 4 | 1.66e-8 | 8.29e-7 | 3.22e-5 |
| High order | 8 | 4.96e-8 | 1.60e-6 | 4.40e-5 |

On the last grid doubling, the FV pressure errors have observed order 2.00 in all modes; the gradient errors have orders 1.91–1.94. The high-order pressure errors have observed orders 5.70, 5.82, and 8.28 (mode dependent), and its gradient errors have orders 5.35–5.39. These are observed orders on this smooth profile, not a theorem about the method or its asymptotic order for every field.

The independently applied solve operator has relative maximum residual below 2.5e-10 in every case. This is much smaller than the divergence-of-pressure-gradient mismatch. Both mismatches decrease under refinement, but neither implementation gives an exactly commuting discrete projection: the FV conservative face Laplacian differs from the separately differentiated pressure gradient, and the high-order composition of vector/scalar derivative matrices differs from its separately stored second-derivative matrix. An algebraic pressure solve residual alone is therefore insufficient to certify incompressibility or a full NS residual.

For FV mode zero the solver changes the source by its recorded compatibility repair. The algebraic solve check includes that recorded change; all physical consistency errors in the table compare against the original analytic source. At the final grid the repair amplitude is 3.96e-8. The high-order solver does not repair the source and reports a computed outer zero-mode flux of 1.36e-11. Neither value is an enclosure.

There is also an axis regularity warning that maximum absolute errors conceal. A smooth scalar angular mode has the form \(p_m=r^{|m|}F_m(r^2,z)\). At the first cell, the exact normalized radial coefficient is about one, whereas FV gives approximately 3.00 for mode 4 and 315.93 for mode 8; high order gives approximately 385.08 for mode 8. High-order mode 4 is accurate there. The absolute fields and the reported absolute near-axis derivative errors still shrink. Thus these tests do not establish an absolute-field nonconvergence, but they do show that even reflection alone does not establish the complete \(r^{|m|}\) regularity needed for a smooth Cartesian reconstruction. A future reconstruction/certificate must explicitly impose and control that regularity rather than infer it from a tiny pressure value at the axis.

Conclusion: the new high-order solve materially improves this manufactured pressure and gradient test. Neither solver has passed a continuum full-residual enclosure. The next numerical validation should keep the actual divergence/gradient mismatch, axis regularity, whole-space tails, and periodic-domain error in its residual budget.

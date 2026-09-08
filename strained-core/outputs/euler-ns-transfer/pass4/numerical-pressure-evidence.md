# Exploratory evidence for the full initial pressure-feedback gate

These are floating-point evaluations of an initial NS derivative, not a time simulation, interval-certified inequality, or blowup result. Both evaluated formulas retain the complete pressure contribution up to the stated numerical resolution.

## Datum and observable

The evaluator uses the exact compact profiles in affine-core-pressure.md and outer-pressure-reservoir.md. Its common C-infinity cutoff is one for s<=1/4, zero for s>=1, and on the transition is

    chi(s)=expit(1/t-1/(1-t)),  t=(s-1/4)/(3/4).

The core cutoff is psi(|x|^2)=chi(|x|^2). The outer pump parameters are sigma=0.7, height=1, distance=4. Its swirl occupies 0.45 times the pump's radial and axial cutoff scales, strictly within the plateau. The outer swirl lies inside r<=|z|/4; both outer fields are separated from the core.

Set b=Omega=1, vary c, and choose A by the exact neutral prescription

    A^2=[38/7+2/5-c^2 C_P]/C_v.

The numerical implementation estimates C_P,C_v by the outer stress formula before forming A; consequently its specified floating-point A is only approximately the mathematically exact tuned amplitude. A proof must bound these constants and tune the exact datum accordingly.

The tested scalar is

    T_nu=p_zz'(0;nu)+32.

At exact neutral tuning, beta=b/Omega satisfies beta''(0)=-T_nu/2. Thus negative T_nu is favorable. This conclusion concerns the initial second derivative and does not supply a persistent or reusable stage.

## Two equivalent whole-space evaluations

The source method solves the decaying whole-space Poisson problem using even Legendre modes, computes the full initial acceleration gradient, and evaluates the contact-plus-degree-two formula for p_zz'. There is no periodic box. Radial Green integrals are exact for the piecewise linear interpolant of each sampled source coefficient.

The split method computes the core contribution through its source and the outer contribution through the equivalent stress integral -2 integral Q_ij u_outer,i u_1,j. This reduces derivative cancellation on the narrow outer packets. Both methods share the approximate pressure solve; they are different integral representations, not independent NS simulations. The exact equality is audited in split-pressure-evaluation-audit.md. A pressure/divergence residual can cause their numerical disagreement.

## Refined results

At b=Omega=1 and c=2:

| Radial nodes | Angular nodes | Largest pressure degree | Source T_0 | Split T_0 | Viscous coefficient dT/dnu |
|---:|---:|---:|---:|---:|---:|
| 1200 | 1024 | 768 | -269.763831 | -271.621045 | 3027.147327 |
| 1600 | 1536 | 1280 | -271.686933 | -271.619999 | 3027.118162 |

The split value changes by approximately 0.00105 between these resolutions. At the finer resolution the two representations differ by approximately 0.06693. These observed differences are diagnostics; neither is a rigorous error bound.

The finer calculation gives C_P approximately -0.42697138636 and C_v approximately 0.00000627156153. For c=2 it gives A approximately 1096.21507043. The independently evaluated source expression for -C_v differs from the stress value by about 5.4 parts per million. Core pressure values agree with -18/7 and 2/5 to approximately 1e-11.

The relative maximum Poisson trace residual falls from about 0.00691 to 0.000559 for this datum. That residual measures angular source truncation. It does not bound radial interpolation error, every individual pressure-Hessian error, or the final scalar error.

At the finer resolution, the split formula gives the following initial gate values:

| Pump amplitude c | T_0 | T_0.01 including viscosity |
|---:|---:|---:|
| 0 | 13.306561 | 42.004603 |
| 1 | -99.702294 | -70.610967 |
| 2 | -271.619999 | -241.348818 |

For c=2, the source formula independently within the same pressure solve gives T_0.01 approximately -241.415752. Both signs are favorable. This is stronger numerical evidence than the earlier isolated swirl-channel gain because it includes the meridional pressure response and the viscous term.

The sufficient swirl-ratio derivative bound is also numerically favorable at c=2, nu=0.01:

    (2381/357)c -4 -nu d_v/C_v approximately 4.415.

The same check is not guaranteed by this sufficient bound for c=1 at nu=0.01. This distinction is retained when selecting the c=2 candidate.

Low angular resolution was misleading: an earlier underresolved source calculation returned the opposite pump sign. The refined results and the stress representation are why the candidate remains under investigation. No sign was accepted from the coarse calculation.

## Independent controls and exact identities

The core-only cubic derivative admits a separate one-dimensional reduction. It gives t300 approximately -6.73446331117 and tWb approximately 2.563122131694 for this cutoff. The general Legendre evaluator approaches their sum with the expected radial refinement; its core-only error is below 0.0001 at 1600 radial nodes. Two independently arranged one-dimensional integrals agree to approximately 5e-11. Those agreements test the implementation, not the full pump sign.

Exactly dS=dW=0 for these radial core profiles. For the outer fields, the viscosity coefficients are evaluated by

    d_j=2 integral Q_ik partial_l j_i partial_l j_k.

The finer numerical values are dP approximately -170.89870411 and dv approximately 0.00308791848. An independently derived second-velocity-derivative formula supplies an additional symbolic and moderate-resolution numerical check; see pressure-quadrature-method.md. Its unresolved coarse core cancellations are reported there rather than treated as nonzero core viscosity coefficients.

## What this changes, and what remains

The proposed meridional pump now has a concrete favorable numerical full-feedback example at positive viscosity. It has not yet passed a rigorous sign test. The next step is to bound source projection, radial interpolation, pressure reconstruction, amplitude tuning, and arithmetic errors tightly enough that an interval enclosure for T_0.01 lies below zero.

Even a certified initial sign would leave persistence, useful receiver gain, inherited profile control, restoration of the exterior environment, and same-solution iteration open. The sign test is not a substitute for those obligations.

Reproduce the finer experiment from this pass's directory with:

    python3 experiments/evaluate_pressure_gate.py --nr 1600 --nmu 1536 --lmax 1280 --c 0 1 2

NumPy, SciPy and SymPy are required. The run uses substantial memory for the resolved profile derivatives. Raw results are in experiments/pressure-split-1200.jsonl and experiments/pressure-split-1600.jsonl.

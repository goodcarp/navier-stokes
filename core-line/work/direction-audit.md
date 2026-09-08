# Independent audit of direction-scout.md

Date: 2026-09-08. Reviewed `work/direction-scout.md` and `work/verify_direction_scout.py`; did not modify either. Ran the verifier with `python3`: all checks passed. Separately checked the time-dependent rescaling derivative and the displayed rational exponents with SymPy.

**Verdict: accept the algebra and scaling calculations as necessary-condition diagnostics. They do not establish a viable nonlinear stage or a Navier–Stokes blowup route. Two force-residual scope clarifications below should accompany downstream use.**

## Accepted calculations

- The circulation and reduced-vorticity diffusion operators in volume coordinates are correct. Their common principal symbol is `q=zeta_1^2+r^2 zeta_2^2=r^2 kappa`. The Fourier-mode eigenvalues are `-nu n^2 N^2 q +/- lambda` for the stated frozen linear system. Resolving the actual nonsinusoidal periodic profile into modes is necessary and is handled explicitly in the scout.
- The physical-coordinate identity for `lambda^2` is correct. Cauchy–Schwarz gives the stated upper bound `2|Gamma| |grad Gamma|/r^3`. Here circulation and viscosity both have units length squared/time, so the thin-ring gate `G0/nu >= C M^2(R/L)^(3/2)` is dimensionally correct. For fixed aspect ratio, reducing the ring radius alone does not improve this gate. This conclusion requires the controlled-gradient and fixed-profile assumptions; it does not exclude all approaches to the axis.
- The nested-core ansatz gives stage time proportional to `rho^(2+alpha)`, core energy to `rho^(1-2alpha)`, core integrated dissipation to `rho^(1-alpha)`, and carrier/strain damping ratio to `rho^(alpha-2beta)`. A carrier of comparable amplitude adds `M^2` to the dissipation, producing exponent `1-alpha-2beta>0`. All are consistent with the stipulated exponent range. For `alpha=1/4,beta=1/16`, that final exponent is `5/8`.
- Geometric scale sequences make these budgets summable. This does not establish suitability, an energy inequality, or same-solution regeneration. Even for overlapping cores, the corresponding geometric sum of square roots of the energy budgets is finite, so overlap does not by itself create an elementary energy contradiction; actual cross terms still require estimates.
- The physical force derivative exponent `3+2alpha+k+(2+alpha)j` is correct for fixed stage scales. Costs of carrier derivatives must remain inside the normalized derivative norms, as the scout states.
- The dynamic rescaling formula is correct when physical pressure is written `p=a(t)^2 P(x/L(t),s(t))` and `s'=a/L`. The amplitude and dilation terms have the displayed signs and coefficients. They cannot generally be absorbed into pressure.

## Clarifications required for a force construction

1. The summability criterion certifies smoothness of a force **series** only when each `f_m` is a compatible smooth global increment (or has controlled smooth extension), and the series really is the residual of the assembled velocity. Residuals of separate stage solutions cannot simply be added: the residual of a velocity sum contains nonlinear cross interactions. Define `f_m` as a difference of complete residuals, or explicitly include every cross-stage and joining term. For Clay admissibility, also establish common bounded spatial support or the required weighted spatial decay.
2. The sentence suggesting exact solution of the principal viscous equation as an alternative to improving residual accuracy needs qualification. Exact principal cancellation alone leaves localization, material-coordinate, nonlinear interaction, and activation residuals. Those remaining physical-force terms must still satisfy the all-order scaled bounds. This is a scope correction, not a failure of the displayed derivative formula.

## Further limitation of the proposed one-stage gate

If the fine carrier has amplitude comparable to the core velocity `U` and wavelength `h=L/M`, its self-advection can be of order `U^2/h`, which is `M` times the assumed older-strain forcing scale `U^2/L`. Thus the linear amplification calculation requires a small carrier amplitude or a proved structural cancellation to remain valid. The scout already labels its Reynolds-number gate as necessary and calls for nonlinear-feedback control; the next lemma must make this quantitative.

There is no algebraic falsification of the nested-core scale budget here. There is also no proof of a nonempty correction-versus-diffusion window, compatible smooth forcing, or an exact PDE orbit. The term “scaling-compatible target” is justified; “viable blowup mechanism” would not be.

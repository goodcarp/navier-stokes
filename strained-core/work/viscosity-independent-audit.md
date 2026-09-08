# Independent audit of viscosity-audit.md

2026-09-08. Reviewed `work/viscosity-audit.md`, ran `python3 work/viscosity-check.py`, and compared the cited statements against the local Euler manuscript extraction, particularly printed pp. 48–49, 60–62, 86–89, 95–97 and 104. Original files were not modified.

**Verdict:** the algebra and the asymptotic conclusion are correct for the explicitly defined auxiliary model using the published Euler coefficient history and scale schedule. I found no blocking mathematical error in that conclusion. It must retain its existing limitation to that model; it does not establish a full viscous-PDE no-go for other histories or correction schemes.

## Fixed phase-symbol lower bound: supported

Section 9 fixes one containing region with `0 < c_- <= c_g <= c_+`, and fixed covector-length endpoints `z_- > 0`, `z_+`. Lemma 9.2 on p. 60 explicitly returns the prescribed covector-length intervals on growth, transition and return; p. 62 again uses the recovered lower endpoint in the central deformation estimate. Section 11 identifies `zeta_± = z_±` and fixes the margins before the scale sequence. Consequently, for the retained Euler history and current layer's growth interval,

`q'_m = zeta_1^2 + c_g zeta_2^2 >= min(1,c_-) z_-^2 = q_* > 0`

has a constant independent of the layer. The audit only needs this central bound; no assertion about arbitrary old layers on all later windows is needed. There is no missing frequency factor in this symbol. The fixed geometric magnification changes viscosity by a constant, as recorded.

This lower bound is conditional on using the already-constructed Euler history. It is not automatically available for the self-consistent flow obtained by actually adding viscosity. The audit's auxiliary-model definition correctly makes that distinction.

## Growth bound and schedule exponent: supported

The potentially easy-to-misread convention is `sigma_i^2`, not `sigma_i`, in equation (11.3):

`c_sigma N_i^beta <= sigma_i^2 <= C_sigma N_i^beta`.

Thus `sigma_i ~ N_i^(beta/2)`, as the audit correctly uses. Equation (12.4) supplies `s_m = Lambda_m sigma_(m-2)/sigma_(m-1)` and `Gamma_m^gr = sigma_(m-1) sin(s_m)`. Lemma 9.2 gives a uniform upper logarithmic-growth bound `g_m(t) <= C Gamma_m^gr` during the growth interval. Hence, for `m >= 3`,

`g_m(t) <= C Lambda_m sigma_(m-2) <= C Lambda_m N_m^[beta/(2 Q_m Q_(m-1))]`.

The two-generation conversion is exact because `N_m = N_(m-2)^[Q_m Q_(m-1)]`. Dividing the uniform principal damping lower bound by this uniform growth upper bound gives exactly the audit's lower bound with exponent

`2 alpha - 1/[16 Q_m Q_(m-1)]`.

For every fixed `alpha > 0`, this exponent is eventually at least `alpha`. Since `Lambda_m = (k_m + 8 - 7/8) x_m`, `k_m <= Q_m/c_Q`, and `x_m = log N_m` grows at least geometrically while `Q_m` grows linearly, `log Lambda_m = o(x_m)`. Therefore the stated ratio tends to infinity, uniformly over the current growing interval. Integrating over the published duration preserves the same conclusion. No numerical limit extrapolation is required.

## Model/PDE boundary: correctly stated, and essential

- The original carrier is nonsinusoidal. Its viscous evolution is not a closed damped version of the original fixed-carrier ansatz. Fourier-resolving the leading phase model, as the audit does, is necessary.
- In that model common scalar damping commutes with the time-dependent matrix. The projective ratio can satisfy the original Riccati equation while the absolute amplitude fails to grow. This distinction is correct.
- The conclusion about every fixed fractional order is a conclusion about the defined principal-symbol model and this particular frequency schedule. A fractional operator's full nonlocal effects, local force recovery, changed transport and nonlinear feedback were not analyzed. The text expressly says so.
- Separate from that model, the same published bounded velocity cannot become an ordinary NS singular solution with admissible smooth force. The displayed `H^1` estimate and continuation argument apply in the stated finite-energy smooth setting, including compactly supported forcing. This is the strongest direct PDE obstruction in the note.

## Verification result

`python3 work/viscosity-check.py` passed all eight checks: the two reduced diffusion operators, physical phase symbol, damped characteristic polynomial, Riccati cancellation, absolute amplitude damping, two-generation exponent, and ordinary-viscosity specialization. These checks verify algebra, not the validity of the Euler paper, its Lean build, or existence/regularity of a redesigned NS solution.

Minor clarity improvement if revising later: use a different letter for Fourier harmonic index and derivative allowance `k_m`, and call `q_*` the lower bound on the **retained Euler phase history** wherever the asymptotic result is summarized. Neither change affects the result.

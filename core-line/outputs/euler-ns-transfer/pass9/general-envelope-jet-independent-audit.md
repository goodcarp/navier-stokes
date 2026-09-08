# Independent audit of the leading-envelope second mean jet

**Accepted as an actual initial-jet enclosure.** I reviewed
`general-envelope-jet.md`, `certify_general_envelope_initial_jet.py`, and
`general-envelope-initial-jet-certificate.json`; checked the pressure and
covariance conventions against the earlier exact identities; and reran
both the symbolic checker and the 2048-panel interval calculation.
Both passed. Exact endpoint comparisons verify, uniformly over the
actual leading family,

\[
 G_{tt}(0,r_*,\pm4)<-13501,
\]
\[
 \frac{d^2}{dt^2}G(t,r_c(t),4)\big|_{0}<-12125,
\]

where `r_c` is the local radial critical branch at **fixed axial height**.
These are not statements about an unrestricted moving global maximum,
nor bounds on a turning time or on the duration of retained gain.

## Exact nonconstant-envelope terms

The covariance `R0=−lambda²mk e²/(2r)` still holds when `e′` is nonzero,
but differentiating its flux and the actual covariance equation changes
the pump and viscosity terms. In particular the partial-time pump term
in `R_t` is `c T0`; the complete pump contribution to `Gtt` is
`−2c B_lambda(2F′+rF″)`. Its sign cannot be borrowed from the flat
envelope. The code retains it and every displayed `J_nu` term, including
`e e‴` and the cylindrical angular-basis contributions.

For the pressure convention in the note, the pressure covariance factors
are `B=e′−ik e` and `−m²e/r²`. Direct differentiation of **both** radial
and axial covariance fluxes gives

\[
 \mathcal T_t^p=\frac{\lambda^2A}{2}
 \operatorname{Re}\{e^{-ikr}[D p_{m,r}+rB s_m]\},
\]

with the stated `D`. Substituting the full Poisson equation cancels
the undifferentiated pressure. The mixed source has the factor two
shown in the note. Its explicit real part, combined with differentiated
mean-shear production at `g′=0`, gives exactly the stated `S_local`.
The real/imaginary multiplication in the interval program has the
correct phase and sign for `k=−20`.

The pump/seed mixed pressure source is zero globally because the seed
is horizontal, divergence free, and wholly within the pump's affine
plateau. Core/seed local products are disjoint. Initial seed-self pressure
has modes `0,2m`, which pair to zero with one `m`-mode seed factor.
Consequently the retained unit swirl/seed mixed pressure supplies all
initial pressure contributions to these covariance derivatives.

## The pressure bound is new for this seed

No radial pressure norm or gradient error from the trailing seed is reused.
The program evaluates the new radial envelope, its negative carrier,
the Green moments, the weighted radial norm `integral r³|P|²dr`, and
the derivative of `P` anew. The local term `−2C b′` includes the
nonconstant envelope derivative.

Prefix and suffix ranges cover each variable Green limit by adding its
partial cell with length in `[0,h]`. Dependence between those quantities
can enlarge the enclosure but does not invalidate it. Outside the new
seed support, the exact pressure forms are `c_in r^m` and `c_out r^(−m)`.
Both corresponding weighted-square tails are included, with denominators
`2m+4` and `2m−4` respectively.

Only the old **axial** interval `integral |(q_s q_v)″|² dz` is reused.
Its defining axial widths, centers, and two-packet factor are unchanged.
I checked its source code: it includes both packets and the product
derivative `q_s″q_v+2q_s′q_v′+q_s q_v″`.

The previously proved harmonic-slab inequality applies structurally to
this new residual. The new physical real-mode residual norm obeys
`n_e²=pi (integral r³|P|²dr)(integral |F″|²dz)`. The factor `pi` cancels
the slab evaluation's `1/pi`, leaving precisely the program's coefficient
`(35/4)r^6/(m²(9/20)^9)` times the odd-series upper bound. The recomputed
upper bound is

\[
 |p_{m,r}(r_*,\pm4)-P'(r_*)|<1.298391<1.3.
\]

This is a bound on the full complex amplitude. No additional factor two
is needed. The pressure contribution is enlarged symmetrically by
`|D|` times this upper bound divided by two, so it contains the actual
pressure correction with either sign.

## Parameter enclosure and scope

The intervals `900<=A<=1020` and `1/64<=lambda²<=1/16` contain every
exactly retuned member of the leading family. Treating them independently
only enlarges the interval; their true correlation is not required.
All pump, viscosity, shear, and full-pressure terms enter that arithmetic.
The exact stored upper endpoints are below `−13501` and `−12125` as
claimed. The final script's `PASS` status requires exact rational endpoint
comparisons for the first derivative, both negative second derivatives,
and the new pressure-error upper-bound expression. The separate checker
replays those comparisons and validates the stored dependency paths.

The nonflat `T0′` is retained in `G_tr`. Since `A g″<0`, the fixed-height
radial branch adds the positive correction `−G_tr²/(A g″)`, with the
correct sign. The original mean radius is reused because the new seed
has zero mean; the pressure-neutral amplitude only multiplies the same
mean radial profile. No old-seed acceleration value is transferred.

These results can coexist with the leading family's positive initial
first derivative and energy gain. Without an evolution/remainder bound,
the negative second derivative does not decide whether a prescribed
finite endpoint gain survives. The H4/H7 validation alternatives in
`evolution-validation-bridge.md` address that separate obligation.

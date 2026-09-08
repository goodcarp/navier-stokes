# Verification scope

The new result is the common-family initial feedback theorem and its qualitative
simultaneous local-time consequence. It does not prove an inherited return, a
fixed useful gain, an infinite sequence or blowup.

Run python3 checks/run_checks.py. It runs seven finite programs:

1. Exact full D tensors, cone signs, and energy constants.
2. Recomputed 1024-panel D viscosity integral. PASS requires R1<120000 and
   nu*d_w<4, using exact endpoint comparisons.
3. Exact E local/source/Green/residual identities, with both radial tails.
4. Recomputed 2048-panel E range integrals. PASS requires E<1/60 and the
   sufficient compatibility bound; amplitude scaling handles either sign
   of a future E upper expression.
5. Independent exact-dyadic replay and common amplitude arithmetic.
6. Both actual certificates joined with the central jets and local-time constants.
7. The auxiliary initial-shear and frozen-model comparisons, explicitly separated from actual-solution persistence.

The interval engine is inherited unchanged from pass5/certificate and uses
mpmath.iv with rational cells and outward elementary-function arithmetic.
The result also depends on the documented analytic integration identities,
local-existence facts and earlier certified profile and retuning bounds.
Displayed intervals for upper-bound expressions are not two-sided intervals
for D, E, or the exact pressure.

Three independent audits inspect the D and E derivations and the range code;
an independent E run reproduced its upper endpoint. Symbolic checks verify
finite identities. No Lean or other formal proof replay was performed.

Exploratory finite-box pressure calculations in work/pass7 helped select the
certificate, but their values and convergence are not premises of the proved
inequalities. No numerical Navier–Stokes evolution was run.

The common existence interval follows qualitatively from strict uniform
initial margins and smooth dependence on a compact parameter family.
No useful numerical duration or return-map bound has been certified.

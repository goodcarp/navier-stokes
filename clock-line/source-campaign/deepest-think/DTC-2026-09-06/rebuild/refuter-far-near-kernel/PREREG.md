# PREREG — REFUTER seat vs GATE H-c (`far-near-kernel-lemma`)
Seat: `refuter-far-near-kernel`.  Written BEFORE any script in this folder was run.  2026-09-06.

Target: ~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/rebuild/far-near-kernel-lemma/
Default posture: REFUTED unless every item below survives.  Numerics falsify, never prove.

## Registered refutation attempts (each with the outcome that would fire it)

R1  PROVENANCE.  PREREG.md mtime must precede every script and every output mtime, and its content
    must contain the kills actually reported.  FIRES if PREREG postdates any run, or if a kill reported
    in the build result is absent from PREREG, or if a registered kill is silently dropped.

R2  REPRODUCTION.  Re-run s1..s6 from a COPY in this folder.  FIRES if any headline number in NOTE.md
    fails to reproduce to the precision it is printed at.

R3  INDEPENDENT DERIVATION.  I re-derive, by my own route and my own code (no import of their kern.py):
    (a) G_5 and K; (b) the interior constant in d_z[rho^l C_l^{3/2}] and the exterior constant;
    (c) alpha_l = g_l/(2l+3) from the single-layer jump; (d) N_l; (e) Phi(0) = -(3/5)g_1 and kappa_0=1/2 sharp;
    (f) C1 and C2in from the stated series; (g) the two collar pieces (bathtub radius R_A, and the pi/6 integral);
    (h) the multipole powers 4 and 5.  FIRES on any disagreement.

R4  INDEPENDENT INSTRUMENT.  An axisymmetric-ring quadrature of a(r,z) written by me from scratch
    (elliptic-form ring kernel, my own panels, no exclusion ball needed away from support) must reproduce
    the series at >=3 points to <1e-3 absolute, and must reproduce c_edge(10 deg) = 0.216773 to <2e-3.
    FIRES otherwise.

R5  FL-043 TWO-SIDEDNESS.  For each registered kill K1,K2,K3,K3',K4,K5,K5',K6 I ask: is there an
    outcome of the experiment as coded that would have fired it?  FIRES for any kill that is
    structurally incapable of firing (e.g. bound and "measurement" both computed from the same series,
    so the comparison only re-checks arithmetic), unless the build already says so.

R6  DECISION-RULE COMPLIANCE.  PREREG's decision rule is "PASSES iff K1,K2,K3,K3',K4,K5,K5' all fail to
    fire".  FIRES if a registered kill fired and the reported verdict is nonetheless PASS.

R7  CONSTANT CONVERGENCE (my own hypothesis, registered before running).  The all-l constant
    D = 124.3991 is an l1 sum  SUM_{l odd} |(l+1)g_l/(2l+3)| ((2^{l+4}-1)/(l+4)) ((l+2)(l+3)/2) 0.5^{l-1}
    truncated at LMAX=120.  The factor 2^{l+4} * 0.5^{l-1} = 32 is l-INDEPENDENT, so the summand behaves
    like ~8 l |g_l|; the bang-bang g_l are expected to decay only like l^{-3/2}, in which case the sum
    DIVERGES like sqrt(LMAX) and 124.40 is a truncation artifact.  FIRES if D grows with LMAX
    (criterion: D(LMAX=2000)/D(LMAX=120) > 1.5).  I will also compute the correct constants:
    (i) the same l1 sum at the stated validity range |x| >= 4 rho_in (factor 0.25^{l-1}, convergent), and
    (ii) a datum-INDEPENDENT Cauchy-Schwarz constant sqrt(2) INT_{2^-j}^{2^{1-j}} S(v) dv/v * 2^{5j}.

R8  STATEMENT AUDIT.  Every hypothesis of the assembled split and of "THE LEMMA" stated; each inequality
    checked for literal truth as written (not as intended).  FIRES on any inequality that is false as
    written for data satisfying its stated hypotheses.

## Verdict rule
refuted = TRUE if any of R1,R2,R3,R4,R6,R7,R8 fires on a load-bearing item, or if R5 finds an
undisclosed one-sided gate.  Severity: FATAL if the shell representation / modal expansions / kappa_0
or the multipole power is wrong; MAJOR if a reported constant or the gate verdict is wrong; MINOR if
only a statement is loose while the operative form is correct.

# PREREG — GATE H-c (far/near split of the 5D-lift strain kernel)
Seat: far-near-kernel-lemma.  Written BEFORE any script was run.  Register time: 2026-09-06.

## Objects
5D lift of axisymmetric no-swirl NS.  x = (y,z) in R^4 x R, r=|y|, rho=|x|, t = z/rho = cos(phi).
eta = omega^theta / r,  -(d_rr + (3/r) d_r + d_zz) psi1 = eta  (= -Laplacian_5 psi1 = eta),
a = u^r/r = -d_z psi1.  Green kernel G5 = 1/(8 pi^2 |x|^3), so
   a(x) = INT K(x-x') eta(x') dx',   K(w) = 3 w_z / (8 pi^2 |w|^5)   (homogeneous of degree -4).

## Claims to be gated (each with its kill)
C-A (shell representation).  With w := omega^theta and the shell operator
   Phi[w](xi) := INT_{S^4} K(xi - sigma) w(sigma)/s(sigma) dOmega_4(sigma),  s(sigma) = sin of polar angle,
   a(x) = INT_{rho0}^{R} Phi[w(rho' .)](x/rho') drho'/rho'.
   KILL K1: if this does not reproduce the validated instrument at 3 test points to < 1e-3 absolute.

C-B (interior modal identity).  Expanding w(sigma)/s(sigma) = SUM_l g_l C_l^{3/2}(t),
   Phi[w](xi) = - SUM_{l>=1} ((2l+1) g_l / (2l+3)) |xi|^{l-1} C_{l-1}^{3/2}(t)  for |xi| < 1.
   In particular the l=1 term is CONSTANT in xi: the exact uniform-strain fact, with
   Phi[w](0) = -(3/5) g_1 = -(3/4) INT_0^pi w cos(phi) sin^2(phi) dphi,  |Phi[w](0)| <= M/2 (sharp, bang-bang).
   KILL K2: if d_z[rho^l C_l^{3/2}(z/rho)] != (2l+1) rho^{l-1} C_{l-1}^{3/2}(z/rho) for any l <= 14 (sympy, exact).

C-C (exterior modal identity / multipole power).
   Phi[w](xi) = SUM_{l>=0} ((l+1) g_l / (2l+3)) |xi|^{-(l+4)} C_{l+1}^{3/2}(t)  for |xi| > 1.
   Leading power: |xi|^{-4} in general (l=0, g_0 propto shell mean of omega^theta);
   |xi|^{-5} when omega^theta is odd in z (g_0 = 0).
   KILL K3: if -d_z[rho^{-l-3} C_l^{3/2}] != (l+1) rho^{-l-4} C_{l+1}^{3/2} for any l <= 14 (sympy, exact).
   KILL K3': if the numerically measured decay exponent of an isolated z-odd octave's strain
   versus distance is not 5.00 +- 0.15 (and 4.00 +- 0.15 for a non-z-odd octave).

C-D (explicit far constant).  |a_far(x) - a_far(0)| <= C1 M with
   C1 = sqrt(2) INT_0^{1/2} Q(u) du/u,  Q(u)^2 = SUM_{l>=2} ((2l+1)/(2l+3))^2 (l(l+1)/2)^2 u^{2(l-1)} / N_l,
   N_l = (l+1)(l+2)/(l+3/2), using Parseval SUM_l g_l^2 N_l = INT_{-1}^1 w^2 dt <= 2 M^2.
   KILL K4: if a direct numerical evaluation of |a_far(x)-a_far(0)| at any test point EXCEEDS C1 M.

C-E (near part is O(M), NO log).  Split rho' < 2 rho into inner (rho' <= rho/2) and collar (rho/2 < rho' < 2 rho).
   Inner:  |a_inner| <= C2in M, C2in = sqrt(2) INT_0^{1/2} S(v) dv/v, S(v)^2 = SUM_{l>=0} ((l+1)/(2l+3))^2 ((l+2)(l+3)/2)^2 v^{2(l+4)}/N_l.
   Collar: |a_collar| <= (4/s + pi/8) M, s = sin(phi) at the evaluation point.
   KILL K5: if either bound is violated by a direct numerical evaluation at any test point.
   KILL K5': if the total near part shows any growth in log(rho/rho0) (fitted slope > 0.02 over 6 octaves).

C-F (the plateau constant).  For the bang-bang cap w = -M sgn(z) on rho0<rho<R, the exact
   inner-edge offset c_edge(phi) := lim_{R/rho0 -> inf} [ a(rho0+, phi) - (M/2) log(R/rho0) ]
   has the closed form  c_edge(phi) = - SUM_{l odd >= 3} ((2l+1) g_l /((2l+3)(l-1))) C_{l-1}^{3/2}(cos phi).
   PREDICTION: c_edge(10 deg) = +0.216773 M (the sharp-clock refuter's number).
   KILL K6: if the closed form differs from 0.216773 by more than 2e-3 -> the refuter's constant is
   CORRECTED and the corrected value is reported (this is a finding either way, not a failure of the gate).

## Decision rule
GATE H-c PASSES iff K1,K2,K3,K3',K4,K5,K5' all fail to fire.  K6 firing does not fail the gate;
it re-prices the plateau constant.  Numerics falsify; they never prove.  Every number printed comes
from a script in this folder.

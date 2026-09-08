#!/usr/bin/env python3
"""G4 — exponent arithmetic (exact Fractions) for the self-stretching consequence inside the estate's record ladder.
Inputs (estate ledger, exponents of M): ell = M^{-1/2}, delta = M^{-2/5}, L_lit = log(delta/ell) = (1/10) log M,
literal-reading interval T_lit = 2/M, staged interval h = M^{-19/20}; G2 constant band |c(phi)| <= 0.27 (numerics).
Self-stretching: a cap-level octave at scale rho inside a multi-octave field sustaining a(0) = (M/2) L(t) is stretched by
   exp( ∫ [ (M/2) L(t) - (M/2) log(rho/ell) - c M ] dt ).
Literal reading: innermost octave (rho ~ ell), L = L_lit constant on T_lit  ->  stretch = M^{1/10} * exp(-2c) >= M^{1/10} e^{-0.54}.
Record gap allowed inside I_j: M^{gamma}, gamma < 1/20 (telescope window).  Kill iff 1/10 > gamma  (always in the window).
Staged reading: only ∫ a(0) dt = (1/10) log M is demanded over h; average L = (1/5) M^{-1/20} log M << 1 : single-octave fields, the
mechanism has no multi-octave structure to bite on -> NO kill (recorded honestly).
Also: the protection regress — sub-structure at scales in (ell, rho) can supply at most (M/2) log(rho/ell) of counter-strain at rho
(extremal identity applied at the shell), so shells with rho < sqrt(R ell) cannot be protected: stretch >= (M/2) log(R ell / rho^2).
Control: the exponent arithmetic must give NO kill when the demanded charge is 2x smaller (deficit-2 sanity), and must give the
estate's gamma = 1/20 edge when T = 1/M.
"""
from fractions import Fraction as F
import math, sys, hashlib
ell, delta = F(-1,2), F(-2,5)
L_lit = delta - ell                    # coefficient of log M in log(delta/ell) = 1/10
assert L_lit == F(1,10)
c_band = 0.27
def stretch_exponent(L_coeff, T_exp, T_pref):  # a = (M/2) L log M sustained for T = T_pref * M^{T_exp}; returns exponent of M in the stretch
    return F(1,2)*L_coeff*T_pref*F(1)*(F(1) if T_exp == -1 else None) if False else F(1,2)*L_coeff*T_pref
# literal: T = 2/M  -> exponent (1/2)(1/10)(2) = 1/10 ; the -cM T term is e^{-2c}, a constant
e_lit = F(1,2)*L_lit*2
print('literal reading: innermost-octave stretch exponent =', e_lit, '| constant factor >= exp(-2*0.27) =', round(math.exp(-2*c_band),3))
gamma_max = F(1,20)
print('telescope window gamma <', gamma_max, '->', 'KILL (record broken mid-interval)' if e_lit > gamma_max else 'no kill')
# T = 1/M: exponent 1/20 = exactly the window edge
e_half = F(1,2)*L_lit*1
print('T = 1/M: stretch exponent =', e_half, '(= window edge gamma = 1/20:', e_half == gamma_max, ')')
# staged reading: average L over h = M^{-19/20}: (1/10) log M = (M/2) <L> h  ->  <L> = (1/5) M^{-1/20} log M -> exponent of M in <L> is -1/20 < 0
Lavg_exp = -(F(1) + F(-19,20))  # <L> = (1/5) log M / (M h) with h = M^{-19/20}: exponent of M is -(1 - 19/20) = -1/20
print('staged reading: <L> ~ M^{', Lavg_exp, '} log M  -> below one octave: mechanism does not bite (NO kill)')
# protection regress threshold: rho* = sqrt(R ell) with R = delta: exponent (delta + ell)/2
rho_star = (delta + ell)/2
print('protection regress: unprotectable shells rho < M^{', rho_star, '} (= sqrt(delta*ell)); innermost stretch coefficient at rho=ell: (1/2)log(R/ell) = (1/2)*(1/10) log M per unit M-time')
# controls
ctrl_deficit2 = F(1,2)*(L_lit/2)*2     # if the demanded charge were half: exponent 1/20 -> equals the edge, no strict kill
print('CTRL half-charge exponent =', ctrl_deficit2, '-> strict kill?', ctrl_deficit2 > gamma_max, '(must be False)')
ok = (e_lit == F(1,10)) and (e_lit > gamma_max) and (e_half == gamma_max) and (Lavg_exp == F(-1,20)) and (rho_star == F(-9,20)) and not (ctrl_deficit2 > gamma_max)
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

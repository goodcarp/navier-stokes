#!/usr/bin/env python3
"""G16 — GENERATION_1 Step 4: can Lei–Zhang's Gronwall be run with a time-dependent form constant δ*(t) that is large only on a
short set of times S (the blob's sustained-compression windows), closing on the time average?
Their integrated inequality (3.6)-(3.7) has the structure  Y' + (1 − K0 δ*(t)) D(t) ≤ g(t),  Y = 3C*²‖J‖²+‖Ω‖² ≥ 0, D = ‖∇Ω‖² ≥ 0, g ∈ L¹.
Claim to test: with δ*(t) ≤ δ0 (K0δ0 < 1/2) off S and δ*(t) = δ1 (K0δ1 > 1) on S, |S| small, the bound sup Y < ∞ still follows.
Exact counterexample (the inequality is an INEQUALITY; take it as an equality with the worst-case D): choose D(t) = D1 on S, 0 off S, g = 0:
  Y(T) = Y(0) + (K0δ1 − 1) D1 |S|  — unbounded in D1 for any |S| > 0. The inequality gives no control of ∫_S D, so averaging fails:
a coefficient that changes sign in front of the dissipation cannot be time-averaged.  Control: if K0δ1 < 1 (no sign change) then Y(T) ≤ Y(0) for every D1.
Consequence: Step 4 must localize in SPACE (the annulus (r_G, r0) carrying the failed 2D Hardy), not average in time. [OPEN, harder than stated]"""
import sympy as sp, sys, hashlib
K0,d0,d1,D1,S,Y0 = sp.symbols('K0 delta0 delta1 D1 S Y0', positive=True)
YT = Y0 + (K0*d1 - 1)*D1*S
grow = sp.limit(YT.subs({K0:1, d1:2, S:sp.Rational(1,100), Y0:1}), D1, sp.oo)
ctrl = sp.simplify(YT.subs({K0:1, d1:sp.Rational(1,2)}) - Y0)   # (1/2 - 1) D1 S = -D1 S/2 <= 0
print('worst-case Y(T) = Y0 + (K0 δ1 − 1) D1 |S| ->', grow, 'as D1 -> oo (sign change: averaging FAILS)')
print('CTRL no sign change (K0 δ1 = 1/2): Y(T) − Y0 =', ctrl, '(≤ 0 for all D1: bounded)')
ok = (grow == sp.oo) and sp.simplify(ctrl + D1*S/2) == 0
print('RESULT', 'PASS' if ok else 'FAIL', '| SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())
sys.exit(0 if ok else 1)

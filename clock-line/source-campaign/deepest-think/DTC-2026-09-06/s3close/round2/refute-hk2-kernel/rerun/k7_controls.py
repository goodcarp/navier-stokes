#!/usr/bin/env python3
"""k7 -- pre-declared falsification controls, fixed before the numbers were read.

R1  The bound must dominate the measurement at every point of N_tau.  If it does not, the
    theorem is refuted.                                                    [checked in k6]
R2  The kernel quadratures must be converged: relative change < 1e-6 between the last two
    refinements of (ns,nT,nc).                                             [k3_conv.json]
R3  The zonal series must satisfy Lap5 a = d_z eta to better than 1e-9 relative, and must
    agree with a direct 5D kernel quadrature of a = K*eta to better than 1e-6 relative.
    Two instruments sharing no code.                                       [k4_results.json]
R4  The ball sups must be grid-converged: doubling the sampling must move sup_B|grad eta| and
    sup_B||Hess eta||_F by less than 1%.                                   [here]
R5  K2 must NOT grow with L.  Doubling L from 10 to 40 (a 4x larger outer radius) must move
    the bound by less than 5%.  If it grows like log R, (H-K2) is false as stated.  [here]
R6  The (D-A) sharp-radial-edge bound must blow up like 1/f: K2hat * f must approach a
    constant as f -> 0.  If it does not, my claim that a sharp radial edge is fatal at f=0
    is wrong.                                                              [k6]
"""
import json, math
import numpy as np
from hk2lib import DatumB, StrainedB, geometry
from k3_bound import sups_on_ball_B

RES = {}
print("=== R4: grid convergence of the ball sups (DatumB, delta=7.5deg, w=0.20, L=10) ===")
DB = DatumB(delta_deg=7.5, w=0.20, L=10.0)
rows = []
for lam in (1.0, 1.25, 1.5):
    g = geometry(30.0, 0.0, 7.5, lam); DL = StrainedB(DB, lam)
    prev = None
    for n in (201, 401, 801):
        s1, s2, sb = sups_on_ball_B(DL, g['r'], g['z'], 0.20, n=n)
        if prev is not None:
            d1 = abs(s1/prev[0]-1); d2 = abs(s2/prev[1]-1)
            print(f"  lam={lam:4.2f} n={n:4d}: sup|grad|={s1:.6f} ({d1:.2e})  "
                  f"sup||Hess||_F={s2:.6f} ({d2:.2e})")
            rows.append(dict(lam=lam, n=n, sup_grad=s1, sup_hess=s2, rel_grad=d1, rel_hess=d2))
        prev = (s1, s2)
RES['R4'] = rows
worst = max(max(r['rel_grad'], r['rel_hess']) for r in rows)
RES['R4_worst_rel'] = worst
print(f"  worst relative move on refinement: {worst:.3e}   (rule: < 1e-2)  ->"
      f" {'NOT VIOLATED' if worst < 1e-2 else 'VIOLATED'}")

print("\n=== R5: L-independence of the bound (k3 rows at L = 10 vs L = 40) ===")
R3 = json.load(open('k3_results.json'))
cmp = []
for lam in (1.0, 1.25, 1.5):
    a = [x for x in R3['DB'] if x['L']==10.0 and abs(x['lam']-lam)<1e-12][0]
    b = [x for x in R3['DB'] if x['L']==40.0 and abs(x['lam']-lam)<1e-12][0]
    cmp.append(dict(lam=lam, K2_L10=a['K2'], K2_L40=b['K2'], rel=abs(b['K2']/a['K2']-1)))
    print(f"  lam={lam:4.2f}:  K2hat(L=10) = {a['K2']:.4f}   K2hat(L=40) = {b['K2']:.4f}"
          f"   relative change = {abs(b['K2']/a['K2']-1):.3e}")
RES['R5'] = cmp
wr = max(c['rel'] for c in cmp); RES['R5_worst_rel'] = wr
print(f"  worst: {wr:.3e}  (rule: < 5e-2; log R would give log40/log10 - 1 = "
      f"{math.log(math.exp(40))/math.log(math.exp(10))-1:.1f})  -> "
      f"{'NOT VIOLATED' if wr < 5e-2 else 'VIOLATED'}")
json.dump(RES, open('k7_results.json','w'), indent=1)
print("\nWROTE k7_results.json")

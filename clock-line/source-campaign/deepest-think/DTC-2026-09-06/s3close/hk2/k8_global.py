#!/usr/bin/env python3
"""k8 -- the GLOBAL K2 question (the refuter's Finding 4).  Theorem V.4's coupling lemma uses
K2 pathwise, i.e. globally, not only on the tube N_tau.  Here I evaluate both the bound and the
independent measurement at the three hard places of the campaign datum (D-B): inside the
equatorial layer, close to the axis taper, and at the inner radial edge."""
import json, math
import numpy as np
from hk2lib import DatumB
from k3_bound import lambda_bulk, assemble, sups_on_ball_B, C_K, C_gradK
from k4_direct import Series, K2_exact

DB = DatumB(delta_deg=7.5, w=0.20, L=10.0)
S  = Series(DB, LMAX=161)
PTS = [("tracked point (phi0=30 deg, rho=rho0)", 0.50000, 0.86603),
       ("inside the equatorial layer (z = 0.05 rho0)", 0.90, 0.05),
       ("equatorial plane exactly (z = 0)", 0.90, 0.0),
       ("near the axis taper (phi = 3 deg)", 0.0523, 1.00),
       ("inner radial edge on the mid-latitude ray", 0.7071, 0.7071)]
rows = []
for lbl, r, z in PTS:
    # measurement
    A = S.a_field(r, z)
    e = [v[0] for v in DB.eta_rz(np.array([r]), np.array([z]))]
    k2m, _ = K2_exact(r, A[1],A[2],A[3],A[4],A[5], e[0], e[1], e[2], ne=161)
    lap = A[3] + 3*A[1]/r + A[5]
    pde = abs(lap - e[2])/max(abs(e[2]), 1e-30)
    # bound, d optimised over a small grid
    best = None
    for frac in (0.05, 0.10, 0.2, 0.3, 0.45, 0.6, 0.8):
        d = frac*min(r, 0.5)
        Ls = lambda_bulk(DB.grad_norm_fast, (r,z), d, 3*DB.R, 300,120,120)
        sg, sh, sb = sups_on_ball_B(DB, r, z, d, n=201)
        B = assemble(sg, sh, sb, d, Ls[4], Ls[5], r, e[0], e[1], e[2])
        B.update(d=d)
        if best is None or B['K2'] < best['K2']: best = B
    rows.append(dict(label=lbl, r=r, z=z, K2_measured=k2m, K2_bound=best['K2'], d=best['d'],
                     grad_a_b=best['grad_a'], hess_a_b=best['hess_a'], pde_residual=pde))
    print(f"  {lbl:44s}  K2hat measured = {k2m:9.4f}   bound = {best['K2']:10.3f}"
          f"  (d*={best['d']:.2f})   [Lap5 a = d_z eta to {pde:.1e}]")
out = dict(points=rows,
           K2_measured_max=max(x['K2_measured'] for x in rows),
           K2_bound_max=max(x['K2_bound'] for x in rows))
print(f"\n  GLOBAL (over these probes): measured K2hat = {out['K2_measured_max']:.4f} ;"
      f"  bound = {out['K2_bound_max']:.3f}")
json.dump(out, open('k8_results.json','w'), indent=1)
print("WROTE k8_results.json")

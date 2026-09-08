#!/usr/bin/env python3
"""r0 -- calibrate my instrument against hk2's own centre numbers (the object under test is
imported by VALUE from its results JSON, read-only; no hk2 code is executed)."""
import json, math, time, numpy as np
from rlib import *
hk2 = json.load(open('../../hk2/k3_results.json'))
rows = [x for x in hk2['DB'] if x['L'] == 10.0]
out = {}
t0 = time.time()
for x in rows:
    lam = x['lam']; F = Field(7.5, 0.20, 10.0, lam)
    r, z = strain_point(30.0, lam)
    # evaluate at hk2's own optimised d, so the comparison is like for like
    b = bound_at(F, r, z, x['d'], dgrid=[x['d']])
    rel = abs(b['K2']/x['K2'] - 1)
    out[str(lam)] = dict(mine=b, hk2=dict(K2=x['K2'], d=x['d'], L4=x['L4'], L5=x['L5'],
                          sup_grad=x['sup_grad'], sup_hess=x['sup_hess']), rel=rel,
                          point_rel=abs(r/x['r']-1)+abs(z/x['z']-1))
    print(f"lam={lam:4.2f} d={x['d']:.4f}: mine K2hat={b['K2']:.4f} hk2={x['K2']:.4f} rel={rel:.2e} | "
          f"L4 {b['L4']:.5f}/{x['L4']:.5f}  L5 {b['L5']:.5f}/{x['L5']:.5f}  supgrad {b['sup_grad']:.5f}/{x['sup_grad']:.5f}"
          f"  suphess {b['sup_hess']:.4f}/{x['sup_hess']:.4f}")
out['worst_rel'] = max(v['rel'] for k,v in out.items() if k != 'worst_rel')
out['seconds'] = time.time()-t0
json.dump(out, open('r0_results.json','w'), indent=1)
print("worst rel", out['worst_rel'], "sec", out['seconds'])

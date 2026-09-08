#!/usr/bin/env python3
"""r1 -- (A) the lambda scan: is the sup over lambda in [1,3/2] at lambda = 3/2, and is it
attained by hk2's three samples?  (B) the tube: hk2 evaluates the majorant only at the
trajectory CENTRE X(lambda); N_tau := U B(X(t), d) is a tube, and the claim is 'on the whole
tube'.  Probe the majorant at points X + rho_t * e for tube radii rho_t in {sqrt(nu tau) at
L=10, 0.10, d*} and four directions (toward the equator, the taper cone, the axis, the origin),
each with its own admissible, optimised ball radius."""
import json, math, time, numpy as np
from rlib import *
Q = dict(ns=200, nT=100, nc=100)
out = {'quad': Q}
t0 = time.time()
# ---------------- (A) lambda scan
print("=== (A) lambda scan, centre of the tube, d optimised ===")
scan = []
for lam in np.linspace(1.0, 1.5, 11):
    F = Field(7.5, 0.20, 10.0, lam); r, z = strain_point(30.0, lam)
    D = distances(r, z, lam)
    dmax = 0.95*min(D['d_eq'], D['d_tap'], 0.95*r)
    b = bound_at(F, r, z, dmax, nd=8, **Q)
    b.update(lam=float(lam), **D); scan.append(b)
    print(f"  lam={lam:5.3f}  d*={b['d']:.4f}  K2hat<={b['K2']:.4f}   (d_eq={D['d_eq']:.4f} d_tap={D['d_tap']:.4f})")
out['scan'] = scan
K = [s['K2'] for s in scan]
out['scan_monotone'] = bool(all(K[i] <= K[i+1] for i in range(len(K)-1)))
out['scan_sup'] = max(K); out['scan_argsup'] = scan[int(np.argmax(K))]['lam']
print(f"  monotone: {out['scan_monotone']}   sup = {out['scan_sup']:.4f} at lam = {out['scan_argsup']}")
# ---------------- (B) tube probes
print("\n=== (B) tube probes: majorant at X(lam) + rho_t * e ===")
nu_tau_L10 = math.sin(math.radians(7.5))**2*THETA_MAX/10.0
out['sqrt_nu_tau_L10'] = math.sqrt(nu_tau_L10)
probes = []
for lam in (1.0, 1.25, 1.5):
    F = Field(7.5, 0.20, 10.0, lam); r0, z0 = strain_point(30.0, lam)
    c = [s for s in scan if abs(s['lam']-lam) < 1e-9][0]
    dstar = c['d']; ang = math.radians(c['cone_deg'])
    dirs = {'equator': (0.0, -1.0), 'taper': (-math.cos(ang), math.sin(ang)),
            'axis': (-1.0, 0.0), 'origin': (-r0/math.hypot(r0,z0), -z0/math.hypot(r0,z0))}
    for rt_lbl, rt in (('sqrt(nu tau)@L10', out['sqrt_nu_tau_L10']), ('0.10', 0.10), ('d*', dstar)):
        for dl, (ur, uz) in dirs.items():
            r, z = r0 + rt*ur, z0 + rt*uz
            D = distances(r, z, lam)
            dmax = 0.95*min(D['d_eq'], D['d_tap'], 0.95*r)
            if dmax <= 0.01:
                row = dict(lam=lam, tube_radius=rt, tube_label=rt_lbl, direction=dl, r=r, z=z,
                           K2=float('inf'), note='ball radius collapses', **D)
            else:
                b = bound_at(F, r, z, dmax, nd=6, **Q)
                b.update(lam=lam, tube_radius=rt, tube_label=rt_lbl, direction=dl, **D)
                row = b
            row['ratio_to_centre'] = row['K2']/c['K2']
            probes.append(row)
            print(f"  lam={lam:4.2f} rho_t={rt:.4f} ({rt_lbl:>16}) {dl:>8}: K2hat<={row['K2']:10.4f}"
                  f"  x{row['ratio_to_centre']:.3f} centre   d_eq={D['d_eq']:.4f} d_tap={D['d_tap']:.4f} d*={row.get('d',0):.4f}")
out['probes'] = probes
for lbl in ('sqrt(nu tau)@L10', '0.10', 'd*'):
    rows = [p for p in probes if p['tube_label'] == lbl]
    out[f'tube_sup_{lbl}'] = max(p['K2'] for p in rows)
    out[f'tube_worst_ratio_{lbl}'] = max(p['ratio_to_centre'] for p in rows)
    print(f"  tube radius {lbl:>16}: sup over probes = {out[f'tube_sup_{lbl}']:.4f}, worst ratio to centre = {out[f'tube_worst_ratio_{lbl}']:.3f}")
out['centre_sup'] = max(c['K2'] for c in scan if abs(c['lam']-1.5)<1e-9)
out['seconds'] = time.time()-t0
json.dump(out, open('r1_results.json','w'), indent=1)
print("WROTE r1_results.json", out['seconds'])

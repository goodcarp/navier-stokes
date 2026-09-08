#!/usr/bin/env python3
"""r2_test -- controls on the refuter's instrument: (i) bulk derivatives vs finite differences,
(ii) a, grad a, Hess a for (D-B) at the tracked point vs hk2's k4/k9 numbers (read from their
results files ONLY for comparison, never as input), (iii) Hessian symmetry and the transverse
control mu = a_r/r, (iv) the PDE control Lap5 a = d_z eta."""
import json, math, time
import numpy as np
from r2lib import *

out = {}
print("=== (i) bulk derivatives vs 4th-order central differences ===")
def fd_check(D, pts, h=1e-3):
    worst = 0.0
    f = lambda R, Z: D.bulk(np.array([R]), np.array([Z]))[0][0]
    d1 = lambda g: (-g(2*h)+8*g(h)-8*g(-h)+g(-2*h))/(12*h)
    d2 = lambda g: (-g(2*h)+16*g(h)-30*g(0.0)+16*g(-h)-g(-2*h))/(12*h*h)
    for (r, z) in pts:
        e, er, ez, err, erz, ezz = [v[0] for v in D.bulk(np.array([r]), np.array([z]))]
        num = [d1(lambda u: f(r+u, z)), d1(lambda u: f(r, z+u)), d2(lambda u: f(r+u, z)),
               d1(lambda u: d1(lambda v: f(r+u, z+v))), d2(lambda u: f(r, z+u))]
        ana = [er, ez, err, erz, ezz]
        rel = max(abs(a_-n_)/max(abs(n_), 1e-2) for a_, n_ in zip(ana, num))
        worst = max(worst, rel)
    return worst
R = math.exp(10.0)
DB = Datum(AngularDB(7.5, 0.20), RadialTanh(R))
DA = Datum(AngularSharp(7.5), RadialSharp(R))
DH = Datum(AngularSharp(7.5), RadialTanh(R))
DBs = Datum(AngularDB(7.5, 0.20), RadialTanh(R), lam=1.5)
pts = [(0.5, 0.866), (0.75, 0.3849), (0.9, 0.2), (0.3, 1.2), (1.5, 2.0), (1.2, -0.7)]
pts2 = [(0.6, 1.0), (0.9, 0.9), (0.3, 1.2), (1.5, 2.0), (1.2, -0.7), (2.0, 0.3)]   # off the sharp sphere
for name, D in (('DB', DB), ('DA', DA), ('hybrid', DH), ('DB_lam1.5', DBs)):
    w = fd_check(D, pts if name.startswith('DB') else pts2); out[f'fd_{name}'] = w; print(f"   {name:10s} worst rel {w:.2e}")
    assert w < 1e-6
# strain identity: eta_lam(T_lam x) = eta_0(x)
r0, z0 = 0.5, 0.866; lam = 1.5
v0 = DB.bulk(np.array([r0]), np.array([z0]))[0][0]; v1 = DBs.bulk(np.array([lam*r0]), np.array([z0/lam**2]))[0][0]
print(f"   material identity eta_lam(T x0) = eta_0(x0): {v0:.12f} vs {v1:.12f}"); assert abs(v0-v1) < 1e-12

print("\n=== (ii) (D-B) at the tracked point (0.5, 0.8660), lam = 1, L = 10 ===")
r, z = 0.5, math.sqrt(3)/2
t0 = time.time()
I = Instrument(DB, r, z, nT=160, nX=120)
m = I.measure(0.25)
print(f"   time {time.time()-t0:.1f}s")
for k in ('a', 'a_r', 'a_z', 'a_rr', 'a_rz', 'a_zz', 'hess_asym', 'mu_integrated', 'mu_from_a_r', 'hess_op', 'grad_a'):
    print(f"   {k:14s} {m[k]:+.9f}")
out['DB_tracked_d0.25'] = {k: float(v) for k, v in m.items()}
# d-independence control: the split radius must not matter
m2 = I.measure(0.12)
print(f"   d-independence (d=0.12 vs 0.25): a_rr {m2['a_rr']:+.9f}  a_rz {m2['a_rz']:+.9f}  a_zz {m2['a_zz']:+.9f}")
out['DB_tracked_d0.12'] = {k: float(v) for k, v in m2.items()}
out['d_independence_rel'] = max(abs(m2[k]-m[k])/abs(m[k]) for k in ('a_rr', 'a_rz', 'a_zz'))
print(f"   worst relative move {out['d_independence_rel']:.2e}")
# PDE control: Lap5 a = a_rr + 3 a_r/r + a_zz = d_z eta
e = DB.bulk(np.array([r]), np.array([z]))
lap = m['a_rr'] + 3*m['a_r']/r + m['a_zz']
print(f"   Lap5 a = {lap:+.9f}   d_z eta = {e[2][0]:+.9f}   rel {abs(lap/e[2][0]-1):.2e}")
out['pde_control_rel'] = float(abs(lap/e[2][0]-1))
# comparison with hk2's numbers (their files read only here, for the comparison)
try:
    k4 = json.load(open('../../hk2/k4_results.json')); k9 = json.load(open('../../hk2/k9_results.json'))
    out['hk2_k9_control_lam1'] = k9['control_lam1']; print("   hk2 k9 control_lam1:", k9['control_lam1'])
    print("   hk2 k9 measured lam=1:", k9['measured_strained'][0])
except Exception as ex:
    print("   (hk2 results not read)", ex)
json.dump(out, open('r2_test_results.json', 'w'), indent=1)
print("WROTE r2_test_results.json")

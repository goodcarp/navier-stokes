#!/usr/bin/env python3
"""Validations for the sharp/viscous-numerics solver. Controls C2, C3, C4 of PREREG.md.
Every control must FIRE (i.e. the deliberately wrong variant must be visibly wrong)."""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from nsring import *

out = {"script_sha256": sha_self()}
print("=" * 78); print("C3 / C2 : exact Hill spherical vortex, full domain, Dirichlet = exact monopole")
print(f"{'h':>7} {'rel err psi1':>14} {'a_nose':>10} {'axis peak/U':>12} {'bad-op rel err':>15}")
rows = []
for h in [0.08, 0.04, 0.02]:
    g = Grid(2.5, -2.5, 2.5, h)
    eta = hill_eta(g, 1.0, 1.0, sub=4, smooth=0.0)
    ex, U = hill_psi1_exact(g)
    RR, ZZ = g.node_mesh(); S = np.sqrt(RR**2 + ZZ**2)
    far = 0.25 * (4.0 / 15.0) / np.maximum(S, 1e-12)**3
    en = cell_to_node(g, eta, "none")
    for bad in (False, True):
        P = Poisson(g, bad=bad)
        psi = P.solve(en, g_bot=far[:, 0], g_top=far[:, -1], g_right=far[-1, :])
        rel = float(np.linalg.norm(psi - ex) / np.linalg.norm(ex))
        if not bad:
            good = rel
            ur, uz, a = node_velocity(g, psi)
            # nose = (r=0, z=rho); a there
            jn = int(round((1.0 - g.z0) / h))
            a_nose = float(a[0, jn])
            peak = float(np.max(uz[0, :]) / U)
        else:
            badrel = rel
    print(f"{h:7.3f} {good:14.3e} {a_nose:10.5f} {peak:12.5f} {badrel:15.3e}")
    rows.append(dict(h=h, rel=good, a_nose=a_nose, peak_over_U=peak, bad_rel=badrel))
out["hill"] = rows
r2, r1 = rows[-1], rows[-2]
rich_a = r2["a_nose"] + (r2["a_nose"] - r1["a_nose"]) / 3.0
print(f"  Richardson (2nd order) a_nose -> {rich_a:.5f}   exact 1/5 = 0.20000")
print(f"  exact axis peak/U = 2.5 ; measured {rows[-1]['peak_over_U']:.5f}")
out["hill_richardson_a_nose"] = rich_a
C3 = rows[-1]["rel"] < 2e-3 and abs(rich_a - 0.2) < 5e-3 and abs(rows[-1]["peak_over_U"] - 2.5) < 5e-3
C2 = all(r["bad_rel"] > 1.0 for r in rows)
print(f"  C3 exact-Hill validation: {'PASS' if C3 else 'FAIL'}")
print(f"  C2 operator control (3/r -> 1/r must exceed 100% error): {'FIRES' if C2 else 'DID NOT FIRE'}")

print("=" * 78); print("C4 : pure 5D diffusion, u = 0, must give d<|X|^2>/dt = 10 nu")
c4 = []
for nu in [0.02, 0.05]:
    h = 0.04; g = Grid(3.0, -3.0, 3.0, h)
    eta = hill_eta(g, 1.0, 1.0, sub=4, smooth=1.5 * h)
    w = eta * g.Rc[:, None]**3
    m2 = lambda e: float(np.sum(e * g.Rc[:, None]**3 * (g.Rc[:, None]**2 + g.Zc[None, :]**2))
                         / np.sum(e * g.Rc[:, None]**3))
    dt = 0.05 * h**2 / nu; T = 0.5; n = int(T / dt)
    m0 = m2(eta)
    for _ in range(n):
        k1 = lap5(g, eta, "none"); e1 = eta + nu * dt * k1
        eta = eta + 0.5 * nu * dt * (k1 + lap5(g, e1, "none"))
    rate = (m2(eta) - m0) / (n * dt)
    print(f"  nu = {nu:.3f}   d<|X|^2>/dt = {rate:.5f}   exact 10 nu = {10*nu:.5f}"
          f"   rel {abs(rate-10*nu)/(10*nu):.2e}")
    c4.append(dict(nu=nu, rate=rate, exact=10 * nu))
out["diffusion5d"] = c4
C4 = all(abs(c["rate"] - c["exact"]) / c["exact"] < 5e-3 for c in c4)
print(f"  C4 5D diffusion law: {'PASS' if C4 else 'FAIL'}")

print("=" * 78); print("odd-symmetry half-domain solver == full-domain solver on an odd datum")
h = 0.05
gf = Grid(2.0, -2.0, 2.0, h); gh = Grid(2.0, 0.0, 2.0, h)
ef, _ = datum_shell(gf, 1); eh, _ = datum_shell(gh, 1)
Pf = Poisson(gf); Ph = Poisson(gh)
zf = np.zeros(gf.Nr + 1); zj = np.zeros(gf.Nz + 1)
pf = Pf.solve(cell_to_node(gf, ef, "none"), g_bot=zf, g_top=zf, g_right=zj)
ph = Ph.solve(cell_to_node(gh, eh, "odd"), g_bot=np.zeros(gh.Nr + 1),
              g_top=np.zeros(gh.Nr + 1), g_right=np.zeros(gh.Nz + 1))
j0 = int(round(2.0 / h))
d = float(np.max(np.abs(pf[:, j0:] - ph)) / np.max(np.abs(ph)))
print(f"  max rel difference psi1 (half vs full) = {d:.3e}")
out["halfdomain_vs_full"] = d
C0 = d < 1e-12
print(f"  half-domain symmetry check: {'PASS' if C0 else 'FAIL'}")

out["verdict"] = dict(C3=bool(C3), C2=bool(C2), C4=bool(C4), halfdomain=bool(C0))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "validate.json"), "w"),
          indent=1, default=float)
print("=" * 78)
print("SCRIPT-SHA256", out["script_sha256"])
sys.exit(0 if (C3 and C2 and C4 and C0) else 2)

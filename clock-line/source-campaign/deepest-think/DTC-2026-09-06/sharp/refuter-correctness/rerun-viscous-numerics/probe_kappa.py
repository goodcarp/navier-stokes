#!/usr/bin/env python3
"""t=0 probe: measure kappa = d a / d log(R/rho0) for datum A and B, and the sign of a at
the argmax of |omega^theta|. Control C6 (single octave shows no log gain) lives here too."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from nsring import *

def setup(N, h, lam=1.5, kind="shell", rho0=1.0):
    R = rho0 * 2.0**(N if kind == "shell" else N - 1)
    L = max(lam * R, 4.0 * rho0)
    g = Grid(L, 0.0, L, h)
    eta, R = (datum_shell(g, N, rho0) if kind == "shell" else datum_rings(g, N, rho0))
    om = eta * g.Rc[:, None]
    A = 1.0 / np.max(np.abs(om))
    eta *= A
    return g, eta, R

print(f"{'kind':>6} {'N':>3} {'R':>7} {'a(0,0)/M':>10} {'a@argmax':>10} {'s*/rho0':>8} "
      f"{'log(R/r0)':>10} {'grid':>12}")
res = []
for kind in ("shell", "rings"):
    for N in [1, 2, 3, 4, 5]:
        h = 1.0 / 8.0
        g, eta, R = setup(N, h, kind=kind)
        P = Poisson(g)
        zr = np.zeros(g.Nr + 1); zj = np.zeros(g.Nz + 1)
        psi = P.solve(cell_to_node(g, eta, "odd"), g_bot=zr, g_top=zr, g_right=zj)
        ur, uz, a = node_velocity(g, psi)
        om = eta * g.Rc[:, None]
        i, j = np.unravel_index(np.argmax(np.abs(om)), om.shape)
        rs, zs = g.Rc[i], g.Zc[j]
        a_at = interp2(g, a, rs, zs)
        a00 = float(a[0, 0])
        lr = math.log(R / 1.0)
        print(f"{kind:>6} {N:3d} {R:7.1f} {a00:10.5f} {a_at:10.5f} {math.hypot(rs,zs):8.3f} "
              f"{lr:10.4f} {g.Nr}x{g.Nz}")
        res.append(dict(kind=kind, N=N, R=R, a00=a00, a_at=a_at, s=math.hypot(rs, zs), logR=lr,
                        rmax_om=float(np.max(np.abs(om)))))
    sub = [x for x in res if x["kind"] == kind]
    x = np.array([s["logR"] for s in sub]); y = np.array([s["a00"] for s in sub])
    k0 = np.polyfit(x, y, 1)
    y2 = np.array([s["a_at"] for s in sub]); k1 = np.polyfit(x, y2, 1)
    print(f"   -> kappa (axis)    d a(0,0)/d log(R/rho0) = {k0[0]:.5f}  intercept {k0[1]:.5f}")
    print(f"   -> kappa (argmax)  d a*/d log(R/rho0)     = {k1[0]:.5f}  intercept {k1[1]:.5f}")
    res.append(dict(kind=kind, kappa_axis=float(k0[0]), kappa_argmax=float(k1[0]),
                    icept_axis=float(k0[1]), icept_argmax=float(k1[1])))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_kappa.json"), "w"),
          indent=1, default=float)
print("SCRIPT-SHA256", hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())

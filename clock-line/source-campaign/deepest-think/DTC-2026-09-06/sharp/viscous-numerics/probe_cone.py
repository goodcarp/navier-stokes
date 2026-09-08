#!/usr/bin/env python3
"""t=0 strain profile a(s) along the 45-degree cone, where |omega^theta| = M.
Claim under test: a(s) = kappa M log(R/s) + O(M), i.e. the innermost ring sees the full log."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from nsring import *
c = 1.0 / math.sqrt(2.0)
res = {}
for kind in ("shell", "rings"):
    for N in [3, 5, 6]:
        R = 2.0**(N if kind == "shell" else N - 1)
        h = 1.0 / 8.0
        L = max(1.5 * R, 4.0)
        g = Grid(L, 0.0, L, h)
        eta, R = (datum_shell(g, N) if kind == "shell" else datum_rings(g, N))
        om = eta * g.Rc[:, None]; eta /= np.max(np.abs(om))
        P = Poisson(g); zr = np.zeros(g.Nr + 1); zj = np.zeros(g.Nz + 1)
        psi = P.solve(cell_to_node(g, eta, "odd"), g_bot=zr, g_top=zr, g_right=zj)
        _, _, a = node_velocity(g, psi)
        om = eta * g.Rc[:, None]
        ss = [1.0, 1.5, 2.0, 3.0, 4.0, 8.0, 16.0, 32.0]
        ss = [s for s in ss if s <= R]
        line = []
        for s in ss:
            av = interp2(g, a, s * c, s * c)
            ov = abs(interp2(g, np.pad(om, ((0, 1), (0, 1)), mode='edge'), s * c, s * c))
            line.append((s, av, ov))
        print(f"[{kind} N={N} R={R:.0f}] a(0,0)={a[0,0]:.4f}")
        print("     s : " + "  ".join(f"{s:7.1f}" for s, _, _ in line))
        print("  a(s) : " + "  ".join(f"{v:7.4f}" for _, v, _ in line))
        print("  |om| : " + "  ".join(f"{v:7.4f}" for _, _, v in line))
        x = np.array([math.log(R / s) for s, _, _ in line if s < R * 0.6])
        y = np.array([v for s, v, _ in line if s < R * 0.6])
        if len(x) > 2:
            p = np.polyfit(x, y, 1)
            print(f"  fit a = {p[0]:.4f} log(R/s) + {p[1]:.4f}")
            res[f"{kind}_N{N}"] = dict(slope=float(p[0]), icept=float(p[1]),
                                       a00=float(a[0, 0]), R=R,
                                       line=[[float(q) for q in t] for t in line])
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_cone.json"), "w"), indent=1)
print("SCRIPT-SHA256", hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())

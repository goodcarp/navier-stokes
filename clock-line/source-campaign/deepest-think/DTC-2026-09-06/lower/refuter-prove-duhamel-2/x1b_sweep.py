#!/usr/bin/env python3
"""x1b -- how much (H2) violation does prove-duhamel's TEST-1 gate tolerate?
Integrate a' = -a^2 + B(t) (a0=1) with B = bp on [0,t1], bn afterwards, and apply the
seat's own gate: HOLDS iff min_{t>0} F_a >= 1-1e-2 and min_{t>0} F_r >= 1-1e-2,
F_a = a(1+t), F_r = exp(int a)/(1+t)."""
import numpy as np, math, json
from scipy.integrate import solve_ivp

def run(bp, t1, bn, tend=0.55):
    def B(t): return bp if t < t1 else bn
    def f(t, y): return [-y[0]**2 + B(t), y[0]]
    s = solve_ivp(f, (0, tend), [1.0, 0.0], rtol=1e-11, atol=1e-13, max_step=1e-3, dense_output=True)
    ts = np.linspace(0, tend, 5001)[1:]
    a, I = s.sol(ts)
    Fa = a*(1+ts); Fr = np.exp(I)/(1+ts)
    return float(np.min(Fa)), float(np.min(Fr))

rows=[]; hold=[]
for bp in (1.5, 2.0, 3.0, 4.0):
    for t1 in (0.25, 0.35, 0.45, 0.55):
        for bn in (-0.25,-0.5,-0.75,-1.0,-1.5,-2.0,-3.0):
            if t1 >= 0.55 and bn < 0: continue
            mFa, mFr = run(bp, t1, bn)
            g = "HOLDS" if (mFa>=0.99 and mFr>=0.99) else "VIOLATED"
            rows.append(dict(bp=bp,t1=t1,bn=bn,minFa=mFa,minFr=mFr,gate=g))
            if g=="HOLDS": hold.append(rows[-1])
print(f"{'bp':>5} {'t1':>5} {'bn':>6} {'minFa':>10} {'minFr':>10} {'gate':>9}")
for r in rows:
    if r['gate']=="HOLDS":
        print(f"{r['bp']:>5.1f} {r['t1']:>5.2f} {r['bn']:>6.2f} {r['minFa']:>10.6f} {r['minFr']:>10.6f} {r['gate']:>9}")
print()
if hold:
    w = min(hold, key=lambda r: r['bn'])
    print("MOST NEGATIVE beta that still passes the gate:", w)
    print(f"  -> the gate certifies HOLDS while (beta+nuV)/a0^2 = {w['bn']:.2f}, i.e. (H2) is violated")
    print(f"     by {abs(w['bn']):.2f} a0^2 over {0.55-w['t1']:.2f}/0.55 of the window.")
json.dump(rows, open("x1b_results.json","w"), indent=1)

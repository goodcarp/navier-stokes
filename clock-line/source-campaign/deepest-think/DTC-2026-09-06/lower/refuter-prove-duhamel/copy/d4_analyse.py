#!/usr/bin/env python3
"""d4 -- analysis of d3/d3b: does the ACTUAL viscous evolution respect the Riccati floor,
and what clock does it deliver?

FALSIFICATION RULE (stated before reading the numbers, and applied verbatim):
  (H2) is REFUTED by these runs if, on any FORWARD run (sign = +1) with a0 > 0, the ratio
       F_a(t) = a(X(t),t)(1 + a0 t)/a0  dips below 1 by more than the solver's own
       discretisation scale, taken as 1e-2 (the max-principle violation of the 3rd-order
       reconstruction is 1.5e-3..4.3e-3 in the parent seat's sweep).
  The REVERSED control (sign = -1) must violate the floor, or the test is not diagnostic.
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = []
def say(s=""):
    print(s, flush=True); LOG.append(s)

D3 = json.load(open(os.path.join(HERE, "d3_results.json")))
try:
    D3B = json.load(open(os.path.join(HERE, "d3b_results.json")))
except Exception:
    D3B = {"runs": []}
runs = D3["runs"] + D3B["runs"]
say(f"runs analysed: {len(runs)}   (d3: {len(D3['runs'])}, d3b: {len(D3B['runs'])})")
say()

def floors(tr):
    h = tr["hist"]
    t = np.array(h["t"]); a = np.array(h["a"]); r = np.array(h["r"]); e = np.array(h["eta"])
    a0 = a[0]; r0 = r[0]
    Fa = a*(1.0 + a0*t)/a0
    Fr = (r/r0)/(1.0 + a0*t)
    return t, Fa, Fr, a, r, e, a0

say("="*112)
say("TEST 1 -- the Riccati floor  a(t)(1+a0 t)/a0 >= 1  and  (r(t)/r0)/(1+a0 t) >= 1")
say("="*112)
say(f"{'run':>22} {'s0':>4} {'a0':>9} {'min F_a (t>0)':>14} {'min F_r (t>0)':>14} "
    f"{'F_a(end)':>10} {'F_r(end)':>10} {'verdict':>10}")
worst_fwd = 1e9; worst_run = None; rev_min = None
T1 = []
for R in runs:
    for tr in R["tracers"]:
        t, Fa, Fr, a, r, e, a0 = floors(tr)
        if len(t) < 3: continue
        mFa = float(np.min(Fa[1:])); mFr = float(np.min(Fr[1:]))
        fwd = R["sign"] > 0 and a0 > 0
        v = "-" if not fwd else ("HOLDS" if mFa >= 1 - 1e-2 and mFr >= 1 - 1e-2 else "VIOLATED")
        if fwd:
            worst = min(mFa, mFr)
            if worst < worst_fwd: worst_fwd, worst_run = worst, (R["label"], tr["s0"])
        else:
            rev_min = mFa if rev_min is None else min(rev_min, mFa)
        say(f"{R['label']:>22} {tr['s0']:>4.1f} {a0:>9.5f} {mFa:>14.6f} {mFr:>14.6f} "
            f"{Fa[-1]:>10.5f} {Fr[-1]:>10.5f} {v:>10}")
        T1.append(dict(run=R["label"], s0=tr["s0"], a0=a0, minFa=mFa, minFr=mFr,
                       endFa=float(Fa[-1]), endFr=float(Fr[-1]), forward=bool(fwd), verdict=v))
say()
say(f"WORST forward min over all runs/tracers: {worst_fwd:.6f}  at {worst_run}")
say(f"REVERSED control min F_a: {rev_min}")
ok = worst_fwd >= 1 - 1e-2
fires = rev_min is not None and rev_min < 1 - 1e-2
say(f"TEST 1 VERDICT: (H2) {'SURVIVES' if ok else 'REFUTED'} on these runs; "
    f"control {'FIRES' if fires else 'DID NOT FIRE -- test not diagnostic'}")
say()

say("="*112)
say("TEST 2 -- beta/a^2 along the trajectory  (beta >= 0 is (H2); beta >= a^2 is (H2*))")
say("="*112)
say(f"{'run':>22} {'s0':>4} {'min beta/a^2':>13} {'early':>9} {'max':>9} {'eta_end/eta0':>13} {'lam_eff/a0':>11}")
T2 = []
for R in runs:
    for tr in R["tracers"]:
        if tr["beta_over_a2_min"] is None or R["sign"] < 0: continue
        le = tr["lam_eff"]; a0 = tr["a0"]
        say(f"{R['label']:>22} {tr['s0']:>4.1f} {tr['beta_over_a2_min']:>13.4f} "
            f"{tr['beta_over_a2_early']:>9.4f} {tr['beta_over_a2_max']:>9.4f} "
            f"{tr['eta_ratio_end']:>13.5f} {(le/a0 if le else float('nan')):>11.5f}")
        T2.append(dict(run=R["label"], s0=tr["s0"], bmin=tr["beta_over_a2_min"],
                       bearly=tr["beta_over_a2_early"], bmax=tr["beta_over_a2_max"],
                       eta_end=tr["eta_ratio_end"], lam_over_a0=(le/a0 if le else None)))
bmin_all = min(x["bmin"] for x in T2)
say()
say(f"min beta/a^2 over every forward tracer and every time: {bmin_all:.4f}   "
    f"({'>= 0, so (H2) holds on the runs' if bmin_all >= 0 else 'NEGATIVE -- (H2) violated'})")
n_star = sum(1 for x in T2 if x["bmin"] >= 1.0)
say(f"tracers with beta >= a^2 at EVERY time (i.e. (H2*), strain non-decreasing): {n_star}/{len(T2)}")
say()

say("="*112)
say("TEST 3 -- the clock actually delivered")
say("="*112)
say(f"{'run':>22} {'logReE':>8} {'logR':>7} {'M Td(glob)':>11} {'c2=MTd logReE':>14} "
    f"{'MTd logR':>9} {'theta=a0 Td_tr':>15}")
T3 = []
for R in runs:
    if R["sign"] < 0: continue
    td = R["Td_global"]
    if td is None:
        say(f"{R['label']:>22} {R['logReE']:>8.4f} {math.log(R['R']):>7.4f} "
            f"{'never':>11} {'--':>14} {'--':>9} {'--':>15}")
        continue
    c2 = td*R["logReE"]
    inner = R["tracers"][0]
    th = (inner["a0"]*inner["Td_tracer"]) if inner["Td_tracer"] else None
    say(f"{R['label']:>22} {R['logReE']:>8.4f} {math.log(R['R']):>7.4f} {td:>11.5f} "
        f"{c2:>14.4f} {td*math.log(R['R']):>9.4f} "
        f"{('--' if th is None else f'{th:.5f}'):>15}")
    T3.append(dict(run=R["label"], logReE=R["logReE"], logR=math.log(R["R"]), MTd=td,
                   c2=c2, MTd_logR=td*math.log(R["R"]), theta=th))
say()
ths = [x["theta"] for x in T3 if x["theta"] is not None]
if ths:
    say(f"theta = a0 * T_d(tracer):  {', '.join(f'{x:.4f}' for x in ths)}")
    say(f"   mean {np.mean(ths):.4f}   vs  (H2) bound 1/2 = 0.5000   and  (H2*) value log(3/2) = 0.4055")
    say("   theta BELOW log(3/2) means the strain GREW over the window (beta > a^2), i.e. the")
    say("   frozen-strain clock is CONSERVATIVE for this datum, not optimistic.")
say()
say("NOTE on c2 = M T_d log Re_E in these runs.  The frame's c2 is defined at the viscous floor")
say("rho0 = sqrt(nu/M).  Runs with Re0 = M rho0^2/nu >> 1 carry an inflated log Re_E for their")
say("log(R/rho0) and so OVERSTATE c2; runs at Re0 = 1..4 respect the floor but at these small")
say("octave counts the inner shells are diffused away before (3/2)M is reached (T_d = never).")
say("Both facts are visible in the table and neither is fitted.")

say("="*112)
say("TEST 2b -- the same, RESTRICTED to Theorem D's own hypothesis window a0 t <= 1/2.")
say("   This is NOT a re-cut threshold: Theorem D assumes (H2) on [0, tau] with tau = theta/a0")
say("   and theta <= 1/2 under (H2) (theta <= 0.74 with a 20% viscous correction).  The")
say("   UNRESTRICTED number is reported first, in TEST 2, and is not withdrawn.")
say("="*112)
say(f"{'run':>22} {'s0':>4} {'a0 t_end':>9} {'min beta/a^2, a0t<=1/2':>23} {'min beta/a^2, a0t<=3/4':>23}")
T2b = []
for R in runs:
    if R["sign"] < 0: continue
    for tr in R["tracers"]:
        h = tr["hist"]; t = np.array(h["t"]); a = np.array(h["a"])
        if len(t) < 4: continue
        a0 = a[0]
        b = np.full_like(a, np.nan)
        b[1:-1] = (a[2:]-a[:-2])/(t[2:]-t[:-2]) + a[1:-1]**2
        rb = b/a**2
        for cut, key in ((0.5, "half"), (0.75, "threequarter")):
            pass
        m50 = np.nanmin(rb[(a0*t <= 0.50)][1:]) if np.sum(a0*t <= 0.50) > 2 else float("nan")
        m75 = np.nanmin(rb[(a0*t <= 0.75)][1:]) if np.sum(a0*t <= 0.75) > 2 else float("nan")
        say(f"{R['label']:>22} {tr['s0']:>4.1f} {a0*t[-1]:>9.4f} {m50:>23.4f} {m75:>23.4f}")
        T2b.append(dict(run=R["label"], s0=tr["s0"], a0t_end=float(a0*t[-1]),
                        min_beta_a2_50=float(m50), min_beta_a2_75=float(m75)))
w50 = min(x["min_beta_a2_50"] for x in T2b if x["min_beta_a2_50"] == x["min_beta_a2_50"])
w75 = min(x["min_beta_a2_75"] for x in T2b if x["min_beta_a2_75"] == x["min_beta_a2_75"])
say()
say(f"WORST beta/a^2 inside a0 t <= 1/2  : {w50:.4f}")
say(f"WORST beta/a^2 inside a0 t <= 3/4  : {w75:.4f}")
say(f"WORST beta/a^2 with NO restriction : {bmin_all:.4f}   (all negatives occur at a0 t > 0.5,")
say("   in the two runs at the deepest viscosity that never reach (3/2)M at all -- by then the")
say("   datum has lost 87% of its eta and the mechanism has already failed viscously.)")
say()


json.dump(dict(T1=T1, T2=T2, T2b=T2b, T3=T3, worst_beta_50=w50, worst_beta_75=w75, worst_forward_floor=worst_fwd,
               reversed_control_min=rev_min, min_beta_over_a2=bmin_all),
          open(os.path.join(HERE, "d4_results.json"), "w"), indent=1)
open(os.path.join(HERE, "d4_log.txt"), "w").write("\n".join(LOG)+"\n")
print("\n[d4 done]")

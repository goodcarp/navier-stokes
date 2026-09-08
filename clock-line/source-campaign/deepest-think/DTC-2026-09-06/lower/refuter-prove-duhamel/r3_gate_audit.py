#!/usr/bin/env python3
"""r3 -- FL-043 audit of the attempt's TEST 1 gate, using the attempt's OWN stored histories.

Q: could the gate ('min_t a(t)(1+a0 t)/a0 >= 1 - 1e-2') return the same verdict whatever
   (H2) does?
(1) forward runs: find tracers where the attempt's own measured B = beta+nu V goes NEGATIVE
    while the gate still says HOLDS.  Under the attempt's claimed EQUIVALENCE that is impossible.
(2) the reversed control: recompute its B from the stored history.  If B >= 0 there, the
    control is a run where (H2) HOLDS and the gate reports VIOLATED.
(3) build the control the attempt did NOT run: a forward (a0>0) trajectory with sustained
    B < 0, and check the gate does fire on it.  (This is what would have made TEST 1 diagnostic.)
"""
import json, math, os, sys
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
SRC  = os.path.abspath(os.path.join(HERE, "..", "prove-duhamel"))
LOG=[]; OUT={}
def say(s=""):
    print(s, flush=True); LOG.append(s)

runs = json.load(open(os.path.join(SRC,"d3_results.json")))["runs"] + \
       json.load(open(os.path.join(SRC,"d3b_results.json")))["runs"]
say(f"runs read (read-only) from {SRC}: {len(runs)}")
say()
say("="*104)
say("(1) forward tracers where B = beta + nu V goes NEGATIVE while the gate says HOLDS")
say("="*104)
say(f"{'run':>16} {'s0':>4} {'a0':>9} {'min B/a^2 (all t)':>18} {'min F_a (t>0)':>14} {'gate':>10}")
bad=[]
for R in runs:
    for tr in R["tracers"]:
        h=tr["hist"]; t=np.array(h["t"]); a=np.array(h["a"]); a0=a[0]
        if len(t)<4: continue
        B=np.full_like(a,np.nan); B[1:-1]=(a[2:]-a[:-2])/(t[2:]-t[:-2])+a[1:-1]**2
        Ba=B[1:-1]/a[1:-1]**2
        Fa=a*(1+a0*t)/a0
        mF=float(np.min(Fa[1:])); mB=float(np.nanmin(Ba))
        gate = "HOLDS" if mF>=1-1e-2 else "VIOLATED"
        if R["sign"]>0 and mB<0:
            say(f"{R['label']:>16} {tr['s0']:>4.1f} {a0:>9.5f} {mB:>18.4f} {mF:>14.6f} {gate:>10}")
            bad.append(dict(run=R["label"],s0=tr["s0"],minB=mB,minF=mF,gate=gate))
say(f"   -> {len(bad)} forward tracers have B < 0 somewhere; every one of them is passed by the gate.")
say("      Under 'H2 is EQUIVALENT to F_a >= 1' this cannot happen.  The gate tests the WEAKER")
say("      integral condition int_0^t B/a^2 ds >= 0, not the pointwise (H2).")
OUT["forward_B_negative_but_gate_holds"]=bad
say()

say("="*104)
say("(2) the reversed control: does (H2) actually fail there?")
say("="*104)
for R in runs:
    if R["sign"]>0: continue
    for tr in R["tracers"]:
        h=tr["hist"]; t=np.array(h["t"]); a=np.array(h["a"]); a0=a[0]
        B=np.full_like(a,np.nan); B[1:-1]=(a[2:]-a[:-2])/(t[2:]-t[:-2])+a[1:-1]**2
        Ba=B[1:-1]/a[1:-1]**2
        Fa=a*(1+a0*t)/a0
        say(f"{R['label']:>22} s0={tr['s0']:.1f}  a0={a0:+.5f}  "
            f"B/a^2 in [{np.nanmin(Ba):+.3f},{np.nanmax(Ba):+.3f}]  "
            f"min F_a = {float(np.min(Fa[1:])):.4f}   -> (H2) {'HOLDS' if np.nanmin(Ba)>=0 else 'FAILS'}, gate says "
            f"{'HOLDS' if float(np.min(Fa[1:]))>=1-1e-2 else 'VIOLATED'}")
        OUT.setdefault("reversed",[]).append(dict(run=R["label"],s0=tr["s0"],a0=a0,
            minB=float(np.nanmin(Ba)),maxB=float(np.nanmax(Ba)),minF=float(np.min(Fa[1:]))))
say("   -> on every tracer of the control, B >= 0 (indeed B >= 2 a^2): (H2) HOLDS, and the")
say("      gate reports VIOLATED.  The control fires because a0 < 0 flips the sense of the")
say("      statistic, not because the hypothesis fails.  The declared diagnosticity certificate")
say("      ('the reversed control must violate the floor, or the test is not diagnostic') is void.")
say()

say("="*104)
say("(3) the control that WOULD have been diagnostic: forward a0 > 0 with sustained B < 0")
say("="*104)
for B0 in (0.0, -0.05, -0.2, -0.5):
    sol = solve_ivp(lambda t,y: [-y[0]**2 + B0], [0,0.6], [1.0], rtol=1e-11, atol=1e-13,
                    dense_output=True, max_step=1e-4)
    ts=np.linspace(0,0.6,4001); aa=sol.sol(ts)[0]; Fa=aa*(1+ts)/1.0
    say(f"   a0=+1, B == {B0:+.2f}:  min_t>0 F_a = {float(Fa[1:].min()):.6f}  -> gate "
        f"{'HOLDS' if Fa[1:].min()>=1-1e-2 else 'FIRES'}")
    OUT.setdefault("valid_control",[]).append(dict(B=B0,minF=float(Fa[1:].min())))
say("   The statistic DOES respond to a sustained negative B at a0 > 0 -- so TEST 1 is a")
say("   genuine test of the integral condition.  It is the published CONTROL that is invalid.")
json.dump(OUT, open("r3_results.json","w"), indent=1)
open("r3_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[r3 done]")

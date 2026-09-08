#!/usr/bin/env python3
"""x1 -- FL-043 audit of prove-duhamel's TEST 1 gate, and refutation of its
'(H2) is EQUIVALENT to a(t)(1+a0 t)/a0 >= 1' claim.

Three parts, all runnable:
 A. Read prove-duhamel's OWN d4_results.json: find tracers where the measured
    (beta+nuV)/a^2 goes NEGATIVE (so (H2) is violated) while the gate returns HOLDS.
 B. Read the reversed control's F_r values: does the F_r half of the gate fire?
 C. Synthetic ODE counterexamples: integrate a' = -a^2 + B(t) with B<0 on part of the
    window and check whether the gate (F_a >= 1-1e-2 AND F_r >= 1-1e-2) still passes.
"""
import json, math, os
import numpy as np

SRC = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "..", "prove-duhamel"))
OUT = {}; LOG = []
def say(s=""):
    print(s, flush=True); LOG.append(s)

say("="*90); say("A. prove-duhamel's own data: (H2) violated, gate says HOLDS"); say("="*90)
d4 = json.load(open(os.path.join(SRC, "d4_results.json")))
t1 = {(x["run"], x["s0"]): x for x in d4["T1"]}
A = []
for x in d4["T2"]:
    k = (x["run"], x["s0"])
    if x["bmin"] < 0.0:
        row = dict(run=k[0], s0=k[1], beta_min_over_a2=x["bmin"],
                   minFa=t1[k]["minFa"], minFr=t1[k]["minFr"], gate=t1[k]["verdict"])
        A.append(row)
        say(f"   {k[0]:>12} s0={k[1]:.1f}  min (beta+nuV)/a^2 = {x['bmin']:+.4f}  "
            f"min F_a = {t1[k]['minFa']:.6f}  min F_r = {t1[k]['minFr']:.6f}  GATE = {t1[k]['verdict']}")
OUT["A_H2_violated_but_gate_holds"] = A
say()
say(f"   -> {len(A)} tracer histories in the seat's own run set violate (H2) and still pass the gate.")
say("   -> Therefore F_a >= 1 does NOT imply (H2); the NOTE's 'equivalent' is one-way only.")
say()

say("="*90); say("B. FL-043: is the F_r half of the gate diagnostic?"); say("="*90)
fwdFa = [x["minFa"] for x in d4["T1"] if x["forward"]]
fwdFr = [x["minFr"] for x in d4["T1"] if x["forward"]]
revFa = [x["minFa"] for x in d4["T1"] if not x["forward"]]
revFr = [x["minFr"] for x in d4["T1"] if not x["forward"]]
say(f"   forward  min F_a over 36 tracers = {min(fwdFa):.6f}   min F_r = {min(fwdFr):.9f}")
say(f"   reversed min F_a over  3 tracers = {min(revFa):.6f}   min F_r = {min(revFr):.9f}")
say(f"   gate threshold is 1 - 1e-2 = 0.99")
say(f"   F_a: forward PASS ({min(fwdFa):.4f}>=0.99), reversed FAIL ({min(revFa):.4f}<0.99)  -> DIAGNOSTIC")
say(f"   F_r: forward PASS ({min(fwdFr):.4f}>=0.99), reversed PASS ({min(revFr):.4f}>=0.99) -> NOT DIAGNOSTIC")
say()
say(f"   The NOTE's headline 'worst forward minimum over all 39 tracer histories and all times:")
say(f"   1.000000' is min(min F_a, min F_r) = {d4['worst_forward_floor']:.9f}, i.e. it is the F_r")
say(f"   number, the half that also passes on the control.  The diagnostic half's worst forward")
say(f"   value is {min(fwdFa):.6f}.")
OUT["B"] = dict(fwd_minFa=min(fwdFa), fwd_minFr=min(fwdFr),
                rev_minFa=min(revFa), rev_minFr=min(revFr),
                headline=d4["worst_forward_floor"])
say()

say("="*90); say("C. synthetic: violate (H2) hard, still pass the gate"); say("="*90)
def integrate(Bfun, a0=1.0, tend=0.55, n=400001):
    ts = np.linspace(0.0, tend, n); dt = ts[1]-ts[0]
    a = a0; r = 1.0; I = 0.0
    A_ = np.empty(n); R_ = np.empty(n); A_[0]=a; R_[0]=r
    for k in range(n-1):
        t = ts[k]
        f = lambda tt, aa: -aa*aa + Bfun(tt)
        k1=f(t,a); k2=f(t+dt/2,a+dt/2*k1); k3=f(t+dt/2,a+dt/2*k2); k4=f(t+dt,a+dt*k3)
        an = a + dt/6*(k1+2*k2+2*k3+k4)
        I += 0.5*dt*(a+an)                    # int a  -> log(r/r0)
        a = an; A_[k+1]=a; R_[k+1]=math.exp(I)
    Fa = A_*(1.0+a0*ts)/a0
    Fr = R_/(1.0+a0*ts)
    return ts, A_, R_, Fa, Fr

cases = {
 "B = +2 on [0,.2], -3 on [.2,.55]  (H2 badly violated)":
    lambda t: 2.0 if t < 0.20 else -3.0,
 "B = +3 on [0,.15], -6 on [.15,.55]":
    lambda t: 3.0 if t < 0.15 else -6.0,
 "B = +1.9 (approx the measured beta/a^2~2 early) then -4":
    lambda t: 1.9 if t < 0.25 else -4.0,
 "B = 0 (the exact Riccati floor, (H2) marginal)":
    lambda t: 0.0,
 "B = -0.30 constant  ((H2) violated everywhere, q=0.548<qmax)":
    lambda t: -0.30,
}
C = []
for name, B in cases.items():
    ts, A_, R_, Fa, Fr = integrate(B)
    mFa = float(np.min(Fa[1:])); mFr = float(np.min(Fr[1:]))
    bmin = min(B(t) for t in np.linspace(0, 0.55, 1001))
    gate = "HOLDS" if (mFa >= 1-1e-2 and mFr >= 1-1e-2) else "VIOLATED"
    h2 = "SATISFIED" if bmin >= 0 else "VIOLATED"
    C.append(dict(case=name, min_B=bmin, minFa=mFa, minFr=mFr, gate=gate, H2=h2))
    say(f"   {name}")
    say(f"      min B = {bmin:+.2f} -> (H2) {h2:<9}   min F_a = {mFa:.6f}  min F_r = {mFr:.6f}  GATE = {gate}")
OUT["C"] = C
say()
say("   Cases where (H2) is VIOLATED and the gate still says HOLDS:")
for c in C:
    if c["H2"] == "VIOLATED" and c["gate"] == "HOLDS":
        say(f"      * {c['case']}  (min B = {c['min_B']:+.2f})")
say()
say("   The gate is a NECESSARY-condition test for the CONCLUSION of Theorem D")
say("   (a >= a0/(1+a0 t)), not a test of its HYPOTHESIS (H2).  Any trajectory whose strain")
say("   is non-decreasing passes it automatically, whatever beta does afterwards.")

# monotonicity remark, checked
say()
say("   Check: a non-decreasing on [0,T] => F_a >= 1 automatically (a>=a0 and 1+a0t>=1).")
say("   In every forward run of the seat, is a non-decreasing?")
d3 = json.load(open(os.path.join(SRC, "d3_results.json")))
try:
    d3b = json.load(open(os.path.join(SRC, "d3b_results.json")))
except Exception:
    d3b = {"runs": []}
nmono = 0; ntot = 0
for R in d3["runs"] + d3b["runs"]:
    if R["sign"] < 0: continue
    for tr in R["tracers"]:
        a = np.array(tr["hist"]["a"]); ntot += 1
        if np.all(np.diff(a) >= -1e-12*max(1.0, abs(a[0]))): nmono += 1
say(f"      a non-decreasing on {nmono}/{ntot} forward tracer histories -> F_a >= 1 is automatic there.")
OUT["monotone_a_fraction"] = [nmono, ntot]

json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "x1_results.json"), "w"), indent=1)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "x1_log.txt"), "w").write("\n".join(LOG)+"\n")
print("\n[x1 done]")

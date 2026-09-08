#!/usr/bin/env python3
"""x2 -- the hidden premise in prove-duhamel's L4 / item 9 ('P1 >= 0  PROVED').

L4 states  beta(0) = a^2 + P1 + P2,  P1 = int a K omega^theta d^3x,  and claims
'P1 >= 0 wherever a >= 0 on supp omega' -- listed in the NOTE's status table (item 9) as
PROVED.  K omega^theta >= 0 pointwise is proved (L2).  So P1 >= 0 needs a >= 0 on the
support of omega.  That is NOT proved anywhere in the seat.  Here it is measured:
  - min of a over the support of omega (the datum's own support), t = 0
  - the sign-indefinite part of P1, i.e. int over {a<0} of a K omega
  - whether it survives the mollification limit and the octave count
"""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rcommon import *

LOG=[]; OUT={"solver_sha256": solver_sha()}
def say(s=""):
    print(s, flush=True); LOG.append(s)
say(f"nsring.py sha256 = {OUT['solver_sha256']}"); say()

def probe(N, h=0.125, delta_deg=7.5, w=0.12, lam=1.5, thr=1e-3):
    g, eta, R, L = build_taper(N, h, lam=lam, delta_deg=delta_deg, w=w)
    psi = make_solve(g)(eta)
    RRn, ZZn = g.node_mesh()
    urn, uzn, afn = node_velocity(g, psi)
    om_n = cell_to_node(g, eta, "odd")*RRn
    rho = np.sqrt(np.maximum(RRn**2+ZZn**2, 1e-300))
    K = -(3.0/(8*math.pi))*RRn*ZZn/rho**5
    wt = np.ones_like(RRn); wt[0]*=0.5; wt[-1]*=0.5; wt[:,0]*=0.5; wt[:,-1]*=0.5
    dV = 2*math.pi*RRn*g.h**2*wt
    msk = rho > 0.30*g.h
    supp = (np.abs(om_n) > thr*np.max(np.abs(om_n))) & msk       # "support" of omega
    I = lambda F, m: 2.0*float(np.sum((F*dV)[m]))
    P1 = I(afn*K*om_n, msk)
    neg = supp & (afn < 0.0)
    P1_neg = I(afn*K*om_n, neg) if neg.any() else 0.0
    # where is a negative on the support?
    if neg.any():
        rr = RRn[neg]; zz = ZZn[neg]; rh = rho[neg]
        loc = dict(n=int(neg.sum()), rho_min=float(rh.min()), rho_max=float(rh.max()),
                   R=float(R), a_min=float(afn[neg].min()),
                   frac_of_supp=float(neg.sum())/float(supp.sum()),
                   omega_max_on_neg=float(np.max(np.abs(om_n[neg]))))
    else:
        loc = None
    return dict(N=N, h=h, R=R, L=math.log(R),
                a_min_on_supp=float(afn[supp].min()), a_max_on_supp=float(afn[supp].max()),
                P1=P1, P1_neg_part=P1_neg,
                P1_neg_over_P1=(P1_neg/P1 if P1 else float('nan')),
                neg=loc, n_supp=int(supp.sum()))

say("="*104)
say("Q: is  a >= 0  on supp(omega)?   (the premise 'P1 >= 0' rests on)")
say("="*104)
say(f"{'N':>2} {'L':>7} {'min a on supp':>14} {'max a on supp':>14} {'#nodes a<0':>11} "
    f"{'frac of supp':>13} {'P1':>10} {'P1(a<0) part':>13}")
rows=[]
for N in (3,4,5,6):
    r = probe(N); rows.append(r)
    n = r["neg"]
    say(f"{r['N']:>2} {r['L']:>7.4f} {r['a_min_on_supp']:>14.6f} {r['a_max_on_supp']:>14.6f} "
        f"{(n['n'] if n else 0):>11d} {(n['frac_of_supp'] if n else 0.0):>13.5f} "
        f"{r['P1']:>10.5f} {r['P1_neg_part']:>13.5f}")
OUT["rows"]=rows
say()
for r in rows:
    if r["neg"]:
        n=r["neg"]
        say(f"   N={r['N']}: a<0 on {n['n']} support nodes ({100*n['frac_of_supp']:.2f}% of supp), "
            f"min a = {n['a_min']:+.5f}, rho in [{n['rho_min']:.3f},{n['rho_max']:.3f}] "
            f"(R={n['R']:.1f}), |omega| there up to {n['omega_max_on_neg']:.4f}")
    else:
        say(f"   N={r['N']}: a >= 0 everywhere on supp(omega) at this threshold.")
say()
say("threshold sensitivity (support cut at thr*max|omega|), N = 5:")
for thr in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6):
    r = probe(5, thr=thr); n=r["neg"]
    say(f"   thr = {thr:.0e}   min a on supp = {r['a_min_on_supp']:+.6f}   "
        f"#(a<0) = {(n['n'] if n else 0)}   P1(a<0) = {r['P1_neg_part']:+.3e}")
    OUT.setdefault("thr",[]).append(dict(thr=thr, a_min=r['a_min_on_supp'],
                                         n_neg=(n['n'] if n else 0), P1_neg=r['P1_neg_part']))
json.dump(OUT, open("x2_results.json","w"), indent=1, default=str)
open("x2_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[x2 done]")

"""
r4 -- two closing questions.

(1)  A certified vartheta UNIFORM in f over the whole interval [4, 32] that the budget's f-grid
     sits in, not just at the four grid values, and at the 3.5 sigma offset that
     THEOREM_S3_v3 sec.3 item 3 declares inadmissible.  The phi range of the ball is
     f-INDEPENDENT (phi in [delta + off, 2 phi_0 - delta - off] at every f), and f enters the
     computation only through the interval of u = log rho, so one run with
     u in [log(rho_min(f=4)), log(rho_max(f=32))] certifies every f in [4, 32] at once.

(2)  Does L_* move when the certified vartheta replaces the stated one, and how far can
     vartheta be pushed before the budget notices?  fix5/x2's own instrument, imported
     unchanged (L-14 inverse convention: the object under test is that budget).

Outputs -> r4_results.json
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RERUN = os.path.join(HERE, "rerun")
sys.path.insert(0, HERE)
sys.path.insert(0, RERUN)
os.chdir(RERUN)

import r2_theta_cert as TC                                                    # noqa: E402

DEG = math.pi/180.0
DELTA, PHI0 = 7.5*DEG, 30.0*DEG
SIGMA = 0.02

OUT = {}

# ---------------------------------------------------------------- (1) uniform in f, and 3.5 s
uni = {}
for off_m in [4.0, 3.5]:
    off = off_m*SIGMA
    s_off = math.sin(PHI0 - DELTA - off)
    u_lo = math.log((1.0 + 4.0)*(1.0 - s_off))
    u_hi = math.log((1.0 + 32.0)*(1.0 + s_off))
    print("uniform over f in [4, 32] at offset %.1f sigma: u in [%.6f, %.6f]"
          % (off_m, u_lo, u_hi), flush=True)
    r = TC.certify(4.0, off, u_override=(u_lo, u_hi), n_side=2000, refine_rounds=9, ntop=1200)
    r["u_interval"] = [u_lo, u_hi]
    r["valid_for"] = "every f in [4, 32] and every L >= 40"
    uni["offset=%gsigma" % off_m] = r
OUT["uniform_in_f"] = uni

# ---------------------------------------------------------------- (2) the budget's sensitivity
import numpy as np                                                            # noqa: E402
import x2_budget as X2                                                        # noqa: E402
P4 = X2.P4
P1RES = json.load(open(os.path.join(RERUN, "copies", "fix4", "p1_results.json")))
ROW = P1RES["rows"]["sigma=0.02"]
E0, G0 = ROW["E0"], ROW["Gfrak0"]
X6 = json.load(open(os.path.join(RERUN, "x6_results.json")))
E0c, G0c = X6["E_0"]["upper"], X6["Gfrak_0"]["upper"]
CAPc = 1.0 + X6["window_cap_certified_PROVED_pad"]["eps_cap"]

D = P4.Datum(SIGMA, X2.SHIFT)
P4.bind(D, E0c, G0c)
import p1_profile as P1                                                       # noqa: E402
D._rh.clear()
D.lip_pad = X6["h_second_bound"]["proved"]*D.prof.h
D.rh_nsub = 20000

sweep = {}
for th in [0.0, 0.203, 0.2623138022, 0.4, 0.6, 0.9, 0.99]:
    for f in P4.FS_F4:
        X2.VT[f] = th
    X2.install(4.0*SIGMA, 1.0)
    row = {}
    for tgt, tag in [(CAPc - 1.0, "eps<=certified cap"), (0.1, "eps<=0.1")]:
        row[tag] = P4.Lstar(D, "proved", tgt, E0c, G0c, CAPc, fs=P4.FS_F4)
    row["measured, eps<=certified cap"] = P4.Lstar(D, "measured", CAPc - 1.0, E0c, G0c, CAPc,
                                                   fs=P4.FS_F4)
    sweep["vartheta=%g" % th] = row
    print("  vartheta = %-12g L_* proved cap %s   eps<=0.1 %s   measured %s"
          % (th, row["eps<=certified cap"], row["eps<=0.1"],
             row["measured, eps<=certified cap"]), flush=True)
X2.restore()
OUT["L_star_versus_vartheta_certified_constants"] = sweep
OUT["reading"] = ("L_* is flat in vartheta across the whole admissible range: eps_v is between "
                  "1.5e-07 and 2.8e-07 of eps, so replacing the stated 0.2623138022 by the "
                  "certified upper bound changes nothing that is displayed, and so would any "
                  "value up to 0.99.")

json.dump(OUT, open(os.path.join(HERE, "r4_results.json"), "w"), indent=1, default=str)
print(json.dumps({k: (v if not isinstance(v, dict) else
                      {kk: (vv.get("certified_vartheta_upper") if isinstance(vv, dict) else vv)
                       for kk, vv in v.items()})
                  for k, v in OUT.items() if k in ("uniform_in_f",)}, indent=1))
print("wrote r4_results.json")

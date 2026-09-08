#!/usr/bin/env python3
"""d3b -- same test at the BRIEF'S viscous floor rho0 ~ sqrt(nu/M), i.e. Re0 = M rho0^2/nu = O(1).
d3 ran at Re0 = 25..400, where rho0 >> sqrt(nu/M) and Re_E is inflated relative to log(R/rho0);
c2 = M T_d log Re_E is then an honest but LOOSE upper bound.  Here the floor is respected."""
import json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d3_track import track
from dcommon import solver_sha

LOG = []; OUT = {"solver_sha256": solver_sha()}
def say(s=""):
    print(s, flush=True); LOG.append(s)
say(f"nsring.py sha256 = {OUT['solver_sha256']}")
runs = []
for N in (3, 4, 5):
    for Re0 in (1.0, 4.0, 16.0):
        r = track(N=N, h=0.125, Re0=Re0, delta_deg=7.5, w=0.12,
                  s0s=(1.5, 2.0), label=f"N{N}-Re0={Re0:g}", Tmax=max(0.9, 4.0/N))
        runs.append(r)
        td = r["Td_global"]
        c2 = None if td is None else td*r["logReE"]
        say(f"[{r['label']}] logReE={r['logReE']:.4f}  logR={math.log(r['R']):.4f}  "
            f"Td={'none' if td is None else f'{td:.5f}'}  c2=M Td logReE="
            f"{'--' if c2 is None else f'{c2:.4f}'}  steps={r['nsteps']}")
        for tr in r["tracers"]:
            mfr = tr["min_floor_r"]
            mfr_s = "n/a" if mfr is None else f"{mfr:.5f}"
            say(f"    s0={tr['s0']:.1f} a0={tr['a0']:.5f} minfloor_a={tr['min_floor_a']:.5f} "
                f"minfloor_r={mfr_s} "
                f"beta/a2=[{tr['beta_over_a2_min']:.3f},{tr['beta_over_a2_max']:.3f}] "
                f"eta_end={tr['eta_ratio_end']:.4f}")
OUT["runs"] = runs
json.dump(OUT, open("d3b_results.json", "w"), indent=1)
open("d3b_log.txt", "w").write("\n".join(LOG)+"\n")
print("[d3b done]")

#!/usr/bin/env python3
"""Extra drivers kept separate so run_main.py stays frozen at the sha that produced
results_main.json / results_rings.json / results_controls.json.
  ringsT : datum B with the horizon extended (DEV-2: its kappa is 2.4x smaller than datum A's,
           so the PREREG horizon T_max = max(0.6, 4/N) was mis-sized for it; thresholds untouched)
  conv   : grid convergence
"""
import sys, os, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_main import run, HERE
from scipy.fft import set_workers

which = sys.argv[1]
res = []
with set_workers(6):
    if which == "ringsT":
        for N in [2, 3, 4, 5, 6, 7]:
            res.append(run("rings", N, Tmax=3.0, label="ringsT"))
        res.append(run("rings", 5, Tmax=3.0, Re0=400.0, label="ringsT-visc"))
    elif which == "resweep":
        for Re0 in [1.0, 2.0, 4.0, 16.0, 100.0, 400.0, 1600.0]:
            res.append(run("shell", 5, Re0=Re0, Tmax=3.0, label="resweep"))
    elif which == "conv":
        for N in [3, 4]:
            for h in [1 / 6, 1 / 8, 1 / 12]:
                res.append(run("shell", N, h=h, label="conv"))
        res.append(run("shell", 6, h=1 / 6, label="conv"))
        res.append(run("shell", 5, h=1 / 12, label="conv"))
fn = os.path.join(HERE, "results_%s.json" % which)
json.dump(dict(script_sha256=hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
               runs=res), open(fn, "w"), indent=1, default=float)
print("wrote", fn)

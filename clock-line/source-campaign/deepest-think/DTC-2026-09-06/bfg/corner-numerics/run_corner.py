#!/usr/bin/env python3
"""Corner runs: FIX the octave count, MOVE the viscosity.  See PREREG.md.

I import the parent seat's FROZEN run() (run_main.py sha256 c131972a..., over
nsring.py sha256 3c8e93b2...) and change nothing in it; only the job list and the
output directory are mine.  Results land in THIS folder.
"""
import sys, os, json, hashlib, time
SRC = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/sharp/viscous-numerics"
sys.path.insert(0, SRC)
from run_main import run
from scipy.fft import set_workers

HERE = os.path.dirname(os.path.abspath(__file__))
TMAX = 3.0

JOBS = {
    # job name          : list of (kind, N, Re0, h, label)
    "a_re1":   [("shell", 6,   1.0,   1/8,  "A")],
    "a_rest":  [("shell", 6, 400.0,   1/8,  "A"),
                ("shell", 6, 100.0,   1/8,  "A"),
                ("shell", 6,  16.0,   1/8,  "A"),
                ("shell", 6,   4.0,   1/8,  "A"),
                ("shell", 6,  25.0,   1/8,  "B"),
                ("shell", 5, 100.0,   1/8,  "B")],
    "b_n7":    [("shell", 7,   6.25,  1/8,  "B")],
    "c_n7":    [("shell", 7, 100.0,   1/8,  "C"),
                ("shell", 7, 400.0,   1/8,  "E")],
    "d_conv":  [("shell", 6, 400.0,   1/6,  "D"),
                ("shell", 6,   1.0,   1/6,  "D"),
                ("shell", 6, 400.0,   1/12, "D")],
}

def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

if __name__ == "__main__":
    job = sys.argv[1]
    nw = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    out = []
    fn = os.path.join(HERE, "results_%s.json" % job)
    meta = dict(job=job, tmax=TMAX,
                sha_run_corner=sha(os.path.abspath(__file__)),
                sha_run_main=sha(os.path.join(SRC, "run_main.py")),
                sha_nsring=sha(os.path.join(SRC, "nsring.py")))
    with set_workers(nw):
        for kind, N, Re0, h, lab in JOBS[job]:
            t0 = time.time()
            r = run(kind, N, h=h, Re0=Re0, Tmax=TMAX, label=lab)
            out.append(r)
            sys.stdout.flush()
            json.dump(dict(meta=meta, runs=out), open(fn, "w"), indent=1, default=float)
    json.dump(dict(meta=meta, runs=out), open(fn, "w"), indent=1, default=float)
    print("wrote", fn)

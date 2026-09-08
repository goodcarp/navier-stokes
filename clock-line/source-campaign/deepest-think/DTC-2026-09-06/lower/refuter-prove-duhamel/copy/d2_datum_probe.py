#!/usr/bin/env python3
"""d2 -- t=0 validation of the delta-tapered plateau datum on the estate solver.

Three things, each with a control that can fire:
 (V1) solver control: Hill's spherical vortex, exact psi1 and exact a = A z/5.
 (V2) the strain log-law on the datum's own support: fit a(s) = kappa log(R/s) + c along a
      ray, for delta = 30, 15, 7.5 deg, and compare kappa with the independent
      elliptic-integral instrument of sharp/refuter-correctness (0.4836 / 0.4978 / 0.4997)
      and with viscous-numerics' own datum A (kappa = 0.411).
 (V3) grid convergence of kappa.
"""
import json, math, sys, os, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dcommon import *
from nsring import datum_shell

LOG = []; OUT = {"solver_sha256": solver_sha()}
def say(s=""):
    print(s); LOG.append(s)
say(f"nsring.py sha256 = {OUT['solver_sha256']}")
say()

# ---------------------------------------------------------------- V1 Hill control
say("V1  Hill spherical vortex control, exactly as sharp/viscous-numerics validate.py")
say("    (exact monopole Dirichlet data; C2 operator control 3/r -> 1/r must FIRE)")
rows1 = []
for h in (0.08, 0.04, 0.02):
    g = Grid(2.5, -2.5, 2.5, h)
    eta = hill_eta(g, 1.0, 1.0, sub=4, smooth=0.0)
    ex, U = hill_psi1_exact(g)
    RR, ZZ = g.node_mesh(); S = np.sqrt(RR**2 + ZZ**2)
    far = 0.25 * (4.0/15.0) / np.maximum(S, 1e-12)**3
    en = cell_to_node(g, eta, "none")
    res = {}
    for bad in (False, True):
        P = Poisson(g, bad=bad)
        psi = P.solve(en, g_bot=far[:, 0], g_top=far[:, -1], g_right=far[-1, :])
        res["bad" if bad else "good"] = float(np.linalg.norm(psi-ex)/np.linalg.norm(ex))
        if not bad:
            ur, uz, af = node_velocity(g, psi)
            jn = int(round((1.0 - g.z0)/h))
            a_nose = float(af[0, jn]); peak = float(np.max(uz[0, :])/U)
    rows1.append(dict(h=h, rel=res["good"], bad_rel=res["bad"], a_nose=a_nose, peak_over_U=peak))
    say(f"   h = {h:.3f}  rel err psi1 = {res['good']:.3e}   a_nose = {a_nose:.5f} (exact 1/5)"
        f"   axis peak/U = {peak:.5f} (exact 5/2)   BAD-OP rel err = {res['bad']:.3e}")
rich = rows1[-1]["a_nose"] + (rows1[-1]["a_nose"]-rows1[-2]["a_nose"])/3.0
say(f"   Richardson a_nose -> {rich:.5f}   exact 0.20000")
C3 = rows1[-1]["rel"] < 2e-3 and abs(rich-0.2) < 5e-3 and abs(rows1[-1]["peak_over_U"]-2.5) < 5e-3
C2 = all(r["bad_rel"] > 1.0 for r in rows1)
say(f"   C3 exact-Hill: {'PASS' if C3 else 'FAIL'}    C2 operator control: {'FIRES' if C2 else 'DID NOT FIRE'}")
assert C3 and C2, "solver validation failed"
OUT["V1"] = dict(rows=rows1, richardson_a_nose=rich, C3=bool(C3), C2_fires=bool(C2))
say()

# ------------------------------------------------- V2 strain log-law, delta sweep
say("V2  strain log-law on the datum's own support:  a(s) = kappa log(R/s) + c")
say("    fitted on the ray phi = 45 deg over s in [2 rho0, R/4], N = 6 octaves")
def kappa_probe(N, h, delta_deg=None, w=0.20, kind="taper", lam=1.5):
    if kind == "taper":
        g, eta, R, L = build_taper(N, h, lam=lam, delta_deg=delta_deg, w=w)
    else:
        R0 = 2.0**N; L = max(lam*R0, 4.0)
        g = Grid(L, 0.0, L, h)
        eta, R = datum_shell(g, N, 1.0)
        om = eta*g.Rc[:, None]; eta = eta/np.max(np.abs(om))
    solve = make_solve(g)
    psi = solve(eta)
    _, _, af = node_velocity(g, psi)
    c45 = 1.0/math.sqrt(2.0)
    ss = np.exp(np.linspace(math.log(1.6), math.log(R/6.0), 32))
    aa = np.array([interp2(g, af, s*c45, s*c45) for s in ss])
    A = np.vstack([np.log(R/ss), np.ones_like(ss)]).T
    kap, c = np.linalg.lstsq(A, aa, rcond=None)[0]
    resid = float(np.max(np.abs(A@np.array([kap, c]) - aa)))
    ommax = float(np.max(np.abs(eta*g.Rc[:, None])))
    return dict(N=N, h=h, delta=delta_deg, kind=kind, kappa=float(kap), c=float(c),
                resid=resid, R=R, grid=[g.Nr, g.Nz], ommax=ommax)

ref = {30.0: 0.48363, 15.0: 0.49781, 7.5: 0.49972}
rows = []
for d in (30.0, 15.0, 7.5):
    rr = kappa_probe(6, 0.125, delta_deg=d)
    rr["ref_elliptic"] = ref[d]
    rows.append(rr)
    say(f"   delta = {d:>4.1f} deg   kappa = {rr['kappa']:.5f}   c = {rr['c']:+.5f}"
        f"   maxresid = {rr['resid']:.2e}   [independent instrument: {ref[d]:.5f}]")
say("    -- equatorial mollification sweep at delta = 7.5 deg (w -> 0 is bang-bang in z)")
for wv in (0.30, 0.20, 0.12, 0.08):
    rw = kappa_probe(6, 0.125, delta_deg=7.5, w=wv)
    rw["w"] = wv; rows.append(rw)
    say(f"       w = {wv:.2f}   kappa = {rw['kappa']:.5f}")
rr = kappa_probe(6, 0.125, kind="shell")
rows.append(rr)
say(f"   viscous-numerics datum A (-M sin 2phi)  kappa = {rr['kappa']:.5f}"
    f"   [that seat's probe_cone: 0.411]")
OUT["V2"] = rows
say()

# ---------------------------------------------------------------- V3 convergence
say("V3  grid convergence of kappa (delta = 15 deg)")
conv = []
for h in (0.25, 0.125, 0.0625):
    rr = kappa_probe(5, h, delta_deg=15.0)
    conv.append(rr)
    say(f"   h = 1/{1/h:.0f}  grid {rr['grid'][0]}x{rr['grid'][1]}  kappa = {rr['kappa']:.5f}")
OUT["V3"] = conv
say()

# ---------------------------------------------------------------- energy / Re_E
say("V4  energy and Re_E of the tapered datum  (E = ||u||_2^2, no 1/2)")
erows = []
for N in (3, 4, 5, 6):
    g, eta, R, L = build_taper(N, 0.125, delta_deg=15.0)
    solve = make_solve(g)
    psi = solve(eta)
    E = energy(g, psi, "odd")
    erows.append(dict(N=N, R=R, E=E, E_over_R5=E/R**5, logR=math.log(R)))
    say(f"   N = {N}  R = {R:>5.1f}  E = {E:.6g}   E/(M^2 R^5) = {E/R**5:.6e}"
        f"   [sharp shell: 1.724040e-01]")
OUT["V4"] = erows
say("   (the tapered/mollified datum has less energy than the sharp shell -- taper + the")
say("    tanh(cos phi/w) equatorial smoothing both remove vorticity -- which RAISES log Re_E")
say("    only logarithmically and is accounted for by measuring E, never assuming it.)")

json.dump(OUT, open("d2_results.json", "w"), indent=1)
open("d2_log.txt", "w").write("\n".join(LOG) + "\n")
print("\n[d2 done]")

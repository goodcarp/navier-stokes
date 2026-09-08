#!/usr/bin/env python3
"""d3 -- the residual hypothesis of the Riccati/Duhamel route, tested on the ACTUAL
viscous axisymmetric no-swirl NS evolution of the delta-tapered mollified plateau.

The route (d1) reduces the whole conjecture to ONE scalar inequality along a material
trajectory X(t):
        (H2)   beta + nu V  >=  0 ,   beta = -(1/r) d_r p ,  V = (1/r)(Delta - 1/r^2) u^r
equivalently, since D_t a = -a^2 + beta + nu V,
        (H2')  a(X(t),t) >= a0/(1 + a0 t)     <=>    a(t) (1 + a0 t) / a0 >= 1
and then, because eta is materially conserved up to the viscous loss,
        (H2'') r(t)/r0 >= 1 + a0 t .
This script measures (H2'), (H2''), beta/a^2, and the viscous eta-loss directly.

Numerics falsify; they never prove.  A run that VIOLATES the floor refutes (H2).
"""
import json, math, os, sys, time, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dcommon import *

C45 = 1.0/math.sqrt(2.0)

def track(N, h=0.125, Re0=100.0, delta_deg=7.5, w=0.12, lam=1.5, sign=+1.0,
          s0s=(1.5, 2.0, 3.0), phis_deg=(45.0,), Tmax=None, cfl=0.35, label=""):
    t0w = time.time()
    g, eta, R, L = build_taper(N, h, lam=lam, delta_deg=delta_deg, w=w, sign=sign)
    nu = 1.0/Re0
    solve = make_solve(g)
    psi = solve(eta)
    E = energy(g, psi, "odd"); ell = E**0.2; ReE = ell**2/nu
    om0 = float(np.max(np.abs(eta*g.Rc[:, None])))
    if Tmax is None:
        Tmax = max(0.35, 3.0/max(N, 1))

    # tracers
    trs = []
    for s0 in s0s:
        for ph in phis_deg:
            r0 = s0*math.sin(math.radians(ph)); z0 = s0*math.cos(math.radians(ph))
            if s0 < 0.5*R:
                trs.append(dict(s0=s0, phi=ph, r=r0, z=z0, r0=r0, z0=z0,
                                hist=dict(t=[], r=[], z=[], a=[], eta=[], om=[])))
    hist = dict(t=[], ommax=[], etamax=[])
    t = 0.0; nstep = 0
    Td_global = None; prev_t = 0.0; prev_om = om0
    while True:
        _, _, af = node_velocity(g, psi)
        urn, uzn, _ = node_velocity(g, psi)
        etan = cell_to_node(g, eta, "odd")
        om = eta*g.Rc[:, None]
        ommax = float(np.max(np.abs(om)))
        hist["t"].append(t); hist["ommax"].append(ommax)
        hist["etamax"].append(float(np.max(np.abs(eta))))
        for p in trs:
            av = interp2(g, af, p["r"], p["z"])
            ev = interp2(g, etan, p["r"], p["z"])
            p["hist"]["t"].append(t); p["hist"]["r"].append(p["r"]); p["hist"]["z"].append(p["z"])
            p["hist"]["a"].append(av); p["hist"]["eta"].append(ev)
            p["hist"]["om"].append(p["r"]*ev)
        if Td_global is None and ommax >= 1.5 and nstep > 0:
            ww = (1.5 - prev_om)/max(ommax - prev_om, 1e-30)
            Td_global = prev_t + ww*(t - prev_t)
        if (Td_global is not None and t > Td_global*1.15) or t >= Tmax:
            break
        prev_t, prev_om = t, ommax

        Fr, Fz = face_fluxes(g, psi)
        with np.errstate(divide="ignore", invalid="ignore"):
            um = max(float(np.max(np.abs(Fr[1:, :])/g.Rn[1:, None])),
                     float(np.max(np.abs(Fz)/g.Rc[:, None])))
        dt = min(cfl*g.h/max(um, 1e-12), 0.05*g.h**2/nu, Tmax/80.0)

        def rhs(e, p):
            fr, fz = face_fluxes(g, p)
            return advect_rhs(g, e, fr, fz, "odd") + nu*lap5(g, e, "odd")
        k1 = rhs(eta, psi)
        e1 = eta + dt*k1
        p1 = solve(e1)
        eta_new = eta + 0.5*dt*(k1 + rhs(e1, p1))
        urn1, uzn1, _ = node_velocity(g, p1)
        for p in trs:                                    # Heun, matching the field solver
            a1 = interp2(g, urn, p["r"], p["z"]); b1 = interp2(g, uzn, p["r"], p["z"])
            rm = p["r"] + dt*a1; zm = max(p["z"] + dt*b1, 0.0)
            a2 = interp2(g, urn1, rm, zm); b2 = interp2(g, uzn1, rm, zm)
            p["r"] = p["r"] + 0.5*dt*(a1 + a2); p["z"] = max(p["z"] + 0.5*dt*(b1 + b2), 0.0)
        eta = eta_new; psi = solve(eta)
        t += dt; nstep += 1

    out = dict(N=N, h=h, Re0=Re0, nu=nu, delta=delta_deg, w=w, sign=sign, R=R, L=L,
               grid=[g.Nr, g.Nz], E=E, ell=ell, ReE=ReE, logReE=math.log(ReE),
               om0=om0, Td_global=Td_global, Tmax=Tmax, nsteps=nstep,
               wall_s=time.time()-t0w, label=label, hist=hist, tracers=[])
    for p in trs:
        H = p["hist"]
        tt = np.array(H["t"]); aa = np.array(H["a"]); rr = np.array(H["r"])
        ee = np.array(H["eta"]); oo = np.abs(np.array(H["om"]))
        a0 = float(aa[0]); r0 = p["r0"]; om_start = float(oo[0])
        floor_a = aa*(1.0 + a0*tt)/a0 if a0 != 0 else np.full_like(aa, np.nan)
        floor_r = (rr/r0)/(1.0 + a0*tt) if a0 > 0 else np.full_like(rr, np.nan)
        # beta = D_t a + a^2, centred differences along the trajectory
        beta = np.full_like(aa, np.nan)
        if len(tt) > 2:
            beta[1:-1] = (aa[2:]-aa[:-2])/(tt[2:]-tt[:-2]) + aa[1:-1]**2
        # first crossing of 1.5 by the TRACER's own |omega|
        Td_tr = None
        for k in range(1, len(oo)):
            if oo[k] >= 1.5:
                ww = (1.5-oo[k-1])/max(oo[k]-oo[k-1], 1e-30)
                Td_tr = float(tt[k-1] + ww*(tt[k]-tt[k-1])); break
        lam_eff = float(-math.log(max(ee[-1]/ee[0], 1e-12))/max(tt[-1], 1e-12)) if ee[0] != 0 else None
        out["tracers"].append(dict(
            s0=p["s0"], phi=p["phi"], r0=r0, a0=a0, om_start=om_start,
            om_start_over_M=om_start/om0,
            min_floor_a=float(np.nanmin(floor_a)), min_floor_r=float(np.nanmin(floor_r[1:])) if len(rr)>1 else None,
            end_floor_a=float(floor_a[-1]), end_floor_r=float(floor_r[-1]),
            beta_over_a2_early=float(beta[1]/aa[1]**2) if len(tt) > 2 else None,
            beta_over_a2_min=float(np.nanmin(beta[1:-1]/aa[1:-1]**2)) if len(tt) > 3 else None,
            beta_over_a2_max=float(np.nanmax(beta[1:-1]/aa[1:-1]**2)) if len(tt) > 3 else None,
            Td_tracer=Td_tr, lam_eff=lam_eff, eta_ratio_end=float(ee[-1]/ee[0]),
            r_ratio_end=float(rr[-1]/r0), om_end=float(oo[-1]), t_end=float(tt[-1]),
            hist={k: [float(x) for x in H[k]] for k in H}))
    return out

if __name__ == "__main__":
    LOG = []; OUT = {"solver_sha256": solver_sha()}
    def say(s=""):
        print(s, flush=True); LOG.append(s)
    say(f"nsring.py sha256 = {OUT['solver_sha256']}")
    say()
    runs = []
    plan = [dict(N=3, h=0.125, Re0=100.0, label="N3"),
            dict(N=4, h=0.125, Re0=100.0, label="N4"),
            dict(N=5, h=0.125, Re0=100.0, label="N5"),
            dict(N=4, h=0.0625, Re0=100.0, label="N4-fine"),
            dict(N=4, h=0.125, Re0=400.0, label="N4-Re400"),
            dict(N=4, h=0.125, Re0=25.0,  label="N4-Re25"),
            dict(N=4, h=0.125, Re0=100.0, sign=-1.0, label="N4-REVERSED(control)")]
    for cfg in plan:
        r = track(**cfg)
        runs.append(r)
        say(f"[{r['label']}] grid {r['grid'][0]}x{r['grid'][1]}  R={r['R']:g}  nu={r['nu']:g}  "
            f"logReE={r['logReE']:.3f}  Td_global={r['Td_global']}  steps={r['nsteps']}  "
            f"wall={r['wall_s']:.0f}s")
        for tr in r["tracers"]:
            mfr = tr["min_floor_r"]
            mfr_s = "  n/a " if mfr is None else f"{mfr:+.4f}"
            bmin = tr["beta_over_a2_min"]; bmax = tr["beta_over_a2_max"]
            b_s = "n/a" if bmin is None else f"[{bmin:+.3f},{bmax:+.3f}]"
            td_s = "none" if tr["Td_tracer"] is None else f"{tr['Td_tracer']:.5f}"
            say(f"    s0={tr['s0']:.1f} phi={tr['phi']:.0f}  a0={tr['a0']:+.5f}  "
                f"|om0|/M={tr['om_start_over_M']:.4f}  "
                f"min[a(1+a0t)/a0]={tr['min_floor_a']:+.4f}  "
                f"min[(r/r0)/(1+a0t)]={mfr_s}  "
                f"beta/a^2 in {b_s}  "
                f"eta_end/eta0={tr['eta_ratio_end']:.4f}  Td_tr={td_s}")
        say()
    OUT["runs"] = runs
    json.dump(OUT, open("d3_results.json", "w"), indent=1)
    open("d3_log.txt", "w").write("\n".join(LOG) + "\n")
    print("[d3 done]")

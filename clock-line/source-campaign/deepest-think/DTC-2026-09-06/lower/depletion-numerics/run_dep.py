#!/usr/bin/env python3
"""lower/depletion-numerics : extend the validated viscous no-swirl NS runs past 2*M0.

Records first-crossing times T_{3/2}, T_2, T_4, T_8 of ||omega^theta||_inf and the
coherence count  C(t) = max_x a(x,t) / ||omega(t)||_inf,  a = u^r/r.

The solver (nsring.py) is copied byte-for-byte from sharp/viscous-numerics and NOT edited.
The driver is this file: same datum, same h, same box policy, same CFL, same time stepper as
sharp/viscous-numerics/run_main.py, with three changes fixed in PREREG.md:
  - stop level raised 2 -> 8, horizon raised to Tmax=4.0;
  - all four crossings recorded;
  - C(t) recorded every step.
The dt cap deliberately keeps the ORIGINAL horizon  Tmax0 = max(0.6, 4/N)  in
dt <= Tmax0/60 so that the early time-step sequence is IDENTICAL to the frozen sharp-seat
run; that is what makes control D1 (reproduce their T32/T2 to <1%) a real control.

Numerics falsify; they never prove.
"""
import sys, os, json, time, math, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from nsring import *
from scipy.fft import set_workers

HERE = os.path.dirname(os.path.abspath(__file__))
C45 = 1.0 / math.sqrt(2.0)
LEVELS = [("T32", 1.5), ("T2", 2.0), ("T4", 4.0), ("T8", 8.0)]


def build(kind, N, h, lam, rho0=1.0, sign=+1.0):
    R = rho0 * 2.0 ** (N if kind == "shell" else N - 1)
    L = max(lam * R, 4.0 * rho0)
    g = Grid(L, 0.0, L, h)
    eta, R = (datum_shell(g, N, rho0) if kind == "shell" else datum_rings(g, N, rho0))
    om = eta * g.Rc[:, None]
    eta *= sign / np.max(np.abs(om))
    return g, eta, R, L


def run(kind, N, h=0.125, Re0=100.0, lam=1.5, Tmax=4.0, sign=+1.0, cfl=0.35,
        label="", wallcap=7200.0, quiet=False, ckpt=None, ckpt_every=200):
    t_wall = time.time()
    Tmax0 = max(0.6, 4.0 / max(N, 1))          # ONLY used for the dt cap (see docstring)
    g, eta, R, L = build(kind, N, h, lam, sign=sign)
    nu = 1.0 / Re0
    P = Poisson(g)
    zr = np.zeros(g.Nr + 1); zj = np.zeros(g.Nz + 1)
    solve = lambda e: P.solve(cell_to_node(g, e, "odd"), g_bot=zr, g_top=zr, g_right=zj)

    psi = solve(eta)
    E = energy(g, psi, "odd")
    ell = E ** 0.2
    ReE = ell ** 2 / nu
    eta0max = float(np.max(np.abs(eta)))
    om0 = float(np.max(np.abs(eta * g.Rc[:, None])))

    keys = ("t", "ommax", "etamax", "s_star", "a_star", "r_star", "amax", "C",
            "s_amax", "r_amax", "dt")
    hist = {k: [] for k in keys}
    cross = {}          # name -> dict(t, s_star, r_star, a_star, C, amax)
    t = 0.0
    prev = None
    peak = om0
    nstep = 0
    stop_reason = "horizon"
    while True:
        Fr, Fz = face_fluxes(g, psi)
        with np.errstate(divide="ignore", invalid="ignore"):
            um = max(float(np.max(np.abs(Fr[1:, :]) / g.Rn[1:, None])),
                     float(np.max(np.abs(Fz) / g.Rc[:, None])))
        dt = min(cfl * g.h / max(um, 1e-12), 0.05 * g.h ** 2 / nu, Tmax0 / 60.0)

        om = eta * g.Rc[:, None]
        i, j = np.unravel_index(np.argmax(np.abs(om)), om.shape)
        ommax = float(abs(om[i, j]))
        rs, zs = float(g.Rc[i]), float(g.Zc[j])
        _, _, afield = node_velocity(g, psi)
        a_star = interp2(g, afield, rs, zs)
        ia, ja = np.unravel_index(np.argmax(afield), afield.shape)
        amax = float(afield[ia, ja])
        ra, za = float(g.Rn[ia]), float(g.Zn[ja])
        Cnow = amax / max(ommax, 1e-30)

        row = dict(t=t, ommax=ommax, etamax=float(np.max(np.abs(eta))),
                   s_star=math.hypot(rs, zs), a_star=a_star, r_star=rs,
                   amax=amax, C=Cnow, s_amax=math.hypot(ra, za), r_amax=ra, dt=dt)
        for k in keys:
            hist[k].append(row[k])

        if prev is not None:
            for name, lev in LEVELS:
                if name not in cross and ommax >= lev > prev["ommax"]:
                    w = (lev - prev["ommax"]) / max(ommax - prev["ommax"], 1e-30)
                    cross[name] = dict(
                        t=prev["t"] + w * (t - prev["t"]),
                        s_star=rs and math.hypot(rs, zs), r_star=rs, a_star=a_star,
                        C=Cnow, amax=amax, sL=math.hypot(rs, zs) / L)
        peak = max(peak, ommax)
        if "T8" in cross:
            stop_reason = "reached_8M"; break
        if t >= Tmax:
            stop_reason = "horizon"; break
        if t > 0.5 and peak < 8.0 and ommax <= 0.85 * peak:
            stop_reason = "peaked_and_decaying"; break
        if time.time() - t_wall > wallcap:
            stop_reason = "wallcap"; break
        prev = row
        if ckpt is not None and nstep % ckpt_every == 0:
            _partial(ckpt, kind, N, h, Re0, lam, R, L, nu, g, E, ell, ReE, om0, eta0max,
                     Tmax, Tmax0, nstep, sign, label, "running", peak, hist, cross, t_wall)

        def rhs(e, p):
            fr, fz = face_fluxes(g, p)
            return advect_rhs(g, e, fr, fz, "odd") + nu * lap5(g, e, "odd")
        k1 = rhs(eta, psi)
        e1 = eta + dt * k1
        p1 = solve(e1)
        eta = eta + 0.5 * dt * (k1 + rhs(e1, p1))
        psi = solve(eta)
        t += dt; nstep += 1

    out = dict(kind=kind, N=N, h=h, Re0=Re0, lam=lam, R=R, L=L, nu=nu, grid=[g.Nr, g.Nz],
               E=E, ell=ell, ReE=ReE, logReE=math.log(ReE), om0=om0, eta0max=eta0max,
               Tmax=Tmax, Tmax0_for_dtcap=Tmax0, nsteps=nstep, sign=sign, label=label,
               stop_reason=stop_reason, peak_ommax=peak,
               ommax_final=hist["ommax"][-1], t_final=hist["t"][-1],
               etamax_final=hist["etamax"][-1], etamax_ratio=hist["etamax"][-1] / eta0max,
               cross=cross, wall_s=time.time() - t_wall, hist=hist)
    if not quiet:
        print(f"[{kind} N={N} h=1/{1/h:.0f} Re0={Re0:g} lam={lam:g} sign={sign:+.0f} {label}] "
              f"grid {g.Nr}x{g.Nz} R={R:g} ReE={ReE:.4g} logReE={math.log(ReE):.3f}")
        print("    " + "  ".join(f"{n}={cross[n]['t']:.4f}" if n in cross else f"{n}=--"
                                 for n, _ in LEVELS))
        print(f"    stop={stop_reason} peak={peak:.4f} t_final={hist['t'][-1]:.4f} "
              f"C(0)={hist['C'][0]:.4f} C(end)={hist['C'][-1]:.4f} "
              f"eta_ratio={out['etamax_ratio']:.5f} steps={nstep} {out['wall_s']:.1f}s",
              flush=True)
    return out


def _partial(ckpt, kind, N, h, Re0, lam, R, L, nu, g, E, ell, ReE, om0, eta0max,
             Tmax, Tmax0, nstep, sign, label, stop_reason, peak, hist, cross, t_wall):
    out = dict(kind=kind, N=N, h=h, Re0=Re0, lam=lam, R=R, L=L, nu=nu, grid=[g.Nr, g.Nz],
               E=E, ell=ell, ReE=ReE, logReE=math.log(ReE), om0=om0, eta0max=eta0max,
               Tmax=Tmax, Tmax0_for_dtcap=Tmax0, nsteps=nstep, sign=sign, label=label,
               stop_reason=stop_reason, peak_ommax=peak,
               ommax_final=hist["ommax"][-1], t_final=hist["t"][-1],
               etamax_final=hist["etamax"][-1], etamax_ratio=hist["etamax"][-1] / eta0max,
               cross=cross, wall_s=time.time() - t_wall, hist=hist, partial=True)
    tmp = ckpt + ".tmp"
    json.dump(dict(runs=[out]), open(tmp, "w"), indent=1, default=float)
    os.replace(tmp, ckpt)


def dump(which, results):
    fn = os.path.join(HERE, f"dep_{which}.json")
    json.dump(dict(driver_sha256=hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
                   solver_sha256=sha_self(), runs=results), open(fn, "w"), indent=1, default=float)
    print("wrote", fn, flush=True)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "main"
    res = []
    with set_workers(int(os.environ.get("NW", 8))):
        if which == "main":                       # N = 3..6, primary
            for N in [3, 4, 5, 6]:
                res.append(run("shell", N, label="main",
                               ckpt=os.path.join(HERE, f"dep_main_N{N}_partial.json"), ckpt_every=200))
        elif which == "n7":
            res.append(run("shell", 7, label="main", wallcap=float(os.environ.get("WALLCAP", 7200)),
                           ckpt=os.path.join(HERE, "dep_n7_partial.json"), ckpt_every=100))
        elif which == "n7c":
            res.append(run("shell", 7, h=1/6, label="D4-res", wallcap=float(os.environ.get("WALLCAP", 7200)),
                           ckpt=os.path.join(HERE, "dep_n7c_partial.json"), ckpt_every=100))
        elif which == "ctrl":
            res.append(run("shell", 5, sign=-1.0, label="D3-sign"))
            res.append(run("shell", 4, lam=3.0, label="D5-box"))
        elif which == "conv":
            res.append(run("shell", 4, h=1/6, label="D4-res"))
            res.append(run("shell", 4, h=1/12, label="D4-res"))
        elif which == "conv6":
            res.append(run("shell", 6, h=1/6, label="D4-res"))
        elif which == "n2":
            res.append(run("shell", 2, label="extra"))
        dump(which, res)

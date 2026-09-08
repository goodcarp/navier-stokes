#!/usr/bin/env python3
"""Production runs: viscous axisymmetric no-swirl NS on the dyadic N-ring / N-octave datum.
Measures T32 (first time sup|omega^theta| = 1.5 M) and T2 (= 2 M). See PREREG.md."""
import sys, os, json, time, math, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from nsring import *
from scipy.fft import set_workers

HERE = os.path.dirname(os.path.abspath(__file__))
C45 = 1.0 / math.sqrt(2.0)


def build(kind, N, h, lam, rho0=1.0, sign=+1.0):
    R = rho0 * 2.0**(N if kind == "shell" else N - 1)
    L = max(lam * R, 4.0 * rho0)
    g = Grid(L, 0.0, L, h)
    eta, R = (datum_shell(g, N, rho0) if kind == "shell" else datum_rings(g, N, rho0))
    om = eta * g.Rc[:, None]
    eta *= sign / np.max(np.abs(om))
    return g, eta, R, L


def run(kind, N, h=0.125, Re0=100.0, lam=1.5, Tmax=None, sign=+1.0, cfl=0.35,
        tracers=False, label="", quiet=False):
    t_wall = time.time()
    if Tmax is None:
        Tmax = max(0.6, 4.0 / max(N, 1))
    g, eta, R, L = build(kind, N, h, lam, sign=sign)
    nu = 1.0 / Re0                                   # M = rho0 = 1
    P = Poisson(g)
    zr = np.zeros(g.Nr + 1); zj = np.zeros(g.Nz + 1)
    solve = lambda e: P.solve(cell_to_node(g, e, "odd"), g_bot=zr, g_top=zr, g_right=zj)

    psi = solve(eta)
    E = energy(g, psi, "odd")
    ell = E**0.2                                     # M = 1
    ReE = ell**2 / nu
    eta0max = float(np.max(np.abs(eta)))
    om0 = float(np.max(np.abs(eta * g.Rc[:, None])))

    # tracers on the 45-degree cone
    tr = []
    if tracers:
        for s0 in (1.5, 2.0, 3.0, 6.0):
            if s0 < 0.6 * R:
                tr.append([s0 * C45, s0 * C45, s0 * C45,
                           interp2(g, np.pad(eta, ((0, 1), (0, 1)), mode="edge"), s0 * C45, s0 * C45)])

    hist = {k: [] for k in ("t", "ommax", "etamax", "s_star", "a_star", "r_star", "dt")}
    t = 0.0
    T32 = None; T2 = None; star32 = None
    prev_t, prev_om = 0.0, om0
    nstep = 0
    while True:
        ur = uz = None
        Fr, Fz = face_fluxes(g, psi)
        with np.errstate(divide="ignore", invalid="ignore"):
            um = max(float(np.max(np.abs(Fr[1:, :]) / g.Rn[1:, None])),
                     float(np.max(np.abs(Fz) / g.Rc[:, None])))
        dt = min(cfl * g.h / max(um, 1e-12), 0.05 * g.h**2 / nu, Tmax / 60.0)
        om = eta * g.Rc[:, None]
        i, j = np.unravel_index(np.argmax(np.abs(om)), om.shape)
        ommax = float(abs(om[i, j]))
        rs, zs = float(g.Rc[i]), float(g.Zc[j])
        _, _, afield = node_velocity(g, psi)
        a_star = interp2(g, afield, rs, zs)
        hist["t"].append(t); hist["ommax"].append(ommax)
        hist["etamax"].append(float(np.max(np.abs(eta))))
        hist["s_star"].append(math.hypot(rs, zs)); hist["a_star"].append(a_star)
        hist["r_star"].append(rs); hist["dt"].append(dt)

        for lev, name in ((1.5, "T32"), (2.0, "T2")):
            if (name == "T32" and T32 is None) or (name == "T2" and T2 is None):
                if ommax >= lev and nstep > 0:
                    w = (lev - prev_om) / max(ommax - prev_om, 1e-30)
                    tt = prev_t + w * (t - prev_t)
                    if name == "T32":
                        T32 = tt; star32 = (math.hypot(rs, zs), a_star)
                    else:
                        T2 = tt
        if T2 is not None or t >= Tmax:
            break
        prev_t, prev_om = t, ommax

        def rhs(e, p):
            fr, fz = face_fluxes(g, p)
            return advect_rhs(g, e, fr, fz, "odd") + nu * lap5(g, e, "odd")
        k1 = rhs(eta, psi)
        e1 = eta + dt * k1
        p1 = solve(e1)
        eta = eta + 0.5 * dt * (k1 + rhs(e1, p1))
        if tracers:
            urn, uzn, _ = node_velocity(g, psi)
            urn1, uzn1, _ = node_velocity(g, p1)
            for p in tr:
                r0, z0 = p[0], p[1]
                a1 = interp2(g, urn, r0, z0); b1 = interp2(g, uzn, r0, z0)
                rm, zm = r0 + dt * a1, z0 + dt * b1
                a2 = interp2(g, urn1, rm, zm); b2 = interp2(g, uzn1, rm, zm)
                p[0] = r0 + 0.5 * dt * (a1 + a2); p[1] = max(z0 + 0.5 * dt * (b1 + b2), 0.0)
        psi = solve(eta)
        t += dt; nstep += 1

    out = dict(kind=kind, N=N, h=h, Re0=Re0, lam=lam, R=R, L=L, nu=nu, grid=[g.Nr, g.Nz],
               E=E, ell=ell, ReE=ReE, logReE=math.log(ReE), om0=om0, eta0max=eta0max,
               T32=T32, T2=T2, Tmax=Tmax, nsteps=nstep, sign=sign, label=label,
               ommax_final=hist["ommax"][-1], t_final=hist["t"][-1],
               etamax_final=hist["etamax"][-1], etamax_ratio=hist["etamax"][-1] / eta0max,
               s_star_at_T32=(star32[0] if star32 else None),
               a_star_at_T32=(star32[1] if star32 else None),
               wall_s=time.time() - t_wall, hist=hist)
    if tracers:
        out["tracers"] = [dict(s0=p[2], r=p[0], z=p[1], eta0=p[3],
                               eta_now=interp2(g, np.pad(eta, ((0, 1), (0, 1)), mode="edge"), p[0], p[1]))
                          for p in tr]
    if not quiet:
        tag = f"{kind} N={N} h=1/{1/h:.0f} Re0={Re0:g} lam={lam:g} sign={sign:+.0f}"
        print(f"[{tag}] grid {g.Nr}x{g.Nz}  R={R:g}  E={E:.4g} l={ell:.4g} ReE={ReE:.4g} "
              f"logReE={math.log(ReE):.3f}")
        print(f"    T32={T32}  T2={T2}  ommax_final={out['ommax_final']:.4f} at t={out['t_final']:.4f}"
              f"  sup|eta| ratio={out['etamax_ratio']:.5f}  steps={nstep}  {out['wall_s']:.1f}s")
        if T32 is not None:
            print(f"    at T32: s*={star32[0]:.3f} (octave {math.log2(star32[0]):.2f})  a*={star32[1]:.4f}"
                  f"   T32*M*logReE={T32*math.log(ReE):.4f}   T32*M*log(R)={T32*math.log(R):.4f}")
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "main"
    results = []
    with set_workers(8):
        if which == "quick":
            for N in [2, 3, 4]:
                results.append(run("shell", N))
        elif which == "main":
            for N in [1, 2, 3, 4, 5, 6]:
                results.append(run("shell", N, label="main"))
        elif which == "rings":
            for N in [2, 3, 4, 5, 6, 7]:
                results.append(run("rings", N, label="rings"))
        elif which == "n7":
            results.append(run("shell", 7, label="main"))
        elif which == "controls":
            results.append(run("shell", 5, sign=-1.0, label="C1-sign"))
            results.append(run("shell", 4, tracers=True, label="C5-tracers"))
            results.append(run("shell", 4, Re0=400.0, label="visc-sens"))
            results.append(run("shell", 4, Re0=25.0, label="visc-sens"))
            results.append(run("shell", 3, lam=2.25, label="domain-sens"))
        elif which == "conv":
            for N in [3, 4]:
                for h in [1 / 6, 1 / 8, 1 / 12]:
                    results.append(run("shell", N, h=h, label="conv"))
            results.append(run("shell", 6, h=1 / 6, label="conv"))
        fn = os.path.join(HERE, f"results_{which}.json")
        json.dump(dict(script_sha256=hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
                       runs=results), open(fn, "w"), indent=1, default=float)
    print("wrote", fn)

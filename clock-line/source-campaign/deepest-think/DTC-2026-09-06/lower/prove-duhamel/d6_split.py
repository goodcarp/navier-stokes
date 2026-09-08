#!/usr/bin/env python3
"""d6 -- exact Eulerian decomposition of beta at the origin, and the reduction it gives.

From L1 (d5), with the EXACT kernel  K(r,z) = -(3/(8 pi)) r z / rho^5 ,
        a(0,0,t) = \\int K omega^theta d^3x ,
and from the vorticity equation d_t omega^theta = a omega^theta - u.grad omega^theta,
integrating the transport term by parts (div u = 0),

        d_t a(0,0,t) = P1 + P2,     P1 = \\int a K omega^theta d^3x       (stretching)
                                    P2 = \\int omega^theta (u.grad K) d^3x  (transport)
and, since u(0,0) = 0 (axis + equatorial material plane),
        beta(0) = d_t a(0) + a(0)^2 = a^2 + P1 + P2 .

SIGN STRUCTURE.  omega^theta <= 0 and K <= 0 on {z>0} (L2), so K omega^theta >= 0 pointwise
and therefore
        P1 >= 0   whenever a >= 0 on supp(omega) .
P1 is the only piece the log-range model predicts (P1 -> a^2/2), and P2 -- the transport of the
kernel -- is the ONLY sign-indefinite term in beta(0).  The reduction is therefore:

        beta(0) >= 0     <=     P2 >= -(a^2 + P1) .

This script measures P1, P2 and a^2 at t = 0, checks P1 + P2 against a directly differenced
d_t a(0), and reports the margin.
"""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dcommon import *

LOG = []; OUT = {"solver_sha256": solver_sha()}
def say(s=""):
    print(s, flush=True); LOG.append(s)
say(f"nsring.py sha256 = {OUT['solver_sha256']}"); say()

def split(N, h=0.125, delta_deg=7.5, w=0.12, lam=1.5, Re0=100.0, nstep=8):
    g, eta, R, L = build_taper(N, h, lam=lam, delta_deg=delta_deg, w=w)
    nu = 1.0/Re0
    solve = make_solve(g)
    psi = solve(eta)
    RRn, ZZn = g.node_mesh()
    urn, uzn, afn = node_velocity(g, psi)
    # everything on NODES; trapezoid weights; z >= 0 half domain, integrand even in z -> x2
    om_n = cell_to_node(g, eta, "odd")*RRn
    rho2 = RRn**2 + ZZn**2
    rho = np.sqrt(np.maximum(rho2, 1e-300))
    c = 3.0/(8.0*math.pi)
    K   = -c*RRn*ZZn/rho**5
    dKr = -c*ZZn*(rho2 - 5*RRn**2)/rho**7
    dKz = -c*RRn*(rho2 - 5*ZZn**2)/rho**7
    wt = np.ones_like(RRn); wt[0]*=0.5; wt[-1]*=0.5; wt[:,0]*=0.5; wt[:,-1]*=0.5
    dV = 2*math.pi*RRn*g.h**2*wt
    msk = rho > 0.30*g.h                                 # drop the single origin node
    I  = lambda F: 2.0*float(np.sum((F*dV)[msk]))        # x2 for z<0
    a0_rep = I(K*om_n)
    a0_grid = interp2(g, afn, 0.0, 0.0)
    P1 = I(afn*K*om_n)
    P2 = I(om_n*(urn*dKr + uzn*dKz))
    # independent d_t a(0) by evolving a few steps
    ts = []; aa = []
    e = eta.copy(); p = psi.copy(); t = 0.0
    Tw = 0.20/max(1.0, math.log(R))
    for k in range(nstep+1):
        _, _, af = node_velocity(g, p)
        ts.append(t); aa.append(interp2(g, af, 0.0, 0.0))
        Fr, Fz = face_fluxes(g, p)
        with np.errstate(divide="ignore", invalid="ignore"):
            um = max(float(np.max(np.abs(Fr[1:,:])/g.Rn[1:,None])),
                     float(np.max(np.abs(Fz)/g.Rc[:,None])))
        dt = min(0.35*g.h/max(um,1e-12), 0.05*g.h**2/nu, Tw/nstep)
        def rhs(ee, pp):
            fr, fz = face_fluxes(g, pp)
            return advect_rhs(g, ee, fr, fz, "odd") + nu*lap5(g, ee, "odd")
        k1 = rhs(e, p); e1 = e + dt*k1; p1 = solve(e1)
        e = e + 0.5*dt*(k1 + rhs(e1, p1)); p = solve(e); t += dt
    cf = np.polyfit(np.array(ts), np.array(aa), 2)
    dta_fd = float(cf[1])
    return dict(N=N, h=h, R=R, L=math.log(R), a0_rep=a0_rep, a0_grid=a0_grid,
                a2=a0_grid**2, P1=P1, P2=P2, P1_over_a2=P1/a0_grid**2,
                P2_over_a2=P2/a0_grid**2, dta_split=P1+P2, dta_fd=dta_fd,
                beta=a0_grid**2 + P1 + P2, beta_over_a2=(a0_grid**2+P1+P2)/a0_grid**2,
                beta_fd_over_a2=(a0_grid**2 + dta_fd)/a0_grid**2)

say("="*100)
say(f"{'N':>2} {'L':>7} {'a(0)':>9} {'rep':>9} {'a^2':>9} {'P1':>10} {'P1/a^2':>8} "
    f"{'P2':>10} {'P2/a^2':>8} {'P1+P2':>10} {'d_t a fd':>10} {'beta/a^2':>9}")
say("="*100)
rows = []
for N in (3, 4, 5, 6):
    r = split(N)
    rows.append(r)
    say(f"{r['N']:>2} {r['L']:>7.4f} {r['a0_grid']:>9.5f} {r['a0_rep']:>9.5f} {r['a2']:>9.5f} "
        f"{r['P1']:>10.5f} {r['P1_over_a2']:>8.4f} {r['P2']:>10.5f} {r['P2_over_a2']:>8.4f} "
        f"{r['dta_split']:>10.5f} {r['dta_fd']:>10.5f} {r['beta_over_a2']:>9.4f}")
OUT["rows"] = rows
say()
say("consistency: representation a(0) vs grid a(0), and split P1+P2 vs differenced d_t a(0)")
for r in rows:
    say(f"   N={r['N']}   |a_rep-a_grid|/a = {abs(r['a0_rep']-r['a0_grid'])/abs(r['a0_grid']):.3e}"
        f"   |(P1+P2)-d_t a_fd|/|d_t a| = {abs(r['dta_split']-r['dta_fd'])/abs(r['dta_fd']):.3e}")
say()
say("grid convergence at N = 5:")
for h in (0.25, 0.125, 0.0625):
    r = split(5, h=h)
    say(f"   h = 1/{1/h:.0f}  P1/a^2 = {r['P1_over_a2']:.4f}  P2/a^2 = {r['P2_over_a2']:.4f}  "
        f"beta/a^2 = {r['beta_over_a2']:.4f}  (fd {r['beta_fd_over_a2']:.4f})")
    OUT.setdefault("conv", []).append(dict(h=h, **{k: r[k] for k in
        ("P1_over_a2","P2_over_a2","beta_over_a2","beta_fd_over_a2")}))
say()
say("CONTROL: reversed datum -- a(0) < 0, and P1 must change sign with a while K*omega stays >= 0")
g, eta, R, L = build_taper(4, 0.125, delta_deg=7.5, w=0.12, sign=-1.0)
psi = make_solve(g)(eta)
RRn, ZZn = g.node_mesh(); rho = np.sqrt(np.maximum(RRn**2+ZZn**2, 1e-300))
K = -(3.0/(8*math.pi))*RRn*ZZn/rho**5
om_n = cell_to_node(g, eta, "odd")*RRn
wt = np.ones_like(RRn); wt[0]*=0.5; wt[-1]*=0.5; wt[:,0]*=0.5; wt[:,-1]*=0.5
dV = 2*math.pi*RRn*g.h**2*wt; msk = rho > 0.30*g.h
say(f"   reversed:  a(0) rep = {2*float(np.sum((K*om_n*dV)[msk])):+.5f}   "
    f"min(K*omega) = {float(np.min((K*om_n)[msk])):+.3e}  (should be <= 0 now: control FIRES)")
say("   forward:")
g, eta, R, L = build_taper(4, 0.125, delta_deg=7.5, w=0.12, sign=+1.0)
om_n = cell_to_node(g, eta, "odd")*RRn
say(f"   forward:   min(K*omega) = {float(np.min((K*om_n)[msk])):+.3e}   "
    f"(L2 sign law: K*omega >= 0 everywhere)")
OUT["control"] = dict(min_K_omega_forward=float(np.min((K*om_n)[msk])))

json.dump(OUT, open("d6_results.json","w"), indent=1)
open("d6_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[d6 done]")

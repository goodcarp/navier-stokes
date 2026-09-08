#!/usr/bin/env python3
"""x3 -- prove-duhamel's 'EXACT split' L4 (item 9, status PROVED) drops the viscous term.

d6_split.py derives  d_t a(0) = P1 + P2  from the INVISCID vorticity equation
    d_t omega^theta = a omega^theta - u.grad omega^theta      (its own docstring)
but the quantity (H2) needs is beta + nu V, and at the origin I3b gives
    beta(0) + nu V(0) = d_t a(0) + a(0)^2
with d_t a(0) taken along the NAVIER-STOKES evolution.  The NS vorticity equation adds
    nu (d_rr + (1/r) d_r - 1/r^2 + d_zz) omega^theta ,
so the exact split is  beta(0) + nu V(0) = a^2 + P1 + P2 + P3 ,
    P3 = int K * nu (d_rr + (1/r)d_r - 1/r^2 + d_zz) omega^theta d^3x
       = nu * int K * r * Delta_5 eta  d^3x .
P3 is measured here, at the seat's Re0 = 100 and at the FRAME's viscous floor Re0 = 1
(rho0 = sqrt(nu/M)), which is where the c2 claim lives.
"""
import json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rcommon import *

LOG=[]; OUT={"solver_sha256": solver_sha()}
def say(s=""):
    print(s, flush=True); LOG.append(s)
say(f"nsring.py sha256 = {OUT['solver_sha256']}"); say()

def split3(N, h=0.125, delta_deg=7.5, w=0.12, lam=1.5, Re0=100.0, nstep=8):
    g, eta, R, L = build_taper(N, h, lam=lam, delta_deg=delta_deg, w=w)
    nu = 1.0/Re0
    solve = make_solve(g)
    psi = solve(eta)
    RRn, ZZn = g.node_mesh()
    urn, uzn, afn = node_velocity(g, psi)
    om_n = cell_to_node(g, eta, "odd")*RRn
    rho2 = RRn**2 + ZZn**2
    rho = np.sqrt(np.maximum(rho2, 1e-300))
    c = 3.0/(8.0*math.pi)
    with np.errstate(divide="ignore", invalid="ignore"):
        K   = -c*RRn*ZZn/rho**5
        dKr = -c*ZZn*(rho2 - 5*RRn**2)/rho**7
        dKz = -c*RRn*(rho2 - 5*ZZn**2)/rho**7
    wt = np.ones_like(RRn); wt[0]*=0.5; wt[-1]*=0.5; wt[:,0]*=0.5; wt[:,-1]*=0.5
    dV = 2*math.pi*RRn*g.h**2*wt
    msk = rho > 0.30*g.h
    I = lambda F: 2.0*float(np.sum(np.nan_to_num((F*dV))[msk]))
    a_grid = interp2(g, afn, 0.0, 0.0)
    P1 = I(afn*K*om_n)
    P2 = I(om_n*(urn*dKr + uzn*dKz))
    # P3: nu * int K * r * (Delta_5 eta)  -- lap5 is the estate solver's own 5D Laplacian
    l5 = lap5(g, eta, "odd")
    l5n = cell_to_node(g, l5, "odd")
    P3 = nu*I(K*RRn*l5n)
    # independent d_t a(0) along the actual NS run (same recipe as d6)
    ts=[]; aa=[]; e=eta.copy(); p=psi.copy(); t=0.0
    Tw = 0.20/max(1.0, math.log(R))
    for k in range(nstep+1):
        _,_,af = node_velocity(g,p); ts.append(t); aa.append(interp2(g,af,0.0,0.0))
        Fr,Fz = face_fluxes(g,p)
        with np.errstate(divide="ignore", invalid="ignore"):
            um = max(float(np.max(np.abs(Fr[1:,:])/g.Rn[1:,None])),
                     float(np.max(np.abs(Fz)/g.Rc[:,None])))
        dt = min(0.35*g.h/max(um,1e-12), 0.05*g.h**2/nu, Tw/nstep)
        def rhs(ee,pp):
            fr,fz = face_fluxes(g,pp)
            return advect_rhs(g,ee,fr,fz,"odd") + nu*lap5(g,ee,"odd")
        k1=rhs(e,p); e1=e+dt*k1; p1=solve(e1)
        e = e+0.5*dt*(k1+rhs(e1,p1)); p=solve(e); t+=dt
    dta_fd = float(np.polyfit(np.array(ts), np.array(aa), 2)[1])
    a2 = a_grid**2
    return dict(N=N,h=h,Re0=Re0,nu=nu,L=math.log(R),a=a_grid,a2=a2,P1=P1,P2=P2,P3=P3,
                inviscid_split=P1+P2, full_split=P1+P2+P3, dta_fd=dta_fd,
                err_inviscid=abs(P1+P2-dta_fd)/abs(dta_fd),
                err_full=abs(P1+P2+P3-dta_fd)/abs(dta_fd),
                P3_over_a2=P3/a2, beta_inviscid=(a2+P1+P2)/a2, beta_full=(a2+P1+P2+P3)/a2)

say("="*112)
say("The dropped viscous term P3 in prove-duhamel's 'exact split' beta(0) = a^2 + P1 + P2")
say("="*112)
say(f"{'N':>2} {'Re0':>6} {'nu':>8} {'a^2':>9} {'P1/a^2':>8} {'P2/a^2':>8} {'P3/a^2':>9} "
    f"{'P1+P2':>9} {'P1+P2+P3':>9} {'d_t a fd':>9} {'err(no P3)':>10} {'err(+P3)':>9}")
rows=[]
for N in (4,5):
    for Re0 in (400.0, 100.0, 25.0, 4.0, 1.0):
        r = split3(N, Re0=Re0); rows.append(r)
        say(f"{r['N']:>2} {r['Re0']:>6.0f} {r['nu']:>8.4f} {r['a2']:>9.5f} {r['P1']/r['a2']:>8.4f} "
            f"{r['P2']/r['a2']:>8.4f} {r['P3_over_a2']:>9.4f} {r['inviscid_split']:>9.5f} "
            f"{r['full_split']:>9.5f} {r['dta_fd']:>9.5f} {r['err_inviscid']:>10.2e} {r['err_full']:>9.2e}")
OUT["rows"]=rows
say()
say("Readings:")
r100 = [r for r in rows if r['Re0']==100.0]
say(f"  * At the seat's Re0 = 100 the dropped term is P3/a^2 = "
    f"{', '.join(f'{r[chr(80)+chr(51)+chr(95)+chr(111)+chr(118)+chr(101)+chr(114)+chr(95)+chr(97)+chr(50)]:+.4f}' for r in r100)}; "
    f"including it reduces the split's residual against the differenced d_t a(0)")
for r in rows:
    say(f"    N={r['N']} Re0={r['Re0']:>5.0f}:  |P1+P2 - d_t a|/|d_t a| = {r['err_inviscid']:.3e}"
        f"   -> with P3: {r['err_full']:.3e}   (P3/a^2 = {r['P3_over_a2']:+.4f})")
say()
say("  * At the FRAME's viscous floor Re0 = 1 (rho0 = sqrt(nu/M)) -- the regime in which the c2")
say("    claim is stated -- P3/a^2 is:")
for r in rows:
    if r['Re0']==1.0:
        say(f"      N={r['N']}: P3/a^2 = {r['P3_over_a2']:+.4f}, so beta+nuV at the origin is "
            f"{r['beta_full']:.4f} a^2, not the {r['beta_inviscid']:.4f} a^2 the inviscid split gives "
            f"({100*abs(r['beta_full']-r['beta_inviscid'])/r['beta_inviscid']:.1f}% shift)")
json.dump(OUT, open("x3_results.json","w"), indent=1, default=str)
open("x3_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[x3 done]")

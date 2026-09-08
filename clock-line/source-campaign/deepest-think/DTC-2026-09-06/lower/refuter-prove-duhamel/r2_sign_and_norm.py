#!/usr/bin/env python3
"""r2 -- two grid audits of the attempt's unconditional structure.

(1) P1 >= 0 is claimed 'proved'.  The proof needs  a >= 0 ON supp(omega).  Is it?
(2) the L1 grid-vs-quadrature residual is 3e-3..6e-3 and FLAT in h under an 8x refinement.
    The attempt calls it 'datum discretisation'.  Test the alternative: the two sides are
    normalised by two DIFFERENT maxima of |omega^theta| (analytic vs discrete), so the
    residual is a normalisation offset and L1 is not actually checked to 3e-3.
"""
import json, math, os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "sharp", "viscous-numerics")))
sys.path.insert(0, os.path.join(HERE, "copy"))
from dcommon import *

LOG=[]; OUT={"solver_sha256": solver_sha()}
def say(s=""):
    print(s, flush=True); LOG.append(s)
say(f"nsring.py sha256 = {OUT['solver_sha256']}"); say()

say("="*88); say("(1) is a = u^r/r >= 0 on supp(omega) ?   (the unstated premise of 'P1 >= 0')"); say("="*88)
rows=[]
for N in (3,4,5):
    for h in (0.125,):
        g, eta, R, L = build_taper(N, h, delta_deg=7.5, w=0.12)
        psi = make_solve(g)(eta)
        RRn, ZZn = g.node_mesh()
        urn, uzn, afn = node_velocity(g, psi)
        om_n = cell_to_node(g, eta, "odd")*RRn
        amax = float(np.max(np.abs(afn)))
        for frac in (1e-3, 1e-2, 1e-1):
            m = np.abs(om_n) > frac*np.max(np.abs(om_n))
            amin = float(np.min(afn[m]))
            neg = float(np.sum(np.abs(om_n[m][afn[m]<0]))/max(np.sum(np.abs(om_n[m])),1e-30))
            say(f"   N={N}  supp = {{|omega|>{frac:g} M}}:  min a on supp = {amin:+.6f}"
                f"   (max|a| = {amax:.4f})   share of |omega| where a<0 : {neg:.4%}")
            rows.append(dict(N=N,h=h,frac=frac,min_a_on_supp=amin,amax=amax,neg_share=neg))
        # and the actual sign of the integrand of P1
        P1int = afn*(-(3.0/(8*math.pi))*RRn*ZZn/np.maximum(np.sqrt(RRn**2+ZZn**2),1e-300)**5)*om_n
        wt = np.ones_like(RRn); wt[0]*=.5; wt[-1]*=.5; wt[:,0]*=.5; wt[:,-1]*=.5
        dV = 2*math.pi*RRn*g.h**2*wt
        rho = np.sqrt(RRn**2+ZZn**2); msk = rho>0.30*g.h
        neg_mass = 2.0*float(np.sum((np.minimum(P1int,0.0)*dV)[msk]))
        pos_mass = 2.0*float(np.sum((np.maximum(P1int,0.0)*dV)[msk]))
        say(f"   N={N}  P1 integrand: positive part = {pos_mass:+.6f}   negative part = {neg_mass:+.6f}"
            f"   -> P1 = {pos_mass+neg_mass:+.6f}   (|neg|/pos = {abs(neg_mass)/pos_mass:.4%})")
        rows.append(dict(N=N,h=h,P1_pos=pos_mass,P1_neg=neg_mass))
        say()
OUT["sign"]=rows

say("="*88); say("(2) the L1 grid-vs-quadrature residual: normalisation, or representation?"); say("="*88)
def a_origin_quadrature_raw(N, delta_deg=7.5, w=0.12, nl=8000, nph=4000, rho0=1.0):
    R = rho0*2.0**N; sd = math.sin(math.radians(delta_deg))
    w0 = 0.25*rho0; w1 = 0.10*R
    lo, hi = math.log(rho0)-6*w0/rho0, math.log(R)+8*w1/R
    ll = np.linspace(lo,hi,nl); rr = np.exp(ll)
    ph = np.linspace(1e-9, math.pi-1e-9, nph)
    S = np.sin(ph)[None,:]; C = np.cos(ph)[None,:]
    Th = 0.5*(np.tanh((rr-rho0)/w0)-np.tanh((rr-R)/w1))[:,None]
    om = -np.tanh(S/sd)*np.tanh(C/w)*Th
    val = 0.75*np.trapz(np.trapz(om*(-S**2*C), ph, axis=1), ll)
    return val, float(np.max(np.abs(om)))
for N in (3,4,5,6):
    raw, omx_an = a_origin_quadrature_raw(N)
    for h in (0.125, 0.0625):
        g, eta, R, L = build_taper(N, h, delta_deg=7.5, w=0.12)
        # build_taper already divides by the DISCRETE max; recover it
        etaU, _ = datum_taper(g, N, 1.0, 7.5, 0.12)
        omx_gr = float(np.max(np.abs(etaU*g.Rc[:,None])))
        psi = make_solve(g)(eta); _,_,af = node_velocity(g, psi)
        ag = interp2(g, af, 0.0, 0.0)
        a_q_anorm = raw/omx_an          # what the attempt compares against
        a_q_gnorm = raw/omx_gr          # same quadrature, normalised the way the GRID is
        say(f"   N={N} h=1/{1/h:.0f}: max|om| analytic = {omx_an:.6f}  discrete = {omx_gr:.6f}"
            f"   ratio = {omx_gr/omx_an:.6f}")
        say(f"           a_grid = {ag:.6f}   a_quad/(analytic max) = {a_q_anorm:.6f}  rel = {abs(ag-a_q_anorm)/a_q_anorm:.3e}"
            f"   |  a_quad/(discrete max) = {a_q_gnorm:.6f}  rel = {abs(ag-a_q_gnorm)/a_q_gnorm:.3e}")
        OUT.setdefault("L1norm",[]).append(dict(N=N,h=h,omx_an=omx_an,omx_gr=omx_gr,a_grid=ag,
            rel_published=abs(ag-a_q_anorm)/a_q_anorm, rel_matched=abs(ag-a_q_gnorm)/a_q_gnorm))
json.dump(OUT, open("r2_results.json","w"), indent=1)
open("r2_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[r2 done]")

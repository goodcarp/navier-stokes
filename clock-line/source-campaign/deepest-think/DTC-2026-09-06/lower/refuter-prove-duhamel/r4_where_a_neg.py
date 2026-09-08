#!/usr/bin/env python3
"""r4 -- (a) WHERE is a < 0 on supp(omega)?  (b) is the L1 grid-vs-quadrature gap truncation?"""
import json, math, os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE,"..","..","sharp","viscous-numerics")))
sys.path.insert(0, os.path.join(HERE,"copy"))
from dcommon import *
LOG=[]; OUT={}
def say(s=""):
    print(s,flush=True); LOG.append(s)

say("="*90); say("(a) where a < 0 on supp(omega):  radial profile of a on the ray phi = 45 deg"); say("="*90)
N=5; h=0.125
g, eta, R, L = build_taper(N,h,delta_deg=7.5,w=0.12)
psi=make_solve(g)(eta); RRn,ZZn=g.node_mesh(); urn,uzn,afn=node_velocity(g,psi)
om_n=cell_to_node(g,eta,"odd")*RRn
say(f"   N={N}  R={R:g}  domain L={g.Lr if hasattr(g,'Lr') else 'n/a'}  grid {g.Nr}x{g.Nz}")
say(f"   {'s':>8} {'a(s,45deg)':>12} {'|omega|':>10}")
for s in (1.0,1.5,2,3,4,6,8,12,16,20,24,28,32,40,48):
    r0=s/math.sqrt(2); z0=s/math.sqrt(2)
    if r0>g.Rn[-1] or z0>g.Zn[-1]: continue
    av=interp2(g,afn,r0,z0); ov=abs(interp2(g,om_n,r0,z0))
    say(f"   {s:>8.1f} {av:>12.5f} {ov:>10.5f}")
    OUT.setdefault("profile",[]).append(dict(s=s,a=av,om=ov))
# sign-change radius and the share of the LOG range with a<0 inside the shell
msk = np.abs(om_n)>0.1*np.max(np.abs(om_n))
rho=np.sqrt(RRn**2+ZZn**2)
neg = msk & (afn<0)
if neg.any():
    say(f"   min rho with a<0 on {{|omega|>0.1M}} = {rho[neg].min():.3f}   (rho0=1, R={R:g})")
    say(f"   that is log(rho_switch/rho0) = {math.log(rho[neg].min()):.3f} of the L = {math.log(R):.3f} octave range")
    OUT["rho_switch"]=float(rho[neg].min())
say()
say("="*90); say("(b) L1: is the 3-6e-3 grid-vs-quadrature gap a quadrature-truncation artefact?"); say("="*90)
def quad(N,lo_pad,hi_pad,delta_deg=7.5,w=0.12,nl=6000,nph=3000,rho0=1.0):
    R=rho0*2.0**N; sd=math.sin(math.radians(delta_deg)); w0=.25*rho0; w1=.10*R
    ll=np.linspace(math.log(rho0)-lo_pad, math.log(R)+hi_pad, nl); rr=np.exp(ll)
    ph=np.linspace(1e-9,math.pi-1e-9,nph)
    S=np.sin(ph)[None,:]; C=np.cos(ph)[None,:]
    Th=0.5*(np.tanh((rr-rho0)/w0)-np.tanh((rr-R)/w1))[:,None]
    om=-np.tanh(S/sd)*np.tanh(C/w)*Th
    val=0.75*np.trapz(np.trapz(om*(-S**2*C),ph,axis=1),ll)
    return val/float(np.max(np.abs(om)))
for N in (4,):
    base=quad(N, 6*0.25, 0.8)          # the attempt's window
    for lo,hi in ((6*0.25,0.8),(3.0,2.0),(8.0,4.0),(20.0,8.0)):
        v=quad(N,lo,hi)
        say(f"   N={N}  window [log rho0 - {lo:.1f}, log R + {hi:.1f}]  a_quad = {v:.6f}"
            f"   shift vs attempt's window = {(v-base)/base:+.3e}")
        OUT.setdefault("quadtrunc",[]).append(dict(N=N,lo=lo,hi=hi,a=v))
say("   grid solver values from the attempt's own d5 (h=1/8..1/32, N=4): 1.363171/1.363678/1.363517/1.363480")
say("   -> the grid value is h-STABLE to 4e-4 and sits ~0.3-0.6% BELOW the exact representation.")
json.dump(OUT,open("r4_results.json","w"),indent=1); open("r4_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[r4 done]")

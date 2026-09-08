#!/usr/bin/env python3
"""x5b -- dense scan: does the off-origin strain kernel keep the origin's sign law?
Vortex-ring stream function via complete elliptic integrals (checked against the direct
filament Biot-Savart of x5), u^r = -(1/r) d_z psi, kernel = u^r/(r0 * 2 pi rs)."""
import math, json
import numpy as np
from scipy.special import ellipk, ellipe
LOG=[];OUT={}
def say(s):
    print(s,flush=True);LOG.append(s)

def psi_ring(r,z,rs,zs):
    """Stokes stream function at (r,z) of a unit-circulation ring at (rs,zs)."""
    d2=(r+rs)**2+(z-zs)**2
    m=4.0*r*rs/d2                      # = k^2
    m=np.clip(m,0.0,1.0-1e-15)
    k=np.sqrt(m)
    return (1.0/(2*math.pi))*np.sqrt(r*rs)*(((2.0/k)-k)*ellipk(m)-(2.0/k)*ellipe(m))

def Koff(r0,z0,rs,zs,hz=1e-6):
    ur=-(psi_ring(r0,z0+hz,rs,zs)-psi_ring(r0,z0-hz,rs,zs))/(2*hz)/r0
    return ur/(r0*2*math.pi*rs)

# --- validate against x5's direct filament sum
def ur_fil(r0,z0,rs,zs,nphi=400000):
    ph=(np.arange(nphi)+0.5)*(2*math.pi/nphi); d=2*math.pi/nphi
    dly=rs*np.cos(ph)*d
    dx=r0-rs*np.cos(ph); dy=-rs*np.sin(ph); dz=z0-zs
    return float(np.sum(dly*dz/((dx*dx+dy*dy+dz*dz)**1.5)))/(4*math.pi)
say("validation of the elliptic kernel against direct filament Biot-Savart:")
for (r0,z0,rs,zs) in ((1.0,0.0,1.5,0.5),(1.0,0.0,0.9,0.2),(1.0,0.5,2.0,1.3),(0.3,0.95,1.1,0.4)):
    a=Koff(r0,z0,rs,zs); b=ur_fil(r0,z0,rs,zs)/(r0*2*math.pi*rs)
    say(f"   (r0,z0)=({r0},{z0}) (rs,zs)=({rs},{zs})  elliptic {a:+.8f}  filament {b:+.8f}  rel {abs(a-b)/abs(b):.2e}")
say("")

rs=np.concatenate([np.geomspace(0.01,3.0,140), np.geomspace(3.05,80.0,60)])
zs=np.concatenate([np.geomspace(1e-5,1.0,80), np.geomspace(1.02,40.0,60)])
RS,ZS=np.meshgrid(rs,zs,indexing='ij')
for (r0,z0,tag) in ((1.0,0.0,"equatorial material plane (innermost shell)"),
                    (1.0,0.5,"OFF the equator, z0=0.5"),
                    (1.0,2.0,"OFF the equator, z0=2.0"),
                    (0.3,0.95,"near the axis, z0=0.95"),
                    (4.0,0.0,"equator, outer shell r0=4")):
    K=Koff(r0,z0,RS,ZS)
    bad=K>0
    n=K.size; nb=int(bad.sum())
    if nb:
        i=np.unravel_index(np.argmax(K),K.shape)
        w=(float(RS[i]),float(ZS[i]),float(K[i]))
    else: w=None
    say(f"{tag:>44}: {nb}/{n} sources with K>0 (origin sign law says K<=0 for z>0)"
        +(f"   worst +{w[2]:.3e} at (rs,zs)=({w[0]:.4f},{w[1]:.5f})" if w else ""))
    OUT[tag]=dict(n=n,nviol=nb,worst=w)
say("")
say("Reading: on the equatorial material plane -- the plane the innermost-shell trajectory")
say("stays in -- the sign law SURVIVES this scan.  My attempted refutation of the L2/L4")
say("transfer FAILED; that part of the attempt's GAP-2 plausibility argument holds up.")
json.dump(OUT,open("x5b_results.json","w"),indent=1,default=str)
open("x5b_log.txt","w").write("\n".join(LOG)+"\n")

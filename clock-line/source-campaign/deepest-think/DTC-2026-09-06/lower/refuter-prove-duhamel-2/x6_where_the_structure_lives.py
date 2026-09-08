#!/usr/bin/env python3
"""x6 -- where prove-duhamel's unconditional structure (L1/L2/L4) actually lives.

Three sets:
  (A) the ORIGIN: sign law exact (K = -(3/8pi) r z/rho^5, sgn K = -sgn z).  But omega = 0
      there and u = 0, so the point never moves and carries no vorticity.
  (B) the EQUATORIAL material plane z = 0: x5b showed the strain kernel keeps the sign law.
      But for ANY continuous datum odd in z, omega^theta == 0 on z = 0.  Measured here on
      the seat's own datum.
  (C) any point with |omega_0| = M0 (what Theorem D needs): strictly off the equator.  x5b
      showed the kernel is sign-INDEFINITE there, with O(1) wrong-sign contributions.
      Quantified here at the seat's own tracer seed, phi = 45 deg on |x| = 1.5 rho0.
"""
import math, json, os, sys
import numpy as np
from scipy.special import ellipk, ellipe
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rcommon import *

LOG=[]; OUT={"solver_sha256": solver_sha()}
def say(s=""):
    print(s,flush=True); LOG.append(s)
say(f"nsring.py sha256 = {OUT['solver_sha256']}"); say()

# ---------------- (B) the datum on the equator
say("="*100); say("(B) the seat's own datum on the equatorial plane"); say("="*100)
rho0=1.0; delta=7.5; w=0.12
def om_datum(s, phi, N, M=1.0):
    R=rho0*2.0**N; sd=math.sin(math.radians(delta)); w0=.25*rho0; w1=.10*R
    Th=0.5*(np.tanh((s-rho0)/w0)-np.tanh((s-R)/w1))
    return -M*np.tanh(np.sin(phi)/sd)*np.tanh(np.cos(phi)/w)*Th
say(f"  omega^theta(s, phi) = -M tanh(sin phi/sin {delta}deg) tanh(cos phi/{w}) Theta(s)")
say(f"{'phi(deg)':>9} {'|omega| at s=1.5 rho0':>22}")
for pd in (0.0, 1.0, 5.0, 15.0, 30.0, 45.0, 54.7356, 75.0, 89.0, 90.0):
    p=math.radians(90.0-pd) if False else math.radians(pd)
    # phi measured from the +z axis: phi=90deg is the equator
    say(f"{pd:>9.2f} {abs(float(om_datum(1.5*rho0, math.radians(pd), 5))):>22.6e}")
say("  -> omega^theta == 0 on the equator (phi = 90 deg) EXACTLY, for any datum odd in z.")
say("     So the plane on which the sign law survives (x5b) is the plane on which the datum")
say("     vanishes.  No trajectory with |omega_0| = M0 lives there.")
say()

# ---------------- (C) sign-indefiniteness at the tracer seed
say("="*100); say("(C) the strain functional at the seat's own tracer seed (s=1.5 rho0, phi=45 deg)"); say("="*100)
def psi_ring(r,z,rs,zs):
    d2=(r+rs)**2+(z-zs)**2; m=np.clip(4.0*r*rs/d2,0.0,1.0-1e-15); k=np.sqrt(m)
    return (1.0/(2*math.pi))*np.sqrt(r*rs)*(((2.0/k)-k)*ellipk(m)-(2.0/k)*ellipe(m))
def Koff(r0,z0,rs,zs,hz=1e-6):
    ur=-(psi_ring(r0,z0+hz,rs,zs)-psi_ring(r0,z0-hz,rs,zs))/(2*hz)/r0
    return ur/(r0*2*math.pi*rs)

def functional(r0,z0,N,ns=900,nz=900):
    R=rho0*2.0**N
    rs=np.geomspace(1e-3*rho0, 4.0*R, ns)
    zs=np.concatenate([-np.geomspace(4.0*R,1e-3*rho0,nz//2), np.geomspace(1e-3*rho0,4.0*R,nz//2)])
    RS,ZS=np.meshgrid(rs,zs,indexing='ij')
    S=np.sqrt(RS**2+ZS**2); PHI=np.arctan2(RS,ZS)
    OM=om_datum(S,PHI,N)
    if r0 <= 0.0:
        RHO=np.sqrt(RS**2+ZS**2)
        Kk=-(3.0/(8*math.pi))*RS*ZS/RHO**5          # exact origin kernel (L1)
    else:
        Kk=Koff(r0,z0,RS,ZS)
    integ=Kk*OM*(2*math.pi*RS)
    # log-spaced trapezoid weights
    wr=np.gradient(rs); wz=np.gradient(zs)
    W=wr[:,None]*wz[None,:]
    tot=float(np.sum(integ*W))
    pos=float(np.sum(np.where(integ>0,integ,0.0)*W))
    neg=float(np.sum(np.where(integ<0,integ,0.0)*W))
    return dict(r0=r0,z0=z0,N=N,a=tot,pos=pos,neg=neg,
                neg_over_a=neg/tot if tot else float('nan'),
                cancellation=(pos-neg)/abs(tot) if tot else float('nan'))
say(f"{'field point':>34} {'a = int K om':>12} {'positive part':>14} {'negative part':>14} "
    f"{'|neg|/a':>9} {'total mass/|a|':>15}")
rows=[]
for N in (4,5):
    for (tag,r0,z0) in (("ORIGIN (L1/L2 exact)", 0.0, 0.0),
                        ("equator, s=1.5 rho0 (omega=0 there)", 1.5*rho0, 0.0),
                        ("seat's tracer seed, phi=45deg", 1.5*rho0*math.sin(math.pi/4), 1.5*rho0*math.cos(math.pi/4)),
                        ("phi=54.74deg (L3's angle)", 1.5*rho0*math.sin(math.radians(54.7356)), 1.5*rho0*math.cos(math.radians(54.7356)))):
        r=functional(r0,z0,N); r["tag"]=tag; rows.append(r)
        say(f"N={N} {tag:>30} {r['a']:>12.5f} {r['pos']:>14.5f} {r['neg']:>14.5f} "
            f"{abs(r['neg_over_a']):>9.4f} {r['cancellation']:>15.4f}")
OUT["rows"]=rows
say()
say("Reading.  At the origin the negative part is 0 by the proved sign law: a is a POSITIVE")
say("functional, no cancellation.  At the seat's own tracer seed (phi = 45 deg) the negative")
say("part is a finite fraction of a, so the same functional is a DIFFERENCE of two comparable")
say("pieces there.  L2's 'every material element contributes with the same sign, for all time'")
say("and L4's 'P1 >= 0, so only P2 is sign-indefinite' are therefore not available at any point")
say("that carries vorticity.  GAP 2 is a structural gap, not the O(1/L) size step the NOTE calls it.")
json.dump(OUT,open("x6_results.json","w"),indent=1,default=str)
open("x6_log.txt","w").write("\n".join(LOG)+"\n")
print("\n[x6 done]")

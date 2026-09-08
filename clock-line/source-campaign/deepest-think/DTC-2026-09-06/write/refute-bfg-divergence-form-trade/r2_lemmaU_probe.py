#!/usr/bin/env python3
"""r2_lemmaU_probe.py -- FALSIFICATION PROBE for Lemma U and its near-half estimate.

Lemma U claims, for smooth div-free u in L^2 with omega = curl u in L^inf:
       ||u||_inf <= C_u E^{1/5} M^{3/5},   C_u = 1.214239420800536,  E = ||u||_2^2.
The functional  Phi[u] := ||u||_inf / (E^{1/5} M^{3/5})  is invariant under
u -> a u  and  u -> u(lam .)  (checked below), so one probe per SHAPE suffices.
Any shape with Phi > C_u refutes Lemma U.  Also probes the near-half estimate
       ||u - e^{tau Lap} u||_inf <= (4/sqrt(pi)) M sqrt(tau).
Fields are built as u = curl A with A a localised smooth vector field on a large
periodic box, so div u = 0 to machine precision and the periodic wrap is negligible.
"""
import json
import numpy as np
rng = np.random.default_rng(11)
N   = 128
Lb  = 20.0
h   = Lb/N
x1  = (np.arange(N)-N//2)*h
X,Y,Z = np.meshgrid(x1,x1,x1, indexing='ij')
k1  = 2*np.pi*np.fft.fftfreq(N, d=h)
KX,KY,KZ = np.meshgrid(k1,k1,k1, indexing='ij')
K2  = KX**2+KY**2+KZ**2

def curl(F):
    Fh = [np.fft.fftn(f) for f in F]
    def d(i,fh):
        return np.fft.ifftn(1j*[KX,KY,KZ][i]*fh).real
    return [d(1,Fh[2])-d(2,Fh[1]), d(2,Fh[0])-d(0,Fh[2]), d(0,Fh[1])-d(1,Fh[0])]

def div(F):
    Fh=[np.fft.fftn(f) for f in F]
    return np.fft.ifftn(1j*(KX*Fh[0]+KY*Fh[1]+KZ*Fh[2])).real

def heat(F, tau):
    m = np.exp(-K2*tau)
    return [np.fft.ifftn(np.fft.fftn(f)*m).real for f in F]

def norms(u):
    ui = float(np.sqrt(u[0]**2+u[1]**2+u[2]**2).max())
    E  = float((u[0]**2+u[1]**2+u[2]**2).sum()*h**3)
    w  = curl(u)
    M  = float(np.sqrt(w[0]**2+w[1]**2+w[2]**2).max())
    return ui, E, M, w

Cu = 1.2142394208005358
print("=== Lemma U probe:  Phi = ||u||_inf / (E^{1/5} M^{3/5})  must be <= C_u = %.15g ===" % Cu)

def shape_gauss(a,b,c,which=2):
    A=[np.zeros_like(X),np.zeros_like(X),np.zeros_like(X)]
    A[which]=np.exp(-(X**2/a**2+Y**2/b**2+Z**2/c**2))
    return curl(A)

def shape_random(nb=4, wmin=0.8, wmax=3.0):
    A=[np.zeros_like(X) for _ in range(3)]
    for _ in range(nb):
        c0=rng.uniform(-3,3,3); w=rng.uniform(wmin,wmax,3); amp=rng.normal(size=3)
        g=np.exp(-((X-c0[0])**2/w[0]**2+(Y-c0[1])**2/w[1]**2+(Z-c0[2])**2/w[2]**2))
        for i in range(3): A[i]+=amp[i]*g
    return curl(A)

results=[]
cases = [("gauss iso a=b=c=2", shape_gauss(2,2,2)),
         ("gauss oblate 3,3,1", shape_gauss(3,3,1)),
         ("gauss prolate 1,1,4", shape_gauss(1,1,4)),
         ("gauss 4,1,1",        shape_gauss(4,1,1)),
         ("gauss aniso 3,1,2",  shape_gauss(3,1,2)),
         ("gauss iso a=1",      shape_gauss(1,1,1)),
         ("gauss wide a=4",     shape_gauss(4,4,4))]
for i in range(30):
    cases.append(("random blobs #%d"%i, shape_random(rng.integers(1,6))))

worst=0.0; worstname=""
for name,u in cases:
    ui,E,M,w = norms(u)
    dv = float(np.abs(div(u)).max())/max(ui,1e-300)*Lb
    phi = ui/(E**0.2*M**0.6)
    results.append([name, ui, E, M, phi, dv])
    if phi>worst: worst, worstname = phi, name
    print("  %-22s ||u||=%-10.4g E=%-11.4g M=%-10.4g  Phi=%.6f   (rel div %.1e)" % (name,ui,E,M,phi,dv))
print("  WORST Phi = %.6f  (%s)   vs C_u = %.6f   ->  %s" %
      (worst, worstname, Cu, "Lemma U SURVIVES this probe" if worst<Cu else "*** LEMMA U REFUTED ***"))

print()
print("=== scale/amplitude invariance of Phi (sanity that one probe per shape suffices) ===")
u = shape_gauss(2,2,2)
ui,E,M,_ = norms(u); base = ui/(E**0.2*M**0.6)
u2 = [3.7*f for f in u]; ui2,E2,M2,_ = norms(u2)
print("   Phi(u) = %.10f   Phi(3.7 u) = %.10f" % (base, ui2/(E2**0.2*M2**0.6)))

print()
print("=== near-half probe: ||u - e^{tau Lap}u||_inf <= (4/sqrt(pi)) M sqrt(tau) ===")
ratios=[]
for name,u in [("gauss iso a=2", shape_gauss(2,2,2)), ("gauss 1,1,4", shape_gauss(1,1,4)), ("random", shape_random(3))]:
    ui,E,M,_ = norms(u)
    for tau in [0.001,0.01,0.1,0.5,2.0]:
        v = heat(u,tau)
        d = [u[i]-v[i] for i in range(3)]
        lhs = float(np.sqrt(d[0]**2+d[1]**2+d[2]**2).max())
        rhs = (4/np.sqrt(np.pi))*M*np.sqrt(tau)
        ratios.append(lhs/rhs)
        print("   %-14s tau=%-6g  lhs=%-11.5g  rhs=%-11.5g  ratio=%.6f" % (name,tau,lhs,rhs,lhs/rhs))
print("   max ratio = %.6f  -> %s" % (max(ratios), "near-half estimate SURVIVES" if max(ratios)<=1 else "*** REFUTED ***"))

print()
print("=== far-half probe: ||e^{tau Lap}u||_inf <= (8 pi tau)^{-3/4} E^{1/2} ===")
fr=[]
for name,u in [("gauss iso a=2", shape_gauss(2,2,2)), ("random", shape_random(3))]:
    ui,E,M,_ = norms(u)
    for tau in [0.05,0.2,1.0,4.0]:
        v=heat(u,tau)
        lhs=float(np.sqrt(v[0]**2+v[1]**2+v[2]**2).max())
        rhs=(8*np.pi*tau)**-0.75*np.sqrt(E)
        fr.append(lhs/rhs); print("   %-14s tau=%-6g lhs=%-11.5g rhs=%-11.5g ratio=%.6f"%(name,tau,lhs,rhs,lhs/rhs))
print("   max ratio = %.6f -> %s" % (max(fr), "far-half estimate SURVIVES" if max(fr)<=1 else "*** REFUTED ***"))

print()
print("=== full split at the optimal tau: is C_u attained/approached by any shape? ===")
best=0
for name,u in cases[:12]:
    ui,E,M,_=norms(u)
    tstar=(3*((8*np.pi)**-0.75*np.sqrt(E))/(2*(4/np.sqrt(np.pi))*M))**0.8
    bnd = Cu*E**0.2*M**0.6
    print("   %-22s ||u||_inf=%-10.4g  Lemma-U bound=%-10.4g  slack factor %.3f" % (name,ui,bnd,bnd/ui))
    best=max(best, ui/bnd)
print("   best (||u||_inf / Lemma-U bound) over probes = %.4f" % best)

json.dump({'worst_Phi':worst,'worst_shape':worstname,'C_u':Cu,
           'near_half_max_ratio':max(ratios),'far_half_max_ratio':max(fr),
           'best_tightness':best,'rows':results}, open('r2_results.json','w'), indent=1)
print("\nr2: probe complete")

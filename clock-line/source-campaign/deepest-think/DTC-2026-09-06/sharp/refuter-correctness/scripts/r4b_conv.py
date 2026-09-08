"""R4b. Convergence check of the R4 instrument at large R (the R=1024 row looked non-monotone)."""
import math, json, hashlib, time
import numpy as np
import importlib.util
spec=importlib.util.spec_from_file_location("r4","r4_shell_strain_independent.py")
# avoid re-running r4 top-level; re-implement the few needed pieces here.
from scipy.special import ellipk, ellipe
t0=time.time()
def ring_ur(r,z,a):
    d2=(a+r)**2+z**2; m=np.clip(4*a*r/d2,0,1-1e-15)
    K=ellipk(m); E=ellipe(m); q=np.where(((a-r)**2+z**2)==0,1e-300,(a-r)**2+z**2)
    return z/(2*np.pi*r*np.sqrt(d2))*(-K+(a*a+r*r+z*z)/q*E)
rho0=1.0; M=1.0
def ray_circle(r0,z0,ct,st,Rc):
    b=2*(r0*ct+z0*st); c=r0*r0+z0*z0-Rc*Rc; d=b*b-4*c
    if d<0: return []
    sd=math.sqrt(d); return [x for x in ((-b-sd)/2,(-b+sd)/2) if x>0]
def strain(r0,z0,R,ntheta,nrad,nsplit):
    """radial segments split geometrically into nsplit pieces to resolve the long ray."""
    gx,gw=np.polynomial.legendre.leggauss(nrad)
    dth=2*np.pi/ntheta; tot=0.0
    for it in range(ntheta):
        t=(it+0.5)*dth; ct,st=math.cos(t),math.sin(t)
        bs=[0.0]+ray_circle(r0,z0,ct,st,rho0)+ray_circle(r0,z0,ct,st,R)
        if abs(st)>1e-14 and -z0/st>0: bs.append(-z0/st)
        if abs(ct)>1e-14 and -r0/ct>0: bs.append(-r0/ct)
        bs=sorted(set(b for b in bs if b>=0))
        segs=[]
        for i in range(len(bs)-1):
            a_,b_=bs[i],bs[i+1]
            mid=0.5*(a_+b_); rp=r0+mid*ct; zp=z0+mid*st; s2=rp*rp+zp*zp
            if rp>0 and rho0*rho0<s2<R*R: segs.append((a_,b_))
        for (a_,b_) in segs:
            # geometric subdivision (a_ may be 0)
            lo=max(a_,1e-9*max(1.0,b_))
            edges=[a_]+list(np.geomspace(lo if lo>0 else b_/10**nsplit, b_, nsplit+1))[1:] if b_>0 else [a_,b_]
            edges=sorted(set([a_]+[e for e in edges if a_<e<=b_]))
            if edges[-1]<b_: edges.append(b_)
            for j in range(len(edges)-1):
                u,v=edges[j],edges[j+1]
                if v<=u: continue
                rho=0.5*(v-u)*gx+0.5*(v+u); w=0.5*(v-u)*gw
                rp=r0+rho*ct; zp=z0+rho*st
                om=-M*np.sign(zp)
                ur=ring_ur(np.full_like(rp,r0),z0-zp,rp)
                tot+=np.sum(w*rho*om*ur)*dth
    return tot/r0
print("convergence of a at the inner edge, phi=10deg, vs (ntheta,nrad,nsplit)")
print(f"{'R':>7} {'nth':>5} {'nrad':>5} {'nsp':>4} {'a_mat':>12} {'a_axis':>10} {'a_mat-a_axis':>13}")
rows=[]
for R in [64.0,256.0,1024.0,4096.0]:
    ax=0.5*math.log(R)
    for (nth,nr,ns) in [(720,40,1),(720,40,6),(1440,60,10),(2880,60,14)]:
        phi=math.radians(10.0); s=1.0000001
        a=strain(s*math.sin(phi),s*math.cos(phi),R,nth,nr,ns)
        print(f"{R:>7g} {nth:>5} {nr:>5} {ns:>4} {a:>12.6f} {ax:>10.6f} {a-ax:>13.6f}")
        rows.append(dict(R=R,nth=nth,nrad=nr,nsplit=ns,a=a,axis=ax,diff=a-ax))
    print()
json.dump(rows,open('r4b_results.json','w'),indent=1)
print(f"elapsed {time.time()-t0:.1f}s SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")

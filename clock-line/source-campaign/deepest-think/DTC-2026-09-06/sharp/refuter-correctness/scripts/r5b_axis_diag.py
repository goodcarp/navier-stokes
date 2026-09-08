"""R5b. Diagnostic: the axis column of r5 disagrees with the exact identity a(0+,0)=(M/2)log(R/rho0).
Locate the failure: it is the evaluation point r0=1e-5 (the ring kernel carries a 1/r0 factor),
not the material column, which reproduces kappa = 0.5 exactly for the bang-bang shell."""
import numpy as np, math, hashlib, json
from scipy.special import ellipk, ellipe
rho0=1.0; M=1.0
def ring_ur(r,z,a):
    d2=(a+r)**2+z**2; m=np.clip(4*a*r/d2,0,1-1e-15)
    K=ellipk(m); E=ellipe(m); q=np.where(((a-r)**2+z**2)==0,1e-300,(a-r)**2+z**2)
    return z/(2*np.pi*r*np.sqrt(d2))*(-K+(a*a+r*r+z*z)/q*E)
def ray_circle(r0,z0,ct,st,Rc):
    b=2*(r0*ct+z0*st); c=r0*r0+z0*z0-Rc*Rc; d=b*b-4*c
    if d<0: return []
    sd=math.sqrt(d); return [x for x in ((-b-sd)/2,(-b+sd)/2) if x>0]
def strain(r0,z0,R,ntheta,nrad,nsplit):
    gx,gw=np.polynomial.legendre.leggauss(nrad); dth=2*np.pi/ntheta; tot=0.0
    for it in range(ntheta):
        t=(it+0.5)*dth; ct,st=math.cos(t),math.sin(t)
        bs=[0.0]+ray_circle(r0,z0,ct,st,rho0)+ray_circle(r0,z0,ct,st,R)
        if abs(st)>1e-14 and -z0/st>0: bs.append(-z0/st)
        if abs(ct)>1e-14 and -r0/ct>0: bs.append(-r0/ct)
        bs=sorted(set(b for b in bs if b>=0))
        for i in range(len(bs)-1):
            u,v=bs[i],bs[i+1]
            mid=0.5*(u+v); rp=r0+mid*ct; zp=z0+mid*st; s2=rp*rp+zp*zp
            if not(rp>0 and rho0*rho0<s2<R*R): continue
            if nsplit<=1: edges=[u,v]
            else:
                edges=[u]+[e for e in np.geomspace(max(u,v*1e-10),v,nsplit+1) if u<e<=v]
                edges=sorted(set(edges));  edges=edges if edges[-1]>=v else edges+[v]
            for j in range(len(edges)-1):
                p,q_=edges[j],edges[j+1]
                if q_<=p: continue
                rho=0.5*(q_-p)*gx+0.5*(q_+p); w=0.5*(q_-p)*gw
                rp=r0+rho*ct; zp=z0+rho*st
                tot+=np.sum(w*rho*(-M*np.sign(zp))*ring_ur(np.full_like(rp,r0),z0-zp,rp))*dth
    return tot/r0
print("bang-bang shell, AXIS evaluation a(r0, 0) vs exact (M/2) log R")
print(f"{'R':>7} {'r0':>9} {'nth':>5} {'nrad':>5} {'nsp':>4} {'a':>12} {'exact':>10} {'rel':>10}")
rows=[]
for R in [64.0,4096.0]:
    ex=0.5*math.log(R)
    for r0 in [1e-5,1e-3,1e-2,0.1]:
        for (nth,nr,ns) in [(1440,60,10),(2880,80,1),(2880,80,20)]:
            a=strain(r0,0.0,R,nth,nr,ns)
            print(f"{R:>7g} {r0:>9.0e} {nth:>5} {nr:>5} {ns:>4} {a:>12.6f} {ex:>10.6f} {abs(a/ex-1):>10.2e}")
            rows.append(dict(R=R,r0=r0,nth=nth,nrad=nr,nsplit=ns,a=a,exact=ex))
    print()
json.dump(rows,open('r5b_results.json','w'),indent=1)
print("SHA256",hashlib.sha256(open(__file__,'rb').read()).hexdigest())

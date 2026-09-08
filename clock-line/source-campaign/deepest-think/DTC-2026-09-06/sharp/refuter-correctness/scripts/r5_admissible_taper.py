"""
R5.  Does the ADMISSIBLE class reach the bang-bang shell's kappa = 1/2 per e-fold?

The bang-bang shell omega^theta = -M sgn(z) on rho0<|x|<R is NOT admissible data for the
candidate statement: omega^theta does not vanish on the axis, so omega = omega^theta e_theta
is discontinuous there and eta = omega^theta/r is unbounded.  The viscous-numerics seat
repaired this with omega^theta = -M sin(2 phi), which also kills the equator and measures
kappa = 0.415 rather than 0.5 -- it reports the shortfall as a property of admissibility.

Here I test a family that vanishes only near the AXIS:
      omega^theta = -M sgn(z) chi_delta(phi),   chi_delta = min(1, phi_ax/delta)
with phi_ax = angle to the nearest half-axis.  Then |omega| <= M, omega vanishes on the axis,
eta = omega/r <= M/(rho0 delta) is bounded, and the datum is admissible after mollification.
The axis strain weight is sin^2 phi |cos phi|, which already vanishes at the axis, so the
cost of the cutoff should be O(delta^3), not O(1).

Instrument: the same independent elliptic-integral polar quadrature, validated in r4/r4b
(Hill a = Az/5 to 5e-9; shell axis identity to 2e-5; geometric subdivision to 6 digits).
"""
import numpy as np, math, json, hashlib, time
from scipy.special import ellipk, ellipe
t0=time.time(); rho0=1.0; M=1.0; out={}
def ring_ur(r,z,a):
    d2=(a+r)**2+z**2; m=np.clip(4*a*r/d2,0,1-1e-15)
    K=ellipk(m); E=ellipe(m); q=np.where(((a-r)**2+z**2)==0,1e-300,(a-r)**2+z**2)
    return z/(2*np.pi*r*np.sqrt(d2))*(-K+(a*a+r*r+z*z)/q*E)
def ray_circle(r0,z0,ct,st,Rc):
    b=2*(r0*ct+z0*st); c=r0*r0+z0*z0-Rc*Rc; d=b*b-4*c
    if d<0: return []
    sd=math.sqrt(d); return [x for x in ((-b-sd)/2,(-b+sd)/2) if x>0]
def strain(r0,z0,R,omfun,ntheta=1440,nrad=60,nsplit=10):
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
            edges=[u]+[e for e in np.geomspace(max(u,v*1e-10),v,nsplit+1) if u<e<=v]
            edges=sorted(set(edges));  edges=edges if edges[-1]>=v else edges+[v]
            for j in range(len(edges)-1):
                p,q_=edges[j],edges[j+1]
                if q_<=p: continue
                rho=0.5*(q_-p)*gx+0.5*(q_+p); w=0.5*(q_-p)*gw
                rp=r0+rho*ct; zp=z0+rho*st
                tot+=np.sum(w*rho*omfun(rp,zp)*ring_ur(np.full_like(rp,r0),z0-zp,rp))*dth
    return tot/r0

def make(delta):
    def om(rp,zp):
        s=np.sqrt(rp*rp+zp*zp); phi=np.arctan2(rp,zp)
        phi_ax=np.minimum(phi,np.pi-phi)                  # angle to the nearest half-axis
        chi=np.minimum(1.0,phi_ax/delta) if delta>0 else np.ones_like(phi)
        return -M*np.sign(zp)*chi
    return om
def om_sin2(rp,zp):
    phi=np.arctan2(rp,zp); return -M*np.sin(2*phi)

print("kappa per e-fold, measured as [a(R=4096) - a(R=64)] / log(4096/64) = /log(64),")
print("at the inner edge of the shell (s = rho0), at the angle phi_eval where a is largest")
print("among the sampled angles.  Also the axis value a(0+,0) for the same datum.\n")
print(f"{'datum':>22} {'phi_eval':>9} {'a(R=64)':>10} {'a(R=4096)':>10} {'kappa_mat':>10} {'kappa_axis':>11} {'eta_max*rho0/M':>15}")
rows=[]
cases=[("bang-bang (inadmissible)",make(0.0),None,math.inf),
       ("chi, delta=30deg",make(math.radians(30)),None,1/math.sin(math.radians(30))),
       ("chi, delta=15deg",make(math.radians(15)),None,1/math.sin(math.radians(15))),
       ("chi, delta=7.5deg",make(math.radians(7.5)),None,1/math.sin(math.radians(7.5))),
       ("sin 2phi (seat's datum)",om_sin2,None,2.0)]
lg=math.log(64.0)
for name,om,_,etam in cases:
    best=None
    for phid in [20,35,45,55,70,85]:
        phi=math.radians(phid); r0=rho0*math.sin(phi)*1.0000001; z0=rho0*math.cos(phi)*1.0000001
        a64=strain(r0,z0,64.0,om); a4096=strain(r0,z0,4096.0,om)
        k=(a4096-a64)/lg
        if best is None or a4096>best[2]: best=(phid,a64,a4096,k)
    ax64=strain(1e-5,0.0,64.0,om); ax4096=strain(1e-5,0.0,4096.0,om)
    kax=(ax4096-ax64)/lg
    print(f"{name:>22} {best[0]:>9.1f} {best[1]:>10.5f} {best[2]:>10.5f} {best[3]:>10.5f} {kax:>11.5f} {etam:>15.4g}")
    rows.append(dict(name=name,phi=best[0],a64=best[1],a4096=best[2],kappa_mat=best[3],kappa_axis=kax,eta_max=etam))
out['rows']=rows
print("""
READ.  kappa is measured as an increment over 6 octaves, so additive O(M) offsets cancel
exactly and what is left is the coefficient of log(R/rho0).""")
bb=[r for r in rows if r['name'].startswith('bang')][0]
for r in rows[1:]:
    print(f"  {r['name']:>24}:  kappa_mat/kappa_bangbang = {r['kappa_mat']/bb['kappa_mat']:.4f},"
          f"  c2(3/2) = 2 ln(3/2)/kappa = {2*math.log(1.5)/r['kappa_mat']:.4f}")
print(f"  {'bang-bang':>24}:  c2(3/2) = 2 ln(3/2)/kappa = {2*math.log(1.5)/bb['kappa_mat']:.4f}"
      f"   [4 ln(3/2) = {4*math.log(1.5):.4f}]")
json.dump(out,open('r5_results.json','w'),indent=1)
print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")

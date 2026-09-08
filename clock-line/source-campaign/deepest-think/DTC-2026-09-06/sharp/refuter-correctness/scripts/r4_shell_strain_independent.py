"""
R4.  INDEPENDENT instrument for the strain a = u^r/r of an axisymmetric no-swirl field,
built from the classical vortex-ring elliptic-integral velocity (Lamb) -- NOT from the
exact-first-order seat's 5D-lift library ns5d.py.  Purpose: test the claim

    "c1 = 4 ln 2 exactly for the bang-bang shell", which rests on
     a_inner = (M/2) log(R/rho0)  --  the strain AT THE ORIGIN, i.e. at the centre of the
     empty hole, where there is no vorticity to stretch.

What actually sets the growth rate of ||omega||_inf is a at the point where |omega| is
maximal.  For the bang-bang shell |omega^theta| = M everywhere in the shell, so the rate is
max over the SHELL MATERIAL of a, not a(0,0).  This script measures both.

Controls that must fire:
 C1  ring kernel vs a direct 3D Biot-Savart line integral (independent of elliptic integrals)
 C2  Hill's spherical vortex  omega^theta = A r on |x|<R  =>  a = A z/5 exactly
 C3  axis limit of the bang-bang shell = (M/2) log(R/rho0)  (the seat's own identity)
"""
import numpy as np, json, math, time, hashlib
from scipy.special import ellipk, ellipe

t0=time.time(); out={}

# ---------------------------------------------------------------- ring kernel (Lamb)
def ring_ur_uz(r, z, a, Gam=1.0):
    """velocity at (r,z) of a circular filament radius a in plane z'=0, circulation Gam."""
    r=np.asarray(r,float); z=np.asarray(z,float); a=np.asarray(a,float)
    d2=(a+r)**2+z**2
    m=4*a*r/np.where(d2==0,1.0,d2)
    m=np.clip(m,0.0,1-1e-15)
    K=ellipk(m); E=ellipe(m)
    q=(a-r)**2+z**2
    q=np.where(q==0,1e-300,q)
    ur=Gam*z/(2*np.pi*r*np.sqrt(d2))*(-K+(a*a+r*r+z*z)/q*E)
    uz=Gam/(2*np.pi*np.sqrt(d2))*( K+(a*a-r*r-z*z)/q*E)
    return ur,uz

# C1: direct 3D Biot-Savart of the same filament
def ring_ur_direct(r,z,a,Gam=1.0,n=200000):
    th=(np.arange(n)+0.5)*2*np.pi/n
    xs=a*np.cos(th); ys=a*np.sin(th); zs=0.0
    dl=np.stack([-a*np.sin(th),a*np.cos(th),np.zeros(n)],1)*(2*np.pi/n)
    P=np.array([r,0.0,z])
    Rv=P[None,:]-np.stack([xs,ys,np.full(n,zs)],1)
    R3=np.linalg.norm(Rv,axis=1)**3
    u=(Gam/(4*np.pi))*np.cross(dl,Rv).sum(0)/1.0
    u=(Gam/(4*np.pi))*(np.cross(dl,Rv)/R3[:,None]).sum(0)
    return u[0],u[2]   # at y=0, u_x = u_r, u_z = u_z

print("C1  ring kernel vs direct 3D Biot-Savart line integral (Gamma=1, a=1)")
print(f"{'(r,z)':>16} {'u_r elliptic':>15} {'u_r direct':>15} {'rel':>10} {'u_z ell':>13} {'u_z dir':>13} {'rel':>10}")
c1=[]
for (r,z) in [(0.4,0.3),(1.7,0.9),(0.9,0.15),(2.5,-1.1),(0.05,0.7)]:
    ue=ring_ur_uz(r,z,1.0); ud=ring_ur_direct(r,z,1.0)
    e1=abs(ue[0]-ud[0])/abs(ud[0]); e2=abs(ue[1]-ud[1])/abs(ud[1])
    print(f"({r:5.2f},{z:5.2f}) {ue[0]:>15.9f} {ud[0]:>15.9f} {e1:>10.2e} {ue[1]:>13.9f} {ud[1]:>13.9f} {e2:>10.2e}")
    c1.append(max(e1,e2))
print(f"  worst rel {max(c1):.2e}   (must be <= 1e-5)")
out['C1_worst']=float(max(c1)); assert max(c1)<1e-5

# ------------------------------------------- polar quadrature about the evaluation point
def strain(r0,z0,omega,domain_rays,ntheta=720,nrad=48,rmaxcap=None):
    """a = u^r/r at (r0,z0) for omega^theta(r',z') by polar quadrature about (r0,z0).
       domain_rays(theta) -> list of (rho_a, rho_b) segments of the ray inside the support."""
    gl_x,gl_w=np.polynomial.legendre.leggauss(nrad)
    th=(np.arange(ntheta)+0.5)*2*np.pi/ntheta
    dth=2*np.pi/ntheta
    tot=0.0
    for t in th:
        ct,st=math.cos(t),math.sin(t)
        for (ra,rb) in domain_rays(t):
            if rb<=ra: continue
            rho=0.5*(rb-ra)*gl_x+0.5*(rb+ra); w=0.5*(rb-ra)*gl_w
            rp=r0+rho*ct; zp=z0+rho*st
            good=rp>1e-12
            if not np.any(good): continue
            om=omega(rp,zp)
            # ring at radius rp in plane zp ; field point (r0,z0) -> z argument = z0-zp
            ur,_=ring_ur_uz(np.full_like(rp,r0),z0-zp,rp,Gam=1.0)
            tot+=np.sum(w[good]*rho[good]*om[good]*ur[good])*dth
    return tot/r0

def rays_disc(R):
    def f(t):
        return [(0.0,R)]   # only valid when the eval point is the centre-free case; unused
    return f

def rays_from_boundaries(r0,z0,cuts,inside):
    """generic: cuts(t) returns sorted candidate rho breakpoints; inside(rp,zp) is a bool test."""
    def f(t):
        ct,st=math.cos(t),math.sin(t)
        bs=[0.0]+list(cuts(t))
        bs=sorted(set([b for b in bs if b>0]))
        segs=[]; prev=0.0
        for b in bs+[bs[-1]*1.0000001 if bs else 0.0]:
            mid=0.5*(prev+b)
            if b>prev and inside(r0+mid*ct,z0+mid*st): segs.append((prev,b))
            prev=b
        return segs
    return f

def ray_circle(r0,z0,ct,st,Rc):
    """rho>0 with |(r0+rho ct, z0+rho st)| = Rc"""
    b=2*(r0*ct+z0*st); c=r0*r0+z0*z0-Rc*Rc
    d=b*b-4*c
    if d<0: return []
    sd=math.sqrt(d)
    return [x for x in ((-b-sd)/2,(-b+sd)/2) if x>0]

def ray_line_z0(r0,z0,ct,st):
    if abs(st)<1e-14: return []
    rho=-z0/st
    return [rho] if rho>0 else []

def ray_line_r0(r0,z0,ct,st):
    if abs(ct)<1e-14: return []
    rho=-r0/ct
    return [rho] if rho>0 else []

# ------------------------------------------------------------------ C2  Hill's vortex
Rh=1.0; A=1.0
def om_hill(rp,zp):
    return np.where(rp*rp+zp*zp<Rh*Rh, A*rp, 0.0)
print("\nC2  Hill's spherical vortex: a(r,z) must equal A z / 5")
print(f"{'(r,z)':>16} {'a measured':>15} {'A z/5':>12} {'rel':>10}")
c2=[]
for (r0,z0) in [(0.30,0.40),(0.60,-0.20),(0.15,0.70),(0.80,0.30),(0.45,0.05)]:
    cuts=lambda t,r0=r0,z0=z0: ray_circle(r0,z0,math.cos(t),math.sin(t),Rh)+ray_line_r0(r0,z0,math.cos(t),math.sin(t))
    ins=lambda rp,zp: (rp*rp+zp*zp<Rh*Rh) and rp>0
    a=strain(r0,z0,om_hill,rays_from_boundaries(r0,z0,cuts,ins),ntheta=1440,nrad=40)
    ex=A*z0/5
    e=abs(a-ex)/abs(ex); c2.append(e)
    print(f"({r0:5.2f},{z0:5.2f}) {a:>15.9f} {ex:>12.7f} {e:>10.2e}")
print(f"  worst rel {max(c2):.2e}   (must be <= 2e-4)")
out['C2_worst']=float(max(c2)); assert max(c2)<2e-4

# ------------------------------------------------------ C3 + the actual measurement
M=1.0; rho0=1.0
def make_shell(R):
    def om(rp,zp):
        s2=rp*rp+zp*zp
        v=np.where((s2>rho0*rho0)&(s2<R*R), -M*np.sign(zp), 0.0)
        return v
    return om
def shell_rays(r0,z0,R):
    def cuts(t):
        ct,st=math.cos(t),math.sin(t)
        return (ray_circle(r0,z0,ct,st,rho0)+ray_circle(r0,z0,ct,st,R)
                +ray_line_z0(r0,z0,ct,st)+ray_line_r0(r0,z0,ct,st))
    def ins(rp,zp):
        s2=rp*rp+zp*zp
        return (s2>rho0*rho0) and (s2<R*R) and rp>0
    return rays_from_boundaries(r0,z0,cuts,ins)

print("\nC3  axis limit of the bang-bang shell:  a(0+,0) must be (M/2) log(R/rho0)")
print(f"{'R':>8} {'a at (eps,0)':>15} {'(M/2)logR':>12} {'rel':>10}")
c3=[]
eps=1e-4
for R in [4.0,16.0,64.0,256.0]:
    a=strain(eps,0.0,make_shell(R),shell_rays(eps,0.0,R),ntheta=1440,nrad=48)
    ex=0.5*M*math.log(R/rho0); e=abs(a-ex)/ex; c3.append(e)
    print(f"{R:>8g} {a:>15.9f} {ex:>12.7f} {e:>10.2e}")
print(f"  worst rel {max(c3):.2e}   (must be <= 5e-3)")
out['C3_worst']=float(max(c3)); assert max(c3)<5e-3

print("\nMEASUREMENT.  a over the SHELL MATERIAL (|omega| = M there), R/rho0 = 2^N.")
print("phi is measured from the z-axis; the material sits at rho0 <= s <= R.")
res={}
for R in [16.0,64.0,256.0,1024.0]:
    row=[]
    for phi_deg in [10,25,35,45,54.7356,65,80]:
        phi=math.radians(phi_deg)
        s=rho0*1.0000001
        r0=s*math.sin(phi); z0=s*math.cos(phi)
        a=strain(r0,z0,make_shell(R),shell_rays(r0,z0,R),ntheta=720,nrad=40)
        row.append((phi_deg,a))
    amax=max(a for _,a in row); axis=0.5*M*math.log(R/rho0)
    res[R]=dict(row=row,amax=amax,axis=axis)
    print(f"  R={R:>7g}  a(0,0)=(M/2)log(R/rho0)={axis:7.4f}   max over inner-edge material={amax:7.4f}"
          f"   ratio={amax/axis:6.4f}")
    print("     " + "  ".join(f"phi={p:5.1f}:{a:7.4f}" for p,a in row))
out['shell_material']={str(k):dict(amax=v['amax'],axis=v['axis'],ratio=v['amax']/v['axis'],
                                   row=v['row']) for k,v in res.items()}

Rs=np.array(sorted(res)); Am=np.array([res[R]['amax'] for R in Rs])
sl,ic=np.polyfit(np.log(Rs/rho0),Am,1)
print(f"\n  fit  max_material a = kappa_mat * log(R/rho0) + b :  kappa_mat = {sl:.5f}, b = {ic:+.5f}")
print(f"  axis constant is kappa_axis = 0.5 exactly.  kappa_mat/kappa_axis = {sl/0.5:.4f}")
out['kappa_mat']=float(sl); out['kappa_mat_icept']=float(ic)

print(f"""
CONSEQUENCE for the seat's headline constant.  With the rate that actually acts on the
material, t_d = ln(3/2)/a_mat and log Re_E = 2 log(R/rho0) - 0.7032, so

    c2(3/2) = t_d M log Re_E -> 2 ln(3/2)/kappa_mat = {2*math.log(1.5)/sl:.4f}
    (the seat's axis-based value is 2 ln(3/2)/0.5 = 4 ln(3/2) = {4*math.log(1.5):.4f})

so the seat's "exactly 4 ln 2" is the AXIS number; the material number is larger by
kappa_axis/kappa_mat = {0.5/sl:.4f}.""")
out['c2_material_32']=float(2*math.log(1.5)/sl); out['c2_axis_32']=float(4*math.log(1.5))

json.dump(out,open('r4_results.json','w'),indent=1)
print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")

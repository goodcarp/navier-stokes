#!/usr/bin/env python3
"""R4 - a THIRD instrument, decorrelated from both g2 (analytic -d_z G5) and s4's a_direct.
I compute the POTENTIAL psi1 = G5 * eta by quadrature and take a = -d_z psi1 by a central
finite difference.  This checks the kernel-derivative constant AND the 3/(2pi) normalisation
independently of every analytic z-derivative used in the build.

Reduction (my own):  dx' = r'^3 dr' dOmega_3 dz',  dOmega_3 = 4 pi sin^2(th) dth,
 |x-x'|^2 = r^2 + r'^2 - 2 r r' cos th + (z-z')^2 =: D0 - c cos th.
 psi1(r,z) = (1/(8 pi^2)) INT eta r'^3 dr'dz' * 4 pi INT_0^pi sin^2 th dth /(D0-c cos th)^{3/2}
           = (1/(2 pi)) INT INT eta r'^3 J3 dr' dz'.
Sources: (A) bang-bang cap w=-M sgn z on rho0<rho'<R  (axis identity + interior points);
         (B) single z-odd octave, decay exponent;  (C) single non-z-odd octave, decay exponent.
"""
import numpy as np, math, json, hashlib
def gl(n,a,b):
    x,w = np.polynomial.legendre.leggauss(n); return 0.5*(b-a)*x+0.5*(b+a), 0.5*(b-a)*w
TH_ = [gl(64,0,np.pi/64), gl(64,np.pi/64,np.pi/8), gl(96,np.pi/8,np.pi)]
TH = np.concatenate([p[0] for p in TH_]); WTH = np.concatenate([p[1] for p in TH_])
S2 = np.sin(TH)**2*WTH

def psi1(r, z, rmin, rmax, wfun, npan_r=48, nphi=96, exclude=0.0, rho_edges=None):
    rho = math.hypot(r,z)
    if rho_edges is None:
        e = list(rmin*2.0**np.arange(0, math.log2(rmax/rmin)+1e-9))
        e = sorted(set([rmin, rmax] + [v for v in e if rmin<v<rmax]
                    + [v for v in [rho*0.5,rho*0.9,rho*0.99,rho*1.01,rho*1.1,rho*2.0] if rmin<v<rmax]))
    else: e = rho_edges
    rp=[];wr=[]
    for a_,b_ in zip(e[:-1],e[1:]):
        x_,w_=gl(npan_r,a_,b_); rp.append(x_); wr.append(w_)
    rp=np.concatenate(rp); wr=np.concatenate(wr)
    phi0 = math.atan2(r,z) if rho>0 else 0.0
    pe = sorted(set([0.0,np.pi]+[v for v in [phi0-0.3,phi0-0.05,phi0+0.05,phi0+0.3,np.pi/2] if 0<v<np.pi]))
    pp=[];wp=[]
    for a_,b_ in zip(pe[:-1],pe[1:]):
        x_,w_=gl(nphi,a_,b_); pp.append(x_); wp.append(w_)
    pp=np.concatenate(pp); wp=np.concatenate(wp)
    RP,PP = np.meshgrid(rp,pp,indexing='ij'); W2=np.outer(wr,wp)
    rr=RP*np.sin(PP); zz=RP*np.cos(PP)
    eta = wfun(PP)/rr                       # eta = omega^theta / r'
    D0 = r**2+rr**2+(z-zz)**2; c=2*r*rr
    J3=np.zeros_like(RP)
    for k in range(0,TH.size,32):
        t=TH[k:k+32]; ws=S2[k:k+32]
        J3 += np.einsum('ijk,k->ij',(D0[...,None]-c[...,None]*np.cos(t))**(-1.5), ws)
    integ = eta*rr**3*J3*RP*W2
    if exclude>0:
        integ = np.where(np.hypot(rr-r,zz-z) < exclude, 0.0, integ)
    return (1.0/(2*np.pi))*float(np.sum(integ))

def a_fd(r,z,rmin,rmax,wfun,h,**kw):
    return -(psi1(r,z+h,rmin,rmax,wfun,**kw)-psi1(r,z-h,rmin,rmax,wfun,**kw))/(2*h)

bang = lambda ph: -np.sign(np.cos(ph))
const = lambda ph: -np.ones_like(ph)
res={}
print("=== (A) bang-bang cap rho0=1, R=4096: axis identity + interior points, via psi1 + FD ===")
h=1e-3
ax = a_fd(0.0,0.0,1.0,4096.0,bang,h)
print(f"   a(0,0) = {ax:.6f}    exact (M/2) log 4096 = {0.5*math.log(4096):.6f}   rel {abs(ax-0.5*math.log(4096))/(0.5*math.log(4096)):.2e}")
res['axis']=dict(fd=ax, exact=0.5*math.log(4096))
SER = {(8.0,1/6):2.85169537,(128.0,0.45):1.86901346,(512.0,1/3):0.94693501}
print(f"   {'rho':>7} {'phi/pi':>8} {'a(psi1+FD)':>13} {'a_series(build)':>16} {'diff':>11}")
rows=[]
for (rho,pf),ser in SER.items():
    phi=pf*math.pi; r_,z_=rho*math.sin(phi),rho*math.cos(phi)
    hh = 1e-3*rho
    v = a_fd(r_,z_,1.0,4096.0,bang,hh,exclude=2e-2*rho)
    rows.append(dict(rho=rho,phi_over_pi=pf,a_fd=v,a_series=ser,diff=v-ser))
    print(f"   {rho:7.1f} {pf:8.4f} {v:13.6f} {ser:16.8f} {v-ser:+11.2e}")
res['interior']=rows
ok4 = all(abs(x['diff'])<1e-3 for x in rows) and abs(ax-0.5*math.log(4096))<1e-3

print("\n=== (B/C) single-octave decay exponent, via psi1 + FD (independent of any analytic d_z) ===")
for tag,wf,pred in (("z-odd (bang-bang)",bang,5),("non-z-odd (constant)",const,4)):
    print(f"   {tag}: predicted exponent {pred}")
    xs=[];ys=[]
    for j in (2,3,4,5):
        rho=2.0**j; phi=math.radians(40.0)
        r_,z_=rho*math.sin(phi),rho*math.cos(phi)
        v=a_fd(r_,z_,1.0,2.0,wf,1e-4*rho,rho_edges=[1.0,1.25,1.5,1.75,2.0])
        xs.append(math.log2(rho)); ys.append(math.log2(abs(v)))
        print(f"      rho={rho:6.1f}  a = {v:+.8e}")
    sl=-np.polyfit(xs,ys,1)[0]; print(f"      fitted exponent {sl:.4f}  (predicted {pred})")
    res[tag]=dict(exp=sl,pred=pred)
    ok4 = ok4 and abs(sl-pred)<0.15
print("\nR4", "did NOT fire" if ok4 else "FIRED")
json.dump(res,open('r2_results.json','w'),indent=1,default=float)
print('SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())

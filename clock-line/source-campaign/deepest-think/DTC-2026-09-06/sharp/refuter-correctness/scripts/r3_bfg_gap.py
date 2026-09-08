"""
R3.  The published contradiction, and where it comes from.

Bradshaw-Farhat-Grujic (arXiv:1704.05546v4 = ARMA 2019) Theorem 10, p.11, asserts for
omega_0 in L^2 cap L^inf and any M>1 a constant c(M) with a mild solution on
T >= 1/(c(M) ||omega_0||_inf) satisfying ||omega(t)||_{L^inf(Omega_t)} <= M ||omega_0||_inf.
Omega_t is a COMPLEX strip {x+iy in C^3 : |y| < sqrt(t/c(M))} (verified in the source text),
so it contains R^3 and the bound is a genuine global sup bound on the real vorticity.
Taking M = 3/2 gives  T_d >= 1/(c(3/2) ||omega_0||_inf), with NO Reynolds number.
That contradicts the candidate's upper bound inf T_d <= c2/(M log Re_E) for large Re.

The load-bearing step of their sketch (p.10) is the weighted BMO inequality
      int |f(x)| / (|x|+1)^4 dx  <=  c ||f||_BMO ,          f = grad u,
which they themselves parenthesise as needing the local average subtracted.

PART 1 shows that inequality FAILS, with an absolute constant, on precisely the family this
campaign is built from -- the bang-bang shell -- and fails by a factor log(R/rho0).
PART 2 restores the average term and shows BFG's own continuation estimate then returns the
logarithmic clock, with constant ~2.
"""
import numpy as np, math, json, hashlib, time
from scipy.special import ellipk, ellipe
from scipy.integrate import quad
from scipy.optimize import brentq
t0=time.time(); out={}
rho0=1.0; M=1.0

def ring_ur(r,z,a):
    d2=(a+r)**2+z**2; m=np.clip(4*a*r/d2,0,1-1e-15)
    K=ellipk(m); E=ellipe(m); q=np.where(((a-r)**2+z**2)==0,1e-300,(a-r)**2+z**2)
    return z/(2*np.pi*r*np.sqrt(d2))*(-K+(a*a+r*r+z*z)/q*E)
def ray_circle(r0,z0,ct,st,Rc):
    b=2*(r0*ct+z0*st); c=r0*r0+z0*z0-Rc*Rc; d=b*b-4*c
    if d<0: return []
    sd=math.sqrt(d); return [x for x in ((-b-sd)/2,(-b+sd)/2) if x>0]
def a_shell(r0,z0,R,ntheta=1440,nrad=48):
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
            rho=0.5*(v-u)*gx+0.5*(v+u); w=0.5*(v-u)*gw
            rp=r0+rho*ct; zp=z0+rho*st
            tot+=np.sum(w*rho*(-M*np.sign(zp))*ring_ur(np.full_like(rp,r0),z0-zp,rp))*dth
    return tot/r0

print("PART 1.  The weighted-BMO step, tested on the bang-bang shell.")
print("  In the empty hole |x| < rho0 the field is smooth and, by axisymmetry +")
print("  incompressibility, grad u has eigenvalues (a, a, -2a) with a = u^r/r, so")
print("  |grad u| = sqrt(6) a  (Frobenius).   ||omega||_inf = M = 1 for every R.")
print(f"\n{'R/rho0':>9} {'a(0,0)':>10} {'(M/2)logR':>10} {'a(0.5,0)':>10} {'a(0,0.5)':>10} {'|grad u|(0)':>12}")
rows=[]
for R in [4.0,16.0,64.0,256.0,1024.0]:
    a0=a_shell(1e-5,0.0,R)
    a1=a_shell(0.5,0.0,R)
    a2=a_shell(1e-5,0.5,R)
    print(f"{R:>9g} {a0:>10.6f} {0.5*math.log(R):>10.6f} {a1:>10.6f} {a2:>10.6f} {math.sqrt(6)*a0:>12.6f}")
    rows.append(dict(R=R,a0=a0,a_r=a1,a_z=a2))
out['hole']=rows
# weight mass of the half-hole
wmass=4*math.pi*quad(lambda r: r*r/(1+r)**4,0,0.5)[0]
print(f"\n  weight mass  int_{{|x|<1/2}} dx/(1+|x|)^4 = {wmass:.6f}")
amin=[min(r['a0'],r['a_r'],r['a_z']) for r in rows]
print(f"{'R/rho0':>9} {'min a on |x|<=1/2':>18} {'lower bd for int |gradu|/(1+|x|)^4':>36}")
for r,am in zip(rows,amin):
    print(f"{r['R']:>9g} {am:>18.6f} {math.sqrt(6)*am*wmass:>36.6f}")
out['bmo_lower']=[math.sqrt(6)*am*wmass for am in amin]
sl,ic=np.polyfit([math.log(r['R']) for r in rows],out['bmo_lower'],1)
print(f"\n  the lower bound grows LINEARLY in log(R/rho0):  slope {sl:.5f} per e-fold, intercept {ic:+.5f}")
print("  while ||grad u||_BMO <= C ||omega||_inf = C M is bounded, uniformly in R.")
print("  => int |f|/(|x|+1)^4 <= c ||f||_BMO is FALSE with an absolute constant.  The")
print("     bang-bang shell -- the campaign's own extremiser -- is a counterexample, and the")
print("     deficit is exactly one power of log(R/rho0), i.e. exactly the clock's logarithm.")
out['bmo_slope']=float(sl)

print("""
PART 2.  BFG's own continuation estimate with the local average restored.
Their string, with the honest weighted inequality
    int |f(z)|/(|z|+1)^4 dz  <=  c ( ||f||_BMO + |f_{B(x,L)}| ),   L = sqrt(nu (t-s)),
and the classical Biot-Savart / telescoping bound
    |f_{B(x,L)}| <= C M (1 + log_+(ell/L)),   ell = E^{1/5} M^{-2/5}
(the L >= ell case is exactly Astra pass 8 eq. (6); the L < ell case telescopes ball means
over log2(ell/L) doublings at cost CM each), gives, with the bootstrap sup||omega|| <= 2M,

    increment <= C M^2 int_0^t [ 1 + (1/2) log_+(Re/(M sigma)) ] d sigma ,  Re = M ell^2/nu.

Holding the record (increment <= M/2) requires, in tau = M t,
    tau [ 3/2 + (1/2) log(Re/tau) ] <= 1/(2C).""")
def solve_tau(Re,C=1.0):
    f=lambda tau: tau*(1.5+0.5*math.log(Re/tau))-1.0/(2*C)
    return brentq(f,1e-12,1.0)
print(f"\n{'Re':>10} {'log Re':>8} {'tau = M t':>11} {'tau log Re':>11}")
tab=[]
for Re in [1e2,1e4,1e6,1e10,1e20,1e40,1e80]:
    tau=solve_tau(Re)
    print(f"{Re:>10.0e} {math.log(Re):>8.3f} {tau:>11.6f} {tau*math.log(Re):>11.5f}")
    tab.append(dict(Re=Re,tau=tau,P=tau*math.log(Re)))
out['repaired_clock']=tab
print("""  tau*log Re tends to a constant of order 1-2: the repaired BFG argument IS a
  logarithmic clock, not a log-free one.  With the average term DROPPED (BFG as written)
  the same string gives C M^2 t <= M/2, i.e. tau <= 1/(2C) -- log-free.  The entire
  discrepancy between BFG Thm 8/10 and Astra pass 8 is that one term.""")
print(f"  For scale: the viscous-numerics seat measures T32*M*log Re_E -> ~1.98 on its family,")
print(f"  and the shell first-order estimate gives 4 ln(3/2) = {4*math.log(1.5):.4f}.")

json.dump(out,open('r3_results.json','w'),indent=1)
print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")

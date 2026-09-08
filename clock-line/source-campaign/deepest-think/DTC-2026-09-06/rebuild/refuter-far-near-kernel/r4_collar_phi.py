#!/usr/bin/env python3
"""R8/additive - test the build's find (3): "the 1/sin(phi) in the proved collar bound is REAL,
not an artifact of crude estimation -- the exact c_edge genuinely blows up at the axis".
That argument substitutes c_edge (the INNER-EDGE plateau offset, which lives in the l>=3 INTERIOR
modes) for the COLLAR piece.  Here I evaluate the exact collar piece itself as phi -> 0, in two
geometries, with my own code:
  (i) deep interior rho >> rho0 (the two-sided collar rho/2 < rho' < 2 rho) -- the object the
      lemma's C_tot(phi) actually bounds;
  (ii) at the inner edge rho = rho0+ (one-sided collar, u in (1/2,1) interior modes only) -- the
      object that actually carries c_edge's divergence.
Bang-bang cap, M=1, rho0=1, R=4096.
"""
import numpy as np, math, json, hashlib
def geg(t,L):
    t=np.atleast_1d(np.asarray(t,float)); C=np.empty((L+1,t.size)); C[0]=1.0
    if L>=1: C[1]=3.0*t
    for l in range(2,L+1): C[l]=((2*l+1)*t*C[l-1]-(l+1)*C[l-2])/l
    return C
Nl=lambda l:(l+1)*(l+2)/(l+1.5)
def g_bang(L,nq=6000):
    x,w=np.polynomial.legendre.leggauss(nq); th=0.5*math.pi*(x+1); w=0.5*math.pi*w
    C=geg(np.cos(th),L); I=C@(w*np.sin(th)**2*(-np.sign(np.cos(th))))
    G=I/Nl(np.arange(L+1)); G[np.arange(L+1)%2==0]=0.0; return G
L=2001; g=g_bang(L); ls=np.arange(1,L+1,2)
def ces(t):
    P=np.cumsum(t); return float(P[-max(4,len(P)//4):].mean())
def collar(rho,phi,rho0=1.0,R=4096.0):
    """exact collar piece: shells rho/2 < rho' < 2 rho, i.e. u = rho/rho' in (1/2,2),
    interior modes on (1/2, min(1,rho/rho0)), exterior modes on (max(1,rho/R), min(2,rho/rho0))."""
    t=math.cos(phi); C=geg(np.array([t]),L+1)[:,0]
    cin=-((ls+2)*g[ls]/(2*ls+3))*C[ls-1]; cout=((ls+1)*g[ls]/(2*ls+3))*C[ls+1]
    eps,U=rho/R,rho/rho0
    hi=min(1.0,U); m=np.maximum(ls-1,1)
    Iin=np.where(ls==1,math.log(hi/0.5),(hi**m-0.5**m)/m); Iin=np.where(hi>0.5,Iin,0.0)
    lo_o,hi_o=max(1.0,eps),min(2.0,U)
    Iout=np.where(hi_o>lo_o,(lo_o**(-(ls+4))-hi_o**(-(ls+4)))/(ls+4),0.0)
    return ces(cin*Iin)+ces(cout*Iout), ces(cin*Iin), ces(cout*Iout)
def far_rem(rho,phi,rho0=1.0,R=4096.0):
    t=math.cos(phi); C=geg(np.array([t]),L+1)[:,0]
    cin=-((ls+2)*g[ls]/(2*ls+3))*C[ls-1]; eps=rho/R; m=np.maximum(ls-1,1)
    I=np.where(ls==1,math.log(0.5/eps),(0.5**m-eps**m)/m)
    return ces(cin*I)-float(cin[0]*I[0])
print("=== (i) DEEP INTERIOR rho=64: exact two-sided collar and far remainder vs phi ===")
print(f"   {'phi(deg)':>9} {'collar_exact':>13} {'bound 4/s+pi/8':>15} {'slack':>8} {'|far-far0|':>11} {'C1=0.292':>9}")
rows=[]
for phid in (45,10,3,1,0.3,0.1,0.03):
    phi=math.radians(phid); s=math.sin(phi)
    c,_,_=collar(64.0,phi); f=far_rem(64.0,phi); b=4.0/s+math.pi/8
    rows.append((phid,c,b,f)); print(f"   {phid:9.2f} {c:13.6f} {b:15.2f} {b/abs(c):8.0f}x {abs(f):11.6f} {0.292:9.3f}")
print("=== (ii) INNER EDGE rho=rho0+: one-sided collar (interior modes u in (1/2,1)) vs phi ===")
print(f"   {'phi(deg)':>9} {'collar_in':>12} {'far_rem':>11} {'sum=c_edge-ish':>15} {'(1/4)ln(1/phi)':>15}")
rows2=[]
for phid in (10,3,1,0.3,0.1,0.03):
    phi=math.radians(phid)
    c,ci,co=collar(1.0000001,phi); f=far_rem(1.0000001,phi)
    rows2.append((phid,ci,f)); print(f"   {phid:9.2f} {ci:12.6f} {f:11.6f} {ci+f:15.6f} {0.25*math.log(1/phi):15.6f}")
A=np.polyfit([-math.log(math.radians(p)) for p,_,_ in rows2[1:]],[ci for _,ci,_ in rows2[1:]],1)[0]
print(f"   fitted growth of the inner-edge collar: {A:.5f} * log(1/phi)   (the BOUND grows like 4/phi)")
json.dump(dict(deep=rows,edge=rows2,edge_lograte=float(A)),open('r4_results.json','w'),indent=1)
print('SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())

#!/usr/bin/env python3
"""R7 - is D = 124.3991 (s4_multipole_decay.py) a convergent constant, or a truncation artifact?

s4 computes
   Dall(LMAX) = SUM_{l odd <= LMAX} |(l+1) g_l/(2l+3)| * ((2^{l+4}-1)/(l+4)) * ((l+2)(l+3)/2) * 0.5^{l-1}
and reports it at LMAX = 120 as "all-l absolute bound valid for rho >= 4 rho_in (|xi|>=2)".
The factor 2^{l+4} * 0.5^{l-1} = 32 is l-INDEPENDENT, so the summand ~ 8 l |g_l|:  the sum converges
only if |g_l| = o(l^-2).  Tests:
  (1) |g_l| decay rate for the bang-bang shell.
  (2) Dall(LMAX) for LMAX = 40..4000.
  (3) The CORRECT l1 constant at the stated validity range rho >= 4 rho_in: factor 0.25^{l-1}.
  (4) A datum-INDEPENDENT Cauchy-Schwarz constant for ANY z-odd |omega^theta| <= M_in:
      |a_octave(x)| <= sqrt(2) M_in INT_{v- }^{v+} S(v) dv/v with v-=rho_in/rho, v+=2 rho_in/rho;
      expressed as D_CS * M_in * (rho_in/rho)^5.
  (5) direct check of both candidate constants against the exactly-evaluated octave field.
"""
import numpy as np, math, json, hashlib
def geg(t, L):
    t=np.atleast_1d(np.asarray(t,float)); C=np.empty((L+1,t.size)); C[0]=1.0
    if L>=1: C[1]=3.0*t
    for l in range(2,L+1): C[l]=((2*l+1)*t*C[l-1]-(l+1)*C[l-2])/l
    return C
Nl=lambda l:(l+1)*(l+2)/(l+1.5)
def g_bang(L,nq=8000):
    x,w=np.polynomial.legendre.leggauss(nq); th=0.5*math.pi*(x+1); w=0.5*math.pi*w
    C=geg(np.cos(th),L); I=C@(w*np.sin(th)**2*(-np.sign(np.cos(th))))
    return I/Nl(np.arange(L+1))

LBIG=4000
g=g_bang(LBIG)
ls=np.arange(1,LBIG+1,2)
print("=== (1) decay of |g_l| (bang-bang shell) ===")
for l in (11,51,101,501,1001,2001,3001):
    print(f"   l={l:5d}  |g_l| = {abs(g[l]):.6e}   l^1.5 |g_l| = {l**1.5*abs(g[l]):.5f}   l^2 |g_l| = {l**2*abs(g[l]):.4f}")
sl=np.polyfit(np.log(ls[200:1200]),np.log(np.abs(g[ls[200:1200]])),1)[0]
print(f"   fitted power-law exponent of |g_l| over l=400..2400 : {sl:.4f}   (l^-1.5 predicted)")

def Dall(LMAX, halving):
    L=np.arange(1,LMAX+1,2)
    # overflow-free:  (2^{L+4}-1) * h^{L-1} = 2^{L+4} h^{L-1} - h^{L-1}
    pw = np.exp2((L+4) + (L-1)*np.log2(halving)) - halving**(L-1.0)
    return float(np.sum(np.abs((L+1)*g[L]/(2*L+3))*(pw/(L+4))*((L+2)*(L+3)/2)))
print("\n=== (2) s4's Dall with factor 0.5^(l-1)  [the number reported as 124.3991] ===")
prev=None
for LM in (40,80,120,240,500,1000,2000,4000):
    v=Dall(LM,0.5); print(f"   LMAX={LM:5d}   Dall = {v:10.3f}   Dall/sqrt(LMAX) = {v/math.sqrt(LM):.4f}")
    prev=v
D120=Dall(120,0.5); D2000=Dall(2000,0.5)
fired = D2000/D120 > 1.5
print(f"   ratio Dall(2000)/Dall(120) = {D2000/D120:.3f}   ->  R7 {'FIRED' if fired else 'did not fire'}  (criterion >1.5)")

print("\n=== (3) the CORRECT l1 constant at the stated validity range rho >= 4 rho_in (factor 0.25^(l-1)) ===")
for LM in (40,120,500,2000,4000):
    print(f"   LMAX={LM:5d}   D_l1(rho>=4rho_in) = {Dall(LM,0.25):10.5f}")
D_l1 = Dall(4000,0.25)

print("\n=== (4) datum-independent Cauchy-Schwarz constant (valid for ANY z-odd |omega^theta|<=M_in) ===")
def S(v,odd=True,L=3000):
    L_=np.arange(1,L,2) if odd else np.arange(0,L)
    c=((L_+1)/(2*L_+3))**2*(((L_+2)*(L_+3)/2)**2)/Nl(L_)
    with np.errstate(under='ignore'): return math.sqrt(float(np.sum(c*v**(2.0*(L_+4)))))
from scipy.integrate import quad
for j in (2,3,4,5,6):
    vm,vp = 2.0**-j, 2.0**(1-j)
    I,_ = quad(lambda v: S(v)/v, vm, vp, limit=200)
    D_cs = math.sqrt(2)*I/ (2.0**(-5*j))
    print(f"   j={j}  |x|={2**j} rho_in :  sqrt(2) INT S dv/v = {math.sqrt(2)*I:.4e}   => D_CS = {D_cs:.4f}")
vm,vp=2.0**-6,2.0**-5
I,_=quad(lambda v: S(v)/v, vm, vp, limit=200); D_CS=math.sqrt(2)*I/(2.0**-30)

print("\n=== (5) both candidate constants vs the exact octave field (bang-bang z-odd octave 1<rho'<2) ===")
def a_oct(rho,phi,L=2000):
    C=geg(np.array([math.cos(phi)]),L+1)[:,0]; L_=np.arange(0,L+1)
    # overflow-free: (2^{L+4}-1)*rho^{-(L+4)} = (2/rho)^{L+4} - rho^{-(L+4)}
    pw = (2.0/rho)**(L_+4.0) - (1.0/rho)**(L_+4.0)
    return float(np.sum(((L_+1)*g[L_]/(2*L_+3))*(pw/(L_+4))*C[L_+1]))
print(f"   {'rho':>7} {'phi':>6} {'|a|/(rho_in/rho)^5':>20}")
worst=0.0
for j in (2,3,4,5):
    for phid in (5,40,90):
        rho=2.0**j; v=abs(a_oct(rho,math.radians(phid)))/(1.0/rho)**5
        worst=max(worst,v); print(f"   {rho:7.1f} {phid:6d} {v:20.5f}")
print(f"   worst measured D_measured = {worst:.4f}")
print(f"\n   reported  D(s4, 0.5^(l-1), LMAX=120) = {D120:.4f}   <-- truncation-dependent")
print(f"   correct   D_l1  (0.25^(l-1), rho>=4rho_in) = {D_l1:.4f}")
print(f"   universal D_CS  (Cauchy-Schwarz, any z-odd datum) = {D_CS:.4f}")
print(f"   exact     D_measured (bang-bang) = {worst:.4f}")
json.dump(dict(D120=D120,D2000=D2000,D_l1=D_l1,D_CS=D_CS,D_measured=worst,g_exponent=sl,fired=bool(fired)),
          open('r3_results.json','w'),indent=1)
print('SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())

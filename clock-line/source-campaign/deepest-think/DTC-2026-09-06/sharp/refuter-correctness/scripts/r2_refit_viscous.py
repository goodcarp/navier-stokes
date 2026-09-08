"""
R2. Independent re-extraction of T32 from the viscous-numerics seat's RAW time series,
and a model-selection test of the load-bearing claim "T32*M falls like 1/log Re_E".

The seat's own summary reports T32 in each run record.  Here T32 is recomputed from
hist['t'], hist['ommax'] by linear interpolation of the first upcrossing of 1.5, so the
extraction is not taken on trust.  Then four models are fitted to (log Re_E, T32*M):

   A   log law with offset       T32*M = c / (logRe - d)          -> 0
   B   saturating (BFG-shaped)   T32*M = A + c / (logRe - d)      -> A >= 0
   C   power law in Re           T32*M = c * Re^-p
   D   pure log power            T32*M = c / (logRe)^q

The question the refuter needs answered is whether the seat's 6 points can DISTINGUISH
A (log clock sharp) from B with A>0 (BFG Thm 10: T32*M bounded below, no log).
"""
import json, os, math, hashlib, time
import numpy as np
from scipy.optimize import curve_fit

t0=time.time()
SRC="~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/sharp/viscous-numerics"
out={}

def load(tag):
    fn=os.path.join(SRC,"results_%s.json"%tag)
    return json.load(open(fn))["runs"] if os.path.exists(fn) else []

def first_cross(t,y,lev):
    for i in range(1,len(y)):
        if y[i-1]<lev<=y[i]:
            return t[i-1]+(lev-y[i-1])*(t[i]-t[i-1])/(y[i]-y[i-1])
    return None

runs=load("main")+load("n7")
runs=[r for r in runs if r["kind"]=="shell" and r.get("label","")=="main" or r.get("label","")=="n7"]
runs=sorted({r["N"]:r for r in runs}.items())
print("RE-EXTRACTION OF T32 FROM RAW HISTORIES (primary datum A, Re0=100)")
print(f"{'N':>2} {'R':>5} {'logReE':>8} {'T32 stored':>11} {'T32 mine':>10} {'rel diff':>9} {'om0':>7} {'max om':>8}")
rows=[]
for N,r in runs:
    h=r["hist"]; t=h["t"]; om=h["ommax"]
    mine=first_cross(t,om,1.5*om[0])
    st=r["T32"]
    rel=(abs(mine-st)/st if (mine and st) else float('nan'))
    print(f"{N:>2} {r['R']:>5g} {r['logReE']:>8.3f} "
          f"{(('%11.5f'%st) if st else '       none')} "
          f"{(('%10.5f'%mine) if mine else '      none')} "
          f"{rel:>9.2e} {om[0]:>7.4f} {max(om):>8.4f}")
    if mine: rows.append((N,r["logReE"],mine,r["R"],r["s_star_at_T32"]))
out['reextract']=[dict(N=a,logRe=b,T32=c) for a,b,c,_,_ in rows]

N=np.array([a for a,_,_,_,_ in rows],float)
x=np.array([b for _,b,_,_,_ in rows])      # log Re_E
y=np.array([c for _,_,c,_,_ in rows])      # T32 * M   (M = 1)
print("\nrecomputed  logRe:",np.round(x,4).tolist())
print("recomputed  T32*M:",np.round(y,6).tolist())
print("1/(T32*M)        :",np.round(1/y,5).tolist())
print("successive diffs of 1/(T32*M) (linear in N would be constant):",
      np.round(np.diff(1/y),5).tolist())

def rms(f,p): return float(np.sqrt(np.mean((f(x,*p)-y)**2)))
def rel(f,p): return float(np.max(np.abs(f(x,*p)/y-1)))

mA=lambda z,c,d: c/(z-d)
mB=lambda z,A,c,d: A+c/(z-d)
mC=lambda z,c,p: c*np.exp(-p*z)
mD=lambda z,c,q: c/z**q

pA,_=curve_fit(mA,x,y,p0=[2,5],maxfev=200000)
pC,_=curve_fit(mC,x,y,p0=[10,0.2],maxfev=200000)
pD,_=curve_fit(mD,x,y,p0=[100,2],maxfev=200000)
try:
    pB,_=curve_fit(mB,x,y,p0=[0.0,pA[0],pA[1]],maxfev=400000)
except Exception as e:
    pB=None; print("model B fit failed:",e)

print("\nMODEL FITS (6 points)")
print(f"  A  T32*M = c/(logRe-d)      c={pA[0]:.4f} d={pA[1]:.4f}      rms={rms(mA,pA):.3e}  maxrel={rel(mA,pA):.2%}")
if pB is not None:
    print(f"  B  T32*M = A+c/(logRe-d)    A={pB[0]:+.5f} c={pB[1]:.4f} d={pB[2]:.4f}  rms={rms(mB,pB):.3e}  maxrel={rel(mB,pB):.2%}")
print(f"  C  T32*M = c Re^-p          c={pC[0]:.4f} p={pC[1]:.5f}     rms={rms(mC,pC):.3e}  maxrel={rel(mC,pC):.2%}")
print(f"  D  T32*M = c/(logRe)^q      c={pD[0]:.4f} q={pD[1]:.4f}     rms={rms(mD,pD):.3e}  maxrel={rel(mD,pD):.2%}")
out['fits']=dict(A=dict(p=pA.tolist(),rms=rms(mA,pA)),
                 C=dict(p=pC.tolist(),rms=rms(mC,pC)),
                 D=dict(p=pD.tolist(),rms=rms(mD,pD)))
if pB is not None: out['fits']['B']=dict(p=pB.tolist(),rms=rms(mB,pB))

# how large an asymptote A can the data tolerate?  profile A over a grid, refitting c,d.
print("\nPROFILE: force a nonzero floor A (BFG-shaped) and refit c,d.")
print(f"{'A (floor)':>10} {'c':>9} {'d':>9} {'rms':>10} {'rms/rmsA':>9}")
prof=[]
base=rms(mA,pA)
for A0 in [0.0,0.02,0.05,0.10,0.15,0.20,0.2272]:
    try:
        f=lambda z,c,d: A0+c/(z-d)
        p,_=curve_fit(f,x,y,p0=[pA[0],pA[1]],maxfev=400000)
        rr=float(np.sqrt(np.mean((f(x,*p)-y)**2)))
        print(f"{A0:>10.4f} {p[0]:>9.4f} {p[1]:>9.4f} {rr:>10.3e} {rr/base:>9.2f}")
        prof.append(dict(A=A0,c=p[0],d=p[1],rms=rr,ratio=rr/base))
    except Exception as e:
        print(f"{A0:>10.4f}  fit failed {e}")
out['floor_profile']=prof
print("""
READ:  a floor as large as the LAST measured value (A = T32*M at N=7) can be forced only
at a large cost in rms; but a floor of a few 0.01 is nearly free.  Six points spanning a
factor 2.07 in log Re cannot exclude a positive asymptote of order 0.05 M^-1, i.e. cannot
falsify a log-free clock whose constant is c(3/2) >= 20.  The numerics BOUND the constant
in the conjectured upper bound; they do not establish that T32*M -> 0.""")

# extrapolate the log model to the reach needed to contradict a log-free clock
c,d=pA
for C in [10,20,50,100]:
    # need T32*M <= 1/C
    z=c*C+d
    print(f"  to reach T32*M = 1/{C:<4d} the log model needs log Re_E = {z:8.2f}  (Re_E = 10^{z/math.log(10):.1f}), "
          f"octaves N = {(z-3.688)/1.3863:6.1f}")

# which ring carries the max, from the raw record
print("\nLOCATION OF THE MAX AT T32 (octave index log2 s*, rho0 = 1):")
for Nn,xx,yy,R,ss in rows:
    print(f"  N={Nn}  R={R:>5g}  s*={ss:.4f}  octave={math.log2(ss):.3f}  log(R/s*)={math.log(R/ss):.4f}  Q'={yy*math.log(R/ss):.4f}")
Qp=[yy*math.log(R/ss) for _,_,yy,R,ss in rows]
print(f"  Q' mean {np.mean(Qp):.4f}  sd {np.std(Qp):.4f}  (seat reports 0.988 +- 0.018)")
out['Qprime']=dict(vals=Qp,mean=float(np.mean(Qp)),sd=float(np.std(Qp)))

json.dump(out,open('r2_results.json','w'),indent=1)
print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")

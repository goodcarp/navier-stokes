"""
R6.  (a) Is the viscous-numerics PREREG gate diagnostic for the claim it tests?
     (b) Does it matter whether Re_E uses the INITIAL energy or the running one?
"""
import json,math,os,hashlib
import numpy as np
from scipy.optimize import curve_fit
out={}
SRC="~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/sharp/viscous-numerics"
runs=[]
for tag in ("main","n7"):
    runs+=json.load(open(os.path.join(SRC,"results_%s.json"%tag)))["runs"]
ok=sorted({r["N"]:r for r in runs if r["T32"] and r.get("label") in ("main","n7")}.items())
N=np.array([k for k,_ in ok],float); T=np.array([v["T32"] for _,v in ok]); X=np.array([v["logReE"] for _,v in ok])
def beta(y):
    p=np.polyfit(np.log(N),np.log(y),1); return -p[0]
b_data=beta(T)
mA=lambda z,c,d: c/(z-d)
pA,_=curve_fit(mA,X,T,p0=[2,5],maxfev=200000)
print("(a) PREREG gate: fit log(T32*M) = const - beta log N ; survive if beta>=0.7, kill if beta<=0.3")
print(f"    beta from the DATA                                    = {b_data:.4f}   (seat reports 1.4059)")
print(f"    beta from model A  T32*M = {pA[0]:.4f}/(logRe-{pA[1]:.4f})   -> 0  = {beta(mA(X,*pA)):.4f}")
rows=[]
for A0 in [0.02,0.05,0.10,0.15]:
    f=lambda z,c,d: A0+c/(z-d)
    p,_=curve_fit(f,X,T,p0=[pA[0],pA[1]],maxfev=400000)
    yb=f(X,*p); bb=beta(yb); mr=float(np.max(np.abs(yb/T-1)))
    print(f"    beta from a SATURATING model with floor A={A0:.2f} (T32*M -> {A0:.2f} > 0) = {bb:.4f}"
          f"   [max rel misfit to the data {mr:.2%}]")
    rows.append(dict(A=A0,beta=bb,maxrel=mr))
out['beta']=dict(data=float(b_data),modelA=float(beta(mA(X,*pA))),saturating=rows)
print("""    Every one of these models fits the six measured points to a few percent, and every
    one of them yields beta well above the 0.7 'survive' line -- including the models whose
    T32*M is BOUNDED BELOW, i.e. the exact negation of the hypothesis under test.  The
    pre-registered beta gate is therefore not diagnostic for 'the log is necessary'; it
    discriminates 'T32 falls appreciably with N' from 'it does not'.  Passing it is real
    evidence that the mechanism is not stopped by viscosity; it is not evidence that the
    clock's logarithm is sharp.""")

print("\n(b) does Re_E drift over [0, T_d] enough to matter?")
print("    Re_E = E^{2/5} M^{1/5} / nu ; E is non-increasing, M grows by at most 3/2 up to T_d.")
f=1.5**0.2
print(f"    so Re_E(t) <= Re_E(0) * (3/2)^(1/5) = {f:.6f} Re_E(0) on [0,T_d],")
print(f"    i.e. log Re_E moves by at most {math.log(f):+.6f} -- an additive O(0.08), absorbed in the constants.")
print("    Measured drift in the seat's runs (E is recorded at t=0 only, so this is the bound, not a measurement).")
out['ReE_drift_log']=float(math.log(f))
json.dump(out,open('r6_results.json','w'),indent=1)
print("\nSHA256",hashlib.sha256(open(__file__,'rb').read()).hexdigest())

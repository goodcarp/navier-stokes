"""
v5_application.py -- gap-V-aronson.  What THEOREM V.1 costs in the lower-bound argument.

Geometry (prove-lagrangian NOTE section 1 and 4(3)):
   rho0 * delta = sqrt(nu/M)   (finest datum feature = dissipation length)
   L = log(R/rho0),  a0 = (M/2) L,  window tau = 2 theta /(M L),  theta = int_0^tau a dt
   grad_5 b for the uniform axisymmetric strain = diag(a,a,a,a,-2a) => Gamma = 2 a0 = M L
   => c = Gamma tau = 2 theta.   theta = log(3/2) (frozen clock) or 2(1-sqrt(2/3))/kappa (accelerated).
   nu tau = (rho0 delta)^2 c / L    =>   q := d^2/(nu tau) = (d/(rho0 delta))^2 * L/c.

Tracked point: rho_* = rho0 (1+f), polar angle phi0.  Everything in units rho0 = 1.
d = min( inner-edge distance f, taper distance (1+f) sin(phi0-delta),
         equatorial-layer distance (1+f) sin(90deg - eps_eq - phi0) ).
Amplitude ratio A = ||eta_0 - eta~_0||_inf / |eta_0(x_*)| <= 2 (1+f) sin(phi0)/sin(delta).
"""
import json, numpy as np
HERE = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/gap-V-aronson"
n = 5; OUT = {}
_ES = np.linspace(1e-3, 1.4, 4001)
def B1(c,q,n=5): return 2*n*np.exp(-q*np.exp(-2*c)/(4*n))
def B2(c,q,n=5):
    k = c/(2*np.exp(2*c)*(np.exp(2*c)-1))
    lv = n*np.log(1+2/_ES)-(1-_ES**2/2)**2*k*q; return float(np.exp(lv.min()))
def BOUND(c,q,n=5): return min(B1(c,q,n), B2(c,q,n))
def BOUND_L(c,q,n=5):   # conditional tier (needs a grad^2 u bound); rate e^{-2c}/4
    k = np.exp(-2*c)/4
    lv = n*np.log(1+2/_ES)-(1-_ES**2/2)**2*k*q; return float(np.exp(lv.min()))
def SHARP(c,q): return float(np.exp(-q*c/(2*(np.exp(2*c)-1))))

theta_frozen = np.log(1.5); theta_acc = 4*(1-np.sqrt(2/3))
c_frozen = 2*theta_frozen; c_acc = theta_acc
print("windows:  theta_frozen = %.6f -> c = %.6f ;  theta_accel = %.6f -> c = %.6f"
      % (theta_frozen, c_frozen, theta_acc, c_acc))
OUT['theta_frozen']=theta_frozen; OUT['c_frozen']=c_frozen
OUT['theta_accel']=theta_acc;     OUT['c_accel']=c_acc
# check the accelerated window really integrates to log(3/2)
kap=0.5; chk = -2*np.log(1-kap*theta_acc/2)
print("check: F(0,theta_accel) = -2 log(1 - kappa theta/2) = %.8f  vs log(3/2) = %.8f"
      % (chk, np.log(1.5)))
OUT['accel_window_check']=[float(chk), float(np.log(1.5))]

# ---- the brief's proposed d = one dissipation length (f such that d = rho0*delta) ----
deg = np.pi/180
print("\n[1] the brief's proposal  d = rho0*delta = sqrt(nu/M)   =>  q = L/c  exactly.")
rows=[]
for delta_deg in (7.5,):
    delta = delta_deg*deg
    for c in (c_acc, c_frozen):
        for L in (5,10,20,40,80,160,320,640):
            q = L/c
            A = 2*np.sin(30*deg)/np.sin(delta)   # phi0=30deg, f~delta so (1+f)~1
            b = BOUND(c,q); print("    delta=%.1f c=%.3f L=%4d  q=%7.1f  A=%.2f  A*BOUND=%.3e   (target 1/L=%.3f)"
                                  % (delta_deg,c,L,q,A,A*b,1.0/L))
            rows.append(dict(delta_deg=delta_deg,c=c,L=L,q=q,A=A,loss=A*b,target=1.0/L,ok=bool(A*b<=1.0/L)))
OUT['brief_d_one_dissipation_length']=rows
Lmin = {}
for c in (c_acc, c_frozen):
    Lm=None
    for L in range(5,4001):
        if 2*np.sin(30*deg)/np.sin(7.5*deg)*BOUND(c,L/c) <= 1.0/L: Lm=L; break
    Lmin[str(round(c,5))]=Lm
    print("    -> smallest L for which d = rho0*delta suffices (loss <= 1/L): L = %s  (log Re ~ %s)"
          % (Lm, None if Lm is None else 2*Lm))
OUT['brief_Lmin']=Lmin

# ---- required inset f ----
def feature_d(f, phi0, delta, eps_eq):
    return min(f, (1+f)*np.sin(phi0-delta), (1+f)*np.sin(np.pi/2-eps_eq-phi0))
def loss(f, L, c, delta, phi0, eps_eq, tier=BOUND):
    d = feature_d(f, phi0, delta, eps_eq)
    if d <= 0: return np.inf, 0.0
    q = (d/delta)**2 * L/c
    A = 2*(1+f)*np.sin(phi0)/np.sin(delta)
    return A*tier(c,q), d
def min_f(L, c, delta, phi0, eps_eq, target, tier=BOUND):
    lo, hi = 1e-4, 20.0
    if loss(hi,L,c,delta,phi0,eps_eq,tier)[0] > target: return None
    for _ in range(80):
        mid=0.5*(lo+hi)
        if loss(mid,L,c,delta,phi0,eps_eq,tier)[0] <= target: hi=mid
        else: lo=mid
    return hi

delta = 7.5*deg; eps_eq = 5*deg
print("\n[2] required inset f = (rho_* - rho0)/rho0 so that the edge loss <= 1/L")
print("    delta=7.5deg, equatorial layer half-width 5deg.  Tier: PROVED min(B1,B2).")
print("    %-6s %-6s | %-7s %-7s %-7s | %-9s %-9s" % ("L","phi0","f_req","d/rho0d","q","c2 infl.","f(V.2)"))
rows2=[]
for c,cname in ((c_acc,'accel'),(c_frozen,'frozen')):
    print("   --- c = %.5f (%s window) ---" % (c,cname))
    for L in (5,10,20,40,80,160):
        for phi0d in (30.,45.,60.):
            phi0=phi0d*deg
            f = min_f(L,c,delta,phi0,eps_eq,1.0/L)
            f2 = min_f(L,c,delta,phi0,eps_eq,1.0/L,tier=BOUND_L)
            if f is None:
                print("    %-6d %-6.0f | %-7s" % (L,phi0d,"none")); 
                rows2.append(dict(c=c,L=L,phi0=phi0d,f=None)); continue
            d = feature_d(f,phi0,delta,eps_eq); q=(d/delta)**2*L/c
            infl = 1.0/(1-np.log(1+f)/L)
            print("    %-6d %-6.0f | %-7.4f %-7.2f %-7.1f | %-9.4f %-9s"
                  % (L,phi0d,f,d/delta,q,infl,("%.4f"%f2) if f2 else "none"))
            rows2.append(dict(c=c,L=L,phi0=phi0d,f=float(f),d_over_rho0delta=float(d/delta),
                              q=float(q),c2_inflation=float(infl),f_tierV2=None if f2 is None else float(f2)))
OUT['required_inset']=rows2

# ---- c2 numbers ----
c2_acc = 8*(1-np.sqrt(2/3)); c2_frozen = 4*np.log(1.5)
print("\n[3] effect on the constant.  c2(accel)=%.7f  c2(frozen)=%.7f" % (c2_acc, c2_frozen))
best = {}
for L in (10,20,40,80,160):
    r=[x for x in rows2 if x['c']==c_acc and x['L']==L and x['f'] is not None]
    if not r: continue
    x=min(r,key=lambda z:z['f'])
    best[str(L)] = dict(phi0=x['phi0'], f=x['f'], infl=x['c2_inflation'],
                        c2=c2_acc*x['c2_inflation'])
    print("    L=%4d  best phi0=%2.0fdeg  f=%.4f  c2 -> %.5f  (asymptotic %.7f)"
          % (L,x['phi0'],x['f'],c2_acc*x['c2_inflation'],c2_acc))
OUT['c2_effect']=best; OUT['c2_accel']=c2_acc; OUT['c2_frozen']=c2_frozen

# ---- V-b: the bulk loss the Gaussian lemma does NOT cover ----
print("\n[4] GAP V-b (the bulk viscous loss on the smooth plateau) -- NOT covered by V.1")
print("    first-moment (grad u only):  rel <= [r0/r_min^2] * delta * sqrt(5(e^{2c}-1)/L)   = O(L^{-1/2})")
print("    second-moment (needs grad^2 u): rel ~ 5 delta^2 (e^{2c}-1)/(L r0^2)              = O(L^{-1})")
rows4=[]
for c in (c_acc,):
    for L in (10,20,40,80,160,320,640):
        phi0=45*deg
        f = min_f(L,c,delta,phi0,eps_eq,1.0/L)
        if f is None: continue
        d = feature_d(f,phi0,delta,eps_eq)
        r0=(1+f)*np.sin(phi0); rmin=r0-d
        b1 = (r0/rmin**2)*delta*np.sqrt(5*(np.exp(2*c)-1)/L)
        b2 = 5*delta**2*(np.exp(2*c)-1)/(L*r0**2)
        exact0 = c*delta**2/(L*np.sin(phi0)**2)   # the note's instantaneous nu/r^2 * tau
        rows4.append(dict(L=L,f=float(f),first_moment=float(b1),second_moment=float(b2),
                          note_nu_over_r2=float(exact0)))
        print("      L=%4d f=%.3f | first-moment bound %.4f | second-moment %.5f | note's nu tau/r^2 %.5f"
              % (L,f,b1,b2,exact0))
OUT['Vb']=rows4

# ---- why the Eulerian route is unavailable ----
print("\n[5] why the Eulerian (drift-in-the-coefficients) route fails quantitatively")
for L in (10,40,160):
    lam=1.5; disp=(lam-1)   # |X(tau)-x0|/rho0 >= lambda-1 for the innermost shell (rho -> lambda rho)
    for c in (c_acc,):
        sq=delta*np.sqrt(c/L)
        print("    L=%4d : drift displacement/rho0 = %.3f ; sqrt(nu tau)/rho0 = %.5f ; ratio = %.1f sigma"
              % (L,disp,sq,disp/sq))
        OUT.setdefault('eulerian',[]).append(dict(L=L,disp=disp,sqrt_nutau=float(sq),sigma=float(disp/sq)))

with open(HERE+"/v5_results.json","w") as fh: json.dump(OUT,fh,indent=1,default=float)
print("\nwrote v5_results.json")

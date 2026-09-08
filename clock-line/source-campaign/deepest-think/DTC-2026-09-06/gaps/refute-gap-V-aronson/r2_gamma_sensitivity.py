#!/usr/bin/env python3
"""
REFUTER r2 -- how much of gap-V-aronson section 5 depends on the UNPROVEN value of Gamma.

Theorem V.1 needs  Gamma := sup over ALL of R^5 x [0,tau] of ||grad_5 b||_op.
Section 5 supplies it as  Gamma = 2 a(0,t),  the ON-AXIS strain, justified only by
the identity ||grad_5 b||_op = 2a for a *uniform* axisymmetric strain.  Neither that
note nor prove-lagrangian bounds the global sup of ||grad u||_op for the actual
(mollified shell) field: prove-lagrangian's L3 bounds the VELOCITY deviation
(sup|a-A(rho)| = 0.7047 M, sup|u^z+2Az|/(M rho) = 0.19143) and a mode sum on psi,
not any derivative of the deviation velocity.

So write Gamma = g * (2 a(0,t)) with an unknown multiplier g >= 1 and sweep g.
tau, and hence nu tau and q, are unchanged (they are physical); only c = Gamma tau moves.
"""
import json, numpy as np
OUT={}; deg=np.pi/180.
kappa=0.5; theta_acc=4*(1-np.sqrt(2/3))
H0=lambda th: 1.0/(1-kappa*th/2)
Gam_max_over_ML = H0(theta_acc)          # sqrt(3/2), from prove-lagrangian's ODE
EPS=np.linspace(1e-4,np.sqrt(2)-1e-4,4000)
def rate1(c,n=5): return np.exp(-2*c)/(4*n)
def rate2(c):     return c/(2*np.exp(2*c)*(np.exp(2*c)-1))
def bound(c,q,n=5):
    b1=2*n*np.exp(-rate1(c,n)*q)
    b2=float(np.exp(np.min(n*np.log1p(2/EPS)-(1-EPS**2/2)**2*rate2(c)*q)))
    return min(b1,b2)
delta=7.5*deg; eps_eq=5*deg
def feat(f,phi0): return min(f,(1+f)*np.sin(phi0-delta),(1+f)*np.sin(np.pi/2-eps_eq-phi0))
def loss(f,L,phi0,c):
    d=feat(f,phi0)
    if d<=0: return np.inf
    return 2*(1+f)*np.sin(phi0)/np.sin(delta)*bound(c,(d/delta)**2*L/theta_acc)
def min_f(L,phi0,c,t):
    lo,hi=1e-4,200.
    if loss(hi,L,phi0,c)>t: return None
    for _ in range(80):
        m=.5*(lo+hi)
        if loss(m,L,phi0,c)<=t: hi=m
        else: lo=m
    return hi
def best_f(L,c,t):
    b=None
    for p in (20.,25.,30.,35.,40.,45.,50.,55.,60.,65.):
        f=min_f(L,p*deg,c,t)
        if f is not None and (b is None or f<b[0]): b=(f,p)
    return b

print("Gamma = g * 2a(0,t);  c = Gamma tau = g * sqrt(3/2) * theta_acc = g * %.6f" % (Gam_max_over_ML*theta_acc))
print("(g = 1 is the theorem applied honestly to the pure-strain part;")
print(" g = %.4f reproduces the value the note actually uses, c = theta_acc.)"
      % (theta_acc/(Gam_max_over_ML*theta_acc)))
print()
print(" %-6s %-9s | %-24s | %-24s | %-24s" % ("g","c","L=10  f (c2)","L=40  f (c2)","L=160 f (c2)"))
c2=8*(1-np.sqrt(2/3)); rows=[]
for g in (0.8165,1.0,1.25,1.5,2.0,3.0,4.0):
    c=g*Gam_max_over_ML*theta_acc
    cells=[]
    for L in (10,40,160):
        b=best_f(L,c,1.0/L)
        if b is None: cells.append(None); continue
        f,p=b; infl=1.0/(1-np.log(1+f)/L); cells.append(dict(f=float(f),phi0=p,c2=float(c2*infl)))
    rows.append(dict(g=g,c=float(c),cells=cells))
    fmt=lambda x: "NONE (no inset works)" if x is None else "f=%7.4f  c2=%.4f"%(x['f'],x['c2'])
    print(" %-6.4f %-9.5f | %-24s | %-24s | %-24s" % (g,c,fmt(cells[0]),fmt(cells[1]),fmt(cells[2])))
OUT['sweep']=rows; OUT['Gamma_max_over_ML']=float(Gam_max_over_ML); OUT['theta_acc']=float(theta_acc)

# where does it break?  a physical inset must satisfy f = O(1): call f > 5 unusable
print("\nlargest g for which the L=10 inset stays below f = 5, 2, 1 :")
lim={}
for cap in (5.0,2.0,1.0):
    glo,ghi=0.5,20.0
    def ok(g):
        b=best_f(10,g*Gam_max_over_ML*theta_acc,1.0/10)
        return b is not None and b[0]<=cap
    if not ok(glo): lim[cap]=None; print("   cap %-4s : none"%cap); continue
    for _ in range(40):
        m=.5*(glo+ghi)
        if ok(m): glo=m
        else: ghi=m
    lim[cap]=glo; print("   cap f<=%-4.1f : g <= %.3f  (c <= %.4f)"%(cap,glo,glo*Gam_max_over_ML*theta_acc))
OUT['g_limits']={str(k):(None if v is None else float(v)) for k,v in lim.items()}
json.dump(OUT,open("r2_results.json","w"),indent=1)
print("\nwrote r2_results.json")

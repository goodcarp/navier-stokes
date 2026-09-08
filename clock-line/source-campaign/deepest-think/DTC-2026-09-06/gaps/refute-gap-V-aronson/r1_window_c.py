#!/usr/bin/env python3
"""
REFUTER r1 -- the window constant c used in gap-V-aronson section 5.

Theorem V.1 as PROVED in that note defines
    Gamma := sup_{t<=tau} ||grad b(.,t)||_{Linf,op},     c := Gamma * tau,
and every appearance of e^{c} in the proof comes from Gronwall with the CONSTANT
majorant Gamma.  Section 5 then sets, for the accelerated doubling window,
    c_accel = theta_accel = 4(1-sqrt(2/3)) = 0.7340105          (v5_application.py line 31)
which is Gamma*tau ONLY IF Gamma = M L.  But in that window the strain accelerates:
prove-lagrangian's own closed-form ODE gives, at the innermost shell sigma=0,
    F(0,theta) = -2 log(1 - kappa theta/2),  lambda = e^F,  a(t) = a0 * lambda(t)   (Lemma 1)
with a0 = (M/2) L and lambda running 1 -> 3/2 over the window.  For the pure
axisymmetric strain the note itself proves ||grad_5 b||_op = 2a, so
    Gamma(t) = 2 a(t) = M L lambda(t),     sup Gamma = (3/2) M L.
Hence c = Gamma tau = (3/2) theta_accel, NOT theta_accel.

This script computes, from the ODE, three internally consistent versions of the
window constant and the resulting B1/B2 rates and required insets f.
Everything is derived here; no number is typed in from the note.
"""
import json, numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

OUT = {}
deg = np.pi/180.0

# ---------- 1. the window, from prove-lagrangian's OWN closed form ----------
# prove-lagrangian (4.1):  d_theta F(sigma,theta) = kappa H(sigma,theta),
#   H(sigma,theta) = int_sigma^1 e^{F(sigma',theta)} dsigma' = (1-sigma)/(1-(1-sigma)kappa theta/2).
# The strain at the tracked shell sigma=0 is a(theta) = M L kappa H(0,theta)   (NOT a0*lambda(0):
# each OUTER shell sigma' is stretched by its own lambda(sigma'), so it is the AVERAGE that enters).
# And ||grad_5 b||_op = 2a for the per-shell pure strain (the note's own section 1), so
#   Gamma(theta) = 2 kappa M L H(0,theta) = M L H(0,theta)  for kappa = 1/2.
kappa      = 0.5
theta_acc  = 4*(1-np.sqrt(2/3))                 # dimensionless window, theta = M L t
H0   = lambda th: 1.0/(1 - kappa*th/2)          # = Gamma/(M L)  when kappa = 1/2
Ffun = lambda th: -2*np.log(1 - kappa*th/2)     # F(0,theta) = log lambda
lam  = lambda th: np.exp(Ffun(th))
assert abs(lam(theta_acc) - 1.5) < 1e-12, lam(theta_acc)
# sanity: 2 * int a dt must equal 2 F = 2 log(3/2)
assert abs(2*kappa*quad(H0,0,theta_acc)[0] - 2*np.log(1.5)) < 1e-12
OUT['theta_accel'] = theta_acc
OUT['lambda_final'] = float(lam(theta_acc))
OUT['Gamma_max_over_ML'] = float(H0(theta_acc))

# (iii) what the note actually uses:  c = theta_acc, i.e. Gamma == M L (the t=0 value)
c_note  = theta_acc
# (i)  the theorem AS PROVED (c := sup_t Gamma * tau)
c_sup   = H0(theta_acc) * theta_acc
# (ii) the theorem with the (easy, but NOT done) time-dependent Gronwall:
#      every e^{c} becomes e^{C(tau)},  C(tau) = int_0^tau Gamma dt = 2 int a dt = 2 log(3/2)
C_tau   = quad(lambda th: H0(th), 0, theta_acc)[0]
OUT['c_theorem_as_proved_sup'] = float(c_sup)
OUT['C_tau_time_dependent']    = float(C_tau)
OUT['c_used_by_note']          = float(c_note)
OUT['C_tau_closed_form_2log32']= float(2*np.log(1.5))
print("window constants (accelerated window, kappa=1/2):")
print("  theta_accel                       = %.7f" % theta_acc)
print("  lambda(tau) (stretch)             = %.7f" % lam(theta_acc))
print("  Gamma(tau)/(M L) = H(0,tau)       = %.7f = sqrt(3/2)" % H0(theta_acc))
print("  c used by the note (= theta)      = %.7f" % c_note)
print("  c = (sup Gamma) tau  [AS PROVED]  = %.7f   (factor %.4f larger)" % (c_sup, c_sup/c_note))
print("  C(tau) = int Gamma dt [repair]    = %.7f   = 2 log(3/2) = %.7f"
      % (C_tau, 2*np.log(1.5)))

# ---------- 2. the three rate laws ----------
def rate_B1(c, n=5):  return np.exp(-2*c)/(4*n)
def rate_B2(c):       return c/(2*np.exp(2*c)*(np.exp(2*c)-1))
def B1(c,q,n=5):      return 2*n*np.exp(-rate_B1(c,n)*q)
EPS = np.linspace(1e-4, np.sqrt(2)-1e-4, 4000)
def _b2(rate,q,n=5):
    return float(np.min(n*np.log1p(2/EPS) - (1-EPS**2/2)**2*rate*q))
def B2(c,q,n=5):      return np.exp(_b2(rate_B2(c),q,n))
def BOUND(c,q,n=5):   return min(B1(c,q,n), B2(c,q,n))

# time-dependent version: replace e^{c}->e^{C(tau)} in the |Z|<=||Phi|| |N| step and
# V0 = 2 nu int_0^tau e^{2C(r)} dr  exactly.
#   C(theta) = int_0^theta lambda   (in units where M L = 1);  tau = theta_acc.
Cfun = lambda th: quad(H0, 0, th)[0]
V0_over_nutau = 2*quad(lambda th: np.exp(2*Cfun(th)), 0, theta_acc)[0]/theta_acc
rate_B1_td = np.exp(-2*C_tau)/(4*5)
rate_B2_td = np.exp(-2*C_tau)/(2*V0_over_nutau)
OUT['V0_over_nu_tau_timedep'] = float(V0_over_nutau)
OUT['rates'] = dict(
    note_B1=float(rate_B1(c_note)), note_B2=float(rate_B2(c_note)),
    asproved_B1=float(rate_B1(c_sup)), asproved_B2=float(rate_B2(c_sup)),
    timedep_B1=float(rate_B1_td), timedep_B2=float(rate_B2_td))
print("\nexponential rates (-log Bound / q):")
print("  %-34s B1 = %.6f   B2 = %.6f" % ("as the note computes", rate_B1(c_note), rate_B2(c_note)))
print("  %-34s B1 = %.6f   B2 = %.6f" % ("theorem as proved (sup Gamma)",   rate_B1(c_sup),  rate_B2(c_sup)))
print("  %-34s B1 = %.6f   B2 = %.6f" % ("time-dependent Gronwall (repair)",rate_B1_td, rate_B2_td))
print("  best-rate degradation vs the note:  as-proved %.3fx   repaired %.3fx"
      % (max(rate_B1(c_note),rate_B2(c_note))/max(rate_B1(c_sup),rate_B2(c_sup)),
         max(rate_B1(c_note),rate_B2(c_note))/max(rate_B1_td,rate_B2_td)))

# ---------- 3. redo section 5(ii): the required inset f ----------
def make_bound(kind):
    if kind=='note':     return lambda q: min(B1(c_note,q), B2(c_note,q))
    if kind=='asproved': return lambda q: min(B1(c_sup,q),  B2(c_sup,q))
    if kind=='timedep':
        def f(q,n=5):
            return min(2*n*np.exp(-rate_B1_td*q), np.exp(_b2(rate_B2_td,q,n)))
        return f
delta=7.5*deg; eps_eq=5*deg
def feature_d(f,phi0): return min(f,(1+f)*np.sin(phi0-delta),(1+f)*np.sin(np.pi/2-eps_eq-phi0))
def loss(f,L,phi0,bnd,c_for_q):
    d=feature_d(f,phi0)
    if d<=0: return np.inf
    q=(d/delta)**2*L/c_for_q
    return 2*(1+f)*np.sin(phi0)/np.sin(delta)*bnd(q)
def min_f(L,phi0,bnd,c_for_q,target):
    lo,hi=1e-4,60.0
    if loss(hi,L,phi0,bnd,c_for_q)>target: return None
    for _ in range(70):
        m=0.5*(lo+hi)
        if loss(m,L,phi0,bnd,c_for_q)<=target: hi=m
        else: lo=m
    return hi

# CAREFUL: nu tau is PHYSICAL and is fixed by the window length theta, not by c:
#   nu = M (rho0 delta)^2,  tau = theta_acc/(M L)  =>  nu tau = (rho0 delta)^2 theta_acc / L,
#   q = (d/(rho0 delta))^2 * L / theta_acc     -- always theta_acc in the denominator.
# The note writes q = (...)^2 L/c and gets the same thing ONLY because it sets c = theta_acc.
# Correcting c must therefore change the RATE but leave q alone.
print("\nrequired inset f so that A*min(B1,B2) <= 1/L   (delta=7.5deg, best phi0 in {30,45,60})")
print("  %-5s | %-22s | %-22s | %-22s" % ("L","note (c=%.4f)"%c_note,"as proved (c=%.4f)"%c_sup,"repaired (time-dep)"))
rows=[]
for L in (10,20,40,80,160):
    cells=[]
    for kind,cq in (('note',theta_acc),('asproved',theta_acc),('timedep',theta_acc)):
        bnd=make_bound(kind)
        best=None
        for phi0d in (30.,45.,60.):
            f=min_f(L,phi0d*deg,bnd,cq,1.0/L)
            if f is not None and (best is None or f<best[0]): best=(f,phi0d)
        if best is None: cells.append((None,None,None))
        else:
            f,phi0d=best
            infl=1.0/(1-np.log(1+f)/L)
            cells.append((f,phi0d,infl))
    rows.append(dict(L=L,note=cells[0],asproved=cells[1],timedep=cells[2]))
    fmt=lambda cc: "none" if cc[0] is None else "f=%.4f phi0=%2.0f c2x%.4f"%(cc[0],cc[1],cc[2])
    print("  %-5d | %-22s | %-22s | %-22s" % (L,fmt(cells[0]),fmt(cells[1]),fmt(cells[2])))
OUT['inset_table']=[{k:(None if v[0] is None else dict(f=float(v[0]),phi0=v[1],c2_inflation=float(v[2])))
                     if k!='L' else v for k,v in r.items()} for r in rows]

# c2 values
c2=8*(1-np.sqrt(2/3))
print("\nc2 = 8(1-sqrt(2/3)) = %.7f ;  inflated c2:" % c2)
for r in rows:
    s="  L=%-4d"%r['L']
    for k in ('note','asproved','timedep'):
        v=r[k]; s+= "  %-8s %s"%(k, "none" if v[0] is None else "%.4f"%(c2*v[2]))
    print(s)
OUT['c2_base']=float(c2)

# ---------- 4. does the O(1/L) survive? ----------
print("\nasymptotics: c2 inflation = L/(L - log(1+f)); f grows only like sqrt(1/rate),")
print("so the inflation is O(1/L) for every one of the three constants -> the")
print("statement c2 = 8(1-sqrt(2/3)) + o(1) is NOT affected; only the table is.")

json.dump(OUT, open("~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/refute-gap-V-aronson/r1_results.json","w"), indent=1)
print("\nwrote r1_results.json")

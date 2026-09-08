#!/usr/bin/env python3
"""x1 -- independent re-derivation of the attempt's LEMMA 1, plus a DEFLATION of it, plus the
truncation status of the constant 0.3958 that the attempt's gate cannot see.

DEFLATION.  eta_0 = -M sgn(z) h(phi)/r is a function of r ALONE (times sgn z and the taper).
T_lam scales r by the single constant lam, so AS A FUNCTION OF THE EULERIAN POINT
      eta_lam(w) = eta_0(T_lam^{-1}w) = -M lam sgn(z_w) h(phi'(w)) / r_w    on T_lam(shell),
i.e. exactly lam times the SAME 1/r profile.  And T_lam(shell) = {rho_0 < |w| G(t) < R} with
G(t) = sqrt((1-t^2)/lam^2 + lam^4 t^2), so the LOG-THICKNESS of the support is log(R/rho_0) = L
for EVERY direction t (both shell radii scale by the same 1/G(t)).  a is linear in eta with a
scale-free kernel and the rho-integral is exactly L in every direction, hence for h == 1
      a_lam(0) = lam a_0(0),
with NO cancellation between "rotation to the equator" and "inward collapse near the axis" to
check.  The attempt's Lagrangian form P_1(lam) = 3 int_0^1 v^2 (Av^2+B)^{-5/2}dv = lam is the
same statement after an unnecessary change of variables.
"""
import numpy as np, sympy as sp, mpmath as mp, json
OUT={}; M=1.0

# --- A. the attempt's Lagrangian identity, symbolically
v,lam = sp.symbols('v lam', positive=True)
A_=lam**2-lam**-4; B_=lam**-4
Ant=v**3/(3*B_*(A_*v**2+B_)**sp.Rational(3,2))
OUT['A_antideriv_residual']=str(sp.simplify(sp.diff(Ant,v)-v**2/(A_*v**2+B_)**sp.Rational(5,2)))
OUT['A_P1_minus_lam']=str(sp.simplify(3*(Ant.subs(v,1)-Ant.subs(v,0))-lam))
mp.mp.dps=30; worst=0
for s in ['1.0','1.1','1.5','2.0','5.0']:
    Lm=mp.mpf(s); Am=Lm**2-Lm**-4; Bm=Lm**-4
    worst=max(worst,abs(3*mp.quad(lambda x:x**2/(Am*x**2+Bm)**mp.mpf(2.5),[0,1])-Lm)/Lm)
OUT['A_mpmath_rel']=float(worst)

# --- B. the Eulerian derivation (mine).  a(0) = (3/4) int int (-cos phi sin^3 phi) eta drho dphi
#     with rho eta_lam = -M lam sgn(t) h(phi'(t))/sqrt(1-t^2)  and rho-support of log-length L
#     for EVERY t.  G(t) therefore cannot appear in the answer for h == 1.
def a0_eulerian(hfun,lamv,n=2000001):
    t=np.linspace(-1+1e-13,1-1e-13,n)
    G=np.sqrt((1-t**2)/lamv**2+lamv**4*t**2)
    tp=np.clip(lamv**2*t/G,-1,1); php=np.arccos(tp)
    f=-M*lamv*np.sign(t)*hfun(php)/np.sqrt(1-t**2)      # = rho*eta_lam
    # (3/4) int dphi (-cos phi) sin^3 phi * f * L ; dphi = -dt/sqrt(1-t^2)
    return 0.75*np.trapz((-t)*(1-t**2)*f, t)
h1=lambda p:np.ones_like(p)
def htap(d): return lambda p: np.minimum(1.0,np.minimum(p,np.pi-p)/d)
OUT['B_eulerian_flat']=[[l,a0_eulerian(h1,l),0.5*l] for l in [1.0,1.1,1.25,1.5,2.0,3.0,5.0]]
d75=np.deg2rad(7.5)
OUT['B_eulerian_taper7.5']=[[l,a0_eulerian(htap(d75),l)] for l in [1.0,1.25,1.5]]
# demonstrate G(t) is absent for h==1: delete it and the answer is unchanged
def a0_noG(lamv,n=2000001):
    t=np.linspace(-1+1e-13,1-1e-13,n)
    return 0.75*np.trapz((-t)*(1-t**2)*(-M*lamv*np.sign(t)/np.sqrt(1-t**2)),t)
OUT['B_flat_without_G']=[[l,a0_noG(l)] for l in [1.0,1.25,1.5,2.0,5.0]]

# --- C. truncation status of "sum_{l>=3} |H_l| ||C_l||_inf/(l(l+3)-4) = 0.3958"
# C_l^{3/2}(t) = P'_{l+1}(t); integrate in theta = arccos t (trapezoid, spectrally accurate)
LM=2001; NT=60001
th=np.linspace(0.0,np.pi,NT); tq=np.cos(th); wq=np.gradient(th)*np.sin(th)   # dt = sin th dth
g=-M*np.sign(tq)/np.sqrt(np.maximum(1-tq**2,1e-300))
wt=wq*(1-tq**2)*g
wt[0]=wt[-1]=0.0
H=np.zeros(LM+1)
Cm2=np.ones_like(tq); Cm1=3*tq
H[1]=np.dot(wt,Cm1)/(2*3/2.5)
for n in range(2,LM+1):
    Cn=(2*tq*(n+0.5)*Cm1-(n+1)*Cm2)/n
    if n%2==1: H[n]=np.dot(wt,Cn)/((n+1)*(n+2)/(n+1.5))
    Cm2,Cm1=Cm1,Cn
l=np.arange(LM+1); Cinf=(l+1)*(l+2)/2.0; den=np.maximum(l*(l+3)-4.0,1.0)
term=np.abs(H)*Cinf/den
OUT['C_H_check']={'H1':float(H[1]),'H3':float(H[3]),'H5':float(H[5]),'exact':[-5/6,3/40,-247/1680]}
OUT['C_partial_sums']={str(k):float(np.sum(term[3:k+1])) for k in (101,401,601,801,1201,1601,2001)}
OUT['C_term_times_l32']={str(k):float(term[k]*k**1.5) for k in (101,401,801,1601)}
# tail estimate beyond k using term ~ c l^{-3/2} over odd l: sum ~ c/sqrt(k)
c=float(np.mean([term[k]*k**1.5 for k in (801,1001,1201,1401,1601)]))
OUT['C_tail_const_c']=c
OUT['C_extrapolated_total_from_601']=float(np.sum(term[3:602])+c/np.sqrt(601))
OUT['C_note_claim']=0.3958; OUT['C_s6_truncation']=601
print(json.dumps(OUT,indent=1,default=str))
json.dump(OUT,open(__file__.replace('.py','_results.json'),'w'),indent=1,default=str)

#!/usr/bin/env python3
"""x1 -- REFUTER instrument: size Theorem V.4's OWN three error terms
(E_hess, E_4, E_tail) at the campaign design point, from the theorem's own
definitions, by a route written independently of b5_sizing.py.

Design point (units rho0 = 1, M = 1):
   nu       = (rho0 sin delta)^2 M          (rho0 delta := sqrt(nu/M), delta = 7.5 deg)
   tau      = theta_max/(M L)
   lam(th)  = (1 - th/4)^{-2},  th = M L t,  th in [0, theta_max], theta_max = 4(1-sqrt(2/3))
   r0       = (1+f) sin(phi0)
   sigma_y  = nu int_0^tau lam^{-2} dt = nu I2/(M L)
   sigma_z  = nu int_0^tau lam^{+4} dt = nu I4/(M L)
   C_tau    = diag(2 sigma_y I_4, 2 sigma_z)
   s_C      = sigma_y/r0^2
   N        = ||eta_0||_inf / |eta_0(x0)| = (M/(rho0 sin delta))/(M/r0) = r0/sin delta
   d        = f rho0 (radial inset to the inner edge)  [or d = rho0 sin delta]
   K2       = K2hat M/rho0
Everything from Theorem V.4 as WRITTEN in write/V-b-bulk-viscous-loss/PROOF.md sec 5:
   E_hess = (r0/R_-^2) (K2/2) e^c tau V1  +  4 N E[W]/d,       E[W] <= (K2/2) e^c tau V1
   E_4    = (r0/R_-^5) [ (tr C)^2 + 2 tr(C^2) ]
   E_tail = 2 N P + P + sum_{k=1..3} r0^{-k} (E|Z^a|^{2k})^{1/2} P^{1/2}
   P      = P(|Z_tau|>d/2) + P(|Z^a_tau|>d/2) + P(|Z^a_tau| >= d)
   V1     = n nu tau (e^{2c}-1)/c,  n = 5, c = sqrt(3/2) theta_max
"""
import json, math
import numpy as np
import sympy as sp
from scipy.stats import chi2
from scipy.special import roots_hermitenorm

R = {}
# ---------- window constants, re-derived from scratch (sympy, exact) ----------
th = sp.symbols('theta', nonnegative=True)
lam = (1 - th/4)**-2                              # kappa = 1/2 : (1 - kappa th/2)^{-2}
theta_max = sp.Rational(4,1)*(1 - sp.sqrt(sp.Rational(2,3)))
assert sp.simplify(lam.subs(th, theta_max) - sp.Rational(3,2)) == 0
I2 = sp.integrate(lam**-2, (th, 0, theta_max))
I4 = sp.integrate(lam**4,  (th, 0, theta_max))
c_w = sp.sqrt(sp.Rational(3,2))*theta_max
R['theta_max'] = float(theta_max); R['I2'] = float(I2); R['I4'] = float(I4)
R['I2_exact'] = str(sp.simplify(I2)); R['I4_exact'] = str(sp.simplify(I4))
R['c_window'] = float(c_w)
R['theta_max_over_I2'] = float(theta_max/I2)
TH, I2f, I4f, C = float(theta_max), float(I2), float(I4), float(c_w)
n = 5
print(f"theta_max={TH:.13f}  I2={I2f:.13f}  I4={I4f:.13f}  c={C:.13f}  th/I2={TH/I2f:.6f}")

sd = math.sin(math.radians(7.5))
R['sin_delta'] = sd

# ---------- exact tail P(|Z^a| > a) for Z^a ~ N(0, diag(2sy I4, 2sz)) ----------
_hx, _hw = roots_hermitenorm(400); _hw = _hw/ _hw.sum()
def P_norm_gt(a, sy, sz):
    """P(|Z|>a), Z=(G_y in R^4 var 2sy each, G_z var 2sz).  |Z|^2 = 2sy X4 + 2sz g^2."""
    if a <= 0: return 1.0
    g = math.sqrt(2*sz)*_hx
    rem = a*a - g*g
    p = np.where(rem <= 0, 1.0, chi2.sf(np.maximum(rem,0.0)/(2*sy), df=4))
    return float((_hw*p).sum())

def moments(sy, sz):
    trC = 8*sy + 2*sz
    trC2 = 4*(2*sy)**2 + (2*sz)**2
    m2 = trC
    m4 = trC**2 + 2*trC2
    # E|Z|^6 for a Gaussian: use the chi-representation numerically
    g = math.sqrt(2*sz)*_hx
    # E[(2sy X4 + g^2)^3] with X4 ~ chi2_4:  E X4=4, EX4^2=24, EX4^3=192
    a = 2*sy
    e = (a**3*192 + 3*a**2*24*(g*g) + 3*a*4*(g**4) + g**6)
    m6 = float((_hw*e).sum())
    return trC, trC2, m2, m4, m6

# ---------- the sizing ----------
def size(L, f=0.0, phi0_deg=30.0, K2hat=1.0, d_mode='inset', Rminus_mode='r0_minus_d'):
    M = 1.0; rho0 = 1.0
    nu  = (rho0*sd)**2 * M
    tau = TH/(M*L)
    r0  = (1+f)*math.sin(math.radians(phi0_deg))
    sy  = nu*I2f/(M*L); sz = nu*I4f/(M*L)
    sC  = sy/r0**2
    N   = r0/sd
    if d_mode == 'inset':
        # gap-V sec 5's geometric distance to the nearest non-plateau set, in the initial
        # configuration: radial inner edge, angular taper edge, equatorial mollification
        rho_star = (1+f)*rho0
        d = min(f*rho0,
                rho_star*math.sin(max(math.radians(phi0_deg)-math.radians(7.5), 1e-9)),
                rho_star*math.sin(max(math.radians(85.0)-math.radians(phi0_deg), 1e-9)))
    else:
        d = rho0*sd
    d = min(d, 0.999*r0)
    Rm  = r0 - d if Rminus_mode == 'r0_minus_d' else r0
    K2  = K2hat*M/rho0
    V1  = n*nu*tau*(math.exp(2*C)-1)/C
    EW  = 0.5*K2*math.exp(C)*tau*V1
    Eh1 = (r0/Rm**2)*EW
    Eh2 = 4*N*EW/d
    trC, trC2, m2, m4, m6 = moments(sy, sz)
    E4  = (r0/Rm**5)*m4
    # P: the theorem's own event set.  P(|Z_tau|>d/2) is bounded below by the affine one
    # (Z and Z^a are coupled and |Z-Z^a| <= W); we use the affine value = OPTIMISTIC.
    Pa_half = P_norm_gt(d/2, sy, sz); Pa_full = P_norm_gt(d, sy, sz)
    Pc = 2*Pa_half + Pa_full
    Et  = 2*N*Pc + Pc + sum(r0**(-k)*math.sqrt([m2,m4,m6][k-1])*math.sqrt(Pc) for k in (1,2,3))
    return dict(L=L, f=f, phi0=phi0_deg, K2hat=K2hat, nu=nu, tau=tau, r0=r0, d=d, N=N,
                sigma_y=sy, sigma_z=sz, s_C=sC, V1=V1, EW=EW,
                E_hess_term1=Eh1, E_hess_term2=Eh2, E_hess=Eh1+Eh2,
                term2_over_term1=(Eh2/Eh1 if Eh1>0 else float("nan")), E_4=E4,
                P_half=Pa_half, P_full=Pa_full, Pcal=Pc, E_tail=Et,
                total=Eh1+Eh2+E4+Et)

# ---------- (1) the E_hess/eps_bulk table, the seat's vs the theorem's ----------
print("\n(1) E_hess/eps_bulk  at K2hat=1, f=0 (d = rho0 sin delta, R_-=r0-d) ")
print(" L    seat's table   term1 only (correct eps)   term1+term2 (theorem)   ratio to seat")
rows=[]
for L in [10,40,160,640]:
    s = size(L, f=0.0, K2hat=1.0, d_mode='onedelta')
    seat = 0.5*math.exp(C)*n*(math.exp(2*C)-1)/C * 1.0*TH*0.5/L     # b5's hess_ratio
    t1 = s['E_hess_term1']/s['s_C']; tot = s['E_hess']/s['s_C']
    rows.append(dict(L=L, seat=seat, term1_over_eps=t1, full_over_eps=tot,
                     full_over_seat=tot/seat, term2_over_term1=s['term2_over_term1']))
    print(f"{L:5d}   {seat:10.4f}     {t1:14.4f}        {tot:14.4f}        {tot/seat:8.2f}")
R['E_hess_table'] = rows

print("\n  (breakdown at L=10, f=0, K2hat=1, d=rho0 sin delta)")
s = size(10, 0.0, K2hat=1.0, d_mode='onedelta')
for k in ['s_C','EW','E_hess_term1','E_hess_term2','term2_over_term1','E_4','Pcal','E_tail','total']:
    print(f"    {k:18s} {s[k]:.6e}")
R['breakdown_L10_f0'] = s

# ---------- (2) corrected K2hat thresholds ----------
print("\n(2) largest K2hat with E_hess <= eps_bulk (d = rho0 sin delta, R_-=r0-d)")
thr = {}
for L in [10,40,160,640,1200]:
    s1 = size(L, 0.0, K2hat=1.0, d_mode='onedelta')
    thr[L] = s1['s_C']/s1['E_hess']
    print(f"   L={L:5d}   seat: {1.0/(0.5*math.exp(C)*n*(math.exp(2*C)-1)/C*TH*0.5/L):9.3f}"
          f"    theorem: {thr[L]:9.5f}")
R['K2hat_threshold_corrected'] = thr
# smallest L at which K2hat = 1 is admissible
lo, hi = 1.0, 1e6
for _ in range(200):
    mid = 0.5*(lo+hi)
    s1 = size(mid, 0.0, K2hat=1.0, d_mode='onedelta')
    if s1['E_hess'] <= s1['s_C']: hi = mid
    else: lo = mid
R['L_min_K2hat_1_theorem'] = hi
print(f"   smallest L with K2hat=1 admissible: seat says ~{0.5*math.exp(C)*n*(math.exp(2*C)-1)/C*TH*0.5:.2f}"
      f" ;  theorem says {hi:.1f}")

# ---------- (3) the theorem's E_tail vs the seat's sec-5 tail table ----------
print("\n(3) tail at d = rho0 sin delta (one dissipation length), f=0")
print(" L     seat A*P(1-D, thr=d)    P(|Z^a|>d/2)   P(|Z^a|>d)   theorem's E_tail   1/L")
rows=[]
for L in [5,10,20,40,80,160,320,640]:
    s = size(L, 0.0, K2hat=0.0, d_mode='onedelta')
    # the seat's own 1-D construction, reproduced
    phi0 = math.radians(30.0)
    var_n = 2*(I2f*math.sin(phi0)**2 + I4f*math.cos(phi0)**2)/L
    seatP = 0.5*math.exp(-1.0/(2*var_n)); A = (1+sd)*math.sin(phi0)/sd + 1.0
    rows.append(dict(L=L, seat_AP=A*seatP, P_half=s['P_half'], P_full=s['P_full'],
                     E_tail=s['E_tail'], target=1.0/L, eps_bulk=s['s_C']))
    print(f"{L:5d}   {A*seatP:.4e}          {s['P_half']:.4e}   {s['P_full']:.4e}   "
          f"{s['E_tail']:.4e}      {1.0/L:.5f}")
R['tail_table'] = rows
lo, hi = 2.0, 1e7
for _ in range(300):
    mid=0.5*(lo+hi); s=size(mid,0.0,K2hat=0.0,d_mode='onedelta')
    if s['E_tail'] <= 1.0/mid: hi=mid
    else: lo=mid
R['L_min_one_dissipation_length_theorem_vs_1L'] = hi
print(f"   smallest L with theorem's E_tail <= 1/L at d=rho0 sin delta: {hi:.4g}   (seat: 24.603)")

# ---------- (4) required inset f, theorem's own bound ----------
print("\n(4) smallest inset f with (E_hess+E_4+E_tail) <= 1/L   [K2hat=0, i.e. tail+E4 only]")
print(" L    f (seat, sec5 vs eps_bulk)  f (theorem, K2hat=0, vs 1/L)   f (gap-V vs 1/L)")
seat_f = {10:0.25972, 20:0.19228, 40:0.14189, 80:0.10441, 160:0.07663}
gapV_f = {10:0.9131, 20:0.6590, 40:0.4694, 80:0.3382, 160:0.2437}
rows=[]
for L in [10,20,40,80,160]:
    lo, hi = 1e-4, 20.0
    for _ in range(200):
        mid=0.5*(lo+hi); s=size(L, mid, K2hat=0.0, d_mode='inset')
        if s['E_tail']+s['E_4'] <= 1.0/L: hi=mid
        else: lo=mid
    rows.append(dict(L=L, f_theorem=hi, f_seat=seat_f[L], f_gapV=gapV_f[L]))
    print(f"{L:5d}   {seat_f[L]:.5f}                    {hi:.5f}                    {gapV_f[L]:.4f}")
R['inset_table'] = rows

json.dump(R, open('x1_results.json','w'), indent=1, default=str)
print("\nWROTE x1_results.json")

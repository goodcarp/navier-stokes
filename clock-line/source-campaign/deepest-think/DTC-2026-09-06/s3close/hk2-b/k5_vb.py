#!/usr/bin/env python3
"""k5 -- what my K2 buys in Theorem V.4 (write/V-b-bulk-viscous-loss/PROOF.md Sec.5),
sized with the refuter's corrected formulae (write/refute-V-b-bulk-viscous-loss Sec.1):

   E_hess = (r0/R_-^2) * (1/2) K2 e^c tau V1   +   4 N E[W]/d ,  E[W] <= (1/2)K2 e^c tau V1
   V1 = n nu tau (e^{2c}-1)/c ,  n = 5
   N  = ||eta_0||_inf r0 / M = (1+f) sin phi0 / sin delta
   eps_bulk = (l_nu/r0)^2 I2 / L ,  l_nu = rho0 sin delta = sqrt(nu/M)

Every constant is rebuilt here from Theorem V.4's own statement; the refuter's published
Eh1/Eh2 at L=10 are used ONLY as an after-the-fact agreement check.
"""
import json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}

th_max = 4.0*(1.0 - math.sqrt(2.0/3.0))
c      = math.sqrt(1.5)*th_max
I2     = 4.0/5.0 - 16.0*math.sqrt(6.0)/135.0
phi0   = math.radians(30.0); delta = math.radians(7.5)
sd     = math.sin(delta)
M = rho0 = 1.0
n5 = 5

RES['theta_max'] = th_max; RES['c'] = c; RES['I2'] = I2; RES['sin_delta'] = sd


def sizing(L, K2hat, f=0.0, d=None):
    r0 = (1.0 + f)*rho0*math.sin(phi0)
    lnu = rho0*sd
    if d is None: d = lnu
    Rm = r0 - d
    tau = th_max/(M*L)
    nutau = lnu**2*th_max/L
    V1 = n5*nutau*(math.exp(2*c) - 1.0)/c
    EW = 0.5*K2hat*(M/rho0)*math.exp(c)*tau*V1
    N = (1.0 + f)*math.sin(phi0)/sd
    Eh1 = (r0/Rm**2)*EW
    Eh2 = 4.0*N*EW/d
    eps = (lnu/r0)**2*I2/L
    E4 = (r0/Rm**5)*0.0            # filled by caller if wanted
    return dict(L=L, K2hat=K2hat, f=f, r0=r0, Rm=Rm, tau=tau, nutau=nutau, V1=V1, EW=EW,
                N=N, Eh1=Eh1, Eh2=Eh2, E_hess=Eh1 + Eh2, eps_bulk=eps,
                ratio_eps=(Eh1 + Eh2)/eps, ratio_1overL=(Eh1 + Eh2)*L)


chk = sizing(10.0, 1.0)
RES['check_vs_refuter_L10'] = dict(Eh1=chk['Eh1'], Eh2=chk['Eh2'],
                                   refuter_Eh1=1.1572e-2, refuter_Eh2=3.7089e-1,
                                   rel_Eh1=abs(chk['Eh1'] - 1.1572e-2)/1.1572e-2,
                                   rel_Eh2=abs(chk['Eh2'] - 3.7089e-1)/3.7089e-1)
print(f"k5  independent rebuild of Theorem V.4's E_hess at L=10, Khat=1, f=0:")
print(f"      Eh1 = {chk['Eh1']:.6e}  (refuter 1.1572e-2, rel {RES['check_vs_refuter_L10']['rel_Eh1']:.1e})")
print(f"      Eh2 = {chk['Eh2']:.6e}  (refuter 3.7089e-1, rel {RES['check_vs_refuter_L10']['rel_Eh2']:.1e})")
print(f"      E_hess/eps_bulk = {chk['ratio_eps']:.4f}  (refuter's 1101.1*Khat/L -> {1101.1/10:.4f})")

# the two thresholds, from the theorem
def khat_max(L, budget):
    s = sizing(L, 1.0)
    if budget == 'eps': return s['eps_bulk']/s['E_hess']
    return (1.0/L)/s['E_hess']

RES['khat_thresholds'] = {str(L): dict(vs_eps=khat_max(L, 'eps'), vs_1overL=khat_max(L, 'L'))
                          for L in (10, 40, 160, 640, 1200)}
print("\nk5  largest admissible Khat = K2 rho0/M")
print("      L        vs eps_bulk      vs 1/L")
for L in (10, 40, 160, 640, 1200):
    t = RES['khat_thresholds'][str(L)]
    print(f"    {L:6d}   {t['vs_eps']:12.6f}  {t['vs_1overL']:12.6f}")

# apply MY measured K2
try:
    k4 = json.load(open(os.path.join(HERE, 'k4_results.json')))['direct']
    k3 = json.load(open(os.path.join(HERE, 'k3_results.json')))['bound_table']
    khat_bound = max(r['K2'] for r in k3 if r['L'] == 10.0)
    khat_exact = max(r['K2_exact'] for r in k4 if r['L'] == 10.0)
    khat_bound_15 = max(r['K2'] for r in k3 if r['L'] == 10.0 and r['rho_star'] >= 1.5)
    khat_exact_15 = max(r['K2_exact'] for r in k4 if r['L'] == 10.0 and r['rho_star'] >= 1.5)
except Exception as e:
    khat_bound = khat_exact = khat_bound_15 = khat_exact_15 = float('nan')
    print("  (k3/k4 results not present yet:", e, ")")

RES['my_khat'] = dict(bound_all=khat_bound, exact_all=khat_exact,
                      bound_inset=khat_bound_15, exact_inset=khat_exact_15)
print(f"\nk5  my K2 (Khat = K2 rho0/M), campaign datum, over N_tau:")
print(f"      proved bound, rho*>=1.0 : {khat_bound:.4f}")
print(f"      measured,     rho*>=1.0 : {khat_exact:.4f}")
print(f"      proved bound, rho*>=1.5 : {khat_bound_15:.4f}")
print(f"      measured,     rho*>=1.5 : {khat_exact_15:.4f}")

tab = []
for L in (10, 40, 160, 640):
    for name, kh in (('bound', khat_bound), ('measured', khat_exact),
                     ('bound_inset', khat_bound_15), ('measured_inset', khat_exact_15)):
        s = sizing(float(L), kh)
        tab.append(dict(L=L, which=name, Khat=kh, E_hess=s['E_hess'],
                        over_eps=s['ratio_eps'], over_1overL=s['ratio_1overL']))
RES['E_hess_table'] = tab
print("\nk5  E_hess with my K2  (E_hess/eps_bulk  and  E_hess*L = E_hess/(1/L))")
print("      L      which             Khat      E_hess     /eps_bulk     *L")
for t in tab:
    print(f"    {t['L']:5d}   {t['which']:15s} {t['Khat']:8.3f}  {t['E_hess']:10.4e}  "
          f"{t['over_eps']:11.2f}  {t['over_1overL']:9.2f}")

# the L at which E_hess drops below each budget, with my Khat
def Lstar(kh, budget):
    lo, hi = 1.0, 1e12
    for _ in range(200):
        mid = math.sqrt(lo*hi)
        s = sizing(mid, kh)
        v = s['E_hess'] - (s['eps_bulk'] if budget == 'eps' else 1.0/mid)
        if v > 0: lo = mid
        else: hi = mid
    return math.sqrt(lo*hi)

RES['Lstar'] = {name: dict(vs_eps=Lstar(kh, 'eps'), vs_1overL=Lstar(kh, 'L'))
                for name, kh in (('bound', khat_bound), ('measured', khat_exact),
                                 ('bound_inset', khat_bound_15), ('measured_inset', khat_exact_15))
                if kh == kh}
print("\nk5  smallest L for which E_hess is below the budget")
for k, v in RES['Lstar'].items():
    print(f"      {k:16s} vs eps_bulk: L >= {v['vs_eps']:12.1f}    vs 1/L: L >= {v['vs_1overL']:10.1f}")

json.dump(RES, open(os.path.join(HERE, 'k5_results.json'), 'w'), indent=1, sort_keys=True)
print("\nwrote k5_results.json")

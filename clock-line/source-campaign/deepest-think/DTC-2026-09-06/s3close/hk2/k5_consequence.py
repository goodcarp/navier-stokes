#!/usr/bin/env python3
"""k5 -- what an explicit K2 buys in Theorem V.4, using the REFUTER's corrected sizing
(write/refute-V-b-bulk-viscous-loss/PROOF.md sec 1):

   E_hess = (r0/R_-^2) * (1/2) K2 e^c tau V1   +   4 N E[W]/d ,   E[W] <= (1/2)K2 e^c tau V1
   V1 = n nu tau (e^{2c}-1)/c ,  n = 5 ,  N = ||eta_0||_inf r0/M
   eps_bulk = (rho0 delta / r0)^2 I2 / L
Window (prove-lagrangian 4.1): theta_max = 4(1-sqrt(2/3)), tau = theta_max/(M L),
   c = sqrt(3/2) theta_max, nu = (rho0 delta)^2 M, I2 = 4/5 - 16 sqrt6/135.
Everything recomputed here from scratch; the refuter's 87.166 is REPRODUCED as a control.
"""
import json, math
import numpy as np
import sympy as sp

RES = {}
th  = sp.symbols('theta', nonnegative=True); kap = sp.Rational(1,2)
lam = (1-kap*th/2)**-2
theta_max = 2*(1-sp.sqrt(sp.Rational(2,3)))/kap
I2 = sp.integrate(lam**-2, (th, 0, theta_max))
TH = float(theta_max); I2f = float(I2); c = math.sqrt(1.5)*TH
RES.update(theta_max=TH, I2=I2f, I2_exact=str(sp.nsimplify(sp.simplify(I2))), c=c)
print(f"theta_max = {TH:.9f}   I2 = {I2f:.10f} ({sp.nsimplify(sp.simplify(I2))})   c = {c:.7f}")

n = 5
def pieces(K2hat, L, delta_deg=7.5, phi0_deg=30.0, f=0.0, d=None, eta0sup=None):
    """all in units rho0 = M = 1."""
    sd = math.sin(math.radians(delta_deg)); s0 = math.sin(math.radians(phi0_deg))
    r0 = (1+f)*s0
    nu = sd**2                      # rho0 delta = sqrt(nu/M) -> nu = delta^2 (rho0=M=1); sin delta used
    tau = TH/L
    if d is None: d = sd
    Rm = r0 - d
    V1 = n*nu*tau*(math.exp(2*c)-1)/c
    EW = 0.5*K2hat*math.exp(c)*tau*V1
    N  = (eta0sup if eta0sup is not None else 1.0/sd)*r0
    Eh1 = (r0/Rm**2)*EW
    Eh2 = 4*N*EW/d
    eps = (sd/r0)**2*I2f/L
    return dict(Eh1=Eh1, Eh2=Eh2, E_hess=Eh1+Eh2, eps_bulk=eps, ratio=(Eh1+Eh2)/eps,
                N=N, d=d, Rminus=Rm, V1=V1, EW=EW, tau=tau, nu=nu, r0=r0)

# ---- control: reproduce the seat's own (withdrawn) sizing and the refuter's factor 87.166
def seat_ratio(K2hat, L, f=0.0, phi0_deg=30.0):
    r0 = (1+f)*math.sin(math.radians(phi0_deg))
    return 0.5*math.exp(c)*n*(math.exp(2*c)-1)/c * K2hat * TH * r0 / L
ctrl = []
for L in (10,40,160,640):
    p = pieces(1.0, L); sr = seat_ratio(1.0, L)
    ctrl.append(dict(L=L, theorem=p['ratio'], seat=sr, factor=p['ratio']/sr,
                     Eh2_over_Eh1=p['Eh2']/p['Eh1']))
    print(f"  L={L:4d}:  TheoremV.4 E_hess/eps = {p['ratio']:.6f}   seat's = {sr:.6f}   "
          f"factor = {p['ratio']/sr:.4f}   Eh2/Eh1 = {p['Eh2']/p['Eh1']:.4f}")
RES['control_refuter_factor'] = ctrl

# ---- with an explicit K2hat -------------------------------------------------------------
print("\n=== E_hess/eps_bulk and E_hess*L, with an explicit K2hat (d = rho0 sin delta) ===")
tab = {}
for K2hat in (1.0, 5.0, 12.0, 25.0, 50.0):
    row = {}
    for L in (10,40,160,640,2560):
        p = pieces(K2hat, L)
        row[L] = dict(ratio=p['ratio'], EhL=p['E_hess']*L)
    tab[K2hat] = row
    print(f"  K2hat={K2hat:6.1f}: " + "  ".join(
        f"L={L}: E_h/eps={row[L]['ratio']:9.2f} E_h*L={row[L]['EhL']:8.3f}" for L in (10,40,160,640)))
RES['table_K2hat'] = {str(k): {str(l): v for l,v in r.items()} for k,r in tab.items()}
# crossings
print("\n=== L needed for E_hess <= eps_bulk and E_hess <= 1/L ===")
cross = {}
for K2hat in (1.0, 5.0, 12.0, 25.0, 50.0):
    # E_hess/eps ~ C1 K2hat / L  ->  L_eps = C1 K2hat ;  E_hess ~ C2 K2hat/L^2 -> L_1overL = C2 K2hat
    p = pieces(K2hat, 100.0)
    L_eps = p['ratio']*100.0
    L_1L  = p['E_hess']*100.0**2
    cross[K2hat] = dict(L_for_eps_bulk=L_eps, L_for_1_over_L=L_1L)
    print(f"  K2hat={K2hat:6.1f}:  E_hess <= eps_bulk from L = {L_eps:10.1f} ;"
          f"   E_hess <= 1/L from L = {L_1L:8.1f}")
RES['crossings'] = {str(k): v for k,v in cross.items()}
json.dump(RES, open('k5_results.json','w'), indent=1)
print("\nWROTE k5_results.json")

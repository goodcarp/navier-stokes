#!/usr/bin/env python3
"""b5 -- sizing V-b for the campaign's datum and window.  Everything computed, nothing typed."""
import json, math
import numpy as np
import sympy as sp

RES = {}
# ---------------------------------------------------------------- window (prove-lagrangian 4.1)
kap = sp.Rational(1,2)
th  = sp.symbols('theta', nonnegative=True)
lam = (1 - kap*th/2)**-2
theta_max = 2*(1-sp.sqrt(sp.Rational(2,3)))/kap                      # lam(theta_max) = 3/2
RES['theta_max'] = float(theta_max)
RES['lam_at_theta_max'] = float(lam.subs(th, theta_max))
I2 = sp.integrate(lam**-2, (th, 0, theta_max))                        # int lam^{-2} dtheta
I4 = sp.integrate(lam**4,  (th, 0, theta_max))                        # int lam^{ 4} dtheta
I0 = theta_max
RES['I2_exact'] = sp.nsimplify(sp.simplify(I2)).__str__(); RES['I2'] = float(I2)
RES['I4_exact'] = sp.nsimplify(sp.simplify(I4)).__str__(); RES['I4'] = float(I4)
# Gamma(t) = ||grad_5 b||_op = 2 a(0,t) = M L lam^{1/2}  (prove-lagrangian (4.1), kappa=1/2)
# so   c := (sup_t Gamma) tau = sqrt(3/2) theta_max ;  C(tau) := int_0^tau Gamma dt = 2 log(3/2)
RES['c_window']  = float(sp.sqrt(sp.Rational(3,2))*theta_max)
RES['C_tau_int'] = float(2*sp.log(sp.Rational(3,2)))
RES['c_frozen_Gamma0'] = float(theta_max)
RES['sigma_z_over_sigma_y'] = float(I4/I2)
RES['I2_over_theta_max'] = float(I2/theta_max)
print(f"theta_max = {RES['theta_max']:.9f}  lam = {RES['lam_at_theta_max']:.9f}")
print(f"I2 = int lam^-2 dtheta = {RES['I2']:.9f}   ({RES['I2_exact']})")
print(f"I4 = int lam^+4 dtheta = {RES['I4']:.9f}   ({RES['I4_exact']})")
print(f"sigma_z/sigma_y = {RES['sigma_z_over_sigma_y']:.6f}   I2/theta_max = {RES['I2_over_theta_max']:.6f}")
print(f"c = (sup Gamma) tau = {RES['c_window']:.7f} ; C(tau) = int Gamma dt = {RES['C_tau_int']:.7f} ;"
      f" Gamma(0) tau = {RES['c_frozen_Gamma0']:.7f}")

# ---------------------------------------------------------------- the bulk term
# design point: rho0 delta = sqrt(nu/M);  tau = c/(M L) with c = theta_max/kappa? NO:
# theta = M L t  =>  t = theta/(M L),  tau = theta_max/(M L).
# nu int_0^tau dt/r(t)^2 = (nu/(r0^2 M L)) I2 = ((rho0 delta)^2/(r0^2 L)) I2 .
delta = math.radians(7.5); sd = math.sin(delta)
def bulk_eps(L, f=0.0, phi0_deg=30.0):
    r0_over_rho0 = (1+f)*math.sin(math.radians(phi0_deg))
    return (sd/r0_over_rho0)**2 * RES['I2']/L
def bulk_eps_lagrangian(L, f=0.0, phi0_deg=30.0):
    """prove-lagrangian sec 4(3)'s pricing: (delta^2/sin^2 phi) * (c/L) with c = theta_max."""
    return (sd/math.sin(math.radians(phi0_deg))/(1+f))**2 * RES['theta_max']/L
rows = []
for L in [10, 20, 40, 80, 160, 320, 640]:
    e0 = bulk_eps(L); el = bulk_eps_lagrangian(L)
    rows.append(dict(L=L, eps_bulk_exact=e0, eps_lagrangian=el, ratio=el/e0,
                     target_1_over_L=1.0/L, eps_over_target=e0*L))
RES['bulk_table'] = rows
print("\nL      eps_bulk (this seat)   prove-lagrangian nu tau/r^2   ratio    (1/L)   eps*L")
for x in rows:
    print(f"{x['L']:4d}   {x['eps_bulk_exact']:.6e}          {x['eps_lagrangian']:.6e}      "
          f"{x['ratio']:.4f}   {x['target_1_over_L']:.5f}  {x['eps_over_target']:.6f}")

# ---------------------------------------------------------------- the tail term (material frame)
# Z ~ N(0, diag(2 sigma_y I_4, 2 sigma_z)) EXACTLY for the affine strain (Lemma V.3).
# sigma_y = (rho0 delta)^2 I2 / L ,  sigma_z = (rho0 delta)^2 I4 / L  (design point rho0 delta = sqrt(nu/M)).
# Along a unit direction n = (sin phi0 in the y-radial slot, cos phi0 in z):
#   Var(Z.n) = 2 sigma_y sin^2 phi0 + 2 sigma_z cos^2 phi0 ,  P(Z.n >= d) <= (1/2) exp(-d^2/(2 Var)).
phi0 = math.radians(30.0); s2 = math.sin(phi0)**2; c2 = math.cos(phi0)**2
def var_n(L, rho0delta=1.0):
    return 2*rho0delta**2*(RES['I2']*s2 + RES['I4']*c2)/L
def tail_edge(f, L):
    """inner edge at rho = rho0, tracked point at rho = (1+f) rho0; d = f rho0 = (f/delta) rho0delta."""
    d = f/sd                                   # in units of rho0 delta
    return 0.5*math.exp(-d**2/(2*var_n(L)))
def amp(f, phi0_deg=30.0):
    """A = ||eta_0||_inf/|eta_0(x_*)| + 1 = (1+f) sin phi0 / sin delta + 1"""
    return (1+f)*math.sin(math.radians(phi0_deg))/sd + 1.0
rows = []
for L in [10, 20, 40, 80, 160]:
    eb = bulk_eps(L)
    lo, hi = 0.0, 5.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if amp(mid)*tail_edge(mid, L) <= eb: hi = mid
        else: lo = mid
    f = hi
    rows.append(dict(L=L, eps_bulk=eb, f_required=f, d_over_rho0delta=f/sd,
                     A=amp(f), tail=amp(f)*tail_edge(f, L),
                     c2_inflation=L/(L-math.log(1+f)),
                     c2=8*(1-math.sqrt(2/3))*L/(L-math.log(1+f))))
RES['inset_table_affine'] = rows
print("\ninset needed so that (edge tail) <= (bulk term), AFFINE strain, exact Gaussian:")
print("L     eps_bulk     f        d/(rho0 delta)   A       tail       c2 infl   c2")
for x in rows:
    print(f"{x['L']:4d}  {x['eps_bulk']:.4e}  {x['f_required']:.5f}   {x['d_over_rho0delta']:7.3f}"
          f"  {x['A']:6.3f}  {x['tail']:.3e}  {x['c2_inflation']:.5f}  {x['c2']:.5f}")
# the brief's d = rho0 delta (f = delta) for reference
rows = []
for L in [5, 10, 20, 40, 80, 160]:
    d = 1.0                                    # one dissipation length, in units of rho0 delta
    P = 0.5*math.exp(-d**2/(2*var_n(L)))
    rows.append(dict(L=L, P=P, A=amp(sd), tail=amp(sd)*P, eps_bulk=bulk_eps(L), target=1.0/L))
RES['brief_one_dissipation_length'] = rows
print("\nthe brief's d = rho0 delta (f = delta = 0.1309):")
print("L      P(edge)      A*P        eps_bulk     1/L")
for x in rows:
    print(f"{x['L']:4d}  {x['P']:.4e}  {x['tail']:.4e}  {x['eps_bulk']:.4e}  {x['target']:.5f}")

# smallest L at which one dissipation length suffices, against the two targets
def AP(L): return amp(sd)*0.5*math.exp(-1.0/(2*var_n(L)))
for name, tgt in [('vs 1/L', lambda L: 1.0/L), ('vs eps_bulk', lambda L: bulk_eps(L))]:
    lo, hi = 1.0, 400.0
    for _ in range(300):
        mid = 0.5*(lo+hi)
        if AP(mid) <= tgt(mid): hi = mid
        else: lo = mid
    RES[f'L_min_one_dissipation_length_{name.replace(" ","_").replace("/","")}'] = hi
    print(f"smallest L with A*P(d=rho0 delta) <= target ({name}):  L = {hi:.3f}")

# ---------------------------------------------------------------- the Hessian (grad^2 u) term
# ratio  E_hess/eps_bulk = (1/2) e^c n (e^{2c}-1)/c * K2 tau r0    (PROOF.md Thm V.4, sec 5)
n = 5
def hess_ratio(K2_over_M_over_rho0, L, f=0.0, phi0_deg=30.0, c=None):
    if c is None: c = RES['c_window']
    r0 = (1+f)*math.sin(math.radians(phi0_deg))            # in units of rho0
    tau_ML = RES['theta_max']                              # tau * M * L
    pref = 0.5*math.exp(c)*n*(math.exp(2*c)-1)/c
    # K2 tau r0 = (K2_hat M/rho0) * (theta_max/(M L)) * (r0 rho0) = K2_hat theta_max r0 / L
    return pref * K2_over_M_over_rho0 * tau_ML * r0 / L
RES['hess_prefactor'] = 0.5*math.exp(RES['c_window'])*n*(math.exp(2*RES['c_window'])-1)/RES['c_window']
rows = []
for K2h in [1.0, 3.0, 10.0]:
    row = dict(K2_hat=K2h)
    for L in [10, 40, 160, 640]:
        row[f'L{L}'] = hess_ratio(K2h, L)
    rows.append(row)
RES['hess_table'] = rows
print(f"\nprefactor (1/2) e^c n (e^2c-1)/c = {RES['hess_prefactor']:.5f}   (c = {RES['c_window']:.7f})")
print("K2 = K2hat * M/rho0 :  E_hess/eps_bulk at L =   10       40      160      640")
for r0_ in rows:
    print(f"   K2hat={r0_['K2_hat']:5.1f}                        {r0_['L10']:8.4f} {r0_['L40']:8.4f}"
          f" {r0_['L160']:8.4f} {r0_['L640']:8.4f}")
# threshold: largest K2hat for which E_hess <= eps_bulk
RES['K2hat_threshold'] = {L: 1.0/hess_ratio(1.0, L) for L in [10, 40, 160, 640]}
print("largest K2hat with E_hess <= eps_bulk :", {k: round(v,3) for k,v in RES['K2hat_threshold'].items()})

json.dump(RES, open('b5_results.json','w'), indent=1)
print("\nWROTE b5_results.json")

#!/usr/bin/env python3
"""r4 -- the (D-B) datum at the origin.  Theta(rho) = (1/2)[tanh((rho-1)/w0) - tanh((rho-R)/w1)]
with w0 = 0.25 does NOT vanish at rho = 0: Theta(0) = (1/2)(1 - tanh 4) > 0.  Then
eta_0 = omega/r with omega -> -Theta(0) tanh(sin phi/sin delta) tanh(cos phi/w) near the origin,
so along phi -> 0, eta_0 ~ -Theta(0)/(rho sin delta): unbounded.  Consequences for hk2's
(D1) '|eta| bounded', Lemma 7.1's ||eta(t)||_inf <= ||eta_0||_inf, and sec 8(d)'s 'every field
is real-analytic and the global sup is finite'.  Also: the L at which Lemma 7.2's tail becomes
small at the d* the tables use."""
import json, math, numpy as np
from scipy.stats import chi2
from rlib import *
out = {}
w0 = 0.25; sd = math.sin(math.radians(7.5)); THETA0 = 0.5*(1.0 - math.tanh(1.0/w0))
out['Theta_at_origin'] = THETA0
F = Field(7.5, 0.20, 10.0, 1.0)
rows = []
for rho in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002):
    phi = math.radians(1.0)                       # near the axis, where tanh(sin phi/sd)/r ~ 1/(rho sd)
    r, z = rho*math.sin(phi), rho*math.cos(phi)
    e, er, ez, err, erz, ezz = [float(v) for v in F.all(np.array([r]), np.array([z]))]
    rows.append(dict(rho=rho, eta=e, eta_times_rho=e*rho, hessF=float(F.hessF(np.array([r]), np.array([z]))[0])))
    print(f"  rho={rho:6.3f} phi=1deg: eta_0 = {e:12.5f}   rho*eta_0 = {e*rho:.6f}   ||Hess||_F = {rows[-1]['hessF']:.3e}")
out['near_origin'] = rows
out['predicted_rho_eta_limit'] = -THETA0/sd*math.tanh(math.cos(math.radians(1))/0.20)
print(f"  Theta(0) = {THETA0:.6e};  predicted rho*eta_0 -> -Theta(0) tanh(1/w)/sin delta = {out['predicted_rho_eta_limit']:.6f}")
# where does the majorant's global claim break: K2 ~ M Theta(0)/(sd rho^2) exceeds 1 for rho < ...
out['rho_where_origin_K2_exceeds_1'] = math.sqrt(THETA0/sd)
print(f"  origin contribution ~ Theta(0)/(sin delta rho^2) exceeds M/rho0 for rho < {out['rho_where_origin_K2_exceeds_1']:.4f}")
# global sup of ||Hess eta_0||_F EXCLUDING the origin region rho < 0.5 (the honest 'global' for the tube question)
rr = np.linspace(0.005, 3.0, 1500); zz = np.linspace(-3.0, 3.0, 3000)
RR, ZZ = np.meshgrid(rr, zz, indexing='ij'); m = np.hypot(RR, ZZ) >= 0.5
H = np.where(m, F.hessF(RR, ZZ), 0.0); G = np.where(m, F.grad(RR, ZZ), 0.0)
i = np.unravel_index(np.argmax(H), H.shape); j = np.unravel_index(np.argmax(G), G.shape)
out['sup_hessF_rho_ge_half'] = dict(value=float(H[i]), r=float(RR[i]), z=float(ZZ[i]), phi_deg=math.degrees(math.atan2(RR[i], ZZ[i])))
out['sup_grad_rho_ge_half'] = dict(value=float(G[j]), r=float(RR[j]), z=float(ZZ[j]), phi_deg=math.degrees(math.atan2(RR[j], ZZ[j])))
print(f"  sup ||Hess eta_0||_F over rho>=0.5: {H[i]:.2f} at (r,z)=({RR[i]:.3f},{ZZ[i]:.3f}), phi={out['sup_hessF_rho_ge_half']['phi_deg']:.1f} deg")
print(f"  sup |grad eta_0| over rho>=0.5: {G[j]:.3f} at phi={out['sup_grad_rho_ge_half']['phi_deg']:.1f} deg")
# Lemma 7.2 tail at the tables' d*: L needed for P(chi2_5 > d^2 L/(8 sd^2 theta)) < 1e-3 and < 1e-6
th = {}
for d in (0.1662, 0.1672, 0.2029, 0.38):
    row = {}
    for p in (1e-2, 1e-3, 1e-6):
        x = chi2.isf(p, 5); row[str(p)] = x*8*sd**2*THETA_MAX/(d*d)
    th[str(d)] = row
    print(f"  d={d}: L needed for tail < 1e-2 / 1e-3 / 1e-6 : {row['0.01']:.1f} / {row['0.001']:.1f} / {row['1e-06']:.1f}")
out['L_needed_for_tail'] = th
# the tail term versus the ball term at L=10 with the tables' d* (lam = 3/2 row: sup_B|grad eta_lam| = 10.5)
d = 0.1672; nutau = sd**2*THETA_MAX/10.0; P = float(chi2.sf(d*d/(8*nutau), 5))
out['tail_vs_ball_L10_lam1.5'] = dict(P=P, grad_tail=P*out['sup_grad_rho_ge_half']['value'], hess_tail=P*out['sup_hessF_rho_ge_half']['value'],
                                       ball_sup_grad_model=10.51, ball_sup_hess_model=71.95)
print(f"  L=10, d*=0.1672: P={P:.3f}; tail terms |grad|: {P*G[j]:.2f} (ball 10.51), ||Hess||: {P*H[i]:.1f} (ball 71.95)")
json.dump(out, open('r4_results.json','w'), indent=1); print("WROTE r4_results.json")

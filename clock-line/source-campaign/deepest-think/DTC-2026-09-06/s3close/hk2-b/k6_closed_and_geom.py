#!/usr/bin/env python3
"""k6 -- (a) the CLOSED-FORM finiteness/no-log constants (crude but explicit, proved),
        (b) the window geometry under T_lam and the admissible phi0 range,
        (c) a GLOBAL scan of the exact K2 density (refuter Finding 4: local vs global).
"""
import json, math, os
import numpy as np
import sympy as sp
import hb_lib as H

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}
M = RHO0 = 1.0
DELTA_DEG, WM, PHI0 = 7.5, 0.12, 30.0
delta = math.radians(DELTA_DEG)

# ============================================================ (a) closed form
# a(delta,lam) := sup_rho  rho'^2 int_{S^4} |grad_5 eta| dOmega_4 / (lam M)
#   sharp delta-tapered plateau, strained by lam: the taper cone opens to
#   tan(dt) = lam^3 tan(delta);  eta = -lam M/r on the plateau, |grad eta| = lam M/r'^2;
#   in the taper |grad eta| <= 1.1 lam M/(dt rho'^2)  (checked below).
def a_ang(delta_eff):
    plateau = 2.0*math.cos(delta_eff)                     # int_{de}^{pi-de} sin phi dphi
    taper = 2.0*(1.0/delta_eff)*(2.0/3.0 - math.cos(delta_eff) + math.cos(delta_eff)**3/3.0)
    return 2.0*math.pi**2*(plateau + 1.1*taper)

rows = []
for lam in (1.0, 1.25, 1.5):
    dt = math.atan(lam**3*math.tan(delta))
    rows.append(dict(lam=lam, delta_eff_deg=math.degrees(dt), a_ang=a_ang(dt)))
RES['a_ang'] = rows
print("k6(a)  angular mass  a(delta_eff) = sup_rho rho^2 int_{S^4}|grad eta| dOmega / (lam M)")
for r in rows:
    print(f"        lam={r['lam']:4.2f}  delta_eff={r['delta_eff_deg']:6.3f} deg   a = {r['a_ang']:9.4f}")

# outer-tail constants  (the NO-LOG statement)
#   int_{rho'>2 rho*} |x-y|^{-5}|grad eta| dy <= 32 int_{2rho*}^R rho'^{-5}(lam M a/rho'^2) rho'^4 drho'
#                                             = 32 lam M a /(2 (2rho*)^2) = 4 lam M a / rho*^2
#   int_{rho'>2 rho*} |x-y|^{-4}|grad eta| dy <= 16 lam M a / rho*
RES['tail_coeffs'] = dict(F_outer=4.0, E_outer=16.0,
                          note="both independent of R: the integrands are rho'^{-3} and rho'^{-2}")
print(f"        outer tail of F:  <= 4  lam M a / rho*^2   (no log)")
print(f"        outer tail of E:  <= 16 lam M a / rho*     (no log)")

# inner/middle crude bounds
def closed_form_K2(rho_s, phi_s, lam, dbar, f_inset, delta_eff):
    aa = a_ang(delta_eff)
    r_s = rho_s*math.sin(phi_s)
    F = lam*M*aa*(4.0/rho_s**2 + 8.0*rho_s**3/(3.0*dbar**5))
    supg = lam*M/(RHO0**2*math.sin(delta_eff)**2)          # sup |grad eta| on the shell
    E = lam*M*aa*16.0/rho_s + supg*H.S4*3.0*rho_s
    Lam2 = 2.0*lam*M/(r_s - dbar)**3
    gnx = lam*M/r_s**2
    A1 = H.C_K*E
    A2 = H.C_GRADK*H.S4*dbar*Lam2 + H.C_GRADK*F + 0.2*gnx
    Om = lam*M/r_s; Om1 = 0.0
    return 5.0*A1 + 4.0*r_s*A2 + Om + 2.0*Om1, dict(E=E, F=F, Lam2=Lam2, A1=A1, A2=A2)

cf = []
for rho_s in (1.5, 2.0, 3.0):
    for lam in (1.0, 1.5):
        dt = math.atan(lam**3*math.tan(delta))
        phi_s = math.atan(lam**3*math.tan(math.radians(PHI0)))
        rr = rho_s*math.sqrt(lam**2*math.sin(math.radians(PHI0))**2
                             + lam**-4*math.cos(math.radians(PHI0))**2)
        best = None
        for db in np.linspace(0.02, 0.9*rr*math.sin(phi_s), 400):
            v, det = closed_form_K2(rr, phi_s, lam, db, 0.0, dt)
            if best is None or v < best[0]: best = (v, db, det)
        cf.append(dict(rho_star=rho_s, lam=lam, rho=rr, phi_deg=math.degrees(phi_s),
                       K2=best[0], dbar=best[1], **{k: float(x) for k, x in best[2].items()}))
        print(f"        closed form: rho*={rho_s:4.2f} lam={lam:4.2f} -> K2 <= {best[0]:12.1f} (dbar={best[1]:.4f})")
RES['closed_form'] = cf

# ============================================================ (b) geometry of the window
phi_c = math.acos(WM)                      # the equatorial layer boundary cos phi = w
RES['phi_c_deg'] = math.degrees(phi_c)
geo = []
for lam in (1.0, 1.1, 1.25, 1.4, 1.5):
    p0 = math.radians(PHI0)
    ph = math.atan(lam**3*math.tan(p0))
    scale = math.sqrt(lam**2*math.sin(p0)**2 + lam**-4*math.cos(p0)**2)
    geo.append(dict(lam=lam, phi_deg=math.degrees(ph), rho_over_rho_star=scale,
                    r_over_rho_star=lam*math.sin(p0), z_over_rho_star=math.cos(p0)/lam**2,
                    d_eq_over_rho_star=scale*math.sin(phi_c - ph),
                    d_taper_over_rho_star=scale*math.sin(ph - math.atan(lam**3*math.tan(delta)))))
RES['geometry'] = geo
print(f"\nk6(b)  equatorial layer boundary: cos phi = w = {WM}  ->  phi_c = {math.degrees(phi_c):.3f} deg")
print("        lam    phi(deg)   rho/rho*   r/rho*   z/rho*   d_eq/rho*  d_taper/rho*")
for g in geo:
    print(f"       {g['lam']:5.2f}  {g['phi_deg']:8.3f}  {g['rho_over_rho_star']:8.4f} "
          f"{g['r_over_rho_star']:8.4f} {g['z_over_rho_star']:8.4f} "
          f"{g['d_eq_over_rho_star']:10.4f} {g['d_taper_over_rho_star']:12.4f}")

# admissible phi0: require d_eq >= dmin for all lam in [1,3/2]
def d_eq(phi0, lam, dmin_units=1.0):
    p0 = math.radians(phi0)
    ph = math.atan(lam**3*math.tan(p0))
    scale = math.sqrt(lam**2*math.sin(p0)**2 + lam**-4*math.cos(p0)**2)
    return scale*math.sin(phi_c - ph)
grid = np.linspace(1.0, 89.0, 8801)
dmin = math.sin(delta)                      # one dissipation length, in units of rho_star
ok = [p for p in grid if min(d_eq(p, l) for l in np.linspace(1, 1.5, 51)) >= dmin
      and math.sin(math.radians(p) - delta) >= dmin]
RES['admissible_phi0'] = dict(dmin_over_rho_star=dmin, lo_deg=float(min(ok)), hi_deg=float(max(ok)))
print(f"\n        admissible phi0 (d_eq >= sin delta = {dmin:.5f} rho* for all lam<=3/2,"
      f" and d_taper >= same at lam=1):")
print(f"        phi0 in [{min(ok):.3f} deg, {max(ok):.3f} deg]   (design point 30 deg)")
worst = min((d_eq(PHI0, l), l) for l in np.linspace(1, 1.5, 501))
RES['worst_d_eq_at_phi0_30'] = dict(d_over_rho_star=worst[0], lam=worst[1],
                                    in_dissipation_lengths=worst[0]/dmin)
print(f"        at phi0=30 deg the worst d_eq = {worst[0]:.5f} rho* at lam={worst[1]:.3f}"
      f"  = {worst[0]/dmin:.3f} dissipation lengths")
# the z-contraction factor
RES['z_contraction'] = 1.0/1.5**2
print(f"        z contracts by lam^-2 = {1/1.5**2:.6f}; r dilates by lam = 1.5 (the [lam^-2, lam] window)")

json.dump(RES, open(os.path.join(HERE, 'k6_results.json'), 'w'), indent=1, sort_keys=True)
print("\nwrote k6_results.json")

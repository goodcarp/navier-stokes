#!/usr/bin/env python3
"""k4 -- DIRECT measurement of grad^2_5 b at points of N_tau for the campaign datum,
and the confrontation with the k3 bound.

The instrument is the polar quadrature validated in k2 (7 quantities to ~1e-13 against an
exact closed form, three internal controls).  Controls re-run here on the real datum:
   C1  a_rr + 3 a_r/r + a_zz = d_z eta          (Lap_5 a = d_z eta)
   C2  d_1 d_5 a = d_5 d_1 a
   C3  the quadrature's d_2 d_2 a equals a_r / r
   C4  independence of the split radius dbar
   C5  grid refinement
"""
import json, math, os
import numpy as np
import hb_lib as H

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}
M, RHO0, DELTA_DEG, WM, PHI0 = 1.0, 1.0, 7.5, 0.12, 30.0
BND = json.load(open(os.path.join(HERE, 'k3_results.json')))['bound_table']


def traj(rho_star, phi0_deg, lam):
    p = math.radians(phi0_deg)
    return lam*rho_star*math.sin(p), rho_star*math.cos(p)/lam**2


def measure(dat, rs, zs, Q, dbar, smax):
    o = Q.eval(dat, rs, zs, dbar, smax)
    er, ez = dat.der(rs, zs)
    err, erz, ezz = dat.hess(rs, zs)
    lap = o['a_rr'] + 3*o['a_r_over_r'] + o['a_zz']
    om = float(dat.omega(rs, zs)); omr, omz = dat.omega_grad(rs, zs)
    T = H.assemble_T(rs, zs, o['a_r'], o['a_z'], o['a_rr'], o['a_rz'], o['a_zz'],
                     om, float(omr), float(omz))
    k2, e = H.K2_of_T(T)
    return dict(quad=o, lap=float(lap), dz_eta=float(ez),
                C1=abs(lap - float(ez))/max(abs(float(ez)), 1e-30),
                C2=abs(o['a_rz'] - o['a_rz_sym'])/max(abs(o['a_rz']), 1e-30),
                C3=abs(o['a_r_over_r'] - o['a_r']/rs)/max(abs(o['a_r']/rs), 1e-30),
                omega=om, om_r=float(omr), om_z=float(omz), K2_exact=float(k2),
                # the hand-derived reduction inequality, evaluated on the measured scalars
                K2_reduction=float(3*abs(o['a_r']) + 2*abs(o['a_z'])
                                   + rs*(abs(o['a_rr']) + abs(o['a_rz']) + abs(o['a_zz'])
                                         + abs(o['a_r'])/rs)
                                   + abs(om)/rs + abs(float(omr)) + abs(float(omz))))


Q = H.PolarQuad(nb=80, ng=80, ns=64, nfar=170)
rows = []
for b in BND:
    L, rho_star, lam = b['L'], b['rho_star'], b['lam']
    dat = H.make_datumB(M=M, rho0=RHO0, L=L, delta_deg=DELTA_DEG, w=WM, lam=lam)
    rs, zs = traj(rho_star, PHI0, lam)
    smax = 3.0*RHO0*math.exp(L)
    m = measure(dat, rs, zs, Q, dbar=0.45*rs, smax=smax)
    m.update(L=L, rho_star=rho_star, lam=lam, r=rs, z=zs,
             K2_bound=b['K2'], slack=b['K2']/m['K2_exact'],
             slack_reduction=b['K2']/m['K2_reduction'])
    rows.append(m)
    print(f"  L={L:5.1f} rho*={rho_star:4.2f} lam={lam:4.2f}  a_r={m['quad']['a_r']:9.4f} "
          f"a_z={m['quad']['a_z']:9.4f} a_rr={m['quad']['a_rr']:9.4f} a_zz={m['quad']['a_zz']:9.4f} "
          f"| K2_exact={m['K2_exact']:8.3f}  K2_red={m['K2_reduction']:8.3f}  "
          f"K2_bound={m['K2_bound']:8.3f}  slack={m['slack']:6.2f}  "
          f"C1={m['C1']:.1e} C2={m['C2']:.1e} C3={m['C3']:.1e}")
RES['direct'] = rows

# ------------------------------------------------------------------ C4, C5
dat = H.make_datumB(M=M, rho0=RHO0, L=10.0, delta_deg=DELTA_DEG, w=WM, lam=1.0)
rs, zs = traj(2.0, PHI0, 1.0)
c4 = []
for db in (0.2, 0.45, 0.8):
    o = Q.eval(dat, rs, zs, db*rs, 3.0*math.exp(10.0))
    c4.append({k: float(o[k]) for k in ('a', 'a_r', 'a_z', 'a_rr', 'a_rz', 'a_zz', 'a_r_over_r')})
RES['C4_dbar'] = c4
spread = {k: (max(c[k] for c in c4) - min(c[k] for c in c4))/max(abs(c4[0][k]), 1e-30)
          for k in c4[0]}
RES['C4_spread'] = spread
print("\nk4  C4  dbar/r in {0.2,0.45,0.8}: relative spread")
print("      " + "  ".join(f"{k}={v:.2e}" for k, v in spread.items()))

c5 = []
for n in (60, 80, 100):
    Qn = H.PolarQuad(nb=n, ng=n, ns=n-16, nfar=n+90)
    o = Qn.eval(dat, rs, zs, 0.45*rs, 3.0*math.exp(10.0))
    c5.append(dict(n=n, **{k: float(o[k]) for k in ('a', 'a_rr', 'a_zz', 'a_r_over_r')}))
RES['C5_grid'] = c5
print("k4  C5  grid refinement at rho*=2, lam=1, L=10")
for c in c5:
    print(f"      n={c['n']:4d} a={c['a']:.10f} a_rr={c['a_rr']:.10f} a_zz={c['a_zz']:.10f} "
          f"a_r/r={c['a_r_over_r']:.10f}")

json.dump(RES, open(os.path.join(HERE, 'k4_results.json'), 'w'), indent=1, sort_keys=True)
print("\nwrote k4_results.json")

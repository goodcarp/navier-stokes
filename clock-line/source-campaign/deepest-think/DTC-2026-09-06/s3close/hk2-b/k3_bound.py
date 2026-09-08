#!/usr/bin/env python3
"""k3 -- the (H-K2) BOUND, assembled from proved ingredients, for the campaign datum.

For x in N_tau and the field eta(.,t) (transported to strain lam):
   ||grad a||        <= C_K * E(x),               E = int |x-y|^{-4}|grad eta| dy
   ||grad^2 a||_op   <= C_gK*|S^4|*dbar*Lam2(dbar) + C_gK*F(x,dbar) + (1/5)|grad eta(x)|
   K2               <= sum_m nu_m(r,z) * |S_m|   (triangle inequality on max_e||T(.,.,e)||_op)
The dbar-minimisation is a legitimate infimum over a family of valid bounds.
"""
import json, math, os
import numpy as np
import hb_lib as H

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}
M, RHO0 = 1.0, 1.0
DELTA_DEG, WM, PHI0 = 7.5, 0.12, 30.0
Q = H.PolarQuad(nb=80, ng=80)


def sgrid(smax, n1=1500, n2=500, scut=20.0):
    a = np.exp(np.linspace(math.log(1e-5), math.log(min(scut, smax)), n1))
    if smax > scut:
        b = np.exp(np.linspace(math.log(scut), math.log(smax), n2))[1:]
        return np.concatenate([a, b])
    return a


def traj(rho_star, phi0_deg, lam):
    p = math.radians(phi0_deg)
    return lam*rho_star*math.sin(p), rho_star*math.cos(p)/lam**2


def bound_at(dat, rs, zs, smax, ngrid_s=None, nh=241):
    er, ez = dat.der(rs, zs)
    gnx = math.hypot(float(er), float(ez))
    sv = ngrid_s if ngrid_s is not None else sgrid(smax)
    g = Q.gfun(dat, rs, zs, sv)
    E, Ffun = H.EF_from_g(sv, g)
    A1 = H.C_K*E
    # cached Hessian sup on the largest disc, then max over sub-discs
    dmax = 0.97*rs
    t = np.linspace(-1.0, 1.0, nh)
    RR, ZZ = np.meshgrid(rs + dmax*t, zs + dmax*t, indexing='ij')
    D = np.sqrt((RR - rs)**2 + (ZZ - zs)**2)
    Vh = np.where(RR > 1e-9, dat.hess_opnorm(np.maximum(RR, 1e-9), ZZ), -np.inf)
    best = None
    for db in np.linspace(0.01*rs, dmax, 120):
        L2 = float(np.max(np.where(D <= db, Vh, -np.inf)))
        F = Ffun(db)
        A2 = H.C_GRADK*H.S4*db*L2 + H.C_GRADK*F + 0.2*gnx
        if best is None or A2 < best[0]:
            best = (A2, db, L2, F)
    A2, db, L2, F = best
    om = float(dat.omega(rs, zs)); omr, omz = dat.omega_grad(rs, zs)
    Om = abs(om)/rs; Om1 = max(abs(float(omr)), abs(float(omz)))
    nu = H.elementary_norms(rs, zs, ngrid=1500)
    K2 = ((nu['a_r'] + nu['a_z'])*A1
          + (nu['a_rr'] + nu['a_rz'] + nu['a_zz'] + nu['a_r_over_r'])*A2
          + nu['om_over_r']*Om + (nu['om_r'] + nu['om_z'])*Om1)
    return dict(r=rs, z=zs, gradeta=gnx, E=E, A1=A1, dbar=float(db), Lam2=L2, F=F, A2=A2,
                omega=om, Om=Om, Om1=Om1, nu=nu, K2=float(K2),
                part_a1=float((nu['a_r'] + nu['a_z'])*A1),
                part_a2=float((nu['a_rr'] + nu['a_rz'] + nu['a_zz'] + nu['a_r_over_r'])*A2),
                part_om=float(nu['om_over_r']*Om),
                part_gom=float((nu['om_r'] + nu['om_z'])*Om1))


nu0 = H.elementary_norms(math.sin(math.radians(PHI0)), math.cos(math.radians(PHI0)), ngrid=6000)
RES['nu_at_rho1_phi30'] = nu0
print("k3  elementary tensor norms nu_m at (r,z) = (sin30, cos30):")
for k, v in nu0.items(): print(f"      nu[{k:12s}] = {v:.6f}")

rows = []
for L in (10.0, 40.0):
    smax = 3.0*RHO0*math.exp(L)
    sv = sgrid(smax)
    for rho_star in (1.0, 1.5, 2.0, 3.0):
        for lam in (1.0, 1.25, 1.5):
            dat = H.make_datumB(M=M, rho0=RHO0, L=L, delta_deg=DELTA_DEG, w=WM, lam=lam)
            rs, zs = traj(rho_star, PHI0, lam)
            b = bound_at(dat, rs, zs, smax, ngrid_s=sv)
            b.update(L=L, rho_star=rho_star, lam=lam,
                     phi_deg=math.degrees(math.atan2(rs, zs)), rho=math.hypot(rs, zs))
            rows.append(b)
            print(f"  L={L:5.1f} rho*={rho_star:4.2f} lam={lam:4.2f} (r,z)=({rs:.4f},{zs:.4f}) "
                  f"phi={b['phi_deg']:5.2f} E={b['E']:8.4f} A1={b['A1']:7.4f} "
                  f"dbar={b['dbar']:5.3f} Lam2={b['Lam2']:8.3f} F={b['F']:8.3f} A2={b['A2']:8.4f} "
                  f"K2={b['K2']:8.3f}")
RES['bound_table'] = rows

pairs = {}
for r in rows:
    pairs.setdefault((r['rho_star'], r['lam']), {})[r['L']] = r['K2']
li = {f"rho*={k[0]},lam={k[1]}": dict(L10=v[10.0], L40=v[40.0], ratio=v[40.0]/v[10.0])
      for k, v in pairs.items()}
RES['L_independence'] = li
print("\nk3  L-independence of the bound (the 'no log' claim):")
for k, v in li.items():
    print(f"      {k:22s} K2(L=10)={v['L10']:8.3f}  K2(L=40)={v['L40']:8.3f}  ratio={v['ratio']:.6f}")

json.dump(RES, open(os.path.join(HERE, 'k3_results.json'), 'w'), indent=1, sort_keys=True)
print("\nwrote k3_results.json")

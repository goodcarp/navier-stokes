#!/usr/bin/env python3
"""k7 -- the REFINED (H-K2) bound and the exact measurement from ONE instrument.

The refined bound keeps the kernel's exact angular cancellation
(int_{S^4} grad K dOmega_4 = 0, proved in k1) and gives up only the cancellation in the
shell radius s:

   |d_j a(x)|      <=  int_0^smax |J_j(s)| ds          J_j = int_{S^4} K(sig)(d_j eta)(x-s sig)
   |d_k d_j a(x)|  <=  C_gK |S^4| dbar Lam2(dbar)                          (near ball, Taylor)
                       + int_dbar^smax |I_kj(s)| ds/s                      (far)
                       + (1/5) delta_{k5}|d_j eta(x)|
   K2 <= 3|a_r| + 2|a_z| + r(|a_rr|+|a_rz|+|a_zz|+|a_r|/r) + |om|/r + |d_r om| + |d_z om|

The EXACT values are the signed sums of the same arrays (near panel subtracted), so bound
and measurement share the instrument and the slack is like-for-like.
"""
import json, math, os
import numpy as np
import hb_lib as H

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}
M = RHO0 = 1.0
DELTA_DEG, WM, PHI0 = 7.5, 0.12, 30.0
Q = H.PolarQuad(nb=64, ng=64)


def sgrid(rs, smax, n1=900, n2=260, scut=25.0):
    a = np.exp(np.linspace(math.log(1e-4*rs), math.log(min(scut, smax)), n1))
    if smax > scut:
        b = np.exp(np.linspace(math.log(scut), math.log(smax), n2))[1:]
        return np.concatenate([a, b])
    return a


def cumint_log(sv, y):
    """returns c[i] = int_{sv[i]}^{sv[-1]} y ds/s  (trapezoid in log s)."""
    ls = np.log(sv)
    inc = 0.5*(y[1:] + y[:-1])*np.diff(ls)
    tail = np.concatenate([np.cumsum(inc[::-1])[::-1], [0.0]])
    return tail


def traj(rho_star, phi0_deg, lam):
    p = math.radians(phi0_deg)
    return lam*rho_star*math.sin(p), rho_star*math.cos(p)/lam**2


rows = []
for L in (10.0, 40.0):
    smax = 3.0*RHO0*math.exp(L)
    for rho_star in (1.0, 1.5, 2.0, 3.0):
        for lam in (1.0, 1.25, 1.5):
            dat = H.make_datumB(M=M, rho0=RHO0, L=L, delta_deg=DELTA_DEG, w=WM, lam=lam)
            rs, zs = traj(rho_star, PHI0, lam)
            sv = sgrid(rs, smax)
            S = Q.shell_integrals(dat, rs, zs, sv)
            ls = np.log(sv)
            # |a_j| bound and exact
            aJ = {'a_r': S['J1'], 'a_z': S['J5']}
            bnd = {}; exa = {}
            for k, y in aJ.items():
                bnd[k] = float(np.trapz(np.abs(y)*sv, ls)) + abs(y[0])*sv[0]
                exa[k] = float(np.trapz(y*sv, ls)) + y[0]*sv[0]
            # far tails of |I| and signed I, as functions of the split radius
            tails_abs = {k: cumint_log(sv, np.abs(S[k])) for k in ('I11', 'I15', 'I55', 'I22')}
            tails_sgn = {k: cumint_log(sv, S[k]) for k in ('I11', 'I15', 'I55', 'I22')}
            # near (subtracted) signed head, cumulative from 0
            def head_sgn(y):
                inc = 0.5*(y[1:] + y[:-1])*np.diff(ls)
                return np.concatenate([[0.0], np.cumsum(inc)]) + y[0]  # + O(s0) head
            heads = {k: head_sgn(S['Is' + k[1:]]) for k in ('I11', 'I15', 'I55')}
            heads['I22'] = head_sgn(S['I22'])          # e_2 component: grad eta(x)_2 = 0
            # exact values (independent of the split radius; take the midpoint index)
            imid = int(np.argmin(np.abs(sv - 0.45*rs)))
            exa['a_rr'] = float(heads['I11'][imid] + tails_sgn['I11'][imid])
            exa['a_rz'] = float(heads['I15'][imid] + tails_sgn['I15'][imid])
            exa['a_zz'] = float(heads['I55'][imid] + tails_sgn['I55'][imid] + 0.2*S['ez_x'])
            exa['a_r_over_r'] = float(heads['I22'][imid] + tails_sgn['I22'][imid])
            om = float(dat.omega(rs, zs)); og = dat.omega_grad(rs, zs)
            omr, omz = float(og[0]), float(og[1])
            # bound: minimise over the split radius
            best = None
            for i in range(len(sv)):
                db = sv[i]
                if db < 0.02*rs or db > 0.9*rs: continue
                L2 = H.Lambda2(dat, rs, zs, db, n=121)
                near = H.C_GRADK*H.S4*db*L2
                b_rr = tails_abs['I11'][i] + near
                b_rz = tails_abs['I15'][i] + near
                b_zz = tails_abs['I55'][i] + near + 0.2*abs(S['ez_x'])
                b_22 = min(tails_abs['I22'][i] + near, bnd['a_r']/rs)
                K2b = (3*bnd['a_r'] + 2*bnd['a_z'] + rs*(b_rr + b_rz + b_zz + b_22)
                       + abs(om)/rs + abs(omr) + abs(omz))
                if best is None or K2b < best[0]:
                    best = (K2b, db, L2, near, dict(a_rr=b_rr, a_rz=b_rz, a_zz=b_zz,
                                                    a_r_over_r=b_22))
            K2b, db, L2, near, bb = best
            T = H.assemble_T(rs, zs, exa['a_r'], exa['a_z'], exa['a_rr'], exa['a_rz'],
                             exa['a_zz'], om, omr, omz)
            k2ex, _ = H.K2_of_T(T)
            K2red = (3*abs(exa['a_r']) + 2*abs(exa['a_z'])
                     + rs*(abs(exa['a_rr']) + abs(exa['a_rz']) + abs(exa['a_zz'])
                           + abs(exa['a_r'])/rs) + abs(om)/rs + abs(omr) + abs(omz))
            scale = max(abs(exa['a_rr']), abs(exa['a_zz']), abs(exa['a_r_over_r']), 1e-12)
            C1 = abs(exa['a_rr'] + 3*exa['a_r_over_r'] + exa['a_zz'] - S['ez_x'])/scale
            C3 = abs(exa['a_r_over_r'] - exa['a_r']/rs)/scale
            # split-radius independence of the exact value
            i2 = int(np.argmin(np.abs(sv - 0.15*rs))); i3 = int(np.argmin(np.abs(sv - 0.8*rs)))
            v = [float(heads['I55'][j] + tails_sgn['I55'][j]) for j in (i2, imid, i3)]
            C4 = (max(v) - min(v))/max(abs(np.mean(v)), 1e-12)
            rows.append(dict(L=L, rho_star=rho_star, lam=lam, r=rs, z=zs,
                             phi_deg=math.degrees(math.atan2(rs, zs)),
                             exact=exa, bound_scalars={**{k: bnd[k] for k in bnd}, **bb},
                             dbar=float(db), Lam2=float(L2), near=float(near),
                             omega=om, om_r=omr, om_z=omz,
                             K2_exact=float(k2ex), K2_reduction=float(K2red),
                             K2_bound=float(K2b), slack_vs_exact=float(K2b/k2ex),
                             slack_vs_reduction=float(K2b/K2red),
                             C1_scaled=float(C1), C3_scaled=float(C3), C4_split=float(C4)))
            print(f"  L={L:5.1f} rho*={rho_star:4.2f} lam={lam:4.2f} phi={rows[-1]['phi_deg']:5.2f}"
                  f" | K2_exact={k2ex:7.4f} reduction={K2red:7.4f} BOUND={K2b:8.4f}"
                  f" slack={K2b/k2ex:6.2f} dbar={db:.4f} C1={C1:.1e} C3={C3:.1e} C4={C4:.1e}")
RES['rows'] = rows
RES['summary'] = dict(
    max_bound_all=max(r['K2_bound'] for r in rows),
    max_exact_all=max(r['K2_exact'] for r in rows),
    max_bound_inset15=max(r['K2_bound'] for r in rows if r['rho_star'] >= 1.5),
    max_exact_inset15=max(r['K2_exact'] for r in rows if r['rho_star'] >= 1.5),
    max_bound_inset20=max(r['K2_bound'] for r in rows if r['rho_star'] >= 2.0),
    max_exact_inset20=max(r['K2_exact'] for r in rows if r['rho_star'] >= 2.0),
    min_slack=min(r['slack_vs_exact'] for r in rows),
    max_slack=max(r['slack_vs_exact'] for r in rows),
    worst_C1=max(r['C1_scaled'] for r in rows), worst_C3=max(r['C3_scaled'] for r in rows),
    worst_C4=max(r['C4_split'] for r in rows))
print("\nk7 summary:", json.dumps(RES['summary'], indent=1))
pairs = {}
for r in rows: pairs.setdefault((r['rho_star'], r['lam']), {})[r['L']] = (r['K2_bound'], r['K2_exact'])
RES['L_independence'] = {f"rho*={k[0]},lam={k[1]}":
                         dict(bound_ratio=v[40.0][0]/v[10.0][0], exact_ratio=v[40.0][1]/v[10.0][1])
                         for k, v in pairs.items()}
print("\nk7  L-independence (L=40 / L=10):")
for k, v in RES['L_independence'].items():
    print(f"      {k:22s} bound {v['bound_ratio']:.6f}   exact {v['exact_ratio']:.6f}")
json.dump(RES, open(os.path.join(HERE, 'k7_results.json'), 'w'), indent=1, sort_keys=True)
print("\nwrote k7_results.json")

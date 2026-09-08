#!/usr/bin/env python3
"""r3 -- the refuter's experiments with the r2lib instrument.
 (a) Lambda_4, Lambda_5 and the assembled hk2 bound for (D-B) at the tracked point, lam = 1 and 3/2,
     at hk2's own d* -- my numbers against their 98.688/135.119, 66.6622 and 158.73/216.47, 161.7735.
 (b) sharp radial edge (D-A) with inset f: the TRUE ||Hess a|| and K2 against hk2's '8/f'.
 (c) the S3 datum of the brief (sharp sgn z, corner taper, tanh radial ramp): hk2's majorant bound and
     the measured K2 at lam = 1 and lam = 3/2.
 (d) delta_m: (D-B) with w = 0.20, 0.10, 0.05, 0.025 measured INSIDE the equatorial layer.
 (e) (D-B) near the origin: a(x) against -3 alpha log|x|, and |grad a| there.
"""
import json, math, time, sys
import numpy as np
from r2lib import *

out = {}
R = math.exp(10.0)
def point(phi0_deg, f, lam):
    phi0 = math.radians(phi0_deg); rho0 = 1.0 + f
    r0, z0 = rho0*math.sin(phi0), rho0*math.cos(phi0)
    return lam*r0, z0/lam**2
def dist_taper(r, z, delta_deg, lam):
    ang = math.atan(lam**3*math.tan(math.radians(delta_deg)))
    return abs(r*math.cos(ang) - z*math.sin(ang))

which = sys.argv[1:] or ['a', 'b', 'c', 'd', 'e']

if 'a' in which:
    print("=== (a) (D-B): Lambda_4, Lambda_5 and the assembled bound at hk2's d* ===")
    rows = []
    for lam, dstar in ((1.0, 0.16619394776998184), (1.5, 0.16715664936537736)):
        D = Datum(AngularDB(7.5, 0.20), RadialTanh(R), lam=lam)
        r, z = point(30.0, 0.0, lam)
        t0 = time.time()
        I = Instrument(D, r, z, nT=200, nX=150)
        m = I.measure(dstar)
        sg, sh, sb = ball_sups(D, r, z, dstar)
        e = [v[0] for v in D.bulk(np.array([r]), np.array([z]))]
        B = hk2_bound(r, dstar, sg, sh, sb, m['L4'], m['L5'], e[0], e[1], e[2])
        K2m = K2_from_fields(r, m['a_r'], m['a_z'], m['a_rr'], m['a_rz'], m['a_zz'], e[0], e[1], e[2])
        row = dict(lam=lam, r=r, z=z, d=dstar, L4=m['L4'], L5=m['L5'], L5_gradK_weighted=m['L5_gradK'],
                   sup_grad=sg, sup_hessF=sh, sup_grad_bdry=sb, bound=B, measured=dict(grad_a=m['grad_a'],
                   hess_a=m['hess_op'], K2=K2m, a=m['a']), time=time.time()-t0)
        rows.append(row)
        print(f"  lam={lam}: L4={m['L4']:.4f} L5={m['L5']:.4f}  sups=({sg:.4f},{sh:.4f},{sb:.4f})"
              f"  bound: |grad a|<={B['grad_a']:.4f} ||Hess||<={B['hess_a']:.4f} K2<={B['K2']:.4f}"
              f"  measured: |grad a|={m['grad_a']:.4f} ||Hess||={m['hess_op']:.4f} K2={K2m:.4f}  a={m['a']:.4f}  [{row['time']:.0f}s]")
        print(f"       with |gradK| weighted instead of C_gradK |w|^-5: C_gradK*L5 = {C_GK*m['L5']:.4f} vs int|gradK||grad eta| = {m['L5_gradK']:.4f}")
    out['a_DB'] = rows

if 'b' in which:
    print("\n=== (b) sharp radial edge (D-A), inset f: TRUE Hess a and K2 vs hk2's majorant ===")
    rows = []
    D = Datum(AngularSharp(7.5), RadialSharp(R))
    for f in (0.4, 0.2, 0.1, 0.05, 0.025, 0.0125):
        r, z = point(30.0, f, 1.0)
        d = 0.5*min(f, z, dist_taper(r, z, 7.5, 1.0))
        t0 = time.time()
        I = Instrument(D, r, z, nT=200, nX=150)
        m = I.measure(d)
        e = [v[0] for v in D.bulk(np.array([r]), np.array([z]))]
        K2m = K2_from_fields(r, m['a_r'], m['a_z'], m['a_rr'], m['a_rz'], m['a_zz'], e[0], e[1], e[2])
        # hk2's assembled majorant at the same d (plateau sups exact on the ball)
        Rm = r - d
        B = hk2_bound(r, d, 1/Rm**2, math.sqrt(7)/Rm**3, 1/Rm**2, m['L4'], m['L5'], e[0], e[1], e[2])
        row = dict(f=f, r=r, z=z, d=d, a=m['a'], grad_a=m['grad_a'], hess_a=m['hess_op'], K2=K2m, K2_times_f=K2m*f,
                   hess_asym=m['hess_asym'], mu_int=m['mu_integrated'], mu_ar=m['mu_from_a_r'], L4=m['L4'], L5=m['L5'],
                   majorant_bound=B, majorant_K2_times_f=B['K2']*f, time=time.time()-t0)
        rows.append(row)
        print(f"  f={f:7.4f}: a={m['a']:.4f} |grad a|={m['grad_a']:.4f} ||Hess a||={m['hess_op']:.4f}  K2={K2m:.4f}  K2*f={K2m*f:.4f}"
              f"   | majorant K2<={B['K2']:.2f} (x f = {B['K2']*f:.2f})  sym {m['hess_asym']:.1e} mu {m['mu_integrated']:+.4f}/{m['mu_from_a_r']:+.4f} [{row['time']:.0f}s]")
    out['b_DA_inset'] = rows

if 'c' in which:
    print("\n=== (c) the S3 datum (sharp sgn z, corner taper 7.5deg, tanh radial w0 = 0.25): bound and measurement ===")
    rows = []
    for lam in (1.0, 1.5):
        D = Datum(AngularSharp(7.5), RadialTanh(R), lam=lam)
        r, z = point(30.0, 0.0, lam)
        dmax = 0.95*min(z, dist_taper(r, z, 7.5, lam), 0.95*r)
        best = None; meas = None
        for d in np.linspace(0.05*dmax, dmax, 8):
            t0 = time.time()
            I = Instrument(D, r, z, nT=200, nX=150)
            m = I.measure(float(d))
            sg, sh, sb = ball_sups(D, r, z, float(d))
            e = [v[0] for v in D.bulk(np.array([r]), np.array([z]))]
            B = hk2_bound(r, float(d), sg, sh, sb, m['L4'], m['L5'], e[0], e[1], e[2])
            B.update(d=float(d), L4=m['L4'], L5=m['L5'], sup_grad=sg, sup_hessF=sh, sup_grad_bdry=sb)
            if meas is None:
                K2m = K2_from_fields(r, m['a_r'], m['a_z'], m['a_rr'], m['a_rz'], m['a_zz'], e[0], e[1], e[2])
                meas = dict(a=m['a'], grad_a=m['grad_a'], hess_a=m['hess_op'], K2=K2m, d=float(d), hess_asym=m['hess_asym'])
            print(f"    lam={lam} d={d:.4f}: L4={m['L4']:.3f} L5={m['L5']:.3f} sups=({sg:.3f},{sh:.3f},{sb:.3f}) K2<={B['K2']:.3f}  [{time.time()-t0:.0f}s]")
            if best is None or B['K2'] < best['K2']: best = B
        rows.append(dict(lam=lam, r=r, z=z, dmax=dmax, best_bound=best, measured=meas))
        print(f"  lam={lam}: BEST bound K2 <= {best['K2']:.4f} at d*={best['d']:.4f};  measured K2 = {meas['K2']:.4f}, a = {meas['a']:.4f}")
    out['c_S3datum'] = rows

if 'd' in which:
    print("\n=== (d) delta_m: (D-B) inside the equatorial layer, w = 0.20 ... 0.025 ===")
    rows = []
    for w in (0.20, 0.10, 0.05, 0.025):
        D = Datum(AngularDB(7.5, w), RadialTanh(R))
        for (r, z) in ((1.5, 0.0), (1.5, 0.3*w*1.5)):
            d = min(0.1, 0.5*w)
            t0 = time.time()
            I = Instrument(D, r, z, nT=320, nX=160, extra_bps=[math.hypot(r, z)])
            m = I.measure(d)
            e = [v[0] for v in D.bulk(np.array([r]), np.array([z]))]
            K2m = K2_from_fields(r, m['a_r'], m['a_z'], m['a_rr'], m['a_rz'], m['a_zz'], e[0], e[1], e[2])
            lap = m['a_rr'] + 3*m['a_r']/r + m['a_zz']
            row = dict(w=w, r=r, z=z, d=d, eta_z=e[2], grad_a=m['grad_a'], hess_a=m['hess_op'], K2=K2m,
                       lap_a=lap, pde_rel=abs(lap/e[2]-1), lower_bound=r*abs(e[2])/5 - 2*m['grad_a'],
                       hess_asym=m['hess_asym'], time=time.time()-t0)
            rows.append(row)
            print(f"  w={w:5.3f} (r,z)=({r},{z:.4f}): eta_z={e[2]:+.4f}  |grad a|={m['grad_a']:.4f} ||Hess a||={m['hess_op']:.4f}"
                  f"  K2={K2m:.4f}   [r|eta_z|/5 - 2|grad a| = {row['lower_bound']:+.4f}]  PDE rel {row['pde_rel']:.1e} sym {m['hess_asym']:.1e} [{row['time']:.0f}s]")
    out['d_delta_m'] = rows

if 'e' in which:
    print("\n=== (e) (D-B) near the origin: a(x) vs -3 alpha log|x|, and |grad a| ===")
    D = Datum(AngularDB(7.5, 0.20), RadialTanh(R))
    r1 = json.load(open('r1_results.json'))
    alpha3 = r1['a(x) ~ -3 alpha log|x| + O(1): coefficient 3 alpha']
    rows = []
    for rho in (0.3, 0.1, 0.03, 0.01, 0.003, 0.001):
        r, z = rho*math.sin(math.radians(45)), rho*math.cos(math.radians(45))
        t0 = time.time()
        I = Instrument(D, r, z, nT=240, nX=160, extra_bps=[rho, 0.5*rho, 2*rho], npan_log=24)
        F = I.bulk_far_and_grad(1e-7*rho)
        rows.append(dict(rho=rho, a=F['a'], grad_a=float(np.hypot(*F['ga'])), a_r=F['ga'][0], a_z=F['ga'][1],
                         predicted_slope_3alpha=alpha3, time=time.time()-t0))
        print(f"  |x|={rho:6.3f}: a={F['a']:.6f}  |grad a|={np.hypot(*F['ga']):.5f}   (3 alpha/|x| = {alpha3/rho:.5f})  [{time.time()-t0:.0f}s]")
    for i in range(1, len(rows)):
        sl = (rows[i]['a']-rows[i-1]['a'])/(math.log(rows[i]['rho'])-math.log(rows[i-1]['rho']))
        rows[i]['slope_da_dlogrho'] = sl
        print(f"     slope da/dlog|x| between {rows[i-1]['rho']} and {rows[i]['rho']}: {sl:+.6f}   predicted -3alpha = {-alpha3:+.6f}")
    out['e_origin'] = rows

json.dump(out, open(f"r3_results_{''.join(which)}.json", 'w'), indent=1)
print(f"\nWROTE r3_results_{''.join(which)}.json")

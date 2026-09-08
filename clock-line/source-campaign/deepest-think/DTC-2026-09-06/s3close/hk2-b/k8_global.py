#!/usr/bin/env python3
"""k8 -- the GLOBAL K2 for the campaign datum (refuter-V-b Finding 4: the coupling lemma
uses K2 globally, while V-b Sec.7 states (H-K2) only on a neighbourhood of the trajectory).

Measures the exact pointwise Lipschitz density of grad_5 b at the points where it is
largest: on the equatorial mollification layer, at the taper, near the axis, and at the
inner radial edge.  Same instrument as k7 (validated in k2).
"""
import json, math, os
import numpy as np
import hb_lib as H

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}
M = RHO0 = 1.0
DELTA_DEG, WM, L = 7.5, 0.12, 10.0
Q = H.PolarQuad(nb=64, ng=64)
smax = 3.0*RHO0*math.exp(L)


def sgrid(rs, n1=700, n2=200, scut=25.0):
    a = np.exp(np.linspace(math.log(1e-4*rs), math.log(min(scut, smax)), n1))
    b = np.exp(np.linspace(math.log(scut), math.log(smax), n2))[1:]
    return np.concatenate([a, b])


def cumint_log(sv, y):
    ls = np.log(sv); inc = 0.5*(y[1:] + y[:-1])*np.diff(ls)
    return np.concatenate([np.cumsum(inc[::-1])[::-1], [0.0]])


def at(dat, rs, zs):
    sv = sgrid(rs)
    S = Q.shell_integrals(dat, rs, zs, sv); ls = np.log(sv)
    ex = {}
    for k, y in (('a_r', S['J1']), ('a_z', S['J5'])):
        ex[k] = float(np.trapz(y*sv, ls)) + y[0]*sv[0]
    tails = {k: cumint_log(sv, S[k]) for k in ('I11', 'I15', 'I55', 'I22')}
    def head(y):
        inc = 0.5*(y[1:] + y[:-1])*np.diff(ls)
        return np.concatenate([[0.0], np.cumsum(inc)]) + y[0]
    hd = {'I11': head(S['Is11']), 'I15': head(S['Is15']), 'I55': head(S['Is55']),
          'I22': head(S['I22'])}
    i = int(np.argmin(np.abs(sv - 0.45*rs)))
    ex['a_rr'] = float(hd['I11'][i] + tails['I11'][i])
    ex['a_rz'] = float(hd['I15'][i] + tails['I15'][i])
    ex['a_zz'] = float(hd['I55'][i] + tails['I55'][i] + 0.2*S['ez_x'])
    ex['a_r_over_r'] = float(hd['I22'][i] + tails['I22'][i])
    om = float(dat.omega(rs, zs)); og = dat.omega_grad(rs, zs)
    omr, omz = float(og[0]), float(og[1])
    T = H.assemble_T(rs, zs, ex['a_r'], ex['a_z'], ex['a_rr'], ex['a_rz'], ex['a_zz'],
                     om, omr, omz)
    k2, _ = H.K2_of_T(T)
    scale = max(abs(ex['a_rr']), abs(ex['a_zz']), abs(ex['a_r_over_r']), 1e-12)
    C1 = abs(ex['a_rr'] + 3*ex['a_r_over_r'] + ex['a_zz'] - S['ez_x'])/scale
    return dict(exact=ex, omega=om, om_r=omr, om_z=omz, K2=float(k2), C1=float(C1))


pts = []
for lam in (1.0, 1.5):
    dat = H.make_datumB(M=M, rho0=RHO0, L=L, delta_deg=DELTA_DEG, w=WM, lam=lam)
    for rho in (1.2, 2.0, 4.0, 10.0):
        for phdeg in (5.0, 10.0, 20.0, 45.0, 75.0, 85.0, 89.5):
            ph = math.radians(phdeg)
            rs, zs = rho*math.sin(ph), rho*math.cos(ph)
            o = at(dat, rs, zs)
            o.update(lam=lam, rho=rho, phi_deg=phdeg, r=rs, z=zs,
                     K2_times_rho=o['K2']*rho)
            pts.append(o)
            print(f"  lam={lam:4.2f} rho={rho:5.2f} phi={phdeg:5.1f}  K2={o['K2']:9.4f}"
                  f"  K2*rho={o['K2']*rho:9.4f}  C1={o['C1']:.1e}")
RES['points'] = pts
RES['max_K2'] = max(p['K2'] for p in pts)
RES['argmax'] = max(pts, key=lambda p: p['K2'])
RES['max_K2_times_rho'] = max(p['K2_times_rho'] for p in pts)
print(f"\nk8  global sup over the scanned points: K2 = {RES['max_K2']:.4f} M/rho0 at "
      f"lam={RES['argmax']['lam']}, rho={RES['argmax']['rho']}, phi={RES['argmax']['phi_deg']} deg")
print(f"    max of K2*rho/M = {RES['max_K2_times_rho']:.4f}  (the M/rho scaling)")
json.dump(RES, open(os.path.join(HERE, 'k8_results.json'), 'w'), indent=1, sort_keys=True)
print("wrote k8_results.json")

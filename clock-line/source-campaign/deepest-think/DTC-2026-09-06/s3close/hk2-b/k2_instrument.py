#!/usr/bin/env python3
"""k2 -- VALIDATION of the polar-quadrature instrument against an exact benchmark,
plus its internal controls.  Numerics falsify; they never prove.

Benchmark: an SO(4)-invariant 5D-radial source centred on the z-axis,
    eta(x) = exp(-|x-x_c|^2/(2 sigma^2)),  x_c = (0,0,0,0,zc).
For a radial source  -Lap_5 psi = eta  gives  psi'(P) = -(1/P^4) int_0^P eta(t) t^4 dt,
hence the CLOSED FORM
    a(x) = -d_z psi = ((z-zc)/P^5) int_0^P e^{-t^2/(2 sigma^2)} t^4 dt ,  P = |x-x_c|,
whose grad a and grad^2 a are obtained by sympy differentiation.  No shared code with
PolarQuad (which never sees this formula).

Controls that must fire on the real datum too:
  C1  Lap_5 a = d_z eta                     (exact identity: a = -d_z psi, -Lap psi = eta)
  C2  symmetry  d_1 d_5 a = d_5 d_1 a       (two different kernel components)
  C3  independence of the split radius dbar (the annulus integral of grad K is exactly 0)
"""
import json, math, os
import numpy as np
import sympy as sp
import hb_lib as H

HERE = os.path.dirname(os.path.abspath(__file__))
RES = {}

# ---------------------------------------------------------------- benchmark
sig, zc = 0.7, 0.4
r, z = sp.symbols('r z', positive=True)
etaB = sp.exp(-((r**2 + (z - zc)**2)/(2*sig**2)))
dat = H._Datum(etaB, r, z, name="gauss")

P = sp.sqrt(r**2 + (z - zc)**2)
t = sp.symbols('t', positive=True)
Iexpr = sp.integrate(sp.exp(-t**2/(2*sig**2))*t**4, (t, 0, P))
aexpr = sp.simplify((z - zc)/P**5*Iexpr)
fa   = sp.lambdify((r, z), aexpr, 'numpy')
fa_r = sp.lambdify((r, z), sp.diff(aexpr, r), 'numpy')
fa_z = sp.lambdify((r, z), sp.diff(aexpr, z), 'numpy')
fa_rr = sp.lambdify((r, z), sp.diff(aexpr, r, 2), 'numpy')
fa_rz = sp.lambdify((r, z), sp.diff(aexpr, r, z), 'numpy')
fa_zz = sp.lambdify((r, z), sp.diff(aexpr, z, 2), 'numpy')

Q = H.PolarQuad(nb=90, ng=90, ns=70, nfar=150)
rows = []
for (rs, zs) in [(0.9, 1.3), (0.5, 0.866), (1.4, 0.9), (0.35, 2.0)]:
    o = Q.eval(dat, rs, zs, dbar=0.6, smax=14.0)
    ex = dict(a=float(fa(rs, zs)), a_r=float(fa_r(rs, zs)), a_z=float(fa_z(rs, zs)),
              a_rr=float(fa_rr(rs, zs)), a_rz=float(fa_rz(rs, zs)), a_zz=float(fa_zz(rs, zs)),
              a_r_over_r=float(fa_r(rs, zs))/rs)
    row = dict(r=rs, z=zs)
    for k in ex:
        row[k] = dict(quad=float(o[k]), exact=ex[k],
                      relerr=abs(o[k]-ex[k])/max(abs(ex[k]), 1e-14), abserr=abs(o[k]-ex[k]))
    row['a_rz_sym_relerr'] = abs(o['a_rz'] - o['a_rz_sym'])/max(abs(o['a_rz']), 1e-14)
    ezr = dat.f['z'](rs, zs)
    lap = o['a_rr'] + 3*o['a_r_over_r'] + o['a_zz']
    dz_eta = dat.f['rz'](rs, zs)*0.0 + sp.lambdify((r, z), sp.diff(etaB, z), 'numpy')(rs, zs)
    row['C1_lap5_a_vs_dz_eta'] = dict(lap=float(lap), dz_eta=float(dz_eta),
                                      relerr=abs(lap-dz_eta)/max(abs(dz_eta), 1e-14))
    rows.append(row)
RES['benchmark'] = rows
print("k2  BENCHMARK (exact 5D-radial Gaussian) -- relative errors of the polar quadrature")
print("     r     z      a         a_r       a_z       a_rr      a_rz      a_zz      a_r/r     C1(Lap5a=dz eta)  C2(sym)")
for w in rows:
    print(f"  {w['r']:5.2f} {w['z']:5.2f} " + " ".join(f"{w[k]['relerr']:9.2e}" for k in
          ['a','a_r','a_z','a_rr','a_rz','a_zz','a_r_over_r'])
          + f"   {w['C1_lap5_a_vs_dz_eta']['relerr']:9.2e}      {w['a_rz_sym_relerr']:8.2e}")

# ---------------------------------------------------------------- C3: dbar independence
o1 = Q.eval(dat, 0.9, 1.3, dbar=0.35, smax=14.0)
o2 = Q.eval(dat, 0.9, 1.3, dbar=0.90, smax=14.0)
o3 = Q.eval(dat, 0.9, 1.3, dbar=1.80, smax=14.0)
sp_ = {}
for k in ['a_rr', 'a_rz', 'a_zz', 'a_r_over_r']:
    v = [o1[k], o2[k], o3[k]]
    sp_[k] = dict(values=[float(x) for x in v],
                  spread=float((max(v)-min(v))/max(abs(np.mean(v)), 1e-14)))
RES['C3_dbar_independence'] = sp_
print("\nk2  C3  dbar in {0.35, 0.90, 1.80}: relative spread of grad^2 a")
for k, v in sp_.items():
    print(f"      {k:11s} {v['values'][0]: .8f} {v['values'][1]: .8f} {v['values'][2]: .8f}   spread {v['spread']:.2e}")

# ---------------------------------------------------------------- convergence in the grids
conv = []
for nb in (50, 70, 90, 110):
    Qc = H.PolarQuad(nb=nb, ng=nb, ns=nb-10, nfar=nb+50)
    o = Qc.eval(dat, 0.9, 1.3, dbar=0.6, smax=14.0)
    conv.append(dict(n=nb, a=float(o['a']), a_zz=float(o['a_zz']), a_rr=float(o['a_rr'])))
RES['convergence'] = conv
print("\nk2  grid refinement at (r,z)=(0.9,1.3)")
for c in conv:
    print(f"      n={c['n']:4d}  a={c['a']:.12f}  a_rr={c['a_rr']:.12f}  a_zz={c['a_zz']:.12f}")
print(f"      exact           a={float(fa(0.9,1.3)):.12f}  a_rr={float(fa_rr(0.9,1.3)):.12f}  a_zz={float(fa_zz(0.9,1.3)):.12f}")

json.dump(RES, open(os.path.join(HERE, 'k2_results.json'), 'w'), indent=1, sort_keys=True)
print("\nwrote k2_results.json")

#!/usr/bin/env python3
"""b3 -- (a) the EXACT bulk rate  nu Delta_5 eta_0 / eta_0  for the campaign's actual datum,
         (b) the exact affine-strain Gaussian average on that datum,
         (c) the confrontation with the measured eta-loss of the campaign runs
             (lower/prove-duhamel/d3_results.json, read-only).

The datum is re-implemented here from the formula in lower/prove-duhamel/dcommon.py's
docstring; agreement with the campaign's own array is checked in T5.
"""
import json, math, os, hashlib
import numpy as np
import sympy as sp
from scipy.special import ive
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite_e import hermegauss

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.abspath(os.path.join(HERE, "..", "..", "lower", "prove-duhamel"))
RES = {}

# ------------------------------------------------------------------ symbolic datum
r, z = sp.symbols('r z', positive=True)
def sym_datum(rho0=1.0, R=16.0, delta_deg=7.5, w=0.12, M=1.0):
    sd = sp.sin(sp.rad(sp.Float(delta_deg)))
    w0 = sp.Float(0.25)*rho0; w1 = sp.Float(0.10)*R
    s = sp.sqrt(r**2+z**2)
    sp_ = r/s; cp = z/s
    x = sp_/sd
    tox = sp.tanh(x)/x
    Th = sp.Rational(1,2)*(sp.tanh((s-rho0)/w0) - sp.tanh((s-R)/w1))
    return -(M/(s*sd))*tox*sp.tanh(cp/w)*Th

def lap5(f):
    return sp.diff(f, r, 2) + 3*sp.diff(f, r)/r + sp.diff(f, z, 2)

ETA = sym_datum()
RATIO = (lap5(ETA)/ETA)                       # = Delta_5 eta_0 / eta_0
f_ratio = sp.lambdify((r, z), RATIO, 'numpy')
f_eta   = sp.lambdify((r, z), ETA, 'numpy')

# plateau reference: Delta_5 eta_P / eta_P = -1/r^2
def bulk_ratio(rr): return -1.0/rr**2

# ------------------------------------------------------------------ T1: the exact rate on the shell
def T1():
    rows = []
    for s0 in [1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]:
        ph = math.radians(45.0)
        r0 = s0*math.sin(ph); z0 = s0*math.cos(ph)
        rat = float(f_ratio(r0, z0)); bl = bulk_ratio(r0)
        rows.append(dict(s0=s0, r0=r0, z0=z0, lap_over_eta=rat, plateau=bl,
                         enhancement=rat/bl))
    return rows
RES['T1_exact_bulk_rate'] = T1()
print("T1  s0   r0     Delta5 eta0/eta0    -1/r0^2     ratio")
for x in RES['T1_exact_bulk_rate']:
    print(f"   {x['s0']:4.2f} {x['r0']:6.4f}  {x['lap_over_eta']:14.6f} {x['plateau']:11.6f}  {x['enhancement']:8.4f}")

# ------------------------------------------------------------------ exact Gaussian average
def gauss_average(fun, r0, z0, sig_y2, sig_z2, nR=1400, nZ=241, nrad=14.0):
    sy = math.sqrt(sig_y2); sz = math.sqrt(sig_z2)
    Rlo = max(1e-10, r0-nrad*sy); Rhi = r0+nrad*sy
    xg, wg = leggauss(nR)
    Rv = 0.5*(Rhi-Rlo)*xg + 0.5*(Rhi+Rlo); wR = 0.5*(Rhi-Rlo)*wg
    pR = (Rv**2/(sig_y2*r0))*np.exp(-(Rv-r0)**2/(2*sig_y2))*ive(1, Rv*r0/sig_y2)
    hx, hw = hermegauss(nZ); Zv = z0+sz*hx; wZ = hw/math.sqrt(2*math.pi)
    RR, ZZ = np.meshgrid(Rv, Zv, indexing='ij')
    val = fun(RR, ZZ)
    return float((wR[:,None]*pR[:,None]*wZ[None,:]*val).sum()), float((wR*pR).sum())

# ------------------------------------------------------------------ the campaign run
D3 = json.load(open(os.path.join(CAMP, 'd3_results.json')))
RES['campaign_d3_sha256'] = hashlib.sha256(open(os.path.join(CAMP,'d3_results.json'),'rb').read()).hexdigest()

def run_rows(label='N4'):
    R = [x for x in D3['runs'] if x['label'] == label][0]
    nu = R['nu']; out = []
    for tr in R['tracers']:
        H = tr['hist']
        t = np.array(H['t']); rr = np.array(H['r']); zz = np.array(H['z']); ee = np.array(H['eta'])
        # measured
        meas_loss = 1.0 - ee[-1]/ee[0]
        # measured INITIAL logarithmic rate along the trajectory (first 5 samples, LSQ)
        k = min(6, len(t))
        A = np.vstack([np.ones(k), t[:k]]).T
        slope = float(np.linalg.lstsq(A, np.log(np.abs(ee[:k])), rcond=None)[0][1])
        # predicted initial rate: nu * Delta5 eta0/eta0 at the tracer's start
        pred_rate0 = nu*float(f_ratio(rr[0], zz[0]))
        plateau_rate0 = -nu/rr[0]**2
        # integrated plateau prediction  nu int dt / r(t)^2   along the MEASURED trajectory
        s_plateau = nu*float(np.trapz(1.0/rr**2, t))
        # integrated EXACT-datum prediction (frozen datum structure, moving point)
        s_datum = -nu*float(np.trapz(f_ratio(rr, zz), t))
        out.append(dict(s0=tr['s0'], r0=float(rr[0]), z0=float(zz[0]), t_end=float(t[-1]),
                        r_end=float(rr[-1]), nu=nu,
                        meas_loss=meas_loss,
                        meas_rate0=slope, pred_rate0=pred_rate0, plateau_rate0=plateau_rate0,
                        rate0_ratio_meas_over_pred=slope/pred_rate0,
                        s_plateau=s_plateau, s_datum=s_datum,
                        meas_over_plateau=meas_loss/s_plateau,
                        meas_over_datum=meas_loss/s_datum))
    return dict(nu=nu, N=R['N'], Re0=R['Re0'], grid=R['grid'], h=R['h'], tracers=out)

ALL = ['N3','N4','N5','N4-fine','N4-Re400','N4-Re25']
for lab in ALL:
    try:
        RES[f'T2_run_{lab}'] = run_rows(lab)
    except IndexError:
        pass
print("\nT2  run  s0   meas_loss   nu*int dt/r^2   exact-datum pred   meas/plateau  meas/datum   rate0 meas/pred")
for lab in ALL:
    if f'T2_run_{lab}' not in RES: continue
    for x in RES[f'T2_run_{lab}']['tracers']:
        print(f"   {lab}  {x['s0']:3.1f}  {x['meas_loss']:.6f}   {x['s_plateau']:.6f}      {x['s_datum']:.6f}"
              f"      {x['meas_over_plateau']:7.3f}    {x['meas_over_datum']:7.3f}     {x['rate0_ratio_meas_over_pred']:7.3f}")

# ------------------------------------------------------------------ T3: exact Gaussian average on the run's datum
def T3(label='N4'):
    R = [x for x in D3['runs'] if x['label'] == label][0]
    nu = R['nu']; rows = []
    for tr in R['tracers']:
        H = tr['hist']; t = np.array(H['t']); rr = np.array(H['r'])
        lam = rr/rr[0]                                   # measured stretch of the tracked shell
        sig_y = nu*float(np.trapz(lam**-2, t)); sig_z = nu*float(np.trapz(lam**4, t))
        r0 = float(rr[0]); z0 = float(H['z'][0])
        val, mass = gauss_average(f_eta, r0, z0, 2*sig_y, 2*sig_z)
        ref = float(f_eta(r0, z0))
        rows.append(dict(s0=tr['s0'], sigma_y=sig_y, sigma_z=sig_z, mass=mass,
                         s=sig_y/r0**2, exact_gauss_rel_loss=1.0-val/ref,
                         meas_loss=1.0-tr['hist']['eta'][-1]/tr['hist']['eta'][0]))
    return rows
print("\nT3  run        s0   nu       s=sigma_y/r0^2  exact-Gauss loss   measured      meas/exact")
for lab in ALL:
    if f'T2_run_{lab}' not in RES: continue
    rows = T3(lab); RES[f'T3_exact_gauss_{lab}'] = rows
    nuv = RES[f'T2_run_{lab}']['nu']
    for x in rows:
        x['meas_over_exact'] = x['meas_loss']/x['exact_gauss_rel_loss']
        print(f"    {lab:11s} {x['s0']:3.1f}  {nuv:7.5f}  {x['s']:.5e}    {x['exact_gauss_rel_loss']:.5e}"
              f"   {x['meas_loss']:.5e}   {x['meas_over_exact']:7.3f}")

# ---------------------------------------------------------------- T4: is the residual numerical?
print("\nT4  absolute EXCESS  (measured - exact affine-Gaussian)  -- a discretisation error is")
print("    nu-independent at fixed h and shrinks with h; a viscous effect does neither.")
print("    run          h        nu       s0    measured    exact      excess")
T4 = []
for lab in ALL:
    if f'T3_exact_gauss_{lab}' not in RES: continue
    hh = RES[f'T2_run_{lab}']['h']; nuv = RES[f'T2_run_{lab}']['nu']
    for x in RES[f'T3_exact_gauss_{lab}']:
        exc = x['meas_loss'] - x['exact_gauss_rel_loss']
        T4.append(dict(run=lab, h=hh, nu=nuv, s0=x['s0'], meas=x['meas_loss'],
                       exact=x['exact_gauss_rel_loss'], excess=exc))
        print(f"    {lab:11s} {hh:7.5f} {nuv:8.5f} {x['s0']:4.1f}  {x['meas_loss']:.5e} "
              f"{x['exact_gauss_rel_loss']:.5e} {exc: .5e}")
RES['T4_excess'] = T4
def exc(lab, s0): return [x['excess'] for x in T4 if x['run'] == lab and x['s0'] == s0][0]
RES['T4_ratios'] = dict(
    h_refinement_1p5 = exc('N4', 1.5)/exc('N4-fine', 1.5),
    h_refinement_2p0 = exc('N4', 2.0)/exc('N4-fine', 2.0),
    nu_16x_1p5_excess = [exc('N4-Re25', 1.5), exc('N4', 1.5), exc('N4-Re400', 1.5)],
    nu_16x_1p5_meas   = [ [x['meas'] for x in T4 if x['run']==r and x['s0']==1.5][0]
                          for r in ('N4-Re25','N4','N4-Re400')],
    nu_values         = [0.04, 0.01, 0.0025])
print("\n    excess(N4)/excess(N4-fine) at s0=1.5 :", round(RES['T4_ratios']['h_refinement_1p5'],4),
      " at s0=2.0 :", round(RES['T4_ratios']['h_refinement_2p0'],4))
print("    excess at s0=1.5 over a 16x range in nu (nu = 0.04, 0.01, 0.0025):",
      [round(v,5) for v in RES['T4_ratios']['nu_16x_1p5_excess']])

# ---------------------------------------------------------------- T5: my datum == the campaign's
# The campaign's build_taper is IMPORTED read-only and compared with the analytic formula used
# above.  Their array is (i) sub-cell averaged (sub=3) and (ii) normalised so max|omega|=1.
import sys
sys.path.insert(0, CAMP)
import dcommon
gg, ea, RR_, LL_ = dcommon.build_taper(4, 0.0625, delta_deg=7.5, w=0.12)
Rc, Zc = gg.cell_mesh()
mine = f_eta(Rc, Zc)
norm = float(np.max(np.abs(mine*Rc)))
mine = mine/norm
sel = (np.sqrt(Rc**2+Zc**2) > 1.5) & (np.sqrt(Rc**2+Zc**2) < 8.0) & (Rc > 0.3) & (Zc > 0.3)
rel = float(np.max(np.abs(mine[sel]-ea[sel]))/np.max(np.abs(ea[sel])))
RES['T5_datum_agreement'] = dict(solver_sha256=dcommon.solver_sha(), grid=[gg.Nr, gg.Nz],
                                 max_rel_diff_bulk=rel, npts=int(sel.sum()))
print(f"\nT5  my analytic datum vs the campaign's build_taper array (h=0.0625, {int(sel.sum())} bulk cells):"
      f"  max rel diff = {rel:.3e}")

json.dump(RES, open('b3_results.json','w'), indent=1)
print("\nWROTE b3_results.json")

#!/usr/bin/env python3
"""r2 -- the majorant at lambda in {1.25, 1.5} with the inputs that hk2's OWN Lemma 7.2 licenses
for the field at time t, instead of the pure-strain model field eta_lam:
   |grad eta(x',t)|     <= e^{int Gamma} |grad eta_0|(T^{-1}x')        (7.2, no diffusion)
   ||Hess eta(x',t)||_F <= e^{2 int Gamma} ||Hess eta_0||_F(T^{-1}x')  (7.2, feedback term dropped)
with e^{int Gamma} = lam^2 (V-b's Gamma = 2a(0,t), i.e. 9/4 at lam = 3/2).  The pure-strain
model instead uses grad eta_lam = T^{-T} grad eta_0, whose norm is <= lam^2 |grad eta_0| with
equality only for a purely z-directed gradient.  Also: the anisotropy ratio, and the Theorem-Gamma
factor e^{C' theta_max / L} that the record's proved Gamma carries in addition."""
import json, math, time, numpy as np
from rlib import *
Q = dict(ns=200, nT=100, nc=100)
out = {'quad': Q}; t0 = time.time()
base = Field(7.5, 0.20, 10.0, 1.0)
for lam in (1.25, 1.5):
    F = Field(7.5, 0.20, 10.0, lam); r, z = strain_point(30.0, lam)
    D = distances(r, z, lam); dmax = 0.95*min(D['d_eq'], D['d_tap'], 0.95*r)
    fg = lam**2; fh = lam**4
    P = PushedField(base, lam, fg, fh)
    # model (pure strain) bound -- same as hk2/r1
    bm = bound_at(F, r, z, dmax, nd=8, **Q)
    # Lemma-7.2-licensed bound: TV density, ball sups, boundary sup, local grad all from P
    bp = bound_at(P, r, z, dmax, nd=8, **Q)
    # anisotropy ratios on the ball used by the model
    d = bm['d']
    g = np.linspace(-d, d, 201); RR, ZZ = np.meshgrid(r+g, z+g, indexing='ij')
    m = (RR-r)**2 + (ZZ-z)**2 <= d*d; RR = np.where(m, RR, r); ZZ = np.where(m, ZZ, z)
    ratio_grad = float(np.max(F.grad(RR, ZZ))/np.max(P.grad(RR, ZZ)))
    ratio_hess = float(np.max(F.hessF(RR, ZZ))/np.max(P.hessF(RR, ZZ)))
    row = dict(lam=lam, model=bm, lemma72=bp, ratio_K2=bp['K2']/bm['K2'],
               ratio_L5=bp['L5']/bm['L5'], ratio_L4=bp['L4']/bm['L4'],
               ratio_supgrad_model_over_licensed=ratio_grad, ratio_suphess_model_over_licensed=ratio_hess)
    out[str(lam)] = row
    print(f"lam={lam}: model K2hat<={bm['K2']:.4f} (d*={bm['d']:.4f});  Lemma-7.2-licensed K2hat<={bp['K2']:.4f} (d*={bp['d']:.4f});"
          f"  ratio {row['ratio_K2']:.3f};  L5 ratio {row['ratio_L5']:.3f};  L4 ratio {row['ratio_L4']:.3f}")
    print(f"        on the model ball: sup|grad eta_lam| / (lam^2 sup|grad eta_0|(label)) = {ratio_grad:.4f};"
          f"  Hess: {ratio_hess:.4f}")
# Theorem Gamma's proved constant, on top: e^{int Gamma} = lam^2 e^{C' theta_max/L}
tg = {}
for Cp in (76.11, 119.33, 151.15):
    tg[str(Cp)] = {str(L): math.exp(Cp*THETA_MAX/L) for L in (10, 40, 160, 640, 2560)}
out['theorem_gamma_extra_factor'] = tg
print("extra factor e^{C' theta_max/L} on the gradient (squared on the Hessian):")
for Cp, row in tg.items():
    print(f"  C'={Cp}: " + "  ".join(f"L={L}: {v:.4g}" for L, v in row.items()))
out['seconds'] = time.time()-t0
json.dump(out, open('r2_results.json','w'), indent=1)
print("WROTE r2_results.json", out['seconds'])

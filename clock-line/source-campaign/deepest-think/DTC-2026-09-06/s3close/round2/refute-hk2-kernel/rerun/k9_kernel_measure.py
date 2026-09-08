#!/usr/bin/env python3
"""k9 -- a THIRD instrument: grad a and Hess a measured by direct 5-D kernel quadrature of the
signed representation (4.1)-(4.2), with the derivatives carried on eta.  Unlike k3 this keeps
the signs (no majorant), and unlike k4 it never expands in zonal harmonics.  It works for the
STRAINED field eta_lam = eta_0 o T_lam^{-1}, which the zonal instrument cannot do cheaply.

Coordinates: w = s uhat, uhat = (sin(Theta) nhat, cos(Theta)), nhat in S^3 at angle chi to yhat,
dw = s^4 ds sin^3(Theta) dTheta 4 pi sin^2(chi) dchi.  Write x = (r ehat_r, z), x' = x - w, so
  r' = |y - s sinTheta nhat| = sqrt(r^2 + s^2 sin^2Theta - 2 r s sinTheta cos chi),  z' = z - s cosTheta,
  c  = yhat' . yhat = (r - s sinTheta cos chi)/r' ,  uhat_r = sinTheta cos chi , uhat_z = cosTheta.
Projections onto the frame at x (perpendicular components integrate to zero over S^3):
  (grad eta)_r = eta_r c ,  (grad eta)_z = eta_z
  (Hess eta)_rr = eta_rr c^2 + (eta_r/r')(1-c^2) ,  (Hess eta)_rz = eta_rz c ,  (Hess eta)_zz = eta_zz
Then, with P := 3/(8 pi^2),
  grad a       = P int ds int uhat_z (grad eta)(x-s uhat) dOmega
  Hess a       = P[ int_0^d ds int uhat_z (Hess eta) dOmega                       (near)
                  + int_d^inf (ds/s) int (delta_jz - 5 uhat_z uhat_j)(grad eta)_k dOmega   (far)
                  + int_{S^4} uhat_z uhat_j (grad eta)_k(x - d uhat) dOmega ]     (surface)
"""
import json, math
import numpy as np
from hk2lib import DatumB, StrainedB, geometry
from k4_direct import Series, K2_exact

P = 3.0/(8*math.pi**2)
def nodes(ns, nT, nc, smin, smax, log=True):
    xg, wg = np.polynomial.legendre.leggauss(ns)
    if log:
        u = 0.5*(math.log(smax)-math.log(smin))*xg + 0.5*(math.log(smax)+math.log(smin))
        s = np.exp(u); ws = 0.5*(math.log(smax)-math.log(smin))*wg*s
    else:
        s = 0.5*(smax-smin)*xg + 0.5*(smax+smin); ws = 0.5*(smax-smin)*wg
    T, wT = np.polynomial.legendre.leggauss(nT); T = 0.5*math.pi*(T+1); wT = 0.5*math.pi*wT
    X, wX = np.polynomial.legendre.leggauss(nc); X = 0.5*math.pi*(X+1); wX = 0.5*math.pi*wX
    return s, ws, T, wT, X, wX

def fields(D, r, z, s, T, X):
    sT = np.sin(T)[:,None]; cT = np.cos(T)[:,None]; cX = np.cos(X)[None,:]
    rp = np.sqrt(np.maximum(r*r + (s*sT)**2 - 2*r*s*sT*cX, 1e-300))
    zp = z - s*cT*np.ones_like(cX)
    e, er, ez, err, erz, ezz = D.eta_rz(rp, zp)
    c = (r - s*sT*cX)/rp
    return rp, er, ez, err, erz, ezz, c, sT*np.ones_like(cX), cT*np.ones_like(cX), cX*np.ones_like(sT)

def grad_a(D, r, z, smax, ns=500, nT=140, nc=140):
    s, ws, T, wT, X, wX = nodes(ns, nT, nc, 1e-7, smax)
    WA = (wT*np.sin(T)**3)[:,None]*(wX*4*math.pi*np.sin(X)**2)[None,:]
    ar = az = 0.0
    for si, wsi in zip(s, ws):
        rp, er, ez, *_ , c, sT, cT, cX = fields(D, r, z, si, T, X)
        ar += wsi*np.sum(WA*cT*er*c); az += wsi*np.sum(WA*cT*ez)
    return P*ar, P*az

def hess_a(D, r, z, d, smax, ns=500, nT=140, nc=140, nsn=200):
    Hrr = Hrz = Hzr = Hzz = 0.0
    # near: int_0^d ds int uhat_z (Hess eta) dOmega
    s, ws, T, wT, X, wX = nodes(nsn, nT, nc, 0.0, d, log=False)
    WA = (wT*np.sin(T)**3)[:,None]*(wX*4*math.pi*np.sin(X)**2)[None,:]
    for si, wsi in zip(s, ws):
        rp, er, ez, err, erz, ezz, c, sT, cT, cX = fields(D, r, z, si, T, X)
        Hrr += wsi*np.sum(WA*cT*(err*c*c + (er/rp)*(1-c*c)))
        Hrz += wsi*np.sum(WA*cT*(erz*c)); Hzr = Hrz
        Hzz += wsi*np.sum(WA*cT*ezz)
    # far: int_d^inf (ds/s) int (delta_jz - 5 uz uj) g_k dOmega
    s, ws, T, wT, X, wX = nodes(ns, nT, nc, d, smax)
    WA = (wT*np.sin(T)**3)[:,None]*(wX*4*math.pi*np.sin(X)**2)[None,:]
    Frr = Frz = Fzr = Fzz = 0.0
    for si, wsi in zip(s, ws):
        rp, er, ez, err, erz, ezz, c, sT, cT, cX = fields(D, r, z, si, T, X)
        ur = sT*cX; uz = cT
        gr = er*c; gz = ez
        Frr += (wsi/si)*np.sum(WA*(-5*uz*ur)*gr)
        Frz += (wsi/si)*np.sum(WA*(-5*uz*ur)*gz)
        Fzr += (wsi/si)*np.sum(WA*(1-5*uz*uz)*gr)
        Fzz += (wsi/si)*np.sum(WA*(1-5*uz*uz)*gz)
    # surface at |w| = d
    rp, er, ez, err, erz, ezz, c, sT, cT, cX = fields(D, r, z, d, T, X)
    ur = sT*cX; uz = cT; gr = er*c; gz = ez
    Srr = np.sum(WA*uz*ur*gr); Srz = np.sum(WA*uz*ur*gz)
    Szr = np.sum(WA*uz*uz*gr); Szz = np.sum(WA*uz*uz*gz)
    return P*np.array([[Hrr+Frr+Srr, Hrz+Frz+Srz],[Hzr+Fzr+Szr, Hzz+Fzz+Szz]])

if __name__ == "__main__":
    out = {}
    DB = DatumB(delta_deg=7.5, w=0.20, L=10.0)
    S  = Series(DB, LMAX=161)
    print("=== control: k9 (signed kernel) vs k4 (zonal series) at lam = 1 ===")
    g = geometry(30.0, 0.0, 7.5, 1.0); r, z = g['r'], g['z']
    A = S.a_field(r, z)
    ar, az = grad_a(DB, r, z, 3*DB.R)
    H = hess_a(DB, r, z, 0.25, 3*DB.R)
    print(f"   a_r  : series {A[1]:+.9f}   kernel {ar:+.9f}   rel {abs(ar/A[1]-1):.2e}")
    print(f"   a_z  : series {A[2]:+.9f}   kernel {az:+.9f}   rel {abs(az/A[2]-1):.2e}")
    print(f"   a_rr : series {A[3]:+.9f}   kernel {H[0,0]:+.9f}   rel {abs(H[0,0]/A[3]-1):.2e}")
    print(f"   a_rz : series {A[4]:+.9f}   kernel {H[0,1]:+.9f}   rel {abs(H[0,1]/A[4]-1):.2e}")
    print(f"        :                       kernel(zr) {H[1,0]:+.9f}  (symmetry control)")
    print(f"   a_zz : series {A[5]:+.9f}   kernel {H[1,1]:+.9f}   rel {abs(H[1,1]/A[5]-1):.2e}")
    out['control_lam1'] = dict(series=[A[1],A[2],A[3],A[4],A[5]],
                               kernel=[ar,az,H[0,0],H[0,1],H[1,1]], kernel_zr=H[1,0])
    out['control_worst_rel'] = max(abs(H[0,0]/A[3]-1), abs(H[0,1]/A[4]-1), abs(H[1,1]/A[5]-1),
                                   abs(ar/A[1]-1), abs(az/A[2]-1))
    out['symmetry_rel'] = abs(H[1,0]/H[0,1]-1)
    print(f"   worst relative disagreement {out['control_worst_rel']:.2e} ;"
          f" Hess symmetry {out['symmetry_rel']:.2e}")

    print("\n=== MEASURED K2hat over the window, on the STRAINED field eta_lam ===")
    rows = []
    for lam in (1.0, 1.1, 1.25, 1.5):
        DL = StrainedB(DB, lam) if lam != 1.0 else DB
        gg = geometry(30.0, 0.0, 7.5, lam); rr, zz = gg['r'], gg['z']
        d = 0.5*min(gg['d_eq'], gg['d_taper'])
        ar, az = grad_a(DL, rr, zz, 3*DB.R)
        H = hess_a(DL, rr, zz, d, 3*DB.R)
        e = [v[0] for v in DL.eta_rz(np.array([rr]), np.array([zz]))]
        k2, _ = K2_exact(rr, ar, az, H[0,0], 0.5*(H[0,1]+H[1,0]), H[1,1], e[0], e[1], e[2])
        ha = max(abs(np.linalg.eigvalsh(0.5*(H+H.T))).max(), abs(ar/rr))
        rows.append(dict(lam=lam, r=rr, z=zz, d=d, a_r=ar, a_z=az, a_rr=H[0,0],
                         a_rz=0.5*(H[0,1]+H[1,0]), a_zz=H[1,1], grad_a=math.hypot(ar,az),
                         hess_a=ha, K2=k2))
        print(f"   lam={lam:4.2f}: |grad a| = {math.hypot(ar,az):.5f}   ||Hess a|| = {ha:.5f}"
              f"    K2hat = {k2:.5f}")
    out['measured_strained'] = rows
    out['K2hat_max_window'] = max(x['K2'] for x in rows)
    print(f"   sup over the window: K2hat = {out['K2hat_max_window']:.5f}")
    json.dump(out, open('k9_results.json','w'), indent=1)
    print("\nWROTE k9_results.json")

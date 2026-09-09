#!/usr/bin/env python3
"""k3 -- the kernel integrals Lambda_4, Lambda_5 and the assembled bound on K2.

  Lambda_p(x,d) := int_{|x-x'|>d} |x-x'|^{-p} d|D eta|(x')      (a majorant integral)

evaluated in 5D spherical coordinates CENTRED AT x (so the excluded ball is a coordinate
domain):  w = s uhat,  uhat = (sin(Theta) nhat, cos(Theta)), nhat in S^3 at angle chi to yhat,
  dw = s^4 ds sin^3(Theta) dTheta * 4 pi sin^2(chi) dchi ,   |S^4| = 8pi^2/3 = 4pi * 2/3 * ... 
  (check: int sin^3 = 4/3, int 4pi sin^2 = 2pi^2  ->  4/3 * 2 pi^2 = 8 pi^2/3  OK)
and  r' = sqrt(r^2 + s^2 sin^2Theta - 2 r s sinTheta cos chi),  z' = z - s cosTheta.

Surface (jump) parts of |D eta| for the idealised datum (D-A) are integrated in SOURCE
coordinates instead, where the jump sets are coordinate surfaces.

Assembly (proved in PROOF.md sec 4):
  |grad a|      <= d * sup_B|grad eta| + C_K Lambda_4
  ||Hess a||_op <= d * sup_B||Hess eta||_F + sup_{dB}|grad eta| + C_gradK Lambda_5
with C_K |S^4| = 1 exactly and C_K = 3/(8 pi^2), C_gradK = 3/(2 pi^2).
"""
import json, math
import numpy as np
from hk2lib import DatumA, DatumB, StrainedB, geometry

C_K     = 3.0/(8*math.pi**2)
C_gradK = 3.0/(2*math.pi**2)
S4      = 8*math.pi**2/3.0
RES = {'C_K': C_K, 'C_gradK': C_gradK, 'S4': S4, 'C_K_times_S4': C_K*S4,
       'C_gradK_times_S4': C_gradK*S4}

from functools import lru_cache
@lru_cache(maxsize=None)
def _lg(n):
    return np.polynomial.legendre.leggauss(n)      # numpy's leggauss is O(n^3): cache it
def gauss(n, a, b):
    x, w = _lg(n)
    return 0.5*(b-a)*x + 0.5*(b+a), 0.5*(b-a)*w

def lambda_bulk(gradfun, x_rz, d, smax, ns=400, nT=160, nc=160, plist=(4,5)):
    """{p: int_{|w|>d} |w|^{-p} |grad eta|(x-w) dw} -- one angular pass per s, reused for all p."""
    r, z = x_rz
    u, wu = gauss(ns, math.log(d), math.log(smax)); s = np.exp(u); ws = wu*s      # ds = s du
    T, wT = gauss(nT, 0.0, math.pi)
    X, wX = gauss(nc, 0.0, math.pi)
    sT = np.sin(T); cT = np.cos(T); cX = np.cos(X)
    WA = (wT*sT**3)[:,None] * (wX*4*math.pi*np.sin(X)**2)[None,:]
    A = np.empty(ns)
    ST = sT[:,None]; CT = cT[:,None]
    for i, si in enumerate(s):
        rp = np.sqrt(np.maximum(r*r + (si*ST)**2 - 2*r*si*ST*cX[None,:], 1e-300))
        zp = z - si*CT
        A[i] = np.sum(WA*gradfun(rp, zp + 0*cX[None,:]))
    return {p: float(np.sum(ws*s**(4-p)*A)) for p in plist}

def lambda_p_surface_DA(D, x_rz, p, nr=1200, nt=800, nc=400):
    """surface (jump) parts of |D eta| for (D-A): equatorial plane {z'=0} with density
    2M h(pi/2)/r' = 2M/r' ; inner sphere {rho'=1} and outer sphere {rho'=R} with density
    M h(phi')/r'.  Integrated in source coordinates."""
    r, z = x_rz; out = {}
    X, wX = gauss(nc, 0.0, math.pi); sX2 = 4*math.pi*np.sin(X)**2
    # --- equatorial plane: dS = d^4 y' = 4 pi r'^3 sin^2(chi) dr' dchi ; density 2M/r'
    lr, wlr = gauss(nr, math.log(1.0), math.log(D.R))   # {z'=0} cap of the support shell
    rp = np.exp(lr); wrp = wlr*rp
    W2 = np.abs(r*r + rp[:,None]**2 - 2*r*rp[:,None]*np.cos(X)[None,:] + z*z)
    val = np.sum((wrp[:,None]*(wX*sX2)[None,:]) * (2.0*rp[:,None]**2) * W2**(-p/2.0))
    out['equator'] = float(val)
    # --- inner sphere rho'=1: dS = 4 pi (1-t'^2) sin^2(chi) dt' dchi ; density M h(phi')/(sqrt(1-t'^2))
    for lbl, rad in (('inner', 1.0), ('outer', D.R)):
        T, wT = gauss(nt, -1.0, 1.0)
        phi = np.arccos(np.clip(T, -1, 1)); h = D.h(phi); sp = np.sqrt(np.maximum(1-T*T, 1e-300))
        rp2 = rad*sp; zp2 = rad*T
        W2 = np.abs(r*r + rp2[:,None]**2 - 2*r*rp2[:,None]*np.cos(X)[None,:] + (z-zp2[:,None])**2)
        dens = (h/sp)*rad**4*sp**2*(1.0/rad)      # M h/r' * rho'^4 * (1-t'^2) with r'=rad*sp
        val = np.sum((wT*dens)[:,None]*(wX*sX2)[None,:] * W2**(-p/2.0))
        out[lbl] = float(val)
    return out

def assemble(sup_grad_B, sup_hess_B, sup_grad_dB, d, L4, L5, r, eta, eta_r, eta_z):
    """PROOF.md sec 4.  Returns dict with the |grad a|, ||Hess a||_op and K2 bounds."""
    ga = d*sup_grad_B + C_K*L4                        # C_K |S^4| = 1 folded into the near term
    Ha = d*sup_hess_B + sup_grad_dB + C_gradK*L5
    q  = ga + abs(eta)                                # |q| = |a_z - eta| <= |grad a| + |eta|
    qr = Ha + abs(eta_r); qz = Ha + abs(eta_z)
    S55 = 2*ga + r*Ha                                 # |S55| = |2 a_z + r a_rz|
    Huz = max(q, q + r*qr + r*qz, r*qz + S55)         # symmetric-matrix row-sum bound
    K2 = r*Ha + 2*ga + Huz
    return dict(grad_a=ga, hess_a=Ha, hess_uz=Huz, K2=K2)

def sups_on_ball_B(D, r0, z0, d, n=201):
    g = np.linspace(-d, d, n)
    RR, ZZ = np.meshgrid(r0+g, z0+g, indexing='ij')
    m = (RR-r0)**2 + (ZZ-z0)**2 <= d*d
    RR = np.where(m, RR, r0); ZZ = np.where(m, ZZ, z0); RR = np.maximum(RR, 1e-6)
    gn = D.grad_norm(RR, ZZ); hf = D.hess_F(RR, ZZ)
    th = np.linspace(0, 2*math.pi, 2001)
    gb = D.grad_norm(np.maximum(r0+d*np.cos(th),1e-6), z0+d*np.sin(th))
    return float(np.max(gn)), float(np.max(hf)), float(np.max(gb))

if __name__ == "__main__":
    import sys
    out = dict(RES)
    NS, NT, NC = 300, 120, 120
    # ---------------- (D-B) the campaign datum ------------------------------------------
    print("=== (D-B) campaign datum: omega = -M tanh(sin phi/sin d) tanh(cos phi/w) Theta(rho) ===")
    rowsB = []
    for L in (10.0, 40.0):
        DB = DatumB(delta_deg=7.5, w=0.20, L=L)
        for lam in (1.0, 1.25, 1.5):
            # the field at time t is eta_lam = eta_0 o T_lam^{-1} (inviscid transport of the
            # scalar eta); evaluating eta_0 at the moved point would be a different object.
            DL = StrainedB(DB, lam)
            g = geometry(30.0, 0.0, 7.5, lam)
            dmax = 0.95*min(g['d_eq'], g['d_taper'], 0.95*g['r'])
            best = None
            for d in np.linspace(0.05*dmax, dmax, 8):
                Ls = lambda_bulk(DL.grad_norm_fast, (g['r'],g['z']), d, 3*DB.R, NS,NT,NC)
                sg, sh, sb = sups_on_ball_B(DL, g['r'], g['z'], d)
                e = [v[0] for v in DL.eta_rz(np.array([g['r']]), np.array([g['z']]))]
                A = assemble(sg, sh, sb, d, Ls[4], Ls[5], g['r'], e[0], e[1], e[2])
                A.update(d=d, L4=Ls[4], L5=Ls[5], sup_grad=sg, sup_hess=sh, sup_grad_bdry=sb)
                if best is None or A['K2'] < best['K2']: best = A
            best.update(L=L, lam=lam, r=g['r'], z=g['z'], dmax=dmax)
            rowsB.append(best)
            print(f"  L={L:5.1f} lam={lam:4.2f}: d*={best['d']:.4f}  |grad a|<={best['grad_a']:.4f}"
                  f"  ||Hess a||<={best['hess_a']:.4f}  ||Hess u^z||<={best['hess_uz']:.4f}"
                  f"   K2hat <= {best['K2']:.4f}")
    out['DB'] = rowsB
    # ---------------- (D-A) idealised, sharp radial edge, inset f -----------------------
    print("\n=== (D-A) idealised datum (sharp sgn(z), corner taper, SHARP radial edges) ===")
    rowsA = []
    DA = DatumA(delta_deg=7.5, L=10.0)
    for f in (0.4, 0.2, 0.1, 0.05, 0.025):
        g = geometry(30.0, f, 7.5, 1.0)
        dmax = 0.95*min(g['d_eq'], g['d_taper'], f)          # ball must avoid the radial edge
        best = None
        surf4 = lambda_p_surface_DA(DA, (g['r'],g['z']), 4)
        surf5 = lambda_p_surface_DA(DA, (g['r'],g['z']), 5)
        for d in np.linspace(0.15*dmax, dmax, 6):
            Ls = lambda_bulk(DA.grad_norm, (g['r'],g['z']), d, 3*DA.R, NS,NT,NC)
            L4 = Ls[4] + sum(surf4.values()); L5 = Ls[5] + sum(surf5.values())
            Rm = g['r'] - d
            sg = 1.0/Rm**2; sh = math.sqrt(7)/Rm**3; sb = 1.0/Rm**2       # exact plateau values
            eta = -1.0/g['r']; eta_r = 1.0/g['r']**2; eta_z = 0.0
            A = assemble(sg, sh, sb, d, L4, L5, g['r'], eta, eta_r, eta_z)
            A.update(d=d, L4=L4, L5=L5, surf4=surf4, surf5=surf5, Rminus=Rm)
            if best is None or A['K2'] < best['K2']: best = A
        best.update(f=f, r=g['r'])
        rowsA.append(best)
        print(f"  f={f:6.3f}: d*={best['d']:.4f}  |grad a|<={best['grad_a']:.4f}"
              f"  ||Hess a||<={best['hess_a']:.4f}   K2hat <= {best['K2']:.4f}"
              f"   [edge surface part of L5 = {best['surf5']['inner']:.3g}]")
    out['DA'] = rowsA
    json.dump(out, open('k3_results.json','w'), indent=1)
    print("\nWROTE k3_results.json")

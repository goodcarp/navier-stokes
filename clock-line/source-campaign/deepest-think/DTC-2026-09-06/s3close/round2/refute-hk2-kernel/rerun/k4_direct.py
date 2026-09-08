#!/usr/bin/env python3
"""k4 -- an INDEPENDENT measurement of grad^2_5 b for the campaign datum (D-B), by the exact
zonal-harmonic solution of -Lap5 psi1 = eta.  Shares no code with k3's kernel quadrature.

  eta(rho,t) = G(rho) Hh(t),  Hh = sum_l hhat_l C_l^{3/2}(t)   (odd l only)
  psi1 = sum_l p_l(rho) C_l ,  p_l'' + (4/rho)p_l' - l(l+3)p_l/rho^2 = -e_l ,  e_l = -M hhat_l Theta/rho
  p_l = A_l (W1 + W2),  A_l = -M hhat_l/(2l+3),
      W1 = rho^{-l-3} int_0^rho s^{l+3}Theta ds ,  W2 = rho^l int_rho^inf s^{-l}Theta ds
  a = -d_z psi1 ,  d_z[p C_l] = (l+2)/(2l+3)[p'+(l+3)p/rho] C_{l-1} + (l+1)/(2l+3)[p'-l p/rho] C_{l+1}
Controls: (i) Lap5 a = d_z eta pointwise; (ii) a against a direct 5D kernel quadrature;
(iii) mode-truncation refinement.
"""
import json, math
import numpy as np
from scipy.special import eval_gegenbauer
from scipy.integrate import quad
from hk2lib import DatumB, geometry

def hhat(DB, LMAX, n=4001):
    """Gegenbauer coefficients of Hh(t) w.r.t. C_l^{3/2}, weight (1-t^2), N_l=(l+1)(l+2)/(l+1.5)."""
    x, w = np.polynomial.legendre.leggauss(n)
    H = DB.Hh(x, 0); ww = w*(1-x*x)
    out = np.zeros(LMAX+1)
    for l in range(LMAX+1):
        Nl = (l+1)*(l+2)/(l+1.5)
        out[l] = np.sum(ww*H*eval_gegenbauer(l, 1.5, x))/Nl
    return out

def W12(DB, l, rho, Rout):
    """W1 = rho int_0^1 v^{l+3}Theta(rho v)dv ,  W2 = rho int_1^inf v^{-l}Theta(rho v)dv,
    both in log variable so the l-dependence is a pure exponential."""
    f1 = lambda xi: math.exp((l+4)*xi)*DB.Th(rho*math.exp(xi), 0)
    lo = -min(60.0/(l+4), math.log(rho/1e-9))
    W1 = rho*quad(f1, lo, 0.0, limit=400, epsabs=1e-14, epsrel=1e-13)[0]
    f2 = lambda xi: math.exp((1-l)*xi)*DB.Th(rho*math.exp(xi), 0)
    hi = math.log(Rout/rho)
    if l >= 2: hi = min(hi, 80.0/(l-1))
    W2 = rho*quad(f2, 0.0, hi, limit=400, epsabs=1e-14, epsrel=1e-13)[0]
    return W1, W2

def p_derivs(DB, l, rho, Rout):
    """p_l and p_l', p_l'', p_l''' at rho."""
    W1, W2 = W12(DB, l, rho, Rout)
    A = -DB.M*0.0  # filled by caller (A_l depends on hhat)
    Th0 = DB.Th(rho,0); Th1 = DB.Th(rho,1)
    p0 = (W1 + W2)
    p1 = (-(l+3)*W1 + l*W2)/rho
    p2 = ((l+3)*(l+4)*W1 + l*(l-1)*W2)/rho**2 - (2*l+3)*Th0/rho
    p3 = ((l+3)*(l+4)*(-(l+5)*W1/rho + Th0) + l*(l-1)*((l-2)*W2/rho - Th0))/rho**2 \
         - (2*l+3)*Th1/rho + (2*l+3)*Th0/rho**2
    return np.array([p0, p1, p2, p3])

class Series:
    def __init__(self, DB, LMAX=161, Rout=None):
        self.DB = DB; self.LMAX = LMAX
        self.Rout = Rout if Rout is not None else 4.0*DB.R
        self.hh = hhat(DB, LMAX+2)
        self.cacheP = {}
    def P(self, l, rho):
        key = (l, round(rho, 14))
        if key not in self.cacheP:
            if l < 0 or abs(self.hh[l]) < 1e-18:
                self.cacheP[key] = np.zeros(4)
            else:
                A = -self.DB.M*self.hh[l]/(2*l+3)
                self.cacheP[key] = A*p_derivs(self.DB, l, rho, self.Rout)
        return self.cacheP[key]
    def alpha(self, m, rho):
        """alpha_m and its first two rho-derivatives, for a = -d_z psi1 = sum_m alpha_m C_m."""
        out = np.zeros(3)
        for src, sgn in ((m+1, 'down'), (m-1, 'up')):
            if src < 0: continue
            l = src; P = self.P(l, rho)
            if src == m+1:
                c = (l+2)/(2*l+3); k = (l+3)
            else:
                c = (l+1)/(2*l+3); k = -l
            # term = c*(p' + k p/rho) ; need term, term', term''
            p0,p1,p2,p3 = P
            t0 = c*(p1 + k*p0/rho)
            t1 = c*(p2 + k*(p1/rho - p0/rho**2))
            t2 = c*(p3 + k*(p2/rho - 2*p1/rho**2 + 2*p0/rho**3))
            out += np.array([t0, t1, t2])
        return -out
    def a_field(self, r, z, mmax=None):
        """returns a, a_r, a_z, a_rr, a_rz, a_zz at the single point (r,z)."""
        rho = math.hypot(r, z); t = z/rho; s = r/rho
        rr_, rz_ = s, t
        tr_, tz_ = -s*t/rho, (1-t*t)/rho
        rrr, rrz, rzz = t*t/rho, -s*t/rho, s*s/rho
        trr = t*(2*s*s - t*t)/rho**2; trz = s*(2*t*t - s*s)/rho**2; tzz = -3*s*s*t/rho**2
        A = np.zeros(6)
        M = self.LMAX if mmax is None else mmax
        for m in range(0, M+1):
            if abs(self.hh[m+1] if m+1 <= self.LMAX+2 else 0.0) < 1e-18 and \
               abs(self.hh[m-1] if 0 <= m-1 else 0.0) < 1e-18:
                continue
            a0, a1, a2 = self.alpha(m, rho)
            if a0 == 0.0 and a1 == 0.0 and a2 == 0.0: continue
            C  = eval_gegenbauer(m, 1.5, t)
            C1 = 3*eval_gegenbauer(m-1, 2.5, t) if m >= 1 else 0.0
            C2 = 15*eval_gegenbauer(m-2, 3.5, t) if m >= 2 else 0.0
            A[0] += a0*C
            A[1] += a1*rr_*C + a0*C1*tr_
            A[2] += a1*rz_*C + a0*C1*tz_
            A[3] += a2*rr_**2*C + a1*rrr*C + 2*a1*rr_*C1*tr_ + a0*C2*tr_**2 + a0*C1*trr
            A[4] += a2*rr_*rz_*C + a1*rrz*C + a1*(rr_*C1*tz_ + rz_*C1*tr_) + a0*C2*tr_*tz_ + a0*C1*trz
            A[5] += a2*rz_**2*C + a1*rzz*C + 2*a1*rz_*C1*tz_ + a0*C2*tz_**2 + a0*C1*tzz
        return A

def K2_exact(r, a_r, a_z, a_rr, a_rz, a_zz, eta, eta_r, eta_z, ne=241):
    """||grad^2_5 b||  = sup_{|e|=1} || d_e grad_5 b ||_op , built from k1's exact tensor."""
    Ha = np.zeros((5,5))
    Ha[0,0]=a_rr; Ha[0,4]=Ha[4,0]=a_rz; Ha[4,4]=a_zz
    Ha[1,1]=Ha[2,2]=Ha[3,3]=a_r/r
    ga = np.zeros(5); ga[0]=a_r; ga[4]=a_z
    q = a_z - eta; qr = a_rz - eta_r; qz = a_zz - eta_z
    Huz = np.zeros((5,5))
    Huz[0,0]=q + r*qr; Huz[1,1]=Huz[2,2]=Huz[3,3]=q
    Huz[0,4]=Huz[4,0]=r*qz; Huz[4,4]=-2*a_z - r*a_rz
    yv = np.zeros(5); yv[0]=r
    best = 0.0; arg = None
    # by axisymmetry e may be taken in span(e_r, e_perp, e_z)
    for th in np.linspace(0, math.pi, ne):
        for ph in np.linspace(0, 2*math.pi, ne//2):
            e = np.zeros(5)
            e[0]=math.sin(th)*math.cos(ph); e[1]=math.sin(th)*math.sin(ph); e[4]=math.cos(th)
            N = np.zeros((5,5))
            He = Ha@e
            for i in range(4):
                for j in range(5):
                    N[i,j] = He[j]*yv[i] + ga[j]*e[i] + (ga@e)*(1.0 if i==j else 0.0)
            N[4,:] = Huz@e
            v = np.linalg.norm(N, 2)
            if v > best: best, arg = v, (th, ph)
    return best, arg

if __name__ == "__main__":
    out = {}
    DB = DatumB(delta_deg=7.5, w=0.20, L=10.0)
    S = Series(DB, LMAX=161)
    out['hhat_decay'] = {str(l): float(S.hh[l]) for l in (1,3,5,11,31,61,101,141,161)}
    print("hhat_l :", {l: f"{S.hh[l]:.3e}" for l in (1,3,5,11,31,61,101,141,161)})

    # ---- control 1: Lap5 a = d_z eta ---------------------------------------------------
    print("\n=== control 1: Lap5 a - d_z eta  (the PDE the series must satisfy) ===")
    ctrl = []
    for (r,z) in [(0.5,0.866),(0.75,0.3849),(0.62,0.554),(0.4,1.1)]:
        A = S.a_field(r,z)
        lap = A[3] + 3*A[1]/r + A[5]
        e = [v[0] for v in DB.eta_rz(np.array([r]),np.array([z]))]
        dz_eta = e[2]
        ctrl.append(dict(r=r,z=z,lap_a=lap,dz_eta=dz_eta,rel=abs(lap-dz_eta)/max(abs(dz_eta),1e-30)))
        print(f"   (r,z)=({r},{z}): Lap5 a = {lap:.10f}   d_z eta = {dz_eta:.10f}   rel = "
              f"{abs(lap-dz_eta)/max(abs(dz_eta),1e-30):.2e}")
    out['control_pde'] = ctrl

    # ---- control 2: a against a direct 5D kernel quadrature ----------------------------
    print("\n=== control 2: a(x) series vs direct 5D kernel quadrature ===")
    def a_kernel(x_rz, ns=900, nT=200, nc=200):
        r,z = x_rz
        xg, wg = np.polynomial.legendre.leggauss(ns)
        smin, smax = 1e-6, 3*DB.R
        u = 0.5*(math.log(smax)-math.log(smin))*xg + 0.5*(math.log(smax)+math.log(smin))
        wu = 0.5*(math.log(smax)-math.log(smin))*wg
        s = np.exp(u)
        T, wT = np.polynomial.legendre.leggauss(nT); T = 0.5*math.pi*(T+1); wT = 0.5*math.pi*wT
        X, wX = np.polynomial.legendre.leggauss(nc); X = 0.5*math.pi*(X+1); wX = 0.5*math.pi*wX
        sT=np.sin(T); cT=np.cos(T); cX=np.cos(X)
        WA = (wT*sT**3)[:,None]*(wX*4*math.pi*np.sin(X)**2)[None,:]
        tot = 0.0
        # K(w) dw = 3 w_z/(8 pi^2 |w|^5) s^4 ds dOmega_4 = (3 cos(Theta)/(8 pi^2)) ds dOmega_4
        for si, wsi in zip(s, wu*s):   # ds = s du
            rp = np.sqrt(np.maximum(r*r + (si*sT[:,None])**2 - 2*r*si*sT[:,None]*cX[None,:],1e-300))
            zp = z - si*cT[:,None]*np.ones_like(cX)[None,:]
            et = DB.eta_rz(rp, zp)[0]
            tot += wsi*(3.0/(8*math.pi**2))*np.sum(WA*cT[:,None]*et)
        return tot
    kq = []
    for (r,z) in [(0.5,0.866),(0.75,0.3849)]:
        av = S.a_field(r,z)[0]; ak = a_kernel((r,z))
        kq.append(dict(r=r,z=z,a_series=av,a_kernel=ak,rel=abs(av-ak)/abs(av)))
        print(f"   (r,z)=({r},{z}): a_series = {av:.9f}   a_kernel = {ak:.9f}   rel = "
              f"{abs(av-ak)/abs(av):.2e}")
    out['control_kernel'] = kq

    # ---- control 3: mode truncation ----------------------------------------------------
    print("\n=== control 3: mode truncation ===")
    tr = []
    for mm in (41, 81, 121, 161):
        A = S.a_field(0.5, 0.866, mmax=mm)
        tr.append(dict(mmax=mm, a=A[0], a_rr=A[3], a_zz=A[5]))
        print(f"   mmax={mm:4d}: a={A[0]:.10f}  a_rr={A[3]:.10f}  a_zz={A[5]:.10f}")
    out['control_truncation'] = tr

    # ---- the measurement ---------------------------------------------------------------
    print("\n=== MEASURED ||grad^2_5 b|| at points of N_tau (L = 10, delta = 7.5 deg, w = 0.20) ===")
    meas = []
    for lam in (1.0, 1.25, 1.5):
        g = geometry(30.0, 0.0, 7.5, lam)
        r, z = g['r'], g['z']
        A = S.a_field(r, z)
        e = [v[0] for v in DB.eta_rz(np.array([r]),np.array([z]))]
        k2, arg = K2_exact(r, A[1],A[2],A[3],A[4],A[5], e[0], e[1], e[2])
        ha_ev = np.linalg.eigvalsh(np.array([[A[3],A[4]],[A[4],A[5]]]))
        hessa = max(abs(ha_ev).max(), abs(A[1]/r))
        meas.append(dict(lam=lam, r=r, z=z, a=A[0], a_r=A[1], a_z=A[2], a_rr=A[3], a_rz=A[4],
                         a_zz=A[5], eta=e[0], eta_r=e[1], eta_z=e[2],
                         grad_a=math.hypot(A[1],A[2]), hess_a=hessa, K2=k2))
        print(f"   lam={lam:4.2f} (r,z)=({r:.5f},{z:.5f}):  a={A[0]:.5f}  |grad a|="
              f"{math.hypot(A[1],A[2]):.5f}  ||Hess a||={hessa:.5f}   ||grad^2 b|| = {k2:.5f}")
    out['measured'] = meas
    json.dump(out, open('k4_results.json','w'), indent=1)
    print("\nWROTE k4_results.json")

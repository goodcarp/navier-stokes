#!/usr/bin/env python3
"""s6_deformed_field.py -- SELF-ATTACK: does Lemma 3 survive on the DEFORMED field?

Lemma 3 (velocity = per-shell pure strain + O(M rho)) was established for the t = 0 plateau.
The argument needs it again at every t, on eta_lam = eta_0 o T_lam^{-1}.

T_lam commutes with isotropic dilation, so eta_lam is again homogeneous of degree -1 with
angular profile  g_lam(phi) = rho eta_lam(rho,phi).  Everything therefore reduces to the
Gegenbauer coefficients H_l(lam) of g_lam:
  * H_1(lam) fixes the log coefficient of a:  kappa(lam) = -3 H_1(lam)/5.  Lemma 1 says this
    must be exactly lam/2 -- a THIRD independent route to P_1(lam) = lam.
  * sum_{l>=3} |H_l(lam)| * ||C_l||_inf / (l(l+3)-4) bounds the Lemma-3 remainder; it must stay
    bounded uniformly for lam in [1, 3/2].
"""
import numpy as np, json
OUT = {}

def geg_all(alpha, L, t):
    C = np.zeros((L+1,)+np.shape(t)); C[0] = 1.0
    if L >= 1: C[1] = 2*alpha*t
    for n in range(2, L+1):
        C[n] = (2*t*(n+alpha-1)*C[n-1] - (n+2*alpha-2)*C[n-2])/n
    return C

def h_taper(delta):
    def h(phi):
        pax = np.minimum(phi, np.pi-phi); return np.minimum(1.0, pax/delta)
    return h
h_flat = lambda phi: np.ones_like(phi)

def g_lam(phi, lam, h):
    """rho * eta_lam at Eulerian polar angle phi (M = 1)."""
    s, c = np.sin(phi), np.cos(phi)
    rp, zp = s/lam, c*lam**2                     # preimage direction (rho = 1)
    rhop = np.hypot(rp, zp); phip = np.arctan2(rp, np.abs(zp))*np.sign(zp)
    phip = np.where(phip < 0, np.pi + phip, phip)
    return -np.sign(c)*h(phip)/(rhop*np.sin(phip))   # omega^theta = -M sgn(z) h

LMAX = 601
xg, wg = np.polynomial.legendre.leggauss(8000)
tq = xg.copy(); wq = wg.copy()                     # t in [-1,1]
phiq = np.arccos(np.clip(tq, -1, 1))
Cq = geg_all(1.5, LMAX, tq)
l = np.arange(LMAX+1); Nl = (l+1)*(l+2)/(l+1.5)
Cinf = (l+1)*(l+2)/2.0                             # ||C_l^{3/2}||_inf = C_l(1)
den = np.maximum(l*(l+3)-4.0, 1.0)

def coeffs(lam, h):
    g = g_lam(phiq, lam, h)
    w = wq*(1-tq**2)                               # Gegenbauer weight (1-t^2)
    Hl = np.einsum('lq,q->l', Cq, w*g)/Nl
    Hl[l % 2 == 0] = 0.0
    return Hl

res = []
for lam in (1.0, 1.1, 1.25, 1.5):
    for name, h in (('flat', h_flat), ('taper7.5', h_taper(np.deg2rad(7.5)))):
        Hl = coeffs(lam, h)
        kap = -3*Hl[1]/5
        tail_worst = float(np.sum(np.abs(Hl[3:])*Cinf[3:]/den[3:]))
        tail_l1 = float(np.sum(np.abs(Hl[3:])))
        res.append(dict(lam=lam, profile=name, kappa=float(kap), kappa_target=lam/2*(1 if name=='flat' else None) if name=='flat' else None,
                        rel_err_vs_lam_over_2=(float(abs(kap-lam/2)/(lam/2)) if name=='flat' else None),
                        sum_absH_tail=tail_l1, worstcase_Lemma3_bound=tail_worst,
                        H_l_decay=[float(Hl[k]) for k in (101,301,601)]))
OUT['deformed_modes'] = res
print(json.dumps(OUT, indent=1, default=str))
json.dump(OUT, open('s6_results.json','w'), indent=1, default=str)

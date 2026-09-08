#!/usr/bin/env python3
"""s5_velocity_structure.py -- LEMMA 3: the velocity of the plateau IS a per-shell pure strain
plus O(M rho), with the O(M rho) constant computed.

5D lift: Delta_5 psi = -eta, Delta_5 = d_rr + (3/r)d_r + d_zz on R^5 = R^4_r x R_z; zonal
harmonics on S^4 are the Gegenbauer C_l^{3/2}(t), t = cos phi, orthogonal with weight (1-t^2),
norm N_l = (l+1)(l+2)/(l+3/2).  Radial operator: Delta_5(rho^g C_l) = [g(g+3)-l(l+3)] rho^{g-2} C_l.

Bang-bang plateau  eta = -M sgn(t)/(rho sqrt(1-t^2))  on rho_0<rho<R  =>  eta = (1/rho) sum_l H_l C_l,
    H_l = -2M int_0^1 sqrt(1-t^2) C_l(t) dt / N_l   (l odd; H_l = 0 for l even).
Source rho^{-1} resonates ONLY with l = 1 (g=1: g(g+3)-l(l+3) = 4-4 = 0), and
    Delta_5(rho log rho C_1) = 5 rho^{-1} C_1,   so   psi_1 = -(H_1/5) rho log rho C_1
                                                        = (M/2) z log rho   (H_1 = -5M/6, C_1 = 3t).
Adding the homogeneous constant:  psi_1 = -A(rho) z,  A(rho) = (M/2) log(R/rho)  --  EXACTLY the
uniform axisymmetric strain (u^r, u^z) = A(rho)(r, -2z) up to O(M rho).
Every l >= 3 mode is NON-resonant: in the bulk f_l = H_l rho /(l(l+3)-4), i.e. psi_dev = M rho S(t),
which contributes O(M) to a and O(M rho) to u -- a factor 1/log(R/rho) below the strain.

This script (i) computes H_l and checks H_1=-5/6, H_3=3/40, H_5=-247/1680 (estate's c6 values),
(ii) sums S(t) and reports sup|a_dev|/M and sup|u^z_dev|/(M rho) -- the explicit Lemma-3 constants,
(iii) does the FULL four-region C^1 matching at rho_0 and R for every l and reproduces, by a route
      completely independent of the campaign's elliptic-integral Biot-Savart, the two numbers
      a(0,0) = (M/2)log(R/rho_0)  and  a(inner material shell) = (M/2)log(R/rho_0) + 0.216773 M.
"""
import numpy as np, json
OUT = {}
M = 1.0

def geg_all(alpha, L, t):
    """C_0..C_L^{alpha}(t) by the standard recurrence, t an array."""
    C = np.zeros((L+1,)+np.shape(t))
    C[0] = 1.0
    if L >= 1: C[1] = 2*alpha*t
    for n in range(2, L+1):
        C[n] = (2*t*(n+alpha-1)*C[n-1] - (n+2*alpha-2)*C[n-2])/n
    return C

LMAX = 801
# ---- H_l by Gauss-Legendre on [0,1] of sqrt(1-t^2) C_l(t)
xg, wg = np.polynomial.legendre.leggauss(6000)
tq = 0.5*(xg+1.0); wq = 0.5*wg
Cq = geg_all(1.5, LMAX, tq)
num = -2*M*np.einsum('lq,q->l', Cq, wq*np.sqrt(1-tq**2))
l = np.arange(LMAX+1)
Nl = (l+1)*(l+2)/(l+1.5)
H = num/Nl
H[l % 2 == 0] = 0.0
OUT['H1'] = float(H[1]); OUT['H1_exact'] = -5/6
OUT['H3'] = float(H[3]); OUT['H3_exact'] = 3/40
OUT['H5'] = float(H[5]); OUT['H5_exact'] = -247/1680
OUT['H_decay_l_51_101_201_401_801'] = [float(H[k]) for k in (51,101,201,401,801)]
OUT['sum_absH_over_l2_tail_from_3'] = float(np.sum(np.abs(H[3:])/np.maximum(l[3:]*(l[3:]+3)-4,1)))

# ---- (ii) bulk deviation  psi_dev = M rho S(t)
tt = np.linspace(-0.999, 0.999, 4001)
C32 = geg_all(1.5, LMAX, tt)
C52 = geg_all(2.5, LMAX, tt)          # dC_l^{3/2}/dt = 3 C_{l-1}^{5/2}
den = l*(l+3)-4.0
odd3 = [k for k in range(3, LMAX+1, 2)]
S  = np.sum([H[k]*C32[k]/den[k] for k in odd3], axis=0)
Sp = np.sum([H[k]*3*C52[k-1]/den[k] for k in odd3], axis=0)
a_dev  = -(tt*S + (1-tt**2)*Sp)
uz_dev = 2*S + (1-tt**2)*(S - tt*Sp)
# l=1 non-strain remainder: a gets -(M/2)t^2 ; u^z gets +(M/2) rho t (1-t^2)
a_dev_tot  = a_dev - 0.5*tt**2
uz_dev_tot = uz_dev + 0.5*tt*(1-tt**2)
# restrict to the region the tracked point lives in (|t| <= cos(20 deg), away from the equator jump)
msk = np.abs(tt) <= np.cos(np.deg2rad(20.0))
OUT['lemma3_constants'] = {
    'sup_a_dev_over_M_allt':      float(np.max(np.abs(a_dev_tot))),
    'sup_a_dev_over_M_bulkband':  float(np.max(np.abs(a_dev_tot[msk]))),
    'sup_uz_dev_over_M_rho_allt': float(np.max(np.abs(uz_dev_tot))),
    'sup_uz_dev_over_M_rho_band': float(np.max(np.abs(uz_dev_tot[msk]))),
}
# partial-sum convergence of S at three sample angles
samples = [np.cos(np.deg2rad(x)) for x in (30.0, 54.7356, 75.0)]
conv = {}
for s0 in samples:
    C32s = geg_all(1.5, LMAX, np.array([s0]))[:,0]
    vals = []
    tot = 0.0
    for k in odd3:
        tot += H[k]*C32s[k]/den[k]
        if k in (51,101,201,401,801): vals.append((k, float(tot)))
    conv[f't={s0:.4f}'] = vals
OUT['S_partial_sums'] = conv

# ---- (iii) full C^1 matching, all four regions (scaled bases: no overflow)
def a_field(rho, t, r0, r1, Lmax=LMAX):
    """a = -d_z psi, psi = sum_l f_l(rho) C_l^{3/2}(t), exact C^1 matching at r0 and r1.
    Bases are scaled so every one is <= 1 on its own region:
      I  : f = A (rho/r0)^l                (rho < r0)
      II : f = p(rho) + B (rho/r1)^l + C (r0/rho)^{l+3}
      III: f = D (r1/rho)^{l+3}            (rho > r1)"""
    C32v = geg_all(1.5, Lmax, np.array([t]))[:,0]
    C52v = geg_all(2.5, Lmax, np.array([t]))[:,0]
    eps = r0/r1
    a = 0.0
    for k in range(1, Lmax+1, 2):
        Hk = H[k]
        if Hk == 0.0: continue
        if k == 1:
            p  = lambda r: -(Hk/5.0)*r*np.log(r)
            pp = lambda r: -(Hk/5.0)*(np.log(r)+1.0)
        else:
            cc = Hk/(k*(k+3)-4.0)
            p  = lambda r, cc=cc: cc*r
            pp = lambda r, cc=cc: cc
        ek  = eps**k if k*np.log(eps) > -700 else 0.0
        ek3 = eps**(k+3) if (k+3)*np.log(eps) > -700 else 0.0
        m = np.array([
            [1.0,      -ek,             -1.0,                 0.0],
            [k/r0,     -k*ek/(eps*r1),   (k+3)/r0,            0.0],
            [0.0,       1.0,             ek3,                -1.0],
            [0.0,       k/r1,           -(k+3)*ek3/r1,       (k+3)/r1],
        ])
        rhs = np.array([p(r0), pp(r0), -p(r1), -pp(r1)])
        A_,B_,C_,D_ = np.linalg.solve(m, rhs)
        if   rho < r0: f, fp = A_*(rho/r0)**k, A_*k/r0*(rho/r0)**(k-1)
        elif rho <= r1:
            wk = (rho/r1)**k; vk = (r0/rho)**(k+3)
            f  = p(rho) + B_*wk + C_*vk
            fp = pp(rho) + B_*k/rho*wk - C_*(k+3)/rho*vk
        else:
            qk = (r1/rho)**(k+3); f, fp = D_*qk, -D_*(k+3)/rho*qk
        a += -(t*fp*C32v[k] + (1-t*t)/rho*f*3*C52v[k-1])
    return float(a)

mm = {}
for Rr in (64.0, 256.0, 1024.0, 4096.0):
    a_ax  = a_field(1e-4*Rr**0, 0.999999, 1.0, Rr)        # deep inside the hole ~ the origin
    a_mat = a_field(1.0+1e-12, np.cos(np.deg2rad(10.0)), 1.0, Rr)
    mm['R=%d'%Rr] = dict(a_origin=a_ax, half_logR=0.5*np.log(Rr),
                         a_origin_minus_halflogR=a_ax-0.5*np.log(Rr),
                         a_material_inner_edge_phi10=a_mat,
                         offset_vs_halflogR=a_mat-0.5*np.log(Rr))
OUT['matched_field'] = mm

# convergence in Lmax of the two headline numbers at R = 1024
cv = {}
for Lm in (101, 201, 401, 801):
    cv['Lmax=%d'%Lm] = dict(a_origin=a_field(1e-4, 0.999999, 1.0, 1024.0, Lm),
                            a_material=a_field(1.0+1e-12, np.cos(np.deg2rad(10.0)), 1.0, 1024.0, Lm))
OUT['matched_field_Lmax_convergence_R1024'] = cv

print(json.dumps(OUT, indent=1, default=str))
json.dump(OUT, open('s5_results.json','w'), indent=1, default=str)

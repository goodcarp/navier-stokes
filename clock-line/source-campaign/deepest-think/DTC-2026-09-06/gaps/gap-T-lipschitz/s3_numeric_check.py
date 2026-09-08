"""s3 - numerical verification of LEMMA T on two map families x two data x three (lambda,mu).

Reports, for each case: the measured sup relative displacement mu, the measured sup relative
Jacobian deviation muJ, the exact-quadrature a for T_lambda and for Phi, |Delta|, the LEMMA T
bound evaluated at the MEASURED (mu,muJ), the slack, and the doubling-convergence residual.
KILL K3 fires if |Delta| > bound.  KILL K4 fires if doubling moves Delta by >1e-6 relative.
KILL K5 fires if a[T_lambda] != (M/2) lambda L to 1e-10 for the bang-bang cap.
"""
import json, numpy as np
from lib_quad import *

M, rho0, R = 1.0, 1.0, 4096.0
L = np.log(R/rho0)
delta = np.deg2rad(7.5)
CASES = [(1.00, 0.05), (1.25, 0.10), (1.50, 0.20)]

def panels_for(data, extra=()):
    p = [0.0, np.pi/2, np.pi]
    if data == 'D2_taper': p += [delta, np.pi-delta]
    p += list(extra)
    return np.array(sorted(set(p)))

def dense_sup(mp, lam, n=2001):
    """dense scan of the two hypothesis constants over the shell (log-radius x angle)"""
    ph = np.linspace(1e-9, np.pi-1e-9, n)
    tau = np.linspace(0, 1, 257)
    PH, TAU = np.meshgrid(ph, tau, indexing='ij')
    RHO = rho0*np.exp(L*TAU)
    r, z = RHO*np.sin(PH), RHO*np.cos(PH)
    _,_,J5,rel = mp(r, z)
    return float(np.max(rel)), float(np.max(np.abs(J5/lam**2-1)))

rows = []
for data, omega in (('D1_bangbang', lambda p: omega_bangbang(p, M)),
                    ('D2_taper',    lambda p: omega_taper(p, delta, M))):
    for lam, mu in CASES:
        mub = 2*np.arcsin(mu/2)                      # so that sup 2 sin(beta/2) = mu
        maps = [MapA(lam, mu, rho0, L, k=1),
                MapAng(lam,
                       lambda p, mub=mub: mub*np.sin(p)**2,
                       lambda p, mub=mub: mub*np.sin(2*p),
                       'B_angular_shear')]
        pT = panels_for(data)
        for mp in maps:
            homog = isinstance(mp, MapAng)
            res = {}
            for n_phi, n_tau in ((60, 60), (120, 120)):
                aT,_,_ = a_of_map(MapT(lam), omega, lam, rho0, R, pT, n_phi, n_tau, True)
                aP, mu_n, muJ_n = a_of_map(mp, omega, lam, rho0, R, pT, n_phi, n_tau, homog)
                res[n_phi] = (aT, aP)
            aT, aP = res[120]
            conv = abs((res[120][1]-res[60][1])/res[120][1])
            mu_s, muJ_s = dense_sup(mp, lam)
            D = abs(aP-aT); B = bound(lam, M, L, mu_s, muJ_s)
            rows.append(dict(data=data, map=mp.name, lam=lam, mu_target=mu,
                             mu_sup=mu_s, muJ_sup=muJ_s, a_T=aT, a_Phi=aP,
                             Delta=aP-aT, absDelta=D, bound=B, slack=B/D if D>0 else float('inf'),
                             K3_fired=bool(D > B), conv_rel=conv, K4_fired=bool(conv > 1e-6)))

# KILL K5: Lemma 1 identity for the bang-bang cap
k5 = []
for lam in (1.0, 1.25, 1.5, 2.0, 0.7):
    aT,_,_ = a_of_map(MapT(lam), lambda p: omega_bangbang(p, M), lam, rho0, R,
                      panels_for('D1_bangbang'), 200, 2, True)
    k5.append(dict(lam=lam, a=aT, exact=M/2*lam*L, resid=abs(aT-M/2*lam*L)))
K5 = any(k['resid'] > 1e-10 for k in k5)

# KILL K1: independent recomputation of the shell integral  I = int_S |eta0| |T_lam x|^-4 dx_5
def shell_I(lam, n=400):
    ph, wph = gl_nodes(np.array([0.0, np.pi/2, np.pi]), n)
    # integrand in (tau,phi): 2 pi^2 rho^4 sin^3 phi * (M/(rho sin phi)) * (rho g5)^-4 * L rho dtau
    g5 = np.sqrt(lam**2*np.sin(ph)**2 + lam**-4*np.cos(ph)**2)
    return 2*np.pi**2*M*L*np.sum(wph*np.sin(ph)**2/g5**4)
K1rows = []
for l in (0.5,1.0,1.25,1.5,3.0):
    q = shell_I(l); c = np.pi**3*M*L/l
    K1rows.append(dict(lam=l, I_quad=q, I_closed=c, rel=abs(q-c)/c))
K1 = any(k['rel'] > 1e-8 for k in K1rows)

out = dict(L=L, M=M, rho0=rho0, R=R, delta_deg=7.5, rows=rows, K5_rows=k5, K5_fired=bool(K5),
           K1_rows=K1rows, K1_fired=bool(K1),
           K3_fired_any=bool(any(r['K3_fired'] for r in rows)),
           K4_fired_any=bool(any(r['K4_fired'] for r in rows)))
json.dump(out, open('s3_results.json','w'), indent=1)

print(f"L = {L:.9f}")
print("\nK1 shell integral  I = pi^3 M L / lambda :")
for k in K1rows: print(f"  lam={k['lam']:<5} quad={k['I_quad']:.10f} closed={k['I_closed']:.10f} rel={k['rel']:.2e}")
print("\nK5 Lemma 1  a[T_lam] = (M/2) lam L :")
for k in k5: print(f"  lam={k['lam']:<5} a={k['a']:.12f} exact={k['exact']:.12f} resid={k['resid']:.2e}")
print("\nLEMMA T check:")
hdr = f"{'data':12s} {'map':18s} {'lam':>5s} {'mu':>8s} {'muJ':>8s} {'a_T':>11s} {'a_Phi':>11s} {'|Delta|':>10s} {'bound':>11s} {'slack':>7s} {'conv':>9s}"
print(hdr)
for r in rows:
    print(f"{r['data']:12s} {r['map']:18s} {r['lam']:5.2f} {r['mu_sup']:8.5f} {r['muJ_sup']:8.5f} "
          f"{r['a_T']:11.6f} {r['a_Phi']:11.6f} {r['absDelta']:10.6f} {r['bound']:11.4f} {r['slack']:7.2f} {r['conv_rel']:9.2e}")
print(f"\nK1 {out['K1_fired']}  K3 {out['K3_fired_any']}  K4 {out['K4_fired_any']}  K5 {out['K5_fired']}")

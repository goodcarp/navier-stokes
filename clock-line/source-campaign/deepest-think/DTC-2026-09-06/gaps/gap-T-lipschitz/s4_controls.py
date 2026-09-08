"""s4 - the registered CONTROLS.

C1 (K6): angular shear with an axis ramp of width phi_c,
    beta(phi) = mu_beta * m(phi/phi_c) * m((pi-phi)/phi_c) * cos(phi),  m(t)=t^2(3-2t) on [0,1], 1 after.
    This is a genuine C^1 SO(4)-equivariant diffeomorphism (1+beta' >= 1-mu_beta > 0), the image-relative
    displacement is sup 2|sin(beta/2)| = 2 sin(mu_beta/2) = mu, INDEPENDENT of phi_c -- but the relative
    Jacobian deviation muJ ~ (sin mu_beta/phi_c)^3 blows up as phi_c -> 0.
    Prediction registered in PREREG: |Delta| ~ 2*(3/4) cos(mu_beta) sin^3(mu_beta) * L * log(1/phi_c),
    hence |Delta|/(mu lambda M L) -> infinity: the LINEAR corollary |Delta| <= C0 mu lambda M L must FAIL.
C2: mu >= 1 (hypothesis (H2) void).
Both are run on the inadmissible bang-bang cap D1 and on the admissible 7.5-degree taper D2.
"""
import json, numpy as np
from lib_quad import *

M, rho0, R = 1.0, 1.0, 4096.0
L = np.log(R/rho0); delta = np.deg2rad(7.5)
C0 = 15*np.pi/8

def m(t):  return np.where(t < 1, t**2*(3-2*t), 1.0)
def dm(t): return np.where(t < 1, 6*t*(1-t), 0.0)

def make_beta(mub, pc):
    def b(p):  return mub*m(p/pc)*m((np.pi-p)/pc)*np.cos(p)
    def db(p):
        A, B = m(p/pc), m((np.pi-p)/pc)
        dA, dB = dm(p/pc)/pc, -dm((np.pi-p)/pc)/pc
        return mub*((dA*B + A*dB)*np.cos(p) - A*B*np.sin(p))
    return b, db

def panels(pc, data):
    lo = max(pc*1e-3, 1e-30)
    left  = np.concatenate(([0.0], np.geomspace(lo, np.pi/2, int(6*np.log10(np.pi/2/lo))+12)))
    p = list(left) + list(np.pi - left[left < np.pi/2][::-1]) + [np.pi]
    if data == 'D2': p += [delta, np.pi-delta]
    return np.array(sorted(set([x for x in p if 0 <= x <= np.pi])))

def run(mub, pc, data, lam=1.0, n=40):
    om = (lambda q: omega_bangbang(q, M)) if data == 'D1' else (lambda q: omega_taper(q, delta, M))
    b, db = make_beta(mub, pc)
    mp = MapAng(lam, b, db, f'C1_axis_ramp_pc={pc:g}')
    P = panels(pc, data)
    aT,_,_ = a_of_map(MapT(lam), om, lam, rho0, R, P, n, 2, True)
    aP, mu_n, muJ_n = a_of_map(mp, om, lam, rho0, R, P, n, 2, True)
    aP2,_,_ = a_of_map(mp, om, lam, rho0, R, P, 2*n, 2, True)
    mu = 2*np.sin(mub/2)
    return dict(mu_beta=mub, phi_c=pc, data=data, lam=lam, mu=mu, muJ_nodes=muJ_n,
                a_T=aT, a_Phi=aP, absDelta=abs(aP-aT),
                linear_rhs=C0*mu*lam*M*L, ratio=abs(aP-aT)/(mu*lam*M*L),
                full_bound=bound(lam, M, L, mu, muJ_n),
                conv=abs((aP2-aP)/aP) if aP != 0 else 0.0)

out = {'C0': C0, 'L': L, 'C1': [], 'C1_taper': [], 'C2': []}
for mub in (0.9, 0.5, 0.3):
    for pc in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-12, 1e-16, 1e-20, 1e-24):
        out['C1'].append(run(mub, pc, 'D1'))
for pc in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
    out['C1_taper'].append(run(0.9, pc, 'D2'))

# C2 : the radial direction, up to and past the diffeomorphism limit of map A (k=1).
#      NOTE (erratum vs PREREG): PREREG's C2 used the k=2 (full-log-period) ripple; s5 shows that
#      family is an EXACT NULL DIRECTION of a(0), so it can never fire.  C2 is run here on the
#      non-null k=1 ripple.  Map A is a diffeomorphism for mu * sqrt(1+(pi/L)^2) < 1, i.e. mu < 0.9342.
for mu in (0.05, 0.2, 0.5, 0.9, 0.93, 1.2, 1.5):
    mp = MapA(1.0, mu, rho0, L, k=1)
    om = lambda q: omega_bangbang(q, M)
    P = np.array([0.0, np.pi/2, np.pi])
    aT,_,_ = a_of_map(MapT(1.0), om, 1.0, rho0, R, P, 120, 120, True)
    aP, mu_n, muJ_n = a_of_map(mp, om, 1.0, rho0, R, P, 160, 240, False)
    out['C2'].append(dict(mu=mu, a_T=aT, a_Phi=aP, absDelta=abs(aP-aT),
                          ratio=abs(aP-aT)/(mu*M*L), muJ_nodes=muJ_n,
                          is_diffeo=bool(mu*np.sqrt(1+(np.pi/L)**2) < 1),
                          note='(1-mu)^-5 undefined for mu>=1; map not injective'))

# empirical growth law: |Delta| = c(mu_beta) * phi_c^{-p}
fits = {}
for mub in (0.9, 0.5, 0.3):
    rs = [r for r in out['C1'] if r['mu_beta'] == mub and 1e-12 <= r['phi_c'] <= 1e-2]
    x = np.log10([r['phi_c'] for r in rs]); y = np.log10([r['absDelta'] for r in rs])
    p, c = np.polyfit(x, y, 1)
    fits[str(mub)] = dict(exponent=float(-p), coeff=float(10**c), n=len(rs))
out['growth_fits'] = fits
mbs = np.array([0.9,0.5,0.3]); cs = np.array([fits[str(m)]['coeff'] for m in mbs])
out['mu_beta_exponent'] = float(np.polyfit(np.log(mbs), np.log(cs), 1)[0])
fired = any((r['ratio'] > C0) for r in out['C1'])
out['K6_fired'] = bool(fired)
json.dump(out, open('s4_results.json','w'), indent=1)

print(f"C0 = {C0:.6f}   L = {L:.6f}")
print("\nCONTROL C1 (bang-bang cap D1): ratio = |Delta|/(mu lambda M L)  -- must exceed C0")
print(f"{'mu_beta':>8s} {'phi_c':>8s} {'mu':>7s} {'|Delta|':>10s} {'C0*mu*L':>9s} {'ratio':>9s} {'muJ':>11s} {'fullbound':>12s} {'conv':>9s}")
for r in out['C1']:
    flag = ' <== VIOLATES linear bound' if r['ratio'] > C0 else ''
    print(f"{r['mu_beta']:8.2f} {r['phi_c']:8.0e} {r['mu']:7.4f} {r['absDelta']:10.4f} {r['linear_rhs']:9.3f} "
          f"{r['ratio']:9.4f} {r['muJ_nodes']:11.3e} {r['full_bound']:12.3e} {r['conv']:9.1e}{flag}")
print("\nCONTROL C1 on the ADMISSIBLE 7.5-degree taper (D2):")
for r in out['C1_taper']:
    print(f"  phi_c={r['phi_c']:8.0e} |Delta|={r['absDelta']:10.5f} ratio={r['ratio']:9.4f} conv={r['conv']:.1e}")
print("\nCONTROL C2 (radial ripple k=1, up to and past the diffeo limit mu=0.9342):")
for r in out['C2']:
    print(f"  mu={r['mu']:4.2f} diffeo={str(r['is_diffeo']):5s} |Delta|={r['absDelta']:10.5f} "
          f"ratio={r['ratio']:8.4f} muJ={r['muJ_nodes']:.3f}")
print("\nempirical growth law |Delta| = c * phi_c^-p  (fit over phi_c in [1e-12,1e-2]):")
for k,v in out['growth_fits'].items(): print(f"  mu_beta={k}: p={v['exponent']:.4f}  c={v['coeff']:.5f}")
print(f"  c(mu_beta) ~ mu_beta^{out['mu_beta_exponent']:.3f}")
print(f"\nK6 fired: {out['K6_fired']}")

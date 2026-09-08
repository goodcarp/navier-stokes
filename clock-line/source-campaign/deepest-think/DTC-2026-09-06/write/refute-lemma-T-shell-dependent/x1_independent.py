"""x1 - REFUTER's independent re-derivation of every load-bearing number of Lemma T'.
Routes deliberately DIFFERENT from the target seat's:
  (A) det D Lambda by finite-difference of the actual 5x5 map (no sympy, no hand matrix).
  (B) the two shell identities by mpmath quadrature at random lambda.
  (C) a[eta_0 o Lambda^{-1}](0) computed ENTIRELY IN THE w-VARIABLES via numerical inversion
      of Lambda (fixed point), which never touches det D Lambda at all.
  (D) the exact constants of the integro-ODE profile from scratch.
"""
import json, math, random
import numpy as np
import mpmath as mp

out = {}
mp.mp.dps = 30

# ---------------- profile ----------------
KT = 2*(1 - math.sqrt(2/3))              # kappa*theta terminal
lam_s   = lambda s: (1 - (1-s)*KT/2)**-2                 # lambda(sigma)
dlogl_s = lambda s: -KT/(1 - (1-s)*KT/2)                 # dlog lambda/dsigma
out['KT'] = KT
out['lam0'] = lam_s(0.0); out['lam1'] = lam_s(1.0)

# lam_bar, int 1/lam, kappa_s : own mpmath quadrature, compared to closed forms
LB   = float(mp.quad(lambda s: (1-(1-s)*mp.mpf(KT)/2)**-2, [0,1]))
LINV = float(mp.quad(lambda s: (1-(1-s)*mp.mpf(KT)/2)**2,  [0,1]))
KS   = max(abs(dlogl_s(s)) for s in np.linspace(0,1,200001))
out['lam_bar_quad'] = LB;   out['lam_bar_closed'] = math.sqrt(6)/2
out['lam_inv_quad'] = LINV; out['lam_inv_closed'] = 5/9 + math.sqrt(6)/9
out['kappa_s_grid'] = KS;   out['kappa_s_closed'] = math.sqrt(6) - 2
out['C_rel_closed'] = 1.5*math.pi*(math.sqrt(6)-2)
out['muJL_closed']  = 2*math.sqrt(6) - 4
out['C0_abs'] = 15*math.pi/8; out['C0_rel'] = 15*math.pi/4
out['C_TLam_limit'] = 15*math.sqrt(6)*math.pi/16
out['muJL_over_c']  = (2*math.sqrt(6)-4)/(2*math.log(1.5))
out['CK'] = 3/(8*math.pi**2); out['CgradK'] = 3/(2*math.pi**2)

# ---------------- (A) det D Lambda by finite differences of the real 5-D map ----------
def Lam(x, lamf):
    rho = np.linalg.norm(x); l = lamf(rho)
    return np.array([l*x[0], l*x[1], l*x[2], l*x[3], x[4]/l**2])
def detFD(x, lamf, h=1e-6):
    Jm = np.zeros((5,5))
    for j in range(5):
        e = np.zeros(5); e[j] = h
        Jm[:,j] = (Lam(x+e,lamf) - Lam(x-e,lamf))/(2*h)
    return np.linalg.det(Jm)
rng = random.Random(20260906)
errs = []
for _ in range(400):
    m = rng.uniform(-0.3, 0.3); a = rng.uniform(0.7, 1.4)
    lamf = lambda r, a=a, m=m: a*r**m
    x = np.array([rng.gauss(0,1) for _ in range(5)])
    x *= rng.uniform(0.6, 3.0)/np.linalg.norm(x)
    rho = np.linalg.norm(x); l = lamf(rho); D = m           # rho lam'/lam = m
    ct = x[4]/rho; st2 = 1-ct**2
    closed = l**2*(1 + D*(st2 - 2*ct**2))
    errs.append(abs(detFD(x,lamf) - closed)/abs(closed))
out['detDLambda_FD_max_rel'] = max(errs)

# ---------------- (B) the identities, own quadrature ----------------
def g2(p, l): return l**2*mp.sin(p)**2 + l**-4*mp.cos(p)**2
def J_id(l):
    l = mp.mpf(l)
    return mp.quad(lambda p: mp.sin(p)**2/g2(p,l)**2, [0, mp.pi/2, mp.pi])
def I2_id(l):
    l = mp.mpf(l)
    return mp.quad(lambda p: mp.cos(p)*mp.sin(p)**2/g2(p,l)**mp.mpf(2.5), [0, mp.pi/2])
lams = [0.37, 1.0, 1.2247448713915890, 1.5, 2.9]
out['J_identity_max_rel'] = max(float(abs(J_id(l) - mp.pi/(2*mp.mpf(l)))/(mp.pi/(2*mp.mpf(l)))) for l in lams)
out['I2_identity_max_rel'] = max(float(abs(I2_id(l) - mp.mpf(l)/3)/(mp.mpf(l)/3)) for l in lams)

# Q(lambda) by quadrature and the seat's closed form
def Q_quad(l):
    l = mp.mpf(l)
    return mp.quad(lambda p: mp.cos(p)*mp.sin(p)**2*(mp.sin(p)**2-2*mp.cos(p)**2)/g2(p,l)**mp.mpf(2.5),
                   [0, mp.pi/2])
def Q_closed(l):
    l = mp.mpf(l)
    if abs(l-1) < mp.mpf('1e-20'): return -mp.mpf(1)/15
    s = mp.sqrt(l**6-1)
    return l*(-2*l**12*s + 9*l**9*mp.asinh(s) - 8*l**6*s + s)/(3*s*(l**12-2*l**6+1))
out['Q_closed_vs_quad_max_rel'] = max(float(abs(Q_closed(l)-Q_quad(l))/abs(Q_quad(l))) for l in [1.05,1.2,1.5,2.0])
out['Q_at_1_quad'] = float(Q_quad(mp.mpf('1.0000000001')))
out['Q_at_1_exact'] = -1/15
out['Q_at_1p5'] = float(Q_quad(1.5))

# ---------------- (C) a[eta_0 o Lambda^{-1}](0) purely in the w-variables --------------
# a = (3M/4) INT dlog|w| INT_0^pi lambda(rho(|w|,psi)) |cos psi| sin^2 psi dpsi   over Lambda(S)
# Change variable to rho at fixed psi:  dlog|w| = [1 - (dlog ghat/dlog lam) D(rho)] dlog rho
# ghat(psi;lam)^2 = lam^-2 sin^2 psi + lam^4 cos^2 psi   ( = |T_lam^{-1}w|^2/|w|^2 )
def Gfun(psi, l):
    s2, c2 = math.sin(psi)**2, math.cos(psi)**2
    gh2 = l**-2*s2 + l**4*c2
    return (-l**-2*s2 + 2*l**4*c2)/gh2               # dlog ghat / dlog lambda
def a_w_route(L, n=1200, m=1200):
    # sigma-Gauss-Legendre x psi-Gauss-Legendre
    xs, ws = np.polynomial.legendre.leggauss(n)
    ss = 0.5*(xs+1); wss = 0.5*ws                       # sigma in [0,1]
    xp, wp = np.polynomial.legendre.leggauss(m)
    ps = 0.5*(xp+1)*math.pi; wps = 0.5*math.pi*wp       # psi in [0,pi]
    tot = 0.0
    for s, wsig in zip(ss, wss):
        l = lam_s(s); D = dlogl_s(s)/L                  # dlog lam/dlog rho
        inner = 0.0
        for p, wpp in zip(ps, wps):
            inner += wpp*l*abs(math.cos(p))*math.sin(p)**2*(1 - Gfun(p,l)*D)
        tot += wsig*L*inner
    return 0.75*tot
def a_x_route(L, n=1200, m=1200):
    # the seat's route: (3M/4) INT dlogrho INT |cos| sin^2 g^-5 (1 + D W) dphi
    xs, ws = np.polynomial.legendre.leggauss(n)
    ss = 0.5*(xs+1); wss = 0.5*ws
    xp, wp = np.polynomial.legendre.leggauss(m)
    ps = 0.5*(xp+1)*math.pi; wps = 0.5*math.pi*wp
    tot = 0.0
    for s, wsig in zip(ss, wss):
        l = lam_s(s); D = dlogl_s(s)/L
        inner = 0.0
        for p, wpp in zip(ps, wps):
            gg = math.sqrt(l**2*math.sin(p)**2 + l**-4*math.cos(p)**2)
            W = math.sin(p)**2 - 2*math.cos(p)**2
            inner += wpp*abs(math.cos(p))*math.sin(p)**2/gg**5*(1 + D*W)
        tot += wsig*L*inner
    return 0.75*tot
rows = []
for L in (8.317766166719343, 20.0, 50.0):
    aw = a_w_route(L); ax = a_x_route(L); aref = 0.5*L*LB
    rows.append(dict(L=L, a_w=aw, a_x=ax, a_ref=aref,
                     offset_w=aw-aref, offset_x=ax-aref,
                     rel_w_x=abs(aw-ax)/abs(ax)))
out['a_rows'] = rows
out['max_rel_w_vs_x'] = max(r['rel_w_x'] for r in rows)
out['offset_spread_over_L'] = (max(r['offset_w'] for r in rows) - min(r['offset_w'] for r in rows))

# closed-form offset (3M/2) INT_0^1 dlog lam/dsigma * Q(lam) dsigma
off_closed = float(1.5*mp.quad(lambda s: mp.mpf(dlogl_s(float(s)))*Q_closed(lam_s(float(s))), [0,1]))
out['offset_closed'] = off_closed
out['offset_rel_times_L'] = 2*off_closed/LB
out['conservatism'] = out['C_rel_closed']/out['offset_rel_times_L']

# ---------------- (D) the (T') bound at Phi = Lambda, rebuilt ----------------
bnd = []
for L in (8.317766166719343, 20.0, 50.0, 100.0, 400.0):
    muJ = (2*math.sqrt(6)-4)/L
    bound = math.pi*L*LB*(3*0*(1+muJ)/(2*(1-0)**5) + 3*muJ/8)
    aref = 0.5*L*LB
    bnd.append(dict(L=L, muJ=muJ, bound=bound, rel=bound/aref, rel_times_L=bound/aref*L,
                    true_rel_times_L=off_closed/(0.5*LB)))
out['bound_rows'] = bnd

# ---------------- (E) P_h(3/2) for the taper, own quadrature ----------------
def P_h(l, deg):
    l = mp.mpf(l); d = mp.mpf(deg)*mp.pi/180
    f = lambda p: 3*min(mp.mpf(1), p/d)*mp.cos(p)*mp.sin(p)**2/g2(p,l)**mp.mpf(2.5)
    return float(mp.quad(f, [0, min(d, mp.pi/2), mp.pi/2]))
out['P_h_1p5'] = {str(d): P_h(1.5, d) for d in (3,5,7.5,10,15,20,30)}
out['P_h_1']   = {str(d): P_h(1.0, d) for d in (3,5,7.5,10,15,20,30)}
out['kappa_delta_7p5'] = P_h(1.0, 7.5)/2

json.dump(out, open(__file__.replace('x1_independent.py','x1_results.json'),'w'), indent=1)
for k,v in out.items():
    if not isinstance(v,(list,dict)): print(f"{k:28s} {v!r}")
print("\n-- a[eta o Lambda^-1](0): w-route (no det D Lambda) vs x-route --")
for r in rows: print("  L=%-20g a_w=%.12f a_x=%.12f offset_w=%.15f rel=%.3e" %
                     (r['L'], r['a_w'], r['a_x'], r['offset_w'], r['rel_w_x']))
print("\n-- (T') at Phi=Lambda --")
for r in bnd: print("  L=%-20g muJ=%.9f bound/a_ref*L=%.12f  true*L=%.12f" %
                    (r['L'], r['muJ'], r['rel_times_L'], r['true_rel_times_L']))
print("\nP_h(3/2):", {k: round(v,6) for k,v in out['P_h_1p5'].items()})

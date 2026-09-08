#!/usr/bin/env python3
"""
u3 -- assembly of the (Gamma-off) constant C'' and the exponent p, plus the
      window fixed point that decides the smallest L at which hypothesis (i)
      of B(t) closes on its own.

Inputs (all produced by u1/u2 in this folder, none typed from memory):
    C_K pi^4/2 = 3 pi^2/16   (u1)
    E0 = sup rho|eta_0|/M , G0 = sup rho^2|grad eta_0|/M   (u2)
    far/near kernel constants 0.291999 (far, z-odd), 0.014754 (inner, z-odd),
    3.999218/sin(phi) + pi/8 (collar) -- recomputed here from their own
    definitions (Cauchy-Schwarz/Parseval series and the bathtub constant), so
    that nothing is imported.

The chain (all on the slab S_s, see PROOF.md sec. 2):
    Gamma <= 2|a| + r|grad a| + |omega^theta|                       (B1, exact)
    |a(x,s)| <= a(0,s) + Chat_a M
    r|grad a| <= Ghat := (3 pi^2/16) e^{p c_G} G0 M
    |omega^theta| <= lambda M
  => Gamma <= 2 a(0,s) + C'' M ,  C'' = 2 Chat_a + lambda + Ghat/M
with p = 1 + 2 Gamma_rad/Gammabar and c_G = Gammabar * tau, tau = c/(M L).
"""
import json, math
import numpy as np
from scipy import integrate

u1 = json.load(open('u1_results.json'))
u2 = json.load(open('u2_results.json'))

OUT = {}
S3 = 2.0*math.pi**2
S4 = 8.0*math.pi**2/3.0
C_K = 3.0/(8.0*math.pi**2)
RIESZ = 3.0*math.pi**2/16.0                      # = C_K * pi^4/2, checked in u1

E0 = u2['E0_axis_closed']
G0 = u2['datum_L10']['G0']
OUT['E0'] = E0; OUT['G0'] = G0; OUT['riesz_C_K_pi4_2'] = RIESZ

# ---- far/near kernel constants, recomputed here from their definitions --------
# N_l = (l+1)(l+2)/(l+3/2) ; ||C_l^{3/2}||_inf = (l+1)(l+2)/2
def N_l(l):    return (l+1.0)*(l+2.0)/(l+1.5)
def Cmax(l):   return (l+1.0)*(l+2.0)/2.0

def Q_zodd(u, lmax=400):
    """far Taylor remainder envelope, z-odd (only odd l, remainder starts at l=3)."""
    s = 0.0
    for l in range(3, lmax+1, 2):
        s += ((l+2.0)/(2*l+3.0))**2*Cmax(l-1)**2*u**(2*(l-1))/N_l(l)
    return math.sqrt(s)

def S_zodd(v, lmax=400):
    """inner exterior-multipole envelope, z-odd (l odd, leading v^5)."""
    s = 0.0
    for l in range(1, lmax+1, 2):
        s += ((l+1.0)/(2*l+3.0))**2*Cmax(l+1)**2*v**(2*(l+4))/N_l(l)
    return math.sqrt(s)

C1_far, _ = integrate.quad(lambda u: math.sqrt(2.0)*Q_zodd(u)/u, 0.0, 0.5,
                           limit=300, epsabs=1e-13, epsrel=1e-13)
C2_in, _ = integrate.quad(lambda v: math.sqrt(2.0)*S_zodd(v)/v, 0.0, 0.5,
                          limit=300, epsabs=1e-13, epsrel=1e-13)
R_A = (2.0**5 - 2.0**-5)**0.2
COLLAR_COEF = 2.0*R_A                            # = 4 R_A/2 = 3.999218...
OUT['C1_far_zodd'] = C1_far
OUT['C2_inner_zodd'] = C2_in
OUT['collar_coef'] = COLLAR_COEF
OUT['collar_axis_piece_pi_over_8'] = math.pi/8
OUT['record_far_near'] = dict(C1=0.291999, C2in=0.014754, collar=3.999218, axis=math.pi/8)

def Chat_a(lam, sigma_star, cG):
    """|a(x,s)| <= a(0,s) + Chat_a * M  on {sin phi >= sigma_star}."""
    bulk = (C1_far + C2_in)*lam
    collar_A = (math.inf if sigma_star <= 0 else
                (COLLAR_COEF/sigma_star + math.pi/8)*lam)         # |omega|<=lam M route
    collar_B = COLLAR_COEF*math.exp(cG)*E0                        # P1 route (axis-safe)
    return bulk + min(collar_A, collar_B), collar_A, collar_B

# ---- the window fixed point --------------------------------------------------
C_WINDOW = 2*math.log(1.5)          # frozen window constant 2 log(3/2)
OUT['c_window_frozen'] = C_WINDOW

def solve(L, lam, sigma_star, c=C_WINDOW, mode='proved', mu=0.0, tol=1e-12, itmax=4000):
    """
    Fixed point in (Ghat, c_G, p).  a(0,s) <= (lam/2) M L  (kappa <= M_s/2).
    mode 'proved'  : radial factor e^{p c_G}, p from the symmetric-part bound
    mode 'p3'      : radial factor e^{3 c_G}  (the crude |b| <= Gamma rho bound)
    mode 'map'     : radial factor ((1+mu) lam)^2 e^{c_G}  (nu = 0, uses the map)
    Returns None if no fixed point (the iteration diverges).
    """
    a0 = 0.5*lam*L                    # a(0,s)/M upper bound
    Ghat = 0.0; cG = lam*c; p = 3.0
    for _ in range(itmax):
        Ca, _, _ = Chat_a(lam, sigma_star, cG)
        if mode == 'p3':
            fac = math.exp(3*cG)
        elif mode == 'map':
            fac = ((1+mu)*lam)**2*math.exp(cG)
        else:
            fac = math.exp(p*cG)
        if fac > 1e250:
            return None
        Gn = RIESZ*fac*G0
        Cpp = 2*Ca + lam + Gn
        Gbar = 2*a0 + Cpp             # in units of M
        cGn = Gbar*c/L
        A = a0 + Ca
        Grad = max(A, A + 2*Gn + lam/2, 2*Ca + 2*Gn + lam/2)
        pn = 1.0 + 2.0*min(Grad, Gbar)/Gbar
        if (abs(Gn-Ghat) <= tol*(1+abs(Gn)) and abs(cGn-cG) <= tol and abs(pn-p) <= tol):
            Ghat, cG, p = Gn, cGn, pn
            break
        Ghat, cG, p = Gn, cGn, pn
        if cG > 200:
            return None
    else:
        return None
    Ca, colA, colB = Chat_a(lam, sigma_star, cG)
    Cpp = 2*Ca + lam + Ghat
    return dict(L=L, lam=lam, sigma_star=sigma_star, mode=mode, mu=mu,
                Chat_a=Ca, collar_A=colA, collar_B=colB,
                Ghat_over_M=Ghat, C2prime=Cpp, c_G=cG, p=p,
                a0_over_M=0.5*lam*L, C1_multiplicative=0.0)

def Lstar(lam, sigma_star, mode='proved', mu=0.0, lo=10.0, hi=1e9):
    """smallest L at which the fixed point exists (bisection on existence)."""
    if solve(hi, lam, sigma_star, mode=mode, mu=mu) is None:
        return None
    for _ in range(200):
        mid = math.sqrt(lo*hi)
        if solve(mid, lam, sigma_star, mode=mode, mu=mu) is None:
            lo = mid
        else:
            hi = mid
        if hi/lo < 1 + 1e-10:
            break
    return hi

# ---- tables ------------------------------------------------------------------
tab = []
for mode, mu in [('p3', 0.0), ('proved', 0.0), ('map', 0.0), ('map', 0.02), ('map', 0.05)]:
    for lam in (1.0, 1.25, 1.5):
        for sig in (0.5, math.sin(math.radians(7.5)), 0.0):
            Ls = Lstar(lam, sig, mode=mode, mu=mu)
            row = dict(mode=mode, mu=mu, lam=lam, sigma_star=sig, Lstar=Ls)
            for L in (1e3, 1e4, 1e5):
                r = solve(L, lam, sig, mode=mode, mu=mu)
                row[f'L{int(L)}'] = None if r is None else dict(
                    C2prime=r['C2prime'], c_G=r['c_G'], p=r['p'], Ghat=r['Ghat_over_M'])
            tab.append(row)
OUT['table'] = tab

# the headline row: mu = 0.05, lambda = 3/2, sigma_* = 1/2
# crossover sigma at which the axis-safe collar (P1) takes over from the |omega|<=lam M route
def sigma_cross(lam, cG):
    return COLLAR_COEF*lam/(COLLAR_COEF*math.exp(cG)*E0 - math.pi*lam/8)
OUT['sigma_crossover_lam1.5_cG1.2164'] = sigma_cross(1.5, 1.5*C_WINDOW)
OUT['sigma_crossover_lam1.0_cG0.8109'] = sigma_cross(1.0, 1.0*C_WINDOW)

head = {}
for mode, mu in [('p3', 0.0), ('proved', 0.0), ('map', 0.05)]:
    Ls = Lstar(1.5, 0.5, mode=mode, mu=mu)
    head[f'{mode}_mu{mu}'] = dict(Lstar=Ls,
                                  at_10Lstar=solve(10*Ls, 1.5, 0.5, mode=mode, mu=mu) if Ls else None,
                                  at_Lstar_x2=solve(2*Ls, 1.5, 0.5, mode=mode, mu=mu) if Ls else None)
OUT['headline'] = head

# the fixed point at the existence boundary (used for the slack quote in PROOF sec 9)
for sg, tag in ((0.5, 'half'), (0.0, 'axis')):
    Ls = Lstar(1.5, sg, mode='proved')
    OUT[f'at_boundary_lam1.5_{tag}'] = dict(Lstar=Ls,
                                            **solve(Ls*(1+1e-7), 1.5, sg, mode='proved'))

# the "no-feedback" reference: what C'' would be if c_G were frozen at its L->inf value
for lam in (1.0, 1.25, 1.5):
  for sg, tag in ((0.5, ''), (math.sin(math.radians(7.5)), '_sin7'), (0.0, '_axis')):
    cG = lam*C_WINDOW
    Ca, colA, colB = Chat_a(lam, sg, cG)
    for nm, fac in [('p3', math.exp(3*cG)), ('p2', math.exp(2*cG)),
                    ('map_mu0.00', ((1.00)*lam)**2*math.exp(cG)),
                    ('map_mu0.02', ((1.02)*lam)**2*math.exp(cG)),
                    ('map_mu0.05', ((1.05)*lam)**2*math.exp(cG))]:
        OUT[f'Linf_lam{lam}{tag}_{nm}'] = dict(c_G=cG, Chat_a=Ca,
                                          Ghat=RIESZ*fac*G0,
                                          C2prime=2*Ca + lam + RIESZ*fac*G0)

with open('u3_results.json', 'w') as f:
    json.dump(OUT, f, indent=1)

print(f"C1_far(z-odd) = {C1_far:.6f}   [record 0.291999]")
print(f"C2_inner      = {C2_in:.6f}   [record 0.014754]")
print(f"collar coeff  = {COLLAR_COEF:.6f}  [record 3.999218]")
print(f"E0 = {E0:.6f}   G0 = {G0:.6f}   3pi^2/16 = {RIESZ:.7f}")
print()
for k in sorted(OUT):
    if k.startswith('Linf_'):
        print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in OUT[k].items()})
print()
print('--- sigma_* = 0 (axis-inclusive, whole shell) ---')
for mode in ('p3', 'proved', 'map'):
    for row in tab:
        if row['mode'] == mode and row['sigma_star'] == 0.0:
            print(f"mode={row['mode']:7s} mu={row['mu']:.2f} lam={row['lam']:.2f} "
                  f"L* = {row['Lstar'] if row['Lstar'] is None else round(row['Lstar'],1)}  "
                  f"L=1e4 -> {None if row['L10000'] is None else {kk: round(vv,4) for kk,vv in row['L10000'].items()}}")
print('--- sigma_* = 1/2 ---')
for mode in ('p3', 'proved', 'map'):
    for row in tab:
        if row['mode'] == mode and abs(row['sigma_star']-0.5) < 1e-9:
            print(f"mode={row['mode']:7s} mu={row['mu']:.2f} lam={row['lam']:.2f} "
                  f"L* = {row['Lstar'] if row['Lstar'] is None else round(row['Lstar'],1)}  "
                  f"L=1e4 -> {None if row['L10000'] is None else {kk: round(vv,4) for kk,vv in row['L10000'].items()}}")

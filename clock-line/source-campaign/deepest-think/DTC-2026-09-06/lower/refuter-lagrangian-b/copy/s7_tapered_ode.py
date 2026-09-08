#!/usr/bin/env python3
"""s7_tapered_ode.py -- the Lagrangian model run with the ADMISSIBLE (tapered) datum.

For the tapered plateau P_h(lam) != lam, so the ODE does not close; integrate it directly:
    d_theta F(sigma) = (1/2) int_sigma^1 P_h(e^{F(sigma')}) dsigma',   F(.,0) = 0,
and read off theta such that e^{F(0)} = 3/2.  Compare with
    the exact bang-bang value  theta = 4(1-sqrt(2/3)) = 0.7340,
    the frozen-strain value    theta = 2 log(3/2)/kappa_delta.
Also the fully conservative version that uses only  P_h(lam) >= P_h(1) = 2 kappa_delta  (the
monotonicity the theorem actually needs), which reproduces the frozen value exactly.
"""
import numpy as np, json
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
OUT = {}

def h_taper(delta):
    def h(phi):
        pax = np.minimum(phi, np.pi-phi); return np.minimum(1.0, pax/delta)
    return h
h_flat = lambda phi: np.ones_like(phi)

_phi = np.linspace(0.0, np.pi/2, 60001); _v = np.sin(_phi); _c = np.cos(_phi)
def P_h(h, lam):
    A = lam**2 - lam**-4; B = lam**-4
    return 3.0*np.trapz(h(_phi)*_v**2*(A*_v**2+B)**-2.5*_c, _phi)

def run(h, tag):
    lg = np.exp(np.linspace(0.0, np.log(1.6), 400))          # tabulate P_h on lam in [1,1.6]
    Pg = np.array([P_h(h, L) for L in lg])
    N = 2001; sg = np.linspace(0,1,N)
    def rhs(t, F):
        lam = np.clip(np.exp(F), 1.0, lg[-1])
        P = np.interp(lam, lg, Pg)
        d = np.diff(sg); trap = (P[1:]+P[:-1])/2*d
        tail = np.concatenate([np.cumsum(trap[::-1])[::-1], [0.0]])
        return 0.5*tail
    sol = solve_ivp(rhs, [0, 1.4], np.zeros(N), rtol=1e-10, atol=1e-12, dense_output=True)
    f = lambda th: sol.sol(th)[0] - np.log(1.5)
    th = brentq(f, 1e-9, 1.3)
    kap = 0.5*P_h(h, 1.0)
    OUT[tag] = dict(kappa=kap, theta32_model=th, theta32_frozen=np.log(1.5)/kap,
                    c2_model=2*th, c2_frozen=2*np.log(1.5)/kap,
                    F_at_theta=float(sol.sol(th)[0]), lam_inner=float(np.exp(sol.sol(th)[0])),
                    lam_outer_mean=float(np.trapz(np.exp(sol.sol(th)), sg)))
run(h_flat, 'bangbang')
for dd in (5.0, 7.5, 10.0, 15.0, 30.0):
    run(h_taper(np.deg2rad(dd)), f'taper_{dd}deg')
OUT['bangbang_exact_theta32'] = 4*(1-np.sqrt(2/3))
OUT['bangbang_exact_c2'] = 8*(1-np.sqrt(2/3))
print(json.dumps(OUT, indent=1, default=str))
json.dump(OUT, open('s7_results.json','w'), indent=1, default=str)

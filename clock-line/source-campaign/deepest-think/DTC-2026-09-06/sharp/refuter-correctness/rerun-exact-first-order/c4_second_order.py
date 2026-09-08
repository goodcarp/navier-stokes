#!/usr/bin/env python3
"""C4 -- the SECOND-order term at t=0 and the window on which first order is valid.

EXACT ALGEBRA (derived here, then checked numerically).  Let a = u^r/r.  Following a fluid particle,
dr/dt = u^r, so
      D_t a = D_t(u^r)/r - (u^r/r)^2 = D_t(u^r)/r - a^2 ,
and for axisymmetric NO-SWIRL flow the radial momentum equation has no centrifugal term:
      D_t u^r = -d_r p + nu (Delta - 1/r^2) u^r .
Hence, at a core (where omega^theta = r eta with eta materially conserved),
      d^2/dt^2 omega^theta = D_t(a omega^theta) = omega^theta (D_t a + a^2) = omega^theta * D_t(u^r)/r
                           = omega^theta * [ -(1/r) d_r p + (nu/r)(Delta - 1/r^2) u^r ].
So the SECOND-order coefficient is beta := D_t a + a^2 = -(1/r) d_r p  (inviscid):
  ** the a^2 "compounding" term is exactly cancelled by the -(u^r)^2/r^2 in D_t a.  The genuine
     second-order driver is the radial pressure gradient, not a^2. **
Numerically we obtain D_t a from the exact evolution of eta under its own flow, self-consistently:
  eta is transported, so each core keeps its peak value, translates with u(x_k), and its covariance
  evolves as C -> (I + tau G)C(I + tau G)^T with G = grad u at the core (exact to first order in tau,
  leading order in sigma/rho).  Then D_t a = [a(tau) - a(-tau)]/(2 tau) at the MOVING inner core.
Internal check: r_0 (D_t a + a^2) must equal D_t u^r computed the same way. Both are printed.
"""
import numpy as np, json, hashlib, time
import ns5d

M = 1.0
PHIS = np.arccos(1.0/np.sqrt(3.0))
NT, NAL = 24, 384
t0 = time.time()

def ring_at(rho, s, sign=-1.0):
    rk, zk = rho*np.sin(PHIS), rho*np.cos(PHIS)
    return ns5d.make_ring(rk, zk, sign*M/rk, (s*rho)**2*np.eye(2))

def a_f(r, z, rings):  return ns5d.a_field(r, z, rings, nt=NT, nal=NAL)
def psi_f(r, z, rings): return ns5d.psi_field(r, z, rings, nt=NT, nal=NAL)
def uz_f(r, z, rings, h):
    p = [psi_f(r + k*h, z, rings) for k in (-2, -1, 1, 2)]
    return 2.0*psi_f(r, z, rings) + r*(p[0] - 8*p[1] + 8*p[2] - p[3])/(12*h)
def ur_f(r, z, rings): return r*a_f(r, z, rings)

def grad_u(r, z, rings, h):
    d = {}
    for name, f in (('ur', ur_f), ('uz', lambda R, Z, RG: uz_f(R, Z, RG, h))):
        d[name+'_r'] = (f(r+h, z, rings) - f(r-h, z, rings))/(2*h)
        d[name+'_z'] = (f(r, z+h, rings) - f(r, z-h, rings))/(2*h)
    return np.array([[d['ur_r'], d['ur_z']], [d['uz_r'], d['uz_z']]])

def evolve(rings, tau, U, G, deform=True):
    out = []
    for R, u, g in zip(rings, U, G):
        Fm = np.eye(2) + tau*g
        C = Fm @ R['C'] @ Fm.T if deform else R['C']
        out.append(ns5d.make_ring(R['r'] + tau*u[0], R['z'] + tau*u[1], R['eta0'], C))
    return out

OUT = {}
S = 0.10
for N in [4, 8, 12]:
    rings = [ring_at(2.0**k, S) for k in range(N)]
    U, G = [], []
    for R in rings:
        h = 2e-3*R['sig']
        U.append(np.array([ur_f(R['r'], R['z'], rings), uz_f(R['r'], R['z'], rings, h)]))
        G.append(grad_u(R['r'], R['z'], rings, 0.05*R['sig']))
    r0, z0 = rings[0]['r'], rings[0]['z']
    a0 = a_f(r0, z0, rings)
    print(f"\nN={N}, s={S}:  a_inner/M = {a0:.8e}   u(core0) = ({U[0][0]:+.5e}, {U[0][1]:+.5e})")
    print(f"   grad u at core0 =\n{G[0]}")
    print(f"   incompressibility residual  d_r u^r + u^r/r + d_z u^z = {G[0][0,0] + U[0][0]/r0 + G[0][1,1]:+.3e}"
          f"   (scale {abs(G[0][0,0]):.3e})")
    res = {}
    for deform in (False, True):
        row = []
        for tau in [4e-3, 2e-3, 1e-3]:
            Rp = evolve(rings, +tau, U, G, deform); Rm = evolve(rings, -tau, U, G, deform)
            xp = np.array([r0, z0]) + tau*U[0]; xm = np.array([r0, z0]) - tau*U[0]
            ap = a_f(xp[0], xp[1], Rp); am = a_f(xm[0], xm[1], Rm)
            Dta = (ap - am)/(2*tau)
            urp = xp[0]*ap; urm = xm[0]*am
            Dtur = (urp - urm)/(2*tau)
            row.append(dict(tau=tau, Dta=Dta, beta=Dta + a0*a0, Dtur_over_r=Dtur/r0))
            print(f"   deform={deform}  tau={tau:.0e}  D_t a = {Dta:+.8e}   beta = D_t a + a^2 = {Dta+a0*a0:+.8e}"
                  f"   [check D_t(u^r)/r = {Dtur/r0:+.8e}, rel {abs(Dtur/r0-(Dta+a0*a0))/abs(Dta+a0*a0):.1e}]")
        res['deform' if deform else 'nodeform'] = row
    beta = res['deform'][-1]['beta']
    print(f"   a^2 = {a0*a0:.6e}   D_t a = {res['deform'][-1]['Dta']:+.6e}   beta = {beta:+.6e}")
    print(f"   -> beta/a^2 = {beta/a0**2:+.5f}    t* := 2a/|beta| = {2*a0/abs(beta):.6f}/M ;"
          f"  t_d = ln2/a = {np.log(2)/a0:.6f}/M ;  t*/t_d = {2*a0**2/(abs(beta)*np.log(2)):.5f}")
    OUT[f'N{N}'] = dict(a=float(a0), res={k: [{kk: float(vv) for kk, vv in r.items()} for r in v] for k, v in res.items()},
                        u0=[float(x) for x in U[0]], G0=[[float(x) for x in r] for r in G[0]],
                        beta=float(beta), t_star=float(2*a0/abs(beta)), t_d=float(np.log(2)/a0))

print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")
json.dump(OUT, open('c4_results.json','w'), indent=1, default=float)

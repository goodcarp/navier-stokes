#!/usr/bin/env python3
"""C5 -- does viscosity kill it?  Exact first-order competition at t = 0, ring by ring.

eta obeys D_t eta = nu Delta_5 eta.  At a core of width sigma << r the 3 sphere directions are flat and
Delta_5 acts as the 2D meridional Laplacian up to a relative O(sigma/r) curvature term, so the peak of a
Gaussian core obeys  sigma^2(t) = sigma_0^2 + 2 nu t,  peak(eta) prop sigma^{-2}, i.e.
      d/dt log eta_peak |_{t=0} = -2 nu / sigma_k^2 .
omega^theta = r eta at the core, and d/dt log r = a, so the NET first-order growth rate of the peak of
|omega^theta| on ring j in an N-ring stack is
      g_j(N) = a_j(N) - 2 nu / sigma_j^2 ,     a_j(N)/M = sum_{m=-j}^{N-1-j} A_m .
With the core-limited viscous floor sigma_0 = sqrt(nu/M)  (nu = M s^2 rho0^2):  2 nu/sigma_j^2 = 2 M 4^{-j}.
Then   d/dt sup|omega| |_{t=0} = M * max_j g_j(N)/M   (all cores carry the same peak at t = 0).
"""
import numpy as np, json, hashlib
C2 = json.load(open('c2_results.json'))
S_LIST = [0.04, 0.06, 0.08, 0.10, 0.125]
OUT = {}
NMAX = 26
print("net first-order growth rate g_j/M = a_j/M - 2*4^{-j}  (core-limited floor sigma_0 = sqrt(nu/M))")
for s in S_LIST:
    A = {int(m): v for m, v in C2['A_table'][str(s)].items()}
    kap = C2['kappa'][str(s)]; Z = C2['a_inner'][str(s)]['Z']
    def Am(m):
        if m > 25: return kap
        if m < -8: return 0.0
        return A[m]
    def a_j(j, N): return sum(Am(m) for m in range(-j, N-j))/Z
    print(f"\n  s={s:6.4f}   kappa={kap/Z:.6e}")
    print(f"   {'N':>3} " + " ".join(f"{'g'+str(j):>11}" for j in range(6)) + f" {'argmax j':>9} {'max g/M':>12} {'a_0/M':>12}")
    rows = []
    for N in [4, 8, 12, 16, 20, 26]:
        g = [a_j(j, N) - 2.0*4.0**(-j) for j in range(min(6, N))]
        gall = [a_j(j, N) - 2.0*4.0**(-j) for j in range(N)]
        jm = int(np.argmax(gall))
        rows.append(dict(N=N, g=g, jmax=jm, gmax=gall[jm], a0=a_j(0, N)))
        print(f"   {N:3d} " + " ".join(f"{x:11.4e}" for x in g) + f" {jm:9d} {gall[jm]:12.5e} {a_j(0,N):12.5e}")
    OUT[str(s)] = rows
    # threshold N for net growth at the innermost core
    Nstar = next((N for N in range(1, 400) if a_j(0, min(N, NMAX)) + kap/Z*max(0, N-NMAX) - 2.0 > 0), None)
    print(f"   inner core has NET growth (g_0>0) once a_0/M > 2, i.e. N >~ {2.0/(kap/Z):.1f}"
          f"  => log Re_E >~ {2*np.log(2)*2.0/(kap/Z):.1f}  (Re_E >~ 10^{2*np.log(2)*2.0/(kap/Z)/np.log(10):.1f})")
    OUT[str(s)+'_Nstar'] = float(2.0/(kap/Z))
print("\nSHA256", hashlib.sha256(open(__file__,'rb').read()).hexdigest())
json.dump(OUT, open('c5_results.json','w'), indent=1, default=float)

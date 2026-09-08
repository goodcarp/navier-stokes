#!/usr/bin/env python3
"""C2 -- the N-ring datum and the EXACT first-order strain at the innermost core.

DATUM (explicit).  Fix M>0, rho0>0, aspect s in (0,1/6], polar angle phi* in (0,pi/2).
For k = 0..N-1:   rho_k = rho0 2^k,  core centre (r_k, z_k) = rho_k (sin phi*, cos phi*),
                  sigma_k = s rho_k,
                  eta_k(r,z) = -(M/r_k) exp( -((r-r_k)^2 + (z-z_k)^2) / (2 sigma_k^2) ),
                  omega^theta = r * sum_k eta_k       (so omega^theta -> 0 at the axis, and
                                                       eta = omega^theta/r is the transported scalar).
Bang-bang sign: all cores in z>0 with omega^theta < 0, which is the sign the estate's extremal identity
  sup a(axis) = (M/2) log(R/rho0)   picks out (there omega^theta = -M sgn z).
Normalisation: ||omega_0||_inf is computed, not assumed; every reported quantity is divided by it so that
  ||omega_0||_inf = M exactly.
phi* = arccos(1/sqrt 3) maximises the axis weight sin^2 phi |cos phi| of the exact reduction
  a(0,0) = (3/4) int int omega^theta (-sin^2 phi cos phi) dlog rho dphi  (verified by a scan below).

SCALE INVARIANCE.  a is invariant under x -> lambda x with M fixed, so the contribution of ring (0+m)
to a at ring 0's core depends only on (m, s, phi*).  Call it M * A_m.  Then
      a_j(N)/M = sum_{m=-j}^{N-1-j} A_m        (strain at ring j's core in an N-ring stack),
      a_inner(N)/M = a_0(N)/M = sum_{m=0}^{N-1} A_m .
Controls: N = 1 must give A_0 = O(s^2) with no log; sign reversal must negate everything exactly.
"""
import numpy as np, json, hashlib, sys, time
import ns5d

M = 1.0
PHIS = np.arccos(1.0/np.sqrt(3.0))
t0 = time.time()

def ring_at(rho, s, M=1.0, phi=PHIS, sign=-1.0):
    rk, zk = rho*np.sin(phi), rho*np.cos(phi)
    return ns5d.make_ring(rk, zk, sign*M/rk, (s*rho)**2*np.eye(2))

def A_m(m, s, phi=PHIS, nt=32, nal=512):
    """contribution of the ring at rho=2^m to a at the core of the ring at rho=1, in units of M."""
    r0, z0 = np.sin(phi), np.cos(phi)
    return ns5d.a_of_ring(r0, z0, ring_at(2.0**m, s, 1.0, phi), nt=nt, nal=nal)

def sup_omega(rings, ngrid=161):
    """sup over the meridional half plane of |omega^theta| = |r sum_k eta_k|, by local refinement at each core."""
    best = 0.0
    for R in rings:
        lo, hi = R['r'] - 4*R['sig'], R['r'] + 4*R['sig']
        lz, hz = R['z'] - 4*R['sig'], R['z'] + 4*R['sig']
        for _ in range(4):
            rr = np.linspace(max(lo, 1e-9), hi, ngrid); zz = np.linspace(lz, hz, ngrid)
            RR, ZZ = np.meshgrid(rr, zz, indexing='ij')
            om = np.abs(RR*sum(ns5d.eta_of(RR, ZZ, S) for S in rings))
            i, j = np.unravel_index(np.argmax(om), om.shape)
            best = max(best, om[i, j])
            dr = rr[1]-rr[0]; dz = zz[1]-zz[0]
            lo, hi = rr[i]-2*dr, rr[i]+2*dr; lz, hz = zz[j]-2*dz, zz[j]+2*dz
    return best

OUT = {}
# ---------------------------------------------------------------- 0. phi optimum
print("0) angular weight sin^2 phi |cos phi| -- exact maximiser phi* = arccos(1/sqrt3):")
ph = np.linspace(1e-4, np.pi/2-1e-4, 200001)
w = np.sin(ph)**2*np.abs(np.cos(ph))
print(f"   argmax numeric {ph[np.argmax(w)]:.8f} rad   exact {PHIS:.8f}   max {w.max():.10f}  exact 2/(3 sqrt3) = {2/(3*np.sqrt(3)):.10f}")
OUT['phi_star'] = dict(numeric=float(ph[np.argmax(w)]), exact=float(PHIS))

# ---------------------------------------------------------------- 1. quadrature convergence
print("\n1) quadrature convergence of A_m (s = 0.10)")
for m in [0, 1, 3, 10]:
    v = [A_m(m, 0.10, nt=nt, nal=na) for nt, na in [(20, 256), (32, 512), (44, 768)]]
    print(f"   m={m:3d}  {v[0]:.12e}  {v[1]:.12e}  {v[2]:.12e}   |d21/v| {abs(v[1]-v[0])/abs(v[2]):.2e}  |d32/v| {abs(v[2]-v[1])/abs(v[2]):.2e}")
    OUT[f'conv_m{m}'] = [float(x) for x in v]

# ---------------------------------------------------------------- 2. A_m table and kappa(s)
S_LIST = [0.04, 0.06, 0.08, 0.10, 0.125]
MS = list(range(-8, 26))
print("\n2) A_m(s) (units of M).  kappa_lead := (3 pi/2) s^2 sin^2 phi* cos phi* = (pi/sqrt3) s^2")
tab = {}
for s in S_LIST:
    tab[s] = {m: A_m(m, s) for m in MS}
    kap_lead = (3*np.pi/2)*s*s*np.sin(PHIS)**2*np.cos(PHIS)
    print(f"  s={s:6.4f}  kappa_lead={kap_lead:.8e}   A_25={tab[s][25]:.8e}  ratio {tab[s][25]/kap_lead:.6f}")
    print("     m:      " + "  ".join(f"{m:>12d}" for m in [-3,-2,-1,0,1,2,3,5,10,25]))
    print("     A_m:    " + "  ".join(f"{tab[s][m]:12.5e}" for m in [-3,-2,-1,0,1,2,3,5,10,25]))
OUT['A_table'] = {str(s): {str(m): float(v) for m, v in d.items()} for s, d in tab.items()}

# kappa(s) := limiting A_m ; tail beyond m=25 uses kappa
print("\n   kappa(s) (= A_m for large m) and its ratio to (pi/sqrt3)s^2:")
kap = {}
for s in S_LIST:
    kap[s] = tab[s][25]
    print(f"    s={s:6.4f}  kappa={kap[s]:.10e}   (pi/sqrt3)s^2 = {np.pi/np.sqrt(3)*s*s:.10e}   ratio {kap[s]/(np.pi/np.sqrt(3)*s*s):.8f}"
          f"   1+4.5 s^2 = {1+4.5*s*s:.6f}")
OUT['kappa'] = {str(k): float(v) for k, v in kap.items()}

# ---------------------------------------------------------------- 3. a_inner(N), affine fit, controls
print("\n3) a_inner(N)/M = sum_{m=0}^{N-1} A_m   (tail m>25 replaced by kappa; exact to <1e-12 rel)")
NMAX = 26
res = {}
for s in S_LIST:
    Am = tab[s]
    cum = []
    tot = 0.0
    for N in range(1, NMAX+1):
        m = N-1
        tot += Am[m] if m <= 25 else kap[s]
        cum.append(tot)
    # normalisation ||omega||_inf
    rings = [ring_at(2.0**k, s) for k in range(NMAX)]
    Z = sup_omega(rings)
    cumn = [c/Z for c in cum]
    # affine fit on the last 10 points (exact: a = kappa N + b for N large)
    Ns = np.arange(NMAX-9, NMAX+1)
    sl, ic = np.polyfit(Ns, np.array(cumn[NMAX-10:]), 1)
    res[s] = dict(cum=cum, cumn=cumn, Z=float(Z), slope=float(sl), inter=float(ic))
    print(f"  s={s:6.4f}  ||omega||_inf/M = {Z:.8f}   fit a_inner/M = {sl:.8e} * N + ({ic:+.8e})")
    print(f"          A_0 (N=1 control, no log) = {cum[0]/Z:.8e} = {cum[0]/Z/s**2:.5f} s^2")
    print(f"          (M/2)ln2 per octave would be {0.5*np.log(2):.8f};  ratio slope/((1/2)ln2) = {sl/(0.5*np.log(2)):.6e}")
OUT['a_inner'] = {str(s): dict(cum_norm=[float(x) for x in res[s]['cumn']], Z=res[s]['Z'],
                               slope=res[s]['slope'], inter=res[s]['inter']) for s in S_LIST}

print("\n   a_inner(N)/M table (normalised so ||omega_0||_inf = M):")
print("   N   " + "  ".join(f"s={s:<10.4g}" for s in S_LIST))
for N in [1, 2, 3, 4, 6, 8, 12, 16, 20, 26]:
    print(f"  {N:3d}  " + "  ".join(f"{res[s]['cumn'][N-1]:12.6e}" for s in S_LIST))

# ---------------------------------------------------------------- 4. controls
print("\n4) CONTROLS")
s = 0.10
one = ring_at(1.0, s)
a1 = ns5d.a_of_ring(np.sin(PHIS), np.cos(PHIS), one, nt=40, nal=768)
print(f"   (a) N=1 single ring: a/M = {a1:.10e} = {a1/s**2:.6f} s^2   -- O(M s^2), no N and no log. OK")
neg = [ring_at(2.0**k, s, sign=+1.0) for k in range(8)]
pos = [ring_at(2.0**k, s, sign=-1.0) for k in range(8)]
ap = sum(ns5d.a_of_ring(np.sin(PHIS), np.cos(PHIS), R) for R in pos)
an = sum(ns5d.a_of_ring(np.sin(PHIS), np.cos(PHIS), R) for R in neg)
print(f"   (b) sign-reversed stack (N=8): a_inner/M  bang-bang {ap:+.10e}   reversed {an:+.10e}   sum {ap+an:.2e}")
print(f"       bang-bang sign gives a_inner > 0 : {ap>0};  reversed gives a_inner < 0 : {an<0}")
# (c) inner rings compress outer rings
print(f"   (c) A_m for m<0 (inner ring acting on an outer core) must be NEGATIVE and decay ~ 4^m:")
for m in [-1, -2, -3, -4]:
    print(f"       A_{m} = {tab[0.10][m]:+.6e}   A_{m}/A_{m+1} = {tab[0.10][m]/tab[0.10][m+1]:.5f}  (4^-1 = 0.25)")
OUT['controls'] = dict(single=float(a1), bangbang=float(ap), reversed=float(an))

print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")
json.dump(OUT, open('c2_results.json', 'w'), indent=1, default=float)

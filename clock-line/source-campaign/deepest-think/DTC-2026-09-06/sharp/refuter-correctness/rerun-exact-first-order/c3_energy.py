#!/usr/bin/env python3
"""C3 -- exact energy, Reynolds number and the doubling-time constant c1 for the N-ring datum.

E := ||u||_2^2 (no 1/2) = 2 pi int int omega^theta Psi dr dz = 2 pi int int eta psi r^3 dr dz
                        = (1/pi) int_{R^5} eta psi  (derived and Hill-verified in c1).
Bilinear:  E = sum_{k,l} E_kl,  E_kl = 2 pi int int eta_k psi_l r^3 dr dz  (symmetric).
Scale invariance: scaling the pair (ring k, ring l) by lambda scales E_kl by lambda^5, so
      E_kl = M^2 rho_min(k,l)^5 h_{|k-l|}(s),
and    E(N) = M^2 rho0^5 [ sum_{k=0}^{N-1} 32^k h_0 + 2 sum_{k} sum_{j=1}^{N-1-k} 32^k h_j ].

VISCOUS FLOOR.  Two conventions, both reported:
  (I)  core-limited (physical): sigma_0 = sqrt(nu/M)  =>  nu = M s^2 rho0^2   -- the core is what
       viscosity smooths, so the SMALLEST length in the datum must sit at the viscous scale;
  (II) brief's:    rho0 = sqrt(nu/M)  =>  nu = M rho0^2.
Re_E := E^{2/5} M^{1/5} / nu  (= M l^2/nu with l = E^{1/5} M^{-2/5}).
First-order doubling time t_d = ln2 / a_inner ;   c1 := t_d * M * log Re_E.
"""
import numpy as np, json, hashlib, time
import ns5d
from ns5d import J3, gl

M = 1.0
PHIS = np.arccos(1.0/np.sqrt(3.0))
t0 = time.time()

def ring_at(rho, s, M=1.0, phi=PHIS, sign=-1.0):
    rk, zk = rho*np.sin(phi), rho*np.cos(phi)
    return ns5d.make_ring(rk, zk, sign*M/rk, (s*rho)**2*np.eye(2))

def gh(n):
    x, w = np.polynomial.hermite_e.hermegauss(n)     # weight exp(-x^2/2), int = sqrt(2 pi)
    return x, w

def outer_nodes(ring, n):
    """tensor Gauss-Hermite nodes for integrating f against ring's Gaussian:
       int f(r,z) exp(-|d|^2/(2 sig^2)) dr dz = sig^2 sum_ij w_i w_j f(r_k+sig x_i, z_k+sig x_j)."""
    x, w = gh(n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    W = np.outer(w, w)*ring['sig']**2
    return (ring['r'] + ring['sig']*X).ravel(), (ring['z'] + ring['sig']*Y).ravel(), W.ravel()

def psi_at_polar(r0, z0, ring, nt=26, nal=256):
    return ns5d.psi_of_ring(r0, z0, ring, nt=nt, nal=nal)

def E_pair(ringA, ringB, nout=18, self_term=False, nt=26, nal=256, ninGH=40):
    """E_AB = 2 pi int int eta_A psi_B r^3 dr dz."""
    rs, zs, ws = outer_nodes(ringA, nout)
    good = rs > 0
    tot = 0.0
    if self_term:
        for r0, z0, w in zip(rs[good], zs[good], ws[good]):
            tot += w*ringA['eta0']*psi_at_polar(r0, z0, ringB, nt=nt, nal=nal)*r0**3
    else:
        # inner integral by tensor Gauss-Hermite on ringB (kernel smooth there: cores are disjoint)
        rb, zb, wb = outer_nodes(ringB, ninGH)
        gb = rb > 0
        rb, zb, wb = rb[gb], zb[gb], wb[gb]
        R0 = rs[good][:, None]; Z0 = zs[good][:, None]
        K = J3(R0, Z0, rb[None, :], zb[None, :])
        psi = ((1/(2*np.pi))*ringB['eta0']*(rb**3*wb)[None, :]*K).sum(axis=1)
        tot = np.sum(ws[good]*ringA['eta0']*psi*rs[good]**3)
    return 2*np.pi*tot

OUT = {}
S_LIST = [0.04, 0.06, 0.08, 0.10, 0.125]
JMAX = 6                       # h_j beyond this is < 1e-4 of h_0 and is added as a bounded remainder
print("h_j(s) = E_{k,k+j} / (M^2 rho_k^5)   [h_0 = self-energy coefficient]")
H = {}
for s in S_LIST:
    r0 = ring_at(1.0, s)
    h = {}
    h[0] = E_pair(r0, r0, nout=18, self_term=True)
    for j in range(1, JMAX+1):
        h[j] = E_pair(r0, ring_at(2.0**j, s), nout=18)
    H[s] = h
    thin = (2*np.pi*(s)**2*M)**2*np.sin(PHIS)*(np.log(8*np.sin(PHIS)/s) - 7/4)   # thin-ring reference (uniform core)
    print(f"  s={s:6.4f}  h0={h[0]:.8e}  (thin-ring uniform-core ref {thin:.3e}, ratio {h[0]/thin:.4f})"
          f"   h1/h0={h[1]/h[0]:+.5f} h2/h0={h[2]/h[0]:+.5f} h3/h0={h[3]/h[0]:+.5f} h4/h0={h[4]/h[0]:+.2e} h6/h0={h[6]/h[0]:+.2e}")
OUT['h'] = {str(s): {str(j): float(v) for j, v in h.items()} for s, h in H.items()}

# convergence check of h_0 and h_1
print("\n  convergence (s=0.10):")
r0 = ring_at(1.0, 0.10)
for nout in [12, 18, 24]:
    v0 = E_pair(r0, r0, nout=nout, self_term=True, nt=26, nal=256)
    v1 = E_pair(r0, ring_at(2.0, 0.10), nout=nout, ninGH=40)
    print(f"    nout={nout:3d}   h0={v0:.10e}   h1={v1:.10e}")
v0b = E_pair(r0, r0, nout=18, self_term=True, nt=38, nal=384)
print(f"    h0 with finer inner polar (nt=38,nal=384): {v0b:.10e}")

def E_of_N(N, s, rho0=1.0, Mv=1.0):
    h = H[s]
    tot = 0.0
    for k in range(N):
        tot += 32.0**k*h[0]
        for j in range(1, min(JMAX, N-1-k)+1):
            tot += 2*32.0**k*h[j]
    return Mv*Mv*rho0**5*tot

# ---------------------------------------------------------------- c1
print("\nc1 := t_d * M * log Re_E  with t_d = ln2/a_inner ;  a_inner from c2 (normalised ||omega||_inf = M)")
C2 = json.load(open('c2_results.json'))
rows = []
for s in S_LIST:
    ai = C2['a_inner'][str(s)]['cum_norm']       # a_inner(N)/M for N=1..26
    Z = C2['a_inner'][str(s)]['Z']
    kap = C2['a_inner'][str(s)]['slope']
    print(f"\n  s={s:6.4f}   kappa (per octave, normalised) = {kap:.8e}    [(1/2)ln2 = {0.5*np.log(2):.6f}]")
    print(f"   {'N':>3} {'a_inner/M':>13} {'E/(M^2 rho0^5)':>15} {'log Re_E (I)':>13} {'c1 (I)':>9} {'log Re_E(II)':>13} {'c1 (II)':>9}")
    for N in [2, 3, 4, 6, 8, 10, 12, 16, 20, 26]:
        a = ai[N-1]
        E = E_of_N(N, s)/Z**2                    # normalisation: omega -> omega/Z  =>  E -> E/Z^2
        for conv, nu in (('I', s*s), ('II', 1.0)):   # nu/(M rho0^2)
            ReE = E**0.4*1.0**0.2/nu             # M=1, rho0=1
            lR = np.log(ReE)
            c1 = (np.log(2)/a)*lR
            if conv == 'I': lR1, c11 = lR, c1
            else: lR2, c12 = lR, c1
        print(f"   {N:3d} {a:13.6e} {E:15.6e} {lR1:13.5f} {c11:9.4f} {lR2:13.5f} {c12:9.4f}")
        rows.append(dict(s=s, N=N, a=a, E=E, logRe_I=lR1, c1_I=c11, logRe_II=lR2, c1_II=c12))
    # asymptotic prediction c1 -> 2 (ln2)^2 / kappa
    print(f"   asymptotic  c1 -> 2(ln2)^2/kappa = {2*np.log(2)**2/kap:.5f}")
OUT['rows'] = rows
OUT['c1_asym'] = {str(s): float(2*np.log(2)**2/C2['a_inner'][str(s)]['slope']) for s in S_LIST}

# the ideal space-filling (bang-bang shell) endpoint, exact
print("\nIDEAL ENDPOINT (exact, no rings): bang-bang shell a = (M/2)log(R/rho0), E ~ M^2 R^5,")
print("  nu = M rho0^2  =>  Re_E = C^{2/5}(R/rho0)^2, log Re_E = 2 log(R/rho0) + O(1),")
print(f"  a = (M/4) log Re_E  =>  t_d = ln2/a = 4 ln2/(M log Re_E)  =>  c1 = 4 ln 2 = {4*np.log(2):.10f}")
OUT['c1_ideal'] = float(4*np.log(2))

print(f"\nelapsed {time.time()-t0:.1f}s  SHA256 {hashlib.sha256(open(__file__,'rb').read()).hexdigest()}")
json.dump(OUT, open('c3_results.json','w'), indent=1, default=float)

#!/usr/bin/env python3
"""C1 -- calibration of the 5D-lift Biot-Savart kernel of ns5d.py against EXACT analytic targets.
  (1) closed-form theta kernels vs mpmath quadrature, including the near-coincidence limit k2 -> 1;
  (2) HILL'S SPHERICAL VORTEX: eta = A on the 5D ball |x| < R  =>  a = u^r/r = A z / 5 at EVERY interior
      point.  Proof (independent of the run): a = -A d_z Phi with Phi the Newtonian potential of the unit
      ball in R^5, -Delta Phi = 1 inside => Phi = C - |x|^2/10 => a = A z/5.  On-axis AND off-axis points.
  (3) Hill's kinetic energy: E = int|u|^2 = (1/pi) int eta psi = 16 pi A^2 R^7/315 = (20/7) pi U^2 R^3
      with U = 2AR^2/15, i.e. (1/2)E = (10/7) pi U^2 R^3, the textbook value.
  (4) BANG-BANG SHELL axis identity a(0,0) = (M/2) log(R/rho0) for omega^theta = -M sgn(z) on rho0<|x|<R.
  (5) sign control.
Numerics falsify, they never prove.  Every target here is an exact closed form derived above.
"""
import numpy as np, mpmath as mp, hashlib, json, sys
import ns5d
from ns5d import J5, J3, gl

mp.mp.dps = 40
OUT = {}
FAIL = []

def chk(name, got, want, tol):
    rel = abs(got - want)/max(abs(want), 1e-300)
    OUT[name] = dict(got=float(got), want=float(want), rel=float(rel), tol=tol, pass_=bool(rel < tol))
    print(f"  {name:44s} got {got: .12e}  want {want: .12e}  rel {rel:.2e}  {'OK' if rel<tol else 'FAIL'}")
    if rel >= tol: FAIL.append(name)

print("(1) closed-form theta kernels, including near-coincidence")
rng = np.random.default_rng(11)
worst5 = worst3 = 0.0
for _ in range(40):
    r, rp = rng.uniform(0.3, 3.0, 2); z = 0.0; zp = rng.uniform(-2, 2)
    A = r*r + rp*rp + (z-zp)**2; B = 2*r*rp
    q5 = mp.quad(lambda th: mp.sin(th)**2/(A - B*mp.cos(th))**mp.mpf('2.5'), [0, mp.pi])
    q3 = mp.quad(lambda th: mp.sin(th)**2/(A - B*mp.cos(th))**mp.mpf('1.5'), [0, mp.pi])
    worst5 = max(worst5, float(abs(J5(r, z, rp, zp)-q5)/q5)); worst3 = max(worst3, float(abs(J3(r,z,rp,zp)-q3)/q3))
print(f"  random points: worst rel err  J5 {worst5:.2e}   J3 {worst3:.2e}")
OUT['kernel_random'] = dict(J5=worst5, J3=worst3)
# near coincidence: r'=r+t, z'=z, t -> 0 ; J5 -> 1/(3 r^3 t^2) (flat-space codim-2 limit)
print("  near-coincidence J5 vs flat-space limit 1/(3 r^3 t^2):")
rows = []
for t in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8]:
    r = 1.0
    num = J5(r, 0.0, r + t, 0.0)
    # reference in EXACT mp arithmetic, written in the well-conditioned form
    # A - B cos th = (r-r')^2 + 2 r r' (1-cos th),  r' = r+t
    tm = mp.mpf(t); rm = mp.mpf(r); rpm = rm + tm
    # integrand is peaked at th ~ t/r of width ~t: guide mp.quad with break points
    bp = [mp.mpf(0)] + [tm*mp.mpf(10)**k for k in range(0, 9) if tm*mp.mpf(10)**k < mp.pi] + [mp.pi]
    q = mp.quad(lambda th: mp.sin(th)**2/(tm**2 + 2*rm*rpm*(1-mp.cos(th)))**mp.mpf('2.5'), bp)
    rows.append(dict(t=t, J5=float(num), mp=float(q), rel=float(abs(num-q)/q), ratio=float(num*3*r**3*t*t)))
    print(f"    t={t:8.1e}  J5 {float(num):.8e}  mp {float(q):.8e}  rel {rows[-1]['rel']:.2e}   J5*3r^3t^2 = {rows[-1]['ratio']:.6f}")
OUT['near_coincidence'] = rows
if max(x['rel'] for x in rows) > 1e-10: FAIL.append('near_coincidence')

# ---------------------------------------------------------------- (2) HILL
print("\n(2) Hill's spherical vortex: eta = A on |x|<R  =>  a = A z/5 exactly inside")
A_H, R_H = 1.0, 1.0
def a_hill(r0, z0, nt=200, nal=1536):
    """polar quadrature about (r0,z0) over the meridional half-disc {r>0, r^2+z^2<R^2}."""
    al = 2.0*np.pi*np.arange(nal)/nal
    def tmax(a):
        ca, sa = np.cos(a), np.sin(a)
        proj = r0*ca + z0*sa
        disc = R_H*R_H - (r0*r0 + z0*z0) + proj*proj
        tc = -proj + np.sqrt(np.maximum(disc, 0.0))          # chord to the sphere
        ta = np.where(ca < -1e-14, r0/np.maximum(-ca, 1e-300), np.inf)   # chord to r'=0
        return np.minimum(tc, ta)
    tm = tmax(al)
    tot = 0.0
    for a, tmx in zip(al, tm):
        if tmx <= 0: continue
        # split the radial panel geometrically near 0 (kernel ~ 1/t^2) and uniformly out to tmx
        edges = [0.0] + [tmx*4.0**(-j) for j in range(6, 0, -1)] + [tmx]
        edges = sorted(set(edges))
        ts, ws = [], []
        for u, v in zip(edges[:-1], edges[1:]):
            x, w = gl(nt//6 + 4, u, v); ts.append(x); ws.append(w)
        t = np.concatenate(ts); w = np.concatenate(ws)
        rp = r0 + t*np.cos(a); zp = z0 + t*np.sin(a)
        m = rp > 0
        K = np.where(m, J5(r0, z0, np.where(m, rp, 1.0), zp), 0.0)
        integ = (3.0/(2*np.pi))*A_H*np.where(m, rp, 0.0)**3*(-t*np.sin(a))*K*t
        tot += (2*np.pi/nal)*np.sum(integ*w)
    return tot
for (r0, z0) in [(0.0, 0.30), (0.0, 0.70), (0.40, 0.30), (0.55, -0.45), (0.80, 0.10)]:
    chk(f"Hill a({r0},{z0}) = A z/5", a_hill(r0, z0), A_H*z0/5.0, 3e-6)

# ---------------------------------------------------------------- (3) Hill potential + energy identity
print("\n(3) Hill potential psi and the energy identity")
# EXACT: -Delta_5 psi = eta = A on |x|<R  =>  psi_in = A(R^2/6 - |x|^2/10), psi_out = A R^5/(15|x|^3)
def psi_hill(r0, z0, nin=160, nal=1024):
    al = 2*np.pi*np.arange(nal)/nal
    ca, sa = np.cos(al), np.sin(al)
    proj = r0*ca + z0*sa
    disc = R_H*R_H - (r0*r0 + z0*z0) + proj*proj
    tc = -proj + np.sqrt(np.maximum(disc, 0.0))
    ta = np.where(ca < -1e-14, r0/np.maximum(-ca, 1e-300), np.inf)
    tm = np.minimum(tc, ta)
    tot = 0.0
    for a, tmx in zip(al, tm):
        if tmx <= 0: continue
        x, w = gl(nin, 0.0, tmx)
        rp = r0 + x*np.cos(a); zp = z0 + x*np.sin(a)
        m = rp > 0
        K = np.where(m, J3(r0, z0, np.where(m, rp, 1.0), zp), 0.0)
        tot += (2*np.pi/nal)*np.sum((1/(2*np.pi))*A_H*np.where(m, rp, 0.0)**3*K*x*w)
    return tot
for (r0, z0) in [(0.30, 0.20), (0.60, -0.30), (0.05, 0.55)]:
    chk(f"Hill psi({r0},{z0}) = A(R^2/6-|x|^2/10)", psi_hill(r0, z0), A_H*(R_H**2/6 - (r0*r0+z0*z0)/10), 3e-6)
# energy identity, symbolically
import sympy as sp
Asy, Rsy, rho = sp.symbols('A R rho', positive=True)
psi_in = Asy*(Rsy**2/6 - rho**2/10)
Esym = sp.simplify((1/sp.pi)*sp.integrate(Asy*psi_in*sp.Rational(8,3)*sp.pi**2*rho**4, (rho, 0, Rsy)))  # |S^4| = 8pi^2/3
chk("E_Hill = 16 pi A^2 R^7/315 (symbolic)", float(Esym.subs({Asy:1, Rsy:1})), 16*np.pi/315, 1e-14)
Usy = 2*Asy*Rsy**2/15
chk("KE_Hill = (10/7) pi U^2 R^3 (symbolic)", float((Esym/2).subs({Asy:1,Rsy:1})), float(((sp.Rational(10,7)*sp.pi*Usy**2*Rsy**3).subs({Asy:1,Rsy:1}))), 1e-14)

# ---------------------------------------------------------------- (4) bang-bang shell axis identity
print("\n(4) bang-bang shell  omega^theta = -M sgn(z) on rho0<|x|<R:  a(0,0) = (M/2) log(R/rho0)")
def a_shell_axis(rho0, R, M=1.0, sign=-1.0, nrho=64, nphi=96):
    """a(0,0) = (3/4) int int omega^theta (-sin^2 phi cos phi) dlog(rho) dphi (exact reduction, r=0 => J=(pi/2)/rho^5)."""
    tot = 0.0
    edges = list(rho0*2.0**np.arange(0, np.log2(R/rho0)+1e-9))
    edges = sorted(set([rho0, R] + [e for e in edges if rho0 < e < R]))
    for u, v in zip(edges[:-1], edges[1:]):
        rg, wg = gl(nrho, np.log(u), np.log(v))
        for ph_a, ph_b in [(0.0, np.pi/2), (np.pi/2, np.pi)]:
            pg, wg2 = gl(nphi, ph_a, ph_b)
            om = sign*M*np.sign(np.cos(pg))
            tot += np.sum(wg)*np.sum(wg2*(3.0/4.0)*om*(-np.sin(pg)**2*np.cos(pg)))
    return tot
for R in [4.0, 4096.0, 2.0**20]:
    chk(f"shell a(0,0), R/rho0={R:g}", a_shell_axis(1.0, R), 0.5*np.log(R), 1e-12)
chk("shell sign control (omega=+M sgn z)", a_shell_axis(1.0, 4096.0, sign=+1.0), -0.5*np.log(4096.0), 1e-12)

# ---------------------------------------------------------------- (5) single Gaussian ring: far field
print("\n(5) single eta-Gaussian ring, far-field check against the exact delta-ring limit")
M = 1.0; phis = np.arccos(1/np.sqrt(3.0))
for s in [0.02, 0.05, 0.10]:
    rk, zk = np.sin(phis), np.cos(phis)
    ring = ns5d.make_ring(rk, zk, -M/rk, (s*s)*np.eye(2))
    got = ns5d.a_of_ring(0.0, 0.0, ring, nt=40, nal=768)
    want = (3*np.pi/2)*s*s*M*np.sin(phis)**2*np.cos(phis)      # leading delta-ring value at the origin
    print(f"  s={s:5.3f}  a(0,0) = {got:.10e}   leading (3pi/2)s^2 sin^2 cos = {want:.10e}   ratio {got/want:.6f}")
    OUT[f'single_ring_s{s}'] = dict(a=float(got), lead=float(want), ratio=float(got/want))

print("\nSHA256", hashlib.sha256(open(__file__,'rb').read()).hexdigest())
json.dump(OUT, open('c1_results.json','w'), indent=1, default=float)
print("FAILURES:", FAIL if FAIL else "none")
sys.exit(1 if FAIL else 0)

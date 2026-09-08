#!/usr/bin/env python3
"""G2 — does the log strain that a multi-scale axisymmetric vorticity field produces AT THE AXIS also act
ON THE FIELD'S OWN SUPPORT?  (the 'self-stretching' supposition, model level, numerics only)

Field: no swirl, omega^theta = -M sign(z) on the meridional annulus rho0 < |x| < R (the bang-bang extremal
field of the estate's identity sup a(axis) = (M/2) log(R/rho0)). Own instrument: 5D-lift Biot-Savart,
a = u^r/r = -d_z psi1, psi1 = G5 * eta, eta = omega^theta/r radial in R^4.
   a(r,z) = (3/(2 pi)) ∫∫ eta(r',z') r'^3 (z-z') J dr'dz',  J = ∫_0^pi sin^2(th) dth / (r^2+r'^2-2 r r' cos th+(z-z')^2)^{5/2}
Calibration (must pass): a(0,0) = (M/2) log(R/rho0) exactly (kernel angular factor 2/3).
Claim under test: for rho0 << |x| = rho << R inside the support, a(x) = (M/2) log(R/rho) + c(phi) M with c = O(1)
independent of rho — i.e. the outer shells stretch the inner shells at the SAME log rate they produce at the axis.
Controls: (i) a field with omega^theta = +M sign(z) must give the NEGATIVE of everything (sign control);
(ii) a single-octave field (R = 2 rho0) must give a(0) = (M/2) log 2 and NO log growth in rho.
Numerics falsify; they never prove. Exclusion ball of relative radius EPS around x (kernel is absolutely
integrable in 5D; dropped mass O(M EPS)). All floats; decision thresholds stated.
"""
import numpy as np, sys, hashlib, time, json
M = 1.0
def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n); return 0.5*(b-a)*x + 0.5*(b+a), 0.5*(b-a)*w
# theta nodes: three panels clustered at 0 (the near-singular direction)
TH = [gl(48, 0, np.pi/64), gl(48, np.pi/64, np.pi/8), gl(64, np.pi/8, np.pi)]
th = np.concatenate([p[0] for p in TH]); wth = np.concatenate([p[1] for p in TH])
s2 = np.sin(th)**2 * wth
def strain_at(r, z, rho0, R, sign=-1.0, eps=2e-2, nr=64, nphi=64):
    """a(r,z) for eta = sign*M*sgn(z')/r' on rho0<rho'<R. Panels refined around the evaluation point."""
    rho = np.hypot(r, z); phi = np.arctan2(r, z)  # phi measured from +z axis, in (0,pi)
    # rho' panels
    geo = list(rho0 * 2.0**np.arange(0, np.log2(R/rho0)+1e-9))  # one panel per octave over the whole support
    edges = sorted(set([rho0, R] + [v for v in geo if rho0 < v < R] + [v for v in [rho*0.5, rho*0.9, rho*0.99, rho*1.01, rho*1.1, rho*2.0] if rho0 < v < R]))
    rp = []; wrp = []
    for a_, b_ in zip(edges[:-1], edges[1:]):
        x_, w_ = gl(nr, a_, b_); rp.append(x_); wrp.append(w_)
    rp = np.concatenate(rp); wrp = np.concatenate(wrp)
    # phi' panels (around phi and around the interface pi/2)
    pedges = sorted(set([0.0, np.pi] + [v for v in [phi-0.3, phi-0.05, phi+0.05, phi+0.3, np.pi/2] if 0 < v < np.pi]))
    pp = []; wpp = []
    for a_, b_ in zip(pedges[:-1], pedges[1:]):
        x_, w_ = gl(nphi, a_, b_); pp.append(x_); wpp.append(w_)
    pp = np.concatenate(pp); wpp = np.concatenate(wpp)
    RP, PP = np.meshgrid(rp, pp, indexing='ij'); W2 = np.outer(wrp, wpp)
    rr = RP*np.sin(PP); zz = RP*np.cos(PP)
    eta = sign*M*np.sign(zz)/rr
    dz = z - zz
    D0 = r**2 + rr**2 + dz**2; c = 2*r*rr
    # exclusion ball
    mask = np.hypot(rr - r, zz - z) < eps*max(rho, rho0)
    out = np.zeros_like(RP)
    # chunk over theta to bound memory
    J = np.zeros_like(RP)
    for k in range(0, th.size, 32):
        t = th[k:k+32]; ws = s2[k:k+32]
        D = D0[..., None] - c[..., None]*np.cos(t)
        J += np.einsum('ijk,k->ij', D**(-2.5), ws)
    integrand = eta * rr**3 * dz * J * RP   # RP = Jacobian of (rho',phi')
    integrand[mask] = 0.0
    return (3.0/(2*np.pi)) * np.sum(integrand * W2)
def run():
    t0 = time.time(); res = {}
    rho0, R = 1.0, 4096.0
    exact_axis = 0.5*M*np.log(R/rho0)
    a00 = strain_at(0.0, 0.0, rho0, R)
    res['axis'] = dict(numeric=a00, exact=exact_axis, relerr=abs(a00-exact_axis)/exact_axis)
    print(f"CAL  a(0,0) numeric {a00:.5f}  exact (M/2)log(R/rho0) = {exact_axis:.5f}  relerr {res['axis']['relerr']:.2e}")
    ok = res['axis']['relerr'] < 5e-3
    # sign control
    a00s = strain_at(0.0, 0.0, rho0, R, sign=+1.0)
    print(f"CTRL sign-flipped field a(0,0) = {a00s:.5f} (must be -{exact_axis:.5f})"); ok &= abs(a00s + exact_axis)/exact_axis < 5e-3
    # single-octave control
    a1 = strain_at(0.0, 0.0, 1.0, 2.0); ex1 = 0.5*M*np.log(2.0)
    print(f"CTRL single octave a(0,0) = {a1:.5f} (exact {ex1:.5f})"); ok &= abs(a1-ex1)/ex1 < 5e-3
    # interior profile
    print("\nInterior strain a(x) for x = rho (sin phi, cos phi); residual c := a/M - (1/2) log(R/rho)")
    print(f"{'rho/rho0':>9} {'phi/pi':>7} {'a/M':>9} {'(1/2)log(R/rho)':>16} {'c':>8}")
    table = []
    for rho in [8.0, 32.0, 128.0, 512.0]:
        for phi in [np.pi/6, np.pi/3, 0.45*np.pi, 0.55*np.pi, 2*np.pi/3, 5*np.pi/6]:
            r, z = rho*np.sin(phi), rho*np.cos(phi)
            a = strain_at(r, z, rho0, R)
            c = a/M - 0.5*np.log(R/rho)
            table.append(dict(rho=rho, phi=phi/np.pi, a=a, c=c))
            print(f"{rho:9.1f} {phi/np.pi:7.3f} {a:9.4f} {0.5*np.log(R/rho):16.4f} {c:8.4f}")
    res['table'] = table
    # decision: is the log coefficient 1/2 at interior points?  fit a vs log(R/rho) per phi
    slopes = {}
    for phi in sorted(set(t['phi'] for t in table)):
        pts = [(np.log(R/t['rho']), t['a']) for t in table if t['phi']==phi]
        x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
        slope = np.polyfit(x, y, 1)[0]; slopes[phi] = slope
    print("\nFitted d a / d log(R/rho) per phi (claim: 1/2):", {f"{k:.3f}": round(v,4) for k,v in slopes.items()})
    res['slopes'] = slopes
    claim_holds = all(abs(v-0.5) < 0.06 for v in slopes.values())
    print("\nVERDICT (numerics):", "log coefficient 1/2 reproduced at interior points — the outer shells stretch the inner shells at the axis rate" if claim_holds else "claim NOT reproduced")
    res['claim_holds'] = bool(claim_holds); res['ok_controls'] = bool(ok); res['elapsed_s'] = time.time()-t0
    json.dump(res, open('g2_results.json','w'), indent=1, default=float)
    print('SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest(), f'elapsed {time.time()-t0:.1f}s')
    sys.exit(0 if ok else 2)
if __name__ == '__main__': run()

"""p7 - a genuinely independent instrument: raw 5-D Monte-Carlo of

        a[eta_0 o Lambda^-1](0) = INT_S Kcal(Lambda x) eta_0(x) J_Lambda(x) dx_5 ,
        a_ref                   = INT_S Kcal(Lambda x) eta_0(x) lambda(rho)^2 dx_5 ,
        I                       = INT_S |eta_0| lambda(rho)^2 |Lambda x|^{-4} dx_5   (Step 3)

with NO meridian reduction and NO Gauss-Legendre: points are drawn uniformly on S^4 (a normalised
5-D Gaussian) and uniformly in log rho.  This tests, in one shot, the 5-D measure
dx_5 = rho^4 drho dOmega_4, the SO(4) reduction of lib5, the closed form of det D(Lambda), and the
kernel -- all of which the p3/p4 route shares with the refuter's r2.

    INT_S f dx_5 = L * |S^4| * E[ f(x) rho^5 ],   |S^4| = 8 pi^2/3,  u = log(rho/rho_0) ~ U[0,L].
For these three integrands f*rho^5 is rho-free, so nothing overflows at large L.
"""
import json, numpy as np, lib5 as G

rng = np.random.default_rng(2026090601)
S4 = 8*np.pi**2/3
M = 1.0
N = 8_000_000
CHUNK = 1_000_000
out = {'N': N, 'rows': []}

for L in (8.317766166719343, 50.0):
    acc_a = acc_a2 = acc_r = acc_r2 = acc_I = acc_I2 = 0.0
    for _ in range(N//CHUNK):
        v = rng.standard_normal((CHUNK, 5))
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        cphi = v[:, 4]
        sphi = np.sqrt(np.maximum(1 - cphi**2, 0.0))
        u = rng.uniform(0.0, L, CHUNK)
        lam = G.lam_sigma(u/L)
        D = G.dloglam_dsig(u/L)/L
        g = np.sqrt(lam**2*sphi**2 + lam**-4*cphi**2)
        Kc = -(3/(8*np.pi**2))*lam**-2*cphi/g**5              # Kcal(Lambda x) * rho^4
        eta = -M*np.sign(cphi)/sphi                            # eta_0 * rho
        JL = lam**2*(1 + D*(sphi**2 - 2*cphi**2))
        fa = Kc*eta*JL
        fr = Kc*eta*lam**2
        fI = (M/sphi)*lam**2/g**4                              # |eta_0| lam^2 |Lambda x|^-4 * rho^5
        acc_a += fa.sum(); acc_a2 += (fa**2).sum()
        acc_r += fr.sum(); acc_r2 += (fr**2).sum()
        acc_I += fI.sum(); acc_I2 += (fI**2).sum()
    sc = L*S4
    def stat(s, s2):
        m = s/N; var = max(s2/N - m*m, 0.0)
        return sc*m, sc*np.sqrt(var/N)
    a_mc, a_se = stat(acc_a, acc_a2)
    r_mc, r_se = stat(acc_r, acc_r2)
    I_mc, I_se = stat(acc_I, acc_I2)
    lam_bar = float(np.sqrt(1.5))
    a_gl = G.a_of_map(L, [('lambda',)])
    r_ex = 0.5*L*lam_bar
    I_ex = np.pi**3*L*lam_bar
    out['rows'].append(dict(
        L=L,
        a_mc=a_mc, a_se=a_se, a_gl=a_gl, a_z=(a_mc - a_gl)/a_se,
        a_ref_mc=r_mc, a_ref_se=r_se, a_ref_exact=r_ex, a_ref_z=(r_mc - r_ex)/r_se,
        I_mc=I_mc, I_se=I_se, I_exact=I_ex, I_z=(I_mc - I_ex)/I_se))

out['max_abs_z'] = max(max(abs(r['a_z']), abs(r['a_ref_z']), abs(r['I_z'])) for r in out['rows'])
json.dump(out, open('p7_results.json', 'w'), indent=1)
print("5-D Monte-Carlo, N = %d per L, no meridian reduction" % N)
for r in out['rows']:
    print(f"L={r['L']:9.4f}")
    print(f"   a[eta0 o Lambda^-1](0) MC = {r['a_mc']:.6f} +- {r['a_se']:.6f}   "
          f"Gauss-Legendre = {r['a_gl']:.6f}   z = {r['a_z']:+.2f}")
    print(f"   a_ref                  MC = {r['a_ref_mc']:.6f} +- {r['a_ref_se']:.6f}   "
          f"exact (M/2)L*sqrt(3/2) = {r['a_ref_exact']:.6f}   z = {r['a_ref_z']:+.2f}")
    print(f"   Step-3 integral        MC = {r['I_mc']:.4f} +- {r['I_se']:.4f}   "
          f"exact pi^3 M A = {r['I_exact']:.4f}   z = {r['I_z']:+.2f}")
print("max |z| =", out['max_abs_z'])

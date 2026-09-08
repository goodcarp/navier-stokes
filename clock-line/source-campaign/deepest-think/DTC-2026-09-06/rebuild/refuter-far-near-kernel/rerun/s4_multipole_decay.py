#!/usr/bin/env python3
"""S4 - KILL K3': the exact multipole power with which ONE inner octave reaches the outer octaves.
Exterior expansion (S1, identity I2):  Phi(xi) = SUM_{l>=0} ((l+1)g_l/(2l+3)) |xi|^{-(l+4)} C_{l+1}^{3/2}(t).
An octave rho' in (rho_in, 2 rho_in) at level M_in therefore gives, at |x| = rho >> rho_in,
   a(x) = SUM_l ((l+1)g_l/(2l+3)) ((2^{l+4}-1)/(l+4)) (rho_in/rho)^{l+4} C_{l+1}^{3/2}(t).
Leading power  rho^{-4}  in general (l=0, g_0 propto shell mean of omega^theta),
               rho^{-5}  when omega^theta is odd in z (g_0 = 0)   <- the near-cap stack.
Both are checked against an INDEPENDENT direct 5D-lift quadrature (not the g2 code)."""
import numpy as np, math, json, hashlib, sys, time
sys.path.insert(0, '.')
from kern import geg_table, Nl

def g_of_w(wfun, LMAX, nquad=600):
    """g_l for a shell profile w(phi):  g_l N_l = INT_0^pi w(phi) C_l(cos phi) sin^2 phi dphi."""
    x, wq = np.polynomial.legendre.leggauss(nquad)
    th = 0.5*math.pi*(x+1.0); wq = 0.5*math.pi*wq
    C = geg_table(np.cos(th), LMAX)
    I = C @ (wq*np.sin(th)**2*wfun(th))
    return I/Nl(np.arange(LMAX+1))

def a_octave_series(rho, phi, rin, g, LMAX):
    t = math.cos(phi); C = geg_table(np.array([t]), LMAX+1)[:, 0]
    ls = np.arange(0, LMAX+1)
    return float(np.sum(((ls+1)*g[ls]/(2*ls+3))*((2.0**(ls+4)-1)/(ls+4))*(rin/rho)**(ls+4)*C[ls+1]))

# ---- independent direct quadrature of the 5D-lift Biot-Savart (own code, own panels) ----
def a_direct(r, z, rin, wfun, nr=140, nph=140, nth=140):
    """a = (3/2pi) INT INT omega^theta r'^2 (z-z') J drho' dphi' * rho',  J = INT_0^pi sin^2 th dth / |..|^{5/2}.
    Octave support rho' in (rin, 2 rin).  Evaluation point is far outside: plain Gauss-Legendre."""
    xg, wg   = np.polynomial.legendre.leggauss(nr);  rp = 0.5*rin*(xg+1)+rin;  wrp = 0.5*rin*wg
    xg2, wg2 = np.polynomial.legendre.leggauss(nph); pp = 0.5*math.pi*(xg2+1); wpp = 0.5*math.pi*wg2
    xg3, wg3 = np.polynomial.legendre.leggauss(nth); th = 0.5*math.pi*(xg3+1); wth = 0.5*math.pi*wg3  # th in (0,pi)
    RP, PP = np.meshgrid(rp, pp, indexing='ij'); W2 = np.outer(wrp, wpp)
    rr = RP*np.sin(PP); zz = RP*np.cos(PP); w = wfun(PP)
    dz = z - zz; D0 = r**2 + rr**2 + dz**2; c = 2*r*rr
    J = np.zeros_like(RP)
    for k in range(0, th.size, 25):
        tt = th[k:k+25]; ws = np.sin(tt)**2*wth[k:k+25]
        J += np.einsum('ijk,k->ij', (D0[..., None]-c[..., None]*np.cos(tt))**(-2.5), ws)
    return (3.0/(2*math.pi))*float(np.sum(w*rr**2*dz*J*RP*W2))

LMAX = 120; rin = 1.0
CASES = {
    "z-odd  (near-cap, omega^theta = -M sgn z)": (lambda ph: -np.sign(np.cos(ph)), 5),
    "not z-odd (omega^theta = -M constant)":     (lambda ph: -np.ones_like(ph),    4),
}
res = {}; k3p = False
for tag, (wf, pred) in CASES.items():
    g = g_of_w(wf, LMAX)
    print(f"\n=== {tag} ===")
    print(f"   g_0 = {g[0]:+.6f}   g_1 = {g[1]:+.6f}   (g_0 = 0 <=> z-odd)")
    phi = math.radians(40.0)
    print(f"{'rho/rho_in':>11} {'a_series':>15} {'a_direct':>15} {'rel diff':>10}")
    rows = []
    for j in range(2, 8):
        rho = 2.0**j
        As = a_octave_series(rho, phi, rin, g, LMAX)
        Ad = a_direct(rho*math.sin(phi), rho*math.cos(phi), rin, wf)
        rows.append((rho, As, Ad)); print(f"{rho:11.1f} {As:15.6e} {Ad:15.6e} {abs(As-Ad)/abs(Ad):10.2e}")
    x = np.log2([r[0] for r in rows]); ys = np.log2(np.abs([r[1] for r in rows])); yd = np.log2(np.abs([r[2] for r in rows]))
    ss = -np.polyfit(x[-4:], ys[-4:], 1)[0]; sd = -np.polyfit(x[-4:], yd[-4:], 1)[0]
    print(f"   fitted decay exponent (last 4 pts): series {ss:.4f}   direct {sd:.4f}   predicted {pred}")
    res[tag] = dict(g0=g[0], exp_series=ss, exp_direct=sd, pred=pred, rows=rows)
    if abs(ss-pred) > 0.15 or abs(sd-pred) > 0.15: k3p = True
print("\nKILL K3'", "FIRED" if k3p else "did not fire", "(threshold 0.15 in the exponent)")

# explicit constant for the near-cap case: |a| <= D * M_in * (rho_in/rho)^5, D from the l=1 term + tail
g = g_of_w(CASES["z-odd  (near-cap, omega^theta = -M sgn z)"][0], LMAX)
l = 1; D1 = abs((l+1)*g[l]/(2*l+3))*((2.0**(l+4)-1)/(l+4))*((l+2)*(l+3)/2)
print(f"\nleading-term constant for one z-odd octave:  |a| <= D * M_in * (rho_in/rho)^5,  D(l=1 term) = {D1:.4f}")
ls = np.arange(1, LMAX+1, 2)
Dall = float(np.sum(np.abs((ls+1)*g[ls]/(2*ls+3))*((2.0**(ls+4)-1)/(ls+4))*((ls+2)*(ls+3)/2)*0.5**(ls-1)))
print(f"   all-l absolute bound valid for rho >= 4 rho_in (|xi|>=2):  D = {Dall:.4f}")
json.dump(res, open('s4_results.json','w'), indent=1, default=float)
print('SCRIPT-SHA256', hashlib.sha256(open(__file__,'rb').read()).hexdigest())

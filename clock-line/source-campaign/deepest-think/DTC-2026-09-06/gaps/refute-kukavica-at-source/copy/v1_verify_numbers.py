"""
Independent re-derivation (decorrelation pass) of the kernel / weighted-BMO numbers behind
the Kukavica-vs-BFG kernel comparison.  Written and run by the gaps/kukavica-at-source
verification seat.  Deliberately different quadrature from the earlier s1_kernel_numbers.py:
Gauss-Legendre on the compactified radius r = s/(1-s) (plus a Cartesian Riemann sum for the
signed gradient integral), not Gauss-Hermite.  Every rule is reported at two resolutions so
the convergence is visible.

[V1] heat-kernel means in R^3 : int G dz, int d_1 G dz, int |grad G| dz  vs  2/sqrt(pi tau)
[V2] weighted-BMO test family  h_R(x) = (1 - log(1+|x|)/log R)_+  on R^3, weight 1/(|x|^4+1):
     unsubtracted J0(R) and locally-average-subtracted J1(R), each times log R
[V3] Duhamel scalings for the two kernels' pointwise bounds (grad-kernel -> sqrt(t), G -> t)
[V4] equivalence of the two weights used, (|x|+1)^4 vs |x|^4+1
"""
import json, math
import numpy as np

_CACHE = {}
def gl(n):
    if n not in _CACHE:
        _CACHE[n] = np.polynomial.legendre.leggauss(n)
    return _CACHE[n]

def radial_nodes(n):
    """nodes/weights for int_0^inf f(r) dr via r = s/(1-s)."""
    x, w = gl(n)
    s = 0.5*(x+1.0); w = 0.5*w
    r = s/(1.0-s)
    return r, w/(1.0-s)**2

def sphere3(f, n):
    r, w = radial_nodes(n)
    return 4.0*math.pi*float(np.sum(w*f(r)*r*r))

def seg(f, a, b, n):
    x, w = gl(n)
    r = 0.5*(b-a)*x + 0.5*(b+a); w = 0.5*(b-a)*w
    return float(np.sum(w*f(r)))

out = {}

# ---------------- [V1] ----------------
v1 = []
for tau in (1.0, 1e-2, 1e-4):
    G      = lambda r: (4*math.pi*tau)**-1.5*np.exp(-r*r/(4*tau))
    absgrd = lambda r: (4*math.pi*tau)**-1.5*np.exp(-r*r/(4*tau))*(r/(2*tau))
    rec = dict(tau=tau)
    for n in (800, 1600):
        rec[f'int_G_n{n}']         = sphere3(G, n)
        rec[f'int_absgradG_n{n}']  = sphere3(absgrd, n)
    rec['two_over_sqrt_pi_tau'] = 2.0/math.sqrt(math.pi*tau)
    # signed int of d_1 G on a symmetric Cartesian grid (must vanish by oddness)
    m = 201; L = 12*math.sqrt(tau)
    x = np.linspace(-L, L, m); h = x[1]-x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    d1G = (4*math.pi*tau)**-1.5*np.exp(-(X*X+Y*Y+Z*Z)/(4*tau))*(-X/(2*tau))
    rec['signed_int_d1G_cartesian'] = float(d1G.sum()*h**3)
    v1.append(rec)
out['V1_heat_kernel_means'] = v1

# ---------------- [V2] ----------------
def hR(r, R):
    return np.maximum(0.0, 1.0 - np.log1p(r)/math.log(R))

v2 = []
for R in (1e2, 1e4, 1e8, 1e16, 1e32, 1e64):
    lg = math.log(R)
    c  = 4.0*math.pi*seg(lambda r: hR(r, R)*r*r, 0.0, 1.0, 800)/(4.0*math.pi/3.0)
    rec = dict(R=R, log_R=lg, avg_B1=c)
    for n in (1200, 2400):
        J0 = sphere3(lambda r: hR(r, R)/(r**4+1.0), n)
        J1 = sphere3(lambda r: np.abs(hR(r, R)-c)/(r**4+1.0), n)
        rec[f'J0_n{n}'] = J0; rec[f'J0_logR_n{n}'] = J0*lg
        rec[f'J1_n{n}'] = J1; rec[f'J1_logR_n{n}'] = J1*lg
    v2.append(rec)
out['V2_weighted_bmo_family'] = v2

# ---------------- [V3] ----------------
# Kukavica's bound  |grad G(z,tau)| <= C/(|z|+sqrt(tau))^{n+1}  -> int_0^t int (|z|+sqrt(t-s))^-4
# BFG's bound       G(z,tau)       <= c sqrt(tau)/(|z|+sqrt(tau))^4 -> extra sqrt(t-s) factor
rg, wg = radial_nodes(2400)
v3 = []
for t in (0.25, 1.0, 4.0):
    xs, ws_ = gl(600)
    s = 0.5*t*(xs+1.0); wq = 0.5*t*ws_
    a = np.sqrt(t-s)
    inner = 4.0*math.pi*np.sum((wg*rg*rg)[None, :]/(rg[None, :]+a[:, None])**4, axis=1)
    A = float(np.sum(wq*inner))        # grad-kernel Duhamel
    B = float(np.sum(wq*a*inner))      # G-kernel Duhamel
    v3.append(dict(t=t, gradG_integral=A, gradG_over_sqrt_t=A/math.sqrt(t),
                   G_integral=B, G_over_t=B/t,
                   closed_form_8pi_over_3=8*math.pi/3,
                   closed_form_4pi_over_3=4*math.pi/3))
out['V3_duhamel_scalings'] = v3

# ---------------- [V4] ----------------
rr = np.linspace(0.0, 60.0, 600_001)
ratio = (rr+1.0)**4/(rr**4+1.0)
out['V4_weight_equivalence'] = dict(sup_ratio=float(ratio.max()),
                                    argmax_r=float(rr[int(ratio.argmax())]),
                                    inf_ratio=float(ratio.min()))
print(json.dumps(out, indent=1))

#!/usr/bin/env python3
"""s1_exact_lemma.py -- THE STRAIN-TRANSPORT LEMMA, three independent ways.

CLAIM (Lemma 1).  Let  eta_0(rho,phi) = -M h(phi) sgn(cos phi) / (rho sin phi)  on the shell
rho_0 < rho < R  (so omega^theta = r eta = -M h(phi) sgn(cos phi), |omega^theta| <= M),
and let T_lam be the volume-preserving axisymmetric strain map on R^3

        T_lam : (r, theta, z) |-> (lam r, theta, lam^{-2} z),

with eta transported as a scalar (eta_lam(T_lam y) = eta_0(y)), i.e. omega^theta |-> lam omega^theta.
Then the axial strain at the ORIGIN is EXACTLY

        a_lam(0) = (M/2) P_h(lam) log(R/rho_0),
        P_h(lam) = 3 int_0^1 h(v) v^2 (A v^2 + B)^{-5/2} dv,  A = lam^2 - lam^{-4}, B = lam^{-4},

with v = sin phi.  For the untapered plateau h == 1 this is  P_1(lam) = lam  EXACTLY.

The rho and phi integrals FACTORISE exactly (no endpoint error) because the strained shell,
read in Lagrangian coordinates, still has the scale-invariant radial profile 1/rho.

Route 1: sympy/mpmath closed form of the v-integral.
Route 2: direct numerical quadrature of the Lagrangian integral.
Route 3: INDEPENDENT -- build eta_lam on an Eulerian (rho,phi) grid by inverting T_lam
         pointwise and integrate  a(0) = (3/4) int int (-cos phi sin^3 phi) eta drho dphi.
"""
import numpy as np, sympy as sp, mpmath as mp, json

OUT = {}

# ---------------------------------------------------------------- Route 1: symbolic
# Verify the antiderivative by DIFFERENTIATION (fast and airtight), not by sp.integrate.
v, lam = sp.symbols('v lam', positive=True)
A = lam**2 - lam**-4
B = lam**-4
# claimed antiderivative:  Ant(v) = v^3 / ( 3 B (A v^2 + B)^{3/2} )
Ant = v**3/(3*B*(A*v**2 + B)**sp.Rational(3,2))
resid = sp.simplify(sp.diff(Ant, v) - v**2/(A*v**2 + B)**sp.Rational(5,2))
OUT['antiderivative_residual'] = str(resid)
I_closed = sp.simplify(Ant.subs(v, 1) - Ant.subs(v, 0))     # = 1/(3 B (A+B)^{3/2}) = lam/3
OUT['I_closed'] = str(sp.simplify(I_closed))
OUT['I_closed_minus_lam_over_3'] = str(sp.simplify(I_closed - lam/3))
OUT['P1_equals_lam'] = str(sp.simplify(3*I_closed - lam))
# high-precision numeric confirmation
mp.mp.dps = 30
worst = 0.0
for L in ['1.0','1.05','1.2','1.5','2.0','3.0','5.0']:
    Lm = mp.mpf(L); Am = Lm**2 - Lm**-4; Bm = Lm**-4
    num = mp.quad(lambda x: x**2/(Am*x**2+Bm)**mp.mpf(2.5), [0,1])
    worst = max(worst, abs(3*num - Lm)/Lm)
OUT['mpmath_P1_minus_lam_relmax'] = float(worst)

# ---------------------------------------------------------------- profiles
def h_flat(phi):    return np.ones_like(phi)
def h_sin2(phi):    return np.abs(np.sin(2*phi))          # viscous-numerics datum A
def h_taper(delta):
    def h(phi):
        pax = np.minimum(phi, np.pi - phi)                 # angle to the nearest axis direction
        return np.minimum(1.0, pax/delta)
    return h

def P_h(h, lam, n=200001):
    """P_h(lam) = 3 int_0^1 h(v) v^2 (A v^2 + B)^{-5/2} dv, v = sin phi, phi in (0,pi/2)."""
    A_ = lam**2 - lam**-4
    B_ = lam**-4
    phi = np.linspace(0.0, np.pi/2, n)
    vv = np.sin(phi)
    integ = h(phi)*vv**2*(A_*vv**2 + B_)**-2.5*np.cos(phi)   # dv = cos phi dphi
    return 3.0*np.trapz(integ, phi)

# kappa = (1/2) P_h(1)
OUT['kappa'] = {
    'flat':      0.5*P_h(h_flat, 1.0),
    'sin2phi':   0.5*P_h(h_sin2, 1.0),
    'taper_30d': 0.5*P_h(h_taper(np.deg2rad(30.0)), 1.0),
    'taper_15d': 0.5*P_h(h_taper(np.deg2rad(15.0)), 1.0),
    'taper_7.5d':0.5*P_h(h_taper(np.deg2rad(7.5)), 1.0),
}

# ---------------------------------------------------------------- Route 2 vs Route 1
lams = [1.0, 1.1, 1.25, 1.5, 1.75, 2.0, 3.0]
OUT['route2_P1_vs_lam'] = [(L, P_h(h_flat, L), abs(P_h(h_flat, L)-L)/L) for L in lams]

# ---------------------------------------------------------------- Route 3: Eulerian, independent
def a_origin_eulerian(h, lam, rho0=1.0, R=None, nr=4000, nphi=4000, Loct=6.0):
    """a(0) for the field  eta_lam = eta_0 o T_lam^{-1}  built pointwise on an EULERIAN grid.
    Returns a(0)/(M L) so it should equal P_h(lam)/2 with L = log(R/rho0)."""
    if R is None: R = rho0*np.exp(Loct)
    # image of the shell lives in rho in [rho0*lam^-2, R*lam] for lam>=1
    rmin, rmax = rho0*min(lam, lam**-2)*0.9, R*max(lam, lam**-2)*1.1
    s = np.linspace(np.log(rmin), np.log(rmax), nr)
    rho = np.exp(s)
    phi = np.linspace(1e-9, np.pi/2-1e-9, nphi)       # z>0 half; double at the end
    RHO, PHI = np.meshgrid(rho, phi, indexing='ij')
    r = RHO*np.sin(PHI); z = RHO*np.cos(PHI)
    # invert T_lam
    rp = r/lam; zp = z*lam**2
    rhop = np.hypot(rp, zp); phip = np.arctan2(rp, zp)
    inside = (rhop > rho0) & (rhop < R)
    eta = np.where(inside, -h(phip)/(rhop*np.sin(phip)), 0.0)   # M = 1
    integ = (-np.cos(PHI)*np.sin(PHI)**3)*eta                    # drho dphi measure
    # integrate: drho = rho ds
    inner = np.trapz(integ, phi, axis=1)
    tot = np.trapz(inner*rho, s)
    return 2.0*0.75*tot/np.log(R/rho0)      # x2 for the z<0 half (odd eta, even contribution)

r3 = []
for L in [1.0, 1.25, 1.5, 2.0]:
    val = a_origin_eulerian(h_flat, L)
    r3.append((L, val, 0.5*P_h(h_flat, L), abs(val - 0.5*L)/(0.5*L)))
OUT['route3_eulerian_flat'] = r3

r3t = []
dl = np.deg2rad(7.5)
for L in [1.0, 1.25, 1.5]:
    val = a_origin_eulerian(h_taper(dl), L)
    r3t.append((L, val, 0.5*P_h(h_taper(dl), L), abs(val - 0.5*P_h(h_taper(dl), L))/(0.5*P_h(h_taper(dl), L))))
OUT['route3_eulerian_taper7.5'] = r3t

# ---------------------------------------------------------------- P'_h(1), the acceleration coefficient
def Pprime1(h, n=400001):
    phi = np.linspace(0.0, np.pi/2, n)
    vv = np.sin(phi)
    num = np.trapz(h(phi)*vv**2*(-5.0)*(3*vv**2-2)*np.cos(phi), phi)
    den = np.trapz(h(phi)*vv**2*np.cos(phi), phi)
    return 3.0*num, num/den          # (dP/dlam at 1, normalised Phi'(1))
OUT['Pprime1'] = {
    'flat':    Pprime1(h_flat),
    'sin2phi': Pprime1(h_sin2),
    'taper7.5':Pprime1(h_taper(dl)),
}

print(json.dumps(OUT, indent=1, default=str))
with open('s1_results.json','w') as f: json.dump(OUT, f, indent=1, default=str)

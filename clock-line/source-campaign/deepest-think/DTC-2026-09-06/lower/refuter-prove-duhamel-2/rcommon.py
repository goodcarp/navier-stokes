"""shared: the delta-tapered mollified plateau datum, on the estate's validated solver.

The solver `nsring.py` is IMPORTED read-only from sharp/viscous-numerics (its sha256 is
recorded by every script that uses it).  Nothing outside this folder is written.

Datum (admissible: omega vanishes on the axis, eta = omega^theta/r is bounded and smooth):

    s = |x|,  sin(phi) = r/s,  cos(phi) = z/s
    omega^theta = -M * tanh( sin(phi)/sin(delta) ) * tanh( cos(phi)/w ) * Theta(s)
    eta = omega^theta / r = -(M/(s sin delta)) * T(x)/x |_{x=sin phi/sin delta}
                              * tanh(cos phi / w) * Theta(s)          (regular on the axis)
    Theta(s) = (1/2)[ tanh((s-rho0)/w0) - tanh((s-R)/w1) ],  w0 = .25 rho0, w1 = .10 R

delta -> 0, w -> 0 recovers the bang-bang plateau omega^theta = -M sgn(z) on rho0<|x|<R.
"""
import os, sys, math, hashlib
import numpy as np

VNUM = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "sharp", "viscous-numerics"))
sys.path.insert(0, VNUM)
import nsring
from nsring import (Grid, Poisson, cell_to_node, face_fluxes, advect_rhs, lap5,
                    node_velocity, energy, interp2, hill_eta, hill_psi1_exact)

def solver_sha():
    return hashlib.sha256(open(os.path.join(VNUM, "nsring.py"), "rb").read()).hexdigest()

def datum_taper(g, N, rho0=1.0, delta_deg=15.0, w=0.20, sub=3):
    """returns (eta on cell centres, R).  eta is odd in z."""
    R = rho0 * 2.0**N
    sd = math.sin(math.radians(delta_deg))
    w0 = 0.25 * rho0; w1 = 0.10 * R
    RR, ZZ = g.cell_mesh(); h = g.h
    off = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros_like(RR)
    for dr in off:
        for dz in off:
            r = RR + dr * h; z = ZZ + dz * h
            s = np.sqrt(r * r + z * z)
            s = np.maximum(s, 1e-300)
            sp = r / s; cp = z / s
            xx = np.maximum(sp / sd, 1e-12)
            tox = np.tanh(xx) / xx                      # -> 1 as xx -> 0
            Th = 0.5 * (np.tanh((s - rho0) / w0) - np.tanh((s - R) / w1))
            acc += -(1.0 / (s * sd)) * tox * np.tanh(cp / w) * Th
    return acc / sub**2, R

def build_taper(N, h, lam=1.5, rho0=1.0, delta_deg=15.0, w=0.20, sign=+1.0, zmin=0.0):
    R = rho0 * 2.0**N
    L = max(lam * R, 4.0 * rho0)
    g = Grid(L, zmin, L, h)
    eta, R = datum_taper(g, N, rho0, delta_deg, w)
    om = eta * g.Rc[:, None]
    eta = eta * sign / np.max(np.abs(om))               # so max|omega^theta| = 1 exactly
    return g, eta, R, L

def make_solve(g):
    P = Poisson(g)
    zr = np.zeros(g.Nr + 1); zj = np.zeros(g.Nz + 1)
    return lambda e: P.solve(cell_to_node(g, e, "odd"), g_bot=zr, g_top=zr, g_right=zj)

def node_field(g, cellF):
    """cell field -> node field, even in r, odd in z (for eta)."""
    return cell_to_node(g, cellF, "odd")

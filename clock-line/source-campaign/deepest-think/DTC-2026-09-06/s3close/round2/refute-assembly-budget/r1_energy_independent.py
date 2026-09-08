"""
r1 -- decorrelated checks of the energy conversion C_E = ||u_0||_2^2/(M^2 R^5).

Three routes that share nothing with a1_datum_energy.py's Gegenbauer / 5D-Green route:

 (A) the 3D identity  ||u||_2^2 = int A . omega dx,  A = (1/4pi) int omega(y)/|x-y| dy
     (u = curl A, div A = 0), specialised to a purely azimuthal omega = w(r,z) e_theta:
        ||u||_2^2 = (1/2) intint w(r,z) w(r',z') r r' G(r,z;r',z') dr dz dr' dz'
        G = int_0^{2pi} cos(t) dt / sqrt(r^2 + r'^2 - 2 r r' cos t + (z-z')^2)
     G is evaluated by its complete-elliptic-integral closed form, which is itself
     checked here against brute-force quadrature.  The 4D integral is done on a
     tensor Gauss grid in (log rho, phi) x (log rho', phi') at finite L, with the
     log-singular diagonal handled by a Duffy-free but symmetric grid; accuracy is
     reported by grid refinement.  Applied to:
        - the bang-bang cap (record 0.172403978 at L -> inf),
        - the runs' datum omega = -M sin 2phi (pure l = 1),
        - the assembly's design datum (delta = 7.5 deg, delta_m = 5 deg),
     each at finite L and compared with the SAME finite-L modal value computed here
     from the assembly's formula (D_l at finite L) -- not from a1's JSON.

 (B) the runs' own Re_E (sharp/viscous-numerics NOTE section 4, corner-numerics
     section 2): C_E = (Re_E nu)^{5/2} / R^5 with nu = 1/Re0, M = rho0 = 1.  If the
     runs used Re_E = E^{2/5} M^{1/5}/nu with E = int |u|^2 (the assembly's convention),
     this must land near the modal value for -sin 2phi (plus the runs' radial smoothing).

 (C) my own Gegenbauer coefficients by an independent quadrature in phi (not t),
     to confirm a1's H_l for the tapered datum.

Every number here is computed; the record values are quoted only as comparison targets
and labelled as such.
"""
import json, math
import numpy as np
from scipy.special import ellipk, ellipe, eval_gegenbauer
from scipy.integrate import quad

deg = math.pi/180.0
OUT = {}

# ---------------------------------------------------------------- profiles
def gprof(phi, delta, dm):
    """|angular profile| in (0, pi); sign = sgn(cos phi) applied outside"""
    phi_ax = min(phi, math.pi-phi)
    h = 1.0 if delta <= 0 else min(1.0, phi_ax/delta)
    m = 1.0 if dm <= 0 else min(1.0, abs(phi-math.pi/2)/dm)
    return h*m

def w_of(kind, delta=0.0, dm=0.0):
    if kind == "sin2phi":
        return lambda phi: -np.sin(2*phi)
    def w(phi):
        phi = np.asarray(phi, dtype=float)
        s = np.sign(np.cos(phi))
        g = np.array([gprof(p, delta, dm) for p in phi.ravel()]).reshape(phi.shape)
        return -s*g
    return w

# ---------------------------------------------------------------- ring kernel
def G_ring(r, rp, dz):
    """int_0^{2pi} cos t / sqrt(r^2+rp^2-2 r rp cos t + dz^2) dt  (vectorised)"""
    a = (r+rp)**2 + dz**2
    k2 = 4.0*r*rp/a
    k2 = np.clip(k2, 0.0, 1.0-1e-15)
    K = ellipk(k2); E = ellipe(k2)
    return (4.0/np.sqrt(a))*(((2.0-k2)*K - 2.0*E)/k2)

def G_brute(r, rp, dz, n=200000):
    t = np.linspace(0, 2*math.pi, n, endpoint=False)
    return float(np.mean(np.cos(t)/np.sqrt(r*r+rp*rp-2*r*rp*np.cos(t)+dz*dz))*2*math.pi)

chk = {}
for (r, rp, dz) in [(1.0, 2.0, 0.5), (1.0, 1.3, 0.1), (0.3, 2.5, 3.0), (1.0, 1.05, 0.02)]:
    chk["G_closed_vs_brute_r=%g_rp=%g_dz=%g" % (r, rp, dz)] = [float(G_ring(r, rp, dz)), G_brute(r, rp, dz)]
OUT["kernel_checks"] = chk

# ---------------------------------------------------------------- route (A): 4D energy
def energy_ring(w, L, n_rho=48, n_phi=48, M=1.0, R=1.0):
    """E = int |u|^2 dx  for omega = w(phi) e_theta on e^{-L} R < rho < R  (sharp radial cutoff).
    Coordinates x = log(rho/R) in [-L, 0], phi in (0, pi).  dx3 = 2 pi r dr dz = 2 pi rho^2 sin phi drho dphi
    but we use the 4D reduced form E = (1/2) intint w w' r r' G dr dz dr' dz' with
    dr dz = rho drho dphi = rho^2 dx dphi."""
    xg, xw = np.polynomial.legendre.leggauss(n_rho)
    pg, pw = np.polynomial.legendre.leggauss(n_phi)
    # split phi at pi/2 (profile kink / sign change) and x panels
    x_nodes = np.concatenate([(-L + L*(xg+1)/2)])
    x_wts = xw*L/2
    ph_nodes = np.concatenate([math.pi/4*(pg+1), math.pi/2 + math.pi/4*(pg+1)])
    ph_wts = np.concatenate([pw*math.pi/4, pw*math.pi/4])
    X, P = np.meshgrid(x_nodes, ph_nodes, indexing="ij")
    WX, WP = np.meshgrid(x_wts, ph_wts, indexing="ij")
    rho = R*np.exp(X)
    r = rho*np.sin(P); z = rho*np.cos(P)
    wv = M*w(P)
    # measure for dr dz : rho^2 dx dphi ; integrand carries r r' w w' G
    f = (wv*r*rho**2*WX*WP).ravel()
    rr = r.ravel(); zz = z.ravel()
    tot = 0.0
    nn = len(rr)
    # symmetric tensor product; the diagonal (log-singular) points are dropped: their
    # measure is O(h^2 log h) and refinement reports the effect.
    for i in range(nn):
        Gi = G_ring(rr[i], rr, zz[i]-zz)
        Gi[i] = 0.0
        tot += f[i]*np.dot(f, Gi)
    return 0.5*tot

# finite-L modal formula (the assembly's formula, my own implementation)
def N_l(l): return (l+1.0)*(l+2.0)/(l+1.5)
def H_l_phi(l, w, M=1.0):
    """H_l = (1/N_l) int_{-1}^{1} W(t) C_l(t) (1-t^2) dt,  W(t) = w(phi)/sin(phi) with t = cos phi
       = (1/N_l) int_0^pi w(phi) C_l(cos phi) sin^2(phi) dphi   (phi-quadrature, not t)"""
    f = lambda p: w(np.array([p]))[0]*eval_gegenbauer(l, 1.5, math.cos(p))*math.sin(p)**2
    pts = sorted({math.pi/2} | {7.5*deg, math.pi-7.5*deg, math.pi/2-5*deg, math.pi/2+5*deg})
    val, _ = quad(f, 0.0, math.pi, points=pts, limit=400, epsabs=1e-13, epsrel=1e-13)
    return M*val/N_l(l)
def D_l_finite(l, L):
    if l == 1:
        return (1.0/5.0)*((1.0-math.exp(-5*L))/5.0 - L*math.exp(-5*L))
    return (1.0/(l+4.0))*((1.0-math.exp(-5*L))/5.0 - (math.exp(-5*L)-math.exp(-(l+4.0)*L))/(l-1.0))
def C_E_modal(w, L, lmax=81, Linf=False):
    tot = 0.0
    for l in range(1, lmax+1, 2):
        H = H_l_phi(l, w)
        D = 1.0/(5.0*(l+4.0)) if Linf else D_l_finite(l, L)
        tot += N_l(l)*H*H*2.0*D/(2*l+3.0)
    return 2.0*math.pi*tot

res = {}
Lfin = 3.0*math.log(2.0)          # R = 8 rho0, the N = 3 run's shell
for name, w in [("bangbang", w_of("plateau")), ("sin2phi", w_of("sin2phi")),
                ("design_d7.5_dm5", w_of("plateau", 7.5*deg, 5*deg))]:
    row = {}
    Rf = math.exp(Lfin)   # in units rho0 = 1, R = e^L
    for n in [24, 36, 48]:
        E = energy_ring(w, Lfin, n_rho=n, n_phi=n, R=Rf)
        row["ring_CE_n=%d" % n] = E/Rf**5
    row["modal_CE_finiteL"] = C_E_modal(w, Lfin)
    row["modal_CE_Linf"] = C_E_modal(w, Lfin, Linf=True)
    res[name] = row
    print(name, row, flush=True)
OUT["route_A_ring_vs_modal_at_L=3ln2"] = res
# H_l check against a1 (read-only) for the design datum
import os
a1 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assembly", "a1_results.json")))
OUT["a1_CE_design_Linf_readonly"] = a1["C_E_table"]["delta=7.5,dm=5,eps=0"]
OUT["a1_CE_bangbang_Linf_readonly"] = a1["C_E_sharp_cap_Linf"]
OUT["my_modal_vs_a1_design_rel"] = abs(res["design_d7.5_dm5"]["modal_CE_Linf"]-OUT["a1_CE_design_Linf_readonly"])/OUT["a1_CE_design_Linf_readonly"]
OUT["record_bangbang_target_quoted"] = 0.172403978
OUT["sin2phi_exact_l1_only"] = 2*math.pi*N_l(1)*(2.0/3.0)**2*2.0/(5.0*5.0)/5.0

# ---------------------------------------------------------------- route (B): runs' Re_E
runs = [  # (N, R, Re_E, Re0) transcribed from sharp/viscous-numerics NOTE section 4 (Re0=100)
    (2, 4, 640.1, 100), (3, 8, 2572.0, 100), (4, 16, 1.029e4, 100), (5, 32, 4.116e4, 100),
    (6, 64, 1.646e5, 100), (7, 128, 6.585e5, 100),
    # corner-numerics section 2 (log Re_E given): (N, R, Re_E, Re0)
    (6, 64, math.exp(7.4063), 1), (6, 64, math.exp(8.7926), 4), (6, 64, math.exp(10.1789), 16),
    (6, 64, math.exp(12.0114), 100), (6, 64, math.exp(13.3977), 400), (7, 128, math.exp(13.3977), 100),
    (7, 128, math.exp(14.7840), 400), (5, 32, math.exp(10.6252), 100), (7, 128, math.exp(10.6252), 6.25)]
rb = []
for N, R, ReE, Re0 in runs:
    nu = 1.0/Re0
    E = (ReE*nu)**2.5            # E^{2/5} M^{1/5}/nu = Re_E with M = 1
    rb.append({"N": N, "Re0": Re0, "C_E_implied_full_energy": E/R**5,
               "C_E_implied_if_half_energy": 2*E/R**5})
OUT["route_B_runs_implied_CE"] = rb
vals = np.array([r["C_E_implied_full_energy"] for r in rb])
OUT["route_B_summary"] = {"mean": float(vals.mean()), "min": float(vals.min()), "max": float(vals.max()),
                          "modal_sin2phi_Linf": OUT["sin2phi_exact_l1_only"],
                          "synthesis_quoted_0.1418961_is_for": "prove-duhamel's own mollified datum d2, not the runs' datum (grep in lower/)"}

# ---------------------------------------------------------------- dictionary constants
s_sin = 1.0/math.sin(7.5*deg); s_lin = 1.0/(7.5*deg)
CEd = res["design_d7.5_dm5"]["modal_CE_Linf"]
OUT["dictionary"] = {
    "s=1/sin(7.5deg)": s_sin, "s=1/(7.5deg in rad)": s_lin,
    "2log s + 0.4 log C_E  (s=1/sin)": 2*math.log(s_sin) + 0.4*math.log(CEd),
    "2log s + 0.4 log C_E  (s=1/delta)": 2*math.log(s_lin) + 0.4*math.log(CEd),
    "difference": 2*math.log(s_sin/s_lin),
    "assembly_quoted_3.3643455": 3.3643455,
    "half_energy_shift": 0.4*math.log(0.5),
    "c2=4log(3/2)": 4*math.log(1.5),
}

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r1_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True)
print(json.dumps(OUT, indent=1, sort_keys=True))

"""
f6 -- the loose ends: controls for f2's shell-density instrument, the three numeric branch
      bounds quoted in FIX2 section 4(c), and L_Gamma^exist at the window L_* actually sits at.

Outputs -> f6_results.json
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "imported"))
import numpy as np                                     # noqa: E402
from f2_datum_s3 import DatumS3                        # noqa: E402
import f5_gamma_exist as F5                            # noqa: E402

DEG = math.pi/180.0
DELTA, DM = 7.5*DEG, 5.0*DEG
SD, CD = math.sin(DELTA), math.cos(DELTA)
S3SPH = 2*math.pi**2
OUT = {}

# ---------------------------------------------------------------- 1. the J instrument, controlled
# hk2 Lemma 5.1:  J_ac = |S^3| int_0^pi sqrt(W^2 + W_phi^2) sin^3 phi dphi , with W = Hh.
def J_ac(Wfun, Wpfun, n=400001):
    phi = np.linspace(1e-9, math.pi-1e-9, n)
    return S3SPH*float(np.trapz(np.sqrt(Wfun(phi)**2 + Wpfun(phi)**2)*np.sin(phi)**3, phi))

# (a) bare plateau, W = 1/sin phi : closed form 2 pi^2 (sqrt2 + arcsinh 1)
Jp = J_ac(lambda p: 1.0/np.sin(p), lambda p: -np.cos(p)/np.sin(p)**2)
Jp_closed = 2*math.pi**2*(math.sqrt(2) + math.asinh(1))
# (b) hk2's (D-A) angular profile, a.c. part only (taper min(1,phi_ax/delta), sharp sgn(z))
def W_DA(phi):
    pax = np.minimum(phi, math.pi - phi)
    return np.minimum(1.0, pax/DELTA)/np.sin(phi)
def Wp_DA(phi):
    pax = np.minimum(phi, math.pi - phi)
    h = np.minimum(1.0, pax/DELTA)
    hp = np.where(pax < DELTA, np.where(phi < math.pi/2, 1.0, -1.0)/DELTA, 0.0)
    s, c = np.sin(phi), np.cos(phi)
    return (hp*s - h*c)/s**2
J_DA = J_ac(W_DA, Wp_DA)
# (c) hk2's (D-B) angular profile, tanh(sin phi/sin delta) tanh(cos phi/0.20)
def W_DB(phi):
    return np.tanh(np.sin(phi)/SD)*np.tanh(np.abs(np.cos(phi))/0.20)/np.sin(phi)
def Wp_DB(phi, h=1e-6):
    return (W_DB(phi+h) - W_DB(phi-h))/(2*h)
J_DB = J_ac(W_DB, Wp_DB)
OUT["J_instrument_controls"] = {
    "bare_plateau_quadrature": Jp, "bare_plateau_closed_form": Jp_closed,
    "bare_plateau_rel": abs(Jp - Jp_closed)/Jp_closed,
    "DA_taper_only_ac_here": J_DA, "DA_taper_only_ac_hk2": 39.162729238636764,
    "DA_rel": abs(J_DA - 39.162729238636764)/39.162729238636764,
    "DB_here": J_DB, "DB_hk2": 65.62590994919597,
    "DB_rel": abs(J_DB - 65.62590994919597)/65.62590994919597}

# ---------------------------------------------------------------- 2. FIX2 sec.4(c) branch bounds
YMAX = (4.5*4.5)**2/64.0
def branch(P, Q):
    return (P + Q) if 81*P <= 8*Q else (YMAX*P + Q)

# bulk phi in (delta, pi/2 - delta_m): P = 1/s^2, Q = c^2/s^4
phis = np.linspace(DELTA, math.pi/2 - DM, 400001)
P = 1.0/np.sin(phis)**2
Q = np.cos(phis)**2/np.sin(phis)**4
crossover = phis[np.argmax(81*P <= 8*Q)] if np.any(81*P <= 8*Q) else None
mask = 81*P > 8*Q
bulk_other = float((YMAX*P[mask] + Q[mask]).max()) if mask.any() else 0.0
bulk_PQ = float((P[~mask] + Q[~mask]).max()) if (~mask).any() else 0.0
# taper phi in (0, delta): P = (phi/(delta sin phi))^2 , Q = ((s - phi c)/(delta s^2))^2
pt = np.linspace(1e-7, DELTA, 400001)
s, c = np.sin(pt), np.cos(pt)
Pt = (pt/(DELTA*s))**2
Qt = ((s - pt*c)/(DELTA*s**2))**2
taper = float(np.max([branch(a, b) for a, b in zip(Pt[::400], Qt[::400])]))
# equator phi in (pi/2 - delta_m, pi/2): P = ((pi/2-phi)/(delta_m sin phi))^2
pe = np.linspace(math.pi/2 - DM, math.pi/2 - 1e-12, 200001)
se, ce = np.sin(pe), np.cos(pe)
A = (math.pi/2 - pe)/DM
Ap = -1.0/DM
Pe = (A/se)**2
Qe = ((Ap*se - A*ce)/se**2)**2
equator = float(np.max([branch(a, b) for a, b in zip(Pe[::200], Qe[::200])]))
OUT["branch_bounds_sec4c"] = {
    "one_over_sin4_delta": 1.0/SD**4,
    "bulk_81P<=8Q_region_max_P+Q": bulk_PQ,
    "bulk_81P>8Q_region_max_6.407P+Q": bulk_other,
    "bulk_crossover_phi_deg": (float(crossover*180/math.pi) if crossover is not None else None),
    "taper_max_branch": taper, "equator_max_branch": equator,
    "taper_Fprime_at_delta_minus": (math.sin(DELTA) - DELTA*math.cos(DELTA))/(DELTA*SD**2),
    "equator_Fprime_at_pi/2": 1.0/DM,
    "all_below_1_over_sin4_delta": bool(max(bulk_PQ, bulk_other, taper, equator) <= 1.0/SD**4)}

# ---------------------------------------------------------------- 3. L_Gamma^exist at the L_* window
KAPPA = 0.49782243818691574
CSTAR = math.log(1.5)/KAPPA
E0, G0 = 1.0/SD, 1.0/SD**2
OUT["L_Gamma_exist_at_L_star_window"] = {}
for x in (1.1360549923591519, 1.0999999955536253):
    c = x*CSTAR
    OUT["L_Gamma_exist_at_L_star_window"]["c=%.7f c_*" % x] = {
        "c": c, "lam_max": math.exp(0.75*c),
        "L_exist": F5.L_exist(c, 1.5, E0, G0, 0.0),
        "L_picard": F5.L_picard(c, 1.5, E0, G0, 0.0)}

with open(os.path.join(HERE, "f6_results.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print(json.dumps(OUT, indent=1, sort_keys=True, default=str))

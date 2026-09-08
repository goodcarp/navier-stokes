"""
t1 -- the datum of the assembled theorem, and every constant that depends only on it.

THE DATUM (stated fully; ASSEMBLY sec.1.1's ANGULAR profile, hk2's RADIAL mollification):

    omega^theta_0(rho,phi) = -M * Theta(rho) * sgn(cos phi)
                                * min(1, phi_ax/delta) * min(1, |phi - pi/2|/delta_m)
    phi_ax = min(phi, pi-phi),   delta = 7.5 deg,  delta_m = 5 deg,
    Theta(rho) = (1/2)[ tanh((u - eps_r)/eps_r) - tanh((u - L + eps_r)/eps_r) ],  u = log(rho/rho0),
    eps_r = 0.25   (this is hk2's tanh radial ramp: rho Theta'(rho0) = 1/(2 eps_r) = 2 = rho0/(2 w0)
                    at w0 = 0.25 rho0, so the inner slope agrees exactly with hk2 (D-B)),
    rho0 = s sqrt(nu/M),  s = 1/sin delta,  R = rho0 e^L,  L = log(R/rho0).

Everything below is computed here, from these definitions, by scripts in this folder.
No code is imported from any other seat.  Cross-checks against other seats are reported
as checks, never used as inputs.

Outputs -> t1_results.json
"""
import json, math
import numpy as np
import mpmath as mp

mp.mp.dps = 25
DEG = mp.pi/180
DELTA = mp.mpf('7.5')*DEG
DM = mp.mpf(5)*DEG
EPS_R = mp.mpf('0.25')

OUT = {}
OUT["datum"] = {"delta_deg": 7.5, "delta_m_deg": 5.0, "eps_r": 0.25,
                "phi0_deg": 30.0, "radial_ramp": "tanh in log rho, hk2 w0 = 0.25 rho0"}

# --------------------------------------------------------------------- angular profiles
def h_axis(phi, delta):
    pax = min(phi, mp.pi - phi)
    return min(mp.mpf(1), pax/delta)

def h_eq(phi, dm):
    return min(mp.mpf(1), abs(phi - mp.pi/2)/dm)

def h_axis_tanh(phi, delta):
    return mp.tanh(mp.sin(phi)/mp.sin(delta))

def h_eq_tanh(phi, w):
    return mp.tanh(abs(mp.cos(phi))/w)

PROFILES = {
    "bangbang":        lambda p: mp.mpf(1),
    "taper_only":      lambda p: h_axis(p, DELTA),
    "THEOREM_datum":   lambda p: h_axis(p, DELTA)*h_eq(p, DM),          # ASSEMBLY sec.1.1
    "campaign_DC":     lambda p: h_axis_tanh(p, DELTA)*h_eq_tanh(p, mp.mpf('0.20')),
}

# --------------------------------------------------------------------- P_h and kappa
KINKS_V = {                      # v = sin phi at the kinks of each angular profile
    "bangbang": [], "taper_only": [mp.sin(DELTA)],
    "THEOREM_datum": [mp.sin(DELTA), mp.cos(DM)], "campaign_DC": [],
}

def P_h(lam, hfun, kinks=()):
    """a_ref = (M/2) P_h(lam);  P_h(lam) = 3 int_0^1 h(arcsin v) v^2 (A v^2 + B)^{-5/2} dv,
       A = lam^2 - lam^-4, B = lam^-4.  Derived in lower/prove-lagrangian sec.2; the
       substitution v = sin phi over the half-range, h symmetric about pi/2.
       The kinks of h are put on panel boundaries so the quadrature is smooth on each panel."""
    A = lam**2 - lam**-4
    B = lam**-4
    f = lambda v: hfun(mp.asin(v))*v**2*(A*v**2 + B)**mp.mpf('-2.5')
    pts = [mp.mpf(0)] + sorted(kinks) + [mp.mpf(1)]
    return 3*mp.quad(f, pts)

kap = {}
for name, hf in PROFILES.items():
    kap[name] = P_h(mp.mpf(1), hf, KINKS_V[name])/2
OUT["kappa"] = {k: float(v) for k, v in kap.items()}
OUT["kappa_str"] = {k: mp.nstr(v, 12) for k, v in kap.items()}

# the clock constants, at the theorem's datum
KAP = kap["THEOREM_datum"]
LOG32 = mp.log(mp.mpf(3)/2)
C_STAR = LOG32/KAP
C2 = 2*LOG32/KAP
LAM_MAX = mp.exp(3*C_STAR/4)
EPS_DELTA = 1 - 2*KAP
FLOOR = EPS_DELTA/(2*KAP)
OUT["clock"] = {"kappa_delta": float(KAP), "c_star": float(C_STAR), "c2": float(C2),
                "lam_max_apriori": float(LAM_MAX), "eps_delta_1m2kappa": float(EPS_DELTA),
                "eps_floor_L_to_inf": float(FLOOR),
                "c2_taper_only": float(2*LOG32/kap["taper_only"]),
                "c2_campaign_DC": float(2*LOG32/kap["campaign_DC"]),
                "floor_taper_only": float((1-2*kap["taper_only"])/(2*kap["taper_only"])),
                "floor_campaign_DC": float((1-2*kap["campaign_DC"])/(2*kap["campaign_DC"]))}

# P_1(lam) = lam exactly (control), and r_h = inf P_h(lam)/lam
OUT["control_P1"] = {"P_1(1.5)/1.5": float(P_h(mp.mpf('1.5'), PROFILES["bangbang"])/mp.mpf('1.5')),
                     "P_1(lam_max)/lam_max": float(P_h(LAM_MAX, PROFILES["bangbang"])/LAM_MAX)}

def r_h(name, lam_hi, n=60):
    hfun, kk = PROFILES[name], KINKS_V[name]
    best = None
    for i in range(n+1):
        lam = 1 + (lam_hi-1)*mp.mpf(i)/n
        val = P_h(lam, hfun, kk)/lam
        if best is None or val < best[0]:
            best = (val, lam)
    return best

rh_lm, arg_lm = r_h("THEOREM_datum", LAM_MAX)
rh_15, arg_15 = r_h("THEOREM_datum", mp.mpf('1.5'))
rh_lm_t, _ = r_h("taper_only", LAM_MAX)
OUT["r_h"] = {"theorem_datum_1_to_lam_max": float(rh_lm), "argmin_lam": float(arg_lm),
              "theorem_datum_1_to_1.5": float(rh_15),
              "taper_only_1_to_lam_max": float(rh_lm_t)}

# monotonicity of P_h on [1, 2.2]: needed for the quasimonotone comparison (u1 sec.5.1)
prev, minc, minc_win, turn = None, None, None, None
for i in range(0, 121):
    lam = 1 + mp.mpf('0.01')*i
    val = P_h(lam, PROFILES["THEOREM_datum"], KINKS_V["THEOREM_datum"])
    if prev is not None:
        inc = val - prev
        if minc is None or inc < minc:
            minc = inc
        if lam <= LAM_MAX and (minc_win is None or inc < minc_win):
            minc_win = inc
        if inc < 0 and turn is None:
            turn = lam
    prev = val
OUT["P_h_monotone"] = {"min_increment_step_0.01_on_[1,2.2]": float(minc),
                       "min_increment_on_[1,lam_max]": float(minc_win),
                       "first_decrease_at_lam": (float(turn) if turn is not None else None),
                       "lam_max": float(LAM_MAX)}

# --------------------------------------------------------------------- E_0 and Gfrak_0
# eta_0 = -(M/rho) F(phi) Theta(u),  F(phi) = sgn(cos phi) * g(phi)/sin phi,  u = log(rho/rho0)
# rho |eta_0|/M      = |F| Theta
# rho^2|grad eta_0|/M = sqrt( F^2 (Theta - dTheta/du)^2 + F'^2 Theta^2 )
def Theta_log(u, L, eps_r=EPS_R):
    return (mp.tanh((u - eps_r)/eps_r) - mp.tanh((u - L + eps_r)/eps_r))/2

def dTheta_log(u, L, eps_r=EPS_R):
    return (mp.sech((u - eps_r)/eps_r)**2 - mp.sech((u - L + eps_r)/eps_r)**2)/(2*eps_r)

D75, DM5 = 7.5*math.pi/180.0, 5.0*math.pi/180.0

def F_np(phi, kind):
    """F(phi) = sgn(cos phi) * g(phi)/sin phi, in float64."""
    s = np.sin(phi)
    if kind == "THEOREM_datum":
        pax = np.minimum(phi, np.pi - phi)
        g = np.minimum(1.0, pax/D75)*np.minimum(1.0, np.abs(phi - np.pi/2)/DM5)
    else:                                    # campaign (D-C): tanh taper x tanh equator
        g = np.tanh(s/math.sin(D75))*np.tanh(np.abs(np.cos(phi))/0.20)
    return np.sign(np.cos(phi))*g/s

def profile_sup(kind, breakpoints, L=1000.0, eps_r=0.25):
    """ess-sup of |F| and |F'|, and the joint sup of sqrt(F^2(Th-Th')^2 + F'^2 Th^2).
       Kinks are handled by evaluating one-sided limits explicitly."""
    hh = 1e-9
    n = 400000
    phis = list(np.pi*np.arange(1, n)/n)
    for bp in breakpoints:
        for sg in (-1.0, 1.0):
            phis.append(bp + sg*1e-8)
    phis = np.array([p for p in phis if 0.0 < p < np.pi])
    F = F_np(phis, kind)
    Fp = (F_np(phis + hh, kind) - F_np(phis - hh, kind))/(2*hh)
    # a central difference straddling a kink is a chord, never larger than the two
    # one-sided slopes, so the explicit kink points above supply the ess-sup.
    us = np.concatenate([np.arange(0.0, 3.0, 0.002), np.array([L/2])])
    supE, supG, arg = 0.0, 0.0, None
    for u in us:
        sech2 = lambda t: 0.0 if abs(t) > 350.0 else 1.0/math.cosh(t)**2
        Th = 0.5*(math.tanh((u - eps_r)/eps_r) - math.tanh((u - L + eps_r)/eps_r))
        dTh = (sech2((u - eps_r)/eps_r) - sech2((u - L + eps_r)/eps_r))/(2*eps_r)
        E = np.abs(F)*Th
        G = np.sqrt(F**2*(Th - dTh)**2 + Fp**2*Th**2)
        if E.max() > supE:
            supE = float(E.max())
        k = int(np.argmax(G))
        if G[k] > supG:
            supG, arg = float(G[k]), (float(phis[k]*180/np.pi), float(u))
    return (mp.mpf(supE), mp.mpf(supG), mp.mpf(float(np.abs(Fp).max())), arg)

bps_theorem = [D75, math.pi/2 - DM5, math.pi/2, math.pi/2 + DM5, math.pi - D75]
E0_t, G0_t, Fp_t, arg_t = profile_sup("THEOREM_datum", bps_theorem)
E0_c, G0_c, Fp_c, arg_c = profile_sup("campaign_DC", [])
OUT["datum_norms"] = {
    "THEOREM_datum": {"E0_sup_rho_eta": float(E0_t), "G0_sup_rho2_grad_eta": float(G0_t),
                      "sup_abs_Fprime": float(Fp_t), "argmax_phi_deg_u": arg_t,
                      "one_over_sin_delta": float(1/mp.sin(DELTA))},
    "campaign_DC":   {"E0_sup_rho_eta": float(E0_c), "G0_sup_rho2_grad_eta": float(G0_c),
                      "sup_abs_Fprime": float(Fp_c), "argmax_phi_deg_u": arg_c},
    "u2_quoted_for_DC": {"E0": 7.66060196, "G0": 20.9197070},
}

# --------------------------------------------------------------------- energy: C_E and the shift
def N_l(l):
    return (l+1.0)*(l+2.0)/(l+1.5)

def C32(l, t):
    """Gegenbauer C_l^{3/2} by the standard recurrence; C_0 = 1, C_1 = 3t."""
    c0 = np.ones_like(t)
    if l == 0:
        return c0
    c1 = 3.0*t
    if l == 1:
        return c1
    for n in range(1, l):
        c2 = (2.0*(n + 1.5)*t*c1 - (n + 2.0*1.5 - 1.0)*c0)/(n + 1.0)
        c0, c1 = c1, c2
    return c1

def g_np(phi, delta, dm):
    pax = np.minimum(phi, np.pi - phi)
    return np.minimum(1.0, pax/delta)*np.minimum(1.0, np.abs(phi - np.pi/2)/dm)

def H_hat(l, delta, dm, npts=400001):
    """Hhat_l = H_l/M = -(2/N_l) int_0^1 g(arccos t) sqrt(1-t^2) C_l^{3/2}(t) dt  (l odd).
       Composite Simpson with the kinks of g placed on nodes."""
    kinks = sorted(set([0.0, 1.0] + [math.cos(x) for x in
                   [delta, math.pi/2 - dm] if 0 < math.cos(x) < 1]))
    tot = 0.0
    for a, b in zip(kinks[:-1], kinks[1:]):
        n = npts//len(kinks) // 2 * 2
        t = np.linspace(a, b, n+1)
        f = g_np(np.arccos(np.clip(t, -1, 1)), delta, dm)*np.sqrt(np.maximum(1-t*t, 0.0))*C32(l, t)
        w = np.ones(n+1); w[1:-1:2] = 4.0; w[2:-1:2] = 2.0
        tot += (b-a)/n/3.0*np.dot(w, f)
    return -(2.0/N_l(l))*tot

def D_l_num(l, L, eps_r, sharp=False, du=2.0e-4):
    """D_l = int int_{x<y} Theta(x)Theta(y) e^{(l+4)x-(l-1)y} dx dy,  x = log(rho/R).
       Substituting b = -y, w = y-x >= 0:
           D_l = int db e^{-5b} Theta(-b) * I(L-b),
           I(u) = int_0^inf Theta(u-w) e^{-(l+4)w} dw,   I' = Theta - (l+4) I .
       The ODE recursion for I is unconditionally stable and never overflows, which the
       direct cumulative form does (e^{-(l-1)x} at x = -L blows up for large l)."""
    pad = 10.0*eps_r if not sharp else 1.0
    u = np.arange(-pad, L + pad + du, du)
    if sharp:
        Th = ((u >= 0.0) & (u <= L)).astype(float)
    else:
        Th = 0.5*(np.tanh((u - eps_r)/eps_r) - np.tanh((u - L + eps_r)/eps_r))
    k = l + 4.0
    E = math.exp(-k*du)
    I = np.empty_like(u)
    I[0] = Th[0]/k                      # equilibrium value for a locally constant Theta
    for i in range(1, len(u)):
        I[i] = I[i-1]*E + 0.5*du*(Th[i-1]*E + Th[i])
    return float(np.trapz(np.exp(5.0*(u - L))*Th*I, u))

def D_l_sharp_exact(l, L):
    if abs(l-1) < 1e-12:
        return (1.0/(l+4.0))*((1-math.exp(-5*L))/5.0 - math.exp(-5*L)*L)
    return (1.0/(l+4.0))*((1-math.exp(-5*L))/5.0
                          - (math.exp(-5*L) - math.exp(-(l+4.0)*L))/(l-1.0))

def C_E(delta, dm, L, eps_r, sharp=False, lmax=41, du=2.0e-4):
    tot = 0.0
    terms = {}
    for l in range(1, lmax+1, 2):
        Hh = H_hat(l, delta, dm)
        D = (D_l_sharp_exact(l, L) if (sharp and eps_r <= 0)
             else D_l_num(l, L, eps_r, sharp=sharp, du=du))
        terms[l] = N_l(l)*Hh*Hh*2.0*D/(2*l+3.0)
        tot += terms[l]
    return 2*math.pi*tot, terms

d75, dm5 = 7.5*math.pi/180, 5.0*math.pi/180
S_REC = 1.0/math.sin(d75)
LL = 40.0
CE_sharp_bb, _ = C_E(1e-300, 1e-300, LL, 0.0, sharp=True, lmax=201)  # bang-bang cap control
CE_sharp, _ = C_E(d75, dm5, LL, 0.0, sharp=True, lmax=201)           # ASSEMBLY sec.1.4 datum
CE_tanh, _ = C_E(d75, dm5, LL, 0.25, sharp=False)                   # THE THEOREM'S datum
CE_tanh_coarse, _ = C_E(d75, dm5, LL, 0.25, sharp=False, du=8.0e-4)
CE_tanh_l21, _ = C_E(d75, dm5, LL, 0.25, sharp=False, lmax=21)
Dctl = {str(l): {"recursion": D_l_num(l, 8.0, 0.0, sharp=True, du=2.0e-5),
                 "closed_form": D_l_sharp_exact(l, 8.0)} for l in [1, 3, 5, 9, 21]}
Dctl_maxrel = max(abs(v["recursion"]-v["closed_form"])/v["closed_form"] for v in Dctl.values())
OUT["energy"] = {
    "C_E_bangbang_control": CE_sharp_bb, "record_bangbang": 0.172403978,
    "C_E_sharp_edges_delta7.5_dm5": CE_sharp, "assembly_value": 0.17032563295410327,
    "C_E_THEOREM_datum_tanh_ramp_eps_r0.25": CE_tanh,
    "D_l_recursion_vs_closed_form_maxrel_sharp_du2e-5": Dctl_maxrel,
    "C_E_tanh_du_convergence_rel": abs(CE_tanh-CE_tanh_coarse)/CE_tanh,
    "C_E_tanh_lmax_convergence_rel": abs(CE_tanh-CE_tanh_l21)/CE_tanh,
    "s_record": S_REC,
    "shift_sharp": 2*math.log(S_REC) + 0.4*math.log(CE_sharp),
    "shift_THEOREM_datum": 2*math.log(S_REC) + 0.4*math.log(CE_tanh),
    "assembly_shift": 3.364345456912714,
    "note": "log Re_E = 2L + 2 log s + (2/5) log C_E ; kappa_delta does NOT enter the shift",
}

# ||omega_0||_2^2 / (M^2 R^3) for the BFG hypothesis check (finite at every L)
def omega_L2sq_over_R3(delta, dm, L, eps_r):
    """int |omega|^2 dx over R^3, with dx = 2 pi r dr dz = 2 pi rho^2 sin phi drho dphi."""
    nu_, nphi = 40001, 4001
    u = np.linspace(-12*eps_r, L + 12*eps_r, nu_)
    Th = 0.5*(np.tanh((u - eps_r)/eps_r) - np.tanh((u - L + eps_r)/eps_r))
    rad = np.trapz(Th**2*np.exp(3*(u - L)), u)          # int Theta^2 rho^3 dlog rho / R^3
    ph = np.linspace(1e-9, math.pi - 1e-9, nphi)
    ang = np.trapz(g_np(ph, delta, dm)**2*np.sin(ph), ph)
    return 2*math.pi*rad*ang

OUT["bfg_hypothesis"] = {
    "omega0_L2sq_over_M2R3_THEOREM_datum": omega_L2sq_over_R3(d75, dm5, LL, 0.25),
    "omega0_Linf_over_M": 1.0,
    "note": "omega_0 in L^2 cap L^inf, decay O(rho^{-8}) at infinity and O(rho^{8}) at 0",
}

with open("t1_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True)

for k in ["kappa_str", "clock", "r_h", "P_h_monotone", "datum_norms", "energy",
          "control_P1", "bfg_hypothesis"]:
    print(k, "=", json.dumps(OUT[k], indent=1, default=str))

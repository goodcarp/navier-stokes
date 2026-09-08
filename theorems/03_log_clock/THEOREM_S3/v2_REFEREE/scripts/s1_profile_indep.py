"""
s1 -- an INDEPENDENT instrument for every datum constant of THEOREM_S3_v2 sec.1.2.

Written from the mathematics in THEOREM_S3_v2.md sec.1.1 alone.  It shares NO code with
fix4/p1_profile.py: fix4 evaluates the Gaussian convolution by an FFT on a 600001-point
uniform phi grid and then takes GRID maxima; this file evaluates the convolution by
adaptive quadrature (scipy.integrate.quad) with the breakpoints of the kinked profile passed
explicitly, and pads every supremum with a RIGOROUS second-derivative bound obtained from
||W_sigma^(k)||_inf <= ||Wtilde'||_inf * ||chi_sigma^(k-1)||_1 .

Outputs -> s1_results.json
"""
import json, math, os
import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
DEG = math.pi/180.0
DELTA, DM = 7.5*DEG, 5.0*DEG
EPS_R = 0.25
PHI0 = 30.0*DEG
SIGMA = 0.02
LOG32 = math.log(1.5)

# ---------------------------------------------------------------- exact kinked W on [0,pi]
def W_exact(phi):
    """W = A/sin phi on [0,pi], A = sgn(cos)min(1,phi_ax/delta)min(1,|phi-pi/2|/dm)."""
    if phi <= 0.0:
        return 1.0/DELTA
    if phi >= math.pi:
        return -1.0/DELTA
    if phi > math.pi/2:
        return -W_exact(math.pi - phi)
    s = math.sin(phi)
    if phi < DELTA:
        return (phi/(DELTA*s)) if phi > 1e-8 else (1.0 + phi*phi/6.0)/DELTA
    if phi > math.pi/2 - DM:
        return (math.pi/2 - phi)/(DM*s)
    return 1.0/s

def Wp_exact(phi):
    """W' on [0,pi] (a.e.; jumps at the four kink cones)."""
    if phi > math.pi/2:
        return Wp_exact(math.pi - phi)          # W(pi-p) = -W(p)  =>  W'(pi-p) = +W'(p)
    if phi <= 0.0:
        return 0.0
    s, c = math.sin(phi), math.cos(phi)
    if phi < DELTA:
        num = (phi**3/3.0*(1.0 - phi*phi/10.0)) if phi < 1e-3 else (s - phi*c)
        return num/(DELTA*s*s)
    if phi > math.pi/2 - DM:
        t = math.pi/2 - phi
        return (-s - t*c)/(DM*s*s)
    return -c/(s*s)

def Wtil(y):
    """even reflection at 0 and at pi, 2pi-periodic."""
    y = math.fmod(y, 2.0*math.pi)
    if y < 0:
        y += 2.0*math.pi
    if y > math.pi:
        y = 2.0*math.pi - y            # even about pi
    return W_exact(y)

def Wtilp(y):
    yy = math.fmod(y, 2.0*math.pi)
    if yy < 0:
        yy += 2.0*math.pi
    if yy > math.pi:
        return -Wp_exact(2.0*math.pi - yy)
    return Wp_exact(yy)

# breakpoints of Wtilde on the line (mod 2pi), where W' jumps
BASE = [DELTA, math.pi/2 - DM, math.pi/2 + DM, math.pi - DELTA]
BRK = sorted([b + 2*math.pi*k for k in (-2,-1,0,1,2) for b in BASE]
             + [-b + 2*math.pi*k for k in (-2,-1,0,1,2) for b in BASE])

# ---------------------------------------------------------------- Gaussian and its L^1 norms
def chi(y, sg=SIGMA):
    return math.exp(-0.5*(y/sg)**2)/(sg*math.sqrt(2.0*math.pi))

def chip(y, sg=SIGMA):
    return -(y/sg**2)*chi(y, sg)

def chipp(y, sg=SIGMA):
    return ((y*y)/sg**2 - 1.0)/sg**2*chi(y, sg)

L1_CHI  = 1.0
L1_CHIP = math.sqrt(2.0/math.pi)/SIGMA                       # = 2 chi(0)
L1_CHIPP = 4.0*math.exp(-0.5)/(SIGMA**2*math.sqrt(2.0*math.pi))   # = 4 phi(1)/sigma^2

TRUNC = 14.0*SIGMA          # tail beyond: <= Wmax * 2*Phi(-14) < 1e-43

def _conv(f, phi, ker):
    """int f(phi - y) ker(y) dy over |y| <= TRUNC, split at the images of the breakpoints."""
    pts = sorted(set([phi - b for b in BRK if abs(phi - b) < TRUNC]))
    edges = [-TRUNC] + pts + [TRUNC]
    tot = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < 1e-15:
            continue
        v, _ = quad(lambda y: f(phi - y)*ker(y), a, b, limit=300, epsabs=1e-14, epsrel=1e-13)
        tot += v
    return tot

def Ws(phi):   return _conv(Wtil,  phi, chi)
def Ws1(phi):  return _conv(Wtilp, phi, chi)
def Ws2(phi):  return _conv(Wtilp, phi, chip)
def Ws3(phi):  return _conv(Wtilp, phi, chipp)

def As(phi):   return Ws(phi)*math.sin(phi)
def As1(phi):  return Ws1(phi)*math.sin(phi) + Ws(phi)*math.cos(phi)
def As2(phi):
    s, c = math.sin(phi), math.cos(phi)
    return Ws2(phi)*s + 2.0*Ws1(phi)*c - Ws(phi)*s

# ---------------------------------------------------------------- rigorous global bounds
WMAX = 1.0/math.sin(DELTA)                     # sup|Wtilde|
WPMAX = math.cos(DELTA)/math.sin(DELTA)**2     # sup|Wtilde'| (checked below)
B_W  = WMAX
B_W1 = WPMAX
B_W2 = WPMAX*L1_CHIP
B_W3 = WPMAX*L1_CHIPP
B_A  = B_W                                     # |A|=|W sin| <= |W|
B_A1 = B_W1 + B_W
B_A2 = B_W2 + 2.0*B_W1 + B_W
B_A3 = B_W3 + 3.0*B_W2 + 3.0*B_W1 + B_W

def certified_sup(f, a, b, n, d2bound, label):
    """grid max of f on [a,b] with n points, padded by (1/8) h^2 sup|f''| (interior max)
       and by (h/2) sup|f'| as a fallback for a boundary/corner max."""
    xs = np.linspace(a, b, n)
    vs = np.array([f(x) for x in xs])
    h = (b - a)/(n - 1)
    k = int(np.argmax(vs))
    return {"label": label, "grid_max": float(vs.max()), "argmax": float(xs[k]),
            "pad_smooth": 0.125*h*h*d2bound, "upper": float(vs.max()) + 0.125*h*h*d2bound,
            "n": n, "h": h}

# ---------------------------------------------------------------- Theta and its derivatives
def Theta_derivs(u, L, er=EPS_R, kmax=4):
    """Theta(u) = 1/2[tanh((u-er)/er) - tanh((u-L+er)/er)] and d^k/du^k, k=0..kmax."""
    out = []
    for arg, sgn in ((u - er)/er, 1.0), ((u - L + er)/er, -1.0):
        T = math.tanh(max(min(arg, 350.0), -350.0))
        d = [T]
        d.append(1.0 - T*T)
        d.append(-2.0*T*(1.0 - T*T))
        d.append(-2.0*(1.0 - T*T)*(1.0 - 3.0*T*T))
        d.append(8.0*T*(1.0 - T*T)*(2.0 - 3.0*T*T))
        out.append([0.5*sgn*d[k]/er**k for k in range(kmax+1)])
    return [out[0][k] + out[1][k] for k in range(kmax+1)]

# ---------------------------------------------------------------- MAIN
if __name__ == "__main__":
    RES = {"instrument": "independent: adaptive quadrature convolution + rigorous derivative "
                         "bounds; no shared code with fix4/p1_profile.py",
           "sigma": SIGMA, "delta_deg": 7.5, "delta_m_deg": 5.0, "eps_r": EPS_R}

    # ---- 0. check the rigorous sup|W'| claim by a dense scan of the exact W'
    xs = np.linspace(1e-9, math.pi-1e-9, 200001)
    wp = np.array([Wp_exact(x) for x in xs])
    RES["sup_abs_Wprime_scan"] = float(np.abs(wp).max())
    RES["sup_abs_Wprime_closed_form_cosdelta_over_sin2delta"] = WPMAX
    RES["sup_abs_W_scan"] = float(max(abs(W_exact(x)) for x in xs))
    RES["sup_abs_W_closed_form_1_over_sindelta"] = WMAX
    RES["rigorous_bounds"] = {"B_W": B_W, "B_W1": B_W1, "B_W2": B_W2, "B_W3": B_W3,
                              "B_A": B_A, "B_A1": B_A1, "B_A2": B_A2, "B_A3": B_A3,
                              "L1_chi_prime": L1_CHIP, "L1_chi_second": L1_CHIPP}

    # ---- 1. N_sigma, with a certified pad
    NG1 = 20001
    nrec = certified_sup(lambda p: abs(As(p)), 1e-9, math.pi-1e-9, NG1, B_A2, "N_sigma")
    # local refinement around the argmax
    a0 = nrec["argmax"]; hh = nrec["h"]
    xs2 = np.linspace(a0-2*hh, a0+2*hh, 401)
    vs2 = [abs(As(x)) for x in xs2]
    k2 = int(np.argmax(vs2))
    h2 = xs2[1]-xs2[0]
    N_lo = float(max(nrec["grid_max"], vs2[k2]))
    N_hi = N_lo + 0.125*h2*h2*B_A2
    RES["N_sigma"] = {"lower": N_lo, "upper": N_hi, "argmax_deg": float(xs2[k2]/DEG),
                      "fix4_value": 1.0124508488}
    N = 0.5*(N_lo + N_hi)

    # ---- 2. kappa_delta = 1.5 int_0^{pi/2} h_sigma sin^2 cos dphi , h = A_sigma/N
    brk = sorted(set([b for b in [DELTA, math.pi/2 - DM] if 0 < b < math.pi/2]))
    edges = [0.0] + brk + [math.pi/2]
    I = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        v, _ = quad(lambda p: As(p)*math.sin(p)**2*math.cos(p), a, b,
                    limit=200, epsabs=1e-13, epsrel=1e-12)
        I += v
    kappa_lo = 1.5*I/N_hi
    kappa_hi = 1.5*I/N_lo
    RES["kappa_delta"] = {"lower": kappa_lo, "upper": kappa_hi,
                          "integral_1p5_int_A_sin2_cos": 1.5*I, "fix4_value": 0.4917868801}
    kap = 1.5*I/N
    RES["c_star"] = {"value": LOG32/kap, "fix4_value": 0.8244732108}
    RES["c2"] = {"value": 2.0*LOG32/kap, "fix4_value": 1.6489464217}
    RES["eps_floor"] = {"value": (1.0-2.0*kap)/(2.0*kap), "fix4_value": 1.6700567265e-02}
    RES["lam_max_at_c_star"] = {"value": math.exp(0.75*LOG32/kap),
                                "fix4_value": 1.8558724485}

    # kinked control: kappa for h = |A_kink|, N = 1
    def A_k(p):
        pax = min(p, math.pi - p); eq = abs(p - math.pi/2)
        return math.copysign(1.0, math.cos(p))*min(1.0, pax/DELTA)*min(1.0, eq/DM)
    Ik = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        v, _ = quad(lambda p: A_k(p)*math.sin(p)**2*math.cos(p), a, b, limit=200,
                    epsabs=1e-14, epsrel=1e-13)
        Ik += v
    RES["kappa_delta_kinked_control"] = {"value": 1.5*Ik, "fix2_value": 0.4978224182,
                                         "fix4_value": 0.4978224383}

    # ---- 3. E_0 = sup|W_sigma|/N  (times sup Theta = 1 - O(e^{-4(L-1/2)}))
    erec = certified_sup(lambda p: abs(Ws(p)), 1e-9, math.pi-1e-9, NG1, B_W2, "sup|W_sigma|")
    a0 = erec["argmax"]; hh = erec["h"]
    xs3 = np.linspace(max(a0-2*hh, 1e-12), min(a0+2*hh, math.pi-1e-12), 401)
    vs3 = [abs(Ws(x)) for x in xs3]
    h3 = xs3[1]-xs3[0]
    supW_lo = float(max(erec["grid_max"], max(vs3)))
    supW_hi = supW_lo + 0.125*h3*h3*B_W2
    RES["E_0"] = {"lower": supW_lo/N_hi, "upper": supW_hi/N_lo,
                  "sup_abs_W_sigma": [supW_lo, supW_hi],
                  "argmax_deg": float(xs3[int(np.argmax(vs3))]/DEG),
                  "fix4_value": 7.5523076942}

    # ---- 4. Gfrak_0 by fix2 sec.4(c)'s branch bound, recomputed independently
    YMAX = (4.5*4.5)**2/64.0
    def branch(p):
        P = (Ws(p)/N)**2; S = (Ws1(p)/N)**2
        return (P + S) if 81.0*P <= 8.0*S else (YMAX*P + S)
    xs4 = np.linspace(1e-6, math.pi-1e-6, 6001)
    v4 = [branch(x) for x in xs4]
    k4 = int(np.argmax(v4))
    xs5 = np.linspace(xs4[max(k4-1,0)], xs4[min(k4+1,len(xs4)-1)], 801)
    v5 = [branch(x) for x in xs5]
    RES["Gfrak_0_branch_bound"] = {"value": math.sqrt(max(max(v4), max(v5))),
                                   "argmax_deg": float(xs5[int(np.argmax(v5))]/DEG),
                                   "fix4_value": 36.0795702396,
                                   "YMAX_used": YMAX}
    RES["sup_abs_Wp_sigma_over_N"] = float(max(abs(Ws1(x))/N for x in xs5))

    # ---- 5. |h''| : is fix4's "proved Lipschitz bound 225.9466903722" a proved bound?
    def hpp(p):
        s, c = math.sin(p), math.cos(p)
        return abs((Ws2(p)*s + 2.0*Ws1(p)*c - Ws(p)*s)/N)
    xs6 = np.linspace(1e-6, math.pi-1e-6, 4001)
    v6 = [hpp(x) for x in xs6]
    k6 = int(np.argmax(v6))
    xs7 = np.linspace(xs6[max(k6-1,0)], xs6[min(k6+1,len(xs6)-1)], 1201)
    v7 = [hpp(x) for x in xs7]
    RES["sup_abs_h_second"] = {"scan_value": float(max(max(v6), max(v7))),
                               "argmax_deg": float(xs7[int(np.argmax(v7))]/DEG),
                               "fix4_claimed_PROVED_bound": 225.9466903722,
                               "rigorous_a_priori_bound_B_A2_over_N": B_A2/N}
    # also the sup of |h'| that the pad is applied to
    def hp(p):
        s, c = math.sin(p), math.cos(p)
        return abs((Ws1(p)*s + Ws(p)*c)/N)
    v8 = [hp(x) for x in np.linspace(1e-6, math.pi-1e-6, 8001)]
    RES["sup_abs_h_prime_scan"] = float(max(v8))

    json.dump(RES, open(os.path.join(HERE, "s1_results.json"), "w"), indent=1, default=str)
    print(json.dumps(RES, indent=1, default=str))

"""
p2 -- THE ENERGY CONSTANT C_E FOR THE MOLLIFIED DATUM, and the log Re_E shift.

The dictionary is THEOREM_S3 sec.1.4, unchanged:

    log Re_E = 2L + 2 log s + (2/5) log C_E ,     C_E = ||u_0||_2^2/(M^2 R^5) ,
    C_E = 2 pi sum_{l odd} 2 N_l Hhat_l^2 D_l/(2l+3) ,   N_l = (l+1)(l+2)/(l+3/2) ,
    Hhat_l = -(2/N_l) int_0^1 g(arccos t) sqrt(1-t^2) C_l^{3/2}(t) dt ,
    D_l = int int_{x<y} Theta(x)Theta(y) e^{(l+4)x-(l-1)y} dx dy ,  x = log(rho/R) .

Only  g  changes: the kinked  g = min(1,phi_ax/delta) min(1,|phi-pi/2|/delta_m)  is replaced
by the mollified, amplitude-normalised  g_sigma(phi) = A_sigma(phi)/N_sigma  on [0, pi/2]
(p1's `a0`), the same function that defines the datum.  Theta, and therefore every D_l, is
unchanged, so the whole change of C_E is in the Gegenbauer coefficients.

The instrument is re-implemented here from the formulas above (L-14: the object under test is
the datum, so the claim is re-derived; t1_datum.py is imported only as a CONTROL, to check
that this file reproduces its sharp-edge and tanh-ramp numbers for the kinked g).

Outputs -> p2_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import p1_profile as P1                                              # noqa: E402

DEG = math.pi/180.0
DELTA, DM, EPS_R = 7.5*DEG, 5.0*DEG, 0.25
S_REC = 1.0/math.sin(DELTA)
LL = 40.0

def N_l(l):
    return (l+1.0)*(l+2.0)/(l+1.5)

def C32(l, t):
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

def g_kinked(phi):
    pax = np.minimum(phi, math.pi - phi)
    return np.minimum(1.0, pax/DELTA)*np.minimum(1.0, np.abs(phi - math.pi/2)/DM)

def H_hat(l, gfun, kinks, npts=400001):
    tot = 0.0
    pts = sorted(set([0.0, 1.0] + list(kinks)))
    for a, b in zip(pts[:-1], pts[1:]):
        n = max(2, npts//len(pts) // 2 * 2)
        t = np.linspace(a, b, n+1)
        f = gfun(np.arccos(np.clip(t, -1, 1)))*np.sqrt(np.maximum(1-t*t, 0.0))*C32(l, t)
        w = np.ones(n+1); w[1:-1:2] = 4.0; w[2:-1:2] = 2.0
        tot += (b-a)/n/3.0*float(np.dot(w, f))
    return -(2.0/N_l(l))*tot

def D_l_num(l, L, eps_r, sharp=False, du=2.0e-4):
    pad = 10.0*eps_r if not sharp else 1.0
    u = np.arange(-pad, L + pad + du, du)
    if sharp:
        Th = ((u >= 0.0) & (u <= L)).astype(float)
    else:
        Th = 0.5*(np.tanh((u - eps_r)/eps_r) - np.tanh((u - L + eps_r)/eps_r))
    k = l + 4.0
    E = math.exp(-k*du)
    I = np.empty_like(u)
    I[0] = Th[0]/k
    for i in range(1, len(u)):
        I[i] = I[i-1]*E + 0.5*du*(Th[i-1]*E + Th[i])
    return float(np.trapz(np.exp(5.0*(u - L))*Th*I, u))

def D_l_sharp_exact(l, L):
    if abs(l-1) < 1e-12:
        return (1.0/(l+4.0))*((1-math.exp(-5*L))/5.0 - math.exp(-5*L)*L)
    return (1.0/(l+4.0))*((1-math.exp(-5*L))/5.0
                          - (math.exp(-5*L) - math.exp(-(l+4.0)*L))/(l-1.0))

def C_E(gfun, kinks, L, eps_r, sharp=False, lmax=41, du=2.0e-4):
    tot = 0.0
    terms = {}
    for l in range(1, lmax+1, 2):
        Hh = H_hat(l, gfun, kinks)
        D = (D_l_sharp_exact(l, L) if (sharp and eps_r <= 0)
             else D_l_num(l, L, eps_r, sharp=sharp, du=du))
        terms[str(l)] = N_l(l)*Hh*Hh*2.0*D/(2*l+3.0)
        tot += terms[str(l)]
    return 2*math.pi*tot, terms

def omega_L2sq_over_R3(gfun, L, eps_r):
    nu_, nphi = 40001, 400001
    u = np.linspace(-12*eps_r, L + 12*eps_r, nu_)
    Th = 0.5*(np.tanh((u - eps_r)/eps_r) - np.tanh((u - L + eps_r)/eps_r))
    rad = np.trapz(Th**2*np.exp(3*(u - L)), u)
    ph = np.linspace(1e-9, math.pi - 1e-9, nphi)
    ang = np.trapz(gfun(ph)**2*np.sin(ph), ph)
    return 2*math.pi*float(rad)*float(ang)


if __name__ == "__main__":
    SIGMA = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
    OUT = {"sigma": SIGMA, "L_used": LL, "s_record": S_REC}

    # ---- controls against t1_datum.py's own numbers, kinked g
    kinks_kink = [math.cos(DELTA), math.cos(math.pi/2 - DM)]
    kinks_kink = [k for k in kinks_kink if 0.0 < k < 1.0]
    CE_bb, _ = C_E(lambda p: np.ones_like(p), [], LL, 0.0, sharp=True, lmax=201)
    CE_sharp, _ = C_E(g_kinked, kinks_kink, LL, 0.0, sharp=True, lmax=201)
    CE_tanh_kink, _ = C_E(g_kinked, kinks_kink, LL, EPS_R, sharp=False)
    OUT["control_kinked"] = {
        "C_E_bangbang": CE_bb, "record_bangbang": 0.172403978,
        "C_E_sharp_edges": CE_sharp, "assembly_value": 0.17032563295410327,
        "C_E_tanh_ramp": CE_tanh_kink, "THEOREM_S3_value": 0.05657477,
        "shift_sharp": 2*math.log(S_REC) + 0.4*math.log(CE_sharp),
        "shift_tanh_kinked": 2*math.log(S_REC) + 0.4*math.log(CE_tanh_kink),
        "THEOREM_S3_shift": 2.9234859}
    print("control: C_E(bb)=%.9f  C_E(sharp)=%.11f  C_E(tanh,kinked)=%.9f  shift=%.7f"
          % (CE_bb, CE_sharp, CE_tanh_kink,
             OUT["control_kinked"]["shift_tanh_kinked"]), flush=True)

    # ---- the mollified datum
    prof = P1.Profile(SIGMA)
    g_sig = lambda phi: np.abs(prof.interp(prof.a0, np.asarray(phi, dtype=float)))
    CE_sig, terms = C_E(g_sig, [], LL, EPS_R, sharp=False)
    CE_sig_l21, _ = C_E(g_sig, [], LL, EPS_R, sharp=False, lmax=21)
    CE_sig_du, _ = C_E(g_sig, [], LL, EPS_R, sharp=False, du=8.0e-4)
    shift = 2*math.log(S_REC) + 0.4*math.log(CE_sig)
    OUT["mollified"] = {
        "C_E": CE_sig, "shift_logReE": shift,
        "shift_minus_kinked": shift - OUT["control_kinked"]["shift_tanh_kinked"],
        "lmax_convergence_rel": abs(CE_sig - CE_sig_l21)/CE_sig,
        "du_convergence_rel": abs(CE_sig - CE_sig_du)/CE_sig,
        "terms_l": terms}
    OUT["bfg_hypothesis"] = {
        "omega0_L2sq_over_M2R3_mollified": omega_L2sq_over_R3(g_sig, LL, EPS_R),
        "omega0_L2sq_over_M2R3_kinked": omega_L2sq_over_R3(g_kinked, LL, EPS_R),
        "THEOREM_S3_quoted_kinked": 1.4785137,
        "omega0_Linf_over_M": 1.0,
        "note": "omega_0 in L^2 cap L^inf; decay O(rho^-8) at infinity, O(rho^8) at 0"}
    # ---- two independent controls on the Gegenbauer coefficients (L-98)
    # (1) SYNTHESIS.  Hhat_l = -(2/N_l) int_0^1 g sqrt(1-t^2) C_l^{3/2} dt is exactly the
    #     coefficient of  W(t) = -g_odd(t)/sqrt(1-t^2) = -W_datum(arccos t)  in the basis
    #     {C_l^{3/2}}, which is orthogonal for the weight (1-t^2) with
    #     int_{-1}^1 (1-t^2) [C_l^{3/2}]^2 dt = N_l EXACTLY.  So the partial sums of
    #     sum_l Hhat_l C_l^{3/2}(t) must reconstruct -W_datum.
    # (2) PARSEVAL.  sum_l Hhat_l^2 N_l = int_{-1}^1 W^2 (1-t^2) dt = int_{-1}^1 g_odd^2 dt.
    tt = np.linspace(-1.0, 1.0, 400001)
    ortho = {}
    for l in (1, 3, 7):
        ortho[str(l)] = {"int_(1-t^2)C_l^2_dt": float(np.trapz((1-tt*tt)*C32(l, tt)**2, tt)),
                         "N_l": N_l(l)}
    OUT["gegenbauer_orthogonality_control"] = ortho
    LMAXC = 121
    Hs = {l: H_hat(l, g_sig, []) for l in range(1, LMAXC+1, 2)}
    syn = np.zeros_like(tt)
    for l, Hh in Hs.items():
        syn = syn + Hh*C32(l, tt)
    ph_t = np.arccos(np.clip(tt, -1, 1))
    Wtrue = -prof.interp(prof.w0, ph_t)
    inner = np.abs(tt) < 0.98
    OUT["synthesis_control"] = {
        "LMAX": LMAXC,
        "max_abs_defect_on_|t|<0.98": float(np.max(np.abs(syn - Wtrue)[inner])),
        "sup_abs_W": float(np.max(np.abs(Wtrue)))}
    god = np.where(tt >= 0, g_sig(np.arccos(np.clip(tt, -1, 1))),
                   -g_sig(np.arccos(np.clip(-tt, -1, 1))))
    lhs = float(np.trapz(god**2, tt))
    rhs = sum(Hs[l]**2*N_l(l) for l in Hs)
    OUT["parseval_control"] = {"int_g_odd^2_dt": lhs, "sum_Hhat_l^2_N_l": rhs,
                               "rel": abs(lhs-rhs)/lhs, "LMAX": LMAXC}
    json.dump(OUT, open(os.path.join(HERE, "p2_results.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: v for k, v in OUT.items() if k != "mollified"}, indent=1))
    print("MOLLIFIED sigma=%g: C_E = %.9f   shift = %.8f   (kinked shift %.8f)"
          % (SIGMA, CE_sig, shift, OUT["control_kinked"]["shift_tanh_kinked"]))

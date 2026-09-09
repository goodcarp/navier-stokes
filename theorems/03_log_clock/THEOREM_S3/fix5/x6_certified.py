"""
x6 -- UNIT F-6.  Every grid-valued entry of THEOREM_S3_v2 sec.3 replaced by a certified bound.

WHAT WAS WRONG.  sec.3 item 2 lists E_0, Gfrak_0, Psi_sigma(0) and N_sigma as ENCLOSED, but
they are suprema read off a 600001-point grid with NO pad.  Item 1's r_h certificate is a
genuine Lipschitz-plus-grid argument, but one of its inputs, the "proved Lipschitz bound
|h_sigma''| <= 225.9466903722", is itself a np.max over the same grid, i.e. a scan value.

WHAT IS DONE HERE.  Every bound below is (value on a grid) + (an explicitly proved pad).  The
proved inputs are the a priori bounds

    ||W_sigma^{(k)}||_inf <= ||Wtil'||_inf ||chi_sigma^{(k-1)}||_1  (k >= 1),
    ||W_sigma||_inf <= ||Wtil||_inf = 1/sin delta ,
    ||A_sigma^{(k)}||_inf <= sum_j C(k,j) ||W_sigma^{(j)}||_inf     (|sin^{(m)}| <= 1) ,

with the L^1 norms of the Gaussian derivatives in closed form (fx_profile.chi_L1, each the total
variation of the previous derivative).  For a SMOOTH interior maximum on a uniform grid of step
h the pad is (1/8) h^2 sup|f''| (the nearest grid point is within h/2 of the argmax and f' there
vanishes); where the maximised function is only piecewise smooth or the maximum can sit at a
region boundary, the pad is (h/2) sup|f'| instead, which needs no interior hypothesis.

r_h and min P_h' are re-run through fix4's OWN certificate (`p1_profile.certified_r_h`,
`Ph_prime_lower`, imported unchanged: L-14's inverse convention, the object under test is that
certificate) with the scan Lipschitz constant replaced by the proved one, at lam_max = e^{3c/4}
for c at c_*, at the certified cap, and at the self-consistent windows the budget actually uses.

Outputs -> x6_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "copies", "fix4"))
import fx_profile as FX
import p1_profile as P1                                     # fix4's own instrument, unchanged

DEG = FX.DEG
SIGMA = 0.02
CONV = FX.Conv(SIGMA)
B = FX.W_bounds(SIGMA)
BW, BA = B["BW"], B["BA"]
LOG32 = math.log(1.5)

# ---------------------------------------------------------------- tanh derivative suprema
_TANH_POLY = {1: [1.0, 0.0, -1.0], 2: [2.0, 0.0, -2.0, 0.0], 3: [-6.0, 0.0, 8.0, 0.0, -2.0],
              4: [24.0, 0.0, -40.0, 0.0, 16.0, 0.0],
              5: [-120.0, 0.0, 240.0, 0.0, -136.0, 0.0, 16.0]}
# coefficients in DESCENDING powers of T = tanh(u); tanh^{(k)}(u) = p_k(T)

def sup_tanh_deriv(k, n=2000001):
    p = np.array(_TANH_POLY[k])
    T = np.linspace(-1.0, 1.0, n)
    v = float(np.max(np.abs(np.polyval(p, T))))
    dp = np.polyder(p)
    return v + (2.0/(n-1))/2.0*float(np.sum(np.abs(dp)))     # grid max + (h/2) sup|p'|

SUP_T = {k: sup_tanh_deriv(k) for k in range(1, 6)}

def Theta_bounds(er=FX.EPS_R):
    """rigorous sup bounds on |Theta^{(k)}|, k = 0..3, and on |rad| = |T2 + T1 - 2 T0|."""
    st = {0: 1.0}
    st.update({k: SUP_T[k] for k in range(1, 6)})
    b = {k: (0.5*(st[k] + st[k])/er**k if k > 0 else 1.0) for k in range(0, 4)}
    b["rad"] = b[2] + b[1] + 2.0*b[0]
    b["rad_prime"] = (0.5*2*st[3]/er**3) + b[2] + 2.0*b[1]
    return b

TH = Theta_bounds()

def Theta_derivs_scalar(u, L, er=FX.EPS_R, kmax=2):
    tot = [np.zeros_like(u) for _ in range(kmax+1)]
    for arg, sgn in (((u - er)/er, 1.0), ((u - L + er)/er, -1.0)):
        T = np.tanh(np.clip(arg, -350.0, 350.0))
        d = [T, 1.0 - T*T, -2.0*T*(1.0 - T*T)]
        for k in range(kmax+1):
            tot[k] = tot[k] + 0.5*sgn*d[k]/er**k
    return tot

# ---------------------------------------------------------------- the certified 1-D suprema
def certify(vals, h, pad_const, mode):
    g = float(np.max(vals))
    pad = 0.125*h*h*pad_const if mode == "d2" else 0.5*h*pad_const
    return {"grid_max": g, "pad": pad, "lower": g, "upper": g + pad, "mode": mode, "h": h}

if __name__ == "__main__":
    OUT = {"sigma": SIGMA, "proved_bounds": {"BW": BW, "BA": BA, "chi_L1": B["chi_L1"],
                                             "Wtil_sup": FX.WMAX, "Wtil_prime_sup": FX.WPMAX},
           "sup_tanh_derivatives": SUP_T, "Theta_bounds": TH}

    # ---------------- the phi grid (my own Gauss-Legendre evaluator, not fix4's FFT)
    NPHI = 100001
    phi = np.linspace(0.0, math.pi, NPHI)
    h = phi[1] - phi[0]
    Wd = CONV.derivs(phi, kmax=3)
    Ad = FX.A_derivs_from_W(phi, Wd + [np.zeros_like(phi)], kmax=2)
    print("grid built, h = %.4e" % h, flush=True)

    # ---- N_sigma = sup |A_sigma| : smooth interior maximum (|A_sigma| ~ 1 there, so |.| is C^2)
    Nrec = certify(np.abs(Ad[0]), h, BA[2], "d2")
    k = int(np.argmax(np.abs(Ad[0])))
    loc = np.linspace(phi[max(k-1, 0)], phi[min(k+1, NPHI-1)], 2001)
    hl = loc[1] - loc[0]
    Aloc = FX.A_derivs_from_W(loc, CONV.derivs(loc, kmax=0) + [np.zeros_like(loc)]*4, kmax=0)[0]
    Nloc = certify(np.abs(Aloc), hl, BA[2], "d2")
    # The LOWER bound may be raised by any local refinement (a grid value is always a lower
    # bound on a supremum).  The UPPER bound must stay the GLOBAL one: a refinement of one
    # window says nothing about the rest of [0,pi].
    N_lo, N_hi = max(Nrec["lower"], Nloc["lower"]), Nrec["upper"]
    OUT["N_sigma"] = {"lower": N_lo, "upper": N_hi, "global": Nrec, "local": Nloc,
                      "argmax_deg": float(loc[int(np.argmax(np.abs(Aloc)))]/DEG),
                      "fix4_value": 1.0124508488,
                      "fix4_inside": bool(N_lo <= 1.0124508488 <= N_hi),
                      "referee_enclosure": [1.0124508487, 1.0124508495]}
    print("N_sigma in [%.12f, %.12f]" % (N_lo, N_hi), flush=True)

    # ---- sup |W_sigma| and E_0 = sup rho|eta_0|/M = sup|W_sigma| sup Theta / N_sigma
    Wrec = certify(np.abs(Wd[0]), h, BW[2], "d2")
    kk = int(np.argmax(np.abs(Wd[0])))
    loc = np.linspace(phi[max(kk-1, 0)], phi[min(kk+1, NPHI-1)], 2001)
    hl = loc[1] - loc[0]
    Wloc = certify(np.abs(CONV.derivs(loc, kmax=0)[0]), hl, BW[2], "d2")
    sW_lo = max(Wrec["lower"], Wloc["lower"])
    sW_hi = Wrec["upper"]                       # global, for the same reason as N_sigma
    OUT["sup_abs_W_sigma"] = {"lower": sW_lo, "upper": sW_hi,
                              "argmax_deg": float(loc[int(np.argmax(np.abs(
                                  CONV.derivs(loc, kmax=0)[0])))]/DEG)}
    OUT["E_0"] = {"lower": sW_lo/N_hi, "upper": sW_hi/N_lo, "fix4_value": 7.5523076942,
                  "fix4_inside": bool(sW_lo/N_hi <= 7.5523076942 <= sW_hi/N_lo),
                  "note": "sup Theta <= 1 is used; the true sup Theta is 1 - O(e^{-2L/eps_r})",
                  "referee_enclosure": [7.5523076888, 7.5523076951]}
    print("E_0 in [%.10f, %.10f]" % (sW_lo/N_hi, sW_hi/N_lo), flush=True)

    # ---- Gfrak_0 by fix2 sec.4(c)'s branch bound, certified cell by cell
    YMAX = P1.YMAX
    w0 = Wd[0]/N_lo
    w1 = Wd[1]/N_lo
    P_ = (Wd[0]/N_hi)**2, (Wd[0]/N_lo)**2         # (lower, upper) since |W|/N is bracketed
    S_ = (Wd[1]/N_hi)**2, (Wd[1]/N_lo)**2
    Plo, Phi_ = P_[0], P_[1]
    Slo, Shi = S_[0], S_[1]
    # the branch condition 81 P <= 8 S is decided on the grid; a cell where it could flip
    # inside the cell is bounded by the LARGER branch, which is safe.
    Lip_cond = (81.0*2.0*BW[0]*BW[1] + 8.0*2.0*BW[1]*BW[2])/N_lo**2
    cond = 8.0*Slo - 81.0*Phi_                                    # > margin => branch 1 for sure
    margin = 0.5*h*Lip_cond
    br1 = Phi_ + Shi
    br2 = YMAX*Phi_ + Shi
    branch_up = np.where(cond > margin, br1, br2)
    Bbr2 = (2.0*YMAX*(BW[1]**2 + BW[0]*BW[2]) + 2.0*(BW[2]**2 + BW[1]*BW[3]))/N_lo**2
    Bbr1 = (2.0*YMAX*BW[0]*BW[1] + 2.0*BW[1]*BW[2])/N_lo**2
    # The maximum of the piecewise branch function over [0,pi] sits either at an interior
    # critical point of ONE smooth branch, where the nearest grid point is within h/2 and the
    # derivative vanishes, so the pad is (1/8) h^2 sup|branch''| ; or at a boundary of the two
    # branch regions, i.e. inside a cell where 81P - 8S can change sign, where the pad is the
    # Lipschitz one (h/2) sup|branch'|.  Cells of the second kind are flagged and padded
    # separately, so the combined bound needs no interior hypothesis anywhere.
    switch = np.abs(cond) <= margin + 0.5*h*Lip_cond
    pad_int = 0.125*h*h*Bbr2
    pad_sw = 0.5*h*Bbr1
    cand_int = float(np.max(np.where(~switch, branch_up, -np.inf))) + pad_int
    cand_sw = (float(np.max(np.where(switch, branch_up, -np.inf))) + pad_sw
               if switch.any() else -np.inf)
    gbr = float(np.max(branch_up))
    kb = int(np.argmax(branch_up))
    G0_up = math.sqrt(max(cand_int, cand_sw))
    G0_lo = math.sqrt(float(np.max(np.where(cond > -margin, Plo + Slo, YMAX*Plo + Slo))))
    OUT["Gfrak_0"] = {"lower": G0_lo, "upper": G0_up, "branch_grid_max": gbr,
                      "pad_interior_d2": pad_int, "pad_switch_d1": pad_sw,
                      "bound_from_interior_cells": cand_int, "bound_from_switch_cells": cand_sw,
                      "n_switch_cells": int(switch.sum()),
                      "argmax_in_switch_cell": bool(switch[kb]),
                      "argmax_deg": float(phi[kb]/DEG),
                      "Lip_branch": Bbr1, "curv_branch": Bbr2,
                      "Lip_condition": Lip_cond, "margin": float(margin),
                      "branch_at_argmax": ("81P<=8S" if cond[kb] > margin else "81P>8S"),
                      "cond_at_argmax": float(cond[kb]),
                      "YMAX": YMAX, "fix4_value": 36.0795702396,
                      "fix4_inside": bool(G0_lo <= 36.0795702396 <= G0_up)}
    print("Gfrak_0 in [%.10f, %.10f]  (fix4 36.0795702396)" % (G0_lo, G0_up), flush=True)

    # ---- sup |h'| and the |h''| Lipschitz bound
    hp = np.abs(Ad[1])/N_lo
    hpp_bound = BA[2]/N_lo
    sup_hp = certify(hp, h, hpp_bound, "d1")
    OUT["sup_abs_h_prime"] = {"grid_max": sup_hp["grid_max"], "pad": sup_hp["pad"],
                              "upper": sup_hp["upper"], "fix4_grid_value": 11.3248852}
    OUT["h_second_bound"] = {"proved": hpp_bound, "fix4_claimed_PROVED": 225.9466903722,
                             "fix4_is_a_grid_scan": True,
                             "referee_rigorous": 2415.5503891825,
                             "ratio_to_fix4": hpp_bound/225.9466903722}
    print("|h''| proved <= %.10f  (fix4 quoted 225.9466903722, a scan)" % hpp_bound, flush=True)

    # ---- Psi_sigma(0), certified by the separated bound
    #      Psi = sup_{u,phi} e^{-2u} | rad(u) w0(phi) + Theta(u) angL(phi) |
    #          <= sup_u e^{-2u} ( |rad| sup|w0| + |Theta| sup|angL| )
    PHI_LO = 0.02
    m = (phi >= PHI_LO) & (phi <= math.pi - PHI_LO)
    s_, c_ = np.sin(phi[m]), np.cos(phi[m])
    angL = (Wd[2][m] + 3.0*(c_/s_)*Wd[1][m])/N_lo
    BangL1 = (BW[3] + 3.0*(1.0/math.sin(PHI_LO))*BW[2] + 3.0*(1.0/math.sin(PHI_LO)**2)*BW[1])/N_lo
    angL_rec = certify(np.abs(angL), h, BangL1, "d1")
    # the two polar caps: |angL| <= 4 sup_{[0,phi_lo]} |w2|   (W_sigma even at the pole, so
    # |w1(phi)| <= phi sup|w2| and 3 cot(phi)|w1| <= 3 sup|w2| since phi cot phi <= 1)
    cap = np.linspace(0.0, PHI_LO, 20001)
    hc = cap[1] - cap[0]
    w2cap = np.abs(CONV.derivs(cap, kmax=2)[2])/N_lo
    sup_w2_cap = float(np.max(w2cap)) + 0.5*hc*BW[3]/N_lo
    angL_cap = 4.0*sup_w2_cap
    supangL = max(angL_rec["upper"], angL_cap)
    # u is split at u = -1.  For u <= -1 the outer tanh contributes e^{-328} at L = 40 and
    # Theta(u) = sigmoid(8(u - 1/4)) <= e^{8(u-1/4)}, |Theta'| <= 8 Theta, |Theta''| <= 64 Theta,
    # so |rad| <= 74 Theta, and f(u) <= e^{6u - 2}(74 E_0 + sup|angL|), decreasing in -u:
    # the whole half line is bounded by its value at u = -1.  For u >= 8, e^{-2u} <= e^{-16}.
    # On [-1, 8] a uniform cell bound is used, with the exponential at its left endpoint and
    # |rad|, |Theta| at their midpoint values plus (hu/2) times the proved bounds on the next
    # derivative.  No global Lipschitz constant multiplies the exponential's own range.
    E0up = sW_hi/N_lo
    uu = np.linspace(-1.0, 8.0, 4500001)
    hu = uu[1] - uu[0]
    T0, T1, T2 = Theta_derivs_scalar(uu, 40.0)
    rad = T2 + T1 - 2.0*T0
    f_u = (np.exp(-2.0*(uu - 0.5*hu))
           * ((np.abs(rad) + 0.5*hu*TH["rad_prime"])*E0up
              + (np.abs(T0) + 0.5*hu*TH[1])*supangL))
    tail_inner = math.exp(-8.0)*(74.0*E0up + supangL)
    tail_outer = math.exp(-16.0)*(TH["rad"]*E0up + supangL)
    Psi_up = max(float(np.max(f_u)), tail_inner, tail_outer)
    OUT["Psi_sigma_0"] = {"certified_upper": Psi_up, "separated_grid_max": float(np.max(f_u)),
                          "u_cell_width": hu, "tail_bound_u_le_minus1": tail_inner,
                          "tail_bound_u_ge_8": tail_outer,
                          "u_grid_max": float(np.max(f_u)),
                          "sup_abs_angL_upper": supangL,
                          "angL_interior": angL_rec, "angL_cap_bound": angL_cap,
                          "sup_w2_on_cap": sup_w2_cap, "phi_lo": PHI_LO,
                          "argmax_u": float(uu[int(np.argmax(f_u))]),
                          "fix4_value": 489.0279393839,
                          "fix4_below_certified": bool(489.0279393839 <= Psi_up),
                          "note": "separated (hence lossy) certified UPPER bound; fix4's tight "
                                  "scan value is the lower end. Psi enters only eps_a, which is "
                                  "3.9e-05 of eps."}
    print("Psi_sigma(0) certified <= %.6f  (fix4 scan 489.0279393839)" % Psi_up, flush=True)
    json.dump(OUT, open(os.path.join(HERE, "x6_results.json"), "w"), indent=1, default=str)

    # ---------------- r_h and min P_h', re-run through fix4's certificate with the PROVED
    #                  Lipschitz constant in place of the scan value
    prof = P1.Profile(SIGMA, nphi=600001)
    inst = P1.PhInstrument(prof)
    hg = prof.h
    kappa = 0.5*inst.P_h(1.0)
    c_star = LOG32/kappa
    lip_scan = P1.lipschitz_hpp(prof)
    OUT["fix4_grid_control"] = {"lipschitz_hpp_scan": lip_scan,
                                "kappa_delta": kappa, "c_star": c_star,
                                "N_sigma_fix4_grid": prof.N,
                                "N_sigma_in_certified_interval":
                                    bool(N_lo <= prof.N <= N_hi)}
    pad_scan = lip_scan*hg
    pad_proved = hpp_bound*hg
    NSUB = 20000                 # fix4 p1.summarise's own value
    rows = {}
    CAP = 1.1664585076
    windows = {"c = c_*": 1.0, "c = 1.0669124225 c_* (L=2e6)": 1.0669124225,
               "c = 1.0238897079 c_* (L=1e7)": 1.0238897079, "c = cap c_*": CAP}
    for tag, fac in windows.items():
        c = fac*c_star
        lam = math.exp(0.75*c)
        r_scan = P1.certified_r_h(inst, prof, lam, nsub=NSUB, lip_pad=pad_scan)
        P1._RHCACHE.clear()
        r_prov = P1.certified_r_h(inst, prof, lam, nsub=NSUB, lip_pad=pad_proved)
        P1._RHCACHE.clear()
        pp_scan = P1.Ph_prime_lower(prof, lam, r_scan["r_h_certified_lower"], pad_scan)
        pp_prov = P1.Ph_prime_lower(prof, lam, r_prov["r_h_certified_lower"], pad_proved)
        rows[tag] = {"c_over_c_star": fac, "c": c, "lam_max": lam,
                     "r_h_grid_min": r_scan["r_h_grid_min"],
                     "r_h_certified_scan_pad": r_scan["r_h_certified_lower"],
                     "r_h_certified_PROVED_pad": r_prov["r_h_certified_lower"],
                     "min_Ph_prime_scan_pad": pp_scan,
                     "min_Ph_prime_PROVED_pad": pp_prov,
                     "nsub": NSUB, "sup_abs_hprime_scan_pad": r_scan["sup_abs_hprime"],
                     "sup_abs_hprime_PROVED_pad": r_prov["sup_abs_hprime"],
                     "delta_r_h": r_prov["r_h_certified_lower"] - r_scan["r_h_certified_lower"]}
        print("  %-32s lam_max=%.10f  r_h: scan-pad %.10f -> proved-pad %.10f   "
              "min P_h': %.10f -> %.10f"
              % (tag, lam, r_scan["r_h_certified_lower"], r_prov["r_h_certified_lower"],
                 pp_scan, pp_prov), flush=True)
    OUT["r_h_and_Ph_prime"] = rows
    OUT["lip_pads"] = {"grid_h": hg, "pad_scan": pad_scan, "pad_proved": pad_proved,
                       "fix4_quoted_pad": 1.1831e-03, "referee_quoted_pad": 1.2648e-02}

    # ---- the certified window cap with the proved pad: largest c/c_* with min P_h' > 0
    lo, hi = 1.0, 1.30
    for _ in range(80):
        mid = 0.5*(lo + hi)
        lam = math.exp(0.75*mid*c_star)
        rr = P1.certified_r_h(inst, prof, lam, nsub=5000, lip_pad=pad_proved)
        P1._RHCACHE.clear()
        if P1.Ph_prime_lower(prof, lam, rr["r_h_certified_lower"], pad_proved) > 0.0:
            lo = mid
        else:
            hi = mid
    OUT["window_cap_certified_PROVED_pad"] = {"c_over_c_star": lo, "eps_cap": lo - 1.0,
                                              "fix4_cap": CAP,
                                              "lam_max_at_cap": math.exp(0.75*lo*c_star)}
    print("certified window cap with the proved pad: c/c_* = %.10f (fix4 %.10f)"
          % (lo, CAP), flush=True)
    json.dump(OUT, open(os.path.join(HERE, "x6_results.json"), "w"), indent=1, default=str)
    print("wrote x6_results.json")

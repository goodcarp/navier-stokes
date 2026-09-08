"""
q2 -- the datum constants that P1 and P2 need, re-derived here from the definitions.

THE DATUM (THEOREM_S3 sec.1.1: ASSEMBLY sec.1.1 angular profile, hk2 tanh radial ramp),
in units M = rho_0 = 1, with u := log(rho/rho_0), phi_ax := min(phi, pi-phi):

    eta_0 = -(1/rho) F(phi) Theta(u) ,
    F(phi)   = sgn(cos phi) g(phi)/sin phi ,
    g(phi)   = min(1, phi_ax/delta) min(1, |phi - pi/2|/delta_m) ,  delta = 7.5 deg, delta_m = 5 deg
    Theta(u) = (1/2)[ tanh((u - eps_r)/eps_r) - tanh((u - L + eps_r)/eps_r) ] ,  eps_r = 0.25 .

Then, with A(u) := Theta - dTheta/du and B(u) := Theta,

    rho |eta_0|        = |F| B ,
    rho^2 |grad eta_0| = sqrt( F^2 A^2 + F'^2 B^2 ) ,
    |eta_0|            = e^{-u} |F| B .

The joint supremum over (phi,u) of F^2A^2 + F'^2B^2 is linear in (A^2,B^2) at fixed phi and
linear in (F^2,F'^2) at fixed u, so it is attained on the PARETO FRONTIER of the point set
{(A(u)^2, B(u)^2)}.  That frontier has a few dozen points and the double supremum becomes an
exact finite maximisation, not a product grid.  Same device for every other joint sup here.

WHAT THIS SCRIPT ESTABLISHES

  1. E_0  := sup rho|eta_0|/M            and its closed form.
  2. G_0  := ess sup rho^2|grad eta_0|/M and its closed form.
  3. ||eta_0||_inf rho_0/M   (a DIFFERENT functional; the brief's expected "4.09" is this one).
  4. ||grad eta_0||_inf rho_0^2/M.
  5. the tail constants of the decay lemma.
  6. TV(eta_0) = int|grad eta_0|dx, finite; needed only for grad a = K * grad eta (u2 sec.6).
  7. the MOLLIFICATION CONTROL that licenses the approximation step in P1/P2.

Nothing is imported from another seat; cross-seat numbers appear only under FOREIGN.

Outputs -> q2_results.json
"""
import json, math, time, sys
import numpy as np
import mpmath as mp

mp.mp.dps = 30
DEG = mp.pi / 180
DELTA = mp.mpf("7.5") * DEG
DM = mp.mpf(5) * DEG
EPS_R = mp.mpf("0.25")
S_REC = 1 / mp.sin(DELTA)

T0 = time.time()


def tick(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg)); sys.stdout.flush()


OUT = {"datum": {"delta_deg": 7.5, "delta_m_deg": 5.0, "eps_r": 0.25,
                 "form": "eta_0 = -(M/rho) F(phi) Theta(u)"}}

d75 = float(DELTA); dm5 = float(DM); epsr = 0.25
hh = 1e-7


def g_np(phi):
    pax = np.minimum(phi, np.pi - phi)
    return np.minimum(1.0, pax / d75) * np.minimum(1.0, np.abs(phi - np.pi / 2) / dm5)


def F_np(phi):
    return np.sign(np.cos(phi)) * g_np(phi) / np.sin(phi)


def Fp_np(phi):
    """F'(phi), EXACT on each of the six smooth arcs of the piecewise profile.  A finite
       difference is not used: F has kinks at delta, pi/2 - delta_m, pi/2, pi/2 + delta_m,
       pi - delta, and a central difference across a kink returns a chord, while at the two
       endpoints it returns a clipped and meaningless value.  Each arc's derivative is written
       out from the closed form of F on that arc; at a kink the grid carries both one-sided
       values, so the ess-sup is the max over arcs."""
    s_, c_ = np.sin(phi), np.cos(phi)
    out = np.empty_like(phi)
    A = phi < d75                                             # taper, z > 0 side
    B = (phi >= d75) & (phi < np.pi / 2 - dm5)                # plateau, z > 0
    C = (phi >= np.pi / 2 - dm5) & (phi < np.pi / 2)          # equatorial ramp, z > 0
    D = (phi >= np.pi / 2) & (phi < np.pi / 2 + dm5)          # equatorial ramp, z < 0
    E = (phi >= np.pi / 2 + dm5) & (phi < np.pi - d75)        # plateau, z < 0
    G = phi >= np.pi - d75                                    # taper, z < 0 side
    out[A] = (s_[A] - phi[A] * c_[A]) / (d75 * s_[A] ** 2)
    out[B] = -c_[B] / s_[B] ** 2
    out[C] = (-s_[C] - (np.pi / 2 - phi[C]) * c_[C]) / (dm5 * s_[C] ** 2)
    out[D] = -(s_[D] - (phi[D] - np.pi / 2) * c_[D]) / (dm5 * s_[D] ** 2)
    out[E] = c_[E] / s_[E] ** 2
    out[G] = (s_[G] + (np.pi - phi[G]) * c_[G]) / (d75 * s_[G] ** 2)
    return out


def Theta_np(u, L):
    return 0.5 * (np.tanh((u - epsr) / epsr) - np.tanh((u - L + epsr) / epsr))


def dTheta_np(u, L):
    return (np.cosh(np.clip((u - epsr) / epsr, -350, 350)) ** -2
            - np.cosh(np.clip((u - L + epsr) / epsr, -350, 350)) ** -2) / (2 * epsr)


def pareto(P, Q):
    """indices of the vertices of the upper-right convex hull of the pairs (P_i, Q_i).
       max_i (alpha P_i + beta Q_i) over alpha, beta >= 0 is attained at one of them, so this
       turns the joint supremum into an exact finite maximisation over a handful of points."""
    order = np.argsort(-P, kind="stable")
    Qo = Q[order]
    run = np.maximum.accumulate(Qo)
    sel = np.empty(len(Qo), dtype=bool)
    sel[0] = True
    sel[1:] = Qo[1:] > run[:-1]
    keep = [int(i) for i in order[sel]]
    pts = np.column_stack([P[order[sel]], Q[order[sel]]]).astype(float)
    # monotone chain on the frontier (P decreasing, Q increasing) -> upper convex hull
    hull = []
    for j in range(len(keep)):
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = pts[hull[-2]], pts[hull[-1]]
            x3, y3 = pts[j]
            if (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) >= 0:
                hull.pop()
            else:
                break
        hull.append(j)
    return np.array([keep[j] for j in hull], dtype=int)


# ------------------------------------------------------------------ angular grid
# The two end caps phi < PHI_CAP and phi > pi - PHI_CAP are excluded from the grid and
# handled analytically instead: there g = t/delta with t = min(phi, pi-phi), so
#     |F| = (1/delta) t/sin t <= (1/delta)(1 + t^2/6 + ...) ,
#     |F'| = |sin t - t cos t|/(delta sin^2 t) = t/(3 delta) (1 + O(t^2)) ,
# both far below the suprema, and a floating-point central difference or even the closed form
# there suffers catastrophic cancellation (sin t - t cos t ~ t^3/3 against sin^2 t ~ t^2).
PHI_CAP = 1e-3
phis = np.linspace(PHI_CAP, np.pi - PHI_CAP, 300001)
for kk in [d75, np.pi / 2 - dm5, np.pi / 2, np.pi / 2 + dm5, np.pi - d75]:
    phis = np.append(phis, [kk - 1e-11, kk + 1e-11])
phis = np.sort(phis[(phis > 0) & (phis < np.pi)])
Fv = F_np(phis)
Fp = Fp_np(phis)
# control: away from the kinks and the endpoints the analytic derivative must agree with a
# central difference.  This is the only use of a finite difference in this script.
_m = np.ones_like(phis, dtype=bool)
for _k in [d75, np.pi / 2 - dm5, np.pi / 2, np.pi / 2 + dm5, np.pi - d75]:
    _m &= np.abs(phis - _k) > 1e-3
_m &= (phis > 1e-3) & (phis < np.pi - 1e-3)
_fd = (F_np(phis[_m] + hh) - F_np(phis[_m] - hh)) / (2 * hh)
FD_CTRL = float(np.max(np.abs(_fd - Fp[_m]) / (1.0 + np.abs(Fp[_m]))))
FF = Fv ** 2
FP = Fp ** 2

# ------------------------------------------------------------------ radial grid
L_BIG = 1000.0
us = np.concatenate([np.linspace(-4.0, 6.0, 200001), np.array([L_BIG / 2]),
                     np.linspace(L_BIG - 6.0, L_BIG + 4.0, 200001)])
Th = Theta_np(us, L_BIG); dTh = dTheta_np(us, L_BIG)
A2 = (Th - dTh) ** 2
B2 = Th ** 2
idx = pareto(A2, B2)
OUT["pareto"] = {"n_u_grid": int(len(us)), "n_pareto": int(len(idx)),
                 "max_A2": float(A2.max()), "max_B2": float(B2.max())}

tick('grids + pareto done')
# ------------------------------------------------------------------ 1-2. E_0 and G_0
E0_closed = S_REC
G0_closed = S_REC ** 2
OUT["closed_forms"] = {
    "E0 = 1/sin(delta)": float(E0_closed), "E0_str": mp.nstr(E0_closed, 18),
    "G0 = 1/sin(delta)^2 = E0^2": float(G0_closed), "G0_str": mp.nstr(G0_closed, 18),
    "derivation": ("on the plateau arc g == 1 so F = sgn(cos phi)/sin phi and "
                   "F^2 + F'^2 = (sin^2 + cos^2)/sin^4 = 1/sin^4 phi; the arc's closest "
                   "approach to the axis is phi = delta.  The taper arc has |F| <= 1/sin delta "
                   "and |F'| small, the equatorial arc has |F'| <= 1/delta_m."),
    "note": "L -> infinity values; at finite L both carry a factor sup Theta <= 1",
}
supE = float(np.max(np.abs(Fv)) * math.sqrt(B2.max()))
i_argE = int(np.argmax(np.abs(Fv)))
best, argG = -1.0, None
for i in idx:
    val = FF * A2[i] + FP * B2[i]
    m = float(val.max())
    if m > best:
        best = m
        argG = (float(phis[int(val.argmax())] * 180 / np.pi), float(us[i]))
supG = math.sqrt(best)
OUT["E0_G0"] = {
    "E0_scan": supE, "E0_closed": float(E0_closed),
    "E0_rel_gap": abs(supE - float(E0_closed)) / float(E0_closed),
    "E0_argmax_phi_deg": float(phis[i_argE] * 180 / np.pi),
    "G0_scan": supG, "G0_closed": float(G0_closed),
    "G0_rel_gap": abs(supG - float(G0_closed)) / float(G0_closed),
    "G0_argmax_phi_deg_u": argG,
    "sup_abs_Fprime": float(np.max(np.abs(Fp))),
    "analytic_vs_central_difference_control": FD_CTRL,
    "cos_delta_over_sin2_delta": float(mp.cos(DELTA) / mp.sin(DELTA) ** 2),
    "check_E0sq_equals_G0": float(abs(float(E0_closed) ** 2 - float(G0_closed))),
    "endcap_phi": PHI_CAP,
    "endcap_bound_absF": float((1.0 / d75) * (PHI_CAP / math.sin(PHI_CAP))),
    "endcap_bound_absFprime": float(PHI_CAP / (3 * d75) * 1.001),
    "endcap_note": ("on phi < 1e-3 and phi > pi - 1e-3 the analytic bounds |F| <= 7.6395 and "
                    "|F'| <= 2.55e-3 hold, both far below the suprema, so excluding the caps "
                    "from the grid cannot move either supremum."),
}

tick('E0/G0 done')
# ------------------------------------------------------------------ 3-4. sup norms
eu = np.exp(-us) * Th
sup_e_u_Theta = float(np.max(eu)); i_u = int(np.argmax(eu))
OUT["sup_norms"] = {
    "sup_u exp(-u)Theta(u)": sup_e_u_Theta,
    "argmax_u": float(us[i_u]), "argmax_rho_over_rho0": float(np.exp(us[i_u])),
    "eta0_Linf_times_rho0_over_M": sup_e_u_Theta * float(E0_closed),
    "note": ("||eta_0||_inf = (M/rho_0) sup_u e^{-u}Theta * sup_phi|F|.  This is NOT E_0. "
             "The brief's expected 4.09 is this functional; E_0 = sup rho|eta_0|/M = 7.6613."),
}
A2g = np.exp(-2 * us) * A2
B2g = np.exp(-2 * us) * B2
idg = pareto(A2g, B2g)
bg = -1.0
for i in idg:
    bg = max(bg, float((FF * A2g[i] + FP * B2g[i]).max()))
OUT["sup_norms"]["grad_eta0_Linf_times_rho0sq_over_M"] = math.sqrt(bg)

tick('sup norms done')
# ------------------------------------------------------------------ 5. tail constants
def tail_consts(L, m_eta=9, m_grad=10):
    uu = np.linspace(-4.0, L + 8.0, 200001)
    Thl = Theta_np(uu, L); dThl = dTheta_np(uu, L)
    Xi_e = (1.0 + np.exp(2 * (uu - L))) ** (-m_eta / 2.0)
    Xi_g = (1.0 + np.exp(2 * (uu - L))) ** (-m_grad / 2.0)
    rad_e = np.exp(-uu) * Thl
    Ae = float(np.max(np.abs(Fv)) * np.max(rad_e / Xi_e))
    ue = float(uu[int(np.argmax(rad_e / Xi_e))])
    a2 = (np.exp(-2 * uu) * (Thl - dThl)) ** 2 / Xi_g ** 2
    b2 = (np.exp(-2 * uu) * Thl) ** 2 / Xi_g ** 2
    ii = pareto(a2, b2)
    bb = -1.0; ug = None
    for i in ii:
        v = FF * a2[i] + FP * b2[i]
        if float(v.max()) > bb:
            bb = float(v.max()); ug = float(uu[i])
    return Ae, ue, math.sqrt(bb), ug


tails = {}
for L in [10.0, 20.0, 40.0]:
    a9, ua9, b10, ub10 = tail_consts(L)
    tails["L=%g" % L] = {"A_9": a9, "argmax_u": ua9, "B_10": b10, "argmax_u_grad": ub10}
OUT["tail_constants"] = tails
OUT["tail_note"] = ("Theta = O(e^{-2(u-L+eps_r)/eps_r}) = O((rho/R)^{-8}) beyond the outer edge, "
                    "so |eta_0| = O(rho^{-9}) and |grad eta_0| = O(rho^{-10}).  With the "
                    "scale-correct envelope (1+(rho/R)^2)^{-m/2} the constants are O(1) and "
                    "L-independent, and A_9 coincides with ||eta_0||_inf rho_0/M.")

tick('tails done')
# ------------------------------------------------------------------ 6. TV(eta_0)
S3 = 2 * math.pi ** 2
step = 40
phig = phis[::step]
Fg = Fv[::step]; Fpg = Fp[::step]
sin3 = np.sin(phig) ** 3
tv = {}
for L in [10.0, 20.0, 40.0]:
    uu = np.linspace(-4.0, L + 8.0, 20001)
    Thl = Theta_np(uu, L); dThl = dTheta_np(uu, L)
    Al = (Thl - dThl)[:, None]; Bl = Thl[:, None]
    tot = 0.0
    for i0 in range(0, len(uu), 4000):
        sl = slice(i0, min(i0 + 4000, len(uu)))
        val = np.sqrt((Fg[None, :] * Al[sl]) ** 2 + (Fpg[None, :] * Bl[sl]) ** 2) * sin3[None, :]
        inner = np.trapz(val, phig, axis=1) * np.exp(3 * uu[sl])
        tot += np.trapz(inner, uu[sl])
    TV = S3 * tot
    tv["L=%g" % L] = {"TV_over_M_rho0cubed": float(TV),
                      "TV_over_M_Rcubed": float(TV / math.exp(3 * L))}
OUT["TV"] = tv

tick('TV done')
# ------------------------------------------------------------------ 7. mollification control
def mollified_F(sigma, npts=100001):
    ph = np.linspace(0.0, np.pi, npts)
    dp = ph[1] - ph[0]
    Fv0 = F_np(np.clip(ph, 1e-12, np.pi - 1e-12))
    Fv0[0] = Fv0[1]; Fv0[-1] = Fv0[-2]
    half = int(math.ceil(sigma / dp))
    t = np.arange(-half, half + 1) * dp / sigma
    ker = np.where(np.abs(t) < 1.0, np.exp(-1.0 / np.maximum(1.0 - t ** 2, 1e-300)), 0.0)
    ker = ker / ker.sum()
    ext = np.concatenate([Fv0[half:0:-1], Fv0, Fv0[-2:-half - 2:-1]])
    Fm = np.convolve(ext, ker, mode="same")[half:half + npts]
    return ph, Fm, np.gradient(Fm, dp)


moll = {}
for sigma in [0.02, 0.01, 0.005, 0.002]:
    ph, Fm, Fmp = mollified_F(sigma)
    e0m = float(np.max(np.abs(Fm)))
    FFm = Fm ** 2; FPm = Fmp ** 2
    bm = -1.0
    for i in idx:
        bm = max(bm, float((FFm * A2[i] + FPm * B2[i]).max()))
    moll["sigma=%g" % sigma] = {
        "E0_mollified": e0m, "E0_ratio_to_E0": e0m / float(E0_closed),
        "G0_mollified": math.sqrt(bm), "G0_ratio_to_G0": math.sqrt(bm) / float(G0_closed),
        "sup_abs_Fprime_mollified": float(np.max(np.abs(Fmp))),
    }
OUT["mollification_control"] = moll
OUT["mollification_claim"] = ("F_sigma = F * chi_sigma has ||F_sigma||_inf <= ||F||_inf and "
                              "||F_sigma'||_inf = ||F' * chi_sigma||_inf <= ||F'||_inf, so "
                              "neither E_0 nor G_0 increases under angular mollification.  The "
                              "table confirms it numerically at four scales.")

tick('mollification done')
# ------------------------------------------------------------------ FOREIGN comparisons
OUT["FOREIGN"] = {
    "THEOREM_S3_t1_E0": {"value": 7.6612975721936944, "origin": "round2/THEOREM_S3/t1_results.json"},
    "THEOREM_S3_t1_G0": {"value": 58.695476310955364, "origin": "round2/THEOREM_S3/t1_results.json"},
    "refute_u2_G0_ASSEMBLY_family": {"value": 58.6948961, "origin": "round2/refute-u2/NOTE.md sec.3"},
    "refute_u2_eta0_Linf_for_DC": {"value": 4.093173, "origin": "round2/refute-u2/NOTE.md MINOR-2"},
    "refute_u2_sup_e_u_Theta_DC": {"value": 0.53431477, "origin": "round2/refute-u2/NOTE.md MINOR-2"},
}
OUT["comparison"] = {
    "E0_here_vs_t1_rel": abs(supE - 7.6612975721936944) / 7.6612975721936944,
    "G0_here_vs_t1_rel": abs(supG - 58.695476310955364) / 58.695476310955364,
    "G0_closed_vs_refute_u2_rel": abs(float(G0_closed) - 58.6948961) / float(G0_closed),
    "eta0_Linf_here_vs_refute_u2_DC_rel":
        abs(OUT["sup_norms"]["eta0_Linf_times_rho0_over_M"] - 4.093173) / 4.093173,
    "sup_e_u_Theta_here_vs_refute_u2_rel":
        abs(sup_e_u_Theta - 0.53431477) / 0.53431477,
}

for k in ["closed_forms", "E0_G0", "sup_norms", "tail_constants", "TV",
          "mollification_control", "comparison"]:
    print("==", k, "==")
    print(json.dumps(OUT[k], indent=1, default=str))
json.dump(OUT, open("q2_results.json", "w"), indent=1, sort_keys=True, default=str)
print("\nwrote q2_results.json")

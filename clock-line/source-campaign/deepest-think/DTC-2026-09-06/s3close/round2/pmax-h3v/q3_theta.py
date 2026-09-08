"""
q3 -- vartheta for the theorem's datum on the ball B(x_*, d):  TASK B's datum input.

Theorem V.4's hypothesis (H3) asks for  eta_0 == eta_P := -M sgn(z)/r  EXACTLY on B(x_*,d).
V.4' replaces it by the quantitative

    (H3'_vartheta)   sup_{|v|=1} | d_v^k ( eta_0 - eta_P ) (x) |  <=  vartheta k! M / r(x)^{k+1}
                     for k = 0,1,2,3,4 and every x in B(x_*, d).

k = 0 is the brief's  |eta_0 + M sgn(z)/r| <= vartheta M/r; the k = 2 clause implies the brief's
Laplacian clause and is what Step 2's Hessian contraction actually needs (see PROOF.md sec.B2:
the Laplacian alone does NOT suffice, because C_tau is anisotropic).  k = 1 is used by Step 1,
k = 3 and 4 by Step 2's Taylor polynomial and its remainder.  The normalisation k! M/r^{k+1} is
the exact value of sup|d_v^k eta_P| (V-b identity (2.4)), so vartheta is a pure relative
deviation, scale free.

WHAT eta_0 - eta_P IS, ON THE BALL.  d is chosen (ASSEMBLY sec.2.5) so that B(x_*,d) misses
both angular kinks, hence g == 1 there and eta_0 = eta_P * Theta(u).  Therefore, on B,

    e := eta_0 - eta_P = - eta_P * psi(rho) = M psi(rho)/r    (z > 0),   psi := 1 - Theta ,

and with eps_r = 1/4, rho_0 = 1,  psi = psi_in + psi_out,
    psi_in(rho)  = 1/(1 + e^{2(u-eps_r)/eps_r}) = e^2/(e^2 + rho^8) ,
    psi_out(rho) = 1/(1 + e^{-2(u-L+eps_r)/eps_r}) <= e^{-8(L - eps_r - u)} ,
so on the ball e is the RATIONAL function e^2/(r(e^2 + rho^8)) up to psi_out, which underflows
at every L the theorem uses and is reported separately.

HOW THE DERIVATIVES ARE TAKEN.  Along a line x + s v with |v| = 1,
    R(s)^2 = r^2 + 2 r a s + q s^2 ,  q := a^2 + b^2 ,
    P(s)^2 = rho^2 + 2 sigma s + s^2 ,  sigma := x.v ,
both quadratics, so
    e(s) = e^2 * (R^2)^{-1/2} * (e^2 + (P^2)^4)^{-1}
and its Taylor coefficients at s = 0 are obtained EXACTLY by formal power series: the
J.C.P. Miller recurrence for a power of a series, polynomial powers for (P^2)^4, the
reciprocal recurrence for the second factor, and one convolution.  d_v^k e(0) = k! c_k.
No difference quotient is used anywhere.  A sympy cross-check on the first two orders is run
at the end.

THE GEOMETRY.  x_* = ((1+f) sin phi_0, 0,0,0, (1+f) cos phi_0) in units rho_0 = 1,
phi_0 = 30 deg; d = min(f, rho_* sin(phi_0 - delta), rho_* sin(pi/2 - delta_m - phi_0)) is
ASSEMBLY sec.2.5's rule as coded in THEOREM_S3/t3_budget.py.  By the SO(4) symmetry a point of
the ball may be taken as x_* + (w1, w2, 0, 0, w5) and a direction as v = a Y-hat + b e_perp +
c e_z, so the search is 3 + 2 dimensional; it is done on a grid whose refinement stability is
reported.

Outputs -> q3_results.json
"""
import json, math
import numpy as np

DEG = math.pi / 180.0
DELTA = 7.5 * DEG
DM = 5.0 * DEG
PHI0 = 30.0 * DEG
EPS_R = 0.25
E2 = math.exp(2.0)
KMAX = 4
FACT = [1.0, 1.0, 2.0, 6.0, 24.0]

OUT = {"setting": {"delta_deg": 7.5, "delta_m_deg": 5.0, "phi0_deg": 30.0, "eps_r": 0.25,
                   "units": "M = rho_0 = 1", "kmax": KMAX}}


def series_pow(g, alpha, n):
    """Taylor coefficients of (sum g_j s^j)^alpha to order n, J.C.P. Miller recurrence.
       g is a list of arrays (g[0] > 0)."""
    a = [g[0] ** alpha]
    for m in range(1, n + 1):
        acc = np.zeros_like(g[0])
        for j in range(1, min(m, len(g) - 1) + 1):
            acc = acc + (j * alpha - (m - j)) * g[j] * a[m - j]
        a.append(acc / (m * g[0]))
    return a


def series_mul(u, v, n):
    """truncated product of two power series, shorter lists padded with zeros."""
    zu, zv = np.zeros_like(u[0]), np.zeros_like(v[0])
    U = list(u) + [zu] * (n + 1 - len(u))
    V = list(v) + [zv] * (n + 1 - len(v))
    return [sum(U[j] * V[m - j] for j in range(0, m + 1)) for m in range(n + 1)]


def series_inv(d, n):
    b = [1.0 / d[0]]
    for m in range(1, n + 1):
        acc = np.zeros_like(d[0])
        for j in range(1, min(m, len(d) - 1) + 1):
            acc = acc + d[j] * b[m - j]
        b.append(-acc / d[0])
    return b


def dv_e(r, z, a, b, c, n=KMAX):
    """d_v^k e at s = 0 for k = 0..n, vectorised over arrays r, z (scalars a, b, c)."""
    rho2 = r * r + z * z
    sigma = r * a + z * c                       # x.v
    q = a * a + b * b
    g = [r * r, 2.0 * r * a, np.full_like(r, q)]              # R^2
    A = series_pow(g, -0.5, n)                                # (R^2)^{-1/2}
    p = [rho2, np.full_like(r, 2.0 * sigma), np.ones_like(r)]  # P^2
    p2 = series_mul(p, p, n)
    p4 = series_mul(p2, p2, n)
    d = list(p4)
    d[0] = d[0] + E2
    B = series_inv(d, n)
    C = series_mul(A, B, n)
    return [E2 * FACT[k] * C[k] for k in range(n + 1)]


def psi_in(rho):
    return E2 / (E2 + rho ** 8)


def psi_out_bound(u_max, L):
    return math.exp(-8.0 * (L - EPS_R - u_max))


def sphere_grid(n_th, n_ph):
    out = []
    for i in range(n_th):
        ct = -1.0 + 2.0 * (i + 0.5) / n_th
        st = math.sqrt(max(0.0, 1.0 - ct * ct))
        for j in range(n_ph):
            ph = 2 * math.pi * j / n_ph
            out.append((st * math.cos(ph), st * math.sin(ph), ct))
    return np.array(out)


def polish(f, d, best, rounds=6, n_local=9):
    """local refinement of each sup, from the grid argmax, by shrinking vectorised local grids.
       A grid search returns a LOWER bound on a supremum; this tightens it.  The gain over the
       coarse grid is reported, so the reader can see how much the grid was missing."""
    rho_s = 1.0 + f
    r_s = rho_s * math.sin(PHI0)
    z_s = rho_s * math.cos(PHI0)
    out = {}
    for k in range(KMAX + 1):
        val0, arg = best[k]
        # recover (w-radius, w-direction, v-direction) is not needed: refine directly in
        # (r, z) around the argmax point, subject to the ball constraint, and in (a,b,c).
        r0, z0, a0, b0, c0 = arg
        cur = (val0, r0, z0, a0, b0, c0)
        span_x, span_v = 0.25 * d, 0.6
        for _ in range(rounds):
            rr = np.linspace(cur[1] - span_x, cur[1] + span_x, n_local)
            zz = np.linspace(cur[2] - span_x, cur[2] + span_x, n_local)
            RR, ZZ = np.meshgrid(rr, zz, indexing="ij")
            RR = RR.ravel(); ZZ = ZZ.ravel()
            ok = ((RR - r_s) ** 2 + (ZZ - z_s) ** 2 <= d * d + 1e-15) & (RR > 0)
            if not ok.any():
                break
            RR, ZZ = RR[ok], ZZ[ok]
            th0 = math.atan2(math.hypot(cur[3], cur[4]), cur[5])
            ph0 = math.atan2(cur[4], cur[3])
            ths = np.linspace(th0 - span_v, th0 + span_v, n_local)
            phs = np.linspace(ph0 - span_v, ph0 + span_v, n_local)
            for th in ths:
                for ph in phs:
                    a = math.sin(th) * math.cos(ph); b = math.sin(th) * math.sin(ph)
                    c = math.cos(th)
                    D = dv_e(RR, ZZ, a, b, c)
                    v = np.abs(D[k]) * RR ** (k + 1) / FACT[k]
                    i = int(np.argmax(v))
                    if v[i] > cur[0]:
                        cur = (float(v[i]), float(RR[i]), float(ZZ[i]), a, b, c)
            span_x *= 0.45; span_v *= 0.45
        out[k] = cur[0]
    return out


def theta_of(f, d, n_rad, n_dir, n_v):
    rho_s = 1.0 + f
    r_s = rho_s * math.sin(PHI0)
    z_s = rho_s * math.cos(PHI0)
    Wd = sphere_grid(*n_dir)
    Vd = sphere_grid(*n_v)
    pts = [(r_s, z_s)]
    for rr in np.linspace(0.0, d, n_rad)[1:]:
        for (u1, u2, u3) in Wd:
            pts.append((math.sqrt((r_s + rr * u1) ** 2 + (rr * u2) ** 2), z_s + rr * u3))
    pts = np.array(pts)
    rv, zv = pts[:, 0].copy(), pts[:, 1].copy()
    best = [(-1.0, None)] * (KMAX + 1)
    for (aa, bb, cc) in Vd:
        D = dv_e(rv, zv, aa, bb, cc)
        for k in range(KMAX + 1):
            val = np.abs(D[k]) * rv ** (k + 1) / FACT[k]
            i = int(np.argmax(val))
            if val[i] > best[k][0]:
                best[k] = (float(val[i]), (float(rv[i]), float(zv[i]),
                                           float(aa), float(bb), float(cc)))
    return best, (r_s, z_s, rho_s)


def theta_for(f, d_scale=1.0, cap_inner=True, fine=False):
    rho_s = 1.0 + f
    d_taper = rho_s * math.sin(PHI0 - DELTA)
    d_eq = rho_s * math.sin(math.pi / 2 - DM - PHI0)
    d_max = min(d_taper, d_eq)
    if cap_inner:
        d_max = min(d_max, f)
    d = d_scale * d_max
    if d <= 0:
        return None
    nr = 15 if fine else 9
    nd = (30, 40) if fine else (16, 22)
    nv = (44, 58) if fine else (24, 32)
    best, geo = theta_of(f, d, nr, nd, nv)
    r_s, z_s, rho_s = geo
    rho_min = rho_s - d
    grid_vals = {k: best[k][0] for k in range(KMAX + 1)}
    pol = polish(f, d, best)
    pol[0] = max(pol[0], psi_in(max(rho_min, 1e-12)))     # exact: vartheta_0 = psi(rho_min)
    return {
        "f": f, "d_scale": d_scale, "d": d, "d_max_rule": d_max,
        "d_taper": d_taper, "d_eq": d_eq, "d_inner_cap": (f if cap_inner else None),
        "rho_star": rho_s, "r_star": r_s, "z_star": z_s,
        "rho_min_on_ball": rho_min, "u_min_on_ball": math.log(max(rho_min, 1e-12)),
        "Theta_at_rho_star": 1.0 - psi_in(rho_s),
        "u_star": math.log(rho_s),
        "R_minus": r_s - d, "psi_at_rho_min": psi_in(max(rho_min, 1e-12)),
        "vartheta_k": pol,
        "vartheta_k_grid_only": grid_vals,
        "grid_to_polish_gain": {k: (pol[k] / max(grid_vals[k], 1e-300)) for k in range(KMAX + 1)},
        "argmax_k_grid": {k: best[k][1] for k in range(KMAX + 1)},
        "vartheta": max(pol.values()),
        "argmax_k_of_max": int(max(range(KMAX + 1), key=lambda k: pol[k])),
    }


# ------------------------------------------------------------------ sympy cross-check
import sympy as sp
_s, _r, _z, _a, _b, _c = sp.symbols("s r z a b c", real=True)
_R = sp.sqrt((_r + _s * _a) ** 2 + (_s * _b) ** 2)
_P2 = _R ** 2 + (_z + _s * _c) ** 2
_e = sp.exp(2) / (_R * (sp.exp(2) + _P2 ** 4))
_chk = {}
_pt = {_r: 1.3, _z: 2.1, _a: 0.3, _b: 0.5, _c: math.sqrt(1 - 0.09 - 0.25)}
for k in range(KMAX + 1):
    sym = float(sp.diff(_e, _s, k).subs(_pt).subs(_s, 0))
    num = float(dv_e(np.array([1.3]), np.array([2.1]), 0.3, 0.5,
                     math.sqrt(1 - 0.09 - 0.25))[k][0])
    _chk["k=%d" % k] = {"sympy": sym, "series": num,
                        "rel": abs(sym - num) / max(abs(sym), 1e-300)}
OUT["series_vs_sympy_control"] = _chk
print("series vs sympy control:", json.dumps(_chk, indent=1))

# ------------------------------------------------------------------ the brief's three f
OUT["brief_f"] = {}
for f in [0.5, 1.0, 2.0]:
    OUT["brief_f"]["f=%g" % f] = theta_for(f, 1.0, cap_inner=True, fine=True)

FS = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
OUT["f_grid_d_from_assembly_rule"] = {}
for f in FS:
    OUT["f_grid_d_from_assembly_rule"]["f=%g" % f] = theta_for(f, 1.0, cap_inner=True)

OUT["f_grid_d_scaled"] = {}
for f in FS:
    for ds in [1.0, 0.5, 0.25]:
        OUT["f_grid_d_scaled"]["f=%g,d_scale=%g" % (f, ds)] = theta_for(f, ds, cap_inner=True)

OUT["f_grid_no_inner_cap"] = {}
for f in FS:
    for ds in [1.0, 0.5]:
        OUT["f_grid_no_inner_cap"]["f=%g,d_scale=%g" % (f, ds)] = theta_for(f, ds,
                                                                           cap_inner=False)

ctrl = {}
for f in [1.0, 8.0]:
    coarse = theta_for(f, 1.0, cap_inner=True)
    fineR = theta_for(f, 1.0, cap_inner=True, fine=True)
    ctrl["f=%g" % f] = {"coarse": coarse["vartheta"], "fine": fineR["vartheta"],
                        "rel_change": abs(fineR["vartheta"] - coarse["vartheta"])
                        / max(fineR["vartheta"], 1e-300)}
OUT["grid_stability"] = ctrl

chk0 = {}
for f in FS:
    row = OUT["f_grid_d_from_assembly_rule"]["f=%g" % f]
    chk0["f=%g" % f] = {"vartheta_0": row["vartheta_k"][0], "psi(rho_min)": row["psi_at_rho_min"],
                        "rel": abs(row["vartheta_k"][0] - row["psi_at_rho_min"])
                        / max(row["psi_at_rho_min"], 1e-300)}
OUT["control_vartheta0_equals_psi"] = chk0
OUT["psi_out_bound"] = {"L=1e4,u_max=4": psi_out_bound(4.0, 1e4),
                        "L=3e5,u_max=4": psi_out_bound(4.0, 3e5),
                        "note": "psi_out underflows at every L the theorem uses"}
# Theta at the tracked point itself: nothing in ell_loss(f) forces this near 1.
OUT["Theta_at_tracked_point"] = {"f=%g" % f: {"rho_star": 1.0 + f, "u_star": math.log(1.0 + f),
                                              "Theta": 1.0 - psi_in(1.0 + f)}
                                 for f in FS}
OUT["asymptotic"] = {
    "psi_in": "e^2/(e^2 + rho^8) ~ e^2 rho^{-8} deep in the plateau",
    "rho_min_needed_for_psi_below": {str(p): math.exp(0.25) * (1.0 / p) ** 0.125
                                     for p in [1e-2, 1e-3, 1e-4, 1e-6]},
}

print("\n   f  |    d      rho_min   psi(rho_min)    vartheta    binding k")
for f in FS:
    row = OUT["f_grid_d_from_assembly_rule"]["f=%g" % f]
    print("%6g | %8.4f %8.4f  %11.4e  %11.4e   k=%d"
          % (f, row["d"], row["rho_min_on_ball"], row["psi_at_rho_min"],
             row["vartheta"], row["argmax_k_of_max"]))
print("\nno inner cap d <= f (tanh datum has no inner kink):")
print("   f  |    d      rho_min    vartheta(d_max)   vartheta(d_max/2)")
for f in FS:
    r1 = OUT["f_grid_no_inner_cap"]["f=%g,d_scale=1" % f]
    r2 = OUT["f_grid_no_inner_cap"]["f=%g,d_scale=0.5" % f]
    print("%6g | %8.4f %8.4f   %12.4e   %12.4e"
          % (f, r1["d"], r1["rho_min_on_ball"], r1["vartheta"], r2["vartheta"]))
print("\ngrid stability:", json.dumps(ctrl, indent=1))
print("vartheta_0 vs psi control (max rel):",
      max(v["rel"] for v in chk0.values()))
json.dump(OUT, open("q3_results.json", "w"), indent=1, sort_keys=True, default=str)
print("wrote q3_results.json")

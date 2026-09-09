"""
x1 -- UNIT F-1 (FATAL).  (H3'_vartheta) recomputed for the sigma-mollified datum on a ball that
stands off the taper cone.

WHAT THE REFEREE FOUND.  THEOREM_S3_v2 row B4 imports vartheta(f=4) = 0.055612 from
round2/pmax-h3v Part B.  That number is a statement about the KINKED angular profile, on a ball
where the angular factor is IDENTICALLY 1, so eta_0 - eta_P is purely the radial tanh ramp.
fix4 replaced the kinked angle by its Gaussian mollification, and at f >= 4 ASSEMBLY's rule
d = rho_* sin(phi_0 - delta) makes the ball EXACTLY TANGENT to the cone phi = delta, so the
smoothing layer |phi - delta| <~ 3 sigma lies inside it.  There the angular derivatives are
O(sigma^{-(k-1)}) and (H3'_vartheta) fails: vartheta >= 8.76 at every admissible f.

THE REPAIR, MEASURED HERE.  Move the ball off the cone:

    d  :=  rho_* min( sin(phi_0 - delta - 4 sigma),
                      sin(pi/2 - delta_m - phi_0 - 4 sigma),
                      f/rho_* ) ,        sigma = 0.02, rho_* = (1+f) rho_0 .

The first branch binds at every f >= 4, and it puts the closest point of the ball at polar angle
exactly delta + 4 sigma, i.e. four mollification widths off the taper cone.

THE INSTRUMENT.  pmax-h3v sec.B3's five-parameter search, re-implemented for the MOLLIFIED
datum: three parameters for the point (a radius in (0,d] and a direction on S^2, using the SO(4)
symmetry to write a ball point as x_* + (w1, w2, 0, 0, w5)), two for the unit direction
v = a Yhat + b e_perp + c e_z.  Along x + s v everything is composed by truncated Taylor
arithmetic to order 4:

    e := eta_0 - eta_P  =  1/r  -  Theta(u) W_sigma(phi) / (rho N_sigma)      (M = rho_0 = 1, z > 0)

with r(s), z(s), rho(s), u(s) = log rho, phi(s) = atan(r/z) all series in s, W_sigma^{(k)}(phi_0)
from fx_profile's Gauss-Legendre convolution and Theta^{(k)}(u_0) in closed form.  For the
KINKED control W_sigma is replaced by the closed-form 1/sin phi and its four derivatives and
N_sigma by 1; the instrument must then return pmax-h3v's 0.055612 at f = 4 on ASSEMBLY's d.

    vartheta_k := sup_{ball, |v|=1} r^{k+1} |d_v^k e| / (k! M) ,   vartheta := max_{k<=4} vartheta_k .

A grid search returns a LOWER bound on a supremum, so the theorem is stated at
SAFETY * vartheta with SAFETY = 2, which is pmax-h3v sec.B5's own convention.

Outputs -> x1_results.json
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fx_profile as FX

DEG = FX.DEG
DELTA, DM, PHI0 = FX.DELTA, FX.DM, FX.PHI0
EPS_R = FX.EPS_R
SIGMA = 0.02
NORD = 5
FACT = [1.0, 1.0, 2.0, 6.0, 24.0]
CONV = FX.Conv(SIGMA)

# ------------------------------------------------------------------ truncated Taylor algebra
def tmul(a, b):
    o = [0.0]*NORD
    for i in range(NORD):
        for j in range(NORD - i):
            o[i+j] = o[i+j] + a[i]*b[j]
    return o

def tcompose(fd, s):
    """f(s(.)) with fd[k] = f^{(k)}(s[0]) (arrays), s a series with s[0] = the base point."""
    du = [np.zeros_like(s[0])] + list(s[1:])
    out = [fd[0]] + [np.zeros_like(s[0]) for _ in range(NORD-1)]
    pw = [np.ones_like(s[0])] + [np.zeros_like(s[0]) for _ in range(NORD-1)]
    fac = 1.0
    for k in range(1, NORD):
        pw = tmul(pw, du)
        fac *= k
        fk = fd[k]/fac
        for i in range(NORD):
            out[i] = out[i] + fk*pw[i]
    return out

def tsqrt(s):
    a = s[0]
    return tcompose([a**0.5, 0.5*a**-0.5, -0.25*a**-1.5, 0.375*a**-2.5, -0.9375*a**-3.5], s)

def tlog(s):
    a = s[0]
    return tcompose([np.log(a), 1.0/a, -1.0/a**2, 2.0/a**3, -6.0/a**4], s)

def tinv(s):
    a = s[0]
    return tcompose([1.0/a, -1.0/a**2, 2.0/a**3, -6.0/a**4, 24.0/a**5], s)

def tatan(s):
    a = s[0]
    d1 = 1.0/(1.0 + a*a)
    return tcompose([np.arctan(a), d1, -2.0*a*d1**2, (-2.0 + 6.0*a*a)*d1**3,
                     (24.0*a - 24.0*a**3)*d1**4], s)

def Theta_derivs(u, L, er=EPS_R):
    """Theta(u) = (1/2)[tanh((u-er)/er) - tanh((u-L+er)/er)] and d^k/du^k, k = 0..4."""
    tot = [np.zeros_like(u) for _ in range(NORD)]
    for arg, sgn in (((u - er)/er, 1.0), ((u - L + er)/er, -1.0)):
        T = np.tanh(np.clip(arg, -350.0, 350.0))
        d = [T, 1.0 - T*T, -2.0*T*(1.0 - T*T), -2.0*(1.0 - T*T)*(1.0 - 3.0*T*T),
             8.0*T*(1.0 - T*T)*(2.0 - 3.0*T*T)]
        for k in range(NORD):
            tot[k] = tot[k] + 0.5*sgn*d[k]/er**k
    return tot

def W_kinked_derivs(phi):
    s, c = np.sin(phi), np.cos(phi)
    return [1.0/s, -c/s**2, (1.0 + c*c)/s**3, -c*(5.0 + c*c)/s**4,
            (5.0 + 18.0*c*c + c**4)/s**5]

# ------------------------------------------------------------------ e and its derivatives
def theta_vals(r, z, Wd, Td, a, b, c, N):
    """r^{k+1}|d_v^k e|/k! at s = 0, k = 0..4, vectorised over the point arrays r, z."""
    q = a*a + b*b
    r2 = [r*r, 2.0*r*a*np.ones_like(r), q*np.ones_like(r),
          np.zeros_like(r), np.zeros_like(r)]
    rs = tsqrt(r2)
    zs = [z, c*np.ones_like(r), np.zeros_like(r), np.zeros_like(r), np.zeros_like(r)]
    rho2 = [r*r + z*z, 2.0*(r*a + z*c), np.ones_like(r), np.zeros_like(r), np.zeros_like(r)]
    rho = tsqrt(rho2)
    u = tlog(rho)
    phi = tatan(tmul(rs, tinv(zs)))
    Th = tcompose(Td, u)
    Wc = tcompose(Wd, phi)
    term = tmul(Th, tmul(Wc, tinv(rho)))
    ri = tinv(rs)
    e = [ri[k] - term[k]/N for k in range(NORD)]
    return [r**(k+1)*np.abs(e[k]) for k in range(NORD)]

def W_and_Theta(r, z, L, kinked):
    phi = np.arctan2(r, z)
    u = np.log(np.hypot(r, z))
    Wd = W_kinked_derivs(phi) if kinked else CONV.derivs(phi)
    return Wd, Theta_derivs(u, L), phi

def sphere_grid(n_th, n_ph):
    out = []
    for i in range(n_th):
        ct = -1.0 + 2.0*(i + 0.5)/n_th
        st = math.sqrt(max(0.0, 1.0 - ct*ct))
        for j in range(n_ph):
            p = 2.0*math.pi*j/n_ph
            out.append((st*math.cos(p), st*math.sin(p), ct))
    return np.array(out)

def d_rule(f, offset):
    """d = rho_* min(sin(phi_0 - delta - off), sin(pi/2 - delta_m - phi_0 - off), f/rho_*)."""
    rho_s = 1.0 + f
    d_taper = rho_s*math.sin(PHI0 - DELTA - offset)
    d_eq = rho_s*math.sin(math.pi/2 - DM - PHI0 - offset)
    d_inner = f
    d = min(d_taper, d_eq, d_inner)
    which = ["taper", "equatorial", "inner"][int(np.argmin([d_taper, d_eq, d_inner]))]
    return d, {"d_taper": d_taper, "d_equatorial": d_eq, "d_inner": d_inner, "binding": which}

def scan(f, d, L, N, kinked, n_rad=9, n_dir=(16, 22), n_v=(24, 32)):
    rho_s = 1.0 + f
    r_s, z_s = rho_s*math.sin(PHI0), rho_s*math.cos(PHI0)
    Wdir = sphere_grid(*n_dir)
    Vdir = sphere_grid(*n_v)
    pts = [(r_s, z_s)]
    for rr in np.linspace(0.0, d, n_rad)[1:]:
        for (u1, u2, u3) in Wdir:
            pts.append((math.sqrt((r_s + rr*u1)**2 + (rr*u2)**2), z_s + rr*u3))
    P = np.array(pts)
    r, z = P[:, 0].copy(), P[:, 1].copy()
    Wd, Td, phi = W_and_Theta(r, z, L, kinked)
    best = [(-1.0, None)]*NORD
    for (a, b, c) in Vdir:
        V = theta_vals(r, z, Wd, Td, a, b, c, N)
        for k in range(NORD):
            i = int(np.argmax(V[k]))
            if V[k][i] > best[k][0]:
                best[k] = (float(V[k][i]), (float(r[i]), float(z[i]), float(a), float(b),
                                            float(c), float(phi[i])))
    return best, (r_s, z_s, rho_s)

def polish(f, d, L, N, kinked, best, rounds=7, n_local=9):
    rho_s = 1.0 + f
    r_s, z_s = rho_s*math.sin(PHI0), rho_s*math.cos(PHI0)
    out = {}
    arg = {}
    for k in range(NORD):
        v0, A = best[k]
        cur = (v0, A[0], A[1], A[2], A[3], A[4])
        sx, sv = 0.25*d, 0.6
        for _ in range(rounds):
            rr = np.linspace(cur[1] - sx, cur[1] + sx, n_local)
            zz = np.linspace(cur[2] - sx, cur[2] + sx, n_local)
            RR, ZZ = np.meshgrid(rr, zz, indexing="ij")
            RR, ZZ = RR.ravel(), ZZ.ravel()
            ok = ((RR - r_s)**2 + (ZZ - z_s)**2 <= d*d + 1e-15) & (RR > 1e-6) & (ZZ > 1e-6)
            if not ok.any():
                break
            RR, ZZ = RR[ok], ZZ[ok]
            Wd, Td, ph = W_and_Theta(RR, ZZ, L, kinked)
            th0 = math.atan2(math.hypot(cur[3], cur[4]), cur[5])
            ph0 = math.atan2(cur[4], cur[3])
            for th in np.linspace(th0 - sv, th0 + sv, n_local):
                for pp in np.linspace(ph0 - sv, ph0 + sv, n_local):
                    a = math.sin(th)*math.cos(pp)
                    b = math.sin(th)*math.sin(pp)
                    c = math.cos(th)
                    V = theta_vals(RR, ZZ, Wd, Td, a, b, c, N)[k]
                    i = int(np.argmax(V))
                    if V[i] > cur[0]:
                        cur = (float(V[i]), float(RR[i]), float(ZZ[i]), a, b, c)
            sx *= 0.45
            sv *= 0.45
        out[k] = cur[0]
        arg[k] = {"r": cur[1], "z": cur[2], "phi_deg": math.atan2(cur[1], cur[2])/DEG,
                  "rho": math.hypot(cur[1], cur[2]),
                  "phi_minus_delta_over_sigma": (math.atan2(cur[1], cur[2]) - DELTA)/SIGMA,
                  "v": [cur[3], cur[4], cur[5]]}
    return out, arg

def run(f, offset, L, N, kinked, tag, fine=False, d_override=None):
    d, dinfo = d_rule(f, offset)
    if d_override is not None:
        d = d_override
    nr, nd, nv = (15, (30, 40), (44, 58)) if fine else (9, (16, 22), (24, 32))
    best, geo = scan(f, d, L, N, kinked, nr, nd, nv)
    grid = {k: best[k][0] for k in range(NORD)}
    pol, arg = polish(f, d, L, N, kinked, best)
    rho_s = geo[2]
    kb = int(max(range(NORD), key=lambda k: pol[k]))
    return {"tag": tag, "f": f, "offset_rad": offset, "offset_over_sigma": offset/SIGMA,
            "d": d, "d_branches": dinfo, "kinked": kinked, "N_sigma": N, "L": L,
            "rho_star": rho_s, "r_star": rho_s*math.sin(PHI0), "R_minus": rho_s*math.sin(PHI0) - d,
            "rho_min_on_ball": rho_s - d, "rho_max_on_ball": rho_s + d,
            "phi_min_on_ball_deg": (PHI0 - math.asin(min(d/rho_s, 1.0)))/DEG,
            "phi_max_on_ball_deg": (PHI0 + math.asin(min(d/rho_s, 1.0)))/DEG,
            "vartheta_k_grid": grid, "vartheta_k": pol, "argmax": arg,
            "grid_to_polish_gain": {k: pol[k]/max(grid[k], 1e-300) for k in range(NORD)},
            "vartheta": max(pol.values()), "binding_k": kb,
            "binding_phi_deg": arg[kb]["phi_deg"],
            "binding_phi_minus_delta_over_sigma": arg[kb]["phi_minus_delta_over_sigma"],
            "SAFETY": 2.0, "safety_vartheta": 2.0*max(pol.values()),
            "admissible_at_SAFETY_2": bool(2.0*max(pol.values()) < 1.0)}

if __name__ == "__main__":
    L = 40.0
    N = 1.0124508488                       # fix4 sec.1.2; x6 certifies it
    OUT = {"definition": "vartheta_k = sup_{B(x_*,d) x S^4} r^{k+1}|d_v^k(eta_0-eta_P)|/(k! M); "
                         "five-parameter search (3 point + 2 direction) with local polish; a "
                         "grid search is a LOWER bound on a supremum, so the theorem is stated "
                         "at SAFETY*vartheta with SAFETY = 2 (pmax-h3v sec.B5)",
           "sigma": SIGMA, "N_sigma_used": N, "L_used": L,
           "delta_deg": 7.5, "delta_m_deg": 5.0, "phi0_deg": 30.0, "eps_r": EPS_R}

    # ---- CONTROL 1 (L-14): the kinked datum on ASSEMBLY's d must give pmax-h3v's 0.055612
    c1 = run(4.0, 0.0, L, 1.0, True, "control_kinked_f4_assembly_d", fine=True)
    c1["pmax_h3v_quoted"] = 0.055612
    c1["rel_to_pmax_h3v"] = abs(c1["vartheta"] - 0.055612)/0.055612
    OUT["control_kinked_f4_assembly_d"] = c1
    print("CONTROL kinked f=4 d=%.6f -> vartheta=%.6g (pmax-h3v 0.055612, rel %.2e), k=%d"
          % (c1["d"], c1["vartheta"], c1["rel_to_pmax_h3v"], c1["binding_k"]), flush=True)

    # ---- CONTROL 2: the mollified datum on ASSEMBLY's d must reproduce the refutation
    c2 = run(4.0, 0.0, L, N, False, "control_mollified_f4_assembly_d", fine=False)
    c2["referee_quoted"] = 8.76183
    OUT["control_mollified_f4_assembly_d"] = c2
    print("CONTROL mollified f=4 stock d=%.6f -> vartheta=%.6g (referee 8.76183), k=%d at phi=%.4f deg"
          % (c2["d"], c2["vartheta"], c2["binding_k"], c2["binding_phi_deg"]), flush=True)

    # ---- THE REPAIR: offset = 4 sigma, f = 4, 8, 16, 32
    rows = {}
    for f in [4.0, 8.0, 16.0, 32.0]:
        r = run(f, 4.0*SIGMA, L, N, False, "repaired_f%g" % f, fine=(f == 4.0))
        rows["f=%g" % f] = r
        print("REPAIRED f=%-4g d=%.6f (%s) -> vartheta=%.6g (k=%d) at phi=%.4f deg "
              "= delta + %.3f sigma ; 2 vartheta = %.6g  admissible=%s"
              % (f, r["d"], r["d_branches"]["binding"], r["vartheta"], r["binding_k"],
                 r["binding_phi_deg"], r["binding_phi_minus_delta_over_sigma"],
                 r["safety_vartheta"], r["admissible_at_SAFETY_2"]), flush=True)
        OUT["repaired"] = rows
        json.dump(OUT, open(os.path.join(HERE, "x1_results.json"), "w"), indent=1, default=str)

    # ---- the offset scan: how far off the cone the ball has to stand
    sc = {}
    for m in [0.0, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]:
        r = run(4.0, m*SIGMA, L, N, False, "offset_%gsigma" % m)
        sc["offset=%g sigma" % m] = {"d": r["d"], "vartheta": r["vartheta"],
                                     "binding_k": r["binding_k"],
                                     "2vartheta": r["safety_vartheta"],
                                     "admissible": r["admissible_at_SAFETY_2"],
                                     "phi_min_deg": r["phi_min_on_ball_deg"]}
        print("  offset %.1f sigma  d=%.6f  vartheta=%.6g  2v=%.6g  admissible=%s"
              % (m, r["d"], r["vartheta"], r["safety_vartheta"],
                 r["admissible_at_SAFETY_2"]), flush=True)
    OUT["offset_scan_f4"] = sc

    # ---- refinement stability of the headline number
    coarse = run(4.0, 4.0*SIGMA, L, N, False, "stab_coarse", fine=False)
    fine = OUT["repaired"]["f=4"]
    OUT["grid_stability_f4"] = {"coarse": coarse["vartheta"], "fine": fine["vartheta"],
                                "rel_change": abs(fine["vartheta"] - coarse["vartheta"])
                                / max(fine["vartheta"], 1e-300)}
    OUT["headline"] = {
        "f": 4.0, "d": fine["d"], "d_formula": "rho_* sin(phi_0 - delta - 4 sigma)",
        "vartheta": fine["vartheta"], "binding_k": fine["binding_k"],
        "binding_phi_deg": fine["binding_phi_deg"],
        "SAFETY": 2.0, "vartheta_stated": fine["safety_vartheta"],
        "Q": (1.0 + fine["safety_vartheta"])/(1.0 - fine["safety_vartheta"]),
        "vt": fine["safety_vartheta"]/(1.0 - fine["safety_vartheta"]),
        "admissible": fine["admissible_at_SAFETY_2"],
        "referee_estimate": 0.1006}
    json.dump(OUT, open(os.path.join(HERE, "x1_results.json"), "w"), indent=1, default=str)
    print("\nHEADLINE:", json.dumps(OUT["headline"], indent=1))
    print("wrote x1_results.json")

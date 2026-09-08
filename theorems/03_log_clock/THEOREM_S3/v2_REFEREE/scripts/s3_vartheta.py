"""
s3 -- (H3'_vartheta) FOR THE sigma-MOLLIFIED DATUM.   (see NOTE.md sec.3 for the argument)

THEOREM_S3_v2 row B4 imports Theorem V.4' "with (H3'_vartheta) for k = 0..4;
vartheta(f = 4) = 0.0556" from round2/pmax-h3v Part B (ADDENDUM_1).  That number was computed
for the KINKED angular profile, where on B(x_*,d) the angular factor is IDENTICALLY 1
(pmax-h3v sec.B1), so eta_0 - eta_P is purely radial.

(H3'_vartheta), pmax-h3v lines 674-680:
    sup_{|v|=1} | d_v^k (eta_0 - eta_P)(x) |  <=  vartheta k! M / r(x)^{k+1},  k = 0..4,
    for every x in B(x_*, d),   vartheta in [0,1) ,   eta_P = -M/r .
    d = min( f , rho_* sin(phi_0 - delta) , rho_* sin(pi/2 - delta_m - phi_0) ) rho_0 ,
    rho_* = (1+f) rho_0 ,  phi_0 = 30 deg .
At f = 4 the second branch binds: d = 5 sin(22.5 deg) = 1.913417 rho_0, so the ball is
EXACTLY TANGENT to the cone phi = delta.  Mollifying at scale sigma puts the smoothing layer
|phi - delta| <~ 3 sigma INSIDE the ball.

This computes a rigorous LOWER bound on vartheta by exhibiting points and meridian directions
and evaluating r^{k+1}|d_v^k e|/k! by truncated Taylor arithmetic along x + t v.
CONTROL: the same instrument on the kinked datum must reproduce pmax-h3v's 0.055612 at f = 4.
Outputs -> s3_results.json
"""
import json, math, os
import numpy as np
from scipy.interpolate import CubicSpline
import s1_profile_indep as P

HERE = os.path.dirname(os.path.abspath(__file__))
DEG = math.pi/180.0
DELTA, DM = 7.5*DEG, 5.0*DEG
PHI0 = 30.0*DEG
SIGMA = P.SIGMA
NORD = 5

def chippp(y, sg=SIGMA):
    return (y/sg**4)*(3.0 - (y*y)/sg**2)*P.chi(y, sg)

# ---------------- W_sigma and its first four derivatives, tabulated then splined
def build_W_tables(lo=0.02, hi=1.40, n=3501):
    ph = np.linspace(lo, hi, n)
    cols = []
    for ker in (P.chi, P.chip, P.chipp, chippp):
        pass
    W0 = np.array([P._conv(P.Wtil,  x, P.chi)    for x in ph])
    W1 = np.array([P._conv(P.Wtilp, x, P.chi)    for x in ph])
    W2 = np.array([P._conv(P.Wtilp, x, P.chip)   for x in ph])
    W3 = np.array([P._conv(P.Wtilp, x, P.chipp)  for x in ph])
    W4 = np.array([P._conv(P.Wtilp, x, chippp)   for x in ph])
    return [CubicSpline(ph, a) for a in (W0, W1, W2, W3, W4)], (lo, hi)

def W_derivs_kinked(phi):
    s, c = math.sin(phi), math.cos(phi)
    return [1.0/s, -c/s**2, (1.0 + c*c)/s**3, -c*(5.0 + c*c)/s**4,
            (5.0 + 18.0*c*c + c**4)/s**5]

# ---------------- truncated Taylor algebra (t^0..t^4)
def tmul(a, b):
    o = [0.0]*NORD
    for i in range(NORD):
        ai = a[i]
        if ai == 0.0:
            continue
        for j in range(NORD - i):
            o[i+j] += ai*b[j]
    return o

def tcompose(fd, s):
    du = [0.0] + list(s[1:])
    out = [fd[0]] + [0.0]*(NORD-1)
    pw = [1.0] + [0.0]*(NORD-1)
    fact = 1.0
    for k in range(1, NORD):
        pw = tmul(pw, du)
        fact *= k
        fk = fd[k]/fact
        for i in range(NORD):
            out[i] += fk*pw[i]
    return out

def tsqrt(s):
    a = s[0]
    return tcompose([a**0.5, 0.5*a**-0.5, -0.25*a**-1.5, 0.375*a**-2.5, -0.9375*a**-3.5], s)

def tlog(s):
    a = s[0]
    return tcompose([math.log(a), 1.0/a, -1.0/a**2, 2.0/a**3, -6.0/a**4], s)

def tinv(s):
    a = s[0]
    return tcompose([1.0/a, -1.0/a**2, 2.0/a**3, -6.0/a**4, 24.0/a**5], s)

def tatan(s):
    a = s[0]
    d1 = 1.0/(1.0+a*a)
    return tcompose([math.atan(a), d1, -2.0*a*d1**2, (-2.0+6.0*a*a)*d1**3,
                     (24.0*a - 24.0*a**3)*d1**4], s)

# ---------------- e = eta_0 - eta_P  along the line x + t v
def vartheta_at(r0, z0, alpha, L, N, spl):
    r = [r0, math.cos(alpha), 0.0, 0.0, 0.0]
    z = [z0, math.sin(alpha), 0.0, 0.0, 0.0]
    rho2 = tmul(r, r)
    zz = tmul(z, z)
    rho2 = [rho2[i] + zz[i] for i in range(NORD)]
    rho = tsqrt(rho2)
    u = tlog(rho)
    phi = tatan(tmul(r, tinv(z)))
    Th = tcompose(P.Theta_derivs(u[0], L), u)
    Wd = W_derivs_kinked(phi[0]) if spl is None else [float(s(phi[0])) for s in spl]
    Wc = tcompose(Wd, phi)
    term = tmul(Th, tmul(Wc, tinv(rho)))
    ri = tinv(r)
    c = [ri[i] - term[i]/N for i in range(NORD)]
    return [r0**(k+1)*abs(c[k]) for k in range(NORD)]

def scan(f, L, N, spl, nrad=48, nang=96, nalpha=36, shrink=1.0):
    rho_s = 1.0 + f
    d = shrink*min(f, rho_s*math.sin(PHI0 - DELTA), rho_s*math.sin(math.pi/2 - DM - PHI0))
    xc = (rho_s*math.sin(PHI0), rho_s*math.cos(PHI0))
    best = [0.0]*NORD
    arg = [None]*NORD
    for i in range(nrad+1):
        rad = d*i/nrad
        for j in range(nang):
            th = 2.0*math.pi*j/nang
            r0 = xc[0] + rad*math.cos(th)
            z0 = xc[1] + rad*math.sin(th)
            if r0 <= 1e-3 or z0 <= 1e-3:
                continue
            ph = math.atan2(r0, z0)
            if spl is not None and not (0.03 < ph < 1.38):
                continue
            for m in range(nalpha):
                al = math.pi*m/nalpha
                v = vartheta_at(r0, z0, al, L, N, spl)
                for k in range(NORD):
                    if v[k] > best[k]:
                        best[k] = v[k]
                        arg[k] = {"r": r0, "z": z0, "alpha_deg": al/DEG,
                                  "phi_deg": ph/DEG, "rho": math.hypot(r0, z0),
                                  "phi_minus_delta_over_sigma": (ph-DELTA)/SIGMA}
    return {"f": f, "d": d, "shrink": shrink, "x_star": list(xc),
            "vartheta_k": best, "vartheta": max(best), "argmax": arg,
            "binding_k": int(np.argmax(best))}


if __name__ == "__main__":
    L = 40.0
    N = 1.0124508490
    print("building W_sigma tables ...", flush=True)
    SPL, rng = build_W_tables()
    print("  done, range", rng, flush=True)
    # spline sanity: W'''' at phi = delta should be ~ Delta * chi''(0)
    dW = P.Wp_exact(DELTA + 1e-7) - P.Wp_exact(DELTA - 1e-7)
    print("  jump of W' at delta = %.6f ; chi''(0) = %.4f ; predicted W'''' = %.5g ; "
          "table = %.5g" % (dW, -1.0/SIGMA**2*P.chi(0.0), dW*(-1.0/SIGMA**2*P.chi(0.0)),
                            float(SPL[4](DELTA))), flush=True)

    OUT = {"definition": "vartheta >= r^{k+1}|d_v^k(eta_0-eta_P)|/(k! M), k=0..4, eta_P=-M/r; "
                         "LOWER bound by scan of B(x_*,d) x meridian directions",
           "L": L, "N_sigma": N, "sigma": SIGMA,
           "W4_at_delta_predicted": dW*(-1.0/SIGMA**2*P.chi(0.0)),
           "W4_at_delta_table": float(SPL[4](DELTA))}

    ctl = scan(4.0, L, 1.0, None, nrad=32, nang=64, nalpha=24)
    ctl["pmax_h3v_quoted_vartheta_f4"] = 0.055612
    OUT["control_kinked_f4"] = ctl
    print("CONTROL kinked  f=4  d=%.6f : " % ctl["d"],
          ["%.6g" % v for v in ctl["vartheta_k"]], "-> vartheta =", "%.6g" % ctl["vartheta"],
          flush=True)

    rows = {}
    for f in [4.0, 8.0, 16.0, 32.0]:
        r = scan(f, L, N, SPL, nrad=48, nang=96, nalpha=36)
        rows["f=%g" % f] = r
        print("MOLLIFIED f=%-4g d=%.6f : " % (f, r["d"]),
              ["%.6g" % v for v in r["vartheta_k"]], "-> vartheta = %.6g (k=%d) at phi=%.4f deg"
              % (r["vartheta"], r["binding_k"], r["argmax"][r["binding_k"]]["phi_deg"]),
              flush=True)
        OUT["mollified"] = rows
        json.dump(OUT, open(os.path.join(HERE, "s3_results.json"), "w"), indent=1, default=str)

    sh = {}
    for s in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]:
        r = scan(4.0, L, N, SPL, nrad=24, nang=48, nalpha=18, shrink=s)
        sh["shrink=%g" % s] = {"d": r["d"], "vartheta": r["vartheta"],
                               "binding_k": r["binding_k"], "vartheta_k": r["vartheta_k"],
                               "min_phi_minus_delta_over_sigma":
                                   r["argmax"][r["binding_k"]]["phi_minus_delta_over_sigma"]}
        print("  shrink %.2f  d=%.6f  vartheta=%.6g  (k=%d)"
              % (s, r["d"], r["vartheta"], r["binding_k"]), flush=True)
    OUT["d_shrink_scan_f4"] = sh
    json.dump(OUT, open(os.path.join(HERE, "s3_results.json"), "w"), indent=1, default=str)

"""
s2 -- THE WINDOW AND THE CAP, re-derived independently for the mollified profile.

    Q(lam) = P_h(lam)/lam = int_0^1 3 lam^9 v^2 h(arcsin v) D^{-5/2} dv ,  D = 1+(lam^6-1)v^2
    (a probability average: int_0^1 3 lam^9 v^2 D^{-5/2} dv = 1 EXACTLY, checked below)
    Q'(lam) = -9 lam^8 int_0^1 h'(arcsin v) v^3 sqrt(1-v^2) D^{-5/2} dv
    kappa_delta = P_h(1)/2 ,  c_* = log(3/2)/kappa ,  lam_max(c) = e^{3c/4} ,
    r_h = inf_{[1,lam_max]} Q ,  P_h' = Q + lam Q' .

Checks:
 (1) r_h and min P_h' on [1, lam_max(c_*)] against fix4's 0.9033435309 / 0.3335513465 ;
 (2) lam_mono and the monotone cap 1.1760705119 / the certified cap 1.1664585076 ;
 (3) fix4/FIX4.md sec.2.1's claim that "the majorant 9 lam^9 sum_i max(M_i,0) I_i is
     INCREASING in lam (lam^9 up, every I_i down), [so] testing it at lam_max certifies the
     whole interval" -- tested directly ;
 (4) lam_max evaluated at the SELF-CONSISTENT c, not at c_*.
Outputs -> s2_results.json
"""
import json, math, os
import numpy as np
from scipy.interpolate import CubicSpline
import s1_profile_indep as P

HERE = os.path.dirname(os.path.abspath(__file__))
DEG = math.pi/180.0
LOG32 = math.log(1.5)

# --------- h and h' on [0, pi/2] from the independent quadrature, splined
NPH = 6001
PH = np.linspace(1e-7, math.pi/2, NPH)
A0 = np.array([P.As(x) for x in PH])
A1 = np.array([P.As1(x) for x in PH])
N_SIG = 1.0124508490
SA = CubicSpline(PH, A0/N_SIG)
SA1 = CubicSpline(PH, A1/N_SIG)

# --------- quadrature in v, panelled and Gauss-Legendre
def panels(edges, ng=200):
    xs, ws = [], []
    x, w = np.polynomial.legendre.leggauss(ng)
    for a, b in zip(edges[:-1], edges[1:]):
        xs.append(0.5*(b-a)*x + 0.5*(a+b))
        ws.append(0.5*(b-a)*w)
    return np.concatenate(xs), np.concatenate(ws)

EDG = [0.0, math.sin(7.5*DEG)*0.5, math.sin(7.5*DEG), math.sin(9.99*DEG),
       math.cos(5.0*DEG), 0.999, 1.0]
VQ, WQ = panels(EDG, ng=300)
PHQ = np.arcsin(np.clip(VQ, 0.0, 1.0))
HQ = SA(PHQ)
HPQ = SA1(PHQ)

def Q(lam):
    D = 1.0 + (lam**6 - 1.0)*VQ**2
    return 3.0*lam**9*float(np.dot(WQ*HQ*VQ**2, D**-2.5))

def norm_check(lam):
    D = 1.0 + (lam**6 - 1.0)*VQ**2
    return 3.0*lam**9*float(np.dot(WQ*VQ**2, D**-2.5))

def Qp(lam):
    D = 1.0 + (lam**6 - 1.0)*VQ**2
    return -9.0*lam**8*float(np.dot(WQ*HPQ*VQ**3*np.sqrt(np.maximum(1-VQ**2, 0.0)), D**-2.5))

def Ph(lam):  return lam*Q(lam)
def Php(lam): return Q(lam) + lam*Qp(lam)

if __name__ == "__main__":
    R = {"instrument": "independent: adaptive-quadrature profile splined, 300-pt Gauss "
                       "panels in v; no shared code with fix4/p1_profile.py"}
    R["probability_normalisation_check"] = {("lam=%g" % l): norm_check(l)
                                            for l in [1.0, 1.2, 1.5, 1.8558724485, 2.0570755604]}
    kap = 0.5*Ph(1.0)
    cst = LOG32/kap
    R["kappa_delta"] = kap
    R["c_star"] = cst
    R["fix4"] = {"kappa_delta": 0.4917868801, "c_star": 0.8244732108,
                 "r_h": 0.9033435309, "min_Ph_prime": 0.3335513465,
                 "Ph_prime_at_lam_max": 0.3549367307, "lam_mono": 2.0693384635,
                 "cap_monotone": 1.1760705119, "cap_certified": 1.1664585076,
                 "lam_max_at_c_star": 1.8558724485, "lam_max_at_cap": 2.0570755604}

    for tag, lmx in [("at c_star", math.exp(0.75*cst)),
                     ("at certified cap 1.1664585076", math.exp(0.75*1.1664585076*cst)),
                     ("at monotone cap 1.1760705119", math.exp(0.75*1.1760705119*cst))]:
        ls = np.linspace(1.0, lmx, 4001)
        qs = np.array([Q(l) for l in ls])
        ps = np.array([Php(l) for l in ls])
        R[tag] = {"lam_max": lmx, "r_h_min": float(qs.min()),
                  "argmin_lam": float(ls[int(np.argmin(qs))]),
                  "min_Ph_prime": float(ps.min()), "argmin_Ph_prime": float(ls[int(np.argmin(ps))]),
                  "Ph_prime_at_lam_max": float(ps[-1]), "Q_at_1": float(qs[0])}
        print(tag, "lam_max=%.10f r_h=%.10f minPh'=%.10f Ph'(lam_max)=%.10f"
              % (lmx, qs.min(), ps.min(), ps[-1]), flush=True)

    lo, hi = 1.5, 3.2
    for _ in range(80):
        mid = 0.5*(lo+hi)
        if Php(mid) > 0: lo = mid
        else: hi = mid
    lam_mono = 0.5*(lo+hi)
    R["lam_mono"] = lam_mono
    R["cap_monotone_c_over_cstar"] = (4.0/3.0)*math.log(lam_mono)/cst
    print("lam_mono = %.10f -> c/c_* cap = %.10f" % (lam_mono, R["cap_monotone_c_over_cstar"]),
          flush=True)

    # ---- (3) is the certificate's majorant really increasing in lam?
    def prim(V, K):
        K = np.asarray(K, dtype=float)
        S = 1.0 + K*V*V
        with np.errstate(divide='ignore', invalid='ignore'):
            val = (2.0/3.0 + S**-1.5/3.0 - S**-0.5)/K**2
        return np.where(K < 1e-9, V**4/4.0, val)
    ncell = 4000
    ed = np.linspace(0.0, 1.0, ncell+1)
    phe = np.arcsin(np.clip(ed, 0.0, 1.0))
    Mi = np.array([float(np.max(SA1(np.linspace(phe[i], phe[i+1], 25))))
                   for i in range(ncell)])
    pos = np.maximum(Mi, 0.0)
    lams = np.linspace(1.0, 2.10, 441)
    maj = []
    for l in lams:
        K = l**6 - 1.0
        I = prim(ed[1:], K) - prim(ed[:-1], K)
        maj.append(9.0*l**9*float(np.dot(pos, I)))
    maj = np.array(maj)
    dif = np.diff(maj)
    R["majorant_monotone_claim"] = {
        "n_lambda": len(lams), "lambda_range": [1.0, 2.10],
        "majorant_at_1": float(maj[0]), "majorant_at_lam_max_cstar": float(np.interp(
            math.exp(0.75*cst), lams, maj)),
        "majorant_max": float(maj.max()), "argmax_lambda": float(lams[int(np.argmax(maj))]),
        "is_increasing_everywhere": bool(np.all(dif >= -1e-12)),
        "n_decreasing_steps": int(np.sum(dif < 0)),
        "first_decreasing_lambda": (float(lams[int(np.argmax(dif < 0))]) if np.any(dif < 0)
                                    else None)}
    print("majorant increasing on [1,2.1]?", R["majorant_monotone_claim"]["is_increasing_everywhere"],
          " max at lam =", R["majorant_monotone_claim"]["argmax_lambda"], flush=True)

    json.dump(R, open(os.path.join(HERE, "s2_results.json"), "w"), indent=1, default=str)
    print(json.dumps(R, indent=1, default=str))

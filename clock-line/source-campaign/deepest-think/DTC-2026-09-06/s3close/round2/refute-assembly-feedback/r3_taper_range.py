"""
r3 -- the taper constants on the EXTENDED strain range the slaved reference needs.

Under the contradiction hypothesis ||omega(t)||_inf < (3/2) M on the window, the l=1 interior
rate of every shell is at most (3/4) M per e-fold (far-near-kernel-lemma Consequence A with
M -> (3/2)M), so the reference profile lambda~(rho,s) = exp int_0^s a_1(rho,.) obeys
1 <= lambda~ <= exp((3/4) c) on the window tau = c/(M L)   (lower bound: Lemma 2, a_1 >= 0).

Lemma T' Corollary 2 needs r_h = inf P_h(lambda)/lambda over the RANGE of the reference
profile, and the conservative strain route needs Phi_h(lambda) = P_h(lambda)/P_h(1) >= 1 there.
prove-lagrangian proved both only on [1, 3/2].  Here: mpmath quadrature of
P_h(lambda) = 3 int_0^1 h(arcsin v) v^2 (A v^2 + B)^{-5/2} dv, h = min(1, phi_ax/delta),
plus the EXACT deficit bound D <= w^3/(B (A w^2+B)^{3/2}), w = sin(delta), of prove-lagrangian
section 2, evaluated on the extended range.
"""
import json, math
import mpmath as mp

mp.mp.dps = 30
DEG = mp.pi/180

def P_h(lam, delta):
    lam = mp.mpf(lam)
    A = lam**2 - lam**(-4); B = lam**(-4)
    def h(v):
        phi = mp.asin(v)
        return min(mp.mpf(1), phi/delta) if delta > 0 else mp.mpf(1)
    f = lambda v: h(v)*v*v*(A*v*v + B)**(-mp.mpf(5)/2)
    w = mp.sin(delta) if delta > 0 else mp.mpf(0)
    return 3*(mp.quad(f, [0, w]) + mp.quad(f, [w, 1]))

def deficit_bound(lam, delta):
    lam = mp.mpf(lam); A = lam**2 - lam**(-4); B = lam**(-4); w = mp.sin(delta)
    return w**3/(B*(A*w*w + B)**(mp.mpf(3)/2))

OUT = {}
for ddeg in [7.5, 15.0]:
    delta = ddeg*DEG
    P1 = P_h(1, delta)
    kappa = P1/2
    row = {"kappa_delta": float(kappa), "P_h(1)": float(P1)}
    for cname, c in [("c_star_0.8113826", 0.8113826), ("c_1", 1.0)]:
        lam_max = math.exp(0.75*c)
        lams = [1 + (lam_max-1)*k/60 for k in range(61)]
        ratios = [float(P_h(l, delta)/l) for l in lams]
        Phis = [float(P_h(l, delta)/P1) for l in lams]
        Pl = [float(P_h(l, delta)) for l in lams]
        # sanity: P_h(lam) = lam - D(lam), D within the exact deficit bound
        Dchk = max(float(l - P) for l, P in zip(lams, Pl))
        Dbnd = max(float(deficit_bound(l, delta)) for l in lams)
        row[cname] = {
            "lambda_max": lam_max,
            "r_h_over_[1,lambda_max]": min(ratios),
            "r_h_at_3over2": float(P_h(1.5, delta)/1.5),
            "argmin_lambda": lams[ratios.index(min(ratios))],
            "min_Phi_h": min(Phis), "Phi_h_ge_1_everywhere": min(Phis) >= 1.0,
            "Phi_h_at_lambda_max": Phis[-1],
            "max_deficit_measured": Dchk, "max_deficit_bound_exact": Dbnd,
            "deficit_within_bound": Dchk <= Dbnd + 1e-12,
            "ratio_monotone_decreasing": all(ratios[i] >= ratios[i+1]-1e-12 for i in range(len(ratios)-1)),
        }
    OUT["delta=%g" % ddeg] = row

with open("r3_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True)
print(json.dumps(OUT, indent=1, sort_keys=True))

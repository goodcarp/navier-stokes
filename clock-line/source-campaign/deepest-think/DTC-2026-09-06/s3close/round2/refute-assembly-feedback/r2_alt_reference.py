"""
r2 -- a2_budget.py re-run FROM A COPY with the reference map redefined.

Variants of the bootstrap (everything not named is exactly as in a2_budget_copy.py, whose
viscous budget, C_a, ell_loss, eps_Tprime, RECORD and columns are imported by exec):

  assembly       : a2 as written.  Reference Lambda fixed a priori (integro-ODE profile,
                   kappa_s = sqrt6-2); the l=1 rate error |a_true - a_ref| <= eps_T' kappa M L
                   is a DRIVE proportional to mu  ->  feedback exponent p c = 34.06.
  slaved_min     : the reference is the flow of the true field's own l=1 (shell-wise uniform
                   strain) part, a_1(|x|,s)(Y,-2Z) at the label radius.  The l=1 rate error is
                   identically zero, so the eps_T' term is removed from the RHS of (2.1)-(2.2)
                   and nothing else changes.  Exponent = c_G = (3/2 + C'/L) c.
  slaved_honest  : slaved_min with the constants the slaved reference actually has:
                   kappa_s = (3/4) c  (|d log lambda~/d log rho| <= (3/4) M tau, Consequence A
                   rate bound at |omega| <= (3/2) M), r_h over [1, e^{3c/4}] (r3), label-radius
                   drive (lambda_omega/2) * 2 log(lambda~_max) instead of kappa c_G, and the
                   lost inner e-folds weighted (3/4)M instead of kappa M.
  conjugated     : the map error measured in the reference's own strained coordinates,
                   Phi = T_lambda~ (x + e)  (r4 (iii)).  Then e is driven ONLY by the remainder
                   velocity and the label-vs-current-radius error, both bounded absolutely;
                   no Lipschitz constant of u and no Jacobian shear enter.  mu <= lambda~^3 eps,
                   mu_J from d/ds log(J_Phi/lambda~^2) = 2[a(Phi x) - a_1(|x|)].  Exponent 0
                   (the growth factor is the anisotropy lambda~^3).
  conjugated_frozen : conjugated, with the vorticity sup taken as the frozen model
                   lambda_omega = e^{kappa theta} <= 3/2 everywhere (a2's own convention) so
                   lambda~(theta) <= exp((lambda_omega - 1)/(2 kappa)) <= e^{1/2}.

Columns 'proved' (C' = 119.33, C_a = 8.6979, Thm V.1 tail, Lemma T' as proved) and 'measured'
(C' = 11.74, C_a = 0.069, exact Gaussian, sensitivity /6.0668) as in a2.
"""
import json, math, os
import numpy as np
from scipy.integrate import solve_ivp, quad

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
src = open("a2_budget_copy.py").read().split('if __name__ == "__main__":')[0]
G = {"__name__": "r2"}
exec(src, G)
RECORD = G["RECORD"]; KAPPA = G["KAPPA"]; R_H = G["R_H"]; KAPPA_S = G["KAPPA_S"]
COLUMNS = G["COLUMNS"]; LAM_MAX = G["LAM_MAX"]; PHI0 = G["PHI0"]; DELTA = G["DELTA"]
eps_Tprime = G["eps_Tprime"]; viscous_budget = G["viscous_budget"]; ell_loss = G["ell_loss"]
C_a_proved = G["C_a_proved"]; lam_frozen = G["lam_frozen"]; SENS = G["SENS_MEASURED_FACTOR"]
CSTAR = math.log(1.5)/KAPPA

# ---- r_h on an arbitrary range, same integrand as a1 (delta_m = 0), scipy quad
def P_h(lam, delta=DELTA):
    A = lam**2 - lam**(-4); B = lam**(-4)
    def f(v):
        phi = math.asin(min(1.0, v))
        return min(1.0, phi/delta)*v*v*(A*v*v + B)**(-2.5)
    w = math.sin(delta)
    return 3.0*(quad(f, 0.0, w, epsabs=1e-13, epsrel=1e-13, limit=200)[0]
                + quad(f, w, 1.0, epsabs=1e-13, epsrel=1e-13, limit=200)[0])
def r_h_range(lam_max):
    lams = np.linspace(1.0, lam_max, 61)
    return float(min(P_h(l)/l for l in lams))
R3 = json.load(open("r3_results.json"))
RH_EXT = R3["delta=7.5"]["c_star_0.8113826"]["r_h_over_[1,lambda_max]"]

# ---------------------------------------------------------------- variants
def bootstrap_variant(L, c, Cprime, C_a, sens, variant):
    GoML = 1.5 + Cprime/L
    c_G = GoML*c
    geom_lo, geom_hi = 4.0/9.0, 1.5
    lam_t_max = math.exp(0.75*c)                      # sup of the reference profile, honest
    if variant == "assembly":
        kappa_s, r_h, label, ell_w = KAPPA_S, R_H, KAPPA*c_G, 1.0
        feedback = True
    elif variant == "slaved_min":
        kappa_s, r_h, label, ell_w = KAPPA_S, R_H, KAPPA*c_G, 1.0
        feedback = False
    elif variant == "slaved_honest":
        kappa_s, r_h, ell_w = 0.75*c, r_h_range(lam_t_max), 0.75/KAPPA
        label = None
        feedback = False
    if variant in ("assembly", "slaved_min", "slaved_honest"):
        def rhs(th, y):
            m, muJL = y
            mu = min(m/geom_lo, 0.9)
            muJ = muJL*(1.0+2.0*kappa_s/L) + 2.0*kappa_s/L
            lam_w = min(lam_frozen(th), LAM_MAX)
            epsT = eps_Tprime(mu, muJ, r_h, sens) if feedback else 0.0
            lab = label if label is not None else (lam_w/2.0)*2.0*math.log(lam_t_max)
            drive = 2.0*(KAPPA*epsT + lab/L) + C_a*lam_w/L
            dm = GoML*m + geom_hi*drive
            dmuJ = 2.0*(KAPPA*epsT + lab/L + C_a*lam_w/L + GoML*mu)
            return [dm, dmuJ]
        def ev(th, y):
            return y[0]/geom_lo - 0.9
        ev.terminal = True; ev.direction = 1
        sol = solve_ivp(rhs, (0.0, c), [0.0, 0.0], rtol=1e-8, atol=1e-14, events=ev,
                        max_step=c/50.0)
        m, muJL = sol.y[0, -1], sol.y[1, -1]
        mu = m/geom_lo
        if sol.t[-1] < c - 1e-12:
            mu = 1.0
        muJ = muJL*(1.0+2.0*kappa_s/L) + 2.0*kappa_s/L
        blew = (not sol.success) or (mu >= 0.899) or (not np.isfinite(mu))
        epsT = eps_Tprime(mu, muJ, r_h, sens) if not blew else float("inf")
        if feedback:
            p_over_ML = GoML + geom_hi*2.0*KAPPA*sens*(15.0*math.pi/4.0)/(r_h*geom_lo)
        else:
            p_over_ML = GoML
        return {"mu": mu, "muJ": muJ, "eps_Tprime": epsT, "c_G": c_G, "closed": not blew,
                "exponent": p_over_ML*c, "p_over_ML": p_over_ML, "r_h": r_h,
                "kappa_s": kappa_s, "ell_weight": ell_w, "lam_ref_max": lam_t_max}
    # ---- conjugated variants
    frozen = (variant == "conjugated_frozen")
    if frozen:
        lam_ref = lambda th: math.exp((min(lam_frozen(th), LAM_MAX) - 1.0)/(2.0*KAPPA))
        lam_ref_max = lam_ref(c)
    else:
        lam_ref = lambda th: math.exp(0.75*th)
        lam_ref_max = lam_t_max
    r_h = r_h_range(lam_ref_max)
    def rhs(th, y):
        e, muJ = y
        e = min(e, 0.9)
        lam_w = min(lam_frozen(th), LAM_MAX)
        lt = lam_ref(th)
        lab = (lam_w/2.0)*(2.0*math.log(lt) + e/(1.0-e))   # |a_1(|X|) - a_1(|x|)| / M
        de = (1.0+e)/L*(2.0*lab + lt**3*C_a*lam_w)         # |D(x+e)| <= 2|x+e| ; |T^-1 R| <= lt^3 C_a M_s |x+e|
        dmuJ = (2.0/L)*(C_a*lam_w + lab)
        return [de, dmuJ]
    def ev(th, y):
        return y[0] - 0.9
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, (0.0, c), [0.0, 0.0], rtol=1e-8, atol=1e-14, events=ev,
                    max_step=c/50.0)
    e, muJ = sol.y[0, -1], sol.y[1, -1]
    if sol.t[-1] < c - 1e-12:
        e = 1.0
    mu = lam_ref_max**3*e
    blew = (not sol.success) or (mu >= 0.999) or (not np.isfinite(mu))
    epsT = eps_Tprime(mu, muJ, r_h, sens) if not blew else float("inf")
    # linearised coefficient of e in de/dtheta, integrated over the window: the exponent
    expo = quad(lambda th: (1.0/L)*(2.0*(min(lam_frozen(th), LAM_MAX)/2.0)*(2.0*math.log(lam_ref(th)) + 1.0)
                                    + lam_ref(th)**3*C_a*min(lam_frozen(th), LAM_MAX)), 0.0, c)[0]
    return {"mu": mu, "muJ": muJ, "eps_Tprime": epsT, "c_G": c_G, "closed": not blew,
            "exponent": expo, "p_over_ML": expo/c, "r_h": r_h, "kappa_s": 0.75*c,
            "ell_weight": 0.75/KAPPA, "lam_ref_max": lam_ref_max, "eps_conj": e,
            "anisotropy_lam3": lam_ref_max**3}

def assemble_variant(L, c, column, variant, C_K=1.0, f=1.0):
    cp_key, ca_kind, use_V1, sens = COLUMNS[column]
    Cprime = RECORD[cp_key]
    C_a = C_a_proved(math.sin(PHI0)) if ca_kind == "far_near" else RECORD["material_offset_M"]
    bs = bootstrap_variant(L, c, Cprime, C_a, sens, variant)
    vb = viscous_budget(L, c, C_K, f, c_G=bs["c_G"], Cprime=Cprime, use_V1_for_Z=use_V1)
    ell = ell_loss(f)
    eps_a = bs["ell_weight"]*ell/L + bs["eps_Tprime"] + LAM_MAX*C_a/(KAPPA*L)
    eps_v = vb["eps_v"]
    eps_delta = 1.0 - 2.0*KAPPA
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        eps = float("inf")
    else:
        eps = (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-eps_delta)) - 1.0
    out = dict(bs)
    out.update({"L": L, "c": c, "column": column, "variant": variant, "C_K": C_K, "f": f,
                "eps_ell": bs["ell_weight"]*ell/L, "eps_Ca": LAM_MAX*C_a/(KAPPA*L),
                "eps_a": eps_a, "eps_v": eps_v, "eps": eps, "Cprime": Cprime, "C_a": C_a})
    return out

FS = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
def eps_best_f(L, column, variant, C_K=1.0):
    return min(assemble_variant(L, CSTAR, column, variant, C_K=C_K, f=f)["eps"] for f in FS)
def L_star(column, variant, target, C_K=1.0):
    lo, hi = 1.0, 1e60
    if eps_best_f(hi, column, variant, C_K) > target:
        return None
    for _ in range(60):
        mid = math.sqrt(lo*hi)
        if eps_best_f(mid, column, variant, C_K) <= target:
            hi = mid
        else:
            lo = mid
        if hi/lo < 1.0001:
            break
    return hi

if __name__ == "__main__":
    VARIANTS = ["assembly", "slaved_min", "slaved_honest", "conjugated", "conjugated_frozen"]
    OUT = {"c_star": CSTAR, "kappa": KAPPA, "r_h_record": R_H, "r_h_ext_r3": RH_EXT,
           "r_h_ext_here": r_h_range(math.exp(0.75*CSTAR)),
           "r_h_frozen_range": r_h_range(math.exp((LAM_MAX-1.0)/(2*KAPPA)))}
    rows = {}
    print("%-18s %-9s %7s %9s %9s %9s %9s %9s %9s" % ("variant", "column", "L", "expo", "mu", "muJ", "epsT", "eps_a", "eps"))
    for variant in VARIANTS:
        for column in ["proved", "measured"]:
            for L in [40.0, 160.0, 640.0, 2560.0, 1e4, 1e5]:
                best = None
                for f in FS:
                    r = assemble_variant(L, CSTAR, column, variant, f=f)
                    if best is None or r["eps"] < best["eps"]:
                        best = r
                key = "%s|%s|L=%g" % (variant, column, L)
                rows[key] = best
                print("%-18s %-9s %7g %9.4g %9.3g %9.3g %9.3g %9.3g %9.4g" %
                      (variant, column, L, best["exponent"], best["mu"], best["muJ"],
                       best["eps_Tprime"], best["eps_a"], best["eps"]))
    OUT["rows"] = rows
    Ls = {}
    print("\n%-18s %-9s %9s %14s %14s" % ("variant", "column", "C_K", "L*(eps<=.5)", "L*(eps<=.1)"))
    for variant in VARIANTS:
        for column in ["proved", "measured"]:
            for C_K in [1.0, 66.6622]:
                a = L_star(column, variant, 0.5, C_K); b = L_star(column, variant, 0.1, C_K)
                Ls["%s|%s|C_K=%g" % (variant, column, C_K)] = {"L_star_eps<=0.5": a, "L_star_eps<=0.1": b}
                print("%-18s %-9s %9g %14.6g %14.6g" % (variant, column, C_K,
                      a if a else float('nan'), b if b else float('nan')))
    OUT["L_star"] = Ls
    # the exponent alone, L -> infinity, at c = c*
    ex = {}
    for variant in VARIANTS:
        for column in ["proved", "measured"]:
            cp_key, ca_kind, use_V1, sens = COLUMNS[column]
            C_a = C_a_proved(math.sin(PHI0)) if ca_kind == "far_near" else RECORD["material_offset_M"]
            for L in [160.0, 1e18]:
                bs = bootstrap_variant(L, CSTAR, RECORD[cp_key], C_a, sens, variant)
                ex["%s|%s|L=%g" % (variant, column, L)] = {"exponent": bs["exponent"],
                                                          "exp(exponent)": math.exp(bs["exponent"]),
                                                          "p_over_ML": bs["p_over_ML"]}
    OUT["exponent"] = ex
    with open("r2_results.json", "w") as fh:
        json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
    print("\nexponents:")
    for k, v in ex.items():
        print("  %-40s expo=%10.5g  e^expo=%12.5g" % (k, v["exponent"], v["exp(exponent)"]))

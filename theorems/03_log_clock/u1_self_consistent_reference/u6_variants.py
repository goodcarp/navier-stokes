"""
u6 -- two variants of the u5 budget that the note has to report but that do not belong in the
headline table.

  (V-strain)  the Theorem V.4 viscous loss charged to the STRAIN lower bound as well as to the
              vorticity (ASSEMBLY's clause (iii) charges it to neither; the brief asks for it).
              eps_a -> eps_a + eps_v .
  (conj)      the budget run with the CONJUGATED bootstrap (2.6) instead of the Groenwall one
              (2.3):  mu = lam_max^3 ehat, no (Gamma-off) in the map error.

Everything else is u5_budget.py's, imported (not re-implemented): this file is a variant
report, not an independent instrument.
"""
import json, math, os
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
src = open("u5_budget.py").read().split('if __name__ == "__main__":')[0]
G = {"__name__": "u6"}
exec(src, G)
assemble = G["assemble"]; COLUMNS = G["COLUMNS"]; CSTAR = G["CSTAR"]
LAM_MAX = G["LAM_MAX"]; LAM_OM = G["LAM_OM"]; KAPPA = G["KAPPA"]
LOGREE_SHIFT = G["LOGREE_SHIFT"]; eps_Tprime = G["eps_Tprime"]
viscous_budget = G["viscous_budget"]; ell_loss = G["ell_loss"]; RECORD = G["RECORD"]
C_a_proved = G["C_a_proved"]; PHI0 = G["PHI0"]; r_h_of = G["r_h_of"]
C_R_PROVED = G["C_R_PROVED"]; C_R_MEASURED = G["C_R_MEASURED"]
SENS = G["SENS_MEASURED_FACTOR"]
FS = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]

def clock(eps_a, eps_v):
    eps_delta = 1.0 - 2.0*KAPPA
    if (not np.isfinite(eps_v)) or eps_v >= 1.0 or eps_a >= 1.0 or not np.isfinite(eps_a):
        return float("inf")
    return (1.0 + math.log(1.0/(1.0-eps_v))/math.log(1.5))/((1.0-eps_a)*(1.0-eps_delta)) - 1.0

# ---------------- (V-strain) ---------------------------------------------------------
def eps_vstrain(L, col, f):
    r = assemble(L, CSTAR, col, f=f)
    return clock(r["eps_a"] + r["eps_v"], r["eps_v"])

# ---------------- (conj) --------------------------------------------------------------
def bootstrap_conj(L, c, C_R, C_a, lam_mx=LAM_MAX, lam_om=LAM_OM):
    K = lam_mx**3*C_R
    def rhs(th, y):
        eh, w = y
        if eh >= 0.999:
            return [0.0, 0.0]
        Lr = 2.0*math.log(lam_mx) + math.log(1.0/(1.0 - eh))
        Lr2 = max(math.log(lam_mx*(1.0 + eh)), Lr)
        return [(lam_om/L)*(Lr + K)*(1.0 + eh), (2.0*lam_om/L)*(C_a + 0.5*Lr2)]
    def ev(th, y):
        return y[0] - 0.99
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(rhs, (0.0, c), [0.0, 0.0], rtol=1e-10, atol=1e-14, events=ev,
                    max_step=c/200.0)
    eh, w = float(sol.y[0, -1]), float(sol.y[1, -1])
    blew = (sol.t[-1] < c - 1e-12) or (not sol.success)
    mu = lam_mx**3*eh if not blew else float("inf")
    muJ = (math.exp(w) - 1.0) if (not blew and w < 50) else float("inf")
    return mu, muJ, (lam_om/L)*(1.0 + 2.0*math.log(lam_mx) + K)*c

def eps_conj(L, col, f):
    cp_key, ca_kind, cr_kind, use_V1, sens, k2key = COLUMNS[col]
    Cpp = RECORD[cp_key]
    C_a = C_a_proved(math.sin(PHI0)) if ca_kind == "far_near" else RECORD["material_offset_M"]
    C_R = C_R_PROVED if cr_kind == "proved" else C_R_MEASURED
    mu, muJ, pc = bootstrap_conj(L, CSTAR, C_R, C_a)
    c_G = (1.5*(1.0 + 1.0/L) + Cpp/L)*CSTAR          # (Gamma-off) still needed for the viscous half
    vb = viscous_budget(L, CSTAR, RECORD[k2key], f, c_G=c_G, Cprime=Cpp, use_V1_for_Z=use_V1)
    epsT = eps_Tprime(mu, muJ, r_h_of(LAM_MAX), sens)
    ell = ell_loss(f, min(mu, 0.999) if np.isfinite(mu) else 0.999)
    eps_a = ell/L + epsT + LAM_OM*C_a/(KAPPA*L)
    return clock(eps_a, vb["eps_v"]), pc, mu

def best(fn, L, col):
    return min(fn(L, col, f) if not isinstance(fn(L, col, f), tuple) else fn(L, col, f)[0]
               for f in FS)

def Lstar(fn, col, target):
    lo, hi = 2.0, 1e18
    if best(fn, hi, col) > target:
        return None
    for _ in range(60):
        mid = math.sqrt(lo*hi)
        if best(fn, mid, col) <= target:
            hi = mid
        else:
            lo = mid
        if hi/lo < 1.0000001:
            break
    return hi

OUT = {}
for col in list(COLUMNS):
    OUT[col] = {}
    for L in [40.0, 160.0, 640.0, 2560.0]:
        OUT[col]["eps_Vstrain_L=%g" % L] = best(eps_vstrain, L, col)
        OUT[col]["eps_conj_L=%g" % L] = best(eps_conj, L, col)
    OUT[col]["L_star_Vstrain_0.5"] = Lstar(eps_vstrain, col, 0.5)
    OUT[col]["L_star_conj_0.5"] = Lstar(eps_conj, col, 0.5)
    OUT[col]["logLambda_star_Vstrain_0.5"] = (2*OUT[col]["L_star_Vstrain_0.5"] + LOGREE_SHIFT
                                              if OUT[col]["L_star_Vstrain_0.5"] else None)
    OUT[col]["logLambda_star_conj_0.5"] = (2*OUT[col]["L_star_conj_0.5"] + LOGREE_SHIFT
                                           if OUT[col]["L_star_conj_0.5"] else None)
    OUT[col]["pc_conj_L=640"] = eps_conj(640.0, col, 1.0)[1]

with open("u6_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
print(json.dumps(OUT, indent=1, sort_keys=True, default=str))

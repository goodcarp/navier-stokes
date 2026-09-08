"""
r1 -- reproduce the assembly's feedback exponent and L_* FROM A COPY of a2_budget.py
(a2_budget_copy.py, sha recorded in SHA256SUMS), and check the linearisation claim of
ASSEMBLY.md section 3.2 numerically against the majorant ODE itself.

Nothing is typed from the record: every number is computed here from the copied script.
"""
import json, math, os
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
src = open("a2_budget_copy.py").read().split('if __name__ == "__main__":')[0]
G = {"__name__": "r1"}
exec(src, G)

assemble = G["assemble"]; bootstrap = G["bootstrap"]; COLUMNS = G["COLUMNS"]
RECORD = G["RECORD"]; KAPPA = G["KAPPA"]; R_H = G["R_H"]; C_a_proved = G["C_a_proved"]
PHI0 = G["PHI0"]; SENS = G["SENS_MEASURED_FACTOR"]; KAPPA_S = G["KAPPA_S"]
CSTAR = math.log(1.5)/KAPPA

OUT = {"c_star": CSTAR, "kappa": KAPPA, "r_h": R_H, "kappa_s_record": KAPPA_S}

# 1. the linearised exponent, as the assembly computes it
Cp = RECORD["Cprime_A1_computedV"]; Ca = C_a_proved(math.sin(PHI0))
bs = bootstrap(1e18, CSTAR, Cp, Ca)
OUT["p_over_ML_proved_Linf"] = bs["p_over_ML"]
OUT["pc_proved_Linf"] = bs["feedback_exponent_pc"]
OUT["exp_pc_proved"] = math.exp(bs["feedback_exponent_pc"])
OUT["p_over_ML_formula"] = 1.5 + 1.5*2*KAPPA*(15*math.pi/4)*(9/4)/R_H
OUT["coefficient_of_m_from_LemmaTprime"] = 1.5*2*KAPPA*(15*math.pi/4)*(9/4)/R_H
bsm = bootstrap(1e18, CSTAR, RECORD["Cprime_measured"], RECORD["material_offset_M"], sens=SENS)
OUT["pc_measured_Linf"] = bsm["feedback_exponent_pc"]

# 2. numerical check of the linearisation: integrate the majorant ODE with a tiny initial
#    m0 and NO drive, and compare the growth factor with exp(p c).  If the linearised
#    coefficient is right, m(c)/m0 -> exp(p c) as m0 -> 0 (muJ also grows, so the check is
#    on the m-equation alone with muJ frozen at 0, which is what "linearising (2.1) about
#    m = 0" means).
#    The assembly's formula uses 15 pi/4 = 2 pi (3/2 + 3/8), i.e. it presumes muJ = mu in the
#    linearisation.  The honest linearisation is the 2x2 system (m, muJ^Lambda) about 0:
#      m'   = G m + (3/2)(2 kappa) (2pi/r_h) [ (3/2) mu + (3/8) muJ ]
#      muJ' = 2 [ kappa (2pi/r_h) ((3/2) mu + (3/8) muJ) + G mu ] ,   mu = (9/4) m ,
#      muJ = muJ^Lambda (1 + 2 kappa_s/L) + 2 kappa_s/L  (the constant part is a drive, not a gain).
def lin_matrix(L, sens=1.0, Cprime=Cp):
    Gm = 1.5 + Cprime/L
    al = sens*(2.0*math.pi/R_H)
    s = 1.0 + 2.0*KAPPA_S/L
    A = np.array([[Gm + 1.5*2*KAPPA*al*1.5*2.25, 1.5*2*KAPPA*al*(3.0/8.0)*s],
                  [2*(KAPPA*al*1.5*2.25 + Gm*2.25), 2*KAPPA*al*(3.0/8.0)*s]])
    return A
def lin_growth(L, c, Cprime, sens=1.0, m0=1e-12):
    A = lin_matrix(L, sens, Cprime)
    def rhs(th, y):
        return list(A @ np.array(y))
    sol = solve_ivp(rhs, (0.0, c), [m0, 0.0], rtol=1e-10, atol=1e-30)
    return sol.y[0, -1]/m0
for L in [1e18, 160.0]:
    A = lin_matrix(L)
    ev = max(np.linalg.eigvals(A).real)
    OUT["lin2x2_dominant_eigenvalue_L%g" % L] = float(ev)
    OUT["lin2x2_pc_L%g" % L] = float(ev*CSTAR)
    OUT["assembly_p_over_ML_L%g" % L] = 1.5 + Cp/L + OUT["coefficient_of_m_from_LemmaTprime"]
    OUT["lin2x2_growth_numeric_L%g" % L] = lin_growth(L, CSTAR, Cp)
    OUT["lin2x2_growth_from_eigenvalue_L%g" % L] = float(math.exp(ev*CSTAR))
# the m-equation alone with muJ frozen at 0 (coefficient 3 pi, not 15 pi/4)
OUT["p_over_ML_muJ_frozen_Linf"] = 1.5 + 1.5*2*KAPPA*(3*math.pi)*(9/4)/R_H

# 3. L_* exactly as a2 does it (same f-set, same bisection), proved and measured, C_K = 1
def eps_at(L, col, C_K=1.0):
    return min(assemble(L, CSTAR, col, C_K=C_K, f=f)["eps"] for f in
               [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0])
Lstar = {}
for col in ["proved", "measured"]:
    for target in [0.5, 0.1]:
        lo, hi = 1.0, 1e60
        if eps_at(hi, col) > target:
            Lstar["%s_eps<=%g" % (col, target)] = None
            continue
        for _ in range(40):
            mid = math.sqrt(lo*hi)
            if eps_at(mid, col) <= target:
                hi = mid
            else:
                lo = mid
            if hi/lo < 1.0001:
                break
        Lstar["%s_eps<=%g" % (col, target)] = hi
OUT["L_star_reproduced"] = Lstar

# 4. eps at the brief's L values, c = c*, best f, proved column
OUT["eps_proved_c_star"] = {str(L): eps_at(L, "proved") for L in [160.0, 640.0, 2560.0]}
OUT["eps_measured_c_star"] = {str(L): eps_at(L, "measured") for L in [160.0, 640.0, 2560.0]}

with open("r1_results.json", "w") as fh:
    json.dump(OUT, fh, indent=1, sort_keys=True, default=str)
for k, v in OUT.items():
    print("%-40s %s" % (k, v))
